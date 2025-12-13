"""
Strukturanalyse mit PyNite
Berechnung der Schnittgrössen N, M, V für den Grenzzustand der Tragfähigkeit
"""

from typing import Dict, List, Tuple, Optional
import numpy as np

from Pynite import FEModel3D

from .bridge_geometry import BridgeGeometry
from .materials import Steel, STEEL_S355
from .loads import Loads
from .load_combinations import LoadCombinations, LoadType


class BridgeAnalysis:
    """Strukturanalyse der Brücke mit PyNite"""

    def __init__(
        self,
        geometry: BridgeGeometry,
        steel: Steel = None,
        nodes_per_span: int = 20
    ):
        self.geometry = geometry
        self.steel = steel or STEEL_S355
        self.nodes_per_span = nodes_per_span
        self.loads = Loads(geometry, steel)
        self.combinations = LoadCombinations()

        # PyNite Modell
        self.model = None
        self.results = {}

    def build_model(self):
        """Erstellt das PyNite FEM-Modell"""
        self.model = FEModel3D()

        # Materialeigenschaften hinzufügen
        # PyNite erwartet: E [ksi oder konsistente Einheiten]
        # Wir verwenden kN und m: E in kN/m²
        E = self.steel.E * 1000  # MPa -> kN/m²
        G = self.steel.G * 1000  # MPa -> kN/m²
        nu = 0.3  # Poissonzahl
        rho = self.steel.rho / 1000  # kg/m³ -> t/m³ (für konsistente Einheiten)

        self.model.add_material("Steel", E, G, nu, rho)

        # Querschnittseigenschaften
        cs = self.geometry.cross_section
        A = cs.area  # m²
        Iy = cs.moment_of_inertia  # m⁴
        Iz = cs.moment_of_inertia  # m⁴ (symmetrisch angenommen)
        J = 2 * Iy  # Torsionsträgheitsmoment (Näherung)

        # Querschnitt hinzufügen
        self.model.add_section("BridgeDeck", A, Iy, Iz, J)

        # Knoten erstellen
        self._create_nodes()

        # Elemente erstellen
        self._create_members()

        # Auflager definieren
        self._create_supports()

    def _create_nodes(self):
        """Erstellt die Knoten des Modells"""
        self.node_positions = []
        node_id = 1

        for i, span in enumerate(self.geometry.span_lengths):
            start_x = self.geometry.support_positions[i]
            end_x = self.geometry.support_positions[i + 1]

            # Knoten im Feld
            x_positions = np.linspace(start_x, end_x, self.nodes_per_span + 1)

            for j, x in enumerate(x_positions):
                # Ersten Knoten nur beim ersten Feld hinzufügen
                if i > 0 and j == 0:
                    continue

                self.model.add_node(f"N{node_id}", x, 0, 0)
                self.node_positions.append((node_id, x))
                node_id += 1

        self.num_nodes = node_id - 1

    def _create_members(self):
        """Erstellt die Balkenelemente"""
        for i in range(1, self.num_nodes):
            member_name = f"M{i}"
            node_i = f"N{i}"
            node_j = f"N{i+1}"

            self.model.add_member(
                member_name,
                node_i,
                node_j,
                "Steel",
                "BridgeDeck"
            )

    def _create_supports(self):
        """Erstellt die Auflager

        Lagersystem aus Maturaarbeit S. 13:
        - Schwimmende Lagerung
        - Widerlager: einseitig fixiert (x-Richtung frei)
        - Mittlere Pfeiler: fixiert oder frei je nach Position
        """
        support_positions = self.geometry.support_positions

        for i, x_pos in enumerate(support_positions):
            # Finde den Knoten an dieser Position
            node_name = None
            for node_id, x in self.node_positions:
                if abs(x - x_pos) < 0.01:
                    node_name = f"N{node_id}"
                    break

            if node_name is None:
                continue

            # Alle Auflager: vertikal (Y) und seitlich (Z) gehalten
            # Rotation um Z erlaubt (Gelenk für Durchlaufträger-Verhalten)

            if i == 0:
                # Erstes Widerlager: fest in x (Festlager)
                self.model.def_support(
                    node_name,
                    support_DX=True,   # x fixiert
                    support_DY=True,   # y fixiert (vertikal)
                    support_DZ=True,   # z fixiert (seitlich)
                    support_RX=False,  # Rotation um x frei
                    support_RY=False,  # Rotation um y frei
                    support_RZ=False   # Rotation um z frei (Gelenk)
                )
            else:
                # Andere Auflager: Loslager (x frei)
                self.model.def_support(
                    node_name,
                    support_DX=False,  # x frei (Verschiebung möglich)
                    support_DY=True,   # y fixiert (vertikal)
                    support_DZ=True,   # z fixiert (seitlich)
                    support_RX=False,
                    support_RY=False,
                    support_RZ=False
                )

    def add_load_case(self, case_name: str):
        """Fügt einen Lastfall hinzu"""
        self.model.add_load_combo(case_name, factors=[(case_name, 1.0)])

    def apply_dead_load(self, case_name: str = "DL"):
        """Wendet das Eigengewicht an"""
        # Streckenlast in -y Richtung (nach unten)
        q = -self.loads.total_permanent_load_line  # kN/m (negativ = nach unten)

        for i in range(1, self.num_nodes):
            member_name = f"M{i}"
            self.model.add_member_dist_load(
                member_name,
                direction="Fy",
                w1=q, w2=q,
                x1=0, x2=None,  # Über gesamte Länge
                case=case_name
            )

    def apply_traffic_load(self, case_name: str = "TL"):
        """Wendet die Verkehrslast an"""
        q = -self.loads.traffic_load_line  # kN/m (negativ = nach unten)

        for i in range(1, self.num_nodes):
            member_name = f"M{i}"
            self.model.add_member_dist_load(
                member_name,
                direction="Fy",
                w1=q, w2=q,
                x1=0, x2=None,
                case=case_name
            )

    def apply_wind_load(self, case_name: str = "WL"):
        """Wendet die Windlast an (horizontal)"""
        q = self.loads.wind_load_line  # kN/m (positiv = in z-Richtung)

        for i in range(1, self.num_nodes):
            member_name = f"M{i}"
            self.model.add_member_dist_load(
                member_name,
                direction="Fz",
                w1=q, w2=q,
                x1=0, x2=None,
                case=case_name
            )

    def create_load_combinations(self):
        """Erstellt die Lastfallkombinationen"""
        # Einzelne Lastfälle
        self.apply_dead_load("DL")
        self.apply_traffic_load("TL")
        self.apply_wind_load("WL")

        # Kombinationen erstellen
        for combo in self.combinations.combinations:
            factors = {}

            # Eigengewicht
            gamma_G = combo.factors.get(LoadType.DEAD, 0)
            if gamma_G > 0:
                factors["DL"] = gamma_G

            # Verkehr
            gamma_T = combo.factors.get(LoadType.TRAFFIC, 0)
            if gamma_T > 0:
                factors["TL"] = gamma_T

            # Wind
            gamma_W = combo.factors.get(LoadType.WIND, 0)
            if gamma_W > 0:
                factors["WL"] = gamma_W

            if factors:
                self.model.add_load_combo(combo.name, factors)

    def analyze(self):
        """Führt die Analyse durch"""
        print("Starte Analyse...")
        self.model.analyze()
        print("Analyse abgeschlossen.")

    def get_member_forces(
        self,
        combo_name: str = "LC1",
        points_per_member: int = 10
    ) -> Dict[str, np.ndarray]:
        """Berechnet die Schnittgrössen entlang der Brücke

        Returns:
            Dict mit:
            - 'x': Positionen [m]
            - 'N': Normalkraft [kN]
            - 'V': Querkraft [kN]
            - 'M': Biegemoment [kNm]
        """
        x_values = []
        N_values = []
        V_values = []
        M_values = []

        for i in range(1, self.num_nodes):
            member_name = f"M{i}"
            member = self.model.members[member_name]
            L = member.L()

            # Position des Elementanfangs
            x_start = self.node_positions[i-1][1]

            # Punkte entlang des Elements
            for j in range(points_per_member + 1):
                xi = j / points_per_member
                x_local = xi * L

                x_global = x_start + x_local

                # Schnittgrössen abrufen
                try:
                    # Querkraft (Fy -> Schub in y-Richtung)
                    V = member.shear("Fy", x_local, combo_name)
                    # Biegemoment (Mz -> um z-Achse)
                    M = member.moment("Mz", x_local, combo_name)
                    # Normalkraft
                    N = member.axial(x_local, combo_name)

                    x_values.append(x_global)
                    N_values.append(N)
                    V_values.append(V)
                    M_values.append(M)
                except Exception:
                    pass

        return {
            'x': np.array(x_values),
            'N': np.array(N_values),
            'V': np.array(V_values),
            'M': np.array(M_values)
        }

    def get_max_forces(self, combo_name: str = "LC1") -> Dict[str, Tuple[float, float]]:
        """Gibt die maximalen Schnittgrössen zurück

        Returns:
            Dict mit (max_value, position) für N, V, M
        """
        forces = self.get_member_forces(combo_name)

        results = {}

        # Maximale Normalkraft (Betrag)
        idx_N = np.argmax(np.abs(forces['N']))
        results['N_Ed'] = (forces['N'][idx_N], forces['x'][idx_N])

        # Maximale Querkraft (Betrag)
        idx_V = np.argmax(np.abs(forces['V']))
        results['V_Ed'] = (forces['V'][idx_V], forces['x'][idx_V])

        # Maximales Biegemoment (Betrag)
        idx_M = np.argmax(np.abs(forces['M']))
        results['M_Ed'] = (forces['M'][idx_M], forces['x'][idx_M])

        return results

    def get_support_reactions(self, combo_name: str = "LC1") -> Dict[str, Dict]:
        """Gibt die Auflagerkräfte zurück"""
        reactions = {}

        for i, x_pos in enumerate(self.geometry.support_positions):
            # Finde den Knoten an dieser Position
            for node_id, x in self.node_positions:
                if abs(x - x_pos) < 0.01:
                    node_name = f"N{node_id}"
                    node = self.model.nodes[node_name]

                    reactions[f"Auflager_{i+1}"] = {
                        'position': x_pos,
                        'Rx': node.RxnFX.get(combo_name, 0),
                        'Ry': node.RxnFY.get(combo_name, 0),
                        'Rz': node.RxnFZ.get(combo_name, 0),
                    }
                    break

        return reactions


def run_analysis(geometry: BridgeGeometry = None) -> BridgeAnalysis:
    """Führt eine vollständige Analyse durch"""
    if geometry is None:
        geometry = BridgeGeometry()

    analysis = BridgeAnalysis(geometry)
    analysis.build_model()
    analysis.create_load_combinations()
    analysis.analyze()

    return analysis


if __name__ == "__main__":
    # Test
    analysis = run_analysis()

    print("\n" + "=" * 60)
    print("ANALYSEERGEBNISSE")
    print("=" * 60)

    for combo in ["LC1", "LC5"]:
        print(f"\n--- {combo} ---")
        max_forces = analysis.get_max_forces(combo)

        for name, (value, pos) in max_forces.items():
            print(f"  {name}: {value:.2f} bei x = {pos:.1f} m")
