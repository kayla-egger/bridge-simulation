"""
Brückengeometrie und Querschnittseigenschaften
Basierend auf der Maturaarbeit "Entwurf einer Brücke"

Zweiter Entwurf: Stahlüberbau mit trapezförmigem Hohlkastenquerschnitt
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class CrossSection:
    """Querschnittseigenschaften des Überbaus (Trapez-Hohlkasten)"""

    # Geometrie aus Maturaarbeit S. 8 (Zweiter Entwurf)
    b1: float = 6.0      # Breite oben [m]
    b2: float = 0.1      # Wandstärke [m]
    b3: float = 3.0      # Breite unten [m]
    h1: float = 0.1      # Deckenstärke [m]
    h2: float = 1.5      # Gesamthöhe [m]

    @property
    def area(self) -> float:
        """Querschnittsfläche A [m²]"""
        # Äusseres Trapez minus inneres Trapez
        A_outer = (self.b1 + self.b3) * self.h2 / 2

        # Innere Abmessungen
        b1_inner = self.b1 - 2 * self.b2
        b3_inner = self.b3 - 2 * self.b2
        h_inner = self.h2 - 2 * self.h1

        A_inner = (b1_inner + b3_inner) * h_inner / 2

        return A_outer - A_inner

    @property
    def moment_of_inertia(self) -> float:
        """Flächenträgheitsmoment I [m⁴]

        Formel für Trapez aus Maturaarbeit S. 22:
        I = h/48 * (b1 + b3) * (b1² + b3²)
        """
        # Vereinfachte Berechnung für Trapez (aus Maturaarbeit)
        I = (self.h2 / 48) * (self.b1 + self.b3) * (self.b1**2 + self.b3**2)
        return I

    @property
    def elastic_section_modulus(self) -> float:
        """Elastisches Widerstandsmoment W_el [m³]

        W_el = I / z_max (Abstand zur äussersten Faser)
        """
        # Schwerpunktabstand zur Oberkante (für Trapez)
        z_s = self.h2 * (self.b3 + 2 * self.b1) / (3 * (self.b1 + self.b3))
        z_max = max(z_s, self.h2 - z_s)

        return self.moment_of_inertia / z_max

    @property
    def shear_area(self) -> float:
        """Schubfläche A_V [m²] - approximiert als Stegfläche"""
        # Zwei geneigte Stege
        # Steghöhe ≈ h2 - 2*h1
        h_steg = self.h2 - 2 * self.h1
        # Steglänge (geneigt)
        delta_b = (self.b1 - self.b3) / 2
        steg_length = np.sqrt(h_steg**2 + delta_b**2)

        return 2 * steg_length * self.b2


@dataclass
class BridgeGeometry:
    """Geometrie der Brücke (Durchlaufträger)"""

    # Gesamtlänge [m]
    total_length: float = 175.0

    # Feldweiten [m] - von links nach rechts
    # 6 Felder mit 5 Pfeilern und 2 Widerlagern
    # Abstände gemäss Maturaarbeit S. 11
    span_lengths: List[float] = None

    # Höhe der Pfeiler [m]
    pier_height: float = 11.5  # h3 + h4 aus S. 11 (ohne Sockel für Berechnung)

    # Querschnitt
    cross_section: CrossSection = None

    def __post_init__(self):
        if self.span_lengths is None:
            # Feldweiten gemäss Maturaarbeit S. 11
            # 21 + 30 + 40 + 33 + 30 + 21 = 175m
            self.span_lengths = [21.0, 30.0, 40.0, 33.0, 30.0, 21.0]

        if self.cross_section is None:
            self.cross_section = CrossSection()

    @property
    def num_spans(self) -> int:
        """Anzahl der Felder"""
        return len(self.span_lengths)

    @property
    def num_supports(self) -> int:
        """Anzahl der Auflager (Widerlager + Pfeiler)"""
        return self.num_spans + 1

    @property
    def support_positions(self) -> List[float]:
        """Positionen der Auflager [m] von links"""
        positions = [0.0]
        cumsum = 0.0
        for span in self.span_lengths:
            cumsum += span
            positions.append(cumsum)
        return positions

    def get_node_positions(self, nodes_per_span: int = 10) -> List[float]:
        """Knotenpositionen für FEM-Modell"""
        positions = []
        for i, span in enumerate(self.span_lengths):
            start = self.support_positions[i]
            end = self.support_positions[i + 1]

            # Knoten im Feld (ohne Endpunkt, ausser beim letzten Feld)
            span_nodes = np.linspace(start, end, nodes_per_span + 1)

            if i == 0:
                positions.extend(span_nodes)
            else:
                positions.extend(span_nodes[1:])  # Ohne Dopplung am Auflager

        return positions


def print_geometry_summary(geometry: BridgeGeometry):
    """Gibt eine Zusammenfassung der Geometrie aus"""
    print("=" * 60)
    print("BRÜCKENGEOMETRIE - Zusammenfassung")
    print("=" * 60)
    print(f"\nGesamtlänge:        {geometry.total_length:.1f} m")
    print(f"Anzahl Felder:      {geometry.num_spans}")
    print(f"Anzahl Auflager:    {geometry.num_supports}")

    print(f"\nFeldweiten [m]:")
    for i, span in enumerate(geometry.span_lengths):
        print(f"  Feld {i+1}: {span:.1f} m")

    print(f"\nAuflagerpositionen [m]: {geometry.support_positions}")

    cs = geometry.cross_section
    print(f"\nQuerschnitt (Stahl-Trapez-Hohlkasten):")
    print(f"  Breite oben b1:   {cs.b1:.2f} m")
    print(f"  Breite unten b3:  {cs.b3:.2f} m")
    print(f"  Höhe h2:          {cs.h2:.2f} m")
    print(f"  Wandstärke:       {cs.b2:.2f} m")
    print(f"\n  Fläche A:         {cs.area:.4f} m²")
    print(f"  Trägheitsmoment I:{cs.moment_of_inertia:.4f} m⁴")
    print(f"  Widerstandsmoment:{cs.elastic_section_modulus:.4f} m³")
    print(f"  Schubfläche A_V:  {cs.shear_area:.4f} m²")
    print("=" * 60)


if __name__ == "__main__":
    # Test
    geometry = BridgeGeometry()
    print_geometry_summary(geometry)
