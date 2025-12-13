"""
Ergebnisauswertung und Visualisierung
Nachweis des Grenzzustands der Tragfähigkeit

R_d ≥ E_d

Formeln aus Maturaarbeit S. 24:
- N_Rd = A * f_yk / γ_M0
- M_Rd = W_el * f_yk / γ_M0
- V_Rd = A_V * f_yk / (√3 * γ_M0)
"""

from dataclasses import dataclass
from typing import Dict, Tuple, Optional
import numpy as np

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

from .bridge_geometry import BridgeGeometry, CrossSection
from .materials import Steel, STEEL_S355
from .analysis import BridgeAnalysis


@dataclass
class ResistanceValues:
    """Bemessungswerte der Widerstände R_d

    Formeln aus Maturaarbeit S. 24 (Ref. [3, S.81, S.70])
    """
    # Querschnitt
    cross_section: CrossSection

    # Material
    steel: Steel

    @property
    def N_Rd(self) -> float:
        """Bemessungswert der Normalkrafttragfähigkeit [kN]

        N_Rd = A * f_yk / γ_M0
        """
        A = self.cross_section.area  # m²
        f_yk = self.steel.f_yk * 1000  # MPa -> kN/m²
        gamma_M0 = self.steel.gamma_M0

        return A * f_yk / gamma_M0

    @property
    def M_Rd(self) -> float:
        """Bemessungswert der Biegetragfähigkeit [kNm]

        M_Rd = W_el * f_yk / γ_M0
        """
        W_el = self.cross_section.elastic_section_modulus  # m³
        f_yk = self.steel.f_yk * 1000  # MPa -> kN/m²
        gamma_M0 = self.steel.gamma_M0

        return W_el * f_yk / gamma_M0

    @property
    def V_Rd(self) -> float:
        """Bemessungswert der Querkrafttragfähigkeit [kN]

        V_Rd = A_V * f_yk / (√3 * γ_M0)
        """
        A_V = self.cross_section.shear_area  # m²
        f_yk = self.steel.f_yk * 1000  # MPa -> kN/m²
        gamma_M0 = self.steel.gamma_M0

        return A_V * f_yk / (np.sqrt(3) * gamma_M0)


@dataclass
class DesignCheck:
    """Ergebnis eines Tragfähigkeitsnachweises"""
    name: str
    E_d: float          # Einwirkung
    R_d: float          # Widerstand
    position: float     # Position [m]
    unit: str

    @property
    def utilization(self) -> float:
        """Ausnutzungsgrad η = E_d / R_d"""
        return abs(self.E_d) / self.R_d

    @property
    def is_ok(self) -> bool:
        """Nachweis erfuellt wenn E_d <= R_d"""
        return abs(self.E_d) <= self.R_d

    def __str__(self) -> str:
        status = "OK" if self.is_ok else "NICHT ERFUELLT"
        return (
            f"{self.name}:\n"
            f"  E_d = {abs(self.E_d):.2f} {self.unit} bei x = {self.position:.1f} m\n"
            f"  R_d = {self.R_d:.2f} {self.unit}\n"
            f"  eta = {self.utilization*100:.1f}% -> {status}"
        )


class ResultsEvaluator:
    """Auswertung der Analyseergebnisse"""

    def __init__(self, analysis: BridgeAnalysis):
        self.analysis = analysis
        self.geometry = analysis.geometry
        self.steel = analysis.steel

        # Bemessungswiderstände
        self.resistance = ResistanceValues(
            cross_section=self.geometry.cross_section,
            steel=self.steel
        )

    def perform_design_checks(
        self,
        combo_name: str = "LC1"
    ) -> Dict[str, DesignCheck]:
        """Führt alle Tragfähigkeitsnachweise durch"""
        max_forces = self.analysis.get_max_forces(combo_name)

        checks = {}

        # Normalkraftnachweis
        N_Ed, x_N = max_forces['N_Ed']
        checks['Normalkraft'] = DesignCheck(
            name="Normalkraft N_Ed <= N_Rd",
            E_d=N_Ed,
            R_d=self.resistance.N_Rd,
            position=x_N,
            unit="kN"
        )

        # Biegemomentnachweis
        M_Ed, x_M = max_forces['M_Ed']
        checks['Biegemoment'] = DesignCheck(
            name="Biegemoment M_Ed <= M_Rd",
            E_d=M_Ed,
            R_d=self.resistance.M_Rd,
            position=x_M,
            unit="kNm"
        )

        # Querkraftnachweis
        V_Ed, x_V = max_forces['V_Ed']
        checks['Querkraft'] = DesignCheck(
            name="Querkraft V_Ed <= V_Rd",
            E_d=V_Ed,
            R_d=self.resistance.V_Rd,
            position=x_V,
            unit="kN"
        )

        return checks

    def print_summary(self, combo_name: str = "LC1"):
        """Gibt eine Zusammenfassung der Ergebnisse aus"""
        print("=" * 70)
        print(f"GRENZZUSTAND DER TRAGFAEHIGKEIT - {combo_name}")
        print("=" * 70)

        print("\n--- Bemessungswidersstaende R_d (Maturaarbeit S. 24) ---")
        print(f"  N_Rd = {self.resistance.N_Rd:,.0f} kN")
        print(f"  M_Rd = {self.resistance.M_Rd:,.0f} kNm")
        print(f"  V_Rd = {self.resistance.V_Rd:,.0f} kN")

        print("\n--- Nachweise ---")
        checks = self.perform_design_checks(combo_name)

        all_ok = True
        for name, check in checks.items():
            print(f"\n{check}")
            if not check.is_ok:
                all_ok = False

        print("\n" + "=" * 70)
        if all_ok:
            print("ERGEBNIS: Alle Nachweise erfuellt!")
        else:
            print("ERGEBNIS: Mindestens ein Nachweis NICHT erfuellt!")
        print("=" * 70)

        return checks

    def get_all_combo_results(self) -> Dict[str, Dict[str, DesignCheck]]:
        """Führt Nachweise für alle Lastfallkombinationen durch"""
        results = {}
        for combo in self.analysis.combinations.combinations:
            results[combo.name] = self.perform_design_checks(combo.name)
        return results

    def get_critical_combo(self) -> Tuple[str, str, float]:
        """Findet die massgebende Lastfallkombination

        Returns:
            (combo_name, check_name, max_utilization)
        """
        all_results = self.get_all_combo_results()

        max_util = 0
        critical_combo = None
        critical_check = None

        for combo_name, checks in all_results.items():
            for check_name, check in checks.items():
                if check.utilization > max_util:
                    max_util = check.utilization
                    critical_combo = combo_name
                    critical_check = check_name

        return critical_combo, critical_check, max_util


def plot_results(
    analysis: BridgeAnalysis,
    combo_name: str = "LC1",
    save_path: Optional[str] = None
):
    """Erstellt Diagramme der Schnittgrössen"""
    if not MATPLOTLIB_AVAILABLE:
        print("Matplotlib nicht installiert. Visualisierung nicht möglich.")
        return

    forces = analysis.get_member_forces(combo_name)
    support_positions = analysis.geometry.support_positions

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

    # Normalkraft N
    ax1 = axes[0]
    ax1.plot(forces['x'], forces['N'], 'b-', linewidth=1.5)
    ax1.fill_between(forces['x'], forces['N'], alpha=0.3)
    ax1.set_ylabel('Normalkraft N [kN]')
    ax1.set_title(f'Schnittgrössen - Lastfallkombination {combo_name}')
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='k', linewidth=0.5)

    # Auflagerpositionen markieren
    for x in support_positions:
        ax1.axvline(x=x, color='gray', linestyle='--', alpha=0.5)

    # Querkraft V
    ax2 = axes[1]
    ax2.plot(forces['x'], forces['V'], 'g-', linewidth=1.5)
    ax2.fill_between(forces['x'], forces['V'], alpha=0.3, color='green')
    ax2.set_ylabel('Querkraft V [kN]')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='k', linewidth=0.5)

    for x in support_positions:
        ax2.axvline(x=x, color='gray', linestyle='--', alpha=0.5)

    # Biegemoment M
    ax3 = axes[2]
    ax3.plot(forces['x'], forces['M'], 'r-', linewidth=1.5)
    ax3.fill_between(forces['x'], forces['M'], alpha=0.3, color='red')
    ax3.set_ylabel('Biegemoment M [kNm]')
    ax3.set_xlabel('Position x [m]')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=0, color='k', linewidth=0.5)

    for x in support_positions:
        ax3.axvline(x=x, color='gray', linestyle='--', alpha=0.5)

    # Biegemoment invertiert darstellen (Konvention: Zug unten)
    ax3.invert_yaxis()

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Diagramm gespeichert: {save_path}")

    plt.show()


def plot_utilization(
    evaluator: ResultsEvaluator,
    save_path: Optional[str] = None
):
    """Erstellt ein Balkendiagramm der Ausnutzungsgrade"""
    if not MATPLOTLIB_AVAILABLE:
        print("Matplotlib nicht installiert. Visualisierung nicht möglich.")
        return

    all_results = evaluator.get_all_combo_results()

    fig, ax = plt.subplots(figsize=(12, 6))

    combo_names = list(all_results.keys())
    check_names = ['Normalkraft', 'Biegemoment', 'Querkraft']
    x = np.arange(len(combo_names))
    width = 0.25

    colors = ['#1f77b4', '#d62728', '#2ca02c']

    for i, check_name in enumerate(check_names):
        utilizations = [
            all_results[combo][check_name].utilization * 100
            for combo in combo_names
        ]
        bars = ax.bar(x + i * width, utilizations, width,
                      label=check_name, color=colors[i], alpha=0.8)

    # 100% Linie (Grenze)
    ax.axhline(y=100, color='red', linestyle='--', linewidth=2,
               label='Grenze (100%)')

    ax.set_xlabel('Lastfallkombination')
    ax.set_ylabel('Ausnutzungsgrad η [%]')
    ax.set_title('Ausnutzungsgrade für alle Lastfallkombinationen')
    ax.set_xticks(x + width)
    ax.set_xticklabels(combo_names)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Diagramm gespeichert: {save_path}")

    plt.show()


if __name__ == "__main__":
    from .analysis import run_analysis

    # Analyse durchführen
    analysis = run_analysis()

    # Ergebnisse auswerten
    evaluator = ResultsEvaluator(analysis)
    evaluator.print_summary("LC1")

    # Kritische Kombination finden
    critical = evaluator.get_critical_combo()
    print(f"\nKritische Kombination: {critical[0]}, {critical[1]}, η = {critical[2]*100:.1f}%")
