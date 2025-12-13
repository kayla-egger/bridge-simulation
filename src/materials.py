"""
Materialeigenschaften für die Brückenberechnung
Basierend auf der Maturaarbeit und Eurocode 3
"""

from dataclasses import dataclass


@dataclass
class Steel:
    """Materialeigenschaften für Baustahl S355"""

    # Bezeichnung
    name: str = "S355"

    # Elastizitätsmodul [N/mm² = MPa]
    # Aus Maturaarbeit S. 22: E = 210'000 N/mm²
    E: float = 210_000.0  # MPa

    # Schubmodul [MPa]
    G: float = 81_000.0  # MPa (≈ E / 2.6)

    # Streckgrenze [MPa]
    # S355: f_yk = 355 MPa für t ≤ 40mm
    f_yk: float = 355.0  # MPa

    # Zugfestigkeit [MPa]
    f_uk: float = 510.0  # MPa

    # Dichte [kg/m³]
    # Aus Maturaarbeit S. 18: ρ = 7.9 * 10³ kg/m³
    rho: float = 7_900.0  # kg/m³

    # Wärmeausdehnungskoeffizient [1/K]
    # Aus Maturaarbeit S. 21: α_T = 1.2 * 10⁻⁵ /K
    alpha_T: float = 1.2e-5  # 1/K

    # Teilsicherheitsbeiwert für Querschnittswiderstand
    gamma_M0: float = 1.0  # Eurocode 3

    # Teilsicherheitsbeiwert für Stabilitätsversagen
    gamma_M1: float = 1.0  # Eurocode 3

    @property
    def E_Pa(self) -> float:
        """E-Modul in Pa"""
        return self.E * 1e6

    @property
    def E_kN_m2(self) -> float:
        """E-Modul in kN/m²"""
        return self.E * 1e3

    @property
    def f_yd(self) -> float:
        """Bemessungswert der Streckgrenze [MPa]"""
        return self.f_yk / self.gamma_M0

    @property
    def weight_per_volume(self) -> float:
        """Wichte [kN/m³]"""
        return self.rho * 9.81 / 1000


@dataclass
class Concrete:
    """Materialeigenschaften für Beton (für Pfeiler und Widerlager)"""

    name: str = "C30/37"

    # Elastizitätsmodul [MPa]
    E: float = 33_000.0  # MPa

    # Charakteristische Druckfestigkeit (Zylinder) [MPa]
    f_ck: float = 30.0  # MPa

    # Dichte [kg/m³]
    # Aus Maturaarbeit S. 17: ρ = 2.2 * 10³ kg/m³
    rho: float = 2_200.0  # kg/m³

    # Teilsicherheitsbeiwert
    gamma_C: float = 1.5

    @property
    def f_cd(self) -> float:
        """Bemessungswert der Druckfestigkeit [MPa] mal Beiwert zur Berücksichtigung der Festigkeitsabnahme unter Dauerlast"""
        return self.f_ck / self.gamma_C * 0.85 

    @property
    def weight_per_volume(self) -> float:
        """Wichte [kN/m³]"""
        return self.rho * 9.81 / 1000


# Standard-Materialien
STEEL_S355 = Steel()
CONCRETE_C30 = Concrete()


def print_material_summary():
    """Gibt eine Zusammenfassung der Materialeigenschaften aus"""
    steel = STEEL_S355
    concrete = CONCRETE_C30

    print("=" * 60)
    print("MATERIALEIGENSCHAFTEN")
    print("=" * 60)

    print(f"\nStahl {steel.name} (Überbau):")
    print(f"  E-Modul:          {steel.E:,.0f} MPa")
    print(f"  Schubmodul G:     {steel.G:,.0f} MPa")
    print(f"  Streckgrenze f_yk:{steel.f_yk:.0f} MPa")
    print(f"  Bemessungswert f_yd: {steel.f_yd:.0f} MPa")
    print(f"  Dichte:           {steel.rho:,.0f} kg/m³")
    print(f"  Wichte:           {steel.weight_per_volume:.2f} kN/m³")
    print(f"  α_T:              {steel.alpha_T:.2e} 1/K")

    print(f"\nBeton {concrete.name} (Pfeiler/Widerlager):")
    print(f"  E-Modul:          {concrete.E:,.0f} MPa")
    print(f"  f_ck:             {concrete.f_ck:.0f} MPa")
    print(f"  f_cd:             {concrete.f_cd:.1f} MPa")
    print(f"  Dichte:           {concrete.rho:,.0f} kg/m³")
    print(f"  Wichte:           {concrete.weight_per_volume:.2f} kN/m³")
    print("=" * 60)


if __name__ == "__main__":
    print_material_summary()
