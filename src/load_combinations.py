"""
Lastfallkombinationen nach Eurocode
Basierend auf Kapitel 3.5.1 der Maturaarbeit

Grenzzustand der Tragfähigkeit (ULS):
E_d = Σ γ_Gj * G_kj + γ_Q1 * Q_k1 + Σ γ_Qi * Ψ_0i * Q_ki
"""

from dataclasses import dataclass, field
from typing import Dict, List
from enum import Enum


class LoadType(Enum):
    """Lastarten"""
    DEAD = "G"           # Eigengewicht
    SURFACING = "G_s"    # Belag
    TRAFFIC = "Q_T"      # Verkehr
    WIND = "Q_W"         # Wind
    TEMPERATURE = "Q_Temp"  # Temperatur


@dataclass
class PartialSafetyFactors:
    """Teilsicherheitsbeiwerte γ nach Eurocode

    Aus Maturaarbeit S. 24
    """
    # Ständige Einwirkungen
    gamma_G_sup: float = 1.35  # ungünstig
    gamma_G_inf: float = 1.00  # günstig

    # Veränderliche Einwirkungen
    gamma_Q_sup: float = 1.50  # ungünstig
    gamma_Q_inf: float = 0.00  # günstig


@dataclass
class CombinationFactors:
    """Kombinationsbeiwerte Ψ nach Eurocode

    Aus Maturaarbeit S. 24
    """
    # Verkehr auf Fussgängerbrücken
    psi_0_traffic: float = 0.40
    psi_1_traffic: float = 0.40
    psi_2_traffic: float = 0.00

    # Wind
    psi_0_wind: float = 0.30  # Aus Maturaarbeit S. 24
    psi_1_wind: float = 0.20
    psi_2_wind: float = 0.00

    # Temperatur
    psi_0_temp: float = 0.80  # Aus Maturaarbeit S. 24
    psi_1_temp: float = 0.60
    psi_2_temp: float = 0.50


@dataclass
class LoadCombination:
    """Eine Lastfallkombination"""
    name: str
    description: str
    factors: Dict[LoadType, float] = field(default_factory=dict)


class LoadCombinations:
    """Lastfallkombinationen für den Grenzzustand der Tragfähigkeit"""

    def __init__(self):
        self.gamma = PartialSafetyFactors()
        self.psi = CombinationFactors()
        self.combinations = self._create_combinations()

    def _create_combinations(self) -> List[LoadCombination]:
        """Erstellt die massgebenden Lastfallkombinationen"""
        combinations = []

        # =====================================================================
        # LC1: Eigengewicht + Verkehr (Leiteinwirkung)
        # =====================================================================
        lc1 = LoadCombination(
            name="LC1",
            description="Eigengewicht + Verkehr (Leiteinwirkung)",
            factors={
                LoadType.DEAD: self.gamma.gamma_G_sup,
                LoadType.SURFACING: self.gamma.gamma_G_sup,
                LoadType.TRAFFIC: self.gamma.gamma_Q_sup,
                LoadType.WIND: self.gamma.gamma_Q_sup * self.psi.psi_0_wind,
                LoadType.TEMPERATURE: self.gamma.gamma_Q_sup * self.psi.psi_0_temp,
            }
        )
        combinations.append(lc1)

        # =====================================================================
        # LC2: Eigengewicht + Wind (Leiteinwirkung)
        # =====================================================================
        lc2 = LoadCombination(
            name="LC2",
            description="Eigengewicht + Wind (Leiteinwirkung)",
            factors={
                LoadType.DEAD: self.gamma.gamma_G_sup,
                LoadType.SURFACING: self.gamma.gamma_G_sup,
                LoadType.TRAFFIC: self.gamma.gamma_Q_sup * self.psi.psi_0_traffic,
                LoadType.WIND: self.gamma.gamma_Q_sup,
                LoadType.TEMPERATURE: self.gamma.gamma_Q_sup * self.psi.psi_0_temp,
            }
        )
        combinations.append(lc2)

        # =====================================================================
        # LC3: Eigengewicht + Temperatur (Leiteinwirkung)
        # =====================================================================
        lc3 = LoadCombination(
            name="LC3",
            description="Eigengewicht + Temperatur (Leiteinwirkung)",
            factors={
                LoadType.DEAD: self.gamma.gamma_G_sup,
                LoadType.SURFACING: self.gamma.gamma_G_sup,
                LoadType.TRAFFIC: self.gamma.gamma_Q_sup * self.psi.psi_0_traffic,
                LoadType.WIND: self.gamma.gamma_Q_sup * self.psi.psi_0_wind,
                LoadType.TEMPERATURE: self.gamma.gamma_Q_sup,
            }
        )
        combinations.append(lc3)

        # =====================================================================
        # LC4: Nur Eigengewicht (für minimale Einwirkung)
        # =====================================================================
        lc4 = LoadCombination(
            name="LC4",
            description="Nur Eigengewicht (minimal)",
            factors={
                LoadType.DEAD: self.gamma.gamma_G_inf,
                LoadType.SURFACING: self.gamma.gamma_G_inf,
                LoadType.TRAFFIC: 0.0,
                LoadType.WIND: 0.0,
                LoadType.TEMPERATURE: 0.0,
            }
        )
        combinations.append(lc4)

        # =====================================================================
        # LC5: Maximale vertikale Last (konservativ)
        # =====================================================================
        lc5 = LoadCombination(
            name="LC5",
            description="Maximale vertikale Last (konservativ)",
            factors={
                LoadType.DEAD: self.gamma.gamma_G_sup,
                LoadType.SURFACING: self.gamma.gamma_G_sup,
                LoadType.TRAFFIC: self.gamma.gamma_Q_sup,
                LoadType.WIND: 0.0,  # Wind reduziert vertikale Last nicht
                LoadType.TEMPERATURE: 0.0,
            }
        )
        combinations.append(lc5)

        return combinations

    def get_combination(self, name: str) -> LoadCombination:
        """Gibt eine spezifische Lastfallkombination zurück"""
        for lc in self.combinations:
            if lc.name == name:
                return lc
        raise ValueError(f"Lastfallkombination '{name}' nicht gefunden")

    def get_design_factors(self, combination_name: str = "LC1") -> Dict[str, float]:
        """Gibt die Bemessungsfaktoren für eine Kombination zurück"""
        lc = self.get_combination(combination_name)
        return {
            "gamma_G": lc.factors[LoadType.DEAD],
            "gamma_Q_traffic": lc.factors[LoadType.TRAFFIC],
            "gamma_Q_wind": lc.factors[LoadType.WIND],
            "gamma_Q_temp": lc.factors[LoadType.TEMPERATURE],
        }


def print_combinations_summary():
    """Gibt eine Zusammenfassung der Lastfallkombinationen aus"""
    lc = LoadCombinations()

    print("=" * 70)
    print("LASTFALLKOMBINATIONEN - Grenzzustand der Tragfähigkeit (ULS)")
    print("=" * 70)

    print("\nTeilsicherheitsbeiwerte γ:")
    print(f"  γ_G,sup = {lc.gamma.gamma_G_sup:.2f} (ständig, ungünstig)")
    print(f"  γ_G,inf = {lc.gamma.gamma_G_inf:.2f} (ständig, günstig)")
    print(f"  γ_Q,sup = {lc.gamma.gamma_Q_sup:.2f} (veränderlich, ungünstig)")

    print("\nKombinationsbeiwerte Ψ_0:")
    print(f"  Ψ_0,Verkehr = {lc.psi.psi_0_traffic:.2f}")
    print(f"  Ψ_0,Wind    = {lc.psi.psi_0_wind:.2f}")
    print(f"  Ψ_0,Temp    = {lc.psi.psi_0_temp:.2f}")

    print("\n" + "-" * 70)
    print("Lastfallkombinationen:")
    print("-" * 70)

    for combo in lc.combinations:
        print(f"\n{combo.name}: {combo.description}")
        print(f"  E_d = ", end="")
        terms = []
        for load_type, factor in combo.factors.items():
            if factor > 0:
                terms.append(f"{factor:.2f}·{load_type.value}")
        print(" + ".join(terms))

    print("\n" + "=" * 70)
    print("Formel aus Maturaarbeit S. 23-24:")
    print("E_d = Σ γ_Gj·G_kj + γ_Q1·Q_k1 + Σ γ_Qi·Ψ_0i·Q_ki")
    print("=" * 70)


if __name__ == "__main__":
    print_combinations_summary()
