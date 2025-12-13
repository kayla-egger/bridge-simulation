"""
Lastdefinitionen für die Brückenberechnung
Basierend auf Kapitel 3 der Maturaarbeit

Einwirkungen:
- Vertikale Lasten: Eigengewicht, Verkehr
- Horizontale Lasten: Wind, Temperatur
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np

from .bridge_geometry import BridgeGeometry, CrossSection
from .materials import Steel, STEEL_S355


@dataclass
class Loads:
    """Lastdefinitionen für die Brücke"""

    # Geometrie-Referenz
    geometry: BridgeGeometry

    # Material-Referenz
    steel: Steel = None

    def __post_init__(self):
        if self.steel is None:
            self.steel = STEEL_S355

    # =========================================================================
    # VERTIKALE LASTEN (Kapitel 3.1)
    # =========================================================================

    @property
    def dead_load_area(self) -> float:
        """Eigengewicht des Überbaus [kN/m²]

        Aus Maturaarbeit S. 18:
        p = F_G / (b1 * l) = 7.6 kN/m²
        """
        # Berechnung aus Querschnittsfläche und Materialdichte
        A = self.geometry.cross_section.area  # m²
        rho = self.steel.rho  # kg/m³
        g = 9.81  # m/s²
        b = self.geometry.cross_section.b1  # m

        # Gewichtskraft pro Meter Brückenlänge
        weight_per_m = A * rho * g / 1000  # kN/m

        # Umrechnung in kN/m²
        return weight_per_m / b

    @property
    def dead_load_line(self) -> float:
        """Eigengewicht als Streckenlast [kN/m]"""
        return self.dead_load_area * self.geometry.cross_section.b1

    @property
    def surfacing_load(self) -> float:
        """Belagsgewicht [kN/m²]

        Aus Maturaarbeit S. 17: 0.25 kN/m²
        """
        return 0.25

    @property
    def surfacing_load_line(self) -> float:
        """Belagsgewicht als Streckenlast [kN/m]"""
        return self.surfacing_load * self.geometry.cross_section.b1

    @property
    def total_permanent_load_line(self) -> float:
        """Gesamte ständige Last als Streckenlast [kN/m]"""
        return self.dead_load_line + self.surfacing_load_line

    # =========================================================================
    # VERKEHRSLASTEN (Kapitel 3.1.1)
    # =========================================================================

    def traffic_load_area(self, span_length: Optional[float] = None) -> float:
        """Verkehrslast [kN/m²]

        Aus Maturaarbeit S. 15-16:
        - Grundwert: q_fk = 5.0 kN/m²
        - Bei Spannweiten > 10m: Formel (3.1)

        q_fk = 2.0 + 120/(L+30) mit 2.5 ≤ q_fk ≤ 5.0 kN/m²
        """
        if span_length is None:
            # Maximale Feldweite verwenden
            span_length = max(self.geometry.span_lengths)

        if span_length <= 10:
            q_fk = 5.0
        else:
            q_fk = 2.0 + 120 / (span_length + 30)
            q_fk = max(2.5, min(5.0, q_fk))

        return q_fk

    @property
    def traffic_load_line(self) -> float:
        """Verkehrslast als Streckenlast [kN/m]

        Konservativ: Maximale Feldweite für Berechnung
        """
        q = self.traffic_load_area()
        return q * self.geometry.cross_section.b1

    @property
    def point_load(self) -> float:
        """Einzellast [kN]

        Aus Maturaarbeit S. 16: Q_fwk = 10 kN
        """
        return 10.0

    # =========================================================================
    # WINDLASTEN (Kapitel 3.2.1)
    # =========================================================================

    @property
    def wind_pressure(self) -> float:
        """Windstaudruck q_p(z_e) [kN/m²]

        Aus Maturaarbeit S. 20: q(z) = 0.86 für z = 20m
        """
        return 0.86 / 1000  # Umrechnung N/m² -> kN/m²

    def wind_load_deck_with_traffic(self) -> float:
        """Windlast auf Überbau mit Verkehr [kN/m²]

        Aus Maturaarbeit S. 21:
        W = q(z) * c_fx,o * ψ_3D = 0.86 * 1.3 * 0.70 = 0.78 N/m²
        """
        return 0.86 * 1.3 * 0.70 / 1000  # kN/m²

    def wind_load_deck_without_traffic(self) -> float:
        """Windlast auf Überbau ohne Verkehr [kN/m²]

        Aus Maturaarbeit S. 21:
        W = q(z) * c_fx,o * ψ_3D = 0.86 * 1.3 * 0.85 = 0.95 N/m²
        """
        return 0.86 * 1.3 * 0.85 / 1000  # kN/m²

    @property
    def wind_reference_area(self) -> float:
        """Windangriffsfläche A_ref,x [m²]

        Aus Maturaarbeit S. 20:
        A_ref,x = (d + 1.2) * L
        d = Höhe des Überbaus = 1.5m
        """
        d = self.geometry.cross_section.h2
        L = self.geometry.total_length
        return (d + 1.2) * L

    @property
    def wind_force_with_traffic(self) -> float:
        """Gesamte Windkraft mit Verkehr [kN]"""
        return self.wind_load_deck_with_traffic() * self.wind_reference_area

    @property
    def wind_load_line(self) -> float:
        """Windlast als horizontale Streckenlast [kN/m]

        Bezogen auf die Überbauhöhe
        """
        d = self.geometry.cross_section.h2
        return self.wind_load_deck_with_traffic() * (d + 1.2)

    # =========================================================================
    # TEMPERATUR (Kapitel 3.2.2)
    # =========================================================================

    @property
    def temperature_range(self) -> tuple:
        """Temperaturbereich [°C]

        Aus Maturaarbeit S. 21:
        ΔT_e,min = -26°C
        ΔT_e,max = +51°C
        """
        return (-26.0, 51.0)

    @property
    def temperature_change(self) -> float:
        """Temperaturschwankung ΔT_N [K]

        Aus Maturaarbeit S. 21:
        ΔT_N = ΔT_e,max - ΔT_e,min = 77K
        """
        T_min, T_max = self.temperature_range
        return T_max - T_min

    @property
    def thermal_expansion(self) -> float:
        """Thermische Längenänderung Δl [m]

        Aus Maturaarbeit S. 21:
        Δl = α_T * L * ΔT_N
        """
        alpha_T = self.steel.alpha_T
        L = self.geometry.total_length
        delta_T = self.temperature_change
        return alpha_T * L * delta_T


def print_loads_summary(loads: Loads):
    """Gibt eine Zusammenfassung der Lasten aus"""
    print("=" * 60)
    print("LASTDEFINITIONEN")
    print("=" * 60)

    print("\n--- Ständige Lasten (Kapitel 3.1.2) ---")
    print(f"  Eigengewicht:     {loads.dead_load_area:.2f} kN/m²")
    print(f"                    {loads.dead_load_line:.2f} kN/m (Streckenlast)")
    print(f"  Belag:            {loads.surfacing_load:.2f} kN/m²")
    print(f"  Gesamt ständig:   {loads.total_permanent_load_line:.2f} kN/m")

    print("\n--- Verkehrslasten (Kapitel 3.1.1) ---")
    for span in loads.geometry.span_lengths:
        q = loads.traffic_load_area(span)
        print(f"  Feld L={span:.0f}m:     {q:.2f} kN/m²")
    print(f"  Streckenlast:     {loads.traffic_load_line:.2f} kN/m")
    print(f"  Einzellast:       {loads.point_load:.1f} kN")

    print("\n--- Windlasten (Kapitel 3.2.1) ---")
    print(f"  Mit Verkehr:      {loads.wind_load_deck_with_traffic()*1000:.3f} N/m²")
    print(f"  Ohne Verkehr:     {loads.wind_load_deck_without_traffic()*1000:.3f} N/m²")
    print(f"  Angriffsfläche:   {loads.wind_reference_area:.1f} m²")
    print(f"  Windkraft:        {loads.wind_force_with_traffic:.1f} kN")

    print("\n--- Temperatur (Kapitel 3.2.2) ---")
    T_min, T_max = loads.temperature_range
    print(f"  T_min:            {T_min:.0f}°C")
    print(f"  T_max:            {T_max:.0f}°C")
    print(f"  ΔT_N:             {loads.temperature_change:.0f} K")
    print(f"  Δl:               {loads.thermal_expansion*1000:.1f} mm")
    print("=" * 60)


if __name__ == "__main__":
    from .bridge_geometry import BridgeGeometry

    geometry = BridgeGeometry()
    loads = Loads(geometry)
    print_loads_summary(loads)
