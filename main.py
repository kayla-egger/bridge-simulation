#!/usr/bin/env python3
"""
Brueckensimulation - Grenzzustand der Tragfaehigkeit
====================================================

Berechnung von N_Ed, M_Ed, V_Ed fuer die Velobruecke
aus der Maturaarbeit "Bessere Veloverbindung zweier Quartiere"

Verwendung:
    python main.py [--plot] [--save-plots]

Autor: Basierend auf der Maturaarbeit von Kayla Egger
"""

import argparse
import sys
from pathlib import Path

# Projektverzeichnis zum Pfad hinzufügen
sys.path.insert(0, str(Path(__file__).parent))

from src.bridge_geometry import BridgeGeometry, print_geometry_summary
from src.materials import print_material_summary, STEEL_S355
from src.loads import Loads, print_loads_summary
from src.load_combinations import print_combinations_summary
from src.analysis import BridgeAnalysis
from src.results import ResultsEvaluator, plot_results, plot_utilization


def main():
    """Hauptprogramm"""
    parser = argparse.ArgumentParser(
        description='Brückensimulation - Grenzzustand der Tragfähigkeit'
    )
    parser.add_argument(
        '--plot', '-p',
        action='store_true',
        help='Zeige Diagramme der Schnittgrössen'
    )
    parser.add_argument(
        '--save-plots', '-s',
        action='store_true',
        help='Speichere Diagramme als PNG'
    )
    parser.add_argument(
        '--combo', '-c',
        default='LC1',
        help='Lastfallkombination für Analyse (Standard: LC1)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Ausführliche Ausgabe'
    )
    args = parser.parse_args()

    print("\n" + "=" * 70)
    print("  BRUECKENSIMULATION - GRENZZUSTAND DER TRAGFAEHIGKEIT")
    print("  Velobruecke Zuerich (Maturaarbeit)")
    print("=" * 70)

    # =========================================================================
    # 1. Geometrie definieren
    # =========================================================================
    print("\n[1/5] Geometrie initialisieren...")

    # Feldweiten gemäss Maturaarbeit S. 11: [21, 30, 40, 33, 30, 21] = 175m
    geometry = BridgeGeometry()

    if args.verbose:
        print_geometry_summary(geometry)

    # =========================================================================
    # 2. Materialien
    # =========================================================================
    print("[2/5] Materialien laden...")

    if args.verbose:
        print_material_summary()

    # =========================================================================
    # 3. Lasten
    # =========================================================================
    print("[3/5] Lasten berechnen...")

    loads = Loads(geometry)

    if args.verbose:
        print_loads_summary(loads)
        print_combinations_summary()

    # =========================================================================
    # 4. FEM-Analyse mit PyNite
    # =========================================================================
    print("[4/5] FEM-Modell erstellen und analysieren...")

    try:
        analysis = BridgeAnalysis(geometry, nodes_per_span=20)
        analysis.build_model()
        analysis.create_load_combinations()
        analysis.analyze()
        print("      Analyse erfolgreich abgeschlossen.")
    except Exception as e:
        print(f"\n FEHLER bei der Analyse: {e}")
        print("\n Stellen Sie sicher, dass PyNite installiert ist:")
        print("   pip install PyNiteFEA")
        return 1

    # =========================================================================
    # 5. Ergebnisauswertung
    # =========================================================================
    print("[5/5] Ergebnisse auswerten...")

    evaluator = ResultsEvaluator(analysis)

    # Nachweise für gewählte Kombination
    print("\n")
    checks = evaluator.print_summary(args.combo)

    # Kritische Kombination finden
    critical_combo, critical_check, max_util = evaluator.get_critical_combo()

    print("\n" + "-" * 70)
    print("MASSGEBENDE LASTFALLKOMBINATION:")
    print(f"  {critical_combo} - {critical_check}")
    print(f"  Maximaler Ausnutzungsgrad: eta = {max_util*100:.1f}%")
    print("-" * 70)

    # Auflagerkraefte
    print("\nAUFLAGERKRAEFTE:")
    reactions = analysis.get_support_reactions(args.combo)
    for name, data in reactions.items():
        print(f"  {name} (x = {data['position']:.1f} m):")
        print(f"    Ry = {data['Ry']:.2f} kN (vertikal)")

    # =========================================================================
    # Diagramme
    # =========================================================================
    if args.plot or args.save_plots:
        print("\nErstelle Diagramme...")

        save_path_forces = "schnittgroessen.png" if args.save_plots else None
        save_path_util = "ausnutzung.png" if args.save_plots else None

        try:
            plot_results(analysis, args.combo, save_path_forces)
            plot_utilization(evaluator, save_path_util)
        except Exception as e:
            print(f"Fehler beim Erstellen der Diagramme: {e}")
            print("Stellen Sie sicher, dass matplotlib installiert ist:")
            print("  pip install matplotlib")

    # =========================================================================
    # Zusammenfassung für Maturaarbeit
    # =========================================================================
    print("\n" + "=" * 70)
    print("ZUSAMMENFASSUNG FÜR MATURAARBEIT (Kapitel 3.5.1)")
    print("=" * 70)

    print("\nBemessungswerte der Einwirkungen E_d (aus Simulation):")
    max_forces = analysis.get_max_forces(args.combo)
    print(f"  N_Ed = {abs(max_forces['N_Ed'][0]):.2f} kN")
    print(f"  M_Ed = {abs(max_forces['M_Ed'][0]):.2f} kNm")
    print(f"  V_Ed = {abs(max_forces['V_Ed'][0]):.2f} kN")

    print("\nBemessungswerte der Widerstände R_d (Formeln S. 24):")
    print(f"  N_Rd = {evaluator.resistance.N_Rd:,.0f} kN")
    print(f"  M_Rd = {evaluator.resistance.M_Rd:,.0f} kNm")
    print(f"  V_Rd = {evaluator.resistance.V_Rd:,.0f} kN")

    print("\nNachweise (R_d >= E_d):")
    all_ok = True
    for name, check in checks.items():
        status = "[OK]" if check.is_ok else "[X]"
        print(f"  {status} {name}: eta = {check.utilization*100:.1f}%")
        if not check.is_ok:
            all_ok = False

    print("\n" + "=" * 70)
    if all_ok:
        print("FAZIT: Der Grenzzustand der Tragfaehigkeit ist ERFUELLT.")
        print("       Die Bruecke ist fuer die angenommenen Lasten tragfaehig.")
    else:
        print("FAZIT: Der Grenzzustand der Tragfaehigkeit ist NICHT ERFUELLT.")
        print("       Der Querschnitt muss verstaerkt werden.")
    print("=" * 70 + "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
