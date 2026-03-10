from pathlib import Path
import argparse
import opengate as gate
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt

from utils.materials import define_all_materials


# ─────────────────────────────────────────────
# Energy presets avec géométrie adaptée
# ─────────────────────────────────────────────
PRESETS = {
    "min": {
        "energy_mev":       62.4,
        "phantom_depth_cm": 3.5,   # range théorique ~3.8 cm → fantôme 3.5 cm
        "pcb_offset_cm":    0.5,   # gap entre sortie fantôme et PCB
    },
    "max": {
        "energy_mev":       252.7,
        "phantom_depth_cm": 30.0,  # range théorique ~38 cm → fantôme 30 cm
        "pcb_offset_cm":    2.0,
    },
}


# ─────────────────────────────────────────────
# Build simulation
# ─────────────────────────────────────────────
def build_simulation(preset: dict, n_primaries: int, output_dir: Path):

    # Units
    m   = gate.g4_units.m
    cm  = gate.g4_units.cm
    mm  = gate.g4_units.mm
    um  = gate.g4_units.um
    MeV = gate.g4_units.MeV

    energy_mev       = preset["energy_mev"]
    phantom_depth_cm = preset["phantom_depth_cm"]
    pcb_offset_cm    = preset["pcb_offset_cm"]

    # ── Simulation ────────────────────────────────────────────────────
    sim = gate.Simulation()
    sim.output_dir        = output_dir
    sim.random_seed       = "auto"
    sim.number_of_threads = 1
    sim.g4_verbose        = False
    sim.visu              = False
    sim.progress_bar      = True

    sim.physics_manager.physics_list_name = "QGSP_BIC_EMZ"
    sim.physics_manager.set_production_cut("world", "all", 1 * mm)

    # ── World ─────────────────────────────────────────────────────────
    world          = sim.world
    world.size     = [1 * m, 1 * m, 1 * m]
    world.material = "G4_AIR"

    # ── Phantom ───────────────────────────────────────────────────────
    # Centré en z=0 → face arrière à +phantom_depth_cm/2
    phantom             = sim.add_volume("Box", "phantom")
    phantom.mother      = world.name
    phantom.size        = [5 * cm, 5 * cm, phantom_depth_cm * cm]
    phantom.translation = [0, 0, 0]
    phantom.material    = "G4_WATER"
    phantom.color       = [0, 1, 1, 0.25]

    # Tous les calculs de position en cm (nombres purs)
    # puis on multiplie par cm uniquement au moment de .translation / .size
    PCB_THICKNESS_CM     = 0.184    # 1.84 mm
    PLANE_IN_THICKNESS_CM = 0.0005  # 5 µm  ← espace réservé explicitement
    SENSOR_THICKNESS_CM  = 0.05     # 0.5 mm

    # Chaîne stricte : phantom | gap | PCB | plane_in | sensor
    phantom_back_z  = phantom_depth_cm / 2
    pcb_center_z    = phantom_back_z + pcb_offset_cm + PCB_THICKNESS_CM / 2
    pcb_back_z      = pcb_center_z + PCB_THICKNESS_CM / 2
    plane_in_z      = pcb_back_z   + PLANE_IN_THICKNESS_CM / 2
    sensor_center_z = pcb_back_z   + PLANE_IN_THICKNESS_CM + SENSOR_THICKNESS_CM / 2

    # ── PCB (FR4) ─────────────────────────────────────────────────────
    pcb             = sim.add_volume("Box", "pcb")
    pcb.mother      = world.name
    pcb.size        = [5 * cm, 5 * cm, PCB_THICKNESS_CM * cm]
    pcb.translation = [0, 0, pcb_center_z * cm]
    pcb.material    = "FR4"
    pcb.color       = [0.8, 0.5, 0.2, 0.7]

    # ── Sensor (Si) ───────────────────────────────────────────────────
    sensor             = sim.add_volume("Box", "sensor")
    sensor.mother      = world.name
    sensor.size        = [2 * cm, 2 * cm, SENSOR_THICKNESS_CM * cm]
    sensor.translation = [0, 0, sensor_center_z * cm]
    sensor.material    = "G4_Si"
    sensor.color       = [0, 0, 1, 1]

    # ── Plan de contrôle (juste avant le sensor) ──────────────────────
    plane_in             = sim.add_volume("Box", "plane_in")
    plane_in.mother      = world.name
    plane_in.size        = [2 * cm, 2 * cm, PLANE_IN_THICKNESS_CM * cm]
    plane_in.translation = [0, 0, plane_in_z * cm]
    plane_in.material    = "G4_AIR"
    plane_in.color       = [1, 0, 0, 1]

    # Cuts et step limits
    sim.physics_manager.set_production_cut("phantom", "all", 0.1 * mm)
    sim.physics_manager.set_production_cut("sensor",  "all", 0.01 * mm)
    phantom.set_max_step_size(0.5 * mm)
    sensor.set_max_step_size(0.01 * mm)
    sim.physics_manager.set_user_limits_particles(["proton", "electron", "positron"])

    # ── Source ────────────────────────────────────────────────────────
    source                      = sim.add_source("GenericSource", "proton_source")
    source.particle             = "proton"
    source.n                    = n_primaries
    source.energy.type          = "mono"
    source.energy.mono          = energy_mev * MeV
    source.position.type        = "point"
    source.position.translation = [0, 0, -(phantom_depth_cm / 2 + 5) * cm]  # 5 cm avant le fantôme
    source.direction.type       = "momentum"
    source.direction.momentum   = [0, 0, 1]

    # ── Actors ────────────────────────────────────────────────────────
    stats = sim.add_actor("SimulationStatisticsActor", "stats")

    # Phase-space au niveau du sensor
    phsp                 = sim.add_actor("PhaseSpaceActor", "phsp")
    phsp.attached_to     = plane_in.name
    phsp.output_filename = Path("phsp") / f"phsp_{energy_mev:.1f}MeV.root"
    phsp.attributes      = [
        "KineticEnergy", "PostPosition", "PostDirection",
        "ParticleName", "PDGCode", "GlobalTime",
    ]

    # Dose dans le fantôme (Bragg peak)
    dose_phantom                       = sim.add_actor("DoseActor", "dose_phantom")
    dose_phantom.attached_to           = phantom.name
    dose_phantom.size                  = [1, 1, 300]
    dose_phantom.spacing               = [
        phantom.size[0],
        phantom.size[1],
        phantom.size[2] / 300,
    ]
    dose_phantom.edep.output_filename  = Path("dose") / f"phantom_edep_{energy_mev:.1f}MeV.mhd"
    dose_phantom.dose.active           = True
    dose_phantom.dose.output_filename  = Path("dose") / f"phantom_dose_{energy_mev:.1f}MeV.mhd"

    # Dose dans le sensor
    dose_sensor                       = sim.add_actor("DoseActor", "dose_sensor")
    dose_sensor.attached_to           = sensor.name
    dose_sensor.size                  = [100, 100, 10]
    dose_sensor.spacing               = [
        sensor.size[0] / 100,
        sensor.size[1] / 100,
        sensor.size[2] / 10,
    ]
    dose_sensor.edep.output_filename  = Path("dose") / f"sensor_edep_{energy_mev:.1f}MeV.mhd"
    dose_sensor.dose.active           = True
    dose_sensor.dose.output_filename  = Path("dose") / f"sensor_dose_{energy_mev:.1f}MeV.mhd"

    print(f"\n{'─'*55}")
    print(f"  Preset        : {preset}")
    print(f"  Phantom depth : {phantom_depth_cm} cm  (face arrière à z={phantom_back_z} cm)")
    print(f"  PCB centre    : z = {pcb_center_z:.3f} cm")
    print(f"  Sensor centre : z = {sensor_center_z:.3f} cm")
    print(f"  Source        : z = {-(phantom_depth_cm/2 + 5):.1f} cm")
    print(f"{'─'*55}")

    return sim, dose_phantom, dose_sensor, phsp, stats


# ─────────────────────────────────────────────
# Analyse Bragg peak
# ─────────────────────────────────────────────
def analyze_depth_dose(dose_file: Path, figure_file: Path, energy_mev: float):
    img         = sitk.ReadImage(str(dose_file))
    arr         = sitk.GetArrayFromImage(img)
    depth_dose  = arr.mean(axis=(1, 2))

    spacing = img.GetSpacing()
    origin  = img.GetOrigin()
    z       = origin[2] + spacing[2] * np.arange(len(depth_dose))
    z_cm    = (z - z[0]) / 10.0   # entrée du fantôme = 0

    peak_idx = int(np.argmax(depth_dose))
    peak_z   = z_cm[peak_idx]

    # Normalisation
    dose_norm = depth_dose / depth_dose.max() if depth_dose.max() > 0 else depth_dose

    plt.figure(figsize=(7, 4))
    plt.plot(z_cm, dose_norm, linewidth=1.5)
    plt.axvline(peak_z, color="red", linestyle="--", alpha=0.6, label=f"Bragg peak : {peak_z:.2f} cm")
    plt.xlabel("Profondeur dans le fantôme (cm)")
    plt.ylabel("Dose normalisée")
    plt.title(f"Bragg peak – {energy_mev:.1f} MeV")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(figure_file, dpi=200)
    plt.close()

    return {
        "peak_index":    peak_idx,
        "peak_depth_cm": float(peak_z),
        "peak_dose":     float(depth_dose[peak_idx]),
    }


# ─────────────────────────────────────────────
# Run
# ─────────────────────────────────────────────
def run_case(mode: str, n_primaries: int = 20000):
    assert mode in PRESETS, f"mode doit être 'min' ou 'max', pas '{mode}'"
    preset     = PRESETS[mode]
    energy_mev = preset["energy_mev"]

    repo_root   = Path(__file__).resolve().parent.parent
    case_dir    = repo_root / "simulation" / "output_first_simu" / f"case_{energy_mev:.1f}MeV"
    figures_dir = repo_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    sim, dose_phantom, dose_sensor, phsp, stats = build_simulation(
        preset     = preset,
        n_primaries = n_primaries,
        output_dir  = case_dir,
    )

    define_all_materials(sim)
    sim.run()

    phantom_dose_file = Path(dose_phantom.dose.get_output_path())
    sensor_dose_file  = Path(dose_sensor.dose.get_output_path())
    phsp_file         = Path(phsp.get_output_path())

    summary = analyze_depth_dose(
        phantom_dose_file,
        figures_dir / f"bragg_peak_{energy_mev:.1f}MeV.png",
        energy_mev,
    )

    print("\n=== Simulation finished ===")
    print(f"Energy         : {energy_mev:.1f} MeV")
    print(f"Phase-space    : {phsp_file}")
    print(f"Phantom dose   : {phantom_dose_file}")
    print(f"Sensor dose    : {sensor_dose_file}")
    print(f"Bragg peak     : {summary}")

    return {
        "energy_mev":        energy_mev,
        "phsp_file":         phsp_file,
        "phantom_dose_file": phantom_dose_file,
        "sensor_dose_file":  sensor_dose_file,
        "summary":           summary,
    }


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="FLASH-BPM simulation")
    parser.add_argument(
        "--mode",
        choices=["min", "max"],
        default="min",
        help="'min' = 62.4 MeV  |  'max' = 252.7 MeV  (default: min)",
    )
    parser.add_argument(
        "--nparticles",
        type=int,
        default=20000,
        help="Nombre de protons primaires (default: 20000)",
    )
    args = parser.parse_args()

    result = run_case(mode=args.mode, n_primaries=args.nparticles)
    print("\n=== Final summary ===")
    print(result)
