from pathlib import Path
import opengate as gate 
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__": 

    #define Units

    m = gate.g4_units.m
    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    um = gate.g4_units.um
    MeV = gate.g4_units.MeV

    #create a simulation

    sim = gate.Simulation()

    # chemin robuste pour le dossier de sortie
    repo_root = Path(__file__).resolve().parent.parent
    output_dir = repo_root / "simulation" / "output_test_gate"
    sim.output_dir = output_dir

    sim.random_seed = "auto"
    sim.number_of_threads = 1
    sim.g4_verbose = False
    sim.visu = False
    sim.progress_bar = True
    sim.visu_type = "vrml"

    #Physics definition

    sim.physics_manager.physics_list_name ="QGSP_BIC_EMZ"

    #World definition

    world = sim.world
    world.size = [1*m,1*m,1*m]
    world.material = "G4_AIR"

    #Detector definition
    detector = sim.add_volume("Box","Detector")
    detector.mother = world.name
    detector.size = [5 * cm, 5 * cm, 1 * mm]
    detector.translation = [0, 0, 20 * cm]
    detector.material = "G4_Si"
    detector.color = [0, 0, 1, 1]

    #Control plane (detect the particles properties without affecting it)
    plane_in = sim.add_volume("Box", "plane_in")
    plane_in.mother = world.name
    plane_in.size = [5 * cm, 5 * cm, 1 * um]
    plane_in.translation = [0, 0, 19.9 * cm]
    plane_in.material = "G4_AIR"
    plane_in.color = [1, 0, 0, 1]

    #Création d'un phantom (qui représente un corps humain)
    #pour voir l'effet de la diffusion multiple sur le profil en profondeur
    phantom = sim.add_volume("Box", "phantom")
    phantom.mother = world.name
    phantom.size = [5 * cm, 5 * cm, 10 * cm]
    phantom.translation = [0, 0, 0]
    phantom.material = "G4_WATER"
    phantom.color = [0,1,1,0.3]

    #Source definition
    source = sim.add_source("GenericSource", "proton_source")
    source.particle = "proton"
    source.n = 20000
    source.energy.type = "mono"
    source.energy.mono = 100 * MeV

    source.position.type = "point" #position of the source
    source.position.translation = [0, 0, -20 * cm]

    # beam along +z
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]

    #Statistics actor 
    stats = sim.add_actor("SimulationStatisticsActor", "stats")

    #Phase space actor (write the particles properties in a root file when they cross the plane_in)
    phsp = sim.add_actor("PhaseSpaceActor", "phsp")
    phsp.attached_to = plane_in.name

    phsp.output_filename = "phsp/protons_on_detector.root"

    phsp.attributes = [
    "KineticEnergy",
    "PostPosition",
    "PostDirection",
    "ParticleName",
    "PDGCode",
    "GlobalTime",
    ]

    # -----------------------------
    # Dose actor on the detector
    # -----------------------------

    dose_detector = sim.add_actor("DoseActor", "dose_detector")

    #Dans quelle volume on veut calculer la dose : detector
    dose_detector.attached_to = detector.name

    #définition du nombre de voxel dans le volume du detector
    dose_detector.size = [1,1,50]

    #définition de la taille des voxels
    dose_detector.spacing = [
        detector.size[0],
        detector.size[1],
        detector.size[2]/50
    ]

    #énergie déposée dans le détecteur
    dose_detector.edep.output_filename = "dose/detector_edep.mhd"

    #active le calcule de dose
    dose_detector.dose.active = True

    #fichier dose détecteur
    dose_detector.dose.output_filename = "dose/detector_dose.mhd"

    # -----------------------------
    # Dose actor on the phantom
    # -----------------------------

    dose_phantom = sim.add_actor("DoseActor", "dose_phantom")

    dose_phantom.attached_to = phantom.name

    dose_phantom.size = [1,1,200]

    dose_phantom.spacing = [
        phantom.size[0],
        phantom.size[1],
        phantom.size[2]/200
    ]

    dose_phantom.edep.output_filename = "dose/phantom_edep.mhd"

    dose_phantom.dose.active = True

    dose_phantom.dose.output_filename = "dose/phantom_dose.mhd"

    # -----------------------------
    # Run
    # -----------------------------

    sim.run()

    # Print paths at the end
    print("\n=== Simulation finished ===")
    print("Phase-space ROOT file:", phsp.get_output_path())
    print("Detector Edep file:", dose_detector.edep.get_output_path())
    print("Detector Dose file:", dose_detector.dose.get_output_path())
    print("Phantom Dose file:", dose_phantom.dose.get_output_path())

    # -----------------------------
    # Read phantom dose to build Bragg peak
    # -----------------------------

    phantom_dose_file = Path(dose_phantom.dose.get_output_path())

    img = sitk.ReadImage(str(phantom_dose_file))

    dose_array = sitk.GetArrayFromImage(img)

    print("Dose array shape:", dose_array.shape)

    # moyenne sur x,y -> profil en profondeur
    depth_dose = dose_array.mean(axis=(1, 2))

    # axe z
    spacing = img.GetSpacing()
    origin = img.GetOrigin()

    z = origin[2] + spacing[2] * np.arange(len(depth_dose))

    # plot Bragg peak
    plt.figure()
    plt.plot(z, depth_dose, marker="o")
    plt.xlabel("Depth (mm)")
    plt.ylabel("Dose")
    plt.title("Bragg Peak")
    plt.grid(True)

    # dossier figures dans le repo
    figures_dir = repo_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    output_plot = figures_dir / "bragg_peak.png"

    plt.savefig(output_plot, dpi=200, bbox_inches="tight")

    print(f"Bragg peak plot saved to: {output_plot}")

    plt.close()