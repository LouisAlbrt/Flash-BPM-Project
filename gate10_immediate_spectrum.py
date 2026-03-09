import opengate as gate
from opengate.utility import g4_units

import uproot
import numpy as np
import matplotlib.pyplot as plt

import os

""" Demonstration of GATE10 for simulation Energy Spectra

Documentation:
https://opengate-python.readthedocs.io

"""

###################### CONTROL PANEL #######################

# ROOT files + plots
output_path = 'output'
world_material = "Air"

# --- Source Parameters
particle = "alpha"                      #'e-', 'gamma', 'proton', 'neutron', 'alpha', 'ion'
Z, A = 6, 12                            # Only necessary for ions
energies = [5.48556, 5.4438, 5.388]       # MeV (MeV/u for ions) - You can also give multiple values as a list
intensities = [0.848, 0.131, 0.0166]
numParticles = 1000    
source_size = 2.5                              # mm (Radius of disc source)

# --- Detector Parameters
size_lat = 1                    #mm  
detector_thickness = 10                  #um
detector_material = "Diamond"   #See data/GateMaterials.db #Silicon, Diamond, SiliconCarbide

# --- Target Parameters
phantom = True
phantom_material = "Water"
phantom_thickness = 1.4e-4      #cm, water equivalent thickness before detector : to simulate the energy loss before reaching the detector. Set to 0 for no phantom.

# --- Plot Parameters (Histogram)
numBins = 100
range = [2500, 4000] #keV - Set this up to your expected range

###################################################################################################

# --- CONSTANTS

EPSILON_DICT = {
    "Silicon": 3.62,
    "Diamond": 13.1,
    "SiliconCarbide": 7.84
}

ELECTRON_CHARGE = 1.602176634e-4 # fC

# GATE Unit definitions
m = g4_units.m
cm = g4_units.cm
keV = g4_units.keV
MeV = g4_units.MeV
mm = g4_units.mm
um = g4_units.um
Bq = g4_units.Bq

###################################################################################################
###################################################################################################
###################################################################################################

def _create_single_source(particle, energy, source_size, i):

    # Create a unique name for each source based on particle type, energy, and index

    SOURCE_NAME = f"{particle}_E{energy:.2f}MeV_{i}"

    if particle == "ion":
        source = sim.add_source("GenericSource", SOURCE_NAME)
        source.particle = f"ion {int(Z)} {int(A)}"
    else:
        source = sim.add_source("GenericSource", SOURCE_NAME)
        source.particle = particle

    source.position.type = "disc"
    source.position.radius = source_size * mm
    source.energy.mono = energy * MeV
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]

    return source

def build_geometry(size_lat, detector_thickness, detector_material, phantom, phantom_material, phantom_thickness):

    detector = sim.add_volume("Box", "Detector")
    detector.size = [size_lat * mm, size_lat * mm, detector_thickness * um]

    if phantom:
        detector.translation = [0 * cm, 0 * cm, phantom_thickness * cm + detector_thickness * um/2]
    else:
        detector.translation = [0 * cm, 0 * cm, detector_thickness * um/2]

    detector.material = detector_material

    # Add a target between source and detector
    if phantom:
        phantom = sim.add_volume("Box", "Phantom")
        phantom.size = [10 * cm, 10 * cm, phantom_thickness * cm]
        phantom.translation = [0 * cm, 0 * cm, phantom_thickness * cm/2]
        phantom.material = phantom_material

def build_source(particle, energies, intensities, numParticles, source_size):

    if len(energies) == 1:
        source = _create_single_source(particle, energies[0], source_size, 0)
        source.n = numParticles

    elif len(energies) > 1:
        for i, e in enumerate(energies):
            source = _create_single_source(particle, e, source_size, i)
            source.n = numParticles * intensities[i] # Scale number of particles by intensity

def build_edep_actor(detector):

    energy_str = "_".join([str(e) for e in energies]) if energies is list else str(energies)
    OUTPUT_NAME = f'{output_path}/{particle}_{energy_str}_{detector_thickness}um'

    # HITS COLLECTION = Every event
    hc = sim.add_actor('DigitizerHitsCollectionActor', 'Hits') #acteur = un truc qui va prendre de mesures
    hc.attached_to = ['Detector'] #can be attached to volumes, sources, or the world (ici le detecteur)
    hc.output_filename = f'{OUTPUT_NAME}_hits.root'
    hc.attributes = ['TotalEnergyDeposit', 'PostPosition', 'GlobalTime', 'PreStepUniqueVolumeID'] 

    # SINGLES = Primary + Secondaries 

    #c'est assez important, en gros les detecteurs vont voir plusieurs interactions (hits) pour une même particule primaire, mais on veut aussi avoir une collection qui regroupe tous les hits d'une même particule primaire en un seul "single". C'est ce que fait le "DigitizerReadoutActor" avec la policy "EnergyWeightedCentroidPosition" : il va faire la somme de l'energie depositée par tous les hits d'une même particule primaire et calculer une position moyenne pondérée par l'energie deposée.
    sc = sim.add_actor("DigitizerReadoutActor", "Singles")
    sc.output_filename = f"{OUTPUT_NAME}_singles.root"
    sc.input_digi_collection = "Hits"
    sc.group_volume = detector.name
    sc.discretize_volume = detector.name
    sc.policy = "EnergyWeightedCentroidPosition"

def immediate_testplot(numBins, range):

    energy_str = "_".join([str(e) for e in energies]) if energies is list else str(energies)
    OUTPUT_NAME = f'{output_path}/{particle}_{energy_str}_{detector_thickness}um'

    # Open the file
    file = uproot.open(f"{OUTPUT_NAME}_singles.root")
    # Get the right branch and convert to numpy array
    data = np.asarray(file['Singles']['TotalEnergyDeposit']) * 1e3 # MeV to keV

    file_hits = uproot.open(f"{OUTPUT_NAME}_hits.root")
    data_hits = np.asarray(file_hits['Hits']['TotalEnergyDeposit']) * 1e3 # MeV to keV

    # Statistics
    print(f"Mean Energy: {np.mean(data):.2f} keV")
    print(f"Most Probable Energy: {(np.histogram(data, bins=numBins, range=range)[0].argmax()) * (range[1]-range[0])/numBins + range[0]:.2f} keV")
    print(f"Max Energy: {np.max(data):.2f} keV")
    print(f"Min Energy: {np.min(data):.2f} keV")
    print(f"Total Counts (Singles): {len(data)}")
    print(f"Total Counts (Hits): {len(data_hits)}")

    # --- Nice plots
    fig, ax = plt.subplots()
    n, bins, _ = ax.hist(data, bins=numBins, label=f'{numParticles} Primaries', range = range, density=False, alpha=0.7)
    n2, bins2, _ = ax.hist(data_hits, bins=numBins, label=f'{numParticles} Primaries (Hits)', range = range, density=False, alpha=0.5)

    diff = np.abs(n - n2)
    ax.step(bins[:-1], diff, where='post', label='Difference', color='red', linestyle='dashed')

    # Second x-axis as charge (fC)
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    epsilon = EPSILON_DICT[detector_material]  # eV
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xticklabels([f"{(tick*1e3)/epsilon * ELECTRON_CHARGE:.2f}" for tick in ax.get_xticks()]) # Q = (E_dep/epsilon) * e_charge, with E in eV and epsilon in eV, we convert E from keV to eV by multiplying by 1e3, and then multiply by electron charge to get charge in fC.
    ax2.set_xlabel('Charge [fC]')
    ax.set_xlabel('Energy [keV]')
    ax.set_ylabel('Counts')

    ax.set_xlim(range)
    ax.set_ylim([0, 1.1*max(ax.get_ylim())])
    
    # Most probably value
    most_probable = (np.histogram(data, bins=numBins, range=range)[0].argmax()) * (range[1]-range[0])/numBins + range[0]
    ax.axvline(most_probable, color='r', linestyle='dashed', linewidth=1)
    
    # Mean value
    mean_value = np.mean(data)
    ax.axvline(mean_value, color='g', linestyle='dashed', linewidth=1)

    # Min/Max value
    ax.axvline(np.min(data), color='b', linestyle='dashed', linewidth=1)
    ax.axvline(np.max(data), color='b', linestyle='dashed', linewidth=1)

    # Add all text annotations to the legend
    handles = []
    handles.append(ax.axvline(most_probable, color='r', linestyle='dashed', linewidth=1, label=f'Most Probable: {most_probable:.2f} keV'))
    handles.append(ax.axvline(mean_value, color='g', linestyle='dashed', linewidth=1, label=f'Mean: {mean_value:.2f} keV'))
    handles.append(ax.axvline(np.min(data), color='b', linestyle='dashed', linewidth=1, label=f'Min: {np.min(data):.2f} keV'))
    handles.append(ax.axvline(np.max(data), color='b', linestyle='dashed', linewidth=1, label=f'Max: {np.max(data):.2f} keV'))
    ax.legend(handles=handles + [ax.patches[0]], loc='best')
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_NAME}.png', dpi=300)
    plt.show()

if __name__ == "__main__":

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Create simulation
    sim = gate.Simulation()
    # Add material database
    sim.volume_manager.add_material_database("data/GateMaterials.db")

    sim.world.material = world_material

    phys = sim.physics_manager
    phys.name = 'QGSP_BERT_EMZ' # Choose approriate physics list: https://opengate-python.readthedocs.io/en/master/user_guide/user_guide_advanced.html#physics-lists-details-label

    # --- Build your simulation
    build_geometry(size_lat, detector_thickness, detector_material, phantom, phantom_material, phantom_thickness)
    build_source(particle, energies, intensities, numParticles, source_size)
    build_edep_actor(detector=sim.volume_manager.volumes["Detector"])

    # Enable visualization (optional)
    # Reduce number of particles!!!
    #sim.visu = True

    # --- Run the simulation
    sim.run(start_new_process=True)

    # Plot results from ROOT file
    immediate_testplot(numBins, range)