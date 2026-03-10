"""
Matériaux personnalisés pour les simulations FLASH-BPM.
Usage : from utils.materials import define_all_materials
"""

def define_fr4(sim):
    sim.volume_manager.material_database.add_material_nb_atoms(
        "FR4",
        ["Si", "O", "C", "H"],
        [1,    2,   11,  12],
        1.85,
    )

def define_all_materials(sim):
    """Appelle tous les matériaux d'un coup."""
    define_fr4(sim)