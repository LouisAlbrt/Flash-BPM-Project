# FLASH-BPM-Project Simulation

Monte Carlo simulation of a proton beam monitor detector system using [OpenGATE](https://opengate.readthedocs.io/) (Geant4).  
The project studies energy deposition and detector response for high-energy proton beams, with applications in medical particle therapy and radiation monitoring.

---

## Project Structure

```
FLASH-BPM-Project/
├── figures/                        # Output plots (Bragg peaks, profiles)
├── report/                         # LaTeX source for the final report
├── simulation/
│   ├── output_first_simu/          # Simulation outputs
│   │   ├── case_62.4MeV/
│   │   │   ├── dose/               # .mhd dose maps
│   │   │   └── phsp/               # .root phase-space files
│   │   └── case_252.7MeV/
│   │       ├── dose/
│   │       └── phsp/
│   └── output_test_gate/           # Early test outputs
├── utils/
│   ├── __init__.py
│   └── materials.py                # Custom Geant4 material definitions (FR4, ...)
├── first_simu.py                   # Main simulation script
├── test_gate.py                    # Early test script
└── README.md
```

---

## Physics Setup

| Parameter | Value |
|---|---|
| Physics list | `QGSP_BIC_EMZ` |
| Particle | Proton |
| Energy range | 62.4 MeV – 252.7 MeV |
| Beam profile | Point source (Gaussian planned) |
| Transport medium | Air |

---

## Geometry

```
Source → [ Water Phantom ] → [ PCB (FR4) ] → [ plane_in ] → [ Si Sensor ]
            depth adapted        1.84 mm          5 µm           0.5 mm
```

| Volume | Material | Purpose |
|---|---|---|
| Water phantom | G4_WATER | Tissue-equivalent medium, Bragg peak study |
| PCB board | FR4 (custom) | Realistic detector support |
| Silicon sensor | G4_Si | Active detector layer |
| Control plane | G4_AIR | Phase-space recording (zero physics impact) |

The phantom thickness is automatically adapted to the beam energy so that protons always reach the detector:

| Mode | Energy | Phantom depth |
|---|---|---|
| `min` | 62.4 MeV | 3.5 cm |
| `max` | 252.7 MeV | 30.0 cm |

---

## Simulation Outputs

| File | Description |
|---|---|
| `phantom_dose_*.mhd` | 3D dose map in water → Bragg peak profile |
| `phantom_edep_*.mhd` | Raw energy deposition in water |
| `sensor_dose_*.mhd` | 2D dose map in silicon sensor |
| `sensor_edep_*.mhd` | Raw energy deposition in sensor |
| `phsp_*.root` | Phase-space: position, energy, direction at sensor entrance |
| `figures/bragg_peak_*.png` | Depth-dose plots |

---

## Usage

### Run a simulation

```bash
# Low energy (62.4 MeV) — thin phantom
python first_simu.py --mode min

# High energy (252.7 MeV) — thick phantom
python first_simu.py --mode max

# Custom number of primaries
python first_simu.py --mode max --nparticles 50000
```

### Environment

```bash
conda activate gate
```

> Requires OpenGATE, SimpleITK, NumPy, Matplotlib.  
> To read `.root` phase-space files: `pip install uproot awkward`

---

## Current Status

- [x] Basic geometry (phantom + PCB + Si sensor)
- [x] Custom FR4 material definition
- [x] Dose actors on phantom and sensor
- [x] Phase-space recording at detector entrance
- [x] Bragg peak analysis and plotting
- [x] Adaptive geometry for min/max energy presets
- [ ] Gaussian beam profile
- [ ] Energy scan across full range (62.4 – 252.7 MeV)
- [ ] Lateral beam profile analysis (uproot)
- [ ] Range–energy curve (R80)
- [ ] Carbon ion beam extension
- [ ] Comparison with experimental measurements

---

## References

- [OpenGATE documentation](https://opengate.readthedocs.io/)
- [Geant4 physics lists](https://geant4.web.cern.ch/support/proc_mod_catalog/physics_lists/)
- ICRU Report 78 – Prescribing, Recording, and Reporting Proton-Beam Therapy