# PowHEDM_simulator
 SImulation engine to generate 2D area detector frames with powder and HEDM data

- **Multi-material simulation**: Supports simulation of multiple crystalline materials simultaneously
- **Hybrid data generation**: Combines HEDM (single crystal spots) and powder diffraction (continuous rings) patterns
- **Realistic backgrounds**: Includes Thomson scattering, Compton scattering, Poisson noise, and Gaussian noise
- **Materials Project integration**: Automatically retrieves crystal structure data via Materials Project API
- **Configurable parameters**: Extensive configuration options for experimental conditions
- **Multiple output formats**: Supports HDF5 and compressed NPZ file formats

<hr />

<ins> **Contents**</ins> :
- dir: PowHEDM_sim: Contains the classes for simulating powder, HEDM and background data
- file: PowHEDM_sim/midas_fwsim.py: Contains the class for running the MIDAS forward simulation
- file: PowHEDM_sim/mpxrd_sim.py: Contains the class for running the materials project XRD simulation
- file: PowHEDM_sim/bkg_sim.py: Contains the class for running the background simulation
- file: PowHEDM_sim/main.py: Main execultion script to run the simulation engine
- file: default_config.json: Contains the default configuration for the simulation engine
- dir: in_grain_files: Contains the grain orientations file for the simulation. It can be downloaded from [https://anl.box.com/s/3s600lzbj0hfoxmdfw1ha671rk9cw4ro]
<hr />

<ins> **Prerequisites**</ins> :

Before using this script, ensure you have the following dependencies installed:
- MIDAS installation: Github link to MIDAS [https://github.com/marinerhemant/MIDAS]
- Python 3.12 or higher
- NumPy
- Pandas
- h5py
- zarr
- pymatgen
- mp_api

<hr />

Run the powHEDM_sim/main.py script to run the simulation engine.

```bash
python PowHEDM_sim/main.py \
    -projName <project_name> \
    -mp_apikey <materials_project_api_key> \
    -nr_frames <number_of_frames>
```

## Command Line Arguments

### Required Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `-projName` | `str` | Project name (used as output directory name). This directory will be created in the current working directory. |
| `-mp_apikey` | `str` | Materials Project API key for crystal structure data |
| `-nr_frames` | `int` | Number of diffraction frames to generate |

### Optional Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `-configFile` | `str` | `"default_config.json"` | Path to configuration file |
| `-grainFile` | `str` | `"path/to/grain/file.csv"` | Path to grain orientations file |
| `-maxFramesPerSim` | `int` | `1` | Maximum omega frames per simulation |
| `-clearZarrDump` | `str` | `"yes"` | Clear zarr dump directory after simulation (`"yes"`/`"no"`) |
| `-clearHKLDump` | `str` | `"yes"` | Clear HKL dump directory after simulation (`"yes"`/`"no"`) |
| `-clearPSDump` | `str` | `"no"` | Clear parameter dump directory after simulation (`"yes"`/`"no"`) |
| `-clearGrainDump` | `str` | `"yes"` | Clear grain dump directory after simulation (`"yes"`/`"no"`) |
| `-saveHdf5` | `str` | `"no"` | Save output as HDF5 file (`"yes"`/`"no"`) |
| `-saveNpz` | `str` | `"yes"` | Save output as compressed NPZ file (`"yes"`/`"no"`) |

<hr />

## Configuration File

The configuration file is a JSON file that contains the configuration for the simulation engine. It is used to specify the experimental conditions for the simulation.

```json
{
    "MIDAS_INSTALL_PATH": "/Users/dbeniwal/MIDAS/",
    "materials_pool": {
        "oxides": [
            "SiO2", "Li2O", "Al2O3", "TiO2", "CeO2", "ZrO2", "NiO", "Fe2O3"
        ],
        "metals": [
            "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Zr", "Nb", "Mo", "Pd", "Ag", "Ta", "W", "Re", "Pt", "Au"
        ]
    },
    "nr_materials": {"min": 1, "max": 2},
    "material_types": ["oxides", "metals"],
    "hedm_always_present": "yes",
    "powder_always_present": "yes",
    "bkg_always_present": "yes",
    "detector": {"NrPixels": 2048, "px": 200},
    "exp_range": {
        "Wavelength": {"min": 0.35, "max": 0.13},
        "Lsd": {"min": 0.7e6, "max": 1.2e6},
        "BC_Ypx_offset": {"min": -100, "max": 100},
        "BC_Zpx_offset": {"min": -100, "max": 100},
        "tx": {"min": -1, "max": 1},
        "ty": {"min": -1, "max": 1},
        "tz": {"min": -1, "max": 1},
        "p0": {"min": -1e-5, "max": 1e-5},
        "p1": {"min": -1e-4, "max": 1e-4},
        "p2": {"min": -1e-3, "max": 1e-3},
        "p3": {"min": 150, "max": 150},
        "RhoD": {"min": 297006.636, "max": 297006.636},
        "EResolution": {"min": 0.0, "max": 0.0001},
        "nScans": {"min": 1, "max": 1},
        "ring_thresh": {"min": 0, "max": 60},
        "beamstop_radius": {"min": 0, "max": 100},
        "beamstop_transmission_frac": {"min": 0, "max": 0.3}
    },
    "mat_specific_range": {
        "GaussWidth": {"min": 1.3, "max": 1.7}
    },
    "hedm_range": {
        "NrGrains": {"min": 500.0, "max": 10000.0},
        "MaxIntensity": {"min": 1e3, "max": 3e3},
        "OmegaStart": {"min": -3, "max": -3},
        "OmegaEnd": {"min": 3, "max": 3},
        "OmegaStep": {"min": 0.25, "max": 0.25}
    },
    "powder_range": {
        "NrGrains": {"min": 1e6, "max": 1e6},
        "MaxIntensity": {"min": 1e3, "max": 3e3},
        "OmegaStart": {"min": -10, "max": 10},
        "OmegaEnd": {"min": -10, "max": 10},
        "OmegaStep": {"min": 2, "max": 2},
        "scale_max_intensity": {"min": 2000, "max": 4000}
    },
    "bkg": {
        "thomson": {
            "scale_Imax": {"min": 10, "max": 100},
            "form_factor_decay": {"min": 0.1, "max": 0.6}
        },
        "compton": {
            "scale_Imax": {"min": 10, "max": 100}
        },
        "poisson": {
            "scale_Imax": {"min": 10, "max": 100}
        },
        "gaussian": {
            "scale_Imax": {"min": 10, "max": 100},
            "loc": {"min": 5, "max": 20},
            "scale": {"min": 2, "max": 6}
        }
    }
}
```

<hr />

## Output

