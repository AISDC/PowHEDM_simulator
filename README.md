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
    -mp_apikey <materials_project_api_key>
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

