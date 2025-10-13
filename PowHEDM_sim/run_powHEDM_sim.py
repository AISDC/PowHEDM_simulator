import numpy as np #type: ignore
import pandas as pd #type: ignore
import h5py #type:ignore
import zarr #type:ignore

from mpxrd_sim import mpxrd
from midas_fwsim import midas_fwsim
from bkg_sim import bkg_sim

import subprocess
import json
import os, argparse, pathlib, shutil, logging, sys, copy

parser = argparse.ArgumentParser(description="PowHEDM simulator")

parser.add_argument("-projName",
                    type=str,
                    help="Project name [used as directory name for saving this model's output in 'trainOutDir' directory]")

parser.add_argument("-mp_apikey",
                    type=str,
                    help="Materials Project API key")

parser.add_argument("-nr_frames",
                    type=int,
                    help="Number of frames to simulate")

parser.add_argument("-configFile",
                    type=str,
                    default="default_config.json",
                    help="Path to configuration file for simulator")

parser.add_argument("-grainFile",
                    type=str,
                    default="in_grain_files/n1000000-id1.tess.orientations.csv",
                    help="Path to grain file")

parser.add_argument("-mpxrdCacheDir",
                    type=str,
                    default="mpxrd_cache",
                    help="Path to mpxrd cache directory")

parser.add_argument("-clearZarrDump",
                    type=str,
                    default="yes",
                    help="(no-don't clear, yes-clear) Whether to clear the zarr dump directory after simulation")

parser.add_argument("-clearPSDump",
                    type=str,
                    default="no",
                    help="(no-don't clear, yes-clear) Whether to clear the ps dump directory after simulation")

parser.add_argument("-clearGrainDump",
                    type=str,
                    default="yes",
                    help="(no-don't clear, yes-clear) Whether to clear the grain dump directory after simulation")

parser.add_argument("-saveHdf5",
                    type=str,
                    default="no",
                    help="(no-don't save, yes-save) Whether to save the data as a hdf5 file")

parser.add_argument("-saveNpz",
                    type=str,
                    default="yes",
                    help="(no-don't save, yes-save) Whether to save the data as a compressed npz file")

args, unparsed = parser.parse_known_args()

#if len(unparsed) > 0:
#    print('Unrecognized argument(s): \n%s \nProgram exiting ... ... ' % '\n'.join(unparsed))
#    exit(0)

print(f"Reading configration file '{args.configFile}'...", end="", flush=True)
sim_config_filepath = args.configFile
with open(sim_config_filepath, 'r') as f:
    powHEDM_sim_config = json.load(f)
print("DONE.")


cwd = pathlib.Path(os.getcwd())
proj_dir = os.path.join(cwd, args.projName)
zarr_dump_dir = os.path.join(proj_dir, "zarr_dump")
ps_dump_dir = os.path.join(proj_dir, "ps_dump")
grain_dump_dir = os.path.join(proj_dir, "grain_dump")
frame_dump_dir = os.path.join(proj_dir, "frame_dump")

if os.path.exists(proj_dir):
    shutil.rmtree(proj_dir)
    print(f"DELETED EXISTING PROJECT DIRECTORY: '{proj_dir}'")

os.makedirs(proj_dir)
os.makedirs(zarr_dump_dir)
os.makedirs(ps_dump_dir)
os.makedirs(grain_dump_dir)
os.makedirs(frame_dump_dir)

    
print("Created output directory: ", proj_dir)
print("Created dump directory for zarr files: ", zarr_dump_dir)
print("Created dump directory for ps files: ", ps_dump_dir)
print("Created dump directory for grains csv files: ", grain_dump_dir)

logFilePath = os.path.join(proj_dir, "simulation.log")
print(f"Setting up logging at: '{logFilePath}'", flush=True)
logging.basicConfig(filename=logFilePath, encoding="utf-8", level=logging.DEBUG)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
logging.info("\n--------------LOGGING STARTED--------------")

powHEDM_config_savepath = os.path.join(proj_dir, "powHEDM_config.json")
logging.info(f"Saving powHEDM configration file at '{powHEDM_config_savepath}'...")
with open(powHEDM_config_savepath, "w") as f:
    json.dump(powHEDM_sim_config, f)

logging.info("Setting up material pool for simulation...")
mat_pool_for_sim = []
for mat_type in powHEDM_sim_config["material_types"]:
    mat_pool_for_sim += (powHEDM_sim_config["materials_pool"][mat_type])
logging.info(f"\t|Materials: {mat_pool_for_sim}\n")

logging.info("Checking if all materials have hkls file in mpxrd cache...")
mpxrd_cache_dir = os.path.join(cwd, args.mpxrdCacheDir)
hkls_files = [i for i in os.listdir(mpxrd_cache_dir) if i.endswith("hkls.txt")]
mats_in_hkl_cache = [i.split("_")[0] for i in hkls_files]

info_files = [i for i in os.listdir(mpxrd_cache_dir) if i.endswith("info.txt")]
mats_in_info_cache = [i.split("_")[0] for i in info_files]

mats_not_in_cache = [i for i in mat_pool_for_sim if i not in mats_in_hkl_cache or i not in mats_in_info_cache]

if len(mats_not_in_cache) >= 1:
    logging.info(f"Creating cache for materials not in cache: {mats_not_in_cache}")
    for mat in mats_not_in_cache:
        mpxrd_mat = mpxrd(mp_api_key=args.mp_apikey, xray_lambda=0.5, mat_formula=mat)
        mpxrd_mat.simulate_xrd()
        mpxrd_mat.create_mp_xrdsim_dict()
        hkls_mat_filepath = os.path.join(mpxrd_cache_dir, f"{mat}_{mpxrd_mat.mp_matid}_hkls.txt")
        mpxrd_mat.save_hkls_text_midas(hkls_mat_filepath)
        sim_info_mat_savepath = os.path.join(mpxrd_cache_dir, f"{mat}_{mpxrd_mat.mp_matid}_info.txt")
        mpxrd_mat.save_sim_info_text(sim_info_mat_savepath)

logging.info(f"Reading master grain file {args.grainFile}")
grains_master = pd.read_csv(args.grainFile, header=None) # Read the grain file


det_shape = (powHEDM_sim_config["detector"]["NrPixels"], powHEDM_sim_config["detector"]["NrPixels"])
logging.info(f"Detector shape: {det_shape}\n")


frame_count = 1

while frame_count <= args.nr_frames:
    frame_files_created = []
    try:
        logging.info("=="*30)
        logging.info(f"------ Simulating frame {frame_count} ------")

        hedm_frame = np.zeros(det_shape)
        pow_frame = np.zeros(det_shape)
        bkg_frame = np.zeros(det_shape)

        nr_materials = np.random.randint(powHEDM_sim_config["nr_materials"]["min"],
                                        powHEDM_sim_config["nr_materials"]["max"] + 1)
        
        materials = np.random.choice(mat_pool_for_sim, nr_materials, replace=False).tolist()
        logging.info(f"Frame {frame_count}: Materials selected: {materials}")

        logging.info(f"Frame {frame_count}: Randomly defining experimental parameters...")
        exp_params_frame = {}
        for i_param in powHEDM_sim_config["exp_range"]:
            min, max = powHEDM_sim_config["exp_range"][i_param]["min"], powHEDM_sim_config["exp_range"][i_param]["max"]
            exp_params_frame[i_param] = np.random.uniform(min, max)

        logging.info(f"\t|Frame {frame_count}: Exp params: {exp_params_frame}")

        # Create list of booleans to determine if each material will be simulated for HEDM or Powder
        if powHEDM_sim_config["hedm_always_present"].lower() in ["y", "yes", "true"]:
            do_hedm_mat_list = nr_materials*[True]
        else:
            do_hedm_mat_list = [np.random.choice([True, False]) for i in range(nr_materials)]

        if powHEDM_sim_config["powder_always_present"].lower() in ["y", "yes", "true"]:
            do_powder_mat_list = nr_materials*[True]
        else:
            do_powder_mat_list = [np.random.choice([True, False]) for i in range(nr_materials)]

        #Ensure that at least one of the "powder" or "hedm" will be simulated for each material
        for i in range(nr_materials):
            if do_powder_mat_list[i] == False and do_hedm_mat_list[i] == False:
                rand_choice = np.random.choice(["powder", "hedm"])
                if rand_choice == "powder": do_powder_mat_list[i] = True
                if rand_choice == "hedm": do_hedm_mat_list[i] = True
        
        # Determine if background will be simulated (once per simulation, not per material like HEDM and Powder)
        if powHEDM_sim_config["bkg_always_present"].lower() in ["y", "yes", "true"]:
            do_bkg = True
        else:
            do_bkg = np.random.choice([True, False])

        # Sample the omega range for HEDM
        omega_start_hedm = np.random.uniform(powHEDM_sim_config["hedm_range"]["OmegaStart"]["min"],
                                                        powHEDM_sim_config["hedm_range"]["OmegaStart"]["max"])
            
        omega_end_hedm = np.random.uniform(powHEDM_sim_config["hedm_range"]["OmegaEnd"]["min"],
                                                    powHEDM_sim_config["hedm_range"]["OmegaEnd"]["max"])
        
        omega_step_hedm = np.random.uniform(powHEDM_sim_config["hedm_range"]["OmegaStep"]["min"],
                                                        powHEDM_sim_config["hedm_range"]["OmegaStep"]["max"])

        # Sample the omega range for Powder
        omega_start_powder = np.random.uniform(powHEDM_sim_config["powder_range"]["OmegaStart"]["min"],
                                                    powHEDM_sim_config["powder_range"]["OmegaStart"]["max"])
            
        omega_end_powder = np.random.uniform(powHEDM_sim_config["powder_range"]["OmegaEnd"]["min"],
                                                        powHEDM_sim_config["powder_range"]["OmegaEnd"]["max"])
            
        omega_step_powder = np.random.uniform(powHEDM_sim_config["powder_range"]["OmegaStep"]["min"],
                                                        powHEDM_sim_config["powder_range"]["OmegaStep"]["max"]) 
        
        frame_BC_Ypx = int(powHEDM_sim_config["detector"]["NrPixels"]/2) + exp_params_frame["BC_Ypx_offset"]
        frame_BC_Zpx = int(powHEDM_sim_config["detector"]["NrPixels"]/2) + exp_params_frame["BC_Zpx_offset"]

        for (mat_i, mat) in enumerate(materials):

            do_hedm_mat = do_hedm_mat_list[mat_i]
            do_powder_mat = do_powder_mat_list[mat_i]
            
            logging.info(f"\t|Frame {frame_count}: For {mat}, do_powder: {do_powder_mat}, do_hedm: {do_hedm_mat}")

            gauss_width_mat = np.random.uniform(powHEDM_sim_config["mat_specific_range"]["GaussWidth"]["min"],
                                                powHEDM_sim_config["mat_specific_range"]["GaussWidth"]["max"])
            
            logging.info(f"\t|Frame {frame_count}: For {mat}: Finding info file in mpxrd cache...")
            mpxrd_info_mat_filename = [i for i in os.listdir(mpxrd_cache_dir) if \
                                i.endswith("info.txt") and i.split("_")[0] == mat][0]

            mpxrd_info_mat_filepath = os.path.join(mpxrd_cache_dir, mpxrd_info_mat_filename)
            mpxrd_info_mat_dict = {}

            with open(mpxrd_info_mat_filepath, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    key, value = line.split(" ", 1)

                    # Try to convert to appropriate data type
                    try:
                        if '.' not in value: mpxrd_info_mat_dict[key] = int(value)
                        else: mpxrd_info_mat_dict[key] = float(value)
                    except ValueError:
                        # Keep as string if conversion fails
                        mpxrd_info_mat_dict[key] = value.strip()

            #--------------------------------
            # Create lattice parameters sting in right format for ps file
            lattice_params = [mpxrd_info_mat_dict["LATTICE_A"], mpxrd_info_mat_dict["LATTICE_B"], mpxrd_info_mat_dict["LATTICE_C"],\
                    mpxrd_info_mat_dict["LATTICE_ALPHA"], mpxrd_info_mat_dict["LATTICE_BETA"], mpxrd_info_mat_dict["LATTICE_GAMMA"]]

            lattice_params = [str(i) for i in lattice_params]
            lattice_params = " ".join(lattice_params)
            
            #--------------------------------
            logging.info(f"\t|Frame {frame_count}: For {mat}: Finding hkls file in mpxrd cache...")
            hkls_mat_filename = [i for i in os.listdir(mpxrd_cache_dir) if \
                                i.endswith("hkls.txt") and i.split("_")[0] == mat][0]

            hkls_mat_filepath = os.path.join(mpxrd_cache_dir, hkls_mat_filename)
            logging.info(f"\t|hkl file found at: {hkls_mat_filepath}")


            # Find max nr. of rings that can fit on the detector
            Rmax = np.sqrt(np.max([frame_BC_Ypx, powHEDM_sim_config["detector"]["NrPixels"]-frame_BC_Ypx])**2 +
                            np.max([frame_BC_Zpx, powHEDM_sim_config["detector"]["NrPixels"]-frame_BC_Zpx])**2)
            Rmax = powHEDM_sim_config["detector"]["px"] * Rmax
            two_theta_max = np.arcsin(Rmax / exp_params_frame["Lsd"]) # in radians
            dmin = exp_params_frame["Wavelength"] / (2 * np.sin(two_theta_max/2))

            with open(hkls_mat_filepath, "r") as f:
                lines = f.readlines()
                for (line_idx, line) in enumerate(lines[1:]):
                    if line.startswith("#"):
                        continue
                    dval = float(line.split(" ")[4])
                    if dval < dmin:
                        nrRings = int(line.split(" ")[5])
                        break
            

            if do_hedm_mat:
                logging.info(f"\t|Frame {frame_count}: For {mat}: Running HEDM simulation...")

                #---------------DO HEDM SIMULATION---------------
                # Define the experimental parameters for the HEDM simulation
                nr_grains_hedm_mat = np.random.randint(powHEDM_sim_config["hedm_range"]["NrGrains"]["min"],
                                                    powHEDM_sim_config["hedm_range"]["NrGrains"]["max"] + 1)
                
                max_int_hedm_mat = np.random.uniform(powHEDM_sim_config["hedm_range"]["MaxIntensity"]["min"],
                                                    powHEDM_sim_config["hedm_range"]["MaxIntensity"]["max"])
            
                #--------------------------------
                # Sample the grains
                grains_use_frac = nr_grains_hedm_mat / len(grains_master[9:])
                grains_hedm_mat = grains_master[9:].sample(frac=grains_use_frac)

                # Concatenate the sampled grains to the new grain file
                grains_hedm_mat = pd.concat([grains_master[:9], grains_hedm_mat], ignore_index=True)
                
                # Update the first row to include the number of sampled grains
                grains_hedm_mat.iloc[0] = f"{grains_hedm_mat.iloc[0].values[0].split(' ')[0]} {str(nr_grains_hedm_mat)}"

                # Save the sampled grains csv file
                grains_hedm_mat_savepath = os.path.join(grain_dump_dir, f"Frame{frame_count}_HEDM_{mat}.csv")
                grains_hedm_mat.to_csv(grains_hedm_mat_savepath, header=False, index=False)

                ps_hedm_mat_filepath = os.path.join(ps_dump_dir, f"Frame{frame_count}_HEDM_{mat}_{nr_grains_hedm_mat}.txt")
                zarr_hedm_mat_filename = os.path.join(zarr_dump_dir, f"Frame{frame_count}_HEDM_{mat}")

                # =================== WRITE HEADER FILE ===================
                with open(ps_hedm_mat_filepath, 'w') as f:
                    f.write("# Parameter file for forward simulation\n")
                    f.write(f"InFileName {grains_hedm_mat_savepath}\n")
                    f.write(f"IntensitiesFile {hkls_mat_filepath}\n")
                    f.write(f"OutFileName {zarr_hedm_mat_filename}\n")
                    f.write(f"LatticeConstant {lattice_params}\n")
                    f.write(f"SpaceGroup {mpxrd_info_mat_dict["SPACE_GROUP"]}\n")
                    f.write(f"Wavelength {exp_params_frame['Wavelength']}\n")
                    f.write(f"Lsd {exp_params_frame['Lsd']}\n")
                    f.write(f"BC {frame_BC_Ypx} {frame_BC_Zpx}\n")
                    f.write(f"tx {exp_params_frame['tx']}\n")
                    f.write(f"ty {exp_params_frame['ty']}\n")
                    f.write(f"tz {exp_params_frame['tz']}\n")
                    f.write(f"p0 {exp_params_frame['p0']}\n")
                    f.write(f"p1 {exp_params_frame['p1']}\n")
                    f.write(f"p2 {exp_params_frame['p2']}\n")
                    f.write(f"p3 {exp_params_frame['p3']}\n")
                    f.write(f"RhoD {exp_params_frame['RhoD']}\n")
                    f.write(f"PeakIntensity {max_int_hedm_mat}\n")
                    f.write(f"GaussWidth {gauss_width_mat}\n")
                    f.write(f"EResolution {exp_params_frame['EResolution']} 10\n")
                    f.write(f"WriteSpots 0\n")
                    f.write(f"nScans {exp_params_frame['nScans']}\n")
                    f.write(f"OmegaStart {omega_start_hedm}\n")
                    f.write(f"OmegaEnd {omega_end_hedm}\n")
                    f.write(f"OmegaStep {omega_step_hedm}\n")
                    f.write(f"NrPixels {powHEDM_sim_config["detector"]["NrPixels"]}\n")
                    f.write(f"px {powHEDM_sim_config["detector"]["px"]}\n")

                    for ring_i in range(1, nrRings+1):
                        ring_thresh = np.random.randint(powHEDM_sim_config["exp_range"]["ring_thresh"]["min"],
                                                        powHEDM_sim_config["exp_range"]["ring_thresh"]["max"] + 1)
                        
                        f.write(f"RingThresh {ring_i} {ring_thresh}\n")

                mat_hedm_fwsim = midas_fwsim(ps_hedm_mat_filepath, powHEDM_sim_config["MIDAS_INSTALL_PATH"], zarr_dump_dir)
                mat_hedm_fwsim.run_midas_forward_sim()

                zarr_hedm_mat_filepath = os.path.join(zarr_dump_dir, f"Frame{frame_count}_HEDM_{mat}_scanNr_0.zip")
                mat_hedm_arr = zarr.open(zarr_hedm_mat_filepath, "r")["exchange"]["data"]
                rand_idx = np.random.randint(0, len(mat_hedm_arr))

                hedm_frame += mat_hedm_arr[rand_idx, :, :]

                logging.info(f"\t|Frame {frame_count}: For {mat}: HEDM simulation complete.")

                if (args.clearGrainDump).lower() in ["y", "yes", "true"]:
                    if os.path.exists(grains_hedm_mat_savepath): os.remove(grains_hedm_mat_savepath)
                else:
                    frame_files_created.append(grains_hedm_mat_savepath)

                if (args.clearPSDump).lower() in ["y", "yes", "true"]:
                    if os.path.exists(ps_hedm_mat_filepath): os.remove(ps_hedm_mat_filepath)
                else:
                    frame_files_created.append(ps_hedm_mat_filepath)

                if (args.clearZarrDump).lower() in ["y", "yes", "true"]:
                    if os.path.exists(zarr_hedm_mat_filepath): os.remove(zarr_hedm_mat_filepath)
                else:
                    frame_files_created.append(zarr_hedm_mat_filepath)

            if do_powder_mat:

                logging.info(f"\t|Frame {frame_count}: For {mat}: Running powder simulation...")
                #---------------DO POWDER SIMULATION---------------
                # Define the experimental parameters for the powder simulation
                nr_grains_powder_mat = np.random.randint(powHEDM_sim_config["powder_range"]["NrGrains"]["min"],
                                                        powHEDM_sim_config["powder_range"]["NrGrains"]["max"] + 1)
                
                max_int_powder_mat = np.random.uniform(powHEDM_sim_config["powder_range"]["MaxIntensity"]["min"],
                                                        powHEDM_sim_config["powder_range"]["MaxIntensity"]["max"])
                
                max_int_powder_mat = np.random.uniform(powHEDM_sim_config["powder_range"]["scale_max_intensity"]["min"],
                                                        powHEDM_sim_config["powder_range"]["scale_max_intensity"]["max"])

                #--------------------------------
                # Sample the grains
                grains_use_frac = nr_grains_powder_mat / len(grains_master[9:])
                grains_powder_mat = grains_master[9:].sample(frac=grains_use_frac)

                # Concatenate the sampled grains to the new grain file
                grains_powder_mat = pd.concat([grains_master[:9], grains_powder_mat], ignore_index=True)
                
                # Update the first row to include the number of sampled grains
                grains_powder_mat.iloc[0] = f"{grains_powder_mat.iloc[0].values[0].split(' ')[0]} {str(nr_grains_powder_mat)}"

                # Save the sampled grains csv file
                grains_powder_mat_savepath = os.path.join(grain_dump_dir, f"Frame{frame_count}_powder_{mat}.csv")
                grains_powder_mat.to_csv(grains_powder_mat_savepath, header=False, index=False)


                ps_powder_mat_filepath = os.path.join(ps_dump_dir, f"Frame{frame_count}_powder_{mat}_{nr_grains_powder_mat}.txt")
                zarr_powder_mat_filename = os.path.join(zarr_dump_dir, f"Frame{frame_count}_powder_{mat}")

                # =================== WRITE HEADER FILE ===================
                with open(ps_powder_mat_filepath, 'w') as f:
                    f.write("# Parameter file for forward simulation\n")
                    f.write(f"InFileName {grains_powder_mat_savepath}\n")
                    f.write(f"IntensitiesFile {hkls_mat_filepath}\n")
                    f.write(f"OutFileName {zarr_powder_mat_filename}\n")
                    f.write(f"LatticeConstant {lattice_params}\n")
                    f.write(f"SpaceGroup {mpxrd_info_mat_dict["SPACE_GROUP"]}\n")
                    f.write(f"Wavelength {exp_params_frame['Wavelength']}\n")
                    f.write(f"Lsd {exp_params_frame['Lsd']}\n")
                    f.write(f"BC {frame_BC_Ypx} {frame_BC_Zpx}\n")
                    f.write(f"tx {exp_params_frame['tx']}\n")
                    f.write(f"ty {exp_params_frame['ty']}\n")
                    f.write(f"tz {exp_params_frame['tz']}\n")
                    f.write(f"p0 {exp_params_frame['p0']}\n")
                    f.write(f"p1 {exp_params_frame['p1']}\n")
                    f.write(f"p2 {exp_params_frame['p2']}\n")
                    f.write(f"p3 {exp_params_frame['p3']}\n")
                    f.write(f"RhoD {exp_params_frame['RhoD']}\n")
                    f.write(f"PeakIntensity {max_int_powder_mat}\n")
                    f.write(f"GaussWidth {gauss_width_mat}\n")
                    f.write(f"EResolution {exp_params_frame['EResolution']} 10\n")
                    f.write(f"WriteSpots 0\n")
                    f.write(f"nScans {exp_params_frame['nScans']}\n")
                    f.write(f"OmegaStart {omega_start_powder}\n")
                    f.write(f"OmegaEnd {omega_end_powder}\n")
                    f.write(f"OmegaStep {omega_step_powder}\n")
                    f.write(f"NrPixels {powHEDM_sim_config["detector"]["NrPixels"]}\n")
                    f.write(f"px {powHEDM_sim_config["detector"]["px"]}\n")
                    for ring_i in range(1, nrRings+1):
                        ring_thresh = np.random.randint(powHEDM_sim_config["exp_range"]["ring_thresh"]["min"],
                                                        powHEDM_sim_config["exp_range"]["ring_thresh"]["max"] + 1)
                        
                        f.write(f"RingThresh {ring_i} {ring_thresh}\n")
                
                mat_powder_fwsim = midas_fwsim(ps_powder_mat_filepath, powHEDM_sim_config["MIDAS_INSTALL_PATH"], zarr_dump_dir)
                mat_powder_fwsim.run_midas_forward_sim()

                zarr_powder_mat_filepath = os.path.join(zarr_dump_dir, f"Frame{frame_count}_powder_{mat}_scanNr_0.zip")
                mat_powder_arr = zarr.open(zarr_powder_mat_filepath, "r")["exchange"]["data"]
                mat_powder_arr = np.sum(mat_powder_arr, axis=0)
                mat_powder_arr = mat_powder_arr * max_int_powder_mat / np.max(mat_powder_arr)

                pow_frame += mat_powder_arr
                logging.info(f"\t|Frame {frame_count}: For {mat}: Powder simulation complete.")

                # Clear the dump files for the powder simulation
                if (args.clearGrainDump).lower() in ["y", "yes", "true"]:
                    if os.path.exists(grains_powder_mat_savepath): os.remove(grains_powder_mat_savepath)
                else:
                    frame_files_created.append(grains_powder_mat_savepath)

                if (args.clearPSDump).lower() in ["y", "yes", "true"]:
                    if os.path.exists(ps_powder_mat_filepath): os.remove(ps_powder_mat_filepath)
                else:
                    frame_files_created.append(ps_powder_mat_filepath)

                if (args.clearZarrDump).lower() in ["y", "yes", "true"]:
                    if os.path.exists(zarr_powder_mat_filepath): os.remove(zarr_powder_mat_filepath)
                else:
                    frame_files_created.append(zarr_powder_mat_filepath)

        if do_bkg:
            logging.info(f"\t|Frame {frame_count}: Running background simulation...")

            # Add background to the frames
            bkg_types = list(powHEDM_sim_config["bkg"].keys())
            bkg = bkg_sim(BC_Ypx=frame_BC_Ypx,
                            BC_Zpx=frame_BC_Zpx,
                            xray_lambda=exp_params_frame["Wavelength"],
                            Lsd=exp_params_frame["Lsd"],
                            px_size=powHEDM_sim_config["detector"]["px"],
                            nr_pixels=powHEDM_sim_config["detector"]["NrPixels"],
                            beamstop_radius=exp_params_frame["beamstop_radius"],
                            beamstop_transmission_frac=exp_params_frame["beamstop_transmission_frac"])


            if "thomson" in bkg_types:
                thomson_ff_decay = np.random.uniform(powHEDM_sim_config["bkg"]["thomson"]["form_factor_decay"]["min"],
                                                powHEDM_sim_config["bkg"]["thomson"]["form_factor_decay"]["max"])

                thomson_scale_Imax = np.random.uniform(powHEDM_sim_config["bkg"]["thomson"]["scale_Imax"]["min"],
                                                    powHEDM_sim_config["bkg"]["thomson"]["scale_Imax"]["max"])

                logging.info(f"\t|Frame {frame_count}: Including thomson scattering.")

                bkg.thomson_scattering(form_factor_decay=thomson_ff_decay, scale_Imax=thomson_scale_Imax)
                bkg_frame += bkg.thomson_bkg["intensity"]

            else:
                thomson_ff_decay = None
                thomson_scale_Imax = None
                

            if "compton" in bkg_types:
                compton_scale_Imax = np.random.uniform(powHEDM_sim_config["bkg"]["compton"]["scale_Imax"]["min"],
                                                    powHEDM_sim_config["bkg"]["compton"]["scale_Imax"]["max"])

                bkg.compton_scattering(scale_Imax=compton_scale_Imax)
                bkg_frame += bkg.compton_bkg["intensity"]

            else:
                compton_scale_Imax = None

            if "poisson" in bkg_types:
                poisson_scale_Imax = np.random.uniform(powHEDM_sim_config["bkg"]["poisson"]["scale_Imax"]["min"],
                                                        powHEDM_sim_config["bkg"]["poisson"]["scale_Imax"]["max"])

                logging.info(f"\t|Frame {frame_count}: Including poisson noise.")

                bkg.poisson_noise(intensity_arr=hedm_frame + pow_frame,
                                scale_Imax=poisson_scale_Imax)
                
                bkg_frame += bkg.poisson_noise["intensity"]

            else:
                poisson_scale_Imax = None

            if "gaussian" in bkg_types:

                gaussian_scale_Imax = np.random.uniform(powHEDM_sim_config["bkg"]["gaussian"]["scale_Imax"]["min"],
                                                            powHEDM_sim_config["bkg"]["gaussian"]["scale_Imax"]["max"])

                gaussian_loc = np.random.uniform(powHEDM_sim_config["bkg"]["gaussian"]["loc"]["min"],
                                                        powHEDM_sim_config["bkg"]["gaussian"]["loc"]["max"])

                gaussian_scale = np.random.uniform(powHEDM_sim_config["bkg"]["gaussian"]["scale"]["min"],
                                                        powHEDM_sim_config["bkg"]["gaussian"]["scale"]["max"])

                logging.info(f"\t|Frame {frame_count}: Including gaussian noise.")

                bkg.gaussian_noise(loc=gaussian_loc, scale=gaussian_scale, scale_Imax=gaussian_scale_Imax)
                
                bkg_frame += bkg.gaussian_noise["intensity"]

            else:
                gaussian_scale_Imax = None
                gaussian_loc = None
                gaussian_scale = None
        

            ps_bkg_filepath = os.path.join(ps_dump_dir, f"Frame{frame_count}_bkg.txt")

            if (args.clearPSDump).lower() in ["y", "yes", "true"]:
                if os.path.exists(ps_bkg_filepath): os.remove(ps_bkg_filepath)
            else:
                frame_files_created.append(ps_bkg_filepath)

                # =================== WRITE HEADER FILE ===================
                with open(ps_bkg_filepath, 'w') as f:
                    f.write("# Parameter file for background simulation\n")
                    f.write(f"thomson_ff_decay_list {thomson_ff_decay}\n")
                    f.write(f"thomson_scale_Imax_list {thomson_scale_Imax}\n")
                    f.write(f"compton_scale_Imax_list {compton_scale_Imax}\n")
                    f.write(f"gaussian_scale_Imax_list {gaussian_scale_Imax}\n")
                    f.write(f"gaussian_loc_list {gaussian_loc}\n")
                    f.write(f"gaussian_scale_list {gaussian_scale}\n")
                    f.write(f"beamstop_radius {exp_params_frame['beamstop_radius']}\n")
                    f.write(f"beamstop_transmission_frac {exp_params_frame['beamstop_transmission_frac']}\n")


        if (args.saveNpz).lower() in ["y", "yes", "true"]:
            frame_out_npz_filename = f"frame_{frame_count}.npz"
            frame_out_npz_filepath = os.path.join(frame_dump_dir, frame_out_npz_filename)
            logging.info(f"Saving frame to npz file at: '{frame_out_npz_filepath}'")
            np.savez_compressed(frame_out_npz_filepath, 
                                HEDM=hedm_frame, 
                                Powder=pow_frame, 
                                Background=bkg_frame)

            frame_files_created.append(frame_out_npz_filepath)

        if (args.saveHdf5).lower() in ["y", "yes", "true"]:
            frame_out_hdf5_filename = f"frame_{frame_count}.h5"
            frame_out_hdf5_filepath = os.path.join(frame_dump_dir, frame_out_hdf5_filename)
            logging.info(f"Saving frame to hdf5 file at: '{frame_out_hdf5_filepath}'")

            with h5py.File(frame_out_hdf5_filepath, "w") as h5file:
                #Saving the HEDM, powder and background frames
                h5file.create_dataset("HEDM", data=hedm_frame)
                h5file.create_dataset("Powder", data=pow_frame)
                h5file.create_dataset("Background", data=bkg_frame)     

            frame_files_created.append(frame_out_hdf5_filepath)


        frame_count += 1
    
    except:
        logging.info("\t" + "x"*30)
        logging.info(f"\tEXCEPTION: Error while simulating frame {frame_count}. Reinitializing...")
        logging.info("\t" + "x"*30)

        if len(frame_files_created) > 0:
            for file in frame_files_created:
                if os.path.exists(file): os.remove(file)

    try:
        if os.path.exists("SpotMatrixGen.csv"):
            os.remove("SpotMatrixGen.csv")
    except:
        pass

    try:
        if os.path.exists("GrainsGen.csv"):
            os.remove("GrainsGen.csv")
    except:
        pass
    
    try:   
        if os.path.exists("GrainsGen.csv"):
            os.remove("GrainsGen.csv")
    except:
        pass

    try:
        if os.path.exists("hkls.csv"):
            os.remove("hkls.csv")
    except:
        pass


logging.info("\n\n")
logging.info("="*30)
logging.info("SIMULATION COMPLETE.")
logging.info("="*30)
