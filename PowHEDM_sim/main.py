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
                    default="/Users/dbeniwal/ANL-research/github/PowHEDM_simulator/in_grain_files/n1000000-id1.tess.orientations.csv",
                    help="Path to grain file")

parser.add_argument("-maxFramesPerSim",
                    type=int,
                    default=1,
                    help="Maximum number of omega frames to use per simulation")

parser.add_argument("-clearZarrDump",
                    type=str,
                    default="yes",
                    help="(no-don't clear, yes-clear) Whether to clear the zarr dump directory after simulation")

parser.add_argument("-clearHKLDump",
                    type=str,
                    default="yes",
                    help="(no-don't clear, yes-clear) Whether to clear the hkl dump directory after simulation")

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


if len(unparsed) > 0:
    print('Unrecognized argument(s): \n%s \nProgram exiting ... ... ' % '\n'.join(unparsed))
    exit(0)


print(f"Reading configration file '{args.configFile}'...", end="", flush=True)
sim_config_filepath = args.configFile
with open(sim_config_filepath, 'r') as f:
    powHEDM_sim_config = json.load(f)
print("DONE.")


cwd = pathlib.Path(os.getcwd())
proj_dir = os.path.join(cwd, args.projName)
hkl_dump_dir = os.path.join(proj_dir, "hkl_dump")
zarr_dump_dir = os.path.join(proj_dir, "zarr_dump")
ps_dump_dir = os.path.join(proj_dir, "ps_dump")
grain_dump_dir = os.path.join(proj_dir, "grain_dump")

if os.path.exists(proj_dir):
    dir_reset = input(f"Output directory '{proj_dir}' already exists. Continuing will reset the directory. Continue? (y/n): ")
    if dir_reset.lower() in ["y", "yes"]:
        shutil.rmtree(proj_dir)
        os.makedirs(proj_dir)
        os.makedirs(hkl_dump_dir)
        os.makedirs(zarr_dump_dir)
        os.makedirs(ps_dump_dir)
        os.makedirs(grain_dump_dir)
    else:
        print("Exiting program.")
        exit()

else:
    os.makedirs(proj_dir)
    os.makedirs(hkl_dump_dir)
    os.makedirs(zarr_dump_dir)
    os.makedirs(ps_dump_dir)
    os.makedirs(grain_dump_dir)
    
print("Created output directory: ", proj_dir)
print("Created dump directory for mpxrd hkls files: ", hkl_dump_dir)
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

logging.info(f"Reading master grain file {args.grainFile}")
grains_master = pd.read_csv(args.grainFile, header=None) # Read the grain file


det_shape = (powHEDM_sim_config["detector"]["NrPixels"], powHEDM_sim_config["detector"]["NrPixels"])
logging.info(f"Detector shape: {det_shape}\n")

hedm_master = np.zeros(shape=(args.nr_frames,) + det_shape)
pow_master = np.zeros(shape=(args.nr_frames,) + det_shape)
bkg_master = np.zeros(shape=(args.nr_frames,) + det_shape)

# Since multiple frames can be simulated from the same simulation, we need to keep track of which frames are simulated from which simulation
# Create a lookup file that tell the simulation for each frame in the final data
frame_sim_lookup = []

sim_count, frame_count = 0, 0

while frame_count < args.nr_frames:
    sim_count += 1
    logging.info("=="*30)
    logging.info(f"------ Starting SIMULATION {sim_count} ------")

    nr_materials = np.random.randint(powHEDM_sim_config["nr_materials"]["min"],
                                     powHEDM_sim_config["nr_materials"]["max"] + 1)
    
    materials = np.random.choice(mat_pool_for_sim, nr_materials, replace=False).tolist()
    logging.info(f"Sim {sim_count}: Materials selected: {materials}")

    logging.info(f"Sim {sim_count}: Randomly defining experimental parameters...")
    exp_params_frame = {}
    for i_param in powHEDM_sim_config["exp_range"]:
        min, max = powHEDM_sim_config["exp_range"][i_param]["min"], powHEDM_sim_config["exp_range"][i_param]["max"]
        exp_params_frame[i_param] = np.random.uniform(min, max)

    logging.info(f"\t|Sim {sim_count}: Exp params: {exp_params_frame}")

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


    # If HEDM is simulated for any material, then we need to calculate the number of frames that will be used from this simulation
    # This is done by calculating the number of frames that will be used from this simulation and then randomly selecting a number of frames
    # to use from the total number of frames that will be used from this simulation
    # If no HEDM is simulated for any material, then we will use 1 frame from this simulation

    if True in do_hedm_mat_list:
        omega_start_hedm = np.random.uniform(powHEDM_sim_config["hedm_range"]["OmegaStart"]["min"],
                                                     powHEDM_sim_config["hedm_range"]["OmegaStart"]["max"])
            
        omega_end_hedm = np.random.uniform(powHEDM_sim_config["hedm_range"]["OmegaEnd"]["min"],
                                                    powHEDM_sim_config["hedm_range"]["OmegaEnd"]["max"])
        
        omega_step_hedm = np.random.uniform(powHEDM_sim_config["hedm_range"]["OmegaStep"]["min"],
                                                     powHEDM_sim_config["hedm_range"]["OmegaStep"]["max"])
        
        # Calculating the number of frames that will be used from this simulation
        zarr_frames = int((omega_end_hedm - omega_start_hedm) / omega_step_hedm)
        frames_from_this_sim = np.min([args.maxFramesPerSim, np.random.randint(1, zarr_frames+1)])
        frames_from_this_sim = np.min([args.nr_frames - frame_count, frames_from_this_sim])
        frame_idx_start, frame_idx_end = frame_count, frame_count + frames_from_this_sim
    
    else:
        frames_from_this_sim = 1
        frame_idx_start, frame_idx_end = frame_count, frame_count + frames_from_this_sim

        
    if True in do_powder_mat_list:
        omega_start_powder = np.random.uniform(powHEDM_sim_config["powder_range"]["OmegaStart"]["min"],
                                                    powHEDM_sim_config["powder_range"]["OmegaStart"]["max"])
            
        omega_end_powder = np.random.uniform(powHEDM_sim_config["powder_range"]["OmegaEnd"]["min"],
                                                        powHEDM_sim_config["powder_range"]["OmegaEnd"]["max"])
            
        omega_step_powder = np.random.uniform(powHEDM_sim_config["powder_range"]["OmegaStep"]["min"],
                                                        powHEDM_sim_config["powder_range"]["OmegaStep"]["max"]) 
    


    frame_sim_lookup.append([[int(sim_count) for i in range(frames_from_this_sim)]])

    frame_BC_Ypx = int(powHEDM_sim_config["detector"]["NrPixels"]/2) + exp_params_frame["BC_Ypx_offset"]
    frame_BC_Zpx = int(powHEDM_sim_config["detector"]["NrPixels"]/2) + exp_params_frame["BC_Zpx_offset"]

    for (mat_i, mat) in enumerate(materials):

        do_hedm_mat = do_hedm_mat_list[mat_i]
        do_powder_mat = do_powder_mat_list[mat_i]
        
        logging.info(f"\t|Sim {sim_count}: For {mat}, do_powder: {do_powder_mat}, do_hedm: {do_hedm_mat}")

        gauss_width_mat = np.random.uniform(powHEDM_sim_config["mat_specific_range"]["GaussWidth"]["min"],
                                            powHEDM_sim_config["mat_specific_range"]["GaussWidth"]["max"])
        
        logging.info(f"\t|Sim {sim_count}: For {mat}: Running materials project xrdsim...")
        mpxrd_mat = mpxrd(mp_api_key=args.mp_apikey, xray_lambda=exp_params_frame["Wavelength"], mat_formula=mat)
        mpxrd_mat.simulate_xrd()
        mpxrd_mat.create_mp_xrdsim_dict()
        hkls_mat_savepath = os.path.join(hkl_dump_dir, f"sim{sim_count}_{mat}.txt")
        mpxrd_mat.save_hkls_text_midas(hkls_mat_savepath)
        

        #--------------------------------
        # Create lattice parameters sting in right format for ps file
        lattice_params = []
        for i in ["a", "b", "c", "alpha", "beta", "gamma"]:
            lattice_params.append(f"{mpxrd_mat.mp_xrdsim_dict['lattice'][i]} ")
        lattice_params = " ".join(lattice_params)
        
        #--------------------------------
        # Find max nr. of rings that can fit on the detector
        Rmax = np.sqrt(np.max([frame_BC_Ypx, powHEDM_sim_config["detector"]["NrPixels"]-frame_BC_Ypx])**2 +
                        np.max([frame_BC_Zpx, powHEDM_sim_config["detector"]["NrPixels"]-frame_BC_Zpx])**2)
        Rmax = powHEDM_sim_config["detector"]["px"] * Rmax
        TwoThetaMax = np.rad2deg(np.arcsin(Rmax / exp_params_frame["Lsd"]))

        # Find the index of the last ring that is less than TwoThetaMax
        max_ring_idx = np.where(mpxrd_mat.mp_xrdsim_dict["xrd_rings"]["twoTheta"] > TwoThetaMax)[0][0] - 1
        nrRings = mpxrd_mat.mp_xrdsim_dict["xrd_rings"]["ringNr"][max_ring_idx]
        

        if do_hedm_mat:
            logging.info(f"\t|Sim {sim_count}: For {mat}: Running HEDM simulation...")

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
            grains_hedm_mat_savepath = os.path.join(grain_dump_dir, f"sim{sim_count}_hedm_{mat}.csv")
            grains_hedm_mat.to_csv(grains_hedm_mat_savepath, header=False, index=False)


            ps_hedm_mat_filepath = os.path.join(ps_dump_dir, f"sim{sim_count}_hedm_{mat}.txt")
            zarr_hedm_mat_filename = os.path.join(zarr_dump_dir, f"sim{sim_count}_hedm_{mat}")

            # =================== WRITE HEADER FILE ===================
            with open(ps_hedm_mat_filepath, 'w') as f:
                f.write("# Parameter file for forward simulation\n")
                f.write(f"InFileName {grains_hedm_mat_savepath}\n")
                f.write(f"IntensitiesFile {hkls_mat_savepath}\n")
                f.write(f"OutFileName {zarr_hedm_mat_filename}\n")
                f.write(f"LatticeConstant {lattice_params}\n")
                f.write(f"SpaceGroup {mpxrd_mat.mp_xrdsim_dict['structure']['SpaceGroup']}\n")
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

            zarr_hedm_mat_filepath = os.path.join(zarr_dump_dir, f"sim{sim_count}_hedm_{mat}_scanNr_0.zip")
            mat_hedm_arr = zarr.open(zarr_hedm_mat_filepath, "r")["exchange"]["data"]
            rand_idx = list(np.random.randint(0, len(mat_hedm_arr), frames_from_this_sim))

            hedm_master[frame_idx_start:frame_idx_end, :, :] += mat_hedm_arr[rand_idx, :, :]

            logging.info(f"\t|Sim {sim_count}: For {mat}: HEDM simulation complete.")

            if (args.clearGrainDump).lower() in ["y", "yes", "true"]:
                if os.path.exists(grains_hedm_mat_savepath): os.remove(grains_hedm_mat_savepath)

            if (args.clearPSDump).lower() in ["y", "yes", "true"]:
                if os.path.exists(ps_hedm_mat_filepath): os.remove(ps_hedm_mat_filepath)

            if (args.clearZarrDump).lower() in ["y", "yes", "true"]:
                if os.path.exists(zarr_hedm_mat_filepath): os.remove(zarr_hedm_mat_filepath)

        if do_powder_mat:

            logging.info(f"\t|Sim {sim_count}: For {mat}: Running powder simulation...")
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
            grains_powder_mat_savepath = os.path.join(grain_dump_dir, f"sim{sim_count}_powder_{mat}.csv")
            grains_powder_mat.to_csv(grains_powder_mat_savepath, header=False, index=False)


            ps_powder_mat_filepath = os.path.join(ps_dump_dir, f"sim{sim_count}_powder_{mat}.txt")
            zarr_powder_mat_filename = os.path.join(zarr_dump_dir, f"sim{sim_count}_powder_{mat}")

            # =================== WRITE HEADER FILE ===================
            with open(ps_powder_mat_filepath, 'w') as f:
                f.write("# Parameter file for forward simulation\n")
                f.write(f"InFileName {grains_powder_mat_savepath}\n")
                f.write(f"IntensitiesFile {hkls_mat_savepath}\n")
                f.write(f"OutFileName {zarr_powder_mat_filename}\n")
                f.write(f"LatticeConstant {lattice_params}\n")
                f.write(f"SpaceGroup {mpxrd_mat.mp_xrdsim_dict['structure']['SpaceGroup']}\n")
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

            zarr_powder_mat_filepath = os.path.join(zarr_dump_dir, f"sim{sim_count}_powder_{mat}_scanNr_0.zip")
            mat_powder_arr = zarr.open(zarr_powder_mat_filepath, "r")["exchange"]["data"]
            mat_powder_arr = np.sum(mat_powder_arr, axis=0)
            mat_powder_arr = mat_powder_arr * max_int_powder_mat / np.max(mat_powder_arr)

            pow_master[frame_idx_start:frame_idx_end, :, :] += mat_powder_arr
            logging.info(f"\t|Sim {sim_count}: For {mat}: Powder simulation complete.")


            # Clear the dump files for the powder simulation
            if (args.clearGrainDump).lower() in ["y", "yes", "true"]:
                if os.path.exists(grains_powder_mat_savepath): os.remove(grains_powder_mat_savepath)

            if (args.clearPSDump).lower() in ["y", "yes", "true"]:
                if os.path.exists(ps_powder_mat_filepath): os.remove(ps_powder_mat_filepath)

            if (args.clearZarrDump).lower() in ["y", "yes", "true"]:
                if os.path.exists(zarr_powder_mat_filepath): os.remove(zarr_powder_mat_filepath)

        # Clear the hkl dump file for material
        if (args.clearHKLDump).lower() in ["y", "yes", "true"]:
            if os.path.exists(hkls_mat_savepath): os.remove(hkls_mat_savepath)


    if do_bkg:
        logging.info(f"\t|Sim {sim_count}: Running background simulation...")

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
            thomson_ff_decay_list = np.random.uniform(powHEDM_sim_config["bkg"]["thomson"]["form_factor_decay"]["min"],
                                            powHEDM_sim_config["bkg"]["thomson"]["form_factor_decay"]["max"],
                                            frames_from_this_sim)
            thomson_scale_Imax_list = np.random.uniform(powHEDM_sim_config["bkg"]["thomson"]["scale_Imax"]["min"],
                                                powHEDM_sim_config["bkg"]["thomson"]["scale_Imax"]["max"],
                                                frames_from_this_sim)
            
            logging.info(f"\t|Sim {sim_count}: Including thomson scattering.")

        else:
            thomson_ff_decay_list = None
            thomson_scale_Imax_list = None

        if "compton" in bkg_types:
            compton_scale_Imax_list = np.random.uniform(powHEDM_sim_config["bkg"]["compton"]["scale_Imax"]["min"],
                                                        powHEDM_sim_config["bkg"]["compton"]["scale_Imax"]["max"],
                                                        frames_from_this_sim)
            logging.info(f"\t|Sim {sim_count}: Including compton scattering.")
        else:
            compton_scale_Imax_list = None


        if "poisson" in bkg_types:
            logging.info(f"\t|Sim {sim_count}: Including poisson noise.")

        if "gaussian" in bkg_types:
            gaussian_scale_Imax_list = np.random.uniform(powHEDM_sim_config["bkg"]["gaussian"]["scale_Imax"]["min"],
                                                        powHEDM_sim_config["bkg"]["gaussian"]["scale_Imax"]["max"],
                                                        frames_from_this_sim)
            gaussian_loc_list = np.random.uniform(powHEDM_sim_config["bkg"]["gaussian"]["loc"]["min"],
                                                        powHEDM_sim_config["bkg"]["gaussian"]["loc"]["max"],
                                                        frames_from_this_sim)
            gaussian_scale_list = np.random.uniform(powHEDM_sim_config["bkg"]["gaussian"]["scale"]["min"],
                                                        powHEDM_sim_config["bkg"]["gaussian"]["scale"]["max"],
                                                        frames_from_this_sim)
            logging.info(f"\t|Sim {sim_count}: Including gaussian noise.")
        else:
            gaussian_scale_Imax_list = None
            gaussian_loc_list = None
            gaussian_scale_list = None
        
        ps_bkg_filepath = os.path.join(ps_dump_dir, f"sim{sim_count}_bkg.txt")

        # =================== WRITE HEADER FILE ===================
        with open(ps_bkg_filepath, 'w') as f:
            f.write("# Parameter file for background simulation\n")
            f.write(f"thomson_ff_decay_list {thomson_ff_decay_list}\n")
            f.write(f"thomson_scale_Imax_list {thomson_scale_Imax_list}\n")
            f.write(f"compton_scale_Imax_list {compton_scale_Imax_list}\n")
            f.write(f"gaussian_scale_Imax_list {gaussian_scale_Imax_list}\n")
            f.write(f"gaussian_loc_list {gaussian_loc_list}\n")
            f.write(f"gaussian_scale_list {gaussian_scale_list}\n")
            f.write(f"beamstop_radius {exp_params_frame['beamstop_radius']}\n")
            f.write(f"beamstop_transmission_frac {exp_params_frame['beamstop_transmission_frac']}\n")

        for (i, frame_i) in enumerate(range(frame_idx_start, frame_idx_end)):

            bkg_sim_i = copy.deepcopy(bkg)
            if "thomson" in bkg_types:
                bkg_sim_i.thomson_scattering(form_factor_decay=thomson_ff_decay_list[i],
                                            scale_Imax=thomson_scale_Imax_list[i])
                
                bkg_master[frame_i, :, :] += bkg_sim_i.thomson_bkg["intensity"]
                
    
            if "compton" in bkg_types:
                bkg_sim_i.compton_scattering(scale_Imax=compton_scale_Imax_list[i])
                bkg_master[frame_i, :, :] += bkg_sim_i.compton_bkg["intensity"]

            if "poisson" in bkg_types:
                frame_img_arr = hedm_master[frame_i] + pow_master[frame_i]
                bkg_sim_i.poisson_noise(intensity_arr=frame_img_arr)
                bkg_master[frame_i, :, :] += bkg_sim_i.poisson_noise["intensity"]

            if "gaussian" in bkg_types:
                bkg_sim_i.gaussian_noise(loc=gaussian_loc_list[i],
                                        scale=gaussian_scale_list[i],
                                        scale_Imax=gaussian_scale_Imax_list[i])
                
                bkg_master[frame_i, :, :] += bkg_sim_i.gaussian_noise["intensity"]

        if (args.clearPSDump).lower() in ["y", "yes", "true"]:
            if os.path.exists(ps_bkg_filepath): os.remove(ps_bkg_filepath)

    frame_count += frames_from_this_sim

logging.info("\n\n")
logging.info("="*30)
logging.info("SIMULATION COMPLETE.")
logging.info("="*30)

if (args.saveHdf5).lower() in ["y", "yes", "true"]:
    sim_out_hdf5_filepath = os.path.join(proj_dir, f"powHEDM_frames.h5")
    logging.info(f"Saving frames to hdf5 file at '{sim_out_hdf5_filepath}'")

    with h5py.File(sim_out_hdf5_filepath, "w") as h5file:
        #Saving the HEDM, powder and background frames
        h5file.create_dataset("HEDM", data=hedm_master)
        h5file.create_dataset("Powder", data=pow_master)
        h5file.create_dataset("Background", data=bkg_master)        
        h5file.create_dataset("frame_sim_lookup", data=frame_sim_lookup)

if (args.saveNpz).lower() in ["y", "yes", "true"]:
    sim_out_npz_filepath = os.path.join(proj_dir, f"powHEDM_frames.npz")
    logging.info(f"Saving frames to npz file at '{sim_out_npz_filepath}'")
    np.savez_compressed(sim_out_npz_filepath, 
                        hedm=hedm_master, 
                        powder=pow_master, 
                        background=bkg_master)

if os.path.exists("SpotMatrixGen.csv"):
    os.remove("SpotMatrixGen.csv")
if os.path.exists("GrainsGen.csv"):
    os.remove("GrainsGen.csv")
if os.path.exists("hkls.csv"):
    os.remove("hkls.csv")