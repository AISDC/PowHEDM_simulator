from mpxrd_sim import mpxrd
import os, argparse, pathlib


mat_list = ["SiO2", "Li2O", "Al2O3", "TiO2", "CeO2", "ZrO2", "NiO", "Fe2O3",\
            "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Zr", "Nb", "Mo", \
            "Pd", "Ag", "Ta", "W", "Re", "Pt", "Au"]

parser = argparse.ArgumentParser(description="Create mpxrd cache")

parser.add_argument("-mp_apikey",
                    type=str,
                    help="Materials Project API key")

parser.add_argument("-cacheDir",
                    default="mpxrd_cache",
                    type=str,
                    help="Directory name for saving the cache")

parser.add_argument("-xrayLambda",
                    default=0.5,
                    type=float,
                    help="X-ray wavelength in Angstrom")

args, unparsed = parser.parse_known_args()

if len(unparsed) > 0:
    print('Unrecognized argument(s): \n%s \nProgram exiting ... ... ' % '\n'.join(unparsed))
    exit(0)

cwd = pathlib.Path(os.getcwd())
cache_dir = os.path.join(cwd, args.cacheDir)
xray_lambda = args.xrayLambda # Angstrom

if not os.path.exists(cache_dir):
    os.makedirs(cache_dir)

for mat in mat_list:
    print("-"*30)
    print(f"Creating cache for {mat}...", flush=True)
    mpxrd_mat = mpxrd(mp_api_key=args.mp_apikey, xray_lambda=args.xrayLambda, mat_formula=mat)
    mpxrd_mat.simulate_xrd()
    mpxrd_mat.create_mp_xrdsim_dict()
    hkls_mat_savepath = os.path.join(cache_dir, f"{mat}_{mpxrd_mat.mp_matid}_hkls.txt")
    mpxrd_mat.save_hkls_text_midas(hkls_mat_savepath)

    sim_info_mat_savepath = os.path.join(cache_dir, f"{mat}_{mpxrd_mat.mp_matid}_info.txt")
    mpxrd_mat.save_sim_info_text(sim_info_mat_savepath)