
import numpy as np #type:ignore

# Import libraries needed for materials project API
from mp_api.client import MPRester #type:ignore
from pymatgen.analysis.diffraction.xrd import XRDCalculator #type:ignore
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer #type:ignore


class mpxrd:
    """Class for simulating XRD patterns of materials using the Materials Project API.
    """

    def __init__(self, mp_api_key=None, mp_matid=None, mat_formula=None, xray_lambda=None, verbose=0):
        self.mp_api_key = mp_api_key
        self.mp_matid = mp_matid
        self.mat_formula = mat_formula
        self.xray_lambda = xray_lambda
        self.verbose = verbose

        if self.verbose > 0:
            print("'mpxrd' object created.", flush=True)
        
        if self.mat_formula is not None and self.mp_matid is None:
            self.assign_mp_matid_from_formula()

        self.check_input_params()


    #--------------------------------
    def check_input_params(self):
        """Check the input parameters for the 'mpxrd' object.
        """
        if self.verbose > 0:
            print("-"*50)
            print("Checking input parameters...", flush=True)
            
        nr_warnings = 0
        if self.mp_api_key is None:
            print("\tWARNING: 'mp_api_key' is not provided."
                  " Please provide a Materials Project API key for the simulation.", flush=True)
            nr_warnings += 1
        if self.mp_matid is None and self.mat_formula is None:
            print("\tWARNING: Both 'mp_matid' and 'mat_formula' are not provided."
                  " Please provide either a Materials Project ID or a chemical formula for the material you want to simulate.", flush=True)
            nr_warnings += 1
        if self.xray_lambda is None:
            print("\tWARNING: 'xray_lambda' is not provided."
                  " Please provide a wavelength (in Angstroms) for the X-ray radiation you want to use for the simulation.", flush=True)
            nr_warnings += 1
        if self.xray_lambda is not None and (self.xray_lambda < 0.06 or self.xray_lambda > 100):
            print(f"\tWARNING: 'xray_lambda={self.xray_lambda}' is not within the valid range of 0.06 to 100 Angstroms."
                  " Please provide a valid wavelength.", flush=True)
            nr_warnings += 1

        if self.verbose > 0:
            if nr_warnings == 0:
                print("All input parameters are valid.", flush=True)
            else:
                print(f"FOUND {nr_warnings} WARNING(S). CHECK BEFORE RUNNING SIMULATION.", flush=True)


    #--------------------------------
    def assign_mp_matid_from_formula(self):
        """Assign the Materials Project ID when material is provided as a chemical formula.
        """

        with MPRester(api_key=self.mp_api_key) as mpr:
            stable_materials = mpr.materials.summary.search(
                energy_above_hull=(0, 1E-10),
                formula=self.mat_formula,
                fields=["material_id", "formula_pretty", "energy_above_hull"]
            )
        
        if len(stable_materials) == 0:
            raise ValueError(f"No stable materials found for {self.mat_formula}."
                             " Please check the chemical formula and try again OR provide a Materials Project ID instead.")
        else:
            self.mp_matid = str(stable_materials[0].material_id[:])
            if self.verbose > 0:
                print(f"Materials Project ID for {self.mat_formula} is set as {self.mp_matid}")


    #--------------------------------
    def simulate_xrd(self):
        """Simulate the XRD pattern of the material.
        """

        if self.verbose > 0:
            print("-"*50)
            print(f"Simulating XRD pattern for {self.mp_matid} ({self.mat_formula})...", flush=True)

        with MPRester(api_key=self.mp_api_key) as mpr:
            # Retrieve the structure of the material
            self.structure = mpr.get_structure_by_material_id(self.mp_matid)

            if self.verbose > 0:
                print("Structure retrieved successfully.", flush=True)

        # important to use the conventional structure to ensure
        # that peaks are labelled with the conventional Miller indices
        sga = SpacegroupAnalyzer(self.structure)
        self.conventional_structure = sga.get_conventional_standard_structure()
        
        if self.verbose > 0:
            print("Conventional structure retrieved successfully.", flush=True)

        self.spacegroup = sga.get_space_group_number()
        self.crystal_system = sga.get_crystal_system()

        calculator = XRDCalculator(wavelength=self.xray_lambda)
        self.xrd_pattern = calculator.get_pattern(self.conventional_structure)

        if self.verbose > 0:
            print("XRD pattern retrieved successfully.", flush=True)


    #--------------------------------
    def create_mp_xrdsim_dict(self):
        """
        Create a dictionary containing the XRD simulation results for a material.
        Args:
            mp_matid (str): Materials Project material ID (e.g., "mp-150" for Fe-fcc)
            xray_wavelength (float): X-ray wavelength in Angstroms
            mp_api_key (str): Materials Project API key
            mp_matname (str; default=None): Name of the material to add an optional label for the material
        Returns:
            mp_xrdsim_dict (dict): dictionary containing the XRD simulation results
        """

        if self.verbose > 0:
            print("-"*50)
            print("Creating XRD simulation results in dictionary format...", flush=True)
        
        if "xrd_pattern" not in vars(self).keys():
            if self.verbose > 0:
                print("XRD pattern not found. Calling 'simulate_xrd()' to generate XRD pattern...", flush=True)
            self.simulate_xrd()
            

        mp_xrdsim_dict = {
            "mp_matid": self.mp_matid,
            "mp_matname": self.mat_formula,
            "structure": {
                "SpaceGroup": self.spacegroup,
                "CrystalSystem": self.crystal_system,
                },
            "lattice": self.conventional_structure._lattice.params_dict,
            "xrd_rings": {
                "xrd_wavelength": self.xray_lambda,
                "ringNr": [i+1 for i in range(len(self.xrd_pattern.x))],
                "twoTheta": list(self.xrd_pattern.x),
                "intensity": list(self.xrd_pattern.y / np.max(self.xrd_pattern.y)),
                "d_hkls": list(self.xrd_pattern.d_hkls)
            }
        }

        rings_hkls = [[ring_i_hkl["hkl"] for ring_i_hkl in self.xrd_pattern.hkls[ring_i]] 
                    for ring_i in range(len(self.xrd_pattern.x))]

        rings_hkls_multiplicity = [[ring_i_hkl["multiplicity"] for ring_i_hkl in self.xrd_pattern.hkls[ring_i]] 
                                for ring_i in range(len(self.xrd_pattern.x))]

        mp_xrdsim_dict["xrd_rings"]["hkls"] = rings_hkls
        mp_xrdsim_dict["xrd_rings"]["hkls_multiplicity"] = rings_hkls_multiplicity


        hkls, hkls_Int, d_hkls, hkls_ringNr = [], [], [], [] # store hkls and their relative intensities

        # List of total multiplicity of all hkls in the ring
        rings_multiplicity = [np.sum(rings_hkls_multiplicity[i]) for i in range(0, len(rings_hkls_multiplicity))]

        # Iterate over all rings
        for ring_i in range(0, len(self.xrd_pattern.x)):
            ringInt = self.xrd_pattern.y[ring_i] # intensity of the ring
            ringM = rings_multiplicity[ring_i] # total multiplicity of all hkls in the ring
            ring_d_hkl = self.xrd_pattern.d_hkls[ring_i] # d-spacings of all hkls in the ring
            ringNr = ring_i + 1 # Miller indices of all hkls in the ring
            
            # Iterate over all hkls in the ring
            for (hkl_i, hkl) in enumerate(rings_hkls[ring_i]):
                hkl_M = rings_hkls_multiplicity[ring_i][hkl_i] # multiplicity of the hkl
                hkl_Int = ringInt * hkl_M / ringM # relative intensity of the hkl
                hkls_Int.append(hkl_Int)
                hkls.append(hkl)
                d_hkls.append(ring_d_hkl)
                hkls_ringNr.append(ringNr)

        mp_xrdsim_dict["xrd_hkls"] = {}
        mp_xrdsim_dict["xrd_hkls"]["xrd_wavelength"] = self.xray_lambda
        mp_xrdsim_dict["xrd_hkls"]["hkls"] = hkls
        mp_xrdsim_dict["xrd_hkls"]["hkls_intensity"] = list(np.array(hkls_Int) / np.max(np.array(hkls_Int)))
        mp_xrdsim_dict["xrd_hkls"]["d_hkls"] = d_hkls
        mp_xrdsim_dict["xrd_hkls"]["ringNr"] = hkls_ringNr

        self.mp_xrdsim_dict = mp_xrdsim_dict
        if self.verbose > 0:
            print("DONE.", flush=True)


    #--------------------------------
    def save_hkls_text_midas(self, hkls_savepath):
        """Save the hkls and their intensities to a text file in MIDAS format.
        This text file can then be used as input for the MIDAS forward simulation code.
        """

        # =================== WRITE TABLE FILE ===================
        with open(hkls_savepath, 'w') as f:
            # Write table header
            f.write("h k l SpotInt D-spacing RingNr\n")
            
            # Write HKL data rows (space-separated)
            xrd_hkls = self.mp_xrdsim_dict['xrd_hkls']
            for i in range(len(xrd_hkls['hkls'])):
                hkl = xrd_hkls['hkls'][i]
                intensity = xrd_hkls['hkls_intensity'][i]
                d_spacing = xrd_hkls['d_hkls'][i]
                ring_nr = xrd_hkls['ringNr'][i]
                
                f.write(f"{hkl[0]} {hkl[1]} {hkl[2]} {intensity:.10f} {d_spacing:.10f} {ring_nr}\n")


    #--------------------------------
    def save_sim_info_text(self, info_savepath):
        """Save the information about the XRD simulation to a text file.
        """

        with open(info_savepath, 'w') as f:
            f.write("# XRD Data File - Header Information\n")
            f.write(f"MATERIAL_CODE {self.mp_xrdsim_dict['mp_matid']}\n")
            f.write(f"MATERIAL_NAME {self.mp_xrdsim_dict['mp_matname']}\n")
            f.write(f"SPACE_GROUP {self.mp_xrdsim_dict['structure']['SpaceGroup']}\n")
            f.write(f"CRYSTAL_SYSTEM {self.mp_xrdsim_dict['structure']['CrystalSystem']}\n")
            
            # Write lattice parameters
            f.write(f"LATTICE_A {self.mp_xrdsim_dict['lattice']['a']}\n")
            f.write(f"LATTICE_B {self.mp_xrdsim_dict['lattice']['b']}\n")
            f.write(f"LATTICE_C {self.mp_xrdsim_dict['lattice']['c']}\n")
            f.write(f"LATTICE_ALPHA {self.mp_xrdsim_dict['lattice']['alpha']}\n")
            f.write(f"LATTICE_BETA {self.mp_xrdsim_dict['lattice']['beta']}\n")
            f.write(f"LATTICE_GAMMA {self.mp_xrdsim_dict['lattice']['gamma']}\n")
            
            # Write XRD parameters
            f.write(f"XRD_WAVELENGTH {self.mp_xrdsim_dict['xrd_hkls']['xrd_wavelength']}\n")
            f.write(f"NUM_HKLS {len(self.mp_xrdsim_dict['xrd_hkls']['hkls'])}\n")
