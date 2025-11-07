import subprocess

import os
SEP = os.sep

class midas_fwsim:
    """Class for running the MIDAS forward simulation..
    """

    def __init__(self, ps_filepath, MIDAS_INSTALL_PATH, out_dir, ncpu=8):
        """Initialize the hedm_fwsim object.
        Args:
            ps_filepath (str): Path to the parameter script txt file.
            MIDAS_INSTALL_PATH (str): Path to the MIDAS installation.
            out_dir (str): Path to the output directory where the simulation results will be saved.
            ncpu (int): Number of CPUs to use for the simulation.
        """

        self.ps_filepath = ps_filepath
        self.MIDAS_INSTALL_PATH = MIDAS_INSTALL_PATH
        self.out_dir = out_dir
        self.ncpu = ncpu


    def run_midas_forward_sim(self):
        cmd = [os.path.join(self.MIDAS_INSTALL_PATH, f"FF_HEDM{SEP}bin{SEP}ForwardSimulationCompressed"), self.ps_filepath, str(self.ncpu)]
        cmd = " ".join(cmd)
        print("Running command: ", cmd)
        outfile_path = os.path.join(self.out_dir, "hedm_fwsim_out.txt")
        errfile_path = os.path.join(self.out_dir, "hedm_fwsim_err.txt")
        with open(outfile_path, 'w') as outfile, open(errfile_path, 'w') as errfile:
            return_code = subprocess.call(cmd, shell=True, stdout=outfile, stderr=errfile)
        
        