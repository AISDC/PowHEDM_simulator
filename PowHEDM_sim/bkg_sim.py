import numpy as np

class bkg_sim:
    """Class for running the background simulation.
    """

    def __init__(self, BC_Ypx, BC_Zpx, xray_lambda, Lsd, px_size, nr_pixels,
                 beamstop_radius=0, beamstop_transmission_frac=0.1):
        
        """Initialize the bkg_sim object.
        """

        self.BC_Ypx = BC_Ypx  # pixels
        self.BC_Zpx = BC_Zpx  # pixels
        self.xray_lambda = xray_lambda  # Angstrom
        self.Lsd = Lsd  # microns
        self.px_size = px_size  # microns
        self.nr_pixels = nr_pixels  # number
        self.beamstop_radius = beamstop_radius  # pixels
        self.beamstop_transmission_frac = beamstop_transmission_frac  # fraction between 0-1
        
        # Initialize common constants
        self.CONSTANTS = {
            "AIR_DENSITY": 1.225,  # kg/m^3 at STP
            "CLASSICAL_ELECTRON_RADIUS": 2.818e-15,  # in meters
            "AVOGADRO_NUMBER": 6.022e23,
            "AIR_MOLAR_MASS": 28.97e-3,  # kg/mol
            "ELECTRONS_PER_AIR_MOLECULE": 14.42,  # weighted average (N2: 14, O2: 16, N2/O2 ratio: 79/21)
            "AIR_PATH_LENGTH": 1.0,  # meters - path length through air
            "H_PLANCK": 6.626e-34,  # J⋅s
            "C_LIGHT": 3e8,  # m/s
            "ELECTRON_REST_ENERGY": 511.0,  # keV
            "AIR_AVERAGE_Z": 7.2,  # Average Z for air (N2, O2)
        }
        
        # Calculate common parameters
        self.Qmap = self.Q_map()


    def YZ_map(self):
        # Create pixel coordinate arrays relative to beam center
        Y = np.arange(self.nr_pixels) - self.BC_Ypx  # Horizontal coordinate
        Z = np.arange(self.nr_pixels) - self.BC_Zpx  # Vertical coordinate
        YY, ZZ = np.meshgrid(Y, Z)  # Y varies along columns, Z varies along rows
        return YY, ZZ
    

    def R_map(self):
        YY, ZZ = self.YZ_map()
        R_map= np.sqrt(YY**2 + ZZ**2)  # physical distance from beam center; pixels
        return R_map


    def scattering_angle(self):
        R_map_microns = self.R_map() * self.px_size # in microns
        scattering_angle = 2 * np.arctan(R_map_microns / self.Lsd)  # scattering angle; radians
        return scattering_angle


    def Q_map(self):
        """Calculate parameters that are common to both Thomson and Compton scattering."""

        scattering_angle = self.scattering_angle()
        scattering_angle = np.maximum(scattering_angle, 1e-6)  # Avoid division by zero
        
        # Momentum transfer
        Qmap = 4 * np.pi * np.sin(scattering_angle / 2) / (self.xray_lambda * 1e-10)  # m^-1
        return Qmap


    def beamstop_factor(self):
        """Create a beamstop transmission factor array that maps the transmission ration at each detector pixel.
        """

        R_map = self.R_map() # in pixels
        beamstop_mask = R_map < self.beamstop_radius
        beamstop_factor = np.where(beamstop_mask, self.beamstop_transmission_frac, 1.0)
        return beamstop_factor


    def thomson_scattering(self, form_factor_decay, scale_Imax=None):
        """
        Simulate Thomson (coherent elastic) scattering.
        
        Args: 
            form_factor_decay (float): Form factor decay constant.
            scale_Imax (float): Maximum intensity of the background.
        """
                
        # Calculate electron density in air [no. of electrons per m3]
        electron_density_air = (self.CONSTANTS["AIR_DENSITY"] / self.CONSTANTS["AIR_MOLAR_MASS"]) *\
            self.CONSTANTS["AVOGADRO_NUMBER"] * self.CONSTANTS["ELECTRONS_PER_AIR_MOLECULE"] # nr. electrons . m^-3
        
        # Thomson cross-section
        thomson_cross_section = self.CONSTANTS["CLASSICAL_ELECTRON_RADIUS"]**2  # m^2
        
        # Form factor for Thomson scattering (accounts for atomic binding)
        # Note: Using q in m^-1, so we need appropriate scaling
        form_factor_air = np.exp(-form_factor_decay * (self.Qmap * 1e-10)**2)  # Dimensionless
        
        # Thomson intensity per electron (CORRECTED: form factor NOT squared)
        thomson_intensity_per_electron = thomson_cross_section * form_factor_air  # m^2
        
        # Intensity per unit volume
        thomson_intensity_per_volume = thomson_intensity_per_electron * electron_density_air  # m^-1
        
        # Total intensity along the beam path
        thomson_intensity = thomson_intensity_per_volume * self.CONSTANTS["AIR_PATH_LENGTH"]  # Dimensionless

        thomson_intensity *= self.beamstop_factor() 
        
        # Scale to desired maximum intensity
        if scale_Imax not in [None, "None", "none"]:
            thomson_intensity *= scale_Imax / np.max(thomson_intensity)

        # Store results
        self.thomson_bkg = {
            "form_factor_decay": form_factor_decay,
            "scale_Imax": scale_Imax,
            "intensity": thomson_intensity,
            "cross_section": thomson_cross_section,
            "form_factor": form_factor_air
        }


    def compton_scattering(self, scale_Imax=None):
        """
        Simulate Compton (incoherent inelastic) scattering using Klein-Nishina formula.
        Args: 
            scale_Imax (float): Maximum intensity of the background.
        """
        
        # Calculate electron density in air [no. of electrons per m3]
        electron_density_air = (self.CONSTANTS["AIR_DENSITY"] / self.CONSTANTS["AIR_MOLAR_MASS"]) *\
            self.CONSTANTS["AVOGADRO_NUMBER"] * self.CONSTANTS["ELECTRONS_PER_AIR_MOLECULE"] # nr. electrons . m^-3
        
        # Calculate incident photon energy
        incident_energy_J = self.CONSTANTS["H_PLANCK"] * self.CONSTANTS["C_LIGHT"] / (self.xray_lambda * 1e-10)
        incident_energy_keV = incident_energy_J / (1.602e-16)  # Convert to keV
        alpha = incident_energy_keV / self.CONSTANTS["ELECTRON_REST_ENERGY"]  # Dimensionless energy

        # Thomson cross-section (base for Klein-Nishina)
        thomson_cross_section = self.CONSTANTS["CLASSICAL_ELECTRON_RADIUS"]**2  # m^2
        
        # Energy of scattered photon (Compton formula)
        cos_theta = np.cos(self.scattering_angle())
        scattered_energy_ratio = 1.0 / (1.0 + alpha * (1.0 - cos_theta))
        
        # Klein-Nishina differential cross-section
        # dσ/dΩ = (r_e²/2) × (E'/E)² × [E/E' + E'/E - sin²(θ)]
        klein_nishina_factor = (scattered_energy_ratio**2 * 
                               (scattered_energy_ratio + 1.0/scattered_energy_ratio - 
                                np.sin(self.scattering_angle())**2))
        
        # Compton cross-section per electron
        compton_cross_section_per_electron = (thomson_cross_section / 2.0) * klein_nishina_factor  # m^2
        
        # Account for atomic number dependence (Compton ∝ Z, not Z²)
        compton_intensity_per_volume = (compton_cross_section_per_electron * electron_density_air)  # m^-1
        
        # Total intensity along the beam path
        compton_intensity = compton_intensity_per_volume * self.CONSTANTS["AIR_PATH_LENGTH"]  # Dimensionless
        
        # Scale to desired maximum intensity
        if scale_Imax not in [None, "None", "none"]:
            compton_intensity *= scale_Imax / np.max(compton_intensity)
        
        # Store results
        self.compton_bkg = {
            "scale_Imax": scale_Imax,
            "intensity": compton_intensity,
            "scattered_energy_ratio": scattered_energy_ratio,
            "klein_nishina_factor": klein_nishina_factor,
            "cross_section": compton_cross_section_per_electron
        }


    def poisson_noise(self, intensity_arr, scale_Imax=None):
        """Create Poisson noise map over the detector.
        Shot noise: Poisson noise is due to the discrete nature of photons.
        """

        poisson_noise = np.random.poisson(intensity_arr) - intensity_arr

        if scale_Imax not in [None, "None", "none"]:
            poisson_noise *= scale_Imax / np.max(np.abs(poisson_noise))

        self.poisson_noise = {
            "intensity": poisson_noise,
        }


    def gaussian_noise(self, loc=10, scale=3, scale_Imax=None):
        """Create Gaussian noise map over the detector.
        Readout noise from electronics
        """

        gaussian_noise = np.random.normal(loc, scale, (self.nr_pixels, self.nr_pixels))
        if scale_Imax not in [None, "None", "none"]:
            gaussian_noise *= scale_Imax / np.max(gaussian_noise)

        self.gaussian_noise = {
            "intensity": gaussian_noise,
            "loc": loc,
            "scale": scale,
            "scale_Imax": scale_Imax
        }
