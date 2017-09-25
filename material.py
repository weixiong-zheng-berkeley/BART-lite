import numpy as np
from math import pi

class material(object):
    def __init__(self):
        '''@brief constructor of material class

        @param self The reference to current instance of the material class
        '''
        self.n_materials = 1#change it
        self.n_group = 1#change it
        self.g_thermal = 0#from which we can see upscattering effect
        # The material properties are stored as dictionaries. The keys are material
        # IDs. Dictionary values are numpy arrays of properties for sigt, nu, sigf
        # and
        # The following is the basic material properties
        self.sigt = dict()
        self.nu = dict()
        self.sigf = dict()
        self.sigs = dict()
        self.chi = dict()

        # The following is the derived properties based on basic properties
        self.sigs_per_str = dict()
        self.chi_nu_sigf = dict()
        self.chi_nu_sigf_per_str = dict()
        self.diff_coef = dict()
        self.inv_sigt = dict()

        # The following is for upscattering acceleration
        self.ksi = dict()
        # TODO: put whatever else necessary parameters if needed

    def read_xsec(self,xsec_filename):
        """@brief read cross sections from xml file

        @param xsec_filename The xml filename
        """
        # TODO: fill in material properties reading. Please read in the following:
        # self.sigt; self.nu; self.sigf; self.sigs
        # data structures are dictionaries of numpy arraies

    def read_material_id(self,id_filename):
        """@brief read mateiral id

        @param id_filename The xml filename
        """
        # TODO: fill in ID reading

    def derived_properties(self):
        '''@brief derive properties after reading in  basic properties
        '''
        # inv_sigt and diff_coef:
        for k, v in self.sigt.items():
            self.inv_sigt[k] = 1. / self.sigt[k]
            self.diff_coef[k] = 1./ (3. * self.sigt[k])
        # sigs_per_str:
        for k, v in self.sigs.items():
            self.sigs_per_str[k] = v / (4.0 * pi)
        # nu_sigf and nu_sigf_per_str:
        for k in range(self.n_materials):
            self.chi_nu_sigf[k] = np.outer(self.chi[k], self.nu[k]*self.sigf[k])
            self.chi_nu_sigf_per_str[k] = self.chi_nu_sigf[k] / (4.0 * pi)

    def derive_scattering_eigenvalue(self):
        for k, v in self.sigs.items():
            # slicing to get thermal group transfer matrices
            thermal_scat = v[self.g_thermal::1, self.g_thermal::1]
            self.ksi[k] = np.linalg.eigvals(thermal_scat)
