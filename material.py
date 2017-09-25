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
        # spectrum
        self.ksi_ua = dict()
        # sigt one group
        self.sigt_ua = dict()
        # diff_coef one group
        self.diff_coef_ua = dict()
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
            thermal_scat = v[self.g_thermal:, self.g_thermal:]
            self.ksi_ua[k] = np.linalg.eigvals(thermal_scat)

    def produce_acceleration_properties(self):
        for k in xrange(n_materials):
            self.sigt_ua[k] = 0.0
            self.diff_coef_ua[k] = 0.0
            self.sigt_ua[k] = np.dot(self.ksi_ua[k], self.sigt[k][self.g_thermal:])
            self.diff_coef_ua[k] = np.dot(self.ksi_ua[k], self.diff_coef[k][self.g_thermal:])

    def estimate_cell_correction_ua(self, cell_corrections_at_qp, material_id):
        '''@brief function used to estimate one-group vector_D at quadrature points

        @param cell_corrections_at_qp vector_D for all groups at the specific quadrature
        point. Use numpy.array
        @param material_id Materail ID for current cell
        @return a one-group vector_D in forms of numpy array
        '''
        vector_D = np.array([0.0, 0.0])
        # x component of vector_D
        vector_D[0] = np.dot(cell_corrections_at_qp[0][self.g_thermal:], self.ksi_ua[material_id])
        # y component of vector_D
        vector_D[1] = np.dot(cell_corrections_at_qp[1][self.g_thermal:], self.ksi_ua[material_id])
        return vector_D

    def estimate_bd_correction_ua(self, bd_corrections_at_qp, material_id):
        '''@brief function used to estimate one-group kappa at quadrature points

        @param bd_corrections_at_qp kappa for all groups at the specific quadrature
        point.
        @param material_id Materail ID for current cell
        @return a one-group kappa
        '''
        return np.dot(bd_corrections_at_qp, self.ksi_ua[material_id])
