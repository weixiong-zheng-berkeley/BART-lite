import numpy as np
from scipy import sparse, sparse.linalg as sla
from itertools import product as pd
from elem import *

class SAAF(object):
    def __init__(self, mat_cls, aq_cls):
        #TODO: address constructor
        # material data
        self.mat_cls = mat_cls
        self.n_grp = mat_cls.get_n_groups()
        # problem type: is problem eigenvalue problem
        self.is_eigen = mat_cls.get_eigen_bool()
        # aq data in forms of dictionary
        self.aq = aq_cls.get_aq_data()
        self.n_dir = self.aq['n_dir']
        # total number of components in HO
        self.n_tot = self.n_grp * self.n_dir
        # get a component indexing mapping
        self.comp = dict()
        # component to group map
        self.comp_grp = dict()
        # component to direction map
        self.comp_dir = dict()
        self._generate_component_map()
        # matrices and vectors
        self.n_dof#
        self.sys_mats = [0 for _ in xrange(self.n_tot)]
        # be very careful about the following
        self.sys_rhses = [np.array(np.ones(self.n_dof)) for _ in xrange(self.n_tot)]
        self.fixed_rhses = [np.array(np.ones(self.n_dof)) for _ in xrange(self.n_tot)]
        self.aflxes = [np.array(np.ones(self.n_dof)) for _ in xrange(self.n_tot)]
        # scalar flux for current calculation
        self.sflxes = [np.array(np.ones(self.n_dof)) for _ in xrange(self.n_grp)]
        # linear solver objects
        self.lu = [0] * self.n_tot
        # source iteration tol
        self._tol = 1.0e-7

    def _generate_component_map(self):
        '''@brief Internal function used to generate mappings between component,
        group and directions
        '''
        ct = 0
        for g in xrange(self.n_grp):
            for d in xrange(self.n_dir):
                self.comp[(g,d)] = ct
                self.comp_grp[ct] = g
                self.comp_dir[ct] = d
                ct += 1

    def assemble_fixed_linear_forms(self, sflxes_prev=None):
        '''@brief a function used to assemble fixed source or fission source on the
        rhs for all components

        Generate numpy arrays and put them in self.
        '''
        # TODO

    def assemble_bilinear_forms(self):
        '''@brief Function used to assemble bilinear forms

        Must be called only once.
        '''
        mass = np.array([[]])
        for i in xrange(self.n_tot):
            # get the group and direction indices
            g,d = self.comp_grp[i],self.comp_dir[i]
            omega = self.aq[d]['omega']
            strm_pre = np.zeros((4,4))
            mass_pre = np.zeros((4,4))
            oxox,oxoy,oyox,oyoy = self.aq[d]['dir_prods']
            local_matrices = []
            for m in xrange(self.n_mat):
                local_sigt = self.mat_cls.get('sig_t',m)[g]#TODO: check with Josh
                local_inv_sigt = self.mat_cls.get('inv_sig_t',m)[g]
                local_matrix = oxox*dxdx + oxoy*(dxdy + dydx) + oyoy*dydy
        # TODO: fill in bilinear forms assembly code

    def assemble_group_linear_forms(self,g):
        '''@brief Function used to assemble linear forms for Group g
        '''
        assert 0<=g<=self.n_grp, 'Group index out of range'
        # TODO: fill in linear form assembly code

    def solve_in_group(self, g):
        '''@brief Called to solve direction by direction inside Group g

        @param g Group index
        '''
        assert 0<=g<self.n_grp, 'Group index out of range'
        # Source iteration
        e = 1.0
        while e>self._tol:
            # assemble group rhses
            self.assemble_group_linear_forms(g)
            sflx_old = self.sflxes[g]
            self.sflxes *= 0
            for d in xrange(self.n_dir):
                # if not factorized, factorize the the HO matrices
                if self.lu[comp[(g,d)]]==0:
                    self.lu[comp[(g,d)]] = sla.splu(self.sys_mats[comp[(g,d)]])
                # solve direction d
                self.aflxes[comp[(g,d)]] = self.lu[comp[(g,d)]].solve(self.sys_rhses[comp[(g,d)]])
                self.sflxes[g] += self.aq[d]['wt'] * self.aflxes[comp[(g,d)]]
            # calculate difference for SI convergence
            e = np.linalg.norm(sflx_old - self.sflxes[g],1) / \
            np.linalg.norm (self.sflxes[g],1)

    def calculate_keff(self):
        assert not self.is_eigen, 'only be called in eigenvalue problems'
        # TODO: fill in this function

    def calculate_sflx_diff(self, sflxes_old, g):
        '''@brief function used to generate ho scalar flux for Group g using
        angular fluxes

        @param sflx_old Scalar flux from previous generation
        @param g The group index
        @return double The relative difference between new and old scalar flux
        '''
        # retrieve old values for scalar flux
        sflxes_old[g] = self.sflxes[g]
        # generate new scalar flux
        self.sflxes[g] = self.aq[0]['wt'] * self.aflxes[comp[(g,0)]]
        for d in xrange(1, self.n_dir):
            self.sflxes[g] += self.aq[d]['wt'] * self.aflxes[comp[(g,d)]]
        # return the l1 norm relative difference
        return np.linalg.norm((self.sflxes[g]-sflxes_old[g]),1) / \
        np.linalg.norm(self.sflxes[g],1)

    def get_sflxes(self, sflxes, g):
        '''@brief Function called outside to retrieve the scalar flux value for Group g
        '''
        sflxes[g] = self.sflxes[g]

    def get_aflxes(self):
        '''@brief A function used to retrieve angular fluxes for NDA use

        Implementation to be determined
        '''
        # TODO: fill this in.
