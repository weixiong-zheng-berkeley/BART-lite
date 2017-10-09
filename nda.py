import numpy as np
from scipy import sparse as sps, sparse.linalg as sla
from itertools import product as pd
import np.linalg.norm as norm
from elem import *

class NDA(object):
    def __init__(self, mat_lib, mat_map, aq_cls):
        # mesh data
        self._mesh = mesh_cls
        self._cell_length = mesh_cls.cell_length()
        # preassembly-interpolation data
        self._elem = Elem(self._cell_length)
        # material lib and map
        self._mlib = mat_lib
        self._mmap = mat_map
        self._n_grp = self._mlib.get('n_groups')
        self._n_thr = self._mlib.get('g_thermal')
        # mesh
        self._mesh = msh_cls
        # problem type
        self._is_eigen = True
        # total number of components: keep consistency with HO
        self._n_tot = self._n_grp
        # all material
        self._dcoefs = self._mlib.get('diff_coef')
        self._sigts = self._mlib.get('sig_t')
        self._sigses = self._mlib.get('sig_s')
        self._fiss_xsecs = self._mlib.get('chi_nu_sig_f')
        # derived material properties
        self._ksi_ua
        self._sigt_ua
        self._diff_coef_ua
        # assistance object
        self._local_dof_pairs = pd(xrange(4),xrange(4))

    def assemble_bilinear_forms(self, ho_cls=None, correction=False):
        '''@brief A function used to assemble bilinear forms of NDA for current
        iterations

        @param correction A boolean used to determine if correction terms are
        assembled. By default, it's not. In this case, the bilinear form is typical
        diffusion
        '''
        # TODO: Boundary is assumed to be reflective so kappa will not be handled
        if correction:
            assert ho_cls is not None, 'ho_cls has to be filled in for correction'
        # basic diffusion Elementary matrices
        diff_mats = {}
        # Elementary correction matrices
        corx,cory,sigt,dcoef = self._elem.corx(),self._elem.cory(),0,0
        for g in xrange(self._n_grp):
            self._sys_mats[g] = sps.lil_matrix((self._mesh.n_node(),self._mesh.n_node()))
            for mid in self._mlib.ids():
                sigt,dcoef = self._sigts[mid][g],self._dcoefs[mid][g]
                diff_mats[(g,mid)] = (dcoef*streaming + sigt*mass)
        # loop over cells for assembly
        for cell in self._mesh.cells():
            # get global dof index and mat id
            idx,mid = cell.id(),cell.global_idx()
            # corrections for all groups in current cell
            corr_vecs = {}
            for g in xrange(self._n_grp):
                cor_mat = np.zeros((4,4))
                corr_at_qp = []
                if correction:
                    # calculate NDA correction in HO class
                    corr_vecs[g]=ho_cls.calculate_nda_cell_correction(
                    g=g, mat_id=mid, idx=idx)
                    for i in xrange(len(corr_vecs[g])):
                        # x-component
                        cor_mat += corr_vecs[g][i][0]*corx[i]
                        # y-component
                        cor_mat += corr_vecs[g][i][1]*cory[i]
                # assemble global system
                for ci,cj in self._local_dof_pairs:
                    self._sys_mats[g][idx[ci],idx[cj]]+=\
                    diff_mats[(g,mid)][ci,cj]+cor_mat[ci,cj]

            if self._is_ua:

        # Transform system matrices to CSC format
        for g in xrange(self._n_grp):
            self._sys_mats[g] = sps.csc_matrix(sys_mats[g])
        if self._is_ua:
            self._sys_mats['ua'] = sps.csc_matrix(sys_mats['ua'])

    def assemble_fixed_linear_forms(self, sflxes_prev=None):
        '''@brief  function used to assemble linear form for fixed source or fission
        source
        '''
        # scale the fission xsec by keff
        scaled_fiss_xsec = {k:v/self._keff for k,v in self._fiss_xsecs}
        for g in xrange(self._n_grp):
            fixed_rhses[g] = np.array(self._mesh.n_node())
            for cell in self._mesh.cells():
                idx,mid,local_fixed = cell.global_idx(),cell.id(),np.zeros(4)
                for gin in filter(lambda: scaled_fiss_xsec[mid][g,x]>1.0e-14,xrange(self._n_grp)):
                    sflx_vtx = self._sflxes[g][idx]
                    local_fixed += scaled_fiss_xsec[mid][g,gin]*np.dot(mass,sflx_vtx)
                self._fixed_rhses[g][idx] += local_fixed

    def assemble_group_linear_forms(self, g):
        '''@brief A function used to assemble linear form for upscattering acceleration
        '''
        # assign fixed source to rhs
        self._sys_rhses[g] = fixed_rhses[g]
        for cell in mesh.cells:
            idx,mid = cell.global_idx(),cell.id()
            sigs = self._sigses[mid][g,:]
            scat_src = np.zeros(4)
            for gin in filter(lambda x: sigs[x]>1.0e-14,xrange(self._n_grp)):
                sflx_vtx = self._sflxes[gin][idx]
                scat_src += sigs[gin] * np.dot(mass, sflx_vtx)
            self._sys_rhses[g][idx] += scat_src

    def assemble_ua_linear_form(self):
        '''@brief A function used to assemble linear form for upscattering acceleration
        '''
        #TODO: fill in this function and put the rhs in _sys_rhses['ua']

    def solve_in_group(self,g):
        assert 0<=g<self._n_grp, 'Group index out of range'
        if g not in self._lu:
            # factorize it if not yet
            self._lu[g] = sla.splu(self._sys_mats[g])
        # direct solve
        self._sflxes[g] = self._lu[g].solve(self._sys_rhses[g])

    #NOTE: this function has to be removed if abstract class is implemented
    def calculate_keff_err(self):
        assert not self._is_eigen, 'only be called in eigenvalue problems'
        # update the previous fission source and previous keff
        self._fiss_src_prev,self._keff_prev = self._fiss_src,self._keff
        # calculate the new fission source
        self._calculate_fiss_src()
        # calculate the new keff
        self._keff = self._keff_prev * self._fiss_src / self._fiss_src_prev
        return abs(self._keff-self._keff_prev)/abs(self._keff)

    #NOTE: this function has to be removed if abstract class is implemented
    def calculate_sflx_diff(self, sflxes_old, g):
        '''@brief function used to generate ho scalar flux for Group g using
        angular fluxes

        @param sflx_old Scalar flux from previous generation
        @param g The group index
        @return double The relative difference between new and old scalar flux
        '''
        # return the l1 norm relative difference
        return norm((self._sflxes[g]-sflxes_old[g]),1) / norm(self._sflxes[g],1)

    #NOTE: this function has to be removed if abstract class is implemented
    def update_sflxes(self, sflxes_old, g):
        '''@brief A function used to update scalar flux for group g

        @param sflxes_old A dictionary
        @param g Group index
        '''
        sflxes_old[g] = self._sflxes[g]

    def update_ua(self):
        '''@brief A function used to update the scalar fluxes after upscattering
        acceleration
        '''
        for g in xrange(self.g_thr,self._n_grp):
            self._sflxes[g] += self._sflxes['ua']

    def clear_factorization(self):
        '''@brief A function used to clear all the factorizations after NDA is dictionaries
        for current iteration

        Every outer iterations, NDA equations are modified so previous factorizations
        are no longer suitable and must be cleared and redone.
        '''
        self._lu.clear()

    def solve_ua(self):
        if 'ua' not in self._lu:
            # factorize it if not yet
            self._lu['ua'] = sla.splu(self._sys_mats['ua'])
        # direct solve
        self._sflxes['ua'] = self._lu['ua'].solve(self._sys_rhses['ua'])
