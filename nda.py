import numpy as np
from scipy import sparse as sps, sparse.linalg as sla
from itertools import product as pd
import np.linalg.norm as norm
from elem import *

class NDA(object):
    def __init__(self, mat_lib, msh_cls, is_ua=False):
        # equation name
        self._name = 'nda'
        # mesh data
        self._mesh = mesh_cls
        self._cell_length = mesh_cls.cell_length()
        # preassembly-interpolation data
        self._elem = Elem(self._cell_length)
        # material lib and map
        self._mlib = mat_lib
        self._n_grp = self._mlib.get('n_groups')
        self._g_thr = self._mlib.get('g_thermal')
        # mesh
        self._mesh = msh_cls
        # problem type
        # TODO: modify problem type read-in
        self._is_eigen = True
        self._is_ua = is_ua
        # total number of components: keep consistency with HO
        self._n_tot = self._n_grp
        # all material
        self._dcoefs = self._mlib.get('diff_coef')
        self._sigts = self._mlib.get('sig_t')
        self._sigses = self._mlib.get('sig_s')
        self._sigrs = self._mlib.get('sig_r')
        self._fiss_xsecs = self._mlib.get('chi_nu_sig_f')
        # derived material properties
        self._ksi_ua = self._mlib.get('ksi_ua')
        self._sigrs_ua = self._mlib.get('sig_r_ua')
        self._dcoefs_ua = self._mlib.get('diff_coef_ua')
        # assistance object
        self._local_dof_pairs = pd(xrange(4),xrange(4))

    def name(self):
        return self._name

    def assemble_bilinear_forms(self, ho_cls=None, correction=False):
        '''@brief A function used to assemble bilinear forms of NDA for current
        iterations

        @param correction A boolean used to determine if correction terms are
        assembled. By default, it's not. In this case, the bilinear form is typical
        diffusion
        '''
        # TODO: Boundary is assumed to be reflective so kappa will not be handled
        streaming,mass = self._elem.streaming(),self._elem.mass()
        if correction:
            assert ho_cls is not None, 'ho_cls has to be filled in for NDA correction'
        # basic diffusion Elementary matrices
        diff_mats = {}
        # Elementary correction matrices
        corx,cory,sigt,dcoef = self._elem.corx(),self._elem.cory(),0,0
        for g in xrange(self._n_grp):
            self._sys_mats[g] = sps.lil_matrix((self._mesh.n_node(),self._mesh.n_node()))
            for mid in self._mlib.ids():
                sigt,sigr,dcoef = self._sigts[mid][g],self._sigrs[mid][g],self._dcoefs[mid][g]
                diff_mats[(g,mid)] = (dcoef*streaming + sigr*mass)
        # preassembled matrices for upscattering acceleration
        if self._is_ua:
            self._sys_mats['ua'] = sps.lil_matrix((self._mesh.n_node(),self._mesh.n_node()))
            for mid in self._mlib.ids():
                dcoef_ua,sigr_ua = self._dcoefs_ua[mid],self._sigrs_ua[mid]
                # basic elementary diffusion matrices for upscattering acceleration
                diff_mats[('ua',mid)] = (dcoef_ua*streaming + sigr_ua*mass)
        # loop over cells for assembly
        for cell in self._mesh.cells():
            # get global dof index and mat id
            idx,mid = cell.id(),cell.global_idx()
            # corrections for all groups in current cell and ua
            corr_vecs,corx_ua,cory_ua = {},{},{}
            for g in xrange(self._n_grp):
                # if correction is asked
                cor_mat = np.zeros((4,4))
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

            # if we do upscattering acceleration
            if self._is_ua:
                # assemble global system of ua matrix without correction
                for ci,cj in self._local_dof_pairs:
                    self._sys_mats['ua'][idx[ci],idx[cj]]+=\
                    diff_mats[('ua',mid)][ci,cj]
                # correction matrix for upscattering acceleration
                cor_mat_ua = np.zeros((4,4))
                if correction:
                    for i in xrange(len(corr_vecs[0])):
                        corx_ua_qp,cory_ua_qp = 0,0
                        for g in xrange(self._g_thr,self._n_grp):
                            ksi = self._ksi_ua[g-self._g_thr]
                            corx_qp,cory_qp = corr_vecs[g][i][0],corr_vecs[g][i][1]
                            corx_ua_qp += ksi*corx_qp
                            cory_ua_qp += ksi*cory_qp
                        corx_ua[i],cory_ua[i] = corx_ua_qp,cory_ua_qp
                        cor_mat_ua+=(corx_ua[i]*corx[i]+cory_ua[i]*cory[i])
                for ci,cj in self._local_dof_pairs:
                    self._sys_mats['ua'][idx[ci],idx[cj]]+=\
                    diff_mats[('ua',mid)][ci,cj]+cor_mat_ua[ci,cj]

        # Transform system matrices to CSC format
        for g in xrange(self._n_grp):
            self._sys_mats[g] = sps.csc_matrix(self._sys_mats[g])
        if self._is_ua:
            self._sys_mats['ua'] = sps.csc_matrix(self._sys_mats['ua'])

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
        self._sys_rhses[g],mass = fixed_rhses[g],self._elem.mass()
        for cell in mesh.cells:
            idx,mid = cell.global_idx(),cell.id()
            sigs = self._sigses[mid][g,:]
            scat_src = np.zeros(4)
            for gin in filter(lambda x: sigs[x]>1.0e-14,xrange(self._n_grp)):
                sflx_vtx = self._sflxes[gin][idx]
                scat_src += sigs[gin] * np.dot(mass, sflx_vtx)
            self._sys_rhses[g][idx] += scat_src

    def assemble_ua_linear_form(self, sflxes_old):
        '''@brief A function used to assemble linear form for upscattering acceleration
        '''
        assert len(sflxes_old)==self._n_grp, \
        'old scalar fluxes should have the same number of groups as current scalar fluxes'
        mass = self._elem.mass()
        self._sys_rhses['ua'] *= 0.0
        for cell in self._mesh.cells():
            idx,mid,scat_src_ua = cell.global_idx(),cell.id(),np.zeros(4)
            for g in xrange(self._g_thr,self._n_grp-1):
                for gin in xrange(g+1,self._n_grp):
                    sigs = self._sigses[mid][g,gin]
                    if sigs>1.0e-14:
                        dsflx_vtx = self._sflxes[gin][idx]-sflxes_old[g][idx]
                        scat_src_ua += sigs*np.dot(mass,dsflx_vtx)
            self._sys_rhses['ua'][idx] += scat_src_ua

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
