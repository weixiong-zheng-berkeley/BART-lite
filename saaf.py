import numpy as np
from scipy import sparse as sps, sparse.linalg as sla
from itertools import product as pd
import np.linalg.norm as norm
from elem import Elem

class SAAF(object):
    def __init__(self, mat_cls, aq_cls, msh_cls):
        # mesh data
        self._mesh = mesh_cls
        self._cell_length = mesh_cls.cell_length()
        # preassembly-interpolation data
        self._elem = Elem(self._cell_length)
        # material data
        self._mlib = mat_lib
        self._n_grp = mat_cls.get_n_groups()
        self._sigts = mat_cls.get('sig_t')
        self._isigts = mat_cls.get('inv_sig_t')
        self._fiss_xsecs = mat_cls.get_per_str('chi_nu_sig_f')
        self._sigses = mat_cls.get_per_str('sig_s')
        self._mids = mat_cls.ids()
        # problem type: is problem eigenvalue problem
        # aq data in forms of dictionary
        self._aq = aq_cls.get_aq_data()
        self._n_dir = self._aq['n_dir']
        # total number of components in HO
        self._n_tot = self._n_grp * self._n_dir
        # get a component indexing mapping
        self._comp = dict()
        # component to group map
        self._comp_grp = dict()
        # component to direction map
        self._comp_dir = dict()
        self._generate_component_map()
        # local vectors
        self._rhs_mats = dict()
        self._preassembly_rhs()
        # related to global matrices and vectors
        self._n_dof = msh_cls.n_node()
        self._sys_mats = {}
        # be very careful about the following
        self._sys_rhses = {k:np.ones(self._n_dof) for k in xrange(self._n_tot)}
        self._fixed_rhses = {k:np.zeros(self._n_dof) for k in xrange(self._n_tot)}
        self._aflxes = {k:np.ones(self._n_dof) for k in xrange(self._n_tot)}
        # scalar flux for current calculation
        self._sflxes = {k:np.ones(self._n_dof) for k in xrange(self._n_grp)}
        # linear solver objects
        self._lu = {}
        # source iteration tol
        self._tol = 1.0e-7
        # fission source
        self._fiss_src = self._calculate_fiss_src(self._sflxes)
        self._fiss_src_prev

    def _generate_component_map(self):
        '''@brief Internal function used to generate mappings between component,
        group and directions
        '''
        ct = 0
        for g in xrange(self._n_grp):
            for d in xrange(self._n_dir):
                self._comp[(g,d)] = ct
                self._comp_grp[ct] = g
                self._comp_dir[ct] = d
                ct += 1

    def _preassembly_rhs(self):
        for mid in self._ids:
            sigts,isigts = self._get('sig_t')[mid],self._get('inv_sig_t')[mid]
            for g in xrange(self._n_grp):
                for d in xrange(self._n_dir):
                    ox,oy = self._aq['omega'][d]
                    # streaming part of rhs
                    rhs_mat = (ox*dxvu+oy*dyvu) * isigts[g]
                    # mass part of rhs
                    rhs_mat += mass
                    self._rhs_mats[mid] = {(g,d):rhs_mat}

    def assemble_bilinear_forms(self):
        '''@brief Function used to assemble bilinear forms

        Must be called only once.
        '''
        # retrieve all the material properties
        for i in xrange(self._n_tot):
            # get the group and direction indices
            g,d = self._comp_grp[i],self._comp_dir[i]
            # get omega_i * omega_j combinations
            oxox,oxoy,oyoy = self._aq['dir_prods'][d].values()
            # dict containing lhs local matrices for all materials for component i
            lhs_mats = dict()
            for mid in self._ids:
                sigt,isigt = self._sigts[mid][g],self._isigts[mid][g]
                # streaming lhs
                matx = isigt * (oxox*self._elem.dxdx() +
                                oxoy*(self._elem.dxdy() + self._elem.dydx()) +
                                oyoy*self._elem.dydy())
                # collision matrix
                matx += sigt * self._elem.mass()
                lhs_mats[mid] = matx
            # loop over cells for assembly
            # sys_mat: temp variable for system matrix for one component
            sys_mat = sps.lil_matrix((self._mesh.n_node(), self._mesh.n_node()))
            for cell in self._mesh.cells():
                idx,mid = cell.global_idx(),cell.id()
                for ci in xrange(4):
                    for cj in xrange(4):
                        sys_mat[idx[ci],idx[cj]] += lhs_mats[mid][ci][cj]
            # transform lil_matrix to csc_matrix for efficient computation
            self._sys_mats[i] = sps.csc_matrix(sys_mat)

    def assemble_fixed_linear_forms(self, sflxes_prev=None, keff=None):
        '''@brief a function used to assemble fixed source or fission source on the
        rhs for all components

        Generate numpy arrays and put them in self._
        '''
        # TODO: add fixed source part. Current program is for eigenvalue Only
        assert sflxes_prev is not None and keff is not None, 'Only eigenvalue is implemented'
        # get properties per str scaled by keff
        for cp in xrange(self._n_tot):
            # re-init fixed rhs. This must be done at the beginning of calling this function
            self._fixed_rhses[cp] = np.array(self._n_dof)
            # get group and direction indices
            g,d = self._comp_grp[cp],self._comp_dir[cp]
            for cell in self._mesh.cells():
                idx,mid = cell.global_idx(),cell.id()
                local_fiss,fiss_xsec = np.zeros(4),self._fiss_xsecs[mid][g]
                # get fission source contribution from ingroups
                for gin in filter(lambda j: fiss_xsec[j]>1.0e-14, xrange(self._n_grp)):
                    local_sflx = sflxes_prev[gin][idx]
                    local_fiss += fiss_xsec[gin]*np.dot(self._rhs_mats[mid][(g,d)],local_sflx)
                self._fixed_rhses[cp][idx] += local_fiss

    def assemble_group_linear_forms(self,g):
        '''@brief Function used to assemble linear forms for Group g
        '''
        assert 0<=g<self._n_grp, 'Group index out of range'
        for d in xrange(self._n_dir):
            cp = self._comp[(g,d)]
            # get fixed/fission source
            sys_rhses[cp] = fixed_rhses[cp]
            for cell in self._mesh.cells():
                # get global dof indices and material ids
                idx,mid = cell.global_idx(),cell.id()
                # get scattering matrix for current cell
                sigs = self._sigses[mid]
                # calculate local scattering source
                scat_src = np.zeros(4)
                for gin in filter(lambda x: sigs[g][x]>1.0e-14, xrange(self._n_grp)):
                    local_sflx = self._sflxes[g][idx]
                    scat_src += sigs[g,gin]*np.dot(self._rhs_mats[mid][(g,d)], local_sflx)
                sys_rhses[cp][idx] += scat_src

    def solve_in_group(self, sflxes_old, g):
        '''@brief Called to solve direction by direction inside Group g

        @param g Group index
        '''
        assert 0<=g<self._n_grp, 'Group index out of range'
        # update old flux
        sflxes_old[g] = self._sflxes[g]
        # Source iteration
        e = 1.0
        while e>self._tol:
            # assemble group rhses
            self._assemble_group_linear_forms(g)
            sflx_old = self._sflxes[g]
            self._sflxes[g] *= 0
            for d in xrange(self._n_dir):
                # if not factorized, factorize the the HO matrices
                if self._lu[comp[(g,d)]]==0:
                    self._lu[comp[(g,d)]] = sla.splu(self._sys_mats[comp[(g,d)]])
                # solve direction d
                self._aflxes[comp[(g,d)]] = self._lu[comp[(g,d)]].solve(self._sys_rhses[comp[(g,d)]])
                self._sflxes[g] += self._aq['wt'][d] * self._aflxes[comp[(g,d)]]
            # calculate difference for SI convergence
            e = norm(sflx_old - self._sflxes[g],1) / norm (self._sflxes[g],1)

    def calculate_keff(self):
        assert not self._is_eigen, 'only be called in eigenvalue problems'
        self._calculate_fiss_src()

    def _calculate_fiss_src(self):
        self._fiss_src = 0
        # loop over cells

    def calculate_sflx_diff(self, sflxes_old, g):
        '''@brief function used to generate ho scalar flux for Group g using
        angular fluxes

        @param sflx_old Scalar flux from previous generation
        @param g The group index
        @return double The relative difference between new and old scalar flux
        '''
        # return the l1 norm relative difference
        return norm((self._sflxes[g]-sflxes_old[g]),1) / norm(self._sflxes[g],1)

    def get_sflxes(self, sflxes, g):
        '''@brief Function called outside to retrieve the scalar flux value for Group g

        @param sflxes A dictionary to be modified to contain all scalar fluxes
        '''
        sflxes[g] = self._sflxes[g]

    def get_keff(self):
        '''@brief A function used to retrieve keff

        @return keff calculated in SAAF class
        '''
        return self._keff

    def get_aflxes_at_qp(self, cell, aflxes_g, g):
        '''@brief A function used to retrieve angular fluxes for NDA use

        @param cell Cell index
        @param aflxes_g A dictionary to be modified to fill in angular fluxes at quadrature points
        @param g Target group
        '''
        for d in xrange(self._n_dir):
            cp,idx = self._comp[(g,d)],cell.global_idx()
            sol_at_vertices = self._aflxes[cp][idx]
            aflxes_g[d] = get_sol_at_qps(sol_at_vertices)

    def get_grad_aflxes_at_qp(self, cell, grad_aflxes_g, g):
        '''@brief A function used to retrieve aflxes gradients for NDA use

        @param cell Cell index
        @param grad_aflxes_g A dictionary to be modified to fill in gradients of
        angular fluxes in terms of tuples at quadrature points
        @param g Group number
        '''
        for d in xrange(self._n_dir):
            cp,idx = self._comp[(g,d)],cell.global_idx()
            sol_at_vertices = self._aflxes[cp][idx]
            grad_aflxes_g[d] = self._elem.get_grad_at_qps(sol_at_vertices)
