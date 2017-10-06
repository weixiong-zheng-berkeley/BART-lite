import numpy as np
from scipy import sparse as sps, sparse.linalg as sla
from itertools import product as pd
import np.linalg.norm as norm
from elem import *

class SAAF(object):
    def __init__(self, mat_lib, mat_map, aq_cls):
        #TODO: address constructor
        # material data
        self.mlib = mat_lib
        self.mmap = mat_map
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
        # local vectors
        self.rhs_mats = dict()
        self._preassembly_rhs()
        # related to global matrices and vectors
        self.n_dof = (1+n_x) * (1+n_x)
        self.sys_mats = {k:v for k,v in zip(xrange(self.n_tot),[0]*self.n_tot)}
        # be very careful about the following
        self.sys_rhses = {k:np.ones(self.n_dof) for k in xrange(self.n_tot)}
        self.fixed_rhses = {k:np.ones(self.n_dof) for k in xrange(self.n_tot)}
        self.aflxes = {k:np.ones(self.n_dof) for k in xrange(self.n_tot)}
        # scalar flux for current calculation
        self.sflxes = {k:np.ones(self.n_dof) for k in xrange(self.n_grp)}
        # linear solver objects
        self.lu = [0] * self.n_tot
        # source iteration tol
        self._tol = 1.0e-7
        # cell mapping about global index, id and reflective sides
        self.cell_info = dict()
        self._init_cell_info ()
        # fission source
        self.fiss_src = self._calculate_fiss_src(self.sflxes)
        self.fiss_src_prev

    #TODO: move the following function to cell related place
    def _init_cell_info (self):
        # TODO: add N_CELL into build_cells
        for c in xrange(N_CELL):
            i,j = c//n_x,c%n_x
            global_idx = np.array([i*n_x+j,i*n_x+j+1,(i+1)*n_x+j,(i+1)*n_x+j+1])
            #TODO: check with Josh for explanation about the cell numbering
            ctr = ((j+0.5)*CELL_LENGTH,(i+0.5)*CELL_LENGTH)
            # list for all reflective boundary sides
            sides =
            #self.cell_info[c] = {'global_idx':global_idx,'id':self.mmap.get('id',c)}
            self.cell_info[c] = {'global_idx':global_idx,'id':self.mmap.get('id',ctr)}

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

    def _preassembly_rhs(self):
        for mid in self.mlib.ids():
            sigts,isigts = self.get('sig_t')[mid],self.get('inv_sig_t')[mid]
            for g in xrange(self.n_grp):
                for d in xrange(self.n_dir):
                    ox,oy = self.aq[d]['dir']
                    # streaming part of rhs
                    rhs_mat = (ox*dxvu+oy*dyvu) * isigts[g]
                    # mass part of rhs
                    rhs_mat += mass
                    self.rhs_mats[mid] = {(g,d):rhs_mat}

    def assemble_bilinear_forms(self):
        '''@brief Function used to assemble bilinear forms

        Must be called only once.
        '''
        # retrieve all the material properties
        for i in xrange(self.n_tot):
            # get the group and direction indices
            g,d = self.comp_grp[i],self.comp_dir[i]
            # get omega_i * omega_j combinations
            oxox,oxoy,oyox,oyoy = self.aq[d]['dir_prods']
            # dict containing lhs local matrices for all materials for component i
            lhs_mats = dict()
            sigts,isigts = self.mmap.get('sig_t'),self.mmap.get('inv_sig_t')
            for mid in self.mlib.ids():
                sigt,isigt = sigts[mid][g],isigts[mid][g]
                # streaming lhs
                matx = isigt * (oxox*dxdx + oxoy*(dxdy + dydx) + oyoy*dydy)
                # collision matrix
                matx += sigt * mass
                lhs_mats[mid] = matx
            # loop over cells for assembly
            # sys_mat: temp variable for system matrix for one component
            sys_mat = sps.lil_matrix((nodes, nodes))
            for cell in mesh.cells:
                idx,mid = cell.global_idx(),cell.id()#TODO: modify after Josh's work
                for ci in xrange(4):
                    for cj in xrange(4):
                        sys_mat[idx[ci]][idx[cj]] += lhs_mats[mid][ci][cj]
            # transform lil_matrix to csc_matrix
            sys_mats[i] = sps.csc_matrix(sys_mat)

    def assemble_fixed_linear_forms(self, sflxes_prev=None, keff=None):
        '''@brief a function used to assemble fixed source or fission source on the
        rhs for all components

        Generate numpy arrays and put them in self.
        '''
        # TODO: add fixed source part. Current program is for eigenvalue Only
        assert sflxes_prev is not None and keff is not None, 'Only eigenvalue is implemented'
        # get properties per str scaled by keff
        xsec_f = {k:(self.mlib.get_per_str('chi_nu_sig_f',mat_id=k))/keff for k in self.mlib.ids()}
        for cp in xrange(self.n_tot):
            # get group and direction indices
            g,d = self.comp_grp[cp],self.comp_dir[cp]
            for c in xrange(N_CELL):
                idx,mid = self.cell_info[c]['global_idx'],self.cell_info[c]['id']
                local_flx = sflxes_prev[idx]
                local_fiss = np.zeros(4)
                # get fission source contribution from ingroups
                for gin in filter(lambda j: xsec_f[mid][g][j]>1.0e-16, xrange(self.n_grp)):
                    local_fiss += np.dot(self.rhs_mats[mid][(g,d)]*xsec_f[mid][g][gin],local_flx)
                # Note: the following set-zero is important
                fixed_rhses[cp][idx] = 0.0
                for ci in xrange(4):
                    fixed_rhses[cp][idx[ci]] += local_fiss[ci]

    def assemble_group_linear_forms(self,g):
        '''@brief Function used to assemble linear forms for Group g
        '''
        assert 0<=g<self.n_grp, 'Group index out of range'
        sigses = {k:self.mlib.get_per_str('sig_s',mat_id=k) for k in self.mlib.ids()}
        for d in xrange(self.n_dir):
            cp = self.comp[(g,d)]
            # get fixed/fission source
            sys_rhses[cp] = fixed_rhses[cp]
            for c in xrange(N_CELL):
                # get global dof indices and material ids
                idx,mid = self.cell_info[c]['global_idx'],self.cell_info[c]['id']
                # get scattering matrix for current cell
                sigs = sigses[mid]
                # calculate local scattering source
                scat_src = np.zeros(4)
                for gin in filter(lambda x: sigs[g][x]>1.0e-14, xrange(self.n_grp)):
                    scat_src += self.rhs_mats[mid][(g,d)] * self.sflxes[g]
                for ci in xrange(4):
                    sys_rhses[cp][idx] += scat_src[ci]

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
            e = norm(sflx_old - self.sflxes[g],1) / norm (self.sflxes[g],1)

    def calculate_keff(self):
        assert not self.is_eigen, 'only be called in eigenvalue problems'
        self._calculate_fiss_src()

    def _calculate_fiss_src(self):


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
        self.sflxes[g] = self.aq['wt'][0] * self.aflxes[comp[(g,0)]]
        for d in xrange(1, self.n_dir):
            self.sflxes[g] += self.aq[d]['wt'] * self.aflxes[comp[(g,d)]]
        # return the l1 norm relative difference
        return norm((self.sflxes[g]-sflxes_old[g]),1) / norm(self.sflxes[g],1)

    def get_sflxes(self, sflxes, g):
        '''@brief Function called outside to retrieve the scalar flux value for Group g

        @param sflxes A dictionary to be modified to contain all scalar fluxes
        '''
        sflxes[g] = self.sflxes[g]

    def get_keff(self):
        '''@brief A function used to retrieve keff

        @return keff calculated in SAAF class
        '''
        return self.keff

    def get_aflxes_at_qp(self, cell, aflxes_g, g):
        '''@brief A function used to retrieve angular fluxes for NDA use

        @param cell Cell index
        @param aflxes_g A dictionary to be modified to fill in angular fluxes at quadrature points
        @param g Target group
        '''
        for d in xrange(self.n_dir):
            cp,idx = self.comp[(g,d)],self.cell_info[cell]['global_idx']
            sol_at_vertices = self.aflxes[cp][idx]
            aflxes_g[d] = get_sol_at_qps(sol_at_vertices)

    def get_grad_aflxes_at_qp(self, cell, grad_aflxes_g, g):
        '''@brief A function used to retrieve aflxes gradients for NDA use

        @param cell Cell index
        @param grad_aflxes_g A dictionary to be modified to fill in gradients of
        angular fluxes in terms of tuples at quadrature points
        @param g Group number
        '''
        for d in xrange(self.n_dir):
            cp,idx = self.comp[(g,d)],self.cell_info[cell]['global_idx']
            sol_at_vertices = self.aflxes[cp][idx]
            grad_aflxes_g[d] = get_grad_at_qps(sol_at_vertices)
