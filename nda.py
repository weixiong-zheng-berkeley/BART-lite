import numpy as np
from scipy import sparse as sps, sparse.linalg as sla
from itertools import product as pd
import np.linalg.norm as norm
from elem import *

class NDA(object):
    def __init__(self, mat_lib, mat_map, aq_cls):
        # material lib and map
        self._mlib = mat_lib
        self._mmap = mat_map
        self._n_grp = self._mlib.get('n_groups')
        self._n_thr = self._mlib.get('g_thermal')
        # mesh
        self._mesh = msh_cls
        # problem type
        self._is_eigen = True
        self._aq = aq_cls.get_aq_data()
        # number of directions in HO
        self._n_dir = self._aq['n_dir']
        # total number of components: keep consistency with HO
        self._n_tot = self._n_grp
        # all material
        self._sigses = {k:self._mmap.get('sig_s',mat_id=k) for k in self._mmap.ids()}

    def assemble_bilinear_forms(self, correction=False):
        '''@brief A function used to assemble bilinear forms of NDA for current
        iterations

        @param correction A boolean used to determine if correction terms are
        assembled. By default, it's not. In this case, the bilinear form is typical
        diffusion
        '''
        for cell in self._mesh.cells():


    def assemble_fixed_linear_forms(self, sflxes_prev=None):
        '''@brief  function used to assemble linear form for fixed source or fission
        source
        '''

    def assemble_group_linear_forms(self, g):
        '''@brief A function used to assemble linear form for upscattering acceleration
        '''
        # assign fixed source to rhs
        self._sys_rhses[g] = fixed_rhses[g]
        for cell in mesh.cells:
            idx,mid = cell.get('global_idx'),cell.get('mat_id')
            sigs = self._sigses[mid][g,:]
            scat_src = np.zeros(4)
            for gin in filter(lambda x: sigs[x]>1.0e-14,xrange(self._n_grp)):
                local_sflx = self._sflxes[gin][idx]
                scat_src += sigs[gin] * np.dot(mass, local_sflx)
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

    def update_sflxes(self):
        '''@brief A function used to update the scalar fluxes after  upscattering
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
