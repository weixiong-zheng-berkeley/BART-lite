import numpy as np
from mg_iterations import MG
'''
class used to perform eigenvalue calculations
'''
class Eigen(object):
    def __init__(self, msh_cls):
        # TODO: address constructor errors
        self._tol = 1e-5
        self._k_tol = 1e-5
        # mesh
        self._mesh = msh_cls
        self._n_dof = self._mesh.n_node()
        self._n_grp
        self._mg = MG()

    def do_iterations(self, ho_cls, nda_cls=None):
        '''@brief Function to be called outside for eigenvalue problems

        @param ho_cls The HO class instance
        @param nda_cls NDA instance
        '''
        if not nda_cls:
            # NDA is not used in this calse
            self.eigen_iterations(ho_cls)
        else:
            # TODO: add NDA
            for g in xrange(self.n_grp)

    def eigen_iterations(self, equ_cls):
        '''@brief Function to be called in do_iterations
        '''
        sflxes_prev = [np.ones(self.n_dof) for _ in xrange(self.n_grp)]
        e,ek,keff = 1.0,1.0,1.0
        while e>self._tol and ek>self._k_tol:
            equ_cls.assemble_fixed_linear_forms()
            self.mg.mg_iterations(equ_cls)
            keff_prev,keff = keff,equ_cls.calculate_keff()
            ek,e = abs((keff-keff_prev)/keff),self._tol*0.1
            for g in xrange(self.n_grp):
                e = max(e, equ_cls.calculate_sflx_diff(sflxes_prev,g))
