from nose.tools import *
from mesh import Cell
import numpy as np

class TestCells:
    # Tests to verify the cell class is working properly

    @classmethod
    def setup_class(cls):
        cls.index = (1,1)
        cls.mesh_params = {'n_cell': 16,
                           'x_cell': 4,
                           'y_cell': 4,
                           'cell_length': 2.5}
        cls.cell = Cell(cls.index, cls.mesh_params)

    def test_area(self):
        """ Cell should be the correct area """
        eq_(self.cell.area(), 2.5**2, "area value")
        eq_(self.cell.length(), 2.5, "length value")

    def test_global_index(self):
        """ Given an index, constructor should generate the correct
        global indexes """
        global_idx = [6, 7, 11, 12]
        ok_(type(self.cell.global_idx()) == list, "returns list")
        eq_(self.cell.global_idx(), global_idx, "correct value")
        
    @raises(AssertionError)
    def test_init_index(self):
        bad_cell = Cell(5, self.mesh_params)

    @raises(AssertionError)
    def test_length_neg(self):
        bad_params = self.mesh_params
        bad_params['cell_length'] = -1
        bad_cell = Cell(self.index, bad_params)

    @raises(TypeError)
    def test_length_wrong_type(self):
        bad_params = self.mesh_params
        bad_params['cell_length'] = 'c'
        bad_cell = Cell(self.index, bad_params)
        
    @raises(KeyError)
    def test_length_missing(self):
        bad_params = self.mesh_params
        del bad_params['cell_length']
        bad_cell = Cell(self.index, bad_params)
