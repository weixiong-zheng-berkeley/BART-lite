from nose.tools import *
from mesh import Cell
import numpy as np

class TestCells:
    # Tests to verify the cell class is working properly

    @classmethod
    def setup_class(cls):
        cls.index = (1,1)
        cls.mesh_params = {'x_cell': 4,
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

    def test_not_boundary(self):
        ok_(not self.cell.bounds(), "non boundary cell")
        
    def test_is_boundary(self):
        for index in [(4,1), (2, 4), (0, 2), (3,0)]:
            b_cell = Cell(index, self.mesh_params)
            ok_(b_cell.bounds(), "boundary cell")

    def test_boundaries(self):
        for param in [('x_max', (4,1)),
                      ('y_max', (2,4)),
                      ('x_min', (0,2)),
                      ('y_min', (3,0))]:
            b_cell = Cell(param[1], self.mesh_params)
            ok_(param[0] in b_cell.bounds())

    def test_multiple_bounds(self):
        index = (4,4)
        b_cell = Cell(index, self.mesh_params)
        ok_('y_max' in b_cell.bounds() and 'x_max' in b_cell.bounds())

    def test_set_bounds(self):
        index = (4,4)
        b_cell = Cell(index, self.mesh_params)
        b_cell.bounds('x_max', 'refl')
        eq_(b_cell.bounds('x_max'), 'refl')
        
    @raises(AssertionError)
    def test_init_index(self):
        bad_cell = Cell(5, self.mesh_params)

    @raises(AssertionError)
    def test_length_neg(self):
        bad_params = {key: value for key, value in
                      self.mesh_params.items()}
        bad_params['cell_length'] = -1
        bad_cell = Cell(self.index, bad_params)

    @raises(TypeError)
    def test_length_wrong_type(self):
        bad_params = {key: value for key, value in
                      self.mesh_params.items()}
        bad_params['cell_length'] = 'c'
        bad_cell = Cell(self.index, bad_params)
        
    @raises(KeyError)
    def test_length_missing(self):
        bad_params = {key: value for key, value in
                      self.mesh_params.items()}
        del bad_params['cell_length']
        bad_cell = Cell(self.index, bad_params)
