from nose.tools import *
from mesh import Cell
from material import mat_lib, mat_map
import numpy as np

class TestIntegration_Cell_Materials:
    # Tests to verify integration of materials and cells
    
    @classmethod
    def setup_class(cls):
        testData_loc = './tests/testData/materials/'
        files = ['test_mat.xml', 'test_mat2.xml']
        filelocs = [testData_loc + f for f in files]
        
        # Materials
        cls.lib = mat_lib(n_grps = 2, files = filelocs)
        cls.mat_dict = {'1': 'test_mat', '2': 'test_mat2'}
        layout = """ 1 1 1 1
                          1 2 2 1
                          1 2 2 1
                          1 1 1 1 """
        cls.testmap = mat_map(lib=cls.lib, layout=layout,
                              layout_dict = cls.mat_dict, x_max = 10, n=20)

        # Cell:
        cls.mesh_params = {'x_cell': 8,
                           'cell_length': 1.25}

    def test_material_at_locations(self):
        for param in [((1,1), 'test_mat'),
                      ((2,2), 'test_mat2'),
                      ((5,4), 'test_mat2'),
                      ((3,0), 'test_mat'),
                      ((7,6), 'test_mat')]:
            cell = Cell(param[0], self.mesh_params, self.testmap)
            eq_(cell.get('id'), param[1],
                "cell at " + str(param[0]) + "correct material")

    @raises(AssertionError)
    def test_misaligned_mesh_and_map(self):
        """ Mesh must be the same size as the mat mapping """
        mesh_params = {'x_cell': 8,
                       'cell_length': 2}
        cell = Cell((1,1), mesh_params, self.testmap)
        
