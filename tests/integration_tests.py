from nose.tools import *
from mesh import Cell, Mesh
from material import mat_lib, mat_map
import numpy as np

class TestIntegration_Mesh_Materials:
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

    def test_mesh_init(self):
        mesh_cells = 4
        domain_upper = 10
        mesh = Mesh(mesh_cells, domain_upper, self.testmap)
        mesh_params = {'x_cell': 4,
                       'cell_length': 2.5}
        eq_(mesh._mesh_params, mesh_params)

    def test_mesh_cells_size(self):
        mesh_cells = 4
        domain_upper = 10
        mesh = Mesh(mesh_cells, domain_upper, self.testmap)
        eq_(len(mesh.cells()), mesh_cells**2)

    def test_mesh_props(self):
        mesh_cells = 4
        domain_upper = 10
        mesh = Mesh(mesh_cells, domain_upper, self.testmap)
        eq_(mesh.n_cell(), 16, "n_cell value")
        eq_(mesh.n_node(), 25, "n_node value")
        eq_(mesh.x_cell(), 4, "x_cell value")
        eq_(mesh.y_cell(), 4, "y_cell value")
        eq_(mesh.x_node(), 5, "x_node value")
        eq_(mesh.y_node(), 5, "y_node value")
        eq_(mesh.cell_length(), 2.5, "cell length value")

    def test_plotting(self):
        mesh_cells = 4
        domain_upper = 10
        mesh = Mesh(mesh_cells, domain_upper, self.testmap)
        mesh_params = {'x_cell': 4,
                       'cell_length': 2.5}
        x_ans   = [0, 2.5, 5.0, 7.5, 10.0,
               0, 2.5, 5.0, 7.5, 10.0,
               0, 2.5, 5.0, 7.5, 10.0,
               0, 2.5, 5.0, 7.5, 10.0,
               0, 2.5, 5.0, 7.5, 10.0]

        y_ans   = [0, 0, 0, 0, 0,
               2.5, 2.5, 2.5, 2.5, 2.5,
               5.0, 5.0, 5.0, 5.0, 5.0,
               7.5, 7.5, 7.5, 7.5, 7.5,
               10.0, 10.0, 10.0, 10.0, 10.0]

        ans = [0.5, 1.0, 1.0, 1.0, 0.5,
               1.0, 2.0, 2.0, 2.0, 1.0,
               1.0, 2.0, 2.0, 2.0, 1.0,
               1.0, 2.0, 2.0, 2.0, 1.0,
               0.5, 1.0, 1.0, 1.0, 0.5]
        
        x, y, z = mesh.test_plot(plot=False)
        ok_(np.array_equal(x, x_ans), "x values")
        ok_(np.array_equal(y, y_ans), "y values")
        ok_(np.array_equal(ans, z), "z values")
        
    @raises(AssertionError)
    def test_mesh_cells_is_int(self):
        mesh_cells = 4.5
        domain_upper = 10
        mesh = Mesh(mesh_cells, domain_upper, self.testmap)
        mesh_params = {'x_cell': 4,
                       'cell_length': 2.5}

    @raises(AssertionError)
    def test_mesh_size_to_map(self):
        """ Mesh size must be the same as the matmap size """
        mesh_cells = 4
        domain_upper = 8
        mesh = Mesh(mesh_cells, domain_upper, self.testmap)
                    

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
        
