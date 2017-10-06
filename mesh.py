import numpy as np

def mesh_gen(cells):
  xx = np.linspace(0, 10, cells+1)#.reshape(cells, 1)
  XY = np.zeros((cells+1, cells+1, 2))
  for i in range(cells+1):
    for j in range(cells+1):
        XY[i, j, 0] = xx[i]
        XY[i, j, 1] = xx[j]
  return XY

class Cell(object):
  """ A single cell in the mesh, holds location and material data """

  def __init__(self, index, mesh_params, bound=None, mat_lib=None):
    """ Cell constructor, give index in a tuple (i,j) """
    
    # Constructor validations
    assert isinstance(index, tuple), "Index must be a tuple"
    assert len(index) == 2, "Index must be a length 2 tuple"
    
    try:
      assert mesh_params['cell_length'] > 0, "cell_length must be greater than 0"
    except KeyError:
      raise KeyError("Missing 'cell_length' parameter in mesh_params")
    
    self._index  = index

    try:
      self._length = float(mesh_params['cell_length'])
      self._area   = np.power(self._length, 2)
    except ValueError:
      raise TypeError("cell_length parameter must be a number")

    # Calculate global_idx
    x_node = mesh_params['x_cell'] + 1
    i,j = index[0], index[1]
    self._global_idx = [x_node*i + j,
                        x_node*i + j + 1,
                        x_node*(i + 1) + j,
                        x_node*(i + 1) + j + 1]
    
    # Get material properties
    if mat_lib:
      self.__material_props__(mat_lib)

    # Determine if on a boundary
    

  # UTILITY FUNCTIONS ================================================

  def __material_props__(self, mat_lib):
    pass
    
  # Getters ==========================================================
    
  def area(self):
    return self._area
   
  def global_idx(self):
    """ Returns global index, a list of the node indices """
    return self._global_idx

  def length(self):
    return self._length
