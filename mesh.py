import numpy as np

def mesh_gen(cells):
  xx = np.linspace(0, 10, cells+1)#.reshape(cells, 1)
  XY = np.zeros((cells+1, cells+1, 2))
  for i in range(cells+1):
    for j in range(cells+1):
        XY[i, j, 0] = xx[i]
        XY[i, j, 1] = xx[j]
  return XY

class Cell():
  """ A single cell in the mesh, holds location and material data """

  def __init__(self):
    pass
