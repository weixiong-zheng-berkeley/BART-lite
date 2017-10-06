from __future__ import division
import matplotlib.pyplot as plt
import build_cells
import fempoi2d
import material
import mesh

#Specify problem here:
from input_kaist_mox_1 import problem

# Build Material Library

try:
  MAT_LIB = material.mat_lib(n_grps = problem['groups'],
                             files = problem['materials'],
                             tr_scatt = problem['tr_scatt'])
except KeyError:
  MAT_LIB = material.mat_lib(n_grps = problem['groups'],
                             files = problem['materials'])

# Build Material Mapping
  
MAT_MAP = material.mat_map(lib = MAT_LIB, layout = problem['layout'],
                           layout_dict = problem['layout_dict'],
                           x_max=problem['domain_upper'],
                           n=problem['mesh_cells'])

# Build Mesh

MESH = mesh.Mesh(problem['mesh_cells'], problem['domain_upper'], MAT_MAP)

def run():
  u = fempoi2d.fempoi2d(CELL_LENGTH, DOMAIN_LENGTH, DATA)
  x = mesh_gen(CELLS)[:, :, 0]
  y = mesh_gen(CELLS)[:, :, 1]
  cset1 = plt.contourf(x, y, u, 10)
  plt.colorbar()
  plt.contour(x, y, u, cset1.levels, hold='on', colors='k')
  plt.axis('equal')
  plt.show()

#run()


