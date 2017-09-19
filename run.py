from input import problem
import build_cells
import fempoi

DOMAIN_LOWER = build_cells.DOMAIN_LOWER
DOMAIN_UPPER = build_cells.DOMAIN_UPPER

CELLS = problem["mesh_cells"]
DOMAIN_LENGTH = DOMAIN_UPPER - DOMAIN_LOWER
CELL_LENGTH = DOMAIN_LENGTH/CELLS

def run():
  u = fempoi.fempoi(CELL_LENGTH, DOMAIN_LENGTH)
  print u

run()


