import numpy as np
import scipy.sparse as sps 
from scipy.sparse.linalg import spsolve

def fempoi2d(h, domain_length, t):
    """Solves 2D Poisson's equation for elements of size h 
    and bilinear basis equations"""
    # Define Global Matrices
    n_x = int(domain_length/h) # number of mesh cells in each direction
    nodes_x = int(n_x + 1) # number of node points in each direction
    nodes = int(nodes_x**2) # number of total node points
    n = int(n_x**2) # number of mesh cells
    A = sps.lil_matrix((nodes, nodes))
    b = sps.lil_matrix((nodes, 1))

    # Find Local Basis Functions
    for cell in range(n):
        cell_i = cell//int(np.sqrt(n))
        cell_j = cell%int(np.sqrt(n))
        
        
        # Retrieve Values for Nodes of Cell
        x_lower = t[cell_i, cell_j][1][0]
        x_upper = t[cell_i, cell_j][1][1]
        y_lower = t[cell_i, cell_j][1][2]
        y_upper = t[cell_i, cell_j][1][3]
        
        V = np.eye(4)
        C = np.array([[1, x_lower, y_lower, x_lower*y_lower],
                      [1, x_lower, y_upper, x_lower*y_upper],
                      [1, x_upper, y_lower, x_upper*y_lower],
                      [1, x_upper, y_upper, x_upper*y_upper]])
        
        # Solve for Coefficients
        c1 = np.linalg.inv(C).dot(V[0])
        c2 = np.linalg.inv(C).dot(V[1])
        c3 = np.linalg.inv(C).dot(V[2])
        c4 = np.linalg.inv(C).dot(V[3])
        c = np.stack((c1, c2, c3, c4), axis=-1)
        
        
        # Assemble Elementary Matrix
        # A_ab = area(partial_ax*partial_bx + partial_ay*partial_by)
        area = h**2
        elem_A = np.zeros((4, 4))
        elem_b = np.zeros(4)
        for ii in range(4):
            for jj in range(4):
                partial_ax = c[1, ii]+c[3, ii]*C[ii, 2]
                partial_bx = c[1, jj]+c[3, jj]*C[jj, 2]
                partial_ay = c[2, ii]+c[3, ii]*C[ii, 1]
                partial_by = c[2, jj]+c[3, jj]*C[jj, 1]
                elem_A[ii,jj] = area*(partial_ax*partial_bx + partial_ay*partial_by)
                elem_b[ii]=area/4.0*(4*c[0, ii] + h*(2*c[1, ii]+2*c[2, ii]+h*c[3, ii]))
        
        # Assemble Global Matrix
        # Global to Local node mapping
        """Nodes are flattened by row with row x=0 going first.
        That means node x0,y0 in cell0 has a global index of 0
        and node x1,y1 in cell0 has a global index of nodes_x+1"""
        mapping = np.array([cell+cell_i, cell+cell_i+1, 
                            cell+cell_i+n_x+1, cell+cell_i+n_x+2])
        xx = 0
        for ii in mapping:
            yy = 0
            for jj in mapping:
                A[ii, jj] += elem_A[xx, yy]
                b[ii, 0] += elem_b[xx]
                yy += 1
            xx += 1

    # Enforce Dirichlet Boundary Conditions u=1 on boundary
    # Identify Boundary Nodes
    boundary_nodes = np.zeros(4*n_x, dtype='int')
    boundary_nodes[:nodes_x] = np.arange(nodes_x) # left
    boundary_nodes[nodes_x:2*n_x] = np.arange(nodes_x, nodes-nodes_x, nodes_x) # bottom 
    boundary_nodes[2*n_x:3*n_x - 1] = np.arange(2*n_x + 1, nodes-nodes_x, nodes_x) # top
    boundary_nodes[3*n_x - 1:] = np.arange(nodes)[nodes - nodes_x:] # right
    
    for i in boundary_nodes:
        A[i, i] = 1
        b[i] = 0 #vacuum boundary
        if 0 < i+1 < nodes:
            A[i+1, i] = 0
            A[i-1, i] = 0
        elif i == 0:
            A[i+1, i] = 0
        else:
            A[i-1, i] = 0

    A = sps.csr_matrix(A)
    b = sps.csr_matrix(b)
    
    u = spsolve(A, b)
    u = u.reshape(nodes_x, nodes_x)
    return u