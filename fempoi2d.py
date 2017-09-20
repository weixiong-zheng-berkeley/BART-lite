def fempoi2D(h, domain_length):
    """Solves 2D Poisson's equation for elements of size h 
    and bilinear basis equations"""
    # Define Global Matrices
    n_x = int(domain_length/h) # number of mesh cells in each direction
    nodes_x = int(n_x + 1) # number of node points in each direction
    nodes = int(nodes_x**2) # number of total node points
    n = int(n_x**2) # number of mesh cells
    A = np.zeros((nodes, nodes))
    b = np.zeros(nodes)

    # Find Local Basis Functions
    for cell in range(n):
        cell_i = cell//int(np.sqrt(n))
        cell_j = cell%int(np.sqrt(n))
        
        lower = 1.0/h
        upper = -1.0/h
        
        x = np.array([lower, lower, upper, upper]) 
        y = np.array([lower, upper, lower, upper]) 
        
        # Assemble Elementary Matrix
        # A_ab = area(partial_ax*partial_bx + partial_ay*partial_by)
        area = h**2
        elem_A = np.zeros((4, 4))
        elem_b = np.zeros(4)
        for ii in range(4):
            for jj in range(4):
                elem_A[ii,jj] = area*(x[ii]*x[jj] + y[ii]*y[jj])
        elem_b[:]=area/4.0
        
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
                b[ii] += elem_b[xx]
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
        b[i] = 1   
        if 0 < i+1 < nodes:
            A[i+1, i] = 0
            A[i-1, i] = 0
        elif i == 0:
            A[i+1, i] = 0
        else:
            A[i-1, i] = 0
            
    u = np.linalg.solve(A, b)
    u = u.reshape(nodes_x, nodes_x)
    return u