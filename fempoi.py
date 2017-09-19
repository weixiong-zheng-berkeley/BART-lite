from __future__ import division
import numpy as np 


def fempoi(h, domain_length):
    """Solves Poisson's equation for elements of size h 
    and bilinear basis equations"""
    n = int(domain_length/h) # number of mesh cells
    A = np.zeros((n+1, n+1))
    b = np.zeros(n+1)
    
    # Assembly
    A_k = n*np.array([[1, -1],[-1, 1]])
    b_k = h/2*np.array([1 , 1])
    for i in range(n):
        A[i:i+2, i:i+2] += A_k
        b[i:i+2] += b_k
        
    # Direchlet Boundary Conditions
    A[0, 0] = 1
    b[0] = 0 
    A[-1, -1] = 1
    b[-1] = 1
    
    # System is singular so eliminate entries above/below 
    # Direchlet variables
    b[1] -= A[1, 0]*b[0]
    A[1, 0] = 0
    b[-2] -= A[-2, -1]*b[-1]
    A[-2, -1] = 0
    
    u = np.linalg.solve(A, b)
    return u