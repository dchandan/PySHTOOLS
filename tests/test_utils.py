
import numpy as np



def periodic2D(m, n, k1, k2):
    grid = np.zeros((m,n))
    
    x = np.sin(np.linspace(0, 2*np.pi, n)*k1)
    y = np.sin(np.linspace(0, np.pi, m)*k2)
    
    for i in range(m):
        grid[i,:] = x
    for i in range(n):
        grid[:,i] = grid[:,i] + y
    return grid

