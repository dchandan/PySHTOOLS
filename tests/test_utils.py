
import numpy as np



def periodic2D(m, n, k1, k2):
    """
    Makes and returns a grid with 2D periodic wave pattern
    ARGUMENTS:
        m - The number of points in output grid along the vertical
        n - The number of points in output grid along the horizontal
        k1 - wavenumber in the vertical
        k2 - wavenumber in the horizontal
    RETURNS
        grid of shape (m,n)
    """
    grid = np.zeros((m,n))
    
    x = np.sin(np.linspace(0, 2*np.pi, n)*k1)
    y = np.sin(np.linspace(0, np.pi, m)*k2)
    
    for i in range(m):
        grid[i,:] = x
    for i in range(n):
        grid[:,i] = grid[:,i] + y
    return grid


if __name__ == "__main__":
    import matplotlib.pylab as plt
    plt.imshow(periodic2D(50, 100, 3, 4))
    plt.show()