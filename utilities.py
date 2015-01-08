

class PySHTOOLSError(Exception):
    pass


class InvalidNormalization(PySHTOOLSError):
    """ This error is raised when the normalization of
    Legendre polynomials is not correct.
    """
    def __init__(self):
        pass


class InvalidPhaseFactor(PySHTOOLSError):
    """ This error is raised when the CSPHASE factor does not
    have the correct value.
    """
    def __init__(self):
        pass


class InvalidSampling(PySHTOOLSError):
    def __init__(self):
        pass


class GridError(PySHTOOLSError):
    def __init__(self):
        pass


class CILMShapeError(PySHTOOLSError):
    def __init__(self):
        pass


class DegreeError(PySHTOOLSError):
    def __init__(self):
        pass
    

def CheckNorm(n):
    if (n > 4) or (n < 1):
        print("ERROR !!!")
        print("Normalization should be: 1 (geodesy), 2 (Schmidt), 3 (unnormalized), or 4 (orthonormalized)")
        print("Received NORM = {0}".format(n))
        raise InvalidNormalization()


def CheckCSphase(l):
    if not ((l == 1) or (l == -1)):
        print("ERROR !!!")
        print("CSPHASE must be 1 (exclude) or -1 (include)")
        print("Received CSPHASE = {0}".format(l))
        raise InvalidPhaseFactor()


def CheckSampling1(n):
    if not ((n == 1) or (n == 2)):
        print("ERROR !!!")
        print("SAMPLING must be 1 (NxN grid) or 2 (Nx2N grid)")
        print("Received SAMPLING = {0}".format(n))
        raise InvalidSampling()


def CheckSampling2(n, tple):
    if n == 1:
        if not (tple[1] == tple[0]):
            print("ERROR !!!")
            print("For SAMPLING=1 grid must be of size NxN")
            print("Received size: {0}x{1}".format(tple[0], tple[1]))
            raise GridError()
    else:
        if not (tple[1] == tple[0]*2):
            print("ERROR !!!")
            print("For SAMPLING=2 grid must be of size Nx2N")
            print("Received size: {0}x{1}".format(tple[0], tple[1]))
            raise GridError()


def CheckCILM(s, l=None):
    """
    Checks that the CILM array has proper dimensions. 
    ARGUMENTS
        s = 3-tuple that is the shape of the CILM array
        l = spherical harmonic degree. 
            If provided, CILM is checked to be of shape (2,l+1,l+1)
            otherwise, CILM is check to be of shape (2,(2d square array))
    """
    if l is not None:
        if not ((s[0] == 2) and (s[1] == l+1) and (s[2] == l+1)):
            print("ERROR !!!")
            print("CILM array should be of shape (2, l+1, l+1) where l is {0}".format(l))
            print("Received shape is ({0}, {1}, {2})".format(s[0],s[1],s[2]))
            raise CILMShapeError()
    else:
        if not ((s[0] == 2) and (s[1] == s[2])):
            print("ERROR !!!")
            print("CILM array should be of shape (2, l+1, l+1)")
            print("Received shape is ({0}, {1}, {2})".format(s[0],s[1],s[2]))
            raise CILMShapeError()


def CheckLmaxNotNegative(l):
    if (l < 0):
        print("ERROR !!!")
        print("LMAX cannot be less than 0. Value received: {0}".format(l))
        raise DegreeError()

def CheckPlmZ(z):
    if (abs(z) > 1.0):
        print("ERROR !!!")
        print("Absolute value of Z argument to Legendre Polynomials must be <= 1.0")
        raise PySHTOOLSError()

# if __name__ == "__main__":
    # CheckSampling2(3, (180,360))
