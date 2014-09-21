from _shtools import plon, plmon, plmbar, plegendrea
from _shtools import plegendre, plmschmidt, plschmidt
from utilities import *

def PlON(lmax, z):
    CheckLmaxNotNegative(lmax)
    CheckPlmZ(z)
    return plon(lmax, z)


def PlmON(lmax, z, csphase=1):
    CheckLmaxNotNegative(lmax)
    CheckCSphase(csphase)
    CheckPlmZ(z)
    return plmon(lmax, z, csphase=csphase, cnorm=0)

def PlmBar(lmax, z, csphase=1):
    CheckLmaxNotNegative(lmax)
    CheckCSphase(csphase)
    CheckPlmZ(z)
    return plmbar(lmax, z, csphase=csphase, cnorm=0)

def PLegendreA(lmax, z, csphase=1):
    CheckLmaxNotNegative(lmax)
    CheckCSphase(csphase)
    CheckPlmZ(z)
    return plegendrea(lmax, z, csphase=csphase)

def PLegendre(lmax, z):
    CheckLmaxNotNegative(lmax)
    CheckPlmZ(z)    
    return plegendre(lmax, z)

def PlmSchmidt(lmax, z, csphase=1):
    CheckLmaxNotNegative(lmax)
    CheckCSphase(csphase)
    CheckPlmZ(z)
    return plmschmidt(lmax, z, csphase=csphase, cnorm=0)

def PlSchmidt(lmax, z):
    CheckLmaxNotNegative(lmax)
    CheckPlmZ(z)
    return plschmidt(lmax, z)