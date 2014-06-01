from _shtools import plon, plmon, plmbar, plegendrea
from _shtools import plegendre, plmschmidt, plschmidt

def PlON(lmax, z):
    return plon(lmax, z)


def PlmON(lmax, z, csphase=1):
    return plmon(lmax, z, csphase=csphase, cnorm=0)

def PlmBar(lmax, z, csphase=1):
    return plmbar(lmax, z, csphase=csphase, cnorm=0)

def PLegendreA(lmax, z, csphase=1):
    return plegendrea(lmax, z, csphase=csphase)

def PLegendre(lmax, z):
    return plegendre(lmax, z)

def PlmSchmidt(lmax, z, csphase=1):
    return plmschmidt(lmax, z, csphase=csphase, cnorm=0)

def PlSchmidt(lmax, z):
    return plschmidt(lmax, z)