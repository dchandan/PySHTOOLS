
SOURCES = DHaj.f95 \
SHExpandDH.f95     \
PlmIndex.f95       \
PlON.f95           \
PlmON.f95          \
PlmBar.f95         \
PlBar.f95          \
PLegendre.f95      \
PLegendreA.f95     \
PlSchmidt.f95      \
PlmSchmidt.f95     \
MakeGrid2D.f95     \
ReadWisdom.f95     \
PreGLQ.f95         \
PreCompute.f95     \
SHExpandLSQ.f95    \
MakeGridPoint.f95  \
MakeGridDH.f95     \
GLQGridCoord.f95   \
SHExpandDHC.f95    \
MakeGridDHC.f95    \
SHConvertCoef.f95

MODS = FFTW3.o SHTOOLS.o

all: $(MODS)
	f2py -c -m _shtools $(SOURCES) -L$(FFTW3_LIB_DIR_GCC) -I$(FFTW3_INC_DIR_GCC) -lfftw3 -llapack -lblas
	mv _shtools.so ../
	


FFTW3.o: FFTW3.f95
	gfortran -c FFTW3.f95

SHTOOLS.o: SHTOOLS.f95
	gfortran -c SHTOOLS.f95

clean:
	rm -f *.o *.mod *.so ../*.so