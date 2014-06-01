! gfortran -03 -I$SHTOOLSMODPATH_GCC -L$SHTOOLSLIBPATH_GCC -L$FFTW3_LIB_DIR_GCC -lSHTOOLS2.8 -llapack -lblas -lfftw3 -lm TEST_RESULTS.f90 -o TEST_RESULTS


PROGRAM TEST_RESULTS
	USE SHTOOLS
	
	REAL*8   ::  DAT(181,361)
	REAL*8   ::  CILM(2,10,10)
	integer  ::  SH_LMAX, I, J
	
	OPEN(34, FILE="/Users/dchandan/Volumes/Elysium/Data1/CommonData/Dynamic_topography/C2E5/dyntopoC2E5_l1-22.txt", STATUS="OLD")
	DO I= 1,181
		READ(34, *) DAT(I,:)
	ENDDO
	CLOSE(34)
	
	OPEN(34, FILE="FieldRead.txt", STATUS="REPLACE")
	DO I=1,181
		WRITE(34, *) DAT(I,:)
	ENDDO
	CLOSE(34)
	
	OPEN(34, FILE="Coeffs.txt", STATUS="REPLACE")
	CALL SHExpandDH(DAT(1:180,1:360), 180, CILM, LMAX=SH_LMAX, NORM=4, SAMPLING=2, LMAX_CALC=9)
	DO I =1,10
		WRITE(34, '(10F6.3)') CILM(1,I,:)
	ENDDO
	CLOSE(34)
	
END PROGRAM TEST_RESULTS