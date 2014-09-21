subroutine PreGLQ(x1, x2, n, zero, w)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will find the zeros and weights that are
!	used in Gauss-Legendre quadrature routines. (Based on routines
!	in Numerical Recipes).
!
!	Calling Parameters:
!		IN
!			x1: Lower bound of integration.
!			x2:	Upper bound of integration.
!			n:	Number of points used in the quadrature. n points
!				will integrate a polynomial of degree 2n-1 exactly.
!		OUT
!			zero:	Array of n Gauss points, which correspond to the zeros
!				of P(n,0).
!			w:	Array of n weights used in the quadrature.
!
!
!	Note 
!		1.	If EPS is less than what is defined, then the do 
!			loop for finding the roots might never terminate for some
!			values of lmax. If the algorithm doesn't converge, consider
!			increasing itermax, or decreasing eps.
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek 2003
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IMPLICIT NONE
	REAL*8, INTENT(IN)  :: 	X1, X2
	REAL*8, INTENT(OUT) ::	ZERO(N), W(N)
	INTEGER, INTENT(IN) ::	N
	
	INTEGER             ::	I, J, M, ITER
	INTEGER, PARAMETER  ::	ITERMAX = 1000
	REAL*8, PARAMETER   ::	EPS     = 1.0D-15
	REAL*8              ::	P1, P2, P3, PP, Z, Z1, XM, XU
	
	
	ZERO = 0.0D0
	W    = 0.0D0
	
	! The roots are symmetric in the interval, so we only have to find half of them. 
	M    = (N+1)/2
	
	XM   = (X2 + X1) / 2.0D0	! midpoint of integration
	XU   = (X2 - X1) / 2.0D0	! Scaling factor between interval of integration, and that of 
				! the -1 to 1 interval for the Gauss-Legendre interval
	
	DO I = 1,M 
		ITER = 0
		Z = COS(3.141592654d0 * (I - .25D0)/(N+.5D0))	! Approximation for the ith root
							
		! Find the true value using newtons method
		DO
			ITER = ITER +1
			
			P1 = 1.0D0
			P2 = 0.0D0
			
			! determine the legendre polynomial evaluated at z (p1) using recurrence relationships
			DO J=1, N 	
				P3 = P2
				P2 = P1
				P1 = (DBLE(2*J-1)*Z*P2-DBLE(J-1)*P3)/DBLE(J)
			ENDDO
			
			! This is the derivative of the legendre polynomial using recurrence relationships
			PP=DBLE(N)*(Z*P1-P2)/(Z*Z-1.0D0)
		
			Z1 = Z
			Z  = Z1-P1/PP 
		
			IF (ABS(Z-Z1) <= EPS) EXIT
			
			IF (ITER > ITERMAX) THEN
				PRINT*, "ERROR --- PreGLQ"
				PRINT*, "Root Finding of PreGLQ not converging."
				PRINT*, "m , n = ", m, n
				STOP
			ENDIF			
		ENDDO
		
		ZERO(I) = XM + XU*Z
		ZERO(N+1-I) = XM - XU*Z
		W(I)= 2.0D0 * XU / ((1.0D0-Z*Z)*PP*PP)
		W(N+1-I)=W(I) 		
	ENDDO
END SUBROUTINE PreGLQ



INTEGER FUNCTION NGLQ(DEGREE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	For a polynomial of order degree, this simple function
!	will determine how many gauss-legendre quadrature points
!	are needed in order to integrate the function exactly.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IMPLICIT NONE
	INTEGER, INTENT(IN) ::	DEGREE
	
	NGLQ = CEILING((DEGREE+1.0D0)/2.0D0) 	
	
END FUNCTION NGLQ


INTEGER FUNCTION NGLQSH(DEGREE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function returns the number of gauss-legendre points that
!	are needed to exactly integrate a spherical harmonic field of 
!	Lmax = degree
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IMPLICIT NONE
	INTEGER, INTENT(IN) ::	DEGREE
	
	NGLQSH = DEGREE + 1.0D0

END FUNCTION NGLQSH


INTEGER FUNCTION NGLQSHN(DEGREE, N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function returns the number of gauss-legendre points that
!	are needed to exactly integrate a spherical harmonic field of 
!	Lmax = degree raised to the nth power. Here, the maximum degree
!	of the integrand is n*lmax + lmax, or (n+1)*lmax
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IMPLICIT NONE
	INTEGER, INTENT(IN) ::	DEGREE, N
	
	NGLQSHN = CEILING( ((N+1.0D0)*DEGREE + 1.0D0)/2.0D0)

END FUNCTION NGLQSHN

