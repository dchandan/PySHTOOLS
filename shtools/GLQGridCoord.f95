subroutine GLQGridCoord(latglq, longlq, lmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Given a maximum spherical harmonic degree lmax, this routine 
!	will determine the latitude and longitude coordinates associated with 
!	grids that are used in the Gauss-Legendre quadratue spherical harmonic 
!	expansion routines. The coordinates are output in DEGREES.
!
!	Calling Parameters
!		IN	
!			lmax		Maximum spherical harmonic degree of the expansion.
!		OUT	
!			latglq		Array of latitude points used in Gauss-Legendre grids, in degrees.
!			longlq		Array of longitude points used in Gauss-Legendre grids, in degrees.
!			nlat, nlong	Number of latidude and longitude points.
!
!	Dependencies:		NGLQSH, PreGLQ
!
!	Written by Mark Wieczorek, 2004
!	Completely rewritten June 6, 2006.
!
!	Copyright (c) 2005-2006 Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: NGLQSH, PreGLQ

	implicit none
	integer, intent(in)  ::	lmax
	real*8, intent(out)  ::	latglq(lmax+1), longlq(2*lmax+1)
	
	real*8               ::	pi, upper, lower, zero(lmax+1), w(lmax+1)
	integer              ::	i, nlat, nlong
	
	pi    = acos(-1.0d0)

	nlat  = NGLQSH(lmax)
	nlong = 2*lmax +1
	
	upper = 1.0d0
	lower = -1.0d0	
						
	call PreGLQ(lower, upper, nlat, zero, w)	! Determine Gauss Points and Weights.

	do i=1, nlong
		longlq(i) = 360.0d0*(i-1)/nlong
	enddo

	do i=1, nlat
		latglq(i) = asin(zero(i))*180.0d0/pi
	enddo

end subroutine GLQGridCoord

