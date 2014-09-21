REAL*8 FUNCTION MakeGridPoint(cilm, lmax, lat, longitude, csphase, norm, dealloc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will determine the value at a given latitude and 
!	longitude corresponding to the given set of spherical harmonics.
!	Latitude and Longitude are assumed to be in DEGREES!
!
!	Calling Parameters:
!		IN
!			cilm:		Spherical harmonic coefficients, with dimensions
!					(2, lmax+1, lmax+1).
!			lmax:		Maximum degree used in the expansion.
!			lat:		latitude (degrees).
!			long:		longitude(degrees).
!		OPTIONAL (IN)
!			norm		Spherical harmonic normalization:
!						(1) "geodesy" (default)
!						(2) Schmidt
!						(3) unnormalized
!						(4) orthonormalized
!			csphase:		1: Do not include the phase factor of (-1)^m
!						-1: Apply the phase factor of (-1)^m.
!			dealloc		If (1) Deallocate saved memory in Legendre function routines.
!					Default (0) is not to deallocate memory.
!
!
!	Dependencies:		PlmBar, PlBar, PlmSchmidt, PlmON, CSPHASE_DEFAULT
!
!	Notes:
!		1.	If lmax is greater than the the maximum spherical harmonic
!			degree of the input file, then this file will be ZERO PADDED!
!			(i.e., those degrees after lmax are assumed to be zero).
!
!	Written by Mark Wieczorek (June 2004)
!
!	August 8, 2012. Modified to precompute sin and cosine terms using a recusion relationship. Also, the default has been modified to 
!			save variables used in the Legendre function routines. To deallocate this space, specify
!			the new optional parameter DEALLOC to be 1.
!
!	Copyright (c) 2005-2012, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: PlmBar, PLegendreA, PlmSchmidt, PlmON

	implicit none
	
	real*8, intent(in)            :: cilm(2,lmax+1,lmax+1), lat, longitude
	integer, intent(in)           :: lmax, csphase
	integer, intent(in), optional :: norm, dealloc
	
	real*8                        :: pi, x, expand, lon
	integer                       :: index, l, m, l1, m1, astat(3)
	real*8, allocatable           :: pl(:), cosm(:), sinm(:)
	
	
 	allocate(pl((lmax+1)*(lmax+2)/2), stat = astat(1))
 	allocate(cosm(lmax+1), stat = astat(2))
 	allocate(sinm(lmax+1), stat = astat(3))
 	
	if (sum(astat(1:3)) /= 0) then
		print*, "Error --- MakeGridPoint"
		print*, "Cannot allocate memory for arrays PL, MCOS and MSIN", astat(1), astat(2), astat(3)
		stop
	endif 

	pi  = acos(-1.0d0)
	x   = sin(lat*pi/180.0d0)
	lon = longitude * pi/180.0d0
	
	
	if (present(norm)) then
		if (norm == 1) call PlmBar(pl, lmax, x, csphase = csphase)
		if (norm == 2) call PlmSchmidt(pl, lmax, x, csphase = csphase)
		if (norm == 3) call PLegendreA(pl, lmax, x, csphase = csphase)
		if (norm == 4) call PlmON(pl, lmax, x, csphase = csphase)
	else
		call PlmBar(pl, lmax, x, csphase = csphase)
	endif
	
	expand = 0.0d0
	
	! Precompute sines and cosines. Use multiple angle identity to minimize number of calls to SIN and COS.
	sinm(1) = 0.0d0
	cosm(1) = 1.0d0
	if (lmax > 0) then
		sinm(2) = sin(lon)
		cosm(2) = cos(lon)
	endif
	do m=2, lmax, 1
		m1 = m+1
		sinm(m1) = 2 * sinm(m)*cosm(2) - sinm(m-1)
		cosm(m1) = 2 * cosm(m)*cosm(2) - cosm(m-1)
	enddo

	do l = lmax, 0, -1
		l1 = l+1
		index = (l+1)*l/2 + 1
		expand = expand + cilm(1,l1,1) * pl(index)

		do m = 1, l, 1
			m1 = m+1
			index = index + 1
			expand = expand + ( cilm(1,l1,m1) * cosm(m1) + &
				cilm(2,l1,m1) * sinm(m1) ) * pl(index)
		enddo
	enddo
	
	MakeGridPoint = expand
	
	! deallocate memory
	if (present(dealloc)) then
		if (dealloc == 1) then
			if (present(norm)) then
				if (norm == 1) call PlmBar(pl, -1, x, csphase = csphase)
				if (norm == 2) call PlmSchmidt(pl, -1, x, csphase = csphase)
				if (norm == 4) call PlmON(pl, -1, x, csphase = csphase)
			else
				call PlmBar(pl, -1, x, csphase = csphase)
			endif
		endif
	endif
	
	deallocate(pl)
	deallocate(cosm)
	deallocate(sinm)
	
END FUNCTION MakeGridPoint

