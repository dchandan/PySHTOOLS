subroutine SHrtoc(rcilm, l1, ccilm, degmax, convention, switchcs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will convert the real "geodesy 4-pi" spherical harmonic coefficients into 
!	complex form with either a 4 pi or unit normalization. The complex 
!	coefficients are only calculated for positive m's (the negative m's are given by 
!
!		c_l,m* = (-1)^m c_l,-m
!
!	This sign convention is designed to be used with the spherical haromic rotation routines
!	taken from Mark Simons' code. If degmax is not specified, then the maximum degree of the 
!	conversion is taken from the size of the input arrays.
!
!
!	Calling Parameters
!		IN
!			rcilm		Real "geodesy" spherical harmonic coefficients with dimensions
!					(2, lmax+1, lmax+1).
!		OUT
!			ccilm		Complex unity-normalized spherical harmonic coefficients, dimensioned
!					as (2, lmax+1, lmax+1). The first index corresponds to the real and
!					complex coefficients, respectively.
!		OPTIONAL
!			degmax		Maximum degree of conversion to be performed.
!			convention	1=output 4-pi normalized coefficients
!					2=output  Varshalovich et al. normalized coefficients
!			switchcs	If 1, Change between different Condon-Shortley phase convenctions.
!					If 0, use consistent phase convention.
!
!		Dependencies:	None
!			
!
!	Written by Mark Wieczorek (August 2003)
!	June 4 2006: Added optional argument to change between Condon-Shortley phase conventions.
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in)           ::    l1, degmax
	real*8, intent(in)            :: 	rcilm(2,l1+1,l1+1)
	real*8, intent(out)           ::	ccilm(2,degmax+1,degmax+1) 
	integer, intent(in), optional ::	convention, switchcs
	integer                       :: 	lmax, l, m
	real*8                        ::	pi
	
	if (present(switchcs)) then
		if (switchcs /= 1 .and. switchcs /=0) then
			print*, "Error --- SHrtoc"
			print*, "switchcs must be equal to either 0 (keep same convention) of 1 (change Condon-Shortley phase)"
			print*, "Input value is ", switchcs
			stop
		endif
	endif
	
	if (present(convention) ) then
		if (convention /=1 .and. convention /=2) then
			print*, "Error --- SHrtoc"
			print*, "CONVENTION must be 1 or 2."
			print*, "Input valuse is ", convention
			stop
		endif
	endif

	lmax = degmax
	if (size(rcilm(:,1,1)) < 2) then
		print*, "Error --- SHrtoc"
		print*, "RCILM must be dimensioned as (2,*,*)."
		print*, "Input array is dimensioned as ",  size(rcilm(:,1,1)), size(rcilm(1,:,1)), size(rcilm(1,1,:))
		stop
	elseif (size(ccilm(:,1,1)) < 2) then
		print*, "Error --- SHrtoc"
		print*, "CCILM must be dimensioned as (2,*,*)."
		print*, "Input array is dimensioned as ",  size(ccilm(:,1,1)), size(ccilm(1,:,1)), size(ccilm(1,1,:))
		stop
	endif

	pi = acos(-1.0d0)
	ccilm = 0.0d0
		
	do l = 0, lmax, 1
		if (convention == 2) then
			ccilm(1,l+1, 1) = sqrt(4.0d0*pi)*rcilm(1,l+1,1)
			ccilm(2,l+1, 1) = 0.0d0
		
			do m=1, l, 1
				if (present(switchcs)) then
					if (switchcs == 1) then
						ccilm(1, l+1, m+1) = sqrt(2.0d0*pi) * rcilm(1,l+1,m+1) * (-1.0d0)**m
						ccilm(2, l+1, m+1) = -sqrt(2.0d0*pi) * rcilm(2,l+1,m+1) * (-1.0d0)**m
					else
						ccilm(1, l+1, m+1) = sqrt(2.0d0*pi) * rcilm(1,l+1,m+1) 
						ccilm(2, l+1, m+1) = -sqrt(2.0d0*pi) * rcilm(2,l+1,m+1)
					endif
				else
					ccilm(1, l+1, m+1) = sqrt(2.0d0*pi) * rcilm(1,l+1,m+1) 
					ccilm(2, l+1, m+1) = -sqrt(2.0d0*pi) * rcilm(2,l+1,m+1)
				endif
			enddo			
		else
			ccilm(1,l+1, 1) = rcilm(1,l+1,1)
			ccilm(2,l+1, 1) = 0.0d0
			
			do m=1, l, 1
				if (present(switchcs)) then
					if (switchcs == 1) then
						ccilm(1, l+1, m+1) = rcilm(1,l+1,m+1)/sqrt(2.0d0) * (-1.0d0)**m
						ccilm(2, l+1, m+1) = -rcilm(2,l+1,m+1)/sqrt(2.0d0) * (-1.0d0)**m
					else
						ccilm(1, l+1, m+1) = rcilm(1,l+1,m+1)/sqrt(2.0d0)
						ccilm(2, l+1, m+1) = -rcilm(2,l+1,m+1)/sqrt(2.0d0)
					endif
				else
					ccilm(1, l+1, m+1) = rcilm(1,l+1,m+1)/sqrt(2.0d0)
					ccilm(2, l+1, m+1) = -rcilm(2,l+1,m+1)/sqrt(2.0d0)
				endif
			enddo
		endif			
	enddo
	
end subroutine SHrtoc


subroutine SHctor(ccilm, l1, rcilm, degmax, convention, switchcs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will convert either "geodesy 4-pi" or "Varshalovich et al." complex spherical 
!	harmonic coefficients into real "geodesy 4-pi" spherical harmonic coefficients.
!
!	If degmax is not specified, then the maximum degree of the 
!	conversion is taken from the size of the input arrays.
!
!
!	Calling Parameters
!		IN
!			ccilm		Complex unity-normalized spherical harmonic coefficients, dimensioned
!					as (2, lmax+1, lmax+1). The first index corresponds to the real and
!					complex coefficients, respectively.
!		OUT
!			rcilm		Real "geodesy" spherical harmonic coefficients with dimensions
!					(2, lmax+1, lmax+1).
!		OPTIONAL
!			degmax		Maximum degree of conversion to be performed.
!			convention	1=input coefficients are 4-pi normalized
!					2=input coefficients are Varshalovich et al. normalized.
!			switchcs	If 1, Change between different Condon-Shortley phase convenctions.
!					If 0, use consistent phase convention.
!
!		Dependencies:	None
!
!	Written by Mark Wieczorek 2004.
!	June 4 2006: Added option to change between Condon-Shortley phase conventions.
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) :: l1, degmax
	real*8, intent(in) :: ccilm(2,l1+1,l1+1)
	real*8, intent(out) :: rcilm(2,degmax+1,degmax+1)
	integer, intent(in), optional :: convention, switchcs
	integer             ::	lmax, l, m
	real*8              ::	pi
	
	if (present(switchcs)) then
		if (switchcs /= 1 .and. switchcs /=0) then
			print*, "Error --- SHrtoc"
			print*, "switchcs must be equal to either 0 (keep same convention) of 1 (change Condon-Shortley phase)"
			print*, "Input value is ", switchcs
			stop
		endif
	endif
	
	if (present(convention) ) then
		if (convention /=1 .and. convention /=2) then
			print*, "Error --- SHrtoc"
			print*, "CONVENTION must be 1 or 2."
			print*, "Input valuse is ", convention
			stop
		endif
	endif

	lmax = degmax

	if (size(rcilm(:,1,1)) < 2) then
		print*, "Error --- SHrtoc"
		print*, "RCILM must be dimensioned as (2,*,*)."
		print*, "Input array is dimensioned as ",  size(rcilm(:,1,1)), size(rcilm(1,:,1)), size(rcilm(1,1,:))
		stop
	elseif (size(ccilm(:,1,1)) < 2) then
		print*, "Error --- SHrtoc"
		print*, "CCILM must be dimensioned as (2,*,*)."
		print*, "Input array is dimensioned as ",  size(ccilm(:,1,1)), size(ccilm(1,:,1)), size(ccilm(1,1,:))
		stop
	endif
			
	pi = acos(-1.0d0)
	rcilm = 0.0d0
		
	do l = 0, lmax, 1
		if (convention == 2) then
				rcilm(1,l+1,1) = ccilm(1,l+1,1)/sqrt(4.0d0*pi)
				rcilm(2, l+1, 1) = 0.0d0
			
				do m=1, l, 1
					if (present(switchcs)) then
						if (switchcs == 1) then
							rcilm(1,l+1, m+1) =  ccilm(1, l+1, m+1) / sqrt(2.0d0*pi) * (-1.0d0)**m
							rcilm(2,l+1, m+1) =  - ccilm(2,l+1,m+1) /sqrt(2.0d0*pi) * (-1.0d0)**m
						else
							rcilm(1,l+1, m+1) =  ccilm(1, l+1, m+1) / sqrt(2.0d0*pi)
							rcilm(2,l+1, m+1) =  - ccilm(2,l+1,m+1) /sqrt(2.0d0*pi)
						endif
					else
						rcilm(1,l+1, m+1) =  ccilm(1, l+1, m+1) / sqrt(2.0d0*pi)
						rcilm(2,l+1, m+1) =  - ccilm(2,l+1,m+1) /sqrt(2.0d0*pi) 
					endif
				enddo			
		else
			rcilm(1,l+1,1) = ccilm(1,l+1,1)
			rcilm(2, l+1, 1) = 0.0d0
			
			do m=1, l, 1
				if (present(switchcs)) then
					if(switchcs == 1) then
						rcilm(1,l+1, m+1) = sqrt(2.0d0) * ccilm(1, l+1, m+1) * (-1.0d0)**m
						rcilm(2,l+1, m+1) = -sqrt(2.0d0) * ccilm(2, l+1, m+1) * (-1.0d0)**m
					else
						rcilm(1,l+1, m+1) = sqrt(2.0d0) * ccilm(1, l+1, m+1)
						rcilm(2,l+1, m+1) = -sqrt(2.0d0) * ccilm(2, l+1, m+1)
					endif
				else
					rcilm(1,l+1, m+1) = sqrt(2.0d0) * ccilm(1, l+1, m+1)
					rcilm(2,l+1, m+1) = -sqrt(2.0d0) * ccilm(2, l+1, m+1)
				endif
			enddo
		
		endif
	enddo
	
end subroutine SHctor


subroutine SHCilmToCindex(degmax, cilm, cindex)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will convert a 3D matrix of spherical harmonics indexed as (i, l+1, m+1)
!	into a 2D matrix that is indexed as (i, index) where index = l(l+1)/2+m+1.
!
!	Calling Parameters:
!		IN
!			cilm	Array of spherical harmonic coefficients with dimensions 
!				(2, lmax+1, lmax+1).
!		OUT
!			cindex	Array of indexed spherical harmonic coefficnets with dimensions
!				(2, (lmax+1)*(lmax+2)/2).
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in)  :: degmax
	real*8, intent(in)   :: cilm(2, degmax+1, degmax+1)
	real*8, intent(out)  :: cindex(2, (degmax+1)*(degmax+2)/2)
	
	integer  :: l, m, index
		
	cindex = 0.0d0
		
	do l=0, degmax
		do m=0, l
			index = (l*(l+1))/2+m+1
			cindex(1, index) = cilm(1,l+1,m+1)
			cindex(2, index) = cilm(2,l+1,m+1)
		enddo
	enddo

end subroutine SHCilmToCindex
		
		
subroutine SHCindexToCilm(degmax, cindex, cilm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will convert a 2D matrix of spherical harmonics indexed as (i, index)
!	into a 2D matrix that is index as (i, l+1, m+1) where index = l(l+1)/2+m+1.
!
!	Calling Parameters:
!		IN
!			cindex	Array of indexed spherical harmonic coefficnets with dimensions
!				(2, (lmax+1)*(lmax+2)/2).
!		OUT
!			cilm	Array of spherical harmonic coefficients with dimensions 
!				(2, lmax+1, lmax+1).
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	degmax
	real*8, intent(out) ::  cilm(2, degmax+1, degmax+1)
	real*8, intent(in)  ::	cindex(2, (degmax+1)*(degmax+2)/2)
	
	integer :: l, m, index
	
	cilm = 0.0d0
		
	do l=0, degmax
		do m=0, l
			index = (l*(l+1))/2+m+1
			cilm(1,l+1,m+1) = cindex(1, index)
			cilm(2,l+1,m+1) = cindex(2, index)	
		enddo
	enddo	

end subroutine SHCindexToCilm

