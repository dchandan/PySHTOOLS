subroutine PLegendreA(p, lmax, z, csphase)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function evalutates all of the unnormalized associated legendre
!	polynomials up to degree lmax. 
!
!	Calling Parameters:
!		IN
!			lmax:		Maximum spherical harmonic degree to compute.
!			z:		[-1, 1], cos(colatitude) or sin(latitude).
!			csphase:	1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!		OUT
!			p:	A vector of all associated Legendgre polynomials evaluated at 
!				z up to lmax. The length must by greater or equal to (lmax+1)*(lmax+2)/2.
!	Notes:
!
!	1.	The integral of plm**2 over (-1,1) is 2 * (l+m)! / (2l+1) / (l-m)!.
!	2.	The index of the array p corresponds to l*(l+1)/2 + m + 1. 
!	3.	The default is to exlude the Condon-Shortley phase of (-1)^m.
!
!
!	This code is based on the stable recursion relationships in numerical recipes.
!
!	August 14, 2012. 	Converted PHASE from REAL*8 to INTEGER*1.
!
!	Written by Mark Wieczorek 2003
!
!	Copyright (c) 2005-2012, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	implicit none
	integer, intent(in) ::	lmax, csphase
	real*8, intent(out) ::	p((lmax+1)*(lmax+2)/2)
   	real*8, intent(in)  ::	z
	
   	real*8              ::	pm2, pm1, pmm, sinsq, sinsqr, fact, plm
  	integer             ::	k, kstart, m, l
	
	if (size(p) < (lmax+1)*(lmax+2)/2) then 
		print*, "Error --- PLegendreA"
 		print*, "P must be dimensioned as (LMAX+1)*(LMAX+2)/2 where LMAX is ", lmax
 		print*, "Input array is dimensioned ", size(p)
 		stop
	endif
     	     	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!	
	!	Calculate P(l,0)
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	sinsq = (1.0d0-z)*(1.0d0+z)
	sinsqr = sqrt(sinsq)

      	pm2  = 1.0d0	
      	p(1) = 1.0d0
      	
      	if (lmax == 0) return
      	
      	pm1  = z	
      	p(2) = pm1
      		
	k = 2

      	do l = 2, lmax, 1
         	k = k+l
         	plm = ( z*(2*l-1)*pm1-(l-1)*pm2 ) / dble(l)
         	p(k) = plm
         	pm2  = pm1
         	pm1  = plm
      	enddo

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate P(m,m), P(m+1,m), and P(l,m)
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
      	pmm  = 1.0d0
      	fact = -1.0d0
      	kstart = 1

      	do m = 1, lmax - 1, 1

		! Calculate P(m,m)
        	kstart = kstart+m+1
         	fact = fact+2.0d0
        	pmm = csphase * pmm * sinsqr * fact
        	p(kstart) = pmm
        	pm2 = pmm

		! Calculate P(m+1,m)
		k = kstart+m+1
	   	pm1 = z*pmm*(2*m+1)
	    	p(k) = pm1

		! Calculate P(l,m)
               	do l = m+2, lmax, 1
               		k = k+l
                  	plm  = ( z*(2*l-1)*pm1-(l+m-1)*pm2 ) / dble(l-m)
                  	p(k) = plm
                  	pm2  = pm1
                  	pm1  = plm
               	enddo

      	enddo
      	
      	! P(lmax, lmax)
      	kstart    = kstart+m+1
        fact      = fact+2.0d0
        pmm       = csphase * pmm * sinsqr * fact
        p(kstart) = pmm
	
end subroutine PLegendreA

