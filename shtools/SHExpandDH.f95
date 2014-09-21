subroutine SHExpandDH(grid, py_m, py_n, lmax_calc, py_lmax, cilm, lmax, norm, sampling, csphase)
	!                 in     in    in   in    in        in      out   out   in    in         in
	use FFTW3
	use SHTOOLS, only: DHaj
		
	implicit none
	
	integer, intent(in)  ::	py_m, py_n, py_lmax, lmax_calc
	real*8, intent(in)   ::	grid(py_m, py_n)
	real*8, intent(out)  ::	cilm(2, py_lmax+1, py_lmax+1)
	integer, intent(out) ::	lmax
	integer, intent(in), optional :: norm, sampling, csphase
	
	complex*16           ::	cc(py_m+1)
	integer              ::	l, m, n, i, l1, m1, i_eq, i_s, astat(4), lmax_comp, nlong
	integer*8            ::	plan
	real*8               ::	pi, gridl(2*py_m), aj(py_m), fcoef1(2, py_m/2+1), fcoef2(2, py_m/2+1)
	real*8               ::	theta, prod, scalef, rescalem, u, p, pmm, pm1, pm2, z, ffc(1:2,-1:1)
	
	real*8,    save, allocatable ::	sqr(:), ff1(:,:), ff2(:,:)
	integer*1, save, allocatable ::	fsymsign(:,:)
	integer, save                ::	lmax_old=0, norm_old = 0
	
	
	n     = py_m
	lmax  = py_m/2 - 1
	nlong = py_n
	
	lmax_comp = min(lmax, lmax_calc)
		
	pi = acos(-1.0d0)
	
	cilm = 0.0d0
	
	scalef = 1.0d-280
	
	call DHaj(n, aj)
	aj(1:n) = aj(1:n)*sqrt(4.0d0*pi) 	! Driscoll and Heally use unity normalized spherical harmonics
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate recursion constants used in computing the Legendre polynomials
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (lmax_comp /= lmax_old .or. norm /= norm_old) then

		if (allocated(sqr)) deallocate(sqr)
		if (allocated(ff1)) deallocate(ff1)
		if (allocated(ff2)) deallocate(ff2)
		if (allocated(fsymsign)) deallocate(fsymsign)
		
		allocate(sqr(2*lmax_comp+1), stat=astat(1))
		allocate(ff1(lmax_comp+1,lmax_comp+1), stat=astat(2))
		allocate(ff2(lmax_comp+1,lmax_comp+1), stat=astat(3))
		allocate(fsymsign(lmax_comp+1,lmax_comp+1), stat=astat(4))
		
		if (sum(astat(1:4)) /= 0) then
			print*, "Error --- SHExpandDH"
			print*, "Problem allocating arrays SQR, FF1, FF2, or FSYMSIGN", &
				astat(1), astat(2), astat(3), astat(4)
			stop
		endif
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! 	Calculate signs used for symmetry of Legendre functions about equator
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		do l = 0, lmax_comp, 1
			do m = 0, l, 1
				if (mod(l-m,2) == 0) then
					fsymsign(l+1,m+1) = 1
				else
					fsymsign(l+1,m+1) = -1
				endif
			enddo
		enddo
			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		!	Precompute square roots of integers that are used several times.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		do l=1, 2*lmax_comp+1
			sqr(l) = sqrt(dble(l))
		enddo

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! 	Precompute multiplicative factors used in recursion relationships
		! 		P(l,m) = x*f1(l,m)*P(l-1,m) - P(l-2,m)*f2(l,m)
		!		k = l*(l+1)/2 + m + 1
		!	Note that prefactors are not used for the case when m=l as a different 
		!	recursion is used. Furthermore, for m=l-1, Plmbar(l-2,m) is assumed to be zero.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		select case(norm)
		
			case(1,4)
	
				if (lmax_comp /= 0) then
					ff1(2,1) = sqr(3)
					ff2(2,1) = 0.0d0
				endif
				
				do l=2, lmax_comp, 1
					ff1(l+1,1) = sqr(2*l-1) * sqr(2*l+1) / dble(l)
					ff2(l+1,1) = dble(l-1) * sqr(2*l+1) / sqr(2*l-3) / dble(l)
					do m=1, l-2, 1
						ff1(l+1,m+1) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                				ff2(l+1,m+1) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                  				 	/ sqr(2*l-3) / sqr(l+m) / sqr(l-m) 
					enddo
					ff1(l+1,l) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                			ff2(l+1,l) = 0.0d0
				enddo
			
			case(2)
			
				if (lmax_comp /= 0) then
					ff1(2,1) = 1.0d0
					ff2(2,1) = 0.0d0
				endif
				
				do l=2, lmax_comp, 1
					ff1(l+1,1) = dble(2*l-1) /dble(l)
					ff2(l+1,1) = dble(l-1) /dble(l)
					do m=1, l-2, 1
						ff1(l+1,m+1) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                  				ff2(l+1,m+1) = sqr(l-m-1) * sqr(l+m-1) / sqr(l+m) / sqr(l-m)
					enddo
					ff1(l+1,l)= dble(2*l-1) / sqr(l+m) / sqr(l-m)
                  			ff2(l+1,l) = 0.0d0
				enddo
			
			case(3)
		
				do l=1, lmax_comp, 1
					ff1(l+1,1) = dble(2*l-1) /dble(l)
					ff2(l+1,1) = dble(l-1) /dble(l)
					do m=1, l-1, 1
						ff1(l+1,m+1) = dble(2*l-1) / dble(l-m)
                  				ff2(l+1,m+1) = dble(l+m-1) / dble(l-m)
					enddo
				enddo

		end select
	
		lmax_old = lmax_comp
		norm_old = norm
	
	endif	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Create generic plan for gridl
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	call dfftw_plan_dft_r2c_1d(plan, nlong, gridl(1:nlong), cc, FFTW_MEASURE)

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!		
	! 	Integrate over all latitudes. Take into account symmetry of the 
	!	Plms about the equator.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	i_eq = n/2 + 1	! Index correspondong to the equator
	
	do i=2, i_eq - 1, 1
	
		theta = (i-1) * pi / dble(n)
		z = cos(theta)
		u = sqrt( (1.0d0-z) * (1.0d0+z) )
		
		gridl(1:nlong) = grid(i,1:nlong)
		call dfftw_execute(plan)	! take fourier transform
		fcoef1(1,1:n/2) = sqrt(2*pi) * aj(i) * dble(cc(1:n/2)) / dble(nlong)
		fcoef1(2,1:n/2) = -sqrt(2*pi) * aj(i) * dimag(cc(1:n/2)) / dble(nlong)
		
		i_s = 2*i_eq -i
		
		gridl(1:nlong) = grid(i_s,1:nlong)
		call dfftw_execute(plan)	! take fourier transform
		fcoef2(1,1:n/2) = sqrt(2*pi) * aj(i_s) * dble(cc(1:n/2)) / dble(nlong)
		fcoef2(2,1:n/2) = -sqrt(2*pi) * aj(i_s) * dimag(cc(1:n/2)) / dble(nlong)
		
		select case(norm)
			case(1,2,3);	pm2 = 1.0d0
			case(4);	pm2 = 1.0d0 / sqrt(4*pi)
		end select

		cilm(1,1,1) = cilm(1,1,1) + pm2 * (fcoef1(1,1) + fcoef2(1,1) )
		! fsymsign = 1
		
		if (lmax_comp == 0) cycle
				
		pm1 =  ff1(2,1) * z * pm2
		cilm(1,2,1) = cilm(1,2,1) + pm1 * ( fcoef1(1,1) - fcoef2(1,1) )
		! fsymsign = -1
		
		ffc(1,-1) = fcoef1(1,1) - fcoef2(1,1)
		ffc(1, 1) = fcoef1(1,1) + fcoef2(1,1)
		do l=2, lmax_comp, 1
			l1 = l+1
			p = ff1(l1,1) * z * pm1 - ff2(l1,1) * pm2
			pm2 = pm1
			pm1 = p
			cilm(1,l1,1) = cilm(1,l1,1) + p * ffc(1,fsymsign(l1,1))
		enddo
				
		select case(norm)
			case(1,2);	pmm = sqr(2) * scalef
			case(3);	pmm = scalef
			case(4);	pmm = sqr(2) * scalef / sqrt(4*pi)
		end select
				
		rescalem = 1.0d0/scalef
	
		do m = 1, lmax_comp-1, 1
				
			m1 = m+1
			rescalem = rescalem * u
					
			select case(norm)
				case(1,4)
					pmm = csphase * pmm * sqr(2*m+1) / sqr(2*m)
					pm2 = pmm
				case(2)
					pmm = csphase * pmm * sqr(2*m+1) / sqr(2*m)
					pm2 = pmm / sqr(2*m+1)
				case(3)
					pmm = csphase * pmm * (2*m-1)
					pm2 = pmm
			end select
					
			fcoef1(1:2,m1) = fcoef1(1:2,m1) * rescalem
			fcoef2(1:2,m1) = fcoef2(1:2,m1) * rescalem
					
			cilm(1:2,m1,m1) = cilm(1:2,m1,m1) + pm2 * &
					( fcoef1(1:2,m1) + fcoef2(1:2,m1) )
			! fsymsign = 1
					
	   		pm1 = z * ff1(m1+1,m1) * pm2
	   			
	   		cilm(1:2,m1+1,m1) = cilm(1:2,m1+1,m1) + pm1 * &
               				( fcoef1(1:2,m1) - fcoef2(1:2,m1) )
               		! fsymsign = -1
					
			ffc(1:2,-1) = fcoef1(1:2,m1) - fcoef2(1:2,m1)
			ffc(1:2, 1) = fcoef1(1:2,m1) + fcoef2(1:2,m1)
			do l = m+2, lmax_comp, 1
				l1 = l+1
                  		p = z * ff1(l1,m1) * pm1-ff2(l1,m1) * pm2
                  		pm2 = pm1
                  		pm1 = p
				cilm(1:2,l1,m1) = cilm(1:2,l1,m1) + p * ffc(1:2,fsymsign(l1,m1))
			enddo
               				
		enddo
				
		rescalem = rescalem * u
			
            	select case(norm)
            		case(1,4);	pmm = csphase * pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem
            		case(2);	pmm = csphase * pmm / sqr(2*lmax_comp) * rescalem
            		case(3);	pmm = csphase * pmm * (2*lmax_comp-1) * rescalem
        	end select
        			
        	cilm(1:2,lmax_comp+1,lmax_comp+1) = cilm(1:2,lmax_comp+1,lmax_comp+1) + pmm * &
        				( fcoef1(1:2,lmax_comp+1) + fcoef2(1:2,lmax_comp+1) )	
        				! fsymsign = 1
	enddo
	
	! Finally, do equator
	
	i = i_eq
	
	z = 0.0d0
	u = 1.0d0
	
	gridl(1:nlong) = grid(i,1:nlong)
	call dfftw_execute(plan)	! take fourier transform
	fcoef1(1,1:n/2) = sqrt(2*pi) * aj(i) * dble(cc(1:n/2)) / dble(nlong)
	fcoef1(2,1:n/2) = -sqrt(2*pi) * aj(i) * dimag(cc(1:n/2)) / dble(nlong)


	select case(norm)
		case(1,2,3);	pm2 = 1.0d0
		case(4);	pm2 = 1.0d0 / sqrt(4*pi)
	end select

	cilm(1,1,1) = cilm(1,1,1) + pm2 * fcoef1(1,1)
	
	if (lmax_comp /= 0) then
				
		do l=2, lmax_comp, 2
			l1 = l+1
			p = - ff2(l1,1) * pm2
			pm2 = p
			cilm(1,l1,1) = cilm(1,l1,1) + p * fcoef1(1,1)
		enddo
				
		select case(norm)
			case(1,2);	pmm = sqr(2) * scalef
			case(3);	pmm = scalef
			case(4);	pmm = sqr(2) * scalef / sqrt(4*pi)
		end select
				
		rescalem = 1.0d0/scalef
	
		do m = 1, lmax_comp-1, 1
				
			m1 = m+1
					
			select case(norm)
				case(1,4)
					pmm = csphase * pmm * sqr(2*m+1) / sqr(2*m)
					pm2 = pmm
				case(2)
					pmm = csphase * pmm * sqr(2*m+1) / sqr(2*m)
					pm2 = pmm / sqr(2*m+1)
				case(3)
					pmm = csphase * pmm * (2*m-1)
					pm2 = pmm
			end select

			fcoef1(1:2,m1) = fcoef1(1:2,m1) * rescalem
					
			cilm(1:2,m1,m1) = cilm(1:2,m1,m1) + pm2 * fcoef1(1:2,m1)
								
			do l = m+2, lmax_comp, 2
				l1 = l+1
                	  	p = - ff2(l1,m1) * pm2
                  		pm2 = p
				cilm(1:2,l1,m1) = cilm(1:2,l1,m1) + p * fcoef1(1:2,m1) 
			enddo

		enddo		
       
       		select case(norm)
            		case(1,4);	pmm = csphase * pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem
            		case(2);	pmm = csphase * pmm / sqr(2*lmax_comp) * rescalem
            		case(3);	pmm = csphase * pmm * (2*lmax_comp-1) * rescalem
       		end select
       			
       		cilm(1:2,lmax_comp+1,lmax_comp+1) = cilm(1:2,lmax_comp+1,lmax_comp+1) + pmm * fcoef1(1:2,lmax_comp+1) 
       		
       	endif

	call dfftw_destroy_plan(plan) 
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Divide by integral of Ylm*Ylm 
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	select case(norm)
	
		case(1)
	
			do l=0, lmax_comp, 1
				cilm(1:2,l+1, 1:l+1) = cilm(1:2,l+1, 1:l+1) / (4*pi)
			enddo
		
		case(2)
	
			do l=0, lmax_comp, 1
				cilm(1:2,l+1, 1:l+1) = cilm(1:2,l+1, 1:l+1) * (2*l+1) / (4*pi)
			enddo
	
		case(3)
		 
			do l=0, lmax_comp, 1
				prod = 4 * pi / dble(2*l+1)
				cilm(1,l+1,1) = cilm(1,l+1,1) / prod
				prod = prod / 2.0d0
				do m=1, l-1, 1
					prod = prod * (l+m) * (l-m+1)
					cilm(1:2,l+1,m+1) = cilm(1:2,l+1,m+1) / prod
				enddo
				!do m=l case
				if (l /= 0) cilm(1:2,l+1,l+1) = cilm(1:2,l+1, l+1)/(prod*2*l)
			enddo
			
	end select
		
end subroutine SHExpandDH
