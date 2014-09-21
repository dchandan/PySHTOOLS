subroutine MakeGrid2D(grid, cilm, lmax, interval, py_m, py_n, north, south, east, west, nlat, nlong, norm, csphase, f, a, dealloc)
	use SHTOOLS, only: PlmBar, PlBar, PlmSchmidt, PlSchmidt, PLegendreA, PLegendre, &
				PlmON, PlON, CSPHASE_DEFAULT

	implicit none

	real*8, intent(in)   :: cilm(2,lmax+1,lmax+1), interval
	real*8, intent(out)  :: grid(py_m, py_n)
	integer, intent(in)  :: lmax, py_m, py_n
	real*8, intent(in)   :: north, south, east, west
	integer, intent(out) :: nlat, nlong
	integer, intent(in), optional :: norm, csphase, dealloc
	real*8, intent(in), optional :: f, a
	
	integer              :: l, m, j, k, index, l1, m1, lmax_comp
	integer              :: temp, astat(4)
	real*8               :: pi, latmax, latmin, longmin, longmax
	real*8               :: lat, longitude, x, intervalrad, r_ref
	real*8, allocatable  ::	pl(:), cosm(:, :), sinm(:, :), cilm2(:,:, :)
	

	latmax = north
	latmin = south
	longmin = west
	longmax = east
		
	nlat = (latmax - latmin) / interval + 1
	nlong = (longmax - longmin) / interval + 1
	

    
	if ((present(f) .and. .not. present(a)) .or. (present(a) .and. .not. present(f)) ) then
		print*, "Error --- MakeGrid2D"
		print*, "Both F and A must be present"
		print*, "F ", present(f)
		print*, "A ", present(a)
		stop
	endif
	
	
	allocate(pl((lmax+1)*(lmax+2)/2), stat = astat(1))
	allocate(cosm(lmax+1, int(360./interval + 1)), stat = astat(2))
	allocate(sinm(lmax+1, int(360./interval + 1)), stat = astat(3))
	allocate(cilm2(2,lmax+1,lmax+1), stat = astat(4))
	if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /=0 .or. astat(4) /= 0) then
		print*, "Error --- MakeGrid2D"
		print*, "Problem allocating arrays PL, COSM, SINM, and CILM2", astat(1), astat(2), astat(3), astat(4)
		stop
	endif
	
	pi = acos(-1.0d0)
	grid = 0.0d0
	
	lmax_comp = min(lmax, size(cilm(1,1,:))-1)

	intervalrad = interval*pi/180.0d0
	
	cilm2(1:2,1:lmax+1,1:lmax+1) = cilm(1:2,1:lmax+1,1:lmax+1)
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Precomputing sines and cosines leads to an increase in speed by a factor
	!	of almost 4 with no optimization, and by a factor of about 15 with normal 
	!	optimizations.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do k=1, nlong
		longitude = longmin*pi/180.0d0 + dble(k-1)*intervalrad
		sinm(1,k) = 0.0d0
		cosm(1,k) = 1.0d0
		if (lmax > 0) then
			sinm(2,k) = sin(longitude)
			cosm(2,k) = cos(longitude)
		endif
		do m=2, lmax, 1
			sinm(m+1,k) = 2 * sinm(m,k)*cosm(2,k) - sinm(m-1,k)
			cosm(m+1,k) = 2 * cosm(m,k)*cosm(2,k) - cosm(m-1,k)
		enddo		
	enddo
	
	do j=1, nlat
		lat = latmax - (j-1)*interval
		x = sin(lat*pi/180.0d0)
		
		if (lat == 90.0d0 .or. lat == -90.0d0) then
			select case(norm)
				case(1); call PlBar(pl, lmax_comp, x)
				case(2); call PlSchmidt(pl, lmax_comp, x)
				case(3); call PLegendre(pl, lmax_comp, x)
				case(4); call PlON(pl, lmax_comp, x)
			end select
			
			do l = lmax_comp, 0, -1
				l1 = l+1
				grid(j, 1) = grid(j, 1) + cilm2(1, l1, 1) * pl(l1)
			enddo
			grid(j, 2:nlong) = grid(j,1)
			if (present(f)) grid(j,1:nlong) = grid(j,1:nlong) - a*(1.0d0-f)
		else
			select case(norm)
				case(1); call PlmBar(pl, lmax_comp, x, csphase = csphase)
				case(2); call PlmSchmidt(pl, lmax_comp, x, csphase = csphase)
				case(3); call PLegendreA(pl, lmax_comp, x, csphase = csphase)
				case(4); call PlmON(pl, lmax_comp, x, csphase = csphase)
			end select
			
			do k = 1, nlong
				! do m = 0 term first
				m1 = 1
				do l=0, lmax_comp, 1
					l1 = l+1
					index = (l+1)*l/2 + 1
					grid(j,k) = grid(j,k) + cilm2(1,l1,1)*cosm(1,k) * pl(index)
				enddo
			
				do m=1, lmax_comp, 1
					m1 = m+1
					do l=m, lmax_comp, 1
						l1 = l+1
						index = (l+1)*l/2 + m + 1
						grid(j,k) = grid(j,k) + ( cilm2(1,l1,m1)*cosm(m1,k) + &
							cilm2(2,l1,m1)*sinm(m1,k) ) * pl(index)
					enddo
				enddo
			enddo
		
			if (present(f)) then
				r_ref = a**2 * (1.0d0 + tan(lat*pi/180.0d0)**2) / &
					(1.0d0  + tan(lat*pi/180.0d0)**2 / (1.0d0 - f)**2 )
				r_ref = sqrt(r_ref)
				grid(j,1:nlong) = grid(j,1:nlong) - r_ref
			endif			
		endif
	enddo
	
	! deallocate memory
	if (present(dealloc)) then
		if (dealloc == 1) then
			select case(norm)
				case(1); call PlmBar(pl, -1, x, csphase = csphase)
				case(2); call PlmSchmidt(pl, -1, x, csphase = csphase)
				case(4); call PlmON(pl, -1, x, csphase = csphase)
			end select
		endif
	endif
	
	deallocate(pl)
	deallocate(cosm)
	deallocate(sinm)
	deallocate(cilm2)
end subroutine MakeGrid2D
