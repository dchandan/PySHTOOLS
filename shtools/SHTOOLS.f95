module SHTOOLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This module contains an interface block defining all the routines
!	used in the archive SHTOOLS. These are necessary in order to use
!	implicitly shaped arrays with most subroutines.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	integer, parameter ::	CSPHASE_DEFAULT = 1	! The default is to EXCLUDE the 
							! CONDON-SHORTLEY phase of (-1)^m
							! in front of the Legendre functions.
							! To use this phase function, set
							! CSPHASE_DEFAULT = -1

	interface
	
		subroutine PlmBar(p, lmax, z, csphase, cnorm)
			integer, intent(in) ::	lmax, csphase
			real*8, intent(out) ::	p((lmax+1)*(lmax+2)/2)
   			real*8, intent(in)  ::	z
   			integer, intent(in), optional :: cnorm
   		end subroutine PlmBar
       		
		subroutine PlmBar_d1(p, dp, lmax, z, csphase, cnorm)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p(:), dp(:)
   			real*8, intent(in) ::	z
   			integer, intent(in), optional :: csphase, cnorm
		end subroutine PlmBar_d1
       		
   		subroutine PlBar(p, lmax, z)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p((lmax+1)*(lmax+2)/2)
   			real*8, intent(in)  ::	z
   		end subroutine PlBar
       		
   		subroutine PlBar_d1(p, dp, lmax, z)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p(:), dp(:)
   			real*8, intent(in) ::	z
   		end subroutine PlBar_d1
       		
   		subroutine PlmSchmidt(p, lmax, z, csphase, cnorm)
			integer, intent(in) ::	lmax, csphase
			real*8, intent(out) ::	p((lmax+1)*(lmax+2)/2)
  		 	real*8, intent(in)  ::	z
  		 	integer, intent(in), optional :: cnorm
  		end subroutine PlmSchmidt
      		
  		subroutine PlSchmidt(p, lmax, z)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p((lmax+1)*(lmax+2)/2)
  		 	real*8, intent(in)  ::	z
  		end subroutine PlSchmidt
      		
  		subroutine PlmSchmidt_d1(p, dp, lmax, z, csphase, cnorm)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p(:), dp(:)
  		 	real*8, intent(in) ::	z
  		 	integer, intent(in), optional :: csphase, cnorm
  		end subroutine PlmSchmidt_d1
      		
  		subroutine PlSchmidt_d1(p, dp, lmax, z)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p(:), dp(:)
  		 	real*8, intent(in) ::	z
  		end subroutine PlSchmidt_d1
      		
  		subroutine PLegendre(p, lmax, z)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p((lmax+1)*(lmax+2)/2)
   			real*8, intent(in)  ::	z
   		end subroutine PLegendre
       		
   		subroutine PLegendreA(p, lmax, z, csphase)
			integer, intent(in) ::	lmax, csphase
			real*8, intent(out) ::	p((lmax+1)*(lmax+2)/2)
   			real*8, intent(in)  ::	z
   		end subroutine PLegendreA
       		
  		subroutine PLegendre_d1(p, dp, lmax, z)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p(:), dp(:)
   			real*8, intent(in) ::	z
   		end subroutine PLegendre_d1
       		
   		subroutine PLegendreA_d1(p, dp, lmax, z, csphase)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p(:), dp(:)
   			real*8, intent(in) ::	z
			integer, intent(in), optional :: csphase
   		end subroutine PLegendreA_d1
       	
   		subroutine CilmPlus(cilm, gridin, lmax, nmax, mass, d, rho, gridtype, w, zero, plx, n, dref)
			real*8, intent(in) :: 	gridin(:,:), mass, rho
			real*8, intent(in), optional :: w(:), zero(:), plx(:,:), dref
			real*8, intent(out) :: 	cilm(:,:,:), d
			integer, intent(in) :: 	lmax, nmax, gridtype
			integer, intent(in), optional :: n
		end subroutine CilmPlus
		
	 	subroutine CilmPlusRhoH(cilm, gridin, lmax, nmax, mass, d, rho, gridtype, w, zero, plx, n, dref)
			real*8, intent(in) :: 	gridin(:,:), mass, rho(:,:)
			real*8, intent(in), optional :: w(:), zero(:), plx(:,:), dref
			real*8, intent(out) :: 	cilm(:,:,:), d
			integer, intent(in) :: 	lmax, nmax, gridtype
			integer, intent(in), optional :: n
		end subroutine CilmPlusRhoH
		
		
		subroutine MakeGrid2d(grid, cilm, lmax, interval, py_m, py_n, north, south, east, west, nlat, nlong, norm, csphase, f, a, dealloc)
			integer, intent(in)  :: lmax, py_m, py_n
			real*8, intent(in)   :: cilm(2,lmax+1,lmax+1), interval
			real*8, intent(out)  :: grid(py_m, py_n)
			real*8, intent(in)   :: north, south, east, west
			integer, intent(out) :: nlat, nlong
			integer, intent(in), optional :: norm, csphase, dealloc
			real*8, intent(in), optional :: f, a
		end subroutine MakeGrid2D

		subroutine GLQGridCoord(latglq, longlq, lmax)
			integer, intent(in)  ::	lmax
			real*8, intent(out)  ::	latglq(lmax+1), longlq(2*lmax+1)
		end subroutine GLQGridCoord
		
		subroutine MakeGridGLQ(gridglq, cilm, lmax, plx, zero, norm, csphase, lmax_calc)
			real*8, intent(in) ::	cilm(:,:,:)
			real*8, intent(in), optional :: plx(:,:), zero(:)
			real*8, intent(out) ::	gridglq(:,:)
			integer, intent(in) ::	lmax
			integer, intent(in), optional :: norm, csphase, lmax_calc
		end subroutine MakeGridGLQ
		
		subroutine SHExpandGLQ(cilm, lmax, gridglq, w, plx, zero, norm, csphase, lmax_calc)	
			real*8, intent(in) ::	w(:), gridglq(:,:)
			real*8, intent(in), optional ::	plx(:,:), zero(:)
			real*8, intent(out) :: 	cilm(:,:,:)
			integer, intent(in) ::	lmax
			integer, intent(in), optional :: norm, csphase, lmax_calc
		end subroutine SHExpandGLQ
		
		subroutine PreCompute(lmax, zero, w, plx, wisdom_file, norm, csphase, cnorm)
			integer, intent(in)                ::  lmax
			real*8, intent(out)                ::  zero(lmax+1), w(lmax+1)
			real*8, intent(out), optional      ::  plx(:,:)
			integer, intent(in), optional      ::  norm, csphase, cnorm
			character(*), intent(in), optional ::  wisdom_file
		end subroutine PreCompute

		subroutine PreGLQ(x1, x2, n, zero, w)
			REAL*8, INTENT(IN)  :: 	X1, X2
			REAL*8, INTENT(OUT) ::	ZERO(N), W(N)
			INTEGER, INTENT(IN) ::	N
		end subroutine PreGLQ
		
		integer function NGLQ(degree)
			integer, intent(in) ::	degree
		end function NGLQ

		integer function NGLQSH(degree)
			integer, intent(in) ::	degree
		end function NGLQSH

		integer function NGLQSHN(degree, n)
			integer, intent(in) ::	degree, n
		end function NGLQSHN

		subroutine SHRead(filename, cilm, lmax, skip, header, error)
			character(*), intent(in) ::		filename
			integer, intent(out) ::			lmax
			real*8, intent(out) ::			cilm(:,:,:)
			real*8, intent(out), optional ::	header(:), error(:,:,:)
			integer, intent(in), optional ::	skip
		end subroutine SHRead
		
		subroutine MakeMagGridDH(cilm, lmax, r0, a, f, rad_grid, theta_grid, phi_grid, total_grid, n, sampling, lmax_calc, pot_grid)
			real*8, intent(in) :: 	cilm(:,:,:), r0, a, f
			real*8, intent(out) ::	rad_grid(:,:), theta_grid(:,:), phi_grid(:,:), total_grid(:,:)
			real*8, intent(out), optional :: pot_grid(:,:)
			integer, intent(in) :: 	lmax
			integer, intent(out) ::	n
			integer, intent(in), optional :: sampling, lmax_calc
		end subroutine MakeMagGridDH
		
		real*8 function SHPowerL(c, l)
			real*8, intent(in) :: c(:,:,:)
			integer, intent(in) :: l
		end function SHPowerL
		
		real*8 function SHPowerDensityL(c, l)
			real*8, intent(in) :: 	c(:,:,:)
			integer, intent(in) :: 	l
		end function SHPowerDensityL
		
		real*8 function SHCrossPowerL(c1, c2, l)
			real*8, intent(in) :: 	c1(:,:,:), c2(:,:,:)
			integer, intent(in) :: 	l
		end function SHCrossPowerL
		
		real*8 function SHCrossPowerDensityL(c1, c2, l)
			real*8, intent(in) :: 	c1(:,:,:), c2(:,:,:)
			integer, intent(in) :: 	l
		end function SHCrossPowerDensityL
		
		subroutine SHPowerSpectrum(c, lmax, spectra)
			real*8, intent(in) :: 	c(:,:,:)
			integer, intent(in) :: 	lmax
			real*8, intent(out) ::	spectra(:)
		end subroutine SHPowerSpectrum
		
		subroutine SHPowerSpectrumDensity(c, lmax, spectra)
			real*8, intent(in) :: 	c(:,:,:)
			integer, intent(in) :: 	lmax
			real*8, intent(out) ::	spectra(:)
		end subroutine SHPowerSpectrumDensity
		
		subroutine SHCrossPowerSpectrum(c1, c2, lmax, cspectra)
			real*8, intent(in) :: 	c1(:,:,:), c2(:,:,:)
			integer, intent(in) :: 	lmax
			real*8, intent(out) ::	cspectra(:)
		end subroutine SHCrossPowerSpectrum
		
		subroutine SHCrossPowerSpectrumDensity(c1, c2, lmax, cspectra)
			real*8, intent(in) :: 	c1(:,:,:), c2(:,:,:)
			integer, intent(in) :: 	lmax
			real*8, intent(out) ::	cspectra(:)
		end subroutine SHCrossPowerSpectrumDensity
		
		subroutine djpi2(dj, lmax)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	dj(:,:,:)
		end subroutine djpi2
		
		subroutine SHrtoc(rcilm, l1, ccilm, degmax, convention, switchcs)
			integer, intent(in)           ::    l1, degmax
			real*8, intent(in)            :: 	rcilm(2,l1+1,l1+1)
			real*8, intent(out)           ::	ccilm(2,degmax+1,degmax+1) 
			integer, intent(in), optional ::	convention, switchcs
		end subroutine SHrtoc
		
		subroutine SHctor(ccilm, l1, rcilm, degmax, convention, switchcs)
			integer, intent(in) :: l1, degmax
			real*8, intent(in) :: ccilm(2,l1+1,l1+1)
			real*8, intent(out) :: rcilm(2,degmax+1,degmax+1)
			integer, intent(in), optional :: convention, switchcs
		end subroutine SHctor
		
		subroutine SHCilmToCindex(degmax, cilm, cindex)
			integer, intent(in)  :: degmax
			real*8, intent(in)   :: cilm(2, degmax+1, degmax+1)
			real*8, intent(out)  :: cindex(2, (degmax+1)*(degmax+2)/2)
		end subroutine SHCilmToCindex
		
		subroutine SHCindexToCilm(degmax, cindex, cilm)
			integer, intent(in) ::	degmax
			real*8, intent(out) ::  cilm(2, degmax+1, degmax+1)
			real*8, intent(in)  ::	cindex(2, (degmax+1)*(degmax+2)/2)
		end subroutine SHCindexToCilm
		
		subroutine SHRotateCoef(x, cof, rcof, dj, lmax)
			real*8, intent(in) :: 	cof(:,:), dj(:,:,:), x(3)
			real*8, intent(out) ::	rcof(:,:)
			integer, intent(in) :: 	lmax
		end subroutine SHRotateCoef
		
		subroutine SHRotateRealCoef(cilmrot, cilm, lmax, x, dj)
			real*8, intent(in) ::	cilm(:,:,:), x(:), dj(:,:,:)
			real*8, intent(out) ::	cilmrot(:,:,:)
			integer, intent(in) ::	lmax
		end subroutine SHRotateRealCoef
		
		subroutine DHaj(n, aj)
			integer, intent(in) ::	n
			real*8, intent(out) ::	aj(:)
		end subroutine DHaj
		
		subroutine SHExpandDH(grid, py_m, py_n, lmax_calc, py_lmax, cilm, lmax, norm, sampling, csphase)
			integer, intent(in)  ::	py_m, py_n, py_lmax, lmax_calc
			real*8, intent(in)   ::	grid(py_m, py_n)
			real*8, intent(out)  ::	cilm(2, py_lmax+1, py_lmax+1)
			integer, intent(out) ::	lmax
			integer, intent(in), optional :: norm, sampling, csphase
		end subroutine SHExpandDH
		
		subroutine MakeGridDH(griddh, py_m, py_n, cilm, lmax, csphase, norm, sampling)
			integer, intent(in)           ::  lmax, py_m, py_n, csphase
			real*8, intent(in)            ::  cilm(2, lmax+1, lmax+1)
			real*8, intent(out)           ::  griddh(py_m, py_n)
			integer, intent(in), optional ::  norm, sampling
		end subroutine MakeGridDH
		
		real*8 function MakeGridPoint(cilm, lmax, lat, longitude, csphase, norm, dealloc)
			real*8, intent(in)            :: cilm(2,lmax+1,lmax+1), lat, longitude
			integer, intent(in)           :: lmax, csphase
			integer, intent(in), optional :: norm, dealloc
		end function MakeGridPoint
		
		real*8 function Wl(l, half, r, d)
			integer, intent(in) ::	l, half
			real*8, intent(in) ::	r, d
		end function Wl
		
		real*8 function WlCurv(l, half, r, d)
			integer, intent(in) ::	l, half
			real*8, intent(in) ::	r, d
		end function WlCurv
		
		subroutine SHExpandLSQ(cilm, chi2, d, lat, lon, nmax, lmax, csphase, norm)
			real*8, intent(in)            ::  d(nmax), lat(nmax), lon(nmax)
			real*8, intent(out)           ::  cilm(2,lmax+1,lmax+1), chi2
			integer, intent(in)           ::  nmax, lmax, csphase
			integer, intent(in), optional ::  norm
		end subroutine SHExpandLSQ

		subroutine SHMultiply(shout, sh1, lmax1, sh2, lmax2, precomp, norm, csphase)
			real*8, intent(out) ::	shout(:,:,:)
			real*8, intent(in) ::	sh1(:,:,:), sh2(:,:,:)
			integer, intent(in) ::	lmax1, lmax2
			integer, intent(in), optional ::	precomp, norm, csphase
		end subroutine SHMultiply
		
		subroutine ComputeD0(D0, lmax, theta0)
			real*8, intent(out) ::	D0(:,:)
			real*8, intent(in) ::	theta0
			integer, intent(in) :: 	lmax
		end subroutine ComputeD0
		
		subroutine ComputeDm(dllm, lmax, m, theta0)
			real*8, intent(out) ::	dllm(:,:)
			real*8, intent(in) ::	theta0
			integer, intent(in) :: 	lmax, m
		end subroutine ComputeDm
		
		subroutine SphericalCapCoef(coef, theta, lmax)
			real*8, intent(out) ::	coef(:)
			real*8, intent(in) ::	theta
			integer, intent(in), optional ::	lmax
		end subroutine SphericalCapCoef
		
		real*8 function SHDeltaL(coef, lmax)
			real*8, intent(in) ::	coef(:)
			integer, intent(in) :: 	lmax
		end function SHDeltaL
		
		real*8 function SHDeltaX(coef, lmax, m)
			real*8, intent(in) ::	coef(:)
			integer, intent(in) :: 	lmax
			integer, intent(in), optional :: m
		end function SHDeltaX

		subroutine EigValVecSym(ain, n, eig, evec, ul, K)
			real*8, intent(in) ::	ain(:,:)
			integer, intent(in) ::	n
			real*8, intent(out) ::	eig(:), evec(:,:)
			character, intent(in), optional ::	ul
			integer, intent(in), optional ::	K
		end subroutine EigValVecSym
		
		subroutine SHReturnTapersM(theta0, lmax, m, tapers, eigenvalues, shannon)
			real*8, intent(in) ::	theta0
			integer, intent(in) ::	lmax, m
			real*8, intent(out) ::	tapers(:,:), eigenvalues(:)
			real*8, intent(out), optional :: shannon
		end subroutine SHReturnTapersM
		
		subroutine EigValSym(ain, n, eval, ul)
			real*8, intent(in) ::	ain(:,:)
			integer, intent(in) ::	n
			real*8, intent(out) ::	eval(:)
			character, intent(in), optional :: ul
		end subroutine EigValSym
		
		integer function SHFindLWin(theta0, m, alpha, taper_number)
			real*8, intent(in) :: 	theta0, alpha
			integer, intent(in) :: 	m
			integer, intent(in), optional ::	taper_number
		end function SHFindLWin
		
		subroutine SHAdmitCorr(G, T, lmax, admit, corr, admit_error)
			real*8, intent(in) ::	G(:,:,:), T(:,:,:)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	admit(:), corr(:)
			real*8, intent(out), optional ::	admit_error(:)
		end subroutine SHAdmitCorr
										
		integer function PlmIndex(l,m)
			integer, intent(in)	:: l, m
		end function PlmIndex
				
		real*8 function RandomN(idum)
			integer, parameter ::	K4B=selected_int_kind(9)
			integer(K4B), intent(inout) ::	idum
		end function RandomN
		
		real*8 function RandomGaussian(idum)
			integer, parameter ::	K4B=selected_int_kind(9)
			integer(K4B), intent(inout) :: 	idum
		end function RandomGaussian
 	 	
 	 	subroutine Wigner3j(w3j, jmin, jmax, j2, j3, m1, m2, m3)
			integer, intent(in) ::	j2, j3, m1, m2, m3
			integer, intent(out) ::	jmin, jmax
			real*8, intent(out) ::	w3j(:)
		end subroutine Wigner3j
 	 	
 	 	subroutine SHBias(Shh, lwin, incspectra, ldata, outcspectra, save_cg)
			real*8, intent(in) ::	Shh(:), incspectra(:)
			real*8, intent(out) ::	outcspectra(:)
			integer, intent(in) ::	lwin, ldata
			integer, intent(in), optional :: save_cg
		end subroutine SHBias
		
		subroutine SHBiasK(tapers, lwin, numk, incspectra, ldata, outcspectra, taper_wt, save_cg)
			real*8, intent(in) ::	tapers(:,:), incspectra(:)
			real*8, intent(out) ::	outcspectra(:)
			integer, intent(in) ::	lwin, ldata, numk
			real*8, intent(in), optional :: taper_wt(:)
			integer, intent(in), optional :: save_cg
		end subroutine SHBiasK
		
		real*8 function SHSjkPG0(incspectra, j, k, l, m, evec, lwin)
			real*8, intent(in) ::	incspectra(:), evec(:,:)
			integer, intent(in) ::	lwin, l, m, j, k
		end function SHSjkPG0
		
		Subroutine SHMTVarOpt0(l, tapers, lwin, kmax, Sff, var_opt, var_unit, weight_opt, unweighted_covar, nocross)
			real*8, intent(in) ::	tapers(:,:), Sff(:)
			real*8, intent(out) ::	var_opt(:), var_unit(:)
			integer, intent(in) ::	l, lwin, kmax
			real*8, intent(out), optional ::	weight_opt(:,:), unweighted_covar(:,:)
			integer, intent(in), optional ::	nocross
		end subroutine SHMTVarOpt0
		
		subroutine SHMultiTaperSE(mtse, sd, sh, lmax, tapers, taper_order, lmaxt, K, alpha, &
			lat, lon, taper_wt, norm, csphase)
			real*8, intent(out) ::	mtse(:), sd(:)
			real*8, intent(in) ::	sh(:,:,:), tapers(:,:)
			integer, intent(in) ::	lmax, lmaxt, K, taper_order(:)
			real*8, intent(in), optional ::	alpha(:), lat, lon, taper_wt(:)
			integer, intent(in), optional :: csphase, norm
		end subroutine SHMultiTaperSE
		
		subroutine SHMultiTaperCSE(mtse, sd, sh1, lmax1, sh2, lmax2, tapers, taper_order, lmaxt, K, &
			alpha, lat, lon, taper_wt, norm, csphase)
			real*8, intent(out) ::	mtse(:), sd(:)
			real*8, intent(in) ::	sh1(:,:,:), sh2(:,:,:), tapers(:,:)
			integer, intent(in) ::	lmax1, lmax2, lmaxt, K, taper_order(:)
			real*8, intent(in), optional ::	alpha(:), lat, lon, taper_wt(:)
			integer, intent(in), optional ::	csphase, norm
		end subroutine SHMultiTaperCSE
		
		subroutine SHReadJPL(filename, cilm, lmax, error, gm, formatstring)
			character(*), intent(in) ::	filename
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	cilm(:,:,:)
			real*8, intent(out), optional :: error(:,:,:), gm(2)
			character, intent(in), optional :: formatstring*6
		end subroutine SHReadJPL
		
		subroutine SHRead2(filename, cilm, lmax, gm, r0_pot, error, dot, doystart, doyend, epoch)
			character(*), intent(in) ::	filename
			integer, intent(out) ::	lmax
			real*8, intent(out) ::	cilm(:,:,:), gm, r0_pot
			real*8, intent(out), optional :: error(:,:,:), dot(:,:,:), doystart, doyend, epoch
		end subroutine SHRead2
				
		subroutine MakeGeoidGrid(geoid, cilm, lmax, r0pot, GM, PotRef, omega, r, gridtype, &
			order, nlat, nlong, interval, lmax_calc, a, f)
			real*8, intent(out) ::	geoid(:,:)
			real*8, intent(in) ::	cilm(:,:,:), r0pot, GM, r, PotRef, omega
			integer, intent(in) :: lmax, order, gridtype
			integer, intent(in), optional :: lmax_calc
			integer, intent(out) :: nlat, nlong
			real*8, intent(in), optional ::	interval, a, f	
		end subroutine MakeGeoidGrid
			
		subroutine PlmON(p, lmax, z, csphase, cnorm)
			integer, intent(in) ::	lmax, csphase
			real*8, intent(out) ::	p((lmax+1)*(lmax+2)/2)
   			real*8, intent(in)  ::	z
   			integer, intent(in), optional :: cnorm
   		end subroutine PlmON
       		
       	subroutine PlON(p, lmax, z)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p((lmax+1)*(lmax+2)/2)
       		real*8, intent(in)  ::	z
       	end subroutine PlON

		subroutine PlmON_d1(p, dp, lmax, z, csphase, cnorm)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p(:), dp(:)
   			real*8, intent(in) ::	z
   			integer, intent(in), optional :: csphase, cnorm
   		end subroutine PlmON_d1

		subroutine PlON_d1(p, dp, lmax, z)
			integer, intent(in) ::	lmax
			real*8, intent(out) ::	p(:), dp(:)
   			real*8, intent(in) ::	z
   		end subroutine PlON_d1
       		
		
		real*8 function SHConfidence(l_conf, r)
			real*8, intent(in) :: r
			integer, intent(in) :: l_conf
		end function SHConfidence
				
		subroutine SHMagPowerSpectrum(c, a, r, lmax, spectra)
			real*8, intent(in) :: c(:,:,:)
			real*8, intent(in) :: a, r
			integer, intent(in) :: lmax
			real*8, intent(out) ::	spectra(:)
		end subroutine SHMagPowerSpectrum
		
		subroutine SHExpandDHC(grid, py_m, py_n, lmax_calc, py_lmax, cilm, lmax, norm, sampling, csphase)
			integer, intent(in)     ::  py_m, py_n, lmax_calc, py_lmax
			complex*16, intent(in)  ::	grid(py_m, py_n)
			complex*16, intent(out) ::	cilm(2, py_lmax+1, py_lmax+1)
			integer, intent(out)    ::	lmax
			integer, intent(in), optional :: norm, sampling, csphase			
		end subroutine SHExpandDHC
		
		subroutine MakeGridDHC(griddh, py_m, py_n, cilm, lmax, csphase, norm, sampling)
			integer, intent(in)     :: 	lmax, py_m, py_n, csphase
			complex*16, intent(in)  :: 	cilm(2, lmax+1, lmax+1)
			complex*16, intent(out) ::	griddh(py_m, py_n)
			integer, intent(in), optional :: norm, sampling
		end subroutine MakeGridDHC
		
		subroutine MakeGridGLQC(gridglq, cilm, lmax, plx, zero, norm, csphase, lmax_calc)
			complex*16, intent(in) ::	cilm(:,:,:)
			real*8, intent(in), optional :: plx(:,:), zero(:)
			complex*16, intent(out) ::	gridglq(:,:)
			integer, intent(in) ::	lmax
			integer, intent(in), optional :: norm, csphase, lmax_calc
		end subroutine MakeGridGLQC
		
		subroutine SHExpandGLQC(cilm, lmax, gridglq, w, plx, zero, norm, csphase, lmax_calc)	
			real*8, intent(in) ::	w(:)
			complex*16, intent(in) ::	gridglq(:,:)
			real*8, intent(in), optional ::	plx(:,:), zero(:)
			complex*16, intent(out) :: 	cilm(:,:,:)
			integer, intent(in) ::	lmax
			integer, intent(in), optional :: norm, csphase, lmax_calc
		end subroutine SHExpandGLQC
		
		real*8 function SHPowerLC(c, l)
			complex*16, intent(in) :: c(:,:,:)
			integer, intent(in) :: l
		end function SHPowerLC
		
		real*8 function SHPowerDensityLC(c, l)
			complex*16, intent(in) :: 	c(:,:,:)
			integer, intent(in) :: 	l
		end function SHPowerDensityLC
		
		complex*16 function SHCrossPowerLC(c1, c2, l)
			complex*16, intent(in) :: 	c1(:,:,:), c2(:,:,:)
			integer, intent(in) :: 	l
		end function SHCrossPowerLC
		
		Complex*16 function SHCrossPowerDensityLC(c1, c2, l)
			complex*16, intent(in) :: 	c1(:,:,:), c2(:,:,:)
			integer, intent(in) :: 	l
		end function SHCrossPowerDensityLC
		
		subroutine SHPowerSpectrumC(c, lmax, spectra)
			complex*16, intent(in) :: 	c(:,:,:)
			integer, intent(in) :: 	lmax
			real*8, intent(out) ::	spectra(:)
		end subroutine SHPowerSpectrumC
		
		subroutine SHPowerSpectrumDensityC(c, lmax, spectra)
			complex*16, intent(in) :: 	c(:,:,:)
			integer, intent(in) :: 	lmax
			real*8, intent(out) ::	spectra(:)
		end subroutine SHPowerSpectrumDensityC
		
		subroutine SHCrossPowerSpectrumC(c1, c2, lmax, cspectra)
			complex*16, intent(in) :: 	c1(:,:,:), c2(:,:,:)
			integer, intent(in) :: 	lmax
			complex*16, intent(out) ::	cspectra(:)
		end subroutine SHCrossPowerSpectrumC
		
		subroutine SHCrossPowerSpectrumDensityC(c1, c2, lmax, cspectra)
			complex*16, intent(in) :: 	c1(:,:,:), c2(:,:,:)
			integer, intent(in) :: 	lmax
			complex*16, intent(out) ::	cspectra(:)
		end subroutine SHCrossPowerSpectrumDensityC

		subroutine SHCilmToVector(cilm, vector, lmax)
			real*8, intent(in) ::	cilm(:,:,:)
			real*8, intent(out) ::	vector(:)
			integer, intent(in) ::	lmax
		end subroutine SHCilmToVector
		
		subroutine SHVectorToCilm(vector, cilm, lmax)
			real*8, intent(out) ::	cilm(:,:,:)
			real*8, intent(in) ::	vector(:)
			integer, intent(in) ::	lmax
		end subroutine SHVectorToCilm

		integer function YilmIndex(i, l, m)
			integer, intent(in) ::	i, l, m
		end function YilmIndex
		
		subroutine ComputeDMap(Dij, dh_mask, n_dh, sampling, lmax)
			real*8, intent(out) ::	Dij(:,:)
			integer, intent(in) ::	dh_mask(:,:), n_dh, sampling, lmax
		end subroutine ComputeDMap
		
	end interface
	
end module SHTOOLS

