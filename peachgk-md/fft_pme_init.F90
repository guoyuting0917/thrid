!**************************************
!*  fft_pme_init.f Ver.1.4 '10.06.30  *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************

!---- original code Copyright
!*-*
!*-* Copyright (C) 2000 Massimo Marchi and Piero Procacci
!*-* Full copyright notice at http://www.chim.unifi.it/orac/copyright4.0.html
!*-* Contact for info M. Marchi, CEA,  Gif Sur Yvette 91191 (FRANCE) 
!*-* Email:marchi@villon.saclay.cea.fr
!*-* 
subroutine fft_pme_init( nfft1, nfft2, nfft3, pme_order, &
     &                   alpha, &
     &                   pot_ewc )

!************************************************************************
!*
!*     To be called from main_md
!*---  ON INPUT
!*     nfft1,nfft2,nfft3: grid points in the k1,k2,k3 directions
!*     order            : order of B-spline interpolation
!*---  ON OUTPUT 
!*     sizfftab is permanent 3d fft table storage
!*     sizffwrk is temporary 3d fft work storage
!*     siztheta is size of arrays theta1-3 dtheta1-3
!*     sizheap is total size of permanent storage
!*     sizstack is total size of temporary storage
!*     bsp_mod1-3 hold the moduli of the inverse DFT of the B splines
!*
!************************************************************************

  use md_global
  use spme_global
#if defined(MPI)
  use mpi_global
#endif

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: nfft1, nfft2, nfft3 ! grid points in SPME
  integer,intent(in):: pme_order       ! B-spline order

  real(8),intent(in):: alpha            ! parameter alpha [1/m]

  real(8),intent(out):: pot_ewc          ! potential of ewald self-energy

! LOCAL:
  integer,parameter:: MAXNFFT = 1000

  integer:: nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw

  real(8):: array(maxord), darray(maxord), w
  real(8):: bsp_arr(MAXNFFT)

  integer:: i,maxn_t

  real*8:: pi

!     - PME MPI
#if defined(MPI)
  integer:: n
#endif

!     +     +     +     +     +     +     +     +     +

!     compute the above output parameters needed for 
!     heap or stack allocation.

  call get_fftdims( nfft1, nfft2, nfft3, &
       &            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
       &            sizfftab, sizffwrk )

  siztheta = natom * pme_order
  siz_Q = 2 * nfftdim1 * nfftdim2 * nfftdim3
  sizheap =  nfft1 + nfft2 + nfft3 + sizfftab
  sizstack = siz_Q+6*siztheta+sizffwrk+3*natom

  call get_fftdims( nfft1, nfft2, nfft3, &
       &            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw)

!     loads the moduli of the inverse DFT of the B splines
!     bsp_mod1-3 hold these values, nfft1-3 are the grid dimensions,
!     Order is the order of the B spline approx.

  if ( pme_order > maxord )then
     write(6,*)'order too large! check on MAXORDER'
     stop
  endif
  maxn_t = max(nfft1,nfft2,nfft3)
  if ( maxn_t > maxn )then 
     write(6,*)'nfft1-3 too large! check on MAXNFFT'
     stop
  endif

  w = 0.d0
  call fill_bspline( w, pme_order, array, darray )
  bsp_arr(1:maxn_t) = 0.d0

  do i = 2,pme_order+1
     bsp_arr(i) = array(i-1)
  end do
  call DFTMOD( bsp_mod1, bsp_arr, nfft1 )
  call DFTMOD( bsp_mod2, bsp_arr, nfft2 )
  call DFTMOD( bsp_mod3, bsp_arr, nfft3 )

#if defined(_FFTW3)
!!! SPME by using FFTW3
#if defined(MPI)
  if (irank == 0) then
#endif
     write(6,*)'using FFTW Ver.3.x'
#if defined(MPI)
  end if
#endif
  call fftw3_init( nfft1, nfft2, nfft3 )

#else
!!! SPME by using public fft

!!! MPI version is not supported for public fft
#if defined(MPI)
  write(6,*) 'Error: Use FFTW for PME MPI calculation'
  stop
#endif

#if defined(MPI)
  if (irank == 0) then
#endif
     write(6,*) 'using public domain fft code'
#if defined(MPI)
  end if
#endif
  call cffti( nfft1, fftable(1) )
  call cffti( nfft2, fftable(nfftable+1) )
  call cffti( nfft3, fftable(nfftable*2+1) )

#endif

!   --- CALCULATE CONSTANT PART OF EWALD POTENTIAL ---
!          -1/sqrt(pi)*alpha * sigmj(qj^2)          
! 

  pi = dacos(-1.0d0)

  pot_ewc = 0.0d0
        
  do i = 1, natom
     pot_ewc = pot_ewc + atmchrg(i)**2
  end do

!---- use at standard ewald
  pot_ewc = -1.0d0/sqrt(pi)*alpha * pot_ewc 


!---- PME MPI initialization
#if defined(MPI)
!     - calculate FFT range for x- and z-direction
  n = nfft1 / nproc
  xlimitmax = n + 1

  do i = 0, nproc - 1
     xlimit1(i) = n * i + max(mod(nfft1, nproc) - nproc + i, 0) + 1
     xlimit2(i) = n * (i + 1) + max(mod(nfft1, nproc) - nproc + i + 1, 0)
  end do

  n = nfft3 / nproc
  zlimitmax = n + 1
      
  do i = 0, nproc - 1
     zlimit1(i) = n * i + max(mod(nfft3, nproc) - nproc + i, 0) + 1
     zlimit2(i) = n * (i + 1) + max(mod(nfft3, nproc) - nproc + i + 1, 0)
  end do
#endif

!     +     +     +     +     +     +     +     +     +

  return
end subroutine fft_pme_init

!   FFT CALLS
!------------------------------------------------------
subroutine get_fftdims( nfft1, nfft2, nfft3, &
     &                  nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
     &                  sizfftab, sizffwrk )
  implicit none
  integer,intent(in):: nfft1,nfft2,nfft3
  integer,intent(out):: nfftdim1,nfftdim2,nfftdim3
  integer,intent(out):: nfftable,nffwork,sizfftab,sizffwrk
  integer:: n
  integer:: nfftmax

  nfftmax = max(nfft1,nfft2,nfft3)
  nfftdim1 = nfft1
  n = nfft1/2
  if ( nfft1 == 2*n )nfftdim1 = nfft1+1
  nfftdim2 = nfft2
  n = nfft2/2
  if ( nfft2 == 2*n )nfftdim2 = nfft2+1
  nfftdim3 = nfft3
  n = nfft3/2
  if ( nfft3 == 2*n )nfftdim3 = nfft3+1

  nfftable = 4*nfftmax + 15
  nffwork = nfftmax
  sizfftab = 3*nfftable
  sizffwrk  = 2*nfftmax

  return
end subroutine get_fftdims

!------------------------------------------------------
subroutine fill_bspline( w, order, array, darray )
!---------- use standard B-spline recursions: see doc file
  implicit none
  integer,intent(in):: order
  real(8),intent(in):: w
  real(8),intent(out):: array(order), darray(order)

  integer:: k
  integer:: j
  real(8):: div


! do linear case
  array(order) = 0.d0
  array(2) = w
  array(1) = 1.d0 - w

! compute standard b-spline recursion
  do k = 3,order-1

     div = 1.d0 / (k-1)
     array(k) = div*w*array(k-1)
     do j = 1, k-2
        array(k-j) = div*((w+j)*array(k-j-1) + (k-j-w)*array(k-j))
     end do
     array(1) = div*(1-w)*array(1)
         
  end do

! perform standard b-spline differentiation
  darray(1) = -array(1)
  do j = 2, order
     darray(j) = array(j-1) - array(j)
  end do

! one more recursion
  div = 1.d0 / (order-1)
  array(order) = div*w*array(order-1)
  do j = 1,order-2
     array(order-j) =  div*((w+j)*array(order-j-1) + (order-j-w)*array(order-j))
  end do
  array(1) = div*(1-w)*array(1)

  return
end subroutine fill_bspline

!---------------------------------------------------
subroutine DFTMOD( bsp_mod, bsp_arr, nfft)
  implicit none
  integer,intent(in):: nfft
  real(8),intent(in):: bsp_arr(nfft)
  real(8),intent(out):: bsp_mod(nfft)
! Computes the modulus of the discrete fourier transform of bsp_arr,
!  storing it into bsp_mod

  integer:: j,k
  real(8):: sum1,sum2,twopi,arg,tiny
  twopi = 2.d0*3.14159265358979323846d0
  tiny = 1.d-7
  do k = 1,nfft
     sum1 = 0.d0
     sum2 = 0.d0
     do j = 1,nfft
        arg = twopi*(k-1)*(j-1)/nfft
        sum1 = sum1 + bsp_arr(j)*dcos(arg)
        sum2 = sum2 + bsp_arr(j)*dsin(arg)
     end do
     bsp_mod(k) = sum1**2 + sum2**2
  end do
  do k = 1,nfft
     if ( bsp_mod(k) < tiny ) bsp_mod(k) = 0.5d0*(bsp_mod(k-1) + bsp_mod(k+1))
  end do
  return
end subroutine DFTMOD

!----------------------------------------------------------
#if defined(_FFTW3)
subroutine fftw3_init( nfft1, nfft2, nfft3 )

  use spme_global

  implicit none

  include 'fftw3.f'

! ARGUMENT:
!     INPUT
  integer,intent(in):: nfft1, nfft2, nfft3 ! grid points in SPME

! LOCAL:

!     +     +     +     +     +     +     +     +     +

!
!     combine 1d routines for 3D-FFT
!

#if defined(_FFTW3_SLOW)
  call DFFTW_PLAN_DFT_1D( planF(1), nfft1, fftWork1, fftWork1, &
       &                  FFTW_FORWARD, FFTW_MEASURE )
  call DFFTW_PLAN_DFT_1D( planF(2), nfft2, fftWork1, fftWork1, &
       &                  FFTW_FORWARD, FFTW_MEASURE)
  call DFFTW_PLAN_DFT_1D( planF(3), nfft3, fftWork1, fftWork1, &
       &                  FFTW_FORWARD, FFTW_MEASURE)

  call DFFTW_PLAN_DFT_1D( planB(1), nfft1, fftWork1, fftWork1, &
       &                  FFTW_BACKWARD, FFTW_MEASURE )
  call DFFTW_PLAN_DFT_1D( planB(2), nfft2, fftWork1, fftWork1, &
       &                       FFTW_BACKWARD, FFTW_MEASURE )
  call DFFTW_PLAN_DFT_1D( planB(3), nfft3, fftWork1, fftWork1, &
       &                  FFTW_BACKWARD, FFTW_MEASURE)

#else
  call DFFTW_PLAN_DFT_1D( planF(1), nfft1, fftWork1, fftWork1, &
       &                  FFTW_FORWARD, FFTW_EXHAUSTIVE)
  call DFFTW_PLAN_DFT_1D( planF(2), nfft2, fftWork1, fftWork1, &
       &                  FFTW_FORWARD, FFTW_EXHAUSTIVE)
  call DFFTW_PLAN_DFT_1D( planF(3), nfft3, fftWork1, fftWork1, &
       &                  FFTW_FORWARD, FFTW_EXHAUSTIVE)

  call DFFTW_PLAN_DFT_1D( planB(1), nfft1, fftWork1, fftWork1, &
       &                  FFTW_BACKWARD, FFTW_EXHAUSTIVE)
  call DFFTW_PLAN_DFT_1D( planB(2), nfft2, fftWork1, fftWork1, &
       &                  FFTW_BACKWARD, FFTW_EXHAUSTIVE)
  call DFFTW_PLAN_DFT_1D( planB(3), nfft3, fftWork1, fftWork1, &
       &                  FFTW_BACKWARD, FFTW_EXHAUSTIVE)

#endif

!     +     +     +     +     +     +     +     +     +

  return
end subroutine fftw3_init
#endif
