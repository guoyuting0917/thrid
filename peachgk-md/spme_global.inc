!**************************************
!*  spme_global.inc Ver.1.0  09.11.11 *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
module spme_global

  implicit none

#include "config.h"

  !---- maximux number of arrays
  !     declare constant parameter
  integer,parameter:: maxfft1 = _MAXFFT1_+1   ! maximum number of pme grid
  integer,parameter:: maxfft2 = _MAXFFT2_+1
  integer,parameter:: maxfft3 = _MAXFFT3_+1
  integer,parameter:: maxord = _MAXORDER_   ! maximum number of B-spline order
  integer,parameter:: maxt = 2*maxfft1*maxfft2*maxfft3
  integer,parameter:: max_atm = _MAXNATOM_
  integer,parameter:: mth = maxord*max_atm
  integer,parameter:: maxn = _MAX_GRID_+1

  !     declare global variables
  real(8),save:: bsp_mod1(maxfft1),bsp_mod2(maxfft2),bsp_mod3(maxfft3)
                                   ! =1/b(1), 1/b(2), 1/b(3)
  integer,save:: sizfftab,sizffwrk
  integer,save:: siztheta
  integer,save:: siz_Q
  integer,save:: sizheap,sizstack

  real(8),save:: rkcut             ! kcut-off

#if defined(_FFTW3)
  complex(8),save:: fftWork1(maxn)   ! work array

  integer(8),save:: planForward    ! for 64 bit addressing mode
  integer(8),save:: planBackward   ! for 64 bit addressing mode
  integer(8),save:: planF(3)       ! for 64 bit addressing mode
  integer(8),save:: planB(3)       ! for 64 bit addressing mode
#else
  real(8),save:: fftable(3*(maxn*4+15)),ffwork(2*maxn)
                                   ! tables for fft
#endif

end module spme_global
