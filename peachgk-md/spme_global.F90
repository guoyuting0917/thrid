!**************************************
!* spme_global.inc Ver.1.0 09.11.11 *
!* for peachgk_md.f *
!* by G.Kikugawa *
!**************************************
module spme_global
  implicit none
  !---- maximux number of arrays
  ! declare constant parameter
  integer,parameter:: maxfft1 = 100 +1 ! maximum number of pme grid
  integer,parameter:: maxfft2 = 100 +1
  integer,parameter:: maxfft3 = 240 +1
  integer,parameter:: maxord = 10 ! maximum number of B-spline order
  integer,parameter:: maxt = 2*maxfft1*maxfft2*maxfft3
  integer,parameter:: max_atm = 30000
  integer,parameter:: mth = maxord*max_atm
  integer,parameter:: maxn = 240 +1
  ! declare global variables
  real(8),save:: bsp_mod1(maxfft1),bsp_mod2(maxfft2),bsp_mod3(maxfft3)
                                   ! =1/b(1), 1/b(2), 1/b(3)
  integer,save:: sizfftab,sizffwrk
  integer,save:: siztheta
  integer,save:: siz_Q
  integer,save:: sizheap,sizstack
  real(8),save:: rkcut ! kcut-off
  real(8),save:: fftable(3*(maxn*4+15)),ffwork(2*maxn)
                                   ! tables for fft
end module spme_global
