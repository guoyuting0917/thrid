!*****************************
!*  spline_interp.f Ver.1.3  *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-02-23 20:47:44 gota>
!
! description:
!   methods for spline interpolation of electric interaction
!
! note: this routine is for spline interpolation of electric interaction
!       Spline order is limited to 3rd order
!
subroutine spline_interp_init( alpha, &
     &                         rrcut )

  use md_global
#if defined(_HF_SKIP_EWALDR)
  use cstmnb, only: rrcut_hfewr
#endif

!C-MPI#if defined(MPI)
  use mpi_global
!C-MPI#endif

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: alpha            ! ewald parameter alpha [non-dimension]
  real(8),intent(in):: rrcut            ! ewald real space cutoff length [non-d]

! LOCAL:
  integer:: i

  real(8):: arg
  real(8):: sqrDist
  real(8):: invDist

  real(8):: ROOT_PI
  real(8):: alphaEwald
  real(8):: piEwald

  real(8):: rrcut_max

!     FUNCTIONS
  real*8:: DERFC
!      intrinsic DERFC

!     +     +     +     +     +     +     +

!---- initialization
  ROOT_PI = DSQRT(DACOS(-1.0d0))
  alphaEwald = alpha
  piEwald = 2.0d0 / ROOT_PI * alphaEwald

  rrcut_max = rrcut
#if defined(_HF_SKIP_EWALDR)
  if (rrcut_max < rrcut_hfewr) then
     rrcut_max = rrcut_hfewr
  end if
#endif

  if (irank == 0) then
     write(6,*) 'Preparing the spline interpolation table'
  end if

  if (nspltbl > max_interpol) then
     write(6,*) 'Error: number of spline interpolation points ', &
          &     'exceeds max_interpol'
     stop
  end if

  spltbl_int = (rrcut_max + 2.0d0) / DBLE(nspltbl)
  if (irank == 0) then
     write(6,*) '  table points: ', nspltbl
     write(6,'(A,2X,F9.6)') '  width [A]: ', spltbl_int
  end if
      
  spl_order = max_spl_order ! fixed to 3
  if (irank == 0) then
     write(6,*) '  spline order: ', spl_order
  end if

  do i = 1, nspltbl
     x_tbl(i) = DBLE(i-1) * spltbl_int
  end do

!---- calculate table value
!     x = 0 terminal forced to be 1.0
  spltbl_real(1) = 1.0d0
  spltbl_realpot(1) = 1.0d0
  spltbl_excs(1) = 1.0d0
  spltbl_excspot(1) = 1.0d0

  do i = 2, nspltbl
     sqrDist = x_tbl(i)
     invDist = 1.0d0 / sqrDist
     arg = alphaEwald * sqrDist

!        - Ewald real
     spltbl_real(i) = invDist * (piEwald * DEXP(-arg * arg) &
          &                    * invDist + DERFC(arg) &
          &                    * invDist * invDist)
     spltbl_realpot(i) = DERFC(arg) * invDist
!        - Ewald excess
     spltbl_excs(i) = invDist * (-piEwald * DEXP(-arg * arg) &
          &                    * invDist + (1.0d0 - DERFC(arg)) &
          &                    * invDist * invDist)
     spltbl_excspot(i) = (1.0d0 - DERFC(arg) ) * invDist
  end do

!---- calculate spline coefficient
!     - Ewald real

  call cal_spl_coeff( spltbl_real, &
       &              spl_b_real, spl_c_real, spl_d_real)

!     - Ewald realpot
  call cal_spl_coeff( spltbl_realpot, &
       &              spl_b_realpot, spl_c_realpot, spl_d_realpot )

!     - Ewald excess
  call cal_spl_coeff( spltbl_excs, &
       &              spl_b_excs,spl_c_excs,spl_d_excs)

!     - Ewald excess pot
  call cal_spl_coeff( spltbl_excspot, &
       &              spl_b_excspot, spl_c_excspot, spl_d_excspot)

  if (irank == 0) then
     write(6,*) '            ...finished'
  end if

  return
end subroutine spline_interp_init

!----------------------------------------------------------
subroutine cal_spl_coeff( spltbl, spl_b, spl_c, spl_d )

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: spltbl(*) ! interpolation data array

!     OUTPUT
  real(8),intent(out):: spl_b(*) ! spline coefficient b
  real(8),intent(out):: spl_c(*) ! spline coefficient c
  real(8),intent(out):: spl_d(*) ! spline coefficient d

! LOCAL:
  real(8):: alpha(nspltbl), l(nspltbl), mu(nspltbl), z(nspltbl)

  integer:: i

!     +     +     +     +     +     +     +

  do i = 2, nspltbl-1
     alpha(i) = 3.0d0 / spltbl_int &
          &   * ((spltbl(i+1) - spltbl(i)) &
          &    - (spltbl(i) - spltbl(i-1)))
  end do

  l(1) = 1.0d0
  mu(1) = 0.0d0
  z(1) = 0.0d0
  do i = 2, nspltbl-1
     l(i) = 2.0d0 * (x_tbl(i+1) - x_tbl(i-1)) - spltbl_int*mu(i-1)
     mu(i) = spltbl_int / l(i)
     z(i) = (alpha(i) - spltbl_int*z(i-1)) / l(i)
  end do

  l(nspltbl) = 1.0d0
  z(nspltbl) = 0.0d0
  spl_c(nspltbl) = 0.0d0
  do i = nspltbl-1, 1, -1
     spl_c(i) = z(i) - mu(i)*spl_c(i+1)
     spl_b(i) = (spltbl(i+1) - spltbl(i)) / spltbl_int &
          &   - spltbl_int * (spl_c(i+1) + 2.0d0*spl_c(i)) / 3.0d0
     spl_d(i) = (spl_c(i+1) - spl_c(i)) &
          &   / (3.0d0 * spltbl_int)
  end do

!     +     +     +     +     +     +     +

  return
end subroutine cal_spl_coeff

!----------------------------------------------------------
!
! functions for spline interpolations
!
function SPLINE_INTERPFUNC( x, spltable, spl_b, spl_c, spl_d ) 

  use md_global

  implicit none

! RESULT:
  real(8):: SPLINE_INTERPFUNC

! ARGUMENT:
!     INPUT
  real(8):: x
  real(8):: spltable(*)
  real(8):: spl_b(*)
  real(8):: spl_c(*)
  real(8):: spl_d(*)

! LOCAL:
  integer:: index_1
!      integer:: index_2
  real(8):: d

!     +     +     +     +     +     +     +

  if (x >= spltbl_int*DBLE(nspltbl-1)) then

     SPLINE_INTERPFUNC = 0.0d0
     return
         
  else

     index_1 = INT(x/spltbl_int) + 1
!     index_2 = index_1 + 1

     d = x - x_tbl(index_1)

     SPLINE_INTERPFUNC = spltable(index_1) + d*(spl_b(index_1) &
          &            + d*(spl_c(index_1) + d*spl_d(index_1)))

     return

  end if

!     +     +     +     +     +     +     +

  return
end function SPLINE_INTERPFUNC


