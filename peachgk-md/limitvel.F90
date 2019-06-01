!**************************
!*  limitvel.f90 Ver.1.0  *
!*     for peachgk_md.f   *
!*         by G.Kikugawa  *
!**************************
! Time-stamp: <2015-01-06 18:32:19 gota>

subroutine limitvel(dt_short_cal, &
     &              limitdist)

  use md_global

  implicit none

!     This subroutine limits the maximum velocity to a specific value.
!     When the atomic distance is very close or molecular distorsion is
!     significant after initial configuration generated,
!     this algorithm works well for the purpose of structural relaxation.
!
! ARGUMENTS:
!     INPUT
  real(8),intent(in):: dt_short_cal     ! time step of short force
  real(8),intent(in):: limitdist  ! maximum atomic displacement when doing
                              ! time integration [non-d] (structure relaxation)

! LOCAL:
  real(8):: atmvel_limit         ! limit of atomic velocity
  real(8):: atmvel_limit2        ! atmvel_limit**2
  real(8):: norm_vel             ! temporary norm of atomic velocity
  real(8):: s_fact               ! factor for velocity rescale

  integer:: i                    ! loop index

!     +     +     +     +     +     +     +

!----- some preparation -----

  atmvel_limit = limitdist / dt_short_cal
  atmvel_limit2 = atmvel_limit * atmvel_limit

!----- impose maximum velocity when velocity is too large

  do i = 1, natom

     norm_vel = atmvel(1,i)*atmvel(1,i) &
          &   + atmvel(2,i)*atmvel(2,i) &
          &   + atmvel(3,i)*atmvel(3,i)

     if (atmvel_limit2 < norm_vel) then
        norm_vel = sqrt(norm_vel)   ! actually calculate norm of the velocity
        s_fact = atmvel_limit / norm_vel
                                ! coefficient of velocity rescaling
        atmvel(1:3,i) = s_fact * atmvel(1:3,i)
     end if

  end do

!     +     +     +     +     +     +     +

end subroutine limitvel
