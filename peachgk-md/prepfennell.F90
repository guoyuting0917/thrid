!**************************************
!*  prepfennell.f Ver.1.0 '08.02.12   *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine prepfennell( pot_ewc, alpha, rrcut, &
     &                  iffennell)

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: alpha            ! parameter alpha [non-d]
  real(8),intent(in):: rrcut            ! ewald real space cutoff length [m]

  logical,intent(in):: iffennell       ! Fennell flag

!     OUTPUT
  real(8),intent(out):: pot_ewc          ! potential of ewald self-energy [non-d]
      
! LOCAL:
  real(8):: rtmp1            ! for temporal use
  real(8):: rtmp2            ! for temporal use
  real(8):: rcut_inv1        ! = 1/rcut
  real(8):: rcut_inv2        ! = 1/rcut^2
  real(8):: rtmp_erfc        ! = erfc(alpha*|rcut|)
  real(8):: const_fewr       ! = 2*alpha / sqrt(pi)

  real(8):: pi               ! = 3.14...

! FUNCTIONS:
  real(8):: derf             ! error function
!      external derf             ! error function
      
!     +     +     +     +     +     +     +

  if (iffennell) then
     fennell_flag = .true.
     pot_ewc = 0.0d0
  else
     fennell_flag = .false.
  end if

!---- some preparation
  pi = dacos(-1.0d0)
  const_fewr = 2.0d0 * alpha / dsqrt(pi) 

!---- calculate shift and damp term
!
  rcut_inv1 = 1.0d0 / rrcut
  rcut_inv2 = rcut_inv1 * rcut_inv1

  rtmp1 = rrcut * alpha     ! = |rcut|*alpha
  rtmp2 = rtmp1 * rtmp1     ! = |rcut|^2*alpha^2

  rtmp_erfc = 1.0d0 - derf(rtmp1) ! = erfc(|rcut|*alpha)

  fennell_shift = rtmp_erfc * rcut_inv1 ! = erfc(|rcut|*alpha) / |rcut|

  fennell_damp = rtmp_erfc * rcut_inv2 &
       &       + const_fewr * rcut_inv1 * dexp(-rtmp2)
                                ! = erfc(|rcut|*alpha) / |rcut|^2
                                ! + 2*alpha/sqrt(pi) * exp(-rtmp2) / |rcut|

!     +     +     +     +     +     +     +

  return
end subroutine prepfennell
