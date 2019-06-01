!**************************************
!*  accforce.f90 Ver.1.4 '11.01.31    *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine accforce( natom, atmmass, force )

  implicit none

!    subroutine to convert force to accel 

! ARGUMENTS: 
!     INPUT
  integer,intent(in):: natom           ! number of atoms
  real(8),intent(in):: atmmass(:)       ! atomic mass

!     OUTPUT
  real(8),intent(inout):: force(:,:)
                                ! input force
                                ! output accel        

! LOCAL:
  integer:: i                   ! do loop index

!     +     +     +     +     +     +     +

  do i=1,natom
     force(1:3,i) = force(1:3,i)/atmmass(i)
  end do

!     +     +     +     +     +     +     +

  return
end subroutine accforce
