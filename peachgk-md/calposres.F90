!**************************************
!*  calposres.f Ver.1.3 '13.12.10     *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine calposres(xcel,ycel,zcel, &
     &               force,pot_posres)

  use md_global
  use mpi_global

  implicit none

!     subroutine to calculate position restraint force & potential 

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! force calculated here

!   OUTPUT
  real(8),intent(out):: pot_posres      ! position restraint potential

! LOCAL:
  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

  real(8):: delta_pot                  ! position restraint potential calculated 
                                       !  in this subroutine

  real(8):: rij(3)                      ! Ri-Rref
  real(8):: rij_abs2                    ! |Rij|^2
  real(8):: fij(3)                      ! temporal storage for 
                                        !  position restraint force(i,ref)
!  real(8):: fij_abs                     ! |Fij|                   

  integer:: i,k                         ! do loop index 
  integer:: iatom                       ! atom indexes 

!      -------------------------------------------------
!      POSITION RESTRAINT FORCE & ENERGY
!         
!          force due to i, ref bond:  Fij = -2*K*(Ri-Rref)
!          and potential energy   :  Eij =    K*(Ri-Rref)**2 
!
!      -------------------------------------------------

!     +     +     +     +     +     +     +

!---- initialization
  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

! --- GRAND LOOP TO CALCULATE BONDED FORCE ---

  pot_posres = 0.0d0


! loop for directions
  DO i = 1, 3
     looplast = nposres(i)

! MPI loop
!      DO k = 1, nposres(i)
     DO k = loopinit, looplast, loopstep
        
!       - calculate Rij -

        iatom   = index_posres(k,i) ! atom index

        rij(i)  = atmcor(i,iatom) - ref_atmcor(i,iatom)

!       --- periodic boundary ---

        rij(i) = rij(i) - box(i) * anint(rij(i)*box_inv(i))

        rij_abs2 = rij(i)*rij(i)

!       - calculate potential -

        delta_pot =  para_posres * rij_abs2
        pot_posres  = pot_posres + delta_pot

!       - calculate Fij -

        fij(i)  = -2.0d0 * para_posres * rij(i)

!       - add to force -

        force(i,iatom) = force(i,iatom) + fij(i)

     END DO

  END DO

!     +     +     +     +     +     +     +

end subroutine calposres
