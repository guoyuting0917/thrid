!**************************************
!*  calbond.f Ver.2.0 '11.01.18       *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine calbond(xcel,ycel,zcel,   &
     &             force,pot_bond)
  
  use md_global
  use mpi_global

  implicit none

!     subroutine to calculate bonded force & potential 

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel           ! x cell length[non-d]
  real(8),intent(in):: ycel           ! y cell length[non-d]
  real(8),intent(in):: zcel           ! z cell length[non-d]

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! force calculated here

!   OUTPUT
  real(8),intent(out):: pot_bond        ! bond potential

! LOCAL:

  real(8):: delta_pot                   ! bonded potential calculated 
                                        !  in this subroutine

  real(8):: rij(3)                      ! Ri-Rj
  real(8):: rij_abs                     ! |Rij| 
  real(8):: fij(3)                      ! temporal storage for 
                                        !  bonding force(i,j)
  real(8):: fij_abs                     ! |Fij|                   

  integer:: k                           ! do loop index 
  integer:: iatom,jatom                 ! atom indexes 

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

!      -------------------------------------------------
!      BOND FORCE & ENERGY
!
!          force due to (i,j) bond:  Fij = -2*K*(Req-Rij)*Rij/|Rij| 
!          and potential energy   :  Eij =    K*(Req-Rij)**2 
!
!      -------------------------------------------------

! --- GRAND LOOP TO CALCULATE BONDED FORCE ---

  pot_bond = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

  looplast = nbond

!c MPI loop
!      DO k = 1, nbond
  DO k = loopinit, looplast, loopstep
        
!    --- if para_cbond = 0.0d0 then cycle ---

     if (abs(para_cbond(indexbondtyp(k))) < 1.0d-16) CYCLE
            
!    - calculate Rij -

     iatom   = ibond(k)                 ! atom index  
     jatom   = jbond(k)                 ! atom index

     rij(1:3)  = atmcor(1:3,iatom) - atmcor(1:3,jatom)
     rij(1:3)  = rij(1:3) - box(1:3)*anint(rij(1:3)*box_inv(1:3))

     rij_abs = rij(1)**2 + rij(2)**2 + rij(3)**2
     rij_abs = sqrt(rij_abs)

!    - calculate potential -

     delta_pot =  para_cbond(indexbondtyp(k))   &
          &    * (rij_abs - para_eqbond(indexbondtyp(k)))**2
     pot_bond  = pot_bond + delta_pot

!    - calculate Fij -

     fij_abs = -2.0d0 * para_cbond(indexbondtyp(k))   &
          &  * (rij_abs - para_eqbond(indexbondtyp(k)))

     fij(1:3)  = fij_abs * rij(1:3)/rij_abs

!    - add to force -

     force(1:3,iatom) = force(1:3,iatom) + fij(1:3)
     force(1:3,jatom) = force(1:3,jatom) - fij(1:3)

  END DO

!     +     +     +     +     +     +     +

end subroutine calbond

! --------------------------------------------------------

subroutine calbondp(xcel,ycel,zcel,   &
     &              force,pot_bond,   &
     &              atm_viri_bond,atm_virit_bond)

  use md_global
  use mpi_global

  implicit none

!     subroutine to calculate bonded force & potential 

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel           ! x cell length[non-d]
  real(8),intent(in):: ycel           ! y cell length[non-d]
  real(8),intent(in):: zcel           ! z cell length[non-d]

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! force calculated here

!   OUTPUT
  real(8),intent(out):: pot_bond        ! bond potential

  real(8),intent(out):: atm_viri_bond   ! virial of each atom
  real(8),intent(out):: atm_virit_bond(:,:) ! virial tensor (bond potential)

! LOCAL:

  real(8):: delta_pot                   ! bonded potential calculated 
                                        !  in this subroutine

  real(8):: rij(3)                      ! Ri-Rj
  real(8):: rij_abs                     ! |Rij| 
  real(8):: fij(3)                      ! temporal storage for 
                                        !  bonding force(i,j)
  real(8):: fij_abs                     ! |Fij|                   
      
  integer:: k,n                         ! do loop index 
  integer:: iatom,jatom                 ! atom indexes 

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

!      -------------------------------------------------
!      BOND FORCE & ENERGY
!         
!          force due to (i,j) bond:  Fij = -2*K*(Req-Rij)*Rij/|Rij| 
!          and potential energy   :  Eij =    K*(Req-Rij)**2 
!
!      -------------------------------------------------

! --- GRAND LOOP TO CALCULATE BONDED FORCE ---

  pot_bond = 0.0d0

  looplast = nbond

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

! MPI loop
!      DO k = 1, nbond
  DO k = loopinit, looplast, loopstep

!    --- if para_cbond = 0.0d0 then cycle ---

     if (abs(para_cbond(indexbondtyp(k))) < 1.0d-16) CYCLE
        
!    - calculate Rij -

     iatom   = ibond(k)                 ! atom index  
     jatom   = jbond(k)                 ! atom index

     rij(1:3)  = atmcor(1:3,iatom) - atmcor(1:3,jatom)
     rij(1:3)  = rij(1:3) - box(1:3)*anint(rij(1:3)*box_inv(1:3))

     rij_abs = rij(1)**2 + rij(2)**2 + rij(3)**2
     rij_abs = sqrt(rij_abs)

!    - calculate potential -

     delta_pot = para_cbond(indexbondtyp(k))   &
          &    * (rij_abs - para_eqbond(indexbondtyp(k)))**2
     pot_bond  = pot_bond + delta_pot

!    - calculate Fij -

     fij_abs = -2.0d0 * para_cbond(indexbondtyp(k))   &
          &  * (rij_abs - para_eqbond(indexbondtyp(k)))

     fij(1:3)  = fij_abs * rij(1:3)/rij_abs

!    - add to force & virial -

     force(1:3,iatom) = force(1:3,iatom) + fij(1:3)
     force(1:3,jatom) = force(1:3,jatom) - fij(1:3)

     atm_viri_bond = atm_viri_bond   &
          &        + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)
     do n=1,3
        atm_virit_bond(1:3,n) = atm_virit_bond(1:3,n)   &
     &                        + fij(1:3)*rij(n)
     end do

  END DO

!     +     +     +     +     +     +     +

end subroutine calbondp
