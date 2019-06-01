!**************************************
!*  calpotbias.f90 Ver.1.4 '11.01.31  *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine calpotbias(xcel,ycel,zcel,   &
     &                npoly,nwater,nmatom, &
     &                force,pot_pbias)

  use interface_mdtech

  use md_global
  use mpi_global

  implicit none

!     subroutine to calculate bias force & potential

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  integer,intent(in):: npoly           ! all number of poly
  integer,intent(in):: nwater          ! number of H2O molecules
  integer,intent(in):: nmatom          ! number of monatomic molecules

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! force calculated here

!   OUTPUT
  real(8),intent(out):: pot_pbias       ! bias potential

! LOCAL:
  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

  real(8):: molecom(3,nmole)            ! center of mass of each molecule

  real(8):: delta_pot                   ! bias potential calculated 
                                        !  in this subroutine

  real(8):: rij(3)                      ! Ri-Rref
!  real(8):: rij_abs2                    ! |Rij|^2
  real(8):: fij(3)                      ! temporal storage for 
                                        !  bias force(i,ref)
  real(8):: fij_tmp
  real(8):: fix_tmp, fiy_tmp

  integer:: k                           ! do loop index 
  integer:: iatom                       ! atom indexes 

  integer:: imole                       ! molecular indexes 

!      -------------------------------------------------
!      POSITION BIAS FORCE & ENERGY
!         
!      (1) Harmonic type
!          force due to i, ref bond: Fiz = -2*K*(Rzi-Rref)
!          and potential energy   :  Ei  =    K*(Rzi-Rref)**2 
!
!      (2) Harmonic type based on COM position
!          force due to i, ref bond: Fiz = -2*K*(Rgzi-Rref)*mi/Mg
!          and potential energy   :  Egi =    K*(Rgzi-Rref)**2 *mi/Mg
!                                                        (only for i)
!
!      (3) Soft wall type
!          force due to i, ref bond: Fiz = | -3K*(z - z+)**2   if z > z+
!                                          | 0                 if z- < z < z+
!                                          | +3K*(z- - z)**2   if z < z-
!
!          and potential energy   :  Ei  = | K*(z - z+)**3   if z > z+
!                                          | 0               if z- < z < z+
!                                          | K*(z- - z)**3   if z < z-
!
!      (4) 2D Harmonic type
!          force due to i, ref bond: Fix = -2*K*(Rxi-Rxref)
!                                    Fiy = -2*K*(Ryi-Ryref)
!          and potential energy   :  Ei  =    K*{(Rxi-Rxref)**2+(Ryi-Ryref)**2}
!
!      (5) 2D Harmonic type based on COM position
!          force due to i, ref bond: Fix = -2*K*(Rgxi-Rxref)*mi/Mg
!                                    Fiy = -2*K*(Rgyi-Ryref)*mi/Mg
!          and potential energy   :  Ei  = K*{(Rgxi-Rxref)**2+(Rgyi-Ryref)**2}
!                                          * mi/Mg   (only for i)
!
!      -------------------------------------------------

!     +     +     +     +     +     +     +

!---- initialization
  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

! --- GRAND LOOP TO CALCULATE BONDED FORCE ---

  pot_pbias = 0.0d0

  looplast = npotbias

! - Harmonic type
  IF (potbias_typ == 'HARM') THEN

! MPI loop
!      DO k = 1, npotbias
     DO k = loopinit, looplast, loopstep
        
!       - calculate Rij -

        iatom   = index_potbias(k)      ! atom index

        rij(3)  = atmcor(3,iatom) - para_potbias(1)

!       --- periodic boundary ---

!        rij(1) = rij(1) - box(1) * anint(rij(1)*box_inv(1))
!        rij(2) = rij(2) - box(2) * anint(rij(2)*box_inv(2))
        rij(3) = rij(3) - box(3) * anint(rij(3)*box_inv(3))

!        rij_abs2 = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

!       - calculate potential -

        fij_tmp = para_potbias(2) * rij(3) ! = K * (rz_i - rz_ref)

        delta_pot = fij_tmp * rij(3)    ! = K * (rz_i - rz_ref)^2

        pot_pbias  = pot_pbias + delta_pot

!       - calculate Fij -

        fij(3)  = -2.0d0 * fij_tmp      ! = -2.0 K * (rz_i - rz_ref)

!       - add to force -

!        force(1,iatom) = force(1,iatom) + fij(1)
!        force(2,iatom) = force(2,iatom) + fij(2)
        force(3,iatom) = force(3,iatom) + fij(3)
         
     END DO

! - Harmonic type based on COM
  ELSE IF (potbias_typ == 'HMCM') THEN

     !---- center of mass of each molecule
     call calcom( npoly,nwater,nmatom, &
     &         xcel,ycel,zcel, &
     &         molecom)

! MPI loop
!      DO k = 1, npotbias
     DO k = loopinit, looplast, loopstep
        
!       - calculate Rij -

        iatom   = index_potbias(k)      ! atom index
        imole   = irmolept_list(iatom)  ! molecular index

        rij(3)  = molecom(3,imole) - para_potbias(1)
!        rij(3)  = atmcor(3,iatom) - para_potbias(1)

!       --- periodic boundary ---

!        rij(1) = rij(1) - box(1) * anint(rij(1)*box_inv(1))
!        rij(2) = rij(2) - box(2) * anint(rij(2)*box_inv(2))
        rij(3) = rij(3) - box(3) * anint(rij(3)*box_inv(3))

!        rij_abs2 = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

!       - calculate potential -

        fij_tmp = para_potbias(2) * rij(3) * atmmass(iatom) / molmass(imole)
                                    ! = K * (rz_g - rz_ref) * mi / Mmol

        delta_pot = fij_tmp * rij(3)   ! = K * (rz_g - rz_ref)^2 * mi / Mmol

        pot_pbias  = pot_pbias + delta_pot

!       - calculate Fij -

        fij(3)  = -2.0d0 * fij_tmp     ! = -2.0 K * (rz_g - rz_ref) * mi / Mmol

!       - add to force -

!        force(1,iatom) = force(1,iatom) + fij(1)
!        force(2,iatom) = force(2,iatom) + fij(2)
        force(3,iatom) = force(3,iatom) + fij(3)
         
     END DO

! - Soft wall type
  ELSE IF (potbias_typ == 'SWAL') THEN

! MPI loop
!      DO k = 1, npotbias
     DO k = loopinit, looplast, loopstep
        
!       - if z is inside the bin -

        iatom   = index_potbias(k)      ! atom index

        if (atmcor(3,iatom) > para_potbias(2)) then ! z > z+

           rij(3) = atmcor(3,iatom) - para_potbias(2) ! = z - z+
           fij_tmp = para_potbias(3) * rij(3)*rij(3) ! = K * (z - z+)^2
               
           delta_pot = fij_tmp * rij(3) ! = K * (z - z+)^3

           fij(3) = -3.0d0 * fij_tmp    ! = -3K * (z - z+)^2

        else if (atmcor(3,iatom) < para_potbias(1)) then ! z < z-

           rij(3) = para_potbias(1) - atmcor(3,iatom) ! = z- - z
           fij_tmp = para_potbias(3) * rij(3)*rij(3) ! = K * (z- - z)^2
               
           delta_pot = fij_tmp * rij(3) ! = K * (z- - z)^3

           fij(3) = 3.0d0 * fij_tmp     ! = 3K * (z- - z)^2

        else                            ! z- < z < z+

           delta_pot = 0.0d0
           fij(3) = 0.0d0

        end if

!       - calculate potential -

        pot_pbias  = pot_pbias + delta_pot

!       - add to force -

!        force(1,iatom) = force(1,iatom) + fij(1)
!        force(2,iatom) = force(2,iatom) + fij(2)
        force(3,iatom) = force(3,iatom) + fij(3)
         
     END DO

! - 2D harmonic type
  ELSE IF (potbias_typ == 'HM2D') THEN

! MPI loop
!      DO k = 1, npotbias
     DO k = loopinit, looplast, loopstep

!       - calculate Rij -

        iatom   = index_potbias(k)      ! atom index

        rij(1)  = atmcor(1,iatom) - para_potbias(1)
        rij(2)  = atmcor(2,iatom) - para_potbias(2)

!       --- periodic boundary ---

        rij(1) = rij(1) - box(1) * anint(rij(1)*box_inv(1))
        rij(2) = rij(2) - box(2) * anint(rij(2)*box_inv(2))
!        rij(3) = rij(3) - box(3) * anint(rij(3)*box_inv(3))

!        rij_abs2 = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

!       - calculate potential -

        fix_tmp = para_potbias(3) * rij(1)   ! = K * (rx_i - rx_ref)
        fiy_tmp = para_potbias(3) * rij(2)   ! = K * (ry_i - ry_ref)

        delta_pot = fix_tmp*rij(1) + fiy_tmp*rij(2)
                                ! = K * {(rx_i - rx_ref)^2 + (ry_i - ry_ref)^2}

        pot_pbias  = pot_pbias + delta_pot

!       - calculate Fij -

        fij(1)  = -2.0d0 * fix_tmp   ! = -2.0 K * (rx_i - rx_ref)
        fij(2)  = -2.0d0 * fiy_tmp   ! = -2.0 K * (ry_i - ry_ref)

!       - add to force -

        force(1,iatom) = force(1,iatom) + fij(1)
        force(2,iatom) = force(2,iatom) + fij(2)
!        force(3,iatom) = force(3,iatom) + fij(3)
         
     END DO

! - 2D harmonic type based on COM
  ELSE IF (potbias_typ == 'HC2D') THEN

     !---- center of mass of each molecule
     call calcom( npoly,nwater,nmatom, &
          &       xcel,ycel,zcel, &
          &       molecom)

! MPI loop
!      DO k = 1, npotbias
     DO k = loopinit, looplast, loopstep

!       - calculate Rij -

        iatom   = index_potbias(k)      ! atom index
        imole   = irmolept_list(iatom)  ! molecular index

        rij(1)  = molecom(1,imole) - para_potbias(1)
        rij(2)  = molecom(2,imole) - para_potbias(2)

!       --- periodic boundary ---

        rij(1) = rij(1) - box(1) * anint(rij(1)*box_inv(1))
        rij(2) = rij(2) - box(2) * anint(rij(2)*box_inv(2))
!        rij(3) = rij(3) - box(3) * anint(rij(3)*box_inv(3))

!        rij_abs2 = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

!       - calculate potential -

        fix_tmp = para_potbias(3) * rij(1) * atmmass(iatom) / molmass(imole)
                                   ! = K * (rx_g - rx_ref) * mi/Mg
        fiy_tmp = para_potbias(3) * rij(2) * atmmass(iatom) / molmass(imole)
                                   ! = K * (ry_g - ry_ref) * mi/Mg


        delta_pot = fix_tmp*rij(1) + fiy_tmp*rij(2)
                        ! = K * {(rx_g - rx_ref)^2 + (ry_g - ry_ref)^2} * mi/Mg

        pot_pbias  = pot_pbias + delta_pot

!       - calculate Fij -

        fij(1)  = -2.0d0 * fix_tmp   ! = -2.0 K * (rx_g - rx_ref) * mi/Mg
        fij(2)  = -2.0d0 * fiy_tmp   ! = -2.0 K * (ry_g - ry_ref) * mi/Mg

!       - add to force -

        force(1,iatom) = force(1,iatom) + fij(1)
        force(2,iatom) = force(2,iatom) + fij(2)
!        force(3,iatom) = force(3,iatom) + fij(3)
         
     END DO

  END IF

!     +     +     +     +     +     +     +

end subroutine calpotbias
