!**************************************
!*  calangl.f Ver.2.0 '11.01.18       *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine calangl(xcel,ycel,zcel,   &
     &             force,pot_angl)

  use md_global
  use mpi_global

  implicit none

!    subroutine to calculate angl force & potential

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel           ! x cell length[non-d]
  real(8),intent(in):: ycel           ! y cell length[non-d]
  real(8),intent(in):: zcel           ! z cell length[non-d]

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! force calculated here

!   OUTPUT
  real(8),intent(out):: pot_angl        ! angle energy 

! LOCAL:

  real(8):: delta_pot                   ! potential calculated 
                                        !  in this subroutine

  real(8):: rij(3)                      ! Rj-Ri 
  real(8):: rij_abs                     ! |Rij| 
  real(8):: rkj(3)                      ! Rj-Rk 
  real(8):: rkj_abs                     ! |Rkj| 

  real(8):: theta                       ! angle (in radian) defined by 
                                        !  atoms i, j, k 
                                        !  0 =< theta < pi 
  real(8):: cos_theta                   ! cos(theta)         
  real(8):: sin_theta                   ! sin(theta)         

  real(8):: fi(3), fj(3), fk(3)         ! force to atoms (i,j,k)
                                        !  calculated here
  real(8):: f_tmp                       ! f = - f_tmp * sin_theta *
                                        !     d(theta)/dr
    

  integer:: n                           ! do loop index 
  integer:: iatom, jatom, katom         ! atom index

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size
      
!      -------------------------------------------------
!      ANGLE FORCE & ENERGY
!         
!          force due to (i,j,k) angl:
!             Fi = - 2*K*(THETA-THETAeq)
!                  *( Rkj/rkj - cosTHETA*Rij/rij)/rij/sinTHETA
!             Fk = - 2*K*(THETA-THETAeq)
!                  *( Rij/rij - cosTHETA*Rkj/rkj)/rkj/sinTHETA
!             Fj = - Fi - Fk
!
!          and potential energy:  
!             E  =  K*(THETA-THETAeq)**2
!
!      -------------------------------------------------


! --- GRAND LOOP TO CALCULATE ANGLE FORCE ---

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)
  
  pot_angl = 0.0d0

  looplast = nangl
      
! MPI loop
!      DO n = 1, nangl 
  DO n = loopinit, looplast, loopstep

!    --- if para_cangl = 0.0d0 then cycle ---

     if (abs(para_cangl(indexangltyp(n))) < 1.0d-32) CYCLE

!    - calculate Rij & Rkj -

     iatom   = iangl(n)                 ! atom index  
     jatom   = jangl(n)                 ! atom index  
     katom   = kangl(n)                 ! atom index  


     rij(1:3)  = atmcor(1:3,jatom) - atmcor(1:3,iatom)
     rij(1:3)  = rij(1:3) - box(1:3)*anint(rij(1:3)*box_inv(1:3))

     rij_abs   = rij(1)**2 + rij(2)**2 + rij(3)**2 
     rij_abs   = dsqrt(rij_abs)


     rkj(1:3)  = atmcor(1:3,jatom) - atmcor(1:3,katom)
     rkj(1:3)  = rkj(1:3) - box(1:3)*anint(rkj(1:3)*box_inv(1:3))

     rkj_abs = rkj(1)**2 + rkj(2)**2 + rkj(3)**2 
     rkj_abs = sqrt(rkj_abs)

!    - calculate theta(ijk) -

     cos_theta = (rij(1)*rkj(1) + rij(2)*rkj(2) + rij(3)*rkj(3))   &
          &    / (rij_abs*rkj_abs)

!    * adjust if cos_theta shows a value such as 1.000000001
     if (cos_theta > 1.0d0) then
        cos_theta =  1.0d0
     else if (cos_theta < -1.0d0) then
        cos_theta = -1.0d0
     end if

     theta     = dacos(cos_theta)       ! shown in radian
     sin_theta = dsin(theta)

!    - calculate potential -

     delta_pot = para_cangl(indexangltyp(n))   &
          &    * (theta - para_eqangl(indexangltyp(n)))**2 
     pot_angl = pot_angl + delta_pot

!    - calculate Fi, Fj, Fk -

     if (abs(sin_theta) > 1.0d-14) then
!       - normal case -
        f_tmp = -2.0d0 * para_cangl(indexangltyp(n))   &
            & * (theta - para_eqangl(indexangltyp(n))) / sin_theta

        fi(1:3) = ( rkj(1:3)/rkj_abs - rij(1:3)*cos_theta/rij_abs ) / rij_abs

        fi(1:3) = f_tmp * fi(1:3)

        fk(1:3) = ( rij(1:3)/rij_abs - rkj(1:3)*cos_theta/rkj_abs ) / rkj_abs

        fk(1:3) = f_tmp * fk(1:3)

        fj(1:3) = - fi(1:3) - fk(1:3)   ! Newton''s 3rd Law of motion

     else 
!       - if theta ~ 0 or pi, then fi,fj,fk are set to 0  -
        fi(1:3) = 0.0d0

        fj(1:3) = 0.0d0

        fk(1:3) = 0.0d0

     end if

!    - add to force -
     force(1:3,iatom) = force(1:3,iatom) + fi(1:3)
     force(1:3,jatom) = force(1:3,jatom) + fj(1:3)
     force(1:3,katom) = force(1:3,katom) + fk(1:3)

  END DO

!     +     +     +     +     +     +     +

end subroutine calangl

! ------------------------------------------------------------

subroutine calanglp(xcel,ycel,zcel,   &
     &              force,pot_angl,   &
     &              atm_viri_angl,atm_virit_angl)

!    subroutine to calculate angl force & potential

  use md_global
  use mpi_global

  implicit none

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel           ! x cell length[non-d]
  real(8),intent(in):: ycel           ! y cell length[non-d]
  real(8),intent(in):: zcel           ! z cell length[non-d]

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! force calculated here

!   OUTPUT
  real(8),intent(out):: pot_angl        ! angle energy 

  real(8),intent(out):: atm_viri_angl   ! virial(angle potential) of each atom
  real(8),intent(out):: atm_virit_angl(:,:) ! virial tensor (angle potential)

! LOCAL:

  real(8):: delta_pot                   ! potential calculated 
                                        !  in this subroutine

  real(8):: rij(3)                      ! Rj-Ri 
  real(8):: rij_abs                     ! |Rij| 
  real(8):: rkj(3)                      ! Rj-Rk 
  real(8):: rkj_abs                     ! |Rkj| 

  real(8):: theta                       ! angle (in radian) defined by 
                                        !  atoms i, j, k 
                                        !  0 =< theta < pi 
  real(8):: cos_theta                   ! cos(theta)         
  real(8):: sin_theta                   ! sin(theta)         

  real(8):: fi(3), fj(3), fk(3)         ! force to atoms (i,j,k)
                                        !  calculated here
  real(8):: f_tmp                       ! f = - f_tmp * sin_theta *
                                        !     d(theta)/dr
    

  integer:: n,l                         ! do loop index 
  integer:: iatom, jatom, katom         ! atom index

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size
      
!      -------------------------------------------------
!      ANGLE FORCE & ENERGY
!         
!          force due to (i,j,k) angl:
!             Fi = - 2*K*(THETA-THETAeq)
!                  *( Rkj/rkj - cosTHETA*Rij/rij)/rij/sinTHETA
!             Fk = - 2*K*(THETA-THETAeq)
!                  *( Rij/rij - cosTHETA*Rkj/rkj)/rkj/sinTHETA
!             Fj = - Fi - Fk
!
!          and potential energy:  
!             E  =  K*(THETA-THETAeq)**2
!
!      -------------------------------------------------


! --- GRAND LOOP TO CALCULATE ANGLE FORCE ---

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

  pot_angl = 0.0d0

  looplast = nangl

! MPI loop
!      DO n = 1, nangl
  DO n = loopinit, looplast, loopstep
        
!    --- if para_cangl = 0.0d0 then cycle ---

     if (abs(para_cangl(indexangltyp(n))) < 1.0d-32) CYCLE

!    - calculate Rij & Rkj -

     iatom   = iangl(n)                 ! atom index  
     jatom   = jangl(n)                 ! atom index  
     katom   = kangl(n)                 ! atom index  

     rij(1:3)  = atmcor(1:3,jatom) - atmcor(1:3,iatom)
     rij(1:3)  = rij(1:3) - box(1:3)*anint(rij(1:3)*box_inv(1:3))

     rij_abs   = rij(1)**2 + rij(2)**2 + rij(3)**2
     rij_abs   = dsqrt(rij_abs)

     rkj(1:3)  = atmcor(1:3,jatom) - atmcor(1:3,katom)
     rkj(1:3)  = rkj(1:3) - box(1:3)*anint(rkj(1:3)*box_inv(1:3))

     rkj_abs = rkj(1)**2 + rkj(2)**2 + rkj(3)**2 
     rkj_abs = sqrt(rkj_abs)

!    - calculate theta(ijk) -

     cos_theta = (rij(1)*rkj(1) + rij(2)*rkj(2) + rij(3)*rkj(3))   &
          &    / (rij_abs*rkj_abs)

!    * adjust if cos_theta shows a value such as 1.000000001
     if (cos_theta > 1.0d0) then
        cos_theta =  1.0d0
     else if (cos_theta < -1.0d0) then
        cos_theta = -1.0d0
     end if

     theta     = dacos(cos_theta)       ! shown in radian
     sin_theta = dsin(theta)

!    - calculate potential -

     delta_pot = para_cangl(indexangltyp(n))   &
          &    * (theta - para_eqangl(indexangltyp(n)))**2
     pot_angl = pot_angl + delta_pot

!    - calculate Fi, Fj, Fk -

     if (abs(sin_theta) > 1.0d-14) then
!       - normal case -
        f_tmp = -2.0d0 * para_cangl(indexangltyp(n))   &
            & * (theta - para_eqangl(indexangltyp(n))) / sin_theta

        fi(1:3) =  ( rkj(1:3)/rkj_abs - rij(1:3)*cos_theta/rij_abs ) / rij_abs

        fi(1:3) = f_tmp * fi(1:3)

        fk(1:3) =  ( rij(1:3)/rij_abs - rkj(1:3)*cos_theta/rkj_abs ) / rkj_abs

        fk(1:3) = f_tmp * fk(1:3)

        fj(1:3) = - fi(1:3) - fk(1:3)   ! Newton''s 3rd Law of motion

     else
!       - if theta ~ 0 or pi, then fi,fj,fk are set to 0  -
        fi(1:3) = 0.0d0

        fj(1:3) = 0.0d0

        fk(1:3) = 0.0d0

     end if

!    - add to force & virial -
     force(1:3,iatom) = force(1:3,iatom) + fi(1:3)

     force(1:3,jatom) = force(1:3,jatom) + fj(1:3)

     force(1:3,katom) = force(1:3,katom) + fk(1:3)

     atm_viri_angl = atm_viri_angl + fi(1)*atmcor(1,iatom)   &
          &        + fi(2)*atmcor(2,iatom) + fi(3)*atmcor(3,iatom)

     atm_viri_angl = atm_viri_angl + fj(1)*atmcor(1,jatom)   &
          &        + fj(2)*atmcor(2,jatom) + fj(3)*atmcor(3,jatom)

     atm_viri_angl = atm_viri_angl + fk(1)*atmcor(1,katom)   &
          &        + fk(2)*atmcor(2,katom) + fk(3)*atmcor(3,katom)

     do l=1,3
        atm_virit_angl(1:3,l) = atm_virit_angl(1:3,l)    &
             &                + fi(1:3)*atmcor(l,iatom)

        atm_virit_angl(1:3,l) = atm_virit_angl(1:3,l)   &
             &                + fj(1:3)*atmcor(l,jatom)

        atm_virit_angl(1:3,l) = atm_virit_angl(1:3,l)   &
             &                + fk(1:3)*atmcor(l,katom)
     end do

  END DO

!     +     +     +     +     +     +     +

end subroutine calanglp
