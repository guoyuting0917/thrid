!******************************
!*  caltorsim_hf.f90 Ver.1.4  *
!*      for peachgk_md.f      *
!*            by G.Kikugawa   *
!******************************
! Time-stamp: <>

subroutine caltorsimp_hf(xcel,ycel,zcel,   &
     &                   force,pot_torsim,   &
     &                   atm_viri_tors,atm_virit_tors,   &
     &                   pot_to_atm,viritot_atm,   &
     &                   ifhfvol,   &
     &                   nhfregion,hfzpos1,hfzpos2,   &
     &                   hftyp_atm,   &
     &                   molecom)

  use interface_interact, only: calijklim, calvirihf

  use md_global
  use mpi_global

  implicit none

!     subroutine to calculate torsion force & potential

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: ifhfvol       ! local volume-based or local surface-based

  integer,intent(in):: nhfregion       ! number of region to calculate heat flux
  real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                        ! z-position of region for heat flux

  integer,intent(in):: hftyp_atm(:)     ! atom- or mole-based heat flux cal.
                                        !   for each atom

  real(8),intent(in):: molecom(:,:)     ! center of mass of molecule

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! force calculated here

!   OUTPUT
  real(8),intent(out):: pot_torsim      ! improper torsion potential

  real(8),intent(out):: atm_viri_tors   ! virial(torsion potential) of each atom
  real(8),intent(out):: atm_virit_tors(:,:) ! virial tensor (torsion potential)

  real(8),intent(out):: pot_to_atm(:)   ! torsion potential of each atom

  real(8),intent(out):: viritot_atm(:,:,:,:)
                                        ! virial tensor of each atom (torsion)

! LOCAL:
  real(8):: delta_pot                   ! potential calculated
                                        !  in calijkl

  real(8):: rij(3)                      ! Rj-Ri
  real(8):: rkj(3)                      ! Rj-Rk
  real(8):: rkl(3)                      ! Rl-Rk
  real(8):: rik(3)                      ! Rk-Ri
  real(8):: rjl(3)                      ! Rl-Rj

  real(8):: fi(3), fj(3), fk(3), fl(3)  ! force to atoms (i,j,k,l)
                                        !  calculated here

  real(8):: phi                         ! torsion (in radian) defined by
                                        !  atoms i, j, k, l

  real(8):: pi                          ! pi

  integer:: n,l                         ! do loop index

  integer:: iatom, jatom, katom, latom  ! atom index

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

  real(8):: nbodycoeff             ! n-body coefficient for virial term of heatf

! ---- initialization

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

  nbodycoeff = 0.25d0       ! = 1/4

! --- GRAND LOOP TO CALCULATE TORSION FORCE ---


  pi = acos(-1.0d0)
  pot_torsim = 0.0d0

  looplast = ntorsim

! MPI loop
!      DO n = 1, ntorsim
  DO n = loopinit, looplast, loopstep

!    - calculate Rij, Rkj, Rkl, etc. -

     iatom   = itorsim(n)               ! atom index
     jatom   = jtorsim(n)               ! atom index
     katom   = ktorsim(n)               ! atom index
     latom   = ltorsim(n)               ! atom index

     rij(1:3)  = atmcor(1:3,jatom) - atmcor(1:3,iatom)
     rkj(1:3)  = atmcor(1:3,jatom) - atmcor(1:3,katom)
     rkl(1:3)  = atmcor(1:3,latom) - atmcor(1:3,katom)
     rik(1:3)  = atmcor(1:3,katom) - atmcor(1:3,iatom)
     rjl(1:3)  = atmcor(1:3,latom) - atmcor(1:3,jatom)

     rij(1:3)  = rij(1:3) - box(1:3)*anint(rij(1:3)*box_inv(1:3))
     rkj(1:3)  = rkj(1:3) - box(1:3)*anint(rkj(1:3)*box_inv(1:3))
     rkl(1:3)  = rkl(1:3) - box(1:3)*anint(rkl(1:3)*box_inv(1:3))
     rik(1:3)  = rik(1:3) - box(1:3)*anint(rik(1:3)*box_inv(1:3))
     rjl(1:3)  = rjl(1:3) - box(1:3)*anint(rjl(1:3)*box_inv(1:3))

!    - calculate torsion force and angle--

     call calijklim(rij,rkj,rkl,rik,rjl,   &
          &         para_barhigim(indextorsimtyp(n)),   &
          &         para_phsangim(indextorsimtyp(n)),   &
          &         delta_pot,fi,fj,fk,fl,phi)

     pot_torsim = pot_torsim + delta_pot

     pot_to_atm(iatom) = pot_to_atm(iatom) + nbodycoeff*delta_pot
     pot_to_atm(jatom) = pot_to_atm(jatom) + nbodycoeff*delta_pot
     pot_to_atm(katom) = pot_to_atm(katom) + nbodycoeff*delta_pot
     pot_to_atm(latom) = pot_to_atm(latom) + nbodycoeff*delta_pot

!    - add to force & potential & virial -

     force(1:3,iatom) = force(1:3,iatom) + fi(1:3)
     force(1:3,jatom) = force(1:3,jatom) + fj(1:3)
     force(1:3,katom) = force(1:3,katom) + fk(1:3)
     force(1:3,latom) = force(1:3,latom) + fl(1:3)

!    - for heat flux calculation

!    pair (i,j)
     call calvirihf(iatom,jatom,fi,fj,   &
          &                  nbodycoeff,   &
          &                  box,box_inv,   &
          &                  ifhfvol,   &
          &                  nhfregion,hfzpos1,hfzpos2,   &
          &                  hftyp_atm,   &
          &                  molecom,   &
          &                  viritot_atm)
!    pair (i,k)
     call calvirihf(iatom,katom,fi,fk,   &
          &                  nbodycoeff,   &
          &                  box,box_inv,   &
          &                  ifhfvol,   &
          &                  nhfregion,hfzpos1,hfzpos2,   &
          &                  hftyp_atm,   &
          &                  molecom,   &
          &                  viritot_atm)
!    pair (i,l)
     call calvirihf(iatom,latom,fi,fl,   &
          &                  nbodycoeff,   &
          &                  box,box_inv,   &
          &                  ifhfvol,   &
          &                  nhfregion,hfzpos1,hfzpos2,   &
          &                  hftyp_atm,   &
          &                  molecom,   &
          &                  viritot_atm)
!    pair (j,k)
     call calvirihf(jatom,katom,fj,fk,   &
          &                  nbodycoeff,   &
          &                  box,box_inv,   &
          &                  ifhfvol,   &
          &                  nhfregion,hfzpos1,hfzpos2,   &
          &                  hftyp_atm,   &
          &                  molecom,   &
          &                  viritot_atm)
!    pair (j,l)
     call calvirihf(jatom,latom,fj,fl,   &
          &                  nbodycoeff,   &
          &                  box,box_inv,   &
          &                  ifhfvol,   &
          &                  nhfregion,hfzpos1,hfzpos2,   &
          &                  hftyp_atm,   &
          &                  molecom,   &
          &                  viritot_atm)
!    pair (k,l)
     call calvirihf(katom,latom,fk,fl,   &
          &                  nbodycoeff,   &
          &                  box,box_inv,   &
          &                  ifhfvol,   &
          &                  nhfregion,hfzpos1,hfzpos2,   &
          &                  hftyp_atm,   &
          &                  molecom,   &
          &                  viritot_atm)

!    - for pressure calculation

     atm_viri_tors = atm_viri_tors + fi(1)*atmcor(1,iatom)   &
          &        + fi(2)*atmcor(2,iatom) + fi(3)*atmcor(3,iatom)
     atm_viri_tors = atm_viri_tors + fj(1)*atmcor(1,jatom)   &
          &        + fj(2)*atmcor(2,jatom) + fj(3)*atmcor(3,jatom)
     atm_viri_tors = atm_viri_tors + fk(1)*atmcor(1,katom)   &
          &        + fk(2)*atmcor(2,katom) + fk(3)*atmcor(3,katom)
     atm_viri_tors = atm_viri_tors + fl(1)*atmcor(1,latom)   &
          &        + fl(2)*atmcor(2,latom) + fl(3)*atmcor(3,latom)

     do l=1,3
        atm_virit_tors(1:3,l) = atm_virit_tors(1:3,l)   &
             &                + fi(1:3)*atmcor(l,iatom)
        atm_virit_tors(1:3,l) = atm_virit_tors(1:3,l)   &
             &                + fj(1:3)*atmcor(l,jatom)
        atm_virit_tors(1:3,l) = atm_virit_tors(1:3,l)   &
             &                + fk(1:3)*atmcor(l,katom)
        atm_virit_tors(1:3,l) = atm_virit_tors(1:3,l)   &
             &                + fl(1:3)*atmcor(l,latom)
     end do

!     phi = phi * 180.0d0 / pi
!     if (phi < 0.0d0) phi = phi + 360.0d0


  END DO

!     +     +     +     +     +     +     +

end subroutine caltorsimp_hf
