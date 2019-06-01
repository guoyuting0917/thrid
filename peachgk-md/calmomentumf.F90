!**********************************
!*  calmomentumf.f90 Ver.1.5      *
!*      for peachgk_md.f          *
!*            by N.Yamamoto       *
!*        modified by G.Kikugawa  *
!**********************************
! Time-stamp: <>

subroutine calmomentumf(npoly,nwater,nmatom, &
     &                  xcel,ycel,zcel, &
     &                  ifhfvol, &
     &                  nhfregion,hfzpos1,hfzpos2, &
     &                  hftyp_atm, &
     &                  molecom, &
     &                  viribot_long_atm,viribot_med_atm, &
     &                  viribot_short_atm, &
     &                  viriant_long_atm,viriant_med_atm, &
     &                  viriant_short_atm, &
     &                  viritot_long_atm,viritot_med_atm, &
     &                  viritot_short_atm, &
     &                  viri14t_long_atm,viri14t_med_atm, &
     &                  viri14t_short_atm, &
     &                  virielint_long_atm,virielint_med_atm, &
     &                  virielint_short_atm, &
     &                  viriljint_long_atm,viriljint_med_atm, &
     &                  viriljint_short_atm, &
     &                  virielct_long_atm,virielct_med_atm, &
     &                  virielct_short_atm, &
     &                  viriljt_long_atm,viriljt_med_atm, &
     &                  viriljt_short_atm, &
     &                  virimort_long_atm,virimort_med_atm, &
     &                  virimort_short_atm, &
     &                  virisht_long_atm,virisht_med_atm, &
     &                  virisht_short_atm, &
     &                  virirfht_long_atm,virirfht_med_atm, &
     &                  virirfht_short_atm, &
     &                  viridout_long_atm,viridout_med_atm, &
     &                  viridout_short_atm, &
     &                  viricstmnbt_long_atm,viricstmnbt_med_atm, &
     &                  viricstmnbt_short_atm, &
     &                  viricstmnbext_long_atm,viricstmnbext_med_atm, &
     &                  viricstmnbext_short_atm, &
     &                  mfkin_atm,mfkinex_atm, &
     &                  mfviribo_atm,mfvirian_atm,mfvirito_atm, &
     &                  mfviri14_atm, &
     &                  mfvirielin_atm,mfviriljin_atm, &
     &                  mfvirielc_atm,mfvirilj_atm, &
     &                  mfvirimor_atm,mfvirish_atm, &
     &                  mfvirirfh_atm,mfviridou_atm, &
     &                  mfviricstmnb_atm, &
     &                  mfviricstmnbex_atm, &
     &                  ncstmnbex, &
     &                  mftot_atm, &
     &                  atmcor_mf,molecom_mf, &
     &                  dt_cross, &
     &                  heatfinterval, &
     &                  current_step)

  use md_global

  implicit none

!
!     momentum flux localized on atomic position
!           PxzV = sigma(atom_i){mi.vx'vz}
!               + 1/2 * sigma(atom_i)simga(atom_j){(rijx*fijz)}
!
!!! In this version, momentum flux (pressure tensor) only without
!!! macroscopic flow can be calculated.
!
!!! Compile macro "_HF_ALL_DIR" is meaningless in this routine
!
! ARGUMENTS:
!     INPUT
  integer,intent(in):: npoly           ! number of polymer1
  integer,intent(in):: nwater          ! number of H2O molecules
  integer,intent(in):: nmatom          ! number of monatomic molecules

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

!---- variables for calculation of heat flux
  logical,intent(in):: ifhfvol        ! local volume-based or local surface-based

  integer,intent(in):: nhfregion
                                ! number of region to calculate momemtum flux
  real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                ! z-position of region for momentum flux

  integer,intent(in):: hftyp_atm(:)    ! atom- or mole-based momentum flux cal.
                                !   for each atom

  real(8),intent(in):: molecom(:,:)     ! center of mass of molecule

  real(8),intent(in):: dt_cross ! interval time to check if cross the surface
                                !   only for surface-based method
                                ! = dt_long_cal * outinterval

  integer,intent(in):: heatfinterval   ! interval of heatf output

  integer,intent(in):: current_step    ! current time step

  integer,intent(in):: ncstmnbex       ! number of extra custom NB output

!     virial tensor for each atom
  real(8),intent(in):: viribot_long_atm(:,:,:,:) ! long-range (bond)
  real(8),intent(in):: viribot_med_atm(:,:,:,:) ! medium-range (bond)
  real(8),intent(in):: viribot_short_atm(:,:,:,:) ! short-range (bond)
  real(8),intent(in):: viriant_long_atm(:,:,:,:) ! long-range (angle)
  real(8),intent(in):: viriant_med_atm(:,:,:,:) ! medium-range (angle)
  real(8),intent(in):: viriant_short_atm(:,:,:,:) ! short-range (angle)
  real(8),intent(in):: viritot_long_atm(:,:,:,:) ! long-range (torsion)
  real(8),intent(in):: viritot_med_atm(:,:,:,:) ! medium-range (torsion)
  real(8),intent(in):: viritot_short_atm(:,:,:,:) ! short-range (torsion)
  real(8),intent(in):: viri14t_long_atm(:,:,:,:) ! long-range (1-4)
  real(8),intent(in):: viri14t_med_atm(:,:,:,:) ! medium-range (1-4)
  real(8),intent(in):: viri14t_short_atm(:,:,:,:) ! short-range (1-4)
  real(8),intent(in):: virielint_long_atm(:,:,:,:) ! long-range (elc intra)
  real(8),intent(in):: virielint_med_atm(:,:,:,:) ! medium-range (elc intra)
  real(8),intent(in):: virielint_short_atm(:,:,:,:) ! short-range (elc intra)
  real(8),intent(in):: viriljint_long_atm(:,:,:,:) ! long-range (L-J intra)
  real(8),intent(in):: viriljint_med_atm(:,:,:,:) ! medium-range (L-J intra)
  real(8),intent(in):: viriljint_short_atm(:,:,:,:) ! short-range (L-J intra)
  real(8),intent(in):: virielct_long_atm(:,:,:,:) ! long-range (Coulomb)
  real(8),intent(in):: virielct_med_atm(:,:,:,:) ! medium-range (Coulomb)
  real(8),intent(in):: virielct_short_atm(:,:,:,:) ! short-range (Coulomb)
  real(8),intent(in):: viriljt_long_atm(:,:,:,:) ! long-range (L-J)
  real(8),intent(in):: viriljt_med_atm(:,:,:,:) ! medium-range (L-J)
  real(8),intent(in):: viriljt_short_atm(:,:,:,:) ! short-range (L-J)
  real(8),intent(in):: virimort_long_atm(:,:,:,:) ! long-range (Morse)
  real(8),intent(in):: virimort_med_atm(:,:,:,:) ! medium-range (Morse)
  real(8),intent(in):: virimort_short_atm(:,:,:,:) ! short-range (Morse)
  real(8),intent(in):: virisht_long_atm(:,:,:,:) ! long-range (SH)
  real(8),intent(in):: virisht_med_atm(:,:,:,:) ! medium-range (SH)
  real(8),intent(in):: virisht_short_atm(:,:,:,:) ! short-range (SH)
  real(8),intent(in):: virirfht_long_atm(:,:,:,:) ! long-range (RFH)
  real(8),intent(in):: virirfht_med_atm(:,:,:,:) ! medium-range (RFH)
  real(8),intent(in):: virirfht_short_atm(:,:,:,:) ! short-range (RFH)
  real(8),intent(in):: viridout_long_atm(:,:,:,:) ! long-range (DOU)
  real(8),intent(in):: viridout_med_atm(:,:,:,:) ! medium-range (DOU)
  real(8),intent(in):: viridout_short_atm(:,:,:,:) ! short-range (DOU)
  real(8),intent(in):: viricstmnbt_long_atm(:,:,:,:) ! long-range (custom NB)
  real(8),intent(in):: viricstmnbt_med_atm(:,:,:,:) ! medium-range (custom NB)
  real(8),intent(in):: viricstmnbt_short_atm(:,:,:,:) ! short-range (custom NB)
  real(8),intent(in):: viricstmnbext_long_atm(:,:,:,:,:)
                                                  ! long-range (custom NB)
  real(8),intent(in):: viricstmnbext_med_atm(:,:,:,:,:)
                                                  ! medium-range (custom NB)
  real(8),intent(in):: viricstmnbext_short_atm(:,:,:,:,:)
                                                  ! short-range (custom NB)

!     OUTPUT
!     for momentum flux sum
  real(8),intent(out):: mfkin_atm(:,:,:) ! momentum flux of kinetic part
  real(8),intent(out):: mfkinex_atm(:,:,:,:)
                              ! momentum flux of kinetic part (extra custom NB)

  real(8),intent(out):: mfviribo_atm(:,:,:) ! momentum flux of virial part (bond)
  real(8),intent(out):: mfvirian_atm(:,:,:) ! momentum flux of virial part (angle)
  real(8),intent(out):: mfvirito_atm(:,:,:) ! momentum flux of virial part (tors)
  real(8),intent(out):: mfviri14_atm(:,:,:) ! momentum flux of virial part (1-4)
  real(8),intent(out):: mfvirielin_atm(:,:,:)
                     ! momentum flux of virial part (elc intra)
  real(8),intent(out):: mfviriljin_atm(:,:,:)
                     ! momentum flux of virial part (L-J intra)

  real(8),intent(out):: mfvirielc_atm(:,:,:) ! momentum flux of virial part (Coulomb)
  real(8),intent(out):: mfvirilj_atm(:,:,:) ! momentum flux of virial part (L-J)
  real(8),intent(out):: mfvirimor_atm(:,:,:) ! momentum flux of virial part (Morse)
  real(8),intent(out):: mfvirish_atm(:,:,:) ! momentum flux of virial part (SH)
  real(8),intent(out):: mfvirirfh_atm(:,:,:) ! momentum flux of virial part (RFH)
  real(8),intent(out):: mfviridou_atm(:,:,:) ! momentum flux of virial part (DOU)
  real(8),intent(out):: mfviricstmnb_atm(:,:,:)
                                ! momentum flux of virial part (custom NB)
  real(8),intent(out):: mfviricstmnbex_atm(:,:,:,:)
                                ! momentum flux of virial part (extra custom NB)

  real(8),intent(out):: mftot_atm(:,:,:)   ! total momentum flux

  real(8),intent(out):: atmcor_mf(:,:)
                     ! old atmcor for momentum flux calculation
  real(8),intent(out):: molecom_mf(:,:)
                     ! old center of mass of molecule (momentum flux)

! LOCAL:
  real(8):: viribot_atm(3,3,natom,nhfregion) ! virial tensor (bond)
  real(8):: viriant_atm(3,3,natom,nhfregion) ! virial tensor (angl)
  real(8):: viritot_atm(3,3,natom,nhfregion) ! virial tensor (torsion)
  real(8):: viri14t_atm(3,3,natom,nhfregion) ! virial tensor (1-4)
  real(8):: virielint_atm(3,3,natom,nhfregion) ! virial tensor (elc intra)
  real(8):: viriljint_atm(3,3,natom,nhfregion) ! virial tensor (L-J intra)
  real(8):: virielct_atm(3,3,natom,nhfregion) ! virial tensor (Coulomb)
  real(8):: viriljt_atm(3,3,natom,nhfregion) ! virial tensor (L-J)
  real(8):: virimort_atm(3,3,natom,nhfregion) ! virial tensor (Morse)
  real(8):: virisht_atm(3,3,natom,nhfregion) ! virial tensor (SH)
  real(8):: virirfht_atm(3,3,natom,nhfregion) ! virial tensor (RFH)
  real(8):: viridout_atm(3,3,natom,nhfregion) ! virial tensor (DOU)
  real(8):: viricstmnbt_atm(3,3,natom,nhfregion) ! virial tensor (custom NB)
  real(8):: viricstmnbext_atm(3,3,natom,nhfregion,1:ncstmnbex)
                                                 ! virial tensor (custom NB)

  real(8):: box(3)
  real(8):: box_inv(3)
  real(8):: area_inv
  real(8):: vol_inv

  integer:: i,n
  integer:: im,i1,i2,ii
  integer:: nex

  integer:: iz

  real(8):: kin_part(3,3)

  real(8):: molemass         ! mass of molecule
  real(8):: inv_molemass     ! = 1/molemass

  real(8):: molevel(3,nmole)  ! velocity of center of mass

  real(8):: parcor(3)        ! center of mass or atom coordinate
  real(8):: parvel(3)        ! velocity of center of mass or atom velocity

  real(8):: parcor_sb(3)     ! COM or atom coordinate (for surface-based)

  integer:: moleindex

  real(8):: parnewold(3)     ! = parcor - parcor_sb
  real(8):: parnewold_abs    ! = |parcor - parcor_sb|
  real(8):: parnew_p1        ! = parnew(3) - hfzpos1
  real(8):: parold_p1        ! = parold(3) - hfzpos1
  real(8):: flagment_fact    ! = 1/dt * sign(vz)

#if defined(_MF_DEBUG) || defined(_MF_DEBUG_CALMOMENF)
  integer:: crosscount(nhfregion) ! counter for crossing the plane
#endif

!     +     +     +     +     +     +     +

#if defined(_MF_DEBUG) || defined(_MF_DEBUG_CALMOMENF)
  do i=1,nhfregion
     crosscount(i) = 0
  end do
#endif

!     --- some preparation ---

!     initialize variables

  ! When volume-based, clear transport terms always
  ! When surface-based, clear these ones after output heat flux data
  if (ifhfvol .or. current_step == 1) then
      mfkin_atm(1:3,1:3,1:nhfregion) = 0.0d0
      mfkinex_atm(1:3,1:3,1:nhfregion,1:ncstmnbex) = 0.0d0
  end if

  mfviribo_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfvirian_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfvirito_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfviri14_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfvirielin_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfviriljin_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfvirielc_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfvirilj_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfvirimor_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfvirish_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfvirirfh_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfviridou_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfviricstmnb_atm(1:3,1:3,1:nhfregion) = 0.0d0
  mfviricstmnbex_atm(1:3,1:3,1:nhfregion,1:ncstmnbex) = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)
  area_inv = box_inv(1) * box_inv(2)

  do iz = 1, nhfregion
     do i = 1, natom

#if defined(_HF_ALL_DIR)
        do n = 1, 3
#else
           n = 3
#endif
           viribot_atm(1:3,n,i,iz) = viribot_long_atm(1:3,n,i,iz) &
!     &                                 + viribot_med_atm(1:3,n,i,iz)
                &                + viribot_short_atm(1:3,n,i,iz)

           viriant_atm(1:3,n,i,iz) = viriant_long_atm(1:3,n,i,iz) &
!     &                                 + viriant_med_atm(1:3,n,i,iz)
                &                + viriant_short_atm(1:3,n,i,iz)

           viritot_atm(1:3,n,i,iz) = viritot_long_atm(1:3,n,i,iz) &
!     &                                 + viritot_med_atm(1:3,n,i,iz)
                &                + viritot_short_atm(1:3,n,i,iz)

           viri14t_atm(1:3,n,i,iz) = viri14t_long_atm(1:3,n,i,iz) &
!     &                                 + viri14t_med_atm(1:3,n,i,iz)
                &                + viri14t_short_atm(1:3,n,i,iz)

           virielint_atm(1:3,n,i,iz) = virielint_long_atm(1:3,n,i,iz) &
!     &                                 + virielint_med_atm(1:3,n,i,iz)
                &                  + virielint_short_atm(1:3,n,i,iz)

           viriljint_atm(1:3,n,i,iz) = viriljint_long_atm(1:3,n,i,iz) &
!     &                                 + viriljint_med_atm(1:3,n,i,iz)
                &                  + viriljint_short_atm(1:3,n,i,iz)

           virielct_atm(1:3,n,i,iz) = virielct_long_atm(1:3,n,i,iz) &
!     &                                  + virielct_med_atm(1:3,n,i,iz)
                &                 + virielct_short_atm(1:3,n,i,iz)

           viriljt_atm(1:3,n,i,iz) = viriljt_long_atm(1:3,n,i,iz) &
!     &                                 + viriljt_med_atm(1:3,n,i,iz)
                &                + viriljt_short_atm(1:3,n,i,iz)

           virimort_atm(1:3,n,i,iz) = virimort_long_atm(1:3,n,i,iz) &
!     &                                  + virimort_med_atm(1:3,n,i,iz)
                &                 + virimort_short_atm(1:3,n,i,iz)

           virisht_atm(1:3,n,i,iz) = virisht_long_atm(1:3,n,i,iz) &
!     &                                 + virisht_med_atm(1:3,n,i,iz)
                &                + virisht_short_atm(1:3,n,i,iz)

           virirfht_atm(1:3,n,i,iz) = virirfht_long_atm(1:3,n,i,iz) &
!     &                                  + virirfht_med_atm(1:3,n,i,iz)
                &                 + virirfht_short_atm(1:3,n,i,iz)

           viridout_atm(1:3,n,i,iz) = viridout_long_atm(1:3,n,i,iz) &
!     &                                  + viridout_med_atm(1:3,n,i,iz)
                &                 + viridout_short_atm(1:3,n,i,iz)

           viricstmnbt_atm(1:3,n,i,iz) = viricstmnbt_long_atm(1:3,n,i,iz) &
!     &                               + viricstmnbt_med_atm(1:3,n,i,iz)
                &                    + viricstmnbt_short_atm(1:3,n,i,iz)

           do nex = 1, ncstmnbex
              viricstmnbext_atm(1:3,n,i,iz,nex) = &
                   &                viricstmnbext_long_atm(1:3,n,i,iz,nex) &
                   ! &              + viricstmnbext_med_atm(1:3,n,i,iz,nex) &
                   &              + viricstmnbext_short_atm(1:3,n,i,iz,nex)
           end do

#if defined(_HF_ALL_DIR)
        end do
#else
#endif

     end do
  end do

!     ---- calculate velocity of center of mass

!     - loop over poly
  DO im = 1, npoly
     i1 = molept_index(im)
     i2 = molept_index(im+1) - 1

     molemass = 0.0d0
     molevel(1:3,im) = 0.0d0

     do ii = i1, i2         ! loop over atom(i)
        i = molept_list(ii)

        molemass = molemass + atmmass(i)

        molevel(1:3,im) = molevel(1:3,im) &
             &          + atmmass(i) * atmvel(1:3,i)

     end do

     inv_molemass = 1.0d0 / molemass

     molevel(1:3,im) = molevel(1:3,im) * inv_molemass

  END DO

!     - loop over water
  DO im = npoly+1, npoly+nwater
     i1 = molept_index(im)
     i2 = molept_index(im+1) - 1

     molemass = 0.0d0
     molevel(1:3,im) = 0.0d0

     do ii = i1, i2         ! loop over atom(i)
        i = molept_list(ii)

        molemass = molemass + atmmass(i)

        molevel(1:3,im) = molevel(1:3,im) &
             &          + atmmass(i) * atmvel(1:3,i)

     end do

     inv_molemass = 1.0d0 / molemass

     molevel(1:3,im) = molevel(1:3,im) * inv_molemass

  END DO

#if defined(_DO_NOT_EXECUTE_THIS)
!     - loop over monatom
  DO im = npoly+nwater+1, npoly+nwater+nmatom
     i1 = molept_index(im)
!         i2 = molept_index(im+1) - 1

!         do ii = i1, i2         ! loop over atom(i)
     ii = i1
     i = molept_list(ii)
     molevel(1:3,im) = atmvel(1:3,i)

!         end do

  END DO
#endif


!     -------- kinetic and potential part (transport by molecules) --------

  DO i = 1, natom

     if (hftyp_atm(i) == HFTYP_MOLE) then ! mole-base
        moleindex = irmolept_list(i)
        parcor(1:3) = molecom(1:3,moleindex)
        parvel(1:3) = molevel(1:3,moleindex)

!           for surface-based method only
        if (.not. ifhfvol) then
           parcor_sb(1:3) = molecom_mf(1:3,moleindex)
        end if

     else                   ! atom-base
        parcor(1:3) = atmcor(1:3,i)
        parvel(1:3) = atmvel(1:3,i)

!           for surface-based method only
        if (.not. ifhfvol) then
           parcor_sb(1:3) = atmcor_mf(1:3,i)
        end if

     end if

!        - P.B.C. (only z-dir)
     if (parcor(3) < 0.0d0) then
        parcor(3) = parcor(3) + box(3)
     else if (parcor(3) >= box(3)) then
        parcor(3) = parcor(3) - box(3)
     end if

     if (.not. ifhfvol) then
        if (parcor_sb(3) < 0.0d0) then
           parcor_sb(3) = parcor_sb(3) + box(3)
        else if (parcor_sb(3) >= box(3)) then
           parcor_sb(3) = parcor_sb(3) - box(3)
        end if
     end if

!        - calculate randam velosity
!     call calrandvel

!        - loop over hfregion
     do iz = 1, nhfregion

!-------- local volume-based method --------
        IF (ifhfvol) THEN

           if ((hfzpos1(iz) <= parcor(3)) .and.  &
                & (parcor(3) <= hfzpos2(iz))) then

!              calculate kinetic part
#if defined(_HF_ALL_DIR)
              do n = 1,3
#else
                 n = 3
#endif
                 kin_part(1,n) = atmmass(i) * parvel(1) * parvel(n)
                 kin_part(2,n) = atmmass(i) * parvel(2) * parvel(n)
                 kin_part(3,n) = atmmass(i) * parvel(3) * parvel(n)

#if defined(_HF_ALL_DIR)
              end do
#else
#endif

#if defined(_HF_ALL_DIR)
              mfkin_atm(1:3,1,iz) = mfkin_atm(1:3,1,iz) + kin_part(1:3,1)
              mfkin_atm(1:3,2,iz) = mfkin_atm(1:3,2,iz) + kin_part(1:3,2)
#endif
              mfkin_atm(1:3,3,iz) = mfkin_atm(1:3,3,iz) + kin_part(1:3,3)

!!! if you want to output some extra data about kinetic part,
!!!  write your procedures here.
! #if defined(_HF_ALL_DIR) || defined(_HF_BULK)
!               nex = 1
!               mfkinex_atm(1:3,1,iz,nex) = mfkinex_atm(1:3,1,iz,nex) &
!                    &                    +  kin_part(1:3,1)
!               mfkinex_atm(1:3,2,iz,nex) = mfkinex_atm(1:3,2,iz,nex) &
!                    &                    +  kin_part(1:3,2)
! #endif
!               mfkinex_atm(1:3,3,iz,nex) = mfkinex_atm(1:3,3,iz,nex) &
!                    &                    +  kin_part(1:3,3)
!!!

           end if

!-------- local surface-based method --------
        ELSE

           !!! old information is not stored
           if (current_step == 1) exit

           parnewold(3) = parcor(3) - parcor_sb(3)
           parnewold_abs = abs(parnewold(3))

           if (parnewold_abs * box_inv(3) >= 0.5d0) exit
                                ! parnewold exceed half cell

           parnew_p1 = parcor(3) - hfzpos1(iz)
           parold_p1 = parcor_sb(3) - hfzpos1(iz)

           if (parnew_p1*parold_p1 < 0.0d0) then

#if defined(_MF_DEBUG) || defined(_MF_DEBUG_CALMOMENF)
              crosscount(iz) = crosscount(iz) + 1
#endif
              if (parnew_p1 > 0.0d0) then ! new parcor is positive
                 flagment_fact = 1.0d0
              else             ! old parcor is positive from the surface
                 flagment_fact = -1.0d0
              endif

              flagment_fact = flagment_fact / dt_cross

!              calculate kinetic part
              kin_part(1,3) = atmmass(i) * parvel(1) * flagment_fact
              kin_part(2,3) = atmmass(i) * parvel(2) * flagment_fact
              kin_part(3,3) = atmmass(i) * parvel(3) * flagment_fact

              mfkin_atm(1:3,3,iz) = mfkin_atm(1:3,3,iz) + kin_part(1:3,3)

!!! if you want to output some extra data about kinetic part,
!!!  write your procedures here.
! #if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              ! nex = 1
              ! mfkinex_atm(1:3,3,iz,nex) = mfkinex_atm(1:3,3,iz,nex) &
              !      &                    +  kin_part(1:3,3)
!!!

           end if

        END IF                 ! local surface-based method end

     end do

  END DO

!---- store current coordinate information (only for surface-based)
#if defined(_MF_DEBUG) || (_MF_DEBUG_CALHEATF)
  write(6,*) '*** momentumf debug info (in calmomentumf) ***'

  write(6,*) '* molecom'
  do i = 1, nmole
     write(6,*) i, (molecom(n,i),n=1,3)
  end do

  write(6,*)
  write(6,*) '* atmcor'
  do i = 1, natom
     write(6,*) i, (atmcor(n,i),n=1,3)
  end do

  write(6,*)
  write(6,*) '* number of crossing'
  do i = 1, nhfregion
     write(6,*) i, crosscount(i)
  end do
#endif
  if (.not. ifhfvol) then
     molecom_mf(1:3,1:nmole) = molecom(1:3,1:nmole)

     atmcor_mf(1:3,1:natom) = atmcor(1:3,1:natom)

  end if

  !!! When surface-based, exit subroutine other than output MD step here
  if (.not. ifhfvol .and. (mod(current_step,heatfinterval) /= 1) .and. &
  &   (heatfinterval /= 1)) then
      return
  end if

!     -------- collision part (transport by interaction) --------
  do iz = 1, nhfregion

     do i = 1, natom
#if defined(_HF_ALL_DIR)
        mfviribo_atm(1:3,1,iz) =  mfviribo_atm(1:3,1,iz) &
             &                  + viribot_atm(1:3,1,i,iz)
        mfviribo_atm(1:3,2,iz) =  mfviribo_atm(1:3,2,iz) &
             &                  + viribot_atm(1:3,2,i,iz)
#endif
        mfviribo_atm(1:3,3,iz) =  mfviribo_atm(1:3,3,iz) &
             &                  + viribot_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfvirian_atm(1:3,1,iz) =  mfvirian_atm(1:3,1,iz) &
             &                  + viriant_atm(1:3,1,i,iz)
        mfvirian_atm(1:3,2,iz) =  mfvirian_atm(1:3,2,iz) &
             &                  + viriant_atm(1:3,2,i,iz)
#endif
        mfvirian_atm(1:3,3,iz) =  mfvirian_atm(1:3,3,iz) &
             &                  + viriant_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfvirito_atm(1:3,1,iz) =  mfvirito_atm(1:3,1,iz) &
             &                  + viritot_atm(1:3,1,i,iz)
        mfvirito_atm(1:3,2,iz) =  mfvirito_atm(1:3,2,iz) &
             &                  + viritot_atm(1:3,2,i,iz)
#endif
        mfvirito_atm(1:3,3,iz) =  mfvirito_atm(1:3,3,iz) &
             &                  + viritot_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfviri14_atm(1:3,1,iz) =  mfviri14_atm(1:3,1,iz) &
             &                  + viri14t_atm(1:3,1,i,iz)
        mfviri14_atm(1:3,2,iz) =  mfviri14_atm(1:3,2,iz) &
             &                  + viri14t_atm(1:3,2,i,iz)
#endif
        mfviri14_atm(1:3,3,iz) =  mfviri14_atm(1:3,3,iz) &
             &                  + viri14t_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfvirielin_atm(1:3,1,iz) =  mfvirielin_atm(1:3,1,iz) &
             &                    + virielint_atm(1:3,1,i,iz)
        mfvirielin_atm(1:3,2,iz) =  mfvirielin_atm(1:3,2,iz) &
             &                    + virielint_atm(1:3,2,i,iz)
#endif
        mfvirielin_atm(1:3,3,iz) =  mfvirielin_atm(1:3,3,iz) &
             &                    + virielint_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfviriljin_atm(1:3,1,iz) =  mfviriljin_atm(1:3,1,iz) &
             &                    + viriljint_atm(1:3,1,i,iz)
        mfviriljin_atm(1:3,2,iz) =  mfviriljin_atm(1:3,2,iz) &
             &                    + viriljint_atm(1:3,2,i,iz)
#endif
        mfviriljin_atm(1:3,3,iz) =  mfviriljin_atm(1:3,3,iz) &
             &                    + viriljint_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfvirielc_atm(1:3,1,iz) =  mfvirielc_atm(1:3,1,iz) &
             &                   + virielct_atm(1:3,1,i,iz)
        mfvirielc_atm(1:3,2,iz) =  mfvirielc_atm(1:3,2,iz) &
             &                   + virielct_atm(1:3,2,i,iz)
#endif
        mfvirielc_atm(1:3,3,iz) =  mfvirielc_atm(1:3,3,iz) &
             &                   + virielct_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfvirilj_atm(1:3,1,iz) =  mfvirilj_atm(1:3,1,iz) &
             &                  + viriljt_atm(1:3,1,i,iz)
        mfvirilj_atm(1:3,2,iz) =  mfvirilj_atm(1:3,2,iz) &
             &                  + viriljt_atm(1:3,2,i,iz)
#endif
        mfvirilj_atm(1:3,3,iz) =  mfvirilj_atm(1:3,3,iz) &
             &                  + viriljt_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfvirimor_atm(1:3,1,iz) =  mfvirimor_atm(1:3,1,iz) &
             &                   + virimort_atm(1:3,1,i,iz)
        mfvirimor_atm(1:3,2,iz) =  mfvirimor_atm(1:3,2,iz) &
             &                   + virimort_atm(1:3,2,i,iz)
#endif
        mfvirimor_atm(1:3,3,iz) =  mfvirimor_atm(1:3,3,iz) &
             &                   + virimort_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfvirish_atm(1:3,1,iz) =  mfvirish_atm(1:3,1,iz) &
             &                  + virisht_atm(1:3,1,i,iz)
        mfvirish_atm(1:3,2,iz) =  mfvirish_atm(1:3,2,iz) &
             &                  + virisht_atm(1:3,2,i,iz)
#endif
        mfvirish_atm(1:3,3,iz) =  mfvirish_atm(1:3,3,iz) &
             &                  + virisht_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfvirirfh_atm(1:3,1,iz) =  mfvirirfh_atm(1:3,1,iz) &
             &                   + virirfht_atm(1:3,1,i,iz)
        mfvirirfh_atm(1:3,2,iz) =  mfvirirfh_atm(1:3,2,iz) &
             &                   + virirfht_atm(1:3,2,i,iz)
#endif
        mfvirirfh_atm(1:3,3,iz) =  mfvirirfh_atm(1:3,3,iz) &
             &                   + virirfht_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfviridou_atm(1:3,1,iz) =  mfviridou_atm(1:3,1,iz) &
             &                   + viridout_atm(1:3,1,i,iz)
        mfviridou_atm(1:3,2,iz) =  mfviridou_atm(1:3,2,iz) &
             &                   + viridout_atm(1:3,2,i,iz)
#endif
        mfviridou_atm(1:3,3,iz) =  mfviridou_atm(1:3,3,iz) &
             &                   + viridout_atm(1:3,3,i,iz)

#if defined(_HF_ALL_DIR)
        mfviricstmnb_atm(1:3,1,iz) =  mfviricstmnb_atm(1:3,1,iz) &
             &                      + viricstmnbt_atm(1:3,1,i,iz)
        mfviricstmnb_atm(1:3,2,iz) =  mfviricstmnb_atm(1:3,2,iz) &
             &                      + viricstmnbt_atm(1:3,2,i,iz)
#endif
        mfviricstmnb_atm(1:3,3,iz) =  mfviricstmnb_atm(1:3,3,iz) &
             &                      + viricstmnbt_atm(1:3,3,i,iz)

        do nex = 1, ncstmnbex
#if defined(_HF_ALL_DIR)
           mfviricstmnbex_atm(1:3,1,iz,nex) = mfviricstmnbex_atm(1:3,1,iz,nex) &
                &                         + viricstmnbext_atm(1:3,1,i,iz,nex)
           mfviricstmnbex_atm(1:3,2,iz,nex) = mfviricstmnbex_atm(1:3,2,iz,nex) &
                &                         + viricstmnbext_atm(1:3,2,i,iz,nex)
#endif
           mfviricstmnbex_atm(1:3,3,iz,nex) = mfviricstmnbex_atm(1:3,3,iz,nex) &
                &                         + viricstmnbext_atm(1:3,3,i,iz,nex)
        end do

     end do
  end do

! --- sum up momentum flux ---

  do iz = 1, nhfregion

     if (ifhfvol) then      ! volume-based (= volume)
        vol_inv = area_inv / (hfzpos2(iz) - hfzpos1(iz))
     else                   ! surface-based (= surface area)
        vol_inv = area_inv
     end if

     mfkin_atm(1:3,1:3,iz) =  mfkin_atm(1:3,1:3,iz) * vol_inv

     do nex = 1, ncstmnbex
        mfkinex_atm(1:3,1:3,iz,nex) = mfkinex_atm(1:3,1:3,iz,nex) * vol_inv
     end do

     mfviribo_atm(1:3,1:3,iz) =  mfviribo_atm(1:3,1:3,iz) * vol_inv

     mfvirian_atm(1:3,1:3,iz) =  mfvirian_atm(1:3,1:3,iz) * vol_inv

     mfvirito_atm(1:3,1:3,iz) =  mfvirito_atm(1:3,1:3,iz) * vol_inv

     mfviri14_atm(1:3,1:3,iz) =  mfviri14_atm(1:3,1:3,iz) * vol_inv

     mfvirielin_atm(1:3,1:3,iz) =  mfvirielin_atm(1:3,1:3,iz) * vol_inv

     mfviriljin_atm(1:3,1:3,iz) =  mfviriljin_atm(1:3,1:3,iz) * vol_inv

     mfvirielc_atm(1:3,1:3,iz) =  mfvirielc_atm(1:3,1:3,iz) * vol_inv

     mfvirilj_atm(1:3,1:3,iz) =  mfvirilj_atm(1:3,1:3,iz) * vol_inv

     mfvirimor_atm(1:3,1:3,iz) =  mfvirimor_atm(1:3,1:3,iz) * vol_inv

     mfvirish_atm(1:3,1:3,iz) =  mfvirish_atm(1:3,1:3,iz) * vol_inv

     mfvirirfh_atm(1:3,1:3,iz) =  mfvirirfh_atm(1:3,1:3,iz) * vol_inv

     mfviridou_atm(1:3,1:3,iz) =  mfviridou_atm(1:3,1:3,iz) * vol_inv

     mfviricstmnb_atm(1:3,1:3,iz) =  mfviricstmnb_atm(1:3,1:3,iz) * vol_inv

     do nex = 1, ncstmnbex
        mfviricstmnbex_atm(1:3,1:3,iz,nex) = &
             &                 mfviricstmnbex_atm(1:3,1:3,iz,nex) * vol_inv
     end do

     mftot_atm(1:3,1:3,iz) = mfkin_atm(1:3,1:3,iz) &
          &                   + mfviribo_atm(1:3,1:3,iz) &
          &                   + mfvirian_atm(1:3,1:3,iz) &
          &                   + mfvirito_atm(1:3,1:3,iz) &
          &                   + mfviri14_atm(1:3,1:3,iz) &
          &                   + mfvirielin_atm(1:3,1:3,iz) &
          &                   + mfviriljin_atm(1:3,1:3,iz) &
          &                   + mfvirielc_atm(1:3,1:3,iz) &
          &                   + mfvirilj_atm(1:3,1:3,iz) &
          &                   + mfvirimor_atm(1:3,1:3,iz) &
          &                   + mfvirish_atm(1:3,1:3,iz) &
          &                   + mfvirirfh_atm(1:3,1:3,iz) &
          &                   + mfviridou_atm(1:3,1:3,iz) &
          &                   + mfviricstmnb_atm(1:3,1:3,iz)

#if defined(_CSTMNB_V2_ADD_ALL)
     !- adding cstmnb extra interaction
     do nex = 1, ncstmnbex
        mftot_atm(1:3,1:3,iz) = mftot_atm(1:3,1:3,iz) &
             &                + mfviricstmnbex_atm(1:3,1:3,iz,nex)
     end do
#endif

  end do

!     +     +     +     +     +     +     +

end subroutine calmomentumf
