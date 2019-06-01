!*****************************
!*  calheatf.f90 Ver.2.3     *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>

subroutine calheatf(npoly,nwater,nmatom, &
     &              xcel,ycel,zcel, &
     &              ifhfvol, &
     &              nhfregion,hfzpos1,hfzpos2, &
     &              hftyp_atm, &
     &              molecom, &
     &              pot_bo_long_atm,pot_bo_med_atm, &
     &              pot_bo_short_atm, &
     &              pot_an_long_atm,pot_an_med_atm, &
     &              pot_an_short_atm, &
     &              pot_to_long_atm,pot_to_med_atm, &
     &              pot_to_short_atm, &
     &              pot_14_long_atm,pot_14_med_atm, &
     &              pot_14_short_atm, &
     &              pot_elin_long_atm,pot_elin_med_atm, &
     &              pot_elin_short_atm, &
     &              pot_ljin_long_atm,pot_ljin_med_atm, &
     &              pot_ljin_short_atm, &
     &              pot_elc_long_atm,pot_elc_med_atm, &
     &              pot_elc_short_atm, &
     &              pot_vdw_long_atm,pot_vdw_med_atm, &
     &              pot_vdw_short_atm, &
     &              pot_mor_long_atm,pot_mor_med_atm, &
     &              pot_mor_short_atm, &
     &              pot_sh_long_atm,pot_sh_med_atm, &
     &              pot_sh_short_atm, &
     &              pot_rfh_long_atm,pot_rfh_med_atm, &
     &              pot_rfh_short_atm, &
     &              pot_dou_long_atm,pot_dou_med_atm, &
     &              pot_dou_short_atm, &
     &              pot_cstmnb_long_atm,pot_cstmnb_med_atm, &
     &              pot_cstmnb_short_atm, &
     &              pot_cstmnbex_long_atm,pot_cstmnbex_med_atm, &
     &              pot_cstmnbex_short_atm, &
     &              viribot_long_atm,viribot_med_atm, &
     &              viribot_short_atm, &
     &              viriant_long_atm,viriant_med_atm, &
     &              viriant_short_atm, &
     &              viritot_long_atm,viritot_med_atm, &
     &              viritot_short_atm, &
     &              viri14t_long_atm,viri14t_med_atm, &
     &              viri14t_short_atm, &
     &              virielint_long_atm,virielint_med_atm, &
     &              virielint_short_atm, &
     &              viriljint_long_atm,viriljint_med_atm, &
     &              viriljint_short_atm, &
     &              virielct_long_atm,virielct_med_atm, &
     &              virielct_short_atm, &
     &              viriljt_long_atm,viriljt_med_atm, &
     &              viriljt_short_atm, &
     &              virimort_long_atm,virimort_med_atm, &
     &              virimort_short_atm, &
     &              virisht_long_atm,virisht_med_atm, &
     &              virisht_short_atm, &
     &              virirfht_long_atm,virirfht_med_atm, &
     &              virirfht_short_atm, &
     &              viridout_long_atm,viridout_med_atm, &
     &              viridout_short_atm, &
     &              viricstmnbt_long_atm,viricstmnbt_med_atm, &
     &              viricstmnbt_short_atm, &
     &              viricstmnbext_long_atm,viricstmnbext_med_atm, &
     &              viricstmnbext_short_atm, &
     &              hfkin_atm,hfkinex_atm, &
     &              hfpotbo_atm,hfpotan_atm,hfpotto_atm, &
     &              hfpot14_atm, &
     &              hfpotelin_atm,hfpotljin_atm, &
     &              hfpotelc_atm,hfpotlj_atm, &
     &              hfpotmor_atm,hfpotsh_atm, &
     &              hfpotrfh_atm,hfpotdou_atm, &
     &              hfpotcstmnb_atm, &
     &              hfpotcstmnbex_atm, &
     &              hfviribo_atm,hfvirian_atm,hfvirito_atm, &
     &              hfviri14_atm, &
     &              hfvirielin_atm,hfviriljin_atm, &
     &              hfvirielc_atm,hfvirilj_atm, &
     &              hfvirimor_atm,hfvirish_atm, &
     &              hfvirirfh_atm,hfviridou_atm, &
     &              hfviricstmnb_atm, &
     &              hfviricstmnbex_atm, &
     &              ncstmnbex, &
     &              hftot_atm, &
     &              atmcor_hf,molecom_hf, &
     &              dt_cross, &
     &              heatfinterval, &
     &              current_step)

  use md_global

  implicit none

!
!     heat flux localized on atomic position
!           JqV = sigma(atom_i){ei.vi}
!               + 1/2 * sigma(atom_i)simga(atom_j){(rij*fij) . vi}
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

  integer,intent(in):: nhfregion       ! number of region to calculate heat flux
  real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                ! z-position of region for heat flux

  integer,intent(in):: hftyp_atm(:)    ! atom- or mole-based heat flux cal.
                                !   for each atom

  real(8),intent(in):: molecom(:,:)     ! center of mass of molecule

  real(8),intent(in):: dt_cross     ! interval time to check if cross the surface
                                !   only for surface-based method
                                ! = dt_long_cal * outinterval

  integer,intent(in):: heatfinterval   ! interval of heatf output

  integer,intent(in):: current_step    ! current time step

  integer,intent(in):: ncstmnbex       ! number of extra custom NB output

!     potential for each atom
  real(8),intent(in):: pot_bo_long_atm(:) ! bond potential of each atom
  real(8),intent(in):: pot_bo_med_atm(:) ! bond potential of each atom
  real(8),intent(in):: pot_bo_short_atm(:) ! bond potential of each atom
  real(8),intent(in):: pot_an_long_atm(:) ! angle potential of each atom
  real(8),intent(in):: pot_an_med_atm(:) ! angle potential of each atom
  real(8),intent(in):: pot_an_short_atm(:) ! angle potential of each atom
  real(8),intent(in):: pot_to_long_atm(:) ! torsion potential of each atom
  real(8),intent(in):: pot_to_med_atm(:) ! torsion potential of each atom
  real(8),intent(in):: pot_to_short_atm(:) ! torsion potential of each atom
  real(8),intent(in):: pot_14_long_atm(:) ! 1-4 int potential of each atom
  real(8),intent(in):: pot_14_med_atm(:) ! 1-4 int potential of each atom
  real(8),intent(in):: pot_14_short_atm(:) ! 1-4 int potential of each atom
  real(8),intent(in):: pot_elin_long_atm(:)
                    ! elc (intra) int potential of each atom
  real(8),intent(in):: pot_elin_med_atm(:)
                    ! elc (intra) int potential of each atom
  real(8),intent(in):: pot_elin_short_atm(:)
                    ! elc (intra) int potential of each atom
  real(8),intent(in):: pot_ljin_long_atm(:)
                    ! VDW (intra) int potential of each atom
  real(8),intent(in):: pot_ljin_med_atm(:)
                    ! VDW (intra) int potential of each atom
  real(8),intent(in):: pot_ljin_short_atm(:)
                    ! VDW (intra) int potential of each atom

  real(8),intent(in):: pot_elc_long_atm(:) ! Coulomb potential of each atom
  real(8),intent(in):: pot_elc_med_atm(:) ! Coulomb potential of each atom
  real(8),intent(in):: pot_elc_short_atm(:) ! Coulomb potential of each atom
  real(8),intent(in):: pot_vdw_long_atm(:) ! VDW potential of each atom
  real(8),intent(in):: pot_vdw_med_atm(:) ! VDW potential of each atom
  real(8),intent(in):: pot_vdw_short_atm(:) ! VDW potential of each atom
  real(8),intent(in):: pot_mor_long_atm(:) ! Morse potential of each atom
  real(8),intent(in):: pot_mor_med_atm(:)   ! Morse potential of each atom
  real(8),intent(in):: pot_mor_short_atm(:) ! Morse potential of each atom
  real(8),intent(in):: pot_sh_long_atm(:) ! SH potential of each atom
  real(8),intent(in):: pot_sh_med_atm(:) ! SH potential of each atom
  real(8),intent(in):: pot_sh_short_atm(:)  ! SH potential of each atom
  real(8),intent(in):: pot_rfh_long_atm(:) ! RFH potential of each atom
  real(8),intent(in):: pot_rfh_med_atm(:) ! RFH potential of each atom
  real(8),intent(in):: pot_rfh_short_atm(:)  ! RFH potential of each atom
  real(8),intent(in):: pot_dou_long_atm(:)  ! DOU potential of each atom
  real(8),intent(in):: pot_dou_med_atm(:)   ! DOU potential of each atom
  real(8),intent(in):: pot_dou_short_atm(:) ! DOU potential of each atom
  real(8),intent(in):: pot_cstmnb_long_atm(:) ! custom NB potential of each atom
  real(8),intent(in):: pot_cstmnb_med_atm(:)  ! custom NB potential of each atom
  real(8),intent(in):: pot_cstmnb_short_atm(:) ! custom NB potential of each atom
  real(8),intent(in):: pot_cstmnbex_long_atm(:,:)
                                        ! custom NB potential of each atom
  real(8),intent(in):: pot_cstmnbex_med_atm(:,:)
                                        ! custom NB potential of each atom
  real(8),intent(in):: pot_cstmnbex_short_atm(:,:)
                                        ! custom NB potential of each atom

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
!     for heat flux sum
  real(8),intent(out):: hfkin_atm(:,:) ! heat flux of kinetic part
  real(8),intent(out):: hfkinex_atm(:,:,:)
                                ! heat flux of kinetic part (extra custom NB)

  real(8),intent(out):: hfpotbo_atm(:,:) ! heat flux of potential part (bond)
  real(8),intent(out):: hfpotan_atm(:,:) ! heat flux of potential part (angle)
  real(8),intent(out):: hfpotto_atm(:,:) ! heat flux of potential part (tors)
  real(8),intent(out):: hfpot14_atm(:,:) ! heat flux of potential part (1-4)
  real(8),intent(out):: hfpotelin_atm(:,:)
                     ! heat flux of potential part (elc intra)
  real(8),intent(out):: hfpotljin_atm(:,:)
                     ! heat flux of potential part (L-J intra)

  real(8),intent(out):: hfpotelc_atm(:,:) ! heat flux of potential part (Coulomb)
  real(8),intent(out):: hfpotlj_atm(:,:) ! heat flux of potential part (L-J)
  real(8),intent(out):: hfpotmor_atm(:,:) ! heat flux of potential part (Morse)
  real(8),intent(out):: hfpotsh_atm(:,:) ! heat flux of potential part (SH)
  real(8),intent(out):: hfpotrfh_atm(:,:) ! heat flux of potential part (RFH)
  real(8),intent(out):: hfpotdou_atm(:,:) ! heat flux of potential part (DOU)
  real(8),intent(out):: hfpotcstmnb_atm(:,:)
                                ! heat flux of potential part (custom NB)
  real(8),intent(out):: hfpotcstmnbex_atm(:,:,:)
                                ! heat flux of potential part (extra custom NB)

  real(8),intent(out):: hfviribo_atm(:,:)   ! heat flux of virial part (bond)
  real(8),intent(out):: hfvirian_atm(:,:)   ! heat flux of virial part (angle)
  real(8),intent(out):: hfvirito_atm(:,:)   ! heat flux of virial part (tors)
  real(8),intent(out):: hfviri14_atm(:,:)   ! heat flux of virial part (1-4)
  real(8),intent(out):: hfvirielin_atm(:,:)
                     ! heat flux of virial part (elc intra)
  real(8),intent(out):: hfviriljin_atm(:,:)
                     ! heat flux of virial part (L-J intra)

  real(8),intent(out):: hfvirielc_atm(:,:) ! heat flux of virial part (Coulomb)
  real(8),intent(out):: hfvirilj_atm(:,:) ! heat flux of virial part (L-J)
  real(8),intent(out):: hfvirimor_atm(:,:) ! heat flux of virial part (Morse)
  real(8),intent(out):: hfvirish_atm(:,:) ! heat flux of virial part (SH)
  real(8),intent(out):: hfvirirfh_atm(:,:) ! heat flux of virial part (RFH)
  real(8),intent(out):: hfviridou_atm(:,:) ! heat flux of virial part (DOU)
  real(8),intent(out):: hfviricstmnb_atm(:,:)
                                        ! heat flux of virial part (custom NB)
  real(8),intent(out):: hfviricstmnbex_atm(:,:,:)
                                  ! heat flux of virial part (extra custom NB)

  real(8),intent(out):: hftot_atm(:,:)   ! total heat flux

  real(8),intent(out):: atmcor_hf(:,:)   ! old atmcor for heat flux calculation
  real(8),intent(out):: molecom_hf(:,:)
                     ! old center of mass of molecule (heat flux)

! LOCAL:
  real(8):: pot_bo_atm(natom) ! bond potential of each atom
  real(8):: pot_an_atm(natom) ! angle potential of each atom
  real(8):: pot_to_atm(natom) ! torsion potential of each atom
  real(8):: pot_14_atm(natom) ! 1-4 potential of each atom
  real(8):: pot_elin_atm(natom) ! elc (intra)  potential of each atom
  real(8):: pot_ljin_atm(natom) ! VDW (intra)  potential of each atom
  real(8):: pot_elc_atm(natom) ! Coulomb potential of each atom
  real(8):: pot_vdw_atm(natom) ! VDW potential of each atom
  real(8):: pot_mor_atm(natom) ! Morse potential of each atom
  real(8):: pot_sh_atm(natom) ! SH potential of each atom
  real(8):: pot_rfh_atm(natom) ! RFH potential of each atom
  real(8):: pot_dou_atm(natom) ! DOU potential of each atom
  real(8):: pot_cstmnb_atm(natom) ! custom NB potential of each atom
  real(8):: pot_cstmnbex_atm(natom,1:ncstmnbex)
                                  ! custom NB extra potential of each atom

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
                                             ! extra virial tensor (custom NB)

  real(8):: box(3)
  real(8):: box_inv(3)
  real(8):: area_inv
  real(8):: vol_inv

  integer:: i,n
  integer:: im,i1,i2,ii
  integer:: nex

  integer:: iz

  real(8):: kin_ene

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

  real(8):: atmvel_tmp(3,natom) ! temporary atom velocity without streaming vel.

#if defined(_HF_DEBUG) || defined(_HF_DEBUG_CALHEATF)
  integer:: crosscount(nhfregion) ! counter for crossing the plane
#endif

!     +     +     +     +     +     +     +

#if defined(_HF_DEBUG) || defined(_HF_DEBUG_CALHEATF)
  do i=1,nhfregion
     crosscount(i) = 0
  end do
#endif

!     --- some preparation ---

  atmvel_tmp(1:3,1:natom) = atmvel(1:3,1:natom)
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) - atmvel_strm(1:3,1:natom)

!     initialize variables

  ! When volume-based, clear transport terms always
  ! When surface-based, clear these ones after output heat flux data
  if (ifhfvol .or. current_step == 1) then
      hfkin_atm(1:3,1:nhfregion) = 0.0d0
      hfkinex_atm(1:3,1:nhfregion,1:ncstmnbex) = 0.0d0

      hfpotbo_atm(1:3,1:nhfregion) = 0.0d0
      hfpotan_atm(1:3,1:nhfregion) = 0.0d0
      hfpotto_atm(1:3,1:nhfregion) = 0.0d0
      hfpot14_atm(1:3,1:nhfregion) = 0.0d0
      hfpotelin_atm(1:3,1:nhfregion) = 0.0d0
      hfpotljin_atm(1:3,1:nhfregion) = 0.0d0
      hfpotelc_atm(1:3,1:nhfregion) = 0.0d0
      hfpotlj_atm(1:3,1:nhfregion) = 0.0d0
      hfpotmor_atm(1:3,1:nhfregion) = 0.0d0
      hfpotsh_atm(1:3,1:nhfregion) = 0.0d0
      hfpotrfh_atm(1:3,1:nhfregion) = 0.0d0
      hfpotdou_atm(1:3,1:nhfregion) = 0.0d0
      hfpotcstmnb_atm(1:3,1:nhfregion) = 0.0d0
      hfpotcstmnbex_atm(1:3,1:nhfregion,1:ncstmnbex) = 0.0d0
  end if

  hfviribo_atm(1:3,1:nhfregion) = 0.0d0
  hfvirian_atm(1:3,1:nhfregion) = 0.0d0
  hfvirito_atm(1:3,1:nhfregion) = 0.0d0
  hfviri14_atm(1:3,1:nhfregion) = 0.0d0
  hfvirielin_atm(1:3,1:nhfregion) = 0.0d0
  hfviriljin_atm(1:3,1:nhfregion) = 0.0d0
  hfvirielc_atm(1:3,1:nhfregion) = 0.0d0
  hfvirilj_atm(1:3,1:nhfregion) = 0.0d0
  hfvirimor_atm(1:3,1:nhfregion) = 0.0d0
  hfvirish_atm(1:3,1:nhfregion) = 0.0d0
  hfvirirfh_atm(1:3,1:nhfregion) = 0.0d0
  hfviridou_atm(1:3,1:nhfregion) = 0.0d0
  hfviricstmnb_atm(1:3,1:nhfregion) = 0.0d0
  hfviricstmnbex_atm(1:3,1:nhfregion,1:ncstmnbex) = 0.0d0

  hfpotlj_atm_ar11(1:3,1:nhfregion) = 0.0d0                       !!!!!!!!!!!!!!!!!ar1-ar2 heat flux (pot + viri) define
  hfpotlj_atm_ar12(1:3,1:nhfregion) = 0.0d0
  hfpotlj_atm_ar22(1:3,1:nhfregion) = 0.0d0 

  hfvirilj_atm_ar11(1:3,1:nhfregion) = 0.0d0
  hfvirilj_atm_ar12(1:3,1:nhfregion) = 0.0d0
  hfvirilj_atm_ar22(1:3,1:nhfregion) = 0.0d0
  hfvirilj_atm_ar1pt(1:3,1:nhfregion) = 0.0d0
  hfvirilj_atm_ar2pt(1:3,1:nhfregion) = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)
  area_inv = box_inv(1) * box_inv(2)

!     sum up each MTS contribution
  pot_bo_atm(1:natom) =  pot_bo_long_atm(1:natom) &
!     &                   + pot_bo_med_atm(i)
       &         + pot_bo_short_atm(1:natom)
  pot_an_atm(1:natom) =  pot_an_long_atm(1:natom) &
!     &                   + pot_an_med_atm(i)
       &         + pot_an_short_atm(1:natom)
  pot_to_atm(1:natom) =  pot_to_long_atm(1:natom) &
!     &                   + pot_to_med_atm(i)
       &         + pot_to_short_atm(1:natom)
  pot_14_atm(1:natom) =  pot_14_long_atm(1:natom) &
!     &                   + pot_14_med_atm(i)
       &         + pot_14_short_atm(1:natom)
  pot_elin_atm(1:natom) =  pot_elin_long_atm(1:natom) &
!     &                   + pot_elin_med_atm(i)
       &           + pot_elin_short_atm(1:natom)
  pot_ljin_atm(1:natom) =  pot_ljin_long_atm(1:natom) &
!     &                   + pot_ljin_med_atm(i)
       &           + pot_ljin_short_atm(1:natom)
  pot_elc_atm(1:natom) =  pot_elc_long_atm(1:natom) &
!     &                   + pot_elc_med_atm(i)
       &          + pot_elc_short_atm(1:natom)
  pot_vdw_atm(1:natom) =  pot_vdw_long_atm(1:natom) &
!     &                   + pot_vdw_med_atm(i)
       &          + pot_vdw_short_atm(1:natom)
  pot_mor_atm(1:natom) =  pot_mor_long_atm(1:natom) &
!     &                   + pot_mor_med_atm(i)
       &          + pot_mor_short_atm(1:natom)
  pot_sh_atm(1:natom) =  pot_sh_long_atm(1:natom) &
!     &                   + pot_sh_med_atm(i)
       &         + pot_sh_short_atm(1:natom)
  pot_rfh_atm(1:natom) =  pot_rfh_long_atm(1:natom) &
!     &                   + pot_rfh_med_atm(i)
       &          + pot_rfh_short_atm(1:natom)
  pot_dou_atm(1:natom) =  pot_dou_long_atm(1:natom) &
!     &                   + pot_dou_med_atm(i)
       &          + pot_dou_short_atm(1:natom)
  pot_cstmnb_atm(1:natom) =  pot_cstmnb_long_atm(1:natom) &
!     &                   + pot_cstmnb_med_atm(i)
       &          + pot_cstmnb_short_atm(1:natom)
  do nex = 1, ncstmnbex
     pot_cstmnbex_atm(1:natom,nex) = &
          &            pot_cstmnbex_long_atm(1:natom,nex) &
       ! &          + pot_cstmnbex_med_atm(1:natom,nex) &
       &          + pot_cstmnbex_short_atm(1:natom,nex)
  end do

  do iz = 1, nhfregion
     do i = 1, natom

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
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
!     &                                  + viricstmnbt_med_atm(1:3,n,i,iz)
                &                 + viricstmnbt_short_atm(1:3,n,i,iz)

           do nex = 1, ncstmnbex
              viricstmnbext_atm(1:3,n,i,iz,nex) = &
                   &                viricstmnbext_long_atm(1:3,n,i,iz,nex) &
                   ! &              + viricstmnbext_med_atm(1:3,n,i,iz,nex) &
                   &              + viricstmnbext_short_atm(1:3,n,i,iz,nex)
           end do

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
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
             &          + atmmass(i) * atmvel_tmp(1:3,i)

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
             &          + atmmass(i) * atmvel_tmp(1:3,i)

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
     molevel(1:3,im) = atmvel_tmp(1:3,i)

!         end do

  END DO
#endif


! -------- kinetic and potential part (transport by molecules) --------

  DO i = 1, natom

                                                             !!!!!!!!!!!!!!!!!ar1 -ar2 kinetic 
     if (atmindex(i) == 110) then
               nex = 1


            else if  (atmindex(i) == 111)then
               nex = 2
            else if  (atmindex(i) == 16)then
               nex = 3
            else
               stop "impossible value"
      end if

     if (hftyp_atm(i) == HFTYP_MOLE) then ! mole-base
        moleindex = irmolept_list(i)
        parcor(1:3) = molecom(1:3,moleindex)
        parvel(1:3) = molevel(1:3,moleindex)

!           for surface-based method only
        if (.not. ifhfvol) then
           parcor_sb(1:3) = molecom_hf(1:3,moleindex)
        end if

     else                   ! atom-base
        parcor(1:3) = atmcor(1:3,i)
        parvel(1:3) = atmvel_tmp(1:3,i)

!           for surface-based method only
        if (.not. ifhfvol) then
           parcor_sb(1:3) = atmcor_hf(1:3,i)
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

!        - loop over hfregion
     do iz = 1, nhfregion

! -------- local volume-based method --------
        IF (ifhfvol) THEN

           if ((hfzpos1(iz) <= parcor(3)) .and.  &
                & (parcor(3) <= hfzpos2(iz))) then

!              calculate kinetic part
              kin_ene =  0.5d0 * atmmass(i) &
                   &   * ( atmvel(1,i)*atmvel(1,i) &
                   &     + atmvel(2,i)*atmvel(2,i) &
                   &     + atmvel(3,i)*atmvel(3,i))

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfkin_atm(1,iz) = hfkin_atm(1,iz) + kin_ene * parvel(1)
              hfkin_atm(2,iz) = hfkin_atm(2,iz) + kin_ene * parvel(2)
#endif
              hfkin_atm(3,iz) = hfkin_atm(3,iz) + kin_ene * parvel(3)

!!! if you want to output some extra data about kinetic part,
!!!  write your procedures here.
#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
               hfkinex_atm(1,iz,nex) = hfkinex_atm(1,iz,nex) &
                    &                +  kin_ene * parvel(1)
               hfkinex_atm(2,iz,nex) = hfkinex_atm(2,iz,nex) &
                    &                +  kin_ene * parvel(2)
#endif
               hfkinex_atm(3,iz,nex) = hfkinex_atm(3,iz,nex) &
                    &                +  kin_ene * parvel(3)
!!!

!              calculate potential part
#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotbo_atm(1,iz) = hfpotbo_atm(1,iz) &
                   &            + pot_bo_atm(i) * parvel(1)
              hfpotbo_atm(2,iz) = hfpotbo_atm(2,iz) &
                   &            + pot_bo_atm(i) * parvel(2)
#endif
              hfpotbo_atm(3,iz) = hfpotbo_atm(3,iz) &
                   &            + pot_bo_atm(i) * parvel(3)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotan_atm(1,iz) = hfpotan_atm(1,iz) &
                   &            + pot_an_atm(i) * parvel(1)
              hfpotan_atm(2,iz) = hfpotan_atm(2,iz) &
                   &            + pot_an_atm(i) * parvel(2)
#endif
              hfpotan_atm(3,iz) = hfpotan_atm(3,iz) &
                   &            + pot_an_atm(i) * parvel(3)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotto_atm(1,iz) = hfpotto_atm(1,iz) &
                   &            + pot_to_atm(i) * parvel(1)
              hfpotto_atm(2,iz) = hfpotto_atm(2,iz) &
                   &            + pot_to_atm(i) * parvel(2)
#endif
              hfpotto_atm(3,iz) = hfpotto_atm(3,iz) &
                   &            + pot_to_atm(i) * parvel(3)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpot14_atm(1,iz) = hfpot14_atm(1,iz) &
                   &            + pot_14_atm(i) * parvel(1)
              hfpot14_atm(2,iz) = hfpot14_atm(2,iz) &
                   &            + pot_14_atm(i) * parvel(2)
#endif
              hfpot14_atm(3,iz) = hfpot14_atm(3,iz) &
                   &            + pot_14_atm(i) * parvel(3)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotelin_atm(1,iz) = hfpotelin_atm(1,iz) &
                   &              + pot_elin_atm(i) * parvel(1)
              hfpotelin_atm(2,iz) = hfpotelin_atm(2,iz) &
                   &              + pot_elin_atm(i) * parvel(2)
#endif
              hfpotelin_atm(3,iz) = hfpotelin_atm(3,iz) &
                   &              + pot_elin_atm(i) * parvel(3)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotljin_atm(1,iz) = hfpotljin_atm(1,iz) &
                   &              + pot_ljin_atm(i) * parvel(1)
              hfpotljin_atm(2,iz) = hfpotljin_atm(2,iz) &
                   &              + pot_ljin_atm(i) * parvel(2)
#endif
              hfpotljin_atm(3,iz) = hfpotljin_atm(3,iz) &
                   &              + pot_ljin_atm(i) * parvel(3)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotelc_atm(1,iz) = hfpotelc_atm(1,iz) &
                   &             + pot_elc_atm(i) * parvel(1)
              hfpotelc_atm(2,iz) = hfpotelc_atm(2,iz) &
                   &             + pot_elc_atm(i) * parvel(2)
#endif
              hfpotelc_atm(3,iz) = hfpotelc_atm(3,iz) &
                   &             + pot_elc_atm(i) * parvel(3)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotlj_atm(1,iz) = hfpotlj_atm(1,iz) &
                   &            + pot_vdw_atm(i) * parvel(1)
              hfpotlj_atm(2,iz) = hfpotlj_atm(2,iz) &
                   &            + pot_vdw_atm(i) * parvel(2)
#endif
              hfpotlj_atm(3,iz) = hfpotlj_atm(3,iz) &
                   &            + pot_vdw_atm(i) * parvel(3)
                                                                                !!!!!!!!!!!!!!! calculate ar1-ar2 potlj heat flux
#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotlj_atm_ar11(1,iz) =  hfpotlj_atm_ar11(1,iz) &
                                     +  pot_vdw_atm_ar11(i) * parvel(1)
            
              hfpotlj_atm_ar11(2,iz) =  hfpotlj_atm_ar11(2,iz) &
                                     +  pot_vdw_atm_ar11(i) * parvel(2)
#endif

              hfpotlj_atm_ar11(3,iz) =  hfpotlj_atm_ar11(3,iz) &
                                     +  pot_vdw_atm_ar11(i) * parvel(3)


#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotlj_atm_ar22(1,iz) =  hfpotlj_atm_ar22(1,iz) &
                                     +  pot_vdw_atm_ar22(i) * parvel(1)
            
              hfpotlj_atm_ar22(2,iz) =  hfpotlj_atm_ar22(2,iz) &
                                     +  pot_vdw_atm_ar22(i) * parvel(2)
#endif

              hfpotlj_atm_ar22(3,iz) =  hfpotlj_atm_ar22(3,iz) &
                                     +  pot_vdw_atm_ar22(i) * parvel(3)


#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotlj_atm_ar12(1,iz) =  hfpotlj_atm_ar12(1,iz) &
                                     +  pot_vdw_atm_ar12(i) * parvel(1)
            
              hfpotlj_atm_ar12(2,iz) =  hfpotlj_atm_ar12(2,iz) &
                                     +  pot_vdw_atm_ar12(i) * parvel(2)
#endif

              hfpotlj_atm_ar12(3,iz) =  hfpotlj_atm_ar12(3,iz) &
                                     +  pot_vdw_atm_ar12(i) * parvel(3)


                                                                                !!!!!!!!!!!!!!! calculate ar1-ar2 potlj heat flux
#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotmor_atm(1,iz) = hfpotmor_atm(1,iz) &
                   &             + pot_mor_atm(i) * parvel(1)
              hfpotmor_atm(2,iz) = hfpotmor_atm(2,iz) &
                   &             + pot_mor_atm(i) * parvel(2)
#endif
              hfpotmor_atm(3,iz) = hfpotmor_atm(3,iz) &
                   &             + pot_mor_atm(i) * parvel(3)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotsh_atm(1,iz) = hfpotsh_atm(1,iz) &
                   &            + pot_sh_atm(i) * parvel(1)
              hfpotsh_atm(2,iz) = hfpotsh_atm(2,iz) &
                   &            + pot_sh_atm(i) * parvel(2)
#endif
              hfpotsh_atm(3,iz) = hfpotsh_atm(3,iz) &
                   &            + pot_sh_atm(i) * parvel(3)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotrfh_atm(1,iz) = hfpotrfh_atm(1,iz) &
                   &             + pot_rfh_atm(i) * parvel(1)
              hfpotrfh_atm(2,iz) = hfpotrfh_atm(2,iz) &
                   &             + pot_rfh_atm(i) * parvel(2)
#endif
              hfpotrfh_atm(3,iz) = hfpotrfh_atm(3,iz) &
                   &             + pot_rfh_atm(i) * parvel(3)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotdou_atm(1,iz) = hfpotdou_atm(1,iz) &
                   &             + pot_dou_atm(i) * parvel(1)
              hfpotdou_atm(2,iz) = hfpotdou_atm(2,iz) &
                   &             + pot_dou_atm(i) * parvel(2)
#endif
              hfpotdou_atm(3,iz) = hfpotdou_atm(3,iz) &
                   &             + pot_dou_atm(i) * parvel(3)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
              hfpotcstmnb_atm(1,iz) = hfpotcstmnb_atm(1,iz) &
                   &                + pot_cstmnb_atm(i) * parvel(1)
              hfpotcstmnb_atm(2,iz) = hfpotcstmnb_atm(2,iz) &
                   &                + pot_cstmnb_atm(i) * parvel(2)
#endif
              hfpotcstmnb_atm(3,iz) = hfpotcstmnb_atm(3,iz) &
                   &                + pot_cstmnb_atm(i) * parvel(3)

              do nex = 1, ncstmnbex
#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
                 hfpotcstmnbex_atm(1,iz,nex) = &
                      &               hfpotcstmnbex_atm(1,iz,nex) &
                      &             + pot_cstmnbex_atm(i,nex) * parvel(1)
                 hfpotcstmnbex_atm(2,iz,nex) = &
                      &               hfpotcstmnbex_atm(2,iz,nex) &
                      &             + pot_cstmnbex_atm(i,nex) * parvel(2)
#endif
                 hfpotcstmnbex_atm(3,iz,nex) = &
                      &               hfpotcstmnbex_atm(3,iz,nex) &
                      &             + pot_cstmnbex_atm(i,nex) * parvel(3)
              end do

           end if

! -------- local surface-based method --------
        ELSE

!           !!! old information is not stored
           if (current_step == 1) exit

           parnewold(3) = parcor(3) - parcor_sb(3)
           parnewold_abs = dabs(parnewold(3))

           if (parnewold_abs * box_inv(3) >= 0.5d0) exit
                                ! parnewold exceed half cell

           parnew_p1 = parcor(3) - hfzpos1(iz)
           parold_p1 = parcor_sb(3) - hfzpos1(iz)

           if (parnew_p1*parold_p1 < 0.0d0) then

#if defined(_HF_DEBUG) || defined(_HF_DEBUG_CALHEATF)
              crosscount(iz) = crosscount(iz) + 1
#endif
              if (parnew_p1 > 0.0d0) then ! new parcor is positive
                 flagment_fact = 1.0d0
              else             ! old parcor is positive from the surface
                 flagment_fact = -1.0d0
              endif

              flagment_fact = flagment_fact / dt_cross

!              calculate kinetic part
              kin_ene =  0.5d0 * atmmass(i) &
                   &   * ( atmvel(1,i)*atmvel(1,i) &
                   &     + atmvel(2,i)*atmvel(2,i) &
                   &     + atmvel(3,i)*atmvel(3,i))

              hfkin_atm(3,iz) = hfkin_atm(3,iz) + kin_ene  &
                   &                            * flagment_fact

!!! if you want to output some extra data about kinetic part,
!!!  write your procedures here.
              ! nex = 1
              ! hfkinex_atm(3,iz,nex) = hfkinex_atm(3,iz,nex) &
              !      &                +  kin_ene * flagment_fact
!!!
                                                                                       !!!!!!!!!kinetic heat flux
!!!
               hfkinex_atm(3,iz,nex) = hfkinex_atm(3,iz,nex) &
                    &                +  kin_ene * flagment_fact
!              calculate potential part
              hfpotbo_atm(3,iz) = hfpotbo_atm(3,iz) &
                   &            + pot_bo_atm(i) * flagment_fact

              hfpotan_atm(3,iz) = hfpotan_atm(3,iz) &
                   &            + pot_an_atm(i) * flagment_fact

              hfpotto_atm(3,iz) = hfpotto_atm(3,iz) &
                   &            + pot_to_atm(i) * flagment_fact

              hfpot14_atm(3,iz) = hfpot14_atm(3,iz) &
                   &            + pot_14_atm(i) * flagment_fact

              hfpotelin_atm(3,iz) = hfpotelin_atm(3,iz) &
                   &              + pot_elin_atm(i) * flagment_fact

              hfpotljin_atm(3,iz) = hfpotljin_atm(3,iz) &
                   &              + pot_ljin_atm(i) * flagment_fact

              hfpotelc_atm(3,iz) = hfpotelc_atm(3,iz) &
                   &             + pot_elc_atm(i) * flagment_fact

              hfpotlj_atm(3,iz) = hfpotlj_atm(3,iz) &
                   &            + pot_vdw_atm(i) * flagment_fact

              hfpotmor_atm(3,iz) = hfpotmor_atm(3,iz) &
                   &             + pot_mor_atm(i) * flagment_fact

              hfpotsh_atm(3,iz) = hfpotsh_atm(3,iz) &
                   &            + pot_sh_atm(i) * flagment_fact

              hfpotrfh_atm(3,iz) = hfpotrfh_atm(3,iz) &
                   &             + pot_rfh_atm(i) * flagment_fact

              hfpotdou_atm(3,iz) = hfpotdou_atm(3,iz) &
                   &             + pot_dou_atm(i) * flagment_fact

              hfpotcstmnb_atm(3,iz) = hfpotcstmnb_atm(3,iz) &
                   &                + pot_cstmnb_atm(i) * flagment_fact

              hfpotcstmnbex_atm(3,iz,1:ncstmnbex) = &
                   &           hfpotcstmnbex_atm(3,iz,1:ncstmnbex) &
                   &         + pot_cstmnbex_atm(i,1:ncstmnbex) * flagment_fact
             
                                                                                  !!!!!!!!!!!!!!! calculate ar1-ar2 potlj heat flux

              hfpotlj_atm_ar11(3,iz) =  hfpotlj_atm_ar11(3,iz) &
                                     +  pot_vdw_atm_ar11(i) * flagment_fact



              hfpotlj_atm_ar22(3,iz) =  hfpotlj_atm_ar22(3,iz) &
                                     +  pot_vdw_atm_ar22(i) * flagment_fact



              hfpotlj_atm_ar12(3,iz) =  hfpotlj_atm_ar12(3,iz) &
                                     +  pot_vdw_atm_ar12(i) * flagment_fact


                                                                                !!!!!!!!!!!!!!! calculate ar1-ar2 potlj heat flux


           end if

        END IF                 ! local surface-based method end

     end do

  END DO

!---- store current coordinate information (only for surface-based)
#if defined(_HF_DEBUG) || (_HF_DEBUG_CALHEATF)
  write(6,*) '*** heatf debug info (in calheatf) ***'

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
     molecom_hf(1:3,1:nmole) = molecom(1:3,1:nmole)

     atmcor_hf(1:3,1:natom) = atmcor(1:3,1:natom)

  end if

  !!! When surface-based, exit subroutine other than output MD step here
  if (.not. ifhfvol .and. (mod(current_step,heatfinterval) /= 1) .and. &
  &   (heatfinterval /= 1)) then
      !--- restore velocity
      atmvel(1:3,1:natom) = atmvel_tmp(1:3,1:natom)
      return
  end if

! -------- collision part (transport by interaction) --------
  do iz = 1, nhfregion

     do i = 1, natom
#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfviribo_atm(1,iz) = hfviribo_atm(1,iz) &
             &             + viribot_atm(1,1,i,iz) * atmvel(1,i) &
             &             + viribot_atm(2,1,i,iz) * atmvel(2,i) &
             &             + viribot_atm(3,1,i,iz) * atmvel(3,i)
        hfviribo_atm(2,iz) = hfviribo_atm(2,iz) &
             &             + viribot_atm(1,2,i,iz) * atmvel(1,i) &
             &             + viribot_atm(2,2,i,iz) * atmvel(2,i) &
             &             + viribot_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfviribo_atm(3,iz) = hfviribo_atm(3,iz) &
             &             + viribot_atm(1,3,i,iz) * atmvel(1,i) &
             &             + viribot_atm(2,3,i,iz) * atmvel(2,i) &
             &             + viribot_atm(3,3,i,iz) * atmvel(3,i)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfvirian_atm(1,iz) = hfvirian_atm(1,iz) &
             &             + viriant_atm(1,1,i,iz) * atmvel(1,i) &
             &             + viriant_atm(2,1,i,iz) * atmvel(2,i) &
             &             + viriant_atm(3,1,i,iz) * atmvel(3,i)
        hfvirian_atm(2,iz) = hfvirian_atm(2,iz) &
             &             + viriant_atm(1,2,i,iz) * atmvel(1,i) &
             &             + viriant_atm(2,2,i,iz) * atmvel(2,i) &
             &             + viriant_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfvirian_atm(3,iz) = hfvirian_atm(3,iz) &
             &             + viriant_atm(1,3,i,iz) * atmvel(1,i) &
             &             + viriant_atm(2,3,i,iz) * atmvel(2,i) &
             &             + viriant_atm(3,3,i,iz) * atmvel(3,i)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfvirito_atm(1,iz) = hfvirito_atm(1,iz) &
             &             + viritot_atm(1,1,i,iz) * atmvel(1,i) &
             &             + viritot_atm(2,1,i,iz) * atmvel(2,i) &
             &             + viritot_atm(3,1,i,iz) * atmvel(3,i)
        hfvirito_atm(2,iz) = hfvirito_atm(2,iz) &
             &             + viritot_atm(1,2,i,iz) * atmvel(1,i) &
             &             + viritot_atm(2,2,i,iz) * atmvel(2,i) &
             &             + viritot_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfvirito_atm(3,iz) = hfvirito_atm(3,iz) &
             &             + viritot_atm(1,3,i,iz) * atmvel(1,i) &
             &             + viritot_atm(2,3,i,iz) * atmvel(2,i) &
             &             + viritot_atm(3,3,i,iz) * atmvel(3,i)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfviri14_atm(1,iz) = hfviri14_atm(1,iz) &
             &             + viri14t_atm(1,1,i,iz) * atmvel(1,i) &
             &             + viri14t_atm(2,1,i,iz) * atmvel(2,i) &
             &             + viri14t_atm(3,1,i,iz) * atmvel(3,i)
        hfviri14_atm(2,iz) = hfviri14_atm(2,iz) &
             &             + viri14t_atm(1,2,i,iz) * atmvel(1,i) &
             &             + viri14t_atm(2,2,i,iz) * atmvel(2,i) &
             &             + viri14t_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfviri14_atm(3,iz) = hfviri14_atm(3,iz) &
             &             + viri14t_atm(1,3,i,iz) * atmvel(1,i) &
             &             + viri14t_atm(2,3,i,iz) * atmvel(2,i) &
             &             + viri14t_atm(3,3,i,iz) * atmvel(3,i)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfvirielin_atm(1,iz) = hfvirielin_atm(1,iz) &
             &               + virielint_atm(1,1,i,iz) * atmvel(1,i) &
             &               + virielint_atm(2,1,i,iz) * atmvel(2,i) &
             &               + virielint_atm(3,1,i,iz) * atmvel(3,i)
        hfvirielin_atm(2,iz) = hfvirielin_atm(2,iz) &
             &               + virielint_atm(1,2,i,iz) * atmvel(1,i) &
             &               + virielint_atm(2,2,i,iz) * atmvel(2,i) &
             &               + virielint_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfvirielin_atm(3,iz) = hfvirielin_atm(3,iz) &
             &               + virielint_atm(1,3,i,iz) * atmvel(1,i) &
             &               + virielint_atm(2,3,i,iz) * atmvel(2,i) &
             &               + virielint_atm(3,3,i,iz) * atmvel(3,i)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfviriljin_atm(1,iz) = hfviriljin_atm(1,iz) &
             &               + viriljint_atm(1,1,i,iz) * atmvel(1,i) &
             &               + viriljint_atm(2,1,i,iz) * atmvel(2,i) &
             &               + viriljint_atm(3,1,i,iz) * atmvel(3,i)
        hfviriljin_atm(2,iz) = hfviriljin_atm(2,iz) &
             &               + viriljint_atm(1,2,i,iz) * atmvel(1,i) &
             &               + viriljint_atm(2,2,i,iz) * atmvel(2,i) &
             &               + viriljint_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfviriljin_atm(3,iz) = hfviriljin_atm(3,iz) &
             &               + viriljint_atm(1,3,i,iz) * atmvel(1,i) &
             &               + viriljint_atm(2,3,i,iz) * atmvel(2,i) &
             &               + viriljint_atm(3,3,i,iz) * atmvel(3,i)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfvirielc_atm(1,iz) = hfvirielc_atm(1,iz) &
             &              + virielct_atm(1,1,i,iz) * atmvel(1,i) &
             &              + virielct_atm(2,1,i,iz) * atmvel(2,i) &
             &              + virielct_atm(3,1,i,iz) * atmvel(3,i)
        hfvirielc_atm(2,iz) = hfvirielc_atm(2,iz) &
             &              + virielct_atm(1,2,i,iz) * atmvel(1,i) &
             &              + virielct_atm(2,2,i,iz) * atmvel(2,i) &
             &              + virielct_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfvirielc_atm(3,iz) = hfvirielc_atm(3,iz) &
             &              + virielct_atm(1,3,i,iz) * atmvel(1,i) &
             &              + virielct_atm(2,3,i,iz) * atmvel(2,i) &
             &              + virielct_atm(3,3,i,iz) * atmvel(3,i)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfvirilj_atm(1,iz) = hfvirilj_atm(1,iz) &
             &             + viriljt_atm(1,1,i,iz) * atmvel(1,i) &
             &             + viriljt_atm(2,1,i,iz) * atmvel(2,i) &
             &             + viriljt_atm(3,1,i,iz) * atmvel(3,i)
        hfvirilj_atm(2,iz) = hfvirilj_atm(2,iz) &
             &             + viriljt_atm(1,2,i,iz) * atmvel(1,i) &
             &             + viriljt_atm(2,2,i,iz) * atmvel(2,i) &
             &             + viriljt_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfvirilj_atm(3,iz) = hfvirilj_atm(3,iz) &
             &             + viriljt_atm(1,3,i,iz) * atmvel(1,i) &
             &             + viriljt_atm(2,3,i,iz) * atmvel(2,i) &
             &             + viriljt_atm(3,3,i,iz) * atmvel(3,i)

                                                                                !!!!!!!!!!!!!!! calculate ar1-ar2 virilj heat flux
       hfvirilj_atm_ar11(1,iz) = hfvirilj_atm_ar11(1,iz) &
             &             + viriljt_atm_ar11(1,1,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar11(2,1,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar11(3,1,i,iz) * atmvel(3,i)

       hfvirilj_atm_ar11(2,iz) = hfvirilj_atm_ar11(2,iz) &
             &             + viriljt_atm_ar11(1,2,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar11(2,2,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar11(3,2,i,iz) * atmvel(3,i)

       hfvirilj_atm_ar11(3,iz) = hfvirilj_atm_ar11(3,iz) &
             &             + viriljt_atm_ar11(1,3,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar11(2,3,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar11(3,3,i,iz) * atmvel(3,i)


       hfvirilj_atm_ar12(1,iz) = hfvirilj_atm_ar12(1,iz) &
             &             + viriljt_atm_ar12(1,1,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar12(2,1,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar12(3,1,i,iz) * atmvel(3,i)

       hfvirilj_atm_ar12(2,iz) = hfvirilj_atm_ar12(2,iz) &
             &             + viriljt_atm_ar12(1,2,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar12(2,2,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar12(3,2,i,iz) * atmvel(3,i)

       hfvirilj_atm_ar12(3,iz) = hfvirilj_atm_ar12(3,iz) &
             &             + viriljt_atm_ar12(1,3,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar12(2,3,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar12(3,3,i,iz) * atmvel(3,i)


       hfvirilj_atm_ar22(1,iz) = hfvirilj_atm_ar22(1,iz) &
             &             + viriljt_atm_ar22(1,1,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar22(2,1,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar22(3,1,i,iz) * atmvel(3,i)

       hfvirilj_atm_ar22(2,iz) = hfvirilj_atm_ar22(2,iz) &
             &             + viriljt_atm_ar22(1,2,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar22(2,2,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar22(3,2,i,iz) * atmvel(3,i)

       hfvirilj_atm_ar22(3,iz) = hfvirilj_atm_ar22(3,iz) &
             &             + viriljt_atm_ar22(1,3,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar22(2,3,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar22(3,3,i,iz) * atmvel(3,i)



       hfvirilj_atm_ar1pt(1,iz) = hfvirilj_atm_ar1pt(1,iz) &
             &             + viriljt_atm_ar1pt(1,1,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar1pt(2,1,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar1pt(3,1,i,iz) * atmvel(3,i)

       hfvirilj_atm_ar1pt(2,iz) = hfvirilj_atm_ar1pt(2,iz) &
             &             + viriljt_atm_ar1pt(1,2,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar1pt(2,2,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar1pt(3,2,i,iz) * atmvel(3,i)

       hfvirilj_atm_ar1pt(3,iz) = hfvirilj_atm_ar1pt(3,iz) &
             &             + viriljt_atm_ar1pt(1,3,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar1pt(2,3,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar1pt(3,3,i,iz) * atmvel(3,i)



       hfvirilj_atm_ar2pt(1,iz) = hfvirilj_atm_ar2pt(1,iz) &
             &             + viriljt_atm_ar2pt(1,1,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar2pt(2,1,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar2pt(3,1,i,iz) * atmvel(3,i)

       hfvirilj_atm_ar2pt(2,iz) = hfvirilj_atm_ar2pt(2,iz) &
             &             + viriljt_atm_ar2pt(1,2,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar2pt(2,2,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar2pt(3,2,i,iz) * atmvel(3,i)

       hfvirilj_atm_ar2pt(3,iz) = hfvirilj_atm_ar2pt(3,iz) &
             &             + viriljt_atm_ar2pt(1,3,i,iz) * atmvel(1,i) &
             &             + viriljt_atm_ar2pt(2,3,i,iz) * atmvel(2,i) &
             &             + viriljt_atm_ar2pt(3,3,i,iz) * atmvel(3,i)
                                                                               !!!!!!!!!!!!!!! calculate ar1-ar2 virilj heat flux

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfvirimor_atm(1,iz) = hfvirimor_atm(1,iz) &
             &              + virimort_atm(1,1,i,iz) * atmvel(1,i) &
             &              + virimort_atm(2,1,i,iz) * atmvel(2,i) &
             &              + virimort_atm(3,1,i,iz) * atmvel(3,i)
        hfvirimor_atm(2,iz) = hfvirimor_atm(2,iz) &
             &              + virimort_atm(1,2,i,iz) * atmvel(1,i) &
             &              + virimort_atm(2,2,i,iz) * atmvel(2,i) &
             &              + virimort_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfvirimor_atm(3,iz) = hfvirimor_atm(3,iz) &
             &              + virimort_atm(1,3,i,iz) * atmvel(1,i) &
             &              + virimort_atm(2,3,i,iz) * atmvel(2,i) &
             &              + virimort_atm(3,3,i,iz) * atmvel(3,i)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfvirish_atm(1,iz) = hfvirish_atm(1,iz) &
             &             + virisht_atm(1,1,i,iz) * atmvel(1,i) &
             &             + virisht_atm(2,1,i,iz) * atmvel(2,i) &
             &             + virisht_atm(3,1,i,iz) * atmvel(3,i)
        hfvirish_atm(2,iz) = hfvirish_atm(2,iz) &
             &             + virisht_atm(1,2,i,iz) * atmvel(1,i) &
             &             + virisht_atm(2,2,i,iz) * atmvel(2,i) &
             &             + virisht_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfvirish_atm(3,iz) = hfvirish_atm(3,iz) &
             &             + virisht_atm(1,3,i,iz) * atmvel(1,i) &
             &             + virisht_atm(2,3,i,iz) * atmvel(2,i) &
             &             + virisht_atm(3,3,i,iz) * atmvel(3,i)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfvirirfh_atm(1,iz) = hfvirirfh_atm(1,iz) &
             &              + virirfht_atm(1,1,i,iz) * atmvel(1,i) &
             &              + virirfht_atm(2,1,i,iz) * atmvel(2,i) &
             &              + virirfht_atm(3,1,i,iz) * atmvel(3,i)
        hfvirirfh_atm(2,iz) = hfvirirfh_atm(2,iz) &
             &              + virirfht_atm(1,2,i,iz) * atmvel(1,i) &
             &              + virirfht_atm(2,2,i,iz) * atmvel(2,i) &
             &              + virirfht_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfvirirfh_atm(3,iz) = hfvirirfh_atm(3,iz) &
             &              + virirfht_atm(1,3,i,iz) * atmvel(1,i) &
             &              + virirfht_atm(2,3,i,iz) * atmvel(2,i) &
             &              + virirfht_atm(3,3,i,iz) * atmvel(3,i)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfviridou_atm(1,iz) = hfviridou_atm(1,iz) &
             &              + viridout_atm(1,1,i,iz) * atmvel(1,i) &
             &              + viridout_atm(2,1,i,iz) * atmvel(2,i) &
             &              + viridout_atm(3,1,i,iz) * atmvel(3,i)
        hfviridou_atm(2,iz) = hfviridou_atm(2,iz) &
             &              + viridout_atm(1,2,i,iz) * atmvel(1,i) &
             &              + viridout_atm(2,2,i,iz) * atmvel(2,i) &
             &              + viridout_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfviridou_atm(3,iz) = hfviridou_atm(3,iz) &
             &              + viridout_atm(1,3,i,iz) * atmvel(1,i) &
             &              + viridout_atm(2,3,i,iz) * atmvel(2,i) &
             &              + viridout_atm(3,3,i,iz) * atmvel(3,i)

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
        hfviricstmnb_atm(1,iz) = hfviricstmnb_atm(1,iz) &
             &                 + viricstmnbt_atm(1,1,i,iz) * atmvel(1,i) &
             &                 + viricstmnbt_atm(2,1,i,iz) * atmvel(2,i) &
             &                 + viricstmnbt_atm(3,1,i,iz) * atmvel(3,i)
        hfviricstmnb_atm(2,iz) = hfviricstmnb_atm(2,iz) &
             &                 + viricstmnbt_atm(1,2,i,iz) * atmvel(1,i) &
             &                 + viricstmnbt_atm(2,2,i,iz) * atmvel(2,i) &
             &                 + viricstmnbt_atm(3,2,i,iz) * atmvel(3,i)
#endif
        hfviricstmnb_atm(3,iz) = hfviricstmnb_atm(3,iz) &
             &                 + viricstmnbt_atm(1,3,i,iz) * atmvel(1,i) &
             &                 + viricstmnbt_atm(2,3,i,iz) * atmvel(2,i) &
             &                 + viricstmnbt_atm(3,3,i,iz) * atmvel(3,i)

        do nex = 1, ncstmnbex
#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
           hfviricstmnbex_atm(1,iz,nex) = &
                &          hfviricstmnbex_atm(1,iz,nex) &
                &        + viricstmnbext_atm(1,1,i,iz,nex) * atmvel(1,i) &
                &        + viricstmnbext_atm(2,1,i,iz,nex) * atmvel(2,i) &
                &        + viricstmnbext_atm(3,1,i,iz,nex) * atmvel(3,i)
           hfviricstmnbex_atm(2,iz,nex) = &
                &          hfviricstmnbex_atm(2,iz,nex) &
                &        + viricstmnbext_atm(1,2,i,iz,nex) * atmvel(1,i) &
                &        + viricstmnbext_atm(2,2,i,iz,nex) * atmvel(2,i) &
                &        + viricstmnbext_atm(3,2,i,iz,nex) * atmvel(3,i)
#endif
           hfviricstmnbex_atm(3,iz,nex) = &
                &          hfviricstmnbex_atm(3,iz,nex) &
                &        + viricstmnbext_atm(1,3,i,iz,nex) * atmvel(1,i) &
                &        + viricstmnbext_atm(2,3,i,iz,nex) * atmvel(2,i) &
                &        + viricstmnbext_atm(3,3,i,iz,nex) * atmvel(3,i)
        end do

     end do
  end do

! --- sum up heat flux ---

  do iz = 1, nhfregion

     if (ifhfvol) then      ! volume-based (= volume)
        vol_inv = area_inv / (hfzpos2(iz) - hfzpos1(iz))
     else                   ! surface-based (= surface area)
        vol_inv = area_inv
     end if

     hfkin_atm(1:3,iz) =  hfkin_atm(1:3,iz) * vol_inv

     do nex = 1, ncstmnbex
        hfkinex_atm(1:3,iz,nex) = hfkinex_atm(1:3,iz,nex) * vol_inv
     end do

     hfpotbo_atm(1:3,iz) =  hfpotbo_atm(1:3,iz) * vol_inv

     hfpotan_atm(1:3,iz) =  hfpotan_atm(1:3,iz) * vol_inv

     hfpotto_atm(1:3,iz) =  hfpotto_atm(1:3,iz) * vol_inv

     hfpot14_atm(1:3,iz) =  hfpot14_atm(1:3,iz) * vol_inv

     hfpotelin_atm(1:3,iz) =  hfpotelin_atm(1:3,iz) * vol_inv

     hfpotljin_atm(1:3,iz) =  hfpotljin_atm(1:3,iz) * vol_inv

     hfpotelc_atm(1:3,iz) =  hfpotelc_atm(1:3,iz) * vol_inv

     hfpotlj_atm(1:3,iz) =  hfpotlj_atm(1:3,iz) * vol_inv

     hfpotmor_atm(1:3,iz) =  hfpotmor_atm(1:3,iz) * vol_inv

     hfpotsh_atm(1:3,iz) =  hfpotsh_atm(1:3,iz) * vol_inv

     hfpotrfh_atm(1:3,iz) =  hfpotrfh_atm(1:3,iz) * vol_inv

     hfpotdou_atm(1:3,iz) =  hfpotdou_atm(1:3,iz) * vol_inv

     hfpotcstmnb_atm(1:3,iz) =  hfpotcstmnb_atm(1:3,iz) * vol_inv

                                                                      !!!!!!!!!!!!!!! calculate ar1-ar2 potlj heat flux
     hfpotlj_atm_ar11(1:3,iz) =  hfpotlj_atm_ar11(1:3,iz) * vol_inv


     hfpotlj_atm_ar22(1:3,iz) =  hfpotlj_atm_ar22(1:3,iz)* vol_inv


     hfpotlj_atm_ar12(1:3,iz) =  hfpotlj_atm_ar12(1:3,iz)* vol_inv

                                                                           !!!!!!!!!!!!!!! calculate ar1-ar2 potlj heat flux   


     do nex = 1, ncstmnbex
        hfpotcstmnbex_atm(1:3,iz,nex) = &
             &       hfpotcstmnbex_atm(1:3,iz,nex) * vol_inv
     end do

     hfviribo_atm(1:3,iz) =  hfviribo_atm(1:3,iz) * vol_inv

     hfvirian_atm(1:3,iz) =  hfvirian_atm(1:3,iz) * vol_inv

     hfvirito_atm(1:3,iz) =  hfvirito_atm(1:3,iz) * vol_inv

     hfviri14_atm(1:3,iz) =  hfviri14_atm(1:3,iz) * vol_inv

     hfvirielin_atm(1:3,iz) =  hfvirielin_atm(1:3,iz) * vol_inv

     hfviriljin_atm(1:3,iz) =  hfviriljin_atm(1:3,iz) * vol_inv

     hfvirielc_atm(1:3,iz) =  hfvirielc_atm(1:3,iz) * vol_inv

     hfvirilj_atm(1:3,iz) =  hfvirilj_atm(1:3,iz) * vol_inv

     hfvirimor_atm(1:3,iz) =  hfvirimor_atm(1:3,iz) * vol_inv

     hfvirish_atm(1:3,iz) =  hfvirish_atm(1:3,iz) * vol_inv

     hfvirirfh_atm(1:3,iz) =  hfvirirfh_atm(1:3,iz) * vol_inv

     hfviridou_atm(1:3,iz) =  hfviridou_atm(1:3,iz) * vol_inv

     hfviricstmnb_atm(1:3,iz) =  hfviricstmnb_atm(1:3,iz) * vol_inv

                                                                      !!!!!!!!!!!!!!! calculate ar1-ar2 viri heat flux
      hfvirilj_atm_ar11(1:3,iz) = hfvirilj_atm_ar11(1:3,iz)  * vol_inv

      hfvirilj_atm_ar22(1:3,iz) = hfvirilj_atm_ar22(1:3,iz)  * vol_inv

      hfvirilj_atm_ar12(1:3,iz) = hfvirilj_atm_ar12(1:3,iz)  * vol_inv

      hfvirilj_atm_ar1pt(1:3,iz) = hfvirilj_atm_ar1pt(1:3,iz)  * vol_inv

      hfvirilj_atm_ar2pt(1:3,iz) = hfvirilj_atm_ar2pt(1:3,iz)  * vol_inv

                                                                          !!!!!!!!!!!!!!! calculate ar1-ar2 viri heat flux


     do nex = 1, ncstmnbex
        hfviricstmnbex_atm(1:3,iz,nex) = &
             &        hfviricstmnbex_atm(1:3,iz,nex) * vol_inv
     end do

     hftot_atm(1:3,iz) = hfkin_atm(1:3,iz) &
          &            + hfpotbo_atm(1:3,iz) + hfpotan_atm(1:3,iz) &
          &            + hfpotto_atm(1:3,iz) + hfpot14_atm(1:3,iz) &
          &            + hfpotelin_atm(1:3,iz) + hfpotljin_atm(1:3,iz) &
          &            + hfpotelc_atm(1:3,iz) &
          &            + hfpotlj_atm(1:3,iz) + hfpotmor_atm(1:3,iz) &
          &            + hfpotsh_atm(1:3,iz) + hfpotrfh_atm(1:3,iz) &
          &            + hfpotdou_atm(1:3,iz) &
          &            + hfpotcstmnb_atm(1:3,iz) &
          &            + hfviribo_atm(1:3,iz) + hfvirian_atm(1:3,iz) &
          &            + hfvirito_atm(1:3,iz) + hfviri14_atm(1:3,iz) &
          &            + hfvirielin_atm(1:3,iz) + hfviriljin_atm(1:3,iz) &
          &            + hfvirielc_atm(1:3,iz) &
          &            + hfvirilj_atm(1:3,iz) + hfvirimor_atm(1:3,iz) &
          &            + hfvirish_atm(1:3,iz) + hfvirirfh_atm(1:3,iz) &
          &            + hfviridou_atm(1:3,iz) &
          &            + hfviricstmnb_atm(1:3,iz)

#if defined(_CSTMNB_V2_ADD_ALL)
     !- adding cstmnb extra interaction
     do nex = 1, ncstmnbex
        hftot_atm(1:3,iz) = hftot_atm(1:3,iz) &
             &            + hfpotcstmnbex_atm(1:3,iz,nex) &
             &            + hfviricstmnbex_atm(1:3,iz,nex)
     end do
#endif

  end do

!--- restore velocity
  atmvel(1:3,1:natom) = atmvel_tmp(1:3,1:natom)

!     +     +     +     +     +     +     +

end subroutine calheatf
