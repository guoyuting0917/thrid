!*****************************
!*  calforce_hf.f90 Ver.3.0  *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-07-17 16:44:51 gota>

subroutine calforce_hf(mts_flag,force, &
     &                 for_viri_coul,pot_viri_coul, &
     &                 for_viri_lj,pot_viri_lj, &
     &                 for_viri_mor,pot_viri_mor, &
     &                 for_viri_sh,pot_viri_sh, &
     &                 for_viri_rfh,pot_viri_rfh, &
     &                 for_viri_dou,pot_viri_dou, &
     &                 for_viri_cstmnb,pot_viri_cstmnb, &
     &                 for_viri_cstmnbex,pot_viri_cstmnbex, &
     &                 atm_viri_bond,atm_viri_angl, &
     &                 atm_viri_tors,atm_viri_14, &
     &                 pot_virit_coul,pot_virit_lj, &
     &                 pot_virit_mor, &
     &                 pot_virit_sh, &
     &                 pot_virit_rfh, &
     &                 pot_virit_dou, &
     &                 pot_virit_cstmnb, &
     &                 pot_virit_cstmnbex, &
     &                 atm_virit_bond,atm_virit_angl, &
     &                 atm_virit_tors,atm_virit_14, &
     &                 pot_bo_atm, &
     &                 pot_an_atm, &
     &                 pot_to_atm, &
     &                 pot_14_atm, &
     &                 pot_elin_atm,pot_ljin_atm, &
     &                 pot_elc_atm, &
     &                 pot_vdw_atm, &
     &                 pot_mor_atm, &
     &                 pot_sh_atm, &
     &                 pot_rfh_atm, &
     &                 pot_dou_atm, &
     &                 pot_cstmnb_atm, &
     &                 pot_cstmnbex_atm, &
     &                 viribot_atm, &
     &                 viriant_atm, &
     &                 viritot_atm, &
     &                 viri14t_atm, &
     &                 virielint_atm,viriljint_atm, &
     &                 virielct_atm, &
     &                 viriljt_atm, &
     &                 virimort_atm, &
     &                 virisht_atm, &
     &                 virirfht_atm, &
     &                 viridout_atm, &
     &                 viricstmnbt_atm, &
     &                 viricstmnbext_atm, &
     &                 mts_bond, mts_angl, mts_anglub, &
     &                 mts_tors,mts_torsrb,mts_torsim, &
     &                 mts_vdw,mts_ewr,mts_ewk, &
     &                 mts_vdw14,mts_elc14, &
     &                 mts_mor,mts_sh,mts_rfh,mts_dou,mts_cnpvw, &
     &                 mts_cstmnb, &
     &                 mts_posres, &
     &                 mts_potbias, &
     &                 npoly,npolytyp,npoly_mole,npoly_atom, &
     &                 nwater,nmatom,nmatyp,nmatomtyp, &
     &                 xcel,ycel,zcel, &
     &                 ifewald,alpha,kmax,rrcut, &
     &                 ifspme,nfft1,nfft2,nfft3,pme_order, &
     &                 eps0,div_factor_14vdw,div_factor_14elc, &
     &                 rcut, &
     &                 rcutmor, &
     &                 rcutsh, &
     &                 rcutrfhfo,rcutrfhoo,rcutrfhoh, &
     &                 rcutdouo,rcutindouo,rcutdouh, &
     &                 rcutrpvw, &
     &                 pot_ewc,pot_ewnq, &
     &                 pot_tot,pot_nonbon,pot_vdw,pot_elc,pot_ewk, &
     &                 pot_vdw14,pot_elc14, &
     &                 pot_bond,pot_angl,pot_anglub, &
     &                 pot_tors,pot_torsrb,pot_torsim, &
     &                 pot_mor, &
     &                 pot_sh, &
     &                 pot_rfh, &
     &                 pot_dou, &
     &                 pot_rpvw,pot_vw, &
     &                 pot_cstmnb, pot_cstmnbex, &
     &                 pot_posres, &
     &                 pot_pbias, &
     &                 ifcalpremole,ifcalpreatom, &
     &                 ifnetqcorrp, &
     &                 netchrgsq, &
     &                 ifcstmnb, &
     &                 ifposres, &
     &                 ifpotbias, &
     &                 ifhfvol, &
     &                 nhfregion,hfzpos1,hfzpos2, &
     &                 hftyp_atm, &
     &                 molecom, &
     &                 outinterval,pressinterval,heatfinterval, &
     &                 recinterval, &
     &                 nstep_book, &
     &                 md_cont, &
     &                 ncstmnbex, &
     &                 current_step, &
     &                 nstep_short,istep_short)

  use interface_interact
#if defined(TIME_M) || defined(TIME_MALL)
  use interface_timer
#endif

  use md_global

#if defined(TIME_M)
  use time_global
#endif

  implicit none

!     subroutine to calculate force

! ARGUMENT:
!     INPUT
  integer,intent(in):: mts_flag        ! flag for MTS integration
                                       ! 1 long-range force mode
                                       ! 2 medium-range force mode
                                       ! 3 short-range force mode

  integer,intent(in):: mts_bond        ! MTS flag for bond
  integer,intent(in):: mts_angl        ! MTS flag for angle
  integer,intent(in):: mts_anglub      ! MTS flag for Urey-Bradley angle
  integer,intent(in):: mts_tors        ! MTS flag for torsion
  integer,intent(in):: mts_torsrb      ! MTS flag for torsionrb
  integer,intent(in):: mts_torsim      ! MTS flag for torsionim
  integer,intent(in):: mts_vdw         ! MTS flag for vdw interaction
  integer,intent(in):: mts_ewr         ! MTS flag for ewald real(=vdw)
  integer,intent(in):: mts_ewk         ! MTS flag for ewald wave
  integer,intent(in):: mts_vdw14       ! MTS flag for 14vdw
  integer,intent(in):: mts_elc14       ! MTS flag for 14elc
  integer,intent(in):: mts_mor         ! MTS flag for Morse interaction
  integer,intent(in):: mts_sh          ! MTS flag for SH interaction
  integer,intent(in):: mts_rfh         ! MTS flag for RFH interaction
  integer,intent(in):: mts_dou         ! MTS flag for DOU interaction
  integer,intent(in):: mts_cnpvw       ! MTS flag for CNP_VW
  integer,intent(in):: mts_cstmnb      ! MTS flag for custom NB interaction
  integer,intent(in):: mts_posres      ! MTS flag for position restraint
  integer,intent(in):: mts_potbias     ! MTS flag for bias potential

  integer,intent(in):: npoly           ! all number of poly
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
  
  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatom          ! number of monatomic molecules
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: ifewald         ! ewald flag
  real(8),intent(in):: alpha            ! parameter alpha [non-d]
  integer,intent(in):: kmax            ! parameter kmax
  real(8),intent(in):: rrcut            ! ewald real space cutoff length [non-d]

  logical,intent(in):: ifspme          ! SPME (Smooth Particle Mesh Ewald) flag
  integer,intent(in):: nfft1, nfft2, nfft3 ! grid points in SPME
  integer,intent(in):: pme_order       ! B-spline order

  real(8),intent(in):: eps0             ! dielectric constant of vacuum [non-d]
  real(8),intent(in):: div_factor_14vdw ! division factor of 14vdw
  real(8),intent(in):: div_factor_14elc ! division factor of 14elc
  real(8),intent(in):: rcut             ! vdw cutoff length [non-d]

  real(8),intent(in):: rcutmor          ! Morse cutoff length [non-d]

  real(8),intent(in):: rcutsh           ! Morse cutoff length [non-d]

  real(8),intent(in):: rcutrfhfo        ! RFH(FO) cutoff length [non-d]
  real(8),intent(in):: rcutrfhoo        ! RFH(OO) cutoff length [non-d]
  real(8),intent(in):: rcutrfhoh        ! RFH(OH) cutoff length [non-d]

  real(8),intent(in):: rcutdouo         ! DOU cutoff length [non-d] for O-Au
  real(8),intent(in):: rcutindouo       ! DOU cutin length [non-d] for O-Au
  real(8),intent(in):: rcutdouh         ! DOU cutoff length [non-d] for H-Au

  real(8),intent(in):: rcutrpvw         ! RP-VW cutoff length

  logical,intent(in):: ifcstmnb         ! flag if using custom NB interaction

!---- potential valiables
  real(8),intent(in):: pot_ewc          ! potential of ewald self-energy

!---- valiables for pressure

  logical,intent(in):: ifcalpremole    ! pressure calculation of molecule
  logical,intent(in):: ifcalpreatom    ! pressure calculation of atom
  logical,intent(in):: ifnetqcorrp     ! net charge correction for pressure

!---- variables for net charge calculation
  real(8),intent(in):: netchrgsq        ! = (sum(qi))**2

!---- variables for position restraint
  logical,intent(in):: ifposres        ! position restraint flag

!---- variable for bias potential
  logical,intent(in):: ifpotbias       ! bias potential flag

!---- variables for calculation of heat flux
  logical,intent(in):: ifhfvol        ! local volume-based or local surface-based

  integer,intent(in):: nhfregion       ! number of region to calculate heat flux
  real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                ! z-position of region for heat flux

  integer,intent(in):: hftyp_atm(:)    ! atom- or mole-based heat flux cal. 
                                !   for each atom

!---- center of mass of each molecule
  real(8),intent(in):: molecom(:,:) ! center of mass of molecule

!---- variables for timestep control
  integer,intent(in):: outinterval     ! output interval of trajectory
  integer,intent(in):: pressinterval   ! interval of pressure output
  integer,intent(in):: heatfinterval   ! interval of heatf output
  integer,intent(in):: recinterval     ! state record interval
  integer,intent(in):: nstep_book      ! bookkeeping interval      
  integer,intent(in):: current_step    ! current time step

  integer,intent(in):: nstep_short     ! number of step for short force
  integer,intent(in):: istep_short     ! number of step for short force (current)

!---- MD control parameter
  integer,intent(in):: md_cont         ! MD control flag

!---- variables for extra custom NB potential
  integer,intent(in):: ncstmnbex       ! number of extra custom NB output

!     OUTPUT
  real(8),intent(inout):: force(:,:)       ! force calculated here

!---- potential
  real(8),intent(inout):: pot_tot          ! all potential
  real(8),intent(inout):: pot_nonbon       ! all non-bonded potential
  real(8),intent(inout):: pot_vdw          ! vdw potential
  real(8),intent(inout):: pot_elc          ! coulomb potential(ewald real)
  real(8),intent(inout):: pot_ewk          ! coulomb potential(ewald wave)
  real(8),intent(inout):: pot_ewnq         ! coulomb potential(ewald netq)
  real(8),intent(inout):: pot_vdw14        ! 1-4vdw potential
  real(8),intent(inout):: pot_elc14        ! 1-4elc potential
  real(8),intent(inout):: pot_bond         ! bond potential
  real(8),intent(inout):: pot_angl         ! angle potential
  real(8),intent(inout):: pot_anglub       ! Urey-Bradley angle potential
  real(8),intent(inout):: pot_tors         ! torsion potential
  real(8),intent(inout):: pot_torsrb       ! RBtorsion potential
  real(8),intent(inout):: pot_torsim       ! improper torsion potential
  real(8),intent(inout):: pot_mor          ! Morse potential
  real(8),intent(inout):: pot_sh           ! SH potential
  real(8),intent(inout):: pot_rfh          ! RFH potential
  real(8),intent(inout):: pot_dou          ! DOU potential
  real(8),intent(inout):: pot_rpvw         ! RP-VW interaction
  real(8),intent(inout):: pot_vw           ! potential of constant force
  real(8),intent(inout):: pot_cstmnb       ! Custom NB potential
  real(8),intent(inout):: pot_cstmnbex(:)  ! extra Custom NB potential
  real(8),intent(inout):: pot_posres       ! position restraint potential
  real(8),intent(inout):: pot_pbias        ! bias potential

!     for molecular pressure
  real(8),intent(inout):: for_viri_coul(:,:) ! virial(coulomb force) of each atom
  real(8),intent(inout):: pot_viri_coul       ! virial(coulomb potential)
  real(8),intent(inout):: pot_virit_coul(:,:) ! virial tensor (coulomb)

  real(8),intent(inout):: for_viri_lj(:,:)  ! virial(L-J force) of each atom
  real(8),intent(inout):: pot_viri_lj       ! virial(L-J potential)
  real(8),intent(inout):: pot_virit_lj(:,:) ! virial tensor (L-J)

  real(8),intent(inout):: for_viri_mor(:,:)  ! virial(Morse force) of each atom
  real(8),intent(inout):: pot_viri_mor       ! virial(Morse potential)
  real(8),intent(inout):: pot_virit_mor(:,:) ! virial tensor (Morse)

  real(8),intent(inout):: for_viri_sh(:,:)   ! virial(SH force) of each atom
  real(8),intent(inout):: pot_viri_sh        ! virial(SH potential)
  real(8),intent(inout):: pot_virit_sh(:,:)  ! virial tensor (SH)

  real(8),intent(inout):: for_viri_rfh(:,:)  ! virial(RFH force) of each atom
  real(8),intent(inout):: pot_viri_rfh       ! virial(RFH potential)
  real(8),intent(inout):: pot_virit_rfh(:,:) ! virial tensor (RFH)

  real(8),intent(inout):: for_viri_dou(:,:)  ! virial(DOU force) of each atom
  real(8),intent(inout):: pot_viri_dou       ! virial(DOU potential)
  real(8),intent(inout):: pot_virit_dou(:,:) ! virial tensor (DOU)

  real(8),intent(inout):: for_viri_cstmnb(:,:)  
                                        ! virial(custom NB force) of each atom
  real(8),intent(inout):: pot_viri_cstmnb       ! virial(custom NB potential)
  real(8),intent(inout):: pot_virit_cstmnb(:,:) ! virial tensor (custom NB)

  real(8),intent(inout):: for_viri_cstmnbex(:,:,:)
                                   ! extra virial(custom NB force) of each atom
  real(8),intent(inout):: pot_viri_cstmnbex(:) 
                                            ! extra virial(custom NB potential)
  real(8),intent(inout):: pot_virit_cstmnbex(:,:,:) 
                                             ! extra virial tensor (custom NB)

!     for atomic pressure
  real(8),intent(inout):: atm_viri_bond       ! virial(bond potential)
  real(8),intent(inout):: atm_virit_bond(:,:) ! virial tensor (bond potential)
  real(8),intent(inout):: atm_viri_angl       ! virial(angle potential)
  real(8),intent(inout):: atm_virit_angl(:,:) ! virial tensor (angle potential)
  real(8),intent(inout):: atm_viri_tors       ! virial(torsion potential)
  real(8),intent(inout):: atm_virit_tors(:,:) ! virial tensor (torsion potential)
  real(8),intent(inout):: atm_viri_14         ! virial(1-4 force potential)
  real(8),intent(inout):: atm_virit_14(:,:) ! virial tensor (1-4 force potential)

!     potential for each atom
  real(8),intent(inout):: pot_bo_atm(:)    ! bond potential of each atom
  real(8),intent(inout):: pot_an_atm(:)    ! angl potential of each atom
  real(8),intent(inout):: pot_to_atm(:)    ! torsion potential of each atom
  real(8),intent(inout):: pot_14_atm(:)    ! 1-4 potential of each atom
  real(8),intent(inout):: pot_elin_atm(:)  ! elc (intra) potential of each atom
  real(8),intent(inout):: pot_ljin_atm(:)  ! L-J (intra) potential of each atom

  real(8),intent(inout):: pot_elc_atm(:)    ! Coulomb potential of each atom
  real(8),intent(inout):: pot_vdw_atm(:)    ! VDW potential of each atom
  real(8),intent(inout):: pot_mor_atm(:)    ! Morse potential of each atom
  real(8),intent(inout):: pot_sh_atm(:)     ! SH potential of each atom
  real(8),intent(inout):: pot_rfh_atm(:)    ! RFH potential of each atom
  real(8),intent(inout):: pot_dou_atm(:)    ! DOU potential of each atom
  real(8),intent(inout):: pot_cstmnb_atm(:) ! custom NB potential of each atom
  real(8),intent(inout):: pot_cstmnbex_atm(:,:)
                                       ! extra custom NB potential of each atom

!     virial tensor term for each atom
  real(8),intent(inout):: viribot_atm(:,:,:,:)
                       ! virial tensor of each atom (bond)
  real(8),intent(inout):: viriant_atm(:,:,:,:)
                       ! virial tensor of each atom (angle)
  real(8),intent(inout):: viritot_atm(:,:,:,:)
                       ! virial tensor of each atom (torsion)
  real(8),intent(inout):: viri14t_atm(:,:,:,:)
                       ! virial tensor of each atom (1-4)
  real(8),intent(inout):: virielint_atm(:,:,:,:)
                                ! virial tensor of each atom (elc intra)
  real(8),intent(inout):: viriljint_atm(:,:,:,:)
                                ! virial tensor of each atom (L-J intra)

  real(8),intent(inout):: virielct_atm(:,:,:,:)
                       ! virial tensor of each atom (Coulomb)
  real(8),intent(inout):: viriljt_atm(:,:,:,:)
                       ! virial tensor of each atom (L-J)
  real(8),intent(inout):: virimort_atm(:,:,:,:)
                       ! virial tensor of each atom (Morse)
  real(8),intent(inout):: virisht_atm(:,:,:,:)
                       ! virial tensor of each atom (SH)
  real(8),intent(inout):: virirfht_atm(:,:,:,:)
                       ! virial tensor of each atom (RFH)
  real(8),intent(inout):: viridout_atm(:,:,:,:)
                       ! virial tensor of each atom (DOU)
  real(8),intent(inout):: viricstmnbt_atm(:,:,:,:)
                       ! virial tensor of each atom (custom NB)
  real(8),intent(inout):: viricstmnbext_atm(:,:,:,:,:)
                       ! virial tensor of each atom (extra custom NB)

! LOCAL:
!      real*8:: pot_tmp,pot_tmp1,pot_tmp2,pot_tmp3,pot_tmp4,pot_tmp5
                                !  temprary potential

  real(8):: pi               ! = 3.14159...
  real(8):: a1               ! = -pi/alpha**2/V

  integer:: i

!     +     +     +     +     +     +     +

!     --- INITIALIZATION ---

  force(1:3,1:natom) = 0.0d0
  for_viri_coul(1:3,1:natom) = 0.0d0
  for_viri_lj(1:3,1:natom) = 0.0d0
  for_viri_mor(1:3,1:natom) = 0.0d0
  for_viri_sh(1:3,1:natom) = 0.0d0
  for_viri_rfh(1:3,1:natom) = 0.0d0
  for_viri_dou(1:3,1:natom) = 0.0d0
  for_viri_cstmnb(1:3,1:natom) = 0.0d0
  for_viri_cstmnbex(1:3,1:natom,1:ncstmnbex) = 0.0d0

  pot_viri_coul   = 0.0d0
  pot_viri_lj     = 0.0d0
  pot_viri_mor    = 0.0d0
  pot_viri_sh     = 0.0d0
  pot_viri_rfh    = 0.0d0
  pot_viri_dou    = 0.0d0
  pot_viri_cstmnb = 0.0d0
  pot_viri_cstmnbex(1:ncstmnbex) = 0.0d0

  atm_viri_bond  = 0.0d0
  atm_viri_angl  = 0.0d0
  atm_viri_tors  = 0.0d0
  atm_viri_14    = 0.0d0

  pot_virit_coul(1:3,1:3)   = 0.0d0
  pot_virit_lj(1:3,1:3)     = 0.0d0
  pot_virit_mor(1:3,1:3)    = 0.0d0
  pot_virit_sh(1:3,1:3)     = 0.0d0
  pot_virit_rfh(1:3,1:3)    = 0.0d0
  pot_virit_dou(1:3,1:3)    = 0.0d0
  pot_virit_cstmnb(1:3,1:3) = 0.0d0
  pot_virit_cstmnbex(1:3,1:3,1:ncstmnbex) = 0.0d0
  atm_virit_bond(1:3,1:3)   = 0.0d0
  atm_virit_angl(1:3,1:3)   = 0.0d0
  atm_virit_tors(1:3,1:3)   = 0.0d0
  atm_virit_14(1:3,1:3)     = 0.0d0

  pot_interface_1 = 0.0d0                                                       !!!!!!!!! denfine interface potential   
  pot_interface_2 = 0.0d0
  pot_interface_no = 0.0d0 

  kin_interface_1 = 0.0d0
  kin_interface_2 = 0.0d0
  kin_interface_no = 0.0d0

  wall_pot_interface_1_ar1 = 0.0d0
  wall_pot_interface_1_ar2 = 0.0d0
  wall_pot_interface_2_ar1 = 0.0d0
  wall_pot_interface_2_ar2 = 0.0d0
  wall_pot_interface_no = 0.0d0

  if (((mod(current_step,heatfinterval) == 1) .or. &
       & (heatfinterval == 1)) .and. &
       & (istep_short == nstep_short)) then ! only when heat flux calculation
     pot_bo_atm(1:natom) = 0.0d0
     pot_an_atm(1:natom) = 0.0d0
     pot_to_atm(1:natom) = 0.0d0
     pot_14_atm(1:natom) = 0.0d0
     pot_elin_atm(1:natom) = 0.0d0
     pot_ljin_atm(1:natom) = 0.0d0
     pot_elc_atm(1:natom) = 0.0d0
     pot_vdw_atm(1:natom) = 0.0d0
     pot_mor_atm(1:natom) = 0.0d0
     pot_sh_atm(1:natom)  = 0.0d0
     pot_rfh_atm(1:natom) = 0.0d0
     pot_dou_atm(1:natom) = 0.0d0
     pot_cstmnb_atm(1:natom) = 0.0d0
     pot_cstmnbex_atm(1:natom,1:ncstmnbex) = 0.0d0

     viribot_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     viriant_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     viritot_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     viri14t_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     virielint_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     viriljint_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     virielct_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     viriljt_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     virimort_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     virisht_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     virirfht_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     viridout_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     viricstmnbt_atm(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
     viricstmnbext_atm(1:3,1:3,1:natom,1:nhfregion,1:ncstmnbex) = 0.0d0


    pot_vdw_atm_ar11(1:natom) = 0.0d0                                         !!!!!!!!! denfine ar1-ar2 potential energy
    pot_vdw_atm_ar12(1:natom) = 0.0d0
    pot_vdw_atm_ar22(1:natom) = 0.0d0

    viriljt_atm_ar11(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
    viriljt_atm_ar12(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
    viriljt_atm_ar22(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
    viriljt_atm_ar1pt(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
    viriljt_atm_ar2pt(1:3,1:3,1:natom,1:nhfregion) = 0.0d0
  end if

!     --- BOND ---


  IF (mts_bond == mts_flag) THEN

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier( bond_time, bond_time_elps, 'START' )
#endif
!!! Time measure

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

        call calbondp_hf( xcel, ycel, zcel, &
             &            force,pot_bond, &
             &            atm_viri_bond,atm_virit_bond, &
             &            pot_bo_atm,viribot_atm, &
             &            ifhfvol, &
             &            nhfregion,hfzpos1,hfzpos2, &
             &            hftyp_atm, &
             &            molecom)
            
     else 

        call calbondp( xcel, ycel, zcel, &
             &         force, pot_bond, &
             &         atm_viri_bond, atm_virit_bond )

     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount(bond_time,bond_time_elps,'STOP')
#endif
!!! Time measure
         
  END IF


!     --- ANGLE ---


  IF (mts_angl == mts_flag) THEN

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(bond_time,bond_time_elps,'START')
#endif
!!! Time measure

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

        call calanglp_hf( xcel, ycel, zcel, &
             &            force, pot_angl, &
             &            atm_viri_angl, atm_virit_angl, &
             &            pot_an_atm, viriant_atm, &
             &            ifhfvol, &
             &            nhfregion, hfzpos1, hfzpos2, &
             &            hftyp_atm, &
             &            molecom )

     else

        call calanglp( xcel, ycel, zcel, &
             &         force, pot_angl, &
             &         atm_viri_angl, atm_virit_angl )

     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount(bond_time,bond_time_elps,'STOP')
#endif
!!! Time measure

  END IF


!     --- ANGLE (Urey-Bradley type) ---


  IF (mts_anglub == mts_flag) THEN

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(bond_time,bond_time_elps,'START')
#endif
!!! Time measure

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

        call calanglubp_hf(xcel,ycel,zcel, &
             &             force,pot_anglub, &
             &             atm_viri_angl,atm_virit_angl, &
             &             pot_an_atm,viriant_atm, &
             &             ifhfvol, &
             &             nhfregion,hfzpos1,hfzpos2, &
             &             hftyp_atm, &
             &             molecom)

     else

        call calanglubp(xcel,ycel,zcel, &
             &          force,pot_anglub, &
             &          atm_viri_angl,atm_virit_angl)

     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount(bond_time,bond_time_elps,'STOP')
#endif
!!! Time measure

  END IF


!     --- TORSION (periodic type) ---


  IF (mts_tors == mts_flag) THEN

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier( bond_time, bond_time_elps, 'START' )
#endif
!!! Time measure

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

        call caltorsp_hf( xcel, ycel, zcel, &
             &            force,pot_tors, &
             &            atm_viri_tors,atm_virit_tors, &
             &            pot_to_atm,viritot_atm, &
             &            ifhfvol, &
             &            nhfregion,hfzpos1,hfzpos2, &
             &            hftyp_atm, &
             &            molecom)

     else

        call caltorsp( xcel, ycel, zcel, &
             &         force, pot_tors, &
             &         atm_viri_tors,atm_virit_tors)

     end if

!           call caltors(nimp_tot,impro,natom_tot,xyz,force,pot_tmp3)

!           call caltrsec(nsec, trssec, ntrsh_tot, ntrso_tot, trso, trsh, &
!           &             natom_tot, xyz, force,                        &
!           &             pot_tmp4, pot_tmp5)

!           pot_impro     = pot_tmp3

!!! Time measure
#if defined(TIME_M)
     call get_timecount(bond_time,bond_time_elps,'STOP')
#endif
!!! Time measure

  END IF


!     --- TORSION_RB (Ryckaert and Bellmans type) ---


  IF (mts_torsrb == mts_flag) THEN

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(bond_time,bond_time_elps,'START')
#endif
!!! Time measure

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

        call caltorsrbp_hf( xcel, ycel, zcel, &
             &              force,pot_torsrb, &
             &              atm_viri_tors,atm_virit_tors, &
             &              pot_to_atm,viritot_atm, &
             &              ifhfvol, &
             &              nhfregion,hfzpos1,hfzpos2, &
             &              hftyp_atm, &
             &              molecom )

     else

        call caltorsrbp( xcel, ycel, zcel, &
             &           force,pot_torsrb, &
             &           atm_viri_tors,atm_virit_tors)

     end if

!           call caltors(nimp_tot,impro,natom_tot,xyz,force,pot_tmp3)

!           call caltrsec(nsec, trssec, ntrsh_tot, ntrso_tot, trso, trsh, &
!           &             natom_tot, xyz, force,                        &
!           &             pot_tmp4, pot_tmp5)

!           pot_impro     = pot_tmp3

!!! Time measure
#if defined(TIME_M)
     call get_timecount(bond_time,bond_time_elps,'STOP')
#endif
!!! Time measure

  END IF


!     --- TORSION_IM ---


  IF (mts_torsim == mts_flag) THEN

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(bond_time,bond_time_elps,'START')
#endif
!!! Time measure

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

        call caltorsimp_hf( xcel,ycel,zcel, &
             &              force,pot_torsim, &
             &              atm_viri_tors, atm_virit_tors, &
             &              pot_to_atm, viritot_atm, &
             &              ifhfvol, &
             &              nhfregion, hfzpos1, hfzpos2, &
             &              hftyp_atm, &
             &              molecom )

     else

        call caltorsimp( xcel, ycel, zcel, &
             &           force, pot_torsim, &
             &           atm_viri_tors, atm_virit_tors )

     end if

!           call caltors(nimp_tot,impro,natom_tot,xyz,force,pot_tmp3)

!           call caltrsec(nsec, trssec, ntrsh_tot, ntrso_tot, trso, trsh, &
!           &             natom_tot, xyz, force,                        &
!           &             pot_tmp4, pot_tmp5)

!           pot_impro     = pot_tmp3

!!! Time measure
#if defined(TIME_M)
     call get_timecount( bond_time, bond_time_elps, 'STOP' )
#endif
!!! Time measure

  END IF


!     --- Morse ---


  IF (mts_mor == mts_flag) THEN 

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier( mor_time, mor_time_elps, 'START' )
#endif
!!! Time measure

!        - Morse interaction -

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

!           --- morsep for heat flux calculation ---
        call morsep_hf( xcel, ycel, zcel, &
             &          rcutmor, &
             &          force,pot_mor, &
             &          for_viri_mor, pot_viri_mor, &
             &          pot_virit_mor, &
             &          pot_mor_atm, virimort_atm, &
             &          ifhfvol, &
             &          nhfregion, hfzpos1, hfzpos2, &
             &          hftyp_atm, &
             &          molecom )

     else

        call morsep( xcel, ycel, zcel, &
             &       rcutmor, &
             &       force,pot_mor, &
             &       for_viri_mor,pot_viri_mor, &
             &       pot_virit_mor)

     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount(mor_time,mor_time_elps,'STOP')
#endif
!!! Time measure

  END IF
      

!     --- SH ---


  IF (mts_sh == mts_flag) THEN 

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(sh_time,sh_time_elps,'START')
#endif
!!! Time measure

!        - SH interaction -

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

!           --- shpotp for heat flux calculation ---
        call shpotp_hf( xcel, ycel, zcel, &
             &          rcutsh, &
             &          force, pot_sh, &
             &          for_viri_sh, pot_viri_sh, &
             &          pot_virit_sh, &
             &          pot_sh_atm, virisht_atm, &
             &          ifhfvol, &
             &          nhfregion, hfzpos1, hfzpos2, &
             &          hftyp_atm, &
             &          molecom )

     else

        call shpotp( xcel, ycel, zcel, &
             &       rcutsh, &
             &       force, pot_sh, &
             &       for_viri_sh, pot_viri_sh, &
             &       pot_virit_sh )

     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount(sh_time,sh_time_elps,'STOP')
#endif
!!! Time measure

  END IF
      

!     --- RFH ---


  IF (mts_rfh == mts_flag) THEN 

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(rfh_time,rfh_time_elps,'START')
#endif
!!! Time measure

!        - RFH interaction -

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

!           --- shpotp for heat flux calculation ---
        call rfhpotp_hf( xcel, ycel, zcel, &
             &           rcutrfhfo,rcutrfhoo,rcutrfhoh, &
             &           force,pot_rfh, &
             &           for_viri_rfh,pot_viri_rfh, &
             &           pot_virit_rfh, &
             &           pot_rfh_atm,virirfht_atm, &
             &           ifhfvol, &
             &           nhfregion,hfzpos1,hfzpos2, &
             &           hftyp_atm, &
             &           molecom )

     else

        call rfhpotp( xcel, ycel, zcel, &
             &        rcutrfhfo,rcutrfhoo,rcutrfhoh, &
             &        force,pot_rfh, &
             &        for_viri_rfh,pot_viri_rfh, &
             &        pot_virit_rfh )

     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount(rfh_time,rfh_time_elps,'STOP')
#endif
!!! Time measure

  END IF
      

!     --- DOU ---


  IF (mts_dou == mts_flag) THEN 

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(dou_time,dou_time_elps,'START')
#endif
!!! Time measure

!        - DOU interaction -

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

        call doupotp_hf( xcel, ycel, zcel, &
             &           rcutdouo, rcutindouo, rcutdouh, &
             &           force, pot_dou, &
             &           for_viri_dou, pot_viri_dou, &
             &           pot_virit_dou, &
             &           pot_dou_atm, viridout_atm, &
             &           ifhfvol, &
             &           nhfregion, hfzpos1, hfzpos2, &
             &           hftyp_atm, &
             &           molecom )

     else

        call doupotp( xcel, ycel, zcel, &
             &        rcutdouo,rcutindouo,rcutdouh, &
             &        force,pot_dou, &
             &        for_viri_dou,pot_viri_dou, &
             &        pot_virit_dou )

     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount(dou_time,dou_time_elps,'STOP')
#endif
!!! Time measure

  END IF


!     --- CSTMNB ---

!!! Time measure
#if defined(TIME_M)
  call get_timecount_barrier(cstmnb_time,cstmnb_time_elps,'START')
#endif
!!! Time measure

!    - custom NB interaction -

  if (ifcstmnb) then

     ! entry to calcstmnb_hf
     call ent_calcstmnb_hf(mts_flag,mts_cstmnb, &
          &                xcel,ycel,zcel, &
          &                npoly,npolytyp,npoly_mole,npoly_atom, &
          &                nwater,nmatom,nmatyp,nmatomtyp, &
          &                force,pot_cstmnb, &
          &                for_viri_cstmnb,pot_viri_cstmnb, &
          &                pot_virit_cstmnb, &
          &                pot_cstmnb_atm,viricstmnbt_atm, &
          &                ifhfvol,   &
          &                nhfregion,hfzpos1,hfzpos2,   &
          &                hftyp_atm,   &
          &                molecom, &
          &                heatfinterval, &
          &                ncstmnbex, &
          &                pot_cstmnbex, &
          &                for_viri_cstmnbex,pot_viri_cstmnbex, &
          &                pot_virit_cstmnbex, &
          &                pot_cstmnbex_atm,viricstmnbext_atm, &
          &                ifcalpreatom,ifcalpremole, &
          &                current_step, &
          &                nstep_short,istep_short)

  end if

!!! Time measure
#if defined(TIME_M)
  call get_timecount(cstmnb_time,cstmnb_time_elps,'STOP')
#endif
!!! Time measure


!     --- NONBONDED INTERACTIONS (ELC & VDW) ---


  IF (mts_vdw == mts_flag) THEN 

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(drct_time,drct_time_elps,'START')
#endif
!!! Time measure

!    - Ewald real and/or VDW -

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

!       --- ewaldrp for heat flux calculation ---
        call ewaldrp_hf( xcel,ycel,zcel, &
             &           alpha,rrcut, &
             &           rcut, &
             &           force,pot_elc,pot_vdw, &
             &           for_viri_coul,pot_viri_coul, &
             &           for_viri_lj,pot_viri_lj, &
             &           pot_virit_coul,pot_virit_lj, &
             &           pot_elin_atm,virielint_atm, &
             &           pot_elc_atm,virielct_atm, &
             &           pot_ljin_atm,viriljint_atm, &
             &           pot_vdw_atm,viriljt_atm, &
             &           ifhfvol, &
             &           nhfregion,hfzpos1,hfzpos2, &
             &           hftyp_atm, &
             &           molecom )

     else

        call ewaldrp( xcel, ycel, zcel, &
             &        alpha,rrcut, &
             &        rcut, &
             &        force,pot_elc,pot_vdw, &
             &        for_viri_coul,pot_viri_coul, &
             &        for_viri_lj,pot_viri_lj, &
             &        pot_virit_coul,pot_virit_lj)

     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount(drct_time,drct_time_elps,'STOP')
#endif
!!! Time measure

  END IF

      
  IF (mts_ewk == mts_flag) THEN

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(rcpr_time,rcpr_time_elps,'START')
#endif
!!! Time measure

!        - Ewald K-space -

!--- If use -D_HF_BULK_EWK, force to use standard Ewald summation method
!---  at a timestep of measuring heat flux.
#if defined(_HF_BULK_EWK)

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

!       --- ewaldkp for heat flux calculation ---
        call ewaldkp_hf(xcel,ycel,zcel,   &
             &          alpha,   &
             &          force,pot_ewk,   &
             &          ifnetqcorrp,   &
             &          netchrgsq,   &
             &          for_viri_coul,pot_viri_coul,pot_virit_coul, &
             &          pot_elc_atm,virielct_atm, &
             &          ifhfvol, &
             &          nhfregion,hfzpos1,hfzpos2, &
             &          hftyp_atm, &
             &          molecom)

!-   At other timesteps, normal procedure is performed.
     else
#endif

!        --- standard Ewald method
     if (ifewald) then

        if (ifcalpremole .or. ifcalpreatom) then

           call ewaldkp( xcel, ycel, zcel, &
                &        alpha, &
                &        force,pot_ewk, &
                &        ifnetqcorrp, &
                &        netchrgsq, &
                &        for_viri_coul, pot_viri_coul, pot_virit_coul )

        else

           call ewaldk( xcel, ycel, zcel, &
                &       force, pot_ewk )

        end if

!        --- SPME method
     else if (ifspme) then

        if (ifcalpremole .or. ifcalpreatom) then

           call fft_pmep( xcel, ycel, zcel, &
                &         alpha, &
                &         nfft1, nfft2, nfft3, pme_order, &
                &         force, pot_ewk, &
                &         ifnetqcorrp, &
                &         netchrgsq, &
                &        for_viri_coul, pot_viri_coul, pot_virit_coul)

        else

           call fft_pme( xcel, ycel, zcel, &
                &        alpha, &
                &        nfft1,nfft2,nfft3,pme_order, &
                &        ifnetqcorrp, &
                &        netchrgsq, &
                &        force,pot_ewk )

        end if

     end if

!--- If use -D_HF_BULK_EWK, force to use standard Ewald summation method.
#if defined(_HF_BULK_EWK)
!-   At other timesteps, normal procedure is performed.
     end if
#endif

!!! Time measure
#if defined(TIME_M)
     call get_timecount(rcpr_time,rcpr_time_elps,'STOP')
#endif
!!! Time measure

  END IF


  IF ((mts_vdw14 == mts_flag) .and. &
       & (mts_elc14 == mts_flag)) THEN 

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(bond_time,bond_time_elps,'START')
#endif
!!! Time measure

!        - 14 VDW & 14 ELC - 

     if (((mod(current_step,heatfinterval) == 1) .or. &
          & (heatfinterval == 1)) .and. &
          & (istep_short == nstep_short)) then

        call calnon14p_hf( xcel, ycel, zcel, &
             &             div_factor_14vdw,div_factor_14elc, &
             &             force,pot_elc14,pot_vdw14, &
             &             atm_viri_14,atm_virit_14, &
             &             pot_14_atm,viri14t_atm, &
             &             ifhfvol, &
             &             nhfregion,hfzpos1,hfzpos2, &
             &             hftyp_atm, &
             &             molecom )

     else

        call calnon14p( xcel, ycel, zcel, &
             &          div_factor_14vdw, div_factor_14elc, &
             &          force, pot_elc14,pot_vdw14, &
             &          atm_viri_14, atm_virit_14 )

     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount(bond_time,bond_time_elps,'STOP')
#endif
!!! Time measure

  END IF


!     --- CNP_VW ---


  IF (mts_cnpvw == mts_flag) THEN

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(cnpvw_time,cnpvw_time_elps,'START')
#endif
!!! Time measure

!        - RP-VW interaction & constant forse -

!     if (((mod(current_step,heatfinterval) == 1) .or. &
!          & (heatfinterval == 1)) .and. &
!          & (istep_short == nstep_short)) then
!!! do not use for pressure calculation and heat flux calculation
!     else

        call cnpvw( rcutrpvw, force, pot_rpvw, pot_vw )

!     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount(cnpvw_time,cnpvw_time_elps,'STOP')
#endif
!!! Time measure

  END IF


!     --- POSITION RESTRAINT ---


  IF (mts_posres == mts_flag) THEN

     call calposres( xcel, ycel, zcel, &
          &          force,pot_posres)

  END IF


!     --- BIAS POTENTIAL ---


  IF (mts_potbias == mts_flag) THEN

     if (ifpotbias) then
        call calpotbias( xcel, ycel, zcel, &
             &           npoly,nwater,nmatom, &
             &           force,pot_pbias)
     end if

  END IF


!---- MPI Sum all potential and force

#if defined(MPI)

!!! Time measure
#if defined(TIME_M)
  call get_timecount_barrier(comm_time,comm_time_elps,'START')
#endif
!!! Time measure

  call MPIsum_all_pot_for(mts_flag, &
       &                  mts_bond,mts_angl,mts_anglub, &
       &                  mts_tors,mts_torsrb,mts_torsim, &
       &                  mts_vdw,mts_ewr,mts_ewk, &
       &                  mts_vdw14,mts_elc14, &
       &                  mts_mor,mts_sh,mts_rfh,mts_dou,mts_cnpvw, &
       &                  mts_cstmnb, &
       &                  mts_posres, &
       &                  mts_potbias, &
       &                  force, &
       &                  pot_vdw,pot_elc,pot_ewk, &
       &                  pot_vdw14,pot_elc14, &
       &                  pot_bond,pot_angl,pot_anglub, &
       &                  pot_tors,pot_torsrb,pot_torsim, &
       &                  pot_mor, &
       &                  pot_sh, &
       &                  pot_rfh, &
       &                  pot_dou, &
       &                  pot_rpvw,pot_vw, &
       &                  pot_cstmnb, &
       &                  ncstmnbex,pot_cstmnbex, &
       &                  pot_posres, &
       &                  pot_pbias)

!!! Time measure
#if defined(TIME_M)
  call get_timecount(comm_time,comm_time_elps,'STOP')
#endif
!!! Time measure
      
#endif

!---- MPI Sum for pressure calculation

#if defined(MPI)

!!! Time measure
#if defined(TIME_M)
  call get_timecount_barrier(comm_time,comm_time_elps,'START')
#endif
!!! Time measure

  if (ifcalpremole .or. ifcalpreatom) then

     if ((mod(current_step,pressinterval) == 1) .or. &
          & (pressinterval == 1)) then

        if (istep_short == nstep_short) then

           call MPIsum_all_viri(mts_flag, &
                &               mts_bond,mts_angl,mts_anglub, &
                &               mts_tors,mts_torsrb,mts_torsim, &
                &               mts_vdw,mts_ewr,mts_ewk, &
                &               mts_vdw14,mts_elc14, &
                &               mts_mor,mts_sh,mts_rfh,mts_dou, &
                &               mts_cstmnb, &
                &               for_viri_coul,pot_viri_coul, &
                &               for_viri_lj,pot_viri_lj, &
                &               for_viri_mor,pot_viri_mor, &
                &               for_viri_sh,pot_viri_sh, &
                &               for_viri_rfh,pot_viri_rfh, &
                &               for_viri_dou,pot_viri_dou, &
                &               for_viri_cstmnb,pot_viri_cstmnb, &
                &               for_viri_cstmnbex,pot_viri_cstmnbex, &
                &               atm_viri_bond,atm_viri_angl, &
                &               atm_viri_tors,atm_viri_14, &
                &               pot_virit_coul,pot_virit_lj, &
                &               pot_virit_mor, &
                &               pot_virit_sh, &
                &               pot_virit_rfh, &
                &               pot_virit_dou, &
                &               pot_virit_cstmnb, &
                &               pot_virit_cstmnbex, &
                &               atm_virit_bond,atm_virit_angl, &
                &               atm_virit_tors,atm_virit_14, &
                &               ifcalpremole,ifcalpreatom, &
                &               ncstmnbex)

        end if

     end if

  end if

!!! Time measure
#if defined(TIME_M)
  call get_timecount(comm_time,comm_time_elps,'STOP')
#endif
!!! Time measure
      
#endif

!---- MPI Sum for heatflux calculation

#if defined(MPI)

!!! Time measure
#if defined(TIME_M)
  call get_timecount_barrier(comm_time,comm_time_elps,'START')
#endif
!!! Time measure
      
  if ((mod(current_step,heatfinterval) == 1) .or. &
       & (heatfinterval == 1)) then

     if (istep_short == nstep_short) then
        call MPIsum_all_atm_heatf(mts_flag, &
             &                    mts_bond,mts_angl,mts_anglub, &
             &                    mts_tors,mts_torsrb,mts_torsim, &
             &                    mts_vdw,mts_ewr,mts_ewk, &
             &                    mts_vdw14,mts_elc14, &
             &                    mts_mor,mts_sh,mts_rfh,mts_dou, &
             &                    mts_cstmnb, &
             &                    nhfregion, &
             &                    pot_bo_atm, &
             &                    pot_an_atm, &
             &                    pot_to_atm, &
             &                    pot_14_atm, &
             &                    pot_elin_atm,pot_ljin_atm, &
             &                    pot_elc_atm, &
             &                    pot_vdw_atm, &
             &                    pot_mor_atm, &
             &                    pot_sh_atm, &
             &                    pot_rfh_atm, &
             &                    pot_dou_atm, &
             &                    pot_cstmnb_atm, &
             &                    pot_cstmnbex_atm, &
             &                    viribot_atm, &
             &                    viriant_atm, &
             &                    viritot_atm, &
             &                    viri14t_atm, &
             &                    virielint_atm,viriljint_atm, &
             &                    virielct_atm, &
             &                    viriljt_atm, &
             &                    virimort_atm, &
             &                    virisht_atm, &
             &                    virirfht_atm, &
             &                    viridout_atm, &
             &                    viricstmnbt_atm, &
             &                    viricstmnbext_atm, &
             &                    ncstmnbex)

     end if

  end if

!!! Time measure
#if defined(TIME_M)
  call get_timecount(comm_time,comm_time_elps,'STOP')
#endif
!!! Time measure

#endif

!     --- calculate potential and virial of net charge MD

  IF (mts_flag == MTS_LONG) THEN   ! = if mts_flag = mts_long(=1)

     if (ifewald .or. ifspme) then

        pi = dacos(-1.0d0)
        a1 = -pi/(alpha*alpha*xcel*ycel*zcel)
        pot_ewnq = 0.5d0*a1*netchrgsq ! = +A1*netchrgsq/2

        if (.not. ifnetqcorrp) then ! no charge correction

           pot_viri_coul = pot_viri_coul + a1*netchrgsq
           pot_virit_coul(1,1) =  pot_virit_coul(1,1) &
                &             + 0.5d0*a1*netchrgsq
           pot_virit_coul(2,2) =  pot_virit_coul(2,2) &
                &             + 0.5d0*a1*netchrgsq
           pot_virit_coul(3,3) =  pot_virit_coul(3,3) &
                &             + 0.5d0*a1*netchrgsq

        else                        ! charge correction

           pot_viri_coul = pot_viri_coul - 0.5d0*a1*netchrgsq

        end if

     end if

  END IF

!     --- ADD UP ALL THE POTENTIAL ENERGIES ---


  pot_nonbon  =   pot_vdw &
       &        + pot_elc   + pot_ewk &
       &        + pot_ewc   + pot_ewnq &
       &        + pot_vdw14 + pot_elc14 &
       &        + pot_sh    + pot_rfh &
       &        + pot_dou &
       &        + pot_rpvw  + pot_vw &
       &        + pot_cstmnb &
       &        + pot_posres + pot_pbias

#if defined(_CSTMNB_V2_ADD_ALL)
  !- adding cstmnb extra potential
  do i = 1, ncstmnbex
     pot_nonbon = pot_nonbon + pot_cstmnbex(i)
  end do
#endif

  pot_tot     =   pot_bond + pot_angl + pot_anglub &
       &        + pot_tors + pot_torsrb + pot_torsim &
       &        + pot_mor &
       &        + pot_nonbon


!     +     +     +     +     +     +     +

end subroutine calforce_hf
