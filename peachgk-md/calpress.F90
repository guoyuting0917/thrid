!*****************************
!*  calpress.f90 Ver.3.0     *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-07-17 12:48:06 gota>

subroutine calpress(for_viric_long,for_viric_med,for_viric_short,   &
     &              pot_viric_long,pot_viric_med,pot_viric_short,   &
     &              pot_virict_long,pot_virict_med,   &
     &              pot_virict_short,   &
     &              for_virilj_long,for_virilj_med,   &
     &              for_virilj_short,   &
     &              pot_virilj_long,pot_virilj_med,   &
     &              pot_virilj_short,   &
     &              pot_viriljt_long,pot_viriljt_med,   &
     &              pot_viriljt_short,   &
     &              for_virimor_long,for_virimor_med,   &
     &              for_virimor_short,   &
     &              pot_virimor_long,pot_virimor_med,   &
     &              pot_virimor_short,   &
     &              pot_virimort_long,pot_virimort_med,   &
     &              pot_virimort_short,   &
     &              for_virish_long,for_virish_med,   &
     &              for_virish_short,   &
     &              pot_virish_long,pot_virish_med,   &
     &              pot_virish_short,   &
     &              pot_virisht_long,pot_virisht_med,   &
     &              pot_virisht_short,   &
     &              for_virirfh_long,for_virirfh_med,   &
     &              for_virirfh_short,   &
     &              pot_virirfh_long,pot_virirfh_med,   &
     &              pot_virirfh_short,   &
     &              pot_virirfht_long,pot_virirfht_med,   &
     &              pot_virirfht_short,   &
     &              for_viridou_long,for_viridou_med,   &
     &              for_viridou_short,   &
     &              pot_viridou_long,pot_viridou_med,   &
     &              pot_viridou_short,   &
     &              pot_viridout_long,pot_viridout_med,   &
     &              pot_viridout_short,   &
     &              for_viricstmnb_long,for_viricstmnb_med,   &
     &              for_viricstmnb_short,   &
     &              pot_viricstmnb_long,pot_viricstmnb_med,   &
     &              pot_viricstmnb_short,   &
     &              pot_viricstmnbt_long,pot_viricstmnbt_med,   &
     &              pot_viricstmnbt_short,   &
     &              for_viricstmnbex_long,for_viricstmnbex_med,   &
     &              for_viricstmnbex_short,   &
     &              pot_viricstmnbex_long,pot_viricstmnbex_med,   &
     &              pot_viricstmnbex_short,   &
     &              pot_viricstmnbext_long,pot_viricstmnbext_med,   &
     &              pot_viricstmnbext_short,   &
     &              atm_viribo_long,atm_viribo_med,   &
     &              atm_viribo_short,   &
     &              atm_viribot_long,atm_viribot_med,   &
     &              atm_viribot_short,   &
     &              atm_virian_long,atm_virian_med,   &
     &              atm_virian_short,   &
     &              atm_viriant_long,atm_viriant_med,   &
     &              atm_viriant_short,   &
     &              atm_virito_long,atm_virito_med,   &
     &              atm_virito_short,   &
     &              atm_viritot_long,atm_viritot_med,   &
     &              atm_viritot_short,   &
     &              atm_viri14_long,atm_viri14_med,   &
     &              atm_viri14_short,   &
     &              atm_viri14t_long,atm_viri14t_med,   &
     &              atm_viri14t_short,   &
     &              atm_viri_const,atm_virit_const,   &
     &              atm_viri_corr,atm_virit_corr,   &
     &              pressmol_ktot,pressmol_vtot,pressmol_tot,   &
     &              pressatm_ktot,pressatm_vtot,pressatm_tot,   &
     &              pot_viric_all,pot_virilj_all,   &
     &              pot_virimor_all,   &
     &              pot_virish_all,   &
     &              pot_virirfh_all,   &
     &              pot_viridou_all,   &
     &              pot_viricstmnb_all,   &
     &              pot_viricstmnbex_all,   &
     &              pot_viri_all,   &
     &              viri_fdotd,   &
     &              atm_viribo_all,atm_virian_all,   &
     &              atm_virito_all,atm_viri14_all,   &
     &              pressmolt_ktot,pressmolt_vtot,   &
     &              pressmolt_tot,   &
     &              pressatmt_ktot,pressatmt_vtot,   &
     &              pressatmt_tot,   &
     &              pot_virict_all,pot_viriljt_all,   &
     &              pot_virimort_all,   &
     &              pot_virisht_all,   &
     &              pot_virirfht_all,   &
     &              pot_viridout_all,   &
     &              pot_viricstmnbt_all,   &
     &              pot_viricstmnbext_all,   &
     &              pot_virit_all,   &
     &              virit_fdotd,   &
     &              atm_viribot_all,atm_viriant_all,   &
     &              atm_viritot_all,atm_viri14t_all,   &
     &              current_step,pressinterval,   &
     &              xcel,ycel,zcel,   &
     &              ifcalpremole,ifcalpreatom,   &
     &              pot_ewc,pref,eref,   &
     &              rcut,   &
     &              ifcalljlong,nsolve,solveindex,   &
     &              pot_ewnq, &
     &              ncstmnbex)

  use md_global

  implicit none

!
!     Summing pressure term of each atom:
!           PV = 1/3 * [ PK + PU]
!        molecular pressure is written by:
!           PmolV = 1/3 * [ sigma(mol){Mmol*vg.vg} +
!                           sigma(mol){rg.Fmol} ]
!
! ARGUMENTS:
!     INPUT
  real(8),intent(in):: for_viric_long(:,:)   ! long-range virial (coulomb force)
  real(8),intent(in):: for_viric_med(:,:)    ! medium-range virial (coulomb force)
  real(8),intent(in):: for_viric_short(:,:)  ! short-range virial (coulomb force)
  real(8),intent(in):: pot_viric_long        ! long-range virial (coulomb pot)
  real(8),intent(in):: pot_viric_med         ! medium-range virial (coulomb pot)
  real(8),intent(in):: pot_viric_short       ! short-range virial (coulomb pot)
  real(8),intent(in):: pot_virict_long(:,:)  ! long-range virial tensor (coulomb)
  real(8),intent(in):: pot_virict_med(:,:)   ! med-range virial tensor (coulomb)
  real(8),intent(in):: pot_virict_short(:,:)
                    ! short-range virial tensor (coulomb)

  real(8),intent(in):: for_virilj_long(:,:)   ! long-range virial (L-J force)
  real(8),intent(in):: for_virilj_med(:,:)    ! medium-range virial (L-J force)
  real(8),intent(in):: for_virilj_short(:,:)  ! short-range virial (L-J force)
  real(8),intent(in):: pot_virilj_long        ! long-range virial (L-J pot)
  real(8),intent(in):: pot_virilj_med         ! medium-range virial (L-J pot)
  real(8),intent(in):: pot_virilj_short       ! short-range virial (L-J pot)
  real(8),intent(in):: pot_viriljt_long(:,:)  ! long-range virial tensor (L-J)
  real(8),intent(in):: pot_viriljt_med(:,:)   ! medium-range virial tensor (L-J)
  real(8),intent(in):: pot_viriljt_short(:,:) ! short-range virial tensor (L-J)

  real(8),intent(in):: for_virimor_long(:,:)   ! long-range virial (Morse force)
  real(8),intent(in):: for_virimor_med(:,:)    ! med-range virial (Morse force)
  real(8),intent(in):: for_virimor_short(:,:)  ! short-range virial (Morse force)
  real(8),intent(in):: pot_virimor_long        ! long-range virial (Morse pot)
  real(8),intent(in):: pot_virimor_med         ! med-range virial (Morse pot)
  real(8),intent(in):: pot_virimor_short       ! short-range virial (Morse pot)
  real(8),intent(in):: pot_virimort_long(:,:)  ! long-range virial tensor (Morse)
  real(8),intent(in):: pot_virimort_med(:,:)   ! med-range virial tensor (Morse)
  real(8),intent(in):: pot_virimort_short(:,:) ! short-range virial tensor(Morse)

  real(8),intent(in):: for_virish_long(:,:)   ! long-range virial (SH force)
  real(8),intent(in):: for_virish_med(:,:)    ! med-range virial (SH force)
  real(8),intent(in):: for_virish_short(:,:)  ! short-range virial (SH force)
  real(8),intent(in):: pot_virish_long        ! long-range virial (SH pot)
  real(8),intent(in):: pot_virish_med         ! med-range virial (SH pot)
  real(8),intent(in):: pot_virish_short       ! short-range virial (SH pot)
  real(8),intent(in):: pot_virisht_long(:,:)  ! long-range virial tensor (SH)
  real(8),intent(in):: pot_virisht_med(:,:)   ! med-range virial tensor (SH)
  real(8),intent(in):: pot_virisht_short(:,:) ! short-range virial tensor(SH)

  real(8),intent(in):: for_virirfh_long(:,:)   ! long-range virial (RFH force)
  real(8),intent(in):: for_virirfh_med(:,:)    ! med-range virial (RFH force)
  real(8),intent(in):: for_virirfh_short(:,:)  ! short-range virial (RFH force)
  real(8),intent(in):: pot_virirfh_long        ! long-range virial (RFH pot)
  real(8),intent(in):: pot_virirfh_med         ! med-range virial (RFH pot)
  real(8),intent(in):: pot_virirfh_short       ! short-range virial (RFH pot)
  real(8),intent(in):: pot_virirfht_long(:,:)  ! long-range virial tensor (RFH)
  real(8),intent(in):: pot_virirfht_med(:,:)   ! med-range virial tensor (RFH)
  real(8),intent(in):: pot_virirfht_short(:,:) ! short-range virial tensor(RFH)

  real(8),intent(in):: for_viridou_long(:,:)   ! long-range virial (DOU force)
  real(8),intent(in):: for_viridou_med(:,:)    ! med-range virial (DOU force)
  real(8),intent(in):: for_viridou_short(:,:)  ! short-range virial (DOU force)
  real(8),intent(in):: pot_viridou_long        ! long-range virial (DOU pot)
  real(8),intent(in):: pot_viridou_med         ! med-range virial (DOU pot)
  real(8),intent(in):: pot_viridou_short       ! short-range virial (DOU pot)
  real(8),intent(in):: pot_viridout_long(:,:)  ! long-range virial tensor (DOU)
  real(8),intent(in):: pot_viridout_med(:,:)   ! med-range virial tensor (DOU)
  real(8),intent(in):: pot_viridout_short(:,:) ! short-range virial tensor(DOU)

  real(8),intent(in):: for_viricstmnb_long(:,:)
                                      ! long-range virial (custom NB force)
  real(8),intent(in):: for_viricstmnb_med(:,:)
                                      ! med-range virial (custom NB force)
  real(8),intent(in):: for_viricstmnb_short(:,:)
                                      ! short-range virial (custom NB force)
  real(8),intent(in):: pot_viricstmnb_long  ! long-range virial (custom NB pot)
  real(8),intent(in):: pot_viricstmnb_med   ! med-range virial (custom NB pot)
  real(8),intent(in):: pot_viricstmnb_short ! short-range virial (custom NB pot)
  real(8),intent(in):: pot_viricstmnbt_long(3,3) 
                                      ! long-range virial tensor (custom NB)
  real(8),intent(in):: pot_viricstmnbt_med(3,3)
                                      ! med-range virial tensor (custom NB)
  real(8),intent(in):: pot_viricstmnbt_short(3,3)
                                      ! short-range virial tensor(custom nB)

  real(8),intent(in):: for_viricstmnbex_long(:,:,:)
                                          ! long-range virial (custom NB force)
  real(8),intent(in):: for_viricstmnbex_med(:,:,:)
                                          ! med-range virial (custom NB force)
  real(8),intent(in):: for_viricstmnbex_short(:,:,:)
                                          ! short-range virial (custom NB force)
  real(8),intent(in):: pot_viricstmnbex_long(:)      
                                          ! long-range virial (custom NB pot)
  real(8),intent(in):: pot_viricstmnbex_med(:)
                                          ! med-range virial (custom NB pot)
  real(8),intent(in):: pot_viricstmnbex_short(:)
                                          ! short-range virial (custom NB pot)
  real(8),intent(in):: pot_viricstmnbext_long(:,:,:)
                                          ! long-range virial tensor (custom NB)
  real(8),intent(in):: pot_viricstmnbext_med(:,:,:)
                                          ! med-range virial tensor (custom NB)
  real(8),intent(in):: pot_viricstmnbext_short(:,:,:)
                                          ! short-range virial tensor(custom NB)

  real(8),intent(in):: atm_viribo_long        ! virial(bond potential)
  real(8),intent(in):: atm_viribo_med         ! virial(bond potential)
  real(8),intent(in):: atm_viribo_short       ! virial(bond potential)
  real(8),intent(in):: atm_viribot_long(:,:)  ! virial tensor (bond)
  real(8),intent(in):: atm_viribot_med(:,:)   ! virial tensor (bond)
  real(8),intent(in):: atm_viribot_short(:,:) ! virial tensor (bond)

  real(8),intent(in):: atm_virian_long        ! virial(angle potential)
  real(8),intent(in):: atm_virian_med         ! virial(angle potential)
  real(8),intent(in):: atm_virian_short       ! virial(angle potential)
  real(8),intent(in):: atm_viriant_long(:,:)  ! virial tensor (angle)
  real(8),intent(in):: atm_viriant_med(:,:)   ! virial tensor (angle)
  real(8),intent(in):: atm_viriant_short(:,:) ! virial tensor (angle)

  real(8),intent(in):: atm_virito_long        ! virial(torsion potential)
  real(8),intent(in):: atm_virito_med         ! virial(torsion potential)
  real(8),intent(in):: atm_virito_short       ! virial(torsion potential)
  real(8),intent(in):: atm_viritot_long(:,:)  ! virial tensor (torsion)
  real(8),intent(in):: atm_viritot_med(:,:)   ! virial tensor (torsion)
  real(8),intent(in):: atm_viritot_short(:,:) ! virial tensor (torsion)

  real(8),intent(in):: atm_viri14_long        ! virial(1-4 force potential)
  real(8),intent(in):: atm_viri14_med         ! virial(1-4 force potential)
  real(8),intent(in):: atm_viri14_short       ! virial(1-4 force potential)
  real(8),intent(in):: atm_viri14t_long(:,:)  ! virial tensor(1-4 force)
  real(8),intent(in):: atm_viri14t_med(:,:)   ! virial tensor(1-4 force)
  real(8),intent(in):: atm_viri14t_short(:,:) ! virial tensor(1-4 force)

  real(8),intent(in):: atm_viri_const         ! virial(constraint force)
  real(8),intent(in):: atm_virit_const(:,:)   ! virial tensor (constraint)

  real(8),intent(inout):: atm_viri_corr          ! virial(L-J correction)
  real(8),intent(inout):: atm_virit_corr(:,:)    ! virial tensor (L-J correction)

  integer,intent(in):: pressinterval   ! interval of pressure output

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  integer,intent(in):: current_step    ! current time step

  logical,intent(in):: ifcalpremole    ! pressure calculation of molecule
  logical,intent(in):: ifcalpreatom    ! pressure calculation of atom

  real(8),intent(in):: pot_ewc          ! coulomb potential(ewald self)

  real(8),intent(in):: pref             ! pressure base value [Pa]
  real(8),intent(in):: eref             ! energy base value [J]

  real(8),intent(in):: pot_ewnq         ! coulomb potential(ewald netq)

  integer,intent(in):: ncstmnbex        ! number of extra custom NB output

!     OUTPUT
  real(8),intent(out):: pressmol_ktot    ! pressure kinetic part of molecule
  real(8),intent(out):: pressmol_vtot    ! pressure virial part of molecule
  real(8),intent(out):: pressmol_tot     ! pressure total of molecule
  real(8),intent(out):: pressatm_ktot    ! pressure kinetic part of atm
  real(8),intent(out):: pressatm_vtot    ! pressure virial part of atm
  real(8),intent(out):: pressatm_tot     ! pressure total of atm

  real(8),intent(out):: pot_viric_all
  real(8),intent(out):: pot_virilj_all
  real(8),intent(out):: pot_virimor_all
  real(8),intent(out):: pot_virish_all
  real(8),intent(out):: pot_virirfh_all
  real(8),intent(out):: pot_viridou_all
  real(8),intent(out):: pot_viricstmnb_all
  real(8),intent(out):: pot_viricstmnbex_all(:)
  real(8),intent(out):: pot_viri_all
  real(8),intent(out):: viri_fdotd

  real(8),intent(out):: atm_viribo_all   ! virial(bond potential) of each atom
  real(8),intent(out):: atm_virian_all   ! virial(angle potential) of each atom
  real(8),intent(out):: atm_virito_all   ! virial(torsion potential) of each atom
  real(8),intent(out):: atm_viri14_all ! virial(1-4 force potential) of each atom

  real(8),intent(out):: pressmolt_ktot(:,:)
                     ! pressure tensor kinetic part of molecule
  real(8),intent(out):: pressmolt_vtot(:,:)
                     ! pressure tensor virial part of molecule
  real(8),intent(out):: pressmolt_tot(:,:)  ! pressure tensor total of molecule
  real(8),intent(out):: pressatmt_ktot(:,:) ! pressure tensor kinetic part of atm
  real(8),intent(out):: pressatmt_vtot(:,:) ! pressure tensor virial part of atm
  real(8),intent(out):: pressatmt_tot(:,:)  ! pressure tensor total of atm

  real(8),intent(out):: pot_virict_all(:,:)  ! virial tensor (coulomb potential)
  real(8),intent(out):: pot_viriljt_all(:,:) ! virial tensor (L-J potential)
  real(8),intent(out):: pot_virimort_all(:,:) ! virial tensor (Morse potential)
  real(8),intent(out):: pot_virisht_all(:,:) ! virial tensor (SH potential)
  real(8),intent(out):: pot_virirfht_all(:,:) ! virial tensor (RFH potential)
  real(8),intent(out):: pot_viridout_all(:,:) ! virial tensor (DOU potential)
  real(8),intent(out):: pot_viricstmnbt_all(:,:)
                                          ! virial tensor (custom NB potential)
  real(8),intent(out):: pot_viricstmnbext_all(:,:,:) 
                                     ! virial tensor (extra custom NB potential)

  real(8),intent(out):: pot_virit_all(:,:) ! all virial tensor of potential term
  real(8),intent(out):: virit_fdotd(:,:) ! virial tensor of correction term F.d

  real(8),intent(out):: atm_viribot_all(:,:)
                     ! virial(bond potential) of each atom
  real(8),intent(out):: atm_viriant_all(:,:)
                     ! virial(angle potential) of each atom
  real(8),intent(out):: atm_viritot_all(:,:)
                     ! virial(torsion potential) of each atom
  real(8),intent(out):: atm_viri14t_all(:,:)
                     ! virial(1-4 force potential) of each atom

  real(8),intent(in):: rcut             ! vdw cutoff length [non-d]

  logical,intent(in):: ifcalljlong      ! long-range correction in pressure
  integer,intent(in):: nsolve           ! number of solvent molecules
!      real*8:: vdw_welij_solve  ! well depth of vdw parameter of solvent
!      real*8:: vdw_radij_solve  ! radius of vdw parameter of solvent
  integer,intent(in):: solveindex       ! atmindex of solvent atom

! LOCAL:

  real(8):: molemass         ! mass of molecule
  real(8):: inv_molemass     ! = 1/molemass
  real(8):: molegrav(3)      ! gravity point of molecule
  real(8):: molegvel(3)      ! gravity velcity of molecule

  real(8):: box(3)
  real(8):: box_inv(3)

  real(8):: for_viri_c_all(3,maxnatom)
  real(8):: for_viri_lj_all(3,maxnatom)
  real(8):: for_viri_mor_all(3,maxnatom)
  real(8):: for_viri_sh_all(3,maxnatom)
  real(8):: for_viri_rfh_all(3,maxnatom)
  real(8):: for_viri_dou_all(3,maxnatom)
  real(8):: for_viri_cstmnb_all(3,maxnatom)
  real(8):: for_viri_cstmnbex_all(3,maxnatom,1:ncstmnbex)
  real(8):: for_viri_all(3,maxnatom)

  real(8):: pressmol_kineall ! pressure kinetic part of system(mole)
  real(8):: pressmol_viriall ! pressure virial part of system(mole)
  real(8):: pressatm_kineall ! pressure kinetic part of system(atom)
  real(8):: pressatm_viriall ! pressure virial part of system(atom)

  real(8):: pressmolt_kineall(3,3) ! pressure tensor kinetic part(mole)
  real(8):: pressmolt_viriall(3,3) ! pressure tensor virial part(mole)
  real(8):: pressatmt_kineall(3,3) ! pressure tensor kinetic part(atom)
  real(8):: pressatmt_viriall(3,3) ! pressure tensor virial part(atom)

  integer:: i,n
  integer:: im,i1,i2,ii

  integer:: itype
  real(8):: rvdw6,rvdw12
  real(8):: rc3,rc9
  real(8):: lj_corr_tmp

  real(8):: pi               ! = 3.141592...

  real(8):: atmvel_tmp(3,natom) ! temporary atom velocity without streaming vel.

!     +     +     +     +     +     +     +

!     --- some preparation ---

  pi = acos(-1.0d0)

  atmvel_tmp(1:3,1:natom) = atmvel(1:3,1:natom)
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) - atmvel_strm(1:3,1:natom)

!     --- initialize atm_viri_corr ---
  atm_viri_corr = 0.0d0
  atm_virit_corr(1:3,1:3) = 0.0d0

!     --- initialize press_kineall & press_viriall ---
  pressmol_kineall = 0.0d0
  pressmol_viriall = 0.0d0
  pressatm_kineall = 0.0d0
  pressatm_viriall = 0.0d0

  pressmolt_kineall(1:3,1:3) = 0.0d0
  pressmolt_viriall(1:3,1:3) = 0.0d0
  pressatmt_kineall(1:3,1:3) = 0.0d0
  pressatmt_viriall(1:3,1:3) = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)


  pot_viric_all  = pot_viric_long   &
       ! &         + pot_viric_med &
       &         + pot_viric_short   &
       &         + pot_ewc   &
       &         + pot_ewnq 
  pot_virilj_all = pot_virilj_long   &
!     &           + pot_virilj_med &
       &         + pot_virilj_short
  pot_virimor_all = pot_virimor_long   &
!     &            + pot_virimor_med &
       &          + pot_virimor_short
  pot_virish_all = pot_virish_long   &
!     &            + pot_virish_med &
       &         + pot_virish_short
  pot_virirfh_all = pot_virirfh_long   &
!     &            + pot_virirfh_med &
       &          + pot_virirfh_short
  pot_viridou_all = pot_viridou_long   &
!     &            + pot_viridou_med &
       &          + pot_viridou_short
  pot_viricstmnb_all = pot_viricstmnb_long   &
!     &               + pot_viricstmnb_med &
       &             + pot_viricstmnb_short
  pot_viricstmnbex_all(1:ncstmnbex) = pot_viricstmnbex_long(1:ncstmnbex) &
       ! &                            + pot_viricstmnbex_med(1:ncstmnbex) &
       &                            + pot_viricstmnbex_short(1:ncstmnbex)

  pot_viri_all = pot_viric_all + pot_virilj_all   &
       &       + pot_virimor_all + pot_virish_all   &
       &       + pot_virirfh_all + pot_viridou_all  &
       &       + pot_viricstmnb_all

#if defined(_CSTMNB_V2_ADD_ALL)
  !- adding cstmnb extra potential for virial
  do i = 1,ncstmnbex
     pot_viri_all = pot_viri_all + pot_viricstmnbex_all(i)
  end do
#endif

  pot_virict_all(1:3,1:3) =  pot_virict_long(1:3,1:3)   &
!     &                           + pot_virict_med(1,n)
       &                   + pot_virict_short(1:3,1:3)
  pot_viriljt_all(1:3,1:3) =  pot_viriljt_long(1:3,1:3)   &
!     &                           + pot_viriljt_med(1,n) 
       &                    + pot_viriljt_short(1:3,1:3)
  pot_virimort_all(1:3,1:3) =  pot_virimort_long(1:3,1:3)   &
!     &                           + pot_virimort_med(1,n) 
       &                     + pot_virimort_short(1:3,1:3)
  pot_virisht_all(1:3,1:3) =  pot_virisht_long(1:3,1:3)   &
!     &                           + pot_virisht_med(1,n) 
       &                    + pot_virisht_short(1:3,1:3)
  pot_virirfht_all(1:3,1:3) =  pot_virirfht_long(1:3,1:3)   &
!     &                           + pot_virirfht_med(1,n) 
       &                     + pot_virirfht_short(1:3,1:3)
  pot_viridout_all(1:3,1:3) =  pot_viridout_long(1:3,1:3)   &
!     &                           + pot_viridout_med(1,n) 
       &                     + pot_viridout_short(1:3,1:3)
  pot_viricstmnbt_all(1:3,1:3) =  pot_viricstmnbt_long(1:3,1:3)   &
!     &                           + pot_viricstmnbt_med(1,n)
       &                        + pot_viricstmnbt_short(1:3,1:3)
  pot_viricstmnbext_all(1:3,1:3,1:ncstmnbex) = &
       &                          pot_viricstmnbext_long(1:3,1:3,1:ncstmnbex) &
       ! &                        + pot_viricstmnbext_med(1:3,1:3,1:ncstmnbex) &
       &                        + pot_viricstmnbext_short(1:3,1:3,1:ncstmnbex)

  pot_virit_all(1:3,1:3) =  pot_virict_all(1:3,1:3) &
       &                  + pot_viriljt_all(1:3,1:3) &
       &                  + pot_virimort_all(1:3,1:3) &
       &                  + pot_virisht_all(1:3,1:3) &
       &                  + pot_virirfht_all(1:3,1:3) &
       &                  + pot_viridout_all(1:3,1:3) &
       &                  + pot_viricstmnbt_all(1:3,1:3)

#if defined(_CSTMNB_V2_ADD_ALL)
  !- adding cstmnb extra potential for virial
  do i = 1,ncstmnbex
     pot_virit_all(1:3,1:3) = pot_virit_all(1:3,1:3) &
          &                 + pot_viricstmnbext_all(1:3,1:3,i)
  end do
#endif


  for_viri_c_all(1:3,1:natom) =  for_viric_long(1:3,1:natom)   &
!     &                          + for_viric_med(1,i)
       &                        + for_viric_short(1:3,1:natom)
  for_viri_lj_all(1:3,1:natom) =  for_virilj_long(1:3,1:natom)   &
!     &                           + for_virilj_med(1,i)
       &                         + for_virilj_short(1:3,1:natom)
  for_viri_mor_all(1:3,1:natom) =  for_virimor_long(1:3,1:natom)   &
!     &                           + for_virimor_med(1,i)
       &                         + for_virimor_short(1:3,1:natom)
  for_viri_sh_all(1:3,1:natom) =  for_virish_long(1:3,1:natom)   &
!     &                           + for_virish_med(1,i)
       &                         + for_virish_short(1:3,1:natom)
  for_viri_rfh_all(1:3,1:natom) =  for_virirfh_long(1:3,1:natom)   &
!     &                           + for_virirfh_med(1,i)
       &                         + for_virirfh_short(1:3,1:natom)
  for_viri_dou_all(1:3,1:natom) =  for_viridou_long(1:3,1:natom)   &
!     &                           + for_viridou_med(1,i)
       &                         + for_viridou_short(1:3,1:natom)
  for_viri_cstmnb_all(1:3,1:natom) =  for_viricstmnb_long(1:3,1:natom)   &
!     &                              + for_viricstmnb_med(1,i)
       &                            + for_viricstmnb_short(1:3,1:natom)
  for_viri_cstmnbex_all(1:3,1:natom,1:ncstmnbex) = &
       &                      for_viricstmnbex_long(1:3,1:natom,1:ncstmnbex) &
       ! &                    + for_viricstmnbex_med(1:3,1:natom,1:ncstmnbex) &
       &                    + for_viricstmnbex_short(1:3,1:natom,1:ncstmnbex)

  for_viri_all(1:3,1:natom) =  for_viri_c_all(1:3,1:natom) &
       &                     + for_viri_lj_all(1:3,1:natom) &
       &                     + for_viri_mor_all(1:3,1:natom) &
       &                     + for_viri_sh_all(1:3,1:natom) &
       &                     + for_viri_rfh_all(1:3,1:natom) &
       &                     + for_viri_dou_all(1:3,1:natom) &
       &                     + for_viri_cstmnb_all(1:3,1:natom)

#if defined(_CSTMNB_V2_ADD_ALL)
  !- adding cstmnb extra force for virial
  do i = 1, ncstmnbex
     for_viri_all(1:3,1:natom) = for_viri_all(1:3,1:natom) &
          &                    + for_viri_cstmnbex_all(1:3,1:natom,i)
  end do
#endif

  atm_viribo_all = atm_viribo_long + atm_viribo_short
  atm_virian_all = atm_virian_long + atm_virian_short
  atm_virito_all = atm_virito_long + atm_virito_short
  atm_viri14_all = atm_viri14_long + atm_viri14_short


  atm_viribot_all(1:3,1:3) =  atm_viribot_long(1:3,1:3)   &
       &                    + atm_viribot_short(1:3,1:3)
  atm_viriant_all(1:3,1:3) =  atm_viriant_long(1:3,1:3)   &
       &                    + atm_viriant_short(1:3,1:3)
  atm_viritot_all(1:3,1:3) =  atm_viritot_long(1:3,1:3)   &
       &                    + atm_viritot_short(1:3,1:3)
  atm_viri14t_all(1:3,1:3) =  atm_viri14t_long(1:3,1:3)   &
       &                    + atm_viri14t_short(1:3,1:3)


!     --- loop over molecule(im) ---

  DO im = 1, nmoleptindex-1
     i1 = molept_index(im)
     i2 = molept_index(im+1) - 1

     molemass = 0.0d0
     molegrav(1:3) = 0.0d0
     molegvel(1:3) = 0.0d0

     DO ii = i1, i2         ! loop over atom(i)
        i = molept_list(ii)

        molemass = molemass + atmmass(i)

!           - calculate gravity point -
        molegrav(1:3) =  molegrav(1:3)   &
             &       + atmmass(i) * atmcor(1:3,i)
!           - calculate gravity velocity -
        molegvel(1:3) =  molegvel(1:3)   &
             &       + atmmass(i) * atmvel(1:3,i)

!           - calculate atomic kinetic energy term -
        if (ifcalpreatom) then
           pressatm_kineall =  pressatm_kineall   &
                &            + atmmass(i)*   &
                &              (  atmvel(1,i)*atmvel(1,i)   &
                &               + atmvel(2,i)*atmvel(2,i)   &
                &               + atmvel(3,i)*atmvel(3,i))

           do n=1,3
              pressatmt_kineall(1:3,n) =  pressatmt_kineall(1:3,n)   &
                   &                  + atmmass(i)   &
                   &                  * atmvel(1:3,i)*atmvel(n,i)
           end do
        end if

     END DO

     inv_molemass = 1.0d0 / molemass

     molegrav(1:3) = molegrav(1:3) * inv_molemass
     molegvel(1:3) = molegvel(1:3) * inv_molemass

     DO ii = i1, i2         ! loop over atom(i)
        i = molept_list(ii)

        pressmol_viriall =  pressmol_viriall    & ! calculate -F.d
             &            - (for_viri_all(1,i)   &
             &               *(atmcor(1,i)-molegrav(1))   &
             &              +for_viri_all(2,i)   &
             &               *(atmcor(2,i)-molegrav(2))   &
             &              +for_viri_all(3,i)   &
             &               *(atmcor(3,i)-molegrav(3)))

        do n=1,3
           pressmolt_viriall(1:3,n) =  pressmolt_viriall(1:3,n)   &
                &                  - for_viri_all(1:3,i)   &
                &                   *(atmcor(n,i)-molegrav(n))
        end do

     END DO

     pressmol_kineall =  pressmol_kineall   &
          &            + molemass *   &
          &             (molegvel(1)**2 +   &
          &              molegvel(2)**2 +   &
          &              molegvel(3)**2)

     do n=1,3
        pressmolt_kineall(1,n) =  pressmolt_kineall(1,n)   &
             &                  + molemass   &
             &                   *molegvel(1)*molegvel(n)
        pressmolt_kineall(2,n) =  pressmolt_kineall(2,n)   &
             &                  + molemass   &
             &                   *molegvel(2)*molegvel(n)
        pressmolt_kineall(3,n) =  pressmolt_kineall(3,n)   &
             &                  + molemass   &
             &                   *molegvel(3)*molegvel(n)
     end do

  END DO

!     --- calculate pressure of L-J long-range correction ---
  if (ifcalljlong) then

     do i=1,natom

        itype = atmindex(i) ! VDW type of atom(i)
        rvdw6 = vdw_radij(itype,solveindex)**6
        rvdw12 = rvdw6*rvdw6
        rc3 = rcut**3
        rc9 = rc3*rc3*rc3

        lj_corr_tmp =  2.0d0 * pi   &
             &       * nsolve*box_inv(1)*box_inv(2)*box_inv(3)   &
             &       * vdw_welij(itype,solveindex)   &
             &       * (rvdw12/(3.0d0*rc9) - 2.0d0*rvdw6/rc3)

        atm_viri_corr = atm_viri_corr + lj_corr_tmp

        atm_virit_corr(1,1) =  atm_virit_corr(1,1)   &
             &               + lj_corr_tmp/3.0d0
        atm_virit_corr(2,2) =  atm_virit_corr(2,2)   &
             &               + lj_corr_tmp/3.0d0
        atm_virit_corr(3,3) =  atm_virit_corr(3,3)   &
             &               + lj_corr_tmp/3.0d0

     end do

  end if

!     --- sum up pressure ---
  viri_fdotd = pressmol_viriall
  pressmol_viriall =  pressmol_viriall + pot_viri_all   & ! add pot term
       &            + atm_viri_corr

  pressatm_viriall = pressatm_viriall + pot_viri_all ! add pot term
  pressatm_viriall =  pressatm_viriall   &           ! add internal term
       &            + atm_viribo_all + atm_virian_all   &
       &            + atm_virito_all + atm_viri14_all   &
       &            + atm_viri_const + atm_viri_corr

  virit_fdotd(1:3,1:3) = pressmolt_viriall(1:3,1:3)
  pressmolt_viriall(1:3,1:3) =  pressmolt_viriall(1:3,1:3)   &
       &                  + pot_virit_all(1:3,1:3)   &  ! add pot term
       &                  + atm_virit_corr(1:3,1:3)

  pressatmt_viriall(1:3,1:3) =  pressatmt_viriall(1:3,1:3)   &
       &                      + pot_virit_all(1:3,1:3)     ! add pot term
  pressatmt_viriall(1:3,1:3) =  pressatmt_viriall(1:3,1:3)  & ! add internal term
       &                      + atm_viribot_all(1:3,1:3)   &
       &                      + atm_viriant_all(1:3,1:3)   &
       &                      + atm_viritot_all(1:3,1:3)   &
       &                      + atm_viri14t_all(1:3,1:3)   &
       &                      + atm_virit_const(1:3,1:3)   &
       &                      + atm_virit_corr(1:3,1:3)

!      if (mod(current_step,pressinterval) .eq. 0) then

!         pressmol_ktot = pressmol_kineall / dble(pressinterval)
!         pressmol_vtot = pressmol_viriall / dble(pressinterval)
  pressmol_ktot = pressmol_kineall
  pressmol_vtot = pressmol_viriall
  pressatm_ktot = pressatm_kineall
  pressatm_vtot = pressatm_viriall

  pressmol_tot =  1.0d0/3.0d0 * box_inv(1)*box_inv(2)*box_inv(3)   &
       &        * (pressmol_ktot + pressmol_vtot)

  pressatm_tot =  1.0d0/3.0d0 * box_inv(1)*box_inv(2)*box_inv(3)   &
       &        * (pressatm_ktot + pressatm_vtot)

  pressmolt_ktot(1:3,1:3) = pressmolt_kineall(1:3,1:3)
  pressmolt_vtot(1:3,1:3) = pressmolt_viriall(1:3,1:3)
  pressatmt_ktot(1:3,1:3) = pressatmt_kineall(1:3,1:3)
  pressatmt_vtot(1:3,1:3) = pressatmt_viriall(1:3,1:3)

  pressmolt_tot(1:3,1:3) =  1.0d0   &
       &                  * box_inv(1)*box_inv(2)*box_inv(3)   &
       &                  * (  pressmolt_ktot(1:3,1:3)   &
       &                     + pressmolt_vtot(1:3,1:3))

  pressatmt_tot(1:3,1:3) =  1.0d0   &
       &                  * box_inv(1)*box_inv(2)*box_inv(3)   &
       &                  * (  pressatmt_ktot(1:3,1:3)   &
       &                     + pressatmt_vtot(1:3,1:3))

  if (.not. ifcalpremole) then
     pressmol_tot  = 0.0d0
     pressmol_ktot = 0.0d0
     pressmol_vtot = 0.0d0
     pressmolt_tot(1:3,1:3)  = 0.0d0
     pressmolt_ktot(1:3,1:3) = 0.0d0
     pressmolt_vtot(1:3,1:3) = 0.0d0
  end if

  if (.not. ifcalpreatom) then
     pressatm_tot  = 0.0d0
     pressatm_ktot = 0.0d0
     pressatm_vtot = 0.0d0
     pressatmt_tot(1:3,1:3)  = 0.0d0
     pressatmt_ktot(1:3,1:3) = 0.0d0
     pressatmt_vtot(1:3,1:3) = 0.0d0
  end if

!--- restore velocity
  atmvel(1:3,1:natom) = atmvel_tmp(1:3,1:natom) 

!     +     +     +     +     +     +     +

end subroutine calpress
