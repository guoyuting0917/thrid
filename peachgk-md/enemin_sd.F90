!********************************
!*  enemin_sd.f90 Ver.1.8       *
!*      for peachgk_md.f        *
!*      EM platform subroutine  *
!*            by G.Kikugawa     *
!********************************
! Time-stamp: <>

subroutine enemin_sd(iostarec,recinterval, &
     &               npoly,npolytyp,npoly_mole,npoly_atom, &
     &               nwater, &
     &               nmatom,nmatyp,nmatomtyp, &
     &               maxnstep,inistep,endstep, &
     &               xcel,ycel,zcel, &
     &               ifcenterfix_all, &
     &               ifcenterfix_poly, &
     &               ifcenterfix_water, &
     &               ifcenterfix_ma, &
     &               cenfix_free, &
     &               ifcenterfix_polytyp, &
     &               ifcenterfix_watertyp, &
     &               ifcenterfix_matyp, &
     &               mts_bond, mts_angl, mts_anglub, &
     &               mts_tors,mts_torsrb,mts_torsim, &
     &               mts_vdw,mts_ewr,mts_ewk,mts_vdw14,mts_elc14, &
     &               mts_mor,mts_sh,mts_rfh,mts_dou,mts_cnpvw, &
     &               mts_cstmnb, &
     &               mts_posres, &
     &               mts_potbias, &
     &               ifewald,alpha,kmax,rrcut, &
     &               ifspme,nfft1,nfft2,nfft3,pme_order, &
     &               eps0,div_factor_14vdw,div_factor_14elc, &
     &               rcut, &
     &               ifcellindex, &
     &               ifbook,nstep_book, &
     &               ifrattle,eps_rattle, &
     &               iflocalfix, &
     &               nlfix,index_nlfix, &
     &               iflocalfixz,iflocalfixzg, &
     &               ifcnp, &
     &               ifposres, &
     &               ifpotbias, &
     &               ifoutpdb,nstep_pdbout, &
     &               resname_poly_pdb, &
     &               resname_water_pdb, &
     &               resname_matom_pdb, &
     &               xref,eref,mref,qref, &
     &               vref,timeref,tempref,pref,fref,eps0ref, &
     &               degfree_poly, degfree_water, &
     &               degfree_ma, degfree_all, &
     &               rcut_book, &
     &               rcutmor,ifbookmor,rcut_bookmor,nstep_bookmor, &
     &               ifcellindex_mor, &
     &               rcutsh,ifbooksh,rcut_booksh,nstep_booksh, &
     &               ifcellindex_sh, &
     &               rcutrfhfo,ifbookrfhfo,rcut_bookrfhfo, &
     &               nstep_bookrfhfo, &
     &               ifcellindex_rfhfo, &
     &               rcutrfhoo,ifbookrfhoo,rcut_bookrfhoo, &
     &               nstep_bookrfhoo, &
     &               ifcellindex_rfhoo, &
     &               rcutrfhoh,ifbookrfhoh,rcut_bookrfhoh, &
     &               nstep_bookrfhoh, &
     &               ifcellindex_rfhoh, &
     &               rcutdouo,rcutindouo,ifbookdouo,rcut_bookdouo, &
     &               nstep_bookdouo, &
     &               ifcellindex_douo, &
     &               rcutdouh,ifbookdouh,rcut_bookdouh, &
     &               nstep_bookdouh, &
     &               ifcellindex_douh, &
     &               rcutrpvw,ifbookrpvw,rcut_bookrpvw, &
     &               nstep_bookrpvw, &
     &               ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
     &               for_long,for_short, &
     &               for_viric_long,for_viric_med,for_viric_short, &
     &               for_virilj_long,for_virilj_med,for_virilj_short, &
     &               for_virimor_long,for_virimor_med, &
     &               for_virimor_short, &
     &               for_virish_long,for_virish_med, &
     &               for_virish_short, &
     &               for_virirfh_long,for_virirfh_med, &
     &               for_virirfh_short, &
     &               for_viridou_long,for_viridou_med, &
     &               for_viridou_short, &
     &               for_viricstmnb_long, for_viricstmnb_med, &
     &               for_viricstmnb_short, &
     &               pot_ewc, &
     &               ouene,oupos,ouvel,oufor, &
     &               outhe,oubar,oupre,oupdb, &
     &               ouumb, &
     &               ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
     &               ifoutbar,ifoutpre, &
     &               outinterval,pressinterval, &
     &               yratio,zratio, &
     &               ifcalpremole,ifcalpreatom, &
     &               ifnetqcorrp, &
     &               mchain, &
     &               pint, pintt, &
     &               ifpatmcont,ifpmolcont, &
     &               ifcalljlong,nsolve,solveindex, &
     &               netchrgsq, &
     &               nstep_expand,r_expand, &
     &               d_rini,d_rmax,d_econv,d_rmsf, &
     &               md_cont)

!     +     +     +     +     +     +     +

  use interface_mdtech
#if defined(TIME_M) || defined(TIME_MALL)
  use interface_timer
#endif

  use md_global

#if defined(MPI)
  use mpi_global
#endif

#if defined(_CSTMNB_V2)
  use cstmnb, only: ncstmnbex_in
#endif

#if defined(TIME_M) || defined(TIME_MALL)
  use time_global
#endif

  implicit none

!---- i/o unit
  integer,intent(in):: iostarec        ! state record file unit

!---- MD control parameters
  integer,intent(in):: npoly           ! all number of poly
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatom          ! number of monatomic molecules
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  integer,intent(in):: maxnstep        ! maximum step of MD
  integer,intent(in):: inistep         ! initial step of MD
  integer,intent(in):: endstep         ! end step of MD

  real(8),intent(inout):: xcel         ! x cell length[non-d]
  real(8),intent(inout):: ycel         ! y cell length[non-d]
  real(8),intent(inout):: zcel         ! z cell length[non-d]

  real(8),intent(in):: yratio           ! y cell ratio of y to x
  real(8),intent(in):: zratio           ! z cell ratio of z to x

  logical,intent(in):: ifcenterfix_all   ! center fix for all
  logical,intent(in):: ifcenterfix_poly  ! center fix for polymer
  logical,intent(in):: ifcenterfix_water ! center fix for water
  logical,intent(in):: ifcenterfix_ma    ! center fix for monatomic mole.
  character(4),intent(in):: cenfix_free  ! COM not fixed in this direction

  logical,intent(in):: ifcenterfix_polytyp(:) ! center fix for each polymer
  logical,intent(in):: ifcenterfix_watertyp   ! center fix for each water
  logical,intent(in):: ifcenterfix_matyp(:)
                                          ! center fix for each monatomic mole.

  integer,intent(inout):: mts_bond        ! MTS flag for bond
  integer,intent(inout):: mts_angl        ! MTS flag for angle
  integer,intent(inout):: mts_anglub      ! MTS flag for Urey-Bradley angle
  integer,intent(inout):: mts_tors        ! MTS flag for torsion
  integer,intent(inout):: mts_torsrb      ! MTS flag for torsionrb
  integer,intent(inout):: mts_torsim      ! MTS flag for torsionim
  integer,intent(inout):: mts_vdw         ! MTS flag for vdw interaction
  integer,intent(inout):: mts_ewr         ! MTS flag for ewald real(=vdw)
  integer,intent(inout):: mts_ewk         ! MTS flag for ewald wave
  integer,intent(inout):: mts_vdw14       ! MTS flag for 14vdw
  integer,intent(inout):: mts_elc14       ! MTS flag for 14elc
  integer,intent(inout):: mts_mor         ! MTS flag for Morse interaction
  integer,intent(inout):: mts_sh          ! MTS flag for SH interaction
  integer,intent(inout):: mts_rfh         ! MTS flag for RFH interaction
  integer,intent(inout):: mts_dou         ! MTS flag for DOU interaction
  integer,intent(inout):: mts_cnpvw       ! MTS flag for CNP_VW
  integer,intent(inout):: mts_cstmnb      ! MTS flag for custom NB interaction
  integer,intent(inout):: mts_posres      ! MTS flag for position restraint
  integer,intent(inout):: mts_potbias     ! MTS flag for bias potential

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

  logical,intent(in):: ifcellindex     ! flag for cell index

  logical,intent(in):: ifbook          ! flag for bookkeeping
  integer,intent(in):: nstep_book      ! bookkeeping interval

  logical,intent(in):: ifrattle        ! rattle flag

  real(8),intent(in):: rcutmor          ! Morse cutoff length [non-d]

  logical,intent(in):: ifcellindex_mor ! flag for cell index (morse)

  logical,intent(in):: ifbookmor      ! flag for bookkeeping of Morse interaction
  real(8),intent(in):: rcut_bookmor
                         ! cut off radius of bookkeeping[non-d] of Morse
  integer,intent(in):: nstep_bookmor  ! bookkeeping interval of Morse interaction

  real(8),intent(in):: rcutsh           ! SH cutoff length [non-d]

  logical,intent(in):: ifcellindex_sh  ! flag for cell index (SH)

  logical,intent(in):: ifbooksh        ! flag for bookkeeping of SH interaction
  real(8),intent(in):: rcut_booksh   ! cut off radius of bookkeeping[non-d] of SH
  integer,intent(in):: nstep_booksh    ! bookkeeping interval of SH interaction

  real(8),intent(in):: rcutrfhfo        ! RFH(FO) cutoff length [non-d]

  logical,intent(in):: ifcellindex_rfhfo ! flag for cell index (RFH(FO))

  logical,intent(in):: ifbookrfhfo  ! flag for bookkeeping of RFH(FO) interaction
  real(8),intent(in):: rcut_bookrfhfo
                         ! cut off radius of bookkeeping[non-d] of RFH(FO)
  integer,intent(in):: nstep_bookrfhfo
                         ! bookkeeping interval of RFH(FO) interaction

  real(8),intent(in):: rcutrfhoo        ! RFH(OO) cutoff length [non-d]

  logical,intent(in):: ifcellindex_rfhoo ! flag for cell index (RFH(OO))

  logical,intent(in):: ifbookrfhoo  ! flag for bookkeeping of RFH(OO) interaction
  real(8),intent(in):: rcut_bookrfhoo
                         ! cut off radius of bookkeeping[non-d] of RFH(OO)
  integer,intent(in):: nstep_bookrfhoo
                         ! bookkeeping interval of RFH(OO) interaction

  real(8),intent(in):: rcutrfhoh        ! RFH(OH) cutoff length [non-d]

  logical,intent(in):: ifcellindex_rfhoh ! flag for cell index (RFH(OH))

  logical,intent(in):: ifbookrfhoh  ! flag for bookkeeping of RFH(OH) interaction
  real(8),intent(in):: rcut_bookrfhoh
                         ! cut off radius of bookkeeping[non-d] of RFH(OH)
  integer,intent(in):: nstep_bookrfhoh
                         ! bookkeeping interval of RFH(OH) interaction

  real(8),intent(in):: rcutdouo         ! DOU cutoff length [non-d] for O-Au
  real(8),intent(in):: rcutindouo       ! DOU cutin length [non-d] for O-Au

  logical,intent(in):: ifcellindex_douo ! flag for cell index (DOU) for O-Au

  logical,intent(in):: ifbookdouo
                         ! flag for bookkeeping of DOU interaction (O-Au)
  real(8),intent(in):: rcut_bookdouo
                         ! cut off radius of bookkeep[non-d] of DOU (O-Au)
  integer,intent(in):: nstep_bookdouo
                         ! bookkeeping interval of DOU interaction (O-Au)

  real(8),intent(in):: rcutdouh         ! DOU cutoff length [non-d] for H-Au

  logical,intent(in):: ifcellindex_douh ! flag for cell index (DOU) for H-Au

  logical,intent(in):: ifbookdouh
                         ! flag for bookkeeping of DOU interaction (H-Au)
  real(8),intent(in):: rcut_bookdouh
                         ! cut off radius of bookkeep[non-d] of DOU (H-Au)
  integer,intent(in):: nstep_bookdouh
                         ! bookkeeping interval of DOU interaction (H-Au)

  real(8),intent(in):: rcutrpvw ! RP-VW cutoff length
  logical,intent(in):: ifbookrpvw ! flag for bookkeeping of RP-VW interaction
  real(8),intent(in):: rcut_bookrpvw
                         ! cut off radius of bookkeep of RP-VW interaction
  integer,intent(in):: nstep_bookrpvw
                         ! bookkeeping interval of RP-VW interaction

  logical,intent(in):: ifcstmnb           ! flag if using custom NB interaction
  logical,intent(in):: ifcellindex_cstmnb ! flag for cell index (custom NB)
  logical,intent(in):: ifbookcstmnb
                                ! flag for bookkeeping of custom NB interaction

!---- variable for fix atoms
  logical,intent(in):: iflocalfix      ! fix atoms flag

  integer,intent(in):: nlfix           ! number of fix atoms
  integer,intent(in):: index_nlfix(:)  ! index of fix atoms

  logical,intent(in):: iflocalfixz     ! flag for fixing z coordinate of atoms

  logical,intent(in):: iflocalfixzg    ! flag for fixing z coordinate of COM of molecules

!---- variable for controlong normal pressure
  logical,intent(in):: ifcnp           ! flag to control normal pressure

!---- variable for position restraint
  logical,intent(in):: ifposres        ! position restraint flag

!---- variable for bias potential
  logical,intent(in):: ifpotbias       ! bias potential flag

!---- PDB output
  logical,intent(in):: ifoutpdb        ! flag for outputting PDB format file
  integer,intent(in):: nstep_pdbout    ! MD step for output of PDB file
  character(4),intent(in):: resname_poly_pdb(:) ! residue name for poly
  character(4),intent(in):: resname_water_pdb ! residue name for water
  character(4),intent(in):: resname_matom_pdb(:) ! residue name for matom

!---- parameters for energy minimization
  real(8),intent(inout):: d_rini   ! initial displacement dr for EM [non-d]
  real(8),intent(inout):: d_rmax   ! maximum displacement dr for EM [non-d]

  real(8),intent(inout):: d_econv  ! convergence condition for energy
                                   !    in EM [non-d]
  real(8),intent(inout):: d_rmsf   ! convergence condition for
                                   !    root mean square force in EM [non-d]

!---- base value for non-dimensionalize
  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: eref             ! energy base value [J]
  real(8),intent(in):: mref             ! mass base value [kg]
  real(8),intent(in):: qref             ! charge base value [C]

  real(8),intent(in):: vref             ! velocity base value [m/s]
  real(8),intent(in):: timeref          ! time base value [sec]
  real(8),intent(in):: tempref          ! temperature base value [K]
  real(8),intent(in):: pref             ! pressure base value [Pa]
  real(8),intent(in):: fref             ! force base value [N]

  real(8),intent(in):: eps0ref          ! dielectric constant base value [c^2/Jm]

!---- i/o unit
  integer,intent(in):: ouene           ! output unit for output energy data
  integer,intent(in):: oupos           ! output unit for output position data
  integer,intent(in):: ouvel           ! output unit for output velocity data
  integer,intent(in):: oufor           ! output unit for output force data
  integer,intent(in):: outhe           ! output unit for output thermostat data
  integer,intent(in):: oubar           ! output unit for output barostat data
  integer,intent(in):: oupre           ! output unit for outpre velocity data
  integer,intent(in):: oupdb           ! output unit for outpdb PDB data
  integer,intent(in):: ouumb         ! output unit for outumb bias potential data

!---- flag for file output
  logical,intent(in):: ifoutene       ! if ouput energy file
  logical,intent(in):: ifoutpos       ! if ouput position file
  logical,intent(in):: ifoutvel       ! if ouput velocity file
  logical,intent(in):: ifoutfor       ! if ouput force file
  logical,intent(in):: ifoutthe       ! if ouput NVT file
  logical,intent(in):: ifoutbar       ! if ouput NPT file
  logical,intent(in):: ifoutpre       ! if ouput pressure file

!---- degree of freedom
  integer,intent(in):: degfree_poly(:) ! degree of freedom of polymer1
  integer,intent(in):: degfree_water   ! degree of freedom of H2O
  integer,intent(in):: degfree_ma(:)   ! degree of freedom of monatomic mole.
  integer,intent(in):: degfree_all     ! degree of freedom of all atoms

!---- radius of bookkeeping
  real(8),intent(in):: rcut_book        ! cut off radius of bookkeeping

!---- force of each atoms
  real(8),intent(inout):: for_long(:,:)  ! long-range force
!  real(8),intent(inout):: for_med(:,:)              ! medium-range force
  real(8),intent(inout):: for_short(:,:) ! short-range force

!---- potential valiables
  real(8):: pot_tot    = 0.0d0 ! all potential
  real(8):: pot_nonbon = 0.0d0 ! all non-bonded potential
  real(8):: pot_vdw    = 0.0d0 ! vdw potential
  real(8):: pot_elc    = 0.0d0 ! coulomb potential(ewald real)
  real(8):: pot_ewk    = 0.0d0 ! coulomb potential(ewald wave)
  real(8),intent(in):: pot_ewc ! coulomb potential(ewald self)
  real(8):: pot_ewnq   = 0.0d0 ! coulomb potential(ewald netq)
  real(8):: pot_vdw14  = 0.0d0 ! 1-4vdw potential
  real(8):: pot_elc14  = 0.0d0 ! 1-4elc potential
  real(8):: pot_bond   = 0.0d0 ! bond potential
  real(8):: pot_angl   = 0.0d0 ! angle potential
  real(8):: pot_anglub = 0.0d0 ! Urey-Bradley angle potential
  real(8):: pot_tors   = 0.0d0 ! torsion potential
  real(8):: pot_torsrb = 0.0d0 ! RBtorsion potential
  real(8):: pot_torsim = 0.0d0 ! improper torsion potential
  real(8):: pot_mor    = 0.0d0 ! Morse potential
  real(8):: pot_sh     = 0.0d0 ! SH potential
  real(8):: pot_rfh    = 0.0d0 ! RFH potential
  real(8):: pot_dou    = 0.0d0 ! DOU potential
  real(8):: pot_rpvw   = 0.0d0 ! RP-VW interaction
  real(8):: pot_vw     = 0.0d0 ! potential of constant force for VW
  real(8):: pot_cstmnb = 0.0d0 ! Custom NB potential
  real(8),allocatable:: pot_cstmnbex(:)   ! extra Custom NB potential
  real(8):: pot_posres = 0.0d0 ! position restraint potential
  real(8):: pot_pbias  = 0.0d0 ! bias potential

!---- kinematic energy

  real(8):: ene_kin_poly(maxnpolytyp) ! kinetic energy of polymer1
  real(8):: ene_kin_water    ! kinetic energy of water
  real(8):: ene_kin_ma(maxnmatyp) ! kinetic energy of monatomic mole.
  real(8):: ene_kin_all      ! kinetic energy

  real(8):: ene_kin_th       ! kinetic energy of NHC thermostat
  real(8):: ene_pot_th       ! potential energy of NHC thermostat

  real(8):: ene_kin_ba       ! kinetic energy of Andersen barostat
  real(8):: ene_pot_ba       ! potential energy of Andersen barostat
  real(8):: extra_pot_ba     ! extra energy of barostat (=kTxi(1))

!---- all energy

  real(8):: ene_tot          ! total energy
  real(8):: ene_conserve     ! total conserved quantity

!---- temperature

  real(8):: temp_poly(maxnpolytyp) ! temperature of polymer1
  real(8):: temp_water       ! temperature of H2O
  real(8):: temp_ma(maxnmatyp) ! temperature of MA
  real(8):: temp_all         ! temperature of all molecules

!---- interval of outputting data
  integer,intent(in):: outinterval     ! output interval of trajectory
  integer,intent(in):: pressinterval   ! interval of pressure output

!---- interval of state recording
  integer,intent(in):: recinterval     ! state record interval

!---- valiables for pressure

  logical,intent(in):: ifcalpremole    ! pressure calculation of molecule
  logical,intent(in):: ifcalpreatom    ! pressure calculation of atom
  logical,intent(in):: ifnetqcorrp     ! net charge correction for pressure

!     for molecular pressure
  real(8),intent(inout):: for_viric_long(:,:)
                            ! long-range virial (coulomb force)
  real(8),intent(inout):: for_viric_med(:,:)
                            ! medium-range virial (coulomb force)
  real(8),intent(inout):: for_viric_short(:,:)
                            ! short-range virial (coulomb force)
  real(8):: pot_viric_long              ! long-range virial (coulomb pot)
  real(8):: pot_viric_med               ! medium-range virial (coulomb pot)
  real(8):: pot_viric_short             ! short-range virial (coulomb pot)
  real(8):: pot_virict_long(3,3)        ! long-range virial tensor (coulomb)
  real(8):: pot_virict_med(3,3)         ! med-range virial tensor (coulomb)
  real(8):: pot_virict_short(3,3)       ! short-range virial tensor (coulomb)

  real(8),intent(inout):: for_virilj_long(:,:)
                            ! long-range virial (L-J force)
  real(8),intent(inout):: for_virilj_med(:,:)
                            ! medium-range virial (L-J force)
  real(8),intent(inout):: for_virilj_short(:,:)
                            ! short-range virial (L-J force)
  real(8):: pot_virilj_long              ! long-range virial (L-J pot)
  real(8):: pot_virilj_med               ! medium-range virial (L-J pot)
  real(8):: pot_virilj_short             ! short-range virial (L-J pot)
  real(8):: pot_viriljt_long(3,3)        ! long-range virial tensor (L-J)
  real(8):: pot_viriljt_med(3,3)         ! medium-range virial tensor (L-J)
  real(8):: pot_viriljt_short(3,3)       ! short-range virial tensor (L-J)

  real(8),intent(inout):: for_virimor_long(:,:)
                            ! long-range virial (Morse force)
  real(8),intent(inout):: for_virimor_med(:,:)
                            ! med-range virial (Morse force)
  real(8),intent(inout):: for_virimor_short(:,:)
                            ! short-range virial (Morse force)
  real(8):: pot_virimor_long              ! long-range virial (Morse pot)
  real(8):: pot_virimor_med               ! med-range virial (Morse pot)
  real(8):: pot_virimor_short             ! short-range virial (Morse pot)
  real(8):: pot_virimort_long(3,3)        ! long-range virial tensor (Morse)
  real(8):: pot_virimort_med(3,3)         ! med-range virial tensor (Morse)
  real(8):: pot_virimort_short(3,3)       ! short-range virial tensor(Morse)

  real(8),intent(inout):: for_virish_long(:,:)
                            ! long-range virial (SH force)
  real(8),intent(inout):: for_virish_med(:,:)
                            ! med-range virial (SH force)
  real(8),intent(inout):: for_virish_short(:,:)
                            ! short-range virial (SH force)
  real(8):: pot_virish_long              ! long-range virial (SH pot)
  real(8):: pot_virish_med               ! med-range virial (SH pot)
  real(8):: pot_virish_short             ! short-range virial (SH pot)
  real(8):: pot_virisht_long(3,3)        ! long-range virial tensor (SH)
  real(8):: pot_virisht_med(3,3)         ! med-range virial tensor (SH)
  real(8):: pot_virisht_short(3,3)       ! short-range virial tensor(SH)

  real(8),intent(inout):: for_virirfh_long(:,:)
                            ! long-range virial (RFH force)
  real(8),intent(inout):: for_virirfh_med(:,:)
                            ! med-range virial (RFH force)
  real(8),intent(inout):: for_virirfh_short(:,:)
                            ! short-range virial (RFH force)
  real(8):: pot_virirfh_long              ! long-range virial (RFH pot)
  real(8):: pot_virirfh_med               ! med-range virial (RFH pot)
  real(8):: pot_virirfh_short             ! short-range virial (RFH pot)
  real(8):: pot_virirfht_long(3,3)        ! long-range virial tensor (RFH)
  real(8):: pot_virirfht_med(3,3)         ! med-range virial tensor (RFH)
  real(8):: pot_virirfht_short(3,3)       ! short-range virial tensor(RFH)

  real(8),intent(inout):: for_viridou_long(:,:)
                            ! long-range virial (DOU force)
  real(8),intent(inout):: for_viridou_med(:,:)
                            ! med-range virial (DOU force)
  real(8),intent(inout):: for_viridou_short(:,:)
                            ! short-range virial (DOU force)
  real(8):: pot_viridou_long              ! long-range virial (DOU pot)
  real(8):: pot_viridou_med               ! med-range virial (DOU pot)
  real(8):: pot_viridou_short             ! short-range virial (DOU pot)
  real(8):: pot_viridout_long(3,3)        ! long-range virial tensor (DOU)
  real(8):: pot_viridout_med(3,3)         ! med-range virial tensor (DOU)
  real(8):: pot_viridout_short(3,3)       ! short-range virial tensor(DOU)

  real(8),intent(inout):: for_viricstmnb_long(:,:)
                            ! long-range virial (custom NB force)
  real(8),intent(inout):: for_viricstmnb_med(:,:)
                            ! med-range virial (custom NB force)
  real(8),intent(inout):: for_viricstmnb_short(:,:)
                            ! short-range virial (custom NB force)
  real(8):: pot_viricstmnb_long           ! long-range virial (custom NB pot)
  real(8):: pot_viricstmnb_med            ! med-range virial (custom NB pot)
  real(8):: pot_viricstmnb_short          ! short-range virial (custom NB pot)
  real(8):: pot_viricstmnbt_long(3,3)     ! long-range virial tensor (custom NB)
  real(8):: pot_viricstmnbt_med(3,3)      ! med-range virial tensor (custom NB)
  real(8):: pot_viricstmnbt_short(3,3)    ! short-range virial tensor(custom nB)

  integer:: ncstmnbex                     ! number of extra custom NB output
  real(8),allocatable:: for_viricstmnbex_long(:,:,:)
                                          ! long-range virial (custom NB force)
  real(8),allocatable:: for_viricstmnbex_med(:,:,:)
                                          ! med-range virial (custom NB force)
  real(8),allocatable:: for_viricstmnbex_short(:,:,:)
                                          ! short-range virial (custom NB force)
  real(8),allocatable:: pot_viricstmnbex_long(:)
                                          ! long-range virial (custom NB pot)
  real(8),allocatable:: pot_viricstmnbex_med(:)
                                          ! med-range virial (custom NB pot)
  real(8),allocatable:: pot_viricstmnbex_short(:)
                                          ! short-range virial (custom NB pot)
  real(8),allocatable:: pot_viricstmnbext_long(:,:,:)
                                          ! long-range virial tensor (custom NB)
  real(8),allocatable:: pot_viricstmnbext_med(:,:,:)
                                          ! med-range virial tensor (custom NB)
  real(8),allocatable:: pot_viricstmnbext_short(:,:,:)
                                          ! short-range virial tensor(custom NB)

!     for atomic pressure
  real(8):: atm_viribo_long        ! virial(bond potential) of each atom
  real(8):: atm_viribo_med         ! virial(bond potential) of each atom
  real(8):: atm_viribo_short       ! virial(bond potential) of each atom
  real(8):: atm_viribot_long(3,3)  ! virial tensor (bond) of each atom
  real(8):: atm_viribot_med(3,3)   ! virial tensor (bond) of each atom
  real(8):: atm_viribot_short(3,3) ! virial tensor (bond) of each atom

  real(8):: atm_virian_long        ! virial(angle potential) of each atom
  real(8):: atm_virian_med         ! virial(angle potential) of each atom
  real(8):: atm_virian_short       ! virial(angle potential) of each atom
  real(8):: atm_viriant_long(3,3)  ! virial tensor (angle) of each atom
  real(8):: atm_viriant_med(3,3)   ! virial tensor (angle) of each atom
  real(8):: atm_viriant_short(3,3) ! virial tensor (angle) of each atom

  real(8):: atm_virito_long        ! virial(torsion potential) of each atom
  real(8):: atm_virito_med         ! virial(torsion potential) of each atom
  real(8):: atm_virito_short       ! virial(torsion potential) of each atom
  real(8):: atm_viritot_long(3,3)  ! virial tensor (torsion) of each atom
  real(8):: atm_viritot_med(3,3)   ! virial tensor (torsion) of each atom
  real(8):: atm_viritot_short(3,3) ! virial tensor (torsion) of each atom

  real(8):: atm_viri14_long        ! virial(1-4 force potential) of each atom
  real(8):: atm_viri14_med         ! virial(1-4 force potential) of each atom
  real(8):: atm_viri14_short       ! virial(1-4 force potential) of each atom
  real(8):: atm_viri14t_long(3,3)  ! virial tensor(1-4 force) of each atom
  real(8):: atm_viri14t_med(3,3)   ! virial tensor(1-4 force) of each atom
  real(8):: atm_viri14t_short(3,3) ! virial tensor(1-4 force) of each atom

  real(8):: atm_viri_const         ! virial(constraint force) of each atom
  real(8):: atm_virit_const(3,3)   ! virial tensor (constraint) of each atom

  real(8):: atm_viri_corr          ! virial(L-J correction)
  real(8):: atm_virit_corr(3,3)    ! virial tensor (L-J correction)

!     for pressure sum
  real(8):: pressmol_ktot       ! pressure kinetic part of molecule
  real(8):: pressmolt_ktot(3,3) ! pressure tensor kinetic part of molecule
  real(8):: pressmol_vtot       ! pressure virial part of molecule
  real(8):: pressmolt_vtot(3,3) ! pressure tensor virial part of molecule
  real(8):: pressmol_tot        ! pressure total of molecule
  real(8):: pressmolt_tot(3,3)  ! pressure tensor total of molecule
  real(8):: pressatm_ktot       ! pressure kinetic part of atm
  real(8):: pressatmt_ktot(3,3) ! pressure tensor kinetic part of atm
  real(8):: pressatm_vtot       ! pressure virial part of atm
  real(8):: pressatmt_vtot(3,3) ! pressure tensor virial part of atm
  real(8):: pressatm_tot        ! pressure total of atm
  real(8):: pressatmt_tot(3,3)  ! pressure tensor total of atm

  real(8):: pot_viric_all          ! virial(coulomb potential)
  real(8):: pot_virict_all(3,3)    ! virial tensor (coulomb potential)
  real(8):: pot_virilj_all         ! virial(L-J potential)
  real(8):: pot_viriljt_all(3,3)   ! virial tensor (L-J potential)
  real(8):: pot_virimor_all        ! virial(Morse potential)
  real(8):: pot_virimort_all(3,3)  ! virial tensor (Morse potential)
  real(8):: pot_virish_all         ! virial(SH potential)
  real(8):: pot_virisht_all(3,3)   ! virial tensor (SH potential)
  real(8):: pot_virirfh_all        ! virial(RFH potential)
  real(8):: pot_virirfht_all(3,3)  ! virial tensor (RFH potential)
  real(8):: pot_viridou_all        ! virial(DOU potential)
  real(8):: pot_viridout_all(3,3)  ! virial tensor (DOU potential)
  real(8):: pot_viricstmnb_all       ! virial(custom NB potential)
  real(8):: pot_viricstmnbt_all(3,3) ! virial tensor (custom NB potential)
  real(8),allocatable:: pot_viricstmnbex_all(:)
                                     ! virial(extra custom NB potential)
  real(8),allocatable:: pot_viricstmnbext_all(:,:,:)
                                     ! virial tensor (extra custom NB potential)
  real(8):: pot_viri_all           ! all virial of potential term
  real(8):: pot_virit_all(3,3)     ! all virial tensor of potential term
  real(8):: viri_fdotd             ! virial of correction term F.d
  real(8):: virit_fdotd(3,3)       ! virial tensor of correction term F.d

  real(8):: atm_viribo_all       ! virial(bond potential) of each atom
  real(8):: atm_viribot_all(3,3) ! virial(bond potential) of each atom
  real(8):: atm_virian_all       ! virial(angle potential) of each atom
  real(8):: atm_viriant_all(3,3) ! virial(angle potential) of each atom
  real(8):: atm_virito_all       ! virial(torsion potential) of each atom
  real(8):: atm_viritot_all(3,3) ! virial(torsion potential) of each atom
  real(8):: atm_viri14_all       ! virial(1-4 force potential) of each atom
  real(8):: atm_viri14t_all(3,3) ! virial(1-4 force potential) of each atom

!---- for RATTLE constraint method

  real(8),intent(in):: eps_rattle       ! tolerance (relative difference)
                                ! for bond length constraint by RATTLE

!---- for Nose-Hoover chain
  integer,intent(in):: mchain          ! Nose-Hoover chain number (>1 in NVT)

!---- variables for Andersen (Hoover type) barostat
  real(8),intent(inout):: pint         ! internal pressure
  real(8),intent(inout):: pintt(:,:)   ! internal pressure tensor
  logical,intent(in):: ifpatmcont      ! atomic pressure control
  logical,intent(in):: ifpmolcont      ! molecular pressure control

!---- pressure calculation of L-J long-range correction
  logical,intent(in):: ifcalljlong     ! long-range correction in pressure
  integer,intent(in):: nsolve          ! number of solvent molecules
!      real*8:: vdw_welij_solve  ! well depth of vdw parameter of solvent
!      real*8:: vdw_radij_solve  ! radius of vdw parameter of solvent
  integer,intent(in):: solveindex      ! atmindex of solvent atom

!---- variables for net charge calculation
  real(8),intent(inout):: netchrgsq        ! = (sum(qi))**2

!---- cell expansion
  integer,intent(in):: nstep_expand    ! time step of cell expansion
  real(8),intent(in):: r_expand         ! expansion ratio of cell (volume)

!---- MD control parameter
  integer,intent(in):: md_cont         ! MD control flag

!     +     +     +     +     +     +     +

! LOCAL:
  integer:: current_step    ! current time step
!      integer:: istep_med
  integer:: istep_short
                                ! do loop index for medium and
                                ! short range forces
  integer:: i, j                ! do loop index

  real(8):: pot_tot_old ! potential energy of old step
  real(8):: atmcor_old(3,maxnatom) ! atmcor at time T for rattle-rRESPA
!      real*8:: atmcor_att(3,maxnatom) ! atmcor at time T for const. force

  real(8):: text_c = 0.0d0   ! external temp. in NVT(NHC) (don't use in NVE)
  real(8):: pext_c = 0.0d0   ! external pressure in NPT(MTK) (")

  integer:: ilfix

  real(8):: atmvel_tmp(3,maxnatom) ! for use of setting velocity of VW to 0

  real(8):: d_t                 ! alpha for coordinate update
  real(8):: d_r                 ! delta r
  real(8):: d_ene               ! delta pot_tot
  real(8):: frms                ! root mean square force
  real(8):: drms                ! root mean square force / natom
  logical:: minflag             ! energy is successfully minimized or not

! store original MTS flags
  integer:: mts_bond_old    ! MTS flag for bond
  integer:: mts_angl_old    ! MTS flag for angle
  integer:: mts_anglub_old  ! MTS flag for Urey-Bradley angle
  integer:: mts_tors_old    ! MTS flag for torsion
  integer:: mts_torsrb_old  ! MTS flag for torsionrb
  integer:: mts_torsim_old  ! MTS flag for torsionim
  integer:: mts_vdw_old     ! MTS flag for vdw interaction
  integer:: mts_ewr_old     ! MTS flag for ewald real(=vdw)
  integer:: mts_ewk_old     ! MTS flag for ewald wave
  integer:: mts_vdw14_old   ! MTS flag for 14vdw
  integer:: mts_elc14_old   ! MTS flag for 14elc
  integer:: mts_mor_old     ! MTS flag for Morse interaction
  integer:: mts_sh_old      ! MTS flag for SH interaction
  integer:: mts_rfh_old     ! MTS flag for RFH interaction
  integer:: mts_dou_old     ! MTS flag for DOU interaction
  integer:: mts_cnpvw_old   ! MTS flag for CNP_VW
  integer:: mts_cstmnb_old  ! MTS flag for custom NB interaction
  integer:: mts_posres_old  ! MTS flag for position restraint
  integer:: mts_potbias_old ! MTS flag for bias potential

!     +     +     +     +     +     +     +

!----- some initialization -----

  atm_viri_const = 0.0d0
  atm_virit_const(1:3,1:3) = 0.0d0

  d_r = d_rini                  ! delta r_ini is copied to local variable
  d_ene = 1.0d+16               ! large dummy value is set
  minflag = .false.

! --- Intitialize extra custom NB output arrays ---

!- CSTMNB function Ver. 1
#if defined(_CSTMNB_V1)
  ncstmnbex = 0

!- CSTMNB function Ver. 2
#elif defined(_CSTMNB_V2)
  ncstmnbex = ncstmnbex_in
#endif

  !- allocate memory
  allocate(pot_cstmnbex(1:ncstmnbex))
  allocate(for_viricstmnbex_long(3,natom,1:ncstmnbex))
  ! allocate(for_viricstmnbex_med(3,natom,1:ncstmnbex))
  allocate(for_viricstmnbex_short(3,natom,1:ncstmnbex))
  allocate(pot_viricstmnbex_long(1:ncstmnbex))
  ! allocate(pot_viricstmnbex_med(1:ncstmnbex))
  allocate(pot_viricstmnbex_short(1:ncstmnbex))
  allocate(pot_viricstmnbext_long(3,3,1:ncstmnbex))
  ! allocate(pot_viricstmnbext_med(3,3,1:ncstmnbex))
  allocate(pot_viricstmnbext_short(3,3,1:ncstmnbex))
  allocate(pot_viricstmnbex_all(1:ncstmnbex))
  allocate(pot_viricstmnbext_all(3,3,1:ncstmnbex))

!----- mts flag is forced to be "LONG"
! store original mts flags
  mts_bond_old = mts_bond
  mts_angl_old = mts_angl
  mts_anglub_old = mts_anglub
  mts_tors_old = mts_tors
  mts_torsrb_old = mts_torsrb
  mts_torsim_old = mts_torsim
  mts_vdw_old = mts_vdw
  mts_ewr_old = mts_ewr
  mts_ewk_old = mts_ewk
  mts_vdw14_old = mts_vdw14
  mts_elc14_old = mts_elc14
  mts_mor_old = mts_mor
  mts_sh_old = mts_sh
  mts_rfh_old = mts_rfh
  mts_dou_old = mts_dou
  mts_cnpvw_old = mts_cnpvw
  mts_cstmnb_old = mts_cstmnb
  mts_posres_old = mts_posres
  mts_potbias_old = mts_potbias

! overwrite LONG flag
  mts_bond = MTS_LONG
  mts_angl = MTS_LONG
  mts_anglub = MTS_LONG
  mts_tors = MTS_LONG
  mts_torsrb = MTS_LONG
  mts_torsim = MTS_LONG
  mts_vdw = MTS_LONG
  mts_ewr = MTS_LONG
  mts_ewk = MTS_LONG
  mts_vdw14 = MTS_LONG
  mts_elc14 = MTS_LONG
  mts_mor = MTS_LONG
  mts_sh = MTS_LONG
  mts_rfh = MTS_LONG
  mts_dou = MTS_LONG
  mts_cnpvw = MTS_LONG
  mts_cstmnb = MTS_LONG
  mts_posres = MTS_LONG
  mts_potbias = MTS_LONG

!-------- EM main loop --------

  DO current_step = inistep, endstep

!!! Time measure initialization
#if defined(TIME_M) || defined(TIME_MALL)
     call init_timecount()
#endif
!!! Time measure initialization

!!! Time measure
#if defined(TIME_M) || defined(TIME_MALL)
     call get_timecount_barrier( total_time, total_time_elps, 'START' )
#endif
!!! Time measure

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier( dyn_time, dyn_time_elps, 'START' )
#endif
!!! Time measure

!---- translate the coordinates (P.B.C.)

     if ((mod(current_step,nstep_book) == 1) .or. &
          & (nstep_book == 1) .or. (current_step == 0)) then

        call transcor( npoly, nwater, nmatom, &
             &         xcel, ycel, zcel )

     end if

!---- Expand cell volume ----
     if (current_step == nstep_expand) then
        call cell_expand( xref, current_step, &
             &            xcel, ycel, zcel, &
             &            yratio, zratio, &
             &            rcut_book, ifbook,rcut, &
             &            npoly, nwater, nmatom, &
             &            r_expand )

        if (ifewald) then
           call upewk( alpha, &
                &      xcel, ycel, zcel )
        end if

        if (ifspme) then
           call erf_corr_cutoff( xcel, ycel, zcel, &
                &                nfft1, nfft2, nfft3 )
        end if

     end if

!---- motion fix of molecules
     if (ifcenterfix_all .or. &
          & ifcenterfix_poly .or. &
          & ifcenterfix_water .or. &
          & ifcenterfix_ma) then
!        if ((mod(current_step,tcontinterval) == 1) .or.   &
!             &               (tcontinterval == 1)) then
        call cenfix(ifcenterfix_all, &
             &      ifcenterfix_poly, &
             &      ifcenterfix_water, &
             &      ifcenterfix_ma, &
             &      cenfix_free, &
             &      ifcenterfix_polytyp, &
             &      ifcenterfix_watertyp, &
             &      ifcenterfix_matyp, &
             &      npoly,npolytyp,npoly_mole,npoly_atom, &
             &      nwater, &
             &      nmatom,nmatyp,nmatomtyp, &
             &      iflocalfix, &
             &      nlfix,index_nlfix,   &
             &      iflocalfixz,iflocalfixzg)
!        end if
     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount( dyn_time, dyn_time_elps, 'STOP')
#endif
!!! Time measure

!---- make list of non-bonded pair list

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier( table_time, table_time_elps, 'START')
#endif
!!! Time measure

     if (ifbook .or. (current_step == 0)) then
        if ((mod(current_step,nstep_book) == 0) .or. &
             & (nstep_book == 1)) then

           if (ifcellindex) then
              call mklist_cell( rcut_book, &
                   &            xcel, ycel, zcel )
           else
              call mklist2a( rcut_book, &
                   &         xcel, ycel, zcel, ifbook )
           end if

        end if
     end if

!---- make list of Morse bonded pair list
     if (ifbookmor .or. (current_step == 0)) then
        if ((mod(current_step,nstep_bookmor) == 0) .or. &
             & (nstep_bookmor == 1)) then

           if (ifcellindex_mor) then
              call mklist_mor_cell( rcut_bookmor, &
                   &                xcel, ycel, zcel )
           else
              call mklist2a_mor( rcut_bookmor, &
                   &             xcel, ycel, zcel, &
                   &             ifbookmor )
           end if

        end if
     end if

!---- make list of SH nonbonded pair list
     if (ifbooksh .or. (current_step == 0)) then
        if ((mod(current_step,nstep_booksh) == 0) .or. &
             & (nstep_booksh == 1)) then

           if (ifcellindex_sh) then
              call mklist_sh_cell( rcut_booksh, &
                   &               xcel, ycel, zcel )
           else
              call mklist2a_sh( rcut_booksh, &
                   &            xcel, ycel, zcel, &
                   &            ifbooksh)
           end if

        end if
     end if

!---- make list of RFH nonbonded pair list
!        RFH(F-O)
     if (ifbookrfhfo .or. (current_step == 0)) then
        if ((mod(current_step,nstep_bookrfhfo) == 0) .or. &
             & (nstep_bookrfhfo == 1)) then

           if (ifcellindex_rfhfo) then
              call mklist_rfh_cell( rcut_bookrfhfo, &
                   &                xcel, ycel, zcel, &
                   &                INTTYPE_RFHFO)
           else
              call mklist2a_rfh( rcut_bookrfhfo, &
                   &             xcel, ycel, zcel, &
                   &             ifbookrfhfo, &
                   &             INTTYPE_RFHFO)
           end if

        end if
     end if

!        RFH(O-O)
     if (ifbookrfhoo .or. (current_step == 0)) then
        if ((mod(current_step,nstep_bookrfhoo) == 0) .or. &
             & (nstep_bookrfhoo == 1)) then

           if (ifcellindex_rfhoo) then
              call mklist_rfh_cell( rcut_bookrfhoo, &
                   &                xcel, ycel, zcel, &
                   &                INTTYPE_RFHOO )
           else
              call mklist2a_rfh( rcut_bookrfhoo, &
                   &             xcel,ycel,zcel, &
                   &             ifbookrfhoo, &
                   &             INTTYPE_RFHOO )
           end if

        end if
     end if

!        RFH(O-H)
     if (ifbookrfhoh .or. (current_step == 0)) then
        if ((mod(current_step,nstep_bookrfhoh) == 0) .or. &
             & (nstep_bookrfhoh == 1)) then

           if (ifcellindex_rfhoh) then
              call mklist_rfh_cell( rcut_bookrfhoh, &
                   &                xcel, ycel, zcel, &
                   &                INTTYPE_RFHOH )
           else
              call mklist2a_rfh( rcut_bookrfhoh, &
                   &             xcel, ycel, zcel, &
                   &             ifbookrfhoh, &
                   &             INTTYPE_RFHOH )
           end if

        end if
     end if

!---- make list of DOU nonbonded pair list
!        DOU(O-Au)
     if (ifbookdouo .or. (current_step == 0)) then
        if ((mod(current_step,nstep_bookdouo) == 0) .or. &
             & (nstep_bookdouo == 1)) then

           if (ifcellindex_douo) then
              call mklist_dou_cell( rcut_bookdouo, &
                   &                xcel, ycel, zcel, &
                   &                INTTYPE_DOUO )
           else
              call mklist2a_dou( rcut_bookdouo, &
                   &             xcel,ycel,zcel, &
                   &             ifbookdouo, &
                   &             INTTYPE_DOUO )
           end if

        end if
     end if

!        DOU(H-Au)
     if (ifbookdouh .or. (current_step == 0)) then
        if ((mod(current_step,nstep_bookdouh) == 0) .or. &
             & (nstep_bookdouh == 1)) then

           if (ifcellindex_douh) then
              call mklist_dou_cell( rcut_bookdouh, &
                   &                xcel, ycel, zcel, &
                   &                INTTYPE_DOUH )
           else
              call mklist2a_dou( rcut_bookdouh, &
                   &             xcel,ycel,zcel, &
                   &             ifbookdouh, &
                   &             INTTYPE_DOUH )
           end if

        end if
     end if

!---- make list of RP-VW interaction pair list
     if (ifbookrpvw .or. (current_step == 0)) then
        if ((mod(current_step,nstep_bookrpvw) == 0) .or. &
             & (nstep_bookrpvw == 1)) then

!!! Cell index omitted
!           if (ifcellindex_rpvw) then
!              call mklist_rpvw_cell( rcut_bookmor, &
!                   &                xcel, ycel, zcel )
!           else
           call mklist2a_rpvw( rcut_bookrpvw, zcel, ifbookrpvw)
!           end if

        end if
     end if

!---- make list of custom NB pair list
     if (ifcstmnb) then

        if (ifbookcstmnb .or. (current_step == 0)) then
           call mklist2a_cstmnb(ifcellindex_cstmnb,ifbookcstmnb, &
                &               xcel,ycel,zcel, &
                &               current_step)
        end if

     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount( table_time, table_time_elps, 'STOP')
#endif
!!! Time measure

!---- Energy minimization (steepest descent method)

!  --- CALCULATE FORCE AT T (CURRENT_TIME) IF CURRENT_STEP = 0 ---
!         * otherwise, the accell at time T has been already calculated

     IF (current_step == 0) THEN

!           - calculate initial force & energy -

        call calforce(MTS_LONG, for_long, &
             &        for_viric_long, pot_viric_long, &
             &        for_virilj_long, pot_virilj_long, &
             &        for_virimor_long, pot_virimor_long, &
             &        for_virish_long, pot_virish_long, &
             &        for_virirfh_long, pot_virirfh_long, &
             &        for_viridou_long, pot_viridou_long, &
             &        for_viricstmnb_long, pot_viricstmnb_long, &
             &        for_viricstmnbex_long, pot_viricstmnbex_long, &
             &        atm_viribo_long, atm_virian_long, &
             &        atm_virito_long, atm_viri14_long, &
             &        pot_virict_long, pot_viriljt_long, &
             &        pot_virimort_long, &
             &        pot_virisht_long, &
             &        pot_virirfht_long, &
             &        pot_viridout_long, &
             &        pot_viricstmnbt_long, &
             &        pot_viricstmnbext_long, &
             &        atm_viribot_long, atm_viriant_long, &
             &        atm_viritot_long, atm_viri14t_long, &
             &        mts_bond, mts_angl, mts_anglub, &
             &        mts_tors, mts_torsrb, mts_torsim, &
             &        mts_vdw, mts_ewr, mts_ewk, &
             &        mts_vdw14, mts_elc14, &
             &        mts_mor, mts_sh, mts_rfh, mts_dou, mts_cnpvw, &
             &        mts_cstmnb, &
             &        mts_posres, &
             &        mts_potbias, &
             &        npoly,npolytyp,npoly_mole,npoly_atom, &
             &        nwater,nmatom,nmatyp,nmatomtyp, &
             &        xcel, ycel, zcel, &
             &        ifewald, alpha, kmax, rrcut, &
             &        ifspme, nfft1, nfft2, nfft3, pme_order, &
             &        eps0, div_factor_14vdw, div_factor_14elc, &
             &        rcut, &
             &        rcutmor, &
             &        rcutsh, &
             &        rcutrfhfo, rcutrfhoo, rcutrfhoh, &
             &        rcutdouo, rcutindouo, rcutdouh, &
             &        rcutrpvw, &
             &        pot_ewc,pot_ewnq, &
             &        pot_tot,pot_nonbon, pot_vdw, pot_elc, pot_ewk, &
             &        pot_vdw14, pot_elc14, &
             &        pot_bond,pot_angl,pot_anglub, &
             &        pot_tors, pot_torsrb, pot_torsim, &
             &        pot_mor, &
             &        pot_sh, &
             &        pot_rfh, &
             &        pot_dou, &
             &        pot_rpvw,pot_vw, &
             &        pot_cstmnb, pot_cstmnbex, &
             &        pot_posres, &
             &        pot_pbias, &
             &        ifcalpremole, ifcalpreatom, &
             &        ifnetqcorrp, &
             &        netchrgsq, &
             &        ifcstmnb, &
             &        ifposres, &
             &        ifpotbias, &
             &        outinterval, pressinterval, recinterval, &
             &        nstep_book, &
             &        md_cont, &
             &        ncstmnbex, &
             &        current_step, &
             &        1,1)
                                ! long-range force

        pot_tot_old = pot_tot

     ELSE

!           --- Store coordinate at T for calculation of constraint force ---
!
!            if (ifrattle) then
!               do i=1,natom
!                  atmcor_att(1,i) = atmcor(1,i)
!                  atmcor_att(2,i) = atmcor(2,i)
!                  atmcor_att(3,i) = atmcor(3,i)
!               end do
!            end if
!
!        --- USUAL CASE (CURRENT_STEP > 0) ---

        call calforce(MTS_LONG, for_long, &
             &        for_viric_long, pot_viric_long, &
             &        for_virilj_long, pot_virilj_long, &
             &        for_virimor_long, pot_virimor_long, &
             &        for_virish_long, pot_virish_long, &
             &        for_virirfh_long, pot_virirfh_long, &
             &        for_viridou_long, pot_viridou_long, &
             &        for_viricstmnb_long, pot_viricstmnb_long, &
             &        for_viricstmnbex_long, pot_viricstmnbex_long, &
             &        atm_viribo_long, atm_virian_long, &
             &        atm_virito_long, atm_viri14_long, &
             &        pot_virict_long, pot_viriljt_long, &
             &        pot_virimort_long, &
             &        pot_virisht_long, &
             &        pot_virirfht_long, &
             &        pot_viridout_long, &
             &        pot_viricstmnbt_long, &
             &        pot_viricstmnbext_long, &
             &        atm_viribot_long, atm_viriant_long, &
             &        atm_viritot_long, atm_viri14t_long, &
             &        mts_bond, mts_angl, mts_anglub, &
             &        mts_tors, mts_torsrb, mts_torsim, &
             &        mts_vdw, mts_ewr, mts_ewk, &
             &        mts_vdw14, mts_elc14, &
             &        mts_mor, mts_sh, mts_rfh, mts_dou, mts_cnpvw, &
             &        mts_cstmnb, &
             &        mts_posres, &
             &        mts_potbias, &
             &        npoly,npolytyp,npoly_mole,npoly_atom, &
             &        nwater,nmatom,nmatyp,nmatomtyp, &
             &        xcel, ycel, zcel, &
             &        ifewald, alpha, kmax, rrcut, &
             &        ifspme, nfft1, nfft2, nfft3, pme_order, &
             &        eps0, div_factor_14vdw, div_factor_14elc, &
             &        rcut, &
             &        rcutmor, &
             &        rcutsh, &
             &        rcutrfhfo, rcutrfhoo, rcutrfhoh, &
             &        rcutdouo, rcutindouo, rcutdouh, &
             &        rcutrpvw, &
             &        pot_ewc,pot_ewnq, &
             &        pot_tot,pot_nonbon, pot_vdw, pot_elc, pot_ewk, &
             &        pot_vdw14, pot_elc14, &
             &        pot_bond,pot_angl,pot_anglub, &
             &        pot_tors, pot_torsrb, pot_torsim, &
             &        pot_mor, &
             &        pot_sh, &
             &        pot_rfh, &
             &        pot_dou, &
             &        pot_rpvw,pot_vw, &
             &        pot_cstmnb, pot_cstmnbex, &
             &        pot_posres, &
             &        pot_pbias, &
             &        ifcalpremole, ifcalpreatom, &
             &        ifnetqcorrp, &
             &        netchrgsq, &
             &        ifcstmnb, &
             &        ifposres, &
             &        ifpotbias, &
             &        outinterval, pressinterval, recinterval, &
             &        nstep_book, &
             &        md_cont, &
             &        ncstmnbex, &
             &        current_step, &
             &        1,1)
                                ! long-range force

!       --- Store old energy data and determine delta r ---

        d_ene = pot_tot - pot_tot_old
        pot_tot_old = pot_tot

        if (d_ene > 0.0d0) then
           d_r = 0.5d0 * d_r   ! d_r is reduced
        else
           d_r = 1.2d0 * d_r   ! d_r is increased
           d_r = min(d_r, d_rmax)
        end if


     END IF

!    --- Calculate root mean square force ---

     frms = 0.0d0

     do i=1,natom
        frms = frms &
             & + for_long(1,i)**2 + for_long(2,i)**2 + for_long(3,i)**2
     end do
     drms = sqrt(frms / dble(natom))  ! drms = sqrt(sum_|f|^2 / natom)
     frms = sqrt(frms)                ! frms = sqrt(sum_|f|^2)

!    --- Check the convergence

     if ((abs(d_ene) < d_econv) .or. (drms < d_rmsf)) then   ! if converged

#if defined(MPI)
        if (irank == 0) then
#endif
           write(6,*) '******************************************'
           write(6,*) 'EM loop: energy is successfully converged.'
           write(6,*) '******************************************'
#if defined(MPI)
        end if
#endif

        minflag = .true.

!!! Time measure
#if defined(TIME_M) || defined(TIME_MALL)
        call get_timecount_barrier(total_time,total_time_elps,'STOP')
#endif
!!! Time measure

!!! Time measure output
#if defined(TIME_M) || defined(TIME_MALL)
        if (current_step /= 0) then
           call out_timecount()
        end if
#endif
!!! Time measure output

        return

     else                                          ! if not converged

!!! Time measure
#if defined(TIME_M)
        call get_timecount_barrier( dyn_time, dyn_time_elps, 'START')
#endif
!!! Time measure

        d_t = d_r / frms

!       --- SHAKE: store old coordinate at T to atmcor_old ---

        if (ifrattle .and. (nconst /= 0)) then
           atmcor_old(1:3,1:natom) = atmcor(1:3,1:natom)
        end if

!       -- UPDATE COORDINATE --

        atmcor(1:3,1:natom) =  atmcor(1:3,1:natom)  &
                &               + d_t * for_long(1:3,1:natom)

!       -- SHAKE_c: reset coordinate --

        if (ifrattle .and. (nconst /= 0)) then

           call shake_c(atmcor_old, &
                &       eps_rattle)

        end if

!       -- SHAKE_c_z: fix z-coordinate --

        if (iflocalfixz) then

           call shake_c_z()

        end if

!       -- SHAKE_c_zg: fix z-coordinate of COM --

        if (iflocalfixzg) then

           call shake_c_zg()

        end if

!!! Time measure
#if defined(TIME_M)
        call get_timecount( dyn_time, dyn_time_elps, 'STOP' )
#endif
!!! Time measure

!       --- setting velocity of VW to 0 when calculating temperature

        if (ifcnp) then
           do j = 1,2
              i = nvwlist(j)
              atmvel_tmp(1:3,i) = atmvel(1:3,i)
              atmvel(1:3,i)     = 0.0d0
           end do
        end if

!       --- CALCULATE PRESSURE OF SYSTEM ---

!!! Time measure
#if defined(TIME_M)
        call get_timecount_barrier( dyn_time, dyn_time_elps, 'START' )
#endif
!!! Time measure

        if ((mod(current_step,pressinterval) == 1) .or. &
             & (pressinterval == 1)) then

           if (ifcalpremole .or. ifcalpreatom) then
              call calpress(for_viric_long,for_viric_med,for_viric_short, &
                   &        pot_viric_long,pot_viric_med,pot_viric_short, &
                   &        pot_virict_long,pot_virict_med, &
                   &        pot_virict_short, &
                   &        for_virilj_long,for_virilj_med, &
                   &        for_virilj_short, &
                   &        pot_virilj_long,pot_virilj_med, &
                   &        pot_virilj_short, &
                   &        pot_viriljt_long,pot_viriljt_med, &
                   &        pot_viriljt_short, &
                   &        for_virimor_long,for_virimor_med, &
                   &        for_virimor_short, &
                   &        pot_virimor_long,pot_virimor_med, &
                   &        pot_virimor_short, &
                   &        pot_virimort_long,pot_virimort_med, &
                   &        pot_virimort_short, &
                   &        for_virish_long,for_virish_med, &
                   &        for_virish_short, &
                   &        pot_virish_long,pot_virish_med, &
                   &        pot_virish_short, &
                   &        pot_virisht_long,pot_virisht_med, &
                   &        pot_virisht_short, &
                   &        for_virirfh_long,for_virirfh_med, &
                   &        for_virirfh_short, &
                   &        pot_virirfh_long,pot_virirfh_med, &
                   &        pot_virirfh_short, &
                   &        pot_virirfht_long,pot_virirfht_med, &
                   &        pot_virirfht_short, &
                   &        for_viridou_long,for_viridou_med, &
                   &        for_viridou_short, &
                   &        pot_viridou_long,pot_viridou_med, &
                   &        pot_viridou_short, &
                   &        pot_viridout_long,pot_viridout_med, &
                   &        pot_viridout_short, &
                   &        for_viricstmnb_long,for_viricstmnb_med,   &
                   &        for_viricstmnb_short,   &
                   &        pot_viricstmnb_long,pot_viricstmnb_med,   &
                   &        pot_viricstmnb_short,   &
                   &        pot_viricstmnbt_long,pot_viricstmnbt_med,   &
                   &        pot_viricstmnbt_short,   &
                   &        for_viricstmnbex_long,for_viricstmnbex_med,   &
                   &        for_viricstmnbex_short,   &
                   &        pot_viricstmnbex_long,pot_viricstmnbex_med,   &
                   &        pot_viricstmnbex_short,   &
                   &        pot_viricstmnbext_long,pot_viricstmnbext_med,   &
                   &        pot_viricstmnbext_short,   &
                   &        atm_viribo_long,atm_viribo_med, &
                   &        atm_viribo_short, &
                   &        atm_viribot_long,atm_viribot_med, &
                   &        atm_viribot_short, &
                   &        atm_virian_long,atm_virian_med, &
                   &        atm_virian_short, &
                   &        atm_viriant_long,atm_viriant_med, &
                   &        atm_viriant_short, &
                   &        atm_virito_long,atm_virito_med, &
                   &        atm_virito_short, &
                   &        atm_viritot_long,atm_viritot_med, &
                   &        atm_viritot_short, &
                   &        atm_viri14_long,atm_viri14_med, &
                   &        atm_viri14_short, &
                   &        atm_viri14t_long,atm_viri14t_med, &
                   &        atm_viri14t_short, &
                   &        atm_viri_const,atm_virit_const, &
                   &        atm_viri_corr,atm_virit_corr, &
                   &        pressmol_ktot,pressmol_vtot,pressmol_tot, &
                   &        pressatm_ktot,pressatm_vtot,pressatm_tot, &
                   &        pot_viric_all,pot_virilj_all, &
                   &        pot_virimor_all, &
                   &        pot_virish_all, &
                   &        pot_virirfh_all, &
                   &        pot_viridou_all, &
                   &        pot_viricstmnb_all,   &
                   &        pot_viricstmnbex_all,   &
                   &        pot_viri_all, &
                   &        viri_fdotd, &
                   &        atm_viribo_all,atm_virian_all, &
                   &        atm_virito_all,atm_viri14_all, &
                   &        pressmolt_ktot,pressmolt_vtot, &
                   &        pressmolt_tot, &
                   &        pressatmt_ktot,pressatmt_vtot, &
                   &        pressatmt_tot, &
                   &        pot_virict_all,pot_viriljt_all, &
                   &        pot_virimort_all, &
                   &        pot_virisht_all, &
                   &        pot_virirfht_all, &
                   &        pot_viridout_all, &
                   &        pot_viricstmnbt_all,   &
                   &        pot_viricstmnbext_all,   &
                   &        pot_virit_all, &
                   &        virit_fdotd, &
                   &        atm_viribot_all,atm_viriant_all, &
                   &        atm_viritot_all,atm_viri14t_all, &
                   &        current_step,pressinterval, &
                   &        xcel,ycel,zcel, &
                   &        ifcalpremole,ifcalpreatom, &
                   &        pot_ewc,pref,eref, &
                   &        rcut, &
                   &        ifcalljlong,nsolve,solveindex, &
                   &        pot_ewnq, &
                   &        ncstmnbex)
           end if

        end if

!       --- Substitute pressmol & pressatm to pint ---

        if (ifpmolcont) then
           pint = pressmol_tot
           pintt(1:3,1:3) = pressmolt_tot(1:3,1:3)
        else
           pint = pressatm_tot
           pintt(1:3,1:3) = pressatmt_tot(1:3,1:3)
        end if

!       - local fix
        if (iflocalfix) then
           do ilfix = 1, nlfix
              i = index_nlfix(ilfix)
              atmvel(1:3,i) = 0.0d0
           end do
        end if

!!! Time measure
#if defined(TIME_M)
        call get_timecount( dyn_time, dyn_time_elps, 'STOP')
#endif
!!! Time measure

!       --- restoration velocity of VW
        if (ifcnp) then
           do j = 1,2
              i = nvwlist(j)
              atmvel(1:3,i) = atmvel_tmp(1:3,i)
           end do
        end if

!       --- Record the state of the end of this EM timestep ---

!!! Time measure
#if defined(TIME_M)
        call get_timecount_barrier( out_time, out_time_elps, 'START' )
#endif
!!! Time measure

#if defined(MPI)
        if (irank == 0) then
#endif
           if (mod(current_step,recinterval) == 0) then
              call wrsta(iostarec, &
                   &     npoly, nwater, nmatom, &
                   &     xcel, ycel, zcel, &
                   &     xref, vref, timeref, pref, &
                   &     current_step, &
                   &     mchain, &
                   &     pint, pintt)
           end if
#if defined(MPI)
        end if
#endif
!

!!! Time measure
#if defined(TIME_M)
        call get_timecount( out_time, out_time_elps, 'STOP')
#endif
!!! Time measure

!       --- CALCULATE KINETIC ENERGY & TOTAL ENERGY AT T + dT ---

!!! Time measure
#if defined(TIME_M)
        call get_timecount_barrier( out_time, out_time_elps, 'START')
#endif
!!! Time measure

!        do iwrite = 1, natom
!           write(6,*) atmvel(1,iwrite),atmvel(2,iwrite),atmvel(3,iwrite)
!        end do

        call calkin( npoly, npolytyp, npoly_mole, &
             &       nwater, &
             &       nmatom, nmatyp, nmatomtyp, &
             &       degfree_poly, degfree_water, degfree_ma, &
             &       ene_kin_poly, ene_kin_water, ene_kin_ma, &
             &       temp_poly, temp_water, temp_ma, &
             &       degfree_all, ene_kin_all, temp_all, &
             &       mchain,text_c, &
             &       ene_kin_th, ene_pot_th, &
             &       pext_c, &
             &       ene_kin_ba, ene_pot_ba, extra_pot_ba)

!            ene_tot = ene_kin_poly + ene_kin_water + ene_kin_ma + pot_tot
        ene_tot = ene_kin_water + pot_tot
        do i = 1, npolytyp
           ene_tot = ene_tot + ene_kin_poly(i)
        end do
        do i = 1, nmatyp
           ene_tot = ene_tot + ene_kin_ma(i)
        end do
        ene_conserve = ene_tot

!           --- Output energy, trajectory, velocity etc. ---
#if defined(MPI)
        if (irank == 0) then
#endif

!           *** in EM, each data is always output
!           if ((mod(current_step,outinterval) == 1) .or. &
!                & (outinterval == 1)) then

           call outene_em(ouene,current_step,xref,eref,tempref,fref, &
                &         npolytyp,nmatyp, &
                &         pot_tot,pot_nonbon,pot_vdw, &
                &         pot_elc,pot_ewk,pot_ewc,pot_ewnq, &
                &         pot_vdw14,pot_elc14, &
                &         pot_bond,pot_angl,pot_anglub, &
                &         pot_tors,pot_torsrb,pot_torsim, &
                &         pot_mor, &
                &         pot_sh, &
                &         pot_rfh, &
                &         pot_dou, &
                &         pot_rpvw,pot_vw, &
                &         pot_cstmnb, &
                &         ncstmnbex,pot_cstmnbex, &
                &         pot_posres, &
                &         pot_pbias, &
                &         ene_kin_poly,ene_kin_water, &
                &         ene_kin_ma,ene_kin_all, &
                &         ene_tot, &
                &         temp_poly,temp_water,temp_ma,temp_all, &
                &         ene_kin_th,ene_pot_th, &
                &         ene_kin_ba,ene_pot_ba, &
                &         ene_conserve, &
                &         d_r,d_ene,drms)

              if (ifoutpos) then
                 call outpos(oupos, current_step, xref)
              end if

              if (ifoutvel) then
                 call outvel(ouvel, current_step, vref)
              end if

              if (ifoutfor) then
                  call outfor(oufor,current_step,fref, &
                  &           for_long,for_short)
              end if

              ! if (ifoutthe) then
              !    call outthe(outhe, current_step, timeref, &
              !         &      mchain, xlogs, vlogs)
              ! end if
              !
              ! if (ifoutbar) then
              !     call outbar(oubar,current_step,xref,timeref, &
              !     &           xcel,ycel,zcel,xlogv,vlogv,xboxh,vboxg, &
              !     &           pcont_axis)
              ! end if

!           end if

           if ((mod(current_step,pressinterval) == 1) .or. &
                & (pressinterval == 1)) then

              if (ifoutpre) then

                 if (ifcalpremole .or. ifcalpreatom) then
                    call outpre(oupre,current_step,pref,eref, &
                         &      pressmol_ktot,pressmolt_ktot, &
                         &      pressmol_vtot, pressmolt_vtot, &
                         &      pressmol_tot,pressmolt_tot, &
                         &      pressatm_ktot,pressatmt_ktot, &
                         &      pressatm_vtot,pressatmt_vtot, &
                         &      pressatm_tot,pressatmt_tot, &
                         &      pot_viric_all,pot_virict_all, &
                         &      pot_virilj_all,pot_viriljt_all, &
                         &      pot_virimor_all,pot_virimort_all, &
                         &      pot_virish_all,pot_virisht_all, &
                         &      pot_virirfh_all,pot_virirfht_all, &
                         &      pot_viridou_all,pot_viridout_all, &
                         &      pot_viricstmnb_all,pot_viricstmnbt_all, &
                         &      ncstmnbex, &
                         &      pot_viricstmnbex_all,pot_viricstmnbext_all, &
                         &      pot_viri_all,pot_virit_all, &
                         &      viri_fdotd,virit_fdotd, &
                         &      atm_viribo_all,atm_viribot_all, &
                         &      atm_virian_all,atm_viriant_all, &
                         &      atm_virito_all,atm_viritot_all, &
                         &      atm_viri14_all,atm_viri14t_all, &
                         &      atm_viri_const,atm_virit_const, &
                         &      atm_viri_corr,atm_virit_corr)
                 end if

              end if

           end if

#if defined(MPI)
        end if
#endif

!!! Time measure
#if defined(TIME_M)
        call get_timecount(out_time,out_time_elps,'STOP')
#endif
!!! Time measure

!!! Time measure
#if defined(TIME_M)
        call get_timecount_barrier(out_time,out_time_elps,'START')
#endif
!!! Time measure

#if defined(MPI)
        if (irank == 0) then
#endif

           if (ifoutpdb .and. (current_step == nstep_pdbout)) then
              call outpdb( oupdb,current_step, &
                   &       npolytyp,npoly_mole,npoly_atom, &
                   &       nwater, &
                   &       nmatyp,nmatomtyp, &
                   &       resname_poly_pdb, &
                   &       resname_water_pdb, &
                   &       resname_matom_pdb )
           end if

           if (ifpotbias .and. (current_step == 0)) then
              call outumb( ouumb, &
                   &       xref, eref )
           endif

#if defined(MPI)
        end if
#endif

!!! Time measure
#if defined(TIME_M)
        call get_timecount(out_time,out_time_elps,'STOP')
#endif
!!! Time measure

!!! Time measure
#if defined(TIME_M) || defined(TIME_MALL)
        call get_timecount_barrier(total_time,total_time_elps,'STOP')
#endif
!!! Time measure

!!! Time measure output
#if defined(TIME_M) || defined(TIME_MALL)
        if (current_step /= 0) then
           call out_timecount()
        end if
#endif
!!! Time measure output

     end if

  END DO              ! end of EM loop

!----- mts flag is restored
  mts_bond = mts_bond_old
  mts_angl = mts_angl_old
  mts_anglub = mts_anglub_old
  mts_tors = mts_tors_old
  mts_torsrb = mts_torsrb_old
  mts_torsim = mts_torsim_old
  mts_vdw = mts_vdw_old
  mts_ewr = mts_ewr_old
  mts_ewk = mts_ewk_old
  mts_vdw14 = mts_vdw14_old
  mts_elc14 = mts_elc14_old
  mts_mor = mts_mor_old
  mts_sh = mts_sh_old
  mts_rfh = mts_rfh_old
  mts_dou = mts_dou_old
  mts_cnpvw = mts_cnpvw_old
  mts_cstmnb = mts_cstmnb_old
  mts_posres = mts_posres_old
  mts_potbias = mts_potbias_old

#if defined(MPI)
  if (irank == 0) then
#endif
     write(6,*) '******************************************'
     write(6,*) 'Warning: In EM, energy is not converged.  '
     write(6,*) '******************************************'
#if defined(MPI)
  end if
#endif


!---- release memory
  deallocate(pot_cstmnbex)
  deallocate(for_viricstmnbex_long)
  ! deallocate(for_viricstmnbex_med)
  deallocate(for_viricstmnbex_short)
  deallocate(pot_viricstmnbex_long)
  ! deallocate(pot_viricstmnbex_med)
  deallocate(pot_viricstmnbex_short)
  deallocate(pot_viricstmnbext_long)
  ! deallocate(pot_viricstmnbext_med)
  deallocate(pot_viricstmnbext_short)
  deallocate(pot_viricstmnbex_all)
  deallocate(pot_viricstmnbext_all)

!     +     +     +     +     +     +     +

end subroutine enemin_sd
