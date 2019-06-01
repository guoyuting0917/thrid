!********************************
!*  moldyn_nhc.f90 Ver.6.7      *
!*      for peachgk_md.f        *
!*      MD pratform subroutine  *
!*            by G.Kikugawa     *
!********************************
! Time-stamp: <>

subroutine moldyn_nhc(iostarec, recinterval, &
     &                npoly, npolytyp, npoly_mole, npoly_atom, &
     &                nwater, &
     &                nmatom, nmatyp, nmatomtyp, &
     &                maxnstep, inistep, endstep, &
     &                xcel, ycel, zcel, &
     &                ifcenterfix_all, &
     &                ifcenterfix_poly, &
     &                ifcenterfix_water, &
     &                ifcenterfix_ma, &
     &                cenfix_free, &
     &                ifcenterfix_polytyp, &
     &                ifcenterfix_watertyp, &
     &                ifcenterfix_matyp, &
     &                dt_long_cal, &
     &                nstep_med, nstep_short, &
     &                nstep_vir, &
     &                iflimitmove,limitdist, &
     &                mts_bond, mts_angl, mts_anglub, &
     &                mts_tors, mts_torsrb, mts_torsim, &
     &                mts_vdw, mts_ewr, mts_ewk, mts_vdw14, mts_elc14, &
     &                mts_mor, mts_sh, mts_rfh, mts_dou, mts_cnpvw, &
     &                mts_cstmnb, &
     &                mts_posres, &
     &                mts_potbias, &
     &                ifewald, alpha, kmax, rrcut, &
     &                ifspme, nfft1, nfft2, nfft3, pme_order, &
     &                eps0, div_factor_14vdw, div_factor_14elc, &
     &                rcut, &
     &                ifcellindex, &
     &                ifbook, nstep_book, &
     &                tcont_poly, tcont_water, tcont_ma, &
     &                ifrattle, eps_rattle, &
     &                iftcratom, &
     &                iflocalfix, &
     &                nlfix, index_nlfix, &
     &                iflocalfixz, iflocalfixzg, &
     &                ifcnp, &
     &                ifposres, &
     &                ifpotbias, &
     &                iflocalvel, ifstrmvel, &
     &                nlvel, index_nlvel, v_nlvel, &
     &                ifoutpdb, nstep_pdbout, &
     &                resname_poly_pdb, &
     &                resname_water_pdb, &
     &                resname_matom_pdb, &
     &                xref, eref, mref, qref, &
     &                vref, timeref, tempref, pref, fref, eps0ref, &
     &                degfree_poly, degfree_water, &
     &                degfree_ma, degfree_all, &
     &                dt_med_cal, dt_short_cal, &
     &                rcut_book, &
     &                rcutmor, ifbookmor, rcut_bookmor, nstep_bookmor, &
     &                ifcellindex_mor, &
     &                rcutsh, ifbooksh, rcut_booksh, nstep_booksh, &
     &                ifcellindex_sh, &
     &                rcutrfhfo, ifbookrfhfo, rcut_bookrfhfo, &
     &                nstep_bookrfhfo, &
     &                ifcellindex_rfhfo, &
     &                rcutrfhoo, ifbookrfhoo, rcut_bookrfhoo, &
     &                nstep_bookrfhoo, &
     &                ifcellindex_rfhoo, &
     &                rcutrfhoh, ifbookrfhoh, rcut_bookrfhoh, &
     &                nstep_bookrfhoh, &
     &                ifcellindex_rfhoh, &
     &                rcutdouo, rcutindouo, ifbookdouo, rcut_bookdouo, &
     &                nstep_bookdouo, &
     &                ifcellindex_douo, &
     &                rcutdouh, ifbookdouh, rcut_bookdouh, &
     &                nstep_bookdouh, &
     &                ifcellindex_douh, &
     &                rcutrpvw, ifbookrpvw, rcut_bookrpvw, &
     &                nstep_bookrpvw, &
     &                ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
     &                for_long, for_short, &
     &                for_viric_long, for_viric_med, for_viric_short, &
     &                for_virilj_long, for_virilj_med, for_virilj_short, &
     &                for_virimor_long, for_virimor_med, &
     &                for_virimor_short, &
     &                for_virish_long, for_virish_med, &
     &                for_virish_short, &
     &                for_virirfh_long, for_virirfh_med, &
     &                for_virirfh_short, &
     &                for_viridou_long, for_viridou_med, &
     &                for_viridou_short, &
     &                for_viricstmnb_long, for_viricstmnb_med, &
     &                for_viricstmnb_short, &
     &                pot_ewc, &
     &                ouene,oupos,ouvel,oufor, &
     &                outhe,oubar,oupre,oupdb, &
     &                ouumb, &
     &                ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
     &                ifoutbar,ifoutpre, &
     &                tcontinterval, outinterval, pressinterval, &
     &                yratio, zratio, &
     &                ifcalpremole, ifcalpreatom, &
     &                ifnetqcorrp, &
     &                mchain, text, &
     &                pint, pintt, &
     &                ifpatmcont, ifpmolcont, &
     &                next, nyosh, &
     &                ifcalljlong, nsolve, solveindex, &
     &                netchrgsq, &
     &                nstep_expand, r_expand, &
     &                md_cont)

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

  real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]

  integer,intent(in):: nstep_med       ! number of step for medium force
  integer,intent(in):: nstep_short     ! number of step for short force

  integer,intent(in):: nstep_vir    ! number of step for virtual time integration

  logical,intent(in):: iflimitmove ! if limit atomic motion to a specific distance
  real(8),intent(in):: limitdist  ! maximum atomic displacement when doing
                              ! time integration [non-d] (structure relaxation)

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

  real(8),intent(in):: tcont_poly(:)    ! poly Temp. [non-d] in NVT
  real(8),intent(in):: tcont_water      ! H2O Temp. [non-d] in NVT
  real(8),intent(in):: tcont_ma(:)  ! monatomic mole. Temp. [K] in NVT (Woodcock)

  logical,intent(in):: ifrattle        ! rattle flag

  real(8),intent(in):: rcutmor          ! Morse cutoff length [m]

  logical,intent(in):: ifcellindex_mor ! flag for cell index (morse)

  logical,intent(in):: ifbookmor      ! flag for bookkeeping of Morse interaction
  real(8),intent(in):: rcut_bookmor   ! cut off radius of bookkeeping[m] of Morse
  integer,intent(in):: nstep_bookmor  ! bookkeeping interval of Morse interaction

  real(8),intent(in):: rcutsh           ! SH cutoff length [m]

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

  logical,intent(in):: ifbookrfhoh
                         ! flag for bookkeeping of RFH(OH) interaction
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

  !---- variable for controling normal pressure
  logical,intent(in):: ifcnp           ! flag to control normal pressure

!---- variable for position restraint
  logical,intent(in):: ifposres        ! position restraint flag

!---- variable for bias potential
  logical,intent(in):: ifpotbias       ! bias potential flag

!---- variable for local velocity
  logical,intent(in):: iflocalvel      ! flag to force local atomic velocity

  integer,intent(in):: nlvel           ! number of atoms for velocity fix
  integer,intent(in):: index_nlvel(:)  ! index of vel-fix atoms
  real(8),intent(in):: v_nlvel(:,:)    ! local velocity values

!---- variable for streaming velocity
  logical,intent(in):: iftcratom    ! region temp. control based on atom or mole.

  logical,intent(in):: ifstrmvel     ! flag to input and use streaming velocity

!---- PDB output
  logical,intent(in):: ifoutpdb        ! flag for outputting PDB format file
  integer,intent(in):: nstep_pdbout    ! MD step for output of PDB file
  character(4),intent(in):: resname_poly_pdb(:) ! residue name for poly
  character(4),intent(in):: resname_water_pdb ! residue name for water
  character(4),intent(in):: resname_matom_pdb(:) ! residue name for matom

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
  integer,intent(in):: oupre           ! output unit for outpre pressure data
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

!---- time step
  real(8),intent(in):: dt_med_cal       ! time step of medium force
  real(8),intent(in):: dt_short_cal     ! time step of short force

!---- radius of bookkeeping
  real(8),intent(in):: rcut_book        ! cut off radius of bookkeeping

!---- force of each atoms
  real(8),intent(inout):: for_long(:,:)  ! long-range force
!      real(8),intent(inout):: for_med(:,:)   ! medium-range force
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

  real(8):: ene_tot          ! total Hamiltonian
  real(8):: ene_conserve     ! total conserved quantity

!---- temperature

  real(8):: temp_poly(maxnpolytyp) ! temperature of polymer1
  real(8):: temp_water       ! temperature of H2O
  real(8):: temp_ma(maxnmatyp) ! temperature of monatomic mole.
  real(8):: temp_all         ! temperature of all molecules

!---- interval of temperature control
  integer,intent(in):: tcontinterval   ! interval of temp. control

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

  real(8):: pot_viric_all            ! virial(coulomb potential)
  real(8):: pot_virict_all(3,3)      ! virial tensor (coulomb potential)
  real(8):: pot_virilj_all           ! virial(L-J potential)
  real(8):: pot_viriljt_all(3,3)     ! virial tensor (L-J potential)
  real(8):: pot_virimor_all          ! virial(Morse potential)
  real(8):: pot_virimort_all(3,3)    ! virial tensor (Morse potential)
  real(8):: pot_virish_all           ! virial(SH potential)
  real(8):: pot_virisht_all(3,3)     ! virial tensor (SH potential)
  real(8):: pot_virirfh_all          ! virial(RFH potential)
  real(8):: pot_virirfht_all(3,3)    ! virial tensor (RFH potential)
  real(8):: pot_viridou_all          ! virial(DOU potential)
  real(8):: pot_viridout_all(3,3)    ! virial tensor (DOU potential)
  real(8):: pot_viricstmnb_all       ! virial(custom NB potential)
  real(8):: pot_viricstmnbt_all(3,3) ! virial tensor (custom NB potential)
  real(8),allocatable:: pot_viricstmnbex_all(:)
                                     ! virial(extra custom NB potential)
  real(8),allocatable:: pot_viricstmnbext_all(:,:,:)
                                     ! virial tensor (extra custom NB potential)
  real(8):: pot_viri_all             ! all virial of potential term
  real(8):: pot_virit_all(3,3)       ! all virial tensor of potential term
  real(8):: viri_fdotd               ! virial of correction term F.d
  real(8):: virit_fdotd(3,3)         ! virial tensor of correction term F.d

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
  real(8),intent(in):: text             ! external temp. [K] (Nose-Hoover chain)

!---- variables for Andersen (Hoover type) barostat
  real(8),intent(inout):: pint         ! internal pressure
  real(8),intent(inout):: pintt(:,:)   ! internal pressure tensor
  logical,intent(in):: ifpatmcont      ! atomic pressure control
  logical,intent(in):: ifpmolcont      ! molecular pressure control

!---- for higher order Trotter expansion
  integer,intent(in):: next            ! iteration number of extended system
  integer,intent(in):: nyosh           ! expansion order of Yoshida-Suzuki method

!---- pressure calculation of L-J long-range correction
  logical,intent(in):: ifcalljlong     ! long-range correction in pressure
  integer,intent(in):: nsolve          ! number of solvent molecules
!      real*8:: vdw_welij_solve ! well depth of vdw parameter of solvent
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
  integer:: i,j             ! do loop index

  real(8):: atmcor_old(3,maxnatom) ! atmcor at time T for rattle-rRESPA
!      real*8:: atmcor_att(3,maxnatom) ! atmcor at time T for const. force

  integer:: ilfix
  integer:: ilvel

  real(8):: atmvel_tmp(3,maxnatom) ! for use of setting velocity of VW to 0

!     temporary temp.
  real(8):: text_c           ! external temp. in NVT(NHC)

  real(8):: pext_c = 0.0d0

  real(8):: expcdt2 = 0.0d0  ! = exp(-dt/2*vxi(1))
  real(8):: expcdt4 = 0.0d0  ! = exp(-dt/4*vxi(1))

  real(8):: e2               ! = 1/3!
  real(8):: e4               ! = 1/5!
  real(8):: e6               ! = 1/7!
  real(8):: e8               ! = 1/9!

  real(8):: argv2            ! = (dt_long/4*vxi(1))^2
  real(8):: polyv            ! = sinh(dt_long/4*vxi(1))/(dt_long/4*vxi(1))
  real(8):: bbv = 0.0d0      ! = polyv*dt_long/2*exp(-dt_long/4*vxi(1))
  real(8):: cf               ! = polyv*exp(-dt_long/4*vxi(1))

!     +     +     +     +     +     +     +

!----- some initialization -----

  atm_viri_const = 0.0d0
  atm_virit_const(1:3,1:3) = 0.0d0

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

!----- some preparation -----

  e2 = 1.0d0/ 6.0d0         ! = 1/3!
  e4 = e2   /20.0d0         ! = 1/5!
  e6 = e4   /42.0d0         ! = 1/7!
  e8 = e6   /72.0d0         ! = 1/9!

  text_c = text

!-------- MD main loop --------

  DO current_step = inistep, endstep

!!! Time measure initialization
#if defined(TIME_M) || defined(TIME_MALL)
     call init_timecount()
#endif
!!! Time measure initialization

!!! Time measure
#if defined(TIME_M) || defined(TIME_MALL)
     call get_timecount_barrier(total_time,total_time_elps,'START')
#endif
!!! Time measure

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier(dyn_time,dyn_time_elps,'START')
#endif
!!! Time measure

!---- translate the coordinates (P.B.C.)

     if ((mod(current_step,nstep_book) == 1) .or. &
          & (nstep_book == 1) .or. (current_step == 0)) then

        call transcor( npoly, nwater, nmatom, &
             &         xcel, ycel, zcel )

     end if

!---- Expand cell volume ---
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
!            if ((mod(current_step,tcontinterval) == 1) .or.
!               &                  (tcontinterval == 1)) then
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
!            end if
     end if

!!! Time measure
#if defined(TIME_M)
     call get_timecount( dyn_time, dyn_time_elps, 'STOP' )
#endif
!!! Time measure

!---- make list of non-bonded pair list

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier( table_time, table_time_elps, 'START' )
#endif
!!! Time measure

#if !defined(MDG3) && !defined(MDG_m3)
     if ((ifbook) .or. (current_step == 0)) then
        if ((mod(current_step,nstep_book) == 0) .or. &
             & (nstep_book == 1)) then

           if (ifcellindex) then
              call mklist_cell( rcut_book, &
                   &            xcel, ycel, zcel )
           else
              call mklist2a( rcut_book, &
                   &         xcel, ycel, zcel, &
                   &         ifbook )
           end if

        end if
     end if
#endif

!---- make list of Morse bonded pair list
     if ((ifbookmor) .or. (current_step == 0)) then
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
     if ((ifbooksh) .or. (current_step == 0)) then
        if ((mod(current_step,nstep_booksh) == 0) .or. &
             & (nstep_booksh == 1)) then

           if (ifcellindex_sh) then
              call mklist_sh_cell( rcut_booksh, &
                   &               xcel, ycel, zcel )
           else
              call mklist2a_sh( rcut_booksh, &
                   &            xcel, ycel, zcel, &
                   &            ifbooksh )
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
                   &                xcel,ycel,zcel, &
                   &                INTTYPE_RFHFO )
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
                   &             xcel, ycel, zcel, &
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
                   &                INTTYPE_RFHOH)
           else
              call mklist2a_rfh( rcut_bookrfhoh, &
                   &             xcel,ycel,zcel, &
                   &             ifbookrfhoh, &
                   &             INTTYPE_RFHOH)
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
                   &                xcel,ycel,zcel, &
                   &                INTTYPE_DOUO)
           else
              call mklist2a_dou( rcut_bookdouo, &
                   &             xcel,ycel,zcel, &
                   &             ifbookdouo, &
                   &             INTTYPE_DOUO)
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
                   &                INTTYPE_DOUH)
           else
              call mklist2a_dou( rcut_bookdouh, &
                   &             xcel, ycel, zcel, &
                   &             ifbookdouh, &
                   &             INTTYPE_DOUH)
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
     call get_timecount(table_time,table_time_elps,'STOP')
#endif
!!! Time measure

!---- NVT-MTS (XO-RESPA)

!  --- CALCULATE FORCE AT T (CURRENT_TIME) IF CURRENT_STEP = 1 ---
!         * otherwise, the accell at time T has been already calculated

     IF (current_step == 0) THEN

!       -- assign streaming velocity
        call ass_strmvel(npoly,npolytyp, &
             &           npoly_mole,npoly_atom, &
             &           nwater, &
             &           nmatom,nmatyp,nmatomtyp, &
             &           xcel,ycel,zcel, &
             &           iftcratom, &
             &           ifstrmvel)

!           - calculate initial force & energy -

        call calforce(MTS_SHORT,for_short, &
             &        for_viric_short, pot_viric_short, &
             &        for_virilj_short, pot_virilj_short, &
             &        for_virimor_short, pot_virimor_short, &
             &        for_virish_short, pot_virish_short, &
             &        for_virirfh_short, pot_virirfh_short, &
             &        for_viridou_short, pot_viridou_short, &
             &        for_viricstmnb_short, pot_viricstmnb_short, &
             &        for_viricstmnbex_short, pot_viricstmnbex_short, &
             &        atm_viribo_short, atm_virian_short, &
             &        atm_virito_short, atm_viri14_short, &
             &        pot_virict_short, pot_viriljt_short, &
             &        pot_virimort_short, &
             &        pot_virisht_short, &
             &        pot_virirfht_short, &
             &        pot_viridout_short, &
             &        pot_viricstmnbt_short, &
             &        pot_viricstmnbext_short, &
             &        atm_viribot_short, atm_viriant_short, &
             &        atm_viritot_short, atm_viri14t_short, &
             &        mts_bond, mts_angl, mts_anglub, &
             &        mts_tors, mts_torsrb, mts_torsim, &
             &        mts_vdw,mts_ewr, mts_ewk, &
             &        mts_vdw14, mts_elc14, &
             &        mts_mor, mts_sh,mts_rfh, mts_dou, mts_cnpvw, &
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
             &        pot_tot,pot_nonbon, pot_vdw,pot_elc, pot_ewk, &
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
             &        nstep_short, istep_short)
                                ! short-range force

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
             &        nstep_short,istep_short)
                                ! long-range force

        call accforce( natom, atmmass, for_short )
        call accforce( natom, atmmass, for_long )


     ELSE

!        --- USUAL CASE (CURRENT_STEP > 0) ---

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
!           --- UPDATE Nose-Hoover chain ---

!!! Time measure
#if defined(TIME_M)
        call get_timecount_barrier( dyn_time, dyn_time_elps, 'START' )
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

        call nhcint( degfree_all, &
             &       dt_long_cal, &
             &       mchain,text_c, &
             &       next,nyosh )

!       --- restoration velocity of VW
        if (ifcnp) then
           do j = 1,2
              i = nvwlist(j)
              atmvel(1:3,i) = atmvel_tmp(1:3,i)
           end do
        end if

!!! Time measure
#if defined(TIME_M)
        call get_timecount( dyn_time, dyn_time_elps, 'STOP' )
#endif
!!! Time measure

!           --- LOOP OVER NVT-XO-RESPA ---

        DO istep_short = 1, nstep_short

!              --- UPDATE VELOCITY BY FOR_LONG ---

!!! Time measure
#if defined(TIME_M)
           call get_timecount_barrier( dyn_time, dyn_time_elps, 'START' )
#endif
!!! Time measure

           if (istep_short == 1) then

!             --- setting velocity of VW to 0
              if (ifcnp) then
                 do j = 1,2
                    i = nvwlist(j)
                    atmvel_tmp(1:3,i) = atmvel(1:3,i)
                    atmvel(1:3,i)     = 0.0d0
                 end do
              end if

!             --- correct velocity with streaming velocity
              atmvel(1:3,1:natom) = atmvel(1:3,1:natom) &
                   &              - atmvel_strm(1:3,1:natom)

!                 calculate coefficient of velocity(for_long) update
              expcdt4 = exp(-0.25d0*dt_long_cal*vlogs(1))
                                ! = exp(-dt/4*vxi(1))
              expcdt2 = expcdt4 * expcdt4
                                ! = exp(-dt/2*vxi(1))
              argv2 =  (0.25d0*dt_long_cal*vlogs(1)) &
                   & * (0.25d0*dt_long_cal*vlogs(1))
                                ! = (dt_long/4*vxi(1))^2
              polyv =  1.0d0 &
                   & + argv2*(e2 + argv2 &
                   &          *(e4 + argv2*(e6 + argv2*e8)))
                                ! = sinh(dt_long/4*vxi(1))/(dt_long/2*vxi(1))
              cf = expcdt4 * polyv
                                ! = polyv*exp(-dt_long/4*vxi(1))
              bbv = cf * 0.5d0 * dt_long_cal
                                ! = polyv*dt_long/2*exp(-dt_long/4*vxi(1))

              atmvel(1:3,1:natom) = atmvel(1:3,1:natom)*expcdt2 &
                   &              + for_long(1:3,1:natom)*bbv
                                ! higher factorized version
                                ! vi' = vi*exp(-dt/2*vxi(1))
                                !      +Fi/mi*bbv

!             --- reset velocity
              atmvel(1:3,1:natom) = atmvel(1:3,1:natom) &
                   &              + atmvel_strm(1:3,1:natom)

!             --- restoration velocity of VW
              if (ifcnp) then
                 do j = 1,2
                    i = nvwlist(j)
                    atmvel(1:3,i) = atmvel_tmp(1:3,i)
                 end do
              end if

!                     do i=1,natom
!                        do j=1,3
!                           atmvel(j,i) =  atmvel(j,i)
!     &                                  * dexp(-0.5d0*dt_long_cal
!     &                                         *vlogs(1))
!     &                                  + 0.5d0*dt_long_cal
!     &                                    *for_long(j,i)
!                                ! lower factorized version
!                                ! vi' = vi*exp(-dt/2*vxi(1))
!                                !      +dt/2*Fi/mi

!                           atmvel(j,i) =  atmvel(j,i)
!     &                                  * dexp(-0.5d0*dt_long_cal
!     &                                         *vlogs(1))
!     &                                  + 0.5d0*dt_long_cal*for_long(j,i)
!     &                                  * dexp(-0.25d0*dt_long_cal
!     &                                         *vlogs(1))
!                                ! lower factorized version
!                                ! vi' = vi*exp(-dt/2*vxi(1))
!                                !      +dt/2*Fi/mi*exp(-dt/4*vxi(1))
!                        end do
!                     end do

           end if


!          -- UPDATE VELOCITY BY FOR_SHORT --

           atmvel(1:3,1:natom) =  atmvel(1:3,1:natom) &
                &               + 0.5d0*dt_short_cal*for_short(1:3,1:natom)

!          - local fix
           if (iflocalfix) then
              do ilfix = 1, nlfix
                 i = index_nlfix(ilfix)
                 atmvel(1:3,i) = 0.0d0
              end do
           end if

!          - local velocity
           if (iflocalvel) then
              do ilvel = 1, nlvel
                 i = index_nlvel(ilvel)
                 atmvel(1:3,i) = v_nlvel(1:3,ilvel)
              end do
           end if

!          - limit velocity
           if (iflimitmove) then
              call limitvel(dt_short_cal, &
                   &        limitdist)
           end if

!          --- RATTLE: store old coordinate at T to atmcor_old ---

           if (ifrattle .and. (nconst /= 0)) then
              atmcor_old(1:3,1:natom) = atmcor(1:3,1:natom)
           end if

!          -- UPDATE COORDINATE --

           atmcor(1:3,1:natom) =  atmcor(1:3,1:natom) &
                &               + dt_short_cal*atmvel(1:3,1:natom)

!          -- RATTLE_c: reset coordinate --

           if (ifrattle .and. (nconst /= 0)) then

              call rattle_c( atmcor_old, &
                   &         eps_rattle, dt_short_cal, dt_long_cal, &
                   &         istep_short, nstep_short )
!                  end if

           end if

!          -- RATTLE_c_z: fix z-coordinate --

           if (iflocalfixz) then

              call rattle_c_z(dt_short_cal,dt_long_cal,   &
                   &          istep_short,nstep_short)

           end if

!          -- RATTLE_c_zg: fix z-coordinate of COM --

           if (iflocalfixzg) then

              call rattle_c_zg(dt_short_cal,dt_long_cal, &
                   &           istep_short,nstep_short)

           end if

!          -- assign streaming velocity
           call ass_strmvel(npoly,npolytyp, &
                &           npoly_mole,npoly_atom, &
                &           nwater, &
                &           nmatom,nmatyp,nmatomtyp, &
                &           xcel,ycel,zcel, &
                &           iftcratom, &
                &           ifstrmvel)

!!! Time measure
#if defined(TIME_M)
           call get_timecount( dyn_time, dyn_time_elps, 'STOP' )
#endif
!!! Time measure

!          -- UPDATE FOR_SHORT & CONVERT TO ACCEL --

           call calforce(MTS_SHORT,for_short, &
                &        for_viric_short, pot_viric_short, &
                &        for_virilj_short, pot_virilj_short, &
                &        for_virimor_short, pot_virimor_short, &
                &        for_virish_short, pot_virish_short, &
                &        for_virirfh_short, pot_virirfh_short, &
                &        for_viridou_short, pot_viridou_short, &
                &        for_viricstmnb_short, pot_viricstmnb_short, &
                &        for_viricstmnbex_short, pot_viricstmnbex_short, &
                &        atm_viribo_short, atm_virian_short, &
                &        atm_virito_short, atm_viri14_short, &
                &        pot_virict_short, pot_viriljt_short, &
                &        pot_virimort_short, &
                &        pot_virisht_short, &
                &        pot_virirfht_short, &
                &        pot_viridout_short, &
                &        pot_viricstmnbt_short, &
                &        pot_viricstmnbext_short, &
                &        atm_viribot_short, atm_viriant_short, &
                &        atm_viritot_short, atm_viri14t_short, &
                &        mts_bond, mts_angl, mts_anglub, &
                &        mts_tors, mts_torsrb, mts_torsim, &
                &        mts_vdw,mts_ewr, mts_ewk, &
                &        mts_vdw14, mts_elc14, &
                &        mts_mor, mts_sh,mts_rfh, mts_dou, mts_cnpvw, &
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
                &        pot_tot,pot_nonbon, pot_vdw,pot_elc, pot_ewk, &
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
                &        nstep_short, istep_short)
                                ! short-range force
           call accforce( natom, atmmass, for_short )

!              --- UPDATE FOR_LONG & CONVERT TO ACCEL ---

           if (istep_short == nstep_short) then

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
                   &        nstep_short,istep_short)
                                ! long-range force
              call accforce( natom, atmmass, for_long )

           end if

!              -- UPDATE VELOCITY BY FOR_SHORT --

!!! Time measure
#if defined(TIME_M)
           call get_timecount_barrier( dyn_time, dyn_time_elps, 'START' )
#endif
!!! Time measure

           atmvel(1:3,1:natom) =  atmvel(1:3,1:natom) &
                &               + 0.5d0*dt_short_cal*for_short(1:3,1:natom)

!          --- UPDATE VELOCITY BY FOR_LONG ---
           if (istep_short == nstep_short) then

!             --- correct velocity with streaming velocity
              atmvel(1:3,1:natom) = atmvel(1:3,1:natom) &
                   &              - atmvel_strm(1:3,1:natom)

              atmvel(1:3,1:natom) =  atmvel(1:3,1:natom)*expcdt2 &
                   &               + for_long(1:3,1:natom)*bbv
                                ! higher factorized version
                                ! vi' = vi*exp(-dt/2*vxi(1))
                                !      +Fi/mi*bbv

!             --- reset velocity
              atmvel(1:3,1:natom) = atmvel(1:3,1:natom) &
                   &              + atmvel_strm(1:3,1:natom)

           end if

!          -- RATTLE_v: reset velocity --

           if (ifrattle .and. (nconst /= 0)) then

              call rattlep_v( eps_rattle, &
                   &          dt_long_cal,dt_short_cal, &
                   &          istep_short,nstep_short, &
                   &          cf, &
                   &          atm_viri_const,atm_virit_const)

           end if

!          - local fix
           if (iflocalfix) then
              do ilfix = 1, nlfix
                 i = index_nlfix(ilfix)
                 atmvel(1:3,i) = 0.0d0
              end do
           end if

!          - local velocity
           if (iflocalvel) then
              do ilvel = 1, nlvel
                 i = index_nlvel(ilvel)
                 atmvel(1:3,i) = v_nlvel(1:3,ilvel)
              end do
           end if

!          - limit velocity
           if (iflimitmove) then
              call limitvel(dt_short_cal, &
                   &        limitdist)
           end if

!          -- RATTLE_v_z: reset velocity --
           if (iflocalfixz) then

              call rattlep_v_z(dt_long_cal, dt_short_cal,   &
                   &           istep_short, nstep_short)

           end if

!          -- RATTLE_v_zg: reset velocity --
           if (iflocalfixzg) then

              call rattlep_v_zg(dt_long_cal,dt_short_cal, &
                   &            istep_short,nstep_short)

           end if

!!! Time measure
#if defined(TIME_M)
           call get_timecount( dyn_time, dyn_time_elps, 'STOP' )
#endif
!!! Time measure

        END DO              ! end of loop over NVT-XO-RESPA

!       --- UPDATE Nose-Hoover chain ---

!!! Time measure
#if defined(TIME_M)
        call get_timecount_barrier( dyn_time, dyn_time_elps, 'START' )
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

        call nhcint( degfree_all, &
             &       dt_long_cal, &
             &       mchain,text_c, &
             &       next,nyosh )

!       --- restoration velocity of VW
        if (ifcnp) then
           do j = 1,2
              i = nvwlist(j)
              atmvel(1:3,i) = atmvel_tmp(1:3,i)
           end do
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

!           --- CALCULATE PRESSURE OF SYSTEM ---

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

!           --- Substitute pressmol & pressatm to pint ---

        if (ifpmolcont) then
           pint = pressmol_tot
           pintt(1:3,1:3) = pressmolt_tot(1:3,1:3)
        else
           pint = pressatm_tot
           pintt(1:3,1:3) = pressatmt_tot(1:3,1:3)
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

!           --- Record the state of the end of this MD timestep ---

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
                   &     npoly,nwater,nmatom, &
                   &     xcel,ycel,zcel, &
                   &     xref,vref,timeref,pref, &
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
        call get_timecount( out_time, out_time_elps, 'STOP' )
#endif
!!! Time measure

!           --- CALCULATE KINETIC ENERGY & TOTAL ENERGY AT T + dT ---

!!! Time measure
#if defined(TIME_M)
        call get_timecount_barrier( out_time, out_time_elps, 'START' )
#endif
!!! Time measure

        call calkin( npoly, npolytyp, npoly_mole, &
             &       nwater, &
             &       nmatom,nmatyp, nmatomtyp, &
             &       degfree_poly, degfree_water, degfree_ma, &
             &       ene_kin_poly, ene_kin_water, ene_kin_ma, &
             &       temp_poly, temp_water, temp_ma, &
             &       degfree_all, ene_kin_all, temp_all, &
             &       mchain,text_c, &
             &       ene_kin_th, ene_pot_th, &
             &       pext_c, &
             &       ene_kin_ba, ene_pot_ba, extra_pot_ba )

!            ene_tot = ene_kin_poly + ene_kin_water + ene_kin_ma + pot_tot
        ene_tot = ene_kin_water + pot_tot
        do i = 1, npolytyp
           ene_tot = ene_tot + ene_kin_poly(i)
        end do
        do i = 1, nmatyp
           ene_tot = ene_tot + ene_kin_ma(i)
        end do
        ene_conserve =  ene_tot + ene_kin_th + ene_pot_th

!           --- Output energy, trajectory, velocity etc. ---
#if defined(MPI)
        if (irank == 0) then
#endif
           if ((mod(current_step,outinterval) == 1) .or. &
                & (outinterval == 1)) then

              if (ifoutene) then
                 call outene(ouene,current_step,eref,tempref, &
                      &      npolytyp,nmatyp, &
                      &      pot_tot,pot_nonbon,pot_vdw, &
                      &      pot_elc,pot_ewk,pot_ewc,pot_ewnq, &
                      &      pot_vdw14,pot_elc14, &
                      &      pot_bond,pot_angl,pot_anglub, &
                      &      pot_tors,pot_torsrb,pot_torsim, &
                      &      pot_mor, &
                      &      pot_sh, &
                      &      pot_rfh, &
                      &      pot_dou, &
                      &      pot_rpvw,pot_vw, &
                      &      pot_cstmnb, &
                      &      ncstmnbex,pot_cstmnbex, &
                      &      pot_posres, &
                      &      pot_pbias, &
                      &      ene_kin_poly,ene_kin_water, &
                      &      ene_kin_ma,ene_kin_all, &
                      &      ene_tot, &
                      &      temp_poly,temp_water,temp_ma,temp_all, &
                      &      ene_kin_th,ene_pot_th, &
                      &      ene_kin_ba,ene_pot_ba, &
                      &      ene_conserve)
              end if

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

              if (ifoutthe) then
                 call outthe(outhe, current_step, timeref, &
                      &      mchain, xlogs, vlogs)
              end if

              ! if (ifoutbar) then
              !     call outbar(oubar,current_step,xref,timeref, &
              !     &           xcel,ycel,zcel,xlogv,vlogv,xboxh,vboxg, &
              !     &           pcont_axis)
              ! end if

           end if

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
        call get_timecount( out_time, out_time_elps, 'STOP' )
#endif
!!! Time measure

     END IF

!!! Time measure
#if defined(TIME_M)
     call get_timecount_barrier( out_time, out_time_elps, 'START' )
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
                &       resname_matom_pdb)
        end if

        if (ifpotbias .and. (current_step == 0)) then
           call outumb( ouumb, &
                &       xref, eref)
        endif

#if defined(MPI)
     end if
#endif

!!! Time measure
#if defined(TIME_M)
     call get_timecount( out_time, out_time_elps, 'STOP' )
#endif
!!! Time measure

!
!---- end of NVT-MTS (XO-RESPA)

!!! Time measure
#if defined(TIME_M) || defined(TIME_MALL)
     call get_timecount_barrier( total_time, total_time_elps, 'STOP' )
#endif
!!! Time measure

!!! Time measure output
#if defined(TIME_M) || defined(TIME_MALL)
     if (current_step /= 0) then
        call out_timecount()
     end if
#endif
!!! Time measure output

  END DO


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

end subroutine moldyn_nhc
