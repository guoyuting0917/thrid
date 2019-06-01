!************************************************************
!*  peachgk_md                        Ver.2.142             *
!*    Molecular Dynamics Program                            *
!*      providing Enhanced Approach to CHallenging problems *
!*      in General Kinds of microscale research field       *
!*               (simplified and augmented PEACH by G. K.)  *
!*                             programmed by Gota Kikugawa  *
!*          Institute of Fluid Science (IFS), Tohoku Univ.  *
!************************************************************
!********     MD run script Ver.5.8     ********
! Time-stamp: <>

program peachgk_md

! include interface modules
  use interface_preproc
  use interface_mddriver
  use interface_mdtech

! include global variables
  use md_global     ! include variables for md calculation
  use spme_global   ! include variables for spme method
!C-MPI#if defined(MPI)
  use mpi_global    ! include mpi variables
!C-MPI#endif

  implicit none

! All variables and arrays:
! MD control parameters
  character(80):: iucorname(maxnpolytyp) ! input file name
  character(80):: iutopname(maxnpolytyp) !       "
  character(80):: iuwtopname             !       "
  character(80):: iuparavdwname          !       "
  character(80):: iuparabondname         !       "
  character(80):: iuparaconstname        !       "
  character(80):: iuparacstmnbname       !       "
  character(80):: iuaddtopname           !       "
  character(80):: iustrmvelname          !       "
  character(80):: iuposresname           !       "

  character(80):: iostarecname    ! state file name

  character(80):: ousumname ! output file name
  character(80):: ouenename ! output file name
  character(80):: ouposname ! output file name
  character(80):: ouvelname ! output file name
  character(80):: ouforname ! output file name
  character(80):: outhename ! output file name
  character(80):: oubarname ! output file name
  character(80):: ouprename ! output file name
  character(80):: outhcname ! output file name
  character(80):: oupdbname ! output file name
  character(80):: ouhtfname ! output file name
  character(80):: oumtfname ! output file name
  character(80):: ouumbname ! output file name

  logical:: ifoutene = .true.  ! if ouput energy file
  logical:: ifoutpos = .true.  ! if ouput position file
  logical:: ifoutvel = .true.  ! if ouput velocity file
  logical:: ifoutfor = .false. ! if ouput force file
  logical:: ifoutthe = .true.  ! if ouput NVT file
  logical:: ifoutbar = .true.  ! if ouput NPT file
  logical:: ifoutpre = .true.  ! if ouput pressure file

  integer:: npoly           ! all number of poly
  integer:: npolytyp        ! number of poly type
  integer:: npoly_mole(maxnpolytyp) ! number of molecules of each poly
  integer:: npoly_atom(maxnpolytyp) ! number of atoms belonging to poly

  integer:: nwater          ! number of H2O molecules

  integer:: nmatom          ! number of monatomic molecules
  integer:: nmatyp          ! number of species of monatomic mole.
  integer:: nmatomtyp(maxnmatyp) ! each number of monatomic mole.

  character(80):: polytyp_free(maxnpolytyp,maxnword)
                                ! use for poly type control
  character(80):: watertyp_free(maxnword) ! use for water type control
  character(80):: matomtyp_free(maxnmatyp,maxnword)
                                ! use for matom type control

  integer:: xmaxpo(maxnpolytyp) ! use for positioning of polymer1
  integer:: ymaxpo(maxnpolytyp) ! use for positioning of polymer1
  integer:: zmaxpo(maxnpolytyp) ! use for positioning of polymer1
  integer:: xmaxw = 0       ! use for positioning of water
  integer:: ymaxw = 0       ! use for positioning of water
  integer:: zmaxw = 0       ! use for positioning of water
  integer:: xmaxma(maxnmatyp) ! use for positioning of monatomic mole.
  integer:: ymaxma(maxnmatyp) ! use for positioning of monatomic mole.
  integer:: zmaxma(maxnmatyp) ! use for positioning of monatomic mole.
  real(8):: inicorpo(maxnpolytyp,maxnword) ! use for positioning of poly
  real(8):: inicorw(maxnword) ! use for positioning of water
  real(8):: inicorma(maxnmatyp,maxnword) ! use for positioning of ma
  integer:: ncrecorpo       ! max number of poly type for createcor
  integer:: index_crecorpo(maxnpolytyp) ! index of polymer for createcor
  integer:: ncrecorw        ! water for createcor
  integer:: index_crecorw   ! index of water type for createcor
  integer:: ncrecorma       ! max number of matom type for createcor
  integer:: index_crecorma(maxnmatyp) ! index of matom type for createcor

  logical:: ifsetcor(maxnatom) ! frag for checking if coordinate has set

!      logical:: ifnve = .false. ! NVE ensemble MD
!      logical:: ifnhc = .false. ! Nose-Hoover chain (NVT) MD
!      logical:: ifmtk = .true.  ! Martyna-Tobias-Klein (NPT) MD

  integer:: maxnstep        ! maximum step of MD (0-maxnstep)
  integer:: nstage          ! stage number of MD
  integer:: nstep_stage(maxnstage) ! time step of each stage

  integer:: mdcont_stage(maxnstage) ! MD control parameter of each stage

  integer:: nstep_maxwell   ! time step of maxwell distribution
  integer:: nstep_expand    ! time step of cell expansion

  real(8):: xcel             ! initial x0 cell length[m]
  real(8):: ycel             ! initial y0 cell length[m]
  real(8):: zcel             ! initial z0 cell length[m]

  real(8):: yratio           ! y cell ratio of y to x
  real(8):: zratio           ! z cell ratio of z to x

  real(8):: r_expand         ! expansion ratio of cell (volume)

  logical:: ifstarec = .false. ! read old state
  logical:: ifcreatecor = .false. ! create new coordinate

  logical:: ifrdaddtop = .false. ! input additional topology information

  logical:: ifcenterfix_all    ! center fix for all
  logical:: ifcenterfix_poly   ! center fix for polymer
  logical:: ifcenterfix_water  ! center fix for water
  logical:: ifcenterfix_ma     ! center fix for monatomic mole.
  character(4):: cenfix_free   ! COM not fixed in this direction
  logical:: ifcenterfix_polytyp(maxnpolytyp)  ! center fix for each polymer
  logical:: ifcenterfix_watertyp              ! center fix for each water
  logical:: ifcenterfix_matyp(maxnmatyp)  ! center fix for each monatomic mole.

  real(8):: dt_long_cal      ! time step of long force [sec]

!!! From Ver.1.74, nstep_med is forced to 1. !!!
  integer:: nstep_med = 1   ! number of step for medium force
  integer:: nstep_short     ! number of step for short force

!!! From Ver.1.74, do not use 2(mts_med). !!!
  integer:: mts_bond = 0    ! MTS flag for bond
  integer:: mts_angl = 0    ! MTS flag for angle
  integer:: mts_anglub = 0  ! MTS flag for Urey-Bradley angle
  integer:: mts_tors = 0    ! MTS flag for torsion
  integer:: mts_torsrb = 0  ! MTS flag for torsionrb
  integer:: mts_torsim = 0  ! MTS flag for torsionim
  integer:: mts_vdw = 0     ! MTS flag for vdw interaction
  integer:: mts_ewr = 0     ! MTS flag for ewald real(=vdw)
  integer:: mts_ewk = 0     ! MTS flag for ewald wave
  integer:: mts_vdw14 = 0   ! MTS flag for 14vdw
  integer:: mts_elc14 = 0   ! MTS flag for 14elc(=mts_vdw14)
  integer:: mts_mor = 0     ! MTS flag for Morse interaction
  integer:: mts_sh = 0      ! MTS flag for SH interaction
  integer:: mts_rfh = 0     ! MTS flag for RFH interaction
  integer:: mts_dou = 0     ! MTS flag for DOU interaction
  integer:: mts_cnpvw = 0   ! MTS flag for CNP_VW
  integer:: mts_cstmnb = 0  ! MTS flag for custom NB interaction
  integer:: mts_posres = 0  ! MTS flag for position restraint
  integer:: mts_potbias = 0 ! MTS flag for bias potential

!---- virtual time step for evaluating constraint force
  integer:: nstep_vir       ! number of step for virtual time integration

!---- limit atomic motion to a specific distance
  logical:: iflimitmove     ! if limit atomic motion to a specific distance
  real(8):: limitdist       ! maximum atomic displacement when doing
                            ! time integration [m] (structure relaxation)

!---- for Nose-Hoover chain
  integer:: mchain          ! Nose-Hoover chain number (>1 in NVT)
  real(8):: tfreq            ! temperature frequency [1/s]
  real(8):: text             ! external temp. [K] (Nose-Hoover chain)

!---- for Andersen (Hoover type) barostat
  real(8):: vfreq           ! volume change frequency [1/s]
  real(8):: pext            ! external pressure [Pa]
  logical:: ifpatmcont      ! atomic pressure control
  logical:: ifpmolcont      ! molecular pressure control
  character(5):: pcont_axis ! axis for pressure control (iso, aniso, etc.)

!---- for higher order Trotter expansion
  integer:: next            ! iteration number of extended system
  integer:: nyosh           ! expansion order of Yoshida-Suzuki method
                                ! choose 1, 3 or 5

!---- some MD flags
  logical:: ifrattle        ! rattle flag

  logical:: ifewald         ! ewald flag
  logical:: ifspme          ! SPME (Smooth Particle Mesh Ewald) flag
  logical:: iffennell       ! Fennell flag

  logical:: ifljari         ! arithmetic mean for LJ cross parameter
  logical:: ifljgeo         ! geometric mean for LJ cross parameter

  logical:: iflocalheat     ! local heating flag
  logical:: ifregionheat    ! region temp control flag
  logical:: ifregionhf      ! region heat flux control flag
  logical:: ifreglange      ! region Langevin thermostat flag
  logical:: iftcratom       ! region temp. or h.f. control
                            ! based on atom or mole.
  logical:: ifoutthc        ! flag for outputting thermal control file
  logical:: iflocalfix      ! fix atoms flag
  logical:: iflocalfixz     ! flag for fixing z coordinate of atoms
  logical:: iflocalfixzg    ! flag for fixing z coordinate of COM of molecules
  logical:: ifcnp           ! flag to control normal pressure
  logical:: ifposres        ! position restraint flag
  logical:: ifpotbias       ! bias potential flag
  logical:: iflocalvel      ! flag to force local atomic velocity
  logical:: ifstrmvel       ! flag to input and use streaming velocity

!!! if you use NPT dynamics, you must choose ifcalpremole or ifcalpreatom !!!
  logical:: ifcalpremole    ! pressure calculation of molecule
  logical:: ifcalpreatom    ! pressure calculation of atom
  logical:: ifnetqcorrp     ! net charge correction for pressure

!---- pressure calculation of L-J long-range correction
  logical:: ifcalljlong     ! long-range correction in pressure
  character(2):: solvetyp ! solvent molecule
  integer:: nsolve          ! number of solvent molecules

!---- parameter for ewald method
  real(8):: alpha            ! parameter alpha [1/m]
!      take care maxwave
  integer:: kmax            ! parameter kmax
  real(8):: rrcut            ! ewald real space cutoff length [m]

!---- parameter for SPME method
  integer:: nfft1, nfft2, nfft3 ! grid points in SPME
  integer:: pme_order       ! B-spline order
                            ! E.g. cubic is order 4
                            ! fifth degree is order 6 etc.
                            ! The order must be an even number
                            !                   and at least 4.

!---- variable for local heating
  integer:: nlheat_poly     ! number of poly type for local heating
  integer:: index_nlheat_poly(maxnpolytyp)
                            ! index of poly type for local heating
  real(8):: tcont_nlheat_poly(maxnpolytyp)
                            ! control temp. of poly type for local heating
  integer:: nlheat_water    ! number of water for local heating
  real(8):: tcont_nlheat_water
                            ! control temp. of water for local heating
  integer:: nlheat_ma       ! number of matom type for local heating
  integer:: index_nlheat_ma(maxnmatyp)
                            ! index of matom type for local heating
  real(8):: tcont_nlheat_ma(maxnmatyp)
                            ! control temp. of matom type for local heating

!---- variable for fix atoms
  integer:: nlfix           ! number of fix atoms
  integer:: index_nlfix(maxnatom) ! index of fix atoms

!---- variable for local velocity
  integer:: nlvel           ! number of atoms for velocity fix
  integer:: index_nlvel(maxnatom) ! index of vel-fix atoms
  real(8):: v_nlvel(1:3,maxnatom) ! local velocity values

!---- other MD parameters
  real(8):: eps0 = 8.854187817d-12 ! dielectric constant of vacuum [C^2/Jm]
  real(8):: div_factor_14vdw = 2.0d0 ! division factor of 14vdw
  real(8):: div_factor_14elc = 1.2d0 ! division factor of 14elc

  real(8):: rcut             ! vdw cutoff length [m]

  logical:: ifcellindex     ! flag for cell index

  logical:: ifbook          ! flag for bookkeeping
  real(8):: rcut_book        ! cut off radius of bookkeeping[m]
  integer:: nstep_book      ! bookkeeping interval

  real(8):: tcont_poly(maxnpolytyp) ! poly Temp. [K] in NVT (Woodcock)
  real(8):: tcont_water      ! H2O Temp. [K] in NVT (Woodcock)
  real(8):: tcont_ma(maxnmatyp) ! monatomic mole. Temp. [K]
                                   ! in NVT (Woodcock)
  real(8):: tcont_poly_ini   ! poly Temp. [K] in NVT (Woodcock) for md_h
  real(8):: tcont_water_ini  ! H2O Temp. [K] in NVT (Woodcock) for md_h
  real(8):: tcont_ma_ini     ! MA Temp. [K] in NVT (Woodcock) for md_h

  integer:: tcontinterval   ! interval of temp. control
  integer:: outinterval     ! interval of outputting data
  integer:: pressinterval   ! interval of pressure output
  integer:: heatfinterval   ! interval of heatf output
  integer:: recinterval     ! state record interval

  character(2):: oatmtyp    ! O atomtype of water model
  character(2):: hatmtyp    ! H atomtype of water model
  character(2):: monoatmtyp(maxnmatyp) ! monatomic mole. type

  integer:: randseed        ! random seed for createcor etc.

  real(8):: compfact         ! compact factor using at poly arrange(<1.0)

  real(8):: eps_rattle       ! tolerance (relative difference)
                             ! for bond length constraint by RATTLE

  real(8):: rcutmor          ! Morse cutoff length [m]

  logical:: ifcellindex_mor ! flag for cell index (morse)

  logical:: ifbookmor       ! flag for bookkeeping of Morse interaction
  real(8):: rcut_bookmor     ! cut off radius of bookkeeping[m] of Morse
  integer:: nstep_bookmor   ! bookkeeping interval of Morse interaction

  real(8):: rcutsh           ! SH cutoff length [m]

  logical:: ifcellindex_sh  ! flag for cell index (SH)

  logical:: ifbooksh        ! flag for bookkeeping of SH interaction
  real(8):: rcut_booksh      ! cut off radius of bookkeeping[m] of SH
  integer:: nstep_booksh    ! bookkeeping interval of SH interaction

  real(8):: rcutrfhfo        ! RFH(FO) cutoff length [m]

  logical:: ifcellindex_rfhfo ! flag for cell index (RFH(FO))

  logical:: ifbookrfhfo     ! flag for bookkeeping of RFH(FO) interaction
  real(8):: rcut_bookrfhfo   ! cut off radius of bookkeeping[m] of RFH(FO)
  integer:: nstep_bookrfhfo ! bookkeeping interval of RFH(FO) interaction

  real(8):: rcutrfhoo        ! RFH(OO) cutoff length [m]

  logical:: ifcellindex_rfhoo ! flag for cell index (RFH(OO))

  logical:: ifbookrfhoo     ! flag for bookkeeping of RFH(OO) interaction
  real(8):: rcut_bookrfhoo   ! cut off radius of bookkeeping[m] of RFH(OO)
  integer:: nstep_bookrfhoo ! bookkeeping interval of RFH(OO) interaction

  real(8):: rcutrfhoh        ! RFH(OH) cutoff length [m]

  logical:: ifcellindex_rfhoh ! flag for cell index (RFH(OH))

  logical:: ifbookrfhoh     ! flag for bookkeeping of RFH(OH) interaction
  real(8):: rcut_bookrfhoh   ! cut off radius of bookkeeping[m] of RFH(OH)
  integer:: nstep_bookrfhoh ! bookkeeping interval of RFH(OH) interaction

  real(8):: rcutdouo         ! DOU cutoff length [m] for O-Au
  real(8):: rcutindouo       ! DOU cutin length [m] for O-Au

  logical:: ifcellindex_douo ! flag for cell index (DOU) for O-Au

  logical:: ifbookdouo      ! flag for bookkeeping of DOU interaction (O-Au)
  real(8):: rcut_bookdouo    ! cut off radius of bookkeeping[m] of DOU (O-Au)
  integer:: nstep_bookdouo  ! bookkeeping interval of DOU interaction (O-Au)

  real(8):: rcutdouh         ! DOU cutoff length [m] for H-Au

  logical:: ifcellindex_douh ! flag for cell index (DOU) for H-Au

  logical:: ifbookdouh      ! flag for bookkeeping of DOU interaction (H-Au)
  real(8):: rcut_bookdouh    ! cut off radius of bookkeeping[m] of DOU (H-Au)
  integer:: nstep_bookdouh  ! bookkeeping interval of DOU interaction (H-Au)

  real(8):: rcutrpvw ! RP-VW cutoff length

  logical:: ifbookrpvw       ! flag for bookkeeping of RP-VW interaction
  real(8):: rcut_bookrpvw    ! cut off radius of bookkeep of RP-VW interaction
  integer:: nstep_bookrpvw   ! bookkeeping interval of RP-VW interaction

!---- parameters for custom interaction ----!
  !!! Other parameters can be found in each custom interaction module program.

  logical:: ifcstmnb          ! flag if using custom NB interaction
  logical:: ifcellindex_cstmnb ! flag for cell index (custom NB)
  logical:: ifbookcstmnb       ! flag for bookkeeping of custom NB interaction

!---- set charge from cor file
  logical:: ifsetchrg(maxnpolytyp) ! if set charge from cor file

!---- PDB output
  logical:: ifoutpdb        ! flag for outputting PDB format file
  integer:: nstep_pdbout    ! MD step for output of PDB file

!---- parameters for energy minimization
  real(8):: d_rini            ! initial displacement dr for EM [m]
  real(8):: d_rmax            ! maximum displacement dr for EM [m]

  real(8):: d_econv           ! convergence condition for energy in EM [J]
  real(8):: d_rmsf            ! convergence condition for
                              !    root mean square force in EM [N]

!---- base value for non-dimensionalize
  real(8):: xref = 1.0d-10   ! distanse base value [m] (adjust to Angstrom)
!      real*8:: eref = 7.60078d-22 ! energy base value [J]
  real(8):: eref = 6.9511042d-21 ! energy base value [J]
                                 ! = 4.18605*1.0e+3/6.0221367e+23
                                 !   (adjust to kcal/mol)
  real(8):: mref = 1.99431d-26 ! mass base value [kg] (C atom)
  real(8):: qref = 1.602117733d-19 ! charge base value [C]

  real(8):: vref             ! velocity base value [m/s]
  real(8):: timeref          ! time base value [sec]
  real(8):: tempref          ! temperature base value [K]
  real(8):: pref             ! pressure base value [Pa]
  real(8):: fref             ! force base value [N]

  real(8):: eps0ref          ! dielectric constant base value [c^2/Jm]

!---- i/o unit
  integer:: iucor(maxnpolytyp) ! input poly coordinate file unit
  integer:: iutop(maxnpolytyp) ! input poly topology file unit
  integer:: iuwtop          ! input water topology file unit
  integer:: iuparavdw       ! input vdw parameter file unit
  integer:: iuparabond      ! input bond parameter file unit
  integer:: iuparaconst     ! input const parameter file unit
  integer:: iuparacstmnb    ! input custom NB parameter file unit
  integer:: iuaddtop        ! input additional topology file unit
  integer:: iustrmvel       ! input streaming velocity file unit
  integer:: iuposres        ! input position restraint ref. file unit

  integer:: iostarec        ! state record file unit

  integer:: ousum           ! output parameter summarization file unit
  integer:: ouene           ! output unit for output energy data
  integer:: oupos           ! output unit for output position data
  integer:: ouvel           ! output unit for output velocity data
  integer:: oufor           ! output unit for output force data
  integer:: outhe           ! output unit for output thermostat data
  integer:: oubar           ! output unit for output barostat data
  integer:: oupre           ! output unit for outpre velocity data
  integer:: outhc           ! output unit for outthc thermal control data
  integer:: oupdb           ! output unit for outpdb PDB data
  integer:: ouhtf           ! output unit for outhtf heat flux data
  integer:: oumtf           ! output unit for outmtf momentum flux data
  integer:: ouumb           ! output unit for outumb bias potential data

!---- degree of freedom
  integer:: degfree_poly(maxnpolytyp) ! degree of freedom of polymer1
  integer:: degfree_water   ! degree of freedom of H2O
  integer:: degfree_ma(maxnmatyp) ! degree of freedom of monatomic mole.
  integer:: degfree_all     ! degree of freedom of all atoms

!---- degree of freedom to fixed
  integer:: nlfix_deg_poly(maxnpolytyp) ! fixed degree of freedom of poly
  integer:: nlfix_deg_water ! fixed degree of freedom of H2O
  integer:: nlfix_deg_ma(maxnmatyp) ! fixed degree of freedom of matom

!---- degree of freedom to fixed (z-coordinate)
  integer:: nlfixz_deg_poly(maxnpolytyp) ! fixed degree of freedom of poly
  integer:: nlfixz_deg_water ! fixed degree of freedom of H2O
  integer:: nlfixz_deg_ma(maxnmatyp) ! fixed degree of freedom of matom

!---- degree of freedom to fixed (z-coordinate) for COM based fixation
  integer:: nlfixzg_deg_poly(maxnpolytyp) ! fixed degree of freedom of poly
  integer:: nlfixzg_deg_water ! fixed degree of freedom of H2O
  integer:: nlfixzg_deg_ma(maxnmatyp) ! fixed degree of freedom of matom

!---- degree of freedom to fixed velocity
  integer:: nlvel_deg_poly(maxnpolytyp) ! fixed degree of freedom of poly
  integer:: nlvel_deg_water ! fixed degree of freedom of H2O
  integer:: nlvel_deg_ma(maxnmatyp) ! fixed degree of freedom of matom

!---- time step
  real(8):: dt_med_cal       ! time step of medium force
  real(8):: dt_short_cal     ! time step of short force

!---- force of each atoms
  real(8):: for_long(3,maxnatom)  ! long-range force
!      real(8):: for_med(3,maxnatom)   ! medium-range force
  real(8):: for_short(3,maxnatom) ! short-range force

!---- force part of molecular pressure
  real(8):: for_viric_long(3,maxnatom)  ! long-range virial (coulomb force)
  real(8):: for_viric_med(3,maxnatom)   ! medium-range virial (coulomb force)
  real(8):: for_viric_short(3,maxnatom) ! short-range virial (coulomb force)

  real(8):: for_virilj_long(3,maxnatom)  ! long-range virial (L-J force)
  real(8):: for_virilj_med(3,maxnatom)   ! medium-range virial (L-J force)
  real(8):: for_virilj_short(3,maxnatom) ! short-range virial (L-J force)

  real(8):: for_virimor_long(3,maxnatom)  ! long-range virial (morse force)
  real(8):: for_virimor_med(3,maxnatom)   ! medium-range virial (morse force)
  real(8):: for_virimor_short(3,maxnatom) ! short-range virial (morse force)

  real(8):: for_virish_long(3,maxnatom)  ! long-range virial (SH force)
  real(8):: for_virish_med(3,maxnatom)   ! medium-range virial (SH force)
  real(8):: for_virish_short(3,maxnatom) ! short-range virial (SH force)

  real(8):: for_virirfh_long(3,maxnatom)  ! long-range virial (RFH force)
  real(8):: for_virirfh_med(3,maxnatom)   ! medium-range virial (RFH force)
  real(8):: for_virirfh_short(3,maxnatom) ! short-range virial (RFH force)

  real(8):: for_viridou_long(3,maxnatom)  ! long-range virial (DOU force)
  real(8):: for_viridou_med(3,maxnatom)   ! medium-range virial (DOU force)
  real(8):: for_viridou_short(3,maxnatom) ! short-range virial (DOU force)

  real(8):: for_viricstmnb_long(3,maxnatom)  ! long-range virial (cstmNB force)
  real(8):: for_viricstmnb_med(3,maxnatom)   ! med-range virial (cstmNB force)
  real(8):: for_viricstmnb_short(3,maxnatom) ! short-range virial (cstmNB force)

!---- potential variables
  real(8):: pot_ewc = 0.0d0  ! potential of ewald self-energy

!---- variables for Andersen (Hoover type) barostat
  real(8):: pint = 0.0d0   ! internal pressure

  real(8):: pintt(3,3) = 0.0d0   ! internal pressure tensor

!---- variables for long-range correction
!      real*8:: vdw_welij_solve  ! well depth of vdw parameter of solvent
!      real*8:: vdw_radij_solve  ! radius of vdw parameter of solvent
  integer:: solveindex      ! atmindex of solvent atom

!!! For old MT random generator
!!---- variables for random number generator
!  integer,parameter:: nperiod = 624   ! period parameter
!  integer:: mt(0:nperiod-1)
!  integer:: mti

!---- timestep control
  integer:: inistep         ! initial step of MD
  integer:: endstep         ! end step of MD
  integer:: md_cont         ! MD control flag

!---- variables for net charge calculation
  real(8):: netchrgsq        ! = (sum(qi))**2

!---- variables for region-based temp. control
  integer:: ntcregion       ! number of region to control temp.
  real(8):: tcxpos1(maxntcregion),tcxpos2(maxntcregion)
                                ! x-position of temp. control region
  real(8):: tcypos1(maxntcregion),tcypos2(maxntcregion)
                                ! y-position of temp. control region
  real(8):: tczpos1(maxntcregion),tczpos2(maxntcregion)
                                ! z-position of temp. control region
  real(8):: r_tcont(maxntcregion) ! control temp. in each region

!---- variables for region-based heat flux control
  integer:: nhfcregion      ! number of region to control heat flux
  real(8):: hfczpos1(maxntcregion),hfczpos2(maxntcregion)
                                ! z-position of heat flux control region
  real(8):: r_hfcont(maxntcregion) ! magnitude of heat flux in each region
                                ! (converted to the input energy)

!---- variables for region-based Langevin thermostat
  integer:: nlangeregion            ! number of region for Langevin thermo.
  real(8):: ltxpos1(maxntcregion),ltxpos2(maxntcregion)
                                    ! x-position of temp. control region
  real(8):: ltypos1(maxntcregion),ltypos2(maxntcregion)
                                    ! y-position of temp. control region
  real(8):: ltzpos1(maxntcregion),ltzpos2(maxntcregion)
                                    ! z-position of temp. control region
  real(8):: r_ltemp(maxntcregion)   ! control temp. in each region
  real(8):: r_ltdamp(maxntcregion)  ! damping factor in each region [1/s]

!---- variables for calculation of heat flux
  logical:: ifhfvol         ! local volume-based or local surface-based

  integer:: nhfregion       ! number of region to calculate heat flux
  real(8):: hfzpos1(maxnhfregion),hfzpos2(maxnhfregion)
                                ! z-position of region for heat flux

  integer:: hftyp_pmole(maxnpolytyp)
                                ! atom- or mole-based heat flux cal. for poly
  integer:: hftyp_water     ! atom- or mole-based heat flux cal. for water

  integer:: hftyp_atm(maxnatom) ! atom- or mole-based heat flux cal.
                                    !   for each atom

!---- flag to calculate & output momentum or heat flux
  logical:: ifcalmtf
  logical:: ifcalhtf

! variables for outputting momentum flux
  character(3):: mtfoutdir ! direction to output data of momentum

!---- variables for setting charge value from cor file
  real(8):: atmchrg_tmp(maxnatom) ! temporary atom charge read from cor file

!---- variables for setting the PDB residue name for each mole type
  character(4):: resname_poly_pdb(maxnpolytyp) ! residue name for poly
  character(4):: resname_water_pdb ! residue name for water
  character(4):: resname_matom_pdb(maxnmatyp) ! residue name for matom

!---- variables for setting mass & charge value from add_top file
  logical:: ifsetatmmass(maxnatom) ! flag for setting mass by add_top
  logical:: ifsetatmchrg(maxnatom) ! flag for setting charge by add_top

!---- local variables
  integer:: istage

!     +     +     +     +     +     +     +

!---- mpi initialization ----

  irank = 0
  nproc = 1

#if defined(MPI)
  call mpi_init( ierror )
  call mpi_comm_size( MPI_COMM_WORLD, nproc, ierror )
  call mpi_comm_rank( MPI_COMM_WORLD, irank, ierror )
  write(6,*) ' process',irank,' starts.'
#endif

!---- program notation ----

  if (irank == 0) then

     write(6,*) '- PEACHGK_MD starts. Ver.2.142 -'
     write(6,*) '-------------------------------------------------------'
     write(6,*) '  Molecular Dynamics Program'
     write(6,*) '    providing Enhanced Approach to CHallenging problems'
     write(6,*) '    in General Kinds of research field'
     write(6,*) '-------------------------------------------------------'
     write(6,*)

  end if

!---- MPI initialization ----

  call MPI_loopinit(0)

!---- Dynamic memory allocation for large arrays ----

  call malloc_larray()

!---- read MD script file ----
  call rdscript(pgkiniver,pgkinidate, &
       &        md_0k,md_h,md_t,md_mtk,md_nhc,md_nve,md_htf,md_ems, &
       &        mts_long,mts_med,mts_short, &
       &        maxnword, &
       &        iucorname,iutopname,iuwtopname, &
       &        iuparavdwname,iuparabondname,iuparaconstname, &
       &        iuparacstmnbname, &
       &        iuaddtopname, &
       &        iustrmvelname, &
       &        iuposresname, &
       &        iostarecname, &
       &        ousumname,ouenename,ouposname,ouvelname,ouforname, &
       &        outhename,oubarname,ouprename,outhcname, &
       &        oupdbname, &
       &        ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
       &        ifoutbar,ifoutpre, &
       &        maxnpolytyp, &
       &        npoly,npolytyp,npoly_mole,npoly_atom, &
       &        nwater, &
       &        maxnmatyp,nmatom,nmatyp,nmatomtyp, &
       &        polytyp_free,watertyp_free,matomtyp_free, &
       &        xmaxpo,ymaxpo,zmaxpo,xmaxw,ymaxw,zmaxw, &
       &        xmaxma,ymaxma,zmaxma, &
       &        inicorpo,inicorw,inicorma, &
       &        ncrecorpo,index_crecorpo, &
       &        ncrecorw,index_crecorw, &
       &        ncrecorma,index_crecorma, &
       &        maxnstage,maxnstep, &
       &        nstage,nstep_stage,mdcont_stage, &
       &        nstep_maxwell,nstep_expand, &
       &        xcel,ycel,zcel,yratio,zratio, &
       &        r_expand, &
       &        ifstarec,ifcreatecor, &
       &        ifrdaddtop, &
       &        ifcenterfix_all, &
       &        ifcenterfix_poly,ifcenterfix_water, &
       &        ifcenterfix_ma, &
       &        cenfix_free, &
       &        ifcenterfix_polytyp, &
       &        ifcenterfix_watertyp, &
       &        ifcenterfix_matyp, &
       &        dt_long_cal, &
       &        nstep_med,nstep_short, &
       &        mts_bond,mts_angl,mts_anglub, &
       &        mts_tors,mts_torsrb,mts_torsim, &
       &        mts_vdw,mts_ewr,mts_ewk, &
       &        mts_vdw14,mts_elc14, &
       &        mts_mor,mts_sh,mts_rfh,mts_dou,mts_cnpvw, &
       &        mts_cstmnb, &
       &        mts_posres, &
       &        mts_potbias, &
       &        nstep_vir, &
       &        iflimitmove,limitdist, &
       &        mchain,tfreq,text, &
       &        vfreq,pext, &
       &        ifpatmcont,ifpmolcont, &
       &        pcont_axis, &
       &        ifnetqcorrp, &
       &        next,nyosh, &
       &        ifrattle, &
       &        ifewald,ifspme,iffennell, &
       &        ifljari,ifljgeo, &
       &        iflocalheat, &
       &        ifregionheat,ifregionhf,ifreglange,iftcratom,ifoutthc, &
       &        iflocalfix,iflocalfixz, &
       &        iflocalfixzg, &
       &        ifcnp, &
       &        ifposres, &
       &        ifpotbias, &
       &        iflocalvel,ifstrmvel, &
       &        ifcalpremole,ifcalpreatom, &
       &        ifcalljlong,solvetyp,nsolve, &
       &        alpha,kmax,rrcut, &
       &        nfft1,nfft2,nfft3,pme_order, &
       &        tcont_poly,tcont_water,tcont_ma, &
       &        tcont_poly_ini,tcont_water_ini,tcont_ma_ini, &
       &        tcontinterval,outinterval, &
       &        pressinterval,heatfinterval, &
       &        recinterval, &
       &        oatmtyp,hatmtyp,monoatmtyp, &
       &        randseed, &
       &        compfact, &
       &        eps_rattle, &
       &        rcut,ifbook,rcut_book,nstep_book, &
       &        ifcellindex, &
       &        rcutmor,ifbookmor,rcut_bookmor,nstep_bookmor, &
       &        ifcellindex_mor, &
       &        rcutsh,ifbooksh,rcut_booksh,nstep_booksh, &
       &        ifcellindex_sh, &
       &        rcutrfhfo,ifbookrfhfo,rcut_bookrfhfo, &
       &        nstep_bookrfhfo, &
       &        ifcellindex_rfhfo, &
       &        rcutrfhoo,ifbookrfhoo,rcut_bookrfhoo, &
       &        nstep_bookrfhoo, &
       &        ifcellindex_rfhoo, &
       &        rcutrfhoh,ifbookrfhoh,rcut_bookrfhoh, &
       &        nstep_bookrfhoh, &
       &        ifcellindex_rfhoh, &
       &        rcutdouo,rcutindouo,ifbookdouo,rcut_bookdouo, &
       &        nstep_bookdouo, &
       &        ifcellindex_douo, &
       &        rcutdouh,ifbookdouh,rcut_bookdouh, &
       &        nstep_bookdouh, &
       &        ifcellindex_douh, &
       &        rcutrpvw, ifbookrpvw, rcut_bookrpvw, &
       &        nstep_bookrpvw, &
       &        ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
       &        ifsetchrg, &
       &        ifoutpdb,nstep_pdbout, &
       &        d_rini,d_rmax,d_econv,d_rmsf, &
       &        nspltbl)

#if defined(HF)
  call rdheatf(ouhtfname, &
       &       oumtfname, &
       &       xref, &
       &       zcel, &
       &       npolytyp, &
       &       ifhfvol, &
       &       nhfregion, hfzpos1, hfzpos2, &
       &       hftyp_pmole, hftyp_water, &
       &       ifcalmtf, ifcalhtf, &
       &       mtfoutdir)
#endif

!---- read script file for bias potential calculation

  if (ifpotbias) then
     call rdpotbias( ouumbname, &
          &          xref, eref )
  endif

!---- intialize random generator ----

!  call sgrnd( randseed, mt, mti )  !!! old MT random generator
  call init_gen_rand(randseed)     ! using SFMT random generator

!-------- Open i/o files --------

  call openfile(iucor,iucorname,iutop,iutopname, &
       &        iuwtop,iuwtopname, &
       &        iuparavdw,iuparavdwname,iuparabond,iuparabondname, &
       &        iuparaconst,iuparaconstname, &
       &        iuparacstmnb,iuparacstmnbname,ifcstmnb, &
       &        iuaddtop,iuaddtopname,ifrdaddtop, &
       &        iustrmvel,iustrmvelname,ifstrmvel, &
       &        iuposres,iuposresname,ifposres, &
       &        iostarec,iostarecname,ifstarec, &
       &        npolytyp, &
       &        ousum,ousumname, &
       &        ouene,ouenename,oupos,ouposname,ouvel,ouvelname, &
       &        oufor,ouforname, &
       &        outhe,outhename,oubar,oubarname, &
       &        oupre,ouprename,outhc,outhcname,ifoutthc, &
       &        oupdb,oupdbname,ifoutpdb, &
       &        ouhtf,ouhtfname, &
       &        oumtf,oumtfname, &
       &        ouumb,ouumbname,ifpotbias, &
       &        ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
       &        ifoutbar,ifoutpre)

!-------- Calculate base value for non-dimensionalization --------

  call calbase(xref,eref,mref,qref, &
       &       vref,timeref,tempref,pref,fref, &
       &       eps0ref, &
       &       npolytyp,nmatyp, &
       &       eps0, &
       &       alpha,rrcut, &
       &       xcel,ycel,zcel, &
       &       rcut, &
       &       tcont_poly,tcont_water,tcont_ma, &
       &       tcont_poly_ini,tcont_water_ini,tcont_ma_ini, &
       &       text,tfreq, &
       &       vfreq,pext, &
       &       rcutmor, &
       &       rcutsh, &
       &       rcutrfhfo,rcutrfhoo,rcutrfhoh, &
       &       rcutdouo,rcutindouo,rcutdouh, &
       &       rcutrpvw, &
       &       d_rini,d_rmax,d_econv,d_rmsf, &
       &       limitdist)

!-------- Read topology & coordinate data of molecules --------

  call rdcor( iucor, xref, &
       &      npoly, npolytyp, npoly_mole, npoly_atom, &
       &      nwater, nmatom, &
       &      nmatyp, nmatomtyp, &
       &      polytyp_free, &
       &      oatmtyp, hatmtyp, monoatmtyp, &
       &      ifsetchrg, atmchrg_tmp )
  call rdtop( iutop, iuwtop, &
       &      npolytyp, npoly_mole, npoly_atom, &
       &      nwater )

!-------- Read bond & vdw parameter --------

  call rdpara(iuparavdw,iuparabond, &
       &      xref,eref,mref, &
       &      xcel,ycel,zcel, &
       &      ifljari,ifljgeo)

!-------- Read const parameter --------

  if (ifrattle) then
     call rdconst( iuparaconst )
  end if

!-------- Make pointer to atoms belonging to molecules --------

  call mkmolept( npoly, npolytyp, npoly_mole, npoly_atom, &
       &         nwater, nmatom )

!-------- Read reference coordinate for position restraint --------

  if (ifposres) then
     call rdrefcor( iuposres, &
          &         npoly, npolytyp, npoly_mole, npoly_atom, &
          &         nwater, nmatom, &
          &         xref )
  end if

!-------- Create coordinate --------

  if (ifstarec) then
     call rdstarec(iostarec, &
          &        npoly, npolytyp, npoly_mole, npoly_atom, &
          &        nwater, nmatom, &
          &        xcel, ycel, zcel, yratio, zratio, &
          &        xref, vref, timeref, pref, &
          &        mchain, &
          &        pint, pintt, &
          &        ifsetcor)
  end if

  if (ifcreatecor) then
     call createcor(npoly,npolytyp,npoly_mole,npoly_atom,   &
          &         nwater,nmatom,nmatyp,nmatomtyp,   &
          &         xmaxpo,ymaxpo,zmaxpo,   &
          &         xmaxw,ymaxw,zmaxw,   &
          &         xmaxma,ymaxma,zmaxma,   &
          &         inicorpo,inicorw,inicorma,   &
          &         ncrecorpo,index_crecorpo,   &
          &         ncrecorw,index_crecorw,   &
          &         ncrecorma,index_crecorma,   &
          &         ifsetcor,   &
          &         xcel,ycel,zcel,   &
          &         xref,   &
          &         compfact,   &
          &         mchain)
  end if

!---- read script file for temp. control

  if (ifregionheat) then
     call rdtcont( xref, tempref, &
          &        xcel, ycel, zcel, &
          &        ntcregion, &
          &        tcxpos1, tcxpos2, &
          &        tcypos1, tcypos2, &
          &        tczpos1, tczpos2, &
          &        r_tcont)
  end if

!---- read script file for heat flux control

  if (ifregionhf) then
     call rdhfcont( xref, eref, timeref, &
          &         xcel, ycel, zcel, &
          &         dt_long_cal, &
          &         tcontinterval, &
          &         nhfcregion, &
          &         hfczpos1, hfczpos2, &
          &         r_hfcont)
  end if

!---- read script file for Langevin thermostat

  if (ifreglange) then
     call rdlangevin(xref, tempref, timeref, &
          &          xcel, ycel, zcel, &
          &          nlangeregion, &
          &          ltxpos1, ltxpos2, &
          &          ltypos1, ltypos2, &
          &          ltzpos1, ltzpos2, &
          &          r_ltemp, r_ltdamp)
  end if

!-------- read and add extra bonded pair --------

  ifsetatmmass(1:natom) = .false.
  ifsetatmchrg(1:natom) = .false.

  if (ifrdaddtop) then

     call rdaddtop(iuaddtop, &
          &        mref, &
          &        ifsetatmmass,ifsetatmchrg)

  end if

!-------- read streaming velocities --------

  if (ifstrmvel) then

     call rdstrmvel(iustrmvel, &
          &         xref,vref, &
          &         xcel,ycel,zcel)

  end if

!-------- Make excluded pair list --------

  call mkexcl()

!-------- Link bond & bond parameter --------

  call linkbond(npolytyp,npoly_mole,npoly_atom, &
       &        polytyp_free, &
       &        atmchrg_tmp, &
       &        ifsetatmmass,ifsetatmchrg)

!-------- Link constraint parameter --------

  if (ifrattle) then
     call mkconst( npolytyp, npoly_mole )
  end if

!-------- Read custom NB parameter --------

  if (ifcstmnb) then
     call rdpara_cstmnb(iuparacstmnb, &
          &             ifcellindex_cstmnb,ifbookcstmnb, &
          &             xref,eref,mref,qref, &
          &             vref,timeref,tempref,pref,fref,eps0ref, &
          &             xcel,ycel,zcel)
  end if

!-------- Insert initial velocity --------

  if (ifcreatecor) then
     call createvel( npoly, npolytyp, npoly_mole, npoly_atom, &
          &          nwater, nmatom, nmatyp, nmatomtyp, &
          &          xmaxpo, ymaxpo, zmaxpo, &
          &          xmaxw, ymaxw, zmaxw, &
          &          xmaxma, ymaxma, zmaxma, &
          &          ncrecorpo, index_crecorpo, &
          &          ncrecorw, index_crecorw, &
          &          ncrecorma, index_crecorma, &
          &          vref, timeref, &
          &          mchain, &
          &          tcont_poly, tcont_water, tcont_ma )
  end if

!-------- Calculate mass of each molecule --------

  call calmolmass()

!-------- Preparation for local heating --------

  if (iflocalheat) then

     call prelocalheat( npolytyp, nmatyp, &
          &             polytyp_free,watertyp_free,matomtyp_free, &
          &             tempref, &
          &             nlheat_poly, index_nlheat_poly, &
          &             tcont_nlheat_poly, &
          &             nlheat_water, tcont_nlheat_water, &
          &             nlheat_ma, index_nlheat_ma, &
          &             tcont_nlheat_ma )

  end if

!-------- Preparation for local fix atoms --------

  if (iflocalfix) then

     call prelocalfix( npolytyp,npoly_mole,npoly_atom, &
          &            nwater, &
          &            nmatyp, nmatomtyp, &
          &            polytyp_free, watertyp_free, matomtyp_free, &
          &            nlfix,index_nlfix, &
          &            nlfix_deg_poly, &
          &            nlfix_deg_water, &
          &            nlfix_deg_ma )

  end if

!-------- Preparation for position fix of z coordinate --------

  if (iflocalfixz) then

     call prelocalfixz(xref,   &
          &            npolytyp, npoly_mole, npoly_atom,   &
          &            nwater,   &
          &            nmatyp,nmatomtyp,   &
          &            polytyp_free, watertyp_free, matomtyp_free,   &
          &            nlfixz_deg_poly, nlfixz_deg_water, nlfixz_deg_ma)

  end if

!-------- Preparation for position fix of COM of molecules
!                                                at a z coordinate --------

  if (iflocalfixzg) then

     call prelocalfixzg(xref,   &
          &             npolytyp,npoly_mole,npoly_atom,   &
          &             nwater,   &
          &             nmatyp,nmatomtyp,   &
          &             polytyp_free,watertyp_free,matomtyp_free,   &
          &             nlfixzg_deg_poly,nlfixzg_deg_water,nlfixzg_deg_ma)

  end if

!-------- Preparation for position restraint --------

  if (ifposres) then

     call preposres( npolytyp, npoly_mole, npoly_atom, &
          &          nwater, &
          &          nmatyp, nmatomtyp, &
          &          polytyp_free, watertyp_free, matomtyp_free )

  end if

!-------- Preparation for PDB output --------

  if (ifoutpdb) then

     call prepdbout( npolytyp,polytyp_free, &
          &          watertyp_free, &
          &          nmatyp,matomtyp_free, monoatmtyp, &
          &          resname_poly_pdb, &
          &          resname_water_pdb, &
          &          resname_matom_pdb)

  end if

!-------- Preparation for heat flux calculation --------

#if defined(HF)
  call preheatf( npoly, npolytyp, npoly_mole, &
       &         nwater, &
       &         nmatom, &
       &         hftyp_pmole, hftyp_water, &
       &         hftyp_atm)
#endif

!-------- Preparation for bias potential --------

  if (ifpotbias) then
     call prepotbias( npolytyp, npoly_mole, npoly_atom, &
          &           nwater, &
          &           nmatyp, nmatomtyp, &
          &           polytyp_free, watertyp_free, matomtyp_free)
  end if

!-------- Preparation for fixing local atomic velocity --------

  if (iflocalvel) then
     call prelocalvel(vref, &
          &           npolytyp,npoly_mole,npoly_atom, &
          &           nwater, &
          &           nmatyp,nmatomtyp, &
          &           polytyp_free,watertyp_free,matomtyp_free, &
          &           nlvel,index_nlvel,v_nlvel, &
          &           nlvel_deg_poly, &
          &           nlvel_deg_water, &
          &           nlvel_deg_ma)
  end if

!-------- Preparation of MD --------

  call prepmd(npoly, npolytyp, npoly_mole, npoly_atom, &
       &      nwater, nmatom, nmatyp, nmatomtyp, &
       &      ifcenterfix_all, &
       &      ifcenterfix_poly, ifcenterfix_water, ifcenterfix_ma, &
       &      cenfix_free, &
       &      ifcenterfix_polytyp, &
       &      ifcenterfix_watertyp, &
       &      ifcenterfix_matyp, &
       &      degfree_poly, degfree_water, &
       &      degfree_ma, degfree_all, &
       &      iflocalfix, iflocalfixz, iflocalfixzg, &
       &      nlfix_deg_poly, &
       &      nlfix_deg_water, &
       &      nlfix_deg_ma, &
       &      nlfixz_deg_poly, nlfixz_deg_water, nlfixz_deg_ma,   &
       &      nlfixzg_deg_poly, nlfixzg_deg_water, nlfixzg_deg_ma, &
       &      iflocalvel, &
       &      nlvel_deg_poly,nlvel_deg_water,nlvel_deg_ma, &
       &      dt_short_cal, dt_med_cal, dt_long_cal, &
       &      nstep_short, nstep_med, &
       &      xref, timeref, eps0, &
       &      rcut, ifbook, rcut_book, nstep_book, &
       &      rcutmor, ifbookmor, rcut_bookmor, nstep_bookmor, &
       &      rcutsh, ifbooksh, rcut_booksh, nstep_booksh, &
       &      rcutrfhfo, ifbookrfhfo, rcut_bookrfhfo, &
       &      nstep_bookrfhfo, &
       &      rcutrfhoo, ifbookrfhoo, rcut_bookrfhoo, &
       &      nstep_bookrfhoo, &
       &      rcutrfhoh, ifbookrfhoh, rcut_bookrfhoh, &
       &      nstep_bookrfhoh, &
       &      rcutdouo, ifbookdouo, rcut_bookdouo, &
       &      nstep_bookdouo, &
       &      rcutdouh, ifbookdouh, rcut_bookdouh, &
       &      nstep_bookdouh, &
       &      rcutrpvw, ifbookrpvw, rcut_bookrpvw, nstep_bookrpvw, &
       &      tcont_poly, tcont_water, tcont_ma, &
       &      xcel, ycel, zcel, &
       &      ifrattle, &
       &      mchain, tfreq, text, &
       &      vfreq, &
       &      pcont_axis, &
       &      nyosh, &
       &      solvetyp, solveindex, &
       &      netchrgsq)

!-------- Preparation of ewald wave --------

!- if -D_HF_BULK_EWK used for heat flux calculation, force to execute prepewk
#if defined(_HF_BULK_EWK)
  ! no execusion
#else
  if (ifewald) then
#endif
     call prepewk(pot_ewc, alpha, kmax, &
          &       xcel, ycel, zcel, &
          &       yratio, zratio)
#if defined(_HF_BULK_EWK)
  ! no execusion
#else
  end if
#endif

!-------- Preparation for SPME --------

  if (ifspme) then
     call fft_pme_init( nfft1, nfft2, nfft3, pme_order, &
          &             alpha, &
          &             pot_ewc )
     call erf_corr_cutoff( xcel, ycel, zcel, &
          &                nfft1, nfft2, nfft3)

  end if

!-------- Preparation for Fennell --------

  call prepfennell( pot_ewc, alpha, rrcut, &
       &            iffennell)

!-------- Preparation for spline table

#if defined(_NOT_SPLINTERP)
  ! no execusion
#else
  if (ifewald .or. ifspme .or. iffennell) then
     call spline_interp_init( alpha, &
          &                   rrcut)
  end if
#endif

!-------- Write summarization of MD --------

  if (irank == 0) then

     call wrsumm(ousum, &
          &      xref,eref,mref,qref, &
          &      vref,timeref,tempref,pref,fref, &
          &      eps0ref, &
          &      xcel,ycel,zcel, &
          &      ifewald,alpha,kmax,rrcut, &
          &      ifspme,nfft1,nfft2,nfft3,pme_order, &
          &      iffennell, &
          &      ifljari,ifljgeo, &
          &      iflocalheat, &
          &      ifregionheat,ifregionhf,ifreglange,iftcratom,ifoutthc, &
          &      iflocalfix,iflocalfixz,iflocalfixzg, &
          &      ifcnp, &
          &      ifposres, &
          &      ifpotbias, &
          &      iflocalvel, ifstrmvel, &
          &      dt_long_cal,dt_med_cal,dt_short_cal, &
          &      nstep_med,nstep_short, &
          &      maxnstep,nstage,nstep_stage,mdcont_stage, &
          &      nstep_maxwell,nstep_expand,r_expand, &
          &      eps0, &
          &      rcut,ifbook,rcut_book,nstep_book, &
          &      ifcellindex, &
          &      iucorname,iutopname, &
          &      npoly,npolytyp,npoly_mole,npoly_atom, &
          &      nwater, &
          &      nmatom,nmatyp,nmatomtyp, &
          &      degfree_poly,degfree_water, &
          &      degfree_ma,degfree_all, &
          &      ifstarec,ifcreatecor, &
          &      ifrdaddtop, &
          &      ifoutpdb,nstep_pdbout, &
          &      tcont_poly,tcont_water,tcont_ma, &
          &      tcont_poly_ini,tcont_water_ini,tcont_ma_ini, &
          &      tcontinterval,outinterval, &
          &      pressinterval,heatfinterval, &
          &      ifcalpremole,ifcalpreatom, &
          &      ifnetqcorrp, &
          &      oatmtyp,hatmtyp,monoatmtyp, &
          &      compfact, &
          &      ifrattle,eps_rattle, &
          &      rcutmor,ifbookmor,rcut_bookmor,nstep_bookmor, &
          &      ifcellindex_mor, &
          &      rcutsh,ifbooksh,rcut_booksh,nstep_booksh, &
          &      ifcellindex_sh, &
          &      rcutrfhfo,ifbookrfhfo,rcut_bookrfhfo, &
          &      nstep_bookrfhfo, &
          &      ifcellindex_rfhfo, &
          &      rcutrfhoo,ifbookrfhoo,rcut_bookrfhoo, &
          &      nstep_bookrfhoo, &
          &      ifcellindex_rfhoo, &
          &      rcutrfhoh,ifbookrfhoh,rcut_bookrfhoh, &
          &      nstep_bookrfhoh, &
          &      ifcellindex_rfhoh, &
          &      rcutdouo,rcutindouo,ifbookdouo,rcut_bookdouo, &
          &      nstep_bookdouo, &
          &      ifcellindex_douo, &
          &      rcutdouh,ifbookdouh,rcut_bookdouh, &
          &      nstep_bookdouh, &
          &      ifcellindex_douh, &
          &      rcutrpvw,ifbookrpvw,rcut_bookrpvw, &
          &      nstep_bookrpvw, &
          &      ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
          &      mchain,tfreq,text, &
          &      vfreq,pext,ifpatmcont,ifpmolcont, &
          &      pcont_axis, &
          &      next,nyosh, &
          &      nstep_vir, &
          &      iflimitmove,limitdist, &
          &      ifcalljlong,solvetyp,nsolve, &
          &      d_rini,d_rmax,d_econv,d_rmsf)

  end if

!     End of Preparation

#if defined(MPI)
  call mpi_barrier( MPI_COMM_WORLD, ierror )
  write(6,*) ' process',irank,' end of preparation'
#endif

!-------- Run MD --------

  inistep = 0
  endstep = 0

  DO istage = 1, nstage

!        end timestep update
     endstep = endstep + nstep_stage(istage)

     md_cont = mdcont_stage(istage)

!        ----- NVE-MD with Woodcock temperature control -----
     if ((md_cont == md_0k) .or. (md_cont == md_h) .or. &
       & (md_cont == md_t)  .or. (md_cont == md_nve)) then
        call moldyn(iostarec, recinterval, &
             &      npoly, npolytyp, npoly_mole, npoly_atom, &
             &      nwater, &
             &      nmatom, nmatyp, nmatomtyp, &
             &      maxnstep, inistep, endstep, &
             &      xcel, ycel, zcel, &
             &      ifcenterfix_all, &
             &      ifcenterfix_poly, &
             &      ifcenterfix_water, &
             &      ifcenterfix_ma, &
             &      cenfix_free, &
             &      ifcenterfix_polytyp, &
             &      ifcenterfix_watertyp, &
             &      ifcenterfix_matyp, &
             &      dt_long_cal, &
             &      nstep_med, nstep_short, &
             &      nstep_vir, &
             &      iflimitmove,limitdist, &
             &      mts_bond, mts_angl, mts_anglub, &
             &      mts_tors, mts_torsrb, mts_torsim, &
             &      mts_vdw, mts_ewr, mts_ewk, mts_vdw14, mts_elc14, &
             &      mts_mor, mts_sh, mts_rfh, mts_dou, mts_cnpvw, &
             &      mts_cstmnb, &
             &      mts_posres, &
             &      mts_potbias, &
             &      ifewald, alpha, kmax, rrcut, &
             &      ifspme, nfft1, nfft2,nfft3, pme_order, &
             &      eps0, div_factor_14vdw, div_factor_14elc, &
             &      rcut, &
             &      ifcellindex, &
             &      ifbook, nstep_book, &
             &      tcont_poly, tcont_water, tcont_ma, &
             &      tcont_poly_ini, tcont_water_ini, tcont_ma_ini, &
             &      ifrattle, eps_rattle, &
             &      iflocalheat, &
             &      nlheat_poly, index_nlheat_poly, &
             &      tcont_nlheat_poly, &
             &      nlheat_water, tcont_nlheat_water, &
             &      nlheat_ma, index_nlheat_ma, &
             &      tcont_nlheat_ma, &
             &      ifregionheat, ifregionhf, ifreglange, &
             &      iftcratom, ifoutthc, &
             &      ntcregion, &
             &      tcxpos1, tcxpos2, &
             &      tcypos1, tcypos2, &
             &      tczpos1, tczpos2, &
             &      r_tcont, &
             &      nhfcregion, &
             &      hfczpos1, hfczpos2, &
             &      r_hfcont, &
             &      nlangeregion, &
             &      ltxpos1, ltxpos2, &
             &      ltypos1, ltypos2, &
             &      ltzpos1, ltzpos2, &
             &      r_ltemp, r_ltdamp, &
             &      iflocalfix, &
             &      nlfix, index_nlfix, &
             &      iflocalfixz, iflocalfixzg, &
             &      ifcnp, &
             &      ifposres, &
             &      ifpotbias, &
             &      iflocalvel, ifstrmvel, &
             &      nlvel, index_nlvel, v_nlvel, &
             &      ifoutpdb, nstep_pdbout, &
             &      resname_poly_pdb, &
             &      resname_water_pdb, &
             &      resname_matom_pdb, &
             &      xref, eref, mref, qref, &
             &      vref, timeref, tempref, pref, fref, eps0ref, &
             &      degfree_poly, degfree_water, &
             &      degfree_ma, degfree_all, &
             &      dt_med_cal, dt_short_cal, &
             &      rcut_book, &
             &      rcutmor,ifbookmor, rcut_bookmor, nstep_bookmor, &
             &      ifcellindex_mor, &
             &      rcutsh, ifbooksh,rcut_booksh, nstep_booksh, &
             &      ifcellindex_sh, &
             &      rcutrfhfo, ifbookrfhfo, rcut_bookrfhfo, &
             &      nstep_bookrfhfo, &
             &      ifcellindex_rfhfo, &
             &      rcutrfhoo, ifbookrfhoo, rcut_bookrfhoo, &
             &      nstep_bookrfhoo, &
             &      ifcellindex_rfhoo, &
             &      rcutrfhoh, ifbookrfhoh, rcut_bookrfhoh, &
             &      nstep_bookrfhoh, &
             &      ifcellindex_rfhoh, &
             &      rcutdouo, rcutindouo, ifbookdouo, rcut_bookdouo, &
             &      nstep_bookdouo, &
             &      ifcellindex_douo, &
             &      rcutdouh, ifbookdouh,rcut_bookdouh, &
             &      nstep_bookdouh, &
             &      ifcellindex_douh, &
             &      rcutrpvw, ifbookrpvw, rcut_bookrpvw, &
             &      nstep_bookrpvw, &
             &      ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
             &      for_long, for_short, &
             &      for_viric_long, for_viric_med, for_viric_short, &
             &      for_virilj_long, for_virilj_med, for_virilj_short, &
             &      for_virimor_long, for_virimor_med, &
             &      for_virimor_short, &
             &      for_virish_long, for_virish_med, &
             &      for_virish_short, &
             &      for_virirfh_long, for_virirfh_med, &
             &      for_virirfh_short, &
             &      for_viridou_long, for_viridou_med, &
             &      for_viridou_short, &
             &      for_viricstmnb_long, for_viricstmnb_med, &
             &      for_viricstmnb_short, &
             &      pot_ewc, &
             &      ouene,oupos,ouvel,oufor,outhe,oubar,oupre,outhc,oupdb, &
             &      ouumb, &
             &      ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
             &      ifoutbar,ifoutpre, &
             &      tcontinterval, outinterval, pressinterval, &
             &      yratio, zratio, &
             &      ifcalpremole, ifcalpreatom, &
             &      ifnetqcorrp, &
             &      mchain, &
             &      pint, pintt, &
             &      ifpatmcont, ifpmolcont, &
             &      ifcalljlong, nsolve, solveindex, &
             &      netchrgsq, &
             &      nstep_expand,r_expand, &
             &      nstep_maxwell, &
             &      md_cont)
     end if

!     ----- NVT-MD with Nose-Hoover chain thermostat -----
     if (md_cont == md_nhc) then
        call moldyn_nhc(iostarec, recinterval, &
             &          npoly, npolytyp, npoly_mole, npoly_atom, &
             &          nwater, &
             &          nmatom, nmatyp, nmatomtyp, &
             &          maxnstep, inistep, endstep, &
             &          xcel, ycel, zcel, &
             &          ifcenterfix_all, &
             &          ifcenterfix_poly, &
             &          ifcenterfix_water, &
             &          ifcenterfix_ma, &
             &          cenfix_free, &
             &          ifcenterfix_polytyp, &
             &          ifcenterfix_watertyp, &
             &          ifcenterfix_matyp, &
             &          dt_long_cal, &
             &          nstep_med, nstep_short, &
             &          nstep_vir, &
             &          iflimitmove,limitdist, &
             &          mts_bond, mts_angl, mts_anglub, &
             &          mts_tors, mts_torsrb, mts_torsim, &
             &          mts_vdw, mts_ewr, mts_ewk, mts_vdw14, mts_elc14, &
             &          mts_mor, mts_sh, mts_rfh, mts_dou, mts_cnpvw, &
             &          mts_cstmnb, &
             &          mts_posres, &
             &          mts_potbias, &
             &          ifewald, alpha, kmax,rrcut, &
             &          ifspme, nfft1, nfft2, nfft3, pme_order, &
             &          eps0, div_factor_14vdw, div_factor_14elc, &
             &          rcut, &
             &          ifcellindex, &
             &          ifbook, nstep_book, &
             &          tcont_poly, tcont_water, tcont_ma, &
             &          ifrattle, eps_rattle, &
             &          iftcratom, &
             &          iflocalfix, &
             &          nlfix, index_nlfix, &
             &          iflocalfixz, iflocalfixzg, &
             &          ifcnp, &
             &          ifposres, &
             &          ifpotbias, &
             &          iflocalvel, ifstrmvel, &
             &          nlvel, index_nlvel, v_nlvel, &
             &          ifoutpdb, nstep_pdbout, &
             &          resname_poly_pdb, &
             &          resname_water_pdb, &
             &          resname_matom_pdb, &
             &          xref, eref, mref, qref, &
             &          vref, timeref, tempref, pref, fref, eps0ref, &
             &          degfree_poly, degfree_water, &
             &          degfree_ma, degfree_all, &
             &          dt_med_cal, dt_short_cal, &
             &          rcut_book, &
             &          rcutmor, ifbookmor, rcut_bookmor, nstep_bookmor, &
             &          ifcellindex_mor, &
             &          rcutsh, ifbooksh, rcut_booksh, nstep_booksh, &
             &          ifcellindex_sh, &
             &          rcutrfhfo, ifbookrfhfo, rcut_bookrfhfo, &
             &          nstep_bookrfhfo, &
             &          ifcellindex_rfhfo, &
             &          rcutrfhoo, ifbookrfhoo, rcut_bookrfhoo, &
             &          nstep_bookrfhoo, &
             &          ifcellindex_rfhoo, &
             &          rcutrfhoh, ifbookrfhoh, rcut_bookrfhoh, &
             &          nstep_bookrfhoh, &
             &          ifcellindex_rfhoh, &
             &          rcutdouo, rcutindouo, ifbookdouo, rcut_bookdouo, &
             &          nstep_bookdouo, &
             &          ifcellindex_douo, &
             &          rcutdouh, ifbookdouh, rcut_bookdouh, &
             &          nstep_bookdouh, &
             &          ifcellindex_douh, &
             &          rcutrpvw, ifbookrpvw, rcut_bookrpvw, &
             &          nstep_bookrpvw, &
             &          ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
             &          for_long, for_short, &
             &          for_viric_long, for_viric_med, for_viric_short, &
             &          for_virilj_long, for_virilj_med, for_virilj_short, &
             &          for_virimor_long, for_virimor_med, &
             &          for_virimor_short, &
             &          for_virish_long, for_virish_med, &
             &          for_virish_short, &
             &          for_virirfh_long, for_virirfh_med, &
             &          for_virirfh_short, &
             &          for_viridou_long, for_viridou_med, &
             &          for_viridou_short, &
             &          for_viricstmnb_long, for_viricstmnb_med, &
             &          for_viricstmnb_short, &
             &          pot_ewc, &
             &          ouene,oupos,ouvel,oufor, &
             &          outhe,oubar,oupre,oupdb, &
             &          ouumb, &
             &          ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
             &          ifoutbar,ifoutpre, &
             &          tcontinterval, outinterval, pressinterval, &
             &          yratio, zratio, &
             &          ifcalpremole, ifcalpreatom, &
             &          ifnetqcorrp, &
             &          mchain, text, &
             &          pint, pintt, &
             &          ifpatmcont, ifpmolcont, &
             &          next, nyosh, &
             &          ifcalljlong, nsolve, solveindex, &
             &          netchrgsq, &
             &          nstep_expand, r_expand, &
             &          md_cont)
     end if

!     ----- NPT-MD by Martyna-Tobias-Klein eqs. of motion -----
     if (md_cont == md_mtk) then
        call moldyn_mtk(iostarec, recinterval, &
             &          npoly,npolytyp, npoly_mole, npoly_atom, &
             &          nwater, &
             &          nmatom, nmatyp, nmatomtyp, &
             &          maxnstep, inistep, endstep, &
             &          xcel, ycel, zcel, &
             &          ifcenterfix_all, &
             &          ifcenterfix_poly, &
             &          ifcenterfix_water, &
             &          ifcenterfix_ma, &
             &          cenfix_free, &
             &          ifcenterfix_polytyp, &
             &          ifcenterfix_watertyp, &
             &          ifcenterfix_matyp, &
             &          dt_long_cal, &
             &          nstep_med, nstep_short, &
             &          nstep_vir, &
             &          iflimitmove,limitdist, &
             &          mts_bond, mts_angl, mts_anglub, &
             &          mts_tors, mts_torsrb, mts_torsim, &
             &          mts_vdw, mts_ewr, mts_ewk, mts_vdw14, mts_elc14, &
             &          mts_mor, mts_sh, mts_rfh, mts_dou, mts_cnpvw, &
             &          mts_cstmnb, &
             &          mts_posres, &
             &          mts_potbias, &
             &          ifewald, alpha, kmax, rrcut, &
             &          ifspme, nfft1, nfft2, nfft3, pme_order, &
             &          eps0, div_factor_14vdw, div_factor_14elc, &
             &          rcut, &
             &          ifcellindex, &
             &          ifbook, nstep_book, &
             &          tcont_poly, tcont_water, tcont_ma, &
             &          ifrattle, eps_rattle, &
             &          iflocalfix, &
             &          nlfix, index_nlfix, &
             &          iflocalfixz, iflocalfixzg, &
             &          ifposres, &
             &          ifpotbias, &
             &          iflocalvel, &
             &          nlvel, index_nlvel, v_nlvel, &
             &          ifoutpdb, nstep_pdbout, &
             &          resname_poly_pdb, &
             &          resname_water_pdb, &
             &          resname_matom_pdb, &
             &          xref,eref,mref,qref, &
             &          vref, timeref, tempref, pref, fref, eps0ref, &
             &          degfree_poly, degfree_water, &
             &          degfree_ma, degfree_all, &
             &          dt_med_cal, dt_short_cal, &
             &          rcut_book, &
             &          rcutmor, ifbookmor, rcut_bookmor, nstep_bookmor, &
             &          ifcellindex_mor, &
             &          rcutsh, ifbooksh, rcut_booksh, nstep_booksh, &
             &          ifcellindex_sh, &
             &          rcutrfhfo, ifbookrfhfo, rcut_bookrfhfo, &
             &          nstep_bookrfhfo, &
             &          ifcellindex_rfhfo, &
             &          rcutrfhoo, ifbookrfhoo, rcut_bookrfhoo, &
             &          nstep_bookrfhoo, &
             &          ifcellindex_rfhoo, &
             &          rcutrfhoh, ifbookrfhoh, rcut_bookrfhoh, &
             &          nstep_bookrfhoh, &
             &          ifcellindex_rfhoh, &
             &          rcutdouo, rcutindouo, ifbookdouo, rcut_bookdouo, &
             &          nstep_bookdouo, &
             &          ifcellindex_douo, &
             &          rcutdouh, ifbookdouh, rcut_bookdouh, &
             &          nstep_bookdouh, &
             &          ifcellindex_douh, &
             &          rcutrpvw, ifbookrpvw, rcut_bookrpvw, &
             &          nstep_bookrpvw, &
             &          ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
             &          for_long, for_short, &
             &          for_viric_long, for_viric_med, for_viric_short, &
             &          for_virilj_long, for_virilj_med, for_virilj_short, &
             &          for_virimor_long, for_virimor_med, &
             &          for_virimor_short, &
             &          for_virish_long, for_virish_med, &
             &          for_virish_short, &
             &          for_virirfh_long, for_virirfh_med, &
             &          for_virirfh_short, &
             &          for_viridou_long, for_viridou_med, &
             &          for_viridou_short, &
             &          for_viricstmnb_long, for_viricstmnb_med, &
             &          for_viricstmnb_short, &
             &          pot_ewc, &
             &          ouene,oupos,ouvel,oufor, &
             &          outhe,oubar,oupre,oupdb, &
             &          ouumb, &
             &          ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
             &          ifoutbar,ifoutpre, &
             &          tcontinterval, outinterval, pressinterval, &
             &          yratio,zratio, &
             &          ifcalpremole, ifcalpreatom, &
             &          ifnetqcorrp, &
             &          mchain,text, &
             &          pext, pint, pintt, &
             &          ifpatmcont, ifpmolcont, &
             &          pcont_axis, &
             &          next,nyosh, &
             &          ifcalljlong, nsolve, solveindex, &
             &          netchrgsq, &
             &          nstep_expand, r_expand, &
             &          md_cont)
     end if

!     ----- Heat flux calculation in NVE ensemble -----

#if defined(HF)
     if (md_cont == md_htf) then
        call moldyn_hf(iostarec, recinterval, &
             &         npoly, npolytyp, npoly_mole, npoly_atom, &
             &         nwater, &
             &         nmatom, nmatyp, nmatomtyp, &
             &         maxnstep, inistep, endstep, &
             &         xcel, ycel, zcel, &
             &         ifcenterfix_all, &
             &         ifcenterfix_poly, &
             &         ifcenterfix_water, &
             &         ifcenterfix_ma, &
             &         cenfix_free, &
             &         ifcenterfix_polytyp, &
             &         ifcenterfix_watertyp, &
             &         ifcenterfix_matyp, &
             &         dt_long_cal, &
             &         nstep_med, nstep_short, &
             &         nstep_vir, &
             &         iflimitmove,limitdist, &
             &         mts_bond, mts_angl, mts_anglub, &
             &         mts_tors, mts_torsrb, mts_torsim, &
             &         mts_vdw, mts_ewr,mts_ewk, mts_vdw14, mts_elc14, &
             &         mts_mor, mts_sh, mts_rfh, mts_dou, mts_cnpvw, &
             &         mts_cstmnb, &
             &         mts_posres, &
             &         mts_potbias, &
             &         ifewald, alpha, kmax, rrcut, &
             &         ifspme, nfft1, nfft2, nfft3, pme_order, &
             &         eps0, div_factor_14vdw, div_factor_14elc, &
             &         rcut, &
             &         ifcellindex, &
             &         ifbook, nstep_book, &
             &         tcont_poly, tcont_water, tcont_ma, &
             &         tcont_poly_ini, tcont_water_ini, tcont_ma_ini, &
             &         ifrattle, eps_rattle, &
             &         iflocalheat, &
             &         nlheat_poly, index_nlheat_poly, &
             &         tcont_nlheat_poly, &
             &         nlheat_water, tcont_nlheat_water, &
             &         nlheat_ma, index_nlheat_ma, &
             &         tcont_nlheat_ma, &
             &         ifregionheat, ifregionhf, ifreglange, &
             &         iftcratom, ifoutthc, &
             &         ntcregion, &
             &         tcxpos1, tcxpos2, &
             &         tcypos1, tcypos2, &
             &         tczpos1, tczpos2, &
             &         r_tcont, &
             &         nhfcregion, &
             &         hfczpos1, hfczpos2, &
             &         r_hfcont, &
             &         nlangeregion, &
             &         ltxpos1, ltxpos2, &
             &         ltypos1, ltypos2, &
             &         ltzpos1, ltzpos2, &
             &         r_ltemp, r_ltdamp, &
             &         iflocalfix, &
             &         nlfix, index_nlfix, &
             &         iflocalfixz, iflocalfixzg, &
             &         ifcnp, &
             &         ifposres, &
             &         ifpotbias, &
             &         iflocalvel, ifstrmvel, &
             &         nlvel, index_nlvel, v_nlvel, &
             &         ifoutpdb, nstep_pdbout, &
             &         resname_poly_pdb, &
             &         resname_water_pdb, &
             &         resname_matom_pdb, &
             &         xref, eref, mref, qref, &
             &         vref, timeref, tempref, pref, fref, eps0ref, &
             &         degfree_poly, degfree_water, &
             &         degfree_ma, degfree_all, &
             &         dt_med_cal, dt_short_cal, &
             &         rcut_book, &
             &         rcutmor, ifbookmor, rcut_bookmor, nstep_bookmor, &
             &         ifcellindex_mor, &
             &         rcutsh, ifbooksh, rcut_booksh, nstep_booksh, &
             &         ifcellindex_sh, &
             &         rcutrfhfo, ifbookrfhfo, rcut_bookrfhfo, &
             &         nstep_bookrfhfo, &
             &         ifcellindex_rfhfo, &
             &         rcutrfhoo, ifbookrfhoo, rcut_bookrfhoo, &
             &         nstep_bookrfhoo, &
             &         ifcellindex_rfhoo, &
             &         rcutrfhoh, ifbookrfhoh, rcut_bookrfhoh, &
             &         nstep_bookrfhoh, &
             &         ifcellindex_rfhoh, &
             &         rcutdouo, rcutindouo, ifbookdouo, rcut_bookdouo, &
             &         nstep_bookdouo, &
             &         ifcellindex_douo, &
             &         rcutdouh, ifbookdouh, rcut_bookdouh, &
             &         nstep_bookdouh, &
             &         ifcellindex_douh, &
             &         rcutrpvw, ifbookrpvw, rcut_bookrpvw, &
             &         nstep_bookrpvw, &
             &         ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
             &         for_long,for_short, &
             &         for_viric_long, for_viric_med, for_viric_short, &
             &         for_virilj_long, for_virilj_med, for_virilj_short, &
             &         for_virimor_long, for_virimor_med, &
             &         for_virimor_short, &
             &         for_virish_long, for_virish_med, &
             &         for_virish_short, &
             &         for_virirfh_long, for_virirfh_med, &
             &         for_virirfh_short, &
             &         for_viridou_long, for_viridou_med, &
             &         for_viridou_short, &
             &         for_viricstmnb_long, for_viricstmnb_med, &
             &         for_viricstmnb_short, &
             &         pot_ewc, &
             &         ouene,oupos,ouvel,oufor, &
             &         outhe,oubar,oupre,outhc,oupdb, &
             &         ouhtf, &
             &         oumtf, &
             &         ouumb, &
             &         ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
             &         ifoutbar,ifoutpre, &
             &         tcontinterval, outinterval, &
             &         pressinterval, heatfinterval, &
             &         yratio,zratio, &
             &         ifcalpremole, ifcalpreatom, &
             &         ifnetqcorrp, &
             &         mchain, &
             &         pint, pintt, &
             &         ifpatmcont, ifpmolcont, &
             &         ifcalljlong, nsolve, solveindex, &
             &         netchrgsq, &
             &         nstep_expand,r_expand, &
             &         nstep_maxwell, &
             &         ifhfvol, &
             &         nhfregion, hfzpos1, hfzpos2, &
             &         hftyp_pmole, hftyp_water, &
             &         hftyp_atm, &
             &         ifcalmtf,ifcalhtf, &
             &         mtfoutdir, &
             &         md_cont)
     end if
#endif

!     ----- Energy minimization by steepest descent (SD) method -----
     if (md_cont == md_ems) then

        call enemin_sd(iostarec,recinterval, &
             &         npoly,npolytyp,npoly_mole,npoly_atom, &
             &         nwater, &
             &         nmatom,nmatyp,nmatomtyp, &
             &         maxnstep,inistep,endstep, &
             &         xcel,ycel,zcel, &
             &         ifcenterfix_all, &
             &         ifcenterfix_poly, &
             &         ifcenterfix_water, &
             &         ifcenterfix_ma, &
             &         cenfix_free, &
             &         ifcenterfix_polytyp, &
             &         ifcenterfix_watertyp, &
             &         ifcenterfix_matyp, &
             &         mts_bond, mts_angl, mts_anglub, &
             &         mts_tors,mts_torsrb,mts_torsim, &
             &         mts_vdw,mts_ewr,mts_ewk,mts_vdw14,mts_elc14, &
             &         mts_mor,mts_sh,mts_rfh,mts_dou,mts_cnpvw, &
             &         mts_cstmnb, &
             &         mts_posres, &
             &         mts_potbias, &
             &         ifewald,alpha,kmax,rrcut, &
             &         ifspme,nfft1,nfft2,nfft3,pme_order, &
             &         eps0,div_factor_14vdw,div_factor_14elc, &
             &         rcut, &
             &         ifcellindex, &
             &         ifbook,nstep_book, &
             &         ifrattle,eps_rattle, &
             &         iflocalfix, &
             &         nlfix,index_nlfix, &
             &         iflocalfixz,iflocalfixzg, &
             &         ifcnp, &
             &         ifposres, &
             &         ifpotbias, &
             &         ifoutpdb,nstep_pdbout, &
             &         resname_poly_pdb, &
             &         resname_water_pdb, &
             &         resname_matom_pdb, &
             &         xref,eref,mref,qref, &
             &         vref,timeref,tempref,pref,fref,eps0ref, &
             &         degfree_poly, degfree_water, &
             &         degfree_ma, degfree_all, &
             &         rcut_book, &
             &         rcutmor,ifbookmor,rcut_bookmor,nstep_bookmor, &
             &         ifcellindex_mor, &
             &         rcutsh,ifbooksh,rcut_booksh,nstep_booksh, &
             &         ifcellindex_sh, &
             &         rcutrfhfo,ifbookrfhfo,rcut_bookrfhfo, &
             &         nstep_bookrfhfo, &
             &         ifcellindex_rfhfo, &
             &         rcutrfhoo,ifbookrfhoo,rcut_bookrfhoo, &
             &         nstep_bookrfhoo, &
             &         ifcellindex_rfhoo, &
             &         rcutrfhoh,ifbookrfhoh,rcut_bookrfhoh, &
             &         nstep_bookrfhoh, &
             &         ifcellindex_rfhoh, &
             &         rcutdouo,rcutindouo,ifbookdouo,rcut_bookdouo, &
             &         nstep_bookdouo, &
             &         ifcellindex_douo, &
             &         rcutdouh,ifbookdouh,rcut_bookdouh, &
             &         nstep_bookdouh, &
             &         ifcellindex_douh, &
             &         rcutrpvw,ifbookrpvw,rcut_bookrpvw, &
             &         nstep_bookrpvw, &
             &         ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
             &         for_long,for_short, &
             &         for_viric_long,for_viric_med,for_viric_short, &
             &         for_virilj_long,for_virilj_med,for_virilj_short, &
             &         for_virimor_long,for_virimor_med, &
             &         for_virimor_short, &
             &         for_virish_long,for_virish_med, &
             &         for_virish_short, &
             &         for_virirfh_long,for_virirfh_med, &
             &         for_virirfh_short, &
             &         for_viridou_long,for_viridou_med, &
             &         for_viridou_short, &
             &         for_viricstmnb_long, for_viricstmnb_med, &
             &         for_viricstmnb_short, &
             &         pot_ewc, &
             &         ouene,oupos,ouvel,oufor, &
             &         outhe,oubar,oupre,oupdb, &
             &         ouumb, &
             &         ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
             &         ifoutbar,ifoutpre, &
             &         outinterval,pressinterval, &
             &         yratio,zratio, &
             &         ifcalpremole,ifcalpreatom, &
             &         ifnetqcorrp, &
             &         mchain, &
             &         pint,pintt, &
             &         ifpatmcont,ifpmolcont, &
             &         ifcalljlong,nsolve,solveindex, &
             &         netchrgsq, &
             &         nstep_expand,r_expand, &
             &         d_rini,d_rmax,d_econv,d_rmsf, &
             &         md_cont)
     end if

!        initial timestep update
     inistep = inistep + nstep_stage(istage)
     if (istage == 1) inistep = inistep + 1

  END DO

!---- Record the state of the end of this MD --------
  if (irank == 0) then

     call wrsta(iostarec, &
          &     npoly,nwater,nmatom, &
          &     xcel,ycel,zcel, &
          &     xref,vref,timeref,pref, &
          &     maxnstep, &
          &     mchain, &
          &     pint, pintt)

!----- Close output files -----
     if (ifoutene) close(ouene)
     if (ifoutpos) close(oupos)
     if (ifoutvel) close(ouvel)
     if (ifoutfor) close(oufor)
     if (ifoutthe) close(outhe)
     if (ifoutbar) close(oubar)
     if (ifoutpre) close(oupre)
     if (ifoutpdb) close(oupdb)
#if defined(HF)
     close(ouhtf)
     close(oumtf)
#endif
     if (ifpotbias) close(ouumb)

     close(iostarec)

     write(6,*) 'program is done'

  end if

!---- MPI finalize ----

#if defined(MPI)
  call mpi_finalize(ierror)
#endif

!     +     +     +     +     +     +     +

end program peachgk_md
