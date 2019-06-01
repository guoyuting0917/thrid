!*****************************
!*  rdscript.f90 Ver.7.0     *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>

subroutine rdscript(pgkiniver,pgkinidate, &
     &              md_0k,md_h,md_t,md_mtk,md_nhc,md_nve,md_htf,md_ems, &
     &              mts_long,mts_med,mts_short, &
     &              maxnword, &
     &              iucorname,iutopname,iuwtopname, &
     &              iuparavdwname,iuparabondname,iuparaconstname, &
     &              iuparacstmnbname, &
     &              iuaddtopname, &
     &              iustrmvelname, &
     &              iuposresname, &
     &              iostarecname, &
     &              ousumname,ouenename,ouposname,ouvelname,ouforname, &
     &              outhename,oubarname,ouprename,outhcname, &
     &              oupdbname, &
     &              ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
     &              ifoutbar,ifoutpre, &
     &              maxnpolytyp, &
     &              npoly,npolytyp,npoly_mole,npoly_atom, &
     &              nwater, &
     &              maxnmatyp,nmatom,nmatyp,nmatomtyp, &
     &              polytyp_free,watertyp_free,matomtyp_free, &
     &              xmaxpo,ymaxpo,zmaxpo,xmaxw,ymaxw,zmaxw, &
     &              xmaxma,ymaxma,zmaxma, &
     &              inicorpo,inicorw,inicorma, &
     &              ncrecorpo,index_crecorpo, &
     &              ncrecorw,index_crecorw, &
     &              ncrecorma,index_crecorma, &
     &              maxnstage,maxnstep, &
     &              nstage,nstep_stage,mdcont_stage, &
     &              nstep_maxwell,nstep_expand, &
     &              xcel,ycel,zcel,yratio,zratio, &
     &              r_expand, &
     &              ifstarec,ifcreatecor, &
     &              ifrdaddtop, &
     &              ifcenterfix_all, &
     &              ifcenterfix_poly,ifcenterfix_water, &
     &              ifcenterfix_ma, &
     &              cenfix_free, &
     &              ifcenterfix_polytyp, &
     &              ifcenterfix_watertyp, &
     &              ifcenterfix_matyp, &
     &              dt_long_cal, &
     &              nstep_med,nstep_short, &
     &              mts_bond,mts_angl,mts_anglub, &
     &              mts_tors,mts_torsrb,mts_torsim, &
     &              mts_vdw,mts_ewr,mts_ewk, &
     &              mts_vdw14,mts_elc14, &
     &              mts_mor,mts_sh,mts_rfh,mts_dou,mts_cnpvw, &
     &              mts_cstmnb, &
     &              mts_posres, &
     &              mts_potbias, &
     &              nstep_vir, &
     &              iflimitmove,limitdist, &
     &              mchain,tfreq,text, &
     &              vfreq,pext, &
     &              ifpatmcont,ifpmolcont, &
     &              pcont_axis, &
     &              ifnetqcorrp, &
     &              next,nyosh, &
     &              ifrattle, &
     &              ifewald,ifspme,iffennell, &
     &              ifljari,ifljgeo, &
     &              iflocalheat, &
     &              ifregionheat,ifregionhf,ifreglange,iftcratom,ifoutthc, &
     &              iflocalfix, iflocalfixz, &
     &              iflocalfixzg, &
     &              ifcnp, &
     &              ifposres, &
     &              ifpotbias, &
     &              iflocalvel,ifstrmvel, &
     &              ifcalpremole,ifcalpreatom, &
     &              ifcalljlong,solvetyp,nsolve, &
     &              alpha,kmax,rrcut, &
     &              nfft1,nfft2,nfft3,pme_order, &
     &              tcont_poly,tcont_water,tcont_ma, &
     &              tcont_poly_ini,tcont_water_ini,tcont_ma_ini, &
     &              tcontinterval,outinterval, &
     &              pressinterval,heatfinterval, &
     &              recinterval, &
     &              oatmtyp,hatmtyp,monoatmtyp, &
     &              randseed, &
     &              compfact, &
     &              eps_rattle, &
     &              rcut,ifbook,rcut_book,nstep_book, &
     &              ifcellindex, &
     &              rcutmor,ifbookmor,rcut_bookmor,nstep_bookmor, &
     &              ifcellindex_mor, &
     &              rcutsh,ifbooksh,rcut_booksh,nstep_booksh, &
     &              ifcellindex_sh, &
     &              rcutrfhfo,ifbookrfhfo,rcut_bookrfhfo, &
     &              nstep_bookrfhfo, &
     &              ifcellindex_rfhfo, &
     &              rcutrfhoo,ifbookrfhoo,rcut_bookrfhoo, &
     &              nstep_bookrfhoo, &
     &              ifcellindex_rfhoo, &
     &              rcutrfhoh,ifbookrfhoh,rcut_bookrfhoh, &
     &              nstep_bookrfhoh, &
     &              ifcellindex_rfhoh, &
     &              rcutdouo,rcutindouo,ifbookdouo,rcut_bookdouo, &
     &              nstep_bookdouo, &
     &              ifcellindex_douo, &
     &              rcutdouh,ifbookdouh,rcut_bookdouh, &
     &              nstep_bookdouh, &
     &              ifcellindex_douh, &
     &              rcutrpvw, ifbookrpvw, rcut_bookrpvw, &
     &              nstep_bookrpvw, &
     &              ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
     &              ifsetchrg, &
     &              ifoutpdb,nstep_pdbout, &
     &              d_rini,d_rmax,d_econv,d_rmsf, &
     &              nspltbl)

  use interface_tools

  use mpi_global
  use spme_global

  implicit none

! ARGUMENT:
!     INPUT
  character(3),intent(in):: pgkiniver  ! peachgk.ini Version num.
  character(8),intent(in):: pgkinidate ! peachgk.ini rev. date

  integer,intent(in):: md_0k           ! 0[K] NVT for clear distorsion
  integer,intent(in):: md_h            ! gradual heating NVT (v-scale)
  integer,intent(in):: md_t            ! target temperature NVT (v-scale)
  integer,intent(in):: md_mtk          ! NPT constant MD (MTK eq.)
  integer,intent(in):: md_nhc          ! NVT constant MD (NHC eq.)
  integer,intent(in):: md_nve          ! NVE constant MD
  integer,intent(in):: md_htf          ! heat flux calculation in NVE MD
  integer,intent(in):: md_ems          ! EM by steepest descent (SD) method

  integer,intent(in):: mts_long        ! MTS flag for long   force
  integer,intent(in):: mts_med         ! MTS flag for medium force
  integer,intent(in):: mts_short       ! MTS flag for short  force

  integer,intent(in):: maxnword        ! maximum number of words from data file

!     OUTPUT
  character(80),intent(out):: iucorname(:)     ! input file name
  character(80),intent(out):: iutopname(:)     !       "
  character(80),intent(out):: iuwtopname       !       "
  character(80),intent(out):: iuparavdwname    !       "
  character(80),intent(out):: iuparabondname   !       "
  character(80),intent(out):: iuparaconstname  !       "
  character(80),intent(out):: iuparacstmnbname !       "
  character(80),intent(out):: iuaddtopname     !       "
  character(80),intent(out):: iustrmvelname    !       "
  character(80),intent(out):: iuposresname     !       "

  character(80),intent(out):: iostarecname ! state file name

  character(80),intent(out):: ousumname ! output file name
  character(80),intent(out):: ouenename ! output file name
  character(80),intent(out):: ouposname ! output file name
  character(80),intent(out):: ouvelname ! output file name
  character(80),intent(out):: ouforname ! output file name
  character(80),intent(out):: outhename ! output file name
  character(80),intent(out):: oubarname ! output file name
  character(80),intent(out):: ouprename ! output file name
  character(80),intent(out):: outhcname ! output file name
  character(80),intent(out):: oupdbname ! output file name

  logical,intent(out):: ifoutene        ! if ouput energy file
  logical,intent(out):: ifoutpos        ! if ouput position file
  logical,intent(out):: ifoutvel        ! if ouput velocity file
  logical,intent(out):: ifoutfor        ! if ouput force file
  logical,intent(out):: ifoutthe        ! if ouput NVT file
  logical,intent(out):: ifoutbar        ! if ouput NPT file
  logical,intent(out):: ifoutpre        ! if ouput pressure file

  integer,intent(in):: maxnpolytyp
  integer,intent(out):: npoly           ! number of polymer1
  integer,intent(out):: npolytyp        ! number of poly type
  integer,intent(out):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(out):: npoly_atom(:)   ! number of atoms belonging to poly

  integer,intent(out):: nwater          ! number of H2O molecules
  integer,intent(out):: nmatom          ! number of monatomic molecules

  integer,intent(in):: maxnmatyp
  integer,intent(out):: nmatyp          ! number of species of monatomic mole.
  integer,intent(out):: nmatomtyp(:)    ! each number of monatomic mole.

  character(80),intent(out):: polytyp_free(:,:)
                                ! use for poly type control
  character(80),intent(out):: watertyp_free(:) ! use for water type control
  character(80),intent(out):: matomtyp_free(:,:)
                                ! use for matom type control

  integer,intent(out):: xmaxpo(:)       ! use for positioning of polymer1
  integer,intent(out):: ymaxpo(:)       ! use for positioning of polymer1
  integer,intent(out):: zmaxpo(:)       ! use for positioning of polymer1
  integer,intent(out):: xmaxw           ! use for positioning of water
  integer,intent(out):: ymaxw           ! use for positioning of water
  integer,intent(out):: zmaxw           ! use for positioning of water
  integer,intent(out):: xmaxma(:)       ! use for positioning of monatomic mole.
  integer,intent(out):: ymaxma(:)       ! use for positioning of monatomic mole.
  integer,intent(out):: zmaxma(:)       ! use for positioning of monatomic mole.
  real(8),intent(out):: inicorpo(:,:)   ! use for positioning of poly
  real(8),intent(out):: inicorw(:)      ! use for positioning of water
  real(8),intent(out):: inicorma(:,:)   ! use for positioning of ma
  integer,intent(out):: ncrecorpo       ! max number of poly type for createcor
  integer,intent(out):: index_crecorpo(:) ! index of polymer type for createcor
  integer,intent(out):: ncrecorw        ! water for createcor
  integer,intent(out):: index_crecorw   ! index of water type for createcor
  integer,intent(out):: ncrecorma       ! max number of matom type for createcor
  integer,intent(out):: index_crecorma(:) ! index of matom type for createcor

  integer,intent(in):: maxnstage        ! maximum number of MD stage
  integer,intent(out):: maxnstep         ! maximum step of MD (0-maxnstep)
  integer,intent(out):: nstage          ! stage number of MD
  integer,intent(inout):: nstep_stage(:) ! time step of each stage

  integer,intent(out):: mdcont_stage(:) ! MD control parameter of each stage

  integer,intent(out):: nstep_maxwell   ! time step of maxwell distribution
  integer,intent(out):: nstep_expand    ! time step of cell expansion

  real(8),intent(out):: xcel             ! initial x0 cell length[m]
  real(8),intent(out):: ycel             ! initial y0 cell length[m]
  real(8),intent(out):: zcel             ! initial z0 cell length[m]

  real(8),intent(out):: yratio           ! y cell ratio of y to x
  real(8),intent(out):: zratio           ! z cell ratio of z to x

  real(8),intent(out):: r_expand         ! expansion ratio of cell (volume)

  logical,intent(out):: ifstarec        ! read old state
  logical,intent(out):: ifcreatecor     ! create new coordinate

  logical,intent(out):: ifrdaddtop      ! input additional topology information

  logical,intent(out):: ifcenterfix_all   ! center fix for all
  logical,intent(out):: ifcenterfix_poly  ! center fix for polymer1
  logical,intent(out):: ifcenterfix_water ! center fix for water
  logical,intent(out):: ifcenterfix_ma    ! center fix for monatomic mole.
  character(4),intent(out):: cenfix_free  ! COM not fixed in this direction

  logical,intent(out):: ifcenterfix_polytyp(:) ! center fix for each polymer
  logical,intent(out):: ifcenterfix_watertyp   ! center fix for each water
  logical,intent(out):: ifcenterfix_matyp(:)
                                          ! center fix for each monatomic mole.

  real(8),intent(out):: dt_long_cal      ! time step of long force [sec]

!!! From Ver.1.74, nstep_med is forced to 1. !!!
  integer,intent(out):: nstep_med       ! number of step for medium force
  integer,intent(out):: nstep_short     ! number of step for short force

!!! From Ver.1.74, do not use 2(mts_med). !!!
  integer,intent(out):: mts_bond        ! MTS flag for bond
  integer,intent(out):: mts_angl        ! MTS flag for angle
  integer,intent(out):: mts_anglub      ! MTS flag for Urey-Bradley angle
  integer,intent(out):: mts_tors        ! MTS flag for torsion
  integer,intent(out):: mts_torsrb      ! MTS flag for torsionrb
  integer,intent(out):: mts_torsim      ! MTS flag for torsionim
  integer,intent(out):: mts_vdw         ! MTS flag for vdw interaction
  integer,intent(out):: mts_ewr         ! MTS flag for ewald real(=vdw)
  integer,intent(out):: mts_ewk         ! MTS flag for ewald wave
  integer,intent(out):: mts_vdw14       ! MTS flag for 14vdw
  integer,intent(out):: mts_elc14       ! MTS flag for 14elc(=mts_vdw14)
  integer,intent(out):: mts_mor         ! MTS flag for Morse interaction
  integer,intent(out):: mts_sh          ! MTS flag for SH interaction
  integer,intent(out):: mts_rfh         ! MTS flag for RFH interaction
  integer,intent(out):: mts_dou         ! MTS flag for DOU interaction
  integer,intent(out):: mts_cnpvw       ! MTS flag for CNP_VW
  integer,intent(out):: mts_cstmnb      ! MTS flag for custom NB interaction
  integer,intent(out):: mts_posres      ! MTS flag for position restraint
  integer,intent(out):: mts_potbias     ! MTS flag for bias potential

!---- virtual time step for evaluating constraint force
  integer,intent(out):: nstep_vir
                          ! number of step for virtual time integration

!---- limit atomic motion to a specific distance
  logical,intent(out):: iflimitmove ! if limit atomic motion to a specific distance
  real(8),intent(out):: limitdist  ! maximum atomic displacement when doing
                                   ! time integration [m] (structure relaxation)

!---- for Nose-Hoover chain
  integer,intent(out):: mchain           ! Nose-Hoover chain number (>1 in NVT)
  real(8),intent(out):: tfreq            ! temperature frequency [1/s]
  real(8),intent(out):: text             ! external temp. [K] (Nose-Hoover chain)

!---- for Andersen (Hoover type) barostat
  real(8),intent(out):: vfreq            ! volume change frequency [1/s]
  real(8),intent(out):: pext             ! external pressure [Pa]
  logical,intent(out):: ifpatmcont       ! atomic pressure control
  logical,intent(out):: ifpmolcont       ! molecular pressure control
  character(5),intent(out):: pcont_axis  ! axis for pressure control (iso, aniso, etc.)

!---- for higher order Trotter expansion
  integer,intent(out):: next           ! iteration number of extended system
  integer,intent(out):: nyosh          ! expansion order of Yoshida-Suzuki method
                                ! choose 1, 3 or 5

!---- some MD flags
  logical,intent(out):: ifrattle        ! rattle flag

  logical,intent(out):: ifewald         ! ewald flag
  logical,intent(out):: ifspme          ! spme (Smooth Particle Mesh Ewald) flag
  logical,intent(out):: iffennell       ! Fennell flag

  logical,intent(out):: ifljari         ! arithmetic mean for LJ cross parameter
  logical,intent(out):: ifljgeo         ! geometric mean for LJ cross parameter

  logical,intent(out):: iflocalheat     ! local heating flag
  logical,intent(out):: ifregionheat    ! region temp. control flag
  logical,intent(out):: ifregionhf      ! region heat flux control flag
  logical,intent(out):: ifreglange      ! region Langevin thermostat flag
  logical,intent(out):: iftcratom
                          ! region temp. control based on atom or mole.
  logical,intent(out):: ifoutthc        ! flag for outputting thermal control file

  logical,intent(out):: iflocalfix      ! fix atoms flag

  logical,intent(out):: iflocalfixz     ! flag for fixing z coordinate of atoms
  logical,intent(out):: iflocalfixzg    ! flag for fixing z coordinate of COM of molecules

  logical,intent(out):: ifcnp           ! flag to control normal pressure

  logical,intent(out):: ifposres        ! position restraint flag
  logical,intent(out):: ifpotbias       ! bias potential flag
  logical,intent(out):: iflocalvel      ! flag to force local atomic velocity
  logical,intent(out):: ifstrmvel       ! flag to input and use streaming velocity

!!! if you use NPT dynamics, you must choose ifcalpremole or ifcalpreatom !!!
  logical,intent(out):: ifcalpremole    ! pressure calculation of molecule
  logical,intent(out):: ifcalpreatom    ! pressure calculation of atom
  logical,intent(out):: ifnetqcorrp     ! net charge correction for pressure

!---- pressure calculation of L-J long-range correction
  logical,intent(out):: ifcalljlong       ! long-range correction in pressure
  character(2),intent(out):: solvetyp ! solvent molecule
  integer,intent(out):: nsolve            ! number of solvent molecules

!---- parameter for ewald method
  real(8),intent(out):: alpha            ! parameter alpha [1/m]
!      take care maxwave
  integer,intent(out):: kmax            ! parameter kmax
  real(8),intent(out):: rrcut            ! ewald real space cutoff length [m]

!---- parameter for SPME method
  integer,intent(out):: nfft1, nfft2, nfft3 ! grid points in SPME
  integer,intent(out):: pme_order       ! B-spline order

!---- other MD parameters
  logical,intent(out):: ifcellindex     ! flag for cell index

  real(8),intent(out):: rcut             ! vdw cutoff length [m]

  logical,intent(out):: ifbook          ! flag for bookkeeping
  real(8),intent(out):: rcut_book        ! cut off radius of bookkeeping[m]
  integer,intent(out):: nstep_book      ! bookkeeping interval

  real(8),intent(out):: tcont_poly(:)    ! poly Temp. [K] in NVT (Woodcock)
  real(8),intent(out):: tcont_water      ! H2O Temp. [K] in NVT (Woodcock)
  real(8),intent(out):: tcont_ma(:)
                          ! monatomic mole. Temp. [K] in NVT (Woodcock)
  real(8),intent(out):: tcont_poly_ini
                          ! poly Temp. [K] in NVT (Woodcock) for md_h
  real(8),intent(out):: tcont_water_ini
                          ! H2O Temp. [K] in NVT (Woodcock) for md_h
  real(8),intent(out):: tcont_ma_ini    ! MA Temp. [K] in NVT (Woodcock) for md_h

  integer,intent(out):: tcontinterval   ! interval of temp. control
  integer,intent(out):: outinterval     ! interval of outputting data
  integer,intent(out):: pressinterval   ! interval of pressure output
  integer,intent(out):: heatfinterval   ! interval of heatf output
  integer,intent(out):: recinterval     ! state record interval

  character(2),intent(out):: oatmtyp    ! O atomtype of water model
  character(2),intent(out):: hatmtyp    ! H atomtype of water model
  character(2),intent(out):: monoatmtyp(:) ! monatomic mole. type

  integer,intent(out):: randseed        ! random seed for createcor etc.

  real(8),intent(out):: compfact     ! compact factor using at poly arrange(<1.0)

  real(8),intent(out):: eps_rattle       ! tolerance (relative difference)
                                ! for bond length constraint by RATTLE

  real(8),intent(out):: rcutmor          ! Morse cutoff length [m]

  logical,intent(out):: ifcellindex_mor ! flag for cell index (morse)

  logical,intent(out):: ifbookmor     ! flag for bookkeeping of Morse interaction
  real(8),intent(out):: rcut_bookmor  ! cut off radius of bookkeeping[m] of Morse
  integer,intent(out):: nstep_bookmor ! bookkeeping interval of Morse interaction

  real(8),intent(out):: rcutsh           ! SH cutoff length [m]

  logical,intent(out):: ifcellindex_sh  ! flag for cell index (SH)

  logical,intent(out):: ifbooksh        ! flag for bookkeeping of SH interaction
  real(8),intent(out):: rcut_booksh      ! cut off radius of bookkeeping[m] of SH
  integer,intent(out):: nstep_booksh    ! bookkeeping interval of SH interaction

  real(8),intent(out):: rcutrfhfo        ! RFH(FO) cutoff length [m]

  logical,intent(out):: ifcellindex_rfhfo ! flag for cell index (RFH(FO))

  logical,intent(out):: ifbookrfhfo
                          ! flag for bookkeeping of RFH(FO) interaction
  real(8),intent(out):: rcut_bookrfhfo
                          ! cut off radius of bookkeeping[m] of RFH(FO)
  integer,intent(out):: nstep_bookrfhfo
                          ! bookkeeping interval of RFH(FO) interaction

  real(8),intent(out):: rcutrfhoo        ! RFH(OO) cutoff length [m]

  logical,intent(out):: ifcellindex_rfhoo ! flag for cell index (RFH(OO))

  logical,intent(out):: ifbookrfhoo ! flag for bookkeeping of RFH(OO) interaction
  real(8),intent(out):: rcut_bookrfhoo
                          ! cut off radius of bookkeeping[m] of RFH(OO)
  integer,intent(out):: nstep_bookrfhoo
                          ! bookkeeping interval of RFH(OO) interaction

  real(8),intent(out):: rcutrfhoh        ! RFH(OH) cutoff length [m]

  logical,intent(out):: ifcellindex_rfhoh ! flag for cell index (RFH(OH))

  logical,intent(out):: ifbookrfhoh ! flag for bookkeeping of RFH(OH) interaction
  real(8),intent(out):: rcut_bookrfhoh
                          ! cut off radius of bookkeeping[m] of RFH(OH)
  integer,intent(out):: nstep_bookrfhoh
                          ! bookkeeping interval of RFH(OH) interaction

  real(8),intent(out):: rcutdouo         ! DOU cutoff length [m] for O-Au
  real(8),intent(out):: rcutindouo       ! DOU cutin length [m] for O-Au

  logical,intent(out):: ifcellindex_douo ! flag for cell index (DOU) for O-Au

  logical,intent(out):: ifbookdouo
                          ! flag for bookkeeping of DOU interaction (O-Au)
  real(8),intent(out):: rcut_bookdouo
                          ! cut off radius of bookkeeping[m] of DOU (O-Au)
  integer,intent(out):: nstep_bookdouo
                          ! bookkeeping interval of DOU interaction (O-Au)

  real(8),intent(out):: rcutdouh         ! DOU cutoff length [m] for H-Au

  logical,intent(out):: ifcellindex_douh ! flag for cell index (DOU) for H-Au

  logical,intent(out):: ifbookdouh
                          ! flag for bookkeeping of DOU interaction (H-Au)
  real(8),intent(out):: rcut_bookdouh
                          ! cut off radius of bookkeeping[m] of DOU (H-Au)
  integer,intent(out):: nstep_bookdouh
                          ! bookkeeping interval of DOU interaction (H-Au)

  real(8),intent(out):: rcutrpvw ! RP-VW cutoff length
  logical,intent(out):: ifbookrpvw ! flag for bookkeeping of RP-VW interaction
  real(8),intent(out):: rcut_bookrpvw
                          ! cut off radius of bookkeep of RP-VW interaction
  integer,intent(out):: nstep_bookrpvw
                          ! bookkeeping interval of RP-VW interaction

!---- parameters for custom interaction ----!
  !!! Other parameters can be found in each custom interaction module program.

  logical,intent(out):: ifcstmnb           ! flag if using custom NB interaction
  logical,intent(out):: ifcellindex_cstmnb ! flag for cell index (custom NB)
  logical,intent(out):: ifbookcstmnb
                               ! flag for bookkeeping of custom NB interaction

!---- set charge from cor file
  logical,intent(out):: ifsetchrg(:)    ! if set charge from cor file

!---- PDB output
  logical,intent(out):: ifoutpdb        ! flag for outputting PDB format file
  integer,intent(out):: nstep_pdbout    ! MD step for output of PDB file

!---- parameters for energy minimization
  real(8),intent(out):: d_rini          ! initial displacement dr for EM [m]
  real(8),intent(out):: d_rmax          ! maximum displacement dr for EM [m]

  real(8),intent(out):: d_econv         ! convergence condition for energy in EM [J]
  real(8),intent(out):: d_rmsf          ! convergence condition for
                                        !    root mean square force in EM [N]

!---- spline interpolation for interaction sum
  integer,intent(out):: nspltbl         ! number of spline interpolation points

! LOCAL:
  character(80):: fredat(maxnword)

  integer:: i               ! do loop index

  integer:: iuini
  integer:: ios

  integer:: ipoly
  integer:: imatom
  integer:: imaxpo
  integer:: imaxw
  integer:: imaxma

  integer:: istage

  character(80):: mdconttyp
  integer:: t_maxnstep

  character(80):: mts_typ

  integer:: nword

  logical:: ifcellength
  logical:: ifcelratio

  logical:: ifcrecorpo(maxnpolytyp)
  logical:: ifcrecorw = .false.
  logical:: ifcrecorma(maxnmatyp)

  logical:: ifmuststarec = .false.

  real(8):: tcont_poly_tmp
  real(8):: tcont_ma_tmp

  integer:: cellcount(3)

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- intialize variable --------

  nstep_stage(1:maxnstage) = 0
  mdcont_stage(1:maxnstage) = 0

  nmatomtyp(1:maxnmatyp) = 0

  ifcellength = .false.
  ifcelratio = .false.

  ncrecorpo = 0
  ncrecorw  = 0
  ncrecorma = 0

  npolytyp = 0
  nwater = 0
  nmatyp = 0

  ifsetchrg(1:maxnpolytyp) = .false.

  ifcenterfix_polytyp(1:maxnpolytyp) = .false.
  ifcenterfix_watertyp = .false.
  ifcenterfix_matyp(1:maxnmatyp) = .false.

  nstep_med = 1

!---- open MD initial script (peachgk.ini)

  iuini = 99
  open(iuini, file='peachgk.ini', status='old', iostat=ios)
  if (ios /= 0) then
     write(6,*) 'Failure in opening file: peachgk.ini'
     stop
  end if

!---- read first line and check peachgk.ini version
  call rdfree_w(iuini, maxnword, fredat, nword)
  if ((fredat(7)(5:7) /= pgkiniver(1:3)) &
       & .or. (fredat(8)(2:9) /= pgkinidate(1:8))) then
     write(6,*) 'Error: Use the correct version of peachgk.ini'
     stop
  end if

!---- read all of the entries in peachgk.ini
  DO
     call rdfree_w(iuini, maxnword, fredat, nword)

     if (fredat(1)(1:3) == 'END') then
        exit

     else if (fredat(1) == ' ') then
        continue

!     input & output file name
     else if (fredat(1) == 'iuwtopname') then
        iuwtopname = fredat(2)

     else if (fredat(1) == 'iuparavdwname') then
        iuparavdwname = fredat(2)

     else if (fredat(1) == 'iuparabondname') then
        iuparabondname = fredat(2)

     else if (fredat(1) == 'iuparaconstname') then
        iuparaconstname = fredat(2)

     else if (fredat(1) == 'iuparacstmnbname') then
        iuparacstmnbname = fredat(2)

     else if (fredat(1) == 'iuaddtopname') then
        iuaddtopname = fredat(2)

     else if (fredat(1) == 'iustrmvelname') then
        iustrmvelname = fredat(2)

     else if (fredat(1) == 'iuposresname') then
        iuposresname = fredat(2)

     else if (fredat(1) == 'iostarecname') then
        iostarecname = fredat(2)

     else if (fredat(1) == 'ousumname') then
        ousumname = fredat(2)

     else if (fredat(1) == 'ouenename') then
        ouenename = fredat(2)
        read(fredat(3),*) ifoutene

     else if (fredat(1) == 'ouposname') then
        ouposname = fredat(2)
        read(fredat(3),*) ifoutpos

     else if (fredat(1) == 'ouvelname') then
        ouvelname = fredat(2)
        read(fredat(3),*) ifoutvel

     else if (fredat(1) == 'ouforname') then
        ouforname = fredat(2)
        read(fredat(3),*) ifoutfor

     else if (fredat(1) == 'outhename') then
        outhename = fredat(2)
        read(fredat(3),*) ifoutthe

     else if (fredat(1) == 'oubarname') then
        oubarname = fredat(2)
        read(fredat(3),*) ifoutbar

     else if (fredat(1) == 'ouprename') then
        ouprename = fredat(2)
        read(fredat(3),*) ifoutpre

     else if (fredat(1) == 'outhcname') then
        outhcname = fredat(2)

     else if (fredat(1) == 'oupdbname') then
        oupdbname = fredat(2)

!     number of particle
     else if (fredat(1) == 'npolytyp') then
        read(fredat(2),*) npolytyp

     else if (fredat(1) == 'npolymoletyp') then
        read(fredat(2),*) ipoly
        if (ipoly > npolytyp) then
           write(6,*) 'Error: npolymoletyp in rdscript'
           stop
        end if
        read(fredat(3),*) npoly_mole(ipoly)
        read(fredat(4),*) npoly_atom(ipoly)
        iucorname(ipoly) = fredat(5)
        iutopname(ipoly) = fredat(6)
        read(fredat(7),*) ifcrecorpo(ipoly)
        do i=8,nword
           polytyp_free(ipoly,i-7) = fredat(i)
        end do

!       - detect setcharge & centerfix
        do i=8,nword
           if (polytyp_free(ipoly,i-7) == 'setcharge') then
              ifsetchrg(ipoly) = .true.
              cycle
           end if

           if (polytyp_free(ipoly,i-7) == 'centerfix') then
              ifcenterfix_polytyp(ipoly) = .true.
              cycle
           end if
        end do

        if (.not. ifcrecorpo(ipoly)) then
           ifmuststarec = .true.
        end if

     else if (fredat(1) == 'nwater') then
        read(fredat(2),*) nwater
        read(fredat(3),*) ifcrecorw
        do i=4,nword
           watertyp_free(i-3) = fredat(i)
        end do

!       - detect centerfix
        do i=4,nword
           if (watertyp_free(i-3) == 'centerfix') then
              ifcenterfix_watertyp = .true.
              cycle
           end if
        end do

        if (.not. ifcrecorw .and. nwater /= 0) then
           ifmuststarec = .true.
        end if

     else if (fredat(1) == 'nmatyp') then
        read(fredat(2),*) nmatyp

     else if (fredat(1) == 'nmatomtyp') then
        read(fredat(2),*) imatom
        if (imatom > nmatyp) then
           write(6,*) 'Error: nmatomtyp in rdscript'
           stop
        end if
        read(fredat(3),*) nmatomtyp(imatom)
        monoatmtyp(imatom) = fredat(4)
        read(fredat(5),*) ifcrecorma(imatom)
        do i=6,nword
           matomtyp_free(imatom,i-5) = fredat(i)
        end do

!       - detect centerfix
        do i=6,nword
           if (matomtyp_free(imatom,i-5) == 'centerfix') then
              ifcenterfix_matyp(imatom) = .true.
              cycle
           end if
        end do

        if (.not. ifcrecorma(imatom)) then
           ifmuststarec = .true.
        end if

     else if (fredat(1) == 'maxpo') then
        read(fredat(2),*) imaxpo
!            if (.not. ifcrecorpo(imaxpo)) then
!               write(6,*) 'Error: ifcrecorpo(',imaxpo,') is .false.,'
!               write(6,*) '       but maxpo is defined.'
!               stop
!            end if
        if (ifcrecorpo(imaxpo)) then
           ncrecorpo = ncrecorpo + 1
           index_crecorpo(ncrecorpo) = imaxpo
        end if

        if (imaxpo > npolytyp) then
           write(6,*) 'Error: maxpo in rdscript'
           stop
        end if
        read(fredat(3),*) xmaxpo(imaxpo)
        read(fredat(4),*) ymaxpo(imaxpo)
        read(fredat(5),*) zmaxpo(imaxpo)
        do i=6,nword
           read(fredat(i),*) inicorpo(imaxpo,i-5)
        end do

     else if (fredat(1) == 'maxw') then
!            if (.not. ifcrecorw) then
!               write(6,*) 'Error: ifcrecorw is .false.,'
!               write(6,*) '       but maxw is defined.'
!               stop
!            end if
        read(fredat(2),*) imaxw
        if (ifcrecorw) then
           ncrecorw = 1
           index_crecorw = 1
        end if

        if (imaxw /= 1) then
           write(6,*) 'Error: maxw in rdscript'
           stop
        end if
        read(fredat(3),*) xmaxw
        read(fredat(4),*) ymaxw
        read(fredat(5),*) zmaxw
        do i=6,nword
           read(fredat(i),*) inicorw(i-5)
        end do

     else if (fredat(1) == 'maxma') then
        read(fredat(2),*) imaxma
!            if (.not. ifcrecorma(imaxma)) then
!               write(6,*) 'Error: ifcrecorma(',imaxma,') is .false.,'
!               write(6,*) '       but maxma is defined.'
!               stop
!            end if
        if (ifcrecorma(imaxma)) then
           ncrecorma = ncrecorma + 1
           index_crecorma(ncrecorma) = imaxma
        end if

        if (imaxma > nmatyp) then
           write(6,*) 'Error: maxma in rdscript'
           stop
        end if
        read(fredat(3),*) xmaxma(imaxma)
        read(fredat(4),*) ymaxma(imaxma)
        read(fredat(5),*) zmaxma(imaxma)
        do i=6,nword
           read(fredat(i),*) inicorma(imaxma,i-5)
        end do

!     parameter of MD stages
     else if (fredat(1) == 'maxnstep') then
        read(fredat(2),*) maxnstep

     else if (fredat(1) == 'nstage') then
        read(fredat(2),*) nstage

     else if (fredat(1) == 'nstep_stage') then
        read(fredat(2),*) istage
        if (istage > nstage) then
           write(6,*) 'Error: nstep_stage in rdscript'
           stop
        end if
        read(fredat(3),*) nstep_stage(istage)

     else if (fredat(1) == 'mdcont_stage') then
        read(fredat(2),*) istage
        if (istage > nstage) then
           write(6,*) 'Error: mdcont_stage in rdscript'
           stop
        end if
        mdconttyp = fredat(3)
        if (mdconttyp == 'md_0k') then
           mdcont_stage(istage) = md_0k
        else if (mdconttyp == 'md_h') then
           mdcont_stage(istage) = md_h
        else if (mdconttyp == 'md_t') then
           mdcont_stage(istage) = md_t
        else if (mdconttyp == 'md_mtk') then
           mdcont_stage(istage) = md_mtk
        else if (mdconttyp == 'md_nhc') then
           mdcont_stage(istage) = md_nhc
        else if (mdconttyp == 'md_nve') then
           mdcont_stage(istage) = md_nve
        else if (mdconttyp == 'md_htf') then
           mdcont_stage(istage) = md_htf
        else if (mdconttyp == 'md_ems') then
           mdcont_stage(istage) = md_ems
        end if

     else if (fredat(1) == 'nstep_maxwell') then
        read(fredat(2),*) nstep_maxwell

     else if (fredat(1) == 'nstep_expand') then
        read(fredat(2),*) nstep_expand

!---- cell dimensions
     else if (fredat(1) == 'xcel') then
        read(fredat(2),*) xcel

     else if (fredat(1) == 'yratio') then
        ifcelratio = .true.
        read(fredat(2),*) yratio

     else if (fredat(1) == 'zratio') then
        ifcelratio = .true.
        read(fredat(2),*) zratio

     else if (fredat(1) == 'ycel') then
        ifcellength = .true.
        read(fredat(2),*) ycel

     else if (fredat(1) == 'zcel') then
        ifcellength = .true.
        read(fredat(2),*) zcel

     else if (fredat(1) == 'r_expand') then
        read(fredat(2),*) r_expand

!---- some important parameters
     else if (fredat(1) == 'ifstarec') then
        read(fredat(2),*) ifstarec

     else if (fredat(1) == 'ifcreatecor') then
        read(fredat(2),*) ifcreatecor

     else if (fredat(1) == 'ifrdaddtop') then
        read(fredat(2),*) ifrdaddtop

     else if (fredat(1) == 'ifcenterfix_all') then
        read(fredat(2),*) ifcenterfix_all

     else if (fredat(1) == 'cenfix_free') then
        cenfix_free = fredat(2)

     else if (fredat(1) == 'ifcenterfix_poly') then
        read(fredat(2),*) ifcenterfix_poly

     else if (fredat(1) == 'ifcenterfix_water') then
        read(fredat(2),*) ifcenterfix_water

     else if (fredat(1) == 'ifcenterfix_ma') then
        read(fredat(2),*) ifcenterfix_ma

!---- PDB output
     else if (fredat(1) == 'ifoutpdb') then
        read(fredat(2),*) ifoutpdb

     else if (fredat(1) == 'nstep_pdbout') then
        read(fredat(2),*) nstep_pdbout

!---- time step and MTS parameters
     else if (fredat(1) == 'dt_long_cal') then
        read(fredat(2),*) dt_long_cal

     else if (fredat(1) == 'nstep_short') then
        read(fredat(2),*) nstep_short

     else if (fredat(1) == 'mts_bond') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_bond = mts_long
        else if (mts_typ == 'med') then
           mts_bond = mts_med
        else if (mts_typ == 'short') then
           mts_bond = mts_short
        end if

     else if (fredat(1) == 'mts_angl') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_angl = mts_long
        else if (mts_typ == 'med') then
           mts_angl = mts_med
        else if (mts_typ == 'short') then
           mts_angl = mts_short
        end if

     else if (fredat(1) == 'mts_anglub') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_anglub = mts_long
        else if (mts_typ == 'med') then
           mts_anglub = mts_med
        else if (mts_typ == 'short') then
           mts_anglub = mts_short
        end if

     else if (fredat(1) == 'mts_tors') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_tors = mts_long
        else if (mts_typ == 'med') then
           mts_tors = mts_med
        else if (mts_typ == 'short') then
           mts_tors = mts_short
        end if

     else if (fredat(1) == 'mts_torsrb') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_torsrb = mts_long
        else if (mts_typ == 'med') then
           mts_torsrb = mts_med
        else if (mts_typ == 'short') then
           mts_torsrb = mts_short
        end if

     else if (fredat(1) == 'mts_torsim') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_torsim = mts_long
        else if (mts_typ == 'med') then
           mts_torsim = mts_med
        else if (mts_typ == 'short') then
           mts_torsim = mts_short
        end if

     else if (fredat(1) == 'mts_vdw') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_vdw = mts_long
        else if (mts_typ == 'med') then
           mts_vdw = mts_med
        else if (mts_typ == 'short') then
           mts_vdw = mts_short
        end if

     else if (fredat(1) == 'mts_ewr') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_ewr = mts_long
        else if (mts_typ == 'med') then
           mts_ewr = mts_med
        else if (mts_typ == 'short') then
           mts_ewr = mts_short
        end if

     else if (fredat(1) == 'mts_ewk') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_ewk = mts_long
        else if (mts_typ == 'med') then
           mts_ewk = mts_med
        else if (mts_typ == 'short') then
           mts_ewk = mts_short
        end if

     else if (fredat(1) == 'mts_vdw14') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_vdw14 = mts_long
        else if (mts_typ == 'med') then
           mts_vdw14 = mts_med
        else if (mts_typ == 'short') then
           mts_vdw14 = mts_short
        end if

     else if (fredat(1) == 'mts_elc14') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_elc14 = mts_long
        else if (mts_typ == 'med') then
           mts_elc14 = mts_med
        else if (mts_typ == 'short') then
           mts_elc14 = mts_short
        end if

     else if (fredat(1) == 'mts_mor') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_mor = mts_long
        else if (mts_typ == 'med') then
           mts_mor = mts_med
        else if (mts_typ == 'short') then
           mts_mor = mts_short
        end if

     else if (fredat(1) == 'mts_sh') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_sh = mts_long
        else if (mts_typ == 'med') then
           mts_sh = mts_med
        else if (mts_typ == 'short') then
           mts_sh = mts_short
        end if

     else if (fredat(1) == 'mts_rfh') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_rfh = mts_long
        else if (mts_typ == 'med') then
           mts_rfh = mts_med
        else if (mts_typ == 'short') then
           mts_rfh = mts_short
        end if

     else if (fredat(1) == 'mts_dou') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_dou = mts_long
        else if (mts_typ == 'med') then
           mts_dou = mts_med
        else if (mts_typ == 'short') then
           mts_dou = mts_short
        end if

     else if (fredat(1) == 'mts_cnpvw') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_cnpvw = mts_long
        else if (mts_typ == 'med') then
           mts_cnpvw = mts_med
        else if (mts_typ == 'short') then
           mts_cnpvw = mts_short
        end if

     else if (fredat(1) == 'mts_cstmnb') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_cstmnb = mts_long
        else if (mts_typ == 'med') then
           mts_cstmnb = mts_med
        else if (mts_typ == 'short') then
           mts_cstmnb = mts_short
        end if

     else if (fredat(1) == 'mts_posres') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_posres = mts_long
        else if (mts_typ == 'med') then
           mts_posres = mts_med
        else if (mts_typ == 'short') then
           mts_posres = mts_short
        end if

     else if (fredat(1) == 'mts_potbias') then
        mts_typ = fredat(2)
        if (mts_typ == 'long') then
           mts_potbias = mts_long
        else if (mts_typ == 'med') then
           mts_potbias = mts_med
        else if (mts_typ == 'short') then
           mts_potbias = mts_short
        end if

     else if (fredat(1) == 'nstep_vir') then
        read(fredat(2),*) nstep_vir

     else if (fredat(1) == 'iflimitmove') then
        read(fredat(2),*) iflimitmove

     else if (fredat(1) == 'limitdist') then
        read(fredat(2),*) limitdist

!---- Nose-Hoover chain and MTK eq. and higher order integration
     else if (fredat(1) == 'mchain') then
        read(fredat(2),*) mchain

     else if (fredat(1) == 'tfreq') then
        read(fredat(2),*) tfreq

     else if (fredat(1) == 'text') then
        read(fredat(2),*) text

     else if (fredat(1) == 'vfreq') then
        read(fredat(2),*) vfreq

     else if (fredat(1) == 'pext') then
        read(fredat(2),*) pext

     else if (fredat(1) == 'ifpatmcont') then
        read(fredat(2),*) ifpatmcont

     else if (fredat(1) == 'ifpmolcont') then
        read(fredat(2),*) ifpmolcont

     else if (fredat(1) == 'pcont_axis') then
        pcont_axis = fredat(2)

     else if (fredat(1) == 'next') then
        read(fredat(2),*) next

     else if (fredat(1) == 'nyosh') then
        read(fredat(2),*) nyosh

!---- some MD technics
     else if (fredat(1) == 'ifrattle') then
        read(fredat(2),*) ifrattle

     else if (fredat(1) == 'ifewald') then
        read(fredat(2),*) ifewald

     else if (fredat(1) == 'ifspme') then
        read(fredat(2),*) ifspme

     else if (fredat(1) == 'iffennell') then
        read(fredat(2),*) iffennell

     else if (fredat(1) == 'ifljari') then
        read(fredat(2),*) ifljari

     else if (fredat(1) == 'ifljgeo') then
        read(fredat(2),*) ifljgeo

     else if (fredat(1) == 'iflocalheat') then
        read(fredat(2),*) iflocalheat

     else if (fredat(1) == 'ifregionheat') then
        read(fredat(2),*) ifregionheat

     else if (fredat(1) == 'ifregionhf') then
        read(fredat(2),*) ifregionhf

     else if (fredat(1) == 'ifreglange') then
        read(fredat(2),*) ifreglange

     else if (fredat(1) == 'iftcratom') then
        read(fredat(2),*) iftcratom

     else if (fredat(1) == 'ifoutthc') then
        read(fredat(2),*) ifoutthc

     else if (fredat(1) == 'iflocalfix') then
        read(fredat(2),*) iflocalfix

     else if (fredat(1) == 'iflocalfixz') then
        read(fredat(2),*) iflocalfixz

     else if (fredat(1) == 'iflocalfixzg') then
        read(fredat(2),*) iflocalfixzg

     else if (fredat(1) == 'ifcnp') then
        read(fredat(2),*) ifcnp

     else if (fredat(1) == 'ifposres') then
        read(fredat(2),*) ifposres

     else if (fredat(1) == 'ifpotbias') then
        read(fredat(2),*) ifpotbias

     else if (fredat(1) == 'iflocalvel') then
        read(fredat(2),*) iflocalvel

     else if (fredat(1) == 'ifstrmvel') then
        read(fredat(2),*) ifstrmvel

     else if (fredat(1) == 'ifcalpremole') then
        read(fredat(2),*) ifcalpremole

     else if (fredat(1) == 'ifcalpreatom') then
        read(fredat(2),*) ifcalpreatom

     else if (fredat(1) == 'ifnetqcorrp') then
        read(fredat(2),*) ifnetqcorrp

!---- pressure calculation of L-J long-range correction
     else if (fredat(1) == 'ifcalljlong') then
        read(fredat(2),*) ifcalljlong

     else if (fredat(1) == 'solvetyp') then
        solvetyp = fredat(2)

     else if (fredat(1) == 'nsolve') then
        read(fredat(2),*) nsolve

!---- parameter for ewald method
     else if (fredat(1) == 'alpha') then
        read(fredat(2),*) alpha

     else if (fredat(1) == 'kmax') then
        read(fredat(2),*) kmax

     else if (fredat(1) == 'rrcut') then
        read(fredat(2),*) rrcut

!---- parameter for SPME method
     else if (fredat(1) == 'nfft1') then
        read(fredat(2),*) nfft1

     else if (fredat(1) == 'nfft2') then
        read(fredat(2),*) nfft2

     else if (fredat(1) == 'nfft3') then
        read(fredat(2),*) nfft3

     else if (fredat(1) == 'pme_order') then
        read(fredat(2),*) pme_order

!---- parameter for energy minimization
     else if (fredat(1) == 'd_rini') then
        read(fredat(2),*) d_rini

     else if (fredat(1) == 'd_rmax') then
        read(fredat(2),*) d_rmax

     else if (fredat(1) == 'd_econv') then
        read(fredat(2),*) d_econv

     else if (fredat(1) == 'd_rmsf') then
        read(fredat(2),*) d_rmsf

!---- other MD parameters
     else if (fredat(1) == 'ifcellindex') then
        read(fredat(2),*) ifcellindex

     else if (fredat(1) == 'rcut') then
        read(fredat(2),*) rcut

     else if (fredat(1) == 'ifbook') then
        read(fredat(2),*) ifbook

     else if (fredat(1) == 'rcut_book') then
        read(fredat(2),*) rcut_book

     else if (fredat(1) == 'nstep_book') then
        read(fredat(2),*) nstep_book

     else if (fredat(1) == 'tcont_poly') then
        read(fredat(2),*) tcont_poly_tmp

     else if (fredat(1) == 'tcont_water') then
        read(fredat(2),*) tcont_water

     else if (fredat(1) == 'tcont_ma') then
        read(fredat(2),*) tcont_ma_tmp

     else if (fredat(1) == 'tcont_poly_ini') then
        read(fredat(2),*) tcont_poly_ini

     else if (fredat(1) == 'tcont_water_ini') then
        read(fredat(2),*) tcont_water_ini

     else if (fredat(1) == 'tcont_ma_ini') then
        read(fredat(2),*) tcont_ma_ini

     else if (fredat(1) == 'tcontinterval') then
        read(fredat(2),*) tcontinterval

     else if (fredat(1) == 'outinterval') then
        read(fredat(2),*) outinterval

     else if (fredat(1) == 'pressinterval') then
        read(fredat(2),*) pressinterval

     else if (fredat(1) == 'heatfinterval') then
        read(fredat(2),*) heatfinterval

     else if (fredat(1) == 'recinterval') then
        read(fredat(2),*) recinterval

     else if (fredat(1) == 'oatmtyp') then
        oatmtyp = fredat(2)

     else if (fredat(1) == 'hatmtyp') then
        hatmtyp = fredat(2)

     else if (fredat(1) == 'randseed') then
        read(fredat(2),*) randseed

     else if (fredat(1) == 'compfact') then
        read(fredat(2),*) compfact

     else if (fredat(1) == 'eps_rattle') then
        read(fredat(2),*) eps_rattle

     else if (fredat(1) == 'rcutmor') then
        read(fredat(2),*) rcutmor

     else if (fredat(1) == 'ifcellindex_mor') then
        read(fredat(2),*) ifcellindex_mor

     else if (fredat(1) == 'ifbookmor') then
        read(fredat(2),*) ifbookmor

     else if (fredat(1) == 'rcut_bookmor') then
        read(fredat(2),*) rcut_bookmor

     else if (fredat(1) == 'nstep_bookmor') then
        read(fredat(2),*) nstep_bookmor

     else if (fredat(1) == 'rcutsh') then
        read(fredat(2),*) rcutsh

     else if (fredat(1) == 'ifcellindex_sh') then
        read(fredat(2),*) ifcellindex_sh

     else if (fredat(1) == 'ifbooksh') then
        read(fredat(2),*) ifbooksh

     else if (fredat(1) == 'rcut_booksh') then
        read(fredat(2),*) rcut_booksh

     else if (fredat(1) == 'nstep_booksh') then
        read(fredat(2),*) nstep_booksh

     else if (fredat(1) == 'rcutrfhfo') then
        read(fredat(2),*) rcutrfhfo

     else if (fredat(1) == 'ifcellindex_rfhfo') then
        read(fredat(2),*) ifcellindex_rfhfo

     else if (fredat(1) == 'ifbookrfhfo') then
        read(fredat(2),*) ifbookrfhfo

     else if (fredat(1) == 'rcut_bookrfhfo') then
        read(fredat(2),*) rcut_bookrfhfo

     else if (fredat(1) == 'nstep_bookrfhfo') then
        read(fredat(2),*) nstep_bookrfhfo

     else if (fredat(1) == 'rcutrfhoo') then
        read(fredat(2),*) rcutrfhoo

     else if (fredat(1) == 'ifcellindex_rfhoo') then
        read(fredat(2),*) ifcellindex_rfhoo

     else if (fredat(1) == 'ifbookrfhoo') then
        read(fredat(2),*) ifbookrfhoo

     else if (fredat(1) == 'rcut_bookrfhoo') then
        read(fredat(2),*) rcut_bookrfhoo

     else if (fredat(1) == 'nstep_bookrfhoo') then
        read(fredat(2),*) nstep_bookrfhoo

     else if (fredat(1) == 'rcutrfhoh') then
        read(fredat(2),*) rcutrfhoh

     else if (fredat(1) == 'ifcellindex_rfhoh') then
        read(fredat(2),*) ifcellindex_rfhoh

     else if (fredat(1) == 'ifbookrfhoh') then
        read(fredat(2),*) ifbookrfhoh

     else if (fredat(1) == 'rcut_bookrfhoh') then
        read(fredat(2),*) rcut_bookrfhoh

     else if (fredat(1) == 'nstep_bookrfhoh') then
        read(fredat(2),*) nstep_bookrfhoh

     else if (fredat(1) == 'rcutdouo') then
        read(fredat(2),*) rcutdouo

     else if (fredat(1) == 'rcutindouo') then
        read(fredat(2),*) rcutindouo

     else if (fredat(1) == 'ifcellindex_douo') then
        read(fredat(2),*) ifcellindex_douo

     else if (fredat(1) == 'ifbookdouo') then
        read(fredat(2),*) ifbookdouo

     else if (fredat(1) == 'rcut_bookdouo') then
        read(fredat(2),*) rcut_bookdouo

     else if (fredat(1) == 'nstep_bookdouo') then
        read(fredat(2),*) nstep_bookdouo

     else if (fredat(1) == 'rcutdouh') then
        read(fredat(2),*) rcutdouh

     else if (fredat(1) == 'ifcellindex_douh') then
        read(fredat(2),*) ifcellindex_douh

     else if (fredat(1) == 'ifbookdouh') then
        read(fredat(2),*) ifbookdouh

     else if (fredat(1) == 'rcut_bookdouh') then
        read(fredat(2),*) rcut_bookdouh

     else if (fredat(1) == 'nstep_bookdouh') then
        read(fredat(2),*) nstep_bookdouh

     else if (fredat(1) == 'rcutrpvw') then
        read(fredat(2),*) rcutrpvw

     else if (fredat(1) == 'ifbookrpvw') then
        read(fredat(2),*) ifbookrpvw

     else if (fredat(1) == 'rcut_bookrpvw') then
        read(fredat(2),*) rcut_bookrpvw

     else if (fredat(1) == 'nstep_bookrpvw') then
        read(fredat(2),*) nstep_bookrpvw

     else if (fredat(1) == 'ifcstmnb') then
        read(fredat(2),*) ifcstmnb

     else if (fredat(1) == 'ifcellindex_cstmnb') then
        read(fredat(2),*) ifcellindex_cstmnb

     else if (fredat(1) == 'ifbookcstmnb') then
        read(fredat(2),*) ifbookcstmnb

     else if (fredat(1) == 'nspltbl') then
        read(fredat(2),*) nspltbl

     end if

  END DO

!---- check the parameter relations and so on
  t_maxnstep = 0
  do i=1,nstage
     t_maxnstep = t_maxnstep + nstep_stage(i)
  end do
  if (t_maxnstep /= maxnstep) then
     write(6,*) 'Error : something is wrong with nstep_stage in rdscript'
     stop
  end if

  npoly = 0
  do i=1,npolytyp
     npoly = npoly + npoly_mole(i)
  end do

  nmatom = 0
  do i=1,nmatyp
     nmatom = nmatom + nmatomtyp(i)
  end do

  if (.not. ifstarec .and. .not. ifcreatecor) then
     write(6,*) 'Error: choose at least ifstarec or ifcreatecor'
     stop
  end if

  if (.not. ifcelratio .and. .not. ifcellength) then
     write(6,*) 'Error: set the cell dimension'
     stop
  else if (ifcelratio .and. .not. ifcellength) then   ! cell ratio is set
     ycel = xcel * yratio
     zcel = xcel * zratio
  else   ! cell length is set
     yratio = ycel / xcel
     zratio = zcel / xcel
  end if

  if (ifmuststarec .and. .not. ifstarec) then
     write(6,*) 'Error: you must use .true. in ifstarec'
     stop
  end if

  do i=1,npolytyp
     if (ifcrecorpo(i) .and. &
          & (npoly_mole(i) /= xmaxpo(i)*ymaxpo(i)*zmaxpo(i))) then
        write(6,*) 'Error: discrepancy between maxpo', &
             &     ' and npoly_mole'
        stop
     end if
  end do

  if (ifcrecorw .and. (nwater /= xmaxw*ymaxw*zmaxw)) then
     write(6,*) 'Error: discrepancy between maxw', &
          &     ' and nwater'
     stop
  end if

  do i=1,nmatyp
     if (ifcrecorma(i) .and. &
          & (nmatomtyp(i) /= xmaxma(i)*ymaxma(i)*zmaxma(i))) then
        write(6,*) 'Error: discrepancy between maxma', &
             &     ' and nmatomtyp'
        stop
     end if
  end do

  if (ifpatmcont .and. ifpmolcont) then
     write(6,*) 'Error: atm and mol-control cannot be implemented together.'
     stop
  else if ((.not.ifpatmcont) .and. (.not.ifpmolcont)) then
     write(6,*) 'Error: use either atm or mol-control in NPT-MD.'
     stop
  end if

!      if (.not.ifcalpremole .and. .not.ifcalpreatom) then
!         write(6,*) 'Error : ifcalpreatom and ifcalpremole are ',
!     &                     ' false together in rdscript'
!         stop
!      end if

  if ((pcont_axis /= 'iso') .and. (pcont_axis /= 'aniso') .and. &
    & (pcont_axis /= 'x') .and. (pcont_axis /= 'y') .and. &
    & (pcont_axis /= 'z')) then
     write(6,*) 'Error: something is wrong with pcont_axis.'
     stop
  end if

  if ((ifewald .and. ifspme) .or. (ifspme .and. iffennell) .or. &
       & (iffennell .and. ifewald)) then
     write(6,*) 'Error : you must choose just one from', &
          &     ' ifewald, ifspme, and iffennell.'
     stop
  end if

#if !defined(_LJ_ONLY)
  if ((.not.ifewald) .and. (.not.ifspme) .and. (.not.iffennell)) then
     write(6,*) 'Error : you must choose at least one from', &
          &     ' ifewald, ifspme, and iffennell.'
     stop
  end if
#endif

  if (ifljari .and. ifljgeo) then
     write(6,*) 'Error : you cannot choose both', &
          &     ' ifljari and ifljgeo.'
     stop
  end if

  if (.not. ifljari .and. .not. ifljgeo) then
     write(6,*) 'Error : you must choose either', &
          &     ' ifljari or ifljgeo.'
     stop
  end if

#if defined(_NOT_FENNELL)
  if (rrcut > rcut) then
     write(6,*) 'Error : rrcut > rcut'
     stop
  end if
#endif

  if (ifcreatecor .and. ifposres) then
     write(6,*) 'Error : you cannot choose both', &
          &     ' ifcreatecor and ifposres'
     stop
  end if

!     - check cellindex method
  if (.not.ifbook .and. ifcellindex) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       choose ifbook'
     stop
  end if

  if (ifcellindex) then
     cellcount(1) = INT(xcel/rcut_book)
     cellcount(2) = INT(ycel/rcut_book)
     cellcount(3) = INT(zcel/rcut_book)

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (vdw & elc)'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex = .false.
     end if
  end if

  if (.not.ifbookmor .and. ifcellindex_mor) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       choose ifbookmor'
     stop
  end if

  if (ifcellindex_mor) then
     cellcount(1) = INT(xcel/rcut_bookmor)
     cellcount(2) = INT(ycel/rcut_bookmor)
     cellcount(3) = INT(zcel/rcut_bookmor)

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (morse)'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_mor = .false.
     end if
  end if

  if (.not.ifbooksh .and. ifcellindex_sh) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       choose ifbooksh'
     stop
  end if

  if (ifcellindex_sh) then
     cellcount(1) = INT(xcel/rcut_booksh)
     cellcount(2) = INT(ycel/rcut_booksh)
     cellcount(3) = INT(zcel/rcut_booksh)

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (SH)'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_sh = .false.
     end if
  end if

  if (.not.ifbookrfhfo .and. ifcellindex_rfhfo) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       choose ifbookrfhfo'
     stop
  end if

  if (ifcellindex_rfhfo) then
     cellcount(1) = INT(xcel/rcut_bookrfhfo)
     cellcount(2) = INT(ycel/rcut_bookrfhfo)
     cellcount(3) = INT(zcel/rcut_bookrfhfo)

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (RFH(F-O))'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_rfhfo = .false.
     end if
  end if

  if (.not.ifbookrfhoo .and. ifcellindex_rfhoo) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       choose ifbookrfhoo'
     stop
  end if

  if (ifcellindex_rfhoo) then
     cellcount(1) = INT(xcel/rcut_bookrfhoo)
     cellcount(2) = INT(ycel/rcut_bookrfhoo)
     cellcount(3) = INT(zcel/rcut_bookrfhoo)

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (RFH(O-O))'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_rfhoo = .false.
     end if
  end if

  if (.not.ifbookrfhoh .and. ifcellindex_rfhoh) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       choose ifbookrfhoh'
     stop
  end if

  if (ifcellindex_rfhoh) then
     cellcount(1) = INT(xcel/rcut_bookrfhoh)
     cellcount(2) = INT(ycel/rcut_bookrfhoh)
     cellcount(3) = INT(zcel/rcut_bookrfhoh)

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (RFH(O-H))'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_rfhoh = .false.
     end if
  end if

  if (.not.ifbookdouo .and. ifcellindex_douo) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       choose ifbookdouo'
     stop
  end if

  if (ifcellindex_douo) then
     cellcount(1) = INT(xcel/rcut_bookdouo)
     cellcount(2) = INT(ycel/rcut_bookdouo)
     cellcount(3) = INT(zcel/rcut_bookdouo)

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (DOU O-Au)'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_douo = .false.
     end if
  end if

  if (.not.ifbookdouh .and. ifcellindex_douh) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       choose ifbookdouh'
     stop
  end if

  if (ifcellindex_douh) then
     cellcount(1) = INT(xcel/rcut_bookdouh)
     cellcount(2) = INT(ycel/rcut_bookdouh)
     cellcount(3) = INT(zcel/rcut_bookdouh)

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (DOU H-Au)'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_douh = .false.
     end if
  end if

  if ((nfft1 > maxfft1) .or. (nfft2 > maxfft2) .or. (nfft3 > maxfft3)) then
     write(6,*) 'Error: nfft* exceeds maximum number of nfft* (MAXFFT*)'
     write(6,*) '       check config.h'
     stop
  end if

  if ((nfft1 < nproc) .or. (nfft2 < nproc) .or. (nfft3 < nproc)) then
     write(6,*) 'Error: nfft* must exceed number of procs in a MPI job.'
     stop
  end if

  if (pme_order > maxord) then
     write(6,*) 'Error: pme_order exceeds maximum pme_order (MAXORDER)'
     write(6,*) '       check config.h'
     stop
  end if

  if (iflocalvel .and. (ifcenterfix_all .or. ifcenterfix_poly &
       &           .or. ifcenterfix_water .or. ifcenterfix_ma)) then
     write(6,*) 'Error: do not use iflocalvel and ifcenterfix_* at once.'
     write(6,*) '       Turn off all of ifcenterfix flags.'
     stop
  end if

  if (cenfix_free /= 'none' .and. .not.ifcenterfix_all) then
     write(6,*) 'Error: use cenfix_free along with ifcenterfix_all.'
     stop
  end if

  if (ifstrmvel .and. (nstep_expand >= 0)) then
     write(6,*) 'Error: do not use ifstrmvel and cell expansion at once.'
     stop
  end if

!---- process for heat flux calculation

#if !defined(HF)
  do i= 1, nstage
     if (mdcont_stage(i) == md_htf) then
        write(6,*) 'Error : You cannot choose md_htf.'
        write(6,*) '        Compile with -DHF switch.'
        stop
     end if
  end do
#endif

!     copy temperature
  tcont_poly(1) = tcont_poly_tmp
  do ipoly = 2, npolytyp
     tcont_poly(ipoly) = tcont_poly_tmp
  end do
  tcont_ma(1) = tcont_ma_tmp
  do imatom = 2, nmatyp
     tcont_ma(imatom) = tcont_ma_tmp
  end do

  close(iuini)

!     +     +     +     +     +     +     +

end subroutine rdscript
