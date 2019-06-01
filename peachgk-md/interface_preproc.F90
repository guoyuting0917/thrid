!***********************************
!*  interface_preproc.f90 Ver.3.0  *
!*      for peachgk_md.f           *
!*            by G.Kikugawa        *
!***********************************
! Time-stamp: <>

!***** This module is interface module for preprocess routines of MD *****

module interface_preproc

  interface

     subroutine MPI_loopinit(startPROC)

       ! ARGUMENT:
       ! INPUT
       integer,intent(in):: startPROC       ! MPI first process ID

     end subroutine MPI_loopinit


     subroutine malloc_larray()

       ! ARGUMENT:
       !   INPUT

     end subroutine malloc_larray

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
       character(80),intent(out):: outhename ! output file name
       character(80),intent(out):: ouforname ! output file name
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

       logical,intent(out):: ifcstmnb          ! flag if using custom NB interaction
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

     end subroutine rdscript

#if defined(HF)
     subroutine rdheatf(ouhtfname, &
          &             oumtfname, &
          &             xref, &
          &             zcel, &
          &             npolytyp, &
          &             ifhfvol, &
          &             nhfregion, hfzpos1, hfzpos2, &
          &             hftyp_pmole, hftyp_water, &
          &             ifcalmtf, ifcalhtf, &
          &             mtfoutdir)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: xref             ! distanse base value [m]

       real(8),intent(in):: zcel             ! initial z0 cell length[m]

       integer,intent(in):: npolytyp        ! number of poly type

       !     OUTPUT
       character(80),intent(out):: ouhtfname ! output file name

       character(80),intent(out):: oumtfname ! output file name

       logical,intent(out):: ifhfvol       ! local volume-based or local surface-based

       integer,intent(out):: nhfregion       ! number of region to calculate heat flux
       real(8),intent(out):: hfzpos1(:),hfzpos2(:)
                          ! z-position of region for heat flux

       integer,intent(out):: hftyp_pmole(:)
                          ! atom- or mole-based heat flux cal. for poly
       integer,intent(out):: hftyp_water
                          ! atom- or mole-based heat flux cal. for water

       logical,intent(out):: ifcalmtf ! flag to calculate & output momentum flux
       logical,intent(out):: ifcalhtf ! flag to calculate & output heat flux

       character(3),intent(out):: mtfoutdir ! direction to output data of momentum

     end subroutine rdheatf
#endif

     subroutine rdpotbias( ouumbname, &
          &                xref, eref )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: eref             ! energy base value [J]

       !     OUTPUT
       character(80),intent(out):: ouumbname ! output file name

     end subroutine rdpotbias

     subroutine openfile(iucor,iucorname,iutop,iutopname, &
          &              iuwtop,iuwtopname, &
          &              iuparavdw,iuparavdwname,iuparabond,iuparabondname, &
          &              iuparaconst,iuparaconstname, &
          &              iuparacstmnb,iuparacstmnbname,ifcstmnb, &
          &              iuaddtop,iuaddtopname,ifrdaddtop, &
          &              iustrmvel,iustrmvelname,ifstrmvel, &
          &              iuposres,iuposresname,ifposres, &
          &              iostarec,iostarecname,ifstarec, &
          &              npolytyp, &
          &              ousum,ousumname, &
          &              ouene,ouenename,oupos,ouposname,ouvel,ouvelname, &
          &              oufor,ouforname, &
          &              outhe,outhename,oubar,oubarname, &
          &              oupre,ouprename,outhc,outhcname,ifoutthc, &
          &              oupdb,oupdbname,ifoutpdb, &
          &              ouhtf,ouhtfname, &
          &              oumtf,oumtfname, &
          &              ouumb,ouumbname,ifpotbias, &
          &              ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
          &              ifoutbar,ifoutpre)

       ! ARGUMENT:
       !     INPUT
       character(80),intent(in):: iucorname(:)    ! input file name
       character(80),intent(in):: iutopname(:)    !       "
       character(80),intent(in):: iuwtopname      !       "
       character(80),intent(in):: iuparavdwname   !       "
       character(80),intent(in):: iuparabondname  !       "
       character(80),intent(in):: iuparacstmnbname !       "
       character(80),intent(in):: iuparaconstname !       "
       character(80),intent(in):: iuaddtopname    !       "
       character(80),intent(in):: iustrmvelname    !       "
       character(80),intent(in):: iuposresname    !       "

       character(80),intent(in):: iostarecname ! state record file name
       logical,intent(in):: ifstarec        ! read old state

       character(80),intent(in):: ousumname ! output file name
       character(80),intent(in):: ouenename ! output file name
       character(80),intent(in):: ouposname ! output file name
       character(80),intent(in):: ouvelname ! output file name
       character(80),intent(in):: ouforname ! output file name
       character(80),intent(in):: outhename ! output file name
       character(80),intent(in):: oubarname ! output file name
       character(80),intent(in):: ouprename ! output file name
       character(80),intent(in):: outhcname ! output file name
       character(80),intent(in):: oupdbname ! output file name
       character(80),intent(in):: ouhtfname ! output file name
       character(80),intent(in):: oumtfname ! output file name
       character(80),intent(in):: ouumbname ! output file name

       integer,intent(in):: npolytyp       ! number of poly type

       logical,intent(in):: ifcstmnb      ! flag if using custom NB interaction

       logical,intent(in):: ifrdaddtop     ! input additional topology information
       logical,intent(in):: ifoutthc       ! flag for outputting thermal control file

       logical,intent(in):: ifoutpdb       ! flag for outputting PDB format file

       logical,intent(in):: ifposres       ! position restraint flag

       logical,intent(in):: ifpotbias      ! bias potential flag
       logical,intent(in):: ifstrmvel      ! flag to input and use streaming velocity

       logical,intent(in):: ifoutene       ! if ouput energy file
       logical,intent(in):: ifoutpos       ! if ouput position file
       logical,intent(in):: ifoutvel       ! if ouput velocity file
       logical,intent(in):: ifoutfor       ! if ouput force file
       logical,intent(in):: ifoutthe       ! if ouput NVT file
       logical,intent(in):: ifoutbar       ! if ouput NPT file
       logical,intent(in):: ifoutpre       ! if ouput pressure file

       !     OUTPUT
       integer,intent(out):: iucor(:)      ! input poly coordinate file unit
       integer,intent(out):: iutop(:)      ! input poly topology file unit
       integer,intent(out):: iuwtop        ! input poly topology file unit (water)
       integer,intent(out):: iuparavdw     ! input vdw parameter file unit
       integer,intent(out):: iuparabond    ! input bond parameter file unit
       integer,intent(out):: iuparaconst   ! input const parameter file unit
       integer,intent(out):: iuparacstmnb  ! input custom NB parameter file unit
       integer,intent(out):: iuaddtop      ! input additional topology file unit
       integer,intent(out):: iustrmvel     ! input streaming velocity file unit
       integer,intent(out):: iuposres      ! input position restraint ref. file unit

       integer,intent(out):: iostarec      ! state record file unit

       integer,intent(out):: ousum         ! output parameter summarization file unit
       integer,intent(out):: ouene         ! output unit for output energy data
       integer,intent(out):: oupos         ! output unit for output position data
       integer,intent(out):: ouvel         ! output unit for output velocity data
       integer,intent(out):: oufor         ! output unit for output force data
       integer,intent(out):: outhe         ! output unit for output thermostat data
       integer,intent(out):: oubar         ! output unit for output barostat data
       integer,intent(out):: oupre         ! output unit for output pressure data
       integer,intent(out):: outhc         ! output unit for outthc thermal control data
       integer,intent(out):: oupdb         ! output unit for outpdb PDB data
       integer,intent(out):: ouhtf         ! output unit for outhtf heat flux data
       integer,intent(out):: oumtf         ! output unit for outmtf momentum flux data
       integer,intent(out):: ouumb         ! output unit for outumb bias potential data

     end subroutine openfile

     subroutine calbase(xref,eref,mref,qref, &
          &             vref,timeref,tempref,pref,fref, &
          &             eps0ref, &
          &             npolytyp,nmatyp, &
          &             eps0, &
          &             alpha,rrcut, &
          &             xcel,ycel,zcel, &
          &             rcut, &
          &             tcont_poly,tcont_water,tcont_ma, &
          &             tcont_poly_ini,tcont_water_ini,tcont_ma_ini, &
          &             text,tfreq, &
          &             vfreq,pext, &
          &             rcutmor, &
          &             rcutsh, &
          &             rcutrfhfo,rcutrfhoo,rcutrfhoh, &
          &             rcutdouo,rcutindouo,rcutdouh, &
          &             rcutrpvw, &
          &             d_rini,d_rmax,d_econv,d_rmsf, &
          &             limitdist)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: eref             ! energy base value [J]
       real(8),intent(in):: mref             ! mass base value [kg]
       real(8),intent(in):: qref             ! charge base value [C]

       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.

       !     OUTPUT
       real(8),intent(out):: vref             ! velocity base value [m/s]
       real(8),intent(out):: timeref          ! time base value [sec]
       real(8),intent(out):: tempref          ! temperature base value [K]
       real(8),intent(out):: pref             ! pressure base value [Pa]
       real(8),intent(out):: fref             ! force base value [N]

       real(8),intent(out):: eps0ref         ! dielectric constant base value [c^2/Jm]

       !     INPUT&OUTPUT
       real(8),intent(inout):: alpha            ! parameter alpha [1/m->non-d]
       real(8),intent(inout):: rrcut
                            ! ewald real space cutoff length [m->non-d]
       real(8),intent(inout):: eps0             ! parameter eps0 [c^2/Jm->non-d]
       real(8),intent(inout):: xcel             ! x cell length [m->non-d]
       real(8),intent(inout):: ycel             ! y cell length [m->non-d]
       real(8),intent(inout):: zcel             ! z cell length [m->non-d]
       real(8),intent(inout):: rcut             ! vdw cutoff length [m->non-d]
       real(8),intent(inout):: tcont_poly(:)    ! poly Temp. [K->non-d]
       real(8),intent(inout):: tcont_water      ! H2O Temp. [K->non-d]
       real(8),intent(inout):: tcont_ma(:)      ! monatomic mole. Temp. [K->non-d]
       real(8),intent(inout):: tcont_poly_ini   ! poly Temp. [K->non-d]
       real(8),intent(inout):: tcont_water_ini  ! H2O Temp. [K->non-d]
       real(8),intent(inout):: tcont_ma_ini     ! MA Temp. [K->non-d]
       real(8),intent(inout):: text             ! external temperature [K->non-d]

       real(8),intent(inout):: tfreq            ! temperature frequency [1/s->non-d]

       real(8),intent(inout):: vfreq            ! volume change frequency [1/s->non-d]
       real(8),intent(inout):: pext             ! external pressure [Pa->non-d]

       real(8),intent(inout):: rcutmor          ! Morse cutoff length [m->non-d]

       real(8),intent(inout):: rcutsh           ! SH cutoff length [m->non-d]

       real(8),intent(inout):: rcutrfhfo        ! RFH(FO) cutoff length [m->non-d]
       real(8),intent(inout):: rcutrfhoo        ! RFH(OO) cutoff length [m->non-d]
       real(8),intent(inout):: rcutrfhoh        ! RFH(OH) cutoff length [m->non-d]

       real(8),intent(inout):: rcutdouo
                            ! DOU cutoff length [m->non-d] for O-Au
       real(8),intent(inout):: rcutindouo       ! DOU cutin length [m->non-d] for O-Au
       real(8),intent(inout):: rcutdouh
                            ! DOU cutoff length [m->non-d] for H-Au

       real(8),intent(inout):: rcutrpvw ! cutoff length for RP-VW interaction

       real(8),intent(inout):: d_rini   ! initial displacement dr for EM [m->non-d]
       real(8),intent(inout):: d_rmax   ! maximum displacement dr for EM [m->non-d]

       real(8),intent(inout):: d_econv  ! convergence condition for energy
                                   !    in EM [J->non-d]
       real(8),intent(inout):: d_rmsf   ! convergence condition for
                                   !    root mean square force in EM [N->non-d]

       real(8),intent(inout):: limitdist   ! maximum atomic displacement [m->non-d]

     end subroutine calbase

     subroutine rdcor( iucor,xref, &
          &            npoly, npolytyp, npoly_mole, npoly_atom, &
          &            nwater, nmatom, &
          &            nmatyp, nmatomtyp, &
          &            polytyp_free, &
          &            oatmtyp, hatmtyp, monoatmtyp, &
          &            ifsetchrg, atmchrg_tmp )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: iucor(:)        ! input poly coordinate file unit
       real(8),intent(in):: xref             ! distanse base value [m]

       integer,intent(in):: npoly           ! all number of poly
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control

       character(2),intent(in):: oatmtyp ! O atomtype of water model
       character(2),intent(in):: hatmtyp ! H atomtype of water model
       character(2),intent(in):: monoatmtyp(:) ! monatomic mole. type

       logical,intent(in):: ifsetchrg(:)    ! if set charge from cor file

       !     OUTPUT
       real(8),intent(out):: atmchrg_tmp(:) ! temporary atom charge read from cor file

     end subroutine rdcor

     subroutine rdtop( iutop, iuwtop, &
          &            npolytyp, npoly_mole, npoly_atom, &
          &            nwater )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: iutop(:)        ! input poly topology file unit
       integer,intent(in):: iuwtop          ! input water topology file unit
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)
                          ! number of atoms belonging to poly
       integer,intent(in):: nwater          ! number of H2O molecules

     end subroutine rdtop

     subroutine rdpara(iuparavdw,iuparabond, &
          &            xref,eref,mref, &
          &            xcel,ycel,zcel, &
          &            ifljari,ifljgeo)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: iuparavdw       ! input vdw parameter file unit
       integer,intent(in):: iuparabond      ! input bond parameter file unit

       real(8),intent(in):: xref            ! distanse base value [m]
       real(8),intent(in):: eref            ! energy base value [J]
       real(8),intent(in):: mref            ! mass base value [kg]

       real(8),intent(in):: xcel            ! x cell length [non-d]
       real(8),intent(in):: ycel            ! y cell length [non-d]
       real(8),intent(in):: zcel            ! z cell length [non-d]

       logical,intent(in):: ifljari         ! arithmetic mean for LJ cross parameter
       logical,intent(in):: ifljgeo         ! geometric mean for LJ cross parameter

     end subroutine rdpara

     subroutine rdpara_cstmnb(iuparacstmnb, &
          &                   ifcellindex_cstmnb,ifbookcstmnb, &
          &                   xref,eref,mref,qref, &
          &                   vref,timeref,tempref,pref,fref,eps0ref, &
          &                   xcel,ycel,zcel)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: iuparacstmnb    ! input custom NB parameter file unit

       logical,intent(inout):: ifcellindex_cstmnb ! flag for cell index (custom NB)
       logical,intent(in):: ifbookcstmnb
                               ! flag for bookkeeping of custom NB interaction

       real(8),intent(in):: xref            ! distanse base value [m]
       real(8),intent(in):: eref            ! energy base value [J]
       real(8),intent(in):: mref            ! mass base value [kg]
       real(8),intent(in):: qref            ! charge base value [C]
       real(8),intent(in):: vref            ! velocity base value [m/s]
       real(8),intent(in):: timeref         ! time base value [sec]
       real(8),intent(in):: tempref         ! temperature base value [K]
       real(8),intent(in):: pref            ! pressure base value [Pa]
       real(8),intent(in):: fref            ! force base value [N]

       real(8),intent(in):: eps0ref         ! dielectric constant base value [c^2/Jm]

       real(8),intent(in):: xcel            ! x cell length [non-d]
       real(8),intent(in):: ycel            ! y cell length [non-d]
       real(8),intent(in):: zcel            ! z cell length [non-d]

     end subroutine rdpara_cstmnb

     subroutine rdconst( iuparaconst )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: iuparaconst     ! input const parameter file unit

     end subroutine rdconst

     subroutine mkmolept( npoly, npolytyp, npoly_mole, npoly_atom, &
          &               nwater, nmatom )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: npoly           ! all number of poly
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
       integer,intent(in):: nwater          ! number of H2O molecules
       integer,intent(in):: nmatom          ! number of monatomic molecules

     end subroutine mkmolept

     subroutine rdrefcor( iuposres, &
          &               npoly, npolytyp, npoly_mole, npoly_atom, &
          &               nwater, nmatom, &
          &               xref )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: iuposres        ! input position restraint ref. file unit

       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
       integer,intent(in):: nwater          ! number of H2O molecules
       integer,intent(in):: nmatom          ! number of monatomic molecules

       real(8),intent(in):: xref             ! distanse base value [m]

     end subroutine rdrefcor

     subroutine rdstarec(iostarec, &
          &              npoly, npolytyp, npoly_mole, npoly_atom, &
          &              nwater, nmatom, &
          &              xcel, ycel, zcel, yratio, zratio, &
          &              xref, vref, timeref, pref, &
          &              mchain, &
          &              pint, pintt, &
          &              ifsetcor)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: iostarec        ! state record file unit

       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
       integer,intent(in):: nwater          ! number of H2O molecules
       integer,intent(in):: nmatom          ! number of monatomic molecules

       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: vref             ! velocity base value [m/s]
       real(8),intent(in):: timeref          ! time base value [sec]
       real(8),intent(in):: pref             ! pressure base value [Pa]

       integer,intent(in):: mchain          ! Nose-Hoover chain number

       !     OUTPUT
       real(8),intent(out):: xcel             ! x cell length
       real(8),intent(out):: ycel             ! y cell length
       real(8),intent(out):: zcel             ! z cell length
       real(8),intent(out):: yratio           ! y cell ratio of y to x
       real(8),intent(out):: zratio           ! z cell ratio of z to x

       real(8),intent(out):: pint             ! internal pressure
       real(8),intent(out):: pintt(:,:)       ! internal pressure tensor

       logical,intent(out):: ifsetcor(:)     ! frag for checking if coordinate has been set

     end subroutine rdstarec

     subroutine createcor(npoly,npolytyp,npoly_mole,npoly_atom,   &
          &               nwater,nmatom,nmatyp,nmatomtyp,   &
          &               xmaxpo,ymaxpo,zmaxpo,   &
          &               xmaxw,ymaxw,zmaxw,   &
          &               xmaxma,ymaxma,zmaxma,   &
          &               inicorpo,inicorw,inicorma,   &
          &               ncrecorpo,index_crecorpo,   &
          &               ncrecorw,index_crecorw,   &
          &               ncrecorma,index_crecorma,   &
          &               ifsetcor,   &
          &               xcel,ycel,zcel,   &
          &               xref,   &
          &               compfact,   &
          &               mchain)

       ! ARGUMENT:
       !   INPUT
       integer,intent(in):: npoly            ! all number of poly
       integer,intent(in):: npolytyp         ! number of poly type
       integer,intent(in):: npoly_mole(:)    ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)    ! number of atoms belonging to poly

       integer,intent(in):: nwater           ! number of H2O molecules

       integer,intent(in):: nmatom           ! number of monatomic molecules
       integer,intent(in):: nmatyp           ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)     ! each number of monatomic mole.

       integer,intent(in):: xmaxpo(:)        ! use for positioning of polymer1
       integer,intent(in):: ymaxpo(:)        ! use for positioning of polymer1
       integer,intent(in):: zmaxpo(:)        ! use for positioning of polymer1
       integer,intent(in):: xmaxw            ! use for positioning of water
       integer,intent(in):: ymaxw            ! use for positioning of water
       integer,intent(in):: zmaxw            ! use for positioning of water
       integer,intent(in):: xmaxma(:)        ! use for positioning of monatomic mole.
       integer,intent(in):: ymaxma(:)        ! use for positioning of monatomic mole.
       integer,intent(in):: zmaxma(:)        ! use for positioning of monatomic mole.
       real(8),intent(in):: inicorpo(:,:) ! use for positioning of poly
       real(8),intent(in):: inicorw(:)       ! use for positioning of water
       real(8),intent(in):: inicorma(:,:) ! use for positioning of ma
       integer,intent(in):: ncrecorpo        ! max number of poly type for createcor
       integer,intent(in):: index_crecorpo(:) ! index of polymer for createcor
       integer,intent(in):: ncrecorw         ! water for createcor
       integer,intent(in):: index_crecorw    ! index of water type for createcor
       integer,intent(in):: ncrecorma        ! max number of matom type for createcor
       integer,intent(in):: index_crecorma(:) ! index of matom type for createcor

       logical,intent(out):: ifsetcor(:)     ! frag for checking if coordinate has set

       real(8),intent(in):: xref             ! distanse base value [m]

       real(8),intent(in):: xcel             ! x cell length
       real(8),intent(in):: ycel             ! y cell length
       real(8),intent(in):: zcel             ! z cell length

       real(8),intent(in):: compfact     ! compact factor using at poly arrange(<1.0)

       integer,intent(in):: mchain           ! Nose-Hoover chain number

     end subroutine createcor

     subroutine rdtcont( xref, tempref, &
          &              xcel, ycel, zcel, &
          &              ntcregion, &
          &              tcxpos1, tcxpos2, &
          &              tcypos1, tcypos2, &
          &              tczpos1, tczpos2, &
          &              r_tcont )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: tempref          ! temperature base value [K]

       real(8),intent(in):: xcel             ! initial x0 cell length[m]
       real(8),intent(in):: ycel             ! initial y0 cell length[m]
       real(8),intent(in):: zcel             ! initial z0 cell length[m]

       !     OUTPUT
       integer,intent(out):: ntcregion       ! number of region to control temp.
       real(8),intent(out):: tcxpos1(:),tcxpos2(:)
                                ! x-position of temp. control region
       real(8),intent(out):: tcypos1(:),tcypos2(:)
                                ! y-position of temp. control region
       real(8),intent(out):: tczpos1(:),tczpos2(:)
                                ! z-position of temp. control region
       real(8),intent(out):: r_tcont(:)       ! control temp. in each region

     end subroutine rdtcont

     subroutine rdhfcont( xref, eref, timeref, &
          &               xcel, ycel, zcel, &
          &               dt_long_cal, &
          &               tcontinterval, &
          &               nhfcregion, &
          &               hfczpos1,hfczpos2, &
          &               r_hfcont )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: eref             ! energy base value [J]
       real(8),intent(in):: timeref          ! time base value [sec]

       real(8),intent(in):: xcel             ! initial x0 cell length[non-d]
       real(8),intent(in):: ycel             ! initial y0 cell length[non-d]
       real(8),intent(in):: zcel             ! initial z0 cell length[non-d]

       real(8),intent(in):: dt_long_cal      ! time step of long force [sec]
       integer,intent(in):: tcontinterval   ! interval of temp. control

       !     OUTPUT
       integer,intent(out):: nhfcregion      ! number of region to control heat flux
       real(8),intent(out):: hfczpos1(:),hfczpos2(:)
                                ! z-position of heat flux control region
       real(8),intent(out):: r_hfcont(:)      ! magnitude of heat flux in each region
                                ! (converted to the input energy)

     end subroutine rdhfcont

     subroutine rdlangevin(xref, tempref, timeref, &
          &                xcel, ycel, zcel, &
          &                nlangeregion, &
          &                ltxpos1, ltxpos2, &
          &                ltypos1, ltypos2, &
          &                ltzpos1, ltzpos2, &
          &                r_ltemp, r_ltdamp)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: tempref          ! temperature base value [K]
       real(8),intent(in):: timeref          ! time base value [sec]

       real(8),intent(in):: xcel             ! initial x0 cell length[m]
       real(8),intent(in):: ycel             ! initial y0 cell length[m]
       real(8),intent(in):: zcel             ! initial z0 cell length[m]

       !     OUTPUT
       integer,intent(out):: nlangeregion    ! number of region for Langevin thermo.
       real(8),intent(out):: ltxpos1(:),ltxpos2(:)
                                        ! x-position of temp. control region
       real(8),intent(out):: ltypos1(:),ltypos2(:)
                                        ! y-position of temp. control region
       real(8),intent(out):: ltzpos1(:),ltzpos2(:)
                                        ! z-position of temp. control region
       real(8),intent(out):: r_ltemp(:)      ! control temp. in each region
       real(8),intent(out):: r_ltdamp(:)
                                 ! damping factor in each region [1/s -> non-d]

     end subroutine rdlangevin

     subroutine rdaddtop(iuaddtop, &
          &              mref, &
          &              ifsetatmmass,ifsetatmchrg)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: iuaddtop        ! input additional topology file unit

       real(8),intent(in):: mref            ! mass base value [kg] (C atom)

       !     OUTPUT
       logical,intent(out):: ifsetatmmass(:) ! flag for setting mass by add_top
       logical,intent(out):: ifsetatmchrg(:) ! flag for setting charge by add_top

     end subroutine rdaddtop

     subroutine rdstrmvel(iustrmvel, &
          &               xref,vref, &
          &               xcel,ycel,zcel)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: iustrmvel       ! input streaming velocity file unit

       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: vref             ! velocity base value [m/s]

       real(8),intent(in):: xcel             ! initial x0 cell length[non-d]
       real(8),intent(in):: ycel             ! initial y0 cell length[non-d]
       real(8),intent(in):: zcel             ! initial z0 cell length[non-d]

     end subroutine rdstrmvel

     subroutine mkexcl()

       ! ARGUMENT:
       !   NONE

     end subroutine mkexcl

     subroutine linkbond(npolytyp,npoly_mole,npoly_atom, &
          &              polytyp_free, &
          &              atmchrg_tmp, &
          &              ifsetatmmass,ifsetatmchrg)

       ! ARGUMENTS:
       ! INPUT
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
       real(8),intent(in):: atmchrg_tmp(:)  ! temporary atom charge read from cor file

       logical,intent(in):: ifsetatmmass(:) ! flag for setting mass by add_top
       logical,intent(in):: ifsetatmchrg(:) ! flag for setting charge by add_top

     end subroutine linkbond

     subroutine mkconst( npolytyp, npoly_mole )

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly

     end subroutine mkconst

     subroutine createvel( npoly, npolytyp, npoly_mole, npoly_atom, &
          &                nwater, nmatom, nmatyp, nmatomtyp, &
          &                xmaxpo, ymaxpo, zmaxpo, &
          &                xmaxw, ymaxw, zmaxw, &
          &                xmaxma, ymaxma, zmaxma, &
          &                ncrecorpo, index_crecorpo, &
          &                ncrecorw, index_crecorw, &
          &                ncrecorma, index_crecorma, &
          &                vref, timeref, &
          &                mchain, &
          &                tcont_poly, tcont_water, tcont_ma )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: npoly           ! all number of poly
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       integer,intent(in):: xmaxpo(:)       ! use for positioning of polymer1
       integer,intent(in):: ymaxpo(:)       ! use for positioning of polymer1
       integer,intent(in):: zmaxpo(:)       ! use for positioning of polymer1
       integer,intent(in):: xmaxw           ! use for positioning of water
       integer,intent(in):: ymaxw           ! use for positioning of water
       integer,intent(in):: zmaxw           ! use for positioning of water
       integer,intent(in):: xmaxma(:)       ! use for positioning of monatomic mole.
       integer,intent(in):: ymaxma(:)       ! use for positioning of monatomic mole.
       integer,intent(in):: zmaxma(:)       ! use for positioning of monatomic mole.

       integer,intent(in):: ncrecorpo       ! max number of poly type for createcor
       integer,intent(in):: index_crecorpo(:) ! index of polymer for createcor
       integer,intent(in):: ncrecorw        ! water for createcor
       integer,intent(in):: index_crecorw   ! index of water type for createcor
       integer,intent(in):: ncrecorma       ! max number of matom type for createcor
       integer,intent(in):: index_crecorma(:) ! index of matom type for createcor

       real(8),intent(in):: vref             ! velocity base value [m/s]
       real(8),intent(in):: timeref          ! time base value [sec]

       integer,intent(in):: mchain          ! Nose-Hoover chain number

       real(8),intent(in):: tcont_poly(:)    ! poly Temp. [non-d]in NVT
       real(8),intent(in):: tcont_water     ! H2O Temp. [non-d] in NVT
       real(8),intent(in):: tcont_ma(:)      ! monatomic mole. Temp. [non-d] in NVT

     end subroutine createvel

     subroutine calmolmass()

       ! ARGUMENTS:

     end subroutine calmolmass

     subroutine prelocalheat( npolytyp, nmatyp, &
          &                   polytyp_free, watertyp_free, matomtyp_free, &
          &                   tempref, &
          &                   nlheat_poly, index_nlheat_poly, &
          &                   tcont_nlheat_poly, &
          &                   nlheat_water, tcont_nlheat_water, &
          &                   nlheat_ma, index_nlheat_ma, &
          &                   tcont_nlheat_ma )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.

       character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
       character(80),intent(in):: watertyp_free(:) ! use for water type control
       character(80),intent(in):: matomtyp_free(:,:)
                                ! use for matom type control

       real(8),intent(in):: tempref          ! temperature base value [K]

       !     OUTPUT
       integer,intent(out):: nlheat_poly     ! number of poly type for local heating
       integer,intent(out):: index_nlheat_poly(:)
                                ! index of poly type for local heating
       real(8),intent(out):: tcont_nlheat_poly(:)
                                ! control temp. of poly type for local heating
       integer,intent(out):: nlheat_water    ! number of water for local heating
       real(8),intent(out):: tcont_nlheat_water
                                ! control temp. of water for local heating
       integer,intent(out):: nlheat_ma       ! number of matom type for local heating
       integer,intent(out):: index_nlheat_ma(:)
                                ! index of matom type for local heating
       real(8),intent(out):: tcont_nlheat_ma(:)
                                ! control temp. of matom type for local heating

     end subroutine prelocalheat

     subroutine prelocalfix( npolytyp, npoly_mole, npoly_atom, &
          &                  nwater, &
          &                  nmatyp, nmatomtyp, &
          &                  polytyp_free, watertyp_free, matomtyp_free, &
          &                  nlfix,index_nlfix, &
          &                  nlfix_deg_poly, &
          &                  nlfix_deg_water, &
          &                  nlfix_deg_ma )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
       character(80),intent(in):: watertyp_free(:) ! use for water type control
       character(80),intent(in):: matomtyp_free(:,:)
                                ! use for matom type control

       !     OUTPUT
       integer,intent(out):: nlfix           ! number of fix atoms
       integer,intent(out):: index_nlfix(:)  ! index of fix atoms

       integer,intent(out):: nlfix_deg_poly(:) ! fixed degree of freedom of poly
       integer,intent(out):: nlfix_deg_water ! fixed degree of freedom of H2O
       integer,intent(out):: nlfix_deg_ma(:) ! fixed degree of freedom of matom

     end subroutine prelocalfix

     subroutine prelocalfixz(xref,   &
          &                  npolytyp, npoly_mole, npoly_atom,   &
          &                  nwater,   &
          &                  nmatyp,nmatomtyp,   &
          &                  polytyp_free, watertyp_free, matomtyp_free,   &
          &                  nlfixz_deg_poly, nlfixz_deg_water, nlfixz_deg_ma)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: xref        ! distanse base value [m] (adjust to Angstrom)

       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
       character(80),intent(in):: watertyp_free(:) ! use for water type control
       character(80),intent(in):: matomtyp_free(:,:)
                                ! use for matom type control

       !    OUTPUT
       integer,intent(out):: nlfixz_deg_poly(:) ! fixed degree of freedom of poly
       integer,intent(out):: nlfixz_deg_water   ! fixed degree of freedom of H2O
       integer,intent(out):: nlfixz_deg_ma(:)   ! fixed degree of freedom of matom

     end subroutine prelocalfixz

     subroutine prelocalfixzg(xref,   &
          &                   npolytyp,npoly_mole,npoly_atom,   &
          &                   nwater,   &
          &                   nmatyp,nmatomtyp,   &
          &                   polytyp_free,watertyp_free,matomtyp_free,   &
          &                   nlfixzg_deg_poly,nlfixzg_deg_water,nlfixzg_deg_ma)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: xref        ! distanse base value [m] (adjust to Angstrom)

       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
       character(80),intent(in):: watertyp_free(:) ! use for water type control
       character(80),intent(in):: matomtyp_free(:,:)
                                ! use for matom type control

       !    OUTPUT
       integer,intent(out):: nlfixzg_deg_poly(:) ! fixed degree of freedom of poly
       integer,intent(out):: nlfixzg_deg_water   ! fixed degree of freedom of H2O
       integer,intent(out):: nlfixzg_deg_ma(:)   ! fixed degree of freedom of matom

     end subroutine prelocalfixzg

     subroutine preposres( npolytyp, npoly_mole, npoly_atom, &
          &                nwater, &
          &                nmatyp,nmatomtyp, &
          &                polytyp_free, watertyp_free, matomtyp_free )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
       character(80),intent(in):: watertyp_free(:) ! use for water type control
       character(80),intent(in):: matomtyp_free(:,:)
                                ! use for matom type control

     end subroutine preposres

     subroutine prepdbout( npolytyp, polytyp_free, &
          &                watertyp_free, &
          &                nmatyp, matomtyp_free, monoatmtyp, &
          &                resname_poly_pdb, &
          &                resname_water_pdb, &
          &                resname_matom_pdb )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: npolytyp        ! number of poly type
       character(80),intent(in):: polytyp_free(:,:)

       character(80),intent(in):: watertyp_free(:) ! use for water type control

       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       character(80),intent(in):: matomtyp_free(:,:)

       character(2),intent(in):: monoatmtyp(:) ! monatomic mole. type

       !     OUTPUT
       character(4),intent(out):: resname_poly_pdb(:)
       character(4),intent(out):: resname_water_pdb ! residue name for poly
       character(4),intent(out):: resname_matom_pdb(:) ! residue name for matom

     end subroutine prepdbout

#if defined(HF)
     subroutine preheatf( npoly, npolytyp, npoly_mole, &
          &               nwater, &
          &               nmatom, &
          &               hftyp_pmole, hftyp_water, &
          &               hftyp_atm )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: npoly           ! all number of poly
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules

       integer,intent(in):: hftyp_pmole(:)
                         ! atom- or mole-based heat flux cal. for poly
       integer,intent(in):: hftyp_water ! atom- or mole-based heat flux cal. for water

       !     OUTPUT
       integer,intent(out):: hftyp_atm(:)    ! atom- or mole-based heat flux cal.
                                !   for each atom

     end subroutine preheatf
#endif

     subroutine prepotbias( npolytyp, npoly_mole, npoly_atom, &
          &                 nwater, &
          &                 nmatyp,nmatomtyp, &
          &                 polytyp_free, watertyp_free, matomtyp_free )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
       character(80),intent(in):: watertyp_free(:) ! use for water type control
       character(80),intent(in):: matomtyp_free(:,:)
                                ! use for matom type control

     end subroutine prepotbias

     subroutine prelocalvel(vref, &
          &                 npolytyp,npoly_mole,npoly_atom, &
          &                 nwater, &
          &                 nmatyp,nmatomtyp, &
          &                 polytyp_free,watertyp_free,matomtyp_free, &
          &                 nlvel,index_nlvel,v_nlvel, &
          &                 nlvel_deg_poly, &
          &                 nlvel_deg_water, &
          &                 nlvel_deg_ma)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: vref            ! velocity base value [m/s]

       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
       character(80),intent(in):: watertyp_free(:) ! use for water type control
       character(80),intent(in):: matomtyp_free(:,:)
                                ! use for matom type control

       !     OUTPUT
       integer,intent(out):: nlvel           ! number of atoms for velocity fix
       integer,intent(out):: index_nlvel(:)  ! index of vel-fix atoms
       real(8),intent(out):: v_nlvel(:,:)    ! local velocity values

       integer,intent(out):: nlvel_deg_poly(:) ! fixed degree of freedom of poly
       integer,intent(out):: nlvel_deg_water ! fixed degree of freedom of H2O
       integer,intent(out):: nlvel_deg_ma(:) ! fixed degree of freedom of matom

     end subroutine prelocalvel

     subroutine prepmd(npoly, npolytyp, npoly_mole, npoly_atom, &
          &            nwater, nmatom,nmatyp, nmatomtyp, &
          &            ifcenterfix_all, &
          &            ifcenterfix_poly, ifcenterfix_water, ifcenterfix_ma, &
          &            cenfix_free, &
          &            ifcenterfix_polytyp, &
          &            ifcenterfix_watertyp, &
          &            ifcenterfix_matyp, &
          &            degfree_poly, degfree_water, &
          &            degfree_ma, degfree_all, &
          &            iflocalfix, iflocalfixz, iflocalfixzg, &
          &            nlfix_deg_poly, &
          &            nlfix_deg_water, &
          &            nlfix_deg_ma, &
          &            nlfixz_deg_poly, nlfixz_deg_water, nlfixz_deg_ma, &
          &            nlfixzg_deg_poly, nlfixzg_deg_water, nlfixzg_deg_ma, &
          &            iflocalvel, &
          &            nlvel_deg_poly,nlvel_deg_water,nlvel_deg_ma, &
          &            dt_short_cal, dt_med_cal, dt_long_cal, &
          &            nstep_short, nstep_med, &
          &            xref, timeref, eps0, &
          &            rcut, ifbook, rcut_book, nstep_book, &
          &            rcutmor, ifbookmor, rcut_bookmor, nstep_bookmor, &
          &            rcutsh, ifbooksh, rcut_booksh, nstep_booksh, &
          &            rcutrfhfo, ifbookrfhfo, rcut_bookrfhfo, &
          &            nstep_bookrfhfo, &
          &            rcutrfhoo, ifbookrfhoo, rcut_bookrfhoo, &
          &            nstep_bookrfhoo, &
          &            rcutrfhoh, ifbookrfhoh, rcut_bookrfhoh, &
          &            nstep_bookrfhoh, &
          &            rcutdouo, ifbookdouo, rcut_bookdouo, &
          &            nstep_bookdouo, &
          &            rcutdouh, ifbookdouh, rcut_bookdouh, &
          &            nstep_bookdouh, &
          &            rcutrpvw, ifbookrpvw, rcut_bookrpvw, nstep_bookrpvw, &
          &            tcont_poly, tcont_water, tcont_ma, &
          &            xcel, ycel, zcel, &
          &            ifrattle, &
          &            mchain, tfreq, text, &
          &            vfreq, &
          &            pcont_axis, &
          &            nyosh, &
          &            solvetyp, solveindex, &
          &            netchrgsq)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: npoly           ! all number of poly
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.


       logical,intent(in):: ifcenterfix_all   ! center fix for all
       logical,intent(in):: ifcenterfix_poly  ! center fix for polymer1
       logical,intent(in):: ifcenterfix_water ! center fix for water
       logical,intent(in):: ifcenterfix_ma    ! center fix for monatomic mole.
       character(4),intent(in):: cenfix_free  ! COM not fixed in this direction

       logical,intent(in):: ifcenterfix_polytyp(:) ! center fix for each polymer
       logical,intent(in):: ifcenterfix_watertyp   ! center fix for each water
       logical,intent(in):: ifcenterfix_matyp(:)
                                          ! center fix for each monatomic mole.

       integer,intent(in):: nstep_med       ! number of step for medium force
       integer,intent(in):: nstep_short     ! number of step for short force

       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: timeref          ! time base value [sec]
       real(8),intent(in):: eps0             ! dielectric constant of vacuum [non-d]

       real(8),intent(in):: tcont_poly(:)    ! poly Temp. [non-d] in NVT
       real(8),intent(in):: tcont_water      ! H2O Temp. [non-d] in NVT
       real(8),intent(in):: tcont_ma(:)      ! monatomic mole. Temp. [non-d] in NVT

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(inout):: rcut             ! vdw cutoff length [non-d]
       logical,intent(in):: ifbook          ! flag for bookkeeping
       real(8),intent(inout):: rcut_book    ! cut off radius of bookkeeping [m->non-d]
       integer,intent(in):: nstep_book      ! bookkeeping interval

       real(8),intent(inout):: rcutmor          ! Morse cutoff length [non-d]
       logical,intent(in):: ifbookmor      ! flag for bookkeeping of Morse interaction
       real(8),intent(inout):: rcut_bookmor
                            ! cut off radius of bookkeeping[m] of Morse
       integer,intent(in):: nstep_bookmor  ! bookkeeping interval of Morse interaction

       real(8),intent(inout):: rcutsh           ! SH cutoff length [non-d]
       logical,intent(in):: ifbooksh        ! flag for bookkeeping of SH interaction
       real(8),intent(inout):: rcut_booksh   ! cut off radius of bookkeeping[m] of SH
       integer,intent(in):: nstep_booksh    ! bookkeeping interval of SH interaction

       real(8),intent(inout):: rcutrfhfo        ! RFH(FO) cutoff length [non-d]
       logical,intent(in):: ifbookrfhfo  ! flag for bookkeeping of RFH(FO) interaction
       real(8),intent(inout):: rcut_bookrfhfo
                         ! cut off radius of bookkeeping[m] of RFH(FO)
       integer,intent(in):: nstep_bookrfhfo
                         ! bookkeeping interval of RFH(FO) interaction

       real(8),intent(inout):: rcutrfhoo        ! RFH(OO) cutoff length [non-d]
       logical,intent(in):: ifbookrfhoo  ! flag for bookkeeping of RFH(OO) interaction
       real(8),intent(inout):: rcut_bookrfhoo
                         ! cut off radius of bookkeeping[m] of RFH(OO)
       integer,intent(in):: nstep_bookrfhoo
                         ! bookkeeping interval of RFH(OO) interaction

       real(8),intent(inout):: rcutrfhoh        ! RFH(OH) cutoff length [non-d]
       logical,intent(in):: ifbookrfhoh  ! flag for bookkeeping of RFH(OH) interaction
       real(8),intent(inout):: rcut_bookrfhoh
                         ! cut off radius of bookkeeping[m] of RFH(OH)
       integer,intent(in):: nstep_bookrfhoh
                         ! bookkeeping interval of RFH(OH) interaction

       real(8),intent(inout):: rcutdouo         ! DOU cutoff length [non-d] for O-Au
       logical,intent(in):: ifbookdouo
                         ! flag for bookkeeping of DOU interaction (O-Au)
       real(8),intent(inout):: rcut_bookdouo
                         ! cut off radius of bookkeep[non-d] of DOU (O-Au)
       integer,intent(in):: nstep_bookdouo
                         ! bookkeeping interval of DOU interaction (O-Au)

       real(8),intent(inout):: rcutdouh         ! DOU cutoff length [non-d] for H-Au
       logical,intent(in):: ifbookdouh
                         ! flag for bookkeeping of DOU interaction (H-Au)
       real(8),intent(inout):: rcut_bookdouh
                         ! cut off radius of bookkeep[non-d] of DOU (H-Au)
       integer,intent(in):: nstep_bookdouh
                         ! bookkeeping interval of DOU interaction (H-Au)

       real(8),intent(inout):: rcutrpvw          ! RP-VW cutoff length [non-d]
       logical,intent(in):: ifbookrpvw      ! flag for bookkeeping of RP-VW interaction
       real(8),intent(inout):: rcut_bookrpvw
                         ! cut off radius of bookkeeping[m] of RP-VW
       integer,intent(in):: nstep_bookrpvw  ! bookkeeping interval of RP-VW interaction

       logical,intent(in):: ifrattle        ! rattle flag

       integer,intent(in):: mchain          ! Nose-Hoover chain number
       real(8),intent(in):: tfreq            ! temperature frequency [non-d]
       real(8),intent(in):: text          ! external temp. [non-d] (Nose-Hoover chain)

       real(8),intent(in):: vfreq            ! volume change frequency [non-d]

       character(5),intent(in):: pcont_axis
                                 ! axis for pressure control (iso, aniso, etc.)

       integer,intent(in):: nyosh           ! expansion order of Yoshida-Suzuki method

       character(2),intent(in):: solvetyp ! solvent molecule

       logical,intent(in):: iflocalfix      ! fix atoms flag
       integer,intent(in):: nlfix_deg_poly(:) ! fixed degree of freedom of poly
       integer,intent(in):: nlfix_deg_water ! fixed degree of freedom of H2O
       integer,intent(in):: nlfix_deg_ma(:) ! fixed degree of freedom of matom

       logical,intent(in):: iflocalfixz     ! flag for fixing z coordinate of atoms
       integer,intent(in):: nlfixz_deg_poly(:) ! fixed degree of freedom of poly
       integer,intent(in):: nlfixz_deg_water   ! fixed degree of freedom of H2O
       integer,intent(in):: nlfixz_deg_ma(:)   ! fixed degree of freedom of matom

       logical,intent(in):: iflocalfixzg    ! flag for fixing z coordinate of COM of molecules
       integer,intent(in):: nlfixzg_deg_poly(:) ! fixed degree of freedom of poly
       integer,intent(in):: nlfixzg_deg_water   ! fixed degree of freedom of H2O
       integer,intent(in):: nlfixzg_deg_ma(:)   ! fixed degree of freedom of matom

       logical,intent(in):: iflocalvel          ! flag to force local atomic velocity

       integer,intent(in):: nlvel_deg_poly(:)   ! fixed degree of freedom of poly
       integer,intent(in):: nlvel_deg_water     ! fixed degree of freedom of H2O
       integer,intent(in):: nlvel_deg_ma(:)     ! fixed degree of freedom of matom

       !     INPUT&OUTPUT
       real(8),intent(inout):: dt_long_cal      ! time step of long force [sec->non-d]

       !     OUTPUT
       integer,intent(out):: degfree_poly(:) ! degree of freedom of polymer1
       integer,intent(out):: degfree_water   ! degree of freedom of H2O
       integer,intent(out):: degfree_ma(:)   ! degree of freedom of monatomic mole.
       integer,intent(out):: degfree_all     ! degree of freedom of all molecules

       real(8),intent(out):: dt_med_cal       ! time step of medium force [non-d]
       real(8),intent(out):: dt_short_cal     ! time step of short force [non-d]

       !      real*8:: vdw_welij_solve  ! well depth of vdw parameter of solvent
       !      real*8:: vdw_radij_solve  ! radius of vdw parameter of solvent
       integer,intent(out):: solveindex      ! atmindex of solvent atom

       real(8),intent(out):: netchrgsq           ! = (sum(qi))**2

     end subroutine prepmd

     subroutine prepewk(pot_ewc, alpha, kmax, &
          &             xcel, ycel, zcel, &
          &             yratio, zratio)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: alpha            ! parameter alpha [non-d]
       integer,intent(in):: kmax           ! parameter kmax

       real(8),intent(in):: xcel             ! x cell length [non-d]
       real(8),intent(in):: ycel             ! y cell length [non-d]
       real(8),intent(in):: zcel             ! z cell length [non-d]

       real(8),intent(in):: yratio           ! y cell ratio of y to x
       real(8),intent(in):: zratio           ! z cell ratio of z to x

       !    OUTPUT
       real(8),intent(out):: pot_ewc     ! potential of ewald self-energy [non-d]


     end subroutine prepewk

     subroutine fft_pme_init( nfft1, nfft2, nfft3, pme_order, &
          &                   alpha, &
          &                   pot_ewc )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: nfft1, nfft2, nfft3 ! grid points in SPME
       integer,intent(in):: pme_order       ! B-spline order

       real(8),intent(in):: alpha            ! parameter alpha [1/m]

       real(8),intent(out):: pot_ewc          ! potential of ewald self-energy

     end subroutine fft_pme_init

    subroutine prepfennell( pot_ewc, alpha, rrcut, &
          &                  iffennell )

       ! ARGUMENT:
      !     INPUT
      real(8),intent(in):: alpha           ! parameter alpha [non-d]
      real(8),intent(in):: rrcut           ! ewald real space cutoff length [m]

      logical,intent(in):: iffennell       ! Fennell flag

      !     OUTPUT
      real(8),intent(out):: pot_ewc      ! potential of ewald self-energy [non-d]

     end subroutine prepfennell

#if defined(_NOT_SPLINTERP)
#else
     subroutine spline_interp_init( alpha, &
          &                         rrcut )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: alpha            ! ewald parameter alpha [non-dimension]
       real(8),intent(in):: rrcut            ! ewald real space cutoff length [non-d]

     end subroutine spline_interp_init
#endif

     subroutine wrsumm(ousum, &
          &            xref,eref,mref,qref, &
          &            vref,timeref,tempref,pref,fref, &
          &            eps0ref, &
          &            xcel,ycel,zcel, &
          &            ifewald,alpha,kmax,rrcut, &
          &            ifspme,nfft1,nfft2,nfft3,pme_order, &
          &            iffennell, &
          &            ifljari,ifljgeo, &
          &            iflocalheat, &
          &            ifregionheat,ifregionhf,ifreglange,iftcratom,ifoutthc, &
          &            iflocalfix,iflocalfixz,iflocalfixzg, &
          &            ifcnp, &
          &            ifposres, &
          &            ifpotbias, &
          &            iflocalvel, ifstrmvel, &
          &            dt_long_cal,dt_med_cal,dt_short_cal, &
          &            nstep_med,nstep_short, &
          &            maxnstep,nstage,nstep_stage,mdcont_stage, &
          &            nstep_maxwell,nstep_expand,r_expand, &
          &            eps0, &
          &            rcut,ifbook,rcut_book,nstep_book, &
          &            ifcellindex, &
          &            iucorname,iutopname, &
          &            npoly,npolytyp,npoly_mole,npoly_atom, &
          &            nwater, &
          &            nmatom,nmatyp,nmatomtyp, &
          &            degfree_poly,degfree_water, &
          &            degfree_ma,degfree_all, &
          &            ifstarec,ifcreatecor, &
          &            ifrdaddtop, &
          &            ifoutpdb,nstep_pdbout, &
          &            tcont_poly,tcont_water,tcont_ma, &
          &            tcont_poly_ini,tcont_water_ini,tcont_ma_ini, &
          &            tcontinterval,outinterval, &
          &            pressinterval,heatfinterval, &
          &            ifcalpremole,ifcalpreatom, &
          &            ifnetqcorrp, &
          &            oatmtyp,hatmtyp,monoatmtyp, &
          &            compfact, &
          &            ifrattle,eps_rattle, &
          &            rcutmor,ifbookmor,rcut_bookmor,nstep_bookmor, &
          &            ifcellindex_mor, &
          &            rcutsh,ifbooksh,rcut_booksh,nstep_booksh, &
          &            ifcellindex_sh, &
          &            rcutrfhfo,ifbookrfhfo,rcut_bookrfhfo, &
          &            nstep_bookrfhfo, &
          &            ifcellindex_rfhfo, &
          &            rcutrfhoo,ifbookrfhoo,rcut_bookrfhoo, &
          &            nstep_bookrfhoo, &
          &            ifcellindex_rfhoo, &
          &            rcutrfhoh,ifbookrfhoh,rcut_bookrfhoh, &
          &            nstep_bookrfhoh, &
          &            ifcellindex_rfhoh, &
          &            rcutdouo,rcutindouo,ifbookdouo,rcut_bookdouo, &
          &            nstep_bookdouo, &
          &            ifcellindex_douo, &
          &            rcutdouh,ifbookdouh,rcut_bookdouh, &
          &            nstep_bookdouh, &
          &            ifcellindex_douh, &
          &            rcutrpvw,ifbookrpvw,rcut_bookrpvw, &
          &            nstep_bookrpvw, &
          &            ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
          &            mchain,tfreq,text, &
          &            vfreq,pext,ifpatmcont,ifpmolcont, &
          &            pcont_axis, &
          &            next,nyosh, &
          &            nstep_vir, &
          &            iflimitmove,limitdist, &
          &            ifcalljlong,solvetyp,nsolve, &
          &            d_rini,d_rmax,d_econv,d_rmsf)

       ! AUGUMENT:
       !     INPUT
       integer,intent(in):: ousum           ! output parameter summarization file unit

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

       logical,intent(in):: iffennell       ! Fennell flag

       logical,intent(in):: ifljari         ! arithmetic mean for LJ cross parameter
       logical,intent(in):: ifljgeo         ! geometric mean for LJ cross parameter

       logical,intent(in):: iflocalheat     ! local heating flag
       logical,intent(in):: ifregionheat    ! region heating flag
       logical,intent(in):: ifregionhf      ! region heat flux control flag
       logical,intent(in):: ifreglange      ! region Langevin thermostat flag
       logical,intent(in):: iftcratom       ! region temp. control based on atom or mole.
       logical,intent(in):: ifoutthc        ! flag for outputting thermal control file

       logical,intent(in):: iflocalfix      ! fix atoms flag
       logical,intent(in):: iflocalfixz     ! flag for fixing z coordinate of atoms
       logical,intent(in):: iflocalfixzg    ! flag for fixing z coordinate of COM of molecules
       logical,intent(in):: ifcnp           ! flag to control normal pressure
       logical,intent(in):: ifposres        ! position restraint flag
       logical,intent(in):: ifpotbias       ! bias potential flag
       logical,intent(in):: iflocalvel      ! flag to force local atomic velocity
       logical,intent(in):: ifstrmvel       ! flag to input and use streaming velocity

       real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]
       real(8),intent(in):: dt_med_cal       ! time step of medium force [non-d]
       real(8),intent(in):: dt_short_cal     ! time step of short force [non-d]
       integer,intent(in):: nstep_med       ! number of step for medium force
       integer,intent(in):: nstep_short     ! number of step for short force

       integer,intent(in):: maxnstep        ! maximum step of MD (0-maxnstep)
       integer,intent(in):: nstage          ! stage number of MD
       integer,intent(in):: nstep_stage(:)  ! time step of each stage

       integer,intent(in):: mdcont_stage(:) ! MD control parameter of each stage

       real(8),intent(in):: eps0             ! dielectric constant of vacuum [non-d]
       real(8),intent(in):: rcut             ! vdw cutoff length [non-d]

       logical,intent(in):: ifcellindex     ! flag for cell index

       integer,intent(in):: nstep_book      ! bookkeeping interval
       real(8),intent(in):: rcut_book        ! cut off radius of bookkeeping
       logical,intent(in):: ifbook          ! flag for bookkeeping

       character(80),intent(in):: iucorname(:) ! input file name
       character(80),intent(in):: iutopname(:) !       "

       integer,intent(in):: npoly           ! all number of poly
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       integer,intent(in):: degfree_poly(:) ! degree of freedom of polymer1
       integer,intent(in):: degfree_water   ! degree of freedom of H2O
       integer,intent(in):: degfree_ma(:)   ! degree of freedom of monatomic mole.
       integer,intent(in):: degfree_all     ! degree of freedom of all atoms

       logical,intent(in):: ifstarec        ! read old state
       logical,intent(in):: ifcreatecor     ! create new coordinate

       logical,intent(in):: ifrdaddtop      ! input additional topology information

       logical,intent(in):: ifoutpdb        ! flag for outputting PDB format file
       integer,intent(in):: nstep_pdbout    ! MD step for output of PDB file

       real(8),intent(in):: tcont_poly(:)    ! poly Temp. [K] in NVT
       real(8),intent(in):: tcont_water      ! H2O Temp. [K] in NVT
       real(8),intent(in):: tcont_ma(:)      ! monatomic mole. Temp. [K] in NVT
       real(8),intent(in):: tcont_poly_ini ! poly Temp. [K] in NVT (Woodcock) for md_h
       real(8),intent(in):: tcont_water_ini ! H2O Temp. [K] in NVT (Woodcock) for md_h
       real(8),intent(in):: tcont_ma_ini     ! MA Temp. [K] in NVT (Woodcock) for md_h

       integer,intent(in):: tcontinterval   ! interval of temp. control
       integer,intent(in):: outinterval     ! interval of outputting data
       integer,intent(in):: pressinterval   ! interval of pressure calculation
       integer,intent(in):: heatfinterval   ! interval of heatf output

       logical,intent(in):: ifcalpremole    ! pressure calculation of molecule
       logical,intent(in):: ifcalpreatom    ! pressure calculation of atom
       logical,intent(in):: ifnetqcorrp     ! net charge correction for pressure

       character(2),intent(in):: oatmtyp ! O atomtype of water model
       character(2),intent(in):: hatmtyp ! H atomtype of water model
       character(2),intent(in):: monoatmtyp(:) ! monatomic mole. type

       real(8),intent(in):: compfact      ! compact factor using at poly arrange(<1.0)

       logical,intent(in):: ifrattle        ! rattle flag
       real(8),intent(in):: eps_rattle     ! tolerance (relative difference) of rattle

       real(8),intent(in):: rcutmor          ! Morse cutoff length [non-d]

       logical,intent(in):: ifcellindex_mor ! flag for cell index (morse)

       logical,intent(in):: ifbookmor      ! flag for bookkeeping of Morse interaction
       real(8),intent(in):: rcut_bookmor
                         ! cut off radius of bookkeeping[non-d] of Morse
       integer,intent(in):: nstep_bookmor  ! bookkeeping interval of Morse interaction

       real(8),intent(in):: rcutsh           ! Sh cutoff length [non-d]

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

       logical,intent(in):: ifbookrfhoo
                         ! flag for bookkeeping of RFH(OO) interaction
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

       logical,intent(in):: ifcstmnb          ! flag if using custom NB interaction
       logical,intent(in):: ifcellindex_cstmnb ! flag for cell index (custom NB)
       logical,intent(in):: ifbookcstmnb
                               ! flag for bookkeeping of custom NB interaction

       integer,intent(in):: mchain          ! Nose-Hoover chain number (>1)
       real(8),intent(in):: tfreq            ! temperature frequency [non-d]
       real(8),intent(in):: text          ! external temp. [non-d] (Nose-Hoover chain)

       real(8),intent(in):: vfreq            ! volume change frequency [non-d]
       real(8),intent(in):: pext             ! external pressure [Pa->non-d]
       logical,intent(in):: ifpatmcont      ! atomic pressure control
       logical,intent(in):: ifpmolcont      ! molecular pressure control
       character(5),intent(in):: pcont_axis  ! axis for pressure control (iso, aniso, etc.)

       integer,intent(in):: next            ! iteration number of extended system
       integer,intent(in):: nyosh           ! expansion order of Yoshida-Suzuki method

       integer,intent(in):: nstep_vir    ! number of step for virtual time integration
       logical,intent(in):: iflimitmove ! if limit atomic motion to a specific distance
       real(8),intent(in):: limitdist  ! maximum atomic displacement when doing
                               ! time integration [non-d] (structure relaxation)

       logical,intent(in):: ifcalljlong     ! long-range correction in pressure
       character(2),intent(in):: solvetyp ! solvent molecule
       integer,intent(in):: nsolve          ! number of solvent molecules

       integer,intent(in):: nstep_maxwell   ! time step of maxwell distribution
       integer,intent(in):: nstep_expand    ! time step of cell expansion
       real(8),intent(in):: r_expand         ! expansion ratio of cell (volume)

       real(8),intent(inout):: d_rini   ! initial displacement dr for EM [non-d]
       real(8),intent(inout):: d_rmax   ! maximum displacement dr for EM [non-d]

       real(8),intent(inout):: d_econv  ! convergence condition for energy
                                   !    in EM [non-d]
       real(8),intent(inout):: d_rmsf   ! convergence condition for
                                   !    root mean square force in EM [non-d]

     end subroutine wrsumm

  end interface

end module interface_preproc
