!*****************************
!*  wrsumm.f90 Ver.6.7       *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2018-02-23 17:51:51 gota>

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

  use md_global

  implicit none

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
  logical,intent(in):: ifpatmcont       ! atomic pressure control
  logical,intent(in):: ifpmolcont       ! molecular pressure control
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

!     LOCAL
  integer:: i
  integer:: degfree_polyall
  integer:: degfree_maall

!     +     +     +     +     +     +

!---- some preparation
  degfree_polyall = 0
  degfree_maall = 0
  do i = 1, npolytyp
     degfree_polyall = degfree_polyall + degfree_poly(i)
  end do
  do i = 1, nmatyp
     degfree_maall = degfree_maall + degfree_ma(i)
  end do

!---- Header
  write(ousum,'(A)') '***** parameter summarization *****'

!---- write base value
  write(ousum,'(A7,E13.5)') 'xref',xref
  write(ousum,'(A7,E13.5)') 'eref',eref
  write(ousum,'(A7,E13.5)') 'mref',mref
  write(ousum,'(A7,E13.5)') 'qref',qref
  write(ousum,'(A7,E13.5)') 'vref',vref
  write(ousum,'(A7,E13.5)') 'timeref',timeref
  write(ousum,'(A7,E13.5)') 'tempref',tempref
  write(ousum,'(A7,E13.5)') 'pref',pref
  write(ousum,'(A7,E13.5)') 'fref',fref
  write(ousum,'(A7,E13.5)') 'eps0ref',eps0ref
  write(ousum,*)

!---- write cell box length
  write(ousum,'(A4,E16.8)') 'xcel',xcel*xref
  write(ousum,'(A4,E16.8)') 'ycel',ycel*xref
  write(ousum,'(A4,E16.8)') 'zcel',zcel*xref
  write(ousum,*)

!---- write parameter for ewald sum
  write(ousum,'(A7,L13)') 'ifewald',ifewald
  write(ousum,'(A7,E16.8)') 'alpha',alpha/xref
  write(ousum,'(A7,I13)') 'kmax',kmax
  write(ousum,'(A7,I13)') 'nwave',nwave
  write(ousum,'(A7,E16.8)') 'rrcut',rrcut*xref
  write(ousum,*)

!---- write parameter for SPME method
  write(ousum,'(A9,L8)') 'ifspme',ifspme
  write(ousum,'(A9,I8)') 'nfft1',nfft1
  write(ousum,'(A9,I8)') 'nfft2',nfft2
  write(ousum,'(A9,I8)') 'nfft3',nfft3
  write(ousum,'(A9,I8)') 'pme_order',pme_order
  write(ousum,*)

!---- write parameter for Fennell method
  write(ousum,'(A9,L8)') 'iffennell',iffennell
  write(ousum,*)

!---- write parameter for LJ parameter
  write(ousum,'(A7,L8)') 'ifljari',ifljari
  write(ousum,'(A7,L8)') 'ifljgeo',ifljgeo
  write(ousum,*)

!---- write parameter for specific control
  write(ousum,'(A12,L8)') 'iflocalheat',iflocalheat
  write(ousum,'(A12,L8)') 'ifregionheat',ifregionheat
  write(ousum,'(A12,L8)') 'ifregionhf',ifregionhf
  write(ousum,'(A12,L8)') 'ifreglange',ifreglange
  write(ousum,'(A12,L8)') 'iftcratom',iftcratom
  write(ousum,'(A12,L8)') 'ifoutthc',ifoutthc
  write(ousum,'(A12,L8)') 'iflocalfix',iflocalfix
  write(ousum,'(A12,L8)') 'iflocalfixz',iflocalfixz
  write(ousum,'(A12,L8)') 'iflocalfixzg',iflocalfixzg
  write(ousum,'(A12,L8)') 'ifcnp',ifcnp
  write(ousum,'(A12,L8)') 'ifposres',ifposres
  write(ousum,'(A12,L8)') 'ifpotbias',ifpotbias
  write(ousum,'(A12,L8)') 'iflocalvel',iflocalvel
  write(ousum,'(A12,L8)') 'ifstrmvel',ifstrmvel
  write(ousum,*)

!---- write time step
  write(ousum,'(A11,E13.5)') 'dt_long',dt_long_cal*timeref
  write(ousum,'(A11,E13.5)') 'dt_med',dt_med_cal*timeref
  write(ousum,'(A11,E13.5)') 'dt_short',dt_short_cal*timeref
  write(ousum,'(A11,I13)') 'nstep_med',nstep_med
  write(ousum,'(A11,I13)') 'nstep_short',nstep_short
  write(ousum,*)

!---- write maximum step
  write(ousum,'(A12,I8)') 'maxnstep',maxnstep
  write(ousum,'(A12,I8)') 'nstage',nstage
  write(ousum,*)
  do i = 1, nstage
     write(ousum,'(A12,I3,I8)') 'nstep_stage',i,nstep_stage(i)
  end do
  do i = 1, nstage
     if (mdcont_stage(i) == md_0k) then
        write(ousum,'(A12,I3,A8)') 'mdcont_stage',i,'md_0k'
     else if (mdcont_stage(i) == md_h) then
        write(ousum,'(A12,I3,A8)') 'mdcont_stage',i,'md_h'
     else if (mdcont_stage(i) == md_t) then
        write(ousum,'(A12,I3,A8)') 'mdcont_stage',i,'md_t'
     else if (mdcont_stage(i) == md_mtk) then
        write(ousum,'(A12,I3,A8)') 'mdcont_stage',i,'md_mtk'
     else if (mdcont_stage(i) == md_nhc) then
        write(ousum,'(A12,I3,A8)') 'mdcont_stage',i,'md_nhc'
     else if (mdcont_stage(i) == md_nve) then
        write(ousum,'(A12,I3,A8)') 'mdcont_stage',i,'md_nve'
     else if (mdcont_stage(i) == md_htf) then
        write(ousum,'(A12,I3,A8)') 'mdcont_stage',i,'md_htf'
     end if
  end do
  write(ousum,*)

!---- write time step of maxwell distribution
  write(ousum,'(A12,I8)') 'nstep_maxwell',nstep_maxwell
!---- write time step of expansion
  write(ousum,'(A12,I8)') 'nstep_expand',nstep_expand
  write(ousum,'(A12,E16.8)') 'r_expand',r_expand
  write(ousum,*)

!---- write parameter for nonbonded interaction
  write(ousum,'(A16,E13.5)') 'eps0',eps0*eps0ref
  write(ousum,'(A16,E13.5)') 'rcut',rcut*xref
  write(ousum,*)

!---- write bookkeeping parameters
  write(ousum,'(A11,L13)') 'ifcellindex',ifcellindex
  write(ousum,'(A11,L13)') 'ifbook',ifbook
  write(ousum,'(A11,I13)') 'nstep_book',nstep_book
  write(ousum,'(A11,E13.5)') 'rcut_book',rcut_book*xref
  write(ousum,*)

!---- write bookkeeping parameters for Morse interaction
  write(ousum,'(A15,E13.5)') 'rcutmor',rcutmor*xref
  write(ousum,*)

  write(ousum,'(A15,L13)') 'ifcellindex_mor',ifcellindex_mor
  write(ousum,'(A15,L13)') 'ifbookmor',ifbookmor
  write(ousum,'(A15,E13.5)') 'rcut_bookmor',rcut_bookmor*xref
  write(ousum,'(A15,I13)') 'nstep_bookmor',nstep_bookmor
  write(ousum,*)

!---- write bookkeeping parameters for SH interaction
  write(ousum,'(A15,E13.5)') 'rcutsh',rcutsh*xref
  write(ousum,*)

  write(ousum,'(A15,L13)') 'ifcellindex_sh',ifcellindex_sh
  write(ousum,'(A15,L13)') 'ifbooksh',ifbooksh
  write(ousum,'(A15,E13.5)') 'rcut_booksh',rcut_booksh*xref
  write(ousum,'(A15,I13)') 'nstep_booksh',nstep_booksh
  write(ousum,*)

!---- write bookkeeping parameters for RFH interaction
  write(ousum,'(A17,E13.5)') 'rcutrfhfo',rcutrfhfo*xref
  write(ousum,*)

  write(ousum,'(A17,L13)') 'ifcellindex_rfhfo',ifcellindex_rfhfo
  write(ousum,'(A17,L13)') 'ifbookrfhfo',ifbookrfhfo
  write(ousum,'(A17,E13.5)') 'rcut_bookrfhfo',rcut_bookrfhfo*xref
  write(ousum,'(A17,I13)') 'nstep_bookrfhfo',nstep_bookrfhfo
  write(ousum,*)

  write(ousum,'(A17,E13.5)') 'rcutrfhoo',rcutrfhoo*xref
  write(ousum,*)

  write(ousum,'(A17,L13)') 'ifcellindex_rfhoo',ifcellindex_rfhoo
  write(ousum,'(A17,L13)') 'ifbookrfhoo',ifbookrfhoo
  write(ousum,'(A17,E13.5)') 'rcut_bookrfhoo',rcut_bookrfhoo*xref
  write(ousum,'(A17,I13)') 'nstep_bookrfhoo',nstep_bookrfhoo
  write(ousum,*)

  write(ousum,'(A17,E13.5)') 'rcutrfhoh',rcutrfhoh*xref
  write(ousum,*)

  write(ousum,'(A17,L13)') 'ifcellindex_rfhoh',ifcellindex_rfhoh
  write(ousum,'(A17,L13)') 'ifbookrfhoh',ifbookrfhoh
  write(ousum,'(A17,E13.5)') 'rcut_bookrfhoh',rcut_bookrfhoh*xref
  write(ousum,'(A17,I13)') 'nstep_bookrfhoh',nstep_bookrfhoh
  write(ousum,*)

!---- write bookkeeping parameters for DOU interaction
  write(ousum,'(A16,E13.5)') 'rcutdouo',rcutdouo*xref
  write(ousum,'(A16,E13.5)') 'rcutindouo',rcutindouo*xref
  write(ousum,*)

  write(ousum,'(A16,L13)') 'ifcellindex_douo',ifcellindex_douo
  write(ousum,'(A16,L13)') 'ifbookdouo',ifbookdouo
  write(ousum,'(A16,E13.5)') 'rcut_bookdouo',rcut_bookdouo*xref
  write(ousum,'(A16,I13)') 'nstep_bookdouo',nstep_bookdouo
  write(ousum,*)

  write(ousum,'(A16,E13.5)') 'rcutdouh',rcutdouh*xref
  write(ousum,*)

  write(ousum,'(A16,L13)') 'ifcellindex_douh',ifcellindex_douh
  write(ousum,'(A16,L13)') 'ifbookdouh',ifbookdouh
  write(ousum,'(A16,E13.5)') 'rcut_bookdouh',rcut_bookdouh*xref
  write(ousum,'(A16,I13)') 'nstep_bookdouh',nstep_bookdouh
  write(ousum,*)

!---- write bookkeeping parameters for RP-VW interaction
  write(ousum,'(A14,E13.5)') 'rcutrpvw',rcutrpvw*xref
  write(ousum,*)

  write(ousum,'(A14,L13)') 'ifbookrpvw',ifbookrpvw
  write(ousum,'(A14,E13.5)') 'rcut_bookrpvw',rcut_bookrpvw*xref
  write(ousum,'(A14,I13)') 'nstep_bookrpvw',nstep_bookrpvw
  write(ousum,*)

!---- write parameters for custom NB interaction
  write(ousum,'(A18,L8)') 'ifcstmnb',ifcstmnb
  write(ousum,'(A18,L8)') 'ifcellindex_cstmnb',ifcellindex_cstmnb
  write(ousum,'(A18,L8)') 'ifbookcstmnb',ifbookcstmnb
  write(ousum,*)

!---- write number of atoms or molucules
  write(ousum,'(A11,I8)') 'natom',natom
  write(ousum,'(A11,I8)') 'nmole',nmole
  write(ousum,'(A11,I8)') 'npoly',npoly
  write(ousum,'(A11,I8)') 'npolytyp',npolytyp
  do i=1,npolytyp
     write(ousum,'(A11,I2,I8,I8,1X,A30,A30)') 'npolymoletyp',i, &
          & npoly_mole(i),npoly_atom(i), &
          & iucorname(i),iutopname(i)
  end do
  write(ousum,'(A11,I8)') 'nwater',nwater
  write(ousum,'(A11,I8)') 'nwater_all',nwater*3
  write(ousum,'(A11,I8)') 'nmatom',nmatom
  write(ousum,'(A11,I8)') 'nmatyp',nmatyp
  do i = 1, nmatyp
     write(ousum,'(A11,I3,I8,A3)') 'nmatomtyp',i,nmatomtyp(i), &
          &                                      monoatmtyp(i)
  end do

  write(ousum,*)

!---- write degree of freedom
  do i=1,npolytyp
     write(ousum,'(A16,I5,I8)') 'degfree_poly',i,degfree_poly(i)
  end do
  write(ousum,'(A16,I5,I8)') 'degfree_water',1,degfree_water
  do i = 1, nmatyp
     write(ousum,'(A16,I5,I8)') 'degfree_ma',i,degfree_ma(i)
  end do
  write(ousum,*)
  write(ousum,'(A16,I8)') 'degfree_polyall',degfree_polyall
  write(ousum,'(A16,I8)') 'degfree_water',degfree_water
  write(ousum,'(A16,I8)') 'degfree_maall',degfree_maall
  write(ousum,'(A16,I8)') 'degfree_all',degfree_all
  write(ousum,*)

!---- write control for creation of coordinate
  write(ousum,'(A11,L2)') 'ifstarec',ifstarec
  write(ousum,'(A11,L2)') 'ifcreatecor',ifcreatecor
  write(ousum,*)

!---- write control for inputting additional topology information
  write(ousum,'(A10,L2)') 'ifrdaddtop',ifrdaddtop
  write(ousum,*)

!---- write control for PDB file output
  write(ousum,'(A12,L2)') 'ifoutpdb',ifoutpdb
  write(ousum,'(A12,I8)') 'nstep_pdbout',nstep_pdbout
  write(ousum,*)

!---- write control temperature
  write(ousum,'(A14,F9.2)') 'temp_poly',tcont_poly(1)*tempref
  write(ousum,'(A14,F9.2)') 'temp_water',tcont_water*tempref
  write(ousum,'(A14,F9.2)') 'temp_ma',tcont_ma(1)*tempref
  write(ousum,'(A14,F9.2)') 'temp_poly_ini',tcont_poly_ini*tempref
  write(ousum,'(A14,F9.2)') 'temp_water_ini',tcont_water_ini*tempref
  write(ousum,'(A14,F9.2)') 'temp_ma_ini',tcont_ma_ini*tempref
  write(ousum,*)

!---- write temp. control interval
  write(ousum,'(A13,I5)') 'tcontinterval',tcontinterval
  write(ousum,*)

!---- write data output interval
  write(ousum,'(A13,I5)') 'outinterval',outinterval
  write(ousum,'(A13,I5)') 'pressinterval',pressinterval
  write(ousum,'(A13,I5)') 'heatfinterval',heatfinterval
  write(ousum,*)

!---- write if pressure calculation or not
  write(ousum,'(A12,L9)') 'ifcalpremole',ifcalpremole
  write(ousum,'(A12,L9)') 'ifcalpreatom',ifcalpreatom
  write(ousum,'(A12,L9)') 'ifnetqcorrp',ifnetqcorrp
  write(ousum,*)

!---- write water model
  write(ousum,'(A9,A3)') 'Owatertyp',oatmtyp
  write(ousum,'(A9,A3)') 'Hwatertyp',hatmtyp
  write(ousum,*)

!---- write compact factor at createcor
  write(ousum,'(A8,F7.3)') 'compfact',compfact
  write(ousum,*)

!---- write if applying rattle method or not
  write(ousum,'(A8,L9)') 'ifrattle',ifrattle
  write(ousum,'(A10,E13.5)') 'eps_rattle',eps_rattle
  write(ousum,*)

!---- write parameter of Nose-Hoover chain
  write(ousum,'(A6,I16)') 'mchain',mchain
  write(ousum,'(A6,E16.8)') 'tfreq',tfreq/timeref
  write(ousum,'(A6,F16.3)') 'text',text*tempref
  write(ousum,*)

!---- write parameter of Andersen (Hoover type) barostat
  write(ousum,'(A6,E16.8)') 'vfreq',vfreq/timeref
  write(ousum,'(A6,E16.8)') 'pext',pext*pref
  write(ousum,'(A11,L2)') 'ifpatmcont',ifpatmcont
  write(ousum,'(A11,L2)') 'ifpmolcont',ifpmolcont
  write(ousum,'(A11,A5)') 'pcont_axis',pcont_axis
  write(ousum,*)

!---- write parameter of Yoshida-Suzuki method
  write(ousum,'(A5,I4)') 'next',next
  write(ousum,'(A5,I4)') 'nyosh',nyosh
  write(ousum,*)

! !---- write parameter of virtual time integration
!   write(ousum,'(A9,I4)') 'nstep_vir',nstep_vir
!   write(ousum,*)

!---- write parameter of limit for atomic motion
  write(ousum,'(A11,L2)') 'iflimitmove',iflimitmove
  write(ousum,'(A11,E16.8)') 'limitdist',limitdist * xref
  write(ousum,*)

!---- write parameter of energy minimization

  write(ousum,'(A7,E16.8)') 'd_rini',d_rini * xref
  write(ousum,'(A7,E16.8)') 'd_rmax',d_rmax * xref
  write(ousum,'(A7,E16.8)') 'd_econv',d_econv * eref
  write(ousum,'(A7,E16.8)') 'd_rmsf',d_rmsf * fref
  write(ousum,*)

!---- write parameter of L-J long-range correction
  write(ousum,'(A11,L8)') 'ifcalljlong',ifcalljlong
  write(ousum,'(A11,A8)') 'solvetyp',solvetyp
  write(ousum,'(A11,I8)') 'nsolve',nsolve
  write(ousum,*)

!---- write parameter of spline interpolation
  write(ousum,'(A7,I8)') 'nspltbl',nspltbl
  write(ousum,*)

!---- close file
  close(ousum)

!     +     +     +     +     +     +

  return
end subroutine wrsumm
