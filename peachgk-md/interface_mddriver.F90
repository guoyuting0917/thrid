!************************************
!*  interface_mddriver.f90 Ver.2.2  *
!*      for peachgk_md.f            *
!*            by G.Kikugawa         *
!************************************
! Time-stamp: <2015-03-30 19:48:03 gota>

!***** This module is interface module for MD driver routines *****

module interface_mddriver

  interface

     subroutine moldyn(iostarec, recinterval, &
          &            npoly, npolytyp, npoly_mole, npoly_atom, &
          &            nwater, &
          &            nmatom, nmatyp, nmatomtyp, &
          &            maxnstep, inistep, endstep, &
          &            xcel,ycel,zcel, &
          &            ifcenterfix_all, &
          &            ifcenterfix_poly, &
          &            ifcenterfix_water, &
          &            ifcenterfix_ma, &
          &            cenfix_free, &
          &            ifcenterfix_polytyp, &
          &            ifcenterfix_watertyp, &
          &            ifcenterfix_matyp, &
          &            dt_long_cal, &
          &            nstep_med, nstep_short, &
          &            nstep_vir, &
          &            iflimitmove,limitdist, &
          &            mts_bond, mts_angl, mts_anglub, &
          &            mts_tors, mts_torsrb, mts_torsim, &
          &            mts_vdw, mts_ewr, mts_ewk, mts_vdw14, mts_elc14, &
          &            mts_mor, mts_sh, mts_rfh, mts_dou, mts_cnpvw, &
          &            mts_cstmnb, &
          &            mts_posres, &
          &            mts_potbias, &
          &            ifewald, alpha, kmax,rrcut, &
          &            ifspme, nfft1, nfft2, nfft3, pme_order, &
          &            eps0, div_factor_14vdw, div_factor_14elc, &
          &            rcut, &
          &            ifcellindex, &
          &            ifbook, nstep_book, &
          &            tcont_poly, tcont_water, tcont_ma, &
          &            tcont_poly_ini, tcont_water_ini, tcont_ma_ini, &
          &            ifrattle, eps_rattle, &
          &            iflocalheat, &
          &            nlheat_poly, index_nlheat_poly, &
          &            tcont_nlheat_poly, &
          &            nlheat_water, tcont_nlheat_water, &
          &            nlheat_ma,index_nlheat_ma, &
          &            tcont_nlheat_ma, &
          &            ifregionheat, ifregionhf, ifreglange, &
          &            iftcratom, ifoutthc, &
          &            ntcregion, &
          &            tcxpos1, tcxpos2, &
          &            tcypos1, tcypos2, &
          &            tczpos1, tczpos2, &
          &            r_tcont, &
          &            nhfcregion, &
          &            hfczpos1, hfczpos2, &
          &            r_hfcont, &
          &            nlangeregion, &
          &            ltxpos1, ltxpos2, &
          &            ltypos1, ltypos2, &
          &            ltzpos1, ltzpos2, &
          &            r_ltemp, r_ltdamp, &
          &            iflocalfix, &
          &            nlfix, index_nlfix, &
          &            iflocalfixz, iflocalfixzg, &
          &            ifcnp, &
          &            ifposres, &
          &            ifpotbias, &
          &            iflocalvel, ifstrmvel, &
          &            nlvel, index_nlvel, v_nlvel, &
          &            ifoutpdb, nstep_pdbout, &
          &            resname_poly_pdb, &
          &            resname_water_pdb, &
          &            resname_matom_pdb, &
          &            xref, eref, mref, qref, &
          &            vref, timeref, tempref, pref, fref, eps0ref, &
          &            degfree_poly, degfree_water, &
          &            degfree_ma, degfree_all, &
          &            dt_med_cal, dt_short_cal, &
          &            rcut_book, &
          &            rcutmor, ifbookmor, rcut_bookmor, nstep_bookmor, &
          &            ifcellindex_mor, &
          &            rcutsh, ifbooksh, rcut_booksh, nstep_booksh, &
          &            ifcellindex_sh, &
          &            rcutrfhfo, ifbookrfhfo, rcut_bookrfhfo, &
          &            nstep_bookrfhfo, &
          &            ifcellindex_rfhfo, &
          &            rcutrfhoo, ifbookrfhoo, rcut_bookrfhoo, &
          &            nstep_bookrfhoo, &
          &            ifcellindex_rfhoo, &
          &            rcutrfhoh, ifbookrfhoh, rcut_bookrfhoh, &
          &            nstep_bookrfhoh, &
          &            ifcellindex_rfhoh, &
          &            rcutdouo, rcutindouo, ifbookdouo, rcut_bookdouo, &
          &            nstep_bookdouo, &
          &            ifcellindex_douo, &
          &            rcutdouh, ifbookdouh, rcut_bookdouh, &
          &            nstep_bookdouh, &
          &            ifcellindex_douh, &
          &            rcutrpvw, ifbookrpvw, rcut_bookrpvw, &
          &            nstep_bookrpvw, &
          &            ifcstmnb,ifcellindex_cstmnb,ifbookcstmnb, &
          &            for_long, for_short, &
          &            for_viric_long, for_viric_med, for_viric_short, &
          &            for_virilj_long, for_virilj_med, for_virilj_short, &
          &            for_virimor_long, for_virimor_med, &
          &            for_virimor_short, &
          &            for_virish_long, for_virish_med, &
          &            for_virish_short, &
          &            for_virirfh_long, for_virirfh_med, &
          &            for_virirfh_short, &
          &            for_viridou_long, for_viridou_med, &
          &            for_viridou_short, &
          &            for_viricstmnb_long, for_viricstmnb_med, &
          &            for_viricstmnb_short, &
          &            pot_ewc, &
          &            ouene,oupos,ouvel,oufor,outhe,oubar,oupre,outhc,oupdb, &
          &            ouumb, &
          &            ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
          &            ifoutbar,ifoutpre, &
          &            tcontinterval, outinterval, pressinterval, &
          &            yratio, zratio, &
          &            ifcalpremole, ifcalpreatom, &
          &            ifnetqcorrp, &
          &            mchain, &
          &            pint, pintt, &
          &            ifpatmcont, ifpmolcont, &
          &            ifcalljlong, nsolve, solveindex, &
          &            netchrgsq, &
          &            nstep_expand, r_expand, &
          &            nstep_maxwell, &
          &            md_cont)

       ! ARGUMENT:
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
       real(8),intent(in):: tcont_poly_ini ! poly Temp. [K] in NVT (Woodcock) for md_h
       real(8),intent(in):: tcont_water_ini ! H2O Temp. [K] in NVT (Woodcock) for md_h
       real(8),intent(in):: tcont_ma_ini     ! MA Temp. [K] in NVT (Woodcock) for md_h

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

       !---- variable for local heating
       logical,intent(in):: iflocalheat     ! local heating flag

       integer,intent(in):: nlheat_poly     ! number of poly type for local heating
       integer,intent(in):: index_nlheat_poly(:)
                                ! index of poly type for local heating
       real(8),intent(in):: tcont_nlheat_poly(:)
                                ! control temp. of poly type for local heating
       integer,intent(in):: nlheat_water    ! number of water for local heating
       real(8),intent(in):: tcont_nlheat_water
                                ! control temp. of water for local heating
       integer,intent(in):: nlheat_ma       ! number of matom type for local heating
       integer,intent(in):: index_nlheat_ma(:)
                                ! index of matom type for local heating
       real(8),intent(in):: tcont_nlheat_ma(:)
                                ! control temp. of matom type for local heating

       !---- variables for region-based temp. control
       logical,intent(in):: ifregionheat    ! region heating flag
       logical,intent(in):: iftcratom    ! region temp. control based on atom or mole.
       logical,intent(in):: ifoutthc     ! flag for outputting thermal control file

       integer,intent(inout):: ntcregion       ! number of region to control temp.
       real(8),intent(in):: tcxpos1(:),tcxpos2(:)
                                ! x-position of temp. control region
       real(8),intent(in):: tcypos1(:),tcypos2(:)
                                ! y-position of temp. control region
       real(8),intent(in):: tczpos1(:),tczpos2(:)
                                ! z-position of temp. control region
       real(8),intent(in):: r_tcont(:)       ! control temp. in each region

       !---- variables for region-based heat flux control
       logical,intent(in):: ifregionhf      ! region heat flux control flag
       integer,intent(inout):: nhfcregion      ! number of region to control heat flux
       real(8),intent(in):: hfczpos1(:),hfczpos2(:)
                                ! z-position of heat flux control region
       real(8),intent(in):: r_hfcont(:)      ! magnitude of heat flux in each region
                                        ! (converted to the input energy)

       !---- variables for region-based Langevin thermostat
       logical,intent(in):: ifreglange      ! region Langevin thermostat flag
       integer,intent(in):: nlangeregion    ! number of region for Langevin thermo.
       real(8),intent(in):: ltxpos1(:),ltxpos2(:)
                                       ! x-position of temp. control region
       real(8),intent(in):: ltypos1(:),ltypos2(:)
                                       ! y-position of temp. control region
       real(8),intent(in):: ltzpos1(:),ltzpos2(:)
                                       ! z-position of temp. control region
       real(8),intent(in):: r_ltemp(:)      ! control temp. in each region
       real(8),intent(in):: r_ltdamp(:)     ! damping factor in each region [non-d]

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
       integer,intent(in):: oupre           ! output unit for outpre velocity data
       integer,intent(in):: outhc           ! output unit for outthc thermal control data
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
       !  real(8),intent(inout):: for_med(:,:)              ! medium-range force
       real(8),intent(inout):: for_short(:,:) ! short-range force

       !---- potential valiables
       real(8),intent(in):: pot_ewc ! coulomb potential(ewald self)

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

       real(8),intent(inout):: for_virilj_long(:,:)
                            ! long-range virial (L-J force)
       real(8),intent(inout):: for_virilj_med(:,:)
                            ! medium-range virial (L-J force)
       real(8),intent(inout):: for_virilj_short(:,:)
                            ! short-range virial (L-J force)

       real(8),intent(inout):: for_virimor_long(:,:)
                            ! long-range virial (Morse force)
       real(8),intent(inout):: for_virimor_med(:,:)
                            ! med-range virial (Morse force)
       real(8),intent(inout):: for_virimor_short(:,:)
                            ! short-range virial (Morse force)

       real(8),intent(inout):: for_virish_long(:,:)
                            ! long-range virial (SH force)
       real(8),intent(inout):: for_virish_med(:,:)
                            ! med-range virial (SH force)
       real(8),intent(inout):: for_virish_short(:,:)
                            ! short-range virial (SH force)

       real(8),intent(inout):: for_virirfh_long(:,:)
                            ! long-range virial (RFH force)
       real(8),intent(inout):: for_virirfh_med(:,:)
                            ! med-range virial (RFH force)
       real(8),intent(inout):: for_virirfh_short(:,:)
                            ! short-range virial (RFH force)

       real(8),intent(inout):: for_viridou_long(:,:)
                            ! long-range virial (DOU force)
       real(8),intent(inout):: for_viridou_med(:,:)
                            ! med-range virial (DOU force)
       real(8),intent(inout):: for_viridou_short(:,:)
                            ! short-range virial (DOU force)

       real(8),intent(inout):: for_viricstmnb_long(:,:)
                            ! long-range virial (custom NB force)
       real(8),intent(inout):: for_viricstmnb_med(:,:)
                            ! med-range virial (custom NB force)
       real(8),intent(inout):: for_viricstmnb_short(:,:)
                            ! short-range virial (custom NB force)

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
       integer,intent(in):: solveindex      ! atmindex of solvent atom

       !---- variables for net charge calculation
       real(8),intent(inout):: netchrgsq        ! = (sum(qi))**2

       !---- cell expansion
       integer,intent(in):: nstep_expand    ! time step of cell expansion
       real(8),intent(in):: r_expand         ! expansion ratio of cell (volume)

       !---- time step of forcing Maxwell distribution
       integer,intent(in):: nstep_maxwell   ! time step of maxwell distribution

       !---- MD control parameter
       integer,intent(in):: md_cont         ! MD control flag

     end subroutine moldyn

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

       ! ARGUMENT:
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
       real(8),intent(in):: pot_ewc ! coulomb potential(ewald self)

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

       real(8),intent(inout):: for_virilj_long(:,:)
                            ! long-range virial (L-J force)
       real(8),intent(inout):: for_virilj_med(:,:)
                            ! medium-range virial (L-J force)
       real(8),intent(inout):: for_virilj_short(:,:)
                            ! short-range virial (L-J force)

       real(8),intent(inout):: for_virimor_long(:,:)
                            ! long-range virial (Morse force)
       real(8),intent(inout):: for_virimor_med(:,:)
                            ! med-range virial (Morse force)
       real(8),intent(inout):: for_virimor_short(:,:)
                            ! short-range virial (Morse force)

       real(8),intent(inout):: for_virish_long(:,:)
                            ! long-range virial (SH force)
       real(8),intent(inout):: for_virish_med(:,:)
                            ! med-range virial (SH force)
       real(8),intent(inout):: for_virish_short(:,:)
                            ! short-range virial (SH force)

       real(8),intent(inout):: for_virirfh_long(:,:)
                         ! long-range virial (RFH force)
       real(8),intent(inout):: for_virirfh_med(:,:)
                         ! med-range virial (RFH force)
       real(8),intent(inout):: for_virirfh_short(:,:)
                         ! short-range virial (RFH force)

       real(8),intent(inout):: for_viridou_long(:,:)
                         ! long-range virial (DOU force)
       real(8),intent(inout):: for_viridou_med(:,:)
                         ! med-range virial (DOU force)
       real(8),intent(inout):: for_viridou_short(:,:)
                         ! short-range virial (DOU force)

       real(8),intent(inout):: for_viricstmnb_long(:,:)
                            ! long-range virial (custom NB force)
       real(8),intent(inout):: for_viricstmnb_med(:,:)
                            ! med-range virial (custom NB force)
       real(8),intent(inout):: for_viricstmnb_short(:,:)
                            ! short-range virial (custom NB force)

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
       integer,intent(in):: solveindex      ! atmindex of solvent atom

       !---- variables for net charge calculation
       real(8),intent(inout):: netchrgsq        ! = (sum(qi))**2

       !---- cell expansion
       integer,intent(in):: nstep_expand    ! time step of cell expansion
       real(8),intent(in):: r_expand         ! expansion ratio of cell (volume)

       !---- MD control parameter
       integer,intent(in):: md_cont         ! MD control flag

     end subroutine moldyn_nhc

     subroutine moldyn_mtk(iostarec, recinterval, &
          &                npoly,npolytyp, npoly_mole, npoly_atom, &
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
          &                iflocalfix, &
          &                nlfix, index_nlfix, &
          &                iflocalfixz, iflocalfixzg, &
          &                ifposres, &
          &                ifpotbias, &
          &                iflocalvel, &
          &                nlvel, index_nlvel, v_nlvel, &
          &                ifoutpdb, nstep_pdbout, &
          &                resname_poly_pdb, &
          &                resname_water_pdb, &
          &                resname_matom_pdb, &
          &                xref,eref,mref,qref, &
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
          &                yratio,zratio, &
          &                ifcalpremole, ifcalpreatom, &
          &                ifnetqcorrp, &
          &                mchain,text, &
          &                pext, pint, pintt, &
          &                ifpatmcont, ifpmolcont, &
          &                pcont_axis, &
          &                next,nyosh, &
          &                ifcalljlong, nsolve, solveindex, &
          &                netchrgsq, &
          &                nstep_expand, r_expand, &
          &                md_cont)

       ! ARGUMENT:
       !---- i/o unit
       integer,intent(in):: iostarec        ! state record file unit

       !---- MD control parameters
       integer,intent(in):: npoly           ! number of polymer1
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

       real(8),intent(inout):: yratio       ! y cell ratio of y to x
       real(8),intent(inout):: zratio       ! z cell ratio of z to x

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
       real(8),intent(in):: tcont_ma(:)      ! monatomic mole. Temp. [non-d] in NVT

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

       logical,intent(in):: iflocalfixzg     ! flag for fixing z coordinate of COM of molecules

       !---- variable for position restraint
       logical,intent(in):: ifposres        ! position restraint flag

       !---- variable for bias potential
       logical,intent(in):: ifpotbias       ! bias potential flag

       !---- variable for local velocity
       logical,intent(in):: iflocalvel      ! flag to force local atomic velocity

       integer,intent(in):: nlvel           ! number of atoms for velocity fix
       integer,intent(in):: index_nlvel(:)  ! index of vel-fix atoms
       real(8),intent(in):: v_nlvel(:,:)    ! local velocity values

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
       !  real(8),intent(inout):: for_med(:,:)   ! medium-range force
       real(8),intent(inout):: for_short(:,:) ! short-range force

       !---- potential valiables
       real(8),intent(in):: pot_ewc ! coulomb potential(ewald self)

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

       real(8),intent(inout):: for_virilj_long(:,:)
                            ! long-range virial (L-J force)
       real(8),intent(inout):: for_virilj_med(:,:)
                            ! medium-range virial (L-J force)
       real(8),intent(inout):: for_virilj_short(:,:)
                            ! short-range virial (L-J force)

       real(8),intent(inout):: for_virimor_long(:,:)
                         ! long-range virial (Morse force)
       real(8),intent(inout):: for_virimor_med(:,:)
                         ! med-range virial (Morse force)
       real(8),intent(inout):: for_virimor_short(:,:)
                         ! short-range virial (Morse force)

       real(8),intent(inout):: for_virish_long(:,:)
                         ! long-range virial (SH force)
       real(8),intent(inout):: for_virish_med(:,:)
                         ! med-range virial (SH force)
       real(8),intent(inout):: for_virish_short(:,:)
                         ! short-range virial (SH force)

       real(8),intent(inout):: for_virirfh_long(:,:)
                            ! long-range virial (RFH force)
       real(8),intent(inout):: for_virirfh_med(:,:)
                            ! med-range virial (RFH force)
       real(8),intent(inout):: for_virirfh_short(:,:)
                            ! short-range virial (RFH force)

       real(8),intent(inout):: for_viridou_long(:,:)
                            ! long-range virial (DOU force)
       real(8),intent(inout):: for_viridou_med(:,:)
                            ! med-range virial (DOU force)
       real(8),intent(inout):: for_viridou_short(:,:)
                            ! short-range virial (DOU force)

       real(8),intent(inout):: for_viricstmnb_long(:,:)
                            ! long-range virial (custom NB force)
       real(8),intent(inout):: for_viricstmnb_med(:,:)
                            ! med-range virial (custom NB force)
       real(8),intent(inout):: for_viricstmnb_short(:,:)
                            ! short-range virial (custom NB force)

       !---- for RATTLE constraint method

       real(8),intent(in):: eps_rattle      ! tolerance (relative difference)
                                       ! for bond length constraint by RATTLE

       !---- for Nose-Hoover chain
       integer,intent(in):: mchain          ! Nose-Hoover chain number (>1 in NVT)
       real(8),intent(in):: text            ! external temp. [K] (Nose-Hoover chain)

       !---- variables for Andersen (Hoover type) barostat
       real(8),intent(in):: pext             ! external pressure [non-d]
       real(8),intent(inout):: pint          ! internal pressure
       real(8),intent(inout):: pintt(:,:)    ! internal pressure tensor
       logical,intent(in):: ifpatmcont       ! atomic pressure control
       logical,intent(in):: ifpmolcont       ! molecular pressure control
       character(5),intent(in):: pcont_axis  ! axis for pressure control (iso, aniso, etc.)

       !---- for higher order Trotter expansion
       integer,intent(inout):: next         ! iteration number of extended system
       integer,intent(inout):: nyosh        ! expansion order of Yoshida-Suzuki method

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

     end subroutine moldyn_mtk

#if defined(HF)
     subroutine moldyn_hf(iostarec,recinterval, &
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
          &               dt_long_cal, &
          &               nstep_med,nstep_short, &
          &               nstep_vir, &
          &               iflimitmove,limitdist, &
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
          &               tcont_poly,tcont_water,tcont_ma, &
          &               tcont_poly_ini,tcont_water_ini,tcont_ma_ini, &
          &               ifrattle,eps_rattle, &
          &               iflocalheat, &
          &               nlheat_poly,index_nlheat_poly, &
          &               tcont_nlheat_poly, &
          &               nlheat_water,tcont_nlheat_water, &
          &               nlheat_ma,index_nlheat_ma, &
          &               tcont_nlheat_ma, &
          &               ifregionheat, ifregionhf, ifreglange, &
          &               iftcratom, ifoutthc, &
          &               ntcregion, &
          &               tcxpos1,tcxpos2, &
          &               tcypos1,tcypos2, &
          &               tczpos1,tczpos2, &
          &               r_tcont, &
          &               nhfcregion, &
          &               hfczpos1,hfczpos2, &
          &               r_hfcont, &
          &               nlangeregion, &
          &               ltxpos1, ltxpos2, &
          &               ltypos1, ltypos2, &
          &               ltzpos1, ltzpos2, &
          &               r_ltemp, r_ltdamp, &
          &               iflocalfix, &
          &               nlfix,index_nlfix, &
          &               iflocalfixz, iflocalfixzg, &
          &               ifcnp, &
          &               ifposres, &
          &               ifpotbias, &
          &               iflocalvel, ifstrmvel, &
          &               nlvel, index_nlvel, v_nlvel, &
          &               ifoutpdb,nstep_pdbout, &
          &               resname_poly_pdb, &
          &               resname_water_pdb, &
          &               resname_matom_pdb, &
          &               xref,eref,mref,qref, &
          &               vref,timeref,tempref,pref,fref,eps0ref, &
          &               degfree_poly,degfree_water, &
          &               degfree_ma,degfree_all, &
          &               dt_med_cal,dt_short_cal, &
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
          &               rcutrpvw, ifbookrpvw, rcut_bookrpvw, &
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
          &               outhe,oubar,oupre,outhc,oupdb, &
          &               ouhtf, &
          &               oumtf, &
          &               ouumb, &
          &               ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
          &               ifoutbar,ifoutpre, &
          &               tcontinterval,outinterval, &
          &               pressinterval,heatfinterval, &
          &               yratio,zratio, &
          &               ifcalpremole,ifcalpreatom, &
          &               ifnetqcorrp, &
          &               mchain, &
          &               pint,pintt, &
          &               ifpatmcont,ifpmolcont, &
          &               ifcalljlong,nsolve,solveindex, &
          &               netchrgsq, &
          &               nstep_expand,r_expand, &
          &               nstep_maxwell, &
          &               ifhfvol, &
          &               nhfregion,hfzpos1,hfzpos2, &
          &               hftyp_pmole,hftyp_water, &
          &               hftyp_atm, &
          &               ifcalmtf,ifcalhtf, &
          &               mtfoutdir, &
          &               md_cont)

       ! ARGUMENT:
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
       real(8),intent(in):: alpha           ! parameter alpha [non-d]
       integer,intent(in):: kmax            ! parameter kmax
       real(8),intent(in):: rrcut           ! ewald real space cutoff length [non-d]

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
       real(8),intent(in):: tcont_poly_ini ! poly Temp. [K] in NVT (Woodcock) for md_h
       real(8),intent(in):: tcont_water_ini ! H2O Temp. [K] in NVT (Woodcock) for md_h
       real(8),intent(in):: tcont_ma_ini     ! MA Temp. [K] in NVT (Woodcock) for md_h

       logical,intent(in):: ifrattle        ! rattle flag

       real(8),intent(in):: rcutmor         ! Morse cutoff length [m]

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

       !---- variable for local heating
       logical,intent(in):: iflocalheat     ! local heating flag

       integer,intent(in):: nlheat_poly     ! number of poly type for local heating
       integer,intent(in):: index_nlheat_poly(:)
                                ! index of poly type for local heating
       real(8),intent(in):: tcont_nlheat_poly(:)
                                ! control temp. of poly type for local heating
       integer,intent(in):: nlheat_water    ! number of water for local heating
       real(8),intent(in):: tcont_nlheat_water
                                ! control temp. of water for local heating
       integer,intent(in):: nlheat_ma       ! number of matom type for local heating
       integer,intent(in):: index_nlheat_ma(:)
                                ! index of matom type for local heating
       real(8),intent(in):: tcont_nlheat_ma(:)
                                ! control temp. of matom type for local heating

       !---- variables for region-based temp. control
       logical,intent(in):: ifregionheat    ! region heating flag
       logical,intent(in):: iftcratom    ! region temp. control based on atom or mole.
       logical,intent(in):: ifoutthc     ! flag for outputting thermal control file

       integer,intent(inout):: ntcregion       ! number of region to control temp.
       real(8),intent(in):: tcxpos1(:),tcxpos2(:)
                                ! x-position of temp. control region
       real(8),intent(in):: tcypos1(:),tcypos2(:)
                                ! y-position of temp. control region
       real(8),intent(in):: tczpos1(:),tczpos2(:)
                                ! z-position of temp. control region
       real(8),intent(in):: r_tcont(:)       ! control temp. in each region

       !---- variables for region-based heat flux control
       logical,intent(in):: ifregionhf      ! region heat flux control flag
       integer,intent(inout):: nhfcregion      ! number of region to control heat flux
       real(8),intent(in):: hfczpos1(:),hfczpos2(:)
                                ! z-position of heat flux control region
       real(8),intent(in):: r_hfcont(:)      ! magnitude of heat flux in each region
                                ! (converted to the input energy)

       !---- variables for region-based Langevin thermostat
       logical,intent(in):: ifreglange      ! region Langevin thermostat flag
       integer,intent(in):: nlangeregion    ! number of region for Langevin thermo.
       real(8),intent(in):: ltxpos1(:),ltxpos2(:)
                                       ! x-position of temp. control region
       real(8),intent(in):: ltypos1(:),ltypos2(:)
                                       ! y-position of temp. control region
       real(8),intent(in):: ltzpos1(:),ltzpos2(:)
                                       ! z-position of temp. control region
       real(8),intent(in):: r_ltemp(:)      ! control temp. in each region
       real(8),intent(in):: r_ltdamp(:)     ! damping factor in each region [non-d]

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
       integer,intent(in):: oupre           ! output unit for outpre velocity data
       integer,intent(in):: outhc           ! output unit for outthc thermal control data
       integer,intent(in):: oupdb           ! output unit for outpdb PDB data
       integer,intent(in):: ouhtf           ! output unit for outhtf heat flux data
       integer,intent(in):: oumtf           ! output unit for outmtf momentum flux data
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
       !  real(8),intent(inout):: for_med(:,:)   ! medium-range force
       real(8),intent(inout):: for_short(:,:) ! short-range force

       !---- potential valiables
       real(8),intent(in):: pot_ewc ! coulomb potential(ewald self)

       !---- interval of temperature control
       integer,intent(in):: tcontinterval   ! interval of temp. control

       !---- interval of outputting data
       integer,intent(in):: outinterval     ! output interval of trajectory
       integer,intent(in):: pressinterval   ! interval of pressure output
       integer,intent(in):: heatfinterval   ! interval of heatf output

       !---- interval of state recording
       integer,intent(in):: recinterval     ! state record interval

       !---- valiables for pressure
       logical,intent(in):: ifcalpremole    ! pressure calculation of molecule
       logical,intent(in):: ifcalpreatom    ! pressure calculation of atom
       logical,intent(in):: ifnetqcorrp     ! net charge correction for pressure

       !---- variables for calculation of heat flux
       logical,intent(in):: ifhfvol        ! local volume-based or local surface-based

       integer,intent(in):: nhfregion     ! number of region to calculate heat flux
       real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                ! z-position of region for heat flux

       integer,intent(in):: hftyp_pmole(:)
                         ! atom- or mole-based heat flux cal. for poly
       integer,intent(in):: hftyp_water
                         ! atom- or mole-based heat flux cal. for water
       integer,intent(in):: hftyp_atm(:)    ! atom- or mole-based heat flux cal.
                                !   for each atom

       !---- flag to calculate & output momentum or heat flux
       logical,intent(in):: ifcalmtf
       logical,intent(in):: ifcalhtf

       !---- variables for outputting momentum flux
       character(3),intent(in):: mtfoutdir ! direction to output data of momentum

       ! for molecular pressure
       real(8),intent(inout):: for_viric_long(:,:)
                            ! long-range virial (coulomb force)
       real(8),intent(inout):: for_viric_med(:,:)
                            ! medium-range virial (coulomb force)
       real(8),intent(inout):: for_viric_short(:,:)
                            ! short-range virial (coulomb force)

       real(8),intent(inout):: for_virilj_long(:,:)
                         ! long-range virial (L-J force)
       real(8),intent(inout):: for_virilj_med(:,:)
                         ! medium-range virial (L-J force)
       real(8),intent(inout):: for_virilj_short(:,:)
                         ! short-range virial (L-J force)

       real(8),intent(inout):: for_virimor_long(:,:)
                         ! long-range virial (Morse force)
       real(8),intent(inout):: for_virimor_med(:,:)
                         ! med-range virial (Morse force)
       real(8),intent(inout):: for_virimor_short(:,:)
                         ! short-range virial (Morse force)

       real(8),intent(inout):: for_virish_long(:,:)
                         ! long-range virial (SH force)
       real(8),intent(inout):: for_virish_med(:,:)
                         ! med-range virial (SH force)
       real(8),intent(inout):: for_virish_short(:,:)
                         ! short-range virial (SH force)

       real(8),intent(inout):: for_virirfh_long(:,:)
                         ! long-range virial (RFH force)
       real(8),intent(inout):: for_virirfh_med(:,:)
                         ! med-range virial (RFH force)
       real(8),intent(inout):: for_virirfh_short(:,:)
                         ! short-range virial (RFH force)

       real(8),intent(inout):: for_viridou_long(:,:)
                         ! long-range virial (DOU force)
       real(8),intent(inout):: for_viridou_med(:,:)
                         ! med-range virial (DOU force)
       real(8),intent(inout):: for_viridou_short(:,:)
                         ! short-range virial (DOU force)

       real(8),intent(inout):: for_viricstmnb_long(:,:)
                            ! long-range virial (custom NB force)
       real(8),intent(inout):: for_viricstmnb_med(:,:)
                            ! med-range virial (custom NB force)
       real(8),intent(inout):: for_viricstmnb_short(:,:)
                            ! short-range virial (custom NB force)

       !---- for RATTLE constraint method

       real(8),intent(in):: eps_rattle      ! tolerance (relative difference)
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

       !---- time step of forcing Maxwell distribution
       integer,intent(in):: nstep_maxwell   ! time step of maxwell distribution

       !---- MD control parameter
       integer,intent(in):: md_cont         ! MD control flag

     end subroutine moldyn_hf
#endif

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
          &               pint,pintt, &
          &               ifpatmcont,ifpmolcont, &
          &               ifcalljlong,nsolve,solveindex, &
          &               netchrgsq, &
          &               nstep_expand,r_expand, &
          &               d_rini,d_rmax,d_econv,d_rmsf, &
          &               md_cont)

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

       real(8),intent(in):: pot_ewc ! coulomb potential(ewald self)

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

       real(8),intent(inout):: for_virilj_long(:,:)
                            ! long-range virial (L-J force)
       real(8),intent(inout):: for_virilj_med(:,:)
                            ! medium-range virial (L-J force)
       real(8),intent(inout):: for_virilj_short(:,:)
                            ! short-range virial (L-J force)

       real(8),intent(inout):: for_virimor_long(:,:)
                            ! long-range virial (Morse force)
       real(8),intent(inout):: for_virimor_med(:,:)
                            ! med-range virial (Morse force)
       real(8),intent(inout):: for_virimor_short(:,:)
                            ! short-range virial (Morse force)

       real(8),intent(inout):: for_virish_long(:,:)
                            ! long-range virial (SH force)
       real(8),intent(inout):: for_virish_med(:,:)
                            ! med-range virial (SH force)
       real(8),intent(inout):: for_virish_short(:,:)
                            ! short-range virial (SH force)

       real(8),intent(inout):: for_virirfh_long(:,:)
                            ! long-range virial (RFH force)
       real(8),intent(inout):: for_virirfh_med(:,:)
                            ! med-range virial (RFH force)
       real(8),intent(inout):: for_virirfh_short(:,:)
                            ! short-range virial (RFH force)

       real(8),intent(inout):: for_viridou_long(:,:)
                            ! long-range virial (DOU force)
       real(8),intent(inout):: for_viridou_med(:,:)
                            ! med-range virial (DOU force)
       real(8),intent(inout):: for_viridou_short(:,:)
                            ! short-range virial (DOU force)

       real(8),intent(inout):: for_viricstmnb_long(:,:)
                            ! long-range virial (custom NB force)
       real(8),intent(inout):: for_viricstmnb_med(:,:)
                            ! med-range virial (custom NB force)
       real(8),intent(inout):: for_viricstmnb_short(:,:)
                            ! short-range virial (custom NB force)

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

     end subroutine enemin_sd

  end interface

end module interface_mddriver
