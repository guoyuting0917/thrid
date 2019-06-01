!**********************************
!*  interface_mdtech.f90 Ver.2.2  *
!*      for peachgk_md.f          *
!*            by G.Kikugawa       *
!**********************************
! Time-stamp: <>

!***** This module is interface module for routines of MD techniques *****

module interface_mdtech

  interface

     subroutine transcor( npoly, nwater, nmatom, &
          &               xcel, ycel, zcel)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: nwater          ! number of H2O molecules
       integer,intent(in):: nmatom          ! number of monatomic molecules

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

     end subroutine transcor

     subroutine mkmaxvel(npoly, nwater, nmatom, &
          &              tcont_poly, tcont_water, tcont_ma, tempref)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: nwater          ! number of H2O molecules
       integer,intent(in):: nmatom          ! number of monatomic molecules

       real(8),intent(in):: tcont_poly(:)    ! poly Temp. [non-d] in NVT
       real(8),intent(in):: tcont_water      ! H2O Temp. [non-d] in NVT
       real(8),intent(in):: tcont_ma(:)      ! monatomic mole. Temp. [non-d] in NVT

       real(8),intent(in):: tempref          ! temperature base value [K]

     end subroutine mkmaxvel

     subroutine cell_expand( xref,current_step, &
          &                  xcel,ycel,zcel, &
          &                  yratio,zratio, &
          &                  rcut_book,ifbook,rcut, &
          &                  npoly,nwater,nmatom, &
          &                  r_expand )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: xref             ! distanse base value [m]
       integer,intent(in):: current_step    ! current time step

       real(8),intent(in):: yratio           ! y cell ratio of y to x
       real(8),intent(in):: zratio           ! z cell ratio of z to x

       real(8),intent(in):: rcut_book        ! cut off radius of bookkeeping
       logical,intent(in):: ifbook          ! flag for bookkeeping

       real(8),intent(in):: rcut             ! vdw cutoff length [non-d]

       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: nwater          ! number of H2O molecules
       integer,intent(in):: nmatom          ! number of monatomic molecules

       real(8),intent(in):: r_expand         ! expansion ratio of cell (volume)

       !     INPUT&OUTPUT
       real(8),intent(inout):: xcel             ! x cell length[non-d]
       real(8),intent(inout):: ycel             ! y cell length[non-d]
       real(8),intent(inout):: zcel             ! z cell length[non-d]

     end subroutine cell_expand

     subroutine upewk( alpha, &
          &            xcel, ycel, zcel )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: alpha            ! parameter alpha [non-d]

       real(8),intent(in):: xcel             ! x cell length [non-d]
       real(8),intent(in):: ycel             ! y cell length [non-d]
       real(8),intent(in):: zcel             ! z cell length [non-d]

     end subroutine upewk

     subroutine erf_corr_cutoff( xcel, ycel, zcel, &
          &                      nfft1, nfft2, nfft3 )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: xcel             ! initial x0 cell length[m]
       real(8),intent(in):: ycel             ! initial y0 cell length[m]
       real(8),intent(in):: zcel             ! initial z0 cell length[m]
       integer,intent(in):: nfft1,nfft2,nfft3

     end subroutine erf_corr_cutoff

     subroutine cenfix(ifcenterfix_all, &
          &            ifcenterfix_poly, &
          &            ifcenterfix_water, &
          &            ifcenterfix_ma, &
          &            cenfix_free, &
          &            ifcenterfix_polytyp, &
          &            ifcenterfix_watertyp, &
          &            ifcenterfix_matyp, &
          &            npoly,npolytyp,npoly_mole,npoly_atom, &
          &            nwater, &
          &            nmatom,nmatyp,nmatomtyp, &
          &            iflocalfix, &
          &            nlfix,index_nlfix,   &
          &            iflocalfixz,iflocalfixzg)

       ! ARGUMENT:
       !     INPUT
       logical,intent(in):: ifcenterfix_all   ! center fix for all
       logical,intent(in):: ifcenterfix_poly  ! center fix for polymer
       logical,intent(in):: ifcenterfix_water ! center fix for water
       logical,intent(in):: ifcenterfix_ma    ! center fix for monatomic mole.
       character(4),intent(in):: cenfix_free  ! COM not fixed in this direction

       logical,intent(in):: ifcenterfix_polytyp(:) ! center fix for each polymer
       logical,intent(in):: ifcenterfix_watertyp   ! center fix for each water
       logical,intent(in):: ifcenterfix_matyp(:)
                                          ! center fix for each monatomic mole.

       integer,intent(in):: npoly           ! all number of poly
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       logical,intent(in):: iflocalfix      ! fix atoms flag

       integer,intent(in):: nlfix           ! number of fix atoms
       integer,intent(in):: index_nlfix(:)  ! index of fix atoms

       logical,intent(in):: iflocalfixz     ! flag for fixing z coordinate of atoms

       logical,intent(in):: iflocalfixzg    ! flag for fixing z coordinate of COM of molecules

     end subroutine cenfix

     subroutine mklist_cell( rcut_book, &
          &                  xcel, ycel, zcel)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: rcut_book        ! cut off radius of bookkeeping

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

     end subroutine mklist_cell

     subroutine mklist2a( rcut_book, &
          &               xcel,ycel,zcel, &
          &               ifbook )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: rcut_book        ! cut off radius of bookkeeping

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       logical,intent(in):: ifbook          ! flag for bookkeeping

     end subroutine mklist2a

     subroutine mklist_mor_cell( rcut_bookmor, &
          &                      xcel, ycel, zcel )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: rcut_bookmor ! cut off radius of bookkeeping of Morse

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

     end subroutine mklist_mor_cell

     subroutine mklist2a_mor( rcut_bookmor, &
          &                   xcel,ycel,zcel, &
          &                   ifbookmor)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: rcut_bookmor ! cut off radius of bookkeeping of Morse

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       logical,intent(in):: ifbookmor ! flag for bookkeeping of Morse interaction


     end subroutine mklist2a_mor

     subroutine mklist_sh_cell( rcut_booksh, &
          &                     xcel, ycel, zcel )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: rcut_booksh     ! cut off radius of bookkeeping of SH

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

     end subroutine mklist_sh_cell

     subroutine mklist2a_sh( rcut_booksh, &
          &                  xcel,ycel,zcel, &
          &                  ifbooksh)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: rcut_booksh     ! cut off radius of bookkeeping of SH

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       logical,intent(in):: ifbooksh     ! flag for bookkeeping of SH interaction

     end subroutine mklist2a_sh

     subroutine mklist_rfh_cell( rcut_bookrfh, &
          &                      xcel, ycel, zcel, &
          &                      inttype_rfh )

       ! ARGUMENT:
       !     INPUT
       real(8):: rcut_bookrfh     ! cut off radius of bookkeeping[non-d] of RFH

       real(8):: xcel             ! x cell length[non-d]
       real(8):: ycel             ! y cell length[non-d]
       real(8):: zcel             ! z cell length[non-d]

       integer:: inttype_rfh     ! flag for RFH interact type
                                 ! = 4, INTTYPE_RFHFO
                                 ! = 5, INTTYPE_RFHOO
                                 ! = 6, INTTYPE_RFHOH

     end subroutine mklist_rfh_cell

     subroutine mklist2a_rfh( rcut_bookrfh, &
          &                   xcel, ycel, zcel, &
          &                   ifbookrfh, &
          &                   inttype_rfh)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: rcut_bookrfh
                         ! cut off radius of bookkeeping[non-d] of RFH

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       logical,intent(in):: ifbookrfh
                         ! flag for bookkeeping of RFH interaction

       integer,intent(in):: inttype_rfh     ! flag for RFH interact type
                                ! = 4, INTTYPE_RFHFO
                                ! = 5, INTTYPE_RFHOO
                                ! = 6, INTTYPE_RFHOH

     end subroutine mklist2a_rfh

     subroutine mklist_dou_cell( rcut_bookdou, &
          &                      xcel, ycel, zcel, &
          &                      inttype_dou )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: rcut_bookdou
                         ! cut off radius of bookkeeping[non-d] of DOU

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       integer,intent(in):: inttype_dou     ! flag for DOU interact type
                                ! = 7, INTTYPE_DOUO
                                ! = 8, INTTYPE_DOUH

     end subroutine mklist_dou_cell

     subroutine mklist2a_dou( rcut_bookdou, &
          &                   xcel, ycel, zcel, &
          &                   ifbookdou, &
          &                   inttype_dou )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: rcut_bookdou
                         ! cut off radius of bookkeeping[non-d] of DOU

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       logical,intent(in):: ifbookdou   ! flag for bookkeeping of DOU interaction

       integer,intent(in):: inttype_dou     ! flag for DOU interact type
                                ! = 7, INTTYPE_DOUO
                                ! = 8, INTTYPE_DOUH

     end subroutine mklist2a_dou

     subroutine mklist2a_rpvw( rcut_bookrpvw, zcel, ifbookrpvw)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: rcut_bookrpvw     ! cut off radius of bookkeeping of RP-VW

       !  real(8),intent(in):: xcel             ! x cell length[non-d]
       !  real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       logical,intent(in):: ifbookrpvw      ! flag for bookkeeping of RP-VW interacti

     end subroutine mklist2a_rpvw

     subroutine mklist2a_cstmnb(ifcellindex_cstmnb,ifbookcstmnb, &
          &                     xcel,ycel,zcel, &
          &                     current_step)

       ! ARGUMENT:
       !     INPUT
       logical,intent(in):: ifcellindex_cstmnb ! flag for cell index (custom NB)
       logical,intent(in):: ifbookcstmnb
                                ! flag for bookkeeping of custom NB interaction

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       integer,intent(in):: current_step     ! current MD step

     end subroutine mklist2a_cstmnb

     subroutine calforce(mts_flag,force, &
          &              for_viri_coul,pot_viri_coul, &
          &              for_viri_lj,pot_viri_lj, &
          &              for_viri_mor,pot_viri_mor, &
          &              for_viri_sh,pot_viri_sh, &
          &              for_viri_rfh,pot_viri_rfh, &
          &              for_viri_dou,pot_viri_dou, &
          &              for_viri_cstmnb,pot_viri_cstmnb, &
          &              for_viri_cstmnbex,pot_viri_cstmnbex, &
          &              atm_viri_bond,atm_viri_angl, &
          &              atm_viri_tors,atm_viri_14, &
          &              pot_virit_coul,pot_virit_lj, &
          &              pot_virit_mor, &
          &              pot_virit_sh, &
          &              pot_virit_rfh, &
          &              pot_virit_dou, &
          &              pot_virit_cstmnb, &
          &              pot_virit_cstmnbex, &
          &              atm_virit_bond,atm_virit_angl, &
          &              atm_virit_tors,atm_virit_14, &
          &              mts_bond, mts_angl, mts_anglub, &
          &              mts_tors,mts_torsrb,mts_torsim, &
          &              mts_vdw,mts_ewr,mts_ewk, &
          &              mts_vdw14,mts_elc14, &
          &              mts_mor,mts_sh,mts_rfh,mts_dou,mts_cnpvw, &
          &              mts_cstmnb, &
          &              mts_posres, &
          &              mts_potbias, &
          &              npoly,npolytyp,npoly_mole,npoly_atom, &
          &              nwater,nmatom,nmatyp,nmatomtyp, &
          &              xcel,ycel,zcel, &
          &              ifewald,alpha,kmax,rrcut, &
          &              ifspme,nfft1,nfft2,nfft3,pme_order, &
          &              eps0,div_factor_14vdw,div_factor_14elc, &
          &              rcut, &
          &              rcutmor, &
          &              rcutsh, &
          &              rcutrfhfo,rcutrfhoo,rcutrfhoh, &
          &              rcutdouo,rcutindouo,rcutdouh, &
          &              rcutrpvw, &
          &              pot_ewc,pot_ewnq, &
          &              pot_tot,pot_nonbon,pot_vdw,pot_elc,pot_ewk, &
          &              pot_vdw14,pot_elc14, &
          &              pot_bond,pot_angl,pot_anglub, &
          &              pot_tors,pot_torsrb,pot_torsim, &
          &              pot_mor, &
          &              pot_sh, &
          &              pot_rfh, &
          &              pot_dou, &
          &              pot_rpvw,pot_vw, &
          &              pot_cstmnb, pot_cstmnbex, &
          &              pot_posres, &
          &              pot_pbias, &
          &              ifcalpremole,ifcalpreatom, &
          &              ifnetqcorrp, &
          &              netchrgsq, &
          &              ifcstmnb, &
          &              ifposres, &
          &              ifpotbias, &
          &              outinterval,pressinterval,recinterval, &
          &              nstep_book, &
          &              md_cont, &
          &              ncstmnbex, &
          &              current_step, &
          &              nstep_short,istep_short)

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
       integer,intent(in):: mts_cnpvw       ! MTS flag for RP-VW interaction
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

       real(8),intent(in):: rcutsh           ! SH cutoff length [non-d]

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

       !---- variables for timestep control
       integer,intent(in):: outinterval     ! output interval of trajectory
       integer,intent(in):: pressinterval   ! interval of pressure output
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

     end subroutine calforce

#if defined(HF)
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

     end subroutine calforce_hf
#endif

     subroutine accforce( natom, atmmass, force )

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: natom           ! number of atoms
       real(8),intent(in):: atmmass(:)       ! atomic mass

       !     OUTPUT
       real(8),intent(inout):: force(:,:)
                                ! input force
                                ! output accel

     end subroutine accforce

     subroutine rattle_c( atmcor_old, &
          &               eps_rattle, dt_short_cal, dt_long_cal, &
          &               istep_short, nstep_short )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: atmcor_old(:,:)  ! atmcor at T

       real(8),intent(in):: eps_rattle       ! tolerance (relative difference)
                                ! for bond length constraint by RATTLE

       real(8),intent(in):: dt_short_cal     ! time step of short force
       real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]

       integer:: istep_short
       integer:: nstep_short     ! number of step for short force

     end subroutine rattle_c

     subroutine rattle_cmtk(atmcor_old, &
         &                  eps_rattle, dt_short_cal, dt_long_cal, &
         &                  istep_short, nstep_short, &
         &                  rv, rvv, &
         &                  pcont_axis)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: atmcor_old(:,:)  ! atmcor at T

       real(8),intent(in):: eps_rattle       ! tolerance (relative difference)
                                ! for bond length constraint by RATTLE

       real(8),intent(in):: dt_short_cal     ! time step of short force
       real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]

       integer,intent(in):: istep_short
       integer,intent(in):: nstep_short     ! number of step for short force

       real(8),intent(in):: rv               ! = poly*exp(dt_short/2*veps)
       real(8),intent(in):: rvv(:)           ! = polyrv(:)*exp(dt_short/2*vboxg(:))

       character(5),intent(in):: pcont_axis  ! axis for pressure control (iso, aniso, etc.)

     end subroutine rattle_cmtk

     subroutine rattlep_v( eps_rattle, &
          &                dt_long_cal, dt_short_cal, &
          &                istep_short, nstep_short, &
          &                expcoeff, &
          &                atm_viri_const, atm_virit_const )

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: eps_rattle       ! tolerance (relative difference)
                                ! for bond length constraint by RATTLE

       real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]
       real(8),intent(in):: dt_short_cal     ! time step of short force

       integer,intent(in):: istep_short     ! do loop index for short range forces
       integer,intent(in):: nstep_short     ! number of step for short force

       real(8),intent(in):: expcoeff         ! = coefficient of exponential function

       !     INPUT&OUTPUT
       real(8),intent(inout):: atm_viri_const   ! virial(constraint force)
       real(8),intent(inout):: atm_virit_const(:,:) ! virial tensor (constraint)

     end subroutine rattlep_v

     subroutine rattle_c_z(dt_short_cal,dt_long_cal,   &
          &                istep_short,nstep_short)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: dt_short_cal    ! time step of short force
       real(8),intent(in):: dt_long_cal     ! time step of long force [non-d]

       integer,intent(in):: istep_short
       integer,intent(in):: nstep_short     ! number of step for short force

     end subroutine rattle_c_z

     subroutine rattlep_v_z(dt_long_cal, dt_short_cal,   &
          &                 istep_short, nstep_short)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]
       real(8),intent(in):: dt_short_cal     ! time step of short force

       integer,intent(in):: istep_short     ! do loop index for short range forces
       integer,intent(in):: nstep_short     ! number of step for short force

     end subroutine rattlep_v_z

     subroutine rattle_c_zg(dt_short_cal,dt_long_cal, &
          &                 istep_short,nstep_short)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: dt_short_cal    ! time step of short force
       real(8),intent(in):: dt_long_cal     ! time step of long force [non-d]

       integer,intent(in):: istep_short
       integer,intent(in):: nstep_short     ! number of step for short force

     end subroutine rattle_c_zg

     subroutine rattlep_v_zg(dt_long_cal,dt_short_cal, &
          &                  istep_short,nstep_short)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]
       real(8),intent(in):: dt_short_cal     ! time step of short force

       integer,intent(in):: istep_short     ! do loop index for short range forces
       integer,intent(in):: nstep_short     ! number of step for short force

     end subroutine rattlep_v_zg

     subroutine shake_c(atmcor_old, &
          &             eps_rattle)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: atmcor_old(:,:)  ! atmcor at T

       real(8),intent(in):: eps_rattle       ! tolerance (relative difference)
                                ! for bond length constraint by RATTLE

     end subroutine shake_c

     subroutine shake_c_z()

       ! ARGUMENT:

     end subroutine shake_c_z

     subroutine shake_c_zg()

       ! ARGUMENT:

     end subroutine shake_c_zg

     subroutine limitvel(dt_short_cal, &
          &              limitdist)

       ! ARGUMENTS:
       !     INPUT
       real(8),intent(in):: dt_short_cal     ! time step of short force
       real(8),intent(in):: limitdist  ! maximum atomic displacement when doing
                              ! time integration [non-d] (structure relaxation)

     end subroutine limitvel

     subroutine ass_strmvel(npoly,npolytyp, &
          &                 npoly_mole,npoly_atom, &
          &                 nwater, &
          &                 nmatom,nmatyp,nmatomtyp, &
          &                 xcel,ycel,zcel, &
          &                 iftcratom, &
          &                 ifstrmvel)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: npoly           ! number of polymer1
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

       logical,intent(in):: iftcratom    ! region temp. control based on atom or mole.

       logical,intent(in):: ifstrmvel     ! flag to input and use streaming velocity
     end subroutine ass_strmvel

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

       real(8),intent(out):: pot_virict_all(:,:)   ! virial tensor (coulomb potential)
       real(8),intent(out):: pot_viriljt_all(:,:)  ! virial tensor (L-J potential)
       real(8),intent(out):: pot_virimort_all(:,:) ! virial tensor (Morse potential)
       real(8),intent(out):: pot_virisht_all(:,:)  ! virial tensor (SH potential)
       real(8),intent(out):: pot_virirfht_all(:,:) ! virial tensor (RFH potential)
       real(8),intent(out):: pot_viridout_all(:,:) ! virial tensor (DOU potential)
       real(8),intent(out):: pot_viricstmnbt_all(:,:)
                                          ! virial tensor (custom NB potential)
       real(8),intent(out):: pot_viricstmnbext_all(:,:,:)
                                     ! virial tensor (extra custom NB potential)

       real(8),intent(out):: pot_virit_all(:,:)  ! all virial tensor of potential term
       real(8),intent(out):: virit_fdotd(:,:)   ! virial tensor of correction term F.d

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

     end subroutine calpress

     subroutine contemp_e( npoly,npolytyp,npoly_mole, &
          &                nwater, &
          &                nmatom,nmatyp,nmatomtyp, &
          &                degfree_poly,degfree_water, &
          &                degfree_ma,degfree_all, &
          &                tcont_polyt,tcont_watert,tcont_mat, &
          &                tcont_polyinit,tcont_waterinit,tcont_mainit, &
          &                tfactor_poly,tfactor_water,tfactor_ma)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       integer,intent(in):: degfree_poly(:) ! degree of freedom of polymer1
       integer,intent(in):: degfree_water   ! degree of freedom of H2O
       integer,intent(in):: degfree_ma(:)   ! degree of freedom of monatomic mole.
       integer,intent(in):: degfree_all     ! degree of freedom of all atoms

       real(8),intent(in):: tcont_polyt(:)   ! poly Temp. [non-d] in NVT
       real(8),intent(in):: tcont_watert     ! H2O Temp. [non-d] in NVT
       real(8),intent(in):: tcont_mat(:)     ! monatomic mole. Temp. [non-d] in NVT

       real(8),intent(in):: tcont_polyinit(:) ! poly Temp. [non-d] in NVT
       real(8),intent(in):: tcont_waterinit  ! H2O Temp. [non-d] in NVT
       real(8),intent(in):: tcont_mainit(:)  ! MA Temp. [non-d] in NVT

       real(8),intent(in):: tfactor_poly(:)  ! temperature control factor of polymer1
       real(8),intent(in):: tfactor_water    ! temperature control factor of H2O
       real(8),intent(in):: tfactor_ma(:)    ! temperature control factor of MA

     end subroutine contemp_e

     subroutine contemp( npoly,npolytyp,npoly_mole, &
          &              nwater, &
          &              nmatom,nmatyp,nmatomtyp, &
          &              degfree_poly,degfree_water, &
          &              degfree_ma,degfree_all, &
          &              tcont_polyt,tcont_watert,tcont_mat, &
          &              tcont_polyinit,tcont_waterinit,tcont_mainit, &
          &              tfactor_poly,tfactor_water,tfactor_ma )

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       integer,intent(in):: degfree_poly(:) ! degree of freedom of polymer1
       integer,intent(in):: degfree_water   ! degree of freedom of H2O
       integer,intent(in):: degfree_ma(:)   ! degree of freedom of monatomic mole.
       integer,intent(in):: degfree_all     ! degree of freedom of all atoms

       real(8),intent(in):: tcont_polyt(:)   ! poly Temp. [non-d] in NVT
       real(8),intent(in):: tcont_watert     ! H2O Temp. [non-d] in NVT
       real(8),intent(in):: tcont_mat(:)     ! monatomic mole. Temp. [non-d] in NVT

       real(8),intent(in):: tcont_polyinit(:) ! poly Temp. [non-d] in NVT
       real(8),intent(in):: tcont_waterinit  ! H2O Temp. [non-d] in NVT
       real(8),intent(in):: tcont_mainit(:)  ! MA Temp. [non-d] in NVT

       real(8),intent(in):: tfactor_poly(:)  ! temperature control factor of polymer1
       real(8),intent(in):: tfactor_water    ! temperature control factor of H2O
       real(8),intent(in):: tfactor_ma(:)    ! temperature control factor of MA

     end subroutine contemp

     subroutine contemp_nve( npoly,npolytyp,npoly_mole,npoly_atom, &
          &                  nwater, &
          &                  nmatom,nmatyp,nmatomtyp, &
          &                  degfree_poly,degfree_water, &
          &                  degfree_ma,degfree_all, &
          &                  tcont_polyt,tcont_watert,tcont_mat, &
          &                  tcont_polyinit,tcont_waterinit, &
          &                  tcont_mainit, &
          &                  tfactor_poly,tfactor_water,tfactor_ma, &
          &                  nlheat_poly,index_nlheat_poly, &
          &                  nlheat_water, &
          &                  nlheat_ma,index_nlheat_ma)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: npoly           ! number of polymer1
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

       real(8),intent(in):: tcont_polyt(:)   ! poly Temp. [non-d] in NVT
       real(8),intent(in):: tcont_watert     ! H2O Temp. [non-d] in NVT
       real(8),intent(in):: tcont_mat(:)     ! monatomic mole. Temp. [non-d] in NVT

       real(8),intent(in):: tcont_polyinit(:) ! poly Temp. [non-d] in NVT
       real(8),intent(in):: tcont_waterinit  ! H2O Temp. [non-d] in NVT
       real(8),intent(in):: tcont_mainit(:)  ! MA Temp. [non-d] in NVT

       real(8),intent(in):: tfactor_poly(:)  ! temperature control factor of polymer1
       real(8),intent(in):: tfactor_water    ! temperature control factor of H2O
       real(8),intent(in):: tfactor_ma(:)    ! temperature control factor of MA

       integer,intent(in):: nlheat_poly     ! number of poly type for local heating
       integer,intent(in):: index_nlheat_poly(:)
                                ! index of poly type for local heating
       integer,intent(in):: nlheat_water    ! number of water for local heating
       integer,intent(in):: nlheat_ma       ! number of matom type for local heating
       integer,intent(in):: index_nlheat_ma(:)
                                ! index of matom type for local heating

     end subroutine contemp_nve

     subroutine contemp_region(npoly,npolytyp, &
          &                    npoly_mole,npoly_atom, &
          &                    nwater, &
          &                    nmatom,nmatyp,nmatomtyp, &
          &                    degfree_poly,degfree_water, &
          &                    degfree_ma, &
          &                    xcel,ycel,zcel, &
          &                    iftcratom, &
          &                    ntcregion, &
          &                    tcxpos1,tcxpos2, &
          &                    tcypos1,tcypos2, &
          &                    tczpos1,tczpos2, &
          &                    r_tcont, &
          &                    ifoutthc, &
          &                    det_ene_kin)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       integer,intent(in):: degfree_poly(:) ! degree of freedom of poly
       integer,intent(in):: degfree_water   ! degree of freedom of H2O
       integer,intent(in):: degfree_ma(:)   ! degree of freedom of monatomic mole.

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       logical,intent(in):: iftcratom    ! region temp. control based on atom or mole.
       integer,intent(in):: ntcregion       ! number of region to control temp.
       real(8),intent(in):: tcxpos1(:),tcxpos2(:)
                                ! x-position of temp. control region
       real(8),intent(in):: tcypos1(:),tcypos2(:)
                                ! y-position of temp. control region
       real(8),intent(in):: tczpos1(:),tczpos2(:)
                                ! z-position of temp. control region
       real(8),intent(in):: r_tcont(:)       ! control temp. in each region

       logical,intent(in):: ifoutthc     ! flag for outputting thermal control file

       !     OUTPUT
       real(8),intent(out):: det_ene_kin(:)          ! for outputting thermal control data

     end subroutine contemp_region

     subroutine conhf_region( npoly,npolytyp, &
          &                   npoly_mole,npoly_atom, &
          &                   nwater, &
          &                   nmatom,nmatyp,nmatomtyp, &
          &                   xcel,ycel,zcel, &
          &                   iftcratom, &
          &                   nhfcregion, &
          &                   hfczpos1,hfczpos2, &
          &                   r_hfcont )

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: npoly           ! number of polymer1
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

       logical,intent(in):: iftcratom    ! region temp. control based on atom or mole.

       integer,intent(in):: nhfcregion      ! number of region to control heat flux
       real(8),intent(in):: hfczpos1(:),hfczpos2(:)
                                ! z-position of heat flux control region
       real(8),intent(in):: r_hfcont(:)      ! magnitude of heat flux in each region
                                ! (converted to the input energy)

     end subroutine conhf_region

     subroutine langevin_regist(npoly,npolytyp, &
          &                     npoly_mole,npoly_atom, &
          &                     nwater, &
          &                     nmatom,nmatyp,nmatomtyp, &
          &                     degfree_poly,degfree_water, &
          &                     degfree_ma, &
          &                     xcel,ycel,zcel, &
          &                     iftcratom, &
          &                     nlangeregion, &
          &                     ltxpos1,ltxpos2, &
          &                     ltypos1,ltypos2, &
          &                     ltzpos1,ltzpos2)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       integer,intent(in):: degfree_poly(:) ! degree of freedom of poly
       integer,intent(in):: degfree_water   ! degree of freedom of H2O
       integer,intent(in):: degfree_ma(:)   ! degree of freedom of monatomic mole.

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       logical,intent(in):: iftcratom    ! region temp. control based on atom or mole.

       integer,intent(in):: nlangeregion    ! number of region for Langevin thermo.
       real(8),intent(in):: ltxpos1(:),ltxpos2(:)
                                       ! x-position of temp. control region
       real(8),intent(in):: ltypos1(:),ltypos2(:)
                                       ! y-position of temp. control region
       real(8),intent(in):: ltzpos1(:),ltzpos2(:)
                                       ! z-position of temp. control region

     end subroutine langevin_regist

     subroutine langevin_region(dt_long_cal, &
          &                     iftcratom, &
          &                     nlangeregion, &
          &                     r_ltemp, r_ltdamp, &
          &                     ifoutthc, &
          &                     det_ene_kin, &
          &                     ifhalf)

       ! ARGUMENTS:
       !     INPUT
       real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]

       logical,intent(in):: iftcratom    ! region temp. control based on atom or mole.

       integer,intent(in):: nlangeregion    ! number of region for Langevin thermo.

       real(8),intent(in):: r_ltemp(:)      ! control temp. in each region
       real(8),intent(in):: r_ltdamp(:)     ! damping factor in each region [non-d]

       logical,intent(in):: ifoutthc     ! flag for outputting thermal control file
       logical,intent(in):: ifhalf   ! flag for former or latter half of integration
                                ! former=.true. or latter=.false.

       !     OUTPUT
       real(8),intent(out):: det_ene_kin(:) ! for outputting thermal control data

     end subroutine langevin_region

     subroutine wrsta(iostarec, &
         &           npoly,nwater,nmatom, &
         &           xcel,ycel,zcel, &
         &           xref,vref,timeref,pref, &
         &           maxnstep, &
         &           mchain, &
         &           pint, pintt)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: iostarec        ! state record file unit

       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: nwater          ! number of H2O molecules
       integer,intent(in):: nmatom          ! number of monatomic molecules

       real(8),intent(in):: xcel             ! x cell length
       real(8),intent(in):: ycel             ! y cell length
       real(8),intent(in):: zcel             ! z cell length

       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: vref             ! velocity base value [m/s]
       real(8),intent(in):: timeref          ! time base value [sec]
       real(8),intent(in):: pref             ! pressure base value [Pa]

       integer,intent(in):: maxnstep        ! current or maximum step of MD

       integer,intent(in):: mchain          ! Nose-Hoover chain number

       real(8),intent(in):: pint             ! internal pressure
       real(8),intent(in):: pintt(:,:)       ! internal pressure tensor

     end subroutine wrsta

     subroutine calkin(npoly,npolytyp,npoly_mole, &
          &            nwater, &
          &            nmatom,nmatyp,nmatomtyp, &
          &            degfree_poly,degfree_water,degfree_ma, &
          &            ene_kin_poly,ene_kin_water,ene_kin_ma, &
          &            temp_poly,temp_water,temp_ma, &
          &            degfree_all,ene_kin_all,temp_all, &
          &            mchain,text_c, &
          &            ene_kin_th,ene_pot_th, &
          &            pext_c, &
          &            ene_kin_ba,ene_pot_ba,extra_pot_ba, &
          &            pcont_axis)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       integer,intent(in):: degfree_poly(:) ! degree of freedom of polymer1
       integer,intent(in):: degfree_water   ! degree of freedom of H2O
       integer,intent(in):: degfree_ma(:)   ! degree of freedom of monatomic mole.
       integer,intent(in):: degfree_all     ! degree of freedom of all atoms

       integer,intent(in):: mchain          ! Nose-Hoover chain number (>1 in NVT)
       real(8),intent(in):: text_c           ! external temp. [K] (Nose-Hoover chain)

       real(8),intent(in):: pext_c           ! external pressure [non-d]

       character(5),intent(in),optional:: pcont_axis  ! axis for pressure control (iso, aniso, etc.)

       !     OUTPUT

       real(8),intent(out):: ene_kin_poly(:)  ! kinetic energy of polymer1
       real(8),intent(out):: ene_kin_water    ! kinetic energy of water
       real(8),intent(out):: ene_kin_ma(:)    ! kinetic energy of monatomic mole.
       real(8),intent(out):: ene_kin_all      ! kinetic energy
       real(8),intent(out):: temp_poly(:)     ! temperature of polymer1
       real(8),intent(out):: temp_water       ! temperature of H2O
       real(8),intent(out):: temp_ma(:)       ! temperature of MA
       real(8),intent(out):: temp_all         ! temperature of all atoms

       real(8),intent(out):: ene_kin_th       ! kinetic energy of NHC thermostat
       real(8),intent(out):: ene_pot_th       ! potential energy of NHC thermostat

       real(8),intent(out):: ene_kin_ba       ! kinetic energy of Andersen barostat
       real(8),intent(out):: ene_pot_ba       ! potential energy of Andersen barostat
       real(8),intent(out):: extra_pot_ba     ! extra potentail (=kTxi(1))

     end subroutine calkin

     subroutine outene(ouene,current_step,eref,tempref, &
          &            npolytyp,nmatyp, &
          &            pot_tot,pot_nonbon,pot_vdw, &
          &            pot_elc,pot_ewk,pot_ewc,pot_ewnq, &
          &            pot_vdw14,pot_elc14, &
          &            pot_bond,pot_angl,pot_anglub, &
          &            pot_tors,pot_torsrb,pot_torsim, &
          &            pot_mor, &
          &            pot_sh, &
          &            pot_rfh, &
          &            pot_dou, &
          &            pot_rpvw,pot_vw, &
          &            pot_cstmnb, &
          &            ncstmnbex,pot_cstmnbex, &
          &            pot_posres, &
          &            pot_pbias, &
          &            ene_kin_poly,ene_kin_water, &
          &            ene_kin_ma,ene_kin_all, &
          &            ene_tot, &
          &            temp_poly,temp_water,temp_ma,temp_all, &
          &            ene_kin_th,ene_pot_th, &
          &            ene_kin_ba,ene_pot_ba, &
          &            ene_conserve)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: ouene           ! output unit for output energy data

       integer,intent(in):: current_step    ! current time step

       !---- base value for non-dimensionalize
       real(8),intent(in):: eref             ! energy base value [J]

       real(8),intent(in):: tempref          ! temperature base value [K]

       !---- number of type of poly and matom
       integer,intent(in):: npolytyp        ! number of poly type

       integer,intent(in):: nmatyp          ! number of species of monatomic mole.

       !---- potential valiables
       real(8),intent(inout):: pot_tot          ! all potential
       real(8),intent(inout):: pot_nonbon       ! all non-bonded potential
       real(8),intent(inout):: pot_vdw          ! vdw potential
       real(8),intent(inout):: pot_elc          ! coulomb potential(ewald real)
       real(8),intent(inout):: pot_ewk          ! coulomb potential(ewald wave)
       real(8),intent(in):: pot_ewc             ! coulomb potential(ewald self)
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

       integer,intent(in):: ncstmnbex           ! number of extra custom NB output

       !---- kinematic energy
       real(8),intent(inout):: ene_kin_poly(:)  ! kinetic energy
       real(8),intent(inout):: ene_kin_water    ! kinetic energy
       real(8),intent(inout):: ene_kin_ma(:)    ! kinetic energy
       real(8),intent(inout):: ene_kin_all      ! kinetic energy

       real(8),intent(inout):: ene_kin_th       ! kinetic energy of NHC thermostat
       real(8),intent(inout):: ene_pot_th       ! potential energy of NHC thermostat

       real(8),intent(inout):: ene_kin_ba       ! kinetic energy of Andersen barostat
       real(8),intent(inout):: ene_pot_ba       ! potential energy of Andersen barostat

       !---- all energy

       real(8),intent(inout):: ene_tot          ! total Hamiltonian
       real(8),intent(inout):: ene_conserve     ! total conserved quantity

       !---- temperature

       real(8),intent(inout):: temp_poly(:)     ! temperature of polymer1
       real(8),intent(inout):: temp_water       ! temperature of H2O
       real(8),intent(inout):: temp_ma(:)       ! temperature of monatomic mole.
       real(8),intent(inout):: temp_all         ! temperature of all atoms

     end subroutine outene

     subroutine outpos(oupos,current_step,xref)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: oupos          ! output unit for output position data

       integer,intent(in):: current_step    ! current time step

       !---- base value for non-dimensionalize
       real(8),intent(in):: xref             ! distanse base value [m]

     end subroutine outpos

     subroutine outvel(ouvel,current_step,vref)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: ouvel          ! output unit for output velocity data

       integer,intent(in):: current_step    ! current time step

       !---- base value for non-dimensionalize
       real(8),intent(in):: vref             ! velocity base value [m/s]

     end subroutine outvel

     subroutine outfor(oufor,current_step,fref, &
         &             for_long,for_short)

         !ARGUMENTS:
         !     INPUT
         integer,intent(in):: oufor           ! output unit for output force data

         integer,intent(in):: current_step    ! current time step

         !---- base value for non-dimensionalize
         real(8),intent(in):: fref             ! force base value [N]

         !---- force of each atoms
         real(8),intent(in):: for_long(:,:)     ! long-range force
         !  real(8),intent(in):: for_med(:,:)   ! medium-range force
         real(8),intent(in):: for_short(:,:)    ! short-range force

     end subroutine outfor

     subroutine outfor_hf(oufor,current_step,fref, &
         &                for_long,for_short, &
         &                ifhfvol, &
         &                nhfregion,hfzpos1,hfzpos2, &
         &                hftyp_atm, &
         &                molecom, &
         &                viribot_long_atm,viribot_short_atm, &
         &                viriant_long_atm,viriant_short_atm, &
         &                viritot_long_atm,viritot_short_atm, &
         &                viri14t_long_atm,viri14t_short_atm, &
         &                virielint_long_atm,virielint_short_atm, &
         &                viriljint_long_atm,viriljint_short_atm, &
         &                virielct_long_atm,virielct_short_atm, &
         &                viriljt_long_atm,viriljt_short_atm, &
         &                virimort_long_atm,virimort_short_atm, &
         &                virisht_long_atm,virisht_short_atm, &
         &                virirfht_long_atm,virirfht_short_atm, &
         &                viridout_long_atm,viridout_short_atm, &
         &                viricstmnbt_long_atm,viricstmnbt_short_atm, &
         &                viricstmnbext_long_atm,viricstmnbext_short_atm)

         !ARGUMENTS:
         !     INPUT
         integer,intent(in):: oufor           ! output unit for output force data

         integer,intent(in):: current_step    ! current time step

         !---- base value for non-dimensionalize
         real(8),intent(in):: fref             ! force base value [N]

         !---- force of each atoms
         real(8),intent(in):: for_long(:,:)     ! long-range force
         !  real(8),intent(in):: for_med(:,:)   ! medium-range force
         real(8),intent(in):: for_short(:,:)    ! short-range force

         !---- variables for calculation of heat flux
         logical,intent(in):: ifhfvol        ! local volume-based or local surface-based

         integer,intent(in):: nhfregion       ! number of region to calculate heat flux
         real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                ! z-position of region for heat flux

         integer,intent(in):: hftyp_atm(:)    ! atom- or mole-based heat flux cal.
                                !   for each atom

         real(8),intent(in):: molecom(:,:)     ! center of mass of molecule

         !     virial tensor for each atom
         real(8),intent(in):: viribot_long_atm(:,:,:,:) ! long-range (bond)
         real(8),intent(in):: viribot_short_atm(:,:,:,:) ! short-range (bond)
         real(8),intent(in):: viriant_long_atm(:,:,:,:) ! long-range (angle)
         real(8),intent(in):: viriant_short_atm(:,:,:,:) ! short-range (angle)
         real(8),intent(in):: viritot_long_atm(:,:,:,:) ! long-range (torsion)
         real(8),intent(in):: viritot_short_atm(:,:,:,:) ! short-range (torsion)
         real(8),intent(in):: viri14t_long_atm(:,:,:,:) ! long-range (1-4)
         real(8),intent(in):: viri14t_short_atm(:,:,:,:) ! short-range (1-4)
         real(8),intent(in):: virielint_long_atm(:,:,:,:) ! long-range (elc intra)
         real(8),intent(in):: virielint_short_atm(:,:,:,:) ! short-range (elc intra)
         real(8),intent(in):: viriljint_long_atm(:,:,:,:) ! long-range (L-J intra)
         real(8),intent(in):: viriljint_short_atm(:,:,:,:) ! short-range (L-J intra)
         real(8),intent(in):: virielct_long_atm(:,:,:,:) ! long-range (Coulomb)
         real(8),intent(in):: virielct_short_atm(:,:,:,:) ! short-range (Coulomb)
         real(8),intent(in):: viriljt_long_atm(:,:,:,:) ! long-range (L-J)
         real(8),intent(in):: viriljt_short_atm(:,:,:,:) ! short-range (L-J)
         real(8),intent(in):: virimort_long_atm(:,:,:,:) ! long-range (Morse)
         real(8),intent(in):: virimort_short_atm(:,:,:,:) ! short-range (Morse)
         real(8),intent(in):: virisht_long_atm(:,:,:,:) ! long-range (SH)
         real(8),intent(in):: virisht_short_atm(:,:,:,:) ! short-range (SH)
         real(8),intent(in):: virirfht_long_atm(:,:,:,:) ! long-range (RFH)
         real(8),intent(in):: virirfht_short_atm(:,:,:,:) ! short-range (RFH)
         real(8),intent(in):: viridout_long_atm(:,:,:,:) ! long-range (DOU)
         real(8),intent(in):: viridout_short_atm(:,:,:,:) ! short-range (DOU)
         real(8),intent(in):: viricstmnbt_long_atm(:,:,:,:) ! long-range (custom NB)
         real(8),intent(in):: viricstmnbt_short_atm(:,:,:,:) ! short-range (custom NB)
         real(8),intent(in):: viricstmnbext_long_atm(:,:,:,:,:)
                                                  ! long-range (custom NB)
         real(8),intent(in):: viricstmnbext_short_atm(:,:,:,:,:)
                                                  ! short-range (custom NB)

     end subroutine outfor_hf

     subroutine outthe(outhe,current_step,timeref, &
          &             mchain,xlogs,vlogs)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: outhe           ! output unit for output thermostat data

       integer,intent(in):: current_step    ! current time step

       real(8),intent(in):: timeref          ! time base value [sec]

       !---- thermostat data
       integer,intent(in):: mchain          ! Nose-Hoover chain number (>1 in NVT)
       real(8),intent(in):: xlogs(:)         ! xi of the thermostat coordinate
       real(8),intent(in):: vlogs(:)         ! vxi of the thermostat velocity

     end subroutine outthe

     subroutine outbar(oubar,current_step,xref,timeref, &
         &             xcel,ycel,zcel,xlogv,vlogv,xboxh,vboxg, &
         &             pcont_axis)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: oubar          ! output unit for output barostat data

       integer,intent(in):: current_step    ! current time step

       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: timeref          ! time base value [sec]

       !---- barostat data
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: xlogv            ! epsilon of the barostat coordinate
       real(8),intent(in):: vlogv            ! vepsilon of the barostat velocity

       real(8),intent(in):: xboxh(:)         ! cell vector h of the barostat coordinate
       real(8),intent(in):: vboxg(:)         ! vg vector of the barostat coordinate

       character(5),intent(in):: pcont_axis  ! axis for pressure control (iso, aniso, etc.)

     end subroutine outbar

     subroutine outpre(oupre,current_step,pref,eref, &
          &            pressmol_ktot,pressmolt_ktot, &
          &            pressmol_vtot, pressmolt_vtot, &
          &            pressmol_tot,pressmolt_tot, &
          &            pressatm_ktot,pressatmt_ktot, &
          &            pressatm_vtot,pressatmt_vtot, &
          &            pressatm_tot,pressatmt_tot, &
          &            pot_viric_all,pot_virict_all, &
          &            pot_virilj_all,pot_viriljt_all, &
          &            pot_virimor_all,pot_virimort_all, &
          &            pot_virish_all,pot_virisht_all, &
          &            pot_virirfh_all,pot_virirfht_all, &
          &            pot_viridou_all,pot_viridout_all, &
          &            pot_viricstmnb_all,pot_viricstmnbt_all, &
          &            ncstmnbex, &
          &            pot_viricstmnbex_all,pot_viricstmnbext_all, &
          &            pot_viri_all,pot_virit_all, &
          &            viri_fdotd,virit_fdotd, &
          &            atm_viribo_all,atm_viribot_all, &
          &            atm_virian_all,atm_viriant_all, &
          &            atm_virito_all,atm_viritot_all, &
          &            atm_viri14_all,atm_viri14t_all, &
          &            atm_viri_const,atm_virit_const, &
          &            atm_viri_corr,atm_virit_corr)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: oupre           ! output unit for outpre velocity data

       integer,intent(in):: current_step    ! current time step

       !---- base value for non-dimensionalize
       real(8),intent(in):: pref             ! pressure base value [Pa]
       real(8),intent(in):: eref             ! energy base value [J]

       !---- valiables for pressure
       real(8),intent(in):: pressmol_ktot       ! pressure kinetic part of molecule
       real(8),intent(in):: pressmolt_ktot(:,:)
                    ! pressure tensor kinetic part of molecule
       real(8),intent(in):: pressmol_vtot       ! pressure virial part of molecule
       real(8),intent(in):: pressmolt_vtot(:,:)
                    ! pressure tensor virial part of molecule
       real(8),intent(in):: pressmol_tot        ! pressure total of molecule
       real(8),intent(in):: pressmolt_tot(:,:)  ! pressure tensor total of molecule
       real(8),intent(in):: pressatm_ktot       ! pressure kinetic part of atm
       real(8),intent(in):: pressatmt_ktot(:,:) ! pressure tensor kinetic part of atm
       real(8),intent(in):: pressatm_vtot       ! pressure virial part of atm
       real(8),intent(in):: pressatmt_vtot(:,:) ! pressure tensor virial part of atm
       real(8),intent(in):: pressatm_tot        ! pressure total of atm
       real(8),intent(in):: pressatmt_tot(:,:)  ! pressure tensor total of atm

       real(8),intent(in):: pot_viric_all        ! virial(coulomb potential)
       real(8),intent(in):: pot_virict_all(:,:)  ! virial tensor (coulomb potential)
       real(8),intent(in):: pot_virilj_all       ! virial(L-J potential)
       real(8),intent(in):: pot_viriljt_all(:,:) ! virial tensor (L-J potential)
       real(8),intent(in):: pot_virimor_all       ! virial(Morse potential)
       real(8),intent(in):: pot_virimort_all(:,:) ! virial tensor (Morse potential)
       real(8),intent(in):: pot_virish_all       ! virial(SH potential)
       real(8),intent(in):: pot_virisht_all(:,:) ! virial tensor (SH potential)
       real(8),intent(in):: pot_virirfh_all       ! virial(RFH potential)
       real(8),intent(in):: pot_virirfht_all(:,:) ! virial tensor (RFH potential)
       real(8),intent(in):: pot_viridou_all       ! virial(DOU potential)
       real(8),intent(in):: pot_viridout_all(:,:) ! virial tensor (DOU potential)
       real(8),intent(in):: pot_viricstmnb_all    ! virial(custom NB potential)
       real(8),intent(in):: pot_viricstmnbt_all(:,:)
                                          ! virial tensor (custom NB potential)
       real(8),intent(in):: pot_viricstmnbex_all(:)
                                     ! virial(extra custom NB potential)
       real(8),intent(in):: pot_viricstmnbext_all(:,:,:)
                                     ! virial tensor (extra custom NB potential)

       real(8),intent(in):: pot_viri_all         ! all virial of potential term
       real(8),intent(in):: pot_virit_all(:,:)   ! all virial tensor of potential term
       real(8),intent(in):: viri_fdotd           ! virial of correction term F.d
       real(8),intent(in):: virit_fdotd(:,:)    ! virial tensor of correction term F.d

       real(8),intent(in):: atm_viribo_all       ! virial(bond potential) of each atom
       real(8),intent(in):: atm_viribot_all(:,:) ! virial(bond potential) of each atom
       real(8),intent(in):: atm_virian_all      ! virial(angle potential) of each atom
       real(8),intent(in):: atm_viriant_all(:,:)
                    ! virial(angle potential) of each atom
       real(8),intent(in):: atm_virito_all    ! virial(torsion potential) of each atom
       real(8),intent(in):: atm_viritot_all(:,:)
                    ! virial(torsion potential) of each atom
       real(8),intent(in):: atm_viri14_all  ! virial(1-4 force potential) of each atom
       real(8),intent(in):: atm_viri14t_all(:,:)
                    ! virial(1-4 force potential) of each atom

       real(8),intent(in):: atm_viri_const     ! virial(constraint force) of each atom
       real(8),intent(in):: atm_virit_const(:,:)
                    ! virial tensor (constraint) of each atom

       real(8),intent(in):: atm_viri_corr          ! virial(L-J correction)
       real(8),intent(in):: atm_virit_corr(:,:)    ! virial tensor (L-J correction)

       integer,intent(in):: ncstmnbex           ! number of extra custom NB output

     end subroutine outpre

     subroutine outthc(outhc,current_step, &
          &            xref,eref,timeref, &
          &            ntcregion, &
          &            det_ene_kin)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: outhc           ! output unit for outthc thermal control data

       integer,intent(in):: current_step    ! current time step

       integer,intent(in):: ntcregion       ! number of region to control temp.

       real(8),intent(in):: det_ene_kin(:)  ! for outputting thermal control data

       !---- base value for non-dimensionalize
       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: eref             ! energy base value [J]

       real(8),intent(in):: timeref          ! time base value [sec]

     end subroutine outthc

     subroutine outpdb(oupdb,current_step,   &
          &            npolytyp,npoly_mole,npoly_atom,   &
          &            nwater,   &
          &            nmatyp,nmatomtyp,   &
          &            resname_poly_pdb,   &
          &            resname_water_pdb,   &
          &            resname_matom_pdb)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: oupdb           ! output unit for outpdb PDB data

       integer,intent(in):: current_step    ! current time step

       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       character(4),intent(in):: resname_poly_pdb(:) ! residue name for poly
       character(4),intent(in):: resname_water_pdb ! residue name for water
       character(4),intent(in):: resname_matom_pdb(:) ! residue name for matom

     end subroutine outpdb

     subroutine outumb( ouumb, &
          &             xref,eref )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: ouumb    ! output unit for outumb bias potential data

       !---- base value for non-dimensionalize
       real(8),intent(in):: xref             ! distanse base value [m]

       real(8),intent(in):: eref             ! energy base value [J]

     end subroutine outumb

     subroutine outene_em(ouene,current_step,xref,eref,tempref,fref, &
          &               npolytyp,nmatyp, &
          &               pot_tot,pot_nonbon,pot_vdw, &
          &               pot_elc,pot_ewk,pot_ewc,pot_ewnq, &
          &               pot_vdw14,pot_elc14, &
          &               pot_bond,pot_angl,pot_anglub, &
          &               pot_tors,pot_torsrb,pot_torsim, &
          &               pot_mor, &
          &               pot_sh, &
          &               pot_rfh, &
          &               pot_dou, &
          &               pot_rpvw,pot_vw, &
          &               pot_cstmnb, &
          &               ncstmnbex,pot_cstmnbex, &
          &               pot_posres, &
          &               pot_pbias, &
          &               ene_kin_poly,ene_kin_water, &
          &               ene_kin_ma,ene_kin_all, &
          &               ene_tot, &
          &               temp_poly,temp_water,temp_ma,temp_all, &
          &               ene_kin_th,ene_pot_th, &
          &               ene_kin_ba,ene_pot_ba, &
          &               ene_conserve, &
          &               d_r,d_ene,drms)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: ouene           ! output unit for output energy data

       integer,intent(in):: current_step    ! current time step

       !---- base value for non-dimensionalize
       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: eref             ! energy base value [J]

       real(8),intent(in):: tempref          ! temperature base value [K]

       real(8),intent(in):: fref             ! force base value [N]

       !---- number of type of poly and matom
       integer,intent(in):: npolytyp        ! number of poly type

       integer,intent(in):: nmatyp          ! number of species of monatomic mole.

       !---- potential valiables
       real(8),intent(inout):: pot_tot          ! all potential
       real(8),intent(inout):: pot_nonbon       ! all non-bonded potential
       real(8),intent(inout):: pot_vdw          ! vdw potential
       real(8),intent(inout):: pot_elc          ! coulomb potential(ewald real)
       real(8),intent(inout):: pot_ewk          ! coulomb potential(ewald wave)
       real(8),intent(in):: pot_ewc          ! coulomb potential(ewald self)
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
       real(8),intent(inout):: pot_rpvw        ! RP-VW interaction
       real(8),intent(inout):: pot_vw          ! potential of constant force
       real(8),intent(inout):: pot_cstmnb       ! Custom NB potential
       real(8),intent(inout):: pot_cstmnbex(:)  ! extra Custom NB potential
       real(8),intent(inout):: pot_posres       ! position restraint potential
       real(8),intent(inout):: pot_pbias        ! bias potential

       integer,intent(in):: ncstmnbex           ! number of extra custom NB output

       !---- kinematic energy
       real(8),intent(inout):: ene_kin_poly(:)  ! kinetic energy
       real(8),intent(inout):: ene_kin_water    ! kinetic energy
       real(8),intent(inout):: ene_kin_ma(:)    ! kinetic energy
       real(8),intent(inout):: ene_kin_all      ! kinetic energy

       real(8),intent(inout):: ene_kin_th       ! kinetic energy of NHC thermostat
       real(8),intent(inout):: ene_pot_th       ! potential energy of NHC thermostat

       real(8),intent(inout):: ene_kin_ba       ! kinetic energy of Andersen barostat
       real(8),intent(inout):: ene_pot_ba      ! potential energy of Andersen barostat

       !---- all energy

       real(8),intent(inout):: ene_tot          ! total Hamiltonian
       real(8),intent(inout):: ene_conserve     ! total conserved quantity

       !---- temperature

       real(8),intent(inout):: temp_poly(:)     ! temperature of polymer1
       real(8),intent(inout):: temp_water       ! temperature of H2O
       real(8),intent(inout):: temp_ma(:)       ! temperature of monatomic mole.
       real(8),intent(inout):: temp_all         ! temperature of all atoms

       !---- variables for EM
       real(8),intent(in):: d_r                 ! delta r
       real(8),intent(in):: d_ene               ! delta pot_tot
       real(8),intent(in):: drms                ! root mean square force / natom

     end subroutine outene_em

#if defined(HF)
     subroutine outhtf(ouhtf,current_step, &
          &            xref,eref,timeref, &
          &            ifhfvol, &
          &            nhfregion,hfzpos1,hfzpos2, &
          &            hfkin_atm,hfkinex_atm, &
          &            hfpotbo_atm,hfpotan_atm, &
          &            hfpotto_atm,hfpot14_atm, &
          &            hfpotelin_atm,hfpotljin_atm, &
          &            hfpotelc_atm,hfpotlj_atm, &
          &            hfpotmor_atm,hfpotsh_atm, &
          &            hfpotrfh_atm,hfpotdou_atm, &
          &            hfpotcstmnb_atm, &
          &            ncstmnbex,hfpotcstmnbex_atm, &
          &            hfviribo_atm,hfvirian_atm, &
          &            hfvirito_atm,hfviri14_atm, &
          &            hfvirielin_atm,hfviriljin_atm, &
          &            hfvirielc_atm,hfvirilj_atm, &
          &            hfvirimor_atm,hfvirish_atm, &
          &            hfvirirfh_atm,hfviridou_atm, &
          &            hfviricstmnb_atm, &
          &            hfviricstmnbex_atm, &
          &            hftot_atm)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: ouhtf           ! output unit for outhtf heat flux data

       integer,intent(in):: current_step    ! current time step

       !---- base value for non-dimensionalize
       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: eref             ! energy base value [J]

       real(8),intent(in):: timeref          ! time base value [sec]

       !---- variables for calculation of heat flux
       logical,intent(in):: ifhfvol        ! local volume-based or local surface-based

       integer,intent(in):: nhfregion       ! number of region to calculate heat flux
       real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                ! z-position of region for heat flux

       !---- valiables for heat flux
       real(8),intent(in):: hfkin_atm(:,:)     ! heat flux of kinetic part
       real(8),intent(in):: hfkinex_atm(:,:,:)
                                 ! heat flux of kinetic part (extra custom NB)

       real(8),intent(in):: hfpotbo_atm(:,:)   ! heat flux of potential part (bond)
       real(8),intent(in):: hfpotan_atm(:,:)   ! heat flux of potential part (angle)
       real(8),intent(in):: hfpotto_atm(:,:)   ! heat flux of potential part (tors)
       real(8),intent(in):: hfpot14_atm(:,:)   ! heat flux of potential part (1-4)
       real(8),intent(in):: hfpotelin_atm(:,:)
                    ! heat flux of potential part (elc intra)
       real(8),intent(in):: hfpotljin_atm(:,:)
                    ! heat flux of potential part (L-J intra)

       real(8),intent(in):: hfpotelc_atm(:,:) ! heat flux of potential part (Coulomb)
       real(8),intent(in):: hfpotlj_atm(:,:) ! heat flux of potential part (L-J)
       real(8),intent(in):: hfpotmor_atm(:,:) ! heat flux of potential part (Morse)
       real(8),intent(in):: hfpotsh_atm(:,:) ! heat flux of potential part (SH)
       real(8),intent(in):: hfpotrfh_atm(:,:) ! heat flux of potential part (RFH)
       real(8),intent(in):: hfpotdou_atm(:,:) ! heat flux of potential part (DOU)
       real(8),intent(in):: hfpotcstmnb_atm(:,:)
                                     ! heat flux of potential part (custom NB)
       real(8),intent(in):: hfpotcstmnbex_atm(:,:,:)
                                ! heat flux of potential part (extra custom NB)

       real(8),intent(in):: hfviribo_atm(:,:)   ! heat flux of virial part (bond)
       real(8),intent(in):: hfvirian_atm(:,:)   ! heat flux of virial part (angle)
       real(8),intent(in):: hfvirito_atm(:,:)   ! heat flux of virial part (tors)
       real(8),intent(in):: hfviri14_atm(:,:)   ! heat flux of virial part (1-4)
       real(8),intent(in):: hfvirielin_atm(:,:) ! heat flux of virial part (elc intra)
       real(8),intent(in):: hfviriljin_atm(:,:) ! heat flux of virial part (L-J intra)

       real(8),intent(in):: hfvirielc_atm(:,:) ! heat flux of virial part (Coulomb)
       real(8),intent(in):: hfvirilj_atm(:,:) ! heat flux of virial part (L-J)
       real(8),intent(in):: hfvirimor_atm(:,:) ! heat flux of virial part (Morse)
       real(8),intent(in):: hfvirish_atm(:,:) ! heat flux of virial part (SH)
       real(8),intent(in):: hfvirirfh_atm(:,:) ! heat flux of virial part (RFH)
       real(8),intent(in):: hfviridou_atm(:,:) ! heat flux of virial part (DOU)
       real(8),intent(in):: hfviricstmnb_atm(:,:)
                                        ! heat flux of virial part (custom NB)
       real(8),intent(in):: hfviricstmnbex_atm(:,:,:)
                                  ! heat flux of virial part (extra custom NB)

       real(8),intent(in):: hftot_atm(:,:) ! total heat flux

       integer,intent(in):: ncstmnbex      ! number of extra custom NB output

     end subroutine outhtf

     subroutine outmtf(oumtf, &
          &            current_step, &
          &            xref,eref,pref,timeref, &
          &            ifhfvol, &
          &            nhfregion,hfzpos1,hfzpos2, &
          &            mfkin_atm,mfkinex_atm, &
          &            mfviribo_atm,mfvirian_atm, &
          &            mfvirito_atm,mfviri14_atm, &
          &            mfvirielin_atm,mfviriljin_atm, &
          &            mfvirielc_atm,mfvirilj_atm, &
          &            mfvirimor_atm,mfvirish_atm, &
          &            mfvirirfh_atm,mfviridou_atm, &
          &            mfviricstmnb_atm, &
          &            ncstmnbex,mfviricstmnbex_atm, &
          &            mftot_atm, &
          &            mtfoutdir)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: oumtf           ! output unit for outmtf momentum flux data

       integer,intent(in):: current_step    ! current time step

       !---- base value for non-dimensionalize
       real(8),intent(in):: xref             ! distanse base value [m]
       real(8),intent(in):: eref             ! energy base value [J]
       real(8),intent(in):: pref             ! pressure base value [J]

       real(8),intent(in):: timeref          ! time base value [sec]

       !---- variables for calculation of transport flux
       logical,intent(in):: ifhfvol        ! local volume-based or local surface-based

       integer,intent(in):: nhfregion       ! number of region to calculate heat flux
       real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                ! z-position of region for heat flux

       !---- variables for heat flux
       real(8),intent(in):: mfkin_atm(:,:,:)     ! momentum flux of kinetic part
       real(8),intent(in):: mfkinex_atm(:,:,:,:)
                              ! momentum flux of kinetic part (extra custom NB)

       real(8),intent(in):: mfviribo_atm(:,:,:)   ! momentum flux of virial part (bond)
       real(8),intent(in):: mfvirian_atm(:,:,:)   ! momentum flux of virial part (angle)
       real(8),intent(in):: mfvirito_atm(:,:,:)
                     ! momentum flux of virial part (tors)
       real(8),intent(in):: mfviri14_atm(:,:,:)
                     ! momentum flux of virial part (1-4)
       real(8),intent(in):: mfvirielin_atm(:,:,:)
                     ! momentum flux of virial part (elc intra)
       real(8),intent(in):: mfviriljin_atm(:,:,:)
                     ! momentum flux of virial part (L-J intra)

       real(8),intent(in):: mfvirielc_atm(:,:,:)
                     ! momentum flux of virial part (Coulomb)
       real(8),intent(in):: mfvirilj_atm(:,:,:)
                     ! momentum flux of virial part (L-J)
       real(8),intent(in):: mfvirimor_atm(:,:,:)
                     ! momentum flux of virial part (Morse)
       real(8),intent(in):: mfvirish_atm(:,:,:)
                     ! momemtum flux of virial part (SH)
       real(8),intent(in):: mfvirirfh_atm(:,:,:)
                     ! momentum flux of virial part (RFH)
       real(8),intent(in):: mfviridou_atm(:,:,:)
                     ! momentum flux of virial part (DOU)
       real(8),intent(in):: mfviricstmnb_atm(:,:,:)
                     ! momentum flux of virial part (custom NB)
       real(8),intent(in):: mfviricstmnbex_atm(:,:,:,:)
                     ! momentum flux of virial part (extra custom NB)

       real(8),intent(in):: mftot_atm(:,:,:) ! total momentum flux

       integer,intent(in):: ncstmnbex        ! number of extra custom NB output

       ! variables for outputting momentum flux
       character(3),intent(in):: mtfoutdir ! direction to output data of momentum

     end subroutine outmtf
#endif

     subroutine nhcint( degfree_all, &
          &             dt_long_cal, &
          &             mchain,text_c, &
          &             next,nyosh )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: degfree_all     ! degree of freedom of all atoms

       real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]

       integer,intent(in):: mchain         ! Nose-Hoover chain number (>1 in NVT)
       real(8),intent(in):: text_c       ! external temp. [K] (Nose-Hoover chain)

       integer,intent(in):: next            ! iteration number of extended system
       integer,intent(in):: nyosh      ! expansion order of Yoshida-Suzuki method

     end subroutine nhcint

     subroutine nhcpisoint( degfree_all, &
          &                 dt_long_cal, &
          &                 mchain,text_c, &
          &                 next,nyosh, &
          &                 pint,pext_c )

       ! ARGUMENT:
       !     INPUT
       integer,intent(in):: degfree_all     ! degree of freedom of all atoms

       real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]

       integer,intent(in):: mchain         ! Nose-Hoover chain number (>1 in NVT)
       real(8),intent(in):: text_c       ! external temp. [K] (Nose-Hoover chain)

       integer,intent(in):: next            ! iteration number of extended system
       integer,intent(in):: nyosh      ! expansion order of Yoshida-Suzuki method

       real(8),intent(in):: pint             ! internal pressure in NPT(MTK)
       real(8),intent(in):: pext_c           ! external pressure in NPT(MTK)

     end subroutine nhcpisoint

     subroutine nhcpanisoint(degfree_all, &
         &                  dt_long_cal, &
         &                  mchain,text_c, &
         &                  next,nyosh, &
         &                  pintt,pext_c,pcont_dir)

         ! ARGUMENT:
         !     INPUT
         integer,intent(in):: degfree_all     ! degree of freedom of all atoms

         real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]

         integer,intent(in):: mchain          ! Nose-Hoover chain number (>1 in NVT)
         real(8),intent(in):: text_c           ! external temp. [K] (Nose-Hoover chain)

         integer,intent(in):: next            ! iteration number of extended system
         integer,intent(in):: nyosh           ! expansion order of Yoshida-Suzuki method

         real(8),intent(in):: pintt(:,:)       ! internal pressure tensor in NPT (MTK)
         real(8),intent(in):: pext_c           ! external pressure in NPT(MTK)
         integer,intent(in):: pcont_dir        ! direction for uniaxial NPT control

     end subroutine nhcpanisoint

     subroutine calvol( xref,current_step, &
          &             xcel,ycel,zcel, &
          &             yratio,zratio, &
          &             xlogv, &
          &             rcut_book,ifbook,rcut, &
          &             rcutmor,ifbookmor,rcut_bookmor)

       ! ARGUMENT:
       !     INPUT
       real(8),intent(in):: xref             ! distanse base value [m]
       integer,intent(in):: current_step    ! current time step

       real(8),intent(in):: yratio           ! y cell ratio of y to x
       real(8),intent(in):: zratio           ! z cell ratio of z to x

       real(8),intent(in):: xlogv            ! epsilon of the barostat coordinate

       real(8),intent(in):: rcut_book        ! cut off radius of bookkeeping
       logical,intent(in):: ifbook          ! flag for bookkeeping

       real(8),intent(in):: rcut             ! vdw cutoff length [non-d]

       real(8),intent(in):: rcutmor          ! Morse cutoff length [m]

       logical,intent(in):: ifbookmor ! flag for bookkeeping of Morse interaction
       real(8),intent(in):: rcut_bookmor
                         ! cut off radius of bookkeeping[m] of Morse

       !     INPUT&OUTPUT
       real(8),intent(inout):: xcel             ! x cell length[non-d]
       real(8),intent(inout):: ycel             ! y cell length[non-d]
       real(8),intent(inout):: zcel             ! z cell length[non-d]

     end subroutine calvol

     subroutine calbox(xref,current_step, &
         &            xcel,ycel,zcel, &
         &            yratio,zratio, &
         &            xboxh, &
         &            rcut_book,ifbook,rcut, &
         &            rcutmor,ifbookmor,rcut_bookmor)

         ! ARGUMENT:
         !     INPUT
         real(8),intent(in):: xref             ! distanse base value [m]
         integer,intent(in):: current_step    ! current time step

         real(8),intent(inout):: yratio        ! y cell ratio of y to x
         real(8),intent(inout):: zratio        ! z cell ratio of z to x

         real(8),intent(in):: xboxh(3)      ! cell vector h of the barostat coordinate

         real(8),intent(in):: rcut_book        ! cut off radius of bookkeeping
         logical,intent(in):: ifbook           ! flag for bookkeeping

         real(8),intent(in):: rcut             ! vdw cutoff length [non-d]

         real(8),intent(in):: rcutmor          ! Morse cutoff length [m]

         logical,intent(in):: ifbookmor        ! flag for bookkeeping of Morse interaction
         real(8),intent(in):: rcut_bookmor     ! cut off radius of bookkeeping[m] of Morse

         !     INPUT&OUTPUT
         real(8),intent(inout):: xcel             ! x cell length[non-d]
         real(8),intent(inout):: ycel             ! y cell length[non-d]
         real(8),intent(inout):: zcel             ! z cell length[non-d]

     end subroutine calbox

     subroutine calcom( npoly,nwater,nmatom, &
          &             xcel,ycel,zcel, &
          &             molecom)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: npoly           ! number of polymer1
       integer,intent(in):: nwater          ! number of H2O molecules
       integer,intent(in):: nmatom          ! number of monatomic molecules

       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       !     OUTPUT
       real(8),intent(out):: molecom(:,:)     ! center of mass of molecule

     end subroutine calcom

#if defined(HF)
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

     end subroutine calheatf

     subroutine calmomentumf(npoly,nwater,nmatom, &
          &                  xcel,ycel,zcel, &
          &                  ifhfvol, &
          &                  nhfregion,hfzpos1,hfzpos2, &
          &                  hftyp_atm, &
          &                  molecom, &
          &                  viribot_long_atm,viribot_med_atm, &
          &                  viribot_short_atm, &
          &                  viriant_long_atm,viriant_med_atm, &
          &                  viriant_short_atm, &
          &                  viritot_long_atm,viritot_med_atm, &
          &                  viritot_short_atm, &
          &                  viri14t_long_atm,viri14t_med_atm, &
          &                  viri14t_short_atm, &
          &                  virielint_long_atm,virielint_med_atm, &
          &                  virielint_short_atm, &
          &                  viriljint_long_atm,viriljint_med_atm, &
          &                  viriljint_short_atm, &
          &                  virielct_long_atm,virielct_med_atm, &
          &                  virielct_short_atm, &
          &                  viriljt_long_atm,viriljt_med_atm, &
          &                  viriljt_short_atm, &
          &                  virimort_long_atm,virimort_med_atm, &
          &                  virimort_short_atm, &
          &                  virisht_long_atm,virisht_med_atm, &
          &                  virisht_short_atm, &
          &                  virirfht_long_atm,virirfht_med_atm, &
          &                  virirfht_short_atm, &
          &                  viridout_long_atm,viridout_med_atm, &
          &                  viridout_short_atm, &
          &                  viricstmnbt_long_atm,viricstmnbt_med_atm, &
          &                  viricstmnbt_short_atm, &
          &                  viricstmnbext_long_atm,viricstmnbext_med_atm, &
          &                  viricstmnbext_short_atm, &
          &                  mfkin_atm,mfkinex_atm, &
          &                  mfviribo_atm,mfvirian_atm,mfvirito_atm, &
          &                  mfviri14_atm, &
          &                  mfvirielin_atm,mfviriljin_atm, &
          &                  mfvirielc_atm,mfvirilj_atm, &
          &                  mfvirimor_atm,mfvirish_atm, &
          &                  mfvirirfh_atm,mfviridou_atm, &
          &                  mfviricstmnb_atm, &
          &                  mfviricstmnbex_atm, &
          &                  ncstmnbex, &
          &                  mftot_atm, &
          &                  atmcor_mf,molecom_mf, &
          &                  dt_cross, &
          &                  heatfinterval, &
          &                  current_step)

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

       integer,intent(in):: nhfregion
                                ! number of region to calculate momemtum flux
       real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                ! z-position of region for momentum flux

       integer,intent(in):: hftyp_atm(:)    ! atom- or mole-based momentum flux cal.
                                !   for each atom

       real(8),intent(in):: molecom(:,:)     ! center of mass of molecule

       real(8),intent(in):: dt_cross ! interval time to check if cross the surface
                                !   only for surface-based method
                                ! = dt_long_cal * outinterval

       integer,intent(in):: heatfinterval   ! interval of heatf output

       integer,intent(in):: current_step    ! current time step

       integer,intent(in):: ncstmnbex       ! number of extra custom NB output

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
       !     for momentum flux sum
       real(8),intent(out):: mfkin_atm(:,:,:) ! momentum flux of kinetic part
       real(8),intent(out):: mfkinex_atm(:,:,:,:)
                              ! momentum flux of kinetic part (extra custom NB)

       real(8),intent(out):: mfviribo_atm(:,:,:)   ! momentum flux of virial part (bond)
       real(8),intent(out):: mfvirian_atm(:,:,:)   ! momentum flux of virial part (angle)
       real(8),intent(out):: mfvirito_atm(:,:,:)   ! momentum flux of virial part (tors)
       real(8),intent(out):: mfviri14_atm(:,:,:)   ! momentum flux of virial part (1-4)
       real(8),intent(out):: mfvirielin_atm(:,:,:)
                     ! momentum flux of virial part (elc intra)
       real(8),intent(out):: mfviriljin_atm(:,:,:)
                     ! momentum flux of virial part (L-J intra)

       real(8),intent(out):: mfvirielc_atm(:,:,:) ! momentum flux of virial part (Coulomb)
       real(8),intent(out):: mfvirilj_atm(:,:,:) ! momentum flux of virial part (L-J)
       real(8),intent(out):: mfvirimor_atm(:,:,:) ! momentum flux of virial part (Morse)
       real(8),intent(out):: mfvirish_atm(:,:,:) ! momentum flux of virial part (SH)
       real(8),intent(out):: mfvirirfh_atm(:,:,:) ! momentum flux of virial part (RFH)
       real(8),intent(out):: mfviridou_atm(:,:,:) ! momentum flux of virial part (DOU)
       real(8),intent(out):: mfviricstmnb_atm(:,:,:)
                                ! momentum flux of virial part (custom NB)
       real(8),intent(out):: mfviricstmnbex_atm(:,:,:,:)
                                ! momentum flux of virial part (extra custom NB)

       real(8),intent(out):: mftot_atm(:,:,:)   ! total momentum flux

       real(8),intent(out):: atmcor_mf(:,:)
                     ! old atmcor for momentum flux calculation
       real(8),intent(out):: molecom_mf(:,:)
                     ! old center of mass of molecule (momentum flux)

     end subroutine calmomentumf
#endif

  end interface

end module interface_mdtech
