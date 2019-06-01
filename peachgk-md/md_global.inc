!*****************************
!*  md_global.inc            *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <Feb 07 2017>

module md_global

  implicit none

#include "config.h"

  !---- peachgk.ini version description
  character(3),parameter:: pgkiniver = _PGKINIVER_  ! peachgk.ini Version num.
  character(8),parameter:: pgkinidate = _PGKINIDATE_   ! peachgk.ini rev. date

  !---- maximum number of arrays
  !     - declare variables for constant parameter
  integer,parameter:: maxnstage = _MAXNSTAGE_   ! maximum number of MD stage
  integer,parameter:: maxnproc = _MAXNPROC_
                                       ! maximum number of processes for MPI run

  integer,parameter:: maxninttyp = _MAXNINTTYP_   ! maximum interaction type

  integer,parameter:: maxnatom = _MAXNATOM_   ! maximum number of atoms

  integer,parameter:: maxnmol = _MAXNMOL_   ! maximum number in each molecule
  integer,parameter:: maxnbond = _MAXNBOND_
  integer,parameter:: maxnangl = _MAXNANGL_
  integer,parameter:: maxntors = _MAXNTORS_

  integer,parameter:: maxnconst = _MAXNCONST_
  integer,parameter:: maxnlist = _MAXNLIST_
                                ! maximum number of atom per one atom in nblist
  integer,parameter:: maxmorlist = _MAXMORLIST_
                                ! maximum number of atom per one atom in morlist
  integer,parameter:: maxshlist = _MAXSHLIST_
                                ! maximum number of atom per one atom in SHlist
  integer,parameter:: maxrfhlist = _MAXRFHLIST_
                                ! maximum number of atom per one atom in RFHlist
  integer,parameter:: maxdoulist = _MAXDOULIST_
                                ! maximum number of atom per one atom in DOUlist
  integer,parameter:: maxrpvwlist = _MAXRPVWLIST_
                                ! maximum number of atom per one atom in RP-VWlist
  integer,parameter:: maxcstmnblist = _MAXCSTMNBLIST_
                                ! maximum number of atom per one atom in cstmnblist

  integer,parameter:: maxnvdwtyp = _MAXNVDWTYP_
                                ! maximum number of types of vdW interaction
  integer,parameter:: maxnatmtyp = _MAXNATMTYP_
                                ! maximum number of atom types
  integer,parameter:: maxnbondtyp = _MAXNBONDTYP_
  integer,parameter:: maxnangltyp = _MAXNANGLTYP_
  integer,parameter:: maxntorstyp = _MAXNTORSTYP_
  integer,parameter:: maxnmortyp = _MAXNMORTYP_
  integer,parameter:: maxnshtyp = _MAXNSHTYP_
  integer,parameter:: maxnrfhtyp = _MAXNRFHTYP_
  integer,parameter:: maxndoutyp = _MAXNDOUTYP_
  integer,parameter:: maxnrpvwtyp = _MAXNRPVWTYP_
  integer,parameter:: maxncstmnbtyp = _MAXNCSTMNBTYP_
  integer,parameter:: maxnconsttyp = _MAXNCONSTTYP_
                                ! maximum number of constraint pair type

  integer,parameter:: maxnpolytyp = _MAXNPOLYTYP_
                                ! maximum number of polymer type species
  integer,parameter:: maxnmatyp = _MAXNMATYP_
                                ! maximum number of monatomic type species

  integer,parameter:: maxwave = _MAXWAVE_   ! =(kmax*2+1)**3/2

  integer,parameter:: maxnhc = _MAXNHC_   ! max num. of Nose-Hoover chain
  integer,parameter:: maxnyosh = _MAXNYOSH_   ! max num. of Yoshida-Suzuki method

  integer,parameter:: max_cellx = _MAXCELLX_
                                        ! maximum number of cell-index cell (x)
  integer,parameter:: max_celly = _MAXCELLY_
                                        ! maximum number of cell-index cell (y)
  integer,parameter:: max_cellz = _MAXCELLZ_
                                        ! maximum number of cell-index cell (z)
  integer,parameter:: max_cellindex = _MAXCELLX_ + _MAXCELLX_*_MAXCELLY_   &
       &                            + _MAXCELLX_*_MAXCELLY_*_MAXCELLZ_
                                        ! maximum number of cell-index cell all
  integer,parameter:: max_atomcell = _MAX_ATOM_CELL_
                                        ! maximum number of atoms for each cell

  integer,parameter:: max_interpol = _MAX_INTERPOL_
                                 ! maximum number of spline interpolation points
  integer,parameter:: max_spl_order = _MAX_SPL_ORDER_
                                 ! maximum order of spline function

  integer,parameter:: maxntcregion = _MAXNTCREGION_
                                   ! maximum number of region for temp. control

  integer,parameter:: maxnword = _MAXNWORD_
                                       ! maximum number of words from data file

  integer,parameter:: maxnhfregion = _MAXNHFREGION_
                                 ! maximum number of region for heat flux calc.

  integer,parameter:: maxnintpbias = _MAXNINTPBIAS_
                                      ! maximum number of sampling according to
                                      ! reaction coordinate with bias potential
  integer,parameter:: maxnvelregion = _MAXNVELREGION_
                                     ! max num. of slabs for streaming velocity

  !---- define constant for MD control flags
  !     - declare variables for constant parameter
  integer,parameter:: MD_0K  = 1          ! 0[K] NVT for clear distorsion
  integer,parameter:: MD_H   = 2          ! gradual heating NVT (v-scale)
  integer,parameter:: MD_T   = 3          ! target temperature NVT (v-scale)
  integer,parameter:: MD_MTK = 4          ! NPT constant MD (MTK eq.)
  integer,parameter:: MD_NHC = 5          ! NVT constant MD (NHC eq.)
  integer,parameter:: MD_NVE = 6          ! NVE constant MD
  integer,parameter:: MD_HTF = 7          ! heat flux calculation in NVE MD
  integer,parameter:: MD_EMS = 8          ! EM by steepest descent (SD) method


  integer,parameter:: INTTYPE_VDW    = 1   ! pair potential type (vdW)
  integer,parameter:: INTTYPE_MOR    = 2   ! pair potential type (Morse)
  integer,parameter:: INTTYPE_SH     = 3   ! pair potential type (SH)
  integer,parameter:: INTTYPE_RFHFO  = 4   ! pair potential type (RFH F-O)
  integer,parameter:: INTTYPE_RFHOO  = 5   ! pair potential type (RFH O-O)
  integer,parameter:: INTTYPE_RFHOH  = 6   ! pair potential type (RFH O-H)
  integer,parameter:: INTTYPE_DOUO   = 7   ! pair potential type (DOU O-Au)
  integer,parameter:: INTTYPE_DOUH   = 8   ! pair potential type (DOU H-Au)
  integer,parameter:: INTTYPE_RPVW   = 9   ! pair potential type (RP-VW)
  integer,parameter:: INTTYPE_CSTMNB = 10  ! pair potential type (CUSTOM-NB)

  integer,parameter:: MTS_LONG  = 1       ! MTS flag for long   force
  integer,parameter:: MTS_MED   = 2       ! MTS flag for medium force
  integer,parameter:: MTS_SHORT = 3       ! MTS flag for short  force

  integer,parameter:: MTK_XAXIS = 1       ! anisotropic NPT in x direction
  integer,parameter:: MTK_YAXIS = 2       ! anisotropic NPT in y direction
  integer,parameter:: MTK_ZAXIS = 3       ! anisotropic NPT in y direction

  integer,parameter:: HFTYP_ATOM = 1      ! atom-based heat flux calculation
  integer,parameter:: HFTYP_MOLE = 2      ! mole-based heat flux calculation

  !***** declare global variables
  !---- coordinate, charge, topology, or some parameters for all atoms
  integer,save:: natom                   ! number of all atoms
  character(2),save:: atmtyp(maxnatom)   ! atom type
  integer,save:: atmindex(maxnatom)      ! index of atom type
  real(8),save:: atmcor(3,maxnatom)      ! coordinate of atoms
  real(8),save:: atmmass(maxnatom)       ! mass of atoms
  real(8),save:: atmchrg(maxnatom)       ! atomic charge

  integer,save:: nbond                             ! number of all bonds
  integer,save:: ibond(maxnbond),jbond(maxnbond)   ! bond topology
  character(5),save:: bondtyp(maxnbond)            ! bond type
  integer,save:: indexbondtyp(maxnbond)            ! bond type index

  integer,save:: nangl                                ! number of all angles
  integer,save:: iangl(maxnangl),jangl(maxnangl),kangl(maxnangl) ! angle top
  character(8),save:: angltyp(maxnangl)               ! angle type
  integer,save:: indexangltyp(maxnangl)               ! angle type index

  integer,save:: nanglub                           ! number of all UB angles
  integer,save:: ianglub(maxnangl),janglub(maxnangl),kanglub(maxnangl)
                                                   ! UB angle top
  character(8),save:: anglubtyp(maxnangl)          ! UB angle type
  integer,save:: indexanglubtyp(maxnangl)          ! UB angle type index

  integer,save:: ntors                             ! number of all torsions
  integer,save:: itors(maxntors),jtors(maxntors),   &
       &         ktors(maxntors),ltors(maxntors)   ! torsion topology
  character(11),save:: torstyp(maxntors)           ! torsion type
  integer,save:: indextorstyp(maxntors)            ! torsion type index
  real(8),save:: ring14(maxntors)       ! correction factor for 1-4 vdW and elc
                                        ! for ring structure
                                        ! 0.0 = within 4- or 5- membered ring
                                        ! 0.5 = within 6-membered ring
                                        ! 1.0 = others

  integer,save:: ntorsrb                              ! number of all RBtorsions
  integer,save:: itorsrb(maxntors),jtorsrb(maxntors),   &
       &         ktorsrb(maxntors),ltorsrb(maxntors)  ! RBtorsion topology
  character(11),save:: torsrbtyp(maxntors)            ! RBtorsion type
  integer,save:: indextorsrbtyp(maxntors)             ! RBtorsion type index
  real(8),save:: ring14_rb(maxntors)    ! correction factor for 1-4 vdW and elc
                                        ! for ring structure
                                        ! 0.0 = within 4- or 5- membered ring
                                        ! 0.5 = within 6-membered ring
                                        ! 1.0 = others

  integer,save:: ntorsim                       ! number of all improper torsions
  integer,save:: itorsim(maxntors),jtorsim(maxntors),   &
       &         ktorsim(maxntors),ltorsim(maxntors) ! improper torsion topology
  character(11),save:: torsimtyp(maxntors)     ! improper torsion type
  integer,save:: indextorsimtyp(maxntors)      ! improper torsion type index

  integer,save:: nmorse	                ! number of all Morse atoms
  integer,save:: nmorselist(maxnatom)   ! list of Morse atoms
  logical,save:: ifmorse(maxnatom)      ! if morse atom or not

  integer,save:: nsho	              ! number of all SH atoms (O atoms)
  integer,save:: nsholist(maxnatom)   ! list of SH atoms (O atoms)
  integer,save:: nshh	              ! number of all SH atoms (H atoms)
  integer,save:: nshhlist(maxnatom)   ! list of SH atoms (H atoms)

  integer,save:: nrfh	              ! number of all RFH atoms
  integer,save:: nrfhlist(maxnatom)   ! list of RFH atoms

  integer,save:: ndou	              ! number of all DOU atoms
  integer,save:: ndoulist(maxnatom)   ! list of DOU atoms

  integer,save:: nrpvw	               ! number of all RP-VW intaraction atoms
  integer,save:: nrpvwlist(maxnatom)   ! list of RP-VW interaction atoms

  integer,save:: nvwlist(2)   ! list of VW atoms (limited to 2 for both sides)

  integer,save:: nconst                              ! number of all contraint
  integer,save:: nconst_poly(maxnpolytyp)          ! number of polymer contraint
  integer,save:: nconst_water                        ! number of water contraint
  integer,save:: iconst(maxnconst),jconst(maxnconst) ! const topology
  real(8),save:: dconst(maxnconst)                   ! constraint length

  integer,save:: nmole                ! number of all molecules

  real(8),save:: molmass(maxnatom)    ! mass of molecule

  !---- variables for intermolecular & intramolecular potential
  integer,save:: inter_inttyp(maxnatmtyp,maxnatmtyp)
                                ! intermolecular interaction type
                                ! = 1,  INTTYPE_VDW
                                ! = 2,  INTTYPE_MOR
                                ! = 3,  INTTYPE_SH
                                ! = 4,  INTTYPE_RFHFO
                                ! = 5,  INTTYPE_RFHOO
                                ! = 6,  INTTYPE_RFHOH
                                ! = 7,  INTTYPE_DOUO
                                ! = 8,  INTTYPE_DOUH
                                ! = 9,  INTTYPE_RPVW
                                ! = 10, INTTYPE_CSTMNB

  integer,save:: nvdwtyp         ! number of vdw type
  character(2),save:: vdw_typ(maxnvdwtyp)
                                 ! atom pair of vdw
  character(1),save:: vdw_comb(maxnvdwtyp)   ! combination rule of vdw
  real(8),save:: vdw_welij(maxnvdwtyp,maxnvdwtyp) ! well depth of vdw parameter
  real(8),save:: vdw_radij(maxnvdwtyp,maxnvdwtyp) ! radius of vdw parameter

  integer,save:: natmtyp                      ! number of atom type
  character(2),save:: para_atmtyp(maxnatmtyp) ! name of atmtyp
  real(8),save:: para_atmmass(maxnatmtyp)     ! atomic mass [g/mol]
  real(8),save:: para_atmchrg(maxnatmtyp)     ! atom charge [non-d]

  integer,save:: nbondtyp                         ! number of bond type
  character(5),save:: para_bondtyp(maxnbondtyp)   ! bondtyp of bond pair
  real(8),save:: para_cbond(maxnbondtyp)       ! bond constant of bond parameter
  real(8),save:: para_eqbond(maxnbondtyp)         ! equilibrium bond length

  integer,save:: nangltyp        ! number of angle type
  character(8),save:: para_angltyp(maxnangltyp) ! angltyp of angle pair
  real(8),save:: para_cangl(maxnangltyp) ! angle constant of angle para
  real(8),save:: para_eqangl(maxnangltyp) ! equilibrium angle of angle para

  integer,save:: nanglubtyp      ! number of Urey-Bradley angle type
  character(8),save:: para_anglubtyp(maxnangltyp) ! angltyp of UB angle pair
  real(8),save:: para_canglub(maxnangltyp)  ! angle constant of UB angle para
  real(8),save:: para_eqanglub(maxnangltyp) ! equilibrium angle of UB angle para
  real(8),save:: para_cbondub(maxnbondtyp)  ! bond constant of UB bond parameter
  real(8),save:: para_eqbondub(maxnbondtyp) ! equilibrium UB bond length

  integer,save:: ntorstyp        ! number of torsion type
  character(11),save:: para_torstyp(maxntorstyp) ! torstyp of tors pair
  integer,save:: para_divfac(maxntorstyp) ! division facter of torsion parameter
  real(8),save:: para_barhig(maxntorstyp) ! barrier height of torsion parameter
  real(8),save:: para_phsang(maxntorstyp) ! phase angle of torsion parameter
  real(8),save:: para_period(maxntorstyp) ! periodicity of torsion parameter
  real(8),save:: para_divfac_vdw(maxntorstyp) ! dividing factor of 1-4 vdw
  real(8),save:: para_divfac_elc(maxntorstyp) ! dividing factor of 1-4 elc
  real(8),save:: para_wfac(maxntorstyp) ! weighting factor of 1-4 vdw & elc for CHARMM !!! by KAWAGOE
  logical,save:: ifcalnonb14(maxntorstyp) ! if calculate 1-4 interaction

  integer,save:: ntorsrbtyp                   ! number of RBtorsion type
  character(11),save:: para_torsrbtyp(maxntorstyp) ! RBtorstyp of tors pair
  real(8),save:: para_barhigrb(6,maxntorstyp)
                                        ! barrier height of RBtorsion parameter
  real(8),save:: para_divfac_vdwrb(maxntorstyp) ! dividing factor of 1-4 vdw (RB)
  real(8),save:: para_divfac_elcrb(maxntorstyp) ! dividing factor of 1-4 elc (RB)
  logical,save:: ifcalnonb14_rb(maxntorstyp)   ! if calculate 1-4 interaction

  integer,save:: ntorsimtyp      ! number of improper torsion type
  character(11),save:: para_torsimtyp(maxntorstyp)
                                 ! improper torstyp of tors pair
  real(8),save:: para_barhigim(maxntorstyp)
                                 ! barrier height of improper torsion parameter
  real(8),save:: para_phsangim(maxntorstyp)
                                 ! phase angle of improper torsion parameter

  integer,save:: nmortyp         ! number of Morse type
  character(5),save:: para_mortyp(maxnmortyp)
                                 ! atom pair of Morse type
  real(8),save:: para_welmor(maxnatmtyp,maxnatmtyp)
	                         ! well depth of Morse parameter
  real(8),save:: para_radmor(maxnatmtyp,maxnatmtyp)
	                         ! radius of Morse parameter
  real(8),save:: para_alphamor(maxnatmtyp,maxnatmtyp)
                                 ! well steepness of Morse parameter

  integer,save:: nshtyp          ! number of SH type (= 2 is expected)
  character(5),save:: para_shtyp(maxnshtyp)
  real(8),save:: para_sh_a(4)    ! a coefficient of SH potential
  real(8),save:: para_sh_b(4)    ! b coefficient of SH potential
  real(8),save:: para_sh_c       ! c coefficient of SH potential

  integer,save:: nrfhtyp          ! number of RFH type
  character(5),save:: para_rfhtyp(maxnrfhtyp)
  real(8),save:: para_rfhfo_a     ! a coefficient of RFH(FO) potential
  real(8),save:: para_rfhfo_b     ! b coefficient of RFH(FO) potential
  real(8),save:: para_rfhfo_c     ! c coefficient of RFH(FO) potential
  real(8),save:: para_rfhfo_d     ! d coefficient of RFH(FO) potential
  real(8),save:: para_rfhoo_a     ! a coefficient of RFH(OO) potential
  real(8),save:: para_rfhoo_b     ! b coefficient of RFH(OO) potential
  real(8),save:: para_rfhoo_c     ! c coefficient of RFH(OO) potential
  real(8),save:: para_rfhoo_d     ! d coefficient of RFH(OO) potential
  real(8),save:: para_rfhoo_lja   ! lja coefficient of RFH(OO) potential
  real(8),save:: para_rfhoo_ljb   ! ljb coefficient of RFH(OO) potential
  real(8),save:: para_rfhoh_a     ! a coefficient of RFH(OH) potential
  real(8),save:: para_rfhoh_b     ! b coefficient of RFH(OH) potential
  real(8),save:: para_rfhoh_c     ! c coefficient of RFH(OH) potential
  real(8),save:: para_rfhoh_d     ! d coefficient of RFH(OH) potential
  real(8),save:: para_rfhoh_e     ! e coefficient of RFH(OH) potential
  real(8),save:: para_rfhoh_req   ! r_eq of RFH(OH) potential

  integer,save:: ndoutyp           ! number of DOU type (= 2 is expected)
  character(5),save:: para_doutyp(maxndoutyp)
  real(8),save:: para_dou_wel_o    ! well depth of DOU potential for O-Au
  real(8),save:: para_dou_wel_h    ! well depth of DOU potential for H-Au
  real(8),save:: para_dou_rad_o    ! radius of DOU potential for O-Au
  real(8),save:: para_dou_rad_h    ! radius of DOU potential for H-Au
  real(8),save:: para_dou_beta_o   ! beta of DOU potential for O-Au
  real(8),save:: para_dou_beta_h   ! beta of DOU potential for H-Au
  real(8),save:: para_dou_gamma    ! gamma of DOU potential for H-Au

  integer,save:: nrpvwtyp         ! number of RP-VW interaction type
  character(5),save:: para_rpvwtyp(maxnrpvwtyp)
                                  ! atom pair of Morse type
  real(8),save:: para_rhorpvw(maxnatmtyp,maxnatmtyp)
                                  ! density of RP-VW interaction
  real(8),save:: para_welrpvw(maxnatmtyp,maxnatmtyp)
                                  ! well depth of RP-VW interaction
  real(8),save:: para_radrpvw(maxnatmtyp,maxnatmtyp)
                                  ! radius of RP-VW interaction
  real(8),save:: para_alpharpvw(maxnatmtyp,maxnatmtyp)
                                  ! well steepness of RP-VW interaction

  real(8),save:: alpha_wall       ! constant force for CNP_VW

  real(8),save:: para_posres            ! position restraint spring constant
  real(8),save:: ref_atmcor(3,maxnatom) ! ref. coordinate of position restraint
  integer,save:: nposres(3)             ! number of atoms for position restraint
  integer,save:: index_posres(maxnatom,3) ! atom index for position restraint

  integer,save:: nlfixz                 ! number of z-coordinate fix atoms
  integer,save:: index_nlfixz(maxnatom) ! index of z-coordinate fix atoms
  real(8),save:: cor_lfixz(maxnatom)    ! position of z-coordinate fix atoms

  integer,save:: nlfixzg                 ! number of z-coordinate fix molecules
  integer,save:: index_nlfixzg(maxnatom) ! index of z-coordinate fix molecules
  real(8),save:: cor_lfixzg(maxnatom)    ! position of z-coordinate fix molecules

  character(4),save:: potbias_typ   ! bias potential type
  integer,save:: nint_pbias         ! sampling number w.r.t. reaction coordinate
  real(8),save:: para_potbias(20)   ! parameters for bias potential
  integer,save:: npotbias          ! number of atoms for imposing bias potential
  integer,save:: index_potbias(maxnatom) ! atom index for bias potential

  !---- variables for constraint pair
  integer,save:: nconsttyp                        ! number of all const type
  integer,save:: nconsttyp_p                      ! number of polymer const type
  integer,save:: nconsttyp_w                      ! number of water const type
  character(5),save:: para_consttyp(maxnconsttyp) ! consttyp of const pair

  !---- list of pointer to atoms belonging to molecules
  integer,save:: molept_index(maxnatom+1)   ! pointer index of molept_list
  integer,save:: molept_list(maxnatom+1)    ! pointer list
  integer,save:: nmoleptindex               ! number of molept index
  integer,save:: nmoleptlist                ! number of molept list
  integer,save:: irmolept_list(maxnatom)    ! reverse list of molept(atom->mole)

  !---- list of excluded list
  integer,save:: index_excl(maxnatom+1)        ! index of excluded list
  integer,save:: list_excl(maxnatom*maxnmol)   ! list of excluded pair
  integer,save:: nlistexcl                     ! number of excluded list

  !---- list of molecular based const list
  integer,save:: nmconstindex               ! number of molconst index
  integer,save:: nmconst                    ! number of molconst list
  integer,save:: index_mconst(maxnatom+1)   ! index list of molconst
  integer,save:: list_mconst(maxnconst)     ! list of molconst

  !---- velocity for all atoms
  real(8),save:: atmvel(3,maxnatom)   ! velocity of all atoms

  !---- list of non-bonded list
  integer,save:: nb_indexall(maxnatom+maxnproc)
                                 ! index for location of the list of atom
  integer,save,allocatable:: nb_listall(:)   ! dynamic allocation
!  integer,save:: nb_listall(maxnatom*maxnlist)
                                 ! list of neighboring atoms
  integer,save:: nlistnball      ! number of non-bonded list

  !---- list of Morse bonded list
  integer,save:: mor_indexall(maxnatom+maxnproc)
                                 ! index for location of the list of atom
  integer,save,allocatable:: mor_listall(:)   ! dynamic allocation
!  integer,save:: mor_listall(maxnatom*maxmorlist)
                                 ! list of neighboring atoms
  integer,save:: nlistmorall     ! number of Morse bonded list

  !---- list of SH interaction list
  integer,save:: sho_indexall(maxnatom+maxnproc)
                                 ! index for location of the list of atom
  integer,save,allocatable:: sho_listall(:)
!  integer,save:: sho_listall(maxnatom*maxshlist)
                                 ! list of neighboring atoms
  integer,save:: nlistshoall     ! number of SH list (O atoms)
  integer,save:: shh_indexall(maxnatom+maxnproc)
                                 ! index for location of the list of atom
  integer,save,allocatable:: shh_listall(:)
!  integer,save:: shh_listall(maxnatom*maxshlist)
                                 ! list of neighboring atoms
  integer,save:: nlistshhall     ! number of SH list (H atoms)

  !---- list of RFH interaction list
  integer,save:: rfhfo_indexall(maxnatom+maxnproc)
                                 ! index for location of the list of atom
  integer,save,allocatable:: rfhfo_listall(:)
!  integer,save:: rfhfo_listall(maxnatom*maxrfhlist)
                                 ! list of neighboring atoms
  integer,save:: nlistrfhfoall   ! number of RFH list (Fe-O pair)
  integer,save:: rfhoo_indexall(maxnatom+maxnproc)
                                 ! index for location of the list of atom
  integer,save,allocatable:: rfhoo_listall(:)
!  integer,save:: rfhoo_listall(maxnatom*maxrfhlist)
                                 ! list of neighboring atoms
  integer,save:: nlistrfhooall   ! number of RFH list (O-O pair)
  integer,save:: rfhoh_indexall(maxnatom+maxnproc)
                                 ! index for location of the list of atom
  integer,save,allocatable:: rfhoh_listall(:)
!  integer,save:: rfhoh_listall(maxnatom*maxrfhlist)
                                 ! list of neighboring atoms
  integer,save:: nlistrfhohall   ! number of RFH list (O-H pair)

  !---- list of DOU interaction list
  integer,save:: douo_indexall(maxnatom+maxnproc)
                                 ! index for location of the list of atom
  integer,save,allocatable:: douo_listall(:)
!  integer,save:: douo_listall(maxnatom*maxshlist)
                                 ! list of neighboring atoms
  integer,save:: nlistdouoall    ! number of DOU list (O atoms)
  integer,save:: douh_indexall(maxnatom+maxnproc)
                                 ! index for location of the list of atom
  integer,save,allocatable:: douh_listall(:)
!  integer,save:: douh_listall(maxnatom*maxshlist)
                                 ! list of neighboring atoms
  integer,save:: nlistdouhall    ! number of DOU list (H atoms)

  !---- list of RP-VW interaction list
  integer,save:: rpvw_indexall(maxnatom+maxnproc)
                                 ! index for location of the list of atom
  integer,save,allocatable:: rpvw_listall(:)
!  integer,save:: rpvw_listall(maxnatom*maxrpvwlist)
                                 ! list of neighboring atoms
  integer,save:: nlistrpvwall    ! number of RP-VW interaction list

  !---- variables for ewald sum
  integer,save:: nwave              ! number of lattice vector
  real(8),save:: nvec(3,maxwave)    ! k vector (k=(kx,ky,kz))
  real(8),save:: kwave(3,maxwave)   ! reciprocal lattice vector
                                    ! (K=(kx/lx,ky/ly,kz/lz))
  real(8),save:: awave(maxwave)     ! = exp(-(pi*K/a)^2)/K^2

  !---- variables for ewald sum
  logical,save:: fennell_flag    ! if use fennell method
  real(8),save:: fennell_shift   ! shitf term of fennell potential
  real(8),save:: fennell_damp    ! damp term of fennell potential

  !---- variables for Nose-Hoover chain
  real(8),save:: xlogs(maxnhc)   ! xi of the thermostat coordinate
  real(8),save:: vlogs(maxnhc)   ! vxi of the thermostat velocity
  real(8),save:: qmass(maxnhc)   ! mass of the NHC thermostat

  !---- variables for Andersen (Hoover type) barostat
  real(8),save:: xlogv           ! epsilon of the barostat coordinate
  real(8),save:: vlogv           ! vepsilon of the barostat velocity
  real(8),save:: vmass           ! mass of the Andersen (Hoover type) barostat

  real(8),save:: xboxh(3)           ! cell vector h of the barostat coordinate
  real(8),save:: vboxg(3)           ! vg vector of the barostat coordinate

  !---- coefficient of Yoshida-Suzuki method
  real(8),save:: yscoeff(maxnyosh)   ! coefficient of Yoshida-Suzuki method

  !---- spline interpolation for interaction sum
  real(8),save:: x_tbl(max_interpol)            ! x position
  real(8),save:: spltbl_real(max_interpol)      ! Ewald real
  real(8),save:: spltbl_realpot(max_interpol)   ! Ewald pot real
  real(8),save:: spltbl_excs(max_interpol)      ! Ewald excess
  real(8),save:: spltbl_excspot(max_interpol)   ! Ewald pot excess
  real(8),save:: spl_b_real(max_interpol)       ! Ewald real
  real(8),save:: spl_c_real(max_interpol)       ! Ewald real
  real(8),save:: spl_d_real(max_interpol)       ! Ewald real
  real(8),save:: spl_b_realpot(max_interpol)    ! Ewald pot real
  real(8),save:: spl_c_realpot(max_interpol)    ! Ewald pot real
  real(8),save:: spl_d_realpot(max_interpol)    ! Ewald pot real
  real(8),save:: spl_b_excs(max_interpol)       ! Ewald excess
  real(8),save:: spl_c_excs(max_interpol)       ! Ewald excess
  real(8),save:: spl_d_excs(max_interpol)       ! Ewald excess
  real(8),save:: spl_b_excspot(max_interpol)    ! Ewald pot excess
  real(8),save:: spl_c_excspot(max_interpol)    ! Ewald pot excess
  real(8),save:: spl_d_excspot(max_interpol)    ! Ewald pot excess

  integer,save:: nspltbl          ! number of spline interpolation points
  real(8),save:: spltbl_int       ! interval of spline points
  integer,save:: spl_order        ! order of spline function

  !---- variables for streaming velocity
  real(8),save:: strmvel(3,maxnvelregion)  ! streaming velocity of each slab
  real(8),save:: atmvel_strm(3,maxnatom)   ! streaming velocity of each atom
  integer,save:: nvelregion                ! number of slabs for streaming vel.
  real(8),save:: strmzpos1(maxnvelregion)  ! slab def. for streaming velocity
  real(8),save:: strmzpos2(maxnvelregion)  ! slab def. for streaming velocity

  !---- variables for region-based Langevin thermostat
  integer,save:: natom_lt(maxntcregion) ! number of atoms in region
  ! integer,save:: degfree_lt(maxntcregion) ! degree of freedom in region
  integer,save:: atmindex_lt(maxnatom,maxntcregion) ! index of atoms in region
  real(8),save:: frand_lt(3,maxnatom)      ! random force

 !---- variables for solid-liquid interface  -by guo

  real(8),save:: pot_interface_1
  real(8),save:: pot_interface_2
  real(8),save:: pot_interface_no

  real(8),save:: kin_interface_1
  real(8),save:: kin_interface_2
  real(8),save:: kin_interface_no

  real(8),save:: wall_pot_interface_1_ar1
  real(8),save:: wall_pot_interface_1_ar2
  real(8),save:: wall_pot_interface_2_ar1
  real(8),save:: wall_pot_interface_2_ar2
  real(8),save:: wall_pot_interface_no

  real(8),save:: pot_vdw_atm_ar12(maxnatom)
  real(8),save:: pot_vdw_atm_ar11(maxnatom)
  real(8),save:: pot_vdw_atm_ar22(maxnatom)

  real(8),save:: hfpotlj_atm_ar11(3,maxnhfregion)
  real(8),save:: hfpotlj_atm_ar12(3,maxnhfregion)
  real(8),save:: hfpotlj_atm_ar22(3,maxnhfregion)

  real(8),save:: viriljt_atm_ar11(3,3,maxnatom,maxnhfregion)
  real(8),save:: viriljt_atm_ar12(3,3,maxnatom,maxnhfregion)
  real(8),save:: viriljt_atm_ar22(3,3,maxnatom,maxnhfregion)
  real(8),save:: viriljt_atm_ar1pt(3,3,maxnatom,maxnhfregion)
  real(8),save:: viriljt_atm_ar2pt(3,3,maxnatom,maxnhfregion)

  real(8),save::  hfvirilj_atm_ar11(3,maxnhfregion)
  real(8),save::  hfvirilj_atm_ar12(3,maxnhfregion)
  real(8),save::  hfvirilj_atm_ar22(3,maxnhfregion)
  real(8),save::  hfvirilj_atm_ar1pt(3,maxnhfregion)
  real(8),save::  hfvirilj_atm_ar2pt(3,maxnhfregion)

end module md_global
