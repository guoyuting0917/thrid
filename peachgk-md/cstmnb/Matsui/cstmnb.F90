!**********************************
!*  cstmnb.f90 Ver.1.3            *
!*      for peachgk_md.f          *
!*            by G.Kikugawa       *
!**********************************
! Time-stamp: <2015-01-22 01:10:54 gota>

module cstmnb

  use md_global
  use mpi_global

  implicit none

  private

!!! define global variables here
  integer,public,save:: ncstmnb	              ! number of all cstm NB atoms
  integer,public,save:: ncstmnblist(maxnatom) ! list of cstm NB atoms
  logical,public,save:: ifcstmnb(maxnatom)    ! if cstm NB atom or not

  integer,public,save:: cstmnb_indexall(maxnatom+maxnproc)
                                    ! index for location of the list of atom
  integer,public,save,allocatable:: cstmnb_listall(:)
!  integer,public,save:: cstmnb_listall(maxnatom*maxcstmnblist)
                                    ! list of neighboring atoms
  integer,public,save:: nlistcstmnball     ! number of custom NB list

  integer,public,save:: ncstmnbtyp  ! number of cstmnb type
  character(2),public,save:: para_cstmnbtyp(maxncstmnbtyp)
                                    ! atom of custom NB type

!  logical,public,save:: ifcalMATSUI(maxnatmtyp,maxnatmtyp)
!                                    ! if this pair is calculated in cstmnb

  real(8),public,save:: rcut_MATSUI  ! cutoff length for Matsui potential
  real(8),public,save:: rcut_bookMATSUI  
                          ! cutoff length for book-keeping of Matsui potential
  integer,public,save:: nstep_bookMATSUI  ! interval for bookkeeping

  real(8),public,save:: AA_MATSUI(maxnatmtyp,maxnatmtyp)
                                     ! Ai + Aj parameter for Matsui potential
  real(8),public,save:: BB_MATSUI(maxnatmtyp,maxnatmtyp)
                                     ! Bi + Bj parameter for Matsui potential
  real(8),public,save:: CC_MATSUI(maxnatmtyp,maxnatmtyp)
                                     ! Ci * Cj parameter for Matsui potential
  real(8),public,save:: D_MATSUI
                                     ! D parameter for Matsui potential

!!! Do not forget to list globally opened procedures here
  public:: calcstmnb, calcstmnbp
#if defined(HF)
  public:: calcstmnbp_hf
#endif

contains

  subroutine calcstmnb(xcel,ycel,zcel, &
       &               npoly,npolytyp,npoly_mole,npoly_atom, &
       &               nwater,nmatom,nmatyp,nmatomtyp, &
       &               force,pot_cstmnb)

    implicit none

!
!
!     Calculate custom NB interactions
!
!
!!! You can implement your own interaction functions by modifying this file,
!!! rdpara_cstmnb.F90, mklist_cstmnb.F90, and corresponding
!!! para_cstmnb.dat files.
!!
!!! This procedure is designed for calculating Matsui potential
!!! which is modeled for metal oxide substrate.
!!! (Matsui, M. Miner. Mag., Vol. 58A, pp. 571-572, (1994).)
!
! ARGUMENTS:
!   INPUT
    real(8),intent(in):: xcel            ! x cell length[non-d]
    real(8),intent(in):: ycel            ! y cell length[non-d]
    real(8),intent(in):: zcel            ! z cell length[non-d]

    integer,intent(in):: npoly           ! all number of poly
    integer,intent(in):: npolytyp        ! number of poly type
    integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
    integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

    integer,intent(in):: nwater          ! number of H2O molecules

    integer,intent(in):: nmatom          ! number of monatomic molecules
    integer,intent(in):: nmatyp          ! number of species of monatomic mole.
    integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

!   INPUT & OUTPUT
    real(8),intent(inout):: force(:,:)   ! force calculated here
                                       ! total force
                                       !  (including vdw, bond, angle etc)

!   OUTPUT
    real(8),intent(out):: pot_cstmnb       ! Custom NB potential

! LOCAL:
    real(8):: fcstmnb(3,natom)              ! custom NB force calculated here

    real(8):: pi                  ! = 3.14
    real(8):: potij               ! temporal potential due to IJ
    real(8):: ftmp                ! temporal force
    real(8):: fij(3)              ! force to be subtracted due to
                                  ! excluded atoms

    real(8):: rij(3)              ! rj-ri
    real(8):: rij_abs             ! |rij|
    real(8):: rij_abs2            ! |rij|**2
    real(8):: rij_abs_inv1        ! 1/|rij|
    real(8):: rij_abs_inv2        ! 1/|rij|**2
    real(8):: rtmp6
    real(8):: rexp

    real(8):: box(3)              ! BOX size
    real(8):: box_inv(3)          ! inverse of BOX size

    real(8):: sqcut               ! rcut*rcut

    integer:: itype, jtype        ! VDW atom types
    integer:: i,j                 ! do loop index
    integer:: j1, j2              ! limits of do loop
    integer:: jj                  ! atom index of the interacting atom

    integer:: cstmnb1

! FUNCTIONS:

!     +     +     +     +     +     +     +


!---- some initialization

    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    box(1) = xcel
    box(2) = ycel
    box(3) = zcel
    box_inv(1:3) = 1.0d0/box(1:3)

    sqcut = rcut_MATSUI * rcut_MATSUI

!   --- CALCULATE NONBONDED FORCE ---

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, ncstmnb - 1
    looplast = ncstmnb - 1

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = cstmnb_indexall(cstmnb1)              ! the first and
       j2 = cstmnb_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting
                                                  ! with i

       i = ncstmnblist(cstmnb1)
       itype = atmindex(i)       ! atom type of atom(i)

       IF (cstmnb_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = cstmnb_listall(j)          ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!         --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3) * anint(rij(1:3)*box_inv(1:3))

!         * box_inv rather than /box accellerates the computation

          rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2
                                ! |rij|^2
!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2
          rij_abs       = sqrt(rij_abs2) ! |rij|
          rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

!         --- LJ term ---
          rtmp6 = rij_abs_inv2 * rij_abs_inv2 * rij_abs_inv2

          potij = - CC_MATSUI(itype,jtype) * rtmp6 ! -CiCj/|rij|^6

          ftmp  = 6.0d0 * potij * rij_abs_inv2 ! -6CiCj/|rij|^8

!         --- exponential term ---
          rexp = D_MATSUI * exp((AA_MATSUI(itype,jtype) - rij_abs)   &
               &               /BB_MATSUI(itype,jtype))
                               ! D*exp((Ai+Aj - |rij|)/(Bi+Bj))
          potij = potij + rexp * BB_MATSUI(itype,jtype)
                               ! += D*(Bi+Bj)*exp((Ai+Aj - |rij|)/(Bi+Bj))
          ftmp = ftmp + rexp * rij_abs_inv1
                               ! += D*exp((Ai+Aj - |rij|)/(Bi+Bj)) / |rij|

!         --- add both terms ---
          pot_cstmnb = pot_cstmnb + potij

          fij(1:3) = ftmp * rij(1:3)
          fcstmnb(1:3,i)  = fcstmnb(1:3,i)  + fij(1:3)
          fcstmnb(1:3,jj) = fcstmnb(1:3,jj) - fij(1:3)

       END DO
    END DO

! --- ADD FCSTMNB TO FORCE ---
    force(1:3,1:natom) = force(1:3,1:natom) + fcstmnb(1:3,1:natom)

!     +     +     +     +     +     +     +

  end subroutine calcstmnb

!-----------------------------------------------------------------

  subroutine calcstmnbp(xcel,ycel,zcel, &
       &                npoly,npolytyp,npoly_mole,npoly_atom, &
       &                nwater,nmatom,nmatyp,nmatomtyp, &
       &                force,pot_cstmnb, &
       &                for_viri_cstmnb,pot_viri_cstmnb, &
       &                pot_virit_cstmnb)

    implicit none

!
!
!     Calculate custom NB interactions
!
!
! ARGUMENTS:
!   INPUT
    real(8),intent(in):: xcel            ! x cell length[non-d]
    real(8),intent(in):: ycel            ! y cell length[non-d]
    real(8),intent(in):: zcel            ! z cell length[non-d]

    integer,intent(in):: npoly           ! all number of poly
    integer,intent(in):: npolytyp        ! number of poly type
    integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
    integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
  
    integer,intent(in):: nwater          ! number of H2O molecules

    integer,intent(in):: nmatom          ! number of monatomic molecules
    integer,intent(in):: nmatyp          ! number of species of monatomic mole.
    integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

!   INPUT & OUTPUT

    real(8),intent(inout):: force(:,:)   ! force calculated here
                                         ! total force
                                         !  (including vdw, bond, angle etc)

    real(8),intent(inout):: for_viri_cstmnb(:,:)
                                         ! virial(custom NB force) of each atom
    real(8),intent(inout):: pot_viri_cstmnb
                                   ! virial(custom NB potential) of each atom
    real(8),intent(inout):: pot_virit_cstmnb(:,:) ! virial tensor (custom NB)

!   OUTPUT
    real(8),intent(out):: pot_cstmnb       ! Custom NB potential

! LOCAL:
    real(8):: fcstmnb(3,natom)              ! custom NB force calculated here

    real(8):: pi                  ! = 3.14
    real(8):: potij               ! temporal potential due to IJ
    real(8):: ftmp                ! temporal force
    real(8):: fij(3)              ! force to be subtracted due to
                                  ! excluded atoms

    real(8):: rij(3)              ! rj-ri
    real(8):: rij_abs             ! |rij|
    real(8):: rij_abs2            ! |rij|**2
    real(8):: rij_abs_inv1        ! 1/|rij|
    real(8):: rij_abs_inv2        ! 1/|rij|**2
    real(8):: rtmp6
    real(8):: rexp

    real(8):: box(3)              ! BOX size
    real(8):: box_inv(3)          ! inverse of BOX size

    real(8):: sqcut               ! rcut*rcut

    integer:: itype, jtype        ! VDW atom types
    integer:: i,j,n               ! do loop index
    integer:: j1, j2              ! limits of do loop
    integer:: jj                  ! atom index of the interacting atom

    integer:: cstmnb1

! FUNCTIONS:

!     +     +     +     +     +     +     +

!---- some initialization

    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    box(1) = xcel
    box(2) = ycel
    box(3) = zcel
    box_inv(1:3) = 1.0d0/box(1:3)

    sqcut = rcut_MATSUI * rcut_MATSUI

!   --- CALCULATE NONBONDED FORCE ---

    pot_virit_cstmnb(1:3,1:3) = 0.0d0

    pot_viri_cstmnb = 0.0d0

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, ncstmnb - 1
    looplast = ncstmnb - 1

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = cstmnb_indexall(cstmnb1)              ! the first and
       j2 = cstmnb_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting
                                                  ! with i
       i = ncstmnblist(cstmnb1)
       itype = atmindex(i)       ! atom type of atom(i)

       IF (cstmnb_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = cstmnb_listall(j)          ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!         --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3) * anint(rij(1:3)*box_inv(1:3))

!         * box_inv rather than /box accellerates the computation

          rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2
                                ! |rij|^2
!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2
          rij_abs       = sqrt(rij_abs2) ! |rij|
          rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

!         --- LJ term ---
          rtmp6 = rij_abs_inv2 * rij_abs_inv2 * rij_abs_inv2

          potij = - CC_MATSUI(itype,jtype) * rtmp6 ! -CiCj/|rij|^6

          ftmp  = 6.0d0 * potij * rij_abs_inv2 ! -6CiCj/|rij|^8

!         --- exponential term ---
          rexp = D_MATSUI * exp((AA_MATSUI(itype,jtype) - rij_abs)   &
               &               /BB_MATSUI(itype,jtype))
                               ! D*exp((Ai+Aj - |rij|)/(Bi+Bj))
          potij = potij + rexp * BB_MATSUI(itype,jtype)
                               ! += D*(Bi+Bj)*exp((Ai+Aj - |rij|)/(Bi+Bj))
          ftmp = ftmp + rexp * rij_abs_inv1
                               ! += D*exp((Ai+Aj - |rij|)/(Bi+Bj)) / |rij|

!         --- add both terms ---
          pot_cstmnb = pot_cstmnb + potij

          fij(1:3) = ftmp * rij(1:3)
          fcstmnb(1:3,i)  = fcstmnb(1:3,i)  + fij(1:3)
          fcstmnb(1:3,jj) = fcstmnb(1:3,jj) - fij(1:3)

!         --- calculate virial tensor ---
          do n=1,3
             pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                  &                  + fij(1:3)*rij(n)
          end do
          pot_viri_cstmnb = pot_viri_cstmnb &
               &          + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

          for_viri_cstmnb(1:3,i)  = for_viri_cstmnb(1:3,i)  + fij(1:3)
          for_viri_cstmnb(1:3,jj) = for_viri_cstmnb(1:3,jj) - fij(1:3)

       END DO
    END DO

!   --- ADD FCSTMNB TO FORCE ---
    force(1:3,1:natom) = force(1:3,1:natom) + fcstmnb(1:3,1:natom)

!     +     +     +     +     +     +     +

  end subroutine calcstmnbp

!-----------------------------------------------------------------

#if defined(HF)
  subroutine calcstmnbp_hf(xcel,ycel,zcel,   &
       &                   npoly,npolytyp,npoly_mole,npoly_atom, &
       &                   nwater,nmatom,nmatyp,nmatomtyp, &
       &                   force,pot_cstmnb, &
       &                   for_viri_cstmnb,pot_viri_cstmnb, &
       &                   pot_virit_cstmnb, &
       &                   pot_cstmnb_atm,viricstmnbt_atm,   &
       &                   ifhfvol,   &
       &                   nhfregion,hfzpos1,hfzpos2,   &
       &                   hftyp_atm,   &
       &                   molecom)

    use interface_interact

    implicit none

!
!
!     Calculate custom NB interactions
!
!
! ARGUMENTS:
!   INPUT
    real(8),intent(in):: xcel            ! x cell length[non-d]
    real(8),intent(in):: ycel            ! y cell length[non-d]
    real(8),intent(in):: zcel            ! z cell length[non-d]

    integer,intent(in):: npoly           ! all number of poly
    integer,intent(in):: npolytyp        ! number of poly type
    integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
    integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
  
    integer,intent(in):: nwater          ! number of H2O molecules

    integer,intent(in):: nmatom          ! number of monatomic molecules
    integer,intent(in):: nmatyp          ! number of species of monatomic mole.
    integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

    logical,intent(in):: ifhfvol    ! local volume-based or local surface-based
    integer,intent(in):: nhfregion  ! number of region to calculate heat flux
    real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                        ! z-position of region for heat flux

    integer,intent(in):: hftyp_atm(:)   ! atom- or mole-based heat flux cal. 
                                        !   for each atom

    real(8),intent(in):: molecom(:,:)   ! center of mass of molecule

!   INPUT & OUTPUT

    real(8),intent(inout):: force(:,:)   ! force calculated here
                                         ! total force
                                         !  (including vdw, bond, angle etc)

    real(8),intent(inout):: for_viri_cstmnb(:,:)
                                         ! virial(custom NB force) of each atom
    real(8),intent(inout):: pot_viri_cstmnb
                                   ! virial(custom NB potential) of each atom
    real(8),intent(inout):: pot_virit_cstmnb(:,:) ! virial tensor (custom NB)

    real(8),intent(inout):: pot_cstmnb_atm(:) 
                                       ! custom NB potential of each atom
    real(8),intent(inout):: viricstmnbt_atm(:,:,:,:)
                                       ! virial tensor of each atom (custom NB)

!   OUTPUT
    real(8),intent(out):: pot_cstmnb       ! Custom NB potential

! LOCAL:
    real(8):: fcstmnb(3,natom)              ! custom NB force calculated here

    real(8):: pi                  ! = 3.14
    real(8):: potij               ! temporal potential due to IJ
    real(8):: ftmp                ! temporal force  
    real(8):: fij(3)              ! force to be subtracted due to 
                                  ! excluded atoms
    real(8):: fji(3)              ! force to be subtracted due to 
                                  ! excluded atoms = -fij

    real(8):: rij(3)              ! rj-ri
    real(8):: rij_abs             ! |Rij|   
    real(8):: rij_abs2            ! |Rij|**2
    real(8):: rij_abs_inv1        ! |Rij|**-1
    real(8):: rij_abs_inv2        ! 1/|rij|**2
    real(8):: rtmp6
    real(8):: rexp

    real(8):: box(3)              ! BOX size
    real(8):: box_inv(3)          ! inverse of BOX size

    real(8):: sqcut               ! rcut*rcut

    integer:: itype, jtype        ! VDW atom types
    integer:: i,j,n               ! do loop index
    integer:: j1, j2              ! limits of do loop
    integer:: jj                  ! atom index of the interacting atom

    integer:: cstmnb1

    real(8):: nbodycoeff          ! n-body coefficient for virial term of heatf

! FUNCTIONS:

!     +     +     +     +     +     +     +

!---- some initialization

    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    box(1) = xcel
    box(2) = ycel
    box(3) = zcel
    box_inv(1:3) = 1.0d0/box(1:3)

    sqcut = rcut_MATSUI * rcut_MATSUI

!     --- CALCULATE NONBONDED FORCE ---

    pot_virit_cstmnb(1:3,1:3) = 0.0d0

    pot_viri_cstmnb = 0.0d0

    nbodycoeff = 0.5d0        ! = 1/2

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, ncstmnb - 1
    looplast = ncstmnb - 1

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = cstmnb_indexall(cstmnb1)              ! the first and
       j2 = cstmnb_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting
                                                  ! with i
       i = ncstmnblist(cstmnb1)
       itype = atmindex(i)       ! atom type of atom(i)

       IF (cstmnb_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = cstmnb_listall(j)          ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!         --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3) * anint(rij(1:3)*box_inv(1:3))

!         * box_inv rather than /box accellerates the computation

          rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2
                                ! |rij|^2
!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2
          rij_abs       = sqrt(rij_abs2) ! |rij|
          rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

!         --- LJ term ---
          rtmp6 = rij_abs_inv2 * rij_abs_inv2 * rij_abs_inv2

          potij = - CC_MATSUI(itype,jtype) * rtmp6 ! -CiCj/|rij|^6

          ftmp  = 6.0d0 * potij * rij_abs_inv2 ! -6CiCj/|rij|^8

!         --- exponential term ---
          rexp = D_MATSUI * exp((AA_MATSUI(itype,jtype) - rij_abs)   &
               &               /BB_MATSUI(itype,jtype))
                               ! D*exp((Ai+Aj - |rij|)/(Bi+Bj))
          potij = potij + rexp * BB_MATSUI(itype,jtype)
                               ! += D*(Bi+Bj)*exp((Ai+Aj - |rij|)/(Bi+Bj))
          ftmp = ftmp + rexp * rij_abs_inv1
                               ! += D*exp((Ai+Aj - |rij|)/(Bi+Bj)) / |rij|

!         --- add both terms ---
          pot_cstmnb = pot_cstmnb + potij

          fij(1:3) = ftmp * rij(1:3)
          fcstmnb(1:3,i)  = fcstmnb(1:3,i)  + fij(1:3)
          fcstmnb(1:3,jj) = fcstmnb(1:3,jj) - fij(1:3)

!         --- calculate virial tensor ---

!         - for heat flux calculation
          fji(1:3) = -fij(1:3)

          pot_cstmnb_atm(i) = pot_cstmnb_atm(i) + 0.5d0 * potij
          pot_cstmnb_atm(jj) = pot_cstmnb_atm(jj) + 0.5d0 * potij

          call calvirihf(i,jj,fij,fji, &
               &         nbodycoeff, &
               &         box,box_inv, &
               &         ifhfvol, &
               &         nhfregion,hfzpos1,hfzpos2, &
               &         hftyp_atm, &
               &         molecom, &
               &         viricstmnbt_atm)

!         - for pressure calculation
          do n=1,3
             pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                  &                  + fij(1:3)*rij(n)
          end do

          pot_viri_cstmnb = pot_viri_cstmnb &
               &          + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

          for_viri_cstmnb(1:3,i)  = for_viri_cstmnb(1:3,i)  + fij(1:3)
          for_viri_cstmnb(1:3,jj) = for_viri_cstmnb(1:3,jj) - fij(1:3)

       END DO
    END DO

!   --- ADD FCSTMNB TO FORCE ---
    force(1:3,1:natom) = force(1:3,1:natom) + fcstmnb(1:3,1:natom)

!     +     +     +     +     +     +     +

  end subroutine calcstmnbp_hf
#endif

end module cstmnb
