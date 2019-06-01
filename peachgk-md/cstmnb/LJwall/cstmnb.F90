!**********************************
!*  cstmnb.f90 Ver.1.2            *
!*      for peachgk_md.f          *
!*            by G.Kikugawa       *
!**********************************
! Time-stamp: <Aug 21 2017>

!!! Virtual wall (continuous media) interaction by LJ potential
!!!   originally coded by Jo Suzuki, IFS, Tohoku University

!!! In this version, you must use calcstmnb() routine, not calcstmnbp(),
!!! because the pressure calculation is not available for this interaction.

module cstmnb
  use md_global
  use mpi_global

  implicit none

  private

!!! define global variables here
  integer,public,save:: ncstmnb                 !   number of all cstm NB atoms
  integer,public,save:: ncstmnblist(maxnatom)   !   list of cstm NB atoms
  logical,public,save:: ifcstmnb(maxnatom)      ! if cstm NB atom or not
  integer,public,save:: cstmnb_listall(maxnatom,2)
                                    ! list for atom interacting with wall-1,2
  integer,public,save:: nlistcstmnball(2)   ! number of custom NB list1,2
  integer,public,save:: nlistcstmnball_lcl(2) ! local variable nlistcstmnball fpr MPI
 
  integer,public,save:: ncstmnbtyp  ! number of cstmnb type
  character(2),public,save:: para_cstmnbtyp(maxncstmnbtyp)  ! atom of custom NB type
  integer,public,save:: index_cstmnbtyp(maxnatom) ! index of  custom NB type
  real(8),public,save:: zcut_WALL       ! cutoff length for wall potential
  real(8),public,save:: zcut_bookWALL   ! cutoff length for book-keeping of wall potential
  integer,public,save:: nstep_bookWALL  ! interval for bookkeeping

  real(8),public,save:: E_WALL(maxncstmnbtyp)
                                     ! epsilon[J] for wall potential
  real(8),public,save:: S_WALL(maxncstmnbtyp)
                                     ! sigma[m] for wall potential
  real(8),public,save:: R_WALL
                                     ! rho[m^-3] for wall potential
  real(8),public,save:: dstnc
                                     ! distance[m]  between pair parallel walls
                                     ! which is the length of cell in z direction with vacuum regions removed

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
!!! This procedure is designed for calculating wall potential
!!! which is modeled for a flat continuum wall.
!!
!!! To privent atoms from being set extremely near walls,
!!! make z-cell size be smaller than walls' dstnc(distance) at first.
!!! After relaxation, enlarge cell to dstnc.
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
    real(8),intent(out):: pot_cstmnb     ! Custom NB potential

!   LOCAL:
    real(8):: pi                   ! = 3.14
    real(8):: poti                 ! potential of atom i
    real(8):: ftmp                 ! temporal force
    real(8):: fz(natom)            ! force acting on atom i by walls

    real(8):: ziw                  ! zi-zw
    real(8):: z_inv                ! 1 / |ziw|
    real(8):: z_inv_sigma          ! sigma / |ziw|
    real(8):: z_inv3_sigma         ! (sigma / |ziw|)**3
    real(8):: z_inv9_sigma         ! (sigma / |ziw|)**9
    real(8):: S_WALL_3             ! S_WALL**3
    real(8):: pfact                ! 4*pi*eps*sigma^3*rho

    integer:: i, j                 ! do loop index

    integer:: cstmnbk

! FUNCTIONS:

!     +     +     +     +     +     +     +

!---- some initialization
    fz(1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    pi = acos(-1.0d0)

!   --- CALCULATE NONBONDED FORCE ---

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!

!   --- wall 1 ---

    DO cstmnbk = 1, nlistcstmnball_lcl(1)
       i = cstmnb_listall(cstmnbk,1)   ! number of atom interacting with wall-k

       ziw = atmcor(3,i)        ! distance between wall-k and particle. ziw > 0

!      --- cut off ---
       if (ziw > zcut_WALL) cycle

       j = index_cstmnbtyp(i)

       z_inv = 1.0d0 / ziw
       z_inv_sigma = S_WALL(j) * z_inv   ! = sigma / |ziw|
       z_inv3_sigma = z_inv_sigma * z_inv_sigma * z_inv_sigma      
                                         ! = sigma^3 / |ziw|^3
       z_inv9_sigma = z_inv3_sigma * z_inv3_sigma * z_inv3_sigma
                                         ! = sigma^9 / |ziw|^9
       S_WALL_3 = S_WALL(j) * S_WALL(j) * S_WALL(j)   ! = sigma^3
       pfact = 4.0d0 * pi * E_WALL(j)* S_WALL_3 * R_WALL
                                         ! = 4*pi*eps*sigma^3*rho

       poti = pfact * (z_inv9_sigma / 45.0d0 - z_inv3_sigma / 6.0d0)
       pot_cstmnb = pot_cstmnb + poti
       ftmp = pfact * (z_inv9_sigma * 0.2d0 - z_inv3_sigma * 0.5d0) * z_inv
       fz(i) = fz(i) + ftmp

    END DO

!   --- wall 2 ---

    DO cstmnbk = 1, nlistcstmnball_lcl(2)
       i = cstmnb_listall(cstmnbk,2)   ! number of atom interacting with wall-k

       ziw = dstnc - atmcor(3,i) ! distance between wall-k and particle. ziw > 0

!      --- cut off ---
       if (ziw > zcut_WALL) cycle

       j = index_cstmnbtyp(i)

       z_inv = 1.0d0 / ziw
       z_inv_sigma = S_WALL(j) * z_inv   ! = sigma / |ziw|
       z_inv3_sigma = z_inv_sigma * z_inv_sigma * z_inv_sigma
                                         ! = sigma^3 / |ziw|^3
       z_inv9_sigma = z_inv3_sigma * z_inv3_sigma * z_inv3_sigma
                                         ! = sigma^9 / |ziw|^9
       S_WALL_3 = S_WALL(j) * S_WALL(j) * S_WALL(j)  ! = sigma^3
       pfact = 4.0d0 * pi * E_WALL(j)* S_WALL_3 * R_WALL
                                         ! = 4*pi*eps*sigma^3*rho

       poti = pfact * (z_inv9_sigma / 45.0d0 - z_inv3_sigma / 6.0d0)
       pot_cstmnb = pot_cstmnb + poti
       ftmp = - pfact * (z_inv9_sigma * 0.2d0 - z_inv3_sigma * 0.5d0) * z_inv
       fz(i) = fz(i) + ftmp

    END DO

!   --- ADD FCSTMNB TO FORCE ---
    force(3,1:natom) = force(3,1:natom) + fz(1:natom)

!     +     +     +     +     +     +     +

  end subroutine calcstmnb

!-----------------------------------------------------------------

!!!
!!!subroutine calcstmnbp below (pressure calculattion) is NOT ready to use.
!!!

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
    real(8):: pi                   ! = 3.14
    real(8):: poti                 ! potential of atom i
    real(8):: ftmp                 ! temporal force
    real(8):: fz(natom)            ! force acting on atom i by walls

    real(8):: ziw                  ! zi-zw
    real(8):: z_inv                ! 1 / |ziw|
    real(8):: z_inv_sigma          ! sigma / |ziw|
    real(8):: z_inv3_sigma         ! (sigma / |ziw|)**3
    real(8):: z_inv9_sigma         ! (sigma / |ziw|)**9
    real(8):: S_WALL_3             ! S_WALL**3
    real(8):: pfact                ! 4*pi*eps*sigma^3*rho

    integer:: i, j                 ! do loop index

    integer:: cstmnbk


! FUNCTIONS:

!     +     +     +     +     +     +     +

!---- some initialization
    fz(1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    pi = 0.0d0   !!! In this version, you must use calcstmnb() routine.
                 !!! DO NOT use pressure calcultion.


!   --- CALCULATE NONBONDED FORCE ---

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!       DO i = 1, nlistcstmnball

!   --- wall 1 ---

    DO cstmnbk = 1, nlistcstmnball_lcl(1)
       i = cstmnb_listall(cstmnbk,1)   ! number of atom interacting with wall-k

       ziw = atmcor(3,i)        ! distance between wall-k and particle. ziw > 0


!      --- cut off ---
       if (ziw > zcut_WALL) cycle

       j = index_cstmnbtyp(i)

       z_inv = 1.0d0 / ziw
       z_inv_sigma = S_WALL(j)* z_inv   ! = sigma / |ziw|
       z_inv3_sigma = z_inv_sigma * z_inv_sigma * z_inv_sigma
                                        ! = sigma^3 / |ziw|^3
       z_inv9_sigma = z_inv3_sigma * z_inv3_sigma * z_inv3_sigma
                                        ! = sigma^9 / |ziw|^9
       S_WALL_3 = S_WALL(j) * S_WALL(j) * S_WALL(j)   ! = sigma^3

       pfact = 4.0d0 * pi * E_WALL(j)* S_WALL_3 * R_WALL
                                        ! = 4*pi*eps*sigma^3*rho

       poti = pfact * (z_inv9_sigma / 45.0d0 - z_inv3_sigma / 6.0d0)
       pot_cstmnb = pot_cstmnb + poti
       ftmp = pfact * (z_inv9_sigma * 0.2d0 - z_inv3_sigma * 0.5d0) * z_inv
       fz(i) = fz(i) + ftmp

!      --- calculate virial tensor and add pot term ---
       pot_virit_cstmnb(3,3) = pot_virit_cstmnb(3,3) &
            &                + fz(i) * atmcor(3,i)

       pot_viri_cstmnb = pot_viri_cstmnb + pot_virit_cstmnb(3,3)
       for_viri_cstmnb(3,1:natom)  = for_viri_cstmnb(3,1:natom) + fz(1:natom) 

    END DO

!   --- wall 2 ---

    DO cstmnbk = 1, nlistcstmnball_lcl(2)
       i = cstmnb_listall(cstmnbk,2)   ! number of atom interacting with wall-k

       ziw = dstnc - atmcor(3,i) ! distance between wall-k and particle. ziw > 0

!      --- cut off ---
       if (ziw > zcut_WALL) cycle

       j = index_cstmnbtyp(i)

       z_inv = 1.0d0 / ziw
       z_inv_sigma = S_WALL(j) * z_inv   ! = sigma / |ziw|
       z_inv3_sigma = z_inv_sigma * z_inv_sigma * z_inv_sigma
                                         ! = sigma^3 / |ziw|^3
       z_inv9_sigma = z_inv3_sigma * z_inv3_sigma * z_inv3_sigma
                                         ! = sigma^9 / |ziw|^9
       S_WALL_3 = S_WALL(j) * S_WALL(j) * S_WALL(j)   ! = sigma^3
       pfact = 4.0d0 * pi * E_WALL(j)* S_WALL_3 * R_WALL
                                      ! = 4*pi*eps*sigma^3*rho

       poti = pfact * (z_inv9_sigma / 45.0d0 - z_inv3_sigma / 6.0d0)
       pot_cstmnb = pot_cstmnb + poti
       ftmp = - pfact * (z_inv9_sigma * 0.2d0 - z_inv3_sigma * 0.5d0) * z_inv
       fz(i) = fz(i) + ftmp

!      --- calculate virial tensor and add pot term ---
       pot_virit_cstmnb(3,3) = pot_virit_cstmnb(3,3) &
            &                + fz(i) * atmcor(3,i)

       pot_viri_cstmnb = pot_viri_cstmnb + pot_virit_cstmnb(3,3)
       for_viri_cstmnb(3,1:natom)  = for_viri_cstmnb(3,1:natom) + fz(1:natom) 

    END DO

!   --- ADD FCSTMNB TO FORCE ---
    force(3,1:natom) = force(3,1:natom) + fz(1:natom)

!     +     +     +     +     +     +     +

  end subroutine calcstmnbp

end module cstmnb
