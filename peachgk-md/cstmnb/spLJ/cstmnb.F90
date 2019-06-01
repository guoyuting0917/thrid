!*****************************
!*  cstmnb.f90 Ver.1.3       *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-01-22 00:45:32 gota>

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
  character(5),public,save:: para_cstmnbtyp(maxncstmnbtyp)
                                    ! atom pair of custom NB type

  logical,public,save:: ifcalspLJ(maxnatmtyp,maxnatmtyp)
                                    ! if this pair is calculated in cstmnb

  real(8),public,save:: rcut_spLJ   ! cutoff length for LJ of spLJ
  real(8),public,save:: rrcut_spLJ  ! cutoff length for Coulomb of spLJ
!  integer,public,save:: nstep_bookspLJ(maxncstmnbtyp)
!                                    ! interval for bookkeeping
  real(8),public,save:: alpha_spLJ  ! Ewald alpha for spLJ

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
!!! This procedure is designed for calculating LJ and Coulomb (real space)
!!! potential for specific pairs of atom type.
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
!    real(8):: fcstmnb(3,natom)              ! custom NB force calculated here

    real(8):: pi                  ! = 3.14
    real(8):: const_fewr          ! = 2/ehta/sqrt(pi), used for Ewald 
                                !   real force
    real(8):: chrgij              ! = qi * qj
    real(8):: potij               ! temporal potential due to IJ
    real(8):: ftmp                ! temporal force  
    real(8):: fij(3)              ! force to be subtracted due to 
                                ! excluded atoms

    real(8):: Eij                 ! well depth 

    real(8):: rij(3)              ! rj-ri
    real(8):: rij_abs             ! |Rij|   
    real(8):: rij_abs2            ! |Rij|**2
    real(8):: rij_abs_inv1        ! |Rij|**-1
    real(8):: rij_abs_inv2        ! 1/|rij|**2
    real(8):: rtmp1
    real(8):: rtmp2, rtmp6, rtmp12           
    real(8):: rtmp_erf, rtmp_erfc           

    real(8):: box(3)              ! BOX size
    real(8):: box_inv(3)          ! inverse of BOX size

    real(8):: sqcut               ! rrcut*rrcut
    real(8):: sqcut2              ! rcut*rcut

    integer:: itype, jtype        ! VDW atom types
    integer:: i,j                 ! do loop index
    integer:: j1, j2              ! limits of do loop
    integer:: jj                  ! atom index of the interacting atom

! FUNCTIONS:
    real(8):: erf                 ! error function
!    real(8),external:: erf             ! error function
    real(8):: SPLINE_INTERPFUNC   ! spline interpolation function

!     +     +     +     +     +     +     +


!---- some initialization

!    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    pi = dacos(-1.0d0)
    const_fewr = 2.0d0 * alpha_spLJ / sqrt(pi)
                                ! coefficient used for calculation of
                                ! Ewald real

    box(1) = xcel
    box(2) = ycel
    box(3) = zcel
    box_inv(1:3) = 1.0d0/box(1:3)

    sqcut = rrcut_spLJ * rrcut_spLJ
    sqcut2 = rcut_spLJ * rcut_spLJ
    if (sqcut2 > sqcut) sqcut = sqcut2

!     --- CALCULATE NONBONDED FORCE ---

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, NATOM - 1
    looplast = natom - 1

    DO i = loopinit, looplast, loopstep

       j1 = nb_indexall(i)              ! the first and
       j2 = nb_indexall(i+loopstep) - 1 ! the last atoms interacting
                                        ! with i

       itype = atmindex(i)       ! VDW type of atom(i)

       IF (nb_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = nb_listall(j)     ! atom index of interacting atom
          jtype      = atmindex(jj) ! VDW type of atom(jj)

          if (.not. ifcalspLJ(itype,jtype)) cycle
                                               ! check if this pair calculated

          rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!         --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!         * box_inv rather than /box accellerates the computation

!         --- EWALD REAL ---

          rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2
                                ! |rij|^2
!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2

!!! if calculate L-J potential only
#if defined(_LJ_ONLY)
        ! no execution
#else
          rij_abs       = sqrt(rij_abs2) ! |rij| 
          rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

!         --- ewald real space cutoff ---

          if (rij_abs <= rrcut_spLJ) then

!---- Coulomb potential by Fennell method
#if defined(_spLJ_FENNELL)

!!! Spline interpolation for ewald real
#if defined(_NOT_SPLINTERP)
             rtmp1 = rij_abs * alpha_spLJ ! |rij|*alpha
             rtmp2 = rtmp1 * rtmp1 ! |rij|^2*alpha^2

             rtmp_erfc  = 1.0d0 - erf(rtmp1) ! erfc(|rij|*alpha)
             chrgij     = atmchrg(i) * atmchrg(jj)
             potij      = chrgij * rij_abs_inv1

!             ftmp       =  potij * (rtmp_erfc  * rij_abs_inv2   &
!                  &      + const_fewr * rij_abs_inv1 * dexp(-rtmp2))

                          ! =   Qi * Qj * [erfc(|rij|*alpha)/|rij|^3 
                          ! + 2*alpha/sqrt(pi)/|rij|^2 exp(-|rij|^2*alpha^2)]

             potij      = potij * rtmp_erfc
                                ! qi qj erfc(|rij|/ehta) / |rij|
#else
             chrgij = atmchrg(i) * atmchrg(jj)
!             ftmp = SPLINE_INTERPFUNC(rij_abs,spltbl_real,   &
!              &                       spl_b_real,spl_c_real,spl_d_real)   &
!              &   * chrgij
             potij = SPLINE_INTERPFUNC(rij_abs,spltbl_realpot,   &
                  &                 spl_b_realpot,spl_c_realpot,spl_d_realpot) &
                  & * chrgij
#endif

!            --- Fennell method
             potij = potij + chrgij   &
                  &     * (-fennell_shift   &
                  &       + fennell_damp*(rij_abs - rrcut_spLJ))
!             ftmp = ftmp - chrgij * fennell_damp * rij_abs_inv1

!             fij(1:3)     = ftmp * rij(1:3)
!
!             felc(1:3,i)  = felc(1:3,i)  + fij(1:3)
!             felc(1:3,jj) = felc(1:3,jj) - fij(1:3)

!---- pure 1/r potential (default)
#else
             potij      = atmchrg(i) * atmchrg(jj) * rij_abs_inv1

!             ftmp = potij * rij_abs_inv2
!
!             fij(1:3) = ftmp * rij(1:3)
#endif

             pot_cstmnb    = pot_cstmnb + potij

          end if

#endif
!!! if calculate L-J potential only

!         --- VDW ---

          if (rij_abs2 <= sqcut2) then

             if (inter_inttyp(itype,jtype) /= INTTYPE_VDW) cycle

             rtmp2      = rij_abs_inv2 * vdw_radij(itype,jtype)**2
             rtmp6      = rtmp2 * rtmp2 * rtmp2
             rtmp12     = rtmp6 * rtmp6
             Eij        = vdw_welij(itype,jtype)
             potij      = Eij * (rtmp12 - 2.0d0 * rtmp6)

             pot_cstmnb    = pot_cstmnb + potij

!             ftmp       = 12.0d0 * Eij * (rtmp12 - rtmp6)   &
!                  &     * rij_abs_inv2
!
!             fij(1:3)     = ftmp * rij(1:3)
!
!             fvdw(1:3,i)  = fvdw(1:3,i)  + fij(1:3)
!             fvdw(1:3,jj) = fvdw(1:3,jj) - fij(1:3)

          end if

       END DO
    END DO

! --- ADD FCSTMNB TO FORCE ---
!    force(1:3,1:natom) = force(1:3,1:natom) + fcstmnb(1:3,1:natom)

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
!    real(8):: fcstmnb(3,natom)              ! custom NB force calculated here

    real(8):: pi                  ! = 3.14
    real(8):: const_fewr          ! = 2/ehta/sqrt(pi), used for Ewald 
                                !   real force
    real(8):: chrgij              ! = qi * qj
    real(8):: potij               ! temporal potential due to IJ
    real(8):: ftmp                ! temporal force  
    real(8):: fij(3)              ! force to be subtracted due to 
                                ! excluded atoms

    real(8):: Eij                 ! well depth 

    real(8):: rij(3)              ! rj-ri
    real(8):: rij_abs             ! |Rij|   
    real(8):: rij_abs2            ! |Rij|**2
    real(8):: rij_abs_inv1        ! |Rij|**-1
    real(8):: rij_abs_inv2        ! 1/|rij|**2
    real(8):: rtmp1
    real(8):: rtmp2, rtmp6, rtmp12           
    real(8):: rtmp_erf, rtmp_erfc           

    real(8):: box(3)              ! BOX size
    real(8):: box_inv(3)          ! inverse of BOX size

    real(8):: sqcut               ! rrcut*rrcut
    real(8):: sqcut2              ! rcut*rcut

    integer:: itype, jtype        ! VDW atom types
    integer:: i,j,n               ! do loop index
    integer:: j1, j2              ! limits of do loop
    integer:: jj                  ! atom index of the interacting atom

! FUNCTIONS:
    real(8):: erf                 ! error function
!    real(8),external:: erf             ! error function
    real(8):: SPLINE_INTERPFUNC   ! spline interpolation function

!     +     +     +     +     +     +     +

!---- some initialization

!    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    pi = dacos(-1.0d0)
    const_fewr = 2.0d0 * alpha_spLJ / sqrt(pi)
                                ! coefficient used for calculation of
                                ! Ewald real

    box(1) = xcel
    box(2) = ycel
    box(3) = zcel
    box_inv(1:3) = 1.0d0/box(1:3)

    sqcut = rrcut_spLJ * rrcut_spLJ
    sqcut2 = rcut_spLJ * rcut_spLJ
    if (sqcut2 > sqcut) sqcut = sqcut2

!     --- CALCULATE NONBONDED FORCE ---

    pot_virit_cstmnb(1:3,1:3) = 0.0d0

    pot_viri_cstmnb = 0.0d0

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, NATOM - 1
    looplast = natom - 1

    DO i = loopinit, looplast, loopstep

       j1 = nb_indexall(i)              ! the first and
       j2 = nb_indexall(i+loopstep) - 1 ! the last atoms interacting
                                        ! with i

       itype = atmindex(i)       ! VDW type of atom(i)

       IF (nb_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = nb_listall(j)     ! atom index of interacting atom
          jtype      = atmindex(jj) ! VDW type of atom(jj)

          if (.not. ifcalspLJ(itype,jtype)) cycle
                                               ! check if this pair calculated

          rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!         --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!         * box_inv rather than /box accellerates the computation

!         --- EWALD REAL ---

          rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2
                                ! |rij|^2
!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2

!!! if calculate L-J potential only
#if defined(_LJ_ONLY)
        ! no execution
#else
          rij_abs       = sqrt(rij_abs2) ! |rij| 
          rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

!         --- ewald real space cutoff ---

          if (rij_abs <= rrcut_spLJ) then

!---- Coulomb potential by Fennell method
#if defined(_spLJ_FENNELL)

!!! Spline interpolation for ewald real
#if defined(_NOT_SPLINTERP)
             rtmp1 = rij_abs * alpha_spLJ ! |rij|*alpha
             rtmp2 = rtmp1 * rtmp1 ! |rij|^2*alpha^2

             rtmp_erfc  = 1.0d0 - erf(rtmp1) ! erfc(|rij|*alpha)
             chrgij     = atmchrg(i) * atmchrg(jj)
             potij      = chrgij * rij_abs_inv1

             ftmp       =  potij * (rtmp_erfc  * rij_abs_inv2   &
                  &      + const_fewr * rij_abs_inv1 * dexp(-rtmp2))

                          ! =   Qi * Qj * [erfc(|rij|*alpha)/|rij|^3 
                          ! + 2*alpha/sqrt(pi)/|rij|^2 exp(-|rij|^2*alpha^2)]

             potij      = potij * rtmp_erfc
                                ! qi qj erfc(|rij|/ehta) / |rij|
#else
             chrgij = atmchrg(i) * atmchrg(jj)
             ftmp = SPLINE_INTERPFUNC(rij_abs,spltbl_real,   &
                  &                   spl_b_real,spl_c_real,spl_d_real)   &
                  &   * chrgij
             potij = SPLINE_INTERPFUNC(rij_abs,spltbl_realpot,   &
                  &                 spl_b_realpot,spl_c_realpot,spl_d_realpot) &
                  &   * chrgij
#endif

!            --- Fennell method
             potij = potij + chrgij   &
                  &     * (-fennell_shift   &
                  &       + fennell_damp*(rij_abs - rrcut_spLJ))
             ftmp = ftmp - chrgij * fennell_damp * rij_abs_inv1

             fij(1:3)     = ftmp * rij(1:3)
!
!             felc(1:3,i)  = felc(1:3,i)  + fij(1:3)
!             felc(1:3,jj) = felc(1:3,jj) - fij(1:3)

!---- pure 1/r potential (default)
#else
             potij      = atmchrg(i) * atmchrg(jj) * rij_abs_inv1

             ftmp = potij * rij_abs_inv2

             fij(1:3) = ftmp * rij(1:3)
#endif

             pot_cstmnb    = pot_cstmnb + potij

!            --- calculate virial tensor ---
             do n=1,3
                pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                     &                  + fij(1:3)*rij(n)
             end do

             pot_viri_cstmnb =  pot_viri_cstmnb + potij

             for_viri_cstmnb(1:3,i)  = for_viri_cstmnb(1:3,i)  + fij(1:3)
             for_viri_cstmnb(1:3,jj) = for_viri_cstmnb(1:3,jj) - fij(1:3)

          end if

#endif
!!! if calculate L-J potential only

!         --- VDW ---

          if (rij_abs2 <= sqcut2) then

             if (inter_inttyp(itype,jtype) /= INTTYPE_VDW) cycle

             rtmp2      = rij_abs_inv2 * vdw_radij(itype,jtype)**2
             rtmp6      = rtmp2 * rtmp2 * rtmp2
             rtmp12     = rtmp6 * rtmp6
             Eij        = vdw_welij(itype,jtype)
             potij      = Eij * (rtmp12 - 2.0d0 * rtmp6)

             pot_cstmnb    = pot_cstmnb + potij

             ftmp       = 12.0d0 * Eij * (rtmp12 - rtmp6)   &
                  &     * rij_abs_inv2

             fij(1:3)     = ftmp * rij(1:3)

!             fvdw(1:3,i)  = fvdw(1:3,i)  + fij(1:3)
!             fvdw(1:3,jj) = fvdw(1:3,jj) - fij(1:3)

!            --- calculate virial tensor ---
             do n=1,3
                pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                     &                  + fij(1:3)*rij(n)
             end do

             pot_viri_cstmnb = pot_viri_cstmnb &
                  &          + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

             for_viri_cstmnb(1:3,i)  = for_viri_cstmnb(1:3,i)  + fij(1:3)
             for_viri_cstmnb(1:3,jj) = for_viri_cstmnb(1:3,jj) - fij(1:3)

          end if

       END DO
    END DO

!   --- ADD FCSTMNB TO FORCE ---
!    force(1:3,1:natom) = force(1:3,1:natom) + fcstmnb(1:3,1:natom)
!    for_viri_cstmnb(1:3,1:natom) = for_viri_cstmnb(1:3,1:natom) &
!         &                       + fcstmnb(1:3,1:natom)

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
!    real(8):: fcstmnb(3,natom)              ! custom NB force calculated here

    real(8):: pi                  ! = 3.14
    real(8):: const_fewr          ! = 2/ehta/sqrt(pi), used for Ewald 
                                  !   real force
    real(8):: chrgij              ! = qi * qj
    real(8):: potij               ! temporal potential due to IJ
    real(8):: ftmp                ! temporal force  
    real(8):: fij(3)              ! force to be subtracted due to 
                                  ! excluded atoms
    real(8):: fji(3)              ! force to be subtracted due to 
                                  ! excluded atoms = -fij

    real(8):: Eij                 ! well depth 

    real(8):: rij(3)              ! rj-ri
    real(8):: rij_abs             ! |Rij|   
    real(8):: rij_abs2            ! |Rij|**2
    real(8):: rij_abs_inv1        ! |Rij|**-1
    real(8):: rij_abs_inv2        ! 1/|rij|**2
    real(8):: rtmp1
    real(8):: rtmp2, rtmp6, rtmp12           
    real(8):: rtmp_erf, rtmp_erfc           

    real(8):: box(3)              ! BOX size
    real(8):: box_inv(3)          ! inverse of BOX size

    real(8):: sqcut               ! rrcut*rrcut
    real(8):: sqcut2              ! rcut*rcut

    integer:: itype, jtype        ! VDW atom types
    integer:: i,j,n               ! do loop index
    integer:: j1, j2              ! limits of do loop
    integer:: jj                  ! atom index of the interacting atom

    real(8):: nbodycoeff          ! n-body coefficient for virial term of heatf

! FUNCTIONS:
    real(8):: erf                 ! error function
!    real(8),external:: erf             ! error function
    real(8):: SPLINE_INTERPFUNC   ! spline interpolation function

!     +     +     +     +     +     +     +

!---- some initialization

!    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    pi = dacos(-1.0d0)
    const_fewr = 2.0d0 * alpha_spLJ / sqrt(pi)
                                ! coefficient used for calculation of
                                ! Ewald real

    box(1) = xcel
    box(2) = ycel
    box(3) = zcel
    box_inv(1:3) = 1.0d0/box(1:3)

    sqcut = rrcut_spLJ * rrcut_spLJ
    sqcut2 = rcut_spLJ * rcut_spLJ
    if (sqcut2 > sqcut) sqcut = sqcut2

!     --- CALCULATE NONBONDED FORCE ---

    pot_virit_cstmnb(1:3,1:3) = 0.0d0

    pot_viri_cstmnb = 0.0d0

    nbodycoeff = 0.5d0        ! = 1/2

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, NATOM - 1
    looplast = natom - 1

    DO i = loopinit, looplast, loopstep

       j1 = nb_indexall(i)              ! the first and
       j2 = nb_indexall(i+loopstep) - 1 ! the last atoms interacting
                                        ! with i

       itype = atmindex(i)       ! VDW type of atom(i)

       IF (nb_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = nb_listall(j)     ! atom index of interacting atom
          jtype      = atmindex(jj) ! VDW type of atom(jj)

          if (.not. ifcalspLJ(itype,jtype)) cycle
                                               ! check if this pair calculated

          rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!         --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!         * box_inv rather than /box accellerates the computation

!         --- EWALD REAL ---

          rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2
                                ! |rij|^2
!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2

!!! if calculate L-J potential only
#if defined(_LJ_ONLY)
        ! no execution
#else
          rij_abs       = sqrt(rij_abs2) ! |rij| 
          rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

!         --- ewald real space cutoff ---

          if (rij_abs <= rrcut_spLJ) then

!---- Coulomb potential by Fennell method
#if defined(_spLJ_FENNELL)

!!! Spline interpolation for ewald real
#if defined(_NOT_SPLINTERP)
             rtmp1 = rij_abs * alpha_spLJ ! |rij|*alpha
             rtmp2 = rtmp1 * rtmp1 ! |rij|^2*alpha^2

             rtmp_erfc  = 1.0d0 - erf(rtmp1) ! erfc(|rij|*alpha)
             chrgij     = atmchrg(i) * atmchrg(jj)
             potij      = chrgij * rij_abs_inv1

             ftmp       =  potij * (rtmp_erfc  * rij_abs_inv2   &
                  &      + const_fewr * rij_abs_inv1 * dexp(-rtmp2))

                          ! =   Qi * Qj * [erfc(|rij|*alpha)/|rij|^3 
                          ! + 2*alpha/sqrt(pi)/|rij|^2 exp(-|rij|^2*alpha^2)]

             potij      = potij * rtmp_erfc
                                ! qi qj erfc(|rij|/ehta) / |rij|
#else
             chrgij = atmchrg(i) * atmchrg(jj)
             ftmp = SPLINE_INTERPFUNC(rij_abs,spltbl_real,   &
                  &                   spl_b_real,spl_c_real,spl_d_real)   &
                  &   * chrgij
             potij = SPLINE_INTERPFUNC(rij_abs,spltbl_realpot,   &
                  &                 spl_b_realpot,spl_c_realpot,spl_d_realpot) &
                  &   * chrgij
#endif

!            --- Fennell method
             potij = potij + chrgij   &
                  &     * (-fennell_shift   &
                  &       + fennell_damp*(rij_abs - rrcut_spLJ))
             ftmp = ftmp - chrgij * fennell_damp * rij_abs_inv1

             fij(1:3)     = ftmp * rij(1:3)
!
!             felc(1:3,i)  = felc(1:3,i)  + fij(1:3)
!             felc(1:3,jj) = felc(1:3,jj) - fij(1:3)

!---- pure 1/r potential (default)
#else
             potij      = atmchrg(i) * atmchrg(jj) * rij_abs_inv1

             ftmp = potij * rij_abs_inv2

             fij(1:3) = ftmp * rij(1:3)
#endif

             pot_cstmnb    = pot_cstmnb + potij

!            --- calculate virial tensor ---

!            - for heat flux calculation
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

!            - for pressure calculation
             do n=1,3
                pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                     &                  + fij(1:3)*rij(n)
             end do

             pot_viri_cstmnb = pot_viri_cstmnb &
                  &          + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

             for_viri_cstmnb(1:3,i)  = for_viri_cstmnb(1:3,i)  + fij(1:3)
             for_viri_cstmnb(1:3,jj) = for_viri_cstmnb(1:3,jj) - fij(1:3)

          end if

#endif
!!! if calculate L-J potential only

!         --- VDW ---

          if (rij_abs2 <= sqcut2) then

             if (inter_inttyp(itype,jtype) /= INTTYPE_VDW) cycle

             rtmp2      = rij_abs_inv2 * vdw_radij(itype,jtype)**2
             rtmp6      = rtmp2 * rtmp2 * rtmp2
             rtmp12     = rtmp6 * rtmp6
             Eij        = vdw_welij(itype,jtype)
             potij      = Eij * (rtmp12 - 2.0d0 * rtmp6)

             pot_cstmnb    = pot_cstmnb + potij

             ftmp       = 12.0d0 * Eij * (rtmp12 - rtmp6)   &
                  &     * rij_abs_inv2

             fij(1:3)     = ftmp * rij(1:3)

!             fvdw(1:3,i)  = fvdw(1:3,i)  + fij(1:3)
!             fvdw(1:3,jj) = fvdw(1:3,jj) - fij(1:3)

!            --- calculate virial tensor ---

!            - for heat flux calculation
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

!            - for pressure calculation
             do n=1,3
                pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                     &                  + fij(1:3)*rij(n)
             end do

             pot_viri_cstmnb = pot_viri_cstmnb &
                  &          + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

             for_viri_cstmnb(1:3,i)  = for_viri_cstmnb(1:3,i)  + fij(1:3)
             for_viri_cstmnb(1:3,jj) = for_viri_cstmnb(1:3,jj) - fij(1:3)

          end if

       END DO
    END DO

!   --- ADD FCSTMNB TO FORCE ---
!    force(1:3,1:natom) = force(1:3,1:natom) + fcstmnb(1:3,1:natom)
!    for_viri_cstmnb(1:3,1:natom) = for_viri_cstmnb(1:3,1:natom) &
!         &                       + fcstmnb(1:3,1:natom)

!     +     +     +     +     +     +     +

  end subroutine calcstmnbp_hf
#endif

end module cstmnb
