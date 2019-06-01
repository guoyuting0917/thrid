!*****************************
!*  cstmnb.f90 Ver.1.3       *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-01-22 00:20:29 gota>

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

  integer,public,save:: nZPcond               ! number of all ZP cond atoms
  integer,public,save:: nZPcondlist(maxnatom) ! list of ZP cond atoms

  integer,public,save:: ZPcond1_indexall(maxnatom+maxnproc)
                              ! list of bookkeeped ZP cond atoms for plane 1
  integer,public,save,allocatable:: ZPcond1_listall(:)
!  integer,public,save:: ZPcond1_listall(maxnatom*maxcstmnblist)
                              ! list of neighboring atoms
  integer,public,save:: nlistZPcond1all     ! number of ZP cond1 NB list

  integer,public,save:: ZPcond2_indexall(maxnatom+maxnproc)
                              ! list of bookkeeped ZP cond atoms for plane 2
  integer,public,save,allocatable:: ZPcond2_listall(:)
!  integer,public,save:: ZPcond2_listall(maxnatom*maxcstmnblist)
                              ! list of neighboring atoms
  integer,public,save:: nlistZPcond2all     ! number of ZP cond1 NB list

  logical,public,save:: ifcellindex_ZPaniso ! cell index flag for ZP aniso term
  logical,public,save:: ifcellindex_ZPiso   ! cell index flag for ZP iso term

  integer,public,save:: ncstmnbtyp  ! number of cstmnb type
  character(5),public,save:: para_cstmnbtyp(maxncstmnbtyp)
                                    ! atom pair of custom NB type
  real(8),public,save:: rcut_ZPcond(maxncstmnbtyp) ! cutoff length for ZP cond
  real(8),public,save:: rcut_bookZPcond(maxncstmnbtyp)
                                    ! cutoff length for bookkeeping of ZP cond
  integer,public,save:: nstep_bookZPcond(maxncstmnbtyp)
                                    ! interval for bookkeeping of ZP cond
  real(8),public,save:: para_ZPint_p1(maxncstmnbtyp)
                                    ! surface position 1 of ZP cond term
  real(8),public,save:: para_ZPint_p2(maxncstmnbtyp)
                                    ! surface position 2 of ZP cond term
  real(8),public,save:: rcut_ZPaniso(maxncstmnbtyp)
                                    ! cutoff length for ZP anisotropic term
  real(8),public,save:: rcut_bookZPaniso(maxncstmnbtyp)
                                    ! cutoff length for bookkeeping of ZP aniso
  integer,public,save:: nstep_bookZPaniso(maxncstmnbtyp)
                                    ! interval for bookkeeping of ZP aniso
  real(8),public,save:: para_alphaZP_aniso(maxnatmtyp,maxnatmtyp)
                                    ! alpha parameter of ZP anisotropic term

  real(8),public,save:: rcut_ZPiso(maxncstmnbtyp) ! cutoff length for ZP iso
  real(8),public,save:: rcut_bookZPiso(maxncstmnbtyp) 
                                    ! cutoff length for bookkeeping of ZP iso
  integer,public,save:: nstep_bookZPiso(maxncstmnbtyp)
                                    ! interval for bookkeeping of ZP iso
  real(8),public,save:: para_cZP_iso(maxnatmtyp,maxnatmtyp)
                                    ! C10 parameter of ZP isotropic term
  real(8),public,save:: para_welZP(maxnatmtyp,maxnatmtyp)
                                    ! vdw well depth of ZP potential
  real(8),public,save:: para_radZP(maxnatmtyp,maxnatmtyp)
                                    ! vdw radius of ZP potential

  logical,public,save:: ifrenew_msurf(maxncstmnbtyp)
                                    ! if renew the mirror surface during MD

!!! Hard coding for indexing surface atoms to renew the mirror surface position
  integer,public,parameter:: matomtyp_ZPint_p1 = 2
                                    ! matom type for the mirror surface 1
  integer,public,parameter:: matomtyp_ZPint_p2 = 3
                                    ! matom type for the mirror surface 1
!!! indexes below starts from 1 and should be counted for each molecular species
  integer,public,parameter:: matom_ZPint_p1_ini = 1
                                    ! first atom index of surface 1
  integer,public,parameter:: matom_ZPint_p1_end = 320
                                    ! last atom index of surface 1
  integer,public,parameter:: matom_ZPint_p2_ini = 1
                                    ! first atom index of surface 2
  integer,public,parameter:: matom_ZPint_p2_end = 320
                                    ! last atom index of surface 2

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
!!! This procedure is designed for calculating ZP potential (Zhu and Philpott,
!!! J. Chem. Phys., 1994) between water and metal surfaces.
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
!      real(8):: fcstmnball(3,natom) ! custom NB force calculated here
!MPI                                       local(in MPI sense) field
    real(8):: potij                      ! temporal potential due to IJ
    real(8):: fcoeff,ftmp                ! temporal force
    real(8):: fij(3)                     ! force to be subtracted due to
                                         ! excluded atoms

    real(8):: rij(3)                     ! rj-ri
    real(8):: rij_abs                    ! |Rij|
    real(8):: xyij_abs2                  ! xij**2 + yij**2
    real(8):: zij_abs2                   ! zij**2
    real(8):: rij_abs2                   ! |Rij|**2
    real(8):: rij_abs_inv2               ! |Rij|**-2
    real(8):: rij_abs_inv1               ! |Rij|**-1

    real(8):: alpha2                     ! alpha**2
    real(8):: alpha2_inv                 ! 1/alpha**2
    real(8):: sigma2                     ! sigma**2

    real(8):: rij_abs2_alpha
    real(8):: rij_abs2_a_inv
    real(8):: rij_abs2_s_a_inv,rij_abs6_s_a_inv,rij_abs12_s_a_inv
    real(8):: rij_abs2_alpha_inv
    real(8):: rij_abs2_ainv_inv
    real(8):: rij_abs2_s_ainv_inv,rij_abs6_s_ainv_inv
    real(8):: rij_abs2_s_inv
    real(8):: rij_abs4_s_inv,rij_abs8_s_inv,rij_abs10_s_inv

    real(8):: rtmp6, rtmp3

    real(8):: box(3)                     ! BOX size
    real(8):: box_inv(3)                 ! inverse of BOX size

    real(8):: sqcut                      ! rcut * rcut
    real(8):: sqcut_aniso,sqcut_iso      ! rcut * rcut

    integer:: itype, jtype               ! VDW atom types
    integer:: i,j,k                      ! do loop index
    integer:: j1, j2                     ! limits of do loop
    integer:: jj                         ! atom index of the interacting atom
    integer:: imole

    real(8):: ZPint_p1_ave, ZPint_p2_ave 
                            ! instantaneous averaged position of mirror surfaces

    integer:: cstmnb1

!     +     +     +     +     +     +     +


!---- some initialization

    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    box(1) = xcel
    box(2) = ycel
    box(3) = zcel
    box_inv(1:3) = 1.0d0/box(1:3)

!---- CALCULATE ZP CONDUCTION FORCE (mirror force) ----

    sqcut = rcut_ZPcond(1) * rcut_ZPcond(1)

!---- update the position of mirror surface #1

    if (ifrenew_msurf(1)) then

       ZPint_p1_ave = 0.0d0

       imole = npoly + nwater
       do i = 1, matomtyp_ZPint_p1 - 1
          imole = imole + nmatomtyp(i)
       end do

       do i = matom_ZPint_p1_ini, matom_ZPint_p1_end
          j1 = molept_index(imole+i)         ! first monatomic mole.
          j2 = molept_index(imole+i+1) - 1

          do j= j1, j2
             k = molept_list(j)

             ZPint_p1_ave = ZPint_p1_ave + atmcor(3,k)

          end do

       end do

       para_ZPint_p1(1) = ZPint_p1_ave &
            &           / dble(matom_ZPint_p1_end - matom_ZPint_p1_ini + 1)

    end if

!---- surface #1
! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO cstmnb1 = 1, nZPcond
    looplast = nZPcond

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = ZPcond1_indexall(cstmnb1)              ! the first and
       j2 = ZPcond1_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting

       i = nZPcondlist(cstmnb1)
       itype = atmindex(i)               ! atmtype of atom(i)

       IF (ZPcond1_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = ZPcond1_listall(j)         ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:2) = atmcor(1:2,i) - atmcor(1:2,jj)
          rij(3) = atmcor(3,i) - (2.0d0*para_ZPint_p1(1) - atmcor(3,jj))
                                          ! = zi - zj' (mirror image)

!         --- periodic boundary ---

          rij(1:2) = rij(1:2) - box(1:2) * anint(rij(1:2)*box_inv(1:2))
                                          ! P.B.C. only for x and y direction

          rij_abs2 = rij(1)**2 + rij(2)**2 + rij(3)**2   ! = |rij|^2

!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2 = 1.0d0 / rij_abs2   ! = 1/|rij|^2
          rij_abs_inv1 = sqrt(rij_abs_inv2) ! = 1/|rij|

          potij = - atmchrg(i) * atmchrg(jj) * rij_abs_inv1 ! = qi*(-qj)/|rij|
          pot_cstmnb = pot_cstmnb + potij * 0.5d0   ! 1/2 is needed for energy

          fcstmnb(1:3,i) = fcstmnb(1:3,i) + potij * rij_abs_inv2 * rij(1:3)
                                            ! = qi*(-qj)/|rij|^3 * rij
       END DO
    END DO

!---- update the position of mirror surface #2

    if (ifrenew_msurf(1)) then

       ZPint_p2_ave = 0.0d0

       imole = npoly + nwater
       do i = 1, matomtyp_ZPint_p2 - 1
          imole = imole + nmatomtyp(i)
       end do

       do i = matom_ZPint_p2_ini, matom_ZPint_p2_end
          j1 = molept_index(imole+i)         ! first monatomic mole.
          j2 = molept_index(imole+i+1) - 1

          do j= j1, j2
             k = molept_list(j)

             ZPint_p2_ave = ZPint_p2_ave + atmcor(3,k)

          end do

       end do

       para_ZPint_p2(1) = ZPint_p2_ave &
            &           / dble(matom_ZPint_p2_end - matom_ZPint_p2_ini + 1)

    end if

!---- surface #2
! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO cstmnb1 = 1, nZPcond
    looplast = nZPcond

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = ZPcond2_indexall(cstmnb1)              ! the first and
       j2 = ZPcond2_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting

       i = nZPcondlist(cstmnb1)
       itype = atmindex(i)               ! atmtype of atom(i)

       IF (ZPcond2_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = ZPcond2_listall(j)         ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:2) = atmcor(1:2,i) - atmcor(1:2,jj)
          rij(3) = atmcor(3,i) - (2.0d0*para_ZPint_p2(1) - atmcor(3,jj))
                                          ! = zi - zj' (mirror image)

!         --- periodic boundary ---

          rij(1:2) = rij(1:2) - box(1:2) * anint(rij(1:2)*box_inv(1:2))
                                          ! P.B.C. only for x and y direction

          rij_abs2 = rij(1)**2 + rij(2)**2 + rij(3)**2   ! = |rij|^2

!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2 = 1.0d0 / rij_abs2   ! = 1/|rij|^2
          rij_abs_inv1 = sqrt(rij_abs_inv2) ! = 1/|rij|

          potij = - atmchrg(i) * atmchrg(jj) * rij_abs_inv1 ! = qi*(-qj)/|rij|
          pot_cstmnb = pot_cstmnb + potij * 0.5d0   ! 1/2 is needed for energy

          fcstmnb(1:3,i) = fcstmnb(1:3,i) + potij * rij_abs_inv2 * rij(1:3)
                                            ! = qi*(-qj)/|rij|^3 * rij
       END DO
    END DO


!---- CALCULATE ZP ANISO- and ISOTROPIC SHORT RANGE FORCE ----

    sqcut_aniso = rcut_ZPaniso(1) * rcut_ZPaniso(1)
    sqcut_iso = rcut_ZPiso(1) * rcut_ZPiso(1)

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO cstmnb1 = 1, ncstmnb - 1
    looplast = ncstmnb - 1

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = cstmnb_indexall(cstmnb1)              ! the first and
       j2 = cstmnb_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting
                                                  ! with i

       i = ncstmnblist(cstmnb1)
       itype = atmindex(i)                  ! atmtype of atom(i)

       IF (cstmnb_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = cstmnb_listall(j)          ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!          --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3) * anint(rij(1:3)*box_inv(1:3))

!         * box_inv rather than /box accellerates the computation

          xyij_abs2 = rij(1)**2 + rij(2)**2
          zij_abs2 = rij(3)**2
          rij_abs2  = xyij_abs2 + zij_abs2 ! |rij|^2

!         --- cut off for anisotropic force ---
          if (rij_abs2 > sqcut_aniso) cycle

          alpha2 = para_alphaZP_aniso(itype,jtype)**2
          alpha2_inv = 1.0d0 / alpha2
          sigma2 = para_radZP(itype,jtype)**2

          rij_abs2_alpha = alpha2*xyij_abs2 + zij_abs2
                                          ! (xij^2+yij^2)*alpha^2 + zij^2
          rij_abs2_a_inv = 1.0d0 / rij_abs2_alpha
                                    ! =1/((xij^2+yij^2)*alpha^2 + zij^2)
          rij_abs2_s_a_inv = sigma2 * rij_abs2_a_inv
                                    ! =sigma^2/((xij^2+yij^2)*alpha^2 + zij^2)
          rij_abs6_s_a_inv = rij_abs2_s_a_inv * rij_abs2_s_a_inv &
               &           * rij_abs2_s_a_inv
                                ! =(sigma^2/((xij^2+yij^2)*alpha^2 + zij^2))^3
          rij_abs12_s_a_inv = rij_abs6_s_a_inv * rij_abs6_s_a_inv
                                ! =(sigma^2/((xij^2+yij^2)*alpha^2 + zij^2))^6

          rij_abs2_alpha_inv = alpha2_inv*xyij_abs2 + zij_abs2
                                          ! (xij^2+yij^2)/alpha^2 + zij^2
          rij_abs2_ainv_inv = 1.0d0 / rij_abs2_alpha_inv
                                    ! =1/((xij^2+yij^2)/alpha^2 + zij^2)
          rij_abs2_s_ainv_inv = sigma2 * rij_abs2_ainv_inv
                                    ! =sigma^2/((xij^2+yij^2)/alpha^2 + zij^2)
          rij_abs6_s_ainv_inv = rij_abs2_s_ainv_inv * rij_abs2_s_ainv_inv &
               &              * rij_abs2_s_ainv_inv
                                ! =(sigma^2/((xij^2+yij^2)/alpha^2 + zij^2))^3

          pot_cstmnb  = pot_cstmnb + 4.0d0*para_welZP(itype,jtype) &
               &                   * (rij_abs12_s_a_inv - rij_abs6_s_ainv_inv)
                                 ! anisotropic potential

          rtmp6 = rij_abs12_s_a_inv * rij_abs2_a_inv  ! f_6 term

          rtmp3 = rij_abs6_s_ainv_inv * rij_abs2_ainv_inv ! f_3 term

          fcoeff = 24.0d0 * para_welZP(itype,jtype)
          ftmp = fcoeff * (2.0d0*alpha2*rtmp6 - alpha2_inv*rtmp3)
                                 ! =24*epsilon*(2*alpha^2*f_6 - 1/alpha^2*f_3)

          fij(1:2)     = ftmp * rij(1:2)

          fij(3) = fcoeff * (2.0d0*rtmp6 - rtmp3) * rij(3)
                                 ! =24*epsilon*(2*f_6 - f_3) * zij

          fcstmnb(1:3,i)  = fcstmnb(1:3,i)  + fij(1:3)
          fcstmnb(1:3,jj) = fcstmnb(1:3,jj) - fij(1:3)

!         --- cut off for isotropic force ---
          if (rij_abs2 > sqcut_iso) cycle

          rij_abs_inv2 = 1.0d0 / rij_abs2                  ! = 1/|rij|^2
          rij_abs2_s_inv = sigma2 * rij_abs_inv2           ! =sigma^2/|rij|^2
          rij_abs4_s_inv = rij_abs2_s_inv * rij_abs2_s_inv ! =sigma^4/|rij|^4
          rij_abs8_s_inv = rij_abs4_s_inv * rij_abs4_s_inv ! =sigma^8/|rij|^8
          rij_abs10_s_inv = rij_abs8_s_inv * rij_abs2_s_inv
                                         ! =sigma^10/|rij|^10

          potij  = -4.0d0*para_welZP(itype,jtype)*para_cZP_iso(itype,jtype) &
               & * rij_abs10_s_inv  ! =-4*epsilon*C10*sigma^10/|rij|^10

          pot_cstmnb = pot_cstmnb + potij

          ftmp = 10.0d0 * potij * rij_abs_inv2 ! =10*phi_ij/|rij|^2
          fij(1:3) = ftmp * rij(1:3) ! =10*phi_ij/|rij|^2*rij

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
!      real(8):: fcstmnball(3,natom) ! custom NB force calculated here
!MPI                                       local(in MPI sense) field
    real(8):: potij                      ! temporal potential due to IJ
    real(8):: fcoeff,ftmp                ! temporal force
    real(8):: fij(3)                     ! force to be subtracted due to
                                         ! excluded atoms

    real(8):: rij(3)                     ! rj-ri
    real(8):: rij_abs                    ! |Rij|
    real(8):: xyij_abs2                  ! xij**2 + yij**2
    real(8):: zij_abs2                   ! zij**2
    real(8):: rij_abs2                   ! |Rij|**2
    real(8):: rij_abs_inv2               ! |Rij|**-2
    real(8):: rij_abs_inv1               ! |Rij|**-1

    real(8):: alpha2                     ! alpha**2
    real(8):: alpha2_inv                 ! 1/alpha**2
    real(8):: sigma2                     ! sigma**2

    real(8):: rij_abs2_alpha
    real(8):: rij_abs2_a_inv
    real(8):: rij_abs2_s_a_inv,rij_abs6_s_a_inv,rij_abs12_s_a_inv
    real(8):: rij_abs2_alpha_inv
    real(8):: rij_abs2_ainv_inv
    real(8):: rij_abs2_s_ainv_inv,rij_abs6_s_ainv_inv
    real(8):: rij_abs2_s_inv
    real(8):: rij_abs4_s_inv,rij_abs8_s_inv,rij_abs10_s_inv

    real(8):: rtmp6, rtmp3

    real(8):: box(3)                      ! BOX size
    real(8):: box_inv(3)                  ! inverse of BOX size

    real(8):: sqcut                      ! rcut * rcut
    real(8):: sqcut_aniso,sqcut_iso      ! rcut * rcut

    integer:: itype, jtype                ! VDW atom types
    integer:: i,j,k,n                     ! do loop index
    integer:: j1, j2                      ! limits of do loop
    integer:: jj                          ! atom index of the interacting atom
    integer:: imole

    real(8):: ZPint_p1_ave, ZPint_p2_ave 
                            ! instantaneous averaged position of mirror surfaces

    integer:: cstmnb1

!     +     +     +     +     +     +     +


!---- some initialization

    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    pot_virit_cstmnb(1:3,1:3) = 0.0d0
    pot_viri_cstmnb = 0.0d0

    box(1) = xcel
    box(2) = ycel
    box(3) = zcel
    box_inv(1:3) = 1.0d0/box(1:3)

!---- CALCULATE ZP CONDUCTION FORCE (mirror force) ----

    sqcut = rcut_ZPcond(1) * rcut_ZPcond(1)

!---- update the position of mirror surface #1

    if (ifrenew_msurf(1)) then

       ZPint_p1_ave = 0.0d0

       imole = npoly + nwater
       do i = 1, matomtyp_ZPint_p1 - 1
          imole = imole + nmatomtyp(i)
       end do

       do i = matom_ZPint_p1_ini, matom_ZPint_p1_end
          j1 = molept_index(imole+i)         ! first monatomic mole.
          j2 = molept_index(imole+i+1) - 1

          do j= j1, j2
             k = molept_list(j)

             ZPint_p1_ave = ZPint_p1_ave + atmcor(3,k)

          end do

       end do

       para_ZPint_p1(1) = ZPint_p1_ave &
            &           / dble(matom_ZPint_p1_end - matom_ZPint_p1_ini + 1)

    end if

!---- surface #1
! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO cstmnb1 = 1, nZPcond
    looplast = nZPcond

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = ZPcond1_indexall(cstmnb1)              ! the first and
       j2 = ZPcond1_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting

       i = nZPcondlist(cstmnb1)
       itype = atmindex(i)               ! atmtype of atom(i)

       IF (ZPcond1_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = ZPcond1_listall(j)         ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:2) = atmcor(1:2,i) - atmcor(1:2,jj)
          rij(3) = atmcor(3,i) - (2.0d0*para_ZPint_p1(1) - atmcor(3,jj))
                                          ! = zi - zj' (mirror image)

!         --- periodic boundary ---

          rij(1:2) = rij(1:2) - box(1:2) * anint(rij(1:2)*box_inv(1:2))
                                          ! P.B.C. only for x and y direction

          rij_abs2 = rij(1)**2 + rij(2)**2 + rij(3)**2   ! = |rij|^2

!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2 = 1.0d0 / rij_abs2   ! = 1/|rij|^2
          rij_abs_inv1 = sqrt(rij_abs_inv2) ! = 1/|rij|

          potij = - atmchrg(i) * atmchrg(jj) * rij_abs_inv1 ! = qi*(-qj)/|rij|
          pot_cstmnb = pot_cstmnb + potij * 0.5d0   ! 1/2 is needed for energy

          fij(1:3) = potij * rij_abs_inv2 * rij(1:3)
          fcstmnb(1:3,i) = fcstmnb(1:3,i) + fij(1:3)
                                            ! = qi*(-qj)/|rij|^3 * rij

!         --- calculate virial tensor ---
          do n=1,3
             pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                  &                  + fij(1:3)*atmcor(n,i)
          end do

          pot_viri_cstmnb = pot_viri_cstmnb + fij(1)*atmcor(1,i) &
               &          + fij(2)*atmcor(2,i) + fij(3)*atmcor(3,i)

       END DO
    END DO

!---- update the position of mirror surface #2

    if (ifrenew_msurf(1)) then

       ZPint_p2_ave = 0.0d0

       imole = npoly + nwater
       do i = 1, matomtyp_ZPint_p2 - 1
          imole = imole + nmatomtyp(i)
       end do

       do i = matom_ZPint_p2_ini, matom_ZPint_p2_end
          j1 = molept_index(imole+i)         ! first monatomic mole.
          j2 = molept_index(imole+i+1) - 1

          do j= j1, j2
             k = molept_list(j)

             ZPint_p2_ave = ZPint_p2_ave + atmcor(3,k)

          end do

       end do

       para_ZPint_p2(1) = ZPint_p2_ave &
            &           / dble(matom_ZPint_p2_end - matom_ZPint_p2_ini + 1)

    end if

!---- surface #2
! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO cstmnb1 = 1, nZPcond
    looplast = nZPcond

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = ZPcond2_indexall(cstmnb1)              ! the first and
       j2 = ZPcond2_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting

       i = nZPcondlist(cstmnb1)
       itype = atmindex(i)               ! atmtype of atom(i)

       IF (ZPcond2_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = ZPcond2_listall(j)         ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:2) = atmcor(1:2,i) - atmcor(1:2,jj)
          rij(3) = atmcor(3,i) - (2.0d0*para_ZPint_p2(1) - atmcor(3,jj))
                                          ! = zi - zj' (mirror image)

!         --- periodic boundary ---

          rij(1:2) = rij(1:2) - box(1:2) * anint(rij(1:2)*box_inv(1:2))
                                          ! P.B.C. only for x and y direction

          rij_abs2 = rij(1)**2 + rij(2)**2 + rij(3)**2   ! = |rij|^2

!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2 = 1.0d0 / rij_abs2   ! = 1/|rij|^2
          rij_abs_inv1 = sqrt(rij_abs_inv2) ! = 1/|rij|

          potij = - atmchrg(i) * atmchrg(jj) * rij_abs_inv1 ! = qi*(-qj)/|rij|
          pot_cstmnb = pot_cstmnb + potij * 0.5d0   ! 1/2 is needed for energy

          fij(1:3) = potij * rij_abs_inv2 * rij(1:3)
          fcstmnb(1:3,i) = fcstmnb(1:3,i) + fij(1:3)
                                            ! = qi*(-qj)/|rij|^3 * rij

!         --- calculate virial tensor ---
          do n=1,3
             pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                  &                  + fij(1:3)*atmcor(n,i)
          end do

          pot_viri_cstmnb = pot_viri_cstmnb + fij(1)*atmcor(1,i) &
               &          + fij(2)*atmcor(2,i) + fij(3)*atmcor(3,i)

       END DO
    END DO


!---- CALCULATE ZP ANISO- and ISOTROPIC SHORT RANGE FORCE ----

    sqcut_aniso = rcut_ZPaniso(1) * rcut_ZPaniso(1)
    sqcut_iso = rcut_ZPiso(1) * rcut_ZPiso(1)

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO cstmnb1 = 1, ncstmnb - 1
    looplast = ncstmnb - 1

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = cstmnb_indexall(cstmnb1)              ! the first and
       j2 = cstmnb_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting
                                                  ! with i

       i = ncstmnblist(cstmnb1)
       itype = atmindex(i)                  ! atmtype of atom(i)

       IF (cstmnb_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = cstmnb_listall(j)          ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!          --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3) * anint(rij(1:3)*box_inv(1:3))

!         * box_inv rather than /box accellerates the computation

          xyij_abs2 = rij(1)**2 + rij(2)**2
          zij_abs2 = rij(3)**2
          rij_abs2  = xyij_abs2 + zij_abs2 ! |rij|^2

!         --- cut off for anisotropic force ---
          if (rij_abs2 > sqcut_aniso) cycle

          alpha2 = para_alphaZP_aniso(itype,jtype)**2
          alpha2_inv = 1.0d0 / alpha2
          sigma2 = para_radZP(itype,jtype)**2

          rij_abs2_alpha = alpha2*xyij_abs2 + zij_abs2
                                          ! (xij^2+yij^2)*alpha^2 + zij^2
          rij_abs2_a_inv = 1.0d0 / rij_abs2_alpha
                                    ! =1/((xij^2+yij^2)*alpha^2 + zij^2)
          rij_abs2_s_a_inv = sigma2 * rij_abs2_a_inv
                                    ! =sigma^2/((xij^2+yij^2)*alpha^2 + zij^2)
          rij_abs6_s_a_inv = rij_abs2_s_a_inv * rij_abs2_s_a_inv &
               &           * rij_abs2_s_a_inv
                                ! =(sigma^2/((xij^2+yij^2)*alpha^2 + zij^2))^3
          rij_abs12_s_a_inv = rij_abs6_s_a_inv * rij_abs6_s_a_inv
                                ! =(sigma^2/((xij^2+yij^2)*alpha^2 + zij^2))^6

          rij_abs2_alpha_inv = alpha2_inv*xyij_abs2 + zij_abs2
                                          ! (xij^2+yij^2)/alpha^2 + zij^2
          rij_abs2_ainv_inv = 1.0d0 / rij_abs2_alpha_inv
                                    ! =1/((xij^2+yij^2)/alpha^2 + zij^2)
          rij_abs2_s_ainv_inv = sigma2 * rij_abs2_ainv_inv
                                    ! =sigma^2/((xij^2+yij^2)/alpha^2 + zij^2)
          rij_abs6_s_ainv_inv = rij_abs2_s_ainv_inv * rij_abs2_s_ainv_inv &
               &              * rij_abs2_s_ainv_inv
                                ! =(sigma^2/((xij^2+yij^2)/alpha^2 + zij^2))^3

          pot_cstmnb  = pot_cstmnb + 4.0d0*para_welZP(itype,jtype) &
               &                   * (rij_abs12_s_a_inv - rij_abs6_s_ainv_inv)
                                 ! anisotropic potential

          rtmp6 = rij_abs12_s_a_inv * rij_abs2_a_inv  ! f_6 term

          rtmp3 = rij_abs6_s_ainv_inv * rij_abs2_ainv_inv ! f_3 term

          fcoeff = 24.0d0 * para_welZP(itype,jtype)
          ftmp = fcoeff * (2.0d0*alpha2*rtmp6 - alpha2_inv*rtmp3)
                                 ! =24*epsilon*(2*alpha^2*f_6 - 1/alpha^2*f_3)

          fij(1:2)     = ftmp * rij(1:2)

          fij(3) = fcoeff * (2.0d0*rtmp6 - rtmp3) * rij(3)
                                 ! =24*epsilon*(2*f_6 - f_3) * zij

!         --- calculate virial tensor ---
          do n=1,3
             pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                  &                  + fij(1:3)*rij(n)
          end do

          pot_viri_cstmnb = pot_viri_cstmnb &
               &          + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

          fcstmnb(1:3,i)  = fcstmnb(1:3,i)  + fij(1:3)
          fcstmnb(1:3,jj) = fcstmnb(1:3,jj) - fij(1:3)

!         --- cut off for isotropic force ---
          if (rij_abs2 > sqcut_iso) cycle

          rij_abs_inv2 = 1.0d0 / rij_abs2                  ! = 1/|rij|^2
          rij_abs2_s_inv = sigma2 * rij_abs_inv2           ! =sigma^2/|rij|^2
          rij_abs4_s_inv = rij_abs2_s_inv * rij_abs2_s_inv ! =sigma^4/|rij|^4
          rij_abs8_s_inv = rij_abs4_s_inv * rij_abs4_s_inv ! =sigma^8/|rij|^8
          rij_abs10_s_inv = rij_abs8_s_inv * rij_abs2_s_inv
                                         ! =sigma^10/|rij|^10

          potij  = -4.0d0*para_welZP(itype,jtype)*para_cZP_iso(itype,jtype) &
               & * rij_abs10_s_inv  ! =-4*epsilon*C10*sigma^10/|rij|^10

          pot_cstmnb = pot_cstmnb + potij

          ftmp = 10.0d0 * potij * rij_abs_inv2 ! =10*phi_ij/|rij|^2
          fij(1:3) = ftmp * rij(1:3) ! =10*phi_ij/|rij|^2*rij

!         --- calculate virial tensor ---
          do n=1,3
             pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                  &                  + fij(1:3)*rij(n)
          end do

          pot_viri_cstmnb = pot_viri_cstmnb &
               &          + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

          fcstmnb(1:3,i)  = fcstmnb(1:3,i)  + fij(1:3)
          fcstmnb(1:3,jj) = fcstmnb(1:3,jj) - fij(1:3)

       END DO
    END DO

!   --- ADD FCSTMNB TO FORCE ---
    force(1:3,1:natom) = force(1:3,1:natom) + fcstmnb(1:3,1:natom)
    for_viri_cstmnb(1:3,1:natom) = for_viri_cstmnb(1:3,1:natom) &
         &                       + fcstmnb(1:3,1:natom)

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
!      real(8):: fcstmnball(3,natom) ! custom NB force calculated here
!MPI                                       local(in MPI sense) field
    real(8):: potij                      ! temporal potential due to IJ
    real(8):: fcoeff,ftmp                ! temporal force
    real(8):: fij(3)                     ! force to be subtracted due to
                                         ! excluded atoms
    real(8):: fji(3)                     ! force to be subtracted due to
                                         ! excluded atoms = -fij

    real(8):: rij(3)                     ! rj-ri
    real(8):: rij_abs                    ! |Rij|
    real(8):: xyij_abs2                  ! xij**2 + yij**2
    real(8):: zij_abs2                   ! zij**2
    real(8):: rij_abs2                   ! |Rij|**2
    real(8):: rij_abs_inv2               ! |Rij|**-2
    real(8):: rij_abs_inv1               ! |Rij|**-1

    real(8):: alpha2                     ! alpha**2
    real(8):: alpha2_inv                 ! 1/alpha**2
    real(8):: sigma2                     ! sigma**2

    real(8):: rij_abs2_alpha
    real(8):: rij_abs2_a_inv
    real(8):: rij_abs2_s_a_inv,rij_abs6_s_a_inv,rij_abs12_s_a_inv
    real(8):: rij_abs2_alpha_inv
    real(8):: rij_abs2_ainv_inv
    real(8):: rij_abs2_s_ainv_inv,rij_abs6_s_ainv_inv
    real(8):: rij_abs2_s_inv
    real(8):: rij_abs4_s_inv,rij_abs8_s_inv,rij_abs10_s_inv

    real(8):: rtmp6, rtmp3

    real(8):: box(3)                      ! BOX size
    real(8):: box_inv(3)                  ! inverse of BOX size

    real(8):: sqcut                      ! rcut * rcut
    real(8):: sqcut_aniso,sqcut_iso      ! rcut * rcut

    integer:: itype, jtype                ! VDW atom types
    integer:: i,j,k,n                     ! do loop index
    integer:: j1, j2                      ! limits of do loop
    integer:: jj                          ! atom index of the interacting atom
    integer:: imole

    real(8):: ZPint_p1_ave, ZPint_p2_ave 
                            ! instantaneous averaged position of mirror surfaces

    integer:: cstmnb1

    real(8):: nbodycoeff          ! n-body coefficient for virial term of heatf

!     +     +     +     +     +     +     +

!---- some initialization

    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    pot_virit_cstmnb(1:3,1:3) = 0.0d0
    pot_viri_cstmnb = 0.0d0

    box(1) = xcel
    box(2) = ycel
    box(3) = zcel
    box_inv(1:3) = 1.0d0/box(1:3)

    nbodycoeff = 0.5d0        ! = 1/2

!---- CALCULATE ZP CONDUCTION FORCE (mirror force) ----

    sqcut = rcut_ZPcond(1) * rcut_ZPcond(1)

!---- update the position of mirror surface #1

    if (ifrenew_msurf(1)) then

       ZPint_p1_ave = 0.0d0

       imole = npoly + nwater
       do i = 1, matomtyp_ZPint_p1 - 1
          imole = imole + nmatomtyp(i)
       end do

       do i = matom_ZPint_p1_ini, matom_ZPint_p1_end
          j1 = molept_index(imole+i)         ! first monatomic mole.
          j2 = molept_index(imole+i+1) - 1

          do j= j1, j2
             k = molept_list(j)

             ZPint_p1_ave = ZPint_p1_ave + atmcor(3,k)

          end do

       end do

       para_ZPint_p1(1) = ZPint_p1_ave &
            &           / dble(matom_ZPint_p1_end - matom_ZPint_p1_ini + 1)

    end if

!---- surface #1
! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO cstmnb1 = 1, nZPcond
    looplast = nZPcond

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = ZPcond1_indexall(cstmnb1)              ! the first and
       j2 = ZPcond1_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting

       i = nZPcondlist(cstmnb1)
       itype = atmindex(i)               ! atmtype of atom(i)

       IF (ZPcond1_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = ZPcond1_listall(j)         ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:2) = atmcor(1:2,i) - atmcor(1:2,jj)
          rij(3) = atmcor(3,i) - (2.0d0*para_ZPint_p1(1) - atmcor(3,jj))
                                          ! = zi - zj' (mirror image)

!         --- periodic boundary ---

          rij(1:2) = rij(1:2) - box(1:2) * anint(rij(1:2)*box_inv(1:2))
                                          ! P.B.C. only for x and y direction

          rij_abs2 = rij(1)**2 + rij(2)**2 + rij(3)**2   ! = |rij|^2

!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2 = 1.0d0 / rij_abs2   ! = 1/|rij|^2
          rij_abs_inv1 = sqrt(rij_abs_inv2) ! = 1/|rij|

          potij = - atmchrg(i) * atmchrg(jj) * rij_abs_inv1 ! = qi*(-qj)/|rij|
          pot_cstmnb = pot_cstmnb + potij * 0.5d0   ! 1/2 is needed for energy

#if defined(_THIS_IS_NOT_CLEAR)
          pot_cstmnb_atm(i) = pot_cstmnb_atm(i) + potij * 0.5d0
#endif

          fij(1:3) = potij * rij_abs_inv2 * rij(1:3)
          fcstmnb(1:3,i) = fcstmnb(1:3,i) + fij(1:3)
                                            ! = qi*(-qj)/|rij|^3 * rij

!         --- calculate virial tensor ---

#if defined(_THIS_IS_NOT_CLEAR)
!      - for heat flux calculation

          fji(1:3) = -fij(1:3)

          call calvirihf(i,jj,fij,fji,   &
               &         nbodycoeff,   &
               &         box,box_inv,   &
               &         ifhfvol,   &
               &         nhfregion,hfzpos1,hfzpos2,   &
               &         hftyp_atm,   &
               &         molecom,   &
               &         viricstmnbt_atm)
#endif

!         - for pressure calculation
          do n=1,3
             pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                  &                  + fij(1:3)*atmcor(n,i)
          end do

          pot_viri_cstmnb = pot_viri_cstmnb + fij(1)*atmcor(1,i) &
               &          + fij(2)*atmcor(2,i) + fij(3)*atmcor(3,i)

       END DO
    END DO

!---- update the position of mirror surface #2

    if (ifrenew_msurf(1)) then

       ZPint_p2_ave = 0.0d0

       imole = npoly + nwater
       do i = 1, matomtyp_ZPint_p2 - 1
          imole = imole + nmatomtyp(i)
       end do

       do i = matom_ZPint_p2_ini, matom_ZPint_p2_end
          j1 = molept_index(imole+i)         ! first monatomic mole.
          j2 = molept_index(imole+i+1) - 1

          do j= j1, j2
             k = molept_list(j)

             ZPint_p2_ave = ZPint_p2_ave + atmcor(3,k)

          end do

       end do

       para_ZPint_p2(1) = ZPint_p2_ave &
            &           / dble(matom_ZPint_p2_end - matom_ZPint_p2_ini + 1)

    end if

!---- surface #2
! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO cstmnb1 = 1, nZPcond
    looplast = nZPcond

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = ZPcond2_indexall(cstmnb1)              ! the first and
       j2 = ZPcond2_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting

       i = nZPcondlist(cstmnb1)
       itype = atmindex(i)               ! atmtype of atom(i)

       IF (ZPcond2_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = ZPcond2_listall(j)         ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:2) = atmcor(1:2,i) - atmcor(1:2,jj)
          rij(3) = atmcor(3,i) - (2.0d0*para_ZPint_p2(1) - atmcor(3,jj))
                                          ! = zi - zj' (mirror image)

!         --- periodic boundary ---

          rij(1:2) = rij(1:2) - box(1:2) * anint(rij(1:2)*box_inv(1:2))
                                          ! P.B.C. only for x and y direction

          rij_abs2 = rij(1)**2 + rij(2)**2 + rij(3)**2   ! = |rij|^2

!         --- cut off ---
          if (rij_abs2 > sqcut) cycle

          rij_abs_inv2 = 1.0d0 / rij_abs2   ! = 1/|rij|^2
          rij_abs_inv1 = sqrt(rij_abs_inv2) ! = 1/|rij|

          potij = - atmchrg(i) * atmchrg(jj) * rij_abs_inv1 ! = qi*(-qj)/|rij|
          pot_cstmnb = pot_cstmnb + potij * 0.5d0   ! 1/2 is needed for energy

#if defined(_THIS_IS_NOT_CLEAR)
          pot_cstmnb_atm(i) = pot_cstmnb_atm(i) + potij * 0.5d0
#endif

          fij(1:3) = potij * rij_abs_inv2 * rij(1:3)
          fcstmnb(1:3,i) = fcstmnb(1:3,i) + fij(1:3)
                                            ! = qi*(-qj)/|rij|^3 * rij

!         --- calculate virial tensor ---

#if defined(_THIS_IS_NOT_CLEAR)
!      - for heat flux calculation

          fji(1:3) = -fij(1:3)

          call calvirihf(i,jj,fij,fji,   &
               &         nbodycoeff,   &
               &         box,box_inv,   &
               &         ifhfvol,   &
               &         nhfregion,hfzpos1,hfzpos2,   &
               &         hftyp_atm,   &
               &         molecom,   &
               &         viricstmnbt_atm)
#endif

!         - for pressure calculation
          do n=1,3
             pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                  &                  + fij(1:3)*atmcor(n,i)
          end do

          pot_viri_cstmnb = pot_viri_cstmnb + fij(1)*atmcor(1,i) &
               &          + fij(2)*atmcor(2,i) + fij(3)*atmcor(3,i)

       END DO
    END DO


!---- CALCULATE ZP ANISO- and ISOTROPIC SHORT RANGE FORCE ----

    sqcut_aniso = rcut_ZPaniso(1) * rcut_ZPaniso(1)
    sqcut_iso = rcut_ZPiso(1) * rcut_ZPiso(1)

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO cstmnb1 = 1, ncstmnb - 1
    looplast = ncstmnb - 1

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = cstmnb_indexall(cstmnb1)              ! the first and
       j2 = cstmnb_indexall(cstmnb1+loopstep) - 1 ! the last atoms interacting
                                                  ! with i

       i = ncstmnblist(cstmnb1)
       itype = atmindex(i)                  ! atmtype of atom(i)

       IF (cstmnb_listall(j1) == 0) CYCLE

       DO j = j1, j2

          jj = cstmnb_listall(j)          ! atom index of interacting atom
          jtype = atmindex(jj)            ! atmtype of atom(jj)

          rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!          --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3) * anint(rij(1:3)*box_inv(1:3))

!         * box_inv rather than /box accellerates the computation

          xyij_abs2 = rij(1)**2 + rij(2)**2
          zij_abs2 = rij(3)**2
          rij_abs2  = xyij_abs2 + zij_abs2 ! |rij|^2

!         --- cut off for anisotropic force ---
          if (rij_abs2 > sqcut_aniso) cycle

          alpha2 = para_alphaZP_aniso(itype,jtype)**2
          alpha2_inv = 1.0d0 / alpha2
          sigma2 = para_radZP(itype,jtype)**2

          rij_abs2_alpha = alpha2*xyij_abs2 + zij_abs2
                                          ! (xij^2+yij^2)*alpha^2 + zij^2
          rij_abs2_a_inv = 1.0d0 / rij_abs2_alpha
                                    ! =1/((xij^2+yij^2)*alpha^2 + zij^2)
          rij_abs2_s_a_inv = sigma2 * rij_abs2_a_inv
                                    ! =sigma^2/((xij^2+yij^2)*alpha^2 + zij^2)
          rij_abs6_s_a_inv = rij_abs2_s_a_inv * rij_abs2_s_a_inv &
               &           * rij_abs2_s_a_inv
                                ! =(sigma^2/((xij^2+yij^2)*alpha^2 + zij^2))^3
          rij_abs12_s_a_inv = rij_abs6_s_a_inv * rij_abs6_s_a_inv
                                ! =(sigma^2/((xij^2+yij^2)*alpha^2 + zij^2))^6

          rij_abs2_alpha_inv = alpha2_inv*xyij_abs2 + zij_abs2
                                          ! (xij^2+yij^2)/alpha^2 + zij^2
          rij_abs2_ainv_inv = 1.0d0 / rij_abs2_alpha_inv
                                    ! =1/((xij^2+yij^2)/alpha^2 + zij^2)
          rij_abs2_s_ainv_inv = sigma2 * rij_abs2_ainv_inv
                                    ! =sigma^2/((xij^2+yij^2)/alpha^2 + zij^2)
          rij_abs6_s_ainv_inv = rij_abs2_s_ainv_inv * rij_abs2_s_ainv_inv &
               &              * rij_abs2_s_ainv_inv
                                ! =(sigma^2/((xij^2+yij^2)/alpha^2 + zij^2))^3

          potij  = 4.0d0*para_welZP(itype,jtype) &
               & * (rij_abs12_s_a_inv - rij_abs6_s_ainv_inv)
          pot_cstmnb = pot_cstmnb + potij
                                 ! anisotropic potential

          pot_cstmnb_atm(i) = pot_cstmnb_atm(i) + 0.5d0*potij
          pot_cstmnb_atm(jj) = pot_cstmnb_atm(jj) + 0.5d0*potij

          rtmp6 = rij_abs12_s_a_inv * rij_abs2_a_inv  ! f_6 term

          rtmp3 = rij_abs6_s_ainv_inv * rij_abs2_ainv_inv ! f_3 term

          fcoeff = 24.0d0 * para_welZP(itype,jtype)
          ftmp = fcoeff * (2.0d0*alpha2*rtmp6 - alpha2_inv*rtmp3)
                                 ! =24*epsilon*(2*alpha^2*f_6 - 1/alpha^2*f_3)

          fij(1:2)     = ftmp * rij(1:2)

          fij(3) = fcoeff * (2.0d0*rtmp6 - rtmp3) * rij(3)
                                 ! =24*epsilon*(2*f_6 - f_3) * zij

!         --- calculate virial tensor ---

!         - for heat flux calculation

          fji(1:3) = -fij(1:3)

          call calvirihf(i,jj,fij,fji,   &
               &         nbodycoeff,   &
               &         box,box_inv,   &
               &         ifhfvol,   &
               &         nhfregion,hfzpos1,hfzpos2,   &
               &         hftyp_atm,   &
               &         molecom,   &
               &         viricstmnbt_atm)

!         - for pressure calculation
          do n=1,3
             pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                  &                  + fij(1:3)*rij(n)
          end do

          pot_viri_cstmnb = pot_viri_cstmnb &
               &          + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

          fcstmnb(1:3,i)  = fcstmnb(1:3,i)  + fij(1:3)
          fcstmnb(1:3,jj) = fcstmnb(1:3,jj) - fij(1:3)

!         --- cut off for isotropic force ---
          if (rij_abs2 > sqcut_iso) cycle

          rij_abs_inv2 = 1.0d0 / rij_abs2                  ! = 1/|rij|^2
          rij_abs2_s_inv = sigma2 * rij_abs_inv2           ! =sigma^2/|rij|^2
          rij_abs4_s_inv = rij_abs2_s_inv * rij_abs2_s_inv ! =sigma^4/|rij|^4
          rij_abs8_s_inv = rij_abs4_s_inv * rij_abs4_s_inv ! =sigma^8/|rij|^8
          rij_abs10_s_inv = rij_abs8_s_inv * rij_abs2_s_inv
                                         ! =sigma^10/|rij|^10

          potij  = -4.0d0*para_welZP(itype,jtype)*para_cZP_iso(itype,jtype) &
               & * rij_abs10_s_inv  ! =-4*epsilon*C10*sigma^10/|rij|^10

          pot_cstmnb = pot_cstmnb + potij

          pot_cstmnb_atm(i) = pot_cstmnb_atm(i) + 0.5d0*potij
          pot_cstmnb_atm(jj) = pot_cstmnb_atm(jj) + 0.5d0*potij

          ftmp = 10.0d0 * potij * rij_abs_inv2 ! =10*phi_ij/|rij|^2
          fij(1:3) = ftmp * rij(1:3) ! =10*phi_ij/|rij|^2*rij

!         --- calculate virial tensor ---

!         - for heat flux calculation

          fji(1:3) = -fij(1:3)

          call calvirihf(i,jj,fij,fji,   &
               &         nbodycoeff,   &
               &         box,box_inv,   &
               &         ifhfvol,   &
               &         nhfregion,hfzpos1,hfzpos2,   &
               &         hftyp_atm,   &
               &         molecom,   &
               &         viricstmnbt_atm)

!         - for pressure calculation
          do n=1,3
             pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                  &                  + fij(1:3)*rij(n)
          end do

          pot_viri_cstmnb = pot_viri_cstmnb &
               &          + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

          fcstmnb(1:3,i)  = fcstmnb(1:3,i)  + fij(1:3)
          fcstmnb(1:3,jj) = fcstmnb(1:3,jj) - fij(1:3)

       END DO
    END DO

!   --- ADD FCSTMNB TO FORCE ---
    force(1:3,1:natom) = force(1:3,1:natom) + fcstmnb(1:3,1:natom)
    for_viri_cstmnb(1:3,1:natom) = for_viri_cstmnb(1:3,1:natom) &
         &                       + fcstmnb(1:3,1:natom)

!     +     +     +     +     +     +     +

  end subroutine calcstmnbp_hf
#endif

end module cstmnb
