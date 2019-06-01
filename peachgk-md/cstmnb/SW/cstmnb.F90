!**********************************
!*  cstmnb.f90 Ver.1.2            *
!*      for peachgk_md.f          *
!*            by G.Kikugawa       *
!*      Coded by T.Nakano         *
!**********************************
! Time-stamp: <2015-03-12 13:52:23 gota>

! This custom NB potential is designed for Stillinger-Weber potential
!   to model a silicon solid material.
! See Stillinger, F.H. and Weber, T.A., Phys. Rev. B, 31, pp. 5262-5271, (1985).

module cstmnb

  use md_global
  use mpi_global

  implicit none

  private

!!! define global variables here
  integer,public,save:: ncstmnb	              ! number of all cstm NB atoms
  integer,public,save:: ncstmnblist(maxnatom) ! list of cstm NB atoms
  logical,public,save:: ifcstmnb(maxnatom)    ! if cstm NB atom or not

  integer,public,save:: cstmnb_indexall2(maxnatom+maxnproc)
                                    ! index for location of the list of atom
  integer,public,save,allocatable:: cstmnb_indexall3(:)
                                    ! index for location of the list of atom
  integer,public,save,allocatable:: cstmnb_listall2(:)
                                    ! list of neighboring atoms
  integer,public,save,allocatable:: cstmnb_listall3(:)
                                    ! list of neighboring atoms
  integer,public,save:: maxcstmnblist2  !
  integer,public,save:: maxcstmnblist3  !

  integer,public,save:: nlistcstmnball2     ! number of custom NB list
  integer,public,save:: nlistcstmnball3     ! number of custom NB list

  integer,public,save:: ncstmnbtyp  ! number of cstmnb type
  character(2),public,save:: para_cstmnbtyp(maxncstmnbtyp)
                                    ! atom of custom NB type
  integer,public,save:: atmindex_ncstmnb(maxnatom)
                                    ! atomindex of custom NB atom

  logical,public,save:: ifcalSW2(0:maxncstmnbtyp,0:maxncstmnbtyp)
  logical,public,save:: ifcalSW3(0:maxncstmnbtyp,0:maxncstmnbtyp,0:maxncstmnbtyp)
                                    ! if this pair is calculated in cstmnb

  real(8),public,save:: rcut_SW(maxnatmtyp)  ! cutoff length for SW potential
  real(8),public,save:: rcut_bookSW  
                          ! cutoff length for book-keeping of SW potential
  integer,public,save:: nstep_bookSW  ! interval for bookkeeping

  real(8),public,save:: A_SW(maxnatmtyp) ! A parameter for SW potential
  real(8),public,save:: B_SW(maxnatmtyp) ! B parameter for SW potential
  integer,public,save:: p_SW(maxnatmtyp) ! p-th power root for SW potential
  integer,public,save:: q_SW(maxnatmtyp) ! q-th power root for SW potential

  real(8),public,save:: gamma_SW(maxnatmtyp) ! gamma parameter for SW potential
  real(8),public,save:: lambda_SW(maxnatmtyp) 
                                             ! lambda parameter for SW potential
  real(8),public,save:: theta0_SW(maxnatmtyp)   ! theta0 for SW potential

  real(8),public,save:: epsilon_SW(maxncstmnbtyp)   ! parameter for SW potential
  real(8),public,save:: sigma2_SW(maxncstmnbtyp)
                                 ! non-deimnsinalized parameter for SW potential
  real(8),public,save:: sigma3_SW(maxncstmnbtyp)
                                 ! non-deimnsinalized parameter for SW potential

  integer,public,save:: cstmnbtypeindex2(maxncstmnbtyp,maxncstmnbtyp)
  integer,public,save:: cstmnbtypeindex3(maxncstmnbtyp,maxncstmnbtyp,maxncstmnbtyp)

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
    real(8):: fcstmnb(3,natom)          ! custom NB force calculated here
    real(8):: atmcor_sig(3,natom)       ! non-dimenstioinal atom coodinate

    real(8):: A_, B_                    ! parameters of SW 2-body potential
    integer:: p_, q_                    !
    integer:: np_, nq_                  !
    integer:: np_1, nq_1                !
    real(8):: gamma, lambda, costheta0  ! parameters of SW 3-body potential

    real(8):: pi                  ! = 3.14
    real(8):: potij               ! temporal potential due to IJ
    real(8):: ftmp                ! temporal force
    real(8):: fi(3), fj(3), fk(3) ! force to be subtracted due to
                                     ! excluded atoms

    real(8):: rij(3),   rik(3),   rjk(3)                ! rj-ri
    real(8):: rij_abs,  rik_abs,  rjk_abs               ! |rij|
    real(8):: rij_abs2, rik_abs2, rjk_abs2              ! |rij|**2
    real(8):: rij_abs_inv1, rik_abs_inv1, rjk_abs_inv1  ! 1/|rij|
    real(8):: rij_abs_inv2, rik_abs_inv2, rjk_abs_inv2  ! 1/|rij|**2

    real(8):: box(3,ncstmnbtyp)      ! BOX size
    real(8):: box_inv(3,ncstmnbtyp)  ! inverse of BOX size

    real(8):: sqcut               ! rcut*rcut

    integer:: icstmnbtype, jcstmnbtype, kcstmnbtype ! SW atom types
    integer:: itype, jtype, ktype ! VDW atom types
    integer:: i,j,k,n             ! do loop index
    integer:: j1, j2              ! limits of do loop
    integer:: k1, k2              ! limits of do loop
    integer:: jj, kk              ! atom index of the interacting atom

    integer:: iclose          ! if calculate force and potential or not

    integer:: cstmnb1, cstmnb2, cstmnb3
    integer:: icmtype, i2btype, i3btype ! cstmnb type of common, 2-body & 3-body

    real(8):: tmp

    real(8):: a
    real(8):: rij_a, rik_a, rjk_a
    real(8):: vcuij,  vcuik,  vcujk        ! vcuij  = 1 / rij_a
    real(8):: evcuij, evcuik, evcujk       ! evcuij = exp(vcuij)
    real(8):: rhoij, rhoik, rhojk

    real(8):: ri_inv,rj_inv,rk_inv
    real(8):: vcui,vcuj,vcuk
    real(8):: ei,ej,ek

    real(8):: hjik,hijk,hikj
    real(8):: hij, hik, hjk

    real(8):: cosjik,  cosijk,  cosikj
    real(8):: cosjik_0,cosijk_0,cosikj_0

    real(8):: eiNjik(3), eiNkij(3)
    real(8):: ejNijk(3), ejNkji(3)
    real(8):: ekNikj(3), ekNjki(3)

    real(8):: hRij(3), hRik(3), hRjk(3)

    real(8):: DELTA = 1.0d-32   ! extremely small value

! FUNCTIONS:

!     +     +     +     +     +     +     +


!---- some initialization

    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

    do cstmnb1 = 1, ncstmnb
       i = ncstmnblist(cstmnb1)
       itype = atmindex_ncstmnb(i)
       atmcor_sig(1:3,i) = atmcor(1:3,i) / sigma2_SW(itype)

       box(1,itype) = xcel / sigma2_SW(itype)
       box(2,itype) = ycel / sigma2_SW(itype)
       box(3,itype) = zcel / sigma2_SW(itype)
       box_inv(1:3,itype) = 1.0d0/box(1:3,itype)
    end do

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, ncstmnb - 1

    looplast = ncstmnb - 1

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = cstmnb_indexall2(cstmnb1)              ! the first and
       j2 = cstmnb_indexall2(cstmnb1+loopstep) - 1 ! the last atoms interacting
                                                  ! with i
       i = ncstmnblist(cstmnb1)
       itype = atmindex_ncstmnb(i)       ! atom type of atom(i)

       IF (cstmnb_listall2(j1) == 0) CYCLE

       DO cstmnb2 = j1, j2

          j = cstmnb_listall2(cstmnb2)    ! atom index of interacting atom
          jtype = atmindex_ncstmnb(j)

          IF (.not. ifcalSW2(itype,jtype)) CYCLE

          i2btype = cstmnbtypeindex2(itype,jtype)

          rij(1:3) = atmcor_sig(1:3,i) - atmcor_sig(1:3,j)

!         --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3,i2btype) &
               &   * anint(rij(1:3)*box_inv(1:3,i2btype))

!         * box_inv rather than /box accellerates the computation

          rij_abs2 = rij(1)**2 + rij(2)**2 + rij(3)**2 ! |rij|^2
          rij_abs  = sqrt(rij_abs2)   ! |rij|

!         --- cut off ---
          sqcut = rcut_SW(i2btype)**2

          if (rij_abs2 < sqcut) then

             A_  = A_SW(i2btype)
             B_  = B_SW(i2btype)
             p_  = p_SW(i2btype)
             q_  = q_SW(i2btype)
             np_ = - p_
             nq_ = - q_
             np_1 = np_ - 1
             nq_1 = nq_ - 1

             rij_abs_inv1  = 1.0d0 / rij_abs  ! 1/|rij|

             vcuij  = 1.0d0 / (rij_abs - rcut_SW(i2btype))
!             evcuij = 0.0d0
!             if (vcuij > -30.0d0) evcuij = exp(vcuij)
             if (vcuij > -30.0d0) then
                evcuij = exp(vcuij)
                tmp = B_ * rij_abs_inv1**p_ - rij_abs_inv1**q_
                potij = A_ * tmp * evcuij

                !---------- 2-body force ----------

                ftmp = potij*(((B_*p_*(rij_abs**np_1)-q_*(rij_abs**nq_1)) &
                     &   / (B_*(rij_abs**np_)-(rij_abs**nq_))) + vcuij*vcuij) &
                     &   * rij_abs_inv1

!               --- add both terms ---
                pot_cstmnb = pot_cstmnb + potij

                fi(1:3) = ftmp * rij(1:3)
                fcstmnb(1:3,i) = fcstmnb(1:3,i) + fi(1:3)
                fcstmnb(1:3,j) = fcstmnb(1:3,j) - fi(1:3)
             end if

          end if

!         --- 3-body potential terms ---

          k1 = cstmnb_indexall3(cstmnb2)              ! the first and
          k2 = cstmnb_indexall3(cstmnb2+1) - 1        ! the last atoms interacting

          IF (cstmnb_listall3(k1) == 0) CYCLE

          DO kk = k1, k2

             k = cstmnb_listall3(kk)     ! atom index of interacting atom
             ktype = atmindex_ncstmnb(k) ! atmtype of atom(k)

             if (.not. ifcalSW3(itype,jtype,ktype)) cycle

             rij_a = 0.0d0
             vcuij = 0.0d0
             rhoij = 0.0d0

             rik_a = 0.0d0
             vcuik = 0.0d0
             rhoik = 0.0d0

             rjk_a = 0.0d0
             vcujk = 0.0d0
             rhojk = 0.0d0

             i3btype = cstmnbtypeindex3(itype,jtype,ktype)

             rik(1:3) = atmcor_sig(1:3,i) - atmcor_sig(1:3,k)
             rik(1:3) = rik(1:3) - box(1:3,i3btype) &
                  &   * anint(rik(1:3)*box_inv(1:3,i3btype))
             rik_abs2 = rik(1)**2 + rik(2)**2 + rik(3)**2
             rik_abs = sqrt(rik_abs2)

             rjk(1:3) = atmcor_sig(1:3,j) - atmcor_sig(1:3,k)
             rjk(1:3) = rjk(1:3) - box(1:3,i3btype) &
                  &   * anint(rjk(1:3)*box_inv(1:3,i3btype))
             rjk_abs2 = rjk(1)**2 + rjk(2)**2 + rjk(3)**2
             rjk_abs = sqrt(rjk_abs2)

             gamma  =  gamma_SW(i3btype)
             lambda = lambda_SW(i3btype)
             costheta0 = cos(theta0_SW(i3btype))
             a = rcut_SW(i3btype)                 ! parameter of 2-body pot.

             iclose = 0

             if (rij_abs2 < sqcut) then
                rij_a = rij_abs - a
                rhoij = rij_abs * rij_a * rij_a
                tmp = gamma / rij_a
                if (tmp > -30.0d0) then
                   vcuij = exp(tmp)
                   iclose = iclose + 1
                end if
             end if

             if (rik_abs2 < sqcut) then
                rik_a = rik_abs - a
                rhoik = rik_abs * rik_a * rik_a
                tmp = gamma / rik_a
                if (tmp > -30.0d0) then
                   vcuik = exp(tmp)
                   iclose = iclose + 1
                end if
             end if

             if (rjk_abs2 < sqcut) then
                rjk_a = rjk_abs - a
                rhojk = rjk_abs * rjk_a * rjk_a
                tmp = gamma / rjk_a
                if (tmp > -30.0d0) then
                   vcujk = exp(tmp)
                   iclose = iclose + 1
                end if
             end if

             if (iclose < 2) cycle   !  Skip the following calculations:
!            Forces and potentials of i,j,k turn out to be zero
!            if not less than two of rxx_abs are in the cutoff distance

             hjik = 0.0d0
             hijk = 0.0d0
             hikj = 0.0d0
             hij  = 0.0d0
             hik  = 0.0d0
             hjk  = 0.0d0
             ei   = 0.0d0
             ej   = 0.0d0
             ek   = 0.0d0
             ri_inv = 0.0d0
             rj_inv = 0.0d0
             rk_inv = 0.0d0
             cosjik = 0.0d0
             cosijk = 0.0d0
             cosikj = 0.0d0
             cosjik_0 = 0.0d0
             cosijk_0 = 0.0d0
             cosikj_0 = 0.0d0

             eiNjik(1:3) = 0.0d0
             eiNkij(1:3) = 0.0d0
             ejNijk(1:3) = 0.0d0
             ejNkji(1:3) = 0.0d0
             ekNikj(1:3) = 0.0d0
             ekNjki(1:3) = 0.0d0

             vcui = vcuij * vcuik
             if (vcui > DELTA) then   ! if (vcui /= 0)
                ri_inv = 1.0d0 / (rij_abs * rik_abs)
                cosjik   = (rij(1)*rik(1) + rij(2)*rik(2) + rij(3)*rik(3)) &
                     &   * ri_inv
                cosjik_0 = cosjik - costheta0
                tmp  = lambda * vcui * cosjik_0
                ei   = 2.0d0 * tmp * ri_inv
                hjik = tmp * cosjik_0
             end if

             vcuj = vcuij * vcujk
             if (vcuj > DELTA) then   ! if (vcuj /= 0)
                rj_inv = 1.0d0 / (rij_abs * rjk_abs)
                cosijk   = -(rij(1)*rjk(1) + rij(2)*rjk(2) + rij(3)*rjk(3)) * rj_inv
                cosijk_0 = cosijk - costheta0
                tmp  = lambda * vcuj * cosijk_0
                ej   = 2.0d0 * tmp * rj_inv
                hijk = tmp * cosijk_0
             end if

             vcuk = vcuik * vcujk
             if (vcuk > DELTA) then   ! if (vcuk /= 0)
                rk_inv = 1.0d0 / (rik_abs * rjk_abs)
                cosikj   = (rjk(1)*rik(1) + rjk(2)*rik(2) + rjk(3)*rik(3)) &
                     &   * rk_inv
                cosikj_0 = cosikj - costheta0
                tmp  = lambda * vcuk * cosikj_0
                ek   = 2.0d0 * tmp * rk_inv
                hikj = tmp * cosikj_0
             end if

!----------- 3-body potential ----------

             pot_cstmnb = pot_cstmnb + (hjik + hijk + hikj)

!----------- 3-body force ----------

             eiNjik(1:3) = ei*( rik(1:3) - rik_abs/rij_abs*cosjik*rij(1:3))
             eiNkij(1:3) = ei*( rij(1:3) - rij_abs/rik_abs*cosjik*rik(1:3))
             ejNijk(1:3) = ej*( rjk(1:3) + rjk_abs/rij_abs*cosijk*rij(1:3))
             ejNkji(1:3) = ej*(-rij(1:3) - rij_abs/rjk_abs*cosijk*rjk(1:3))
             ekNikj(1:3) = ek*(-rjk(1:3) + rjk_abs/rik_abs*cosikj*rik(1:3))
             ekNjki(1:3) = ek*(-rik(1:3) + rik_abs/rjk_abs*cosikj*rjk(1:3))

             if (rhoij > DELTA) hij = gamma * (hijk + hjik)/rhoij
             if (rhoik > DELTA) hik = gamma * (hjik + hikj)/rhoik
             if (rhojk > DELTA) hjk = gamma * (hikj + hijk)/rhojk

             hRij(1:3) = hij * rij(1:3)
             hRik(1:3) = hik * rik(1:3)
             hRjk(1:3) = hjk * rjk(1:3)

             fi(1:3) =   hRij(1:3) + hRik(1:3) &
                  &  - eiNjik(1:3) - eiNkij(1:3) + ejNijk(1:3) + ekNikj(1:3)
             fj(1:3) = - hRij(1:3) + hRjk(1:3) &
                  &  - ejNijk(1:3) - ejNkji(1:3) + eiNjik(1:3) + ekNjki(1:3)
             fk(1:3) = - hRjk(1:3) - hRik(1:3) &
                  &  - ekNikj(1:3) - ekNjki(1:3) + eiNkij(1:3) + ejNkji(1:3)

             fcstmnb(1:3,i) = fcstmnb(1:3,i) + fi(1:3)
             fcstmnb(1:3,j) = fcstmnb(1:3,j) + fj(1:3)
             fcstmnb(1:3,k) = fcstmnb(1:3,k) + fk(1:3)

          END DO
       END DO
    END DO

!   --- ADD FCSTMNB TO FORCE ---
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
    real(8):: fcstmnb(3,natom)          ! custom NB force calculated here
    real(8):: atmcor_sig(3,natom)        ! non-dimenstioinal atom coodinate

    real(8):: A_, B_                    ! parameters of SW 2-body potential
    integer:: p_, q_                    ! 
    integer:: np_, nq_                  ! 
    integer:: np_1, nq_1                ! 

    real(8):: gamma, lambda, costheta0 ! parameters of SW 3-body potential

    real(8):: pi                  ! = 3.14
    real(8):: potij               ! temporal potential due to IJ
    real(8):: ftmp                ! temporal force
    real(8):: fi(3), fj(3), fk(3) ! force to be subtracted due to
                                     ! excluded atoms

    real(8):: rij(3),   rik(3),   rjk(3)                ! rj-ri
    real(8):: rij_abs,  rik_abs,  rjk_abs               ! |rij|
    real(8):: rij_abs2, rik_abs2, rjk_abs2              ! |rij|**2
    real(8):: rij_abs_inv1, rik_abs_inv1, rjk_abs_inv1  ! 1/|rij|
    real(8):: rij_abs_inv2, rik_abs_inv2, rjk_abs_inv2  ! 1/|rij|**2

    real(8):: box(3,ncstmnbtyp)      ! BOX size
    real(8):: box_inv(3,ncstmnbtyp)  ! inverse of BOX size

    real(8):: sqcut               ! rcut*rcut

    integer:: icstmnbtype, jcstmnbtype, kcstmnbtype ! SW atom types
    integer:: itype, jtype, ktype ! VDW atom types
    integer:: i,j,k,n             ! do loop index
    integer:: j1, j2              ! limits of do loop
    integer:: k1, k2              ! limits of do loop
    integer:: jj, kk              ! atom index of the interacting atom

    integer:: iclose              ! judge calculation or not for force and potential

    integer:: cstmnb1, cstmnb2, cstmnb3
    integer:: icmtype, i2btype, i3btype ! cstmnb type of common, 2-body & 3-body

    real(8):: tmp

    real(8):: a
    real(8):: rij_a, rik_a, rjk_a
    real(8):: vcuij,  vcuik,  vcujk        ! vcuij  = 1 / rij_a  
    real(8):: evcuij, evcuik, evcujk       ! evcuij = exp(vcuij)  
    real(8):: rhoij, rhoik, rhojk

    real(8):: ri_inv,rj_inv,rk_inv
    real(8):: vcui,vcuj,vcuk
    real(8):: ei,ej,ek

    real(8):: hjik,hijk,hikj
    real(8):: hij, hik, hjk

    real(8):: cosjik,  cosijk,  cosikj
    real(8):: cosjik_0,cosijk_0,cosikj_0

    real(8):: eiNjik(3), eiNkij(3)
    real(8):: ejNijk(3), ejNkji(3)
    real(8):: ekNikj(3), ekNjki(3)

    real(8):: hRij(3), hRik(3), hRjk(3)

    real(8):: DELTA = 1.0d-32   ! extremely small value

! FUNCTIONS:

!     +     +     +     +     +     +     +

!---- some initialization

    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

!   --- CALCULATE NONBONDED FORCE ---

    pot_virit_cstmnb(1:3,1:3) = 0.0d0

    pot_viri_cstmnb = 0.0d0

    do cstmnb1 = 1, ncstmnb
       i = ncstmnblist(cstmnb1)
       itype = atmindex_ncstmnb(i)
       atmcor_sig(1:3,i) = atmcor(1:3,i) / sigma2_SW(itype)

       box(1,itype) = xcel / sigma2_SW(itype)
       box(2,itype) = ycel / sigma2_SW(itype)
       box(3,itype) = zcel / sigma2_SW(itype)
       box_inv(1:3,itype) = 1.0d0/box(1:3,itype)
    end do

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, ncstmnb - 1

    looplast = ncstmnb - 1

    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = cstmnb_indexall2(cstmnb1)              ! the first and
       j2 = cstmnb_indexall2(cstmnb1+loopstep) - 1 ! the last atoms interacting
                                                  ! with i
       i = ncstmnblist(cstmnb1)
       itype = atmindex_ncstmnb(i)       ! atom type of atom(i)

       IF (cstmnb_listall2(j1) == 0) CYCLE

       DO cstmnb2 = j1, j2

          j = cstmnb_listall2(cstmnb2)    ! atom index of interacting atom
          jtype = atmindex_ncstmnb(j)

          IF (.not. ifcalSW2(itype,jtype)) CYCLE

          i2btype = cstmnbtypeindex2(itype,jtype)

          rij(1:3) = atmcor_sig(1:3,i) - atmcor_sig(1:3,j)

!         --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3,i2btype) &
               &   * anint(rij(1:3)*box_inv(1:3,i2btype))

!         * box_inv rather than /box accellerates the computation

          rij_abs2 = rij(1)**2 + rij(2)**2 + rij(3)**2 ! |rij|^2
          rij_abs  = sqrt(rij_abs2)   ! |rij|

!         --- cut off ---
          sqcut = rcut_SW(i2btype) **2

          if (rij_abs2 < sqcut) then

             A_  = A_SW(i2btype)
             B_  = B_SW(i2btype)
             p_  = p_SW(i2btype)
             q_  = q_SW(i2btype)
             np_ = - p_
             nq_ = - q_
             np_1 = np_ - 1
             nq_1 = nq_ - 1

             rij_abs_inv1  = 1.0d0 / rij_abs  ! 1/|rij|

             vcuij  = 1.0d0 / (rij_abs - rcut_SW(i2btype))
!             evcuij = 0.0d0
!             if (vcuij > -30.0d0) evcuij = exp(vcuij)
             if (vcuij > -30.0d0) then
                evcuij = exp(vcuij)
                tmp = B_ * rij_abs_inv1**p_ - rij_abs_inv1**q_
                potij = A_ * tmp * evcuij

                !---------- 2-body force ----------

                ftmp = potij*(((B_*p_*(rij_abs**np_1)-q_*(rij_abs**nq_1)) &
                     &   / (B_*(rij_abs**np_)-(rij_abs**nq_))) + vcuij*vcuij) &
                     &   * rij_abs_inv1

!               --- add both terms ---
                pot_cstmnb = pot_cstmnb + potij

                fi(1:3) = ftmp * rij(1:3)
                fcstmnb(1:3,i) = fcstmnb(1:3,i) + fi(1:3)
                fcstmnb(1:3,j) = fcstmnb(1:3,j) - fi(1:3)

                do n=1,3
                   pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                        &                  + fi(1:3)*rij(n)
                end do

                pot_viri_cstmnb = pot_viri_cstmnb   &
                     &          + fi(1)*atmcor(1,i) &
                     &          + fi(2)*atmcor(2,i) &
                     &          + fi(3)*atmcor(3,i)

                for_viri_cstmnb(1:3,i) = for_viri_cstmnb(1:3,i) + fi(1:3)
                for_viri_cstmnb(1:3,j) = for_viri_cstmnb(1:3,j) - fi(1:3)
             end if

          end if

!         --- 3-body potential terms ---

          k1 = cstmnb_indexall3(cstmnb2)           ! the first and
          k2 = cstmnb_indexall3(cstmnb2+1) - 1     ! the last atoms interacting

          IF (cstmnb_listall3(k1) == 0) CYCLE

          DO kk = k1, k2

             k = cstmnb_listall3(kk)     ! atom index of interacting atom
             ktype = atmindex_ncstmnb(k) ! atmtype of atom(k)

             if (.not. ifcalSW3(itype,jtype,ktype)) cycle

             rij_a = 0.0d0
             vcuij = 0.0d0
             rhoij = 0.0d0

             rik_a = 0.0d0
             vcuik = 0.0d0
             rhoik = 0.0d0

             rjk_a = 0.0d0
             vcujk = 0.0d0
             rhojk = 0.0d0

             i3btype = cstmnbtypeindex3(itype,jtype,ktype)

             rik(1:3) = atmcor_sig(1:3,i) - atmcor_sig(1:3,k)
             rik(1:3) = rik(1:3) - box(1:3,i3btype) &
                  &   * anint(rik(1:3)*box_inv(1:3,i3btype))
             rik_abs2 = rik(1)**2 + rik(2)**2 + rik(3)**2
             rik_abs = sqrt(rik_abs2)

             rjk(1:3) = atmcor_sig(1:3,j) - atmcor_sig(1:3,k)
             rjk(1:3) = rjk(1:3) - box(1:3,i3btype) &
                  &   * anint(rjk(1:3)*box_inv(1:3,i3btype))
             rjk_abs2 = rjk(1)**2 + rjk(2)**2 + rjk(3)**2
             rjk_abs = sqrt(rjk_abs2)

             gamma  =  gamma_SW(i3btype)
             lambda = lambda_SW(i3btype)
             costheta0 = cos(theta0_SW(i3btype))
             a = rcut_SW(i3btype)                 ! parameter of 2-body pot.

             iclose = 0

             if (rij_abs2 < sqcut) then

!                rij_abs = sqrt(rij_abs2)
                rij_a = rij_abs - a
                rhoij = rij_abs * rij_a * rij_a
                tmp = gamma / rij_a
                if (tmp > -30.0d0) then
                   vcuij = exp(tmp)
                   iclose = iclose + 1
                end if
             end if

             if (rik_abs2 < sqcut) then

                rik_a = rik_abs - a
                rhoik = rik_abs * rik_a * rik_a
                tmp = gamma / rik_a
                if (tmp > -30.0d0) then
                   vcuik = exp(tmp)
                   iclose = iclose + 1
                end if
             end if

             if (rjk_abs2 < sqcut) then

                rjk_a = rjk_abs - a
                rhojk = rjk_abs * rjk_a * rjk_a
                tmp = gamma / rjk_a
                if (tmp > -30.0d0) then
                   vcujk = exp(tmp)
                   iclose = iclose + 1
                end if
             end if

             if (iclose < 2) cycle   !  Skip the following calculations: 
!            Forces and potentials of i,j,k turn out to be zero
!            if not less than two of rxx_abs are in the cutoff distance

             hjik = 0.0d0
             hijk = 0.0d0
             hikj = 0.0d0
             hij  = 0.0d0
             hik  = 0.0d0
             hjk  = 0.0d0
             ei   = 0.0d0
             ej   = 0.0d0
             ek   = 0.0d0
             ri_inv = 0.0d0
             rj_inv = 0.0d0
             rk_inv = 0.0d0
             cosjik = 0.0d0
             cosijk = 0.0d0
             cosikj = 0.0d0
             cosjik_0 = 0.0d0
             cosijk_0 = 0.0d0
             cosikj_0 = 0.0d0

             eiNjik(1:3) = 0.0d0
             eiNkij(1:3) = 0.0d0
             ejNijk(1:3) = 0.0d0
             ejNkji(1:3) = 0.0d0
             ekNikj(1:3) = 0.0d0
             ekNjki(1:3) = 0.0d0

             vcui = vcuij * vcuik
             if (vcui > DELTA) then   ! if (vcui /= 0)
                ri_inv = 1.0d0 / (rij_abs * rik_abs)
                cosjik   = (rij(1)*rik(1) + rij(2)*rik(2) + rij(3)*rik(3)) &
                     &   * ri_inv
                cosjik_0 = cosjik - costheta0
                tmp  = lambda * vcui * cosjik_0
                ei   = 2.0d0 * tmp * ri_inv
                hjik = tmp * cosjik_0
             end if

             vcuj = vcuij * vcujk
             if (vcuj > DELTA) then   ! if (vcuj /= 0)
                rj_inv = 1.0d0 / (rij_abs * rjk_abs)
                cosijk   = -(rij(1)*rjk(1) + rij(2)*rjk(2) + rij(3)*rjk(3)) &
                     &   * rj_inv
                cosijk_0 = cosijk - costheta0
                tmp  = lambda * vcuj * cosijk_0
                ej   = 2.0d0 * tmp * rj_inv
                hijk = tmp * cosijk_0
             end if

             vcuk = vcuik * vcujk
             if (vcuk > DELTA) then   ! if (vcuk /= 0)
                rk_inv = 1.0d0 / (rik_abs * rjk_abs)
                cosikj   = (rjk(1)*rik(1) + rjk(2)*rik(2) + rjk(3)*rik(3)) &
                     &   * rk_inv
                cosikj_0 = cosikj - costheta0
                tmp  = lambda * vcuk * cosikj_0
                ek   = 2.0d0 * tmp * rk_inv
                hikj = tmp * cosikj_0
             end if

!----------- 3-body potential ----------

             pot_cstmnb = pot_cstmnb + (hjik + hijk + hikj)

!----------- 3-body force ----------

             eiNjik(1:3) = ei*( rik(1:3) - rik_abs/rij_abs*cosjik*rij(1:3))
             eiNkij(1:3) = ei*( rij(1:3) - rij_abs/rik_abs*cosjik*rik(1:3))
             ejNijk(1:3) = ej*( rjk(1:3) + rjk_abs/rij_abs*cosijk*rij(1:3))
             ejNkji(1:3) = ej*(-rij(1:3) - rij_abs/rjk_abs*cosijk*rjk(1:3))
             ekNikj(1:3) = ek*(-rjk(1:3) + rjk_abs/rik_abs*cosikj*rik(1:3))
             ekNjki(1:3) = ek*(-rik(1:3) + rik_abs/rjk_abs*cosikj*rjk(1:3))

             if (rhoij > DELTA) hij = gamma * (hijk + hjik)/rhoij
             if (rhoik > DELTA) hik = gamma * (hjik + hikj)/rhoik
             if (rhojk > DELTA) hjk = gamma * (hikj + hijk)/rhojk

             hRij(1:3) = hij * rij(1:3)
             hRik(1:3) = hik * rik(1:3)
             hRjk(1:3) = hjk * rjk(1:3)

             fi(1:3) =   hRij(1:3) + hRik(1:3) &
                  &  - eiNjik(1:3) - eiNkij(1:3) + ejNijk(1:3) + ekNikj(1:3)
             fj(1:3) = - hRij(1:3) + hRjk(1:3) &
                  &  - ejNijk(1:3) - ejNkji(1:3) + eiNjik(1:3) + ekNjki(1:3)
             fk(1:3) = - hRjk(1:3) - hRik(1:3) &
                  &  - ekNikj(1:3) - ekNjki(1:3) + eiNkij(1:3) + ejNkji(1:3)

             fcstmnb(1:3,i) = fcstmnb(1:3,i) + fi(1:3)
             fcstmnb(1:3,j) = fcstmnb(1:3,j) + fj(1:3)
             fcstmnb(1:3,k) = fcstmnb(1:3,k) + fk(1:3)

!         --- calculate virial tensor ---
             do n=1,3
                pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                     &                  + fi(1:3)*atmcor(n,i)
                pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                     &                  + fj(1:3)*atmcor(n,j)
                pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                     &                  + fk(1:3)*atmcor(n,k)
             end do

             pot_viri_cstmnb = pot_viri_cstmnb &
                  & + fi(1)*atmcor(1,i) + fi(2)*atmcor(2,i) + fi(3)*atmcor(3,i)
             pot_viri_cstmnb = pot_viri_cstmnb &
                  & + fj(1)*atmcor(1,j) + fj(2)*atmcor(2,j) + fj(3)*atmcor(3,j)
             pot_viri_cstmnb = pot_viri_cstmnb &
                  & + fk(1)*atmcor(1,k) + fk(2)*atmcor(2,k) + fk(3)*atmcor(3,k)

             for_viri_cstmnb(1:3,i) = for_viri_cstmnb(1:3,i) + fi(1:3)
             for_viri_cstmnb(1:3,j) = for_viri_cstmnb(1:3,j) + fj(1:3)
             for_viri_cstmnb(1:3,k) = for_viri_cstmnb(1:3,k) + fk(1:3)

          END DO
       END DO
    END DO

!   --- ADD FCSTMNB TO FORCE ---
    force(1:3,1:natom) = force(1:3,1:natom) + fcstmnb(1:3,1:natom)

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
    real(8):: fcstmnb(3,natom)          ! custom NB force calculated here
    real(8):: atmcor_sig(3,natom)        ! non-dimenstioinal atom coodinate

    real(8):: A_, B_                    ! parameters of SW 2-body potential
    integer:: p_, q_                    ! 
    integer:: np_, nq_                  ! 
    integer:: np_1, nq_1                ! 
    real(8):: gamma, lambda, costheta0 ! parameters of SW 3-body potential

    real(8):: pi                  ! = 3.14
    real(8):: potij               ! temporal potential due to IJ
    real(8):: ftmp                ! temporal force
    real(8):: fi(3), fj(3), fk(3) ! force to be subtracted due to
                                     ! excluded atoms

    real(8):: rij(3),   rik(3),   rjk(3)                ! rj-ri
    real(8):: rij_abs,  rik_abs,  rjk_abs               ! |rij|
    real(8):: rij_abs2, rik_abs2, rjk_abs2              ! |rij|**2
    real(8):: rij_abs_inv1, rik_abs_inv1, rjk_abs_inv1  ! 1/|rij|
    real(8):: rij_abs_inv2, rik_abs_inv2, rjk_abs_inv2  ! 1/|rij|**2

    real(8):: box(3,ncstmnbtyp)      ! BOX size
    real(8):: box_inv(3,ncstmnbtyp)  ! inverse of BOX size

    real(8):: sqcut               ! rcut*rcut

    integer:: icstmnbtype, jcstmnbtype, kcstmnbtype ! SW atom types
    integer:: itype, jtype, ktype ! VDW atom types
    integer:: i,j,k,n             ! do loop index
    integer:: j1, j2              ! limits of do loop
    integer:: k1, k2              ! limits of do loop
    integer:: jj, kk              ! atom index of the interacting atom

    integer:: iclose              ! judge calculation or not for force and potential

    integer:: cstmnb1, cstmnb2, cstmnb3
    integer:: icmtype, i2btype, i3btype ! cstmnb type of common, 2-body & 3-body

    real(8):: tmp

    real(8):: a
    real(8):: rij_a, rik_a, rjk_a
    real(8):: vcuij,  vcuik,  vcujk        ! vcuij  = 1 / rij_a  
    real(8):: evcuij, evcuik, evcujk       ! evcuij = exp(vcuij)  
    real(8):: rhoij, rhoik, rhojk

    real(8):: ri_inv,rj_inv,rk_inv
    real(8):: vcui,vcuj,vcuk
    real(8):: ei,ej,ek

    real(8):: hjik,hijk,hikj
    real(8):: hij, hik, hjk

    real(8):: cosjik,  cosijk,  cosikj
    real(8):: cosjik_0,cosijk_0,cosikj_0

    real(8):: eiNjik(3), eiNkij(3)
    real(8):: ejNijk(3), ejNkji(3)
    real(8):: ekNikj(3), ekNjki(3)

    real(8):: hRij(3), hRik(3), hRjk(3)

    real(8):: nbodycoeff2          ! 2-body coefficient for virial term of heatf
    real(8):: nbodycoeff3          ! 3-body coefficient for virial term of heatf

    real(8):: box_hf(3)            ! BOX size
    real(8):: box_inv_hf(3)        ! inverse of BOX size

    real(8):: DELTA = 1.0d-32   ! extremely small value

! FUNCTIONS:

!     +     +     +     +     +     +     +

!---- some initialization

    fcstmnb(1:3,1:natom) = 0.0d0
    pot_cstmnb = 0.0d0

!     --- CALCULATE NONBONDED FORCE ---

    pot_virit_cstmnb(1:3,1:3) = 0.0d0

    pot_viri_cstmnb = 0.0d0

    nbodycoeff2 = 0.5d0          ! = 1/2
    nbodycoeff3 = 1.0d0/3.0d0    ! = 1/3

    do cstmnb1 = 1, ncstmnb
       i = ncstmnblist(cstmnb1)
       itype = atmindex_ncstmnb(i)
       atmcor_sig(1:3,i) = atmcor(1:3,i) / sigma2_SW(itype)

       box(1,itype) = xcel / sigma2_SW(itype)
       box(2,itype) = ycel / sigma2_SW(itype)
       box(3,itype) = zcel / sigma2_SW(itype)
       box_inv(1:3,itype) = 1.0d0/box(1:3,itype)
    end do

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, ncstmnb - 1
    looplast = ncstmnb - 1
    DO cstmnb1 = loopinit, looplast, loopstep

       j1 = cstmnb_indexall2(cstmnb1)              ! the first and
       j2 = cstmnb_indexall2(cstmnb1+loopstep) - 1 ! the last atoms interacting
                                                  ! with i
       i = ncstmnblist(cstmnb1)
       itype = atmindex_ncstmnb(i)       ! atom type of atom(i)

       IF (cstmnb_listall2(j1) == 0) CYCLE

       DO cstmnb2 = j1, j2

          j = cstmnb_listall2(cstmnb2)    ! atom index of interacting atom
          jtype = atmindex_ncstmnb(j)

          IF (.not. ifcalSW2(itype,jtype)) CYCLE

          i2btype = cstmnbtypeindex2(itype,jtype)

          rij(1:3) = atmcor_sig(1:3,i) - atmcor_sig(1:3,j)

!         --- periodic boundary ---

          rij(1:3) = rij(1:3) - box(1:3,i2btype) &
               &   * anint(rij(1:3)*box_inv(1:3,i2btype))

!         * box_inv rather than /box accellerates the computation

          rij_abs2 = rij(1)**2 + rij(2)**2 + rij(3)**2 ! |rij|^2
          rij_abs  = sqrt(rij_abs2)   ! |rij|

!         --- cut off ---
          sqcut = rcut_SW(i2btype) **2

          if (rij_abs2 < sqcut) then

             A_  = A_SW(i2btype)
             B_  = B_SW(i2btype)
             p_  = p_SW(i2btype)
             q_  = q_SW(i2btype)
             np_ = - p_
             nq_ = - q_
             np_1 = np_ - 1
             nq_1 = nq_ - 1

             rij_abs_inv1  = 1.0d0 / rij_abs  ! 1/|rij|

             vcuij  = 1.0d0 / (rij_abs - rcut_SW(i2btype))
!             evcuij = 0.0d0
!             if (vcuij > -30.0d0) evcuij = exp(vcuij)
             if (vcuij > -30.0d0) then
                evcuij = exp(vcuij)
                tmp = B_ * rij_abs_inv1**p_ - rij_abs_inv1**q_
                potij = A_ * tmp * evcuij

                !---------- 2-body force ----------

                ftmp = potij*(((B_*p_*(rij_abs**np_1)-q_*(rij_abs**nq_1)) &
                     &   / (B_*(rij_abs**np_)-(rij_abs**nq_))) + vcuij*vcuij) &
                     &   * rij_abs_inv1

!               --- add both terms ---
                pot_cstmnb = pot_cstmnb + potij

                fi(1:3) = ftmp * rij(1:3)
                fcstmnb(1:3,i) = fcstmnb(1:3,i) + fi(1:3)
                fcstmnb(1:3,j) = fcstmnb(1:3,j) - fi(1:3)

!               - for heat flux calculation
                fj(1:3) = -fi(1:3)

                box_hf(1:3) = box(1:3,i2btype)
                box_inv_hf(1:3) = box_inv(1:3,i2btype)

                call calvirihf(i,j,fi,fj, &
                     &         nbodycoeff2, &
                     &         box_hf,box_inv_hf, &
                     &         ifhfvol, &
                     &         nhfregion,hfzpos1,hfzpos2, &
                     &         hftyp_atm, &
                     &         molecom, &
                     &         viricstmnbt_atm)

                do n=1,3
                   pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                        &                  + fi(1:3)*rij(n)
                end do

                pot_viri_cstmnb = pot_viri_cstmnb   &
                     &          + fi(1)*atmcor(1,i) &
                     &          + fi(2)*atmcor(2,i) &
                     &          + fi(3)*atmcor(3,i)

                for_viri_cstmnb(1:3,i) = for_viri_cstmnb(1:3,i) + fi(1:3)
                for_viri_cstmnb(1:3,j) = for_viri_cstmnb(1:3,j) - fi(1:3)
             end if

          end if

!         --- 3-body potential terms ---

          k1 = cstmnb_indexall3(cstmnb2)              ! the first and
          k2 = cstmnb_indexall3(cstmnb2+1) - 1        ! the last atoms interacting

          IF (cstmnb_listall3(k1) == 0) CYCLE

          DO kk = k1, k2

             k = cstmnb_listall3(kk)     ! atom index of interacting atom
             ktype = atmindex_ncstmnb(k) ! atmtype of atom(k)

             if (.not. ifcalSW3(itype,jtype,ktype)) cycle

             rij_a = 0.0d0
             vcuij = 0.0d0
             rhoij = 0.0d0

             rik_a = 0.0d0
             vcuik = 0.0d0
             rhoik = 0.0d0

             rjk_a = 0.0d0
             vcujk = 0.0d0
             rhojk = 0.0d0

             i3btype = cstmnbtypeindex3(itype,jtype,ktype)

             rik(1:3) = atmcor_sig(1:3,i) - atmcor_sig(1:3,k)
             rik(1:3) = rik(1:3) - box(1:3,i3btype) &
                  &   * anint(rik(1:3)*box_inv(1:3,i3btype))
             rik_abs2 = rik(1)**2 + rik(2)**2 + rik(3)**2
             rik_abs = sqrt(rik_abs2)

             rjk(1:3) = atmcor_sig(1:3,j) - atmcor_sig(1:3,k)
             rjk(1:3) = rjk(1:3) - box(1:3,i3btype) &
                  &   * anint(rjk(1:3)*box_inv(1:3,i3btype))
             rjk_abs2 = rjk(1)**2 + rjk(2)**2 + rjk(3)**2
             rjk_abs = sqrt(rjk_abs2)

             gamma  =  gamma_SW(i3btype)
             lambda = lambda_SW(i3btype)
             costheta0 = cos(theta0_SW(i3btype))
             a = rcut_SW(i3btype)                 ! parameter of 2-body pot.

             if (rij_abs2 < sqcut) then
!                rij_abs = sqrt(rij_abs2)
                rij_a = rij_abs - a
                rhoij = rij_abs * rij_a * rij_a
                tmp = gamma / rij_a
                if (tmp > -30.0d0) then
                   vcuij = exp(tmp)
                   iclose = iclose + 1
                end if
             end if

             if (rik_abs2 < sqcut) then
                rik_a = rik_abs - a
                rhoik = rik_abs * rik_a * rik_a
                tmp = gamma / rik_a
                if (tmp > -30.0d0) then
                   vcuik = exp(tmp)
                   iclose = iclose + 1
                end if
             end if

             if (rjk_abs2 < sqcut) then
                rjk_a = rjk_abs - a
                rhojk = rjk_abs * rjk_a * rjk_a
                tmp = gamma / rjk_a
                if (tmp > -30.0d0) then
                   vcujk = exp(tmp)
                   iclose = iclose + 1
                end if
             end if

             if (iclose < 2) cycle   !  Skip the following calculations: 
!            Forces and potentials of i,j,k turn out to be zero
!            if not less than two of rxx_abs are in the cutoff distance

             hjik = 0.0d0
             hijk = 0.0d0
             hikj = 0.0d0
             hij  = 0.0d0
             hik  = 0.0d0
             hjk  = 0.0d0
             ei   = 0.0d0
             ej   = 0.0d0
             ek   = 0.0d0
             ri_inv = 0.0d0
             rj_inv = 0.0d0
             rk_inv = 0.0d0
             cosjik = 0.0d0
             cosijk = 0.0d0
             cosikj = 0.0d0
             cosjik_0 = 0.0d0
             cosijk_0 = 0.0d0
             cosikj_0 = 0.0d0

             eiNjik(1:3) = 0.0d0
             eiNkij(1:3) = 0.0d0
             ejNijk(1:3) = 0.0d0
             ejNkji(1:3) = 0.0d0
             ekNikj(1:3) = 0.0d0
             ekNjki(1:3) = 0.0d0

             vcui = vcuij * vcuik
             if (vcui > DELTA) then   ! if (vcui /= 0)
                ri_inv = 1.0d0 / (rij_abs * rik_abs)
                cosjik   = (rij(1)*rik(1) + rij(2)*rik(2) + rij(3)*rik(3)) &
                     &   * ri_inv
                cosjik_0 = cosjik - costheta0
                tmp  = lambda * vcui * cosjik_0
                ei   = 2.0d0 * tmp * ri_inv
                hjik = tmp * cosjik_0
             end if

             vcuj = vcuij * vcujk
             if (vcuj > DELTA) then   ! if (vcuj /= 0)
                rj_inv = 1.0d0 / (rij_abs * rjk_abs)
                cosijk   = -(rij(1)*rjk(1) + rij(2)*rjk(2) + rij(3)*rjk(3)) &
                     &   * rj_inv
                cosijk_0 = cosijk - costheta0
                tmp  = lambda * vcuj * cosijk_0
                ej   = 2.0d0 * tmp * rj_inv
                hijk = tmp * cosijk_0
             end if

             vcuk = vcuik * vcujk
             if (vcuk > DELTA) then   ! if (vcuk /= 0)
                rk_inv = 1.0d0 / (rik_abs * rjk_abs)
                cosikj   = (rjk(1)*rik(1) + rjk(2)*rik(2) + rjk(3)*rik(3)) &
                     &   * rk_inv
                cosikj_0 = cosikj - costheta0
                tmp  = lambda * vcuk * cosikj_0
                ek   = 2.0d0 * tmp * rk_inv
                hikj = tmp * cosikj_0
             end if

!----------- 3-body potential ----------

             pot_cstmnb = pot_cstmnb + (hjik + hijk + hikj)

!----------- 3-body force ----------

             eiNjik(1:3) = ei*( rik(1:3) - rik_abs/rij_abs*cosjik*rij(1:3))
             eiNkij(1:3) = ei*( rij(1:3) - rij_abs/rik_abs*cosjik*rik(1:3))
             ejNijk(1:3) = ej*( rjk(1:3) + rjk_abs/rij_abs*cosijk*rij(1:3))
             ejNkji(1:3) = ej*(-rij(1:3) - rij_abs/rjk_abs*cosijk*rjk(1:3))
             ekNikj(1:3) = ek*(-rjk(1:3) + rjk_abs/rik_abs*cosikj*rik(1:3))
             ekNjki(1:3) = ek*(-rik(1:3) + rik_abs/rjk_abs*cosikj*rjk(1:3))

             if (rhoij > DELTA) hij = gamma * (hijk + hjik)/rhoij
             if (rhoik > DELTA) hik = gamma * (hjik + hikj)/rhoik
             if (rhojk > DELTA) hjk = gamma * (hikj + hijk)/rhojk

             hRij(1:3) = hij * rij(1:3)
             hRik(1:3) = hik * rik(1:3)
             hRjk(1:3) = hjk * rjk(1:3)

             fi(1:3) =   hRij(1:3) + hRik(1:3) &
                  &  - eiNjik(1:3) - eiNkij(1:3) + ejNijk(1:3) + ekNikj(1:3)
             fj(1:3) = - hRij(1:3) + hRjk(1:3) &
                  &  - ejNijk(1:3) - ejNkji(1:3) + eiNjik(1:3) + ekNjki(1:3)
             fk(1:3) = - hRjk(1:3) - hRik(1:3) &
                  &  - ekNikj(1:3) - ekNjki(1:3) + eiNkij(1:3) + ejNkji(1:3)

             fcstmnb(1:3,i) = fcstmnb(1:3,i) + fi(1:3)
             fcstmnb(1:3,j) = fcstmnb(1:3,j) + fj(1:3)
             fcstmnb(1:3,k) = fcstmnb(1:3,k) + fk(1:3)
!
!            --- calculate virial tensor ---

             box_hf(1:3) = box(1:3,i3btype)
             box_inv_hf(1:3) = box_inv(1:3,i3btype)

!            - for heat flux calculation
!            pair (i,j)
             call calvirihf(i,j,fi,fj, &
                  &         nbodycoeff3, &
                  &         box_hf,box_inv_hf, &
                  &         ifhfvol, &
                  &         nhfregion,hfzpos1,hfzpos2, &
                  &         hftyp_atm, &
                  &         molecom, &
                  &         viricstmnbt_atm)

!            pair (j,k)
             call calvirihf(j,k,fj,fk, &
                  &         nbodycoeff3, &
                  &         box_hf,box_inv_hf, &
                  &         ifhfvol, &
                  &         nhfregion,hfzpos1,hfzpos2, &
                  &         hftyp_atm, &
                  &         molecom, &
                  &         viricstmnbt_atm)

!            pair (k,i)
             call calvirihf(k,i,fk,fi, &
                  &         nbodycoeff3, &
                  &         box_hf,box_inv_hf, &
                  &         ifhfvol, &
                  &         nhfregion,hfzpos1,hfzpos2, &
                  &         hftyp_atm, &
                  &         molecom, &
                  &         viricstmnbt_atm)

!            - for pressure calculation
             do n=1,3
                pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                     &                  + fi(1:3)*atmcor(n,i)
                pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                     &                  + fj(1:3)*atmcor(n,j)
                pot_virit_cstmnb(1:3,n) = pot_virit_cstmnb(1:3,n) &
                     &                  + fk(1:3)*atmcor(n,k)
             end do

             pot_viri_cstmnb = pot_viri_cstmnb &
                  & + fi(1)*atmcor(1,i) + fi(2)*atmcor(2,i) + fi(3)*atmcor(3,i)
             pot_viri_cstmnb = pot_viri_cstmnb &
                  & + fj(1)*atmcor(1,j) + fj(2)*atmcor(2,j) + fj(3)*atmcor(3,j)
             pot_viri_cstmnb = pot_viri_cstmnb &
                  & + fk(1)*atmcor(1,k) + fk(2)*atmcor(2,k) + fk(3)*atmcor(3,k)

             for_viri_cstmnb(1:3,i) = for_viri_cstmnb(1:3,i) + fi(1:3)
             for_viri_cstmnb(1:3,j) = for_viri_cstmnb(1:3,j) + fj(1:3)
             for_viri_cstmnb(1:3,k) = for_viri_cstmnb(1:3,k) + fk(1:3)

          END DO
       END DO
    END DO

!   --- ADD FCSTMNB TO FORCE ---
    force(1:3,1:natom) = force(1:3,1:natom) + fcstmnb(1:3,1:natom)

!     +     +     +     +     +     +     +

  end subroutine calcstmnbp_hf
#endif

end module cstmnb
