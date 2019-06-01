!**************************************
!*  shpot.f Ver.1.1 '09.11.26         *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine shpot(xcel,ycel,zcel,   &
     &           rcutsh,   &
     &           force,pot_sh)

  use md_global
  use mpi_global

  implicit none

!
!    
!     Calculate SH interactions 
!
!
! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  real(8),intent(in):: rcutsh           ! SH cutoff length [non-d]

!   INPUT & OUTPUT

  real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc) 

!   OUTPUT

  real(8),intent(out):: pot_sh          ! SH potential

! LOCAL:
  real(8):: fsh(3,natom)                ! SH force calculated here
!      real(8):: fshall(3,natom) ! SH force calculated here
!MPI                                       local(in MPI sense) field
  real(8):: potij                       ! temporal potential due to IJ  
  real(8):: ftmp                        ! temporal force  
  real(8):: fij(3)                      ! force to be subtracted due to 
                                        ! excluded atoms

  real(8):: rij(3)                      ! rj-ri
  real(8):: rij_abs                     ! |Rij|   
  real(8):: rij_abs2                    ! |Rij|**2
  real(8):: rij_abs_inv1                ! |Rij|**-1

  real(8):: rho_abs2                    ! rho**2 = rijx**2 + rijy**2
  real(8):: f_rho                       ! f(rho)
  real(8):: sh_1exp                     ! = a1*exp(-b1*|rij|)
  real(8):: sh_2exp                     ! = a2*exp(-b2*|rij|)
  real(8):: sh_3exp                     ! = a3*exp(-b3*|rij|)
  real(8):: sh_1expall                  ! = a1*exp(-b1*|rij|) * f(rho)
  real(8):: sh_2expall                  ! = a2*exp(-b2*|rij|) * f(rho)
  real(8):: sh_3expall                  ! = a3*exp(-b3*|rij|) * (1 - f(rho))
  real(8):: ftmp_z

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

  real(8):: sqcut                       ! rcutsh * rcutsh

!
! MPI  following variable are declared in MPI version
!
!      real*8:: pot_sh_lcl      !local SH potential
!      real*8:: local,global     !tmps for global summation
! MPI end

  integer:: itype, jtype                ! VDW atom types
  integer:: i,j                         ! do loop index
  integer:: j1, j2                      ! limits of do loop
  integer:: jj                          ! atom index of the interacting atom

  integer:: sh1

!      real*8:: xij,xij0
!      real*8:: yij,yij0
!      real*8:: zij,zij0
!      real*8:: r_ij

!      integer:: im,ip,jm,jp,km,kp
!      real*8:: xm,xp,ym,yp,zm,zp

!      real*8:: fxct,fyct,fzct,pott

!           real(kind=dp), external::  ferfc, ferf   
!          * whether erfc and erf are external or intrinsic depends
!            on systems. so, ferfc and ferf are defined in 
!            MACHINE/(system)/functions.f90.
!            by Bestsystems.

!     +     +     +     +     +     +     +


! --- ALLOCATE ARRAY FOR FORCE CALCULATED IN THIS SUBROUTINE --- 

  fsh(1:3,1:natom) = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

  sqcut = rcutsh * rcutsh

! --- CALCULATE NONBONDED FORCE ---

  pot_sh = 0.0d0

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, NATOM - 1

!
! SH O atoms
!
  looplast = nsho - 1

  DO sh1 = loopinit, looplast, loopstep

     j1 = sho_indexall(sh1)              ! the first and
     j2 = sho_indexall(sh1+loopstep) - 1 ! the last atoms interacting
                                         ! with i

     i = nsholist(sh1)
     itype = atmindex(i)                 ! atmtype of atom(i)

     IF (sho_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = sho_listall(j)              ! atom index of interacting atom

        jtype = atmindex(jj)             ! atmtype of atom(jj)

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rho_abs2      = rij(1)**2 + rij(2)**2 ! rho**2

        rij_abs2      = rho_abs2 + rij(3)**2
                                        ! |rij|^2
!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs       = dsqrt(rij_abs2) ! |rij| 

!        rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2
        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

        f_rho = dexp(-para_sh_c * rho_abs2) ! = f(rho)

        sh_1exp = para_sh_a(1) * dexp(-para_sh_b(1)*rij_abs)
                                        ! = a1*exp(-b1*|rij|)
        sh_2exp = para_sh_a(2) * dexp(-para_sh_b(2)*rij_abs)
                                        ! = a2*exp(-b2*|rij|)
        sh_3exp = para_sh_a(3) * dexp(-para_sh_b(3)*rij_abs)
                                        ! = a3*exp(-b3*|rij|)

        sh_1expall = sh_1exp * f_rho    ! = a1*exp(-b1*|rij|) * f(rho)
        sh_2expall = sh_2exp * f_rho    ! = a2*exp(-b2*|rij|) * f(rho)
        sh_3expall = sh_3exp * (1.0d0 - f_rho)
                                        ! = a3*exp(-b3*|rij|) * (1 - f(rho))

        potij = sh_1expall - sh_2expall + sh_3expall

        ftmp_z = (sh_1expall * para_sh_b(1)   &
             &  - sh_2expall * para_sh_b(2)   &
             &  + sh_3expall * para_sh_b(3))   &
             & * rij_abs_inv1

        ftmp = ftmp_z   &
           & + (sh_1expall - sh_2expall - sh_3exp*f_rho)   &
           & * 2.0d0 * para_sh_c

        pot_sh    = pot_sh + potij

        fij(1)     = ftmp   * rij(1)
        fij(2)     = ftmp   * rij(2)
        fij(3)     = ftmp_z * rij(3)

        fsh(1:3,i)  = fsh(1:3,i)  + fij(1:3)
        fsh(1:3,jj) = fsh(1:3,jj) - fij(1:3)

     END DO
  END DO

!
! SH H atoms
!
  looplast = nshh - 1

  DO sh1 = loopinit, looplast, loopstep

     j1 = shh_indexall(sh1)              ! the first and
     j2 = shh_indexall(sh1+loopstep) - 1 ! the last atoms interacting
                                         ! with i

     i = nshhlist(sh1)
     itype = atmindex(i)                ! atmtype of atom(i)

     IF (shh_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = shh_listall(j)             ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2
                                        ! |rij|^2

!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs       = dsqrt(rij_abs2) ! |rij| 

!        rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2
        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

        potij = para_sh_a(4) * dexp(-para_sh_b(4)*rij_abs)
                                        ! = a4*exp(-b4*|rij|)

        ftmp = potij * para_sh_b(4) * rij_abs_inv1

        pot_sh    = pot_sh + potij

        fij(1:3)     = ftmp * rij(1:3)

        fsh(1:3,i)  = fsh(1:3,i)  + fij(1:3)
        fsh(1:3,jj) = fsh(1:3,jj) - fij(1:3)

     END DO
  END DO


! --- ADD FSH TO FORCE ---

  force(1:3,1:natom) = force(1:3,1:natom) + fsh(1:3,1:natom)

!     write(20,'(6f12.3)') fshall(:,:)

!     +     +     +     +     +     +     +

end subroutine shpot

!----------------------------------------------------------
subroutine shpotp(xcel,ycel,zcel,   &
     &            rcutsh,   &
     &            force,pot_sh,   &
     &            for_viri_sh,pot_viri_sh,   &
     &            pot_virit_sh)

  use md_global
  use mpi_global

  implicit none
      

!
!    
!     Calculate SH interactions 
!
!
! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  real(8),intent(in):: rcutsh           ! SH cutoff length [non-d]

!   INPUT & OUTPUT

  real(8),intent(inout):: force(:,:)   ! force calculated here
                                       ! total force
                                       !  (including vdw, bond, angle etc)

  real(8),intent(inout):: for_viri_sh(:,:)   ! virial(SH force) of each atom
  real(8),intent(inout):: pot_viri_sh        ! virial(SH potential)
  real(8),intent(inout):: pot_virit_sh(:,:)  ! virial tensor (SH)

!   OUTPUT

  real(8),intent(out):: pot_sh         ! SH potential

! LOCAL:
  real(8):: fsh(3,natom)               ! SH force calculated here
!  real*8:: fshall(3,natom) ! Sh force calculated here
!MPI                                       local(in MPI sense) field
  real(8):: potij                      ! temporal potential due to IJ  
  real(8):: ftmp                       ! temporal force  
  real(8):: fij(3)                     ! force to be subtracted due to 
                                       ! excluded atoms

  real(8):: rij(3)                     ! rj-ri
  real(8):: rij_abs                    ! |Rij|   
  real(8):: rij_abs2                   ! |Rij|**2
  real(8):: rij_abs_inv1               ! |Rij|**-1

  real(8):: rho_abs2                   ! rho**2 = rijx**2 + rijy**2
  real(8):: f_rho                      ! f(rho)
  real(8):: sh_1exp                    ! = a1*exp(-b1*|rij|)
  real(8):: sh_2exp                    ! = a2*exp(-b2*|rij|)
  real(8):: sh_3exp                    ! = a3*exp(-b3*|rij|)
  real(8):: sh_1expall                 ! = a1*exp(-b1*|rij|) * f(rho)
  real(8):: sh_2expall                 ! = a2*exp(-b2*|rij|) * f(rho)
  real(8):: sh_3expall                 ! = a3*exp(-b3*|rij|) * (1 - f(rho))
  real(8):: ftmp_z

  real(8):: box(3)                     ! BOX size
  real(8):: box_inv(3)                 ! inverse of BOX size

  real(8):: sqcut                      ! rcutsh * rcutsh

!
! MPI  following variable are declared in MPI version
!
!      real*8:: pot_sh_lcl      !local elc potential
!      real*8:: local(3*3+2),global(3*3+2) !tmps for global summation
!
!      real*8:: pvtsh_lcl(3,3)   ! pot_virit_sh   local array
!
!      real*8:: pvsh_lcl         ! pot_viri_sh   local array
!
! MPI end

  integer:: itype, jtype               ! VDW atom types
  integer:: i,j,n                      ! do loop index
  integer:: j1, j2                     ! limits of do loop
  integer:: jj                         ! atom index of the interacting atom

  integer:: sh1

!      integer:: im,i1,i2,ii,jm  ! loop index of molecule

!      logical:: exclatom        ! flag for if atom belongs list_excl

!      real*8:: xij,xij0
!      real*8:: yij,yij0
!      real*8:: zij,zij0
!      real*8:: r_ij

!      integer:: im,ip,jm,jp,km,kp
!      real*8:: xm,xp,ym,yp,zm,zp

!      real*8:: fxct,fyct,fzct,pott

!           real(kind=dp), external::  ferfc, ferf   
!          * whether erfc and erf are external or intrinsic depends
!            on systems. so, ferfc and ferf are defined in 
!            MACHINE/(system)/functions.f90.
!            by Bestsystems.

!     +     +     +     +     +     +     +


! --- ALLOCATE ARRAY FOR FORCE CALCULATED IN THIS SUBROUTINE --- 

  fsh(1:3,1:natom) = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

  sqcut = rcutsh * rcutsh
      
! --- CALCULATE NONBONDED FORCE ---

  pot_sh = 0.0d0

  pot_virit_sh(1:3,1:3) = 0.0d0

  pot_viri_sh   = 0.0d0

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, NATOM - 1

!
! SH O atoms
!      
  looplast = nsho - 1

  DO sh1 = loopinit, looplast, loopstep

     j1 = sho_indexall(sh1)              ! the first and
     j2 = sho_indexall(sh1+loopstep) - 1 ! the last atoms interacting
                                         ! with i

     i = nsholist(sh1)
     itype = atmindex(i)                ! atmtype of atom(i)

!     im = irmolept_list(i)  ! what molecule atom(i) belongs in

     IF (sho_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = sho_listall(j)             ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

!            jm = irmolept_list(jj) ! what molecule atom(i) belongs in

!       --- CALCULATE VIRIAL & FORCE ---

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rho_abs2      = rij(1)**2 + rij(2)**2 ! rho**2

        rij_abs2      = rho_abs2 + rij(3)**2 ! |rij|^2
 
!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs       = dsqrt(rij_abs2) ! |rij| 

!        rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2
        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

        f_rho = dexp(-para_sh_c * rho_abs2) ! = f(rho)

        sh_1exp = para_sh_a(1) * dexp(-para_sh_b(1)*rij_abs)
                                        ! = a1*exp(-b1*|rij|)
        sh_2exp = para_sh_a(2) * dexp(-para_sh_b(2)*rij_abs)
                                        ! = a2*exp(-b2*|rij|)
        sh_3exp = para_sh_a(3) * dexp(-para_sh_b(3)*rij_abs)
                                        ! = a3*exp(-b3*|rij|)

        sh_1expall = sh_1exp * f_rho    ! = a1*exp(-b1*|rij|) * f(rho)
        sh_2expall = sh_2exp * f_rho    ! = a2*exp(-b2*|rij|) * f(rho)
        sh_3expall = sh_3exp * (1.0d0 - f_rho)
                                        ! = a3*exp(-b3*|rij|) * (1 - f(rho))

        potij = sh_1expall - sh_2expall + sh_3expall

        ftmp_z = (sh_1expall * para_sh_b(1)   &
             &  - sh_2expall * para_sh_b(2)   &
             &  + sh_3expall * para_sh_b(3))   &
             & * rij_abs_inv1

        ftmp = ftmp_z   &
           & + (sh_1expall - sh_2expall - sh_3exp*f_rho)   &
           & * 2.0d0 * para_sh_c

        pot_sh    = pot_sh + potij

        fij(1)     = ftmp   * rij(1)
        fij(2)     = ftmp   * rij(2)
        fij(3)     = ftmp_z * rij(3)

        fsh(1:3,i)  = fsh(1:3,i)  + fij(1:3)
        fsh(1:3,jj) = fsh(1:3,jj) - fij(1:3)

!       --- calculate virial tensor ---
        do n=1,3
           pot_virit_sh(1:3,n) = pot_virit_sh(1:3,n)   &
                &              + fij(1:3)*rij(n)
        end do
               
        pot_viri_sh = pot_viri_sh   &
             &      + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

     END DO

  END DO

!
! SH H atoms
!      
  looplast = nshh - 1

  DO sh1 = loopinit, looplast, loopstep

     j1 = shh_indexall(sh1)              ! the first and
     j2 = shh_indexall(sh1+loopstep) - 1 ! the last atoms interacting
                                         ! with i

     i = nshhlist(sh1)
     itype = atmindex(i)                ! atmtype of atom(i)

!     im = irmolept_list(i)  ! what molecule atom(i) belongs in

     IF (shh_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = shh_listall(j)             ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

!            jm = irmolept_list(jj) ! what molecule atom(i) belongs in

!       --- CALCULATE VIRIAL & FORCE ---

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2 ! |rij|^2
 
!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs       = dsqrt(rij_abs2) ! |rij| 

!        rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2
        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

        potij = para_sh_a(4) * dexp(-para_sh_b(4)*rij_abs)
                                        ! = a4*exp(-b4*|rij|)

        ftmp = potij * para_sh_b(4) * rij_abs_inv1

        pot_sh    = pot_sh + potij

        fij(1:3)     = ftmp * rij(1:3)

        fsh(1:3,i)  = fsh(1:3,i)  + fij(1:3)
        fsh(1:3,jj) = fsh(1:3,jj) - fij(1:3)

!       --- calculate virial tensor ---
        do n=1,3
           pot_virit_sh(1:3,n) = pot_virit_sh(1:3,n)   &
                &              + fij(1:3)*rij(n)
        end do
               
        pot_viri_sh = pot_viri_sh   &
             &      + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

     END DO

  END DO


! --- ADD FSH TO FORCE ---

  force(1:3,1:natom) = force(1:3,1:natom) + fsh(1:3,1:natom)
  for_viri_sh(1:3,1:natom) = for_viri_sh(1:3,1:natom) + fsh(1:3,1:natom)

!     +     +     +     +     +     +     +

end subroutine shpotp
