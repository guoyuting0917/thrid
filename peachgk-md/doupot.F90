!**************************************
!*  doupot.f Ver.1.2 '10.06.28        *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine doupot(xcel,ycel,zcel,   &
     &            rcutdouo,rcutindouo,rcutdouh,   &
     &            force,pot_dou)

  use md_global
  use mpi_global

  implicit none


!
!    
!     Calculate DOU interactions
!       for Au-water interaction proposed by Dou et al.
!       ref. Dou, Y., Zhigilei, L.V., Winograd, N., Garrison, B.J.,
!       J. Phys. Chem. A, 105 (2001), pp. 2748-2755.
!
!
! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  real(8),intent(in):: rcutdouo         ! DOU cutoff length [non-d] for O-Au
  real(8),intent(in):: rcutindouo       ! DOU cutin length [non-d] for O-Au
  real(8),intent(in):: rcutdouh         ! DOU cutoff length [non-d] for H-Au

!   INPUT & OUTPUT

  real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc) 

!   OUTPUT

  real(8),intent(out):: pot_dou         ! DOU potential

! LOCAL:
  real(8):: fdou(3,natom)               ! DOU force calculated here
!      real*8:: fdouall(3,natom) ! DOU force calculated here
!MPI                                       local(in MPI sense) field
  real(8):: potij                       ! temporal potential due to IJ
  real(8):: ftmp                        ! temporal force
  real(8):: fij(3)                      ! force to be subtracted due to
                                        ! excluded atoms

  real(8):: rij(3)                      ! rj-ri
  real(8):: rij_abs                     ! |Rij|   
  real(8):: rij_abs2                    ! |Rij|**2   
!  real(8):: rij_abs_inv1     ! |Rij|**-1

  real(8):: rtmp1                       ! exp(-beta*(|rij| - r0))
  real(8):: rtmp2                       ! exp(-2*beta*(|rij| - r0))

  real(8):: roff_r2                     ! r_off**2 - |rij|**2
  real(8):: r2_ron                      ! |rij|**2 - r_on**2


  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

  real(8):: sqcut                       ! rcutdoux * rcutdoux
  real(8):: sqcutin                     ! rcutindouo * rcutindouo

  real(8):: denom_rcut                  ! 1/(r_off**2 - r_on**2)**3

  real(8):: S2                          ! switching function

  integer:: itype, jtype                ! VDW atom types
  integer:: i,j                         ! do loop index
  integer:: j1, j2                      ! limits of do loop
  integer:: jj                          ! atom index of the interacting atom

  integer:: dou1

!      real*8:: xij,xij0
!      real*8:: yij,yij0
!      real*8:: zij,zij0
!      real*8:: r_ij

!      integer:: im,ip,jm,jp,km,kp
!      real*8:: xm,xp,ym,yp,zm,zp

!      real*8:: fxct,fyct,fzct,pott

!     +     +     +     +     +     +     +


! --- ALLOCATE ARRAY FOR FORCE CALCULATED IN THIS SUBROUTINE --- 

  fdou(1:3,1:natom) = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

! --- CALCULATE NONBONDED FORCE ---

  pot_dou = 0.0d0

!
! - DOU (O) interaction
!

  sqcut = rcutdouo * rcutdouo
  sqcutin = rcutindouo * rcutindouo

  denom_rcut = 1.0d0 / (sqcut - sqcutin) ! = 1 / (r_off**2 - r_on**2)
  denom_rcut = denom_rcut * denom_rcut * denom_rcut
                                        ! = 1 / (r_off**2 - r_on**2)**3

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, NATOM - 1
  looplast = ndou - 1

  DO dou1 = loopinit, looplast, loopstep

     j1 = douo_indexall(dou1)              ! the first and
     j2 = douo_indexall(dou1+loopstep) - 1 ! the last atoms interacting
                                           ! with i

     i = ndoulist(dou1)
     itype = atmindex(i)      ! atmtype of atom(i)

     IF (douo_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = douo_listall(j)            ! atom index of interacting atom

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

!        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

        rtmp1 = dexp(-para_dou_beta_o*(rij_abs - para_dou_rad_o))
                                        ! exp(-beta*(|rij|-r0))
        rtmp2 = rtmp1 * rtmp1           ! exp(-2*beta*(|rij|-r0))

        ftmp  = 2.0d0 * para_dou_beta_o   &
            & * para_dou_wel_o   &
            & * (rtmp2 - rtmp1) / rij_abs
                                        ! = 2*beta*D*(rtmp2-rtmp1) / |rij|

        potij = para_dou_wel_o * (rtmp2 - 2.0d0*rtmp1)
                                        ! = D * (rtmp2 - 2*rtmp1)

        if (rij_abs2 > sqcutin) then

           roff_r2 = sqcut - rij_abs2   ! r_off**2 - |rij|**2
           r2_ron  = rij_abs2 - sqcutin ! |rij|**2 - r_on**2

           S2 = roff_r2*roff_r2 * (roff_r2 + 3.0d0*r2_ron)
                              ! = (r_off**2 - |rij|**2)**2 
                              ! * ((roff**2 - |rij|**2) + 3*(|rij|**2 - ron**2))

           ftmp  = ftmp * S2   &
               & + potij * 12.0d0 * roff_r2 * r2_ron
           ftmp  = ftmp * denom_rcut
               
           potij = potij * S2 * denom_rcut

        end if

        pot_dou    = pot_dou + potij

        fij(1:3)     = ftmp * rij(1:3)

        fdou(1:3,i)  = fdou(1:3,i)  + fij(1:3)
        fdou(1:3,jj) = fdou(1:3,jj) - fij(1:3)

     END DO
  END DO

!
! - DOU (H) interaction
!

  sqcut = rcutdouh * rcutdouh

  looplast = ndou - 1

  DO dou1 = loopinit, looplast, loopstep

     j1 = douh_indexall(dou1)              ! the first and
     j2 = douh_indexall(dou1+loopstep) - 1 ! the last atoms interacting
                                           ! with i

     i = ndoulist(dou1)
     itype = atmindex(i)                ! atmtype of atom(i)

     IF (douh_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = douh_listall(j)            ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2 ! |rij|^2

!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs       = dsqrt(rij_abs2) ! |rij|

!        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|
            
        potij = para_dou_gamma * para_dou_wel_h   &
            & * dexp(-2.0d0 * para_dou_beta_h   &
            &      * (rij_abs - para_dou_rad_h))
                                        ! = gamma*D * exp(-2*beta*(|rij|-r0))

        ftmp = 2.0d0 * para_dou_beta_h * potij / rij_abs
                                    ! = 2*beta*gamma*D * exp(-2*beta*(|rij|-r0))
                                    !   / |rij|

        pot_dou      = pot_dou + potij

        fij(1:3)     = ftmp * rij(1:3)

        fdou(1:3,i)  = fdou(1:3,i)  + fij(1:3)
        fdou(1:3,jj) = fdou(1:3,jj) - fij(1:3)

     END DO
  END DO

! --- ADD FDOU TO FORCE ---

  force(1:3,1:natom) = force(1:3,1:natom) + fdou(1:3,1:natom)

!     write(20,'(6f12.3)') fdouall(:,:)

!     +     +     +     +     +     +     +

end subroutine doupot

!----------------------------------------------------------
subroutine doupotp(xcel,ycel,zcel,   &
     &             rcutdouo,rcutindouo,rcutdouh,   &
     &             force,pot_dou,   &
     &             for_viri_dou,pot_viri_dou,   &
     &             pot_virit_dou)
  
  use md_global
  use mpi_global

  implicit none
      

!
!    
!     Calculate DOU interactions 
!
!
! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  real(8),intent(in):: rcutdouo         ! DOU cutoff length [non-d] for O-Au
  real(8),intent(in):: rcutindouo       ! DOU cutin length [non-d] for O-Au
  real(8),intent(in):: rcutdouh         ! DOU cutoff length [non-d] for H-Au

!   INPUT & OUTPUT

  real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc)

  real(8),intent(inout):: for_viri_dou(:,:)  ! virial(DOU force) of each atom
  real(8),intent(inout):: pot_viri_dou       ! virial(DOU potential)
  real(8),intent(inout):: pot_virit_dou(:,:) ! virial tensor (DOU)

!   OUTPUT

  real(8),intent(out):: pot_dou         ! DOU potential

! LOCAL:
  real(8):: fdou(3,natom)    ! DOU force calculated here
!      real*8:: fdouall(3,natom) ! DOU force calculated here
!MPI                                       local(in MPI sense) field
  real(8):: potij                       ! temporal potential due to IJ  
  real(8):: ftmp                        ! temporal force  
  real(8):: fij(3)                      ! force to be subtracted due to 
                                        ! excluded atoms

  real(8):: rij(3)                      ! rj-ri
  real(8):: rij_abs                     ! |Rij|   
  real(8):: rij_abs2                    ! |Rij|**2
!  real(8):: rij_abs_inv1     ! |Rij|**-1

  real(8):: rtmp1                       ! exp(-beta*(|rij| - r0))
  real(8):: rtmp2                       ! exp(-2*beta*(|rij| - r0))

  real(8):: roff_r2                     ! r_off**2 - |rij|**2
  real(8):: r2_ron                      ! |rij|**2 - r_on**2

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

  real(8):: sqcut                       ! rcutdoux * rcutdoux
  real(8):: sqcutin                     ! rcutindouo * rcutindouo

  real(8):: denom_rcut                  ! 1/(r_off**2 - r_on**2)**3

  real(8):: S2                          ! switching function

  integer:: itype, jtype                ! VDW atom types
  integer:: i,j,n                       ! do loop index
  integer:: j1, j2                      ! limits of do loop
  integer:: jj                          ! atom index of the interacting atom

  integer:: dou1

!      integer:: im,i1,i2,ii,jm  ! loop index of molecule

!      logical:: exclatom        ! flag for if atom belongs list_excl

!      real*8:: xij,xij0
!      real*8:: yij,yij0
!      real*8:: zij,zij0
!      real*8:: r_ij

!      integer:: im,ip,jm,jp,km,kp
!      real*8:: xm,xp,ym,yp,zm,zp

!      real*8:: fxct,fyct,fzct,pott

!     +     +     +     +     +     +     +


! --- ALLOCATE ARRAY FOR FORCE CALCULATED IN THIS SUBROUTINE --- 

  fdou(1:3,1:natom) = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

! --- CALCULATE NONBONDED FORCE ---

  pot_dou = 0.0d0

  pot_virit_dou(1:3,1:3) = 0.0d0

  pot_viri_dou   = 0.0d0

!
! - DOU (O) interaction
!

  sqcut = rcutdouo * rcutdouo
  sqcutin = rcutindouo * rcutindouo

  denom_rcut = 1.0d0 / (sqcut - sqcutin) ! = 1 / (r_off**2 - r_on**2)
  denom_rcut = denom_rcut * denom_rcut * denom_rcut
                                        ! = 1 / (r_off**2 - r_on**2)**3

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, NATOM - 1
  looplast = ndou - 1

  DO dou1 = loopinit, looplast, loopstep

     j1 = douo_indexall(dou1)              ! the first and
     j2 = douo_indexall(dou1+loopstep) - 1 ! the last atoms interacting
                                           ! with i

     i = ndoulist(dou1)
     itype = atmindex(i)                ! atmtype of atom(i)

!     im = irmolept_list(i)  ! what molecule atom(i) belongs in

     IF (douo_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = douo_listall(j)            ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

!        jm = irmolept_list(jj) ! what molecule atom(i) belongs in

!       --- CALCULATE VIRIAL & FORCE ---

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2
                                        ! |rij|^2
 
!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs       = dsqrt(rij_abs2) ! |rij| 

!        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

        rtmp1 = dexp(-para_dou_beta_o*(rij_abs - para_dou_rad_o))
                                         ! exp(-beta*(|rij|-r0))
        rtmp2 = rtmp1 * rtmp1            ! exp(-2*beta*(|rij|-r0))

        ftmp  = 2.0d0 * para_dou_beta_o   &
            & * para_dou_wel_o   &
            & * (rtmp2 - rtmp1) / rij_abs
                                         ! = 2*beta*D*(rtmp2-rtmp1) / |rij|

        potij = para_dou_wel_o * (rtmp2 - 2.0d0*rtmp1)
                                         ! = D * (rtmp2 - 2*rtmp1)

        if (rij_abs2 > sqcutin) then

           roff_r2 = sqcut - rij_abs2    ! r_off**2 - |rij|**2
           r2_ron  = rij_abs2 - sqcutin  ! |rij|**2 - r_on**2

           S2 = roff_r2*roff_r2 * (roff_r2 + 3.0d0*r2_ron)
                              ! = (r_off**2 - |rij|**2)**2 
                              ! * ((roff**2 - |rij|**2) + 3*(|rij|**2 - ron**2))

           ftmp  = ftmp * S2   &
               & + potij * 12.0d0 * roff_r2 * r2_ron
           ftmp  = ftmp * denom_rcut
               
           potij = potij * S2 * denom_rcut

        end if

        pot_dou      = pot_dou + potij

        fij(1:3)     = ftmp * rij(1:3)

        fdou(1:3,i)  = fdou(1:3,i)  + fij(1:3)
        fdou(1:3,jj) = fdou(1:3,jj) - fij(1:3)

!       --- calculate virial tensor ---
        do n=1,3
           pot_virit_dou(1:3,n) = pot_virit_dou(1:3,n)   &
                &               + fij(1:3)*rij(n)
        end do
               
        pot_viri_dou = pot_viri_dou   &
             &       + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

     END DO

  END DO

!
! - DOU (H) interaction
!

  sqcut = rcutdouh * rcutdouh

  looplast = ndou - 1

  DO dou1 = loopinit, looplast, loopstep

     j1 = douh_indexall(dou1)              ! the first and
     j2 = douh_indexall(dou1+loopstep) - 1 ! the last atoms interacting
                                           ! with i

     i = ndoulist(dou1)
     itype = atmindex(i)                 ! atmtype of atom(i)

!     im = irmolept_list(i)  ! what molecule atom(i) belongs in

     IF (douh_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = douh_listall(j)             ! atom index of interacting atom

        jtype = atmindex(jj)             ! atmtype of atom(jj)

!        jm = irmolept_list(jj) ! what molecule atom(i) belongs in

!       --- CALCULATE VIRIAL & FORCE ---

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2 ! |rij|^2
 
!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs       = dsqrt(rij_abs2)  ! |rij|

!        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|
            
        potij = para_dou_gamma * para_dou_wel_h   &
            & * dexp(-2.0d0 * para_dou_beta_h   &
            & * (rij_abs - para_dou_rad_h))
                                         ! = gamma*D * exp(-2*beta*(|rij|-r0))

        ftmp = 2.0d0 * para_dou_beta_h * potij / rij_abs
                                    ! = 2*beta*gamma*D * exp(-2*beta*(|rij|-r0))
                                    !   / |rij|

        pot_dou      = pot_dou + potij

        fij(1:3)     = ftmp * rij(1:3)

        fdou(1:3,i)  = fdou(1:3,i)  + fij(1:3)
        fdou(1:3,jj) = fdou(1:3,jj) - fij(1:3)

!       --- calculate virial tensor ---
        do n=1,3
           pot_virit_dou(1:3,n) = pot_virit_dou(1:3,n)   &
                &               + fij(1:3)*rij(n)
        end do
               
        pot_viri_dou = pot_viri_dou   &
             &       + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

     END DO

  END DO

! --- ADD FDOU TO FORCE ---

  force(1:3,1:natom) = force(1:3,1:natom) + fdou(1:3,1:natom)
  for_viri_dou(1:3,1:natom) = for_viri_dou(1:3,1:natom) + fdou(1:3,1:natom)

!     +     +     +     +     +     +     +

end subroutine doupotp
