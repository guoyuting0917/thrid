!**************************************
!*  rfhpot.f Ver.1.2 '10.06.28        *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine rfhpot(xcel,ycel,zcel,   &
     &            rcutrfhfo,rcutrfhoo,rcutrfhoh,   &
     &            force,pot_rfh)

  use md_global
  use mpi_global

  implicit none


!
!    
!     Calculate RFH interactions
!       for iron oxide materials proposed by Rustad, Felmy, and Hay.
!       ref. Rustad, J.R., Felmy, A.R., Hay, B.P., Geochim. Cosmochim. Acta,
!              60 (1996), pp. 1553-1562.
!
!
! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  real(8),intent(in):: rcutrfhfo        ! RFH(FO) cutoff length [non-d]
  real(8),intent(in):: rcutrfhoo        ! RFH(OO) cutoff length [non-d]
  real(8),intent(in):: rcutrfhoh        ! RFH(OH) cutoff length [non-d]

!  INPUT & OUTPUT

  real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc) 

!  OUTPUT

  real(8),intent(out):: pot_rfh         ! RFH potential

! LOCAL:
  real(8):: frfh(3,natom)               ! RFH force calculated here
!  real(8):: frfhall(3,natom) ! RFH force calculated here
!MPI                                       local(in MPI sense) field
  real(8):: potij                       ! temporal potential due to IJ
  real(8):: ftmp                        ! temporal force
  real(8):: fij(3)                      ! force to be subtracted due to
                                        ! excluded atoms

  real(8):: rij(3)                      ! rj-ri
  real(8):: rij_abs                     ! |Rij|   
  real(8):: rij_abs2                    ! |Rij|**2
  real(8):: rij_abs_inv1                ! |Rij|**-1
  real(8):: rij_abs_inv2                ! |Rij|**-2

  real(8):: rfhfo_exp                   ! = a*exp(-b*|rij|)
  real(8):: rtmp6                       ! = c*|Rij|**-6
  real(8):: rtmp12                      ! = d*|Rij|**-12

  real(8):: rfhoo_exp1                  ! = a*exp(-b*|rij|)
  real(8):: rfhoo_exp2                  ! = exp(c*(|rij|-d))
  real(8):: rfhoo_exp3                  ! = 1/[1+exp(c*(|rij|-d))]
  real(8):: rfhoo_enecoeff = 3.3192d+2  ! = C0 [kcal/mol]
  real(8):: rfhoo_rinv12                ! = lja*|Rij|**-12
  real(8):: rfhoo_rinv6                 ! = ljb*|Rij|**-6

  real(8):: rfhoh_r_req                 ! = |rij| - req
  real(8):: rfhoh_exp1                  ! = a*exp(-b*|rij|) / |rij|
  real(8):: rfhoh_exp2                  ! = exp(-e*(|rij| - req)^2)
  real(8):: rfhoh_c_req                 ! = c*(|rij| - req)
  real(8):: rfhoh_e_req                 ! = e*(|rij| - req)
  real(8):: rfhoh_c_e                   ! = c*(|rij| - req)^2 - d*(|rij| - req)

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

  real(8):: sqcut                       ! rcutrshxx * rcutrshxx

  integer:: itype, jtype                ! VDW atom types
  integer:: i,j                         ! do loop index
  integer:: j1, j2                      ! limits of do loop
  integer:: jj                          ! atom index of the interacting atom

  integer:: rfh1

!      real*8:: xij,xij0
!      real*8:: yij,yij0
!      real*8:: zij,zij0
!      real*8:: r_ij

!      integer:: im,ip,jm,jp,km,kp
!      real*8:: xm,xp,ym,yp,zm,zp

!      real*8:: fxct,fyct,fzct,pott

!     +     +     +     +     +     +     +


! --- ALLOCATE ARRAY FOR FORCE CALCULATED IN THIS SUBROUTINE --- 

  frfh(1:3,1:natom) = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

! --- CALCULATE NONBONDED FORCE ---

  pot_rfh = 0.0d0

!
! - RSH (F-O) interaction
!

  sqcut = rcutrfhfo * rcutrfhfo

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, NATOM - 1
  looplast = nrfh - 1

  DO rfh1 = loopinit, looplast, loopstep

     j1 = rfhfo_indexall(rfh1)              ! the first and
     j2 = rfhfo_indexall(rfh1+loopstep) - 1 ! the last atoms interacting
                                            ! with i

     i = nrfhlist(rfh1)
     itype = atmindex(i)                ! atmtype of atom(i)

     IF (rfhfo_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = rfhfo_listall(j)           ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2 ! |rij|^2
!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs       = dsqrt(rij_abs2) ! |rij|

        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|
        rij_abs_inv2  = rij_abs_inv1 * rij_abs_inv1 ! 1/|rij|^2

        rfhfo_exp = para_rfhfo_a * dexp(-para_rfhfo_b * rij_abs) 
                                        ! A exp(-B|rij|)

        rtmp6 = rij_abs_inv2 * rij_abs_inv2 * rij_abs_inv2 ! 1/|rij|^6
        rtmp12 = rtmp6 * rtmp6          ! 1/|rij|^12
            
        rtmp6 = rtmp6 * para_rfhfo_c    ! C/|rij|^6
        rtmp12 = rtmp12 * para_rfhfo_d  ! D/|rij|^12

#if defined(_RFH_FOINT_ORG)
        potij = rfhfo_exp + rtmp6 + rtmp12
#else
        potij = rfhfo_exp + rtmp6 - rtmp12
#endif
                                       ! A exp(-B|rij|) + C/|rij|^6 + D/|rij|^12

#if defined(_RFH_FOINT_ORG)
        ftmp = para_rfhfo_b*rfhfo_exp*rij_abs_inv1   &
           & + (rtmp6 * 6.0d0 + rtmp12 * 12.0d0) * rij_abs_inv2
#else
        ftmp = para_rfhfo_b*rfhfo_exp*rij_abs_inv1   &
           & + (rtmp6 * 6.0d0 - rtmp12 * 12.0d0) * rij_abs_inv2
#endif
                                        ! BA exp(-B|rij|)/|rij| + 6C/|rij|^8
                                        !                       + 12D/|rij|^14

        pot_rfh    = pot_rfh + potij

        fij(1:3)     = ftmp * rij(1:3)

        frfh(1:3,i)  = frfh(1:3,i)  + fij(1:3)
        frfh(1:3,jj) = frfh(1:3,jj) - fij(1:3)

     END DO
  END DO

!
! - RSH (O-O) interaction
!

  sqcut = rcutrfhoo * rcutrfhoo

  looplast = nrfh - 1

  DO rfh1 = loopinit, looplast, loopstep

     j1 = rfhoo_indexall(rfh1)              ! the first and
     j2 = rfhoo_indexall(rfh1+loopstep) - 1 ! the last atoms interacting
                                            ! with i

     i = nrfhlist(rfh1)
     itype = atmindex(i)                ! atmtype of atom(i)

     IF (rfhoo_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = rfhoo_listall(j)           ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2   ! |rij|^2

!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

#if defined(_RFH_OOINT_ORG) || 1
        rij_abs       = dsqrt(rij_abs2) ! |rij|

        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|
!        rij_abs_inv2  = rij_abs_inv1 * rij_abs_inv1 ! 1/|rij|^2

        rfhoo_exp1 = para_rfhoo_a * dexp(-para_rfhoo_b * rij_abs)
                                        ! A exp(-B|rij|)

        rfhoo_exp2 = dexp(para_rfhoo_c * (rij_abs - para_rfhoo_d))
                                        ! exp(C*(|rij|-D))

        rfhoo_exp3 = 1.0d0 / (1.0d0 + rfhoo_exp2)
                                        ! 1/[1 + exp(C*(|rij|-D))]
            
        potij = rfhoo_exp1 + rfhoo_exp3 * rfhoo_enecoeff
                                  ! A exp(-B|rij|) + 1/[1 + exp(C*(|rij|-D))]*C0

        ftmp = (para_rfhoo_b * rfhoo_exp1   &
           &  + rfhoo_enecoeff * para_rfhoo_c * rfhoo_exp2   &
           &  * rfhoo_exp3 * rfhoo_exp3)  &
           & * rij_abs_inv1
                               ! 1/r * {BA exp(-B|rij|)
                               ! + C0 C exp(C*(|rij|-D))/[1+exp(C*(|rij|-D))]^2}

#else
        rij_abs_inv2 = 1.0d0 / rij_abs2 ! |rij|^-2
            
        rfhoo_rinv6 = rij_abs_inv2 * rij_abs_inv2 * rij_abs_inv2   ! |rij|^-6
        rfhoo_rinv12 = rfhoo_rinv6 * rfhoo_rinv6   ! |rij|^-12
            
        rfhoo_rinv12 = rfhoo_rinv12 * para_rfhoo_lja   ! A / |rij|^12
        rfhoo_rinv6 = rfhoo_rinv6 * para_rfhoo_ljb   ! B / |rij|^6

        potij = rfhoo_rinv12 + rfhoo_rinv6
                                        ! A / |rij|^12 + B / |rij|^6

        ftmp = (12.0d0 * rfhoo_rinv12 + 6.0d0 * rfhoo_rinv6)   &
           & * rij_abs_inv2             ! 12A / |rij|^14 + 6B / |rij|^8
#endif

        pot_rfh    = pot_rfh + potij

        fij(1:3)     = ftmp * rij(1:3)

        frfh(1:3,i)  = frfh(1:3,i)  + fij(1:3)
        frfh(1:3,jj) = frfh(1:3,jj) - fij(1:3)

     END DO
  END DO

!
! - RSH (O-H) interaction
!

  sqcut = rcutrfhoh * rcutrfhoh

  looplast = nrfh - 1

  DO rfh1 = loopinit, looplast, loopstep

     j1 = rfhoh_indexall(rfh1)              ! the first and
     j2 = rfhoh_indexall(rfh1+loopstep) - 1 ! the last atoms interacting
                                            ! with i

     i = nrfhlist(rfh1)
     itype = atmindex(i)                ! atmtype of atom(i)

     IF (rfhoh_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = rfhoh_listall(j)           ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2   ! |rij|^2

!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs       = dsqrt(rij_abs2) ! |rij|

        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|
        rij_abs_inv2  = rij_abs_inv1 * rij_abs_inv1   ! 1/|rij|^2

        rfhoh_r_req = rij_abs - para_rfhoh_req   ! |rij| - req

        rfhoh_exp1 = para_rfhoh_a * dexp(-para_rfhoh_b * rij_abs)   &
             &     * rij_abs_inv1
                                        ! a*exp(-b*|rij|) / |rij|

        rfhoh_c_req = para_rfhoh_c * rfhoh_r_req   ! c*(|rij| - req)
        rfhoh_e_req = para_rfhoh_e * rfhoh_r_req   ! e*(|rij| - req)

        rfhoh_exp2 = dexp(-rfhoh_e_req * rfhoh_r_req)
                                        ! exp(-e*(|rij| - req)^2)

        rfhoh_c_e = rfhoh_c_req * rfhoh_r_req   &
             &    - para_rfhoh_d * rfhoh_r_req
                                        ! c*(|rij| - req)^2 - d*(|rij| - req)

        potij = rfhoh_exp1 + rfhoh_c_e * rfhoh_exp2
                                       ! a*exp(-b*|rij|) / |rij|
                                       ! + [c*(|rij| - req)^2 - d*(|rij| - req)]
                                       ! * exp(-e*(|rij| - req)^2)

        ftmp = (rfhoh_exp1 * (para_rfhoh_b + rij_abs_inv1)   &
           &  + (2.0d0*rfhoh_e_req * rfhoh_c_e   &
           &   + para_rfhoh_d - 2.0d0*rfhoh_c_req) * rfhoh_exp2)   &
           & * rij_abs_inv1
                                     ! a*exp(-b*|rij|) *(b + 1/|rij|) / |rij|^2
                                     ! + [ 2e*(|rij|-req) 
                                     !   * {c*(|rij| - req)^2 - d*(|rij| - req)}
                                     !   + d - 2c*(|rij| - req)]
                                     !  * exp(-e*(|rij| - req)^2) / |rij|
            
        pot_rfh    = pot_rfh + potij

        fij(1:3)     = ftmp * rij(1:3)

        frfh(1:3,i)  = frfh(1:3,i)  + fij(1:3)
        frfh(1:3,jj) = frfh(1:3,jj) - fij(1:3)

     END DO
  END DO

! --- ADD FSH TO FORCE ---

  force(1:3,1:natom) = force(1:3,1:natom) + frfh(1:3,1:natom)

!     write(20,'(6f12.3)') frfhall(:,:)

!     +     +     +     +     +     +     +

end subroutine rfhpot

!----------------------------------------------------------
subroutine rfhpotp(xcel,ycel,zcel,   &
     &             rcutrfhfo,rcutrfhoo,rcutrfhoh,   &
     &             force,pot_rfh,   &
     &             for_viri_rfh,pot_viri_rfh,   &
     &             pot_virit_rfh)

  use md_global
  use mpi_global
  
  implicit none
      

!
!    
!     Calculate RFH interactions 
!
!
! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  real(8),intent(in):: rcutrfhfo        ! RFH(FO) cutoff length [non-d]
  real(8),intent(in):: rcutrfhoo        ! RFH(OO) cutoff length [non-d]
  real(8),intent(in):: rcutrfhoh        ! RFH(OH) cutoff length [non-d]

!   INPUT & OUTPUT

  real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc)

  real(8),intent(inout):: for_viri_rfh(:,:)  ! virial(RFH force) of each atom
  real(8),intent(inout):: pot_viri_rfh       ! virial(RFH potential)
  real(8),intent(inout):: pot_virit_rfh(:,:) ! virial tensor (RFH)

!   OUTPUT

  real(8),intent(out):: pot_rfh         ! RFH potential

! LOCAL:
  real(8):: frfh(3,natom)               ! RFH force calculated here
!  real(8):: frfhall(3,natom) ! RFH force calculated here
!MPI                                       local(in MPI sense) field
  real(8):: potij                       ! temporal potential due to IJ  
  real(8):: ftmp                        ! temporal force  
  real(8):: fij(3)                      ! force to be subtracted due to 
                                        ! excluded atoms

  real(8):: rij(3)                      ! rj-ri
  real(8):: rij_abs                     ! |Rij|   
  real(8):: rij_abs2                    ! |Rij|**2
  real(8):: rij_abs_inv1                ! |Rij|**-1
  real(8):: rij_abs_inv2                ! |Rij|**-2

  real(8):: rfhfo_exp                   ! = a*exp(-b*|rij|)
  real(8):: rtmp6                       ! = c*|Rij|**-6
  real(8):: rtmp12                      ! = d*|Rij|**-12

  real(8):: rfhoo_exp1                  ! = a*exp(-b*|rij|)
  real(8):: rfhoo_exp2                  ! = exp(c*(|rij|-d))
  real(8):: rfhoo_exp3                  ! = 1/[1+exp(c*(|rij|-d))]
  real(8):: rfhoo_enecoeff = 3.3192d+2  ! = C0 [kcal/mol]
  real(8):: rfhoo_rinv12                ! = lja*|Rij|**-12
  real(8):: rfhoo_rinv6                 ! = ljb*|Rij|**-6

  real(8):: rfhoh_r_req                 ! = |rij| - req
  real(8):: rfhoh_exp1                  ! = a*exp(-b*|rij|) / |rij|
  real(8):: rfhoh_exp2                  ! = exp(-e*(|rij| - req)^2)
  real(8):: rfhoh_c_req                 ! = c*(|rij| - req)
  real(8):: rfhoh_e_req                 ! = e*(|rij| - req)
  real(8):: rfhoh_c_e                   ! = c*(|rij| - req)^2 - d*(|rij| - req)

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

  real(8):: sqcut                       ! rcutrfhxx * rcutrfhxx

  integer:: itype, jtype                ! VDW atom types
  integer:: i,j,n                       ! do loop index
  integer:: j1, j2                      ! limits of do loop
  integer:: jj                          ! atom index of the interacting atom

  integer:: rfh1

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

  frfh(1:3,1:natom) = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

! --- CALCULATE NONBONDED FORCE ---

  pot_rfh = 0.0d0

  pot_virit_rfh(1:3,1:3) = 0.0d0

  pot_viri_rfh = 0.0d0

!
! - RSH (F-O) interaction
!

  sqcut = rcutrfhfo * rcutrfhfo

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, NATOM - 1
  looplast = nrfh - 1

  DO rfh1 = loopinit, looplast, loopstep

     j1 = rfhfo_indexall(rfh1)              ! the first and
     j2 = rfhfo_indexall(rfh1+loopstep) - 1 ! the last atoms interacting
                                            ! with i

     i = nrfhlist(rfh1)
     itype = atmindex(i)                ! atmtype of atom(i)

!     im = irmolept_list(i)  ! what molecule atom(i) belongs in

     IF (rfhfo_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = rfhfo_listall(j)           ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

!        jm = irmolept_list(jj) ! what molecule atom(i) belongs in

!       --- CALCULATE VIRIAL & FORCE ---

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2 ! |rij|^2
 
!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs       = dsqrt(rij_abs2) ! |rij| 

        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|
        rij_abs_inv2  = rij_abs_inv1 * rij_abs_inv1 ! 1/|rij|^2

        rfhfo_exp = para_rfhfo_a * dexp(-para_rfhfo_b * rij_abs) 
                                        ! A exp(-B|rij|)

        rtmp6 = rij_abs_inv2 * rij_abs_inv2 * rij_abs_inv2 ! 1/|rij|^6
        rtmp12 = rtmp6 * rtmp6          ! 1/|rij|^12
            
        rtmp6 = rtmp6 * para_rfhfo_c    ! C/|rij|^6
        rtmp12 = rtmp12 * para_rfhfo_d  ! D/|rij|^12

#if defined(_RFH_FOINT_ORG)
        potij = rfhfo_exp + rtmp6 + rtmp12
#else
        potij = rfhfo_exp + rtmp6 - rtmp12
#endif

#if defined(_RFH_FOINT_ORG)
        ftmp = para_rfhfo_b*rfhfo_exp*rij_abs_inv1   &
           & + (rtmp6 * 6.0d0 + rtmp12 * 12.0d0) * rij_abs_inv2
#else
        ftmp = para_rfhfo_b*rfhfo_exp*rij_abs_inv1   &
           & + (rtmp6 * 6.0d0 - rtmp12 * 12.0d0) * rij_abs_inv2
#endif
                                        ! BA exp(-B|rij|)/|rij| + 6C/|rij|^8
                                        !                       + 12D/|rij|^14

        pot_rfh    = pot_rfh + potij

        fij(1:3)     = ftmp * rij(1:3)

        frfh(1:3,i)  = frfh(1:3,i)  + fij(1:3)
        frfh(1:3,jj) = frfh(1:3,jj) - fij(1:3)

!       --- calculate virial tensor ---
        do n=1,3
           pot_virit_rfh(1:3,n) = pot_virit_rfh(1:3,n)   &
                &               + fij(1:3)*rij(n)
        end do
               
        pot_viri_rfh = pot_viri_rfh   &
             &       + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

     END DO

  END DO

!
! - RSH (O-O) interaction
!

  sqcut = rcutrfhoo * rcutrfhoo

  looplast = nrfh - 1

  DO rfh1 = loopinit, looplast, loopstep

     j1 = rfhoo_indexall(rfh1)              ! the first and
     j2 = rfhoo_indexall(rfh1+loopstep) - 1 ! the last atoms interacting
                                            ! with i

     i = nrfhlist(rfh1)
     itype = atmindex(i)                ! atmtype of atom(i)

!     im = irmolept_list(i)  ! what molecule atom(i) belongs in

     IF (rfhoo_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = rfhoo_listall(j)           ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

!        jm = irmolept_list(jj) ! what molecule atom(i) belongs in

!       --- CALCULATE VIRIAL & FORCE ---

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2 ! |rij|^2
 
!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

#if defined(_RFH_OOINT_ORG) || 1
        rij_abs       = dsqrt(rij_abs2) ! |rij| 

        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|
!        rij_abs_inv2  = rij_abs_inv1 * rij_abs_inv1 ! 1/|rij|^2

        rfhoo_exp1 = para_rfhoo_a * dexp(-para_rfhoo_b * rij_abs)
                                        ! A exp(-B|rij|)

        rfhoo_exp2 = dexp(para_rfhoo_c * (rij_abs - para_rfhoo_d))
                                        ! exp(C*(|rij|-D))

        rfhoo_exp3 = 1.0d0 / (1.0d0 + rfhoo_exp2)
                                        ! 1/[1 + exp(C*(|rij|-D))]
            
        potij = rfhoo_exp1 + rfhoo_exp3 * rfhoo_enecoeff
                                 ! A exp(-B|rij|) + 1/[1 + exp(C*(|rij|-D))]*C0

        ftmp = (para_rfhoo_b * rfhoo_exp1   &
     &        + rfhoo_enecoeff * para_rfhoo_c * rfhoo_exp2   &
     &        * rfhoo_exp3 * rfhoo_exp3)   &
     &       * rij_abs_inv1
                               ! 1/r * {BA exp(-B|rij|)
                               ! + C0 C exp(C*(|rij|-D))/[1+exp(C*(|rij|-D))]^2}
#else
        rij_abs_inv2 = 1.0d0 / rij_abs2 ! |rij|^-2
            
        rfhoo_rinv6 = rij_abs_inv2 * rij_abs_inv2 * rij_abs_inv2 ! |rij|^-6
        rfhoo_rinv12 = rfhoo_rinv6 * rfhoo_rinv6 ! |rij|^-12
            
        rfhoo_rinv12 = rfhoo_rinv12 * para_rfhoo_lja ! A / |rij|^12
        rfhoo_rinv6 = rfhoo_rinv6 * para_rfhoo_ljb ! B / |rij|^6

        potij = rfhoo_rinv12 + rfhoo_rinv6
                                        ! A / |rij|^12 + B / |rij|^6

        ftmp = (12.0d0 * rfhoo_rinv12 + 6.0d0 * rfhoo_rinv6)   &
           & * rij_abs_inv2             ! 12A / |rij|^14 + 6B / |rij|^8
#endif

        pot_rfh    = pot_rfh + potij

        fij(1:3)   = ftmp * rij(1:3)

        frfh(1:3,i)  = frfh(1:3,i)  + fij(1:3)
        frfh(1:3,jj) = frfh(1:3,jj) - fij(1:3)

!       --- calculate virial tensor ---
        do n=1,3
           pot_virit_rfh(1:3,n) = pot_virit_rfh(1:3,n)   &
                &               + fij(1:3)*rij(n)
        end do
               
        pot_viri_rfh = pot_viri_rfh   &
             &       + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

     END DO

  END DO

!
! - RSH (O-H) interaction
!

  sqcut = rcutrfhoh * rcutrfhoh

  looplast = nrfh - 1

  DO rfh1 = loopinit, looplast, loopstep

     j1 = rfhoh_indexall(rfh1)              ! the first and
     j2 = rfhoh_indexall(rfh1+loopstep) - 1 ! the last atoms interacting
                                            ! with i

     i = nrfhlist(rfh1)
     itype = atmindex(i)                ! atmtype of atom(i)

     IF (rfhoh_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = rfhoh_listall(j)           ! atom index of interacting atom

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

        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|
        rij_abs_inv2  = rij_abs_inv1 * rij_abs_inv1 ! 1/|rij|^2

        rfhoh_r_req = rij_abs - para_rfhoh_req ! |rij| - req

        rfhoh_exp1 = para_rfhoh_a * dexp(-para_rfhoh_b * rij_abs)   &
             &     * rij_abs_inv1
                                        ! a*exp(-b*|rij|) / |rij|

        rfhoh_c_req = para_rfhoh_c * rfhoh_r_req ! c*(|rij| - req)
        rfhoh_e_req = para_rfhoh_e * rfhoh_r_req ! e*(|rij| - req)

        rfhoh_exp2 = dexp(-rfhoh_e_req * rfhoh_r_req) ! exp(-e*(|rij| - req)^2)

        rfhoh_c_e = rfhoh_c_req * rfhoh_r_req   &
             &    - para_rfhoh_d * rfhoh_r_req
                                        ! c*(|rij| - req)^2 - d*(|rij| - req)

        potij = rfhoh_exp1 + rfhoh_c_e * rfhoh_exp2
                                       ! a*exp(-b*|rij|) / |rij|
                                       ! + [c*(|rij| - req)^2 - d*(|rij| - req)]
                                       ! * exp(-e*(|rij| - req)^2)

        ftmp = (rfhoh_exp1 * (para_rfhoh_b + rij_abs_inv1)   &
           &  + (2.0d0*rfhoh_e_req * rfhoh_c_e   &
           &   + para_rfhoh_d - 2.0d0*rfhoh_c_req) * rfhoh_exp2)   &
           & * rij_abs_inv1
                                     ! a*exp(-b*|rij|) *(b + 1/|rij|) / |rij|^2
                                     ! + [ 2e*(|rij|-req) 
                                     !   * {c*(|rij| - req)^2 - d*(|rij| - req)}
                                     !   + d - 2c*(|rij| - req)]
                                     !  * exp(-e*(|rij| - req)^2) / |rij|

        pot_rfh    = pot_rfh + potij

        fij(1:3)   = ftmp * rij(1:3)

        frfh(1:3,i)  = frfh(1:3,i)  + fij(1:3)
        frfh(1:3,jj) = frfh(1:3,jj) - fij(1:3)

!       --- calculate virial tensor ---
        do n=1,3
           pot_virit_rfh(1:3,n) = pot_virit_rfh(1:3,n)   &
                &               + fij(1:3)*rij(n)
        end do
               
        pot_viri_rfh = pot_viri_rfh   &
             &       + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

     END DO

  END DO

! --- ADD FRFH TO FORCE ---

  force(1:3,1:natom) = force(1:3,1:natom) + frfh(1:3,1:natom)
  for_viri_rfh(1:3,1:natom) = for_viri_rfh(1:3,1:natom) + frfh(1:3,1:natom)

!     +     +     +     +     +     +     +

end subroutine rfhpotp
