!*****************************
!*  rattle_mtk.f90 Ver.1.6   *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-08-24 12:10:17 gota>

subroutine rattle_cmtk(atmcor_old, &
    &                  eps_rattle, dt_short_cal, dt_long_cal, &
    &                  istep_short, nstep_short, &
    &                  rv, rvv, &
    &                  pcont_axis)

!    subroutine to constrain bond-lengh
!      by RATTLE (Andersen, H. C., 1983, J. Comput. Phys. 52, 24-34.)
!
!      'rattle_C' is for Coordinate resetting
!       for the first half of RATTLE.
!
!      All the bonds are constraint within a precision of
!      eps_rattle (A).
!      Namely, coord(:,:) will be moved so that the following
!      equation be satisfied for all constraint atom pairs (i,j):
!
!         [|coord(:,i)-coord(:,j)| - dij]/dij < eps_rattle .
!
!
!      Thus the resultant coord(:,:) will contain the constraint
!      coordinates at time = T+dT.
!
!
!      INPUT
!        coord(:,:)     :  coordinate at T+dT before constraint
!        veloc(:,:)     :  velocity at T+dT/2 before constraint
!      OUTPUT
!        coord(:,:)     :  coordinate at T+dT after constraint
!        veloc(:,:)     :  velocity at T+dT/2 after constraint
!

  use md_global

  implicit none

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

! LOCAL:
  integer:: maxiter = 300
                                ! maximum allowable
                                ! number of iterations

  integer:: iter            ! iteration index

  integer:: ib, jb          ! atom index making bonds


  real(8):: dij              ! the bond length between ib & jb
                                ! to be satisfied.
  real(8):: dij2             ! dij**2

  real(8):: Rij(3)           ! = coord(:,i)-coord(:,j)
  real(8):: Rij2             ! Rij**2
  real(8):: Rij_old(3)       ! = coord_old(:,i)-coord_old(:,j)

  real(8):: amassi_inv       ! 1/atmmass(i)
  real(8):: amassj_inv       ! 1/atmmass(j)

  real(8):: gij, GRij(3), RR, delta
                                ! constraint force etc.

  logical:: mvlast(maxnatom)
                                ! set true if the coordinate of atom(i)
                                ! was moved in the last iteration.
                                ! * if both ib & jb were not moved in the last
                                !   iteration,
                                !   then they need not be
                                !   moved in this iteration.

  logical:: mvnow(maxnatom)
                                ! set true if the coordinate of atom(i)
                                ! was moved in the current iteration.
                                ! this information will be passed to
                                ! mv last at the end of the iteration.

  logical:: finished = .false.  ! set true if no more coordinate resetting
                                ! is needed.

  real(8):: one = 1.0d0
  real(8):: two = 2.0d0

  real(8):: eps_rattle2         ! = eps_rattle^2

  integer:: i, m                ! do loop indexes
  integer:: i1,i2,ii
  integer:: nc

  real(8):: dt_cal              ! time step

  real(8):: rv_inv(3)           ! = 1/rv or 1/rvv(:)

!     +     +     +     +     +     +     +     +     +     +

!     --- initialization ---
  eps_rattle2 = eps_rattle * eps_rattle

  ! set rvinv
  if (pcont_axis == 'iso') then
      rv_inv(1:3) = 1.0d0 / rv
  else
      rv_inv(1:3) = 1.0d0 / rvv(1:3)
  end if

!     --- Set timestep ---

  if (istep_short == 1) then
     dt_cal = dt_short_cal
  else
     dt_cal = dt_short_cal
  end if

!     --- SET FLAGS BEFORE ITERATION-LOOP ETC ---

  do m=1,natom
     mvlast(m) = .true.     ! before iteration,
                                ! no atoms were updated.
     mvnow(m)  = .true.
  end do

!     --- GRAND LOOP TO APPLY BOND LENGHT CONSTRAINT ---

!      * this loop is iterated until all the constraints
!        are satisfied.

!     -- LOOP OVER CONSTRAINTS  --

  DOconstr: DO i = 1, nconst

     i1 = index_mconst(i)
     i2 = index_mconst(i+1) - 1

     DOiter: DO iter = 1, maxiter ! loop over iterations

!           -- set finished flag --

        finished = .true.   ! this flag will be set false
                                ! if at least one coordinate was reset
                                ! in the iteration

        do ii = i1, i2

           nc = list_mconst(ii)

!              - pick up constraint to examine -

           ib  = iconst(nc)  ! bond length between ib & jb
           jb  = jconst(nc)  ! will be constraint
           dij = dconst(nc)  ! to dij

!              - CYCLE DOcostr -
!              if neither ib nor jb have been updated last iteration
!                                          or in the current iteration.
!              * if both ib & jb were not moved in the last or current
!                iterations , then they need not be moved this iteration.

           IF (.not. (mvlast(ib) .or. mvlast(jb) .or. &
                & mvnow(ib) .or. mvnow(jb))) then
              cycle
           END IF

!              - SEE if Rij satisfy the constraint -

           Rij(1:3) = atmcor(1:3,ib) - atmcor(1:3,jb)

           Rij2   = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)

           dij2 = dij * dij

           delta  = dabs(Rij2-dij2)/dij2


           IF (delta < eps_rattle) CYCLE
!                 *  Rij satisfy the constraint

!              - NOW atoms ib & jb should be moved to satisfy the constraint -


           finished   = .false. ! DOinter is not finished
                                    ! in the current iteration.

           mvnow(ib)  = .true. ! ib & jb are moved now
           mvnow(jb)  = .true. !


!              - Coordinate update -

           Rij_old(1:3) = atmcor_old(1:3,ib) - atmcor_old(1:3,jb)
                                ! Rij at time T

           RR = Rij(1)*Rij_old(1) + Rij(2)*Rij_old(2) + &
                & Rij(3)*Rij_old(3)

           amassi_inv = one/atmmass(ib) ! 1/mi
           amassj_inv = one/atmmass(jb) ! 1/mj

!               dij2 = dij*dij   ! dij**2

           gij  = (Rij2-dij2)/(two * RR * (amassi_inv + amassj_inv) &
                & * dt_cal)      ! additional dividing timestep

           GRij(1:3) = gij*Rij_old(1:3)

           atmcor(1:3,ib) =  atmcor(1:3,ib) - amassi_inv * GRij(1:3) &
                &          * dt_short_cal
           atmcor(1:3,jb) =  atmcor(1:3,jb) + amassj_inv * GRij(1:3) &
                &          * dt_short_cal

!              - Velocity update -
!              divided by rv or rvv
           GRij(1:3) = GRij(1:3) * rv_inv(1:3)

           atmvel(1:3,ib) = atmvel(1:3,ib) - amassi_inv * GRij(1:3)
           atmvel(1:3,jb) = atmvel(1:3,jb) + amassj_inv * GRij(1:3)

        end do

!           -- EXIT IF no atom coordinate have been updated --

        if (finished) EXIT DOiter

!           -- Transfer mvnow to mvlast and go to next iteration -

!            do m = 1, natom
!               mvlast(m) = mvnow(m)
!               mvnow(m)  = .false.
!            end do
        mvlast(ib) = mvnow(ib)
        mvnow(ib)  = .false.
        mvlast(jb) = mvnow(jb)
        mvnow(jb)  = .false.

     END DO DOiter          ! end of loop over iterations

!        --- WRITE NUMBER OF ITERATIONS IF FAILED OR IF DEBUG MODE -

     if (.not. finished) then

!           - RATTLE failed to converge, so quitting -

        write(6,999)           maxiter
999     format(/3x, &
             & 'ERROR: RATTLE(C) failed to converge within MAXITER = ', &
             & i3)
        STOP
     end if

  END DO DOconstr

!     +     +     +     +     +     +     +     +     +     +

end subroutine rattle_cmtk

!-----------------------------------------------------------------------
subroutine rattlep_v(eps_rattle, &
    &                dt_long_cal, dt_short_cal, &
    &                istep_short, nstep_short, &
    &                expcoeff, &
    &                atm_viri_const, atm_virit_const)

!    subroutine to constrain veloctity
!      by RATTLE (Andersen, H. C., 1983, J. Comput. Phys. 52, 24-34.)
!
!      'rattle_V' is for Velocity resetting
!       for the second half of RATTLE.
!
!      All the velocities are constraint within a precision of
!       eps_rattle/dt (/sqrt(g/kcal)).
!      Namely, veloc(:,:) will be moved so that the following
!      equation be satisfied for all constraint atom pairs (i,j):
!
!      Rij.Vij/|Rij||Vij| < eps_rattle .
!
!      Thus the resultant veloc(:,:) will contain the constraint
!      velocities at time = T+dT.
!
!
!      INPUT
!        coord(:,:)     :  coordinate at T+dT after constraint
!        veloc(:,:)     :  velocity at T+dT before constraint
!      OUTPUT
!        veloc(:,:)     :  velocity at T+dT/2 after constraint
!

  use md_global

  implicit none

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

! LOCAL:
  integer:: maxiter = 300
                                ! maximum allowable
                                ! number of iterations

  integer:: iter            ! iteration index

  integer:: ib, jb          ! atom index making bonds

  real(8):: dij              ! the bond length between ib & jb
                                ! to be satisfied.

  real(8):: Rij(3)           ! = coord(:,i)-coord(:,j)
  real(8):: Vij(3)           ! = veloc(:,i)-veloc(:,j)
  real(8):: RVij             ! Rij.Vij
  real(8):: Rij2             ! |Rij|**2
  real(8):: Vij2             ! |Vij|**2

  real(8):: phi              ! (pi/2-phi) is the angle between
                                !               Vij & Rij
  real(8):: amassi_inv       ! 1/atmmass(i)
  real(8):: amassj_inv       ! 1/atmmass(j)

  real(8):: kij, KRij(3)
                                ! constraint force etc.

  logical:: mvlast(maxnatom)
                                ! set true if the coordinate of atom(i)
                                ! was moved in the last iteration.
                                ! * if both ib & jb were not moved in the last
                                !   iteration,
                                !   then they need not be moved
                                !   in this iteration.

  logical:: mvnow(maxnatom)
                                ! set true if the coordinate of atom(i)
                                ! was moved in the current iteration.
                                ! this information will be passed to
                                ! mv last at the end of the iteration.

  logical:: finished = .false. ! set true if no more coordinate resetting
                                ! is needed.

  real(8):: one = 1.0d0

  real(8):: vtimestep        ! set timestep for virial calculation
  real(8):: inv_vtime        ! = 1/vtimestep

  real(8):: eps_rattle2      ! = eps_rattle^2

  integer:: i,m,n           ! do loop indexes
  integer:: i1,i2,ii
  integer:: nc

!     +     +     +     +     +     +     +     +     +     +

!     --- initialization ---
  eps_rattle2 = eps_rattle * eps_rattle

!     --- Initilize virial for constraint force ---

  atm_viri_const = 0.0d0
  atm_virit_const(1:3,1:3) = 0.0d0

!     --- Set timestep for calculating virial of constraint force ---

  if (istep_short == nstep_short) then
     vtimestep = dt_long_cal*expcoeff
  else
     vtimestep = dt_short_cal
  end if

  inv_vtime = 1.0d0 / vtimestep

!     --- SET FLAGS BEFORE ITERATION-LOOP ETC ---

  do m=1,natom
     mvlast(m) = .true.     ! before iteration,
                                ! no atoms were updated.
     mvnow(m)  = .true.
  end do

!     --- GRAND LOOP TO APPLY VELOCITY CONSTRAINT ---

!      * this loop is iterated until all the constraints
!        are satisfied.

!     -- LOOP OVER CONSTRAINTS  --

  DOconstr: DO i = 1, nconst

     i1 = index_mconst(i)
     i2 = index_mconst(i+1) - 1

!        -- iteration for convergence --
     DOiter: DO iter = 1, maxiter ! loop over iterations

!           -- set finished flag --

        finished = .true.   ! this flag will be set false
                                ! if at least one coordinate was reset
                                ! in the iteration

        do ii = i1, i2

           nc = list_mconst(ii)

!              - pick up constraint to examine -

           ib  = iconst(nc) ! bond length between ib & jb
           jb  = jconst(nc) ! will be constraint
           dij = dconst(nc) ! to dij

!              - CYCLE DOcostr -
!              if neither ib nor jb have been updated last iteration
!                                          or in the current iteration.
!              * if both ib & jb were not moved in the last or current
!                iterations , then they need not be moved this iteration.

           IF (.not. (mvlast(ib) .or. mvlast(jb) .or. &
                & mvnow(ib) .or. mvnow(jb))) then
              cycle
           END IF

!              - SEE if (Rij.Vij) satisfy the constraint -

           Rij(1:3) = atmcor(1:3,ib) - atmcor(1:3,jb)
           Vij(1:3) = atmvel(1:3,ib) - atmvel(1:3,jb)

           RVij   = Rij(1)*Vij(1) + Rij(2)*Vij(2) + Rij(3)*Vij(3)
           Rij2   = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
           Vij2   = Vij(1)*Vij(1) + Vij(2)*Vij(2) + Vij(3)*Vij(3)

           phi    = RVij*RVij/Rij2*Vij2
                                ! (pi/2-phi) is the angle between
                                !               Vij & Rij
                                ! phi is squared for efficiency

           IF (phi < eps_rattle2) CYCLE
                                ! Vij satisfy the constraint

!              - NOW atoms ib & jb should be moved to satisfy the constraint -

           finished   = .false. ! DOinter is not finished
                                    ! in the current iteration.

           mvnow(ib)  = .true. ! ib & jb are moved now
           mvnow(jb)  = .true. !


!              - Velocity update -

           amassi_inv = one/atmmass(ib) ! 1/mi
           amassj_inv = one/atmmass(jb) ! 1/mj

           kij  = RVij/(dij*dij*(amassi_inv + amassj_inv))

           KRij(1:3) = kij* Rij(1:3)

           atmvel(1:3,ib) = atmvel(1:3,ib) - amassi_inv * KRij(1:3)
           atmvel(1:3,jb) = atmvel(1:3,jb) + amassj_inv * KRij(1:3)

#if defined(_DO_NOT_USE_THIS)
!!!        pressure calculation with constraint force does not work properly.
!          - Virial update -
           KRij(1:3) = 2.0d0*KRij(1:3)*inv_vtime

           atm_viri_const = atm_viri_const - KRij(1)*Rij(1) &
                &                          - KRij(2)*Rij(2) &
                &                          - KRij(3)*Rij(3)

           do n=1, 3
              atm_virit_const(1:3,n) = atm_virit_const(1:3,n) &
                   &                 - KRij(1:3)*Rij(n)
           end do
#endif

        end do

!           -- EXIT IF no atom velocities have been updated --

        if (finished) EXIT DOiter

!           -- Transfer mvnow to mvlast and go to next iteration -

!            do m = 1, natom
!               mvlast(m) = mvnow(m)
!               mvnow(m)  = .false.
!            end do
        mvlast(ib) = mvnow(ib)
        mvnow(ib)  = .false.
        mvlast(jb) = mvnow(jb)
        mvnow(jb)  = .false.

     END DO DOiter          ! end of loop over iterations

!        --- WRITE NUMBER OF ITERATIONS IF FAILED OR IF DEBUG MODE -

     if (.not. finished) then

!           - RATTLE failed to converge, so quitting -

        write(6,999)           maxiter
999     format(/3x, &
             & 'ERROR: RATTLE(V) failed to converge within MAXITER = ', &
             & i3)
        STOP
     end if

  END DO DOconstr

!     +     +     +     +     +     +     +     +     +     +

end subroutine rattlep_v
