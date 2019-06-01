!***************************************
!*  rattle_mtk_mpi.f Ver.1.4 '10.07.09 *
!*      for peachgk_md.f               *
!*            by G.Kikugawa            *
!***************************************
subroutine rattle_cmtk( atmcor_old, &
     &                  eps_rattle,dt_short_cal,dt_long_cal, &
     &                  istep_short,nstep_short, &
     &                  rv )

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
  use mpi_global

  implicit none

! ARGUMENT:
!     INPUT
  real(8):: atmcor_old(:,:)  ! atmcor at T

  real(8):: eps_rattle       ! tolerance (relative difference)
                                ! for bond length constraint by RATTLE

  real(8):: dt_short_cal     ! time step of short force
  real(8):: dt_long_cal      ! time step of long force [non-d]

  integer:: istep_short
  integer:: nstep_short     ! number of step for short force

  real(8):: rv               ! = poly*exp(dt_short/2*veps)

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
             
  logical:: finished = .false. ! set true if no more coordinate resetting
                                ! is needed.

  real(8):: one = 1.0d0  
  real(8):: two = 2.0d0

  real(8):: eps_rattle2      ! = eps_rattle^2

  integer:: i, m      ! do loop indexes

  real(8):: dt_cal           ! time step

! MPI  The following arrays and variables are declared in MPI version.
  real(8):: atmcor_tmp(3,natom) !temporary atmcor
  real(8):: atmvel_tmp(3,natom) !temporary atmvel
  real(8):: atmcor_gbl(3,natom) !temporary global atmcor
  real(8):: atmvel_gbl(3,natom) !temporary global atmvel

  integer:: count_mc(nproc)
  integer:: istamc(nproc),iendmc(nproc)
  integer:: mod_mc
  integer:: mcconst_ll,mcconst_l

  integer:: i1,i2,ii
  integer:: nc
  integer:: ir
! MPI end

!     +     +     +     +     +     +     +     +     +     +

!     --- initialization ---
  eps_rattle2 = eps_rattle * eps_rattle

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

!     --- some preparation ---

! MPI      count the number of ncount for each preocess.
!             c for nwave
  mcconst_l = int(nmconstindex/nproc)
  mod_mc = mod(nmconstindex,nproc)
  do i=1,nproc
     mcconst_ll = mcconst_l
     if(i <= mod_mc) mcconst_ll = mcconst_l + 1
     count_mc(i) = mcconst_ll
  end do

! MPI     start and end indices for each process
!             w and w3 for nwave
  istamc(1) = 1
  iendmc(1) = count_mc(1)
  do i=2,nproc
     istamc(i) = iendmc(i-1) + 1
     iendmc(i) = iendmc(i-1) + count_mc(i)
  end do

! MPI     set the temporary atmcor and atmvel
  atmcor_tmp(1:3,1:natom) = 0.0d0
  atmvel_tmp(1:3,1:natom) = 0.0d0
  atmcor_gbl(1:3,1:natom) = 0.0d0
  atmvel_gbl(1:3,1:natom) = 0.0d0

  do i = istamc(irank+1), iendmc(irank+1)
     i1 = index_mconst(i)
     i2 = index_mconst(i+1) - 1
     do ii = i1, i2
        nc = list_mconst(ii)
        atmcor_tmp(1:3,iconst(nc)) = atmcor(1:3,iconst(nc))
        atmcor_tmp(1:3,jconst(nc)) = atmcor(1:3,jconst(nc))
        atmvel_tmp(1:3,iconst(nc)) = atmvel(1:3,iconst(nc))
        atmvel_tmp(1:3,jconst(nc)) = atmvel(1:3,jconst(nc))
     end do
  end do

!     --- GRAND LOOP TO APPLY BOND LENGHT CONSTRAINT ---

!      * this loop is iterated until all the constraints
!        are satisfied.

!        -- LOOP OVER CONSTRAINTS  --
! MPI 
!         DOconstr: DO i = 1, nconst
  DOconstr: DO i = istamc(irank+1), iendmc(irank+1)

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

           ib  = iconst(nc)  ! bond length between ib & jb  
           jb  = jconst(nc)  ! will be constraint
           dij = dconst(nc)  ! to dij

!           - CYCLE DOcostr -
!              if neither ib nor jb have been updated last iteration 
!                                          or in the current iteration.      
!              * if both ib & jb were not moved in the last or current 
!                iterations , then they need not be moved this iteration.
                
           IF (.not. (mvlast(ib) .or. mvlast(jb) .or. &
                & mvnow(ib) .or. mvnow(jb))) then
              cycle
           END IF

!             - SEE if Rij satisfy the constraint -

           Rij(1:3) = atmcor_tmp(1:3,ib) - atmcor_tmp(1:3,jb)

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

           gij  = (Rij2-dij2)/(two * RR * (amassi_inv + amassj_inv)) &
                & / dt_cal      ! additional dividing timestep

           GRij(1:3) = gij*Rij_old(1:3)

           atmcor_tmp(1:3,ib) =  atmcor_tmp(1:3,ib) &
                &              - amassi_inv * GRij(1:3) * dt_short_cal
           atmcor_tmp(1:3,jb) =  atmcor_tmp(1:3,jb) &
                &              + amassj_inv * GRij(1:3) * dt_short_cal

!              - Velocity update -
!              divided by rv
           GRij(1:3) = GRij(1:3) / rv

           atmvel_tmp(1:3,ib) = atmvel_tmp(1:3,ib) &
                &              - amassi_inv * GRij(1:3) 
           atmvel_tmp(1:3,jb) = atmvel_tmp(1:3,jb) &
                &              + amassj_inv * GRij(1:3)

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

! MPI    Then sum-up all atmcor and atmvel
!
  call mpi_allreduce( atmcor_tmp, atmcor_gbl, 3*natom, &
       &              MPI_DOUBLE_PRECISION, &
       &              MPI_SUM, MPI_COMM_WORLD, ierror)

  call mpi_allreduce( atmvel_tmp,atmvel_gbl, 3*natom, &
       &              MPI_DOUBLE_PRECISION, &
       &              MPI_SUM, MPI_COMM_WORLD, ierror)

  do i = 1, nconst

     ib  = iconst(i)
     jb  = jconst(i)

     atmcor(1:3,ib) = atmcor_gbl(1:3,ib)
     atmcor(1:3,jb) = atmcor_gbl(1:3,jb)
     atmvel(1:3,ib) = atmvel_gbl(1:3,ib)
     atmvel(1:3,jb) = atmvel_gbl(1:3,jb)

  end do

!     +     +     +     +     +     +     +     +     +     +

  return
end subroutine rattle_cmtk

!-----------------------------------------------------------------------
subroutine rattlep_v( eps_rattle, &
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
  use mpi_global

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

! MPI  The following arrays and variables are declared in MPI version.
  real(8):: atmvel_tmp(3,natom) !temporaty atmvel
  real(8):: atmvel_gbl(3,natom) !temporary global atmvel
  real(8):: av_const_lcl     !local atm_viri_const
  real(8):: avt_const_lcl(3,3) !local atm_virit_const

  integer:: count_mc(nproc)
  integer:: istamc(nproc),iendmc(nproc)
  integer:: mod_mc
  integer:: mcconst_ll,mcconst_l

  integer:: i1,i2,ii
  integer:: nc
  integer:: ir

  real(8):: local(7),global(7)
! MPI end

!     +     +     +     +     +     +     +     +     +     +

!     --- initialization ---
  eps_rattle2 = eps_rattle * eps_rattle

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

!     --- some preparation ---

! MPI      count the number of ncount for each preocess.
!             c for nwave
  mcconst_l = int(nmconstindex/nproc)
  mod_mc = mod(nmconstindex,nproc)
  do i=1,nproc
     mcconst_ll = mcconst_l
     if(i <= mod_mc) mcconst_ll = mcconst_l + 1
     count_mc(i) = mcconst_ll
  end do

! MPI     start and end indices for each process
!             w and w3 for nwave
  istamc(1) = 1
  iendmc(1) = count_mc(1)
  do i=2,nproc
     istamc(i) = iendmc(i-1) + 1
     iendmc(i) = iendmc(i-1) + count_mc(i)
  end do

! MPI     set the temporary atmcor and atmvel
  atmvel_tmp(1:3,1:natom) = 0.0d0

  do i = istamc(irank+1), iendmc(irank+1)
     i1 = index_mconst(i)
     i2 = index_mconst(i+1) - 1
     do ii = i1, i2
        nc = list_mconst(ii)
        atmvel_tmp(1:3,iconst(nc)) = atmvel(1:3,iconst(nc))
        atmvel_tmp(1:3,jconst(nc)) = atmvel(1:3,jconst(nc))
     end do
  end do

!     --- GRAND LOOP TO APPLY VELOCITY CONSTRAINT ---

!     --- Initilize virial for constraint force ---

  av_const_lcl = 0.0d0
  avt_const_lcl(1,1) = 0.0d0
  avt_const_lcl(1,2) = 0.0d0
  avt_const_lcl(1,3) = 0.0d0
  avt_const_lcl(2,2) = 0.0d0
  avt_const_lcl(2,3) = 0.0d0
  avt_const_lcl(3,3) = 0.0d0

!      * this loop is iterated until all the constraints
!        are satisfied.

!        -- LOOP OVER CONSTRAINTS  --
! MPI 
!         DOconstr: DO i = 1, nconst
  DOconstr: DO i = istamc(irank+1), iendmc(irank+1)

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

!              - SEE if (Rij.Vij) satisfy the constraint -

           Rij(1:3) = atmcor(1:3,ib) - atmcor(1:3,jb)
           Vij(1:3) = atmvel_tmp(1:3,ib) - atmvel_tmp(1:3,jb)

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

           atmvel_tmp(1:3,ib) = atmvel_tmp(1:3,ib) &
                &             - amassi_inv * KRij(1:3)
           atmvel_tmp(1:3,jb) = atmvel_tmp(1:3,jb) &
                &             + amassj_inv * KRij(1:3)

!              - Virial update -

           KRij(1:3) = 2.0d0*KRij(1:3)*inv_vtime

           av_const_lcl = av_const_lcl - KRij(1)*Rij(1) &
                &                      - KRij(2)*Rij(2) &
                &                      - KRij(3)*Rij(3)

           avt_const_lcl(1,1) =  avt_const_lcl(1,1) &
                &              - KRij(1)*Rij(1)
           avt_const_lcl(1,2) =  avt_const_lcl(1,2) &
                &              - KRij(1)*Rij(2)
           avt_const_lcl(1,3) =  avt_const_lcl(1,3) &
                &              - KRij(1)*Rij(3)
           avt_const_lcl(2,2) =  avt_const_lcl(2,2) &
                &              - KRij(2)*Rij(2)
           avt_const_lcl(2,3) =  avt_const_lcl(2,3) &
                &              - KRij(2)*Rij(3)
           avt_const_lcl(3,3) =  avt_const_lcl(3,3) &
                &              - KRij(3)*Rij(3)

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

! MPI     sum-up all atmvel
!
  call mpi_allreduce( atmvel_tmp,atmvel_gbl,3*natom, &
       &              MPI_DOUBLE_PRECISION, &
       &              MPI_SUM,MPI_COMM_WORLD,ierror)

  do i = 1, nconst

     ib  = iconst(i)
     jb  = jconst(i)

     atmvel(1:3,ib) = atmvel_gbl(1:3,ib)
     atmvel(1:3,jb) = atmvel_gbl(1:3,jb)

  end do

! MPI    sum-up all virial
  local(1) = av_const_lcl
      
  local(2) = avt_const_lcl(1,1)
  local(3) = avt_const_lcl(1,2)
  local(4) = avt_const_lcl(1,3)
  local(5) = avt_const_lcl(2,2)
  local(6) = avt_const_lcl(2,3)
  local(7) = avt_const_lcl(3,3)

  call mpi_allreduce( local,global,7, &
       &              MPI_DOUBLE_PRECISION, &
       &              MPI_SUM,MPI_COMM_WORLD,ierror)
      
  atm_viri_const = global(1)

  atm_virit_const(1,1) = global(2)
  atm_virit_const(1,2) = global(3)
  atm_virit_const(1,3) = global(4)
  atm_virit_const(2,2) = global(5)
  atm_virit_const(2,3) = global(6)
  atm_virit_const(3,3) = global(7)

  atm_virit_const(2,1) = atm_virit_const(1,2)
  atm_virit_const(3,1) = atm_virit_const(1,3)
  atm_virit_const(3,2) = atm_virit_const(2,3)
      
!     +     +     +     +     +     +     +     +     +     +

  return
end subroutine rattlep_v
