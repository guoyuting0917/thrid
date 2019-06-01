!*****************************
!*  ewaldk_mpi.f Ver.1.6     *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-06-18 10:49:35 gota>

subroutine ewaldk(xcel,ycel,zcel,   &
     &            force,pot_ewk)

  use md_global
  use mpi_global

  implicit none

!
!     Calculate Ewald K-space force & potential
!                                   & pressure virial
!

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! total force
                                        !  (including vdw, bond, angle etc) 

!   OUTPUT

  real(8),intent(out):: pot_ewk         ! electrostatic potential

! LOCAL:

  real(8):: box(3)                      ! BOX size
  real(8):: poti                        ! temporal potential due to atom I  

  real(8):: vol                         ! volume = box(1)*box(2)*box(3)
  real(8):: inv_vol                     ! = 1/vol

  real(8):: pot_fact                    ! factor to multiply the potentail
                                        ! =  1/(2piV)

  real(8):: bsk, bck                    ! = sigmaj qj*sin(2pi * K.rj)
                                        ! = sigmaj qj*cos(2pi * K.rj)

  real(8):: csk(3),cck(3)               ! = bsk(bck) * A(K) * K

  real(8):: fi(3,maxnatom)              ! = sigmak {cck(3,k)*sin(2pi * K.rj)
                                        !         + csk(3,k)*cos(2pi * K.rj)}

  real(8):: pi                          ! 3.1415
  real(8):: pi2                         ! 2 * pi

  real(8):: theta
  real(8):: sini(natom),cosi(natom)

  integer::  i,k,m                      ! do loop index

  integer:: count_w(nproc),displ_w(nproc)
                          ! local nwave count and displacement of each process
  integer:: count_w3(nproc),displ_w3(nproc)
                          ! data transfer count and displacement of each process
  integer:: nwave_l,nwave_ll,mod_w
  integer:: ii,kk                       ! new loop indices
! MPI end

!     +     +     +     +     +     +     +

! MPI      count the number of nwaves and natom for each preocess.
  nwave_l = int(nwave/nproc)
  mod_w   = mod(nwave,nproc)

  do i=1,nproc
     nwave_ll = nwave_l
     if(i <= mod_w) nwave_ll = nwave_l + 1
     count_w(i) = nwave_ll
  end do

! MPI     data transfer number of nwave for each process
  do i=1,nproc
     count_w3(i) = count_w(i)*3
  end do

! MPI     diplacement indices for each process
!             w and w3 for nwave
  displ_w (1) = 0
  displ_w3(1) = 0
  do i=2,nproc
     displ_w (i) = displ_w (i-1)+count_w (i-1)
     displ_w3(i) = displ_w3(i-1)+count_w3(i-1)
  end do

! --- SOME PREPARATIONS ---

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel

  pi  = acos(-1.0d0)
  vol = box(1) * box(2) * box(3)

  pot_fact = 0.5d0 / pi / vol

  pi2 = 2.0d0 * pi

  inv_vol = 1.0d0 / vol

  pot_ewk = 0.0d0

  fi(1:3,1:natom) = 0.0d0

! MPI     only a part of nwave contributions are calculated in this process.
!        irank is a process number of this process.
! --- STEP 1 SCALAR SIN -  bsk = sigmaj qj*sin(2pi * K.rj)  and
!     STEP 2 SCALAR COS -  bck = sigmaj qj*cos(2pi * K.rj) ---
!  DO k = 1, nwave
  DO kk = 1,count_w(irank+1)

     k  = kk + displ_w(irank+1)

     bsk = 0.0d0
     bck = 0.0d0

     DO i = 1, natom                    ! summation over i

        theta  = kwave(1,k) * atmcor(1,i)   &
             & + kwave(2,k) * atmcor(2,i)   &
             & + kwave(3,k) * atmcor(3,i)
        theta  = pi2 * theta
        sini(i) = sin(theta)
        cosi(i) = cos(theta)
        bsk = bsk + atmchrg(i) * sini(i)
        bck = bck + atmchrg(i) * cosi(i)

     END DO

!    --- STEP 3 CALCULATE POTENTIAL = 1/2piV * sigmak A(K)*
!                                    {(bsk)^2 + (bck)^2}

     poti  = awave(k) * (bsk**2 + bck**2)
     pot_ewk = pot_ewk + poti

!    --- STEP 4 MAKE csk(:) = bsk * Ak * K
!                    cck(:) = bck * Ak * K

     csk(1:3) = awave(k) * bsk * kwave(1:3,k)
     cck(1:3) = awave(k) * bck * kwave(1:3,k)

! --- STEP 5 VECTOR SIN fsi(:) = sigmak csk(:) * sin(2pi * K.Ri(:)) and
!     STEP 6 VECTOR COS fci(:) = sigmak cck(:) * cos(2pi * K.Ri(:)) ---

     DO i = 1, natom

        fi(1:3,i) = fi(1:3,i) + cck(1:3) * sini(i) - csk(1:3) * cosi(i)

     END DO

  END DO

  pot_ewk = pot_fact * pot_ewk          ! pot_fact = 1/(2piV)

! --- STEP 7 CALUCLATE FINAL FORCE TO ADD ---

  DO i = 1, natom

     fi(1:3,i) = 2.0d0 * atmchrg(i) * inv_vol * fi(1:3,i)

  END DO

  do i = 1, natom
     force(1:3,i) = force(1:3,i) + fi(1:3,i)
  end do

!     +     +     +     +     +     +     +

end subroutine ewaldk

! -------------------------------------------------------------

subroutine ewaldkp(xcel,ycel,zcel,   &
     &             alpha,   &
     &             force,pot_ewk,   &
     &             ifnetqcorrp,   &
     &             netchrgsq,   &
     &             for_viri_coul,pot_viri_coul,pot_virit_coul)

  use md_global
  use mpi_global

  implicit none

!
!     Calculate Ewald K-space force & potential
!                                   & pressure virial
!

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  real(8),intent(in):: alpha            ! parameter alpha [non-d]

  logical,intent(in):: ifnetqcorrp      ! net charge correction for pressure
  real(8),intent(in):: netchrgsq        ! = (sum(qi))**2

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! total force
                                        !  (including vdw, bond, angle etc) 

  real(8),intent(inout):: for_viri_coul(:,:) ! virial(coulomb force) of each atom
  real(8),intent(inout):: pot_viri_coul ! virial(coulomb potential) of each atom
  real(8),intent(inout):: pot_virit_coul(:,:) ! virial tensor (coulomb)

!   OUTPUT
  real(8),intent(out):: pot_ewk         ! electrostatic potential

! LOCAL:

  real(8):: box(3)                      ! BOX size
  real(8):: poti                        ! temporal potential due to atom I  

  real(8):: vol                         ! volume = box(1)*box(2)*box(3)
  real(8):: inv_vol                     ! = 1/vol
  real(8):: pot_fact                    ! factor to multiply the potentail
                                        ! =  1/(2piV)

  real(8):: bsk, bck                    ! = sigmaj qj*sin(2pi * K.rj)
                                        ! = sigmaj qj*cos(2pi * K.rj)

  real(8):: csk(3),cck(3)               ! = bsk(bck) * A(K) * K

  real(8):: fi(3,maxnatom)              ! = sigmak {cck(3,k)*sin(2pi * K.rj)
                                        !         + csk(3,k)*cos(2pi * K.rj)}

  real(8):: pi                          ! 3.1415
  real(8):: pi2                         ! 2 * pi 

  real(8):: theta
  real(8):: sini(natom),cosi(natom)

  integer:: i,k,m,n                     ! do loop index

  real(8):: ktmp2
  real(8):: k_fact(3,3)

  real(8):: inv2_alpha                  ! = 1/(alpha*alpha)

  real(8):: fac_vir                     ! = 2*pi^2/alpha^2
  real(8):: fac_vir2                    ! = pot_fact*awave(i)*(bsk(i)^2+bck(i)^2)
  real(8):: fac_vir3                    ! = pot_fact*awave(i)*netchrgsq
  real(8):: fac_vir4                    ! = pot_fact*awave(i)*netchrgsq
                                        !  *(-1 + fac_vir*ktmp2)

  real(8):: pot_ewk_cc                  ! ewk potentioal (charge corr)

  integer:: count_w(nproc),displ_w(nproc)
                          ! local nwave count and displacement of each process
  integer:: count_w3(nproc),displ_w3(nproc)
                          ! data transfer count and displacement of each process
  integer:: nwave_l,nwave_ll,mod_w  
  integer:: ii,kk                       ! new loop indices

!     +     +     +     +     +     +     +

! MPI      count the number of nwaves and natom for each preocess.
!             w for nwave
  nwave_l = int(nwave/nproc)
  mod_w   = mod(nwave,nproc)
  do i=1,nproc
     nwave_ll = nwave_l
     if(i <= mod_w) nwave_ll = nwave_l + 1
     count_w(i) = nwave_ll
  end do

! MPI     data transfer number of nwave for each process
  do i=1,nproc
     count_w3(i) = count_w(i)*3
  end do

! MPI     diplacement indices for each process
!             w and w3 for nwave
  displ_w (1) = 0
  displ_w3(1) = 0
  do i=2,nproc
     displ_w (i) = displ_w (i-1)+count_w (i-1)
     displ_w3(i) = displ_w3(i-1)+count_w3(i-1)
  end do

! --- SOME PREPARATIONS ---

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel

  pi  = acos(-1.0d0)
  vol = box(1) * box(2) * box(3)

  pot_fact = 0.5d0 / pi / vol

  pi2 = 2.0d0 * pi

  inv_vol = 1.0d0 / vol

  inv2_alpha = 1.0d0 / (alpha*alpha)

  fac_vir = 2.0d0*pi*pi*inv2_alpha

  pot_ewk    = 0.0d0
  pot_ewk_cc = 0.0d0

  fi(1:3,1:natom) = 0.0d0

! MPI     only a part of nwave contributions are calculated in this process.
!        irank is a process number of this process.
!  --- STEP 1 SCALAR SIN -  bsk = sigmaj qj*sin(2pi * K.rj)  and
!      STEP 2 SCALAR COS -  bck = sigmaj qj*cos(2pi * K.rj) ---
!  DO k = 1, nwave
  DO kk = 1,count_w(irank+1)

     k  = kk + displ_w(irank+1)

     bsk = 0.0d0
     bck = 0.0d0

     DO i = 1, natom                    ! summation over i

        theta  = kwave(1,k) * atmcor(1,i)   &
             & + kwave(2,k) * atmcor(2,i)   &
             & + kwave(3,k) * atmcor(3,i)
        theta  = pi2 * theta
        sini(i) = sin(theta)
        cosi(i) = cos(theta)
        bsk = bsk + atmchrg(i) * sini(i)
        bck = bck + atmchrg(i) * cosi(i)

     END DO

!    --- STEP 3 CALCULATE POTENTIAL = 1/2piV * sigmak A(K)*
!                                      {(bsk)^2 + (bck)^2}

     poti  = awave(k) * (bsk**2 + bck**2)
     pot_ewk = pot_ewk + poti
     fac_vir3 = pot_fact * awave(k) * netchrgsq

!    --- STEP 4 MAKE csk(:,k) = bsk * Ak * K
!                    cck(:,k) = bck * Ak * K

     csk(1:3) = awave(k) * bsk * kwave(1:3,k)
     cck(1:3) = awave(k) * bck * kwave(1:3,k)

! --- STEP 5 VECTOR SIN fsi(:) = sigmak csk(:) * sin(2pi * K.Ri(:)) and
!     STEP 6 VECTOR COS fci(:) = sigmak cck(:) * cos(2pi * K.Ri(:)) ---

     DO i = 1, natom

        fi(1:3,i) = fi(1:3,i) + cck(1:3) * sini(i) - csk(1:3) * cosi(i)

     END DO

!    --- Calculate virial tensor  ---
     if (k == 1) cycle

     ktmp2 = kwave(1,k)*kwave(1,k) + kwave(2,k)*kwave(2,k)   &
         & + kwave(3,k)*kwave(3,k)
     fac_vir2 = pot_fact * poti

     if (.not. ifnetqcorrp) then        ! no charge correction

        do n=1,3

           k_fact(1:3,n) = -1.0d0*(2.0d0/ktmp2 + fac_vir)   &
                &        * kwave(1:3,k)*kwave(n,k)

           k_fact(n,n) = k_fact(n,n) + 1.0d0

           k_fact(1:3,n) =  k_fact(1:3,n) * fac_vir2

           pot_virit_coul(1:3,n) = pot_virit_coul(1:3,n)   &
                &                + k_fact(1:3,n)

        end do

     else                               ! charge correction

        fac_vir4 = fac_vir3 * (-1.0d0 + fac_vir*ktmp2)
        pot_ewk_cc = pot_ewk_cc + fac_vir4

        do n=1,3

           k_fact(1:3,n) = -1.0d0*(2.0d0/ktmp2 + fac_vir)   &
                &        * kwave(1:3,k)*kwave(n,k)

           k_fact(n,n) = k_fact(n,n) + 1.0d0

           k_fact(1:3,n) = k_fact(1:3,n) * (fac_vir2 - fac_vir3)

           pot_virit_coul(1:3,n) = pot_virit_coul(1:3,n)   &
                &                + k_fact(1:3,n)

        end do

     end if

  END DO

  pot_ewk = pot_fact * pot_ewk          ! pot_fact = 1/(2piV)

  pot_viri_coul =  pot_viri_coul + pot_ewk + pot_ewk_cc


! --- STEP 7 CALUCLATE FINAL FORCE TO ADD ---
!         fi =  2Qi/V * (fsi - fci)
!         fi will be stored in fsi

  DO i = 1, natom

     fi(1:3,i) = 2.0d0 * atmchrg(i) * inv_vol * fi(1:3,i)

  END DO

  do i = 1, natom
     force(1:3,i) = force(1:3,i) + fi(1:3,i)
     for_viri_coul(1:3,i) = for_viri_coul(1:3,i) + fi(1:3,i)
  end do

!     +     +     +     +     +     +     +

end subroutine ewaldkp
