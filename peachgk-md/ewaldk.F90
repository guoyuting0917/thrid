!********************************
!*  ewaldk.f Ver.2.0 '10.06.28  *
!*      for peachgk_md.f        *
!*            by G.Kikugawa     *
!********************************      
subroutine ewaldk(xcel,ycel,zcel,   &
     &            force,pot_ewk)

  use md_global

  implicit none

!    
!     Calculate Ewald K-space force & potential   
!                                   & pressure virial
!
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

  real(8):: bsk(maxwave),bck(maxwave)   ! = sigmaj qj*sin(2pi * K.rj)
                                        ! = sigmaj qj*cos(2pi * K.rj) 
  real(8):: csk(3,maxwave),cck(3,maxwave) ! = bsk(bck) * A(K) * K
  real(8):: fsi(3,maxnatom)             ! = sigmak cck(3,k)*sin(2pi * K.rj)
  real(8):: fci(3,maxnatom)             ! = sigmak csk(3,k)*cos(2pi * K.rj)

  real(8):: pi                          ! 3.1415
  real(8):: pi2                         ! 2 * pi 

  real(8)::  theta

  integer::  i,j,k                      ! do loop index

!     +     +     +     +     +     +     +


! --- SOME PREPARATIONS ---

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel

  pi  = acos(-1.0d0)
  vol = box(1) * box(2) * box(3)

  pot_fact = 0.5d0 / pi / vol

  pi2 = 2.0d0 * pi

  inv_vol = 1.0d0 / vol

! --- INITIALIZE ARRAYS FOR FORCE CALCULATED IN THE SUBROUTINE ---

  bsk(1:maxwave)   = 0.0d0
  bck(1:maxwave)   = 0.0d0

  csk(1:3,1:maxwave) = 0.0d0
  cck(1:3,1:maxwave) = 0.0d0

  fsi(1:3,1:maxnatom) = 0.0d0
  fci(1:3,1:maxnatom) = 0.0d0

! --- STEP 1 SCALAR SIN -  bsk = sigmaj qj*sin(2pi * K.rj)  and
!     STEP 2 SCALAR COS -  bck = sigmaj qj*cos(2pi * K.rj) ---


  DO k = 1, nwave

     DO j = 1, natom                    ! summation over J

        theta  = kwave(1,k) * atmcor(1,j)   &
             & + kwave(2,k) * atmcor(2,j)   &
             & + kwave(3,k) * atmcor(3,j)
        theta  = pi2 * theta
        bsk(k) = bsk(k) + atmchrg(j) * dsin(theta)
        bck(k) = bck(k) + atmchrg(j) * dcos(theta)

     END DO

  END DO


! --- STEP 3 CALCULATE POTENTIAL = 1/2piV * sigmak A(K)*
!                                    {(bsk)^2 + (bck)^2}

  pot_ewk = 0.0d0
  do i = 1, nwave
     poti  = awave(i) * (bsk(i)**2 + bck(i)**2)
     pot_ewk = pot_ewk + poti
  end do
  pot_ewk = pot_fact * pot_ewk ! pot_fact = 1/(2piV)


! --- STEP 4 MAKE csk(:,k) = bsk * Ak * K
!                 cck(:,k) = bck * Ak * K

  do k = 1, nwave
     csk(1:3,k) = awave(k) * bsk(k) * kwave(1:3,k)
     cck(1:3,k) = awave(k) * bck(k) * kwave(1:3,k)
  end do


! --- STEP 5 VECTOR SIN fsi(:) = sigmak csk(:) * sin(2pi * K.Ri(:)) and
!     STEP 6 VECTOR COS fci(:) = sigmak cck(:) * cos(2pi * K.Ri(:)) ---


  DO i = 1, natom

     DO k = 1, nwave        ! summation over K

        theta  = kwave(1,k) * atmcor(1,i)   &
             & + kwave(2,k) * atmcor(2,i)   &
             & + kwave(3,k) * atmcor(3,i)
        theta  =   pi2 * theta

        fsi(1:3,i) = fsi(1:3,i) + cck(1:3,k) * dsin(theta)
        fci(1:3,i) = fci(1:3,i) + csk(1:3,k) * dcos(theta)

     END DO

  END DO


! --- STEP 7 CALUCLATE FINAL FORCE TO ADD ---
!         fi =  2Qi/V * (fsi - fci)
!         fi will be stored in fsi

  DO i = 1, natom
!     fsi(1,i) = fsi(1,i) - fci(1,i)
!     fsi(1,i) = 2.0d0 * atmchrg(i) * inv_vol * fsi(1,i)
!     fsi(2,i) = fsi(2,i) - fci(2,i)
!     fsi(2,i) = 2.0d0 * atmchrg(i) * inv_vol * fsi(2,i)
!     fsi(3,i) = fsi(3,i) - fci(3,i)
!     fsi(3,i) = 2.0d0 * atmchrg(i) * inv_vol * fsi(3,i)     
     fsi(1:3,i) = 2.0d0 * atmchrg(i) * inv_vol * (fsi(1:3,i) - fci(1:3,i))

  END DO


! --- ADD FSI TO FORCE  ---

  force(1:3,1:natom) = force(1:3,1:natom) + fsi(1:3,1:natom)

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

  implicit none

!    
!     Calculate Ewald K-space force & potential   
!                                   & pressure virial
!
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

  real(8):: bsk(maxwave),bck(maxwave)   ! = sigmaj qj*sin(2pi * K.rj)
                                        ! = sigmaj qj*cos(2pi * K.rj) 
  real(8):: csk(3,maxwave),cck(3,maxwave) ! = bsk(bck) * A(K) * K
  real(8):: fsi(3,maxnatom)             ! = sigmak cck(3,k)*sin(2pi * K.rj)
  real(8):: fci(3,maxnatom)             ! = sigmak csk(3,k)*cos(2pi * K.rj)
                                
  real(8):: pi                          ! 3.1415
  real(8):: pi2                         ! 2 * pi 

  real(8)::  theta

  integer::  i,j,k,n                    ! do loop index

  real(8):: ktmp2                       ! |kwave|**2
  real(8):: k_fact(3,3)

  real(8):: inv2_alpha                  ! = 1/(alpha*alpha)

  real(8):: fac_vir                    ! = 2*pi^2/alpha^2
  real(8):: fac_vir2                   ! = pot_fact*awave(i)*(bsk(i)^2+bck(i)^2)
  real(8):: fac_vir3                   ! = pot_fact*awave(i)*netchrgsq
  real(8):: fac_vir4                   ! = pot_fact*awave(i)*netchrgsq
                                       !  *(-1 + fac_vir*ktmp2)
      
!     +     +     +     +     +     +     +


! --- SOME PREPARATIONS ---

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel

  pi  = dacos(-1.0d0)
  vol = box(1) * box(2) * box(3)

  pot_fact = 0.5d0 / pi / vol

  pi2 = 2.0d0 * pi

  inv_vol = 1.0d0 / vol

  inv2_alpha = 1.0d0 / (alpha*alpha)

! --- INITIALIZE ARRAYS FOR FORCE CALCULATED IN THE SUBROUTINE ---

  bsk(1:maxwave)   = 0.0d0
  bck(1:maxwave)   = 0.0d0

  csk(1:3,1:maxwave) = 0.0d0
  cck(1:3,1:maxwave) = 0.0d0

  fsi(1:3,1:maxnatom) = 0.0d0
  fci(1:3,1:maxnatom) = 0.0d0

! --- STEP 1 SCALAR SIN -  bsk = sigmaj qj*sin(2pi * K.rj)  and
!     STEP 2 SCALAR COS -  bck = sigmaj qj*cos(2pi * K.rj) ---


  DO k = 1, nwave

     DO j = 1, natom                    ! summation over J

        theta  = kwave(1,k) * atmcor(1,j)   &
             & + kwave(2,k) * atmcor(2,j)   &
             & + kwave(3,k) * atmcor(3,j)
        theta  = pi2 * theta
        bsk(k) = bsk(k) + atmchrg(j) * dsin(theta)
        bck(k) = bck(k) + atmchrg(j) * dcos(theta)

     END DO

  END DO


! --- STEP 3 CALCULATE POTENTIAL = 1/2piV * sigmak A(K)*
!                                     {(bsk)^2 + (bck)^2}

  pot_ewk = 0.0d0
  do i = 1, nwave
     poti  = awave(i) * (bsk(i)**2 + bck(i)**2)
     pot_ewk = pot_ewk + poti
  end do
  pot_ewk = pot_fact * pot_ewk          ! pot_fact = 1/(2piV)

  pot_viri_coul =  pot_viri_coul + pot_ewk

! --- STEP 4 MAKE csk(:,k) = bsk * Ak * K
!                 cck(:,k) = bck * Ak * K

  do k = 1, nwave
     csk(1:3,k) = awave(k) * bsk(k) * kwave(1:3,k)
     cck(1:3,k) = awave(k) * bck(k) * kwave(1:3,k)
  end do


! --- STEP 5 VECTOR SIN fsi(:) = sigmak csk(:) * sin(2pi * K.Ri(:)) and
!     STEP 6 VECTOR COS fci(:) = sigmak cck(:) * cos(2pi * K.Ri(:)) ---


  DO i = 1, natom

     DO k = 1, nwave                    ! summation over K

        theta  = kwave(1,k) * atmcor(1,i)   &
             & + kwave(2,k) * atmcor(2,i)   &
             & + kwave(3,k) * atmcor(3,i)
        theta  =   pi2 * theta

        fsi(1:3,i) = fsi(1:3,i) + cck(1:3,k) * dsin(theta)
        fci(1:3,i) = fci(1:3,i) + csk(1:3,k) * dcos(theta)

     END DO

  END DO


! --- STEP 7 CALUCLATE FINAL FORCE TO ADD ---
!         fi =  2Qi/V * (fsi - fci)
!         fi will be stored in fsi

  DO i = 1, natom
!     fsi(1,i) = fsi(1,i) - fci(1,i)
!     fsi(1,i) = 2.0d0 * atmchrg(i) * inv_vol * fsi(1,i)
!     fsi(2,i) = fsi(2,i) - fci(2,i)
!     fsi(2,i) = 2.0d0 * atmchrg(i) * inv_vol * fsi(2,i)
!     fsi(3,i) = fsi(3,i) - fci(3,i)
!     fsi(3,i) = 2.0d0 * atmchrg(i) * inv_vol * fsi(3,i)     
     fsi(1:3,i) = 2.0d0 * atmchrg(i) * inv_vol * (fsi(1:3,i) - fci(1:3,i))

  END DO


! --- ADD FSI TO FORCE  ---

  force(1:3,1:natom) = force(1:3,1:natom) + fsi(1:3,1:natom)
         
  for_viri_coul(1:3,1:natom) = for_viri_coul(1:3,1:natom) + fsi(1:3,1:natom)

! --- Calculate virial tensor  ---

  fac_vir = 2.0d0*pi*pi*inv2_alpha

  if (.not. ifnetqcorrp) then           ! no charge correction

     do i = 2, nwave                    ! i=1(n=(0,0,0)) should be excluded
        ktmp2 = kwave(1,i)*kwave(1,i) + kwave(2,i)*kwave(2,i)   &
            & + kwave(3,i)*kwave(3,i)

        fac_vir2 = pot_fact * awave(i) * (bsk(i)**2 + bck(i)**2)
         
        do n=1,3

           k_fact(1:3,n) = -1.0d0*(2.0d0/ktmp2 + fac_vir)   &
                &        * kwave(1:3,i)*kwave(n,i)

           k_fact(n,n) = k_fact(n,n) + 1.0d0

           k_fact(1:3,n) = k_fact(1:3,n) * fac_vir2

           pot_virit_coul(1:3,n) = pot_virit_coul(1:3,n) + k_fact(1:3,n)

        end do

     end do

  else                                  ! charge correction

     do i = 2, nwave                    ! i=1(n=(0,0,0)) should be excluded
        ktmp2 = kwave(1,i)*kwave(1,i) + kwave(2,i)*kwave(2,i)   &
     &        + kwave(3,i)*kwave(3,i)

        fac_vir3 = pot_fact * awave(i)
        fac_vir2 = fac_vir3 * (bsk(i)**2 + bck(i)**2)
        fac_vir3 = fac_vir3 * netchrgsq

        fac_vir4 = fac_vir3 * (-1.0d0 + fac_vir*ktmp2)
        pot_viri_coul = pot_viri_coul + fac_vir4

        do n=1,3

           k_fact(1:3,n) = -1.0d0*(2.0d0/ktmp2 + fac_vir)   &
                &        * kwave(1:3,i)*kwave(n,i)

           k_fact(n,n) = k_fact(n,n) + 1.0d0

           k_fact(1:3,n) =  k_fact(1:3,n) * (fac_vir2 - fac_vir3)

           pot_virit_coul(1:3,n) = pot_virit_coul(1:3,n) + k_fact(1:3,n)

        end do

     end do
     
  end if

!     +     +     +     +     +     +     +

end subroutine ewaldkp
