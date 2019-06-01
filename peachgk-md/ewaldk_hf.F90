!*****************************
!*  ewaldk_hf.f90 Ver.1.0    *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-06-18 01:08:43 gota>

subroutine ewaldkp_hf(xcel,ycel,zcel, &
     &                alpha, &
     &                force,pot_ewk, &
     &                ifnetqcorrp, &
     &                netchrgsq, &
     &                for_viri_coul,pot_viri_coul,pot_virit_coul, &
     &                pot_elc_atm,virielct_atm, &
     &                ifhfvol, &
     &                nhfregion,hfzpos1,hfzpos2, &
     &                hftyp_atm, &
     &                molecom)

  use md_global

  implicit none

!
!     Calculate Ewald K-space force & potential (per-atom)
!                                   & pressure virial (per-atom)
!      for heat flux in the bulk region
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

  logical,intent(in):: ifhfvol       ! local volume-based or local surface-based

  integer,intent(in):: nhfregion     ! number of region to calculate heat flux
  real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                     ! z-position of region for heat flux

  integer,intent(in):: hftyp_atm(:)  ! atom- or mole-based heat flux cal.
                                     !   for each atom

  real(8),intent(in):: molecom(:,:)  ! center of mass of molecule

!   OUTPUT
  real(8),intent(out):: pot_ewk         ! electrostatic potential

  real(8),intent(inout):: pot_elc_atm(:)   ! Coulomb potential of each atom
  real(8),intent(inout):: virielct_atm(:,:,:,:)
                                     ! virial tensor of each atom (Coulomb)

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

  integer::  i,k,n                      ! do loop index

  real(8):: ktmp2                       ! |kwave|**2
  real(8):: k_fact(3,3)

  real(8):: inv2_alpha                  ! = 1/(alpha*alpha)

  real(8):: fac_vir                    ! = 2*pi^2/alpha^2
  real(8):: fac_vir2                   ! = pot_fact*awave(i)*(bsk(i)^2+bck(i)^2)
  real(8):: fac_vir3                   ! = pot_fact*awave(i)*netchrgsq
  real(8):: fac_vir4                   ! = pot_fact*awave(i)*netchrgsq
                                       !  *(-1 + fac_vir*ktmp2)
  real(8):: pot_fact_awave             ! = pot_fact*awave(i)

  real(8):: pot_ewk_cc                  ! ewk potentioal (charge corr)

  real(8):: alpha_sqrtpi               ! = alpha / sqrt(pi)

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

  inv2_alpha = 1.0d0 / (alpha*alpha)

  fac_vir = 2.0d0*pi*pi*inv2_alpha

  alpha_sqrtpi = alpha / sqrt(pi)

  pot_ewk    = 0.0d0
  pot_ewk_cc = 0.0d0

  fi(1:3,1:natom) = 0.0d0

! --- STEP 1 SCALAR SIN -  bsk = sigmaj qj*sin(2pi * K.rj)  and
!     STEP 2 SCALAR COS -  bck = sigmaj qj*cos(2pi * K.rj) ---

  DO k = 1, nwave

     bsk = 0.0d0
     bck = 0.0d0

     DO i = 1, natom                    ! summation over i

        theta  = kwave(1,k) * atmcor(1,i) &
             & + kwave(2,k) * atmcor(2,i) &
             & + kwave(3,k) * atmcor(3,i)
        theta  = pi2 * theta
        sini(i) = sin(theta)
        cosi(i) = cos(theta)
        bsk = bsk + atmchrg(i) * sini(i)
        bck = bck + atmchrg(i) * cosi(i)

     END DO

! --- STEP 3 CALCULATE POTENTIAL = 1/2piV * sigmak A(K)*
!                                     {(bsk)^2 + (bck)^2}

     pot_fact_awave = pot_fact * awave(k)
     poti  = pot_fact_awave * (bsk**2 + bck**2)
     pot_ewk = pot_ewk + poti

     !- calculate potential of each atom (_HF_BULK_EWK)
     !        pot_ewk(i) =  1/2piV * sigmak A(K)* {(bski*(bsk) + bcki*(bck)}
     do i = 1, natom

        pot_elc_atm(i) = pot_elc_atm(i) &
             &         + pot_fact_awave * atmchrg(i) &
             &         * (sini(i)*bsk + cosi(i)*bck)

     end do

! --- STEP 4 MAKE csk(:,k) = bsk * Ak * K
!                 cck(:,k) = bck * Ak * K

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
     fac_vir2 = poti

     if (.not. ifnetqcorrp) then        ! no charge correction

        do n=1,3

           k_fact(1:3,n) = -1.0d0*(2.0d0/ktmp2 + fac_vir)   &
                &        * kwave(1:3,k)*kwave(n,k)

           k_fact(n,n) = k_fact(n,n) + 1.0d0

           k_fact(1:3,n) =  k_fact(1:3,n) * fac_vir2

           pot_virit_coul(1:3,n) = pot_virit_coul(1:3,n)   &
                &                + k_fact(1:3,n)

        end do

        !- calculate virial tensor of each atom (_HF_BULK_EWK)
        do i = 1, natom

           fac_vir2 = pot_fact_awave &
                &   * atmchrg(i) * (sini(i)*bsk + cosi(i)*bck)

           do n=1,3
              k_fact(1:3,n) = -1.0d0*(2.0d0/ktmp2 + fac_vir) &
                   &        * kwave(1:3,k)*kwave(n,k)

              k_fact(n,n) = k_fact(n,n) + 1.0d0

              k_fact(1:3,n) = k_fact(1:3,n) * fac_vir2

              virielct_atm(1:3,n,i,1) = virielct_atm(1:3,n,i,1) &
                   &                  + k_fact(1:3,n)

           end do

        end do


     else                               ! charge correction

        fac_vir3 = pot_fact_awave * netchrgsq
        fac_vir4 = fac_vir3 * (-1.0d0 + fac_vir*ktmp2)
        pot_ewk_cc = pot_ewk_cc + fac_vir4

        do n=1,3

           k_fact(1:3,n) = -1.0d0*(2.0d0/ktmp2 + fac_vir) &
                &        * kwave(1:3,k)*kwave(n,k)

           k_fact(n,n) = k_fact(n,n) + 1.0d0

           k_fact(1:3,n) = k_fact(1:3,n) * (fac_vir2 - fac_vir3)

           pot_virit_coul(1:3,n) = pot_virit_coul(1:3,n) &
                &                + k_fact(1:3,n)

        end do

        !- calculate virial tensor of each atom (_HF_BULK_EWK)
        do i = 1, natom

           fac_vir2 = pot_fact_awave &
                &   * atmchrg(i) * (sini(i)*bsk + cosi(i)*bck)

           do n=1,3

              k_fact(1:3,n) = -1.0d0*(2.0d0/ktmp2 + fac_vir)   &
                   &        * kwave(1:3,k)*kwave(n,k)

              k_fact(n,n) = k_fact(n,n) + 1.0d0

              k_fact(1:3,n) = k_fact(1:3,n) * (fac_vir2 - fac_vir3)

              virielct_atm(1:3,n,i,1) = virielct_atm(1:3,n,i,1) &
                   &                  + k_fact(1:3,n)

           end do

        end do

     end if

  END DO

  pot_viri_coul =  pot_viri_coul + pot_ewk + pot_ewk_cc

! --- STEP 7 CALUCLATE FINAL FORCE TO ADD ---
!         fi =  2Qi/V * (fsi - fci)

  DO i = 1, natom
     fi(1:3,i) = 2.0d0 * atmchrg(i) * inv_vol * fi(1:3,i)
  END DO

! --- ADD FSI TO FORCE  ---

  force(1:3,1:natom) = force(1:3,1:natom) + fi(1:3,1:natom)
  for_viri_coul(1:3,1:natom) = for_viri_coul(1:3,1:natom) + fi(1:3,1:natom)

!--- Calculate self-interaction term of each atom (_HF_BULK_EWK)
  do i = 1, natom

     pot_elc_atm(i) = pot_elc_atm(i) - alpha_sqrtpi * atmchrg(i)**2

  end do

!     +     +     +     +     +     +     +

end subroutine ewaldkp_hf
