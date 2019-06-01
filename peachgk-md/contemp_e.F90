!*****************************
!*  contemp_e.f90 Ver.2.5    *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-02-17 18:58:39 gota>

subroutine contemp_e(npoly,npolytyp,npoly_mole, &
     &               nwater, &
     &               nmatom,nmatyp,nmatomtyp, &
     &               degfree_poly,degfree_water, &
     &               degfree_ma,degfree_all, &
     &               tcont_polyt,tcont_watert,tcont_mat, &
     &               tcont_polyinit,tcont_waterinit,tcont_mainit, &
     &               tfactor_poly,tfactor_water,tfactor_ma)

  use md_global

  implicit none

!     subroutine to control temperature by WOODCOCK's method
!
!     temperature T = 2<K>/(Df.R)
!     where <K> is the total kinetic energy, 
!           Df  is the degree of freedom
!--------------------------------------------------------

! ARGUMENTS:
!     INPUT
  integer,intent(in):: npoly           ! number of polymer1
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatom          ! number of monatomic molecules
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  integer,intent(in):: degfree_poly(:) ! degree of freedom of polymer1
  integer,intent(in):: degfree_water   ! degree of freedom of H2O
  integer,intent(in):: degfree_ma(:)   ! degree of freedom of monatomic mole.
  integer,intent(in):: degfree_all     ! degree of freedom of all atoms

  real(8),intent(in):: tcont_polyt(:)   ! poly Temp. [non-d] in NVT
  real(8),intent(in):: tcont_watert     ! H2O Temp. [non-d] in NVT
  real(8),intent(in):: tcont_mat(:)     ! monatomic mole. Temp. [non-d] in NVT

  real(8),intent(in):: tcont_polyinit(:) ! poly Temp. [non-d] in NVT
  real(8),intent(in):: tcont_waterinit  ! H2O Temp. [non-d] in NVT
  real(8),intent(in):: tcont_mainit(:)  ! MA Temp. [non-d] in NVT

  real(8),intent(in):: tfactor_poly(:)  ! temperature control factor of polymer1
  real(8),intent(in):: tfactor_water    ! temperature control factor of H2O
  real(8),intent(in):: tfactor_ma(:)    ! temperature control factor of MA

! LOCAL:       

  real(8):: factor           ! factor to multipy to Ekin 
  real(8):: veloc2           ! |veloc|**2                 
  real(8):: ene_kin          ! kinetic energy                 
  real(8):: real_temp        ! real temperature               
  real(8):: sfact            ! factor to multiply velocity   
  integer:: i,j             ! do loop index
  integer:: j1,j2
  integer:: jj
  integer:: imole1,imole2
  integer:: ipolytyp
  integer:: imatyp

  real(8):: atmvel_tmp(3,natom) ! temporary atom velocity without streaming vel.

!     +     +     +     +     +     +     +

! - some preparation -

  atmvel_tmp(1:3,1:natom) = atmvel(1:3,1:natom)
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) - atmvel_strm(1:3,1:natom)

!-------- temperature control about polymer1 --------

!     -- CALCULATE KINETIC ENERGY --

  imole1 = 0
  imole2 = 0

  DO ipolytyp = 1, npolytyp

     ene_kin = 0.0d0

     DO i = 1, npoly_mole(ipolytyp)
        imole1 = imole1 + 1
        j1 = molept_index(imole1)
        j2 = molept_index(imole1+1) - 1
        do j = j1, j2
           jj = molept_list(j)
           veloc2 =  atmvel(1,jj)**2 + atmvel(2,jj)**2 &
                &  + atmvel(3,jj)**2
           ene_kin = ene_kin + atmmass(jj) * veloc2
        end do
     END DO

!     -- CALCULATE TEMPERATURE ---      

     if (degfree_poly(ipolytyp) /= 0) then
        factor = 1.0d0 / dble(degfree_poly(ipolytyp))   ! ene_kin = 2.0 * K
     else
        factor = 0.0d0
     end if
     real_temp = factor*ene_kin

!     -- RESET VELOCITY --

     if (real_temp > 1.0d-5) then
        sfact = sqrt((tfactor_poly(ipolytyp)  &
            & * (tcont_polyt(ipolytyp)-tcont_polyinit(ipolytyp)) &
            & + tcont_polyinit(ipolytyp)) / real_temp)
     else
        sfact = 1.0d0
     end if

     DO i = 1, npoly_mole(ipolytyp)
        imole2 = imole2 + 1
        j1 = molept_index(imole2)
        j2 = molept_index(imole2+1) - 1
        do j = j1, j2
           jj = molept_list(j)

           atmvel(1:3,jj) = 1.0d0 * (sfact-1.0d0) * atmvel(1:3,jj) &
                &         + atmvel(1:3,jj)
        end do
     END DO

  END DO

!-------- temperature control about water --------

!     -- CALCULATE KINETIC ENERGY --

  ene_kin = 0.0d0

  DO i = 1+npoly, npoly+nwater
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1, j2
        jj = molept_list(j)
        veloc2 = atmvel(1,jj)**2 + atmvel(2,jj)**2 + atmvel(3,jj)**2
        ene_kin = ene_kin + atmmass(jj) * veloc2
     end do
  END DO

!     -- CALCULATE TEMPERATURE ---      

  if (degfree_water .ne. 0) then
     factor = 1.0d0 / dble(degfree_water)   ! ene_kin = 2.0 * K
  else
     factor = 0.0d0
  end if
  real_temp = factor*ene_kin

!     -- RESET VELOCITY --

  if (real_temp > 1.0d-5) then
     sfact = sqrt((tfactor_water * (tcont_watert-tcont_waterinit) &
          &       + tcont_waterinit) / real_temp)
  else
     sfact = 1.0d0
  end if

  DO i = 1+npoly, npoly+nwater
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1, j2
        jj = molept_list(j)

        atmvel(1:3,jj) = 1.0d0 * (sfact-1.0d0) * atmvel(1:3,jj) &
             &          + atmvel(1:3,jj)

     end do
  END DO

!-------- temperature control about MATOM --------

!     -- CALCULATE KINETIC ENERGY --

  imole1 = 0
  imole2 = 0

  DO imatyp = 1, nmatyp

     ene_kin = 0.0d0

     DO i = 1, nmatomtyp(imatyp)
        imole1 = imole1 + 1
        j1 = molept_index(npoly+nwater+imole1)
        j2 = molept_index(npoly+nwater+imole1+1) - 1
        do j = j1, j2
           jj = molept_list(j)
           veloc2 = atmvel(1,jj)**2 + atmvel(2,jj)**2 &
                & + atmvel(3,jj)**2
           ene_kin = ene_kin + atmmass(jj) * veloc2
        end do
     END DO

!     -- CALCULATE TEMPERATURE ---      

     if (degfree_ma(imatyp) .ne. 0) then
        factor = 1.0d0 / dble(degfree_ma(imatyp))   ! ene_kin = 2.0 * K
     else
        factor = 0.0d0
     end if
     real_temp = factor*ene_kin

!     -- RESET VELOCITY --

     if (real_temp > 1.0d-5) then
        sfact = sqrt((tfactor_ma(imatyp)  &
             &      * (tcont_mat(imatyp)-tcont_mainit(imatyp)) &
             &      + tcont_mainit(imatyp)) / real_temp)
     else
        sfact = 1.0d0
     end if

     DO i = 1, nmatomtyp(imatyp)
        imole2 = imole2 + 1
        j1 = molept_index(npoly+nwater+imole2)
        j2 = molept_index(npoly+nwater+imole2+1) - 1
        do j = j1, j2
           jj = molept_list(j)

           atmvel(1:3,jj) = 1.0d0 * (sfact-1.0d0) * atmvel(1:3,jj) &
                &         + atmvel(1:3,jj)

        end do
     END DO

  END DO

!--- reset velocity
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) + atmvel_strm(1:3,1:natom)

!     +     +     +     +     +     +     +

end subroutine contemp_e
