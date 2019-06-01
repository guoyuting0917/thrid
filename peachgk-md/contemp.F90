!*****************************
!*  contemp.f90 Ver.2.5      *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-02-17 18:55:44 gota>

subroutine contemp(npoly,npolytyp,npoly_mole, &
     &             nwater, &
     &             nmatom,nmatyp,nmatomtyp, &
     &             degfree_poly,degfree_water, &
     &             degfree_ma,degfree_all, &
     &             tcont_polyt,tcont_watert,tcont_mat, &
     &             tcont_polyinit,tcont_waterinit,tcont_mainit, &
     &             tfactor_poly,tfactor_water,tfactor_ma)

  use md_global

  implicit none

!     subroutine to control temperature (for all system) by WOODCOCK's method
!     !!! The temperature is set to the value of tcont_polyt(1) !!!
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

  real(8):: atmvel_tmp(3,natom) ! temporary atom velocity without streaming vel.

!     +     +     +     +     +     +     +

! - some preparation -

  atmvel_tmp(1:3,1:natom) = atmvel(1:3,1:natom)
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) - atmvel_strm(1:3,1:natom)

!-------- temperature control to all --------

!     -- CALCULATE KINETIC ENERGY --

  ene_kin = 0.0d0

  DO i = 1, npoly
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1, j2
        jj = molept_list(j)
        veloc2 = atmvel(1,jj)**2 + atmvel(2,jj)**2 + atmvel(3,jj)**2
        ene_kin = ene_kin + atmmass(jj) * veloc2
     end do
  END DO

  DO i = 1+npoly, npoly+nwater
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1, j2
        jj = molept_list(j)
        veloc2 = atmvel(1,jj)**2 + atmvel(2,jj)**2 + atmvel(3,jj)**2
        ene_kin = ene_kin + atmmass(jj) * veloc2
     end do
  END DO

  DO i = 1+npoly+nwater, npoly+nwater+nmatom
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1, j2
        jj = molept_list(j)
        veloc2 = atmvel(1,jj)**2 + atmvel(2,jj)**2 + atmvel(3,jj)**2
        ene_kin = ene_kin + atmmass(jj) * veloc2
     end do
  END DO

!     -- CALCULATE TEMPERATURE ---      

  if (degfree_all /= 0) then
     factor = 1.0d0 / dble(degfree_all)   ! ene_kin = 2.0 * K
  else
     factor = 0.0d0
  end if
  real_temp = factor*ene_kin

!     -- RESET VELOCITY --
!!! tfactor_poly is representative for tfactor_all in this version
  if (real_temp > 1.0d-5) then
     sfact = sqrt((tfactor_poly(1)  &
          &     * (tcont_polyt(1)-tcont_polyinit(1)) &
          &      + tcont_polyinit(1)) / real_temp)
  else
     sfact = 1.0d0
  end if

  DO i = 1, npoly
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1, j2
        jj = molept_list(j)

        atmvel(1:3,jj) = 1.0d0 * (sfact-1.0d0) * atmvel(1:3,jj) &
             &         + atmvel(1:3,jj)
     end do
  END DO

  DO i = 1+npoly, npoly+nwater
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1, j2
        jj = molept_list(j)

        atmvel(1:3,jj) = 1.0d0 * (sfact-1.0d0) * atmvel(1:3,jj) &
             &         + atmvel(1:3,jj)

     end do
  END DO

  DO i = 1+npoly+nwater, npoly+nwater+nmatom
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1, j2
        jj = molept_list(j)

        atmvel(1:3,jj) = 1.0d0 * (sfact-1.0d0) * atmvel(1:3,jj) &
             &         + atmvel(1:3,jj)
 
     end do
  END DO

!--- reset velocity
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) + atmvel_strm(1:3,1:natom)

!     +     +     +     +     +     +     +

end subroutine contemp
