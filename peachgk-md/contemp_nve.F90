!*****************************
!*  contemp_nve.f90 Ver.1.4  *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-02-17 19:00:57 gota>

subroutine contemp_nve(npoly,npolytyp,npoly_mole,npoly_atom, &
     &                 nwater, &
     &                 nmatom,nmatyp,nmatomtyp, &
     &                 degfree_poly,degfree_water, &
     &                 degfree_ma,degfree_all, &
     &                 tcont_polyt,tcont_watert,tcont_mat, &
     &                 tcont_polyinit,tcont_waterinit, &
     &                 tcont_mainit, &
     &                 tfactor_poly,tfactor_water,tfactor_ma, &
     &                 nlheat_poly,index_nlheat_poly, &
     &                 nlheat_water, &
     &                 nlheat_ma,index_nlheat_ma)

  use md_global

  implicit none

!     subroutine to control temperature by WOODCOCK's method
!       for specific molecules
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
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

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

  integer,intent(in):: nlheat_poly     ! number of poly type for local heating
  integer,intent(in):: index_nlheat_poly(:) 
                                ! index of poly type for local heating
  integer,intent(in):: nlheat_water    ! number of water for local heating
  integer,intent(in):: nlheat_ma       ! number of matom type for local heating
  integer,intent(in):: index_nlheat_ma(:) 
                                ! index of matom type for local heating

! LOCAL:       

  real(8):: factor           ! factor to multipy to Ekin 
  real(8):: veloc2           ! |veloc|**2                 
  real(8):: ene_kin          ! kinetic energy                 
  real(8):: real_temp        ! real temperature               
  real(8):: sfact            ! factor to multiply velocity   
  integer:: i,j             ! do loop index
  integer:: j1,j2
  integer:: jj

  integer:: ipolytyp
  integer:: imatyp
      
  integer:: index_atom, index_atom2
  integer:: index_polytyp
  integer:: index_matyp

  integer:: npolyall
  integer:: npolywaterall

  real(8):: atmvel_tmp(3,natom) ! temporary atom velocity without streaming vel.

!     +     +     +     +     +     +     +

!---- some preparation

  atmvel_tmp(1:3,1:natom) = atmvel(1:3,1:natom)
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) - atmvel_strm(1:3,1:natom)

  npolyall = 0
  do i = 1, npolytyp
     npolyall = npolyall + npoly_mole(i)*npoly_atom(i)
  end do
  npolywaterall = npolyall + nwater*3

!-------- temperature control about polymer1 --------

!     -- CALCULATE KINETIC ENERGY --

  DO ipolytyp = 1, nlheat_poly

     ene_kin = 0.0d0
     index_atom = 0
     index_polytyp = index_nlheat_poly(ipolytyp)

     do i = 1, index_polytyp - 1
        index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
     end do
     index_atom2 = index_atom

     do i = 1, npoly_mole(index_polytyp)
        do j = 1, npoly_atom(index_polytyp)
           index_atom = index_atom + 1
           veloc2 =  atmvel(1,index_atom)**2  &
                &  + atmvel(2,index_atom)**2 &
                &  + atmvel(3,index_atom)**2
           ene_kin = ene_kin + atmmass(index_atom) * veloc2
        end do
     end do

!     -- CALCULATE TEMPERATURE ---      

     if (degfree_poly(index_polytyp) .ne. 0) then
        factor = 1.0d0 / dble(degfree_poly(index_polytyp))
                                                           ! ene_kin = 2.0 * K
     else
        factor = 0.0d0
     end if
     real_temp = factor*ene_kin

!     -- RESET VELOCITY --

     if (real_temp > 1.0d-5) then
        sfact = sqrt((tfactor_poly(index_polytyp)  &
             &     * (tcont_polyt(index_polytyp) &
             &       - tcont_polyinit(index_polytyp)) &
             &      + tcont_polyinit(index_polytyp)) / real_temp)
     else
        sfact = 1.0d0
     end if

     do i = 1, npoly_mole(index_polytyp)

        do j = 1, npoly_atom(index_polytyp)
           index_atom2 = index_atom2 + 1
           atmvel(1:3,index_atom2) = 1.0d0 * (sfact-1.0d0)   &
                &                  * atmvel(1:3,index_atom2) &
                &                  + atmvel(1:3,index_atom2)
        end do
     end do

  END DO

!-------- temperature control about water --------

!     -- CALCULATE KINETIC ENERGY --

  ene_kin = 0.0d0

  if (nlheat_water == 1) then

     DO i = 1+npoly, npoly+nwater
        j1 = molept_index(i)
        j2 = molept_index(i+1) - 1
        do j = j1, j2
           jj = molept_list(j)
           veloc2 = atmvel(1,jj)**2 + atmvel(2,jj)**2 &
                & + atmvel(3,jj)**2
           ene_kin = ene_kin + atmmass(jj) * veloc2
        end do
     END DO

!     -- CALCULATE TEMPERATURE ---      

     if (degfree_water /= 0) then
        factor = 1.0d0 / dble(degfree_water)   ! ene_kin = 2.0 * K
     else
        factor = 0.0d0
     end if
     real_temp = factor*ene_kin

!     -- RESET VELOCITY --

     if (real_temp > 1.0d-5) then
        sfact = sqrt((tfactor_water  &
             &     * (tcont_watert-tcont_waterinit) &
             &      + tcont_waterinit) / real_temp)
     else
        sfact = 1.0d0
     end if

     DO i = 1+npoly, npoly+nwater
        j1 = molept_index(i)
        j2 = molept_index(i+1) - 1
        do j = j1, j2
           jj = molept_list(j)

           atmvel(1:3,jj) = 1.0d0 * (sfact-1.0d0) * atmvel(1:3,jj) &
                &         + atmvel(1:3,jj)

        end do
     END DO

  end if

!-------- temperature control about MATOM --------

!     -- CALCULATE KINETIC ENERGY --

  DO imatyp = 1, nlheat_ma

     ene_kin = 0.0d0
     index_atom = npolywaterall
     index_matyp = index_nlheat_ma(imatyp)

     do i = 1, index_matyp-1
        index_atom = index_atom + nmatomtyp(i)
     end do
     index_atom2 = index_atom

     do i = 1, nmatomtyp(index_matyp)
        index_atom = index_atom + 1
        veloc2 = atmvel(1,index_atom)**2 + atmvel(2,index_atom)**2 &
             & + atmvel(3,index_atom)**2
        ene_kin = ene_kin + atmmass(index_atom) * veloc2
     end do

!     -- CALCULATE TEMPERATURE ---      

     if (degfree_ma(index_matyp) /= 0) then
        factor = 1.0d0 / dble(degfree_ma(index_matyp))   ! ene_kin = 2.0 * K
     else
        factor = 0.0d0
     end if
     real_temp = factor*ene_kin

!     -- RESET VELOCITY --

     if (real_temp > 1.0d-5) then
        sfact = sqrt((tfactor_ma(index_matyp) &
             &      * (tcont_mat(index_matyp) &
             &        -tcont_mainit(index_matyp)) &
             &      + tcont_mainit(index_matyp)) / real_temp)
     else
        sfact = 1.0d0
     end if

     do i = 1, nmatomtyp(index_matyp)
        index_atom2 = index_atom2 + 1
        atmvel(1:3,index_atom2) = 1.0d0 * (sfact-1.0d0) &
             &                  * atmvel(1:3,index_atom2) &
             &                  + atmvel(1:3,index_atom2)
     end do

  END DO

!--- reset velocity
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) + atmvel_strm(1:3,1:natom)

!     +     +     +     +     +     +     +

end subroutine contemp_nve
