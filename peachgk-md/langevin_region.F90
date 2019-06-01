!*********************************
!*  langevin_region.f90 Ver.1.2  *
!*      for peachgk_md.f         *
!*            by G.Kikugawa      *
!*********************************
! Time-stamp: <2015-03-19 11:30:01 gota>

subroutine langevin_region(dt_long_cal, &
     &                     iftcratom, &
     &                     nlangeregion, &
     &                     r_ltemp, r_ltdamp, &
     &                     ifoutthc, &
     &                     det_ene_kin, &
     &                     ifhalf)

  use md_global

  implicit none

!     subroutine to control temperature by Langevin thermostat
!       in specific regions described in "langecont.ini"
!
!--------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!! Caution !!!!!!!!!!!!!!!!!!!!!!!!!
!  In this version, time descritization scheme of Langevin equation
!  with holonomic constraint is not accurate one.
!  make sure temperature constancy if some molecules have contraint bonding.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ARGUMENTS:
!     INPUT
  real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]

  logical,intent(in):: iftcratom    ! region temp. control based on atom or mole.

  integer,intent(in):: nlangeregion    ! number of region for Langevin thermo.

  real(8),intent(in):: r_ltemp(:)      ! control temp. in each region
  real(8),intent(in):: r_ltdamp(:)     ! damping factor in each region [non-d]

  logical,intent(in):: ifoutthc     ! flag for outputting thermal control file

  logical,intent(in):: ifhalf   ! flag for former or latter half of integration
                                ! former=.true. or latter=.false.

!     OUTPUT
  real(8),intent(out):: det_ene_kin(:) ! for outputting thermal control data

! LOCAL:
  integer:: i,j             ! do loop index
  integer:: j1,j2           ! do loop index
  integer:: ipolytyp
  integer:: ipoly
  integer:: iwater
  integer:: imatyp
  integer:: atm_index
  integer:: atm_index_tmp

  integer:: iltr

  real(8):: ene_kin(nlangeregion) ! kinetic energy in region
  real(8):: ene_kin_new(nlangeregion) ! kinetic energy in region (updated)

  real(8):: atmvel_tmp(3,natom) ! temporary atom velocity without streaming vel.

  real(8):: var_randf             ! variance of random force
  real(8):: var_randf_tmp

  real(8):: u                     ! random number from SFMT

  real(8):: pi

! FUNCTIONS:
  real(8):: genrand_res53        ! function creating random number

!     +     +     +     +     +     +     +

!---- some preparation

  atmvel_tmp(1:3,1:natom) = atmvel(1:3,1:natom)
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) - atmvel_strm(1:3,1:natom)

  atm_index = 0

  pi = acos(-1.0d0)

!-------- region-based Langevin thermostat --------

  !- calculate imposed (or extracted) kinetic energy
  if (ifoutthc) then
     ene_kin(1:nlangeregion) = 0.0d0

     do iltr = 1, nlangeregion
        do i = 1, natom_lt(iltr)
           atm_index = atmindex_lt(i,iltr)
           ene_kin(iltr) = ene_kin(iltr) &
                &        + atmmass(atm_index) &
                &        * (atmvel(1,atm_index)**2 &
                &         + atmvel(2,atm_index)**2 &
                &         + atmvel(3,atm_index)**2)
        end do
     end do

  end if

  !- calculate random force by Box-Muller method
  if (ifhalf) then    ! only former half of integration

     frand_lt(1:3,1:natom) = 0.0d0

     do iltr = 1, nlangeregion

        var_randf_tmp = 2.0d0 * r_ltdamp(iltr) * r_ltemp(iltr) / dt_long_cal

        do i = 1, natom_lt(iltr)
           atm_index = atmindex_lt(i,iltr)

           var_randf = var_randf_tmp * atmmass(atm_index)   ! variance

           u = genrand_res53()
           frand_lt(1,atm_index) = sqrt(-2.0d0*log(u)*var_randf)
           frand_lt(2,atm_index) = sqrt(-2.0d0*log(u)*var_randf)
           u = genrand_res53()
           frand_lt(1,atm_index) = frand_lt(1,atm_index) * cos(2.0d0*pi*u)
           frand_lt(2,atm_index) = frand_lt(2,atm_index) * sin(2.0d0*pi*u)
           u = genrand_res53()
           frand_lt(3,atm_index) = sqrt(-2.0d0*log(u)*var_randf)
           u = genrand_res53()
           frand_lt(3,atm_index) = frand_lt(3,atm_index) * cos(2.0d0*pi*u)

           ! acceleration
           frand_lt(1:3,atm_index) = frand_lt(1:3,atm_index) &
                &                  / atmmass(atm_index) 

        end do

     end do
  end if

!---- time integration of Langevin equation of motion
  if (ifhalf) then   ! former half of integration
     do iltr = 1, nlangeregion

        do i = 1, natom_lt(iltr)
           atm_index = atmindex_lt(i,iltr)

           ! velocity update
           atmvel(1:3,atm_index) = atmvel(1:3,atm_index) &
                &                + 0.5d0*dt_long_cal * frand_lt(1:3,atm_index)
                                                                ! random term
           atmvel(1:3,atm_index) = atmvel(1:3,atm_index) &
                &                * exp(-0.5d0*dt_long_cal * r_ltdamp(iltr))
                                                                ! friction term
        end do

     end do

  else   ! latter half of integration

     do iltr = 1, nlangeregion

        do i = 1, natom_lt(iltr)
           atm_index = atmindex_lt(i,iltr)

           ! velocity update
           atmvel(1:3,atm_index) = atmvel(1:3,atm_index) &
                &                * exp(-0.5d0*dt_long_cal * r_ltdamp(iltr))
                                                                ! friction term
           atmvel(1:3,atm_index) = atmvel(1:3,atm_index) &
                &                + 0.5d0*dt_long_cal * frand_lt(1:3,atm_index)
                                                                ! random term
        end do

     end do

  end if

!---- Calculate imposed (or extracted) kinetic energy
  if (ifoutthc) then
     ene_kin_new(1:nlangeregion) = 0.0d0

     do iltr = 1, nlangeregion
        do i = 1, natom_lt(iltr)
           atm_index = atmindex_lt(i,iltr)
           ene_kin_new(iltr) = ene_kin_new(iltr) &
                &        + atmmass(atm_index) &
                &        * (atmvel(1,atm_index)**2 &
                &         + atmvel(2,atm_index)**2 &
                &         + atmvel(3,atm_index)**2)
        end do

        det_ene_kin(iltr) = det_ene_kin(iltr) &
             &            + 0.5d0 * (ene_kin_new(iltr) - ene_kin(iltr))
     end do

  end if

!--- reset velocity
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) + atmvel_strm(1:3,1:natom)

!     +     +     +     +     +     +     +

end subroutine langevin_region
