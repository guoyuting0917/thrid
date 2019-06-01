!*****************************
!*  output_for.f90 Ver.1.0   *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>

!----------------------------------------------------------
subroutine outfor(oufor,current_step,fref, &
    &             for_long,for_short)

  use md_global

  implicit none

!     subroutine for outputting force-related data

!ARGUMENTS:
!     INPUT
  integer,intent(in):: oufor           ! output unit for output force data

  integer,intent(in):: current_step    ! current time step

!---- base value for non-dimensionalize
  real(8),intent(in):: fref             ! force base value [N]

!---- force of each atoms
  real(8),intent(in):: for_long(:,:)     ! long-range force
!  real(8),intent(in):: for_med(:,:)   ! medium-range force
  real(8),intent(in):: for_short(:,:)    ! short-range force

! LOCAL:
  integer:: i,m             ! do loop index

#if defined(_DOUBLE_OUTFOR)
  real(8),allocatable:: force_out(:,:)   ! output type for double precision
#else
  real,allocatable:: force_out(:,:)   ! single precision is default
#endif

!     +     +     +     +     +     +     +     +

  !---- allocalte memory
  allocate(force_out(3,natom))

  !---- force are copied and summed to arrays for output

  do i = 1, natom
      force_out(1:3,i) = (for_long(1:3,i)+for_short(1:3,i)) &
      &                      * atmmass(i) * fref
  end do

  !---- output velocity
  write(oufor) current_step

  do i=1,natom
     write(oufor) i,(force_out(m,i),m=1,3),atmtyp(i)
  end do

  !---- release memory
  deallocate(force_out)

!     +     +     +     +     +     +     +     +

end subroutine outfor

!----------------------------------------------------------
#if defined(HF)
subroutine outfor_hf(oufor,current_step,fref, &
    &                for_long,for_short, &
    &                ifhfvol, &
    &                nhfregion,hfzpos1,hfzpos2, &
    &                hftyp_atm, &
    &                molecom, &
    &                viribot_long_atm,viribot_short_atm, &
    &                viriant_long_atm,viriant_short_atm, &
    &                viritot_long_atm,viritot_short_atm, &
    &                viri14t_long_atm,viri14t_short_atm, &
    &                virielint_long_atm,virielint_short_atm, &
    &                viriljint_long_atm,viriljint_short_atm, &
    &                virielct_long_atm,virielct_short_atm, &
    &                viriljt_long_atm,viriljt_short_atm, &
    &                virimort_long_atm,virimort_short_atm, &
    &                virisht_long_atm,virisht_short_atm, &
    &                virirfht_long_atm,virirfht_short_atm, &
    &                viridout_long_atm,viridout_short_atm, &
    &                viricstmnbt_long_atm,viricstmnbt_short_atm, &
    &                viricstmnbext_long_atm,viricstmnbext_short_atm)

  use md_global

  implicit none

!     subroutine for outputting force- or heat flux- related data

!ARGUMENTS:
!     INPUT
  integer,intent(in):: oufor           ! output unit for output force data

  integer,intent(in):: current_step    ! current time step

!---- base value for non-dimensionalize
  real(8),intent(in):: fref             ! force base value [N]

!---- force of each atoms
  real(8),intent(in):: for_long(:,:)     ! long-range force
!  real(8),intent(in):: for_med(:,:)   ! medium-range force
  real(8),intent(in):: for_short(:,:)    ! short-range force

!---- variables for calculation of heat flux
  logical,intent(in):: ifhfvol        ! local volume-based or local surface-based

  integer,intent(in):: nhfregion       ! number of region to calculate heat flux
  real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                ! z-position of region for heat flux

  integer,intent(in):: hftyp_atm(:)    ! atom- or mole-based heat flux cal.
                                !   for each atom

  real(8),intent(in):: molecom(:,:)     ! center of mass of molecule

!     virial tensor for each atom
  real(8),intent(in):: viribot_long_atm(:,:,:,:) ! long-range (bond)
  real(8),intent(in):: viribot_short_atm(:,:,:,:) ! short-range (bond)
  real(8),intent(in):: viriant_long_atm(:,:,:,:) ! long-range (angle)
  real(8),intent(in):: viriant_short_atm(:,:,:,:) ! short-range (angle)
  real(8),intent(in):: viritot_long_atm(:,:,:,:) ! long-range (torsion)
  real(8),intent(in):: viritot_short_atm(:,:,:,:) ! short-range (torsion)
  real(8),intent(in):: viri14t_long_atm(:,:,:,:) ! long-range (1-4)
  real(8),intent(in):: viri14t_short_atm(:,:,:,:) ! short-range (1-4)
  real(8),intent(in):: virielint_long_atm(:,:,:,:) ! long-range (elc intra)
  real(8),intent(in):: virielint_short_atm(:,:,:,:) ! short-range (elc intra)
  real(8),intent(in):: viriljint_long_atm(:,:,:,:) ! long-range (L-J intra)
  real(8),intent(in):: viriljint_short_atm(:,:,:,:) ! short-range (L-J intra)
  real(8),intent(in):: virielct_long_atm(:,:,:,:) ! long-range (Coulomb)
  real(8),intent(in):: virielct_short_atm(:,:,:,:) ! short-range (Coulomb)
  real(8),intent(in):: viriljt_long_atm(:,:,:,:) ! long-range (L-J)
  real(8),intent(in):: viriljt_short_atm(:,:,:,:) ! short-range (L-J)
  real(8),intent(in):: virimort_long_atm(:,:,:,:) ! long-range (Morse)
  real(8),intent(in):: virimort_short_atm(:,:,:,:) ! short-range (Morse)
  real(8),intent(in):: virisht_long_atm(:,:,:,:) ! long-range (SH)
  real(8),intent(in):: virisht_short_atm(:,:,:,:) ! short-range (SH)
  real(8),intent(in):: virirfht_long_atm(:,:,:,:) ! long-range (RFH)
  real(8),intent(in):: virirfht_short_atm(:,:,:,:) ! short-range (RFH)
  real(8),intent(in):: viridout_long_atm(:,:,:,:) ! long-range (DOU)
  real(8),intent(in):: viridout_short_atm(:,:,:,:) ! short-range (DOU)
  real(8),intent(in):: viricstmnbt_long_atm(:,:,:,:) ! long-range (custom NB)
  real(8),intent(in):: viricstmnbt_short_atm(:,:,:,:) ! short-range (custom NB)
  real(8),intent(in):: viricstmnbext_long_atm(:,:,:,:,:)
                                                  ! long-range (custom NB)
  real(8),intent(in):: viricstmnbext_short_atm(:,:,:,:,:)
                                                  ! short-range (custom NB)

! LOCAL:
  integer:: i,m             ! do loop index
  integer:: ihfr            ! do loop index

#if defined(_DOUBLE_OUTFOR)
  real(8),allocatable:: force_out(:,:,:)   ! output type for double precision
#else
  real,allocatable:: force_out(:,:,:)   ! single precision is default
#endif

!     +     +     +     +     +     +     +     +

  !---- check ifhfvol
#if defined(_OUT_HFFOR)
  if (ifhfvol) then
      write(6,*) 'Error: outfor for heat flux is only compatible'
      write(6,*) '       with surface-based heat flux.'
      stop
  end if
#endif

  !---- allocalte memory
#if defined(_OUT_HFFOR)
! output pointed force
  allocate(force_out(3,natom,nhfregion))
#else
! normal output (just force)
  allocate(force_out(3,natom,1))
#endif

  !---- force are copied and summed to arrays for output

#if defined(_OUT_HFFOR)
! output pointed force
  do ihfr = 1, nhfregion
     do i = 1, natom
        force_out(1:3,i,ihfr) = fref &
        & * (viribot_long_atm(1:3,3,i,ihfr) + viribot_short_atm(1:3,3,i,ihfr) &
        &  + viriant_long_atm(1:3,3,i,ihfr) + viriant_short_atm(1:3,3,i,ihfr) &
        &  + viritot_long_atm(1:3,3,i,ihfr) + viritot_short_atm(1:3,3,i,ihfr) &
        &  + viri14t_long_atm(1:3,3,i,ihfr) + viri14t_short_atm(1:3,3,i,ihfr) &
        &  + virielint_long_atm(1:3,3,i,ihfr) &
        &  + virielint_short_atm(1:3,3,i,ihfr) &
        &  + viriljint_long_atm(1:3,3,i,ihfr) &
        &  + viriljint_short_atm(1:3,3,i,ihfr) &
        &  + virielct_long_atm(1:3,3,i,ihfr) &
        &  + virielct_short_atm(1:3,3,i,ihfr) &
        &  + viriljt_long_atm(1:3,3,i,ihfr) + viriljt_short_atm(1:3,3,i,ihfr) &
        &  + virimort_long_atm(1:3,3,i,ihfr) &
        &  + virimort_short_atm(1:3,3,i,ihfr) &
        &  + virisht_long_atm(1:3,3,i,ihfr) + virisht_short_atm(1:3,3,i,ihfr) &
        &  + virirfht_long_atm(1:3,3,i,ihfr) &
        &  + virirfht_short_atm(1:3,3,i,ihfr) &
        &  + viridout_long_atm(1:3,3,i,ihfr) &
        &  + viridout_short_atm(1:3,3,i,ihfr) &
        &  + viricstmnbt_long_atm(1:3,3,i,ihfr) &
        &  + viricstmnbt_short_atm(1:3,3,i,ihfr))
     end do
  end do

#else
! normal output (just force)
  do i = 1, natom
      force_out(1:3,i,1) = (for_long(1:3,i)+for_short(1:3,i)) &
      &                      * atmmass(i) * fref
  end do
#endif

  !---- output velocity
  write(oufor) current_step

#if defined(_OUT_HFFOR)
! output pointed force
  do ihfr = 1, nhfregion
     do i = 1, natom
        write(oufor) i,(force_out(m,i,ihfr),m=1,3),atmtyp(i)
     end do
  end do

#else
! normal output (just force)
  do i=1,natom
     write(oufor) i,(force_out(m,i,1),m=1,3),atmtyp(i)
  end do
#endif

  !---- release memory
  deallocate(force_out)

!     +     +     +     +     +     +     +     +

end subroutine outfor_hf
#endif
!----------------------------------------------------------
