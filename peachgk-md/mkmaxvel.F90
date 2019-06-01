!**************************************
!*  mkmaxvel_.f Ver.1.9 '10.04.12     *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine mkmaxvel(npoly, nwater, nmatom, &
     &              tcont_poly, tcont_water, tcont_ma, tempref)

  use md_global

#if defined(MPI)
  use mpi_global
#endif

  implicit none

!     This subroutine is that Maxwell-Boltzmann distribution are given to
!       all molecules.
!     !!! The temperature is set to the value of tcont_poly(1),
!         tcont_water, and tcont_ma(1). !!!
!
! ARGUMENT:
!     INPUT
  integer,intent(in):: npoly           ! number of polymer1
  integer,intent(in):: nwater          ! number of H2O molecules
  integer,intent(in):: nmatom          ! number of monatomic molecules

  real(8),intent(in):: tcont_poly(:)    ! poly Temp. [non-d] in NVT
  real(8),intent(in):: tcont_water      ! H2O Temp. [non-d] in NVT
  real(8),intent(in):: tcont_ma(:)      ! monatomic mole. Temp. [non-d] in NVT

  real(8),intent(in):: tempref          ! temperature base value [K]

! LOCAL:
  real(8):: u                ! temporary valiables
!      real*8:: v                ! temporary valiables
  real(8):: vcal,vmp         ! temporary valiables      
  real(8):: pi               ! =3.14159265...
  integer:: i,j,k
  integer:: j1,j2

! FUNCTIONS:
  real(8):: genrand_res53        ! function creating random number

!     +     +     +     +     +     +     +

#if defined(MPI)
  if (irank == 0) then
#endif
     write(6,*) 'Give Maxwell-Boltzmann distribution: tcont= ', &
          &     tcont_poly(1)*tempref
#if defined(MPI)
  end if
#endif

!---- some preparation
  pi = dacos(-1.0d0)

!-------- Velocity of Polymer1 --------

  do i = 1,npoly
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1,j2
        k = molept_list(j)
 
10      u = genrand_res53()
        atmvel(1,k) = dsqrt(-2.0d0*tcont_poly(1)/atmmass(k)*dlog(u))
        atmvel(2,k) = dsqrt(-2.0d0*tcont_poly(1)/atmmass(k)*dlog(u))
        u = genrand_res53()
        atmvel(1,k) = atmvel(1,k) * dcos(2.0d0*pi*u)
        atmvel(2,k) = atmvel(2,k) * dsin(2.0d0*pi*u)

        u = genrand_res53()
        atmvel(3,k) = dsqrt(-2.0d0*tcont_poly(1)/atmmass(k)*dlog(u))
        u = genrand_res53()
        atmvel(3,k) = atmvel(3,k) * dcos(2.0d0*pi*u)

        vcal = atmvel(1,k)**2 + atmvel(2,k)**2 + atmvel(3,k)**2
        vmp  = 2.0d0*tcont_poly(1)/atmmass(k)*3.5d0**2
        if (vcal > vmp) goto 10

     end do
  end do

!-------- Velocity of H2O --------

  do i = npoly+1, npoly+nwater
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1,j2
        k = molept_list(j)

20      u = genrand_res53()
        atmvel(1,k) = dsqrt(-2.0d0*tcont_water/atmmass(k)*dlog(u))
        atmvel(2,k) = dsqrt(-2.0d0*tcont_water/atmmass(k)*dlog(u))
        u = genrand_res53()
        atmvel(1,k) = atmvel(1,k) * dcos(2.0d0*pi*u)
        atmvel(2,k) = atmvel(2,k) * dsin(2.0d0*pi*u)

        u = genrand_res53()
        atmvel(3,k) = dsqrt(-2.0d0*tcont_water/atmmass(k)*dlog(u))
        u = genrand_res53()
        atmvel(3,k) = atmvel(3,k) * dcos(2.0d0*pi*u)

        vcal = atmvel(1,k)**2 + atmvel(2,k)**2 + atmvel(3,k)**2
        vmp  = 2.0d0*tcont_water/atmmass(k)*3.5d0**2
        if (vcal > vmp) goto 20

     end do
  end do

!-------- Velocity of MA --------

  do i = npoly+nwater+1,npoly+nwater+nmatom
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1,j2
        k = molept_list(j)

30      u = genrand_res53()
        atmvel(1,k) = dsqrt(-2.0d0*tcont_ma(1)/atmmass(k)*dlog(u))
        atmvel(2,k) = dsqrt(-2.0d0*tcont_ma(1)/atmmass(k)*dlog(u))
        u = genrand_res53()
        atmvel(1,k) = atmvel(1,k) * dcos(2.0d0*pi*u)
        atmvel(2,k) = atmvel(2,k) * dsin(2.0d0*pi*u)

        u = genrand_res53()
        atmvel(3,k) = dsqrt(-2.0d0*tcont_ma(1)/atmmass(k)*dlog(u))
        u = genrand_res53()
        atmvel(3,k) = atmvel(3,k) * dcos(2.0d0*pi*u)

        vcal = atmvel(1,k)**2 + atmvel(2,k)**2 + atmvel(3,k)**2
        vmp  = 2.0d0*tcont_ma(1)/atmmass(k)*3.5d0**2
        if (vcal > vmp) goto 30

     end do
  end do

!     +     +     +     +     +     +     +

  return
end subroutine mkmaxvel
