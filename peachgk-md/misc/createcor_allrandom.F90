!*************************************
!*  createcor.f Ver.2.7 '10.06.30    *
!*      for peachgk_md.f             *
!*            by G.Kikugawa          *
!*************************************
subroutine createcor(npoly,npolytyp,npoly_mole,npoly_atom,   &
     &               nwater,nmatom,nmatyp,nmatomtyp,   &
     &               xmaxpo,ymaxpo,zmaxpo,   &
     &               xmaxw,ymaxw,zmaxw,   &
     &               xmaxma,ymaxma,zmaxma,   &
     &               inicorpo,inicorw,inicorma,   &
     &               ncrecorpo,index_crecorpo,   &
     &               ncrecorw,index_crecorw,   &
     &               ncrecorma,index_crecorma,   &
     &               ifsetcor,   &
     &               xcel,ycel,zcel,   &
     &               xref,   &
     &               compfact,   &
     &               mchain)

!#################################################################
!# Usage of the inicorpo(), inicorw(), inicorma()                #
!#   inicorpo(), inicorw(), inicorma() are the real*8 type.      #
!#   Therefore, you have to cast the type which you want to use. #
!#   For example, int(inicorpo()) to cast the int type           #
!#################################################################

  use md_global
#if defined(MPI)
  use mpi_global
#endif

  implicit none

! ARGUMENT:
!   INPUT
  integer,intent(in):: npoly            ! all number of poly
  integer,intent(in):: npolytyp         ! number of poly type
  integer,intent(in):: npoly_mole(:)    ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)    ! number of atoms belonging to poly

  integer,intent(in):: nwater           ! number of H2O molecules

  integer,intent(in):: nmatom           ! number of monatomic molecules
  integer,intent(in):: nmatyp           ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)     ! each number of monatomic mole.

  integer,intent(in):: xmaxpo(:)        ! use for positioning of polymer1
  integer,intent(in):: ymaxpo(:)        ! use for positioning of polymer1
  integer,intent(in):: zmaxpo(:)        ! use for positioning of polymer1
  integer,intent(in):: xmaxw            ! use for positioning of water
  integer,intent(in):: ymaxw            ! use for positioning of water
  integer,intent(in):: zmaxw            ! use for positioning of water
  integer,intent(in):: xmaxma(:)        ! use for positioning of monatomic mole.
  integer,intent(in):: ymaxma(:)        ! use for positioning of monatomic mole.
  integer,intent(in):: zmaxma(:)        ! use for positioning of monatomic mole.
  real(8),intent(in):: inicorpo(:,:) ! use for positioning of poly
  real(8),intent(in):: inicorw(:)       ! use for positioning of water
  real(8),intent(in):: inicorma(:,:) ! use for positioning of ma
  integer,intent(in):: ncrecorpo        ! max number of poly type for createcor
  integer,intent(in):: index_crecorpo(:) ! index of polymer for createcor
  integer,intent(in):: ncrecorw         ! water for createcor
  integer,intent(in):: index_crecorw    ! index of water type for createcor
  integer,intent(in):: ncrecorma        ! max number of matom type for createcor
  integer,intent(in):: index_crecorma(:) ! index of matom type for createcor

  logical,intent(out):: ifsetcor(:)     ! frag for checking if coordinate has set

  real(8),intent(in):: xref             ! distanse base value [m]

  real(8),intent(in):: xcel             ! x cell length
  real(8),intent(in):: ycel             ! y cell length
  real(8),intent(in):: zcel             ! z cell length

!!! old MT random generator
!      real*8:: randseed         ! random seed for createcor etc.
!      integer:: nperiod         ! period parameter
!      parameter(nperiod = 624)
!      integer:: mt(0:nperiod-1)
!      integer:: mti

  real(8),intent(in):: compfact     ! compact factor using at poly arrange(<1.0)

  integer,intent(in):: mchain           ! Nose-Hoover chain number

! LOCAL:
  integer:: i
  integer:: j,jj                        ! do loop index
  integer:: m,n                         ! do loop index
  integer:: n1,n2,nn
  integer:: j1,j2,k,kk
  integer:: ix,iy,iz
  integer:: imole
  integer:: icount

  real(8):: rot_m(3,3)                  ! rotational matrix
  real(8):: psi,phi,theta               ! Euler angle
  real(8):: rxp,ryp,rzp                 ! coordinate of poly1
  real(8):: rxw,ryw,rzw                 ! coordinate of water
  real(8):: rxma,ryma,rzma              ! coordinate of monatomic mole.

  real(8):: pi                          ! = 3.14159...

  real(8):: atmcor_t(3)                 ! atmcor tmp

  real(8):: r_ij(3)                     ! xi - xj etc.
  real(8):: rij                         ! distance of r_ij

  integer:: ipoly
  integer:: imatom

  integer:: icrepo, icrewa, icrema

! FUNCTIONS:
!  real(8):: grnd                        ! function creating random number
  real(8):: genrand_res53                ! SFMT random generator
      
!     +     +     +     +     +     +     +

!---- notation ----
#if defined(MPI)
  if (irank == 0) then
#endif

     write(6,*) 'Create atmcor from createcor'

#if defined(MPI)
  end if
#endif

  pi = acos(-1.0d0)

!-------- Locate molecules of poly --------

  do icrepo = 1, ncrecorpo
     ipoly = index_crecorpo(icrepo)

     imole = 0
     do i = 1, ipoly-1
        imole = imole + npoly_mole(i)
     end do

     do ix=0,xmaxpo(ipoly)-1
        do iy=0,ymaxpo(ipoly)-1
           do iz=0,zmaxpo(ipoly)-1

              icount = 0
              imole = imole + 1
              ifsetcor(imole) = .true.

              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

!             --- make gravity coordinate ---
              rxp = compfact*(genrand_res53()-0.5d0)*xcel + 0.5d0*xcel
              ryp = compfact*(genrand_res53()-0.5d0)*ycel + 0.5d0*ycel
              rzp = compfact*(genrand_res53()-0.5d0)*zcel + 0.5d0*zcel

!             --- make rotational matrix ---
              psi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
              phi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
              theta = (genrand_res53() - 0.5d0) * 2.0d0 * pi
!              psi   = 0.5d0 * pi
!              phi   = 0.0d0
!              theta = 0.5d0 * pi

              rot_m(1,1) =   cos(psi)*cos(phi)   &
                   &       - cos(theta)*sin(phi)*sin(psi)
              rot_m(1,2) =   cos(psi)*sin(phi)   &
                   &       + cos(theta)*cos(phi)*sin(psi)
              rot_m(1,3) =   sin(psi)*sin(theta)
              rot_m(2,1) = - sin(psi)*cos(phi)   &
                   &       - cos(theta)*sin(phi)*cos(psi)
              rot_m(2,2) = - sin(psi)*sin(phi)   &
                   &       + cos(theta)*cos(phi)*cos(psi)
              rot_m(2,3) =   cos(psi)*sin(theta)
              rot_m(3,1) =   sin(theta)*sin(phi)
              rot_m(3,2) = - sin(theta)*cos(phi)
              rot_m(3,3) =   cos(theta)

              !--- make coordinate ---
              do j=j1,j2

                 k = molept_list(j)

                 atmcor_t(1:3) = atmcor(1:3,k)

                 atmcor(1,k) = rot_m(1,1)*atmcor_t(1)   &
                      &      + rot_m(2,1)*atmcor_t(2)   &
                      &      + rot_m(3,1)*atmcor_t(3)   &
                      &      + rxp
                 atmcor(2,k) = rot_m(1,2)*atmcor_t(1)   &
                      &      + rot_m(2,2)*atmcor_t(2)   &
                      &      + rot_m(3,2)*atmcor_t(3)   &
                      &      + ryp
                 atmcor(3,k) = rot_m(1,3)*atmcor_t(1)   &
                      &      + rot_m(2,3)*atmcor_t(2)   &
                      &      + rot_m(3,3)*atmcor_t(3)   &
                      &      + rzp

              end do
#if defined(MPI)
              if (irank == 0) then
#endif
                 write(6,*) 'Set atmcor of mole No. ',imole,',',icount
#if defined(MPI)
              end if
#endif
           end do
        end do
     end do

  end do

!-------- Locate molecules of H2O --------

  do icrewa = 1, ncrecorw

     imole = npoly

     !---- arrange water molecules
     do ix=0,xmaxw-1
        do iy=0,ymaxw-1
           do iz=0,zmaxw-1

              icount = 0
              imole = imole + 1
              ifsetcor(imole) = .true.

              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

              !--- make rotational matrix ---
              psi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
              phi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
              theta = (genrand_res53() - 0.5d0) * 2.0d0 * pi

              rxw = genrand_res53()*xcel
              ryw = genrand_res53()*ycel
              rzw = genrand_res53()*zcel

              rot_m(1,1) =   cos(psi)*cos(phi)   &
                   &       - cos(theta)*sin(phi)*sin(psi)
              rot_m(1,2) =   cos(psi)*sin(phi)   &
              &            + cos(theta)*cos(phi)*sin(psi)
              rot_m(1,3) =   sin(psi)*sin(theta)
              rot_m(2,1) = - sin(psi)*cos(phi)   &
                   &       - cos(theta)*sin(phi)*cos(psi)
              rot_m(2,2) = - sin(psi)*sin(phi)   &
                   &       + cos(theta)*cos(phi)*cos(psi)
              rot_m(2,3) =   cos(psi)*sin(theta)
              rot_m(3,1) =   sin(theta)*sin(phi)
              rot_m(3,2) = - sin(theta)*cos(phi)
              rot_m(3,3) =   cos(theta)

              do j= j1, j2
                 k = molept_list(j)

                 atmcor_t(1:3) = atmcor(1:3,k)

                 atmcor(1,k) = rot_m(1,1)*atmcor_t(1)   &
                      &      + rot_m(2,1)*atmcor_t(2)   &
                      &      + rot_m(3,1)*atmcor_t(3)   &
                      &      + rxw
                 atmcor(2,k) = rot_m(1,2)*atmcor_t(1)   &
                      &      + rot_m(2,2)*atmcor_t(2)   &
                      &      + rot_m(3,2)*atmcor_t(3)   &
                      &      + ryw
                 atmcor(3,k) = rot_m(1,3)*atmcor_t(1)   &
                      &      + rot_m(2,3)*atmcor_t(2)   &
                      &      + rot_m(3,3)*atmcor_t(3)   &
                      &      + rzw

              end do

#if defined(MPI)
              if (irank == 0) then
#endif
                 write(6,*) 'Set atmcor of mole No. ',imole,',',icount
#if defined(MPI)
              end if
#endif
           end do
        end do
     end do

  end do

!-------- Locate molecules of monatomic mole. --------

  do icrema = 1, ncrecorma
     imatom = index_crecorma(icrema)

     imole = npoly + nwater
     do i = 1, imatom-1
        imole = imole + nmatomtyp(i)
     end do

     !---- arrange MA molecules
     do ix=0,xmaxma(imatom)-1
        do iy=0,ymaxma(imatom)-1
           do iz=0,zmaxma(imatom)-1

              icount = 0
              imole = imole + 1
              ifsetcor(imole) = .true.

              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

              rxma = genrand_res53()*xcel
              ryma = genrand_res53()*ycel
              rzma = genrand_res53()*zcel

              do j= j1, j2
                 k = molept_list(j)

                 atmcor_t(1:3) = atmcor(1:3,k)

                 atmcor(1,k) = rxma
                 atmcor(2,k) = ryma
                 atmcor(3,k) = rzma

              end do

#if defined(MPI)
              if (irank == 0) then
#endif
                 write(6,*) 'Set atmcor of mole No. ',imole,',',icount
#if defined(MPI)
              end if
#endif
           end do
        end do
     end do

  end do

!----- Preset the Nose-Hoover chain thermostat -----

#if defined(MPI)
  if (irank == 0) then
#endif
     write(6,*) 'Set Nose-Hoover chain thermostat'
     write(6,*) 'chain number= ',mchain
#if defined(MPI)
  end if
#endif
  xlogs(1:mchain) = 0.0d0

!     +     +     +     +     +     +     +

end subroutine createcor
