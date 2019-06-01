!*****************************
!*  createvel.f90 Ver.2.1    *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>

subroutine createvel(npoly, npolytyp, npoly_mole, npoly_atom, &
     &               nwater, nmatom, nmatyp, nmatomtyp, &
     &               xmaxpo, ymaxpo, zmaxpo, &
     &               xmaxw, ymaxw, zmaxw, &
     &               xmaxma, ymaxma, zmaxma, &
     &               ncrecorpo, index_crecorpo, &
     &               ncrecorw, index_crecorw, &
     &               ncrecorma, index_crecorma, &
     &               vref, timeref, &
     &               mchain, &
     &               tcont_poly, tcont_water, tcont_ma)

  use md_global
#if defined(MPI)
  use mpi_global
#endif

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: npoly           ! all number of poly
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatom          ! number of monatomic molecules
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  integer,intent(in):: xmaxpo(:)       ! use for positioning of polymer1
  integer,intent(in):: ymaxpo(:)       ! use for positioning of polymer1
  integer,intent(in):: zmaxpo(:)       ! use for positioning of polymer1
  integer,intent(in):: xmaxw           ! use for positioning of water
  integer,intent(in):: ymaxw           ! use for positioning of water
  integer,intent(in):: zmaxw           ! use for positioning of water
  integer,intent(in):: xmaxma(:)       ! use for positioning of monatomic mole.
  integer,intent(in):: ymaxma(:)       ! use for positioning of monatomic mole.
  integer,intent(in):: zmaxma(:)       ! use for positioning of monatomic mole.

  integer,intent(in):: ncrecorpo       ! max number of poly type for createcor
  integer,intent(in):: index_crecorpo(:) ! index of polymer for createcor
  integer,intent(in):: ncrecorw        ! water for createcor
  integer,intent(in):: index_crecorw   ! index of water type for createcor
  integer,intent(in):: ncrecorma       ! max number of matom type for createcor
  integer,intent(in):: index_crecorma(:) ! index of matom type for createcor

  real(8),intent(in):: vref             ! velocity base value [m/s]
  real(8),intent(in):: timeref          ! time base value [sec]

  integer,intent(in):: mchain          ! Nose-Hoover chain number

  real(8),intent(in):: tcont_poly(:)    ! poly Temp. [non-d]in NVT
  real(8),intent(in):: tcont_water     ! H2O Temp. [non-d] in NVT
  real(8),intent(in):: tcont_ma(:)      ! monatomic mole. Temp. [non-d] in NVT

! LOCAL:

  integer:: i,j,k             ! do loop index
  integer:: j1,j2

  integer:: icrepo, icrew, icrema

  integer:: ipoly
  integer:: imatom

  integer:: imole
  integer:: ix,iy,iz

!     +     +     +     +     +     +     +

!---- notation ----

#if defined(MPI)
  if (irank == 0) then
#endif
     write(6,*) 'Create atmvel from createvel'
#if defined(MPI)
  end if
#endif

!-------- Input velocity for polymer1 --------
!     poly velocity is substituted ( 0[K] control version )

  do icrepo = 1, ncrecorpo
     ipoly = index_crecorpo(icrepo)

     imole = 0
     do i = 1, ipoly-1
        imole = imole + npoly_mole(i)
     end do

     do ix = 0, xmaxpo(ipoly)-1
        do iy = 0, ymaxpo(ipoly)-1
           do iz = 0, zmaxpo(ipoly)-1

              imole = imole + 1
              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

              do j = j1, j2

                 k = molept_list(j)
                    atmvel(1:3,k) = 0.0d0

              end do

           end do
        end do
     end do

  end do

!---- non-dimensionalize
  do icrepo = 1, ncrecorpo
     ipoly = index_crecorpo(icrepo)

     imole = 0
     do i = 1, ipoly-1
        imole = imole + npoly_mole(i)
     end do

     do ix = 0, xmaxpo(ipoly)-1
        do iy = 0, ymaxpo(ipoly)-1
           do iz = 0, zmaxpo(ipoly)-1

              imole = imole + 1
              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

              do j = j1, j2

                 k = molept_list(j)
                 atmvel(1:3,k) = atmvel(1:3,k) / vref

              end do

           end do
        end do
     end do

  end do

!-------- Input velocity for H2O --------
!     water molecule velocity is substituted ( 0[K] control version )
  do icrew = 1, ncrecorw

     imole = npoly

     do ix = 0, xmaxw-1
        do iy = 0, ymaxw-1
           do iz = 0, zmaxw-1

              imole = imole + 1
              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

              do j = j1, j2

                 k = molept_list(j)
                 atmvel(1:3,k) = 0.0d0

              end do

           end do
        end do
     end do

  end do

!---- non-dimensionalize
  do icrew = 1, ncrecorw

     imole = npoly

     do ix = 0, xmaxw-1
        do iy = 0, ymaxw-1
           do iz = 0, zmaxw-1

              imole = imole + 1
              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

              do j = j1, j2

                 k = molept_list(j)
                 atmvel(1:3,k) = atmvel(1:3,k) / vref

              end do

           end do
        end do
     end do

  end do

!-------- Input velocity for MA --------
!     monatomic velocity is substituted ( 0[K] control version )

  do icrema = 1, ncrecorma
     imatom = index_crecorma(icrema)

     imole = npoly + nwater
     do i = 1, imatom-1
        imole = imole + nmatomtyp(i)
     end do

     do ix = 0, xmaxma(imatom)-1
        do iy = 0, ymaxma(imatom)-1
           do iz = 0, zmaxma(imatom)-1

              imole = imole + 1
              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

              do j = j1, j2

                 k = molept_list(j)
                 atmvel(1:3,k) = 0.0d0

              end do

           end do
        end do
     end do

  end do

!---- non-dimensionalize
  do icrema = 1, ncrecorma
     imatom = index_crecorma(icrema)

     imole = npoly + nwater
     do i = 1, imatom-1
        imole = imole + nmatomtyp(i)
     end do

     do ix = 0 ,xmaxma(imatom)-1
        do iy = 0, ymaxma(imatom)-1
           do iz = 0, zmaxma(imatom)-1

              imole = imole + 1
              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

              do j = j1, j2

                 k = molept_list(j)
                 atmvel(1:3,k) = atmvel(1:3,k) / vref

              end do

           end do
        end do
     end do

  end do

!-------- Input velocity for Nose-Hoover chain --------
!     vlogs[1/s] is set to 0

  vlogs(1:mchain) = 0.0d0

!-------- Input velocity for Andersen (Hoover type) barostat --------
!     vlogv[1/s] is set to 0

  vlogv = 0.0d0
  vboxg(1:3) = 0.0d0

!---- non-dimensionalize

  vlogs(1:mchain) = vlogs(1:mchain) * timeref

  vlogv = vlogv * timeref
  vboxg(1:3) = vboxg(1:3) * timeref

!     +     +     +     +     +     +     +

end subroutine createvel
