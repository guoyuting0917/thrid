!*****************************
!*  rdstarec.f90 Ver.2.2     *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>

subroutine rdstarec(iostarec, &
     &              npoly, npolytyp, npoly_mole, npoly_atom, &
     &              nwater, nmatom, &
     &              xcel, ycel, zcel, yratio, zratio, &
     &              xref, vref, timeref, pref, &
     &              mchain, &
     &              pint, pintt, &
     &              ifsetcor)

  use interface_tools

  use md_global
#if defined(MPI)
  use mpi_global
#endif

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: iostarec        ! state record file unit

  integer,intent(in):: npoly           ! number of polymer1
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
  integer,intent(in):: nwater          ! number of H2O molecules
  integer,intent(in):: nmatom          ! number of monatomic molecules

  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: vref             ! velocity base value [m/s]
  real(8),intent(in):: timeref          ! time base value [sec]
  real(8),intent(in):: pref             ! pressure base value [Pa]

  integer,intent(in):: mchain          ! Nose-Hoover chain number

!     OUTPUT
  real(8),intent(out):: xcel             ! x cell length
  real(8),intent(out):: ycel             ! y cell length
  real(8),intent(out):: zcel             ! z cell length
  real(8),intent(out):: yratio           ! y cell ratio of y to x
  real(8),intent(out):: zratio           ! z cell ratio of z to x

  real(8),intent(out):: pint             ! internal pressure
  real(8),intent(out):: pintt(:,:)       ! internal pressure tensor

  logical,intent(out):: ifsetcor(:)     ! frag for checking if coordinate has been set

! LOCAL:
!      integer:: npoly_s         ! number of polymer1
!      integer:: nwater_s        ! number of H2O molecules
!      integer:: nmatom_s        ! number of monatomic molecules

!      integer:: natom_s         ! number of all atoms

  integer:: mchain_s        ! Nose-Hoover chain number

!      real*8:: xcel_s           ! x cell length
!      real*8:: ycel_s           ! y cell length
!      real*8:: zcel_s           ! z cell length

  integer:: nstep_s         ! maxnstep of last MD

  character(2):: atmtyp_s ! atom type

  integer:: imole
  integer:: inhc

  character(80):: fredat(20)

  integer:: i

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

#if defined(MPI)
  if (irank == 0) then
#endif
     write(6,*) 'Read atmcor&vel from state_file'
#if defined(MPI)
  end if
#endif

!---- Initialization of valiables
  imole = 0
  inhc = 0

!     clear the ifsetcor
  do i=1,nmole
     ifsetcor(i) = .false.
  end do

!---- Read from state file

  DOREAD:do
     call rdfree( iostarec, 20, fredat )

     if (fredat(1) == '#end') exit

!         if (fredat(1) .eq. '#npoly') then
!            read(fredat(2),*) npoly_s
!            if (npoly_s .ne. npoly) then
!               write(6,*) 'Failure: npoly in rdstarec'
!               stop
!            end if
!         end if
!
!         if (fredat(1) .eq. '#nwater') then
!            read(fredat(2),*) nwater_s
!            if (nwater_s .ne. nwater) then
!               write(6,*) 'Failure: nwater in rdstarec'
!               stop
!            end if
!         end if
!
!         if (fredat(1) .eq. '#nmatom') then
!            read(fredat(2),*) nmatom_s
!            if (nmatom_s .ne. nmatom) then
!               write(6,*) 'Failure: nmatom in rdstarec'
!               stop
!            end if
!         end if
!
!         if (fredat(1) .eq. '#natom') then
!            read(fredat(2),*) natom_s
!            if (natom_s .ne. natom) then
!               write(6,*) 'Failure: natom in rdstarec'
!               stop
!            end if
!         end if

     if (fredat(1) == '#mchain') then
        read(fredat(2),*) mchain_s
        if (mchain_s /= mchain) then
           write(6,*) 'Failure: mchain in rdstarec'
           stop
        end if
     end if

     if (fredat(1) == '#xcel') then
        read(fredat(2),*) xcel
#if defined(MPI)
        if (irank == 0) then
#endif
           write(6,'(A9,E27.18)') 'old xcel=',xcel
#if defined(MPI)
        end if
#endif
        xcel = xcel / xref
     end if

     if (fredat(1) == '#ycel') then
        read(fredat(2),*) ycel
#if defined(MPI)
        if (irank == 0) then
#endif
           write(6,'(A9,E27.18)') 'old ycel=',ycel
#if defined(MPI)
        end if
#endif
        ycel = ycel / xref
     end if

     if (fredat(1) .eq. '#zcel') then
        read(fredat(2),*) zcel
#if defined(MPI)
        if (irank == 0) then
#endif
           write(6,'(A9,E27.18)') 'old zcel=',zcel
#if defined(MPI)
        end if
#endif
        zcel = zcel / xref
     end if

     if (fredat(1) .eq. '#nstep') then
        read(fredat(2),*) nstep_s
#if defined(MPI)
        if (irank == 0) then
#endif
           write(6,*) 'Last MD maxnstep: ',nstep_s
#if defined(MPI)
        end if
#endif
     end if

     if ((fredat(1) == '#poly1') .or. &
          & (fredat(1) == '#poly')) then
        do
           call rdfree(iostarec,20,fredat)

           if (fredat(1) == ' ') cycle DOREAD

           imole = imole + 1
           read(fredat(2),*) atmcor(1,imole)
           read(fredat(3),*) atmcor(2,imole)
           read(fredat(4),*) atmcor(3,imole)
           read(fredat(5),*) atmvel(1,imole)
           read(fredat(6),*) atmvel(2,imole)
           read(fredat(7),*) atmvel(3,imole)
           atmtyp(imole) = fredat(8)(1:2)
           ifsetcor(irmolept_list(imole)) = .true.
#if defined(_DO_NOT_CHECK_THIS)
           !!! This error check is abbreviated from Ver.2.135
           atmtyp_s = fredat(8)(1:2)
           if (atmtyp_s /= atmtyp(imole)) then
              write(6,*) 'Failure: atmtyp in rdstarec, atom No.', imole
              stop
           end if
#endif

!              - nondimensionalize
           atmcor(1,imole) = atmcor(1,imole) / xref
           atmcor(2,imole) = atmcor(2,imole) / xref
           atmcor(3,imole) = atmcor(3,imole) / xref
           atmvel(1,imole) = atmvel(1,imole) / vref
           atmvel(2,imole) = atmvel(2,imole) / vref
           atmvel(3,imole) = atmvel(3,imole) / vref

        end do
     end if

     if (fredat(1) == '#water') then
        imole = 0
        do i = 1, npolytyp
           imole = imole + npoly_mole(i) * npoly_atom(i)
        end do

        do
           call rdfree( iostarec, 20, fredat )

           if (fredat(1) == ' ') cycle DOREAD

           imole = imole + 1
           read(fredat(2),*) atmcor(1,imole)
           read(fredat(3),*) atmcor(2,imole)
           read(fredat(4),*) atmcor(3,imole)
           read(fredat(5),*) atmvel(1,imole)
           read(fredat(6),*) atmvel(2,imole)
           read(fredat(7),*) atmvel(3,imole)
           atmtyp(imole) = fredat(8)(1:2)
           ifsetcor(irmolept_list(imole)) = .true.
#if defined(_DO_NOT_CHECK_THIS)
           !!! This error check is abbreviated from Ver.2.135
           atmtyp_s = fredat(8)(1:2)
           if (atmtyp_s .ne. atmtyp(imole)) then
              write(6,*) 'Failure: atmtyp in rdstarec, atom No.', &
                   &     imole,atmtyp(imole),atmtyp_s
              stop
           end if
#endif

!              - nondimensionalize
           atmcor(1,imole) = atmcor(1,imole) / xref
           atmcor(2,imole) = atmcor(2,imole) / xref
           atmcor(3,imole) = atmcor(3,imole) / xref
           atmvel(1,imole) = atmvel(1,imole) / vref
           atmvel(2,imole) = atmvel(2,imole) / vref
           atmvel(3,imole) = atmvel(3,imole) / vref

        end do
     end if

     if (fredat(1) .eq. '#matom') then
        imole = 0
        do i = 1, npolytyp
           imole = imole + npoly_mole(i) * npoly_atom(i)
        end do
        imole = imole + nwater*3

        do
           call rdfree( iostarec, 20, fredat)

           if (fredat(1) == ' ') cycle DOREAD

           imole = imole + 1
           read(fredat(2),*) atmcor(1,imole)
           read(fredat(3),*) atmcor(2,imole)
           read(fredat(4),*) atmcor(3,imole)
           read(fredat(5),*) atmvel(1,imole)
           read(fredat(6),*) atmvel(2,imole)
           read(fredat(7),*) atmvel(3,imole)
           atmtyp(imole) = fredat(8)(1:2)
           ifsetcor(irmolept_list(imole)) = .true.
#if defined(_DO_NOT_CHECK_THIS)
           !!! This error check is abbreviated from Ver.2.135
           atmtyp_s = fredat(8)(1:2)
           if (atmtyp_s .ne. atmtyp(imole)) then
              write(6,*) 'Failure: atmtyp in rdstarec, atom No.', &
                   &     imole,atmtyp(imole),atmtyp_s
              stop
           end if
#endif

!              - nondimensionalize
           atmcor(1,imole) = atmcor(1,imole) / xref
           atmcor(2,imole) = atmcor(2,imole) / xref
           atmcor(3,imole) = atmcor(3,imole) / xref
           atmvel(1,imole) = atmvel(1,imole) / vref
           atmvel(2,imole) = atmvel(2,imole) / vref
           atmvel(3,imole) = atmvel(3,imole) / vref

        end do
     end if

     if (fredat(1) == '#nhc') then
        do
           call rdfree( iostarec, 20, fredat )

           if (fredat(1) == ' ') cycle DOREAD

           inhc = inhc + 1
           read(fredat(2),*) xlogs(inhc)
           read(fredat(3),*) vlogs(inhc)

!              - nondimensionalize
           vlogs(inhc) = vlogs(inhc) * timeref

        end do
     end if

     if (fredat(1) == '#andersen') then
        do
           call rdfree( iostarec, 20, fredat )

           if (fredat(1) == ' ') cycle DOREAD

           read(fredat(2),*) vlogv
           read(fredat(3),*) pint

!              - nondimensionalize
           vlogv = vlogv * timeref
           pint = pint / pref

        end do
     end if

     if (fredat(1) == '#andersen_aniso') then
        do
           call rdfree( iostarec, 20, fredat )

           if (fredat(1) == ' ') cycle DOREAD

           read(fredat(2),*) vboxg(1)
           read(fredat(3),*) vboxg(2)
           read(fredat(4),*) vboxg(3)
           read(fredat(5),*)  pintt(1,1)
           read(fredat(6),*)  pintt(2,1)
           read(fredat(7),*)  pintt(3,1)
           read(fredat(8),*)  pintt(1,2)
           read(fredat(9),*)  pintt(2,2)
           read(fredat(10),*) pintt(3,2)
           read(fredat(11),*) pintt(1,3)
           read(fredat(12),*) pintt(2,3)
           read(fredat(13),*) pintt(3,3)

!              - nondimensionalize
           vboxg(1:3) = vboxg(1:3) * timeref
           pintt(1:3,1:3) = pintt(1:3,1:3) / pref

        end do
     end if

  end do DOREAD

!      if (imole .ne. natom) then
!         write(6,*) 'Failure: atm table is wrong in state file'
!         stop
!      end if

  if (inhc /= mchain) then
     write(6,*) 'Failure: nhc table is wrong in state file'
     stop
  end if

! - calculate yratio and xratio
  yratio = ycel / xcel
  zratio = zcel / xcel

!     +     +     +     +     +     +     +

  return
end subroutine rdstarec
