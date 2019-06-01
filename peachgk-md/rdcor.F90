!*****************************
!*  rdcor.f Ver.2.1          *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2014-06-13 16:34:56 gota>

subroutine rdcor(iucor,xref, &
     &           npoly, npolytyp, npoly_mole, npoly_atom, &
     &           nwater, nmatom, &
     &           nmatyp, nmatomtyp, &
     &           polytyp_free, &
     &           oatmtyp, hatmtyp, monoatmtyp, &
     &           ifsetchrg, atmchrg_tmp)

  use interface_tools

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: iucor(:)        ! input poly coordinate file unit
  real(8),intent(in):: xref             ! distanse base value [m]

  integer,intent(in):: npoly           ! all number of poly
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatom          ! number of monatomic molecules
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control

  character(2),intent(in):: oatmtyp ! O atomtype of water model
  character(2),intent(in):: hatmtyp ! H atomtype of water model
  character(2),intent(in):: monoatmtyp(:) ! monatomic mole. type

  logical,intent(in):: ifsetchrg(:)    ! if set charge from cor file

!     OUTPUT
  real(8),intent(out):: atmchrg_tmp(:) ! temporary atom charge read from cor file

! LOCAL:
  character(80):: fredat(20)
  integer:: i,m,n           ! do loop index

  real(8):: h1cor(3)      ! coordinate of atoms
  real(8):: h2cor(3)      ! coordinate of atoms
  real(8):: oxcor(3)      ! coordinate of atoms

  integer:: polyindex
  integer:: maindex

  integer:: ipoly

  integer:: npoly_atom_tmp  ! for checking the npoly_atom(*)

!      integer:: iword

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- some preparation

  do i = 1, maxnatom
     atmchrg_tmp(i) = 0.0d0
  end do

!---- read from coordinate

  natom = 0

  DO ipoly=1,npolytyp

     polyindex = natom
     npoly_atom_tmp = 0

     DOREAD:do
        call rdfree(iucor(ipoly),20,fredat)

        IF (fredat(1) == 'END') THEN
           exit
        END IF

        if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                         ! comment line

        if (fredat(1) /= ' ') then

           natom = natom + 1
           if (natom > maxnatom) then
              write(6,*) 'Error : natom exceeds maxnatom= ',maxnatom
              stop
           end if

           npoly_atom_tmp = npoly_atom_tmp + 1
           atmtyp(natom) = fredat(2)(1:2)
           read(fredat(3),*) atmcor(1,natom)
           read(fredat(4),*) atmcor(2,natom)
           read(fredat(5),*) atmcor(3,natom)

!---- fixed '07.01.22
!               do iword = 1, maxnword
!
!                  if (polytyp_free(ipoly,iword) .eq. 'setcharge') then
!                     read(fredat(6),*) atmchrg_tmp(natom)
!                     exit
!                  end if
!
!               end do
           if (ifsetchrg(ipoly)) then
              read(fredat(6),*) atmchrg_tmp(natom)
           end if
!---- fixed '07.01.22

        end if

     END DO DOREAD

!----    check the npoly_atom

     if (npoly_atom_tmp /= npoly_atom(ipoly)) then
        write(6,*) 'Error in rdcor:', &
             &     ' descripancy npoly_atom and that in cor file'

        write(6,*) '                npolymoletyp No.',ipoly
        stop
     end if

!----    non-dimensionalize coordinate
     do i=1,npoly_atom(ipoly)
        atmcor(1,polyindex+i) = atmcor(1,polyindex+i) / xref
        atmcor(2,polyindex+i) = atmcor(2,polyindex+i) / xref
        atmcor(3,polyindex+i) = atmcor(3,polyindex+i) / xref
     end do

!-------- duplicate atmcor to all polymer1
     do n=2,npoly_mole(ipoly)
        do m=1,npoly_atom(ipoly)
           natom = natom + 1
           if (natom > maxnatom) then
              write(6,*) 'Error : natom exceeds maxnatom= ',maxnatom
              stop
           end if

           atmcor(1,natom) = atmcor(1,polyindex+m)
           atmcor(2,natom) = atmcor(2,polyindex+m)
           atmcor(3,natom) = atmcor(3,polyindex+m)
           atmtyp(natom) = atmtyp(polyindex+m)
           atmchrg_tmp(natom) = atmchrg_tmp(polyindex+m)
        end do
     end do

     close(iucor(ipoly))

  END DO


!-------- create H2O coordinate --------

  h1cor(1) =  0.0d0             / xref
  h1cor(2) = -0.8166415552d-10  / xref
  h1cor(3) = -0.5125623357d-10  / xref
  h2cor(1) =  0.0d0             / xref
  h2cor(2) =  0.8166415552d-10  / xref
  h2cor(3) = -0.5125623357d-10  / xref
  oxcor(1) =  0.0d0             / xref
  oxcor(2) =  0.0d0             / xref
  oxcor(3) =  0.0645828543d-10  / xref

  do n=1,nwater
     natom = natom + 1
     if (natom > maxnatom) then
        write(6,*) 'Error : natom exceeds maxnatom= ',maxnatom
        stop
     end if

     atmcor(1,natom) = h1cor(1)
     atmcor(2,natom) = h1cor(2)
     atmcor(3,natom) = h1cor(3)
     atmtyp(natom) = hatmtyp

     natom = natom + 1
     if (natom > maxnatom) then
        write(6,*) 'Error : natom exceeds maxnatom= ',maxnatom
        stop
     end if

     atmcor(1,natom) = oxcor(1)
     atmcor(2,natom) = oxcor(2)
     atmcor(3,natom) = oxcor(3)
     atmtyp(natom) = oatmtyp

     natom = natom + 1
     if (natom > maxnatom) then
        write(6,*) 'Error : natom exceeds maxnatom= ',maxnatom
        stop
     end if

     atmcor(1,natom) = h2cor(1)
     atmcor(2,natom) = h2cor(2)
     atmcor(3,natom) = h2cor(3)
     atmtyp(natom) = hatmtyp
  end do

!-------- create monatomic coordinate --------

  maindex = natom
  do n=1,nmatom
     natom = natom + 1
     if (natom > maxnatom) then
        write(6,*) 'Error : natom exceeds maxnatom= ',maxnatom
        stop
     end if

     atmcor(1,natom) = 0.0d0
     atmcor(2,natom) = 0.0d0
     atmcor(3,natom) = 0.0d0
  end do

  do n=1,nmatyp
     do m=1,nmatomtyp(n)
        maindex = maindex + 1
        atmtyp(maindex) = monoatmtyp(n)
     end do
  end do

!---- print natom & nmole

  nmole = npoly + nwater + nmatom

!      write(6,'(3X,A6,I8)') 'natom: ',natom
!      write(6,'(3X,A6,I8)') 'nmole: ',nmole

!     +     +     +     +     +     +     +

  return
end subroutine rdcor
