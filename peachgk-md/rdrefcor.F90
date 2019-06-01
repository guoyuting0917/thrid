!**************************************
!*  rdrefcor.f Ver.1.1 '10.07.01      *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine rdrefcor( iuposres, &
     &               npoly, npolytyp, npoly_mole, npoly_atom, &
     &               nwater, nmatom, &
     &               xref )

  use interface_tools

  use md_global
#if defined(MPI)
  use mpi_global
#endif

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: iuposres        ! input position restraint ref. file unit

  integer,intent(in):: npoly           ! number of polymer1
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
  integer,intent(in):: nwater          ! number of H2O molecules
  integer,intent(in):: nmatom          ! number of monatomic molecules

  real(8),intent(in):: xref             ! distanse base value [m]

! LOCAL:
  integer:: imole

  character(80):: fredat(20)

  integer:: i
      
! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

#if defined(MPI)
  if (irank == 0) then
#endif
     write(6,*) 'Read reference coordinate for position restraint'
#if defined(MPI)
  end if
#endif

!---- Initialization of valiables
  imole = 0

!---- Read from state file

  DOREAD:do
     call rdfree(iuposres,20,fredat)

     if (fredat(1) == '#end') exit

     if ((fredat(1) == '#poly1') .or. (fredat(1) == '#poly')) then
        do
           call rdfree(iuposres,20,fredat)

           if (fredat(1) == ' ') cycle DOREAD
               
           imole = imole + 1
           read(fredat(2),*) ref_atmcor(1,imole)
           read(fredat(3),*) ref_atmcor(2,imole)
           read(fredat(4),*) ref_atmcor(3,imole)

!              - nondimensionalize
           ref_atmcor(1:3,imole) = ref_atmcor(1:3,imole) / xref

        end do
     end if

     if (fredat(1) == '#water') then
        imole = 0
        do i = 1, npolytyp
           imole = imole + npoly_mole(i) * npoly_atom(i)
        end do

        do
           call rdfree( iuposres, 20, fredat )

           if (fredat(1) == ' ') cycle DOREAD
               
           imole = imole + 1
           read(fredat(2),*) ref_atmcor(1,imole)
           read(fredat(3),*) ref_atmcor(2,imole)
           read(fredat(4),*) ref_atmcor(3,imole)

!              - nondimensionalize
           ref_atmcor(1:3,imole) = ref_atmcor(1:3,imole) / xref

        end do
     end if

     if (fredat(1) == '#matom') then
        imole = 0
        do i = 1, npolytyp
           imole = imole + npoly_mole(i) * npoly_atom(i)
        end do
        imole = imole + nwater*3

        do
           call rdfree( iuposres, 20, fredat )

           if (fredat(1) == ' ') cycle DOREAD
               
           imole = imole + 1
           read(fredat(2),*) ref_atmcor(1,imole)
           read(fredat(3),*) ref_atmcor(2,imole)
           read(fredat(4),*) ref_atmcor(3,imole)

!              - nondimensionalize
           ref_atmcor(1:3,imole) = ref_atmcor(1:3,imole) / xref

        end do
     end if

  end do DOREAD

!---- close file

  close(iuposres)

!     +     +     +     +     +     +     +

  return
end subroutine rdrefcor
