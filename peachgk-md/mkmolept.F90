!*************************************
!*  mkmolept.f Ver.1.6 '10.06.30     *
!*      for peachgk_md.f             *
!*            by G.Kikugawa          *
!*************************************
subroutine mkmolept( npoly, npolytyp, npoly_mole, npoly_atom, &
     &               nwater, nmatom )

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: npoly           ! all number of poly
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
  integer,intent(in):: nwater          ! number of H2O molecules
  integer,intent(in):: nmatom          ! number of monatomic molecules

! LOCAL:
  integer:: i,j             ! do loop index
  integer:: j1,j2,jj

  integer:: ipoly

!     +     +     +     +     +     +     +

!---- initialize number of pointer
  nmoleptindex = 0
  nmoleptlist = 0

!-------- pointer for poly --------

  do ipoly = 1, npolytyp

     do i = 1, npoly_mole(ipoly)

        nmoleptindex = nmoleptindex + 1
        molept_index(nmoleptindex) = nmoleptlist + 1

        do j= 1, npoly_atom(ipoly)
           nmoleptlist = nmoleptlist + 1
           molept_list(nmoleptlist) = nmoleptlist
        end do

     end do

  end do

!-------- pointer for H2O --------

  do i = 1, nwater

     nmoleptindex = nmoleptindex + 1
     molept_index(nmoleptindex) = nmoleptlist + 1

     do j = 1, 3
        nmoleptlist = nmoleptlist + 1
        molept_list(nmoleptlist) = nmoleptlist
     end do

  end do

!-------- pointer for PT --------

  do i = 1, nmatom

     nmoleptindex = nmoleptindex + 1
     molept_index(nmoleptindex) = nmoleptlist + 1

     nmoleptlist = nmoleptlist + 1
     molept_list(nmoleptlist) = nmoleptlist

  end do

!---- additional pointer

  nmoleptindex = nmoleptindex + 1
  molept_index(nmoleptindex) = nmoleptlist + 1

!---- Error check for nmoleptindex

  if (nmoleptindex /= nmole+1) then
     write(6,*) 'Failure in sum of molecules(nmoleptindex)'
     stop
  end if

!     - make irmolept_list -

  DOATOM: do i = 1, natom

     do j = 1, nmoleptindex-1
        j1 = molept_index(j)
        j2 = molept_index(j+1) - 1

        do jj = j1, j2
           if (i == molept_list(jj)) then
              irmolept_list(i) = j
              cycle DOATOM
           end if
        end do

     end do

  end do DOATOM

!     +     +     +     +     +     +     +

  return
end subroutine mkmolept
