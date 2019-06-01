!**************************************
!*  rdheatf.f Ver.1.1 '10.06.30       *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine preheatf( npoly, npolytyp, npoly_mole, &
     &               nwater, &
     &               nmatom, &
     &               hftyp_pmole, hftyp_water, &
     &               hftyp_atm )

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: npoly           ! all number of poly
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatom          ! number of monatomic molecules

  integer,intent(in):: hftyp_pmole(:)  
                         ! atom- or mole-based heat flux cal. for poly
  integer,intent(in):: hftyp_water ! atom- or mole-based heat flux cal. for water

!     OUTPUT
  integer,intent(out):: hftyp_atm(:)    ! atom- or mole-based heat flux cal. 
                                !   for each atom

! LOCAL:
  integer:: i               ! do loop index
  integer:: j
  integer:: jj,j1,j2
  integer:: ipoly, ipolyindex

!     +     +     +     +     +     +     +

!---- intialize variable --------

  do i = 1, natom
     hftyp_atm(i) = 0
  end do

!---- assign hftyp_atom or hftyp_mole to each atom

!     polytyp
  ipolyindex = 0
  do i = 1, npolytyp
     do ipoly = 1, npoly_mole(i)
        ipolyindex = ipolyindex + 1
        j1 = molept_index(ipolyindex)
        j2 = molept_index(ipolyindex+1) - 1
        do j = j1, j2
           jj = molept_list(j)
               
           hftyp_atm(jj) = hftyp_pmole(i)
        end do
     end do
  end do

!     water
  do i = npoly+1, npoly+nwater
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1, j2
        jj = molept_list(j)

        hftyp_atm(jj) = hftyp_water
     end do
  end do

!     monatom
  do i = npoly+nwater+1, npoly+nwater+nmatom
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1
     do j = j1, j2
        jj = molept_list(j)

        hftyp_atm(jj) = HFTYP_ATOM
     end do
  end do

!     - error check
  do i = 1, natom
     if (hftyp_atm(i) == 0) then
        write(6,*) 'Error: hftyp_atm is not set in hftyp_atm(', i,').'
        stop
     end if
  end do

#if defined(_HF_DEBUG)
  write(6,*) '**************'
  write(6,*) '* Debug info *'
  write(6,*) '**************'

  write(6,*) 'HTF_ATOM=1, HTF_MOLE=2'
  write(6,*)
  write(6,*) '- polytyp'
  do i = 1, npolytyp
     write(6,*) i, hftyp_pmole(i)
  end do

  write(6,*) '- water'
  write(6,*) hftyp_water

  write(6,*)
  write(6,*) '- atom'
  do i = 1, natom
     write(6,*) i, hftyp_atm(i)
  end do
#endif

!     +     +     +     +     +     +     +

  return
end subroutine preheatf
