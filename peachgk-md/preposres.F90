!**************************************
!*  preposres.f Ver.1.3 '13.12.10     *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine preposres(npolytyp,npoly_mole,npoly_atom, &
     &               nwater, &
     &               nmatyp,nmatomtyp, &
     &               polytyp_free,watertyp_free,matomtyp_free)

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
  character(80),intent(in):: watertyp_free(:) ! use for water type control
  character(80),intent(in):: matomtyp_free(:,:)
                                ! use for matom type control

! LOCAL:
  integer:: ipoly
  integer:: imatom
  integer:: iword
      
  integer:: index_atom
  integer:: index_atom_tmp

  integer:: i,j

  integer:: nresatm
  integer:: index_resatm(maxnmol)

!     +     +     +     +     +     +

!---- initialization

  nposres = 0

!---- make list for restraint atoms of poly type

  do ipoly = 1, npolytyp
     do iword = 1, maxnword

!       constraint for a whole molecule
        if (polytyp_free(ipoly,iword) == 'posres') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, npoly_atom(ipoly)
                 nposres(1:3) = nposres(1:3) + 1
                 index_atom = index_atom + 1
                 index_posres(nposres(1),1) = index_atom
                 index_posres(nposres(2),2) = index_atom
                 index_posres(nposres(3),3) = index_atom
              end do
           end do

           exit

        else if (polytyp_free(ipoly,iword) == 'posresx') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, npoly_atom(ipoly)
                 nposres(1) = nposres(1) + 1
                 index_atom = index_atom + 1
                 index_posres(nposres(1),1) = index_atom
              end do
           end do

        else if (polytyp_free(ipoly,iword) == 'posresy') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, npoly_atom(ipoly)
                 nposres(2) = nposres(2) + 1
                 index_atom = index_atom + 1
                 index_posres(nposres(2),2) = index_atom
              end do
           end do

        else if (polytyp_free(ipoly,iword) == 'posresz') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, npoly_atom(ipoly)
                 nposres(3) = nposres(3) + 1
                 index_atom = index_atom + 1
                 index_posres(nposres(3),3) = index_atom
              end do
           end do

!       constraint for specific atoms
        else if (polytyp_free(ipoly,iword) == 'posresatm') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           read(polytyp_free(ipoly,iword+1),*) nresatm
           do i = 1, nresatm
              read(polytyp_free(ipoly,iword+1+i),*) index_resatm(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, nresatm
                 nposres(1:3) = nposres(1:3) + 1
                 index_atom_tmp = index_atom + index_resatm(j)
                 index_posres(nposres(1),1) = index_atom_tmp
                 index_posres(nposres(2),2) = index_atom_tmp
                 index_posres(nposres(3),3) = index_atom_tmp
              end do
              index_atom = index_atom + npoly_atom(ipoly)
           end do

           exit

        else if (polytyp_free(ipoly,iword) == 'posresatmx') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           read(polytyp_free(ipoly,iword+1),*) nresatm
           do i = 1, nresatm
              read(polytyp_free(ipoly,iword+1+i),*) index_resatm(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, nresatm
                 nposres(1) = nposres(1) + 1
                 index_atom_tmp = index_atom + index_resatm(j)
                 index_posres(nposres(1),1) = index_atom_tmp
              end do
              index_atom = index_atom + npoly_atom(ipoly)
           end do

        else if (polytyp_free(ipoly,iword) == 'posresatmy') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           read(polytyp_free(ipoly,iword+1),*) nresatm
           do i = 1, nresatm
              read(polytyp_free(ipoly,iword+1+i),*) index_resatm(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, nresatm
                 nposres(2) = nposres(2) + 1
                 index_atom_tmp = index_atom + index_resatm(j)
                 index_posres(nposres(2),2) = index_atom_tmp
              end do
              index_atom = index_atom + npoly_atom(ipoly)
           end do

        else if (polytyp_free(ipoly,iword) == 'posresatmz') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           read(polytyp_free(ipoly,iword+1),*) nresatm
           do i = 1, nresatm
              read(polytyp_free(ipoly,iword+1+i),*) index_resatm(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, nresatm
                 nposres(3) = nposres(3) + 1
                 index_atom_tmp = index_atom + index_resatm(j)
                 index_posres(nposres(3),3) = index_atom_tmp
              end do
              index_atom = index_atom + npoly_atom(ipoly)
           end do

        end if

     end do
  end do

!---- make list for restraint atoms of water type

  do iword = 1, maxnword

!    constraint for a whole molecule
     if (watertyp_free(iword) == 'posres') then
        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        do i = 1, nwater
           do j = 1, 3
              nposres(1:3) = nposres(1:3) + 1
              index_atom = index_atom + 1
              index_posres(nposres(1),1) = index_atom
              index_posres(nposres(2),2) = index_atom
              index_posres(nposres(3),3) = index_atom
           end do
        end do

        exit

     else if (watertyp_free(iword) == 'posresx') then
        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        do i = 1, nwater
           do j = 1, 3
              nposres(1) = nposres(1) + 1
              index_atom = index_atom + 1
              index_posres(nposres(1),1) = index_atom
           end do
        end do

     else if (watertyp_free(iword) == 'posresy') then
        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        do i = 1, nwater
           do j = 1, 3
              nposres(2) = nposres(2) + 1
              index_atom = index_atom + 1
              index_posres(nposres(2),2) = index_atom
           end do
        end do

     else if (watertyp_free(iword) == 'posresz') then
        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        do i = 1, nwater
           do j = 1, 3
              nposres(3) = nposres(3) + 1
              index_atom = index_atom + 1
              index_posres(nposres(3),3) = index_atom
           end do
        end do

!    constraint for specific atoms
     else if (watertyp_free(iword) == 'posresatm') then
        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        read(watertyp_free(iword+1),*) nresatm
        do i = 1, nresatm
           read(watertyp_free(iword+1+i),*) index_resatm(i)
        end do

        do i = 1, nwater
           do j = 1, nresatm
              nposres(1:3) = nposres(1:3) + 1
              index_atom_tmp = index_atom + index_resatm(j)
              index_posres(nposres(1),1) = index_atom_tmp
              index_posres(nposres(2),2) = index_atom_tmp
              index_posres(nposres(3),3) = index_atom_tmp
           end do
           index_atom = index_atom + 3
        end do

        exit

     else if (watertyp_free(iword) == 'posresatmx') then
        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        read(watertyp_free(iword+1),*) nresatm
        do i = 1, nresatm
           read(watertyp_free(iword+1+i),*) index_resatm(i)
        end do

        do i = 1, nwater
           do j = 1, nresatm
              nposres(1) = nposres(1) + 1
              index_atom_tmp = index_atom + index_resatm(j)
              index_posres(nposres(1),1) = index_atom_tmp
           end do
           index_atom = index_atom + 3
        end do

     else if (watertyp_free(iword) == 'posresatmy') then
        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        read(watertyp_free(iword+1),*) nresatm
        do i = 1, nresatm
           read(watertyp_free(iword+1+i),*) index_resatm(i)
        end do

        do i = 1, nwater
           do j = 1, nresatm
              nposres(2) = nposres(2) + 1
              index_atom_tmp = index_atom + index_resatm(j)
              index_posres(nposres(2),2) = index_atom_tmp
           end do
           index_atom = index_atom + 3
        end do

     else if (watertyp_free(iword) == 'posresatmz') then
        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        read(watertyp_free(iword+1),*) nresatm
        do i = 1, nresatm
           read(watertyp_free(iword+1+i),*) index_resatm(i)
        end do

        do i = 1, nwater
           do j = 1, nresatm
              nposres(3) = nposres(3) + 1
              index_atom_tmp = index_atom + index_resatm(j)
              index_posres(nposres(3),3) = index_atom_tmp
           end do
           index_atom = index_atom + 3
        end do

     end if

  end do

!---- make list for restraint atoms of matom type

  do imatom = 1, nmatyp
     do iword = 1, maxnword

!       constraint for all atoms
        if (matomtyp_free(imatom,iword) == 'posres') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do
           index_atom = index_atom + nwater*3
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           do i = 1, nmatomtyp(imatom)
              nposres(1:3) = nposres(1:3) + 1
              index_atom = index_atom + 1
              index_posres(nposres(1),1) = index_atom
              index_posres(nposres(2),2) = index_atom
              index_posres(nposres(3),3) = index_atom
           end do

           exit

        else if (matomtyp_free(imatom,iword) == 'posresx') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do
           index_atom = index_atom + nwater*3
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           do i = 1, nmatomtyp(imatom)
              nposres(1) = nposres(1) + 1
              index_atom = index_atom + 1
              index_posres(nposres(1),1) = index_atom
           end do

        else if (matomtyp_free(imatom,iword) == 'posresy') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do
           index_atom = index_atom + nwater*3
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           do i = 1, nmatomtyp(imatom)
              nposres(2) = nposres(2) + 1
              index_atom = index_atom + 1
              index_posres(nposres(2),2) = index_atom
           end do

        else if (matomtyp_free(imatom,iword) == 'posresz') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do
           index_atom = index_atom + nwater*3
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           do i = 1, nmatomtyp(imatom)
              nposres(3) = nposres(3) + 1
              index_atom = index_atom + 1
              index_posres(nposres(3),3) = index_atom
           end do

        else if (matomtyp_free(imatom,iword) == 'posresatm') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do
           index_atom = index_atom + nwater*3
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           read(matomtyp_free(imatom,iword+1),*) nresatm
           do i = 1, nresatm
              read(matomtyp_free(imatom,iword+1+i),*) index_resatm(i)
           end do

           do j = 1, nresatm
              nposres(1:3) = nposres(1:3) + 1
              index_atom_tmp = index_atom + index_resatm(j)
              index_posres(nposres(1),1) = index_atom_tmp
              index_posres(nposres(2),2) = index_atom_tmp
              index_posres(nposres(3),3) = index_atom_tmp
           end do

           exit

        else if (matomtyp_free(imatom,iword) == 'posresatmx') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do
           index_atom = index_atom + nwater*3
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           read(matomtyp_free(imatom,iword+1),*) nresatm
           do i = 1, nresatm
              read(matomtyp_free(imatom,iword+1+i),*) index_resatm(i)
           end do

           do j = 1, nresatm
              nposres(1) = nposres(1) + 1
              index_atom_tmp = index_atom + index_resatm(j)
              index_posres(nposres(1),1) = index_atom_tmp
           end do

        else if (matomtyp_free(imatom,iword) == 'posresatmy') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do
           index_atom = index_atom + nwater*3
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           read(matomtyp_free(imatom,iword+1),*) nresatm
           do i = 1, nresatm
              read(matomtyp_free(imatom,iword+1+i),*) index_resatm(i)
           end do

           do j = 1, nresatm
              nposres(2) = nposres(2) + 1
              index_atom_tmp = index_atom + index_resatm(j)
              index_posres(nposres(2),2) = index_atom_tmp
           end do

        else if (matomtyp_free(imatom,iword) == 'posresatmz') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do
           index_atom = index_atom + nwater*3
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           read(matomtyp_free(imatom,iword+1),*) nresatm
           do i = 1, nresatm
              read(matomtyp_free(imatom,iword+1+i),*) index_resatm(i)
           end do

           do j = 1, nresatm
              nposres(3) = nposres(3) + 1
              index_atom_tmp = index_atom + index_resatm(j)
              index_posres(nposres(3),3) = index_atom_tmp
           end do

        end if

     end do
  end do

!     +     +     +     +     +     +     +

end subroutine preposres
