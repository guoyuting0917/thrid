!**************************************
!*  prepotbias.f Ver.1.1 '10.06.30    *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine prepotbias( npolytyp, npoly_mole, npoly_atom, &
     &                 nwater, &
     &                 nmatyp,nmatomtyp, &
     &                 polytyp_free, watertyp_free, matomtyp_free)

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

  integer:: nbiasatm
  integer:: index_biasatm(maxnmol)

!     +     +     +     +     +     +

!---- initialization

  npotbias = 0

!---- make list for atoms imposed by bias potential of poly type

  do ipoly = 1, npolytyp
     do iword = 1, maxnword

        if (polytyp_free(ipoly,iword) == 'potbias') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, npoly_atom(ipoly)
                 npotbias = npotbias + 1
                 index_atom = index_atom + 1
                 index_potbias(npotbias) = index_atom
              end do
           end do

           exit

        else if (polytyp_free(ipoly,iword) == 'potbiasatm') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           read(polytyp_free(ipoly,iword+1),*) nbiasatm
           do i = 1, nbiasatm
              read(polytyp_free(ipoly,iword+1+i),*) index_biasatm(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, nbiasatm
                 npotbias = npotbias + 1
                 index_atom_tmp = index_atom + index_biasatm(j)
                 index_potbias(npotbias) = index_atom_tmp
              end do
              index_atom = index_atom + npoly_atom(ipoly)
           end do

           exit

        end if

     end do
  end do

!---- make list for atoms imposed by bias potential of water type

  do iword = 1, maxnword

     if (watertyp_free(iword) .eq. 'potbias') then
        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        do i = 1, nwater
           do j = 1, 3
              npotbias = npotbias + 1
              index_atom = index_atom + 1
              index_potbias(npotbias) = index_atom
           end do
        end do
                     
        exit

     else if (watertyp_free(iword) .eq. 'potbiasatm') then
        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        read(watertyp_free(iword+1),*) nbiasatm
        do i = 1, nbiasatm
           read(watertyp_free(iword+1+i),*) index_biasatm(i)
        end do

        do i = 1, nwater
           do j = 1, nbiasatm
              npotbias = npotbias + 1
              index_atom_tmp = index_atom + index_biasatm(j)
              index_potbias(npotbias) = index_atom_tmp
           end do
           index_atom = index_atom + 3
        end do
                     
        exit

     end if

  end do

!---- make list for atoms imposed by bias potential of matom type

  do imatom = 1, nmatyp
     do iword = 1, maxnword

        if (matomtyp_free(imatom,iword) .eq. 'potbias') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do
           index_atom = index_atom + nwater*3
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           do i = 1, nmatomtyp(imatom)
              npotbias = npotbias + 1
              index_atom = index_atom + 1
              index_potbias(npotbias) = index_atom
           end do
                  
           exit

        else if (matomtyp_free(imatom,iword) .eq. 'potbiasatm') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do
           index_atom = index_atom + nwater*3
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           read(matomtyp_free(imatom,iword+1),*) nbiasatm
           do i = 1, nbiasatm
              read(matomtyp_free(imatom,iword+1+i),*) index_biasatm(i)
           end do

           do j = 1, nbiasatm
              npotbias = npotbias + 1
              index_atom_tmp = index_atom + index_biasatm(j)
              index_potbias(npotbias) = index_atom_tmp
           end do

           exit

        end if

     end do
  end do

#if defined(_POTBIAS_DEBUG)
  write(6,*) '***** potbias debug info *****'
  write(6,*) 'number of atoms imposed by bias potential = ', npotbias
  write(6,*) '*** atom indexes'
  do i = 1, npotbias
     write(6,*) index_potbias(i)
  end do
#endif

!     +     +     +     +     +     +     +

  return
end subroutine prepotbias
