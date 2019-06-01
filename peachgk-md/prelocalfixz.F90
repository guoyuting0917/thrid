!**************************************
!*  prelocalfixz.f Ver.1.1 '11.01.31  *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine prelocalfixz(xref,   &
     &                  npolytyp, npoly_mole, npoly_atom,   &
     &                  nwater,   &
     &                  nmatyp,nmatomtyp,   &
     &                  polytyp_free, watertyp_free, matomtyp_free,   &
     &                  nlfixz_deg_poly, nlfixz_deg_water, nlfixz_deg_ma)

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: xref        ! distanse base value [m] (adjust to Angstrom)

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

!    OUTPUT
  integer,intent(out):: nlfixz_deg_poly(:) ! fixed degree of freedom of poly
  integer,intent(out):: nlfixz_deg_water   ! fixed degree of freedom of H2O
  integer,intent(out):: nlfixz_deg_ma(:)   ! fixed degree of freedom of matom

! LOCAL:
  integer:: ipoly
  integer:: imatom
  integer:: iword
      
  integer:: index_atom

  integer:: i

  integer:: ilfixz_atom

  real(8):: cor_lfixz_tmp

!     +     +     +     +     +     +

!---- initialization

  nlfixz = 0
  nlfixz_deg_poly(1:npolytyp) = 0
  nlfixz_deg_water = 0
  nlfixz_deg_ma(1:nmatyp) = 0

!---- make list for fixed atoms and z-coordinate of poly type

  do ipoly = 1, npolytyp
     do iword = 1, maxnword

        if (polytyp_free(ipoly,iword) == 'localfixz') then

           read(polytyp_free(ipoly,iword+1),*) ilfixz_atom
           read(polytyp_free(ipoly,iword+2),*) cor_lfixz_tmp

           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           do i = 1, npoly_mole(ipoly)
              nlfixz = nlfixz + 1
              index_nlfixz(nlfixz) = index_atom + npoly_atom(ipoly)*(i-1)   &
                   &               + ilfixz_atom
              cor_lfixz(nlfixz) = cor_lfixz_tmp / xref

              nlfixz_deg_poly(ipoly) = nlfixz_deg_poly(ipoly) + 1
           end do

           exit
           
        end if

     end do
  end do

!---- make list for fixed atoms and z-coordinate of water type

  do iword = 1, maxnword

     if (watertyp_free(iword) == 'localfixz') then

        read(watertyp_free(iword+1),*) ilfixz_atom
        read(watertyp_free(iword+2),*) cor_lfixz_tmp

        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do
        
        do i = 1, nwater
           nlfixz = nlfixz + 1
           index_nlfixz(nlfixz) = index_atom + 3*(i-1) + ilfixz_atom
           cor_lfixz(nlfixz) = cor_lfixz_tmp / xref

           nlfixz_deg_water = nlfixz_deg_water + 1
        end do
                     
        exit
        
     end if

  end do

!---- make list for fixed atoms and z-coordinate of matom type

  do imatom = 1, nmatyp
     do iword = 1, maxnword

        if (matomtyp_free(imatom,iword) == 'localfixz') then

           read(matomtyp_free(imatom,iword+1),*) cor_lfixz_tmp

           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do
           index_atom = index_atom + nwater*3
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           do i = 1, nmatomtyp(imatom)
              nlfixz = nlfixz + 1
              index_nlfixz(nlfixz) = index_atom + i
              cor_lfixz(nlfixz) = cor_lfixz_tmp / xref

              nlfixz_deg_ma(imatom) = nlfixz_deg_ma(imatom) + 1
           end do
                  
           exit
           
        end if

     end do
  end do

#if defined(_PRELFIXZ_DEBUG)
  write(6,*) '***** localfixz debug info *****'
  write(6,*) 'number of atoms of which z-coordinate are fixed = ', nlfixz
  write(6,*) '*** atom indexes and position'
  do i = 1, nlfixz
     write(6,*) index_nlfixz(i), cor_lfixz(i)
  end do
#endif

!     +     +     +     +     +     +     +

  return
end subroutine prelocalfixz
