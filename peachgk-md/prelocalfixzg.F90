!**************************************
!*  prelocalfixzg.f Ver.1.2 '12.07.18 *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!*    customized by T.Nakano          *
!**************************************
subroutine prelocalfixzg(xref,   &
     &                   npolytyp,npoly_mole,npoly_atom,   &
     &                   nwater,   &
     &                   nmatyp,nmatomtyp,   &
     &                   polytyp_free,watertyp_free,matomtyp_free,   &
     &                   nlfixzg_deg_poly,nlfixzg_deg_water,nlfixzg_deg_ma)

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
  integer,intent(out):: nlfixzg_deg_poly(:) ! fixed degree of freedom of poly
  integer,intent(out):: nlfixzg_deg_water   ! fixed degree of freedom of H2O
  integer,intent(out):: nlfixzg_deg_ma(:)   ! fixed degree of freedom of matom

! LOCAL:
  integer:: ipoly
  integer:: imatom
  integer:: iword
      
  integer:: index_mole

  integer:: i

  real(8):: cor_lfixzg_tmp

!     +     +     +     +     +     +

!---- initialization

  nlfixzg = 0
  nlfixzg_deg_poly(1:npolytyp) = 0
  nlfixzg_deg_water = 0
  nlfixzg_deg_ma(1:nmatyp) = 0

!---- make list for fixed molecules and z-coordinate of poly type

  do ipoly = 1, npolytyp
     do iword = 1, maxnword

        if (polytyp_free(ipoly,iword) == 'localfixzg') then

           read(polytyp_free(ipoly,iword+1),*) cor_lfixzg_tmp

           index_mole = 0
           do i = 1, ipoly-1
              index_mole = index_mole + npoly_mole(i)
           end do

           do i = 1, npoly_mole(ipoly)
              nlfixzg = nlfixzg + 1
              index_nlfixzg(nlfixzg) = index_mole + i
              cor_lfixzg(nlfixzg) = cor_lfixzg_tmp / xref

              nlfixzg_deg_poly(ipoly) = nlfixzg_deg_poly(ipoly) + 1
           end do

           exit
           
        end if

     end do
  end do

!---- make list for fixed molecules and z-coordinate of water type

  do iword = 1, maxnword

     if (watertyp_free(iword) == 'localfixzg') then

        read(watertyp_free(iword+1),*) cor_lfixzg_tmp

        index_mole = 0
        do i = 1, npolytyp
           index_mole = index_mole + npoly_mole(i)
        end do

        do i = 1, nwater
           nlfixzg = nlfixzg + 1
           index_nlfixzg(nlfixzg) = index_mole + i
           cor_lfixzg(nlfixzg) = cor_lfixzg_tmp / xref

           nlfixzg_deg_water = nlfixzg_deg_water + 1
        end do

        exit
        
     end if

  end do

!---- make list for fixed atoms and z-coordinate of matom type

  do imatom = 1, nmatyp
     do iword = 1, maxnword

        if (matomtyp_free(imatom,iword) == 'localfixzg') then

           read(matomtyp_free(imatom,iword+1),*) cor_lfixzg_tmp

           index_mole = 0
           do i = 1, npolytyp
              index_mole = index_mole + npoly_mole(i)
           end do
           index_mole = index_mole + nwater
           do i = 1, imatom - 1
              index_mole = index_mole + nmatomtyp(i)
           end do

           do i = 1, nmatomtyp(imatom)
              nlfixzg = nlfixzg + 1
              index_nlfixzg(nlfixzg) = index_mole + i
              cor_lfixzg(nlfixzg) = cor_lfixzg_tmp / xref

              nlfixzg_deg_ma(imatom) = nlfixzg_deg_ma(imatom) + 1
           end do
                  
           exit
           
        end if

     end do
  end do

#if defined(_PRELFIXZG_DEBUG)
  write(6,*) '***** localfixzg debug info *****'
  write(6,*) 'number of atoms of which z-coordinate are fixed = ', nlfixzg
  write(6,*) '*** atom indexes and position'
  do i = 1, nlfixzg
     write(6,*) index_nlfixzg(i), cor_lfixzg(i)
  end do
#endif

!     +     +     +     +     +     +     +

end subroutine prelocalfixzg
