!**************************************
!*  prelocalfix.f Ver.1.2 '10.06.30   *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine prelocalfix( npolytyp, npoly_mole, npoly_atom, &
     &                  nwater, &
     &                  nmatyp, nmatomtyp, &
     &                  polytyp_free, watertyp_free, matomtyp_free, &
     &                  nlfix,index_nlfix, &
     &                  nlfix_deg_poly, &
     &                  nlfix_deg_water, &
     &                  nlfix_deg_ma )

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

!     OUTPUT
  integer,intent(out):: nlfix           ! number of fix atoms
  integer,intent(out):: index_nlfix(:)  ! index of fix atoms
 
  integer,intent(out):: nlfix_deg_poly(:) ! fixed degree of freedom of poly
  integer,intent(out):: nlfix_deg_water ! fixed degree of freedom of H2O
  integer,intent(out):: nlfix_deg_ma(:) ! fixed degree of freedom of matom

! LOCAL:
  integer:: ipoly
  integer:: imatom
  integer:: iword
      
  integer:: index_atom

  integer:: i,j

!     +     +     +     +     +     +

!---- initialization
  nlfix_deg_poly(1:npolytyp) = 0
  nlfix_deg_water = 0
  nlfix_deg_ma(1:nmatyp) = 0

!---- make list for fix atoms of poly type

  nlfix = 0
  do ipoly = 1, npolytyp
     do iword = 1, maxnword

        if (polytyp_free(ipoly,iword) == 'localfix') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, npoly_atom(ipoly)
                 nlfix = nlfix + 1
                 index_atom = index_atom + 1
                 index_nlfix(nlfix) = index_atom

                 nlfix_deg_poly(ipoly) = nlfix_deg_poly(ipoly) + 3
              end do
           end do

           exit
        end if

     end do
  end do

!---- make list for fix atoms of water type

  do iword = 1, maxnword

     if (watertyp_free(iword) == 'localfix') then
        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        do i = 1, nwater
           do j = 1, 3
              nlfix = nlfix + 1
              index_atom = index_atom + 1
              index_nlfix(nlfix) = index_atom
              nlfix_deg_water = nlfix_deg_water + 3
           end do
        end do
                     
        exit
     end if

  end do

!---- make list for fix atoms of matom type

  do imatom = 1, nmatyp
     do iword = 1, maxnword

        if (matomtyp_free(imatom,iword) == 'localfix') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           index_atom = index_atom + nwater*3
               
           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           do i = 1, nmatomtyp(imatom)
              nlfix = nlfix + 1
              index_atom = index_atom + 1
              index_nlfix(nlfix) = index_atom

              nlfix_deg_ma(imatom) = nlfix_deg_ma(imatom) + 3
           end do
                  
           exit
        end if

     end do
  end do

!     +     +     +     +     +     +     +

  return
end subroutine prelocalfix
