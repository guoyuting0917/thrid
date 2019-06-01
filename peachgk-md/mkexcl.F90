!**************************************
!*  mkexcl.f90 Ver.1.7 '13.11.05      *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine mkexcl()

  use interface_tools

  use md_global

  implicit none

!     This subroutine makes excluded neighbor list
!     (1-2, 1-3, 1-4 pairs)
!     by using BOND, ANGLE, TORSION PAIRS
 
! ARGUMENTS

!     Excluded atom index of atom(i) 
!      begins at list_excl(index_excl(i)), and  
!      ends at list_excl(index_excl(i+1)-1)  
!      If atom(i) has no excluded neighbors, 
!      0 is put in the list_excl()

!     FUNCTION:
!      isort(idata,ndata)
!      iskip(idata,ndata)

!     LOCAL
  integer:: i,j,j1,j2       ! do loop index

  integer:: iexcltmp(maxnatom) ! temporary storage for excluded atoms 
                                !  for atom (i)
  integer:: nexcltmp        ! temporal number of excluded atoms 
                                !  for atom (i)

  integer:: itmp
  logical:: foundi
  integer:: list_tmp(maxnatom*maxnmol)
  integer:: index_tmp(maxnatom)

!     +     +     +     +     +     +

!   --- INITIALIZATION ---

  list_excl(1:maxnatom*maxnmol) = 0

  index_excl(1:natom+1) = 0
  nlistexcl = 0

!   --- LOOP OVER ATOMS --- 

  DOATOM: DO i = 1, natom

     iexcltmp(1:maxnatom) = 0     ! temporary storage for excluded atoms
     nexcltmp = 0 ! temporal number of excluded atoms for atom (i)

!       - loop over BOND -

     do j = 1, nbond 
        if (i == ibond(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = jbond(j)
        else if (i == jbond(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = ibond(j)
        end if
     end do

!       - loop over ANGLE -

     do j = 1, nangl
        if (i == iangl(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = kangl(j)
        else if (i == kangl(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = iangl(j)
        end if
     end do

!       - loop over Urey-Bradley ANGLE -

     do j = 1, nanglub
        if (i == ianglub(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = kanglub(j)
        else if (i == kanglub(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = ianglub(j)
        end if
     end do

!       - loop over TORSION -

     do j = 1, ntors 
        if (i == itors(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = ltors(j)
        else if (i == ltors(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = itors(j)
        end if
     end do

!       - loop over TORSION_RB -

     do j = 1, ntorsrb
        if (i == itorsrb(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = ltorsrb(j)
        else if (i == ltorsrb(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = itorsrb(j)
        end if
     end do

!       - loop over TORSION_IM -

     do j = 1, ntorsim
        if (i == itorsim(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = ltorsim(j)
        else if (i == ltorsim(j)) then
           nexcltmp = nexcltmp + 1
           iexcltmp(nexcltmp) = itorsim(j)
        end if
     end do

!       - even if no excluded atom has been found for atom(i),
!          nexcltmp = 1 and iexcltmp(nexcltmp) = 0  

     if (nexcltmp == 0) then
        nexcltmp = 1
        iexcltmp(nexcltmp) = 0
     end if

!       - sort iexclemtp and skip duplicated ones -


     if (nexcltmp > 1) then
        call isort( iexcltmp, nexcltmp ) 
        call iskip( iexcltmp, nexcltmp ) 
            ! iskip debugged  
     end if


!       - put excluded list for atom(I) into list_excluded 
!         and set index_excluded - 

     index_excl(i) = nlistexcl + 1
     do j = 1, nexcltmp
        list_excl(nlistexcl+j) = iexcltmp(j)
     end do
     nlistexcl = nlistexcl + nexcltmp


  END DO DOATOM

  index_excl(natom+1) = nlistexcl + 1

! --- REMAKE EXCLUDED LIST ---
!
!       the excluded list from partop (list_excl) contains
!       data for both i>j and i<j for use by GRAPE . But i<j is not 
!       needed for GENERAL COMPUTERS.
!       so, list_excl is reconstructed so that it contains only
!       j >i
!       
!       use list_excl(nexcl_tot)
!           atom(natom_tot)%index_excl
!       make list_tmp(nexcl_tot)
!            index_tmp(natom_tot)
!       list_tmp and index_tmp are finally transferred 
!            to list_excl and index_excl



!    -- make index_tmp and list_tmp --

  itmp = 0    ! number excluded atoms (new)

  DO i = 1, natom
     j1 = index_excl(i)     ! position at which excle for i begins
     j2 = index_excl(i+1) - 1 ! position at which excle for i ends
     foundi = .false.       ! set true if at least one excluded
                                ! atom for i is found
     DO j = j1, j2
        if (list_excl(j) .gt. i) then
           itmp           = itmp + 1
           list_tmp(itmp) = list_excl(j)
           if (.not.foundi) then
!     - this is the first excluded atom for i -
              index_tmp(i) = itmp           
              foundi       = .true.
           end if
        end if
     END DO
         
!     - if no j was found for i, still 0 is input -
         
     if (.not.foundi) then
        itmp           = itmp + 1
        index_tmp(i)   = itmp
        list_tmp(itmp) = 0
     end if
         
  END DO

  nlistexcl = itmp 
  index_tmp(natom + 1) = nlistexcl + 1 

!    -- transfer index_tmp and list_tmp to
!        index_excl and list_excl --

  index_excl(1:natom+1) = index_tmp(1:natom+1)

  list_excl(1:nlistexcl) = list_tmp(1:nlistexcl)

!     +     +     +     +     +     +

end subroutine mkexcl
