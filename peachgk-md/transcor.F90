!**************************************
!*  transcor.f Ver.1.3 '09.12.01      *
!*     for peachgk_md.f               *
!*            by G.Kikugawa           *
!**************************************      
subroutine transcor( npoly, nwater, nmatom, &
     &               xcel, ycel, zcel )

  use md_global

  implicit none

!     This subroutine translate the coordinates in the case of 
!     periodic boundry condition (PBC).
!     PBC is applied according to the first atom of each molecule;
!     namely, all of the atoms belonging to the molecule will be 
!     translated if the first atom is outside the BOX.
!     Attention: the BOX is supposed to be centered at the origin.  
!
! ARGUMENTS: 
!     INPUT
  integer,intent(in):: npoly           ! number of polymer1
  integer,intent(in):: nwater          ! number of H2O molecules
  integer,intent(in):: nmatom          ! number of monatomic molecules

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

! LOCAL:
  real(8):: box(3)           ! box cell length
      
  integer:: i,j,k           ! do loop index
  integer:: j1, j2

!     +     +     +     +     +     +     +


!     --- SET BOXMIN & BOXMAX SO THAT THE BOX SHOULD BE CENTERED  
!             AT THE ORIGIN ---

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel

!     --- TRANSLATE THE COORDINATES ---
 
!---- molecule polymer1
  DO i = 1, npoly

!        - find the first (j1) & last (j2) atoms  of mol(i) ---

     j1 = molept_index(i)       ! the first atom of mol(i)  
     j2 = molept_index(i+1) - 1 ! the first atom of mol(i+1)
         
!        - translate the coordinates -

!        loop over X, Y, Z
!        - x
     if (atmcor(1,molept_list(j1)) < 0.0d0 ) then

        do j=j1,j2
           k = molept_list(j)
           atmcor(1,k) = atmcor(1,k) + box(1)
        end do

     else if (atmcor(1,molept_list(j1)) > box(1) ) then

        do j=j1,j2
           k = molept_list(j)                  
           atmcor(1,k) = atmcor(1,k) - box(1)
        end do

     end if
!        - y
     if (atmcor(2,molept_list(j1)) < 0.0d0 ) then

        do j=j1,j2
           k = molept_list(j)
           atmcor(2,k) = atmcor(2,k) + box(2)
        end do

     else if (atmcor(2,molept_list(j1)) > box(2) ) then

        do j=j1,j2
           k = molept_list(j)                  
           atmcor(2,k) = atmcor(2,k) - box(2)
        end do

     end if
!        - z
     if (atmcor(3,molept_list(j1)) < 0.0d0 ) then

        do j=j1,j2
           k = molept_list(j)
           atmcor(3,k) = atmcor(3,k) + box(3)
        end do

     else if (atmcor(3,molept_list(j1)) > box(3) ) then

        do j=j1,j2
           k = molept_list(j)                  
           atmcor(3,k) = atmcor(3,k) - box(3)
        end do

     end if

  END DO

!---- molecule water
  DO i = npoly+1, npoly+nwater

!        - find the first (j1) & last (j2) atoms  of mol(i) ---

     j1 = molept_index(i)       ! the first atom of mol(i)  
     j2 = molept_index(i+1) - 1 ! the first atom of mol(i+1)
         
!        - translate the coordinates -

!        loop over X, Y, Z
!        - x
     if (atmcor(1,molept_list(j1)) < 0.0d0 ) then
            
        do j=j1,j2
           k = molept_list(j)
           atmcor(1,k) = atmcor(1,k) + box(1)
        end do

     else if (atmcor(1,molept_list(j1)) > box(1) ) then

        do j=j1,j2
           k = molept_list(j)                  
           atmcor(1,k) = atmcor(1,k) - box(1)
        end do
               
     end if
!        - y
     if (atmcor(2,molept_list(j1)) < 0.0d0 ) then
            
        do j=j1,j2
           k = molept_list(j)
           atmcor(2,k) = atmcor(2,k) + box(2)
        end do

     else if (atmcor(2,molept_list(j1)) > box(2) ) then

        do j=j1,j2
           k = molept_list(j)                  
           atmcor(2,k) = atmcor(2,k) - box(2)
        end do
               
     end if
!        - z
     if (atmcor(3,molept_list(j1)) < 0.0d0 ) then
            
        do j=j1,j2
           k = molept_list(j)
           atmcor(3,k) = atmcor(3,k) + box(3)
        end do

     else if (atmcor(3,molept_list(j1)) > box(3) ) then

        do j=j1,j2
           k = molept_list(j)                  
           atmcor(3,k) = atmcor(3,k) - box(3)
        end do
               
     end if
         
  END DO

!---- atom pt
  DO i = npoly+nwater+1, npoly+nwater+nmatom

!        - find the first (k1) & last (k2) atoms  of mol(i) ---

     j1 = molept_index(i)       ! the first atom of mol(i)  
     j2 = molept_index(i+1) - 1 ! the first atom of mol(i+1)
         
!        - translate the coordinates -

!        loop over X, Y, Z
!        - x
     if (atmcor(1,molept_list(j1)) < 0.0d0 ) then

        do j=j1,j2
           k = molept_list(j)
           atmcor(1,k) = atmcor(1,k) + box(1)
        end do

     else if (atmcor(1,molept_list(j1)) > box(1) ) then

        do j=j1,j2
           k = molept_list(j)                  
           atmcor(1,k) = atmcor(1,k) - box(1)
        end do

     end if
!        - y
     if (atmcor(2,molept_list(j1)) < 0.0d0 ) then

        do j=j1,j2
           k = molept_list(j)
           atmcor(2,k) = atmcor(2,k) + box(2)
        end do

     else if (atmcor(2,molept_list(j1)) > box(2) ) then

        do j=j1,j2
           k = molept_list(j)                  
           atmcor(2,k) = atmcor(2,k) - box(2)
        end do

     end if
!        - z
     if (atmcor(3,molept_list(j1)) < 0.0d0 ) then

        do j=j1,j2
           k = molept_list(j)
           atmcor(3,k) = atmcor(3,k) + box(3)
        end do

     else if (atmcor(3,molept_list(j1)) > box(3) ) then

        do j=j1,j2
           k = molept_list(j)                  
           atmcor(3,k) = atmcor(3,k) - box(3)
        end do

     end if

  END DO

!     +     +     +     +     +     +     +     +

  return
end subroutine transcor
