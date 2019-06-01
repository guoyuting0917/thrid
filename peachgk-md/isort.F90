subroutine isort(idata, ndata)
!
!     subroutine to sort IDATA. 
!     so that idata(1) =< idata(2) =< idata(3) .... 
!     by BUBBLE SORT

  implicit none

! ARGUMENT:

!     INPUT
  integer,intent(in)::   ndata         ! number of data 

!     INPUT & OUTPUT 
  integer,intent(inout)::   idata(ndata)  ! data array to sort

!     LOCAL

  logical::   swapped       ! set true if swapped

  integer::   itmp          ! temporal storage for swapping
  integer::   i             ! do loop index


!     --------------------------------------------------------------

      
  DO 

     swapped = .false.
     do i = 1, ndata - 1
        if (idata(i)>idata(i+1)) then
           itmp       = idata(i+1)
           idata(i+1) = idata(i)
           idata(i)   = itmp
           swapped    = .true.
        end if
     end do
     if (.not.swapped) EXIT

  END DO

  return
end subroutine isort
