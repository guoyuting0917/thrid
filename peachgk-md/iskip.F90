subroutine iskip(idata, ndata)
!
!     subroutine to find equal data from  IDATA. 
!     For instance,
!       1, 2, 2, 4, 5 becomes
!       1, 2, 4, 5
!
!     ATTENTION! idata is supposed to have been sorted.
!     debugged: Jan 15, 1995 for last sequence such as 1,2,5,5,5

  implicit none

! ARGUMENT:

!     INPUT & OUTPUT 
  integer,intent(inout)::   ndata         ! number of data 
  integer,intent(inout)::   idata(ndata)  ! data array to sort

! LOCAL:

  integer::   i,j           ! do loop index


!     --------------------------------------------------------------

  i = 0 ! current index of data
  DO  
     i = i + 1
     if (i>=ndata) exit

     if (idata(i)==idata(i+1)) then
        ndata = ndata - 1
        do j = i+1, ndata 
           idata(j) = idata(j+1)
        end do
        i = i -1 ! work around for such sequence as 7,7,7
     end if
  END DO

  return
end subroutine iskip
