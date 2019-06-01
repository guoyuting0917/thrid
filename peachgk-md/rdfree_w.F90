!-----------------------------------------------------------------------
subroutine rdfree_w(iuin, ndata, fredat, nword)
!     
!     subroutine to read in a free format 
!
!     This subroutine reads in a free format a line containing various
!     types of data separated by a space ' ' or a comma ','.
!     For instance,
!         'ATOM 500 RES 101  50.2 30.3 1.99 10e-10'  
!     DATA ARE RETURNED AS CHARACTERS
!     If ndata > actual number of data, then spaces are input in the
!     rest of fredat.

  implicit none

!     ARGUMENT

!        INPUT: 
  integer,intent(in)::     iuin          ! unit to read
  integer,intent(in)::     ndata         ! number of data to be read.

!        OUTPUT:
  character(80),intent(inout):: fredat(ndata) ! data read as characters
!                   length of each datum is assumed =<80
!                   each datum begins at the left end of data(i)
!
  integer,intent(inout)::     nword         ! actual number of data

!     LOCAL

!        * LINE 
  integer,parameter:: lnleng = 2000   ! line length
  character(len=lnleng):: line ! line is read as a character
 
!        * COUNTS
  integer::     idata         ! temporal count of data
  integer::     ldata         ! length of the datum
 
!        * FLAGS
  logical::     spcom0, spcom1 ! true  if i-1 (or i) is space
!                                                  or comma

!        * OHTERS
  character(len=1):: space,  comma,  tab
                                ! separators
  integer::       i           ! do loop index 

!          +    +    +    +    +


!     --- INITIALIZE THE DATA ---

  space=' '
  comma=','
  tab='	'
  do i=1,ndata
     fredat(i) = space
  end do

!     --- READ THE LINE FROM IUIN AS A LONG CHARACTER ---

  read(iuin,'(a)') line

!     --- FIND DATA SEPARATED BY SPACE OR COMMA ---

  idata = 0       ! temporal number of data
  ldata = 0
  spcom0 = .true.

  do i = 1, lnleng 
     spcom1 = (line(i:i)==space.or.line(i:i)==comma &
          &    .or.line(i:i)==tab) 
     if (spcom0) then
!           * the last character was space or comma
        if (.not.spcom1) then
!              * i is the begining of a data
           idata = idata + 1       ! temporal number of data
           if (idata.gt.ndata) then
              nword = idata
              return
           end if
           ldata = 1               ! length of the datum
           fredat(idata)(ldata:ldata) = line(i:i)
        end if
     else 
!           * the last character was not space or comma
        if (.not.spcom1) then
           ldata = ldata + 1                               
           fredat(idata)(ldata:ldata) = line(i:i)
        end if
     end if
     spcom0 = spcom1
  end do

  nword = idata
            
!          +    +    +    +    +

  return
end subroutine rdfree_w
!-----------------------------------------------------------------------
