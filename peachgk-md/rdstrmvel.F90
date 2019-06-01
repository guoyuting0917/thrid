!*****************************
!*  rdstrmvel.f90 Ver.1.0    *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-02-17 12:04:32 gota>

subroutine rdstrmvel(iustrmvel, &
     &               xref,vref, &
     &               xcel,ycel,zcel)

  use interface_tools
  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: iustrmvel       ! input streaming velocity file unit

  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: vref             ! velocity base value [m/s]

  real(8),intent(in):: xcel             ! initial x0 cell length[non-d]
  real(8),intent(in):: ycel             ! initial y0 cell length[non-d]
  real(8),intent(in):: zcel             ! initial z0 cell length[non-d]

!     OUTPUT

! LOCAL:
  character(80):: fredat(20)

  integer:: nword

  integer:: i               ! do loop index

  integer:: index_reg

  integer:: read_sta        ! indication of reading position (one of following)
  integer:: rd_head = 1     ! header part 1
  integer:: rd_head2 = 2    ! header part 2
  integer:: rd_bhead = 3    ! block header
  integer:: rd_data = 4     ! data block

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- intialize variable --------

  read_sta = rd_head

!---- read parameter from streming velocity file (out_strmvel.dat)

  DO
     call rdfree_w(iustrmvel,20,fredat,nword)

     !- header part 1
     if (read_sta == rd_head) then

        if (fredat(1) == ' ') then   ! blank line
           read_sta = rd_head2
           index_reg = 0
           cycle

        else if (fredat(1) == '#nint_lvs') then  ! find nint_lvs parameter
           read(fredat(2),*) nvelregion

        end if

     !- header part 2
     else if (read_sta == rd_head2) then

        if (fredat(1) == ' ') then   ! blank line
           read_sta = rd_bhead
           cycle

        else if (fredat(1) == '#region') then   ! first line of header 2
           cycle

        end if

        index_reg = index_reg + 1
        read(fredat(2),*) strmzpos1(index_reg)
        read(fredat(4),*) strmzpos2(index_reg)
        strmzpos1(index_reg) = strmzpos1(index_reg) / xref
        strmzpos2(index_reg) = strmzpos2(index_reg) / xref

     !- block header
     else if (read_sta == rd_bhead) then

        call rdfree_w(iustrmvel,20,fredat,nword) ! read extra line

        read_sta = rd_data
        index_reg = 0
        cycle

     !- data block
     else if (read_sta == rd_data) then

        if (fredat(1) == ' ') then   ! blank line
           exit

        else if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) then
           cycle                     ! comment line

        end if

        index_reg = index_reg + 1
        read(fredat(3),*) strmvel(1,index_reg)
        read(fredat(4),*) strmvel(2,index_reg)
        read(fredat(5),*) strmvel(3,index_reg)
        strmvel(1:3,index_reg) = strmvel(1:3,index_reg) / vref

     end if

  END DO

!---- check the parameter relations and so on

  if (index_reg /= nvelregion) then
     write(6,*) 'Error: number of slabs in header block and', &
          &     ' in actual data block does not match.'
     stop
  end if

  do i = 1, nvelregion

     if ((strmzpos1(i) < 0.0d0) .or. (strmzpos1(i) > zcel)) then
        write(6,*) 'Error: strmzpos1(',i,') exceeds cell dimension'
        stop
     end if
     if ((strmzpos2(i) < 0.0d0) .or. (strmzpos2(i) > zcel)) then
        write(6,*) 'Error: strmzpos2(',i,') exceeds cell dimension'
        stop
     end if
     if (strmzpos1(i) >= strmzpos2(i)) then
        write(6,*) 'Error: strmzpos1(',i,') exceeds strmzpos2(',i,')'
        stop
     end if

  end do

#if defined(_RDSTRMVEL_DEBUG)
  write(6,*) '***** rdstrmvel debug info *****'
  do i=1,nvelregion
     write(6,*) 'region No. ',i,': ',strmzpos1(i)*xref,' <=> ', &
          &     strmzpos2(i)*xref,' : ', &
          &     strmvel(1,i)*vref,strmvel(2,i)*vref,strmvel(3,i)*vref
  end do
#endif

  close(iustrmvel)

!     +     +     +     +     +     +     +

end subroutine rdstrmvel
