!*****************************
!*  wrsta.f90 Ver.2.1        *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>

subroutine wrsta(iostarec, &
     &           npoly,nwater,nmatom, &
     &           xcel,ycel,zcel, &
     &           xref,vref,timeref,pref, &
     &           maxnstep, &
     &           mchain, &
     &           pint, pintt)

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: iostarec        ! state record file unit

  integer,intent(in):: npoly           ! number of polymer1
  integer,intent(in):: nwater          ! number of H2O molecules
  integer,intent(in):: nmatom          ! number of monatomic molecules

  real(8),intent(in):: xcel             ! x cell length
  real(8),intent(in):: ycel             ! y cell length
  real(8),intent(in):: zcel             ! z cell length

  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: vref             ! velocity base value [m/s]
  real(8),intent(in):: timeref          ! time base value [sec]
  real(8),intent(in):: pref             ! pressure base value [Pa]

  integer,intent(in):: maxnstep        ! current or maximum step of MD

  integer,intent(in):: mchain          ! Nose-Hoover chain number

  real(8),intent(in):: pint             ! internal pressure
  real(8),intent(in):: pintt(:,:)       ! internal pressure tensor

! LOCAL:
  integer:: i,j             ! do loop index
  integer:: j1,j2,jj
  integer:: imole

!     +     +     +     +     +     +     +     +

  write(6,*) 'Write atmcor&vel to state_file'

!---- writing some parameter
  rewind(iostarec)

  write(iostarec,'(A6,I8)') '#npoly',npoly
  write(iostarec,*)
  write(iostarec,'(A7,I8)') '#nwater',nwater
  write(iostarec,*)
  write(iostarec,'(A7,I8)') '#nmatom',nmatom
  write(iostarec,*)
  write(iostarec,'(A6,I10)') '#natom',natom
  write(iostarec,*)
  write(iostarec,'(A7,I3)') '#mchain',mchain
  write(iostarec,*)
  write(iostarec,'(A5,E24.16)') '#xcel',xcel*xref
  write(iostarec,*)
  write(iostarec,'(A5,E24.16)') '#ycel',ycel*xref
  write(iostarec,*)
  write(iostarec,'(A5,E24.16)') '#zcel',zcel*xref
  write(iostarec,*)
  write(iostarec,'(A6,I12)') '#nstep',maxnstep
  write(iostarec,*)

!---- Initialization of valiables
  imole = 0

!---- writing atmcor&vel of poly1
  write(iostarec,'(A6)') '#poly'

  do i = 1, npoly
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1

     do j = j1,j2
        jj = molept_list(j)

        imole = imole + 1
        write(iostarec,'(I10,6E24.16,A3)') &
             & imole, &
             & atmcor(1,jj)*xref,atmcor(2,jj)*xref,atmcor(3,jj)*xref, &
             & atmvel(1,jj)*vref,atmvel(2,jj)*vref,atmvel(3,jj)*vref, &
             & atmtyp(jj)

     end do

  end do

  write(iostarec,*)

!---- writing atmcor&vel of water
  write(iostarec,'(A6)') '#water'

  do i = npoly+1, npoly+nwater
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1

     do j = j1,j2
        jj = molept_list(j)

        imole = imole + 1
        write(iostarec,'(I10,6E24.16,A3)') &
             & imole, &
             & atmcor(1,jj)*xref,atmcor(2,jj)*xref,atmcor(3,jj)*xref, &
             & atmvel(1,jj)*vref,atmvel(2,jj)*vref,atmvel(3,jj)*vref, &
             & atmtyp(jj)

     end do

  end do

  write(iostarec,*)

!---- writing atmcor&vel of MA
  write(iostarec,'(A6)') '#matom'

  do i = npoly+nwater+1, npoly+nwater+nmatom
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1

     do j = j1,j2
        jj = molept_list(j)

        imole = imole + 1
        write(iostarec,'(I10,6E24.16,A3)') &
             & imole, &
             & atmcor(1,jj)*xref,atmcor(2,jj)*xref,atmcor(3,jj)*xref, &
             & atmvel(1,jj)*vref,atmvel(2,jj)*vref,atmvel(3,jj)*vref, &
             & atmtyp(jj)

     end do

  end do

  write(iostarec,*)

!---- writing coordinate&velocity of Nose-Hoover chain
  write(iostarec,'(A4)') '#nhc'

  do i = 1,mchain
     write(iostarec,'(I10,2E24.16)') &
          & i,xlogs(i),vlogs(i)/timeref

  end do

  write(iostarec,*)

!---- writing velocity of Andersen barostat
  write(iostarec,'(A9)') '#andersen'

  write(iostarec,'(I10,2E24.16)') &
       & 1,vlogv/timeref,pint*pref

  write(iostarec,*)

!---- writing velocity of Andersen barostat
  write(iostarec,'(A15)') '#andersen_aniso'

  write(iostarec,'(I10,3E24.16,9E24.16)') &
       & 1,vboxg(1:3)/timeref,pintt(1:3,1:3)*pref

  write(iostarec,*)

!---- writing end mark
  write(iostarec,'(A4)') '#end'

!     +     +     +     +     +     +     +     +

  return
end subroutine wrsta
