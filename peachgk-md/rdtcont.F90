!*****************************
!*  rdtcont.f Ver.1.3        *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2014-06-13 16:53:21 gota>

subroutine rdtcont(xref, tempref, &
     &             xcel, ycel, zcel, &
     &             ntcregion, &
     &             tcxpos1, tcxpos2, &
     &             tcypos1, tcypos2, &
     &             tczpos1, tczpos2, &
     &             r_tcont)

  use interface_tools

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: tempref          ! temperature base value [K]

  real(8),intent(in):: xcel             ! initial x0 cell length[m]
  real(8),intent(in):: ycel             ! initial y0 cell length[m]
  real(8),intent(in):: zcel             ! initial z0 cell length[m]

!     OUTPUT
  integer,intent(out):: ntcregion       ! number of region to control temp.
  real(8),intent(out):: tcxpos1(:),tcxpos2(:)
                                ! x-position of temp. control region
  real(8),intent(out):: tcypos1(:),tcypos2(:)
                                ! y-position of temp. control region
  real(8),intent(out):: tczpos1(:),tczpos2(:)
                                ! z-position of temp. control region
  real(8),intent(out):: r_tcont(:)       ! control temp. in each region

! LOCAL:
  character(80):: fredat(20)

  integer:: nword

  integer:: i               ! do loop index

  integer:: iuini
  integer:: ios

  integer:: itcregion

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- intialize variable --------

!---- read parameter from MD script (peachgk.ini)

  iuini = 99
  open(iuini,file='tempcont.ini',status='old',iostat=ios)
  if (ios /= 0) then
     write(6,*) 'Failure in opening file: tempcont.ini'
     stop
  end if

  DO
     call rdfree_w( iuini, 20, fredat, nword)

     if (fredat(1)(1:3) == 'END') then
        exit

     else if (fredat(1) == ' ') then
        cycle

     else if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) then
        cycle                     ! comment line

!        definition of temp. control region
     else if (fredat(1) == 'ntcregion') then
        read(fredat(2),*) ntcregion

     else if (fredat(1) == 'tcregion') then
        read(fredat(2),*) itcregion
        if (itcregion > ntcregion) then
           write(6,*) 'Error: tcregion in rdtcont'
           stop
        end if

        if (fredat(3) == 'xcel') then
           tcxpos1(itcregion) = 0.0d0
        else
           read(fredat(3),*) tcxpos1(itcregion)
           tcxpos1(itcregion) = tcxpos1(itcregion) / xref
        end if

        if (fredat(4) == 'xcel') then
           tcxpos2(itcregion) = xcel
        else
           read(fredat(4),*) tcxpos2(itcregion)
           tcxpos2(itcregion) = tcxpos2(itcregion) / xref
        end if

        if (fredat(5) == 'ycel') then
           tcypos1(itcregion) = 0.0d0
        else
           read(fredat(5),*) tcypos1(itcregion)
           tcypos1(itcregion) = tcypos1(itcregion) / xref
        end if

        if (fredat(6) == 'ycel') then
           tcypos2(itcregion) = ycel
        else
           read(fredat(6),*) tcypos2(itcregion)
           tcypos2(itcregion) = tcypos2(itcregion) / xref
        end if

        if (fredat(7) == 'zcel') then
           tczpos1(itcregion) = 0.0d0
        else
           read(fredat(7),*) tczpos1(itcregion)
           tczpos1(itcregion) = tczpos1(itcregion) / xref
        end if

        if (fredat(8) == 'zcel') then
           tczpos2(itcregion) = zcel
        else
           read(fredat(8),*) tczpos2(itcregion)
           tczpos2(itcregion) = tczpos2(itcregion) / xref
        end if

        read(fredat(9),*) r_tcont(itcregion)

     end if

  END DO

!---- check the parameter relations and so on

  do i=1,ntcregion
     if ((tcxpos1(i) < 0.0d0) .or. (tcxpos1(i) > xcel)) then
        write(6,*) 'Error: tcxpos1(',i,') exceeds cell dimension'
        stop
     end if
     if ((tcxpos2(i) < 0.0d0) .or. (tcxpos2(i) > xcel)) then
        write(6,*) 'Error: tcxpos2(',i,') exceeds cell dimension'
        stop
     end if
     if (tcxpos1(i) >= tcxpos2(i)) then
!            write(6,*) 'Error: tcxpos1(',i,') exceeds tcxpos2(',i,')'
!            stop
        tcxpos1(i) = tcxpos1(i) - xcel
     end if

     if ((tcypos1(i) < 0.0d0) .or. (tcypos1(i) > ycel)) then
        write(6,*) 'Error: tcypos1(',i,') exceeds cell dimension'
        stop
     end if
     if ((tcypos2(i) < 0.0d0) .or. (tcypos2(i) > ycel)) then
        write(6,*) 'Error: tcypos2(',i,') exceeds cell dimension'
        stop
     end if
     if (tcypos1(i) >= tcypos2(i)) then
!            write(6,*) 'Error: tcypos1(',i,') exceeds tcypos2(',i,')'
!            stop
        tcypos1(i) = tcypos1(i) - ycel
     end if

     if ((tczpos1(i) < 0.0d0) .or. (tczpos1(i) > zcel)) then
        write(6,*) 'Error: tczpos1(',i,') exceeds cell dimension'
        stop
     end if
     if ((tczpos2(i) < 0.0d0) .or. (tczpos2(i) > zcel)) then
        write(6,*) 'Error: tczpos2(',i,') exceeds cell dimension'
        stop
     end if
     if (tczpos1(i) >= tczpos2(i)) then
!            write(6,*) 'Error: tczpos1(',i,') exceeds tczpos2(',i,')'
!            stop
        tczpos1(i) = tczpos1(i) - zcel
     end if

  end do

!---- non-dimensionalize parameters

  do i=1,ntcregion
     r_tcont(i) = r_tcont(i) / tempref
  end do

  close(iuini)

!     +     +     +     +     +     +     +

end subroutine rdtcont
