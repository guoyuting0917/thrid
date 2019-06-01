!*****************************
!*  rdhfcont.f Ver.1.3       *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2014-06-13 16:47:08 gota>

subroutine rdhfcont(xref, eref, timeref, &
     &              xcel, ycel, zcel, &
     &              dt_long_cal, &
     &              tcontinterval, &
     &              nhfcregion, &
     &              hfczpos1,hfczpos2, &
     &              r_hfcont)

  use interface_tools

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: eref             ! energy base value [J]
  real(8),intent(in):: timeref          ! time base value [sec]

  real(8),intent(in):: xcel             ! initial x0 cell length[non-d]
  real(8),intent(in):: ycel             ! initial y0 cell length[non-d]
  real(8),intent(in):: zcel             ! initial z0 cell length[non-d]

  real(8),intent(in):: dt_long_cal      ! time step of long force [sec]
  integer,intent(in):: tcontinterval   ! interval of temp. control

!     OUTPUT
  integer,intent(out):: nhfcregion      ! number of region to control heat flux
  real(8),intent(out):: hfczpos1(:),hfczpos2(:)
                                ! z-position of heat flux control region
  real(8),intent(out):: r_hfcont(:)      ! magnitude of heat flux in each region
                                ! (converted to the input energy)

! LOCAL:
  character(80):: fredat(20)

  integer:: nword

  integer:: i               ! do loop index

  integer:: iuini
  integer:: ios

  integer:: ihfcregion

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- intialize variable --------

!---- read parameter from MD script (hfcont.ini)

  iuini = 99
  open(iuini,file='hfcont.ini',status='old',iostat=ios)
  if (ios /= 0) then
     write(6,*) 'Failure in opening file: hfcont.ini'
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

!        definition of h.f. control region
     else if (fredat(1) == 'nhfcregion') then
        read(fredat(2),*) nhfcregion

     else if (fredat(1) == 'hfcregion') then
        read(fredat(2),*) ihfcregion
        if (ihfcregion > nhfcregion) then
           write(6,*) 'Error: hfcregion in rdhfcont'
           stop
        end if

        read(fredat(3),*) hfczpos1(ihfcregion)
        hfczpos1(ihfcregion) = hfczpos1(ihfcregion) / xref

        read(fredat(4),*) hfczpos2(ihfcregion)
        hfczpos2(ihfcregion) = hfczpos2(ihfcregion) / xref

        read(fredat(5),*) r_hfcont(ihfcregion)
        r_hfcont(ihfcregion) = r_hfcont(ihfcregion) &
             &               / (eref/(xref*xref*timeref))

     end if

  END DO

!---- check the parameter relations and so on

  do i = 1, nhfcregion

     if ((hfczpos1(i) < 0.0d0) .or. (hfczpos1(i) > zcel)) then
        write(6,*) 'Error: hfczpos1(',i,') exceeds cell dimension'
        stop
     end if
     if ((hfczpos2(i) < 0.0d0) .or. (hfczpos2(i) > zcel)) then
        write(6,*) 'Error: hfczpos2(',i,') exceeds cell dimension'
        stop
     end if
     if (hfczpos1(i) >= hfczpos2(i)) then
!            write(6,*) 'Error: hfczpos1(',i,') exceeds hfczpos2(',i,')'
!            stop
        hfczpos1(i) = hfczpos1(i) - zcel
     end if

  end do

!---- convert heat flux to input (output) energy

  do i=1,nhfcregion
     r_hfcont(i) = r_hfcont(i) * 2.0d0 * xcel * ycel &
          &      * dt_long_cal * tcontinterval / timeref
                                ! dt_long_cal is not non-dimensionalized here
  end do

#if defined(_HFCONT_DEBUG)
  write(6,*) '***** h.f. control debug info *****'
  do i=1,nhfcregion
     write(6,*) 'region No. ',i,', input energy[J]: ', &
          &     r_hfcont(i) * eref
  end do
#endif

  close(iuini)

!     +     +     +     +     +     +     +

end subroutine rdhfcont
