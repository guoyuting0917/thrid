!*****************************
!*  rdlangevin.f90 Ver.1.0   *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-02-28 16:31:15 gota>

subroutine rdlangevin(xref, tempref, timeref, &
     &                xcel, ycel, zcel, &
     &                nlangeregion, &
     &                ltxpos1, ltxpos2, &
     &                ltypos1, ltypos2, &
     &                ltzpos1, ltzpos2, &
     &                r_ltemp, r_ltdamp)

  use interface_tools

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: tempref          ! temperature base value [K]
  real(8),intent(in):: timeref          ! time base value [sec]

  real(8),intent(in):: xcel             ! initial x0 cell length[m]
  real(8),intent(in):: ycel             ! initial y0 cell length[m]
  real(8),intent(in):: zcel             ! initial z0 cell length[m]

!     OUTPUT
  integer,intent(out):: nlangeregion    ! number of region for Langevin thermo.
  real(8),intent(out):: ltxpos1(:),ltxpos2(:)
                                        ! x-position of temp. control region
  real(8),intent(out):: ltypos1(:),ltypos2(:)
                                        ! y-position of temp. control region
  real(8),intent(out):: ltzpos1(:),ltzpos2(:)
                                        ! z-position of temp. control region
  real(8),intent(out):: r_ltemp(:)      ! control temp. in each region
  real(8),intent(out):: r_ltdamp(:)     
                                 ! damping factor in each region [1/s -> non-d]

! LOCAL:
  character(80):: fredat(20)

  integer:: nword

  integer:: i               ! do loop index

  integer:: iuini
  integer:: ios

  integer:: iltregion

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- intialize variable --------

!---- read parameter from MD script (peachgk.ini)

  iuini = 99
  open(iuini,file='langecont.ini',status='old',iostat=ios)
  if (ios /= 0) then
     write(6,*) 'Failure in opening file: langecont.ini'
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
     else if (fredat(1) == 'nlangeregion') then
        read(fredat(2),*) nlangeregion

     else if (fredat(1) == 'ltregion') then
        read(fredat(2),*) iltregion
        if (iltregion > nlangeregion) then
           write(6,*) 'Error: ltregion in rdlangevin'
           stop
        end if

        if (fredat(3) == 'xcel') then
           ltxpos1(iltregion) = 0.0d0
        else
           read(fredat(3),*) ltxpos1(iltregion)
           ltxpos1(iltregion) = ltxpos1(iltregion) / xref
        end if

        if (fredat(4) == 'xcel') then
           ltxpos2(iltregion) = xcel
        else
           read(fredat(4),*) ltxpos2(iltregion)
           ltxpos2(iltregion) = ltxpos2(iltregion) / xref
        end if

        if (fredat(5) == 'ycel') then
           ltypos1(iltregion) = 0.0d0
        else
           read(fredat(5),*) ltypos1(iltregion)
           ltypos1(iltregion) = ltypos1(iltregion) / xref
        end if

        if (fredat(6) == 'ycel') then
           ltypos2(iltregion) = ycel
        else
           read(fredat(6),*) ltypos2(iltregion)
           ltypos2(iltregion) = ltypos2(iltregion) / xref
        end if

        if (fredat(7) == 'zcel') then
           ltzpos1(iltregion) = 0.0d0
        else
           read(fredat(7),*) ltzpos1(iltregion)
           ltzpos1(iltregion) = ltzpos1(iltregion) / xref
        end if

        if (fredat(8) == 'zcel') then
           ltzpos2(iltregion) = zcel
        else
           read(fredat(8),*) ltzpos2(iltregion)
           ltzpos2(iltregion) = ltzpos2(iltregion) / xref
        end if

        read(fredat(9),*) r_ltemp(iltregion)

        read(fredat(10),*) r_ltdamp(iltregion)

     end if

  END DO

!---- check the parameter relations and so on

  do i=1,nlangeregion
     if ((ltxpos1(i) < 0.0d0) .or. (ltxpos1(i) > xcel)) then
        write(6,*) 'Error: ltxpos1(',i,') exceeds cell dimension'
        stop
     end if
     if ((ltxpos2(i) < 0.0d0) .or. (ltxpos2(i) > xcel)) then
        write(6,*) 'Error: ltxpos2(',i,') exceeds cell dimension'
        stop
     end if
     if (ltxpos1(i) >= ltxpos2(i)) then
!            write(6,*) 'Error: tcxpos1(',i,') exceeds ltxpos2(',i,')'
!            stop
        ltxpos1(i) = ltxpos1(i) - xcel
     end if

     if ((ltypos1(i) < 0.0d0) .or. (ltypos1(i) > ycel)) then
        write(6,*) 'Error: ltypos1(',i,') exceeds cell dimension'
        stop
     end if
     if ((ltypos2(i) < 0.0d0) .or. (ltypos2(i) > ycel)) then
        write(6,*) 'Error: ltypos2(',i,') exceeds cell dimension'
        stop
     end if
     if (ltypos1(i) >= ltypos2(i)) then
!            write(6,*) 'Error: ltypos1(',i,') exceeds ltypos2(',i,')'
!            stop
        ltypos1(i) = ltypos1(i) - ycel
     end if

     if ((ltzpos1(i) < 0.0d0) .or. (ltzpos1(i) > zcel)) then
        write(6,*) 'Error: ltzpos1(',i,') exceeds cell dimension'
        stop
     end if
     if ((ltzpos2(i) < 0.0d0) .or. (ltzpos2(i) > zcel)) then
        write(6,*) 'Error: ltzpos2(',i,') exceeds cell dimension'
        stop
     end if
     if (ltzpos1(i) >= ltzpos2(i)) then
!            write(6,*) 'Error: ltzpos1(',i,') exceeds ltzpos2(',i,')'
!            stop
        ltzpos1(i) = ltzpos1(i) - zcel
     end if

  end do

!---- non-dimensionalize parameters

  r_ltemp(1:nlangeregion) = r_ltemp(1:nlangeregion) / tempref
  r_ltdamp(1:nlangeregion) = r_ltdamp(1:nlangeregion) * timeref

  close(iuini)

!     +     +     +     +     +     +     +

end subroutine rdlangevin
