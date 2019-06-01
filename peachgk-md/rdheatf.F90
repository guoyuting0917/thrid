!*****************************
!*  rdheatf.f Ver.1.6        *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-05-12 17:04:44 gota>

subroutine rdheatf(ouhtfname, &
     &             oumtfname, &
     &             xref, &
     &             zcel, &
     &             npolytyp, &
     &             ifhfvol, &
     &             nhfregion, hfzpos1, hfzpos2, &
     &             hftyp_pmole, hftyp_water, &
     &             ifcalmtf, ifcalhtf, &
     &             mtfoutdir)

  use interface_tools

  use md_global

  implicit none
  
! ARGUMENT:
!     INPUT
  real(8),intent(in):: xref             ! distanse base value [m]

  real(8),intent(in):: zcel             ! initial z0 cell length[m]

  integer,intent(in):: npolytyp        ! number of poly type

!     OUTPUT
  character(80),intent(out):: ouhtfname ! output file name

  character(80),intent(out):: oumtfname ! output file name

  logical,intent(out):: ifhfvol       ! local volume-based or local surface-based

  integer,intent(out):: nhfregion       ! number of region to calculate heat flux
  real(8),intent(out):: hfzpos1(:),hfzpos2(:)
                          ! z-position of region for heat flux

  integer,intent(out):: hftyp_pmole(:)
                          ! atom- or mole-based heat flux cal. for poly
  integer,intent(out):: hftyp_water
                          ! atom- or mole-based heat flux cal. for water

  logical,intent(out):: ifcalmtf ! flag to calculate & output momentum flux
  logical,intent(out):: ifcalhtf ! flag to calculate & output heat flux

  character(3),intent(out):: mtfoutdir ! direction to output data of momentum

! LOCAL:
  character(80):: fredat(20)

  integer:: nword

  integer:: i               ! do loop index

  integer:: iuini
  integer:: ios

  integer:: ihfzpos

  integer:: ihftyp_pmole

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- error check for HF calculation

#if defined(_HF_ALL_DIR) && defined(_HF_BULK)
  write(6,*) 'Error: Do not use -D_HF_ALL_DIR and -D_HF_BULK switches'
  write(6,*) '       at once when compiling peachgk_md.'
  stop
#endif

!---- intialize variable --------

  do i = 1, npolytyp
     hftyp_pmole(i) = 0
  end do
  hftyp_water = 0

  ifhfvol = .true.

!---- read parameter from MD script (peachgk.ini)

  iuini = 99
  open(iuini,file='transflux.ini',status='old',iostat=ios)
  if (ios /= 0) then
     write(6,*) 'Failure in opening file: transflux.ini'

     open(iuini,file='heatflux.ini',status='old',iostat=ios)
     if (ios /= 0) then
        write(6,*) 'Failure in opening file: heatflux.ini'

        open(iuini,file='momentumflux.ini',status='old',iostat=ios)
        if (ios /= 0) then
           write(6,*) 'Failure in opening file: momentumflux.ini'
           write(6,*) 'Failure in opening any initial files'
           write(6,*) ' for transport flux calculation'
           stop
        end if

     end if

  end if

  DO
     call rdfree_w( iuini,20, fredat, nword )

     if (fredat(1)(1:3) == 'END') then
        exit

     else if (fredat(1) == ' ') then
        cycle

     else if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) then
        cycle                     ! comment line

!        output file name
     else if (fredat(1) == 'ouhtfname') then
        ouhtfname = fredat(2)

     else if (fredat(1) == 'oumtfname') then
        oumtfname = fredat(2)

!        local volume-based or surface-based method
     else if (fredat(1) == 'ifhfvol') then
        read(fredat(2),*) ifhfvol

!        definition of calculation region
     else if (fredat(1) == 'nhfregion') then
        read(fredat(2),*) nhfregion

     else if (fredat(1) == 'hfzpos') then
        read(fredat(2),*) ihfzpos
        if (ihfzpos > nhfregion) then
           write(6,*) 'Error: hfzpos in rdheatf'
           stop
        end if
        read(fredat(3),*) hfzpos1(ihfzpos)
        if (ifhfvol) then   ! volume-based only
           read(fredat(4),*) hfzpos2(ihfzpos)
        end if

!        declare atom-base or molecule-base calculation
     else if (fredat(1) == 'hftyp_pmole') then
        read(fredat(2),*) ihftyp_pmole
        if (ihftyp_pmole > npolytyp) then
           write(6,*) 'Error: hftyp_pmole in rdheatf'
           stop
        end if
        if (fredat(3) == 'ATOM') then
           hftyp_pmole(ihftyp_pmole) = HFTYP_ATOM
        else if (fredat(3) == 'MOLE') then
           hftyp_pmole(ihftyp_pmole) = HFTYP_MOLE
        end if

     else if (fredat(1) == 'hftyp_water') then
        if (fredat(2) == 'ATOM') then
           hftyp_water = HFTYP_ATOM
        else if (fredat(2) == 'MOLE') then
           hftyp_water = HFTYP_MOLE
        end if

!        flag to calculate & output momentum or heat flux
     else if (fredat(1) == 'ifcalmtf') then
        read(fredat(2),*) ifcalmtf
     else if (fredat(1) == 'ifcalhtf') then
        read(fredat(2),*) ifcalhtf

     else if (fredat(1) == 'mtfoutdir') then
        read(fredat(2),*) mtfoutdir

     end if

  END DO

  close(iuini)

!---- check the parameter relations and so on

#if defined(_HF_BULK)
! when using bulk mode, force to put zcel dimension to hfzpos1(1) and hfzpos1(2)
  hfzpos1(1) = 0.0d0
  hfzpos2(1) = zcel
  if (.not. ifhfvol) then      ! surface-based method is prohibited
     write(6,*) 'Error: do not use surface-based method in HF_BULK mode'
     stop
  end if

#else
  do i=1,nhfregion
     if ((hfzpos1(i) < 0.0d0) .or. (hfzpos1(i) >= zcel)) then
        write(6,*) 'Error: hfzpos1(',i,') exceeds cell dimension.'
        stop
     end if
     if ((hfzpos2(i) < 0.0d0) .or. (hfzpos2(i) >= zcel)) then
        write(6,*) 'Error: hfzpos2(',i,') exceeds cell dimension.'
        stop
     end if
     if (ifhfvol) then      ! only for volume-based method
        if (hfzpos1(i) > hfzpos2(i)) then
           write(6,*) 'Error: hfzpos1 must be smaller then hfzpos2.'
           stop
        end if
     end if
  end do
#endif

  do i = 1, npolytyp
     if (hftyp_pmole(i) == 0) then
        write(6,*) 'Error: hftyp_pmole is not set in hftyp_pmole(', i, ').'
        stop
     end if
  end do

  if ((mtfoutdir /= 'X') .and. (mtfoutdir /= 'Y') .and. (mtfoutdir /= 'Z')   &
       &                                       .and. (mtfoutdir /= 'ALL')) then
     write(6,*) 'Error: mtfoutdir must be "X", "Y", "Z", or "ALL"'
     stop
  end if

!---- non-dimensionalize parameters

  hfzpos1(1:nhfregion) = hfzpos1(1:nhfregion) / xref
  hfzpos2(1:nhfregion) = hfzpos2(1:nhfregion) / xref

!     +     +     +     +     +     +     +

end subroutine rdheatf
