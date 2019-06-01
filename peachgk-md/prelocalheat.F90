!**************************************
!*  prelocalheat.f Ver.1.1 '10.06.30  *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine prelocalheat( npolytyp, nmatyp, &
     &                   polytyp_free, watertyp_free, matomtyp_free, &
     &                   tempref, &
     &                   nlheat_poly, index_nlheat_poly, &
     &                   tcont_nlheat_poly, &
     &                   nlheat_water, tcont_nlheat_water, &
     &                   nlheat_ma, index_nlheat_ma, &
     &                   tcont_nlheat_ma )

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.

  character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
  character(80),intent(in):: watertyp_free(:) ! use for water type control
  character(80),intent(in):: matomtyp_free(:,:)
                                ! use for matom type control

  real(8),intent(in):: tempref          ! temperature base value [K]

!     OUTPUT
  integer,intent(out):: nlheat_poly     ! number of poly type for local heating
  integer,intent(out):: index_nlheat_poly(:) 
                                ! index of poly type for local heating
  real(8),intent(out):: tcont_nlheat_poly(:)
                                ! control temp. of poly type for local heating
  integer,intent(out):: nlheat_water    ! number of water for local heating
  real(8),intent(out):: tcont_nlheat_water
                                ! control temp. of water for local heating
  integer,intent(out):: nlheat_ma       ! number of matom type for local heating
  integer,intent(out):: index_nlheat_ma(:) 
                                ! index of matom type for local heating
  real(8),intent(out):: tcont_nlheat_ma(:)
                                ! control temp. of matom type for local heating
 
! LOCAL:
  integer:: ipoly
  integer:: imatom
  integer:: iword

!     +     +     +     +     +     +

!---- make list for local heating of poly type

  nlheat_poly = 0
  do ipoly = 1, npolytyp
     do iword = 1, maxnword

        if (polytyp_free(ipoly,iword) .eq. 'localheat') then
           nlheat_poly = nlheat_poly + 1
           index_nlheat_poly(nlheat_poly) = ipoly
           read(polytyp_free(ipoly,iword+1),*) &
                & tcont_nlheat_poly(nlheat_poly)
           tcont_nlheat_poly(nlheat_poly) = tcont_nlheat_poly(nlheat_poly) &
                &                         / tempref
           exit
        end if

     end do
  end do

!---- make list for local heating of water type

  nlheat_water = 0
  do iword = 1, maxnword

     if (watertyp_free(iword) == 'localheat') then
        nlheat_water = nlheat_water + 1
        read(watertyp_free(iword+1),*) tcont_nlheat_water
        tcont_nlheat_water = tcont_nlheat_water / tempref
        exit
     end if

  end do

!---- make list for local heating of matom type

  nlheat_ma = 0
  do imatom = 1, nmatyp
     do iword = 1, maxnword

        if (matomtyp_free(imatom,iword) == 'localheat') then
           nlheat_ma = nlheat_ma + 1
           index_nlheat_ma(nlheat_ma) = imatom
           read(matomtyp_free(imatom,iword+1),*) tcont_nlheat_ma(nlheat_ma)
           tcont_nlheat_ma(nlheat_ma) = tcont_nlheat_ma(nlheat_ma) &
                &                     / tempref
           exit
        end if

     end do
  end do

!     +     +     +     +     +     +     +

  return
end subroutine prelocalheat
