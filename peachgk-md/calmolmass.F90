!**************************************
!*  calmolmass.f90 Ver.1.0 '11.01.20  *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine calmolmass()

  use md_global

  implicit none

!
!     Calculate mass of each molecule
!
! ARGUMENTS:
!     OUTPUT

! LOCAL:
  integer:: i
  integer:: im,i1,i2,ii

!     +     +     +     +     +     +     +

!     --- some preparation ---

!     --- loop over molecule(im) ---

  DO im = 1, nmoleptindex-1
     i1 = molept_index(im)
     i2 = molept_index(im+1) - 1

     molmass(im) = 0.0d0

     DO ii = i1, i2         ! loop over atom(i)
        i = molept_list(ii)

        molmass(im) = molmass(im) + atmmass(i)

     END DO

  END DO

!     +     +     +     +     +     +     +

end subroutine calmolmass
