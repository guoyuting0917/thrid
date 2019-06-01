!***********************************
!*  calcom.f Ver.1.1 '10.06.24     *
!*      for peachgk_md.f           *
!*            by G.Kikugawa        *
!***********************************
subroutine calcom( npoly,nwater,nmatom, &
     &             xcel,ycel,zcel, &
     &             molecom)

  use md_global

  implicit none

!
!     Summing pressure term of each atom:
!           PV = 1/3 * [ PK + PU]
!        molecular pressure is written by:
!           PmolV = 1/3 * [ sigma(mol){Mmol*vg.vg} +
!                           sigma(mol){rg.Fmol} ]
!
! ARGUMENTS:
!     INPUT
  integer,intent(in):: npoly           ! number of polymer1
  integer,intent(in):: nwater          ! number of H2O molecules
  integer,intent(in):: nmatom          ! number of monatomic molecules

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

!     OUTPUT
  real(8),intent(out):: molecom(:,:)     ! center of mass of molecule

! LOCAL:

  real(8):: molemass         ! mass of molecule
  real(8):: inv_molemass     ! = 1/molemass

  real(8):: box(3)
  real(8):: box_inv(3)

  integer:: i
  integer:: im
  integer:: ii,i1,i2

!     +     +     +     +     +     +     +

!     ---- initialization ----

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

!---- calculate center of mass
!     - loop over poly

  DO im = 1, npoly
     i1 = molept_index(im)
     i2 = molept_index(im+1) - 1

     molemass = 0.0d0
     molecom(1:3,im) = 0.0d0

     do ii = i1, i2         ! loop over atom(i)
        i = molept_list(ii)

        molemass = molemass + atmmass(i)

        molecom(1:3,im) = molecom(1:3,im) &
             &          + atmmass(i) * atmcor(1:3,i)

     end do

     inv_molemass = 1.0d0 / molemass

     molecom(1:3,im) = molecom(1:3,im) * inv_molemass

  END DO

!     - loop over water

  DO im = npoly+1, npoly+nwater
     i1 = molept_index(im)
     i2 = molept_index(im+1) - 1

     molemass = 0.0d0
     molecom(1:3,im) = 0.0d0

     do ii = i1, i2         ! loop over atom(i)
        i = molept_list(ii)

        molemass = molemass + atmmass(i)

        molecom(1:3,im) = molecom(1:3,im) &
             &          + atmmass(i) * atmcor(1:3,i)

     end do

     inv_molemass = 1.0d0 / molemass

     molecom(1:3,im) = molecom(1:3,im) * inv_molemass

  END DO

!     - loop over monatom

  DO im = npoly+nwater+1, npoly+nwater+nmatom
     i1 = molept_index(im)
!         i2 = molept_index(im+1) - 1

!         do ii = i1, i2         ! loop over atom(i)
     ii = i1
     i = molept_list(ii)
     molecom(1:3,im) = atmcor(1:3,i)

!         end do

  END DO

!---- P.B.C.
      
  do im = 1, nmole
     if (molecom(1,im) < 0.0d0) then
        molecom(1,im) = molecom(1,im) + xcel
     else if (molecom(1,im) >= xcel) then
        molecom(1,im) = molecom(1,im) - xcel
     end if

     if (molecom(2,im) < 0.0d0) then
        molecom(2,im) = molecom(2,im) + ycel
     else if (molecom(2,im) >= ycel) then
        molecom(2,im) = molecom(2,im) - ycel
     end if

     if (molecom(3,im) < 0.0d0) then
        molecom(3,im) = molecom(3,im) + zcel
     else if (molecom(3,im) >= zcel) then
        molecom(3,im) = molecom(3,im) - zcel
     end if

  end do

!     +     +     +     +     +     +     +

  return
end subroutine calcom
