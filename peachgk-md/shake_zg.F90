!*****************************************
!*  shake_zg.f90 Ver.1.0 '12.07.19       *
!*      for peachgk_md.f                 *
!*            by G.Kikugawa              *
!*****************************************
subroutine shake_c_zg()

!    subroutine to constraint z coordinate of COM of selected molecules
!      by SHAKE (Rickaert, J. et al., 1977, J. Comput. Phys. 23, 327-341.)
!
  use md_global

  implicit none
      
! ARGUMENT:

! LOCAL:
  integer:: iatom, imole          
  integer:: i, ii, i1, i2

  real(8):: corcent(3)            ! center of mass of molecule

  real(8):: molemass
  real(8):: delta_cor(3)

#if defined(_DO_NOT_USE_THIS)
  real(8):: gi               ! constraint force constant
#endif

!     +     +     +     +     +     +     +     +     +     +

! --- Set timestep ---

  DO i = 1, nlfixzg

     imole = index_nlfixzg(i)      ! molecular index

     i1 = molept_index(imole)
     i2 = molept_index(imole+1) - 1

     molemass = 0.0d0
     corcent(3) = 0.0d0

!    - Coordinate update -

     do iatom = i1, i2         ! loop over atom(ii)
        ii = molept_list(iatom)

        molemass = molemass + atmmass(ii)

        corcent(3) = corcent(3) + atmmass(ii) * atmcor(3,ii)
     end do

     corcent(3) = corcent(3) / molemass

#if defined(_DO_NOT_USE_THIS)
     ! strict RATTLE procedure (not used)
     gi = molemass * (corcent(3) - cor_lfixzg(i))

     delta_cor(3) = gi/molemass

     do iatom = i1, i2         ! loop over atoms of a fix molecule
        ii = molept_list(iatom)
        atmcor(3,ii) = atmcor(3,ii) - delta_cor(3)
     end do
#else
     delta_cor(3) = corcent(3) - cor_lfixzg(i)

     do iatom = i1, i2         ! loop over atoms of a fix molecule
        ii = molept_list(iatom)
        atmcor(3,ii) = atmcor(3,ii) - delta_cor(3)
     end do
#endif

  END DO

!     +     +     +     +     +     +     +     +     +     +

end subroutine shake_c_zg
