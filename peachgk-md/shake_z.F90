!*****************************************
!*  shake_z.f90 Ver.1.0 '12.05.23        *
!*      for peachgk_md.f                 *
!*            by G.Kikugawa              *
!*****************************************
subroutine shake_c_z()

!    subroutine to constraint z coordinate of selected atoms
!      by SHAKE (Rickaert, J. et al., 1977, J. Comput. Phys. 23, 327-341.)
!
  use md_global

  implicit none
      
! ARGUMENT:

! LOCAL:
  integer:: iatom          
  integer:: i

#if defined(_DO_NOT_USE_THIS)      
  real(8):: gi               ! constraint force constant
#endif

!     +     +     +     +     +     +     +     +     +     +

! --- Set timestep ---

  DO i = 1, nlfixz

     iatom = index_nlfixz(i)

!    - Coordinate update -

#if defined(_DO_NOT_USE_THIS)
     ! strict RATTLE procedure (not used)
     gi = atmmass(iatom) * (atmcor(3,iatom) - cor_lfixz(i))

     atmcor(3,iatom) = atmcor(3,iatom) - gi/atmmass(iatom)
#else
     atmcor(3,iatom) = cor_lfixz(i)

#endif

  END DO

!     +     +     +     +     +     +     +     +     +     +

end subroutine shake_c_z
