!*****************************************
!*  rattle_zg.f90 Ver.1.0 '12.07.19      *
!*      for peachgk_md.f                 *
!*            by G.Kikugawa              *
!*   (originally programmed by J. Kato)  *
!*              customized by T.Nakano   *
!*****************************************
subroutine rattle_c_zg(dt_short_cal,dt_long_cal, &
     &                 istep_short,nstep_short)

!    subroutine to constraint z coordinate of COM of selected molecules
!      by RATTLE (Andersen, H. C., 1983, J. Comput. Phys. 52, 24-34.)
!
!      'rattle_C' is for Coordinate resetting
!       for the first half of RATTLE.
!
  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: dt_short_cal    ! time step of short force
  real(8),intent(in):: dt_long_cal     ! time step of long force [non-d]

  integer,intent(in):: istep_short
  integer,intent(in):: nstep_short     ! number of step for short force

! LOCAL:
  integer:: iatom, imole          
  integer:: i, ii, i1, i2

  real(8):: corcent(3)            ! center of mass of molecule

  real(8):: molemass
  real(8):: delta_cor(3)

#if defined(_DO_NOT_USE_THIS)
  real(8):: gi               ! constraint force constant
#endif

  real(8):: dt_cal           ! time step

!     +     +     +     +     +     +     +     +     +     +

! ---- initialization ----

! --- Set timestep ---

#if defined(_DO_NOT_USE_THIS)
  if (istep_short == 1) then
     dt_cal = dt_short_cal
  else
     dt_cal = dt_short_cal
  end if
#endif

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

end subroutine rattle_c_zg

!----------------------------------------------------
subroutine rattlep_v_zg(dt_long_cal,dt_short_cal, &
     &                  istep_short,nstep_short)

!    subroutine to constraint z coordinate of COM of selected molecles
!      by RATTLE (Andersen, H. C., 1983, J. Comput. Phys. 52, 24-34.)
!
!      'rattle_V' is for Velocity resetting
!       for the secondx half of RATTLE.
!
  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]
  real(8),intent(in):: dt_short_cal     ! time step of short force

  integer,intent(in):: istep_short     ! do loop index for short range forces
  integer,intent(in):: nstep_short     ! number of step for short force

! LOCAL:
#if defined(_DO_NOT_USE_THIS)
  real(8):: amassi_inv       ! 1/atmmass(i)

  real(8):: ki               ! constraint force etc.

  real(8):: one = 1.0d0  

  real(8):: vtimestep        ! set timestep for virial calculation
  real(8):: inv_vtime        ! = 1/vtimestep

  real(8):: delta_vel(3)     !
#endif

  integer:: i, ii, i1, i2
  integer:: iatom, imole          ! atom & molecule index

  real(8):: molemass
  real(8):: velcent(3)         ! vel of centroid 

!     +     +     +     +     +     +     +     +     +     +

! --- GRAND LOOP TO APPLY VELOCITY CONSTRAINT ---

  DO i = 1, nlfixzg

     imole = index_nlfixzg(i)      ! molecular index

     i1 = molept_index(imole)
     i2 = molept_index(imole+1) - 1

     molemass = 0.0d0
     velcent(3) = 0.0d0

     do iatom = i1, i2         ! loop over atom(ii)
        ii = molept_list(iatom)

        molemass = molemass + atmmass(ii)

        velcent(3) = velcent(3) + atmmass(ii) * atmvel(3,ii)

     end do

     velcent(3) = velcent(3) / molemass

     ! - Velocity update -

#if defined(_DO_NOT_USE_THIS)
     ! strict RATTLE procedure
     amassi_inv = one/molemass(imole) ! 1/mi
     ki  = molemass(imole)*velcent(3)

     delta_vel(3) = amass_inv * ki - 0.0d0
#else
     do iatom = i1, i2         ! loop over atoms of a fix molecule
        ii = molept_list(iatom)

!        delta_vel(3) = velcent(3)
!        atmvel(3,ii) = atmvel(3,ii) - delta_vel(3)

        atmvel(3,ii) = atmvel(3,ii) - velcent(3)
     end do
#endif

  END DO

!     +     +     +     +     +     +     +     +     +     +

end subroutine rattlep_v_zg
