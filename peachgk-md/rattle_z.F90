!*****************************************
!*  rattle_z.f90 Ver.1.1 '11.01.31       *
!*      for peachgk_md.f                 *
!*            by G.Kikugawa              *
!*   (originally programmed by J. Kato)  *
!*****************************************
subroutine rattle_c_z(dt_short_cal, dt_long_cal,   &
     &                istep_short, nstep_short)

!    subroutine to constraint z coordinate of selected atoms
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
  integer:: iatom          
  integer:: i

#if defined(_DO_NOT_USE_THIS)      
  real(8):: gi               ! constraint force constant
#endif
            
  real(8):: dt_cal           ! time step
            
!     +     +     +     +     +     +     +     +     +     +
      
! --- Set timestep ---

#if defined(_DO_NOT_USE_THIS)
  if (istep_short == 1) then
     dt_cal = dt_short_cal
  else 
     dt_cal = dt_short_cal
  end if
#endif
         
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

  return      
end subroutine rattle_c_z

!----------------------------------------------------
subroutine rattlep_v_z(dt_long_cal, dt_short_cal,   &
     &                 istep_short, nstep_short)

!    subroutine to constraint z coordinate of selected atoms
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
#endif
      
  integer:: i                ! do loop indexes
  integer:: iatom
      
!     +     +     +     +     +     +     +     +     +     +
            
! --- GRAND LOOP TO APPLY VELOCITY CONSTRAINT ---
         
  DO i = 1, nlfixz
                     
     iatom = index_nlfixz(i)

     ! - Velocity update -

#if defined(_DO_NOT_USE_THIS)
     ! strict RATTLE procedure
     amassi_inv = one/atmmass(iatom) ! 1/mi
     ki  = atmmass(iatom)*atmvel(3,iatom)
            
     atmvel(3,iatom) = atmvel(3,iatom) - amassi_inv * ki
#else
     atmvel(3,iatom) = 0.0d0
#endif
                     
  END DO
      
!     +     +     +     +     +     +     +     +     +     +
      
  return
end subroutine rattlep_v_z
