!*****************************
!*  caltorsim.f90 Ver.1.5    *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>

subroutine caltorsim(xcel,ycel,zcel,   &
     &               force,pot_torsim)

  use interface_interact, only: calijklim

  use md_global
  use mpi_global

  implicit none

!     subroutine to calculate improper torsion force & potential

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel           ! x cell length[non-d]
  real(8),intent(in):: ycel           ! y cell length[non-d]
  real(8),intent(in):: zcel           ! z cell length[non-d]

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! force calculated here

!   OUTPUT
  real(8),intent(out):: pot_torsim      ! improper torsion potential

! LOCAL:
  real(8):: delta_pot                   ! potential calculated
                                        !  in calijkl

  real(8):: rij(3)                      ! Rj-Ri
  real(8):: rkj(3)                      ! Rj-Rk
  real(8):: rkl(3)                      ! Rl-Rk
  real(8):: rik(3)                      ! Rk-Ri
  real(8):: rjl(3)                      ! Rl-Rj

  real(8):: fi(3), fj(3), fk(3), fl(3)  ! force to atoms (i,j,k,l)
                                        !  calculated here

  real(8):: phi                         ! torsion (in radian) defined by
                                        !  atoms i, j, k, l

  real(8):: pi                          ! pi

  integer:: n                           ! do loop index

  integer:: iatom, jatom, katom, latom  ! atom index

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

! --- GRAND LOOP TO CALCULATE TORSION FORCE ---


  pi = dacos(-1.0d0)
  pot_torsim = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

  looplast = ntorsim

! MPI loop
!      DO n = 1, ntorsim
  DO n = loopinit, looplast, loopstep

!    - calculate Rij, Rkj, Rkl, etc. -

     iatom   = itorsim(n)               ! atom index
     jatom   = jtorsim(n)               ! atom index
     katom   = ktorsim(n)               ! atom index
     latom   = ltorsim(n)               ! atom index

     rij(1:3)  = atmcor(1:3,jatom) - atmcor(1:3,iatom)
     rkj(1:3)  = atmcor(1:3,jatom) - atmcor(1:3,katom)
     rkl(1:3)  = atmcor(1:3,latom) - atmcor(1:3,katom)
     rik(1:3)  = atmcor(1:3,katom) - atmcor(1:3,iatom)
     rjl(1:3)  = atmcor(1:3,latom) - atmcor(1:3,jatom)

     rij(1:3)  = rij(1:3) - box(1:3)*anint(rij(1:3)*box_inv(1:3))
     rkj(1:3)  = rkj(1:3) - box(1:3)*anint(rkj(1:3)*box_inv(1:3))
     rkl(1:3)  = rkl(1:3) - box(1:3)*anint(rkl(1:3)*box_inv(1:3))
     rik(1:3)  = rik(1:3) - box(1:3)*anint(rik(1:3)*box_inv(1:3))
     rjl(1:3)  = rjl(1:3) - box(1:3)*anint(rjl(1:3)*box_inv(1:3))

!    - calculate torsion force and angle--

     call calijklim(rij,rkj,rkl,rik,rjl,   &
          &         para_barhigim(indextorsimtyp(n)),   &
          &         para_phsangim(indextorsimtyp(n)),   &
          &         delta_pot,fi,fj,fk,fl,phi)

!    - add to force & potential -

     force(1:3,iatom) = force(1:3,iatom) + fi(1:3)
     force(1:3,jatom) = force(1:3,jatom) + fj(1:3)
     force(1:3,katom) = force(1:3,katom) + fk(1:3)
     force(1:3,latom) = force(1:3,latom) + fl(1:3)

     pot_torsim = pot_torsim + delta_pot

!     phi = phi * 180.0d0 / pi
!     if (phi < 0.0d0) phi = phi + 360.0d0


  END DO

!     +     +     +     +     +     +     +

end subroutine caltorsim

! -----------------------------------------------------------------

subroutine caltorsimp(xcel,ycel,zcel,   &
     &                force,pot_torsim,   &
     &                atm_viri_tors,atm_virit_tors)

  use interface_interact, only: calijklim

  use md_global
  use mpi_global

  implicit none

!     subroutine to calculate improper torsion force & potential

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel           ! x cell length[non-d]
  real(8),intent(in):: ycel           ! y cell length[non-d]
  real(8),intent(in):: zcel           ! z cell length[non-d]

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! force calculated here

!   OUTPUT
  real(8),intent(out):: pot_torsim      ! improper torsion potential

  real(8),intent(out):: atm_viri_tors   ! virial(torsion potential) of each atom
  real(8),intent(out):: atm_virit_tors(:,:)  ! virial tensor (torsion potential)

! LOCAL:
  real(8):: delta_pot                   ! potential calculated
                                        !  in calijkl

  real(8):: rij(3)                      ! Rj-Ri
  real(8):: rkj(3)                      ! Rj-Rk
  real(8):: rkl(3)                      ! Rl-Rk
  real(8):: rik(3)                      ! Rk-Ri
  real(8):: rjl(3)                      ! Rl-Rj

  real(8):: fi(3), fj(3), fk(3), fl(3)  ! force to atoms (i,j,k,l)
                                        !  calculated here

  real(8):: phi                         ! torsion (in radian) defined by
                                        !  atoms i, j, k, l

  real(8):: pi                          ! pi

  integer:: n,l                         ! do loop index

  integer:: iatom, jatom, katom, latom  ! atom index

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

! --- GRAND LOOP TO CALCULATE TORSION FORCE ---


  pi = acos(-1.0d0)
  pot_torsim = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

  looplast = ntorsim

! MPI loop
!      DO n = 1, ntorsim
  DO n = loopinit, looplast, loopstep

!    - calculate Rij, Rkj, Rkl, etc. -

     iatom   = itorsim(n)               ! atom index
     jatom   = jtorsim(n)               ! atom index
     katom   = ktorsim(n)               ! atom index
     latom   = ltorsim(n)               ! atom index

     rij(1:3)  = atmcor(1:3,jatom) - atmcor(1:3,iatom)
     rkj(1:3)  = atmcor(1:3,jatom) - atmcor(1:3,katom)
     rkl(1:3)  = atmcor(1:3,latom) - atmcor(1:3,katom)
     rik(1:3)  = atmcor(1:3,katom) - atmcor(1:3,iatom)
     rjl(1:3)  = atmcor(1:3,latom) - atmcor(1:3,jatom)

     rij(1:3)  = rij(1:3) - box(1:3)*anint(rij(1:3)*box_inv(1:3))
     rkj(1:3)  = rkj(1:3) - box(1:3)*anint(rkj(1:3)*box_inv(1:3))
     rkl(1:3)  = rkl(1:3) - box(1:3)*anint(rkl(1:3)*box_inv(1:3))
     rik(1:3)  = rik(1:3) - box(1:3)*anint(rik(1:3)*box_inv(1:3))
     rjl(1:3)  = rjl(1:3) - box(1:3)*anint(rjl(1:3)*box_inv(1:3))

!    - calculate torsion force and angle--

     call calijklim(rij,rkj,rkl,rik,rjl,   &
          &         para_barhigim(indextorsimtyp(n)),   &
          &         para_phsangim(indextorsimtyp(n)),   &
          &         delta_pot,fi,fj,fk,fl,phi)

!    - add to force & potential & virial -

     force(1:3,iatom) = force(1:3,iatom) + fi(1:3)
     force(1:3,jatom) = force(1:3,jatom) + fj(1:3)
     force(1:3,katom) = force(1:3,katom) + fk(1:3)
     force(1:3,latom) = force(1:3,latom) + fl(1:3)

     atm_viri_tors = atm_viri_tors + fi(1)*atmcor(1,iatom)   &
          &        + fi(2)*atmcor(2,iatom) + fi(3)*atmcor(3,iatom)
     atm_viri_tors = atm_viri_tors + fj(1)*atmcor(1,jatom)   &
          &        + fj(2)*atmcor(2,jatom) + fj(3)*atmcor(3,jatom)
     atm_viri_tors = atm_viri_tors + fk(1)*atmcor(1,katom)   &
          &        + fk(2)*atmcor(2,katom) + fk(3)*atmcor(3,katom)
     atm_viri_tors = atm_viri_tors + fl(1)*atmcor(1,latom)   &
          &        + fl(2)*atmcor(2,latom) + fl(3)*atmcor(3,latom)

     do l=1,3
        atm_virit_tors(1:3,l) = atm_virit_tors(1:3,l)   &
             &                + fi(1:3)*atmcor(l,iatom)
        atm_virit_tors(1:3,l) = atm_virit_tors(1:3,l)   &
             &                + fj(1:3)*atmcor(l,jatom)
        atm_virit_tors(1:3,l) = atm_virit_tors(1:3,l)   &
             &                + fk(1:3)*atmcor(l,katom)
        atm_virit_tors(1:3,l) = atm_virit_tors(1:3,l)   &
             &                + fl(1:3)*atmcor(l,latom)
     end do

     pot_torsim = pot_torsim + delta_pot

!     phi = phi * 180.0d0 / pi
!     if (phi < 0.0d0) phi = phi + 360.0d0


  END DO

!     +     +     +     +     +     +     +

end subroutine caltorsimp

!------------------------------------------------------------------
subroutine calijklim(rij, rkj, rkl, rik, rjl,   &
     &               barhig, phsang,   &
     &               delta_pot, fi,fj,fk,fl, phi)

  implicit none

!     subroutine to calculate improper torsion force & potential
!     Harmonic and umbrella (used in DREIDING) type improper torsion

! ARGUMENTS:
!   INPUT
!   -- vectors ---

  real(8),intent(in)::  rij(3)          ! Rj-Ri
  real(8),intent(in)::  rkj(3)          ! Rj-Rk
  real(8),intent(in)::  rkl(3)          ! Rl-Rk
  real(8),intent(in)::  rik(3)          ! Rk-Ri
  real(8),intent(in)::  rjl(3)          ! Rl-Rj

! -- force field parameters ---

  real(8),intent(in):: barhig           ! barrier height (barhig)
  real(8),intent(in):: phsang           ! phase angle    (phsang)

!   OUTPUT

  real(8),intent(out):: delta_pot       ! potential calculated
                                        !  in this subroutine

  real(8),intent(out):: fi(3), fj(3), fk(3), fl(3) ! force to atoms (i,j,k,l)
                                                   !  calculated here
  real(8),intent(out):: phi             ! torsion (in radian) defined by
                                        !  atoms i, j, k, l
                                        !  -pi =< phi < pi

! LOCAL:

  real(8):: pi                          ! 3.1415

  real(8):: rij_abs                     ! |Rij|
  real(8):: rkj_abs                     ! |Rkj|
  real(8):: rkl_abs                     ! |Rkl|

  real(8):: norm_ijk(3)                 ! normal vector of plane (i,j,k)
  real(8):: norm_jkl(3)                 ! normal vector of plane (j,k,l)
  real(8):: norm_ijk_abs                ! absolute value of norm_ijk
  real(8):: norm_jkl_abs                ! absolute value of norm_jkl

  real(8):: cos_phi                     ! cos(phi)
!  real(8):: sin_phi                     ! sin(phi)
  real(8):: cos_phsang                  ! cos(phsang)
  real(8):: barhig_tmp                  ! C_umb

  real(8):: phi_tmp

  real(8):: f_tmp_a                     ! force for temporal usage
  real(8):: f_tmp_b(3)                  !
  real(8):: f_tmp_c(3)                  !

  real(8):: r_tmp1(3)                   ! vectors for temporal usage
  real(8):: r_tmp2(3)                   !

! FUNCTIONS:

  real(8),external:: dot_product

  real(8),external:: dsinc

!  -------------------------------------------------
!      IMPROPER TORSION FORCE & ENERGY
!
!          force due to (i,j,k,l) improper torsion:
!
!          and potential energy:
!          <Harmonic type>
!             U = V * (phi - phi_0)^2
!          <Umbrella type>
!             U = V * (cos(phi) - cos(phi_0))
!             where V = -K (when phi_0 = 0),
!                     = 1/2 K * (1+cos(phi_0))^2/sin(phi_0)^2
!  -------------------------------------------------


! --- SOME INITIAL EVALUATION ---

  pi = acos(-1.0d0)


! --- CALCULATION OF SIZE OF RIJ ETC ---

  rij_abs = sqrt(dot_product(rij,rij))
  rkj_abs = sqrt(dot_product(rkj,rkj))
  rkl_abs = sqrt(dot_product(rkl,rkl))


! --- CALCULATION OF NORMAL VECTORS ---

  call vec_product(rij, rkj, norm_ijk)
  call vec_product(rkj, rkl, norm_jkl)

  norm_ijk_abs = sqrt(dot_product(norm_ijk,norm_ijk))
  norm_jkl_abs = sqrt(dot_product(norm_jkl,norm_jkl))


! --- CALCULATION OF PHI AND RELATED VARIABLES ---

! - calculate cos_phi --

  cos_phi = dot_product(norm_ijk,norm_jkl)   &
       &  / (norm_ijk_abs*norm_jkl_abs)

! - adjust if cos_phi shows a value such as 1.000000001 -

  if (cos_phi > 1.0d0) then
     cos_phi =  1.0d0
  else if (cos_phi < -1.0d0) then
     cos_phi = -1.0d0
  end if

! - calculate phi(ijkl) -

  call vec_product(norm_ijk, norm_jkl, r_tmp1)
                               ! r_tmp1(:) and rkj(:) determines the sign of phi


! --- CALCULATION OF POTENTIAL ---

#if defined(_IMPROP_UMBR)
! Umbrella type
  cos_phsang = cos(phsang)

  if (abs(phsang) < 1.0d-8) then
     barhig_tmp = -barhig   ! C_umb = -K

  else
     barhig_tmp = 0.5d0 * barhig * (1.0d0 + cos_phsang)**2 / (sin(phsang))**2
                               ! C_umb = 1/2 K * (1+cos(phi_0))^2/sin(phi_0)^2

  endif

  delta_pot = barhig_tmp * (cos_phi - cos_phsang)

#else
! Harmonic type
  phi = acos(cos_phi)                   !  0 <phi<pi
!!!!! phi must be positive !!!!!
!      phi = -dsign(phi,dot_product(rkj,r_tmp1)) ! -pi<phi<pi
!
!      sin_phi = dsin(phi)

  phi_tmp = phi - phsang

  delta_pot = barhig * phi_tmp * phi_tmp

#endif

! --- CALCULATION OF FORCE ---

! - set so that |norm| = 1 -

  norm_ijk(1:3) = norm_ijk(1:3)/norm_ijk_abs
  norm_jkl(1:3) = norm_jkl(1:3)/norm_jkl_abs

! - calculate f_tmp_a -

!      * by definition f_tmp_a is evaluated using the formula below:
!           f_tmp_a  = 2.0 * barhig *
!                      * (phi - phi_0) /sin(phi)
!         However, it causes problems when sin(phi) ~ 0.

#if defined(_IMPROP_UMBR)
! f_a for umbrella type
  f_tmp_a = -barhig_tmp   ! f_a = -C_umb

#else
! f_a for harmonic type

  if (abs(phsang) < 1.0d-8) then

     f_tmp_a = 2.0d0 * barhig / dsinc(phi)

  else if (abs(phsang - pi) < 1.0d-8) then

     f_tmp_a = - 2.0d0 * barhig / dsinc(phi_tmp)

  else

     phi = phi + 1.0d-16
     f_tmp_a = 2.0d0 * barhig * phi_tmp / dsin(phi)

  end if

#endif

! - calculate f_tmp_b, f_tmp_c -

  f_tmp_b(1:3) = (norm_jkl(1:3) - cos_phi*norm_ijk(1:3))/norm_ijk_abs
  f_tmp_c(1:3) = (norm_ijk(1:3) - cos_phi*norm_jkl(1:3))/norm_jkl_abs

! - calculate Fi, Fj, Fk, Fl -

  call vec_product(f_tmp_b, rkj, r_tmp1)
  fi(1:3) = f_tmp_a * r_tmp1(1:3)

  call vec_product(f_tmp_c, rkl, r_tmp1)
  call vec_product(f_tmp_b, rik, r_tmp2)
  fj(1:3) = f_tmp_a * ( - r_tmp1(1:3) + r_tmp2(1:3) )

  call vec_product(f_tmp_b, rij, r_tmp1)
  call vec_product(f_tmp_c, rjl, r_tmp2)
  fk(1:3) = f_tmp_a * ( - r_tmp1(1:3) + r_tmp2(1:3) )

  call vec_product(f_tmp_c, rkj, r_tmp1)
  fl(1:3) = f_tmp_a * r_tmp1(1:3)

!     +     +     +     +     +     +     +

end subroutine calijklim

!-----------------------------------------------------------------
function dsinc(x)

  implicit none

!     this function is to calculate the sinc function, namely, sin(x)/x

  real(8):: dsinc                       ! = sin(x)/x
  real(8),intent(in):: x                ! input

  real(8):: taylor_0_bound = 1.0d-32
  real(8):: taylor_2_bound
  real(8):: taylor_n_bound

  real(8):: x2                          ! = x*x

  taylor_2_bound = sqrt(taylor_0_bound)
  taylor_n_bound = sqrt(taylor_2_bound)

  if (abs(x) >= taylor_n_bound) then

     dsinc = sin(x)/x
     return

  else
!---- approximation by taylor series in x at 0 up tp order 0
     dsinc = 1.0d0

     if (dabs(x) >= taylor_0_bound) then
        x2 = x*x
!----   approximation by taylor series in x at 0 up tp order 2
        dsinc = dsinc - x2 / 6.0d0

        if (abs(x) >= taylor_2_bound) then
!----      approximation by taylor series in x at 0 up tp order 4
           dsinc = dsinc + x2*x2 / 120.0d0
        end if

     end if

  end if

end function dsinc
