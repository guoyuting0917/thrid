!**************************************
!*  cnpvw.f90 Ver.1.3                 *
!*      for peachgk_md.f              *
!*            by N.Yamamoto           *
!*            modified by G.Kikugawa  *
!**************************************
! Time-stamp: <>

subroutine cnpvw(rcutrpvw,force,pot_rpvw,pot_vw)

  use md_global
  use mpi_global

  implicit none

!
!     Calculate interactions to control normal pressure
!       including normal force to VW
!                 and interaction between VW and physical system
!
!        Lupkowski, M., et al. J.Chem.Phys. 93(1), 1990
!
!
! ARGUMENTS:
!   INPUT
  real(8),intent(in):: rcutrpvw    ! RP-VW interaction cutoff length [non-d]

!   INPUT & OUTPUT

  real(8),intent(inout):: force(:,:)   ! force calculated here
                                       ! total force
                                       !  (including vdw, bond, angle etc)

!   OUTPUT

  real(8),intent(out):: pot_rpvw        ! potential of RP-VW interaction
  real(8),intent(out):: pot_vw          ! potential of constant force


! LOCAL:
  real(8):: frpvw(3,natom)              ! Morse force calculated here
!      real(8):: fmorall(3,natom) ! Morse force calculated here
!MPI                                       local(in MPI sense) field
!  real(8):: ene_co(maxnatmtyp,maxnatmtyp) ! coefficient of potential
  real(8):: potij                      ! temporal potential due to IJ
  real(8):: ftmp                       ! temporal force
  real(8):: fij                        ! force to be subtracted due to
                                       ! excluded atoms

  real(8):: zij                        ! rj(3)-ri(3)
  real(8):: zij_abs                    ! |Zij|
  real(8):: zij_abs_inv1               ! |Zij|**-1
  real(8):: alpha_inv                  ! alpha**-1
  real(8):: rtmp1
  real(8):: rtmp2
  real(8):: wel_alpha_inv1             ! (pi*rho*D)/alpha
  real(8):: wel_alpha_inv2             ! (pi*rho*D)/(2*alpha**2)
  real(8):: zij_abs_alpha_inv_2        ! Z+1/2alpha
  real(8):: zij_abs_alpha_inv          ! Z+1/alpha
  real(8):: zij_abs_2_alpha_inv        ! Z+2/alpha

!  real(8):: box(3)                     ! BOX size
!  real(8):: box_inv(3)                 ! inverse of BOX size

!  real(8):: sqcut                      ! rcutmor * rcutmor

!
! MPI  following variable are declared in MPI version
!
!      real*8:: pot_mor_lcl      !local Morse potential
!      real*8:: local,global     !tmps for global summation
! MPI end

  integer:: itype, jtype               ! VDW atom types
  integer:: i,j                        ! do loop index
  integer:: j1, j2                     ! limits of do loop
  integer:: jj                         ! atom index of the interacting atom

  integer:: i_rpvw

  real(8):: pi  ! = 3.1415...

!      real*8:: xij,xij0
!      real*8:: yij,yij0
!      real*8:: zij,zij0
!      real*8:: r_ij

!      integer:: im,ip,jm,jp,km,kp
!      real*8:: xm,xp,ym,yp,zm,zp

!      real*8:: fxct,fyct,fzct,pott

!           real(kind=dp), external::  ferfc, ferf
!          * whether erfc and erf are external or intrinsic depends
!            on systems. so, ferfc and ferf are defined in
!            MACHINE/(system)/functions.f90.
!            by Bestsystems.

!     +     +     +     +     +     +     +


! --- some preparations

  frpvw(1:3,1:natom) = 0.0d0

  pi = dacos(-1.0d0)

!  box(1) = xcel
!  box(2) = ycel
!  box(3) = zcel
!  box_inv(1:3) = 1.0d0/box(1:3)

!  sqcut = rcutrpvw * rcutrpvw

! --- CALCULATE NONBONDED FORCE BETWEEN VW AND PHYSICAL SYSTEM ---
!!!!!   Note that this routine can be used for Morse type interaction only.

  pot_rpvw = 0.0d0

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO mor1 = 1, NMORSE - 1

  looplast = nrpvw - 1

  DO i_rpvw = loopinit, looplast, loopstep

     j1 = rpvw_indexall(i_rpvw)              ! the first and
     j2 = rpvw_indexall(i_rpvw+loopstep) - 1 ! the last atoms interacting
                                          ! with i

     i = nrpvwlist(i_rpvw)
     itype = atmindex(i)                  ! atmtype of atom(i)

     IF (rpvw_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = rpvw_listall(j)             ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

        zij     = atmcor(3,i) - atmcor(3,jj)

        zij_abs = abs( zij )           ! |zij|

!       --- periodic boundary ---

!        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

!        zij_2      = zij(3)**2          ! |zij|^2

!       --- cut off ---
        if (zij_abs > rcutrpvw) cycle

!        rij_abs       = dsqrt(rij_abs2) ! |rij|

!            rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2
        zij_abs_inv1  = 1.0d0 / zij_abs ! 1/|zij|
        alpha_inv     = 1.0d0 / para_alpharpvw(itype,jtype)

        rtmp1 = exp(-1.0d0*(zij_abs - para_radrpvw(itype,jtype))   &
             &     * para_alpharpvw(itype,jtype)) ! exp(-(|rij|-r0)*alpha)

        rtmp2 = rtmp1 * rtmp1           ! exp(-2(|rij|-r0)*alpha)

        wel_alpha_inv1 =  pi * para_rhorpvw(itype,jtype) &
             &          * para_welrpvw(itype,jtype) * alpha_inv
!        wel_alpha_inv1 =  ene_co(itype,jtype) * para_welrpvw(itype,jtype) &
!             &          * alpha_inv
                     ! = (pi*rho*D)/alpha
        wel_alpha_inv2 = 0.5d0 * wel_alpha_inv1 * alpha_inv
                     ! = (pi*rho*D)/(2*alpha**2)

        zij_abs_alpha_inv_2 = zij_abs + alpha_inv *0.5d0
                          ! = Z+1/(2*alpha)
        zij_abs_alpha_inv   = zij_abs + alpha_inv
                          ! = Z+1/alpha
        zij_abs_2_alpha_inv = zij_abs_alpha_inv + alpha_inv
                          ! = Z+2/alpha

        ftmp       = wel_alpha_inv1 &
             &     * ( zij_abs_alpha_inv_2 * rtmp2 &
             &       - zij_abs_alpha_inv   * rtmp1 * 4.0d0 )
                 ! = (pi*rho*D)/alpha
                 ! * ((Z+1/(2*alpha))*rtmp2
                 !   -(Z+1/alpha)*4*rtmp1)

        potij      = wel_alpha_inv2 &
             &     * ( zij_abs_alpha_inv   * rtmp2 &
             &       - zij_abs_2_alpha_inv * rtmp1 * 8.0d0 )
                 ! = (pi*rho*D)/(2*alpha**2)
                 ! * ((Z+1/alpha)*rtmp2
                 !   -(Z+2/alpha)*8*rtmp1)

        pot_rpvw   = pot_rpvw + potij

        fij        = ftmp * zij * zij_abs_inv1
                 ! = (pi*rho*D)/alpha
                 ! * ((Z+1/(2*alpha))*rtmp2
                 !   -(Z+1/alpha)*4*rtmp1) * zij / |zij|


        frpvw(3,i)  = frpvw(3,i)  + fij
        frpvw(3,jj) = frpvw(3,jj) - fij

     END DO
  END DO

  force(1:3,1:natom) = force(1:3,1:natom) + frpvw(1:3,1:natom)

! --- CALCULATE CONSTANT FORCE EXERTING VIRTUAL WALL ---

  pot_vw = 0.0d0

  if( irank == 0 ) then

     i = nvwlist(1)
     j = nvwlist(2)

     if( (i/=0) .and. (j/=0) ) then

#if defined(_CNPVW_DEBUG)
        write(19,'(2E16.8)')force(3,i),force(3,j)
#endif
        if( atmcor(3,i) < atmcor(3,j) ) then

           force(3,i) = force(3,i) + alpha_wall
                                   ! add constant force "+alpha_wall"
           force(3,j) = force(3,j) - alpha_wall
                                   ! add constant force "-alpha_wall"

           pot_vw = - alpha_wall * ( atmcor(3,i) - atmcor(3,j) )
              ! U = - alpha_wall * zwi
        else

           force(3,i) = force(3,i) - alpha_wall
           force(3,j) = force(3,j) + alpha_wall

           pot_vw =   alpha_wall * ( atmcor(3,i) - atmcor(3,j) )

        end if

     end if

  end if

!     +     +     +     +     +     +     +

end subroutine cnpvw
