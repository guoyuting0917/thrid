!*****************************
!*  calnon14_hf.f90 Ver.1.6  *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>

subroutine calnon14p_hf(xcel,ycel,zcel,   &
     &                  div_factor_14vdw,div_factor_14elc,   &
     &                  force,pot_elc14,pot_vdw14,   &
     &                  atm_viri_14,atm_virit_14,   &
     &                  pot_14_atm,viri14t_atm,   &
     &                  ifhfvol,   &
     &                  nhfregion,hfzpos1,hfzpos2,   &
     &                  hftyp_atm,   &
     &                  molecom)

  use interface_interact, only: calvirihf

  use md_global
  use mpi_global

  implicit none

!    subroutine to calculate 14 nonbonded force & potential
!    (namely 1-4 electorstatic & vdw interactions)
!
!    CHARMM type 1-4 interaction (multiple dihedral) is treated here
!    when -D switch (-D_CHARMM_NONB14) is activated.
!    Partly modified by Y. Kawagoe, Tohoku University.

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: div_factor_14vdw ! division factor of 14vdw
  real(8),intent(in):: div_factor_14elc ! division factor of 14elc

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: ifhfvol       ! local volume-based or local surface-based

  integer,intent(in):: nhfregion       ! number of region to calculate heat flux
  real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                        ! z-position of region for heat flux

  integer,intent(in):: hftyp_atm(:)     ! atom- or mole-based heat flux cal.
                                        !   for each atom

  real(8),intent(in):: molecom(:,:)     ! center of mass of molecule

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! atomic force

!   OUTPUT
  real(8),intent(out):: atm_viri_14   ! virial(1-4 force potential) of each atom
  real(8),intent(out):: atm_virit_14(:,:) ! virial tensor (1-4 force potential)

  real(8),intent(out):: pot_elc14       ! 14 electrostatic energy
  real(8),intent(out):: pot_vdw14       ! 14 vdw  energy

  real(8),intent(out):: pot_14_atm(:)   ! 1-4 potential of each atom

  real(8),intent(out):: viri14t_atm(:,:,:,:)
                                        ! virial tensor of each atom (1-4)

! LOCAL:

  real(8):: rij(3)                      ! rj-ri
  real(8):: rij2                        ! |rij|**2
  real(8):: rij_abs_inv1                ! |rij|**-1
  real(8):: rij_abs_inv2                ! 1/|rij|**2
  real(8):: rtmp2                       ! |Rij/rij|**2
  real(8):: rtmp6                       ! |Rij/rij|**6
  real(8):: rtmp12                      ! |Rij/rij|**12

  real(8):: felc14(3)                   ! 14 electrostatic force calculated here
  real(8):: fvdw14(3)                   ! 14 vdw force calculated here
  real(8):: felc41(3)                   ! = -felc14
  real(8):: fvdw41(3)                   ! = -fvdw41
  real(8):: delta_pot_elc               ! 14 elc potential calculated here
  real(8):: delta_pot_vdw               ! 14 vdw potential calculated here
  real(8):: Eij                         ! epsilon
  real(8):: div_factor_vdw              ! 1/vdw_scale_14
  real(8):: div_factor_elc              ! 1/elc_scale_14
#if defined(_CHARMM_NONB14)
  real(8):: vdw_welij_charmm            ! epsilon for LJ14 of CHARMM F.F.
  real(8):: vdw_radij_charmm            ! radius for LJ14 of CHARMM F.F.
#endif

  real(8):: ftmp                        ! temporal force

  integer:: iatom, jatom                ! atom index
  integer:: itype, jtype                ! atom type index

  integer:: i,n                         ! do loop index

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

  real(8):: nbodycoeff             ! n-body coefficient for virial term of heatf

! ---- initialization

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

  nbodycoeff = 0.5d0        ! = 1/2

! --- GRAND LOOP TO CALCULATE 14 ELC & VDW INTERACTIONS ---

  pot_elc14 = 0.0d0
  pot_vdw14 = 0.0d0
  felc14(1:3) = 0.0d0
  fvdw14(1:3) = 0.0d0

!---- for periodic torsion

  looplast = ntors

! MPI loop
!      DO i = 1, ntors
  DO i = loopinit, looplast, loopstep

     if (.not. ifcalnonb14(indextorstyp(i))) cycle

!    - set iatom & jatom -

     iatom = itors(i)
     jatom = ltors(i)
     if (iatom == jatom) cycle   ! three-membered ring like epoxy group

     itype = atmindex(iatom) ! atom type of iatom
     jtype = atmindex(jatom) ! atom type of jatom

!    - set the dividing factor -
    ! correction factor for 1-4 vdW and elc for ring structure
    ! 0.0 = within 4- or 5- membered ring
    ! 0.5 = within 6-membered ring
    ! 1.0 = others
#if defined(_CHARMM_NONB14)
     vdw_welij_charmm = para_divfac_vdw(indextorstyp(i))   &
     &                * para_wfac(indextorstyp(i)) !!! by KAWAGOE
     vdw_radij_charmm = para_divfac_elc(indextorstyp(i))
     div_factor_elc = para_wfac(indextorstyp(i)) !!! by KAWAGOE
#else
     div_factor_vdw = 1.0d0 / para_divfac_vdw(indextorstyp(i)) * ring14(i)
     div_factor_elc = 1.0d0 / para_divfac_elc(indextorstyp(i)) * ring14(i)
#endif

!    - calculate Rij etc. -

     rij(1:3)      = atmcor(1:3,iatom) - atmcor(1:3,jatom)
     rij(1:3)      = rij(1:3) - box(1:3)*anint(rij(1:3)*box_inv(1:3))

     rij2          = rij(1)**2 + rij(2)**2 + rij(3)**2
     rij_abs_inv2  = 1.0d0 / rij2
     rij_abs_inv1  = sqrt(rij_abs_inv2)

!---- calculate 14 electrostatic

     delta_pot_elc = atmchrg(iatom) * atmchrg(jatom) * rij_abs_inv1
     delta_pot_elc = delta_pot_elc * div_factor_elc

     felc14(1:3)   = delta_pot_elc * rij_abs_inv2 * rij(1:3)

!    - add to force & potential & virial -

     pot_elc14 = pot_elc14 + delta_pot_elc

     pot_14_atm(iatom) = pot_14_atm(iatom)   &
          &            + nbodycoeff*delta_pot_elc
     pot_14_atm(jatom) = pot_14_atm(jatom)   &
          &            + nbodycoeff*delta_pot_elc

     force(1:3,iatom) = force(1:3,iatom) + felc14(1:3)
     force(1:3,jatom) = force(1:3,jatom) - felc14(1:3)

!    - for heat flux calculation
     felc41(1:3) = -felc14(1:3)

     call calvirihf(iatom,jatom,felc14,felc41,   &
          &         nbodycoeff,   &
          &         box,box_inv,   &
          &         ifhfvol,   &
          &         nhfregion,hfzpos1,hfzpos2,   &
          &         hftyp_atm,   &
          &         molecom,   &
          &         viri14t_atm)

!    - for pressure calculation
     atm_viri_14 = atm_viri_14 + felc14(1)*rij(1)   &
          &      + felc14(2)*rij(2) + felc14(3)*rij(3)

     do n=1,3
        atm_virit_14(1:3,n) = atm_virit_14(1:3,n)   &
             &              + felc14(1:3)*rij(n)
     end do

!---- calculate 14 VDW Interaction

#if defined(_CHARMM_NONB14)
     ! 1-4 LJ of CHARMM F.F.
     rtmp2         = rij_abs_inv2 * vdw_radij_charmm**2
     rtmp6         = rtmp2 * rtmp2 * rtmp2
     rtmp12        = rtmp6 * rtmp6
     Eij           = vdw_welij_charmm
#else
     ! normal 1-4 LJ
     rtmp2         = rij_abs_inv2 * vdw_radij(itype,jtype)**2
     rtmp6         = rtmp2 * rtmp2 * rtmp2
     rtmp12        = rtmp6 * rtmp6
     Eij           = vdw_welij(itype,jtype) * div_factor_vdw
#endif
     delta_pot_vdw = Eij * (rtmp12 - 2.0d0 * rtmp6)
     ftmp          = 12.0d0 * Eij * (rtmp12 - rtmp6)   &
          &        * rij_abs_inv2

     fvdw14(1:3) = ftmp * rij(1:3)

!    - add to force & potential & virial -

     pot_vdw14 = pot_vdw14 + delta_pot_vdw

     pot_14_atm(iatom) = pot_14_atm(iatom)   &
          &            + nbodycoeff*delta_pot_vdw
     pot_14_atm(jatom) = pot_14_atm(jatom)   &
          &            + nbodycoeff*delta_pot_vdw

     force(1:3,iatom) = force(1:3,iatom) + fvdw14(1:3)
     force(1:3,jatom) = force(1:3,jatom) - fvdw14(1:3)

!    - for heat flux calculation
     fvdw41(1:3) = -fvdw14(1:3)

     call calvirihf(iatom,jatom,fvdw14,fvdw41,   &
          &         nbodycoeff,   &
          &         box,box_inv,   &
          &         ifhfvol,   &
          &         nhfregion,hfzpos1,hfzpos2,   &
          &         hftyp_atm,   &
          &         molecom,   &
          &         viri14t_atm)

!    - for pressure calculation
     atm_viri_14 = atm_viri_14 + fvdw14(1)*rij(1)   &
          &      + fvdw14(2)*rij(2) + fvdw14(3)*rij(3)

     do n=1,3
        atm_virit_14(1:3,n) = atm_virit_14(1:3,n)   &
             &              + fvdw14(1:3)*rij(n)
     end do

  END DO

!---- for RB torsion (OPLS FF)

  looplast = ntorsrb

! MPI loop
!      DO i = 1, ntorsrb
  DO i = loopinit, looplast, loopstep

     if (.not. ifcalnonb14_rb(indextorsrbtyp(i))) cycle

!    - set iatom & jatom -

     iatom = itorsrb(i)
     jatom = ltorsrb(i)
     if (iatom == jatom) cycle   ! three-membered ring like epoxy group

     itype = atmindex(iatom)            ! atom type of iatom
     jtype = atmindex(jatom)            ! atom type of jatom

!    - set the dividing factor -
    ! correction factor for 1-4 vdW and elc for ring structure
    ! 0.0 = within 4- or 5- membered ring
    ! 0.5 = within 6-membered ring
    ! 1.0 = others
     div_factor_vdw = 1.0d0 / para_divfac_vdwrb(indextorsrbtyp(i)) &
     &              * ring14_rb(i)
     div_factor_elc = 1.0d0 / para_divfac_elcrb(indextorsrbtyp(i)) &
     &              * ring14_rb(i)

!    - calculate Rij etc. -

     rij(1:3)      = atmcor(1:3,iatom) - atmcor(1:3,jatom)
     rij(1:3)      = rij(1:3) - box(1:3)*anint(rij(1:3)*box_inv(1:3))

     rij2          = rij(1)**2 + rij(2)**2 + rij(3)**2
     rij_abs_inv2  = 1.0d0 / rij2
     rij_abs_inv1  = sqrt(rij_abs_inv2)

!---- calculate 14 electrostatic

     delta_pot_elc = atmchrg(iatom) * atmchrg(jatom) * rij_abs_inv1
     delta_pot_elc = delta_pot_elc * div_factor_elc

     felc14(1:3)     = delta_pot_elc * rij_abs_inv2 * rij(1:3)

!    - add to force & potential & virial -

     pot_elc14 = pot_elc14 + delta_pot_elc

     pot_14_atm(iatom) = pot_14_atm(iatom)   &
          &            + nbodycoeff*delta_pot_elc
     pot_14_atm(jatom) = pot_14_atm(jatom)   &
          &            + nbodycoeff*delta_pot_elc

     force(1:3,iatom) = force(1:3,iatom) + felc14(1:3)
     force(1:3,jatom) = force(1:3,jatom) - felc14(1:3)

!    - for heat flux calculation
     felc41(1:3) = -felc14(1:3)

     call calvirihf(iatom,jatom,felc14,felc41,   &
          &         nbodycoeff,   &
          &         box,box_inv,   &
          &         ifhfvol,   &
          &         nhfregion,hfzpos1,hfzpos2,   &
          &         hftyp_atm,   &
          &         molecom,   &
          &         viri14t_atm)

!    - for pressure calculation
     atm_viri_14 = atm_viri_14 + felc14(1)*rij(1)   &
          &      + felc14(2)*rij(2) + felc14(3)*rij(3)

     do n=1,3
        atm_virit_14(1:3,n) = atm_virit_14(1:3,n)   &
             &              + felc14(1:3)*rij(n)
     end do

!---- calculate 14 VDW Interaction

     rtmp2         = rij_abs_inv2 * vdw_radij(itype,jtype)**2
     rtmp6         = rtmp2 * rtmp2 * rtmp2
     rtmp12        = rtmp6 * rtmp6
     Eij           = vdw_welij(itype,jtype) * div_factor_vdw

     delta_pot_vdw = Eij * (rtmp12 - 2.0d0 * rtmp6)
     ftmp          = 12.0d0 * Eij * (rtmp12 - rtmp6)   &
          &        * rij_abs_inv2

     fvdw14(1:3) = ftmp * rij(1:3)

!    - add to force & potential & virial -

     pot_vdw14 = pot_vdw14 + delta_pot_vdw

     pot_14_atm(iatom) = pot_14_atm(iatom)   &
          &            + nbodycoeff*delta_pot_vdw
     pot_14_atm(jatom) = pot_14_atm(jatom)   &
          &            + nbodycoeff*delta_pot_vdw

     force(1:3,iatom) = force(1:3,iatom) + fvdw14(1:3)
     force(1:3,jatom) = force(1:3,jatom) - fvdw14(1:3)

!    - for heat flux calculation
     fvdw41(1:3) = -fvdw14(1:3)

     call calvirihf(iatom,jatom,fvdw14,fvdw41,   &
          &         nbodycoeff,   &
          &         box,box_inv,   &
          &         ifhfvol,   &
          &         nhfregion,hfzpos1,hfzpos2,   &
          &         hftyp_atm,   &
          &         molecom,   &
          &         viri14t_atm)

!    - for pressure calculation
     atm_viri_14 = atm_viri_14 + fvdw14(1)*rij(1)   &
          &      + fvdw14(2)*rij(2) + fvdw14(3)*rij(3)

     do n=1,3
        atm_virit_14(1:3,n) = atm_virit_14(1:3,n)   &
             &              + fvdw14(1:3)*rij(n)
     end do

  END DO

!     +     +     +     +     +     +     +

end subroutine calnon14p_hf
