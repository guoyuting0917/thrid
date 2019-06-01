!*****************************
!*  ewaldr_hf.f Ver.2.0      *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-06-02 12:59:28 gota>

subroutine ewaldrp_hf(xcel,ycel,zcel,   &
     &                alpha,rrcut,   &
     &                rcut,   &
     &                force,pot_elc,pot_vdw,   &
     &                for_viri_coul,pot_viri_coul,   &
     &                for_viri_lj,pot_viri_lj,   &
     &                pot_virit_coul,pot_virit_lj,   &
     &                pot_elin_atm,virielint_atm,   &
     &                pot_elc_atm,virielct_atm,    &
     &                pot_ljin_atm,viriljint_atm,   &
     &                pot_vdw_atm,viriljt_atm,   &
     &                ifhfvol,   &
     &                nhfregion,hfzpos1,hfzpos2,   &
     &                hftyp_atm,   &
     &                molecom)

  use interface_interact, only: calvirihf

  use md_global
  use mpi_global

  implicit none
      
!
!    
!     Calculate VDW and Ewald-REAl interactions & PRESSURE VIRIAL
!        
!        VDW:
!           U   = sigma(i>j) [ Eij {(Rij/rij)^12 - 2 (Rij/rij)^6} ]
!           Fix = sigma(j) [ 12Eij {(Rij/rij)^12 - (Rij/rij)^6} xij/rij^2 ]
!
!        Ewald-real: 
!           U   = sigma(i>j) [Qi Qj/rij erfc(rij/ehta)] 
!           Fix = sigma(j) [- Qi Qj/rij^2 
!              {erfc(rij/ehta)/rij + 2/ehta/sqrt(pi) exp(-rij^2/ehta^2}] xij 
!
!     Comment on Excluded Neighbors
!           The nonbonded list does not contain excluded neighbors.
!           Hence U and F calculated here do not contain data of 
!           excluded neighbors.  However, in Ewald summation,
!               Utot = Uewald(all pairs) - Ucoulomb(excluded pairs)
!           Hence, Uewr - Ucoulomb due to excluded pairs are 
!               added to the potential and force calculated here.
!
!
! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel          ! x cell length[non-d]
  real(8),intent(in):: ycel          ! y cell length[non-d]
  real(8),intent(in):: zcel          ! z cell length[non-d]

  real(8),intent(in):: alpha         ! parameter alpha [non-d]
  real(8),intent(in):: rrcut         ! ewald real space cutoff length [non-d]

  real(8),intent(in):: rcut          ! vdw cutoff length [non-d]

  logical,intent(in):: ifhfvol       ! local volume-based or local surface-based

  integer,intent(in):: nhfregion     ! number of region to calculate heat flux
  real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                     ! z-position of region for heat flux

  integer,intent(in):: hftyp_atm(:)  ! atom- or mole-based heat flux cal. 
                                     !   for each atom

  real(8),intent(in):: molecom(:,:)  ! center of mass of molecule

!   INPUT & OUTPUT

  real(8),intent(inout):: force(:,:) ! force calculated here
                                     ! total force
                                     !  (including vdw, bond, angle etc) 

  real(8),intent(inout):: for_viri_coul(:,:) ! virial(coulomb force) of each atom
  real(8),intent(inout):: pot_viri_coul  ! virial(coulomb potential) of each atom
  real(8),intent(inout):: pot_virit_coul(:,:) ! virial tensor (coulomb)

  real(8),intent(inout):: for_viri_lj(:,:)   ! virial(L-J force) of each atom
  real(8),intent(inout):: pot_viri_lj        ! virial(L-J potential) of each atom
  real(8),intent(inout):: pot_virit_lj(:,:)  ! virial tensor (L-J)

  real(8),intent(inout):: pot_elin_atm(:)  ! elc (intra) potential of each atom
  real(8),intent(inout):: pot_ljin_atm(:)  ! L-J (intra) potential of each atom
  real(8),intent(inout):: pot_elc_atm(:)   ! Coulomb potential of each atom
  real(8),intent(inout):: pot_vdw_atm(:)   ! VDW potential of each atom
  
  real(8),intent(inout):: virielint_atm(:,:,:,:)
                                     ! virial tensor of each atom (elc intra)
  real(8),intent(inout):: viriljint_atm(:,:,:,:)
                                     ! virial tensor of each atom (L-J intra)
  real(8),intent(inout):: virielct_atm(:,:,:,:)
                                     ! virial tensor of each atom (Coulomb)
  real(8),intent(inout):: viriljt_atm(:,:,:,:)
                                     ! virial tensor of each atom (L-J)

!   OUTPUT
  real(8),intent(out):: pot_elc      ! electrostatic potential
  real(8),intent(out):: pot_vdw      ! vdw potential

! LOCAL:
  real(8):: felc(3,maxnatom)         ! electrostatic force calculated here
  real(8):: fvdw(3,maxnatom)         ! VDW force calculated here

  real(8):: pi                       ! = 3.14...

  real(8):: const_fewr               ! = 2/ehta/sqrt(pi), used for Ewald 
                                     !   real force   
  real(8):: chrgij                   ! = qi * qj
  real(8):: potij                    ! temporal potential due to IJ  
  real(8):: ftmp                     ! temporal force  
  real(8):: fij(3)                   ! force to be subtracted due to 
                                     ! excluded atoms
  real(8):: fji(3)                   ! force to be subtracted due to 
                                     ! excluded atoms = -fij

  real(8):: Eij                      ! well depth 

  real(8):: rij(3)                   ! rj-ri
  real(8):: rij_abs                  ! |Rij|   
  real(8):: rij_abs2                 ! |Rij|**2
  real(8):: rij_abs_inv1             ! |Rij|**-1
  real(8):: rij_abs_inv2             ! 1/|rij|**2
  real(8):: rtmp1
  real(8):: rtmp2, rtmp6, rtmp12           
  real(8):: rtmp_erf, rtmp_erfc           

  real(8):: box(3)                   ! BOX size
  real(8):: box_inv(3)               ! inverse of BOX size

  real(8):: sqcut                    ! rrcut * rrcut
  real(8):: sqcut2                   ! rcut * rcut

  integer:: itype, jtype             ! VDW atom types
  integer:: i,j,n                    ! do loop index
  integer:: j1, j2                   ! limits of do loop
  integer:: jj                       ! atom index of the interacting atom
  
!  real(8):: rij_fij(3)               ! rij dyad fij (tensor product)

  real(8):: nbodycoeff             ! n-body coefficient for virial term of heatf

  integer:: im, jm                   ! index of molecule
!      integer:: im,i1,i2,ii,jm  ! loop index of molecule

!      logical:: exclatom        ! flag for if atom belongs list_excl

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

! FUNCTIONS:
  real(8):: derf                     ! error function
!  real(8),external:: derf             ! error function
  real(8):: SPLINE_INTERPFUNC        ! spline interpolation function

!     +     +     +     +     +     +     +

! --- ALLOCATE ARRAY FOR FORCE CALCULATED IN THIS SUBROUTINE --- 

  felc(1:3,1:natom) = 0.0d0
  fvdw(1:3,1:natom) = 0.0d0

  pi         = dacos(-1.0d0)
  const_fewr = 2.0d0 * alpha / dsqrt(pi) 
                                     ! coefficient used for calculation of
                                     ! Ewald real

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

  nbodycoeff = 0.5d0                 ! = 1/2

  sqcut = rrcut * rrcut
  sqcut2 = rcut * rcut
  if (sqcut2 > sqcut) sqcut = sqcut2

! --- CALCULATE NONBONDED FORCE ---


  pot_elc = 0.0d0
  pot_vdw = 0.0d0

! MPI loop
! DO i = 1, natom - 1

  looplast = natom - 1

  do i = loopinit, looplast, loopstep

     j1 = nb_indexall(i)              ! the first and
     j2 = nb_indexall(i+loopstep) - 1 ! the last atoms interacting
                                      ! with i

     itype = atmindex(i)              ! VDW type of atom(i)

     im = irmolept_list(i)            ! what molecule atom(i) belongs in

     IF (nb_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = nb_listall(j)            ! atom index of interacting atom

        jm = irmolept_list(jj)        ! what molecule atom(jj) belongs in

!       --- CALCULATE VIRIAL & FORCE ---

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

!       --- EWALD REAL ---

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2 ! |rij|^2

!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2

!!! if calculate L-J potential only
#if defined(_LJ_ONLY)
        ! no execution
#else
        rij_abs       = dsqrt(rij_abs2) ! |rij| 
        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

        if (rij_abs <= rrcut) then

!!! Spline interpolation for ewald real
#if defined(_NOT_SPLINTERP)
           rtmp1 = rij_abs * alpha   ! |rij|*alpha
           rtmp2 = rtmp1 * rtmp1     ! |rij|^2*alpha^2

           rtmp_erfc  = 1.0d0 - derf(rtmp1) ! erfc(|rij|*alpha)
           chrgij     = atmchrg(i) * atmchrg(jj)
           potij      = chrgij * rij_abs_inv1 

           ftmp       = potij * (rtmp_erfc  * rij_abs_inv2   &
                &     + const_fewr * rij_abs_inv1 * dexp(-rtmp2))
                            ! =   Qi * Qj * [erfc(|rij|*alpha)/|rij|^3 
                            ! + 2*alpha/sqrt(pi)/|rij|^2 exp(-|rij|^2*alpha^2)]

           potij      = potij * rtmp_erfc 
                                     ! qi qj erfc(|rij|/ehta) / |rij|

#else
           chrgij = atmchrg(i) * atmchrg(jj)
           ftmp = SPLINE_INTERPFUNC(rij_abs,spltbl_real,   &
              &                     spl_b_real,spl_c_real,spl_d_real)   &
              & * chrgij
           potij = SPLINE_INTERPFUNC(rij_abs,spltbl_realpot,   &
               &                     spl_b_realpot,spl_c_realpot,spl_d_realpot) &
               & * chrgij
#endif

!          --- Fennell method
           if (fennell_flag) then

              potij = potij + chrgij   &
                  & * (-fennell_shift   &
                  &   + fennell_damp*(rij_abs - rrcut))
              ftmp = ftmp - chrgij * fennell_damp * rij_abs_inv1

           end if

           pot_elc    = pot_elc + potij

           fij(1:3)   = ftmp * rij(1:3)

           felc(1:3,i)  = felc(1:3,i)  + fij(1:3)
           felc(1:3,jj) = felc(1:3,jj) - fij(1:3)

           for_viri_coul(1:3,i)  = for_viri_coul(1:3,i)  + fij(1:3)
           for_viri_coul(1:3,jj) = for_viri_coul(1:3,jj) - fij(1:3)

!          --- calculate virial tensor ---

!          - for pressure calculation
           do n=1,3
              pot_virit_coul(1:3,n) = pot_virit_coul(1:3,n)   &
              &                     + fij(1:3)*rij(n)
           end do

           pot_viri_coul =  pot_viri_coul + potij

!          - for heat flux calculation

!          --- SKIP transport flux calculation for ewaldr
#if defined(_HF_SKIP_EWALDR)
           ! no execution
#else

!          --- Bulk mode with Ewald calculation (no execution)
#if defined(_HF_BULK_EWK)
           ! no execution

!          --- Fennell method
#elif defined(_HF_FENNELL)
           if (fennell_flag) then
                                  ! no execution (fennell potential calculated)
           else
              potij = potij + chrgij   &
                  & * (-fennell_shift   &
                  &   + fennell_damp*(rij_abs - rrcut))
              ftmp = ftmp - chrgij * fennell_damp * rij_abs_inv1

              fij(1:3) = ftmp * rij(1:3)
           end if

!          --- pure 1/r potential (default)
#else
           potij      = atmchrg(i) * atmchrg(jj) * rij_abs_inv1

           ftmp = potij * rij_abs_inv2

           fij(1:3) = ftmp * rij(1:3)
#endif

           fji(1:3) = -fij(1:3)

           if (im == jm) then         ! intramolecular interaction
              pot_elin_atm(i) = pot_elin_atm(i) + 0.5d0 * potij
              pot_elin_atm(jj) = pot_elin_atm(jj) + 0.5d0 * potij

              call calvirihf(i,jj,fij,fji,   &
                   &         nbodycoeff,   &
                   &         box,box_inv,   &
                   &         ifhfvol,   &
                   &         nhfregion,hfzpos1,hfzpos2,   &
                   &         hftyp_atm,   &
                   &         molecom,   &
                   &         virielint_atm)
           else                       ! intermolecular interaction
              pot_elc_atm(i) = pot_elc_atm(i) + 0.5d0 * potij
              pot_elc_atm(jj) = pot_elc_atm(jj) + 0.5d0 * potij

              call calvirihf(i,jj,fij,fji,   &
                   &         nbodycoeff,   &
                   &         box,box_inv,   &
                   &         ifhfvol,   &
                   &         nhfregion,hfzpos1,hfzpos2,   &
                   &         hftyp_atm,   &
                   &         molecom,   &
                   &         virielct_atm)
           end if

#endif
!!! end if SKIP transport flux calculation for ewaldr

        end if

#endif
!!! end if calculate L-J potential only

!       --- VDW ---

        if (rij_abs2 <= sqcut2) then

           jtype      = atmindex(jj) ! VDW type of atom(jj)

           if (inter_inttyp(itype,jtype) /= INTTYPE_VDW) cycle

           rtmp2      = rij_abs_inv2 * vdw_radij(itype,jtype)**2
           rtmp6      = rtmp2 * rtmp2 * rtmp2
           rtmp12     = rtmp6 * rtmp6
           Eij        = vdw_welij(itype,jtype)
           potij      = Eij * (rtmp12 - 2.0d0 * rtmp6)
           pot_vdw    = pot_vdw + potij

  if ( i <= 1792 .and. jtype == 110 ) then                                                      !!!!!!!!! wall-liquid potential energy

     wall_pot_interface_1_ar1 = wall_pot_interface_1_ar1 + potij

    else if ( i <= 1792 .and. jtype == 111 ) then  

     wall_pot_interface_1_ar2 = wall_pot_interface_1_ar2 + potij



    else if ( i >1792 .and. i <= 3584 .and. jtype == 110 ) then                                        

     wall_pot_interface_2_ar1 = wall_pot_interface_2_ar1 + potij

    else if ( i >1792 .and. i <= 3584 .and. jtype == 111 ) then                           

     wall_pot_interface_2_ar2 = wall_pot_interface_2_ar2 + potij

    else

     wall_pot_interface_no =  wall_pot_interface_no + potij

    end if                                                                                           !!!!!!!!! wall-liquid potential energy


           ftmp       = 12.0d0 * Eij * (rtmp12 - rtmp6)   &
                &     * rij_abs_inv2

           fij(1:3)     = ftmp * rij(1:3)

           fvdw(1:3,i)  = fvdw(1:3,i)  + fij(1:3)
           fvdw(1:3,jj) = fvdw(1:3,jj) - fij(1:3)

           for_viri_lj(1:3,i)  = for_viri_lj(1:3,i)  + fij(1:3)
           for_viri_lj(1:3,jj) = for_viri_lj(1:3,jj) - fij(1:3)

!          --- calculate virial tensor ---

!          - for heat flux calculation

           fji(1:3) = -fij(1:3)

           if (im == jm) then        ! intramolecular interaction
              pot_ljin_atm(i) = pot_ljin_atm(i) + 0.5d0*potij
              pot_ljin_atm(jj) = pot_ljin_atm(jj) + 0.5d0*potij

              call calvirihf(i,jj,fij,fji,   &
                   &         nbodycoeff,   &
                   &         box,box_inv,   &
                   &         ifhfvol,   &
                   &         nhfregion,hfzpos1,hfzpos2,   &
                   &         hftyp_atm,   &
                   &         molecom,   &
                   &         viriljint_atm)
           else                      ! intermolecular interaction
              pot_vdw_atm(i) = pot_vdw_atm(i) + 0.5d0*potij
              pot_vdw_atm(jj) = pot_vdw_atm(jj) + 0.5d0*potij

              call calvirihf(i,jj,fij,fji,   &
                   &         nbodycoeff,   &
                   &         box,box_inv,   &
                   &         ifhfvol,   &
                   &         nhfregion,hfzpos1,hfzpos2,   &
                   &         hftyp_atm,   &
                   &         molecom,   &
                   &         viriljt_atm)

  if ( itype ==110 .and. jtype == 111 ) then                                        !!!!!!!!!calculate ar1-ar2 pot LJ energy

    pot_vdw_atm_ar12(i) = pot_vdw_atm_ar12(i) + 0.5d0*potij
    pot_vdw_atm_ar12(jj) = pot_vdw_atm_ar12(jj) + 0.5d0*potij

            call calvirihf(i,jj,fij,fji,   &
                   &         nbodycoeff,   &
                   &         box,box_inv,   &
                   &         ifhfvol,   &
                   &         nhfregion,hfzpos1,hfzpos2,   &
                   &         hftyp_atm,   &
                   &         molecom,   &
                   &         viriljt_atm_ar12)

    else if ( itype ==110 .and. jtype == 110 ) then  

    pot_vdw_atm_ar11(i) = pot_vdw_atm_ar11(i) + 0.5d0*potij
    pot_vdw_atm_ar11(jj) = pot_vdw_atm_ar11(jj) + 0.5d0*potij

            call calvirihf(i,jj,fij,fji,   &
                   &         nbodycoeff,   &
                   &         box,box_inv,   &
                   &         ifhfvol,   &
                   &         nhfregion,hfzpos1,hfzpos2,   &
                   &         hftyp_atm,   &
                   &         molecom,   &
                   &         viriljt_atm_ar11)

    else if ( itype ==111 .and. jtype == 111 ) then  

    pot_vdw_atm_ar22(i) = pot_vdw_atm_ar22(i) +  0.5d0*potij
    pot_vdw_atm_ar22(jj) = pot_vdw_atm_ar22(jj) +  0.5d0*potij

             call calvirihf(i,jj,fij,fji,   &
                   &         nbodycoeff,   &
                   &         box,box_inv,   &
                   &         ifhfvol,   &
                   &         nhfregion,hfzpos1,hfzpos2,   &
                   &         hftyp_atm,   &
                   &         molecom,   &
                   &         viriljt_atm_ar22)


   else if ( itype ==16 .and. jtype == 110 ) then  

             call calvirihf(i,jj,fij,fji,   &
                   &         nbodycoeff,   &
                   &         box,box_inv,   &
                   &         ifhfvol,   &
                   &         nhfregion,hfzpos1,hfzpos2,   &
                   &         hftyp_atm,   &
                   &         molecom,   &
                   &         viriljt_atm_ar1pt)

   else if ( itype ==16 .and. jtype == 111 ) then  

             call calvirihf(i,jj,fij,fji,   &
                   &         nbodycoeff,   &
                   &         box,box_inv,   &
                   &         ifhfvol,   &
                   &         nhfregion,hfzpos1,hfzpos2,   &
                   &         hftyp_atm,   &
                   &         molecom,   &
                   &         viriljt_atm_ar2pt)


    end if                                                                             !!!!!!!!!calculate ar1-ar2 pot LJ energy

 end if
!          - for pressure calculation
           do n=1,3
              pot_virit_lj(1:3,n) = pot_virit_lj(1:3,n)   &
                   &              + fij(1:3)*rij(n)
           end do

           pot_viri_lj = pot_viri_lj   &
                &      + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

        end if

     END DO

  END DO

! ---- Fennnel method: return
  if (fennell_flag) then
     force(1:3,1:natom) = force(1:3,1:natom)   &
          &             + fvdw(1:3,1:natom) + felc(1:3,1:natom)
     return
  end if


!!! if calculate L-J potential only
#if defined(_LJ_ONLY)
  ! no execution
#else

! --- ADD AND SUBTRACT CONTRIBUTION OF EXCLUDED NEIGHBORS ---

! MPI loop
!      DO i = 1, natom - 1

  looplast = natom - 1
      
  DO i = loopinit, looplast, loopstep

     j1 = index_excl(i)
     j2 = index_excl(i+1) - 1
     IF (list_excl(j1) == 0) CYCLE

     DO j = j1, j2                    ! loop over excluded atoms for atom(i)

        jj = list_excl(j)

        rij(1:3)        = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       --- ADD EWALD REAL AND SUBTRACT COULOMB
!                          DUE TO EXCLUDED NEIGHBORS ---

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2
                                        ! |rij|^2
        rij_abs       = dsqrt(rij_abs2) ! |rij| 

!!! Spline interpolation for ewald real
#if defined(_NOT_SPLINTERP)
        rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2
        rij_abs_inv1  = 1.0d0 / rij_abs  ! 1/|rij|

        rtmp1 = rij_abs * alpha          ! |rij|*alpha
        rtmp2 = rtmp1 * rtmp1            ! |rij|^2/ehta^2
        rtmp_erf   = derf(rtmp1)         ! erf(|rij|*alpha)

        potij      = atmchrg(i) * atmchrg(jj) * rij_abs_inv1 

        ftmp       =  potij * ( -rtmp_erf * rij_abs_inv2   &
             &                + const_fewr * rij_abs_inv1 * dexp(-rtmp2))

                           ! =   Qi * Qj * [-erf(|rij|*alpha)/|rij|^3 
                           ! + 2*alpha/sqrt(pi)/|rij|^2 exp(-|rij|^2*alpha^2)]

        potij      = potij * rtmp_erf
                                     ! qi qj erf(|rij|*alpha) / |rij|
#else
        potij = atmchrg(i) * atmchrg(jj)
        ftmp = -SPLINE_INTERPFUNC(rij_abs,spltbl_excs,   &
           &                      spl_b_excs,spl_c_excs,spl_d_excs)   &
           & * potij
        potij = SPLINE_INTERPFUNC(rij_abs,spltbl_excspot,   &
            &                     spl_b_excspot,spl_c_excspot,spl_d_excspot)   &
            & * potij
#endif

        pot_elc    = pot_elc - potij
        pot_viri_coul =  pot_viri_coul - potij

        fij(1:3)     = ftmp * rij(1:3)

        felc(1:3,i)  = felc(1:3,i)  + fij(1:3)
        felc(1:3,jj) = felc(1:3,jj) - fij(1:3)

        for_viri_coul(1:3,i)  = for_viri_coul(1:3,i)  + fij(1:3)
        for_viri_coul(1:3,jj) = for_viri_coul(1:3,jj) - fij(1:3)

!       --- calculate virial tensor ---
        do n=1,3
           pot_virit_coul(1:3,n) = pot_virit_coul(1:3,n)   &
                &                + fij(1:3)*rij(n)
        end do

!       - for heat flux calculation
!       --- Bulk mode with Ewald calculation
#if defined(_HF_BULK_EWK)
        fji(1:3) = -fij(1:3)

        ! all contribution is assigned to intramolecular interaction
        pot_elin_atm(i) = pot_elin_atm(i) - 0.5d0 * potij
        pot_elin_atm(jj) = pot_elin_atm(jj) - 0.5d0 * potij

        call calvirihf(i,jj,fij,fji,   &
             &         nbodycoeff,   &
             &         box,box_inv,   &
             &         ifhfvol,   &
             &         nhfregion,hfzpos1,hfzpos2,   &
             &         hftyp_atm,   &
             &         molecom,   &
             &         virielint_atm)
#endif

     END DO
  END DO

#endif
!!! end if calculate L-J potential only

! --- ADD FELC and/or FVDW TO FORCE ---

  force(1:3,1:natom) = force(1:3,1:natom)   &
       &             + felc(1:3,1:natom) + fvdw(1:3,1:natom)

!     write(20,'(6f12.3)') felc(:,:)

!     +     +     +     +     +     +     +

end subroutine ewaldrp_hf
