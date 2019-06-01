!**************************************
!*  morse_hf.f90 Ver.1.5 '10.06.28    *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine morsep_hf(xcel,ycel,zcel,   &
     &               rcutmor,   &
     &               force,pot_mor,   &
     &               for_viri_mor,pot_viri_mor,   &
     &               pot_virit_mor,   &
     &               pot_mor_atm,virimort_atm,   &
     &               ifhfvol,   &
     &               nhfregion,hfzpos1,hfzpos2,   &
     &               hftyp_atm,   &
     &               molecom)

  use interface_interact, only: calvirihf

  use md_global
  use mpi_global

  implicit none
      
!
!    
!     Calculate MORSE interactions & PRESSURE VIRIAL
!        

! ARGUMENTS:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: ifhfvol      ! local volume-based or local surface-based

  real(8),intent(in):: rcutmor          ! Morse cutoff length [non-d]

  integer,intent(in):: nhfregion        ! number of region to calculate heat flux
  real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                        ! z-position of region for heat flux

  integer,intent(in):: hftyp_atm(:)     ! atom- or mole-based heat flux cal. 
                                        !   for each atom

  real(8),intent(in):: molecom(:,:)     ! center of mass of molecule

!   INPUT & OUTPUT

  real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc) 

  real(8),intent(inout):: for_viri_mor(:,:)  ! virial(Morse force) of each atom
  real(8),intent(inout):: pot_viri_mor   ! virial(Morse potential) of each atom
  real(8),intent(inout):: pot_virit_mor(:,:) ! virial tensor (Morse)

  real(8),intent(inout):: pot_mor_atm(:) ! Morse potential of each atom
  real(8),intent(inout):: virimort_atm(:,:,:,:)
                                         ! virial tensor of each atom (Morse)

!   OUTPUT

  real(8),intent(out):: pot_mor         ! Morse potential

! LOCAL:
  real(8):: fmor(3,maxnatom)            ! Morce force calculated here

  real(8):: potij                       ! temporal potential due to IJ  
  real(8):: ftmp                        ! temporal force  
  real(8):: fij(3)                      ! force to be subtracted due to 
                                        ! excluded atoms
  real(8):: fji(3)                      ! force to be subtracted due to 
                                        ! excluded atoms = -fij

  real(8):: rij(3)                      ! rj-ri
  real(8):: rij_abs                     ! |Rij|   
  real(8):: rij_abs2                    ! |Rij|**2
  real(8):: rij_abs_inv1                ! |Rij|**-1
  real(8):: rtmp1
  real(8):: rtmp2

  real(8):: box(3)                      ! BOX size
  real(8):: box_inv(3)                  ! inverse of BOX size

  real(8):: sqcut                       ! rcutmor * rcutmor

  integer:: itype, jtype                ! VDW atom types
  integer:: i,j,n                       ! do loop index
  integer:: j1, j2                      ! limits of do loop
  integer:: jj                          ! atom index of the interacting atom

  integer:: mor1

!  real(8):: rij_fij                     ! rij dyad fij (tensor product)

  real(8):: nbodycoeff             ! n-body coefficient for virial term of heatf

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

!     +     +     +     +     +     +     +


! --- ALLOCATE ARRAY FOR FORCE CALCULATED IN THIS SUBROUTINE --- 

  fmor(1:3,1:natom) = 0.0d0

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1:3) = 1.0d0/box(1:3)

  nbodycoeff = 0.5d0        ! = 1/2

  sqcut = rcutmor * rcutmor

! --- CALCULATE NONBONDED FORCE ---

  pot_mor = 0.0d0

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO mor1 = 1, NMORSE - 1

  looplast = nmorse - 1

  DO mor1 = loopinit, looplast, loopstep

     j1 = mor_indexall(mor1)              ! the first and
     j2 = mor_indexall(mor1+loopstep) - 1 ! the last atoms interacting
                                          ! with i

     i = nmorselist(mor1)
     itype = atmindex(i)                ! atmtype of atom(i)

!         im = irmolept_list(i)  ! what molecule atom(i) belongs in

     IF (mor_listall(j1) == 0) CYCLE

     DO j = j1, j2

        jj = mor_listall(j)             ! atom index of interacting atom

        jtype = atmindex(jj)            ! atmtype of atom(jj)

!            jm = irmolept_list(jj) ! what molecule atom(i) belongs in

!       --- CALCULATE VIRIAL & FORCE ---

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,jj)

!       --- periodic boundary ---

        rij(1:3) = rij(1:3) - box(1:3) * dnint(rij(1:3)*box_inv(1:3))

!       * box_inv rather than /box accellerates the computation

        rij_abs2      = rij(1)**2 + rij(2)**2 + rij(3)**2 ! |rij|^2
 
!       --- cut off ---
        if (rij_abs2 > sqcut) cycle

        rij_abs       = dsqrt(rij_abs2) ! |rij| 

!            rij_abs_inv2  = 1.0d0 / rij_abs2 ! 1/|rij|^2
        rij_abs_inv1  = 1.0d0 / rij_abs ! 1/|rij|

        rtmp1 = dexp(-1.0d0*(rij_abs - para_radmor(itype,jtype))   &
             &     * para_alphamor(itype,jtype)) ! exp(-(|rij|-r0)*alpha)
        rtmp2 = rtmp1 * rtmp1           ! exp(-2(|rij|-r0)*alpha)

        ftmp       = 2.0d0 * para_alphamor(itype,jtype)   &
             &     * para_welmor(itype,jtype)   &
             &     * (rtmp2 - rtmp1) * rij_abs_inv1 
                                        ! = 2*alpha*D*(rtmp2-rtmp1) / |rij|

        potij      = para_welmor(itype,jtype)   &
             &     * (rtmp2 - 2.0d0*rtmp1)
                                        ! = D * (rtmp2 - 2*rtmp1)

        pot_mor    = pot_mor + potij

        pot_mor_atm(i) = pot_mor_atm(i) + 0.5d0*potij
        pot_mor_atm(jj) = pot_mor_atm(jj) + 0.5d0*potij

        fij(1:3)     = ftmp * rij(1:3)

        fmor(1:3,i)  = fmor(1:3,i)  + fij(1:3)
        fmor(1:3,jj) = fmor(1:3,jj) - fij(1:3)

        for_viri_mor(1:3,i)  = for_viri_mor(1:3,i)  + fij(1:3)
        for_viri_mor(1:3,jj) = for_viri_mor(1:3,jj) - fij(1:3)

!       --- calculate virial tensor ---

!       - for heat flux calculation

        fji(1:3) = -fij(1:3)

        call calvirihf(i,jj,fij,fji,   &
             &         nbodycoeff,   &
             &         box,box_inv,   &
             &         ifhfvol,   &
             &         nhfregion,hfzpos1,hfzpos2,   &
             &         hftyp_atm,   &
             &         molecom,   &
             &         virimort_atm)

!       - for pressure calculation
        do n=1,3
           pot_virit_mor(1:3,n) = pot_virit_mor(1:3,n)   &
                &               + fij(1:3)*rij(n)
        end do
               
        pot_viri_mor = pot_viri_mor   &
             &       + fij(1)*rij(1) + fij(2)*rij(2) + fij(3)*rij(3)

     END DO

  END DO


! --- ADD FELC and/or FVDW TO FORCE ---

  force(1:3,1:natom) = force(1:3,1:natom) + fmor(1:3,1:natom)

!     write(20,'(6f12.3)') felc(:,:)

!     +     +     +     +     +     +     +

end subroutine morsep_hf
