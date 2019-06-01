!************************************
!*  interface_interact.f90 Ver.1.4  *
!*      for peachgk_md.f            *
!*            by G.Kikugawa         *
!************************************
! Time-stamp: <2015-07-16 17:49:57 gota>

!***** This module is interface module for routines of interactions *****

module interface_interact

  interface

     subroutine ewaldr(xcel,ycel,zcel,   &
          &            alpha,rrcut,   &
          &            rcut,   &
          &            force,pot_elc,pot_vdw)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       real(8),intent(in):: alpha          ! parameter alpha [non-d]
       real(8),intent(in):: rrcut          ! ewald real space cutoff length [non-d]
       
       real(8),intent(in):: rcut           ! vdw cutoff length [non-d]

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)  ! force calculated here
                                      ! total force
                                      !  (including vdw, bond, angle etc) 

       !   OUTPUT

       real(8),intent(out):: pot_elc       ! electrostatic potential
       real(8),intent(out):: pot_vdw       ! vdw potential

     end subroutine ewaldr

     subroutine ewaldrp(xcel,ycel,zcel,   &
          &             alpha,rrcut,   &
          &             rcut,   &
          &             force,pot_elc,pot_vdw,   &
          &             for_viri_coul,pot_viri_coul,   &
          &             for_viri_lj,pot_viri_lj,   &
          &             pot_virit_coul,pot_virit_lj)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       real(8),intent(in):: alpha          ! parameter alpha [non-d]
       real(8),intent(in):: rrcut          ! ewald real space cutoff length [non-d]

       real(8),intent(in):: rcut           ! vdw cutoff length [non-d]

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)  ! force calculated here
                                      ! total force
                                      ! (including vdw, bond, angle etc) 

       real(8),intent(inout):: for_viri_coul(:,:) ! virial(coulomb force) of each atom
       real(8),intent(inout):: pot_viri_coul ! virial(coulomb potential) of each atom
       real(8),intent(inout):: pot_virit_coul(:,:) ! virial tensor (coulomb)

       real(8),intent(inout):: for_viri_lj(:,:)  ! virial(L-J force) of each atom
       real(8),intent(inout):: pot_viri_lj       ! virial(L-J potential) of each atom
       real(8),intent(inout):: pot_virit_lj(:,:) ! virial tensor (L-J)

       !   OUTPUT

       real(8),intent(out):: pot_elc       ! electrostatic potential
       real(8),intent(out):: pot_vdw       ! vdw potential

     end subroutine ewaldrp

#if defined(HF)
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

     end subroutine ewaldrp_hf
#endif

     subroutine ewaldk(xcel,ycel,zcel,   &
          &            force,pot_ewk)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! total force
                                        !  (including vdw, bond, angle etc)

       !   OUTPUT

       real(8),intent(out):: pot_ewk         ! electrostatic potential

     end subroutine ewaldk

     subroutine ewaldkp(xcel,ycel,zcel,   &
          &             alpha,   &
          &             force,pot_ewk,   &
          &             ifnetqcorrp,   &
          &             netchrgsq,   &
          &             for_viri_coul,pot_viri_coul,pot_virit_coul)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: alpha            ! parameter alpha [non-d]

       logical,intent(in):: ifnetqcorrp      ! net charge correction for pressure
       real(8),intent(in):: netchrgsq        ! = (sum(qi))**2

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! total force
                                        !  (including vdw, bond, angle etc) 

       real(8),intent(inout):: for_viri_coul(:,:) ! virial(coulomb force) of each atom
       real(8),intent(inout):: pot_viri_coul ! virial(coulomb potential) of each atom
       real(8),intent(inout):: pot_virit_coul(:,:) ! virial tensor (coulomb)

       !   OUTPUT

       real(8),intent(out):: pot_ewk         ! electrostatic potential

     end subroutine ewaldkp

#if defined(HF)
     subroutine ewaldkp_hf(xcel,ycel,zcel, &
          &                alpha, &
          &                force,pot_ewk, &
          &                ifnetqcorrp, &
          &                netchrgsq, &
          &                for_viri_coul,pot_viri_coul,pot_virit_coul, &
          &                pot_elc_atm,virielct_atm, &
          &                ifhfvol, &
          &                nhfregion,hfzpos1,hfzpos2, &
          &                hftyp_atm, &
          &                molecom)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: alpha            ! parameter alpha [non-d]

       logical,intent(in):: ifnetqcorrp      ! net charge correction for pressure
       real(8),intent(in):: netchrgsq        ! = (sum(qi))**2

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! total force
                                        !  (including vdw, bond, angle etc)

       real(8),intent(inout):: for_viri_coul(:,:) ! virial(coulomb force) of each atom
       real(8),intent(inout):: pot_viri_coul ! virial(coulomb potential) of each atom
       real(8),intent(inout):: pot_virit_coul(:,:) ! virial tensor (coulomb)

       logical,intent(in):: ifhfvol       ! local volume-based or local surface-based

       integer,intent(in):: nhfregion     ! number of region to calculate heat flux
       real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                     ! z-position of region for heat flux

       integer,intent(in):: hftyp_atm(:)  ! atom- or mole-based heat flux cal.
                                     !   for each atom

       real(8),intent(in):: molecom(:,:)  ! center of mass of molecule

       !   OUTPUT

       real(8),intent(out):: pot_ewk         ! electrostatic potential

       real(8),intent(inout):: pot_elc_atm(:)   ! Coulomb potential of each atom
       real(8),intent(inout):: virielct_atm(:,:,:,:)
                                     ! virial tensor of each atom (Coulomb)

     end subroutine ewaldkp_hf

#endif

     subroutine fft_pme(xcel,ycel,zcel,   &
          &             alpha,   &
          &             nfft1,nfft2,nfft3,pme_order,   &
          &             ifnetqcorrp,   &
          &             netchrgsq,   &
          &             force,pot_ewk)

       ! ARGUMENT:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: alpha            ! parameter alpha [non-d]

       integer,intent(in):: nfft1, nfft2, nfft3 ! grid points in SPME
       integer,intent(in):: pme_order        ! B-spline order

       logical,intent(in):: ifnetqcorrp      ! net charge correction for pressure
       real(8),intent(in):: netchrgsq        ! = (sum(qi))**2

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! total force

       !   OUTPUT
       real(8),intent(out):: pot_ewk         ! electrostatic potential

     end subroutine fft_pme

     subroutine fft_pmep(xcel,ycel,zcel,   &
          &              alpha,   &
          &              nfft1,nfft2,nfft3,pme_order,   &
          &              force,pot_ewk,   &
          &              ifnetqcorrp,   &
          &              netchrgsq,   &
          &              for_viri_coul,pot_viri_coul,pot_virit_coul)

       ! ARGUMENT:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: alpha            ! parameter alpha [non-d]

       integer,intent(in):: nfft1, nfft2, nfft3 ! grid points in SPME
       integer,intent(in):: pme_order        ! B-spline order

       logical,intent(in):: ifnetqcorrp      ! net charge correction for pressure
       real(8),intent(in):: netchrgsq        ! = (sum(qi))**2

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! total force

       !   OUTPUT
       real(8),intent(out):: pot_ewk         ! electrostatic potential

       real(8),intent(out):: for_viri_coul(:,:) ! virial(coulomb force) of each atom
       real(8),intent(out):: pot_viri_coul   ! virial(coulomb potential) of each atom
       real(8),intent(out):: pot_virit_coul(:,:) ! virial tensor (coulomb)

     end subroutine fft_pmep

     subroutine get_scaled_fractionals(nfft1,nfft2,nfft3,   &
          &                            recip,   &
          &                            fr1,fr2,fr3)

       ! ARGUMENT:
       !   INPUT
       integer,intent(in):: nfft1, nfft2, nfft3 ! number of grid points in SPME
       real(8),intent(in):: recip(:)

       !   OUTPUT
       real(8),intent(out):: fr1(:),fr2(:),fr3(:)

     end subroutine get_scaled_fractionals

     subroutine get_bspline_coeffs(pme_order,   &
          &                        fr1,fr2,fr3,   &
          &                        theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)

       use md_global

       integer,intent(in):: pme_order       ! B-spline order
       real(8),intent(in):: fr1(:),fr2(:),fr3(:)

       real(8),intent(out):: theta1(pme_order,natom),theta2(pme_order,natom),   &
            &                theta3(pme_order,natom)
       real(8),intent(out):: dtheta1(pme_order,natom),   &
            &                dtheta2(pme_order,natom),dtheta3(pme_order,natom)

     end subroutine get_bspline_coeffs

     subroutine fill_charge_grid(theta1,theta2,theta3,fr1,fr2,fr3,   &
          &                      pme_order,nfft1,nfft2,nfft3,   &
          &                      nfftdim1,nfftdim2,nfftdim3,Q)

       use md_global

       ! ARGUMENT:
       integer,intent(in):: pme_order,nfft1,nfft2,nfft3
       integer,intent(in):: nfftdim1,nfftdim2,nfftdim3
       real(8),intent(in):: fr1(:),fr2(:),fr3(:)
       real(8),intent(out):: Q(2,nfftdim1,nfftdim2,nfftdim3)

       real(8),intent(in):: theta1(pme_order,natom),theta2(pme_order,natom),   &
            &               theta3(pme_order,natom)

     end subroutine fill_charge_grid

#if defined(_FFTW3)
     subroutine fft_back(array,   &
          &              nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3)

       ! ARGUMENT:
       real(8),intent(out):: array(2,nfftdim1,nfftdim2,nfftdim3)

       integer,intent(in):: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

     end subroutine fft_back
#else
     subroutine fft_back(array,   &
          &              nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,   &
          &              nfftable,nffwork)

       ! ARGUMENT:
       real(8),intent(out):: array(:)
       integer,intent(in):: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

       integer,intent(out):: nfftable,nffwork

     end subroutine fft_back
#endif

#if defined(_FFTW3)
     subroutine fft_forward(array,   &
          &                 nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3)

       ! ARGUMENT:
       real(8),intent(out):: array(2,nfftdim1,nfftdim2,nfftdim3)
       integer,intent(in):: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

     end subroutine fft_forward
#else
     subroutine fft_forward(array,   &
          &                 nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,   &
          &                 nfftable,nffwork)

       ! ARGUMENT:
       real(8),intent(out):: array(:)
       integer,intent(in):: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
       integer,intent(out):: nfftable,nffwork

     end subroutine fft_forward
#endif

     subroutine scalar_sum(alpha,pot_ewk,   &
          &     Q,vol,recip,   &
          &     nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,vir,   &
          &     vir_sca,   &
          &     ifnetqcorrp,   &
          &     netchrgsq)

       ! ARGUMENT:
       real(8),intent(in):: alpha           ! parameter alpha [non-d]

       integer,intent(in):: nfft1, nfft2, nfft3 ! grid points in SPME
       integer,intent(in):: nfftdim1,nfftdim2,nfftdim3

       real(8),intent(out):: pot_ewk        ! electrostatic potential

       real(8),intent(inout):: Q(2,nfftdim1,nfftdim2,nfftdim3)

       real(8),intent(in):: vol
       real(8),intent(in):: recip(:)

       real(8),intent(out):: vir(:,:)       ! virial tensor
       real(8),intent(out):: vir_sca        ! virial scalar for charge correction

       logical,intent(in):: ifnetqcorrp     ! net charge correction for pressure
       real(8),intent(in):: netchrgsq       ! = (sum(qi))**2

     end subroutine scalar_sum

     subroutine grad_sum(recip,   &
          &              theta1,theta2,theta3,   &
          &              dtheta1,dtheta2,dtheta3,   &
          &              force,   &
          &              fr1,fr2,fr3,   &
          &              pme_order,nfft1,nfft2,nfft3,   &
          &              nfftdim1,nfftdim2,nfftdim3,   &
          &              Q,   &
          &              vir_for)

       use md_global

       ! ARGUMENT:
       integer,intent(in):: pme_order,nfft1,nfft2,nfft3
       integer,intent(in):: nfftdim1,nfftdim2,nfftdim3
       real(8),intent(in):: recip(:)
       real(8),intent(in):: fr1(:),fr2(:),fr3(:)
       
       real(8),intent(inout):: force(:,:)

       real(8),intent(in):: theta1(pme_order,natom),theta2(pme_order,natom),   &
            &               theta3(pme_order,natom)
       real(8),intent(in):: dtheta1(pme_order,natom),dtheta2(pme_order,natom),   &
            &               dtheta3(pme_order,natom)
       real(8),intent(in):: Q(2,nfftdim1,nfftdim2,nfftdim3)

       real(8),intent(inout):: vir_for(:,:)

     end subroutine grad_sum

#if defined(MPI)
     subroutine FFTdata_trans_xz(send_array,array,   &
          &                      nfft1,nfft2,nfft3,   &
          &                      nfftdim1,nfftdim2,nfftdim3)

       use mpi_global

       ! ARGUMENT:
       !   INPUT
       real(8),intent(in):: send_array(2,nfft2,zlimitmax,xlimitmax*nproc)

       integer,intent(in):: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

       !   OUTPUT
       real(8),intent(out):: array(2,nfftdim1,nfftdim2,nfftdim3)

     end subroutine FFTdata_trans_xz

     subroutine FFTdata_trans_zx(send_array,array,   &
          &                      nfft1,nfft2,nfft3,   &
          &                      nfftdim1,nfftdim2,nfftdim3)

       use mpi_global

       ! ARGUMENT:
       !   INPUT
       real(8),intent(in):: send_array(2,nfft2,xlimitmax,zlimitmax*nproc)

       integer,intent(in):: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

       !   OUTPUT
       real(8),intent(out):: array(2,nfftdim1,nfftdim2,nfftdim3)

     end subroutine FFTdata_trans_zx
#endif

     subroutine morse(xcel,ycel,zcel,   &
          &           rcutmor,   &
          &           force,pot_mor)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel            ! x cell length[non-d]
       real(8),intent(in):: ycel            ! y cell length[non-d]
       real(8),intent(in):: zcel            ! z cell length[non-d]

       real(8),intent(in):: rcutmor         ! Morse cutoff length [non-d]

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)   ! force calculated here
                                       ! total force
                                       !  (including vdw, bond, angle etc) 

       !   OUTPUT

       real(8),intent(out):: pot_mor        ! Morse potential

     end subroutine morse

     subroutine morsep(xcel,ycel,zcel,   &
          &            rcutmor,   &
          &            force,pot_mor,   &
          &            for_viri_mor,pot_viri_mor,   &
          &            pot_virit_mor)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: rcutmor          ! Morse cutoff length [non-d]

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc) 

       real(8),intent(inout):: for_viri_mor(:,:)  ! virial(Morse force) of each atom
       real(8),intent(inout):: pot_viri_mor   ! virial(Morse potential) of each atom
       real(8),intent(inout):: pot_virit_mor(:,:) ! virial tensor (Morse)

       !   OUTPUT

       real(8),intent(out):: pot_mor         ! Morse potential

     end subroutine morsep

#if defined(HF)
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

     end subroutine morsep_hf
#endif

     subroutine shpot(xcel,ycel,zcel,   &
          &           rcutsh,   &
          &           force,pot_sh)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: rcutsh           ! SH cutoff length [non-d]

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc) 

       !   OUTPUT

       real(8),intent(out):: pot_sh          ! SH potential

     end subroutine shpot

     subroutine shpotp(xcel,ycel,zcel,   &
          &            rcutsh,   &
          &            force,pot_sh,   &
          &            for_viri_sh,pot_viri_sh,   &
          &            pot_virit_sh)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: rcutsh           ! SH cutoff length [non-d]

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)   ! force calculated here
                                       ! total force
                                       !  (including vdw, bond, angle etc)

       real(8),intent(inout):: for_viri_sh(:,:)   ! virial(SH force) of each atom
       real(8),intent(inout):: pot_viri_sh        ! virial(SH potential)
       real(8),intent(inout):: pot_virit_sh(:,:)  ! virial tensor (SH)

       !   OUTPUT

       real(8),intent(out):: pot_sh         ! SH potential

     end subroutine shpotp

#if defined(HF)
     subroutine shpotp_hf(xcel,ycel,zcel,   &
          &               rcutsh,   &
          &               force,pot_sh,   &
          &               for_viri_sh,pot_viri_sh,   &
          &               pot_virit_sh,   &
          &               pot_sh_atm,virisht_atm,   &
          &               ifhfvol,   &
          &               nhfregion,hfzpos1,hfzpos2,   &
          &               hftyp_atm,   &
          &               molecom)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel            ! x cell length[non-d]
       real(8),intent(in):: ycel            ! y cell length[non-d]
       real(8),intent(in):: zcel            ! z cell length[non-d]

       real(8),intent(in):: rcutsh          ! SH cutoff length [non-d]

       logical,intent(in):: ifhfvol      ! local volume-based or local surface-based

       integer,intent(in):: nhfregion       ! number of region to calculate heat flux
       real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                       ! z-position of region for heat flux

       integer,intent(in):: hftyp_atm(:)    ! atom- or mole-based heat flux cal. 
                                       !   for each atom

       real(8),intent(in):: molecom(:,:)    ! center of mass of molecule

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)   ! force calculated here
                                       ! total force
                                       !  (including vdw, bond, angle etc)

       real(8),intent(inout):: for_viri_sh(:,:)  ! virial(SH force) of each atom
       real(8),intent(inout):: pot_viri_sh       ! virial(SH potential)
       real(8),intent(inout):: pot_virit_sh(:,:) ! virial tensor (SH)

       real(8),intent(inout):: pot_sh_atm(:) ! SH potential of each atom
       real(8),intent(inout):: virisht_atm(:,:,:,:) 
                                       ! virial tensor of each atom (SH)

       !   OUTPUT

       real(8),intent(out):: pot_sh         ! SH potential

     end subroutine shpotp_hf
#endif

     subroutine rfhpot(xcel,ycel,zcel,   &
          &            rcutrfhfo,rcutrfhoo,rcutrfhoh,   &
          &            force,pot_rfh)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: rcutrfhfo        ! RFH(FO) cutoff length [non-d]
       real(8),intent(in):: rcutrfhoo        ! RFH(OO) cutoff length [non-d]
       real(8),intent(in):: rcutrfhoh        ! RFH(OH) cutoff length [non-d]

       !  INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc) 

       !  OUTPUT

       real(8),intent(out):: pot_rfh         ! RFH potential

     end subroutine rfhpot

     subroutine rfhpotp(xcel,ycel,zcel,   &
          &             rcutrfhfo,rcutrfhoo,rcutrfhoh,   &
          &             force,pot_rfh,   &
          &             for_viri_rfh,pot_viri_rfh,   &
          &             pot_virit_rfh)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: rcutrfhfo        ! RFH(FO) cutoff length [non-d]
       real(8),intent(in):: rcutrfhoo        ! RFH(OO) cutoff length [non-d]
       real(8),intent(in):: rcutrfhoh        ! RFH(OH) cutoff length [non-d]

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc)

       real(8),intent(inout):: for_viri_rfh(:,:)  ! virial(RFH force) of each atom
       real(8),intent(inout):: pot_viri_rfh       ! virial(RFH potential)
       real(8),intent(inout):: pot_virit_rfh(:,:) ! virial tensor (RFH)

       !   OUTPUT

       real(8),intent(out):: pot_rfh         ! RFH potential

     end subroutine rfhpotp

#if defined(HF)
     subroutine rfhpotp_hf(xcel,ycel,zcel,   &
          &                rcutrfhfo,rcutrfhoo,rcutrfhoh,   &
          &                force,pot_rfh,   &
          &                for_viri_rfh,pot_viri_rfh,   &
          &                pot_virit_rfh,   &
          &                pot_rfh_atm,virirfht_atm,   &
          &                ifhfvol,   &
          &                nhfregion,hfzpos1,hfzpos2,   &
          &                hftyp_atm,   &
          &                molecom)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: rcutrfhfo        ! RFH(FO) cutoff length [non-d]
       real(8),intent(in):: rcutrfhoo        ! RFH(OO) cutoff length [non-d]
       real(8),intent(in):: rcutrfhoh        ! RFH(OH) cutoff length [non-d]

       logical,intent(in):: ifhfvol       ! local volume-based or local surface-based

       integer,intent(in):: nhfregion       ! number of region to calculate heat flux
       real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                        ! z-position of region for heat flux

       integer,intent(in):: hftyp_atm(:)     ! atom- or mole-based heat flux cal. 
                                        !   for each atom

       real(8),intent(in):: molecom(:,:)     ! center of mass of molecule

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc)

       real(8),intent(inout):: for_viri_rfh(:,:)  ! virial(RFH force) of each atom
       real(8),intent(inout):: pot_viri_rfh       ! virial(RFH potential)
       real(8),intent(inout):: pot_virit_rfh(:,:) ! virial tensor (RFH)

       real(8),intent(inout):: pot_rfh_atm(:) ! RFH potential of each atom
       real(8),intent(inout):: virirfht_atm(:,:,:,:)
                                        ! virial tensor of each atom (RFH)

       !   OUTPUT
       
       real(8),intent(out):: pot_rfh         ! RFH potential

     end subroutine rfhpotp_hf
#endif

     subroutine doupot(xcel,ycel,zcel,   &
          &            rcutdouo,rcutindouo,rcutdouh,   &
          &            force,pot_dou)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: rcutdouo         ! DOU cutoff length [non-d] for O-Au
       real(8),intent(in):: rcutindouo       ! DOU cutin length [non-d] for O-Au
       real(8),intent(in):: rcutdouh         ! DOU cutoff length [non-d] for H-Au

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc) 

       !   OUTPUT

       real(8),intent(out):: pot_dou         ! DOU potential

     end subroutine doupot

     subroutine doupotp(xcel,ycel,zcel,   &
          &             rcutdouo,rcutindouo,rcutdouh,   &
          &             force,pot_dou,   &
          &             for_viri_dou,pot_viri_dou,   &
          &             pot_virit_dou)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: rcutdouo         ! DOU cutoff length [non-d] for O-Au
       real(8),intent(in):: rcutindouo       ! DOU cutin length [non-d] for O-Au
       real(8),intent(in):: rcutdouh         ! DOU cutoff length [non-d] for H-Au

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc)

       real(8),intent(inout):: for_viri_dou(:,:)  ! virial(DOU force) of each atom
       real(8),intent(inout):: pot_viri_dou       ! virial(DOU potential)
       real(8),intent(inout):: pot_virit_dou(:,:) ! virial tensor (DOU)

       !   OUTPUT

       real(8),intent(out):: pot_dou         ! DOU potential

     end subroutine doupotp

#if defined(HF)
     subroutine doupotp_hf(xcel,ycel,zcel,   &
          &                rcutdouo,rcutindouo,rcutdouh,   &
          &                force,pot_dou,   &
          &                for_viri_dou,pot_viri_dou,   &
          &                pot_virit_dou,   &
          &                pot_dou_atm,viridout_atm,   &
          &                ifhfvol,   &
          &                nhfregion,hfzpos1,hfzpos2,   &
          &                hftyp_atm,   &
          &                molecom)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       real(8),intent(in):: rcutdouo         ! DOU cutoff length [non-d] for O-Au
       real(8),intent(in):: rcutindouo       ! DOU cutin length [non-d] for O-Au
       real(8),intent(in):: rcutdouh         ! DOU cutoff length [non-d] for H-Au

       logical,intent(in):: ifhfvol       ! local volume-based or local surface-based

       integer,intent(in):: nhfregion       ! number of region to calculate heat flux
       real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                        ! z-position of region for heat flux

       integer,intent(in):: hftyp_atm(:)     ! atom- or mole-based heat flux cal. 
                                        !   for each atom

       real(8),intent(in):: molecom(:,:)     ! center of mass of molecule

       !   INPUT & OUTPUT

       real(8),intent(inout):: force(:,:)    ! force calculated here
                                        ! total force
                                        !  (including vdw, bond, angle etc)

       real(8),intent(inout):: for_viri_dou(:,:)  ! virial(DOU force) of each atom
       real(8),intent(inout):: pot_viri_dou       ! virial(DOU potential)
       real(8),intent(inout):: pot_virit_dou(:,:) ! virial tensor (DOU)

       real(8),intent(inout):: pot_dou_atm(:) ! DOU potential of each atom
       real(8),intent(inout):: viridout_atm(:,:,:,:) 
                                        ! virial tensor of each atom (DOU)

       !   OUTPUT

       real(8),intent(out):: pot_dou         ! DOU potential

     end subroutine doupotp_hf
#endif

     subroutine calbond(xcel,ycel,zcel,   &
          &             force,pot_bond)
  
       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_bond        ! bond potential


     end subroutine calbond

     subroutine calbondp(xcel,ycel,zcel,   &
          &              force,pot_bond,   &
          &              atm_viri_bond,atm_virit_bond)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_bond        ! bond potential

       real(8),intent(out):: atm_viri_bond   ! virial of each atom
       real(8),intent(out):: atm_virit_bond(:,:) ! virial tensor (bond potential)

     end subroutine calbondp

#if defined(HF)
     subroutine calbondp_hf(xcel,ycel,zcel,   &
          &                 force,pot_bond,   &
          &                 atm_viri_bond,atm_virit_bond,   &
          &                 pot_bo_atm,viribot_atm,   &
          &                 ifhfvol,   &
          &                 nhfregion,hfzpos1,hfzpos2,   &
          &                 hftyp_atm,   &
          &                 molecom)

       ! ARGUMENTS:
       !   INPUT
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
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_bond        ! bond potential

       real(8),intent(out):: atm_viri_bond   ! virial of each atom
       real(8),intent(out):: atm_virit_bond(:,:) ! virial tensor (bond potential)

       real(8),intent(out):: pot_bo_atm(:)   ! bond potential of each atom

       real(8),intent(out):: viribot_atm(:,:,:,:)
                                        ! virial tensor of each atom (bond)

     end subroutine calbondp_hf
#endif

     subroutine calangl(xcel,ycel,zcel,   &
          &             force,pot_angl)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_angl        ! angle energy 

     end subroutine calangl

     subroutine calanglp(xcel,ycel,zcel,   &
          &              force,pot_angl,   &
          &              atm_viri_angl,atm_virit_angl)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_angl        ! angle energy 

       real(8),intent(out):: atm_viri_angl   ! virial(angle potential) of each atom
       real(8),intent(out):: atm_virit_angl(:,:) ! virial tensor (angle potential)

     end subroutine calanglp

#if defined(HF)
     subroutine calanglp_hf(xcel,ycel,zcel,   &
          &                 force,pot_angl,   &
          &                 atm_viri_angl,atm_virit_angl,   &
          &                 pot_an_atm,viriant_atm,   &
          &                 ifhfvol,   &
          &                 nhfregion,hfzpos1,hfzpos2,   &
          &                 hftyp_atm,   &
          &                 molecom)

       ! ARGUMENTS:
       !   INPUT
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
       real(8),intent(inout):: force(:,:)       ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_angl        ! angle energy 

       real(8),intent(out):: atm_viri_angl   ! virial(angle potential) of each atom
       real(8),intent(out):: atm_virit_angl(:,:) ! virial tensor (angle potential)

       real(8),intent(out):: pot_an_atm(:)   ! angl potential of each atom

       real(8),intent(out):: viriant_atm(:,:,:,:)
                                        ! virial tensor of each atom (angle)

     end subroutine calanglp_hf
#endif

     subroutine calanglub(xcel,ycel,zcel, &
          &               force,pot_anglub)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_anglub      ! Urey-Bradley angle potential

     end subroutine calanglub

     subroutine calanglubp(xcel,ycel,zcel, &
          &                force,pot_anglub, &
          &                atm_viri_angl,atm_virit_angl)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_anglub      ! Urey-Bradley angle potential

       real(8),intent(out):: atm_viri_angl   ! virial(angle potential) of each atom
       real(8),intent(out):: atm_virit_angl(:,:) ! virial tensor (angle potential)

     end subroutine calanglubp

#if defined(HF)
     subroutine calanglubp_hf(xcel,ycel,zcel, &
          &                   force,pot_anglub, &
          &                   atm_viri_angl,atm_virit_angl, &
          &                   pot_an_atm,viriant_atm, &
          &                   ifhfvol, &
          &                   nhfregion,hfzpos1,hfzpos2, &
          &                   hftyp_atm, &
          &                   molecom)

       ! ARGUMENTS:
       !   INPUT
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
       real(8),intent(inout):: force(:,:)       ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_anglub      ! Urey-Bradley angle potential

       real(8),intent(out):: atm_viri_angl   ! virial(angle potential) of each atom
       real(8),intent(out):: atm_virit_angl(:,:) ! virial tensor (angle potential)

       real(8),intent(out):: pot_an_atm(:)   ! angl potential of each atom

       real(8),intent(out):: viriant_atm(:,:,:,:)
                                        ! virial tensor of each atom (angle)

     end subroutine calanglubp_hf
#endif

     subroutine caltors(xcel,ycel,zcel,   &
          &             force,pot_tors)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_tors        ! torsion potential

     end subroutine caltors

     subroutine caltorsp(xcel,ycel,zcel,   &
          &              force,pot_tors,   &
          &              atm_viri_tors,atm_virit_tors)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_tors        ! torsion potential

       real(8),intent(out):: atm_viri_tors   ! virial(torsion potential) of each atom
       real(8),intent(out):: atm_virit_tors(:,:) ! virial tensor (torsion potential)

     end subroutine caltorsp

#if defined(HF)
     subroutine caltorsp_hf(xcel,ycel,zcel,   &
          &                 force,pot_tors,   &
          &                 atm_viri_tors,atm_virit_tors,   &
          &                 pot_to_atm,viritot_atm,   &
          &                 ifhfvol,   &
          &                 nhfregion,hfzpos1,hfzpos2,   &
          &                 hftyp_atm,   &
          &                 molecom)

       ! ARGUMENTS:
       !   INPUT
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
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_tors        ! torsion potential

       real(8),intent(out):: atm_viri_tors   ! virial(torsion potential) of each atom
       real(8),intent(out):: atm_virit_tors(:,:) ! virial tensor (torsion potential)
       
       real(8),intent(out):: pot_to_atm(:)   ! torsion potential of each atom

       real(8),intent(out):: viritot_atm(:,:,:,:) 
                                        ! virial tensor of each atom (torsion)

     end subroutine caltorsp_hf
#endif

     subroutine calijkl(rij, rkj, rkl, rik, rjl,   &
          &             divfac, barhig, phsang, period,   &
          &             delta_pot, fi,fj,fk,fl, phi)

       ! ARGUMENTS:
       !   INPUT
       !   -- vectors ---

       real(8),intent(in)::  rij(3)          ! Rj-Ri 
       real(8),intent(in)::  rkj(3)          ! Rj-Rk 
       real(8),intent(in)::  rkl(3)          ! Rl-Rk 
       real(8),intent(in)::  rik(3)          ! Rk-Ri 
       real(8),intent(in)::  rjl(3)          ! Rl-Rj 

       ! -- force field parameters ---

       integer,intent(in):: divfac           ! division factor
       real(8),intent(in):: barhig           ! barrier height (barhig)
       real(8),intent(in):: phsang           ! phase angle    (phsang)
       real(8),intent(in):: period           ! periodicity    (period)

       !   OUTPUT

       real(8),intent(out):: delta_pot       ! potential calculated 
                                        !  in this subroutine  

       real(8),intent(out):: fi(3), fj(3), fk(3), fl(3) ! force to atoms (i,j,k,l)
                                                   !  calculated here
       real(8),intent(out):: phi             ! torsion (in radian) defined by
                                        !  atoms i, j, k, l 
                                        !  -pi =< phi < pi

     end subroutine calijkl

     subroutine caltorsrb(xcel,ycel,zcel,   &
          &               force,pot_torsrb)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_torsrb      ! RBtorsion potential

     end subroutine caltorsrb

     subroutine caltorsrbp(xcel,ycel,zcel,   &
          &                force,pot_torsrb,   &
          &                atm_viri_tors,atm_virit_tors)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_torsrb      ! RBtorsion potential

       real(8),intent(out):: atm_viri_tors   ! virial(torsion potential) of each atom
       real(8),intent(out):: atm_virit_tors(:,:) ! virial tensor (torsion potential)

     end subroutine caltorsrbp

#if defined(HF)
     subroutine caltorsrbp_hf(xcel,ycel,zcel,   &
          &                   force,pot_torsrb,   &
          &                   atm_viri_tors,atm_virit_tors,   &
          &                   pot_to_atm,viritot_atm,   &
          &                   ifhfvol,   &
          &                   nhfregion,hfzpos1,hfzpos2,   &
          &                   hftyp_atm,   &
          &                   molecom)

       ! ARGUMENTS:
       !   INPUT
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
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_torsrb      ! RBtorsion potential

       real(8),intent(out):: atm_viri_tors   ! virial(torsion potential) of each atom
       real(8),intent(out):: atm_virit_tors(:,:) ! virial tensor (torsion potential)

       real(8),intent(out):: pot_to_atm(:)   ! torsion potential of each atom

       real(8),intent(out):: viritot_atm(:,:,:,:)
                                        ! virial tensor of each atom (torsion)

     end subroutine caltorsrbp_hf
#endif

     subroutine calijklrb(rij, rkj, rkl, rik, rjl,   &
          &               barhig,   &
          &               delta_pot, fi,fj,fk,fl, phi)

       ! ARGUMENTS:
       !   INPUT
       !   -- vectors ---

       real(8),intent(in)::  rij(3)          ! Rj-Ri 
       real(8),intent(in)::  rkj(3)          ! Rj-Rk 
       real(8),intent(in)::  rkl(3)          ! Rl-Rk 
       real(8),intent(in)::  rik(3)          ! Rk-Ri 
       real(8),intent(in)::  rjl(3)          ! Rl-Rj 

       !   -- force field parameters ---

       real(8),intent(in):: barhig(:)        ! barrier height (barhig)

       !   OUTPUT

       real(8),intent(out):: delta_pot       ! potential calculated
                                        !  in this subroutine  

       real(8),intent(out):: fi(3), fj(3), fk(3), fl(3) ! force to atoms (i,j,k,l)
                                                   !  calculated here
       real(8),intent(out):: phi             ! torsion (in radian) defined by
                                        !  atoms i, j, k, l 
                                        !  -pi =< phi < pi

     end subroutine calijklrb

     subroutine caltorsim(xcel,ycel,zcel,   &
          &               force,pot_torsim)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_torsim      ! improper torsion potential

     end subroutine caltorsim

     subroutine caltorsimp(xcel,ycel,zcel,   &
          &                force,pot_torsim,   &
          &                atm_viri_tors,atm_virit_tors)

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

     end subroutine caltorsimp

#if defined(HF)
     subroutine caltorsimp_hf(xcel,ycel,zcel,   &
          &                   force,pot_torsim,   &
          &                   atm_viri_tors,atm_virit_tors,   &
          &                   pot_to_atm,viritot_atm,   &
          &                   ifhfvol,   &
          &                   nhfregion,hfzpos1,hfzpos2,   &
          &                   hftyp_atm,   &
          &                   molecom)

       ! ARGUMENTS:
       !   INPUT
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
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_torsim      ! improper torsion potential

       real(8),intent(out):: atm_viri_tors   ! virial(torsion potential) of each atom
       real(8),intent(out):: atm_virit_tors(:,:) ! virial tensor (torsion potential)

       real(8),intent(out):: pot_to_atm(:)   ! torsion potential of each atom

       real(8),intent(out):: viritot_atm(:,:,:,:)
                                        ! virial tensor of each atom (torsion)

     end subroutine caltorsimp_hf
#endif

     subroutine calijklim(rij, rkj, rkl, rik, rjl,   &
          &               barhig, phsang,   &
          &               delta_pot, fi,fj,fk,fl, phi)

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

     end subroutine calijklim

     subroutine calnon14(xcel,ycel,zcel, &
          &      div_factor_14vdw,div_factor_14elc,   &
          &      force,pot_elc14,pot_vdw14)

       ! ARGUMENTS:
       !   INPUT

       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       real(8),intent(in):: div_factor_14vdw ! division factor of 14vdw
       real(8),intent(in):: div_factor_14elc ! division factor of 14elc

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! atomic force

       !   OUTPUT
       real(8),intent(out):: pot_elc14       ! 14 electrostatic energy
       real(8),intent(out):: pot_vdw14       ! 14 vdw  energy 

     end subroutine calnon14

     subroutine calnon14p(xcel,ycel,zcel, &
          &               div_factor_14vdw,div_factor_14elc,   &
          &               force,pot_elc14,pot_vdw14,   &
          &               atm_viri_14,atm_virit_14)

       ! ARGUMENTS:
       !   INPUT

       real(8),intent(in):: xcel           ! x cell length[non-d]
       real(8),intent(in):: ycel           ! y cell length[non-d]
       real(8),intent(in):: zcel           ! z cell length[non-d]

       real(8),intent(in):: div_factor_14vdw ! division factor of 14vdw
       real(8),intent(in):: div_factor_14elc ! division factor of 14elc

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! atomic force

       !   OUTPUT
       real(8),intent(out):: atm_viri_14   ! virial(1-4 force potential) of each atom
       real(8),intent(out):: atm_virit_14(:,:) ! virial tensor (1-4 force potential)

       real(8),intent(out):: pot_elc14       ! 14 electrostatic energy
       real(8),intent(out):: pot_vdw14       ! 14 vdw  energy 

     end subroutine calnon14p

#if defined(HF)
     subroutine calnon14p_hf(xcel,ycel,zcel,   &
          &                  div_factor_14vdw,div_factor_14elc,   &
          &                  force,pot_elc14,pot_vdw14,   &
          &                  atm_viri_14,atm_virit_14,   &
          &                  pot_14_atm,viri14t_atm,   &
          &                  ifhfvol,   &
          &                  nhfregion,hfzpos1,hfzpos2,   &
          &                  hftyp_atm,   &
          &                  molecom)

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

     end subroutine calnon14p_hf
#endif

     subroutine cnpvw( rcutrpvw, force, pot_rpvw, pot_vw )

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

     end subroutine cnpvw

     subroutine calposres(xcel,ycel,zcel,   &
          &               force,pot_posres)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_posres      ! position restraint potential

     end subroutine calposres

     subroutine calpotbias(xcel,ycel,zcel,   &
          &                npoly,nwater,nmatom, &
          &                force,pot_pbias)

       ! ARGUMENTS:
       !   INPUT
       real(8),intent(in):: xcel             ! x cell length[non-d]
       real(8),intent(in):: ycel             ! y cell length[non-d]
       real(8),intent(in):: zcel             ! z cell length[non-d]

       integer,intent(in):: npoly           ! all number of poly
       integer,intent(in):: nwater          ! number of H2O molecules
       integer,intent(in):: nmatom          ! number of monatomic molecules

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       real(8),intent(out):: pot_pbias       ! bias potential

     end subroutine calpotbias

     subroutine ent_calcstmnb(mts_flag,mts_cstmnb, &
          &                   xcel,ycel,zcel, &
          &                   npoly,npolytyp,npoly_mole,npoly_atom, &
          &                   nwater,nmatom,nmatyp,nmatomtyp, &
          &                   force,pot_cstmnb, &
          &                   for_viri_cstmnb,pot_viri_cstmnb, &
          &                   pot_virit_cstmnb, &
          &                   ncstmnbex, &
          &                   pot_cstmnbex, &
          &                   for_viri_cstmnbex,pot_viri_cstmnbex, &
          &                   pot_virit_cstmnbex, &
          &                   ifcalpreatom,ifcalpremole, &
          &                   current_step, &
          &                   nstep_short,istep_short)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: mts_flag        ! flag for MTS integration
                                       ! 1 long-range force mode
                                       ! 2 medium-range force mode
                                       ! 3 short-range force mode
       integer,intent(in):: mts_cstmnb      ! MTS flag for custom NB interaction

       integer,intent(in):: npoly           ! all number of poly
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
  
       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       real(8),intent(in):: xcel            ! x cell length[non-d]
       real(8),intent(in):: ycel            ! y cell length[non-d]
       real(8),intent(in):: zcel            ! z cell length[non-d]

       logical,intent(in):: ifcalpremole    ! pressure calculation of molecule
       logical,intent(in):: ifcalpreatom    ! pressure calculation of atom

       integer,intent(in):: current_step    ! current time step
       integer,intent(in):: nstep_short     ! number of step for short force
       integer,intent(in):: istep_short

       integer,intent(in):: ncstmnbex       ! number of extra custom NB output

       !     OUTPUT
       real(8),intent(inout):: force(:,:)   ! force calculated here

       real(8),intent(inout):: pot_cstmnb   ! Custom NB potential

       real(8),intent(inout):: for_viri_cstmnb(:,:)
                                         ! virial(custom NB force) of each atom
       real(8),intent(inout):: pot_viri_cstmnb
                                   ! virial(custom NB potential) of each atom
       real(8),intent(inout):: pot_virit_cstmnb(:,:) ! virial tensor (custom NB)

       real(8),intent(inout):: pot_cstmnbex(:) ! Custom NB extra potential

       real(8),intent(inout):: for_viri_cstmnbex(:,:,:)
                                   ! extra virial(custom NB force) of each atom
       real(8),intent(inout):: pot_viri_cstmnbex(:)
                               ! extra virial(custom NB potential) of each atom
       real(8),intent(inout):: pot_virit_cstmnbex(:,:,:) 
                                          ! extra virial tensor (custom NB)

     end subroutine ent_calcstmnb

     subroutine ent_calcstmnb_hf(mts_flag,mts_cstmnb, &
          &                      xcel,ycel,zcel, &
          &                      npoly,npolytyp,npoly_mole,npoly_atom, &
          &                      nwater,nmatom,nmatyp,nmatomtyp, &
          &                      force,pot_cstmnb, &
          &                      for_viri_cstmnb,pot_viri_cstmnb, &
          &                      pot_virit_cstmnb, &
          &                      pot_cstmnb_atm,viricstmnbt_atm, &
          &                      ifhfvol,   &
          &                      nhfregion,hfzpos1,hfzpos2,   &
          &                      hftyp_atm,   &
          &                      molecom, &
          &                      heatfinterval, &
          &                      ncstmnbex, &
          &                      pot_cstmnbex, &
          &                      for_viri_cstmnbex,pot_viri_cstmnbex, &
          &                      pot_virit_cstmnbex, &
          &                      pot_cstmnbex_atm,viricstmnbext_atm, &
          &                      ifcalpreatom,ifcalpremole, &
          &                      current_step, &
          &                      nstep_short,istep_short)

       ! ARGUMENTS:
       !     INPUT
       integer,intent(in):: mts_flag        ! flag for MTS integration
                                       ! 1 long-range force mode
                                       ! 2 medium-range force mode
                                       ! 3 short-range force mode
       integer,intent(in):: mts_cstmnb      ! MTS flag for custom NB interaction

       integer,intent(in):: npoly           ! all number of poly
       integer,intent(in):: npolytyp        ! number of poly type
       integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
       integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
  
       integer,intent(in):: nwater          ! number of H2O molecules

       integer,intent(in):: nmatom          ! number of monatomic molecules
       integer,intent(in):: nmatyp          ! number of species of monatomic mole.
       integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

       real(8),intent(in):: xcel            ! x cell length[non-d]
       real(8),intent(in):: ycel            ! y cell length[non-d]
       real(8),intent(in):: zcel            ! z cell length[non-d]

       logical,intent(in):: ifhfvol    ! local volume-based or local surface-based
       integer,intent(in):: nhfregion  ! number of region to calculate heat flux
       real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                       ! z-position of region for heat flux

       integer,intent(in):: hftyp_atm(:)    ! atom- or mole-based heat flux cal. 
                                       !   for each atom

       real(8),intent(in):: molecom(:,:)    ! center of mass of molecule

       integer,intent(in):: heatfinterval   ! interval of heatf output

       logical,intent(in):: ifcalpremole    ! pressure calculation of molecule
       logical,intent(in):: ifcalpreatom    ! pressure calculation of atom

       integer,intent(in):: current_step    ! current time step
       integer,intent(in):: nstep_short     ! number of step for short force
       integer,intent(in):: istep_short

       integer,intent(in):: ncstmnbex       ! number of extra custom NB output

       !     OUTPUT
       real(8),intent(inout):: force(:,:)   ! force calculated here

       real(8),intent(inout):: pot_cstmnb   ! Custom NB potential

       real(8),intent(inout):: for_viri_cstmnb(:,:)
                                       ! virial(custom NB force) of each atom
       real(8),intent(inout):: pot_viri_cstmnb
                                   ! virial(custom NB potential) of each atom
       real(8),intent(inout):: pot_virit_cstmnb(:,:) ! virial tensor (custom NB)

       real(8),intent(inout):: pot_cstmnb_atm(:) 
                                       ! custom NB potential of each atom
       real(8),intent(inout):: viricstmnbt_atm(:,:,:,:)
                                       ! virial tensor of each atom (custom NB)

       real(8),intent(inout):: pot_cstmnbex(:) ! Custom NB extra potential

       real(8),intent(inout):: for_viri_cstmnbex(:,:,:)
                                   ! extra virial(custom NB force) of each atom
       real(8),intent(inout):: pot_viri_cstmnbex(:)
                               ! extra virial(custom NB potential) of each atom
       real(8),intent(inout):: pot_virit_cstmnbex(:,:,:) 
                                          ! extra virial tensor (custom NB)

       real(8),intent(inout):: pot_cstmnbex_atm(:,:)
                                       ! extra custom NB potential of each atom
       real(8),intent(inout):: viricstmnbext_atm(:,:,:,:,:)
                                 ! virial tensor of each atom (extra custom NB)

     end subroutine ent_calcstmnb_hf

#if defined(MPI)
     subroutine MPIsum_all_pot_for(mts_flag, &
          &                        mts_bond,mts_angl,mts_anglub, &
          &                        mts_tors,mts_torsrb,mts_torsim, &
          &                        mts_vdw,mts_ewr,mts_ewk, &
          &                        mts_vdw14,mts_elc14, &
          &                        mts_mor,mts_sh,mts_rfh,mts_dou,mts_cnpvw, &
          &                        mts_cstmnb, &
          &                        mts_posres, &
          &                        mts_potbias, &
          &                        force, &
          &                        pot_vdw,pot_elc,pot_ewk, &
          &                        pot_vdw14,pot_elc14, &
          &                        pot_bond,pot_angl,pot_anglub, &
          &                        pot_tors,pot_torsrb,pot_torsim, &
          &                        pot_mor, &
          &                        pot_sh, &
          &                        pot_rfh, &
          &                        pot_dou, &
          &                        pot_rpvw,pot_vw, &
          &                        pot_cstmnb, &
          &                        ncstmnbex,pot_cstmnbex, &
          &                        pot_posres, &
          &                        pot_pbias)

       ! ARGUMENT:
       !   INPUT
       integer,intent(in):: mts_flag         ! flag for MTS integration
                                        ! 1 long-range force mode
                                        ! 2 medium-range force mode
                                        ! 3 short-range force mode

       integer,intent(in):: mts_bond         ! MTS flag for bond
       integer,intent(in):: mts_angl         ! MTS flag for angle
       integer,intent(in):: mts_anglub       ! MTS flag for Urey-Bradley angle
       integer,intent(in):: mts_tors         ! MTS flag for torsion
       integer,intent(in):: mts_torsrb       ! MTS flag for torsionrb
       integer,intent(in):: mts_torsim       ! MTS flag for torsionim
       integer,intent(in):: mts_vdw          ! MTS flag for vdw interaction
       integer,intent(in):: mts_ewr          ! MTS flag for ewald real(=vdw)
       integer,intent(in):: mts_ewk          ! MTS flag for ewald wave
       integer,intent(in):: mts_vdw14        ! MTS flag for 14vdw
       integer,intent(in):: mts_elc14        ! MTS flag for 14elc
       integer,intent(in):: mts_mor          ! MTS flag for Morse interaction
       integer,intent(in):: mts_sh           ! MTS flag for SH interaction
       integer,intent(in):: mts_rfh          ! MTS flag for RFH interaction
       integer,intent(in):: mts_dou          ! MTS flag for DOU interaction
       integer,intent(in):: mts_cnpvw        ! MTS flag for CNP_VW
       integer,intent(in):: mts_cstmnb       ! MTS flag for custom NB interaction
       integer,intent(in):: mts_posres       ! MTS flag for position restraint
       integer,intent(in):: mts_potbias      ! MTS flag for bias potential

       integer,intent(in):: ncstmnbex        ! number of extra custom NB output

       !   INPUT & OUTPUT
       real(8),intent(inout):: force(:,:)    ! force calculated here

       !   OUTPUT
       !---- potential 
       real(8),intent(out):: pot_vdw         ! vdw potential
       real(8),intent(out):: pot_elc         ! coulomb potential(ewald real)
       real(8),intent(out):: pot_ewk         ! coulomb potential(ewald wave)
       real(8),intent(out):: pot_vdw14       ! 1-4vdw potential
       real(8),intent(out):: pot_elc14       ! 1-4elc potential
       real(8),intent(out):: pot_bond        ! bond potential
       real(8),intent(out):: pot_angl        ! angle potential
       real(8),intent(out):: pot_anglub      ! Urey-Bradley angle potential
       real(8),intent(out):: pot_tors        ! torsion potential
       real(8),intent(out):: pot_torsrb      ! RBtorsion potential
       real(8),intent(out):: pot_torsim      ! improper torsion potential
       real(8),intent(out):: pot_mor         ! Morse potential
       real(8),intent(out):: pot_sh          ! SH potential
       real(8),intent(out):: pot_rfh         ! RFH potential
       real(8),intent(out):: pot_dou         ! DOU potential
       real(8),intent(out):: pot_rpvw        ! RP-VW interaction
       real(8),intent(out):: pot_vw          ! potential of constant force
       real(8),intent(out):: pot_cstmnb      ! Custom NB potential
       real(8),intent(out):: pot_cstmnbex(:) ! extra Custom NB potential
       real(8),intent(out):: pot_posres      ! position restraint potential
       real(8),intent(out):: pot_pbias       ! bias potential

     end subroutine MPIsum_all_pot_for

     subroutine MPIsum_all_viri(mts_flag, &
          &                     mts_bond,mts_angl,mts_anglub, &
          &                     mts_tors,mts_torsrb,mts_torsim, &
          &                     mts_vdw,mts_ewr,mts_ewk, &
          &                     mts_vdw14,mts_elc14, &
          &                     mts_mor,mts_sh,mts_rfh,mts_dou, &
          &                     mts_cstmnb, &
          &                     for_viri_coul,pot_viri_coul, &
          &                     for_viri_lj,pot_viri_lj, &
          &                     for_viri_mor,pot_viri_mor, &
          &                     for_viri_sh,pot_viri_sh, &
          &                     for_viri_rfh,pot_viri_rfh, &
          &                     for_viri_dou,pot_viri_dou, &
          &                     for_viri_cstmnb,pot_viri_cstmnb, &
          &                     for_viri_cstmnbex,pot_viri_cstmnbex, &
          &                     atm_viri_bond,atm_viri_angl, &
          &                     atm_viri_tors,atm_viri_14, &
          &                     pot_virit_coul,pot_virit_lj, &
          &                     pot_virit_mor, &
          &                     pot_virit_sh, &
          &                     pot_virit_rfh, &
          &                     pot_virit_dou, &
          &                     pot_virit_cstmnb, &
          &                     pot_virit_cstmnbex, &
          &                     atm_virit_bond,atm_virit_angl, &
          &                     atm_virit_tors,atm_virit_14, &
          &                     ifcalpremole,ifcalpreatom, &
          &                     ncstmnbex)

       ! ARGUMENT:
       !   INPUT
       integer,intent(in):: mts_flag         ! flag for MTS integration
                                        ! 1 long-range force mode    
                                        ! 2 medium-range force mode
                                        ! 3 short-range force mode

       integer,intent(in):: mts_bond         ! MTS flag for bond
       integer,intent(in):: mts_angl         ! MTS flag for angle
       integer,intent(in):: mts_anglub       ! MTS flag for Urey-Bradley angle
       integer,intent(in):: mts_tors         ! MTS flag for torsion
       integer,intent(in):: mts_torsrb       ! MTS flag for torsionrb
       integer,intent(in):: mts_torsim       ! MTS flag for torsionim
       integer,intent(in):: mts_vdw          ! MTS flag for vdw interaction
       integer,intent(in):: mts_ewr          ! MTS flag for ewald real(=vdw)
       integer,intent(in):: mts_ewk          ! MTS flag for ewald wave
       integer,intent(in):: mts_vdw14        ! MTS flag for 14vdw
       integer,intent(in):: mts_elc14        ! MTS flag for 14elc
       integer,intent(in):: mts_mor          ! MTS flag for Morse interaction
       integer,intent(in):: mts_sh           ! MTS flag for SH interaction
       integer,intent(in):: mts_rfh          ! MTS flag for RFH interaction
       integer,intent(in):: mts_dou          ! MTS flag for DOU interaction
       integer,intent(in):: mts_cstmnb       ! MTS flag for custom NB interaction

       !---- valiables for pressure

       logical,intent(in):: ifcalpremole     ! pressure calculation of molecule
       logical,intent(in):: ifcalpreatom     ! pressure calculation of atom

       !---- variables for extra custom NB potential
       integer,intent(in):: ncstmnbex       ! number of extra custom NB output

       !   OUTPUT
       !   for molecular pressure
       real(8),intent(out):: for_viri_coul(:,:)  ! virial(coulomb force) of each atom
       real(8),intent(out):: pot_viri_coul       ! virial(coulomb potential)
       real(8),intent(out):: pot_virit_coul(:,:) ! virial tensor (coulomb)

       real(8),intent(out):: for_viri_lj(:,:)  ! virial(L-J force) of each atom
       real(8),intent(out):: pot_viri_lj       ! virial(L-J potential)
       real(8),intent(out):: pot_virit_lj(:,:) ! virial tensor (L-J)

       real(8),intent(out):: for_viri_mor(:,:)  ! virial(Morse force) of each atom
       real(8),intent(out):: pot_viri_mor       ! virial(Morse potential)
       real(8),intent(out):: pot_virit_mor(:,:) ! virial tensor (Morse)

       real(8),intent(out):: for_viri_sh(:,:)   ! virial(SH force) of each atom
       real(8),intent(out):: pot_viri_sh        ! virial(SH potential)
       real(8),intent(out):: pot_virit_sh(:,:)  ! virial tensor (SH)

       real(8),intent(out):: for_viri_rfh(:,:)   ! virial(RFH force) of each atom
       real(8),intent(out):: pot_viri_rfh        ! virial(RFH potential)
       real(8),intent(out):: pot_virit_rfh(:,:)  ! virial tensor (RFH)

       real(8),intent(out):: for_viri_dou(:,:)  ! virial(DOU force) of each atom
       real(8),intent(out):: pot_viri_dou       ! virial(DOU potential)
       real(8),intent(out):: pot_virit_dou(:,:) ! virial tensor (DOU)

       real(8),intent(out):: for_viri_cstmnb(:,:)
                                        ! virial(custom NB force) of each atom
       real(8),intent(out):: pot_viri_cstmnb       ! virial(custom NB potential)
       real(8),intent(out):: pot_virit_cstmnb(:,:) ! virial tensor (custom NB)

       real(8),intent(out):: for_viri_cstmnbex(:,:,:)
                                   ! extra virial(custom NB force) of each atom
       real(8),intent(out):: pot_viri_cstmnbex(:)
                                            ! extra virial(custom NB potential)
       real(8),intent(out):: pot_virit_cstmnbex(:,:,:)
                                             ! extra virial tensor (custom NB)

       !   for atomic pressure
       real(8),intent(out):: atm_viri_bond       ! virial(bond potential)
       real(8),intent(out):: atm_virit_bond(:,:) ! virial tensor (bond potential)
       real(8),intent(out):: atm_viri_angl       ! virial(angle potential)
       real(8),intent(out):: atm_virit_angl(:,:) ! virial tensor (angle potential)
       real(8),intent(out):: atm_viri_tors       ! virial(torsion potential)
       real(8),intent(out):: atm_virit_tors(:,:) ! virial tensor (torsion potential)
       real(8),intent(out):: atm_viri_14         ! virial(1-4 force potential)
       real(8),intent(out):: atm_virit_14(:,:)   ! virial tensor (1-4 force potential)

     end subroutine MPIsum_all_viri

#if defined(HF)
     subroutine MPIsum_all_atm_heatf(mts_flag, &
          &                          mts_bond,mts_angl,mts_anglub, &
          &                          mts_tors,mts_torsrb,mts_torsim, &
          &                          mts_vdw,mts_ewr,mts_ewk, &
          &                          mts_vdw14,mts_elc14, &
          &                          mts_mor,mts_sh,mts_rfh,mts_dou, &
          &                          mts_cstmnb, &
          &                          nhfregion, &
          &                          pot_bo_atm, &
          &                          pot_an_atm, &
          &                          pot_to_atm, &
          &                          pot_14_atm, &
          &                          pot_elin_atm,pot_ljin_atm, &
          &                          pot_elc_atm, &
          &                          pot_vdw_atm, &
          &                          pot_mor_atm, &
          &                          pot_sh_atm, &
          &                          pot_rfh_atm, &
          &                          pot_dou_atm, &
          &                          pot_cstmnb_atm, &
          &                          pot_cstmnbex_atm, &
          &                          viribot_atm, &
          &                          viriant_atm, &
          &                          viritot_atm, &
          &                          viri14t_atm, &
          &                          virielint_atm,viriljint_atm, &
          &                          virielct_atm, &
          &                          viriljt_atm, &
          &                          virimort_atm, &
          &                          virisht_atm, &
          &                          virirfht_atm, &
          &                          viridout_atm, &
          &                          viricstmnbt_atm, &
          &                          viricstmnbext_atm, &
          &                          ncstmnbex)

       ! ARGUMENT:
       !   INPUT
       integer,intent(in):: mts_flag         ! flag for MTS integration
                                        ! 1 long-range force mode    
                                        ! 2 medium-range force mode
                                        ! 3 short-range force mode

       integer,intent(in):: mts_bond         ! MTS flag for bond
       integer,intent(in):: mts_angl         ! MTS flag for angle
       integer,intent(in):: mts_anglub       ! MTS flag for Urey-Bradley angle
       integer,intent(in):: mts_tors         ! MTS flag for torsion
       integer,intent(in):: mts_torsrb       ! MTS flag for torsionrb
       integer,intent(in):: mts_torsim       ! MTS flag for torsionim
       integer,intent(in):: mts_vdw          ! MTS flag for vdw interaction
       integer,intent(in):: mts_ewr          ! MTS flag for ewald real(=vdw)
       integer,intent(in):: mts_ewk          ! MTS flag for ewald wave
       integer,intent(in):: mts_vdw14        ! MTS flag for 14vdw
       integer,intent(in):: mts_elc14        ! MTS flag for 14elc
       integer,intent(in):: mts_mor          ! MTS flag for Morse interaction
       integer,intent(in):: mts_sh           ! MTS flag for SH interaction
       integer,intent(in):: mts_rfh          ! MTS flag for RFH interaction
       integer,intent(in):: mts_dou          ! MTS flag for DOU interaction
       integer,intent(in):: mts_cstmnb       ! MTS flag for custom NB interaction

       integer,intent(in):: nhfregion        ! number of region to calculate heat flux

       integer,intent(in):: ncstmnbex        ! number of extra custom NB output

       !   OUTPUT
       !   potential for each atom
       real(8),intent(out):: pot_bo_atm(:)   ! bond potential of each atom
       real(8),intent(out):: pot_an_atm(:)   ! angl potential of each atom
       real(8),intent(out):: pot_to_atm(:)   ! torsion potential of each atom
       real(8),intent(out):: pot_14_atm(:)   ! 1-4 potential of each atom
       real(8),intent(out):: pot_elin_atm(:) ! elc (intra) potential of each atom
       real(8),intent(out):: pot_ljin_atm(:) ! L-J (intra) potential of each atom

       real(8),intent(out):: pot_elc_atm(:)  ! Coulomb potential of each atom
       real(8),intent(out):: pot_vdw_atm(:)  ! VDW potential of each atom
       real(8),intent(out):: pot_mor_atm(:)  ! Morse potential of each atom
       real(8),intent(out):: pot_sh_atm(:)   ! SH potential of each atom
       real(8),intent(out):: pot_rfh_atm(:)  ! RFH potential of each atom
       real(8),intent(out):: pot_dou_atm(:)  ! DOU potential of each atom
       real(8),intent(out):: pot_cstmnb_atm(:) ! custom NB potential of each atom
       real(8),intent(out):: pot_cstmnbex_atm(:,:)
                                       ! extra custom NB potential of each atom

       !   virial tensor for each atom
       real(8),intent(out):: viribot_atm(:,:,:,:)
                                        ! virial tensor of each atom (bond)
       real(8),intent(out):: viriant_atm(:,:,:,:)
                                        ! virial tensor of each atom (angle)
       real(8),intent(out):: viritot_atm(:,:,:,:)
                                        ! virial tensor of each atom (torsion)
       real(8),intent(out):: viri14t_atm(:,:,:,:)
                                        ! virial tensor of each atom (1-4)
       real(8),intent(out):: virielint_atm(:,:,:,:)
                                        ! virial tensor of each atom (elc intra)
       real(8),intent(out):: viriljint_atm(:,:,:,:)
                                        ! virial tensor of each atom (L-J intra)

       real(8),intent(out):: virielct_atm(:,:,:,:)
                                        ! virial tensor of each atom (Coulomb)
       real(8),intent(out):: viriljt_atm(:,:,:,:)
                                        ! virial tensor of each atom (L-J)
       real(8),intent(out):: virimort_atm(:,:,:,:)
                                        ! virial tensor of each atom (Morse)
       real(8),intent(out):: virisht_atm(:,:,:,:)
                                        ! virial tensor of each atom (SH)
       real(8),intent(out):: virirfht_atm(:,:,:,:)
                                        ! virial tensor of each atom (RFH)
       real(8),intent(out):: viridout_atm(:,:,:,:)
                                        ! virial tensor of each atom (DOU)
       real(8),intent(out):: viricstmnbt_atm(:,:,:,:)
                                        ! virial tensor of each atom (custom NB)
       real(8),intent(out):: viricstmnbext_atm(:,:,:,:,:)
                                  ! virial tensor of each atom (extra custom NB)

     end subroutine MPIsum_all_atm_heatf
#endif
#endif

#if defined(HF)
     subroutine calvirihf(i,j,fi,fj,   &
          &               nbodycoeff,   &
          &               box,box_inv,   &
          &               ifhfvol,   &
          &               nhfregion,hfzpos1,hfzpos2,   &
          &               hftyp_atm,   &
          &               molecom,   &
          &               virihft_atm)

       ! ARGUMENTS:
       !   INPUT
       integer,intent(in):: i,j              ! particle number
       !  real(8),intent(in):: rij(:)           ! Ri - Rj
       real(8),intent(in):: fi(:),fj(:)      ! force to i and j
       real(8),intent(in):: nbodycoeff  ! n-body coefficient for virial term of heatf

       real(8),intent(in):: box(:)           ! BOX size
       real(8),intent(in):: box_inv(:)       ! inverse of BOX size

       logical,intent(in):: ifhfvol       ! local volume-based or local surface-based

       integer,intent(in):: nhfregion       ! number of region to calculate heat flux
       real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                        ! z-position of region for heat flux

       integer,intent(in):: hftyp_atm(:)     ! atom- or mole-based heat flux cal. 
                                        !   for each atom

       real(8),intent(in):: molecom(:,:)     ! center of mass of molecule

       !   OUTPUT
       real(8),intent(out):: virihft_atm(:,:,:,:)
                                        ! virial tensor of each atom (L-J)

     end subroutine calvirihf
#endif

  end interface

end module interface_interact
