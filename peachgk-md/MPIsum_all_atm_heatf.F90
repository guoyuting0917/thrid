!**************************************
!*  MPIsum_all_atm_heatf.f90 Ver.1.9  *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
! Time-stamp: <2015-07-18 00:13:52 gota>

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

  use md_global
  use mpi_global

  implicit none

! subroutine to calculate force

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

  real(8),intent(out):: pot_elc_atm(:)    ! Coulomb potential of each atom
  real(8),intent(out):: pot_vdw_atm(:)    ! VDW potential of each atom
  real(8),intent(out):: pot_mor_atm(:)    ! Morse potential of each atom
  real(8),intent(out):: pot_sh_atm(:)     ! SH potential of each atom
  real(8),intent(out):: pot_rfh_atm(:)    ! RFH potential of each atom
  real(8),intent(out):: pot_dou_atm(:)    ! DOU potential of each atom
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

! LOCAL:
  integer:: i                            ! do loop index
  integer:: ihfr

! MPI local variable for message passing
  integer:: msg_length_pot           ! message length for potential of each atom
  real(8),allocatable:: pot_local(:,:)  ! local potential for MPI
  real(8),allocatable:: pot_global(:,:) ! global potential for MPI

  integer:: msg_length_virt      ! message length for virial tensor of each atom
  real(8),allocatable:: virt_local(:,:,:,:,:)
                                         ! local virial tensor for MPI
  real(8),allocatable:: virt_global(:,:,:,:,:)
                                         ! global virial tensor for MPI

  integer:: maxnintrct

!     +     +     +     +     +     +     +

!---- memory allocation

  maxnintrct = 13   ! fixed number of interations
  maxnintrct = maxnintrct + ncstmnbex

  allocate(pot_local(natom,maxnintrct))
  allocate(pot_global(natom,maxnintrct))

  allocate(virt_local(3,3,natom,nhfregion,maxnintrct))
  allocate(virt_global(3,3,natom,nhfregion,maxnintrct))

! --- INITIALIZATION ---

!  msg_length_pot  = 0
!  msg_length_virt = 0

! abbreviate the initialization of local parameters

!---- MPI for virial tensor of each atom

! store to local variables
  msg_length_virt = 0

  if (mts_bond == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                             viribot_atm(1:3,1:3,1:natom,1:nhfregion)
  end if

  if ((mts_angl == mts_flag) .or. (mts_anglub == mts_flag)) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                             viriant_atm(1:3,1:3,1:natom,1:nhfregion)
  end if

  if ((mts_tors == mts_flag) .or. (mts_torsrb == mts_flag) &
       &                     .or. (mts_torsim == mts_flag)) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                             viritot_atm(1:3,1:3,1:natom,1:nhfregion)
  end if

  if ((mts_vdw14 == mts_flag) .and. (mts_elc14 == mts_flag)) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                             viri14t_atm(1:3,1:3,1:natom,1:nhfregion)
  end if

#if defined(_DO_NOT_USE_THIS_)
#else
  if ((mts_ewr == mts_flag) .and. (mts_vdw == mts_flag)) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                           virielint_atm(1:3,1:3,1:natom,1:nhfregion)
  end if

  if ((mts_ewr == mts_flag) .and. (mts_vdw == mts_flag)) then      
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                            virielct_atm(1:3,1:3,1:natom,1:nhfregion)
  end if
#endif

  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                           viriljint_atm(1:3,1:3,1:natom,1:nhfregion)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                             viriljt_atm(1:3,1:3,1:natom,1:nhfregion)
  end if

                                                                                                     !!!!!!!! mts_flag viri
  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                             viriljt_atm_ar11(1:3,1:3,1:natom,1:nhfregion)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                             viriljt_atm_ar12(1:3,1:3,1:natom,1:nhfregion)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                             viriljt_atm_ar22(1:3,1:3,1:natom,1:nhfregion)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                             viriljt_atm_ar1pt(1:3,1:3,1:natom,1:nhfregion)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                             viriljt_atm_ar2pt(1:3,1:3,1:natom,1:nhfregion)
  end if

                                                                                                     !!!!!!!! mts_flag viri

  if (mts_mor == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                            virimort_atm(1:3,1:3,1:natom,1:nhfregion)
  endif

  if (mts_sh == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                             virisht_atm(1:3,1:3,1:natom,1:nhfregion)
  end if

  if (mts_rfh == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                            virirfht_atm(1:3,1:3,1:natom,1:nhfregion)
  end if

  if (mts_dou == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                            viridout_atm(1:3,1:3,1:natom,1:nhfregion)
  end if

  if (mts_cstmnb == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                         viricstmnbt_atm(1:3,1:3,1:natom,1:nhfregion)
  end if

  ! mts_flag is ignored for cstmnbex
  do i = 1, ncstmnbex
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,1:natom,1:nhfregion,msg_length_virt) =   &
          &                   viricstmnbext_atm(1:3,1:3,1:natom,1:nhfregion,i)
  end do

! MPI allreduce for virial tensor of each atom
#if defined(_DO_NOT_USE_THIS)
  call MPI_ALLREDUCE(virt_local,virt_global,   &
       &             msg_length_virt*3*3*natom*nhfregion,   &
       &             MPI_DOUBLE_PRECISION,   &
       &             MPI_SUM,MPI_COMM_WORLD,ierror)
#else
  call MPI_REDUCE(virt_local,virt_global,   &
       &          msg_length_virt*3*3*natom*nhfregion,   &
       &          MPI_DOUBLE_PRECISION,   &
       &          MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

! restore from global variables
  msg_length_virt = 0

  if (mts_bond == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     viribot_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if ((mts_angl == mts_flag) .or. (mts_anglub == mts_flag)) then
     msg_length_virt = msg_length_virt + 1
     viriant_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if ((mts_tors == mts_flag) .or. (mts_torsrb == mts_flag) &
       &                     .or. (mts_torsim == mts_flag)) then
     msg_length_virt = msg_length_virt + 1
     viritot_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if ((mts_vdw14 == mts_flag) .and. (mts_elc14 == mts_flag)) then
     msg_length_virt = msg_length_virt + 1
     viri14t_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

#if defined(_DO_NOT_USE_THIS_)
#else
  if ((mts_ewr == mts_flag) .and. (mts_vdw == mts_flag)) then
     msg_length_virt = msg_length_virt + 1
     virielint_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if ((mts_ewr == mts_flag) .and. (mts_vdw == mts_flag)) then
     msg_length_virt = msg_length_virt + 1
     virielct_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if
#endif

  if (mts_vdw == mts_flag) then      
     msg_length_virt = msg_length_virt + 1
     viriljint_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     viriljt_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if
                                                                                                     !!!!!!!! mts_flag viri
  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     viriljt_atm_ar11(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     viriljt_atm_ar12(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     viriljt_atm_ar22(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     viriljt_atm_ar1pt(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     viriljt_atm_ar2pt(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if
                                                                                                     !!!!!!!! mts_flag viri

  if (mts_mor == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virimort_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if (mts_sh == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virisht_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if
      
  if (mts_rfh == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     virirfht_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if (mts_dou == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     viridout_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  if (mts_cstmnb == mts_flag) then
     msg_length_virt = msg_length_virt + 1
     viricstmnbt_atm(1:3,1:3,1:natom,1:nhfregion) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end if

  ! mts_flag is ignored for cstmnbex
  do i = 1, ncstmnbex
     msg_length_virt = msg_length_virt + 1
     viricstmnbext_atm(1:3,1:3,1:natom,1:nhfregion,i) =   &
          &             virt_global(1:3,1:3,1:natom,1:nhfregion,msg_length_virt)
  end do

!---- MPI for potential of each atom

! store to local variables
  msg_length_pot = 0

  if (mts_bond == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_bo_atm(1:natom)
  end if

  if ((mts_angl == mts_flag) .or. (mts_anglub == mts_flag)) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_an_atm(1:natom)
  end if

  if ((mts_tors == mts_flag) .or. (mts_torsrb == mts_flag) &
       &                     .or. (mts_torsim == mts_flag)) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_to_atm(1:natom)
  end if

  if ((mts_vdw14 == mts_flag) .and. (mts_elc14 == mts_flag)) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_14_atm(1:natom)
  end if

#if defined(_DO_NOT_USE_THIS)
#else
  if ((mts_ewr == mts_flag) .and. (mts_vdw == mts_flag)) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_elin_atm(1:natom)
  end if

  if ((mts_ewr == mts_flag) .and. (mts_vdw == mts_flag)) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_elc_atm(1:natom)
  end if
#endif

  if (mts_vdw == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_ljin_atm(1:natom)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_vdw_atm(1:natom)
  end if
                                                                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (mts_vdw == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_vdw_atm_ar11(1:natom)
  end if
  if (mts_vdw == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_vdw_atm_ar12(1:natom)
  end if
  if (mts_vdw == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_vdw_atm_ar22(1:natom)
  end if
                                                                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                    
  if (mts_mor == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_mor_atm(1:natom)
  end if

  if (mts_sh == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_sh_atm(1:natom)
  end if

  if (mts_rfh == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_rfh_atm(1:natom)
  end if

  if (mts_dou == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_dou_atm(1:natom)
  end if

  if (mts_cstmnb == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_cstmnb_atm(1:natom)
  end if

  ! mts_flag is ignored for cstmnbex
  do i = 1, ncstmnbex
     msg_length_pot = msg_length_pot + 1
     pot_local(1:natom,msg_length_pot) = pot_cstmnbex_atm(1:natom,i)
  end do

! MPI allreduce for potential of each atom
#if defined(_DO_NOT_USE_THIS)
  call MPI_ALLREDUCE(pot_local,pot_global,msg_length_pot*natom,   &
       &             MPI_DOUBLE_PRECISION,   &
       &             MPI_SUM,MPI_COMM_WORLD,ierror)
#else
  call MPI_REDUCE(pot_local,pot_global,msg_length_pot*natom,   &
       &          MPI_DOUBLE_PRECISION,   &
       &          MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

! restore from global variables
  msg_length_pot = 0

  if (mts_bond == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_bo_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  if ((mts_angl == mts_flag) .or. (mts_anglub == mts_flag)) then
     msg_length_pot = msg_length_pot + 1
     pot_an_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  if ((mts_tors == mts_flag) .or. (mts_torsrb == mts_flag) &
       &                     .or. (mts_torsim == mts_flag)) then
     msg_length_pot = msg_length_pot + 1
     pot_to_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  if ((mts_vdw14 == mts_flag) .and. (mts_elc14 == mts_flag)) then
     msg_length_pot = msg_length_pot + 1
     pot_14_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

#if defined(_DO_NOT_USE_THIS_)
#else
  if ((mts_ewr == mts_flag) .and. (mts_vdw == mts_flag)) then
     msg_length_pot = msg_length_pot + 1
     pot_elin_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  if ((mts_ewr == mts_flag) .and. (mts_vdw == mts_flag)) then
     msg_length_pot = msg_length_pot + 1
     pot_elc_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if
#endif

  if (mts_vdw == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_ljin_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_vdw_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if
                                                                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pot flag
  if (mts_vdw == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_vdw_atm_ar11(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_vdw_atm_ar12(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  if (mts_vdw == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_vdw_atm_ar22(1:natom) = pot_global(1:natom,msg_length_pot)
  end if
                                                                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pot flag
  if (mts_mor == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_mor_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  if (mts_sh == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_sh_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  if (mts_rfh == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_rfh_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  if (mts_dou == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_dou_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  if (mts_cstmnb == mts_flag) then
     msg_length_pot = msg_length_pot + 1
     pot_cstmnb_atm(1:natom) = pot_global(1:natom,msg_length_pot)
  end if

  ! mts_flag is ignored for cstmnbex
  do i = 1, ncstmnbex
     msg_length_pot = msg_length_pot + 1
     pot_cstmnbex_atm(1:natom,i) = pot_global(1:natom,msg_length_pot)
  end do

!---- memory release
  deallocate(pot_local)
  deallocate(pot_global)

  deallocate(virt_local)
  deallocate(virt_global)

!     +     +     +     +     +     +     +

end subroutine MPIsum_all_atm_heatf
