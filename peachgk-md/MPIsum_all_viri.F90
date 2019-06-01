!*********************************
!*  MPIsum_all_viri.f90 Ver.1.8  *
!*      for peachgk_md.f         *
!*            by G.Kikugawa      *
!*********************************
! Time-stamp: <2015-07-18 00:09:00 gota>

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

  use md_global
  use mpi_global

  implicit none

!     subroutine to calculate force  

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

! LOCAL:
  integer:: i                           ! do loop index

! MPI local variable for message passing
  integer:: msg_length_vir              ! message length for virial
  real(8),allocatable:: vir_local(:)       ! local virial for MPI
  real(8),allocatable:: vir_global(:)      ! global virial for MPI

  integer:: msg_length_virt             ! message length for virial tensor
  real(8),allocatable:: virt_local(:,:,:)  ! local virial tensor for MPI
  real(8),allocatable:: virt_global(:,:,:) ! global virial tensor for MPI

  integer:: msg_length_fc         ! message length for correction force (virial)
  real(8),allocatable:: fc_local(:,:,:)
                                      ! local correction force (virial) for MPI
  real(8),allocatable:: fc_global(:,:,:)
                                      ! global correction force (virial) for MPI

  integer:: maxnintrct

!     +     +     +     +     +     +     +

!---- memory allocation

  maxnintrct = 7   ! fixed number of interations
  maxnintrct = maxnintrct + ncstmnbex

  allocate(vir_local(maxnintrct))
  allocate(vir_global(maxnintrct))

  allocate(virt_local(3,3,maxnintrct))
  allocate(virt_global(3,3,maxnintrct))

  allocate(fc_local(3,natom,maxnintrct))
  allocate(fc_global(3,natom,maxnintrct))

! --- INITIALIZATION ---

  msg_length_vir  = 0
  msg_length_virt = 0
  msg_length_fc   = 0

! abbreviate the initialization of local parameters

!---- MPI for intra-molecular virial

  if (ifcalpreatom) then

!    store to local variables
     msg_length_vir = 0

     msg_length_vir = msg_length_vir + 1
     vir_local(msg_length_vir) = atm_viri_bond

     msg_length_vir = msg_length_vir + 1
     vir_local(msg_length_vir) = atm_viri_angl

     msg_length_vir = msg_length_vir + 1
     vir_local(msg_length_vir) = atm_viri_tors

     msg_length_vir = msg_length_vir + 1
     vir_local(msg_length_vir) = atm_viri_14

     msg_length_virt = 0

     msg_length_virt = msg_length_virt + 1
!     do i = 1,3
!        virt_local(1:3,i,msg_length_virt) = atm_virit_bond(1:3,i)
!     end do
     virt_local(1:3,1:3,msg_length_virt) = atm_virit_bond(1:3,1:3)

     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,msg_length_virt) = atm_virit_angl(1:3,1:3)

     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,msg_length_virt) = atm_virit_tors(1:3,1:3)

     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,msg_length_virt) = atm_virit_14(1:3,1:3)

!    MPI allreduce for atomic virial
     call MPI_ALLREDUCE(vir_local,vir_global,msg_length_vir,   &
          &             MPI_DOUBLE_PRECISION,   &
          &             MPI_SUM,MPI_COMM_WORLD,ierror)

!    MPI allreduce for atomic virial tensor
     call MPI_ALLREDUCE(virt_local,virt_global,msg_length_virt*3*3,   &
          &             MPI_DOUBLE_PRECISION,   &
          &             MPI_SUM,MPI_COMM_WORLD,ierror)

!    restore from global variables
     msg_length_vir = 0

     msg_length_vir = msg_length_vir + 1
     atm_viri_bond = vir_global(msg_length_vir)

     msg_length_vir = msg_length_vir + 1
     atm_viri_angl = vir_global(msg_length_vir)

     msg_length_vir = msg_length_vir + 1
     atm_viri_tors = vir_global(msg_length_vir)

     msg_length_vir = msg_length_vir + 1
     atm_viri_14 = vir_global(msg_length_vir)

     msg_length_virt = 0

     msg_length_virt = msg_length_virt + 1
     atm_virit_bond(1:3,1:3) = virt_global(1:3,1:3,msg_length_virt)

     msg_length_virt = msg_length_virt + 1
     atm_virit_angl(1:3,1:3) = virt_global(1:3,1:3,msg_length_virt)

     msg_length_virt = msg_length_virt + 1
     atm_virit_tors(1:3,1:3) = virt_global(1:3,1:3,msg_length_virt)

     msg_length_virt = msg_length_virt + 1
     atm_virit_14(1:3,1:3) = virt_global(1:3,1:3,msg_length_virt)

  end if

!---- MPI for inter-molecular virial

! store to local variables
  msg_length_vir = 0

  msg_length_vir = msg_length_vir + 1
  vir_local(msg_length_vir) = pot_viri_coul

  msg_length_vir = msg_length_vir + 1
  vir_local(msg_length_vir) = pot_viri_lj

  msg_length_vir = msg_length_vir + 1
  vir_local(msg_length_vir) = pot_viri_mor

  msg_length_vir = msg_length_vir + 1
  vir_local(msg_length_vir) = pot_viri_sh

  msg_length_vir = msg_length_vir + 1
  vir_local(msg_length_vir) = pot_viri_rfh

  msg_length_vir = msg_length_vir + 1
  vir_local(msg_length_vir) = pot_viri_dou

  msg_length_vir = msg_length_vir + 1
  vir_local(msg_length_vir) = pot_viri_cstmnb

  do i = 1, ncstmnbex
     msg_length_vir = msg_length_vir + 1
     vir_local(msg_length_vir) = pot_viri_cstmnbex(i)
  end do

  msg_length_virt = 0

  msg_length_virt = msg_length_virt + 1
  virt_local(1:3,1:3,msg_length_virt) = pot_virit_coul(1:3,1:3)

  msg_length_virt = msg_length_virt + 1
  virt_local(1:3,1:3,msg_length_virt) = pot_virit_lj(1:3,1:3)

  msg_length_virt = msg_length_virt + 1
  virt_local(1:3,1:3,msg_length_virt) = pot_virit_mor(1:3,1:3)

  msg_length_virt = msg_length_virt + 1
  virt_local(1:3,1:3,msg_length_virt) = pot_virit_sh(1:3,1:3)

  msg_length_virt = msg_length_virt + 1
  virt_local(1:3,1:3,msg_length_virt) = pot_virit_rfh(1:3,1:3)

  msg_length_virt = msg_length_virt + 1
  virt_local(1:3,1:3,msg_length_virt) = pot_virit_dou(1:3,1:3)

  msg_length_virt = msg_length_virt + 1
  virt_local(1:3,1:3,msg_length_virt) = pot_virit_cstmnb(1:3,1:3)

  do i = 1, ncstmnbex
     msg_length_virt = msg_length_virt + 1
     virt_local(1:3,1:3,msg_length_virt) = pot_virit_cstmnbex(1:3,1:3,i)
  end do

! MPI allreduce for molecular virial
  call MPI_ALLREDUCE(vir_local,vir_global,msg_length_vir,   &
       &             MPI_DOUBLE_PRECISION,   &
       &             MPI_SUM,MPI_COMM_WORLD,ierror)

! MPI allreduce for molecular virial tensor
  call MPI_ALLREDUCE(virt_local,virt_global,msg_length_virt*3*3,   &
       &             MPI_DOUBLE_PRECISION,   &
       &             MPI_SUM,MPI_COMM_WORLD,ierror)

! restore from global variables
  msg_length_vir = 0

  msg_length_vir = msg_length_vir + 1
  pot_viri_coul = vir_global(msg_length_vir)

  msg_length_vir = msg_length_vir + 1
  pot_viri_lj = vir_global(msg_length_vir)

  msg_length_vir = msg_length_vir + 1
  pot_viri_mor = vir_global(msg_length_vir)

  msg_length_vir = msg_length_vir + 1
  pot_viri_sh = vir_global(msg_length_vir)

  msg_length_vir = msg_length_vir + 1
  pot_viri_rfh = vir_global(msg_length_vir)

  msg_length_vir = msg_length_vir + 1
  pot_viri_dou = vir_global(msg_length_vir)

  msg_length_vir = msg_length_vir + 1
  pot_viri_cstmnb = vir_global(msg_length_vir)

  do i = 1, ncstmnbex
     msg_length_vir = msg_length_vir + 1
     pot_viri_cstmnbex(i) = vir_global(msg_length_vir)
  end do

  msg_length_virt = 0

  msg_length_virt = msg_length_virt + 1
  pot_virit_coul(1:3,1:3) = virt_global(1:3,1:3,msg_length_virt)

  msg_length_virt = msg_length_virt + 1
  pot_virit_lj(1:3,1:3) = virt_global(1:3,1:3,msg_length_virt)

  msg_length_virt = msg_length_virt + 1
  pot_virit_mor(1:3,1:3) = virt_global(1:3,1:3,msg_length_virt)

  msg_length_virt = msg_length_virt + 1
  pot_virit_sh(1:3,1:3) = virt_global(1:3,1:3,msg_length_virt)

  msg_length_virt = msg_length_virt + 1
  pot_virit_rfh(1:3,1:3) = virt_global(1:3,1:3,msg_length_virt)

  msg_length_virt = msg_length_virt + 1
  pot_virit_dou(1:3,1:3) = virt_global(1:3,1:3,msg_length_virt)

  msg_length_virt = msg_length_virt + 1
  pot_virit_cstmnb(1:3,1:3) = virt_global(1:3,1:3,msg_length_virt)

  do i = 1, ncstmnbex
     msg_length_virt = msg_length_virt + 1
     pot_virit_cstmnbex(1:3,1:3,i) = virt_global(1:3,1:3,msg_length_virt)
  end do

!---- MPI for correction force (molecular virial)

  if (ifcalpremole) then

!    store to local variables
     msg_length_fc = 0

     if ((mts_ewr == mts_flag) .and. (mts_vdw == mts_flag)) then
        msg_length_fc = msg_length_fc + 1
        fc_local(1:3,1:natom,msg_length_fc) = for_viri_coul(1:3,1:natom)
     end if

     if (mts_vdw == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        fc_local(1:3,1:natom,msg_length_fc) = for_viri_lj(1:3,1:natom)
     end if

     if (mts_mor == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        fc_local(1:3,1:natom,msg_length_fc) = for_viri_mor(1:3,1:natom)
     end if

     if (mts_sh == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        fc_local(1:3,1:natom,msg_length_fc) = for_viri_sh(1:3,1:natom)
     end if

     if (mts_rfh == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        fc_local(1:3,1:natom,msg_length_fc) = for_viri_rfh(1:3,1:natom)
     end if

     if (mts_dou == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        fc_local(1:3,1:natom,msg_length_fc) = for_viri_dou(1:3,1:natom)
     end if

     if (mts_cstmnb == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        fc_local(1:3,1:natom,msg_length_fc) = for_viri_cstmnb(1:3,1:natom)
     end if

     ! mts_flag is ignored for cstmnbex
     do i = 1, ncstmnbex
        msg_length_fc = msg_length_fc + 1
        fc_local(1:3,1:natom,msg_length_fc) = for_viri_cstmnbex(1:3,1:natom,i)
     end do

!    MPI allreduce for correction force (virial)
     call MPI_ALLREDUCE(fc_local,fc_global,msg_length_fc*natom*3,   &
          &             MPI_DOUBLE_PRECISION,   &
          &             MPI_SUM,MPI_COMM_WORLD,ierror)

!    restore from global variables
     msg_length_fc = 0

     if ((mts_ewr == mts_flag) .and. (mts_vdw == mts_flag)) then
        msg_length_fc = msg_length_fc + 1
        for_viri_coul(1:3,1:natom) = fc_global(1:3,1:natom,msg_length_fc)
     end if

     if (mts_vdw == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        for_viri_lj(1:3,1:natom) = fc_global(1:3,1:natom,msg_length_fc)
     end if

     if (mts_mor == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        for_viri_mor(1:3,1:natom) = fc_global(1:3,1:natom,msg_length_fc)
     end if

     if (mts_sh == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        for_viri_sh(1:3,1:natom) = fc_global(1:3,1:natom,msg_length_fc)
     end if

     if (mts_rfh == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        for_viri_rfh(1:3,1:natom) = fc_global(1:3,1:natom,msg_length_fc)
     end if

     if (mts_dou == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        for_viri_dou(1:3,1:natom) = fc_global(1:3,1:natom,msg_length_fc)
     end if

     if (mts_cstmnb == mts_flag) then
        msg_length_fc = msg_length_fc + 1
        for_viri_cstmnb(1:3,1:natom) = fc_global(1:3,1:natom,msg_length_fc)
     end if

     ! mts_flag is ignored for cstmnbex
     do i = 1, ncstmnbex
        msg_length_fc = msg_length_fc + 1
        for_viri_cstmnbex(1:3,1:natom,i) = fc_global(1:3,1:natom,msg_length_fc)
     end do

  end if

!---- memory release
  deallocate(vir_local)
  deallocate(vir_global)

  deallocate(virt_local)
  deallocate(virt_global)

  deallocate(fc_local)
  deallocate(fc_global)

!     +     +     +     +     +     +     +

end subroutine MPIsum_all_viri
