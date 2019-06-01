!**************************************
!*  output_mf.f90 Ver.1.5             *
!*      for peachgk_md.f              *
!*            by N.Yamamoto           *
!*            modified by G.Kikugawa  *
!**************************************
! Time-stamp: <2015-08-28 13:42:22 gota>

subroutine outmtf(oumtf, &
     &            current_step, &
     &            xref,eref,pref,timeref, &
     &            ifhfvol, &
     &            nhfregion,hfzpos1,hfzpos2, &
     &            mfkin_atm,mfkinex_atm, &
     &            mfviribo_atm,mfvirian_atm, &
     &            mfvirito_atm,mfviri14_atm, &
     &            mfvirielin_atm,mfviriljin_atm, &
     &            mfvirielc_atm,mfvirilj_atm, &
     &            mfvirimor_atm,mfvirish_atm, &
     &            mfvirirfh_atm,mfviridou_atm, &
     &            mfviricstmnb_atm, &
     &            ncstmnbex,mfviricstmnbex_atm, &
     &            mftot_atm, &
     &            mtfoutdir)

  use interface_tools

  use md_global

  implicit none

!     subroutine for outputting momentum flux data
!!! Compile macro "_HF_ALL_DIR" is meaningless in this routine

! ARGUMENTS:
!     INPUT
  integer,intent(in):: oumtf           ! output unit for outmtf momentum flux data

  integer,intent(in):: current_step    ! current time step

!---- base value for non-dimensionalize
  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: eref             ! energy base value [J]
  real(8),intent(in):: pref             ! pressure base value [J]

  real(8),intent(in):: timeref          ! time base value [sec]

!---- variables for calculation of transport flux
  logical,intent(in):: ifhfvol        ! local volume-based or local surface-based

  integer,intent(in):: nhfregion       ! number of region to calculate heat flux
  real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                ! z-position of region for heat flux

!---- variables for heat flux
  real(8),intent(in):: mfkin_atm(:,:,:)     ! momentum flux of kinetic part
  real(8),intent(in):: mfkinex_atm(:,:,:,:)
                              ! momentum flux of kinetic part (extra custom NB)

  real(8),intent(in):: mfviribo_atm(:,:,:)   ! momentum flux of virial part (bond)
  real(8),intent(in):: mfvirian_atm(:,:,:)   ! momentum flux of virial part (angle)
  real(8),intent(in):: mfvirito_atm(:,:,:)
                     ! momentum flux of virial part (tors)
  real(8),intent(in):: mfviri14_atm(:,:,:)
                     ! momentum flux of virial part (1-4)
  real(8),intent(in):: mfvirielin_atm(:,:,:)
                     ! momentum flux of virial part (elc intra)
  real(8),intent(in):: mfviriljin_atm(:,:,:)
                     ! momentum flux of virial part (L-J intra)

  real(8),intent(in):: mfvirielc_atm(:,:,:)
                     ! momentum flux of virial part (Coulomb)
  real(8),intent(in):: mfvirilj_atm(:,:,:)
                     ! momentum flux of virial part (L-J)
  real(8),intent(in):: mfvirimor_atm(:,:,:)
                     ! momentum flux of virial part (Morse)
  real(8),intent(in):: mfvirish_atm(:,:,:)
                     ! momemtum flux of virial part (SH)
  real(8),intent(in):: mfvirirfh_atm(:,:,:)
                     ! momentum flux of virial part (RFH)
  real(8),intent(in):: mfviridou_atm(:,:,:)
                     ! momentum flux of virial part (DOU)
  real(8),intent(in):: mfviricstmnb_atm(:,:,:)
                     ! momentum flux of virial part (custom NB)
  real(8),intent(in):: mfviricstmnbex_atm(:,:,:,:)
                     ! momentum flux of virial part (extra custom NB)

  real(8),intent(in):: mftot_atm(:,:,:) ! total momentum flux

  integer,intent(in):: ncstmnbex        ! number of extra custom NB output

! variables for outputting momentum flux
  character(3),intent(in):: mtfoutdir ! direction to output data of momentum

! LOCAL:
  character(80):: sub_tmp
  integer:: sub_len

  character(80):: sub_mftot_atm(3,nhfregion)
  character(80):: sub_mfkin_atm(3,nhfregion)
  character(80):: sub_mfkinex_atm(3,nhfregion,1:ncstmnbex)
  character(80):: sub_mfviribo_atm(3,nhfregion)
  character(80):: sub_mfvirian_atm(3,nhfregion)
  character(80):: sub_mfvirito_atm(3,nhfregion)
  character(80):: sub_mfviri14_atm(3,nhfregion)
  character(80):: sub_mfvirielin_atm(3,nhfregion)
  character(80):: sub_mfviriljin_atm(3,nhfregion)
  character(80):: sub_mfvirielc_atm(3,nhfregion)
  character(80):: sub_mfvirilj_atm(3,nhfregion)
  character(80):: sub_mfvirimor_atm(3,nhfregion)
  character(80):: sub_mfvirish_atm(3,nhfregion)
  character(80):: sub_mfvirirfh_atm(3,nhfregion)
  character(80):: sub_mfviridou_atm(3,nhfregion)
  character(80):: sub_mfviricstmnb_atm(3,nhfregion)
  character(80):: sub_mfviricstmnbex_atm(3,nhfregion,1:ncstmnbex)

  integer:: n,i
  integer:: nex

  character(25):: fmt
  integer:: nmffield

!     +     +     +     +     +     +     +     +

!---- some preparation

  !- count number of output fields used below
  nmffield = 15                ! fixed number of fields (in one direction)
  nmffield = (nmffield + ncstmnbex*2) * nhfregion

!     --- Output heat flux data ---

  if (current_step == 1) then

!        make subject
     sub_mftot_atm(1:3,1:nhfregion) = ' '
     sub_mfkin_atm(1:3,1:nhfregion) = ' '
     sub_mfkinex_atm(1:3,1:nhfregion,1:ncstmnbex) = ' '
     sub_mfviribo_atm(1:3,1:nhfregion) = ' '
     sub_mfvirian_atm(1:3,1:nhfregion) = ' '
     sub_mfvirito_atm(1:3,1:nhfregion) = ' '
     sub_mfviri14_atm(1:3,1:nhfregion) = ' '
     sub_mfvirielin_atm(1:3,1:nhfregion) = ' '
     sub_mfviriljin_atm(1:3,1:nhfregion) = ' '
     sub_mfvirielc_atm(1:3,1:nhfregion) = ' '
     sub_mfvirilj_atm(1:3,1:nhfregion) = ' '
     sub_mfvirimor_atm(1:3,1:nhfregion) = ' '
     sub_mfvirish_atm(1:3,1:nhfregion) = ' '
     sub_mfvirirfh_atm(1:3,1:nhfregion) = ' '
     sub_mfviridou_atm(1:3,1:nhfregion) = ' '
     sub_mfviricstmnb_atm(1:3,1:nhfregion) = ' '
     sub_mfviricstmnbex_atm(1:3,1:nhfregion,1:ncstmnbex) = ' '

     do i = 1, nhfregion
        write(sub_tmp,*) 'mftot_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mftot_atm(1,i))
        write(sub_tmp,*) 'mftot_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mftot_atm(2,i))
        write(sub_tmp,*) 'mftot_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mftot_atm(3,i))
        write(sub_tmp,*) 'mfkin_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfkin_atm(1,i))
        write(sub_tmp,*) 'mfkin_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfkin_atm(2,i))
        write(sub_tmp,*) 'mfkin_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfkin_atm(3,i))
        do nex = 1, ncstmnbex
           write(sub_tmp,*) 'mfkinex_atm_',nex,'_xz_',i
           call excl_sp(sub_tmp,sub_len,sub_mfkinex_atm(1,i,nex))
           write(sub_tmp,*) 'mfkinex_atm_',nex,'_yz_',i
           call excl_sp(sub_tmp,sub_len,sub_mfkinex_atm(2,i,nex))
           write(sub_tmp,*) 'mfkinex_atm_',nex,'_zz_',i
           call excl_sp(sub_tmp,sub_len,sub_mfkinex_atm(3,i,nex))
        end do

        write(sub_tmp,*) 'mfviribo_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviribo_atm(1,i))
        write(sub_tmp,*) 'mfviribo_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviribo_atm(2,i))
        write(sub_tmp,*) 'mfviribo_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviribo_atm(3,i))

        write(sub_tmp,*) 'mfvirian_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirian_atm(1,i))
        write(sub_tmp,*) 'mfvirian_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirian_atm(2,i))
        write(sub_tmp,*) 'mfvirian_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirian_atm(3,i))

        write(sub_tmp,*) 'mfvirito_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirito_atm(1,i))
        write(sub_tmp,*) 'mfvirito_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirito_atm(2,i))
        write(sub_tmp,*) 'mfvirito_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirito_atm(3,i))

        write(sub_tmp,*) 'mfviri14_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviri14_atm(1,i))
        write(sub_tmp,*) 'mfviri14_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviri14_atm(2,i))
        write(sub_tmp,*) 'mfviri14_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviri14_atm(3,i))

        write(sub_tmp,*) 'mfvirielin_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirielin_atm(1,i))
        write(sub_tmp,*) 'mfvirielin_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirielin_atm(2,i))
        write(sub_tmp,*) 'mfvirielin_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirielin_atm(3,i))

        write(sub_tmp,*) 'mfviriljin_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviriljin_atm(1,i))
        write(sub_tmp,*) 'mfviriljin_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviriljin_atm(2,i))
        write(sub_tmp,*) 'mfviriljin_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviriljin_atm(3,i))

        write(sub_tmp,*) 'mfvirielc_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirielc_atm(1,i))
        write(sub_tmp,*) 'mfvirielc_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirielc_atm(2,i))
        write(sub_tmp,*) 'mfvirielc_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirielc_atm(3,i))

        write(sub_tmp,*) 'mfvirilj_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirilj_atm(1,i))
        write(sub_tmp,*) 'mfvirilj_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirilj_atm(2,i))
        write(sub_tmp,*) 'mfvirilj_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirilj_atm(3,i))

        write(sub_tmp,*) 'mfvirimor_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirimor_atm(1,i))
        write(sub_tmp,*) 'mfvirimor_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirimor_atm(2,i))
        write(sub_tmp,*) 'mfvirimor_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirimor_atm(3,i))

        write(sub_tmp,*) 'mfvirish_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirish_atm(1,i))
        write(sub_tmp,*) 'mfvirish_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirish_atm(2,i))
        write(sub_tmp,*) 'mfvirish_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirish_atm(3,i))

        write(sub_tmp,*) 'mfvirirfh_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirirfh_atm(1,i))
        write(sub_tmp,*) 'mfvirirfh_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirirfh_atm(2,i))
        write(sub_tmp,*) 'mfvirirfh_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfvirirfh_atm(3,i))

        write(sub_tmp,*) 'mfviridou_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviridou_atm(1,i))
        write(sub_tmp,*) 'mfviridou_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviridou_atm(2,i))
        write(sub_tmp,*) 'mfviridou_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviridou_atm(3,i))

        write(sub_tmp,*) 'mfviricstmnb_atm_xz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviricstmnb_atm(1,i))
        write(sub_tmp,*) 'mfviricstmnb_atm_yz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviricstmnb_atm(2,i))
        write(sub_tmp,*) 'mfviricstmnb_atm_zz_',i
        call excl_sp(sub_tmp,sub_len,sub_mfviricstmnb_atm(3,i))

        do nex = 1, ncstmnbex
           write(sub_tmp,*) 'mfviricsex_atm_',nex,'_xz_',i
           call excl_sp(sub_tmp,sub_len,sub_mfviricstmnbex_atm(1,i,nex))
           write(sub_tmp,*) 'mfviricsex_atm_',nex,'_yz_',i
           call excl_sp(sub_tmp,sub_len,sub_mfviricstmnbex_atm(2,i,nex))
           write(sub_tmp,*) 'mfviricsex_atm_',nex,'_zz_',i
           call excl_sp(sub_tmp,sub_len,sub_mfviricstmnbex_atm(3,i,nex))
        end do

     end do

#if defined(_HF_ALL_DIR)
#else
     !- creating output format
     write(fmt,'(''(A9,1X,''I8''(1X,A23))'')') nmffield
     ! 101  format(A9,1X,<n>(1X,A23))

     if (mtfoutdir == 'X') then
        write(oumtf,fmt) 'nstep', &
             &           (sub_mftot_atm(1,i),i=1,nhfregion), &
             &           (sub_mfkin_atm(1,i),i=1,nhfregion), &
             &           ((sub_mfkinex_atm(1,i,nex),i=1,nhfregion),nex=1,ncstmnbex), &
             &           (sub_mfviribo_atm(1,i),i=1,nhfregion), &
             &           (sub_mfvirian_atm(1,i),i=1,nhfregion), &
             &           (sub_mfvirito_atm(1,i),i=1,nhfregion), &
             &           (sub_mfviri14_atm(1,i),i=1,nhfregion), &
             &           (sub_mfvirielin_atm(1,i),i=1,nhfregion), &
             &           (sub_mfviriljin_atm(1,i),i=1,nhfregion), &
             &           (sub_mfvirielc_atm(1,i),i=1,nhfregion), &
             &           (sub_mfvirilj_atm(1,i),i=1,nhfregion), &
             &           (sub_mfvirimor_atm(1,i),i=1,nhfregion), &
             &           (sub_mfvirish_atm(1,i),i=1,nhfregion), &
             &           (sub_mfvirirfh_atm(1,i),i=1,nhfregion), &
             &           (sub_mfviridou_atm(1,i),i=1,nhfregion), &
             &           (sub_mfviricstmnb_atm(1,i),i=1,nhfregion), &
             &           ((sub_mfviricstmnbex_atm(1,i,nex),i=1,nhfregion),nex=1,ncstmnbex)
     else if (mtfoutdir == 'Y') then
        write(oumtf,fmt) 'nstep', &
             &           (sub_mftot_atm(2,i),i=1,nhfregion), &
             &           (sub_mfkin_atm(2,i),i=1,nhfregion), &
             &           ((sub_mfkinex_atm(2,i,nex),i=1,nhfregion),nex=1,ncstmnbex), &
             &           (sub_mfviribo_atm(2,i),i=1,nhfregion), &
             &           (sub_mfvirian_atm(2,i),i=1,nhfregion), &
             &           (sub_mfvirito_atm(2,i),i=1,nhfregion), &
             &           (sub_mfviri14_atm(2,i),i=1,nhfregion), &
             &           (sub_mfvirielin_atm(2,i),i=1,nhfregion), &
             &           (sub_mfviriljin_atm(2,i),i=1,nhfregion), &
             &           (sub_mfvirielc_atm(2,i),i=1,nhfregion), &
             &           (sub_mfvirilj_atm(2,i),i=1,nhfregion), &
             &           (sub_mfvirimor_atm(2,i),i=1,nhfregion), &
             &           (sub_mfvirish_atm(2,i),i=1,nhfregion), &
             &           (sub_mfvirirfh_atm(2,i),i=1,nhfregion), &
             &           (sub_mfviridou_atm(2,i),i=1,nhfregion), &
             &           (sub_mfviricstmnb_atm(2,i),i=1,nhfregion), &
             &           ((sub_mfviricstmnbex_atm(2,i,nex),i=1,nhfregion),nex=1,ncstmnbex)
     else if (mtfoutdir == 'Z') then
        write(oumtf,fmt) 'nstep', &
             &           (sub_mftot_atm(3,i),i=1,nhfregion), &
             &           (sub_mfkin_atm(3,i),i=1,nhfregion), &
             &           ((sub_mfkinex_atm(3,i,nex),i=1,nhfregion),nex=1,ncstmnbex), &
             &           (sub_mfviribo_atm(3,i),i=1,nhfregion), &
             &           (sub_mfvirian_atm(3,i),i=1,nhfregion), &
             &           (sub_mfvirito_atm(3,i),i=1,nhfregion), &
             &           (sub_mfviri14_atm(3,i),i=1,nhfregion), &
             &           (sub_mfvirielin_atm(3,i),i=1,nhfregion), &
             &           (sub_mfviriljin_atm(3,i),i=1,nhfregion), &
             &           (sub_mfvirielc_atm(3,i),i=1,nhfregion), &
             &           (sub_mfvirilj_atm(3,i),i=1,nhfregion), &
             &           (sub_mfvirimor_atm(3,i),i=1,nhfregion), &
             &           (sub_mfvirish_atm(3,i),i=1,nhfregion), &
             &           (sub_mfvirirfh_atm(3,i),i=1,nhfregion), &
             &           (sub_mfviridou_atm(3,i),i=1,nhfregion), &
             &           (sub_mfviricstmnb_atm(3,i),i=1,nhfregion), &
             &           ((sub_mfviricstmnbex_atm(3,i,nex),i=1,nhfregion),nex=1,ncstmnbex)
     else if (mtfoutdir == 'ALL') then
        !- creating output format
        write(fmt,'(''(A9,1X,''I8''(1X,A23))'')') nmffield*3
        ! 101  format(A9,1X,<n>(1X,A23))

        write(oumtf,fmt) 'nstep', &
             &           ((sub_mftot_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfkin_atm(n,i),n=1,3),i=1,nhfregion), &
             &           (((sub_mfkinex_atm(n,i,nex),n=1,3),i=1,nhfregion),nex=1,ncstmnbex), &
             &           ((sub_mfviribo_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfvirian_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfvirito_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfviri14_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfvirielin_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfviriljin_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfvirielc_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfvirilj_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfvirimor_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfvirish_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfvirirfh_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfviridou_atm(n,i),n=1,3),i=1,nhfregion), &
             &           ((sub_mfviricstmnb_atm(n,i),n=1,3),i=1,nhfregion), &
             &           (((sub_mfviricstmnbex_atm(n,i,nex),n=1,3),i=1,nhfregion),nex=1,ncstmnbex)
     end if
#endif

  end if

!     - if surface-based method and current step = 1, do not output data
  if ((.not. ifhfvol) .and. (current_step == 1)) return

#if defined(_HF_ALL_DIR)
#else
  !- creating output format
#if defined(_OUTHF_HIGH)
  write(fmt,'(''(I10,''I8''(1X,E23.16))'')') nmffield
  ! 102 format(I10,<n>(1X,E23.16)
#else
  write(fmt,'(''(I10,''I8''(1X,E23.8))'')') nmffield
  ! 102 format(I10,<n>(1X,E23.8))
#endif

  if (mtfoutdir == 'X') then
     write(oumtf,fmt) current_step, &
          &           (mftot_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfkin_atm(1,3,i)*pref,i=1,nhfregion), &
          &           ((mfkinex_atm(1,3,i,nex)*pref,i=1,nhfregion),nex=1,ncstmnbex), &
          &           (mfviribo_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfvirian_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfvirito_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfviri14_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfvirielin_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfviriljin_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfvirielc_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfvirilj_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfvirimor_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfvirish_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfvirirfh_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfviridou_atm(1,3,i)*pref,i=1,nhfregion), &
          &           (mfviricstmnb_atm(1,3,i)*pref,i=1,nhfregion), &
          &           ((mfviricstmnbex_atm(1,3,i,nex)*pref,i=1,nhfregion),nex=1,ncstmnbex)
  else if (mtfoutdir == 'Y') then
     write(oumtf,fmt) current_step, &
          &           (mftot_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfkin_atm(2,3,i)*pref,i=1,nhfregion), &
          &           ((mfkinex_atm(2,3,i,nex)*pref,i=1,nhfregion),nex=1,ncstmnbex), &
          &           (mfviribo_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfvirian_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfvirito_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfviri14_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfvirielin_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfviriljin_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfvirielc_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfvirilj_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfvirimor_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfvirish_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfvirirfh_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfviridou_atm(2,3,i)*pref,i=1,nhfregion), &
          &           (mfviricstmnb_atm(2,3,i)*pref,i=1,nhfregion), &
          &           ((mfviricstmnbex_atm(2,3,i,nex)*pref,i=1,nhfregion),nex=1,ncstmnbex)
  else if (mtfoutdir == 'Z') then
     write(oumtf,fmt) current_step, &
          &           (mftot_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfkin_atm(3,3,i)*pref,i=1,nhfregion), &
          &           ((mfkinex_atm(3,3,i,nex)*pref,i=1,nhfregion),nex=1,ncstmnbex), &
          &           (mfviribo_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfvirian_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfvirito_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfviri14_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfvirielin_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfviriljin_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfvirielc_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfvirilj_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfvirimor_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfvirish_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfvirirfh_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfviridou_atm(3,3,i)*pref,i=1,nhfregion), &
          &           (mfviricstmnb_atm(3,3,i)*pref,i=1,nhfregion), &
          &           ((mfviricstmnbex_atm(3,3,i,nex)*pref,i=1,nhfregion),nex=1,ncstmnbex)
  else if (mtfoutdir == 'ALL') then
     !- creating output format
#if defined(_OUTHF_HIGH)
     write(fmt,'(''(I10,''I8''(1X,E23.16))'')') nmffield*3
     ! 102 format(I10,<n>(1X,E23.16)
#else
     write(fmt,'(''(I10,''I8''(1X,E23.8))'')') nmffield*3
     ! 102 format(I10,<n>(1X,E23.8))
#endif

     write(oumtf,fmt) current_step, &
          &           ((mftot_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfkin_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           (((mfkinex_atm(n,3,i,nex)*pref,n=1,3),i=1,nhfregion),nex=1,ncstmnbex), &
          &           ((mfviribo_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfvirian_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfvirito_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfviri14_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfvirielin_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfviriljin_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfvirielc_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfvirilj_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfvirimor_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfvirish_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfvirirfh_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfviridou_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           ((mfviricstmnb_atm(n,3,i)*pref,n=1,3),i=1,nhfregion), &
          &           (((mfviricstmnbex_atm(n,3,i,nex)*pref,n=1,3),i=1,nhfregion),nex=1,ncstmnbex)
  end if
#endif

!     +     +     +     +     +     +     +     +

end subroutine outmtf
