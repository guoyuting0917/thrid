!*****************************
!*  output_hf.f Ver.2.3      *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-08-28 13:21:18 gota>

subroutine outhtf(ouhtf,current_step, &
     &            xref,eref,timeref, &
     &            ifhfvol, &
     &            nhfregion,hfzpos1,hfzpos2, &
     &            hfkin_atm,hfkinex_atm, &
     &            hfpotbo_atm,hfpotan_atm, &
     &            hfpotto_atm,hfpot14_atm, &
     &            hfpotelin_atm,hfpotljin_atm, &
     &            hfpotelc_atm,hfpotlj_atm, &
     &            hfpotmor_atm,hfpotsh_atm, &
     &            hfpotrfh_atm,hfpotdou_atm, &
     &            hfpotcstmnb_atm, &
     &            ncstmnbex,hfpotcstmnbex_atm, &
     &            hfviribo_atm,hfvirian_atm, &
     &            hfvirito_atm,hfviri14_atm, &
     &            hfvirielin_atm,hfviriljin_atm, &
     &            hfvirielc_atm,hfvirilj_atm, &
     &            hfvirimor_atm,hfvirish_atm, &
     &            hfvirirfh_atm,hfviridou_atm, &
     &            hfviricstmnb_atm, &
     &            hfviricstmnbex_atm, &
     &            hftot_atm)

  use interface_tools

  use md_global

  implicit none

!     subroutine for outputting heat flux data

! ARGUMENTS:
!     INPUT
  integer,intent(in):: ouhtf           ! output unit for outhtf heat flux data

  integer,intent(in):: current_step    ! current time step

!---- base value for non-dimensionalize
  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: eref             ! energy base value [J]

  real(8),intent(in):: timeref          ! time base value [sec]

!---- variables for calculation of heat flux
  logical,intent(in):: ifhfvol        ! local volume-based or local surface-based

  integer,intent(in):: nhfregion       ! number of region to calculate heat flux
  real(8),intent(in):: hfzpos1(:),hfzpos2(:)
                                ! z-position of region for heat flux

!---- valiables for heat flux
  real(8),intent(in):: hfkin_atm(:,:)     ! heat flux of kinetic part
  real(8),intent(in):: hfkinex_atm(:,:,:)
                                 ! heat flux of kinetic part (extra custom NB)

  real(8),intent(in):: hfpotbo_atm(:,:)   ! heat flux of potential part (bond)
  real(8),intent(in):: hfpotan_atm(:,:)   ! heat flux of potential part (angle)
  real(8),intent(in):: hfpotto_atm(:,:)   ! heat flux of potential part (tors)
  real(8),intent(in):: hfpot14_atm(:,:)   ! heat flux of potential part (1-4)
  real(8),intent(in):: hfpotelin_atm(:,:) 
                    ! heat flux of potential part (elc intra)
  real(8),intent(in):: hfpotljin_atm(:,:)
                    ! heat flux of potential part (L-J intra)

  real(8),intent(in):: hfpotelc_atm(:,:) ! heat flux of potential part (Coulomb)
  real(8),intent(in):: hfpotlj_atm(:,:) ! heat flux of potential part (L-J)
  real(8),intent(in):: hfpotmor_atm(:,:) ! heat flux of potential part (Morse)
  real(8),intent(in):: hfpotsh_atm(:,:) ! heat flux of potential part (SH)
  real(8),intent(in):: hfpotrfh_atm(:,:) ! heat flux of potential part (RFH)
  real(8),intent(in):: hfpotdou_atm(:,:) ! heat flux of potential part (DOU)
  real(8),intent(in):: hfpotcstmnb_atm(:,:)
                                     ! heat flux of potential part (custom NB)
  real(8),intent(in):: hfpotcstmnbex_atm(:,:,:)
                                ! heat flux of potential part (extra custom NB)

  real(8),intent(in):: hfviribo_atm(:,:)   ! heat flux of virial part (bond)
  real(8),intent(in):: hfvirian_atm(:,:)   ! heat flux of virial part (angle)
  real(8),intent(in):: hfvirito_atm(:,:)   ! heat flux of virial part (tors)
  real(8),intent(in):: hfviri14_atm(:,:)   ! heat flux of virial part (1-4)
  real(8),intent(in):: hfvirielin_atm(:,:) ! heat flux of virial part (elc intra)
  real(8),intent(in):: hfviriljin_atm(:,:) ! heat flux of virial part (L-J intra)

  real(8),intent(in):: hfvirielc_atm(:,:) ! heat flux of virial part (Coulomb)
  real(8),intent(in):: hfvirilj_atm(:,:) ! heat flux of virial part (L-J)
  real(8),intent(in):: hfvirimor_atm(:,:) ! heat flux of virial part (Morse)
  real(8),intent(in):: hfvirish_atm(:,:) ! heat flux of virial part (SH)
  real(8),intent(in):: hfvirirfh_atm(:,:) ! heat flux of virial part (RFH)
  real(8),intent(in):: hfviridou_atm(:,:) ! heat flux of virial part (DOU)
  real(8),intent(in):: hfviricstmnb_atm(:,:) 
                                        ! heat flux of virial part (custom NB)
  real(8),intent(in):: hfviricstmnbex_atm(:,:,:)
                                  ! heat flux of virial part (extra custom NB)

  real(8),intent(in):: hftot_atm(:,:) ! total heat flux

  integer,intent(in):: ncstmnbex      ! number of extra custom NB output

! LOCAL:
  character(80):: sub_tmp
  integer:: sub_len

  character(80):: sub_hftot_atm(3,nhfregion)
  character(80):: sub_hfkin_atm(3,nhfregion)
  character(80):: sub_hfkinex_atm(3,nhfregion,1:ncstmnbex)
  character(80):: sub_hfpotbo_atm(3,nhfregion)
  character(80):: sub_hfpotan_atm(3,nhfregion)
  character(80):: sub_hfpotto_atm(3,nhfregion)
  character(80):: sub_hfpot14_atm(3,nhfregion)
  character(80):: sub_hfpotelin_atm(3,nhfregion)
  character(80):: sub_hfpotljin_atm(3,nhfregion)
  character(80):: sub_hfpotelc_atm(3,nhfregion)
  character(80):: sub_hfpotlj_atm(3,nhfregion)
  character(80):: sub_hfpotmor_atm(3,nhfregion)
  character(80):: sub_hfpotsh_atm(3,nhfregion)
  character(80):: sub_hfpotrfh_atm(3,nhfregion)
  character(80):: sub_hfpotdou_atm(3,nhfregion)
  character(80):: sub_hfpotcstmnb_atm(3,nhfregion)
  character(80):: sub_hfpotcstmnbex_atm(3,nhfregion,1:ncstmnbex)
  character(80):: sub_hfviribo_atm(3,nhfregion)
  character(80):: sub_hfvirian_atm(3,nhfregion)
  character(80):: sub_hfvirito_atm(3,nhfregion)
  character(80):: sub_hfviri14_atm(3,nhfregion)
  character(80):: sub_hfvirielin_atm(3,nhfregion)
  character(80):: sub_hfviriljin_atm(3,nhfregion)
  character(80):: sub_hfvirielc_atm(3,nhfregion)
  character(80):: sub_hfvirilj_atm(3,nhfregion)
  character(80):: sub_hfvirimor_atm(3,nhfregion)
  character(80):: sub_hfvirish_atm(3,nhfregion)
  character(80):: sub_hfvirirfh_atm(3,nhfregion)
  character(80):: sub_hfviridou_atm(3,nhfregion)
  character(80):: sub_hfviricstmnb_atm(3,nhfregion)
  character(80):: sub_hfviricstmnbex_atm(3,nhfregion,1:ncstmnbex)


                                                                       !!!!!!!!!!!!!!!!!!!!!!define ar1-ar2 hfsub
  character(80):: sub_hfpotlj_atm_ar11(3,nhfregion)                                         
  character(80):: sub_hfpotlj_atm_ar12(3,nhfregion)
  character(80):: sub_hfpotlj_atm_ar22(3,nhfregion)

  character(80):: sub_hfvirilj_atm_ar11(3,nhfregion)
  character(80):: sub_hfvirilj_atm_ar12(3,nhfregion)
  character(80):: sub_hfvirilj_atm_ar22(3,nhfregion)
  character(80):: sub_hfvirilj_atm_ar1pt(3,nhfregion)
  character(80):: sub_hfvirilj_atm_ar2pt(3,nhfregion)

  real(8):: jref             ! heat flux base value [W/m2]

  integer:: i
  integer:: n
  integer:: nex

  character(25):: fmt
  integer:: nhffield

!     +     +     +     +     +     +     +     +

!---- some preparation

  !- count number of output fields used below
  nhffield = 28 +6 +2               ! fixed number of fields (in one direction)
  nhffield = (nhffield + ncstmnbex*1) * nhfregion

!     --- calculate base value of heat flux ---

  jref = eref / (xref*xref * timeref)

!     --- Output heat flux data ---

  if (current_step == 1) then

!        make subject
     sub_hftot_atm(1:3,1:nhfregion) = ' '
     sub_hfkin_atm(1:3,1:nhfregion) = ' '
     sub_hfkinex_atm(1:3,1:nhfregion,1:ncstmnbex) = ' '
     sub_hfpotbo_atm(1:3,1:nhfregion) = ' '
     sub_hfpotan_atm(1:3,1:nhfregion) = ' '
     sub_hfpotto_atm(1:3,1:nhfregion) = ' '
     sub_hfpot14_atm(1:3,1:nhfregion) = ' '
     sub_hfpotelin_atm(1:3,1:nhfregion) = ' '
     sub_hfpotljin_atm(1:3,1:nhfregion) = ' '
     sub_hfpotelc_atm(1:3,1:nhfregion) = ' '
     sub_hfpotlj_atm(1:3,1:nhfregion) = ' '
     sub_hfpotmor_atm(1:3,1:nhfregion) = ' '
     sub_hfpotsh_atm(1:3,1:nhfregion) = ' '
     sub_hfpotrfh_atm(1:3,1:nhfregion) = ' '
     sub_hfpotdou_atm(1:3,1:nhfregion) = ' '
     sub_hfpotcstmnb_atm(1:3,1:nhfregion) = ' '
     sub_hfpotcstmnbex_atm(1:3,1:nhfregion,1:ncstmnbex) = ' '
     sub_hfviribo_atm(1:3,1:nhfregion) = ' '
     sub_hfvirian_atm(1:3,1:nhfregion) = ' '
     sub_hfvirito_atm(1:3,1:nhfregion) = ' '
     sub_hfviri14_atm(1:3,1:nhfregion) = ' '
     sub_hfvirielin_atm(1:3,1:nhfregion) = ' '
     sub_hfviriljin_atm(1:3,1:nhfregion) = ' '
     sub_hfvirielc_atm(1:3,1:nhfregion) = ' '
     sub_hfvirilj_atm(1:3,1:nhfregion) = ' '
     sub_hfvirimor_atm(1:3,1:nhfregion) = ' '
     sub_hfvirish_atm(1:3,1:nhfregion) = ' '
     sub_hfvirirfh_atm(1:3,1:nhfregion) = ' '
     sub_hfviridou_atm(1:3,1:nhfregion) = ' '
     sub_hfviricstmnb_atm(1:3,1:nhfregion) = ' '
     sub_hfviricstmnbex_atm(1:3,1:nhfregion,1:ncstmnbex) = ' '

                                                                    !!!!!!!!!!!!!!!!!!!!!!value ar1-ar2 sub_hf(pot + viri)          
     sub_hfpotlj_atm_ar11(1:3,1:nhfregion) = ' '
     sub_hfpotlj_atm_ar12(1:3,1:nhfregion) = ' '
     sub_hfpotlj_atm_ar22(1:3,1:nhfregion) = ' '
    
     sub_hfvirilj_atm_ar11(1:3,1:nhfregion) = ' '
     sub_hfvirilj_atm_ar12(1:3,1:nhfregion) = ' '
     sub_hfvirilj_atm_ar22(1:3,1:nhfregion) = ' '
     sub_hfvirilj_atm_ar1pt(1:3,1:nhfregion) = ' '
     sub_hfvirilj_atm_ar2pt(1:3,1:nhfregion) = ' '


     do i = 1, nhfregion
        write(sub_tmp,*) 'hftot_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hftot_atm(1,i))
        write(sub_tmp,*) 'hftot_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hftot_atm(2,i))
        write(sub_tmp,*) 'hftot_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hftot_atm(3,i))
        write(sub_tmp,*) 'hfkin_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfkin_atm(1,i))
        write(sub_tmp,*) 'hfkin_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfkin_atm(2,i))
        write(sub_tmp,*) 'hfkin_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfkin_atm(3,i))
        do nex = 1, ncstmnbex
           write(sub_tmp,*) 'hfkinex_atm_',nex,'_x_',i
           call excl_sp(sub_tmp,sub_len,sub_hfkinex_atm(1,i,nex))
           write(sub_tmp,*) 'hfkinex_atm_',nex,'_y_',i
           call excl_sp(sub_tmp,sub_len,sub_hfkinex_atm(2,i,nex))
           write(sub_tmp,*) 'hfkinex_atm_',nex,'_z_',i
           call excl_sp(sub_tmp,sub_len,sub_hfkinex_atm(3,i,nex))
        end do

        write(sub_tmp,*) 'hfpotbo_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotbo_atm(1,i))
        write(sub_tmp,*) 'hfpotbo_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotbo_atm(2,i))
        write(sub_tmp,*) 'hfpotbo_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotbo_atm(3,i))
        write(sub_tmp,*) 'hfpotan_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotan_atm(1,i))
        write(sub_tmp,*) 'hfpotan_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotan_atm(2,i))
        write(sub_tmp,*) 'hfpotan_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotan_atm(3,i))
        write(sub_tmp,*) 'hfpotto_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotto_atm(1,i))
        write(sub_tmp,*) 'hfpotto_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotto_atm(2,i))
        write(sub_tmp,*) 'hfpotto_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotto_atm(3,i))
        write(sub_tmp,*) 'hfpot14_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpot14_atm(1,i))
        write(sub_tmp,*) 'hfpot14_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpot14_atm(2,i))
        write(sub_tmp,*) 'hfpot14_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpot14_atm(3,i))
        write(sub_tmp,*) 'hfpotelin_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotelin_atm(1,i))
        write(sub_tmp,*) 'hfpotelin_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotelin_atm(2,i))
        write(sub_tmp,*) 'hfpotelin_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotelin_atm(3,i))
        write(sub_tmp,*) 'hfpotljin_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotljin_atm(1,i))
        write(sub_tmp,*) 'hfpotljin_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotljin_atm(2,i))
        write(sub_tmp,*) 'hfpotljin_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotljin_atm(3,i))
        write(sub_tmp,*) 'hfpotelc_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotelc_atm(1,i))
        write(sub_tmp,*) 'hfpotelc_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotelc_atm(2,i))
        write(sub_tmp,*) 'hfpotelc_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotelc_atm(3,i))
        write(sub_tmp,*) 'hfpotlj_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm(1,i))
        write(sub_tmp,*) 'hfpotlj_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm(2,i))
        write(sub_tmp,*) 'hfpotlj_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm(3,i))
        write(sub_tmp,*) 'hfpotmor_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotmor_atm(1,i))
        write(sub_tmp,*) 'hfpotmor_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotmor_atm(2,i))
        write(sub_tmp,*) 'hfpotmor_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotmor_atm(3,i))
        write(sub_tmp,*) 'hfpotsh_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotsh_atm(1,i))
        write(sub_tmp,*) 'hfpotsh_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotsh_atm(2,i))
        write(sub_tmp,*) 'hfpotsh_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotsh_atm(3,i))
        write(sub_tmp,*) 'hfpotrfh_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotrfh_atm(1,i))
        write(sub_tmp,*) 'hfpotrfh_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotrfh_atm(2,i))
        write(sub_tmp,*) 'hfpotrfh_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotrfh_atm(3,i))
        write(sub_tmp,*) 'hfpotdou_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotdou_atm(1,i))
        write(sub_tmp,*) 'hfpotdou_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotdou_atm(2,i))
        write(sub_tmp,*) 'hfpotdou_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotdou_atm(3,i))
        write(sub_tmp,*) 'hfpotcstmnb_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotcstmnb_atm(1,i))
        write(sub_tmp,*) 'hfpotcstmnb_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotcstmnb_atm(2,i))
        write(sub_tmp,*) 'hfpotcstmnb_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotcstmnb_atm(3,i))

                                                                 !!!!!!!!!!!!!!!!!!!!!!!!write ar1-ar2 hfsub  
        write(sub_tmp,*) 'hfpotlj_ar11_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm_ar11(1,i))
        write(sub_tmp,*) 'hfpotlj_ar11_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm_ar11(2,i))
        write(sub_tmp,*) 'hfpotlj_ar11_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm_ar11(3,i))

        write(sub_tmp,*) 'hfpotlj_ar12_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm_ar12(1,i))
        write(sub_tmp,*) 'hfpotlj_ar12_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm_ar12(2,i))
        write(sub_tmp,*) 'hfpotlj_ar12_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm_ar12(3,i))

        write(sub_tmp,*) 'hfpotlj_ar22_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm_ar22(1,i))
        write(sub_tmp,*) 'hfpotlj_ar22_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm_ar22(2,i))
        write(sub_tmp,*) 'hfpotlj_ar22_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfpotlj_atm_ar22(3,i))    




        write(sub_tmp,*) 'hfvirilj_ar11_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar11(1,i))
        write(sub_tmp,*) 'hfvirilj_ar11_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar11(2,i))
        write(sub_tmp,*) 'hfvirilj_ar11_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar11(3,i))

        write(sub_tmp,*) 'hfvirilj_ar12_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar12(1,i))
        write(sub_tmp,*) 'hfvirilj_ar12_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar12(2,i))
        write(sub_tmp,*) 'hfvirilj_ar12_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar12(3,i))

        write(sub_tmp,*) 'hfvirilj_ar22_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar22(1,i))
        write(sub_tmp,*) 'hfvirilj_ar22_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar22(2,i))
        write(sub_tmp,*) 'hfvirilj_ar22_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar22(3,i))

        write(sub_tmp,*) 'hfvirilj_ar1pt_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar1pt(1,i))
        write(sub_tmp,*) 'hfvirilj_ar1pt_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar1pt(2,i))
        write(sub_tmp,*) 'hfvirilj_ar1pt_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar1pt(3,i))

        write(sub_tmp,*) 'hfvirilj_ar2pt_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar2pt(1,i))
        write(sub_tmp,*) 'hfvirilj_ar2pt_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar2pt(2,i))
        write(sub_tmp,*) 'hfvirilj_ar2pt_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm_ar2pt(3,i))

                                                                    !!!!!!!!!!!!!!!!!!!!!!!!write ar1-ar2 hfsub 
        do nex = 1, ncstmnbex
           write(sub_tmp,*) 'hfpotcsex_atm_',nex,'_x_',i
           call excl_sp(sub_tmp,sub_len,sub_hfpotcstmnbex_atm(1,i,nex))
           write(sub_tmp,*) 'hfpotcsex_atm_',nex,'_y_',i
           call excl_sp(sub_tmp,sub_len,sub_hfpotcstmnbex_atm(2,i,nex))
           write(sub_tmp,*) 'hfpotcsex_atm_',nex,'_z_',i
           call excl_sp(sub_tmp,sub_len,sub_hfpotcstmnbex_atm(3,i,nex))
        end do

        write(sub_tmp,*) 'hfviribo_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviribo_atm(1,i))
        write(sub_tmp,*) 'hfviribo_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviribo_atm(2,i))
        write(sub_tmp,*) 'hfviribo_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviribo_atm(3,i))
        write(sub_tmp,*) 'hfvirian_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirian_atm(1,i))
        write(sub_tmp,*) 'hfvirian_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirian_atm(2,i))
        write(sub_tmp,*) 'hfvirian_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirian_atm(3,i))
        write(sub_tmp,*) 'hfvirito_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirito_atm(1,i))
        write(sub_tmp,*) 'hfvirito_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirito_atm(2,i))
        write(sub_tmp,*) 'hfvirito_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirito_atm(3,i))
        write(sub_tmp,*) 'hfviri14_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviri14_atm(1,i))
        write(sub_tmp,*) 'hfviri14_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviri14_atm(2,i))
        write(sub_tmp,*) 'hfviri14_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviri14_atm(3,i))
        write(sub_tmp,*) 'hfvirielin_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirielin_atm(1,i))
        write(sub_tmp,*) 'hfvirielin_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirielin_atm(2,i))
        write(sub_tmp,*) 'hfvirielin_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirielin_atm(3,i))
        write(sub_tmp,*) 'hfviriljin_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviriljin_atm(1,i))
        write(sub_tmp,*) 'hfviriljin_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviriljin_atm(2,i))
        write(sub_tmp,*) 'hfviriljin_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviriljin_atm(3,i))
        write(sub_tmp,*) 'hfvirielc_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirielc_atm(1,i))
        write(sub_tmp,*) 'hfvirielc_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirielc_atm(2,i))
        write(sub_tmp,*) 'hfvirielc_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirielc_atm(3,i))
        write(sub_tmp,*) 'hfvirilj_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm(1,i))
        write(sub_tmp,*) 'hfvirilj_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm(2,i))
        write(sub_tmp,*) 'hfvirilj_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirilj_atm(3,i))
        write(sub_tmp,*) 'hfvirimor_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirimor_atm(1,i))
        write(sub_tmp,*) 'hfvirimor_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirimor_atm(2,i))
        write(sub_tmp,*) 'hfvirimor_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirimor_atm(3,i))
        write(sub_tmp,*) 'hfvirish_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirish_atm(1,i))
        write(sub_tmp,*) 'hfvirish_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirish_atm(2,i))
        write(sub_tmp,*) 'hfvirish_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirish_atm(3,i))
        write(sub_tmp,*) 'hfvirirfh_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirirfh_atm(1,i))
        write(sub_tmp,*) 'hfvirirfh_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirirfh_atm(2,i))
        write(sub_tmp,*) 'hfvirirfh_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfvirirfh_atm(3,i))
        write(sub_tmp,*) 'hfviridou_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviridou_atm(1,i))
        write(sub_tmp,*) 'hfviridou_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviridou_atm(2,i))
        write(sub_tmp,*) 'hfviridou_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviridou_atm(3,i))
        write(sub_tmp,*) 'hfviricstmnb_atm_x_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviricstmnb_atm(1,i))
        write(sub_tmp,*) 'hfviricstmnb_atm_y_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviricstmnb_atm(2,i))
        write(sub_tmp,*) 'hfviricstmnb_atm_z_',i
        call excl_sp(sub_tmp,sub_len,sub_hfviricstmnb_atm(3,i))
        do nex = 1, ncstmnbex
           write(sub_tmp,*) 'hfviricsex_atm_',nex,'_x_',i
           call excl_sp(sub_tmp,sub_len,sub_hfviricstmnbex_atm(1,i,nex))
           write(sub_tmp,*) 'hfviricsex_atm_',nex,'_y_',i
           call excl_sp(sub_tmp,sub_len,sub_hfviricstmnbex_atm(2,i,nex))
           write(sub_tmp,*) 'hfviricsex_atm_',nex,'_z_',i
           call excl_sp(sub_tmp,sub_len,sub_hfviricstmnbex_atm(3,i,nex))
        end do

     end do

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
     !- creating output format
     write(fmt,'(''(A9,1X,''I8''(1X,A23))'')') nhffield*3
     ! 101  format(A9,1X,<n>(1X,A23))

     write(ouhtf,fmt) 'nstep', &
          &           ((sub_hftot_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfkin_atm(n,i),n=1,3),i=1,nhfregion), &
          &           (((sub_hfkinex_atm(n,i,nex),n=1,3),i=1,nhfregion),nex=1,ncstmnbex), &
          &           ((sub_hfpotbo_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotan_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotto_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpot14_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotelin_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotljin_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotelc_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotlj_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotmor_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotsh_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotrfh_atm(n,i),n=1,3),i=1,nhfregion), & 
          &           ((sub_hfpotdou_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotcstmnb_atm(n,i),n=1,3),i=1,nhfregion), &
  !        &           (((sub_hfpotcstmnbex_atm(n,i,nex),n=1,3),i=1,nhfregion),nex=1,ncstmnbex), &
          &           ((sub_hfviribo_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirian_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirito_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfviri14_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirielin_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfviriljin_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirielc_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirilj_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirimor_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirish_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirirfh_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfviridou_atm(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfviricstmnb_atm(n,i),n=1,3),i=1,nhfregion), &
    !      &           (((sub_hfviricstmnbex_atm(n,i,nex),n=1,3),i=1,nhfregion),nex=1,ncstmnbex), &
          &           ((sub_hfpotlj_atm_ar11(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotlj_atm_ar22(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfpotlj_atm_ar12(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirilj_atm_ar11(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirilj_atm_ar22(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirilj_atm_ar12(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirilj_atm_ar1pt(n,i),n=1,3),i=1,nhfregion), &
          &           ((sub_hfvirilj_atm_ar2pt(n,i),n=1,3),i=1,nhfregion)

#else
     !- creating output format
     write(fmt,'(''(A9,1X,''I8''(1X,A23))'')') nhffield
     ! 101  format(A9,1X,<n>(1X,A23))

     write(ouhtf,fmt) 'nstep', &
          &           (sub_hftot_atm(3,i),i=1,nhfregion), &
          &           (sub_hfkin_atm(3,i),i=1,nhfregion), &
          &           ((sub_hfkinex_atm(3,i,nex),i=1,nhfregion),nex=1,ncstmnbex), &
          &           (sub_hfpotbo_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpotan_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpotto_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpot14_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpotelin_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpotljin_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpotelc_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpotlj_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpotmor_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpotsh_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpotrfh_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpotdou_atm(3,i),i=1,nhfregion), &
          &           (sub_hfpotcstmnb_atm(3,i),i=1,nhfregion), &
  !        &           ((sub_hfpotcstmnbex_atm(3,i,nex),i=1,nhfregion),nex=1,ncstmnbex), &
          &           (sub_hfviribo_atm(3,i),i=1,nhfregion), &
          &           (sub_hfvirian_atm(3,i),i=1,nhfregion), &
          &           (sub_hfvirito_atm(3,i),i=1,nhfregion), &
          &           (sub_hfviri14_atm(3,i),i=1,nhfregion), &
          &           (sub_hfvirielin_atm(3,i),i=1,nhfregion), &
          &           (sub_hfviriljin_atm(3,i),i=1,nhfregion), &
          &           (sub_hfvirielc_atm(3,i),i=1,nhfregion), &
          &           (sub_hfvirilj_atm(3,i),i=1,nhfregion), &
          &           (sub_hfvirimor_atm(3,i),i=1,nhfregion), &
          &           (sub_hfvirish_atm(3,i),i=1,nhfregion), &
          &           (sub_hfvirirfh_atm(3,i),i=1,nhfregion), &
          &           (sub_hfviridou_atm(3,i),i=1,nhfregion), &
          &           (sub_hfviricstmnb_atm(3,i),i=1,nhfregion), &
   !       &           ((sub_hfviricstmnbex_atm(3,i,nex),i=1,nhfregion),nex=1,ncstmnbex), &
          &           (sub_hfpotlj_atm_ar11(3,i),i=1,nhfregion), &
          &           (sub_hfpotlj_atm_ar22(3,i),i=1,nhfregion), &
          &           (sub_hfpotlj_atm_ar12(3,i),i=1,nhfregion), &
          &           (sub_hfvirilj_atm_ar11(3,i),i=1,nhfregion), &
          &           (sub_hfvirilj_atm_ar22(3,i),i=1,nhfregion), &
          &           (sub_hfvirilj_atm_ar12(3,i),i=1,nhfregion), &
          &           (sub_hfvirilj_atm_ar1pt(3,i),i=1,nhfregion), &
          &           (sub_hfvirilj_atm_ar2pt(3,i),i=1,nhfregion)

#endif

  end if

!     - if surface-based method and current step = 1, do not output data
  if ((.not. ifhfvol) .and. (current_step == 1)) return

#if defined(_HF_ALL_DIR) || defined(_HF_BULK)
  !- creating output format
#if defined(_OUTHF_HIGH)
  write(fmt,'(''(I10,''I8''(1X,E23.16))'')') nhffield*3
  ! 102 format(I10,<n>(1X,E23.16))
#else
  write(fmt,'(''(I10,''I8''(1X,E23.8))'')') nhffield*3
  ! 102 format(I10,<n>(1X,E23.8))
#endif

  write(ouhtf,fmt) current_step, &
       &           ((hftot_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfkin_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           (((hfkinex_atm(n,i,nex)*jref,n=1,3),i=1,nhfregion),nex=1,ncstmnbex), &
       &           ((hfpotbo_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotan_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotto_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpot14_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotelin_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotljin_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotelc_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotlj_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotmor_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotsh_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotrfh_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotdou_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotcstmnb_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
!       &           (((hfpotcstmnbex_atm(n,i,nex)*jref,n=1,3),i=1,nhfregion),nex=1,ncstmnbex), &
       &           ((hfviribo_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirian_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirito_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfviri14_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirielin_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfviriljin_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirielc_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirilj_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirimor_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirish_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirirfh_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfviridou_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfviricstmnb_atm(n,i)*jref,n=1,3),i=1,nhfregion), &
   !    &           (((hfviricstmnbex_atm(n,i,nex)*jref,n=1,3),i=1,nhfregion),nex=1,ncstmnbex), &
       &           ((hfpotlj_atm_ar11(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotlj_atm_ar22(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfpotlj_atm_ar12(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirilj_atm_ar11(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirilj_atm_ar22(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirilj_atm_ar12(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirilj_atm_ar1pt(n,i)*jref,n=1,3),i=1,nhfregion), &
       &           ((hfvirilj_atm_ar2pt(n,i)*jref,n=1,3),i=1,nhfregion)
#else
  !- creating output format
#if defined(_OUTHF_HIGH)
  write(fmt,'(''(I10,''I8''(1X,E23.16))'')') nhffield
  ! 102 format(I10,<n>(1X,E23.16))
#else
  write(fmt,'(''(I10,''I8''(1X,E23.8))'')') nhffield
  ! 102 format(I10,<n>(1X,E23.8))
#endif

  write(ouhtf,fmt) current_step, &
       &           (hftot_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfkin_atm(3,i)*jref,i=1,nhfregion), &
       &           ((hfkinex_atm(3,i,nex)*jref,i=1,nhfregion),nex=1,ncstmnbex), &
       &           (hfpotbo_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpotan_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpotto_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpot14_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpotelin_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpotljin_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpotelc_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpotlj_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpotmor_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpotsh_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpotrfh_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpotdou_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfpotcstmnb_atm(3,i)*jref,i=1,nhfregion), &
 !      &           ((hfpotcstmnbex_atm(3,i,nex)*jref,i=1,nhfregion),nex=1,ncstmnbex), &
       &           (hfviribo_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfvirian_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfvirito_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfviri14_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfvirielin_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfviriljin_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfvirielc_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfvirilj_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfvirimor_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfvirish_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfvirirfh_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfviridou_atm(3,i)*jref,i=1,nhfregion), &
       &           (hfviricstmnb_atm(3,i)*jref,i=1,nhfregion), &
!       &           ((hfviricstmnbex_atm(3,i,nex)*jref,i=1,nhfregion),nex=1,ncstmnbex), &
       &           (hfpotlj_atm_ar11(3,i)*jref,i=1,nhfregion), &
       &           (hfpotlj_atm_ar22(3,i)*jref,i=1,nhfregion), &
       &           (hfpotlj_atm_ar12(3,i)*jref,i=1,nhfregion), &
       &           (hfvirilj_atm_ar11(3,i)*jref,i=1,nhfregion), &
       &           (hfvirilj_atm_ar22(3,i)*jref,i=1,nhfregion), &
       &           (hfvirilj_atm_ar12(3,i)*jref,i=1,nhfregion), &
       &           (hfvirilj_atm_ar1pt(3,i)*jref,i=1,nhfregion), &
       &           (hfvirilj_atm_ar2pt(3,i)*jref,i=1,nhfregion)
#endif

!     +     +     +     +     +     +     +     +

end subroutine outhtf
