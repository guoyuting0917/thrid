!*****************************
!*  output_thc.f Ver.1.1     *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-05-01 18:23:04 gota>

subroutine outthc(outhc,current_step, &
     &            xref,eref,timeref, &
     &            ntcregion, &
     &            det_ene_kin)

  use interface_tools

  use md_global

  implicit none

!     subroutine for outputting heat flux data

! ARGUMENTS:
!     INPUT
  integer,intent(in):: outhc           ! output unit for outthc thermal control data

  integer,intent(in):: current_step    ! current time step

  integer,intent(in):: ntcregion       ! number of region to control temp.

  real(8),intent(in):: det_ene_kin(:)  ! for outputting thermal control data

!---- base value for non-dimensionalize
  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: eref             ! energy base value [J]

  real(8),intent(in):: timeref          ! time base value [sec]

! LOCAL:
  character(80):: sub_tmp
  integer:: sub_len

  character(80):: sub_dene_atm(ntcregion)

  integer:: i

!     +     +     +     +     +     +     +     +

!     --- Output heat flux data ---

  if (current_step == 1) then

!    make subject
     sub_dene_atm(1:ntcregion) = ' '

     do i = 1, ntcregion
        write(sub_tmp,*) 'D_ene_',i
        call excl_sp(sub_tmp,sub_len,sub_dene_atm(i))
     end do

     write(outhc,101) 'nstep', &
          &           (sub_dene_atm(i),i=1,ntcregion)

101  format(A9,1X,2000(A15,1X))

  end if

  write(outhc,102) current_step, &
       &           (det_ene_kin(i)*eref,i=1,ntcregion)

102 format(I10,2000E16.8)

!     +     +     +     +     +     +     +     +

end subroutine outthc
