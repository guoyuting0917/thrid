!********************************************
!*  interface_timer.f90 Ver.1.0  '09.11.16  *
!*      for peachgk_md.f                    *
!*            by G.Kikugawa                 *
!********************************************

!***** This module is interface module for timer routines *****

module interface_timer

  interface

     subroutine get_timecount(sum_time,elps_time,mode_time)

       ! ARGUMENT:
       !     INPUT 
       character(5),intent(in):: mode_time ! mode for time count

       !     INPUT & OUTPUT
       real(8),intent(inout):: sum_time         ! cumulative time
       real(8),intent(inout):: elps_time        ! elapsed time

     end subroutine get_timecount

     subroutine get_timecount_barrier(sum_time,elps_time,mode_time)

       ! ARGUMENT:
       !     INPUT 
       character(5),intent(in):: mode_time ! mode for time count

       !     INPUT & OUTPUT
       real(8),intent(inout):: sum_time         ! cumulative time
       real(8),intent(inout):: elps_time        ! elapsed time

     end subroutine get_timecount_barrier

     subroutine init_timecount()

       ! ARGUMENT:
       !   NONE

     end subroutine init_timecount

     subroutine out_timecount()

       ! ARGUMENT:
       !   NONE

     end subroutine out_timecount
     
  end interface

end module interface_timer
