!**************************************
!*  time_method.f Ver.1.6 '13.12.20   *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine get_timecount(sum_time,elps_time,mode_time)

!     subroutine for time count (no barrier function of MPI)

  implicit none

#if defined(MPI)
  include 'mpif.h'          ! MPI interface declaration
#endif

! ARGUMENTS: 
!     INPUT 
  character(5),intent(in):: mode_time ! mode for time count

!     INPUT & OUTPUT
  real(8),intent(inout):: sum_time         ! cumulative time
  real(8),intent(inout):: elps_time        ! elapsed time

! LOCAL:
  real(4):: time(2)          ! argument for etime function
  real(4):: dummy_time       ! dummy

! FUNCTION:
  real(4):: etime            ! etime function

!     +     +     +     +     +     +     +

#if defined(MPI)
  if (mode_time == 'START') then ! Initialize of sum_time
     elps_time = MPI_WTIME()
  else
     sum_time = sum_time + MPI_WTIME() - elps_time
  end if
#else
  if (mode_time == 'START') then ! Initialize of sum_time
     dummy_time = etime(time)
     elps_time = DBLE(time(1))
  else
     dummy_time = etime(time)         
     sum_time = sum_time + DBLE(time(1)) - elps_time
  end if
#endif

!     +     +     +     +     +     +     +
    
  return
end subroutine get_timecount

!----------------------------------------------------------
subroutine get_timecount_barrier(sum_time,elps_time,mode_time)

!     subroutine for time count (with barrier function of MPI)

  implicit none

#if defined(MPI)
  include 'mpif.h'          ! MPI interface declaration
#endif

! ARGUMENTS: 
!     INPUT 
  character(5),intent(in):: mode_time ! mode for time count

!     INPUT & OUTPUT
  real(8),intent(inout):: sum_time         ! cumulative time
  real(8),intent(inout):: elps_time        ! elapsed time

! LOCAL:
  integer:: err             ! MPI error cord
  real(4):: time(2)          ! argument for etime function
  real(4):: dummy_time       ! dummy

! FUNCTION:
  real(4):: etime            ! etime function

!     +     +     +     +     +     +     +

#if defined(MPI)
  if (mode_time == 'START') then ! Initialize of sum_time
     call MPI_BARRIER(MPI_COMM_WORLD,err)
     elps_time = MPI_WTIME()
  else
     call MPI_BARRIER(MPI_COMM_WORLD,err)
     sum_time = sum_time + MPI_WTime() - elps_time
  end if
#else
  if (mode_time == 'START') then ! Initialize of sum_time
     dummy_time = etime(time)
     elps_time = DBLE(time(1))
  else
     dummy_time = etime(time)         
     sum_time = sum_time + DBLE(time(1)) - elps_time
  end if
#endif

!     +     +     +     +     +     +     +
    
  return
end subroutine get_timecount_barrier

!----------------------------------------------------------
subroutine init_timecount()

!     subroutine for initialize of timer variables

  use time_global

  implicit none

! ARGUMENTS:
!     INPUT

! LOCAL:

!     +     +     +     +     +     +     +
      
  total_time = 0.0d0

  dyn_time = 0.0d0
  int_time = 0.0d0

  comm_time = 0.0d0
      
  out_time = 0.0d0

  table_time = 0.0d0

  bond_time = 0.0d0
  drct_time = 0.0d0
  rcpr_time = 0.0d0
  mor_time = 0.0d0
  sh_time = 0.0d0
  rfh_time = 0.0d0
  dou_time = 0.0d0
  cnpvw_time = 0.0d0
  cstmnb_time = 0.0d0

  mdgcell_time = 0.0d0
  mdgboard_time = 0.0d0
  mdgexcs_time = 0.0d0
  mdgsum_time = 0.0d0

#if defined(HF)
  heatf_time = 0.0d0
#endif

!     +     +     +     +     +     +     +
  return
end subroutine init_timecount

!----------------------------------------------------------
subroutine out_timecount()

!     subroutine for output of timer result

#if defined(MPI)
  use mpi_global
#endif
  use time_global

  implicit none

! ARGUMENTS: 
!     INPUT

! LOCAL:
  integer:: unit            ! file unit to output time

!     +     +     +     +     +     +     +

#if defined(MPI)
  unit = 99 - irank
#else 
  unit = 99
#endif

!---- Calculate sum of certain component

  int_time = bond_time + drct_time + rcpr_time &
       &   + mor_time + sh_time + rfh_time + dou_time + cstmnb_time

!---- Calculate parcentage of cumulative time
  dyn_ratio = dyn_time / total_time
  int_ratio = int_time / total_time

  comm_ratio = comm_time / total_time
      
  out_ratio = out_time / total_time

  table_ratio = table_time / total_time

  bond_ratio = bond_time / total_time
  drct_ratio = drct_time / total_time
  rcpr_ratio = rcpr_time / total_time
  mor_ratio = mor_time / total_time
  sh_ratio = sh_time / total_time
  rfh_ratio = rfh_time / total_time
  dou_ratio = dou_time / total_time
  cnpvw_ratio = cnpvw_time / total_time
  cstmnb_ratio = cstmnb_time / total_time

#if defined(MDG3) || defined(MDG_m3)
  mdgcell_ratio = mdgcell_time / total_time
  mdgboard_ratio = mdgboard_time / total_time
  mdgexcs_ratio = mdgexcs_time / total_time
  mdgsum_ratio = mdgsum_time / total_time
#endif

#if defined(HF)
  heatf_ratio = heatf_time / total_time
#endif

!---- Output timer result
  write(unit,'(A12,E15.7,F10.4)') 'Total_Time',total_time,1.0
#if defined(TIME_MALL)
!     do not ouput other timer result
#else
  write(unit,'(A12,E15.7,F10.4)') 'Dyn_Time',dyn_time,dyn_ratio
  write(unit,'(A12,E15.7,F10.4)') 'Int_Time',int_time,int_ratio

  write(unit,'(A12,E15.7,F10.4)') 'Comm_Time',comm_time,comm_ratio

  write(unit,'(A12,E15.7,F10.4)') 'Out_Time',out_time,out_ratio

  write(unit,'(A12,E15.7,F10.4)') 'Table_Time',table_time, &
       &                           table_ratio

  write(unit,'(A12,E15.7,F10.4)') 'Bond_Time',bond_time,bond_ratio
  write(unit,'(A12,E15.7,F10.4)') 'Drct_Time',drct_time,drct_ratio
  write(unit,'(A12,E15.7,F10.4)') 'Rcpr_Time',rcpr_time,rcpr_ratio
  write(unit,'(A12,E15.7,F10.4)') 'Mor_Time',mor_time,mor_ratio
  write(unit,'(A12,E15.7,F10.4)') 'SH_Time',sh_time,sh_ratio
  write(unit,'(A12,E15.7,F10.4)') 'RFH_Time',rfh_time,rfh_ratio
  write(unit,'(A12,E15.7,F10.4)') 'DOU_Time',dou_time,dou_ratio
  write(unit,'(A12,E15.7,F10.4)') 'CNPVW_Time',cnpvw_time,cnpvw_ratio
  write(unit,'(A12,E15.7,F10.4)') 'CSTMNB_Time',cstmnb_time,cstmnb_ratio

#if defined(MDG3) || defined(MDG_m3)
  write(unit,'(A12,E15.7,F10.4)') 'MDG_Cell',mdgcell_time, &
       &                          mdgcell_ratio
  write(unit,'(A12,E15.7,F10.4)') 'MDG_Board',mdgboard_time, &
       &                          mdgboard_ratio
  write(unit,'(A12,E15.7,F10.4)') 'MDG_Excs',mdgexcs_time, &
       &                          mdgexcs_ratio
  write(unit,'(A12,E15.7,F10.4)') 'MDG_Sum',mdgsum_time, &
       &                          mdgsum_ratio
#endif

#if defined(HF)
  write(unit,'(A12,E15.7,F10.4)') 'Heatf_Sum',heatf_time, &
       &                          heatf_ratio
#endif

#endif
  write(unit,*)

!     +     +     +     +     +     +     +

  return
end subroutine out_timecount
