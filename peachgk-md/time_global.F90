!*******************************
!*  time_global.f90  13.12.20  *
!*      for peachgk_md.f       *
!*            by G.Kikugawa    *
!*******************************
module time_global

  implicit none

  !---- cumulative time
  real(8),save:: total_time        ! total time

  real(8),save:: dyn_time          ! dynamics time
  real(8),save:: int_time          ! interaction time

  real(8),save:: comm_time         ! communication time

  real(8),save:: out_time          ! output time

  real(8),save:: table_time        ! table time

  real(8),save:: bond_time         ! bond time
  real(8),save:: drct_time         ! direct sum time
  real(8),save:: rcpr_time         ! reciprocal time
  real(8),save:: mor_time          ! morse time
  real(8),save:: sh_time           ! SH time
  real(8),save:: rfh_time          ! RFH time
  real(8),save:: dou_time          ! DOU time
  real(8),save:: cnpvw_time        ! CNP_VW time
  real(8),save:: cstmnb_time       ! CSTMNB time

  real(8),save:: mdgcell_time      ! MDGRAPE cell-index time
  real(8),save:: mdgboard_time     ! MDGRAPE board time
  real(8),save:: mdgexcs_time      ! MDGRAPE excs time
  real(8),save:: mdgsum_time       ! MDGRAPE sum time

  real(8),save:: heatf_time        ! heatflux sum time

  !---- elapsed time
  real(8),save:: total_time_elps   ! total elapsed time

  real(8),save:: dyn_time_elps     ! dynamics elapsed time
  real(8),save:: int_time_elps     ! interaction elapsed time

  real(8),save:: comm_time_elps    ! communication elapsed time

  real(8),save:: out_time_elps     ! output elapsed time

  real(8),save:: table_time_elps   ! table elapsed time

  real(8),save:: bond_time_elps    ! bond elapsed time
  real(8),save:: drct_time_elps    ! direct sum elapsed time
  real(8),save:: rcpr_time_elps    ! reciprocal elapsed time
  real(8),save:: mor_time_elps     ! morse elapsed time
  real(8),save:: sh_time_elps      ! SH elapsed time
  real(8),save:: rfh_time_elps     ! RFH elapsed time
  real(8),save:: dou_time_elps     ! DOU elapsed time
  real(8),save:: cnpvw_time_elps   ! CNP_VW elapsed time
  real(8),save:: cstmnb_time_elps  ! CSTMNB elapsed time

  real(8),save:: mdgcell_time_elps  ! MDGRAPE cell-index elapsed time
  real(8),save:: mdgboard_time_elps ! MDGRAPE board elapsed time
  real(8),save:: mdgexcs_time_elps  ! MDGRAPE excs elapsed time
  real(8),save:: mdgsum_time_elps   ! MDGRAPE sum elapsed time

  real(8),save:: heatf_time_elps   ! heatflux sum elapsed time

  !---- time ratio
  real(8),save:: dyn_ratio         ! dynamics time ratio
  real(8),save:: int_ratio         ! interaction time ratio

  real(8),save:: comm_ratio        ! communication time ratio

  real(8),save:: out_ratio         ! output time ratio

  real(8),save:: table_ratio       ! table time ratio

  real(8),save:: bond_ratio        ! bond time ratio
  real(8),save:: drct_ratio        ! direct sum time ratio
  real(8),save:: rcpr_ratio        ! reciprocal time ratio
  real(8),save:: mor_ratio         ! morse time ratio
  real(8),save:: sh_ratio          ! SH time ratio
  real(8),save:: rfh_ratio         ! RFH time ratio
  real(8),save:: dou_ratio         ! DOU time ratio
  real(8),save:: cnpvw_ratio       ! CNP_VW time ratio
  real(8),save:: cstmnb_ratio       ! CSTMNB time ratio

  real(8),save:: mdgcell_ratio     ! MDGRAPE cell-index time ratio
  real(8),save:: mdgboard_ratio    ! MDGRAPE board time ratio
  real(8),save:: mdgexcs_ratio     ! MDGRAPE excs time ratio
  real(8),save:: mdgsum_ratio      ! MDGRAPE sum time ratio

  real(8),save:: heatf_ratio      ! heatflux sum time ratio

end module time_global
