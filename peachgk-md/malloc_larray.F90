!*******************************
!*  malloc_larray.f90 Ver.1.0  *
!*      for peachgk_md.f       *
!*            by G.Kikugawa    *
!*******************************
! Time-stamp: <2015-01-22 00:14:32 gota>

subroutine malloc_larray()

  use md_global
  use mpi_global

  implicit none

! ARGUMENT:
!   INPUT

! LOCAL:

!     +     +     +     +     +     +     +

!---- dynamic memory allocation for NB lists
  allocate(nb_listall(maxnatom*maxnlist / nproc))
  allocate(mor_listall(maxnatom*maxmorlist / nproc))
  allocate(sho_listall(maxnatom*maxshlist / nproc))
  allocate(shh_listall(maxnatom*maxshlist / nproc))
  allocate(rfhfo_listall(maxnatom*maxrfhlist / nproc))
  allocate(rfhoo_listall(maxnatom*maxrfhlist / nproc))
  allocate(rfhoh_listall(maxnatom*maxrfhlist / nproc))
  allocate(douo_listall(maxnatom*maxshlist / nproc))
  allocate(douh_listall(maxnatom*maxshlist / nproc))
  allocate(rpvw_listall(maxnatom*maxrpvwlist / nproc))

  write(6,*) 'Dynamically allocate memory for non-bond lists at process: ',irank

!     +     +     +     +     +     +     +

end subroutine malloc_larray
