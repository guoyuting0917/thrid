!**************************************
!*  MPI_loopinit.f Ver.1.1 '09.12.10  *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************

subroutine MPI_loopinit(startPROC)

  use mpi_global

  implicit none

! ARGUMENT:
!   INPUT
  integer,intent(in):: startPROC         ! MPI first process ID

! LOCAL:

!     +     +     +     +     +     +     +

!---- initialization
  loopinit = 1
  loopstep = 1
  looplast = 0

!---- distribute loop bound
#if defined(MPI)
  if (irank >= startPROC) then
     loopinit = loopinit + irank - startPROC
     loopstep = nproc - startPROC

  else

!       to change init to not execute loop
!
!        loopstep = 1
!        loopinit = looplast + loopstep
     
  endif
#endif

#if defined(MPI_DEBUG)
  write(6,*) 'irank= ',irank,' ',loopinit,loopstep
#endif

!     +     +     +     +     +     +     +

end subroutine MPI_loopinit
