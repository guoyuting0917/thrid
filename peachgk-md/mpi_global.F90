!**************************************
!* mpi_global.inc Ver.1.0 09.10.15 *
!* for peachgk_md.f *
!* by G.Kikugawa *
!**************************************
module mpi_global
  implicit none
  integer,save:: nproc,irank
                       ! nproc : number of concurrent processes
                       ! irank : process number of this process(0 to nproc-1)
  integer,save:: ierror
                       ! ierror : return code of MPI routines
  !---- loop variables
  integer,save:: loopinit,looplast ! do loop start and end point
  integer,save:: loopstep ! do loop step
  !---- PME MPI
  integer,parameter:: maxnproc_fft = 1000
  integer,save:: xlimit1(0:maxnproc_fft), xlimit2(0:maxnproc_fft)
  integer,save:: zlimit1(0:maxnproc_fft), zlimit2(0:maxnproc_fft)
  integer,save:: xlimitmax
  integer,save:: zlimitmax
end module mpi_global
