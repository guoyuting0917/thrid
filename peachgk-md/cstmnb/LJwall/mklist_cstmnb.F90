!*****************************
!*  mklist_cstmnb.f Ver.1.1  *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-01-28 12:40:47 gota>

!!! Virtual wall (continuous media) interaction by LJ potential
!!!   originally coded by Jo Suzuki, IFS, Tohoku University

subroutine mklist2a_cstmnb(ifcellindex_cstmnb,ifbookcstmnb, &
     &                     xcel,ycel,zcel, &
     &                     current_step)

  use md_global
  use mpi_global

  use cstmnb

  implicit none

!    subroutine to generate ATOM-BASED custom NB ATOM LIST
!
!    ARGUMENT:
!    INPUT
  logical,intent(in):: ifcellindex_cstmnb     ! flag for cell index (custom NB)
  logical,intent(in):: ifbookcstmnb           ! flag for bookkeeping of custom NB interaction

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  integer,intent(in):: current_step     ! current MD step 

! MPI version by G.Kikugawa

! LOCAL:
  real(8):: ziw                  ! zi-zw
!---- for general use --

!     -- misc --

  integer:: i, k            ! do loop index for atoms
  integer:: itmp
  integer:: cstmnbk


!     +     +     +     +     +     +     +


!---- Generate NB list for Matsui interaction

  IF ((mod(current_step,nstep_bookWALL) == 0) .or. &
       &  (nstep_bookWALL == 1)) THEN

     if (.not. ifbookcstmnb) then          ! if bookkeeping is not used
        zcut_bookWALL =  zcel * 0.5d0
     end if

!---- Initialize nb_list or other ----

     cstmnb_listall(1:maxnatom,1:2) = 0

!--- ordinary bookkeeping

!---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!       - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!

     nlistcstmnball_lcl(1:2) = 0

!    --- wall 1 ---
     looplast = ncstmnb
     itmp = 0

     DO cstmnbk = loopinit, looplast, loopstep   ! loop over ncstmnblist
 
        i = ncstmnblist(cstmnbk)

!       - CHECK IF ATOM(i) is close to wall-1 -
        ziw = atmcor(3,i)       ! distance between wall-k and particles (ziw > 0)

        if (ziw > zcut_bookWALL) cycle
        itmp = itmp + 1
        cstmnb_listall(itmp,1) = i

     END DO                      ! end of loop over ncstmnblist

     nlistcstmnball_lcl(1) = itmp


!    --- wall 2 ---
     looplast = ncstmnb
     itmp = 0

     DO cstmnbk = loopinit, looplast, loopstep   ! loop over ncstmnblist

        i = ncstmnblist(cstmnbk)

!       - CHECK IF ATOM(i) is close to wall-2 -
        ziw = dstnc - atmcor(3,i)   ! distance between wall-k and particles (ziw > 0)

        if (ziw > zcut_bookWALL) cycle
        itmp = itmp + 1
        cstmnb_listall(itmp,2) = i

     END DO                 ! end of loop over ncstmnblist

     nlistcstmnball_lcl(2) = itmp

#if defined(MPI)
     call mpi_reduce(nlistcstmnball_lcl,nlistcstmnball,2,MPI_INTEGER, &
          &          MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

     if (irank == 0) then
        write(6,998) nlistcstmnball(1),nlistcstmnball(2)
     end if

     if (maxnatom < nlistcstmnball_lcl(1)) then
        write(6,*) 'Error: nlistcstmnball(1) (local) exceed number', &
             &     ' of custom NB list'
        stop
     end if
     if (maxnatom < nlistcstmnball_lcl(2)) then
        write(6,*) 'Error: nlistcstmnball(2) (local) exceed number', &
             &     ' of custom NB list'
        stop
     end if

998  format(/3x,'Custom nonbonded pairs updated. ','Wall-1:',i4,'Wall-2:',i4)

  END IF

!      +      +      +      +      +      +

end subroutine mklist2a_cstmnb
