!*****************************
!*  mklist_mor.f Ver.1.3     *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-01-22 00:49:19 gota>

subroutine mklist2a_mor( rcut_bookmor, &
     &                   xcel,ycel,zcel, &
     &                   ifbookmor)

  use md_global
  use mpi_global

  implicit none

!    subroutine to generate ATOM-BASED Morse bonded ATOM LIST   
!
! ARGUMENT:
!     INPUT
  real(8),intent(in):: rcut_bookmor     ! cut off radius of bookkeeping of Morse

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: ifbookmor      ! flag for bookkeeping of Morse interaction

!        ATTENTION: list_excl and index_ecxl should have been 
!                   reconstructed so that only i < j are contained.
        
!               i.e.  list for atom(i) begins at nb_indexall(i) of nb_listall
!        If atom(i) does not have any neighbor, then 0 is input in nb_listall.
!
! MPI version by G.Kikugawa

! LOCAL:

!---- for general use --

  real*8:: sqcut            ! rcut2**2
  real*8:: rij(3), dr
                                ! cooridnate of cent of res i & j 
                                ! difference, etc.

  integer:: list_tmp_i(maxnatom)
                                !  temporal storage for list for atom(i)

  integer:: icut            ! current size of lists  

  real*8:: box(3)           ! BOX size
  real*8:: box_inv(3)       ! inverse of BOX size

!     -- misc --

  integer:: i, j            ! do loop index for atoms
  integer:: itmp, jtmp,iitmp
  integer:: mor1, mor2
  integer:: itype, jtype

!      logical:: ifmorsepair
!      logical:: ifmorse_former, ifmorse_latter

! MPI
  integer:: nlistmorall_lcl ! local variable nlistmorall

!     +     +     +     +     +     +     +

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1) = 1.0d0/xcel
  box_inv(2) = 1.0d0/ycel
  box_inv(3) = 1.0d0/zcel

  sqcut = rcut_bookmor * rcut_bookmor

  if (.not. ifbookmor) then          ! if bookkeeping is not used
     sqcut =  (xcel/2.0d0)*(xcel/2.0d0) &
          & + (ycel/2.0d0)*(ycel/2.0d0) &
          & + (zcel/2.0d0)*(zcel/2.0d0)
  end if

!---- Initialize nb_list or other ----

  do i = 1, maxnatom+maxnproc
     mor_indexall(i) = 0
  end do

  do i = 1, maxnatom*maxmorlist / nproc
     mor_listall(i) = 0
  end do


!---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!     - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, natom

  icut = 0
  nlistmorall_lcl = 0

  looplast = nmorse - 1

  DO mor1 = loopinit, looplast, loopstep ! loop over nmorselist(mor1)
     i = nmorselist(mor1)
     itype = atmindex(i)

     itmp = 0               ! temporal number of neighboring atoms
                                ! of atom(i)
            
     DO mor2 = mor1+1, nmorse     ! loop over nmorselist(mor2)
        j = nmorselist(mor2)
        jtype = atmindex(j)

!----       if not morse interaction type, skip the loop

        if (inter_inttyp(itype,jtype) /= INTTYPE_MOR) cycle

!       - CHECK IF ATOM(i) is close to ATOM(j) -

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,j)
 
!           - periodic boundary -
            
        rij(1:3) = rij(1:3) - box(1:3) * &
             &                dnint(rij(1:3)*box_inv(1:3))

!                 * quite large box has been set in case of
!                   NON-periodicity to avoid imaging

        dr = rij(1)**2 + rij(2)**2 + rij(3)**2

        if (dr < sqcut) then

           itmp = itmp + 1
           list_tmp_i(itmp) = j

        end if

     END DO                 ! end of loop over nmorselist(mor2)
                  

     iitmp = itmp           ! = itmp - excluded atoms

!        - transfer data of list_tmp_i to mor_list

     mor_indexall(mor1) = icut + 1

     IF (iitmp == 0) then

!        - atom(i) has no neighbor but still 0 is put to mor_list -

        mor_listall(icut+1) = 0
        icut = icut + 1
        nlistmorall_lcl = nlistmorall_lcl + 1

     ELSE

!           - atom(i) has at least one neighbor -

        jtmp = 0
        do j = 1, itmp
           if (list_tmp_i(j) > 0) then
              jtmp = jtmp + 1
              mor_listall(icut+jtmp) = list_tmp_i(j)
           end if
        end do
        icut = icut + iitmp
        nlistmorall_lcl = nlistmorall_lcl + iitmp

     END IF

  END DO                    ! end of loop over nmorselist(mor1)


!     - end of loop -

  nlistmorall = icut
!!! non-existence last loop index mor1
  mor_indexall(mor1) = nlistmorall + 1

#if defined(MPI)
  call mpi_reduce(nlistmorall_lcl,nlistmorall,1,MPI_INTEGER, &
       &                MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

  if (irank == 0) then
     write(6,999) nlistmorall
  end if

  if (maxmorlist*maxnatom / nproc < nlistmorall_lcl) then
     write(6,*) 'Error: nlistmorall(local) exceed number', &
             &     ' of Morse list at process: ',irank
     stop
  end if

999 format(/3x,'Morse bonded pairs updated: ', i18) 

!      +      +      +      +      +      +

  return
end subroutine mklist2a_mor
