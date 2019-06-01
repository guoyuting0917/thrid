!*****************************
!*  mklist_rfh.f90 Ver.1.1   *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-01-22 00:51:28 gota>

subroutine mklist2a_rfh( rcut_bookrfh, &
     &                   xcel, ycel, zcel, &
     &                   ifbookrfh, &
     &                   inttype_rfh)

  use md_global
  use mpi_global

  implicit none

!    subroutine to generate ATOM-BASED RFH interaction ATOM LIST   
!
! ARGUMENT:
!     INPUT
  real(8),intent(in):: rcut_bookrfh ! cut off radius of bookkeeping[non-d] of RFH

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: ifbookrfh       ! flag for bookkeeping of RFH interaction

  integer,intent(in):: inttype_rfh     ! flag for RFH interact type
                                ! = 4, INTTYPE_RFHFO
                                ! = 5, INTTYPE_RFHOO
                                ! = 6, INTTYPE_RFHOH

!        ATTENTION: list_excl and index_excl should have been 
!                   reconstructed so that only i < j are contained.
        
!               i.e.  list for atom(i) begins at nb_indexall(i) of nb_listall
!        If atom(i) does not have any neighbor, then 0 is input in nb_listall.
!
! MPI version by G.Kikugawa

! LOCAL:
!---- for core list array --
  integer:: rfh_indexall(maxnatom+maxnproc)
  integer:: rfh_listall(maxnatom*maxrfhlist/nproc)
  integer:: nlistrfhall

!---- for general use --

  real(8):: sqcut            ! rcut_rfhxx**2
  real(8):: rij(3), dr
                                ! cooridnate of cent of res i & j 
                                ! difference, etc.

  integer:: list_tmp_i(maxnatom)
                                !  temporal storage for list for atom(i)

  integer:: icut            ! current size of lists

  real(8):: box(3)           ! BOX size
  real(8):: box_inv(3)       ! inverse of BOX size

!     -- misc --

  integer:: i, j            ! do loop index for atoms
  integer:: itmp
  integer:: iitmp
  integer:: jtmp
  integer:: rfh1, rfh2
  integer:: itype, jtype

! MPI
  integer:: nlistrfhall_lcl ! local variable nlistrfhxxall

!     +     +     +     +     +     +     +

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1) = 1.0d0/xcel
  box_inv(2) = 1.0d0/ycel
  box_inv(3) = 1.0d0/zcel

  sqcut = rcut_bookrfh * rcut_bookrfh

!     for save the memory usage
#if defined(_DO_NOT_USE_THIS)
  if (.not. ifbookrfh) then          ! if bookkeeping is not used
     sqcut =  (xcel/2.0d0)*(xcel/2.0d0) &
          & + (ycel/2.0d0)*(ycel/2.0d0) &
          & + (zcel/2.0d0)*(zcel/2.0d0)
  end if
#endif

!---- Initialize nb_list or other ----

  do i = 1, maxnatom+maxnproc
     rfh_indexall(i) = 0
  end do

  do i = 1, maxnatom*maxrfhlist / nproc
     rfh_listall(i) = 0
  end do

!---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!     - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, natom

  icut = 0
  nlistrfhall_lcl = 0

  looplast = nrfh - 1

  DO rfh1 = loopinit, looplast, loopstep ! loop over nrfhlist(rfh1)
     i = nrfhlist(rfh1)
     itype = atmindex(i)

     itmp = 0               ! temporal number of neighboring atoms
                                ! of atom(i)
            
     DO rfh2 = rfh1+1, nrfh     ! loop over nrfhlist(rfh2)
        j = nrfhlist(rfh2)
        jtype = atmindex(j)

!----       if RFH(FO) interaction type
        if (inter_inttyp(itype,jtype) /= inttype_rfh) cycle

!           - CHECK IF ATOM(i) is close to ATOM(j) -

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,j)
 
!           - periodic boundary -
            
        rij(1:3) = rij(1:3) - box(1:3) * &
             &     dnint(rij(1:3)*box_inv(1:3))

!                 * quite large box has been set in case of
!                   NON-periodicity to avoid imaging

        dr = rij(1)**2 + rij(2)**2 + rij(3)**2

        if (dr <= sqcut) then

           itmp = itmp + 1
           list_tmp_i(itmp) = j

        end if

     END DO                 ! end of loop over nrfhlist(rfh2)
                  

     iitmp = itmp           ! = itmp - excluded atoms

!        - transfer data of list_tmp_i to nrfh_list

     rfh_indexall(rfh1) = icut + 1

     IF (iitmp == 0) then

!        - atom(i) has no neighbor but still 0 is put to rfh_list -

        rfh_listall(icut+1) = 0
        icut = icut + 1
        nlistrfhall_lcl = nlistrfhall_lcl + 1

     ELSE

!           - atom(i) has at least one neighbor -

        jtmp = 0
        do j = 1, itmp
           if (list_tmp_i(j) > 0) then
              jtmp = jtmp + 1
              rfh_listall(icut+jtmp) = list_tmp_i(j)
           end if
        end do
        icut = icut + iitmp
        nlistrfhall_lcl = nlistrfhall_lcl + iitmp

     END IF

  END DO                    ! end of loop over nrfhlist(rfh1)


!     - end of loop -

  nlistrfhall = icut
!!! non-existence last loop index rfh1
  rfh_indexall(rfh1) = nlistrfhall + 1

!     - store to selected array
  if (inttype_rfh == INTTYPE_RFHFO) then
     do i = 1, nrfh+nproc
        rfhfo_indexall(i) = rfh_indexall(i)
     end do

     do i = 1, nlistrfhall
        rfhfo_listall(i) =  rfh_listall(i)
     end do

     nlistrfhfoall = nlistrfhall

  else if (inttype_rfh == INTTYPE_RFHOO) then
     do i = 1, nrfh+nproc
        rfhoo_indexall(i) = rfh_indexall(i)
     end do

     do i = 1, nlistrfhall
        rfhoo_listall(i) =  rfh_listall(i)
     end do

     nlistrfhooall = nlistrfhall

  else if (inttype_rfh == INTTYPE_RFHOH) then
     do i = 1, nrfh+nproc
        rfhoh_indexall(i) = rfh_indexall(i)
     end do

     do i = 1, nlistrfhall
        rfhoh_listall(i) =  rfh_listall(i)
     end do

     nlistrfhohall = nlistrfhall

  endif

#if defined(MPI)
  call mpi_reduce( nlistrfhall_lcl, nlistrfhall, 1, MPI_INTEGER, &
       &           MPI_SUM, 0, MPI_COMM_WORLD, ierror)
#endif

  if (inttype_rfh == INTTYPE_RFHFO) then
     if (irank == 0) then
        write(6,999) nlistrfhall
     end if

     if (maxrfhlist*maxnatom / nproc < nlistrfhall_lcl) then
        write(6,*) 'Error: nlistrfhfoall(local) exceed number', &
             &     ' of RFH list at process ',irank
        stop
     end if

  else if (inttype_rfh == INTTYPE_RFHOO) then
     if (irank == 0) then
        write(6,998) nlistrfhall
     end if

     if (maxrfhlist*maxnatom / nproc < nlistrfhall_lcl) then
        write(6,*) 'Error: nlistrfhooall(local) exceed number', &
             &     ' of RFH list at process ',irank
        stop
     end if

  else if (inttype_rfh == INTTYPE_RFHOH) then
     if (irank == 0) then
        write(6,997) nlistrfhall
     end if

     if (maxrfhlist*maxnatom / nproc < nlistrfhall_lcl) then
        write(6,*) 'Error: nlistrfhohall(local) exceed number', &
             &     ' of RFH list at process ',irank
        stop
     end if
  end if

999 format(/3x,'RFH(FO) nonbonded pairs updated: ', i18)
998 format(/3x,'RFH(OO) nonbonded pairs updated: ', i18)
997 format(/3x,'RFH(OH) nonbonded pairs updated: ', i18)

!      +      +      +      +      +      +

end subroutine mklist2a_rfh
