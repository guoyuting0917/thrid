!*****************************
!*  mklist_dou.f90 Ver.1.1   *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-01-22 00:52:26 gota>

subroutine mklist2a_dou( rcut_bookdou, &
     &                   xcel, ycel, zcel, &
     &                   ifbookdou, &
     &                   inttype_dou)

  use md_global
  use mpi_global

  implicit none

!    subroutine to generate ATOM-BASED DOU interaction ATOM LIST
!
! ARGUMENT:
!     INPUT
  real(8),intent(in):: rcut_bookdou ! cut off radius of bookkeeping[non-d] of DOU

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: ifbookdou       ! flag for bookkeeping of DOU interaction

  integer,intent(in):: inttype_dou     ! flag for DOU interact type
                                ! = 7, INTTYPE_DOUO
                                ! = 8, INTTYPE_DOUH

!        ATTENTION: list_excl and index_excl should have been 
!                   reconstructed so that only i < j are contained.
        
!               i.e.  list for atom(i) begins at nb_indexall(i) of nb_listall
!        If atom(i) does not have any neighbor, then 0 is input in nb_listall.
!
! MPI version by G.Kikugawa

! LOCAL:
!---- for core list array --
  integer:: dou_indexall(maxnatom+maxnproc)
  integer:: dou_listall(maxnatom*maxdoulist/nproc)
  integer:: nlistdouall

!---- for general use --

  real(8):: sqcut            ! rcut_doux**2
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
  integer:: dou1, dou2
  integer:: itype, jtype

! MPI
  integer:: nlistdouall_lcl ! local variable nlistdouxall

!     +     +     +     +     +     +     +

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1) = 1.0d0/xcel
  box_inv(2) = 1.0d0/ycel
  box_inv(3) = 1.0d0/zcel

  sqcut = rcut_bookdou * rcut_bookdou

!     for save the memory usage
#if defined(_DO_NOT_USE_THIS)
  if (.not. ifbookdou) then          ! if bookkeeping is not used
     sqcut =  (xcel/2.0d0)*(xcel/2.0d0) &
          & + (ycel/2.0d0)*(ycel/2.0d0) &
          & + (zcel/2.0d0)*(zcel/2.0d0)
  end if
#endif

!---- Initialize nb_list or other ----

  do i = 1, maxnatom+maxnproc
     dou_indexall(i) = 0
  end do

  do i = 1, maxnatom*maxdoulist / nproc
     dou_listall(i) = 0
  end do

!---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!     - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, natom

  icut = 0
  nlistdouall_lcl = 0

  looplast = ndou - 1

  DO dou1 = loopinit, looplast, loopstep ! loop over ndoulist(dou1)
     i = ndoulist(dou1)
     itype = atmindex(i)

     itmp = 0               ! temporal number of neighboring atoms
                                ! of atom(i)
            
     DO dou2 = dou1+1, ndou ! loop over ndoulist(dou2)
        j = ndoulist(dou2)
        jtype = atmindex(j)

!----       if DOU interaction type
        if (inter_inttyp(itype,jtype) /= inttype_dou) cycle

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

     END DO                 ! end of loop over ndoulist(dou2)
                  

     iitmp = itmp           ! = itmp - excluded atoms

!        - transfer data of list_tmp_i to ndou_list

     dou_indexall(dou1) = icut + 1

     IF (iitmp == 0) then

!        - atom(i) has no neighbor but still 0 is put to dou_list -

        dou_listall(icut+1) = 0
        icut = icut + 1
        nlistdouall_lcl = nlistdouall_lcl + 1

     ELSE

!           - atom(i) has at least one neighbor -

        jtmp = 0
        do j = 1, itmp
           if (list_tmp_i(j) > 0) then
              jtmp = jtmp + 1
              dou_listall(icut+jtmp) = list_tmp_i(j)
           end if
        end do
        icut = icut + iitmp
        nlistdouall_lcl = nlistdouall_lcl + iitmp

     END IF

  END DO                    ! end of loop over ndoulist(dou1)


!     - end of loop -

  nlistdouall = icut
!!! non-existence last loop index dou1
  dou_indexall(dou1) = nlistdouall + 1

!     - store to selected array
  if (inttype_dou == INTTYPE_DOUO) then
     do i = 1, ndou+nproc
        douo_indexall(i) = dou_indexall(i)
     end do

     do i = 1, nlistdouall
        douo_listall(i) =  dou_listall(i)
     end do

     nlistdouoall = nlistdouall

  else if (inttype_dou == INTTYPE_DOUH) then
     do i = 1, ndou+nproc
        douh_indexall(i) = dou_indexall(i)
     end do

     do i = 1, nlistdouall
        douh_listall(i) =  dou_listall(i)
     end do

     nlistdouhall = nlistdouall

  endif

#if defined(MPI)
  call mpi_reduce( nlistdouall_lcl, nlistdouall, 1, MPI_INTEGER, &
       &           MPI_SUM, 0, MPI_COMM_WORLD, ierror )
#endif

  if (inttype_dou == INTTYPE_DOUO) then
     if (irank == 0) then
        write(6,999) nlistdouall
     end if

     if (maxdoulist*maxnatom / nproc < nlistdouall_lcl) then
        write(6,*) 'Error: nlistdouoall(local) exceed number', &
             &     ' of DOU list at process ',irank
        stop
     end if
  else if (inttype_dou == INTTYPE_DOUH) then
     if (irank == 0) then
        write(6,998) nlistdouall
     end if

     if (maxdoulist*maxnatom / nproc < nlistdouall_lcl) then
        write(6,*) 'Error: nlistdouhall(local) exceed number', &
             &     ' of DOU list at process ',irank
        stop
     end if
  end if

999 format(/3x,'DOU(O) nonbonded pairs updated: ', i18)
998 format(/3x,'DOU(H) nonbonded pairs updated: ', i18)

!      +      +      +      +      +      +

end subroutine mklist2a_dou
