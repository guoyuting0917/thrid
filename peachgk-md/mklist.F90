!*****************************
!*  mklist.f Ver.2.5         *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-01-21 23:13:49 gota>

subroutine mklist2a( rcut_book, &
     &               xcel,ycel,zcel, &
     &               ifbook )

  use md_global
  use mpi_global

  implicit none

!    subroutine to generate ATOM-BASED nonbonded ATOM LIST   
!    intended for use with EWALD option
!
! ARGUMENT:
!     INPUT
  real(8),intent(in):: rcut_book        ! cut off radius of bookkeeping

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: ifbook          ! flag for bookkeeping

!        ATTENTION: list_excl and index_ecxl should have been 
!                   reconstructed so that only i < j are contained.
        
!               i.e.  list for atom(i) begins at nb_indexall(i) of nb_listall
!        If atom(i) does not have any neighbor, then 0 is input in nb_listall.
!
! MPI version by G.Kikugawa

! LOCAL:

!---- for general use --

  real(8):: sqcut            ! rcut2**2
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
  integer:: ix, ix1, ix2
  integer:: itmp, jtmp,iitmp

  integer:: itype, jtype

!      logical:: ifmorsepair
!      logical:: ifmorse_former, ifmorse_latter

! MPI
  integer:: nlistnball_lcl ! local variable nlistnball

!     +     +     +     +     +     +     +

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1) = 1.0d0/xcel
  box_inv(2) = 1.0d0/ycel
  box_inv(3) = 1.0d0/zcel

  sqcut = rcut_book * rcut_book

  if (.not. ifbook) then          ! if bookkeeping is not used
     sqcut =  (xcel/2.0d0)*(xcel/2.0d0) &
          & + (ycel/2.0d0)*(ycel/2.0d0) &
          & + (zcel/2.0d0)*(zcel/2.0d0)
  end if

!---- Initialize nb_list or other ----

  do i = 1, maxnatom+maxnproc
     nb_indexall(i) = 0
  end do

  do i = 1, maxnatom*maxnlist/nproc
     nb_listall(i) = 0
  end do


!---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!     - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, natom

  icut = 0
  nlistnball_lcl = 0

  looplast = natom-1

  DO i = loopinit, looplast, loopstep ! loop over atom(i)

     itmp = 0               ! temporal number of neighboring atoms
                                ! of atom(i)
     itype = atmindex(i)

     DO j = i+1, natom          ! loop over atom(j)

        jtype = atmindex(j)

#if defined(_DO_NOT_USE_THIS)
!----       if morse interaction type, skip the loop
        if (inter_inttyp(itype,jtype) == INTTYPE_MOR) cycle

!----       if SH interaction type, skip the loop
        if (inter_inttyp(itype,jtype) == INTTYPE_SH) cycle

!----       if RFH interaction type, skip the loop
        if (inter_inttyp(itype,jtype) == INTTYPE_RFH) cycle
#endif
!            if (inter_inttyp(itype,jtype) /= INTTYPE_VDW) cycle

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

     END DO                 ! end of loop over atom(j)
                  

!        - delete excluded atoms from list_tmp_i -

     iitmp = itmp           ! = itmp - excluded atoms

     ix1 = index_excl(i)    ! the beginning and    
     ix2 = index_excl(i+1) - 1 ! the end of excluded list
                                ! for atom(i)

     DO j = 1, itmp         ! loop over neighbors of atom (i)

        DO ix = ix1, ix2    ! loop over excluded of atom(i)
           if (list_tmp_i(j) == list_excl(ix)) then
              list_tmp_i(j) = -999
              iitmp = iitmp - 1
              EXIT
           end if
        END DO

     END DO

!        - transfer data of list_tmp_i to nb_list

     nb_indexall(i) = icut + 1

     IF (iitmp == 0) then

!        - atom(i) has no neighbor but still 0 is put to nb_list -

        nb_listall(icut+1) = 0
        icut = icut + 1
        nlistnball_lcl = nlistnball_lcl + 1

     ELSE

!           - atom(i) has at least one neighbor -

        jtmp = 0
        do j = 1, itmp
           if (list_tmp_i(j) > 0) then
              jtmp = jtmp + 1
              nb_listall(icut+jtmp) = list_tmp_i(j)
           end if
        end do
        icut = icut + iitmp
        nlistnball_lcl = nlistnball_lcl + iitmp

     END IF

  END DO                    ! end of loop over atom(i)

!     - end of loop -

  nlistnball = icut
!!! non-existence last loop index i
  nb_indexall(i) = nlistnball + 1

! MPI sum up the nlistnball
#if defined(MPI)
  call mpi_reduce( nlistnball_lcl,nlistnball,1,MPI_INTEGER, &
       &           MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

  if (irank == 0) then
     write(6,999) nlistnball
  end if

! check the array bound of nb_listall in the local process
  if (maxnlist*maxnatom / nproc < nlistnball_lcl) then
     write(6,*) 'Error: nlistnball(local) exceed number', &
          &     ' of nonbonded list at process: ',irank
     stop
  end if

999 format(/3x,'vdW and Coulomb nonbonded pairs updated: ', i18)

!      +      +      +      +      +      +

end subroutine mklist2a
