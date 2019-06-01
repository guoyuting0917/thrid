!*****************************
!*  mklist_sh.f90 Ver.1.1    *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-01-22 00:50:30 gota>

subroutine mklist2a_sh( rcut_booksh, &
     &                  xcel,ycel,zcel, &
     &                  ifbooksh)

  use md_global
  use mpi_global

  implicit none

!    subroutine to generate ATOM-BASED SH interaction ATOM LIST   
!
! ARGUMENT:
!     INPUT
  real(8),intent(in):: rcut_booksh      ! cut off radius of bookkeeping of SH

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: ifbooksh        ! flag for bookkeeping of SH interaction

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
  integer:: itmp, jtmp,iitmp
  integer:: sh1, sh2
  integer:: itype, jtype

! MPI
  integer:: nlistshoall_lcl ! local variable nlistshoall
  integer:: nlistshhall_lcl ! local variable nlistshhall

!     +     +     +     +     +     +     +

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1) = 1.0d0/xcel
  box_inv(2) = 1.0d0/ycel
  box_inv(3) = 1.0d0/zcel

  sqcut = rcut_booksh * rcut_booksh

!     for save the memory usage
#if defined(_DO_NOT_USE_THIS)
  if (.not. ifbooksh) then          ! if bookkeeping is not used
     sqcut =  (xcel/2.0d0)*(xcel/2.0d0) &
          & + (ycel/2.0d0)*(ycel/2.0d0) &
          & + (zcel/2.0d0)*(zcel/2.0d0)
  end if
#endif

!---- Initialize nb_list or other ----

  do i = 1, maxnatom+maxnproc
     sho_indexall(i) = 0
     shh_indexall(i) = 0
  end do

  do i = 1, maxnatom*maxshlist / nproc
     sho_listall(i) = 0
     shh_listall(i) = 0
  end do


!---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!     - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, natom

!---- SH O atoms
  icut = 0
  nlistshoall_lcl = 0

  looplast = nsho - 1

  DO sh1 = loopinit, looplast, loopstep ! loop over nsholist(sh1)
     i = nsholist(sh1)
     itype = atmindex(i)

     itmp = 0               ! temporal number of neighboring atoms
                                ! of atom(i)
            
     DO sh2 = sh1+1, nsho     ! loop over nsholist(sh2)
        j = nsholist(sh2)
        jtype = atmindex(j)

!----       if not SH interaction type, skip the loop

        if (inter_inttyp(itype,jtype) /= INTTYPE_SH) cycle

!       - CHECK IF ATOM(i) is close to ATOM(j) -

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

     END DO                 ! end of loop over nsholist(sh2)
                  

     iitmp = itmp           ! = itmp - excluded atoms

!        - transfer data of list_tmp_i to sho_list

     sho_indexall(sh1) = icut + 1

     IF (iitmp == 0) then

!        - atom(i) has no neighbor but still 0 is put to sho_list -

        sho_listall(icut+1) = 0
        icut = icut + 1
        nlistshoall_lcl = nlistshoall_lcl + 1

     ELSE

!           - atom(i) has at least one neighbor -

        jtmp = 0
        do j = 1, itmp
           if (list_tmp_i(j) > 0) then
              jtmp = jtmp + 1
              sho_listall(icut+jtmp) = list_tmp_i(j)
           end if
        end do
        icut = icut + iitmp
        nlistshoall_lcl = nlistshoall_lcl + iitmp

     END IF

  END DO                    ! end of loop over nsholist(sh1)


!     - end of loop -

  nlistshoall = icut
!!! non-existence last loop index sh1
  sho_indexall(sh1) = nlistshoall + 1

#if defined(MPI)
  call mpi_reduce( nlistshoall_lcl, nlistshoall, 1, MPI_INTEGER, &
       &           MPI_SUM, 0, MPI_COMM_WORLD, ierror )
#endif

  if (irank == 0) then
     write(6,999) nlistshoall
  end if

  if (maxshlist*maxnatom / nproc < nlistshoall_lcl) then
     write(6,*) 'Error: nlistshoall(local) exceed number', &
          &     ' of SH list at process ',irank
     stop
  end if

999 format(/3x,'SH(O) nonbonded pairs updated: ', i18)


!----- SH H atoms
  icut = 0
  nlistshhall_lcl = 0

  looplast = nshh - 1

  DO sh1 = loopinit, looplast, loopstep ! loop over nshhlist(sh1)
     i = nshhlist(sh1)
     itype = atmindex(i)

     itmp = 0               ! temporal number of neighboring atoms
                                ! of atom(i)
            
     DO sh2 = sh1+1, nshh     ! loop over nshhlist(sh2)
        j = nshhlist(sh2)
        jtype = atmindex(j)

!----       if not SH interaction type, skip the loop

        if (inter_inttyp(itype,jtype) /= INTTYPE_SH) cycle

!       - CHECK IF ATOM(i) is close to ATOM(j) -

        rij(1:3) = atmcor(1:3,i) - atmcor(1:3,j)

!           - periodic boundary -
            
        rij(1:3) = rij(1:3) - box(1:3) * &
             &     dnint(rij(1:3)*box_inv(1:3))

!                 * quite large box has been set in case of
!                   NON-periodicity to avoid imaging

        dr = rij(1)**2 + rij(2)**2 + rij(3)**2

        if (dr < sqcut) then

           itmp = itmp + 1
           list_tmp_i(itmp) = j

        end if

     END DO                 ! end of loop over nshhlist(sh2)
                  

     iitmp = itmp           ! = itmp - excluded atoms

!        - transfer data of list_tmp_i to shh_list

     shh_indexall(sh1) = icut + 1

     IF (iitmp == 0) then

!        - atom(i) has no neighbor but still 0 is put to sh_list -

        shh_listall(icut+1) = 0
        icut = icut + 1
        nlistshhall_lcl = nlistshhall_lcl + 1

     ELSE

!           - atom(i) has at least one neighbor -

        jtmp = 0
        do j = 1, itmp
           if (list_tmp_i(j) > 0) then
              jtmp = jtmp + 1
              shh_listall(icut+jtmp) = list_tmp_i(j)
           end if
        end do
        icut = icut + iitmp
        nlistshhall_lcl = nlistshhall_lcl + iitmp

     END IF

  END DO                    ! end of loop over nshhlist(sh1)


!     - end of loop -

  nlistshhall = icut
!!! non-existence last loop index sh1
  shh_indexall(sh1) = nlistshhall + 1

#if defined(MPI)
  call mpi_reduce( nlistshhall_lcl, nlistshhall, 1, MPI_INTEGER, &
       &           MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

  if (irank == 0) then
     write(6,998) nlistshhall
  end if

  if (maxshlist*maxnatom / nproc < nlistshhall_lcl) then
     write(6,*) 'Error: nlistshhall(local) exceed number', &
          &     ' of SH list at process ',irank
     stop
  end if

998 format(/3x,'SH(H) nonbonded pairs updated: ', i18)

!      +      +      +      +      +      +

end subroutine mklist2a_sh
