!********************************
!*  mklist_sh_cell.f90 Ver.1.2  *
!*      for peachgk_md.f        *
!*            by G.Kikugawa     *
!********************************
! Time-stamp: <2015-01-22 00:50:56 gota>

subroutine mklist_sh_cell( rcut_booksh, &
     &                     xcel, ycel, zcel )

  use md_global
  use mpi_global

  implicit none

!    subroutine to generate ATOM-BASED SH interaction ATOM LIST
!
! ARGUMENT:
!     INPUT
  real(8),intent(in):: rcut_booksh     ! cut off radius of bookkeeping of SH

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

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

!---- for cell-index method
  integer:: cellcount(3)    ! number of cell for each direction
  real(8):: cellsize(3)      ! cell length for each direction
  integer:: atom_to_cell(3,natom) ! 3-dimensional cell-index
  integer:: cellindex       ! cell-index

  integer:: cell_natom(0:max_cellindex) ! number of atoms in each cell
  integer:: cell_atomid(max_atomcell,0:max_cellindex)
                                ! atom ID for each cell

  integer:: cellindex_last

  integer:: cell_exist(3)
  integer:: x1,y1,z1
  integer:: x, y, z
  integer:: posyz

!     -- misc --

  integer:: i, j            ! do loop index for atoms
  integer:: jj
  integer:: itmp, jtmp,iitmp
  integer:: sh1
!      integer:: sh2
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

!      if (.not. ifbooksh) then          ! if bookkeeping is not used
!         sqcut =  (xcel/2.0d0)*(xcel/2.0d0)
!     &          + (ycel/2.0d0)*(ycel/2.0d0)
!     &          + (zcel/2.0d0)*(zcel/2.0d0)
!      end if

!---- Initialize nb_list or other ----

  do i = 1, maxnatom+maxnproc
     sho_indexall(i) = 0
     shh_indexall(i) = 0
  end do

  do i = 1, maxnatom*maxshlist / nproc
     sho_listall(i) = 0
     shh_listall(i) = 0
  end do

!---- Calculation cell length and number of cells

!     - cell count      
  cellcount(1) = INT(box(1) / rcut_booksh)
  cellcount(2) = INT(box(2) / rcut_booksh)
  cellcount(3) = INT(box(3) / rcut_booksh)

!      if (cellcount(1) <= 3) then
!         write(6,*) 'Warning: xcel is too small',
!     &              ' to execute cell-index'
!      else if (cellcount(1) > max_cellx) then
  if (cellcount(1) > max_cellx) then
     write(6,*) 'Error: number of cells (x)', &
          &     ' exceeds max num. of cells'
     stop
  end if

!      if (cellcount(2) <= 3) then
!         write(6,*) 'Warning: ycel is too small',
!     &              ' to execute cell-index'
!      else if (cellcount(2) > max_celly) then
  if (cellcount(2) > max_celly) then
     write(6,*) 'Error: number of cells (y)', &
          &     ' exceeds max num. of cells'
     stop
  end if

!      if (cellcount(3) <= 3) then
!         write(6,*) 'Warning: zcel is too small',
!     &              ' to execute cell-index'
!      else if (cellcount(3) > max_cellz) then
  if (cellcount(3) > max_cellz) then
     write(6,*) 'Error: number of cells (z)', &
          &     ' exceeds max num. of cells'
     stop
  end if

  cellindex_last = cellcount(1) - 1  &
       &         + cellcount(1) &
       &         * ((cellcount(2) - 1) &
       &          + cellcount(2) * (cellcount(3) - 1))

!     - cell size
  cellsize(1:3) = box(1:3) / DBLE(cellcount(1:3))

!
! SH (O) list
!

!---- Initialize cell-index variables  ----
  cell_natom(0:cellindex_last) = 0

!---- Cell-index registration

  do sh1 = 1, nsho
     i = nsholist(sh1)

     atom_to_cell(1:3,i) = MOD(INT((atmcor(1:3,i) + box(1:3)) &
          &                    / cellsize(1:3)), cellcount(1:3))

     cellindex = atom_to_cell(1,i) &
          &    + cellcount(1) &
          &    * (atom_to_cell(2,i) &
          &     + cellcount(2) * atom_to_cell(3,i))

     cell_natom(cellindex) = cell_natom(cellindex) + 1
     if (cell_natom(cellindex) > max_atomcell) then
        write(6,*) 'Error: number of atoms in cell exceeds max', &
             &     ' (SH(O))'
        write(6,*) '       cell index= ',atom_to_cell(1,i),',', &
             &     atom_to_cell(2,i),',',atom_to_cell(3,i)
        stop
     end if

     cell_atomid(cell_natom(cellindex),cellindex) = i

  end do

#if defined(_CELL_DEBUG)
  write(6,*) 'DEBUG: cellindex for SH(O) interaction'
  write(6,*) 'Number of cells= ',cellcount(1),'x',cellcount(2),'x', &
       &                         cellcount(3)
  write(6,*) 'Number of cellindex= ',cellindex_last
      
  do i = 0, cellindex_last
     write(6,*) i, cell_natom(i)
  end do
#endif

!---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!     - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, natom

  icut = 0
  nlistshoall_lcl = 0

  looplast = nsho - 1

  DO sh1 = loopinit, looplast, loopstep ! loop over nsholist(sh1)
     i = nsholist(sh1)
     itype = atmindex(i)

     itmp = 0               ! temporal number of neighboring atoms
                            ! of atom(i)
            
!         DO sh2 = sh1+1, nsho     ! loop over nsholist(sh2)
!            j = nsholist(sh2)
!            jtype = atmindex(j)

     cell_exist(1:3) = atom_to_cell(1:3,i)

     do z1 = cell_exist(3)-1, cell_exist(3)+1
        z = MOD(z1 + cellcount(3), cellcount(3))
        do y1 = cell_exist(2)-1, cell_exist(2)+1
           y = MOD(y1 + cellcount(2), cellcount(2))
           posyz = y * cellcount(1) + z * cellcount(1)*cellcount(2)
           do x1 = cell_exist(1)-1, cell_exist(1)+1
              x = MOD(x1 + cellcount(1), cellcount(1))
              cellindex = x + posyz
                  
              do jj = 1, cell_natom(cellindex)

                 j = cell_atomid(jj,cellindex) ! atom(j)
                 if (i >= j) cycle
                 jtype = atmindex(j)

!----                if not SH interaction type, skip the loop

                 if (inter_inttyp(itype,jtype) /= INTTYPE_SH) cycle

!                    - CHECK IF ATOM(i) is close to ATOM(j) -

                 rij(1:3) = atmcor(1:3,i) - atmcor(1:3,j)
 
!                    - periodic boundary -
            
                 rij(1:3) = rij(1:3) - box(1:3) &
                      &   * dnint(rij(1:3)*box_inv(1:3))

!                 * quite large box has been set in case of
!                   NON-periodicity to avoid imaging

                 dr = rij(1)**2 + rij(2)**2 + rij(3)**2

                 if (dr <= sqcut) then

                    itmp = itmp + 1
                    list_tmp_i(itmp) = j

                 end if

!         END DO                 ! end of loop over nsholist(sh2)
              end do

           end do
        end do
     end do


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
  call mpi_reduce( nlistshoall_lcl, nlistshoall,1, MPI_INTEGER, &
       &           MPI_SUM, 0, MPI_COMM_WORLD, ierror )
#endif

  if (irank == 0) then
     write(6,999) nlistshoall
  end if

  if (maxshlist*maxnatom / nproc < nlistshoall_lcl) then
     write(6,*) 'Error: nlistshoall(local) exceed number',&
          &     ' of SH list at process ',irank
     stop
  end if

999 format(/3x,'SH(O) nonbonded pairs updated: ', i18) 


!
! SH (H) list
!

!---- Initialize cell-index variables  ----
  cell_natom(0:cellindex_last) = 0

!---- Cell-index registration

  do sh1 = 1, nshh
     i = nshhlist(sh1)

     atom_to_cell(1:3,i) = MOD(INT((atmcor(1:3,i) + box(1:3)) &
          &                    / cellsize(1:3)), cellcount(1:3))

     cellindex = atom_to_cell(1,i) &
          &    + cellcount(1) &
          &    * (atom_to_cell(2,i) &
          &     + cellcount(2) * atom_to_cell(3,i))

     cell_natom(cellindex) = cell_natom(cellindex) + 1
     if (cell_natom(cellindex) > max_atomcell) then
        write(6,*) 'Error: number of atoms in cell exceeds max', &
             &     ' (SH(H))'
        write(6,*) '       cell index= ',atom_to_cell(1,i),',', &
             &     atom_to_cell(2,i),',',atom_to_cell(3,i)
        stop
     end if

     cell_atomid(cell_natom(cellindex),cellindex) = i

  end do

#if defined(_CELL_DEBUG)
  write(6,*) 'DEBUG: cellindex for SH(H) interaction'
  write(6,*) 'Number of cells= ',cellcount(1),'x',cellcount(2),'x', &
       &                         cellcount(3)
  write(6,*) 'Number of cellindex= ',cellindex_last
      
  do i = 0, cellindex_last
     write(6,*) i, cell_natom(i)
  end do
#endif

!---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!     - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, natom

  icut = 0
  nlistshhall_lcl = 0

  looplast = nshh - 1

  DO sh1 = loopinit, looplast, loopstep ! loop over nshhlist(sh1)
     i = nshhlist(sh1)
     itype = atmindex(i)

     itmp = 0               ! temporal number of neighboring atoms
                            ! of atom(i)
            
!         DO sh2 = sh1+1, nshh     ! loop over nshhlist(sh2)
!            j = nshhlist(sh2)
!            jtype = atmindex(j)

     cell_exist(1:3) = atom_to_cell(1:3,i)

     do z1 = cell_exist(3)-1, cell_exist(3)+1
        z = MOD(z1 + cellcount(3), cellcount(3))
        do y1 = cell_exist(2)-1, cell_exist(2)+1
           y = MOD(y1 + cellcount(2), cellcount(2))
           posyz = y * cellcount(1) + z * cellcount(1)*cellcount(2)
           do x1 = cell_exist(1)-1, cell_exist(1)+1
              x = MOD(x1 + cellcount(1), cellcount(1))
              cellindex = x + posyz
                  
              do jj = 1, cell_natom(cellindex)

                 j = cell_atomid(jj,cellindex) ! atom(j)
                 if (i >= j) cycle
                 jtype = atmindex(j)

!----                if not SH interaction type, skip the loop

                 if (inter_inttyp(itype,jtype) /= INTTYPE_SH) cycle

!                    - CHECK IF ATOM(i) is close to ATOM(j) -

                 rij(1:3) = atmcor(1:3,i) - atmcor(1:3,j)

!                    - periodic boundary -
            
                 rij(1:3) = rij(1:3) - box(1:3) &
                      &   * dnint(rij(1)*box_inv(1))

!                 * quite large box has been set in case of
!                   NON-periodicity to avoid imaging

                 dr = rij(1)**2 + rij(2)**2 + rij(3)**2

                 if (dr <= sqcut) then

                    itmp = itmp + 1
                    list_tmp_i(itmp) = j

                 end if

!         END DO                 ! end of loop over nshhlist(sh2)
              end do

           end do
        end do
     end do


     iitmp = itmp           ! = itmp - excluded atoms

!        - transfer data of list_tmp_i to shh_list

     shh_indexall(sh1) = icut + 1

     IF (iitmp == 0) then

!        - atom(i) has no neighbor but still 0 is put to shh_list -

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
       &           MPI_SUM, 0, MPI_COMM_WORLD, ierror )
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

end subroutine mklist_sh_cell
