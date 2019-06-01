!*****************************
!*  mklist_cstmnb.f Ver.1.2  *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-01-22 01:06:51 gota>

subroutine mklist2a_cstmnb(ifcellindex_cstmnb,ifbookcstmnb, &
     &                     xcel,ycel,zcel, &
     &                     current_step)

  use md_global
  use mpi_global

  use cstmnb

  implicit none

!    subroutine to generate ATOM-BASED custom NB ATOM LIST
!
! ARGUMENT:
!     INPUT
  logical,intent(in):: ifcellindex_cstmnb ! flag for cell index (custom NB)
  logical,intent(in):: ifbookcstmnb       
                                ! flag for bookkeeping of custom NB interaction

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  integer,intent(in):: current_step     ! current MD step 

!        ATTENTION: list_excl and index_excl should have been
!                   reconstructed so that only i < j are contained.
!
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

!--- for cell-index method
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
  integer:: jj              ! do loop index for atoms
  integer:: ix, ix1, ix2
  integer:: itmp, jtmp,iitmp
  integer:: itype, jtype

! MPI
  integer:: nlistcstmnball_lcl ! local variable nlistcstmnball

!     +     +     +     +     +     +     +

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel
  box_inv(1) = 1.0d0/xcel
  box_inv(2) = 1.0d0/ycel
  box_inv(3) = 1.0d0/zcel

!---- Generate NB list for cstmnb interaction

  IF ((mod(current_step,nstep_bookhfewr) == 0) .or. &
       & (nstep_bookhfewr == 1)) THEN

     sqcut = rrcutbook_hfewr * rrcutbook_hfewr

     if (.not. ifbookcstmnb) then          ! if bookkeeping is not used
        sqcut =  (xcel/2.0d0)*(xcel/2.0d0) &
             & + (ycel/2.0d0)*(ycel/2.0d0) &
             & + (zcel/2.0d0)*(zcel/2.0d0)
     end if

!---- Initialize nb_list or other ----

     do i = 1, maxnatom+maxnproc
        cstmnb_indexall(i) = 0
     end do

     do i = 1, maxnatom*maxcstmnblist / nproc
        cstmnb_listall(i) = 0
     end do

!    --- ordinary bookkeeping
     IF (.not.ifcellindex_cstmnb) THEN

!---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!       - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, natom

        icut = 0
        nlistcstmnball_lcl = 0

        looplast = natom - 1

        DO i = loopinit, looplast, loopstep ! loop over natom
           itype = atmindex(i)

           itmp = 0               ! temporal number of neighboring atoms
                                  ! of atom(i)

           DO j = i+1, natom     ! loop over natom
              jtype = atmindex(j)

! !----         if not custom NB interaction type, skip the loop

!               if (inter_inttyp(itype,jtype) /= INTTYPE_CSTMNB) cycle

!             - CHECK IF ATOM(i) is close to ATOM(j) -

              rij(1:3) = atmcor(1:3,i) - atmcor(1:3,j)

!             - periodic boundary -

              rij(1:3) = rij(1:3) - box(1:3) * &
                   &                anint(rij(1:3)*box_inv(1:3))

!                 * quite large box has been set in case of
!                   NON-periodicity to avoid imaging

              dr = rij(1)**2 + rij(2)**2 + rij(3)**2

              if (dr < sqcut) then

                 itmp = itmp + 1
                 list_tmp_i(itmp) = j

              end if

           END DO                 ! end of loop over j

!          - delete excluded atoms from list_tmp_i -

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

!          - transfer data of list_tmp_i to cstmnb_list

           cstmnb_indexall(i) = icut + 1

           IF (iitmp == 0) then

!             - atom(i) has no neighbor but still 0 is put to cstmnb_list -

              cstmnb_listall(icut+1) = 0
              icut = icut + 1
              nlistcstmnball_lcl = nlistcstmnball_lcl + 1

           ELSE

!             - atom(i) has at least one neighbor -

              jtmp = 0
              do j = 1, itmp
                 if (list_tmp_i(j) > 0) then
                    jtmp = jtmp + 1
                    cstmnb_listall(icut+jtmp) = list_tmp_i(j)
                 end if
              end do
              icut = icut + iitmp
              nlistcstmnball_lcl = nlistcstmnball_lcl + iitmp

           END IF

        END DO                    ! end of loop over i


!       - end of loop -

        nlistcstmnball = icut
!!! non-existence last loop index cstmnb1
        cstmnb_indexall(i) = nlistcstmnball + 1

#if defined(MPI)
        call mpi_reduce(nlistcstmnball_lcl,nlistcstmnball,1,MPI_INTEGER, &
             &          MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

        if (irank == 0) then
           write(6,998) nlistcstmnball
        end if
        
        if (maxcstmnblist*maxnatom / nproc < nlistcstmnball_lcl) then
           write(6,*) 'Error: nlistcstmnball(local) exceed number', &
                &     ' of custom NB list at process',irank
           stop
        end if

998     format(/3x,'Custom nonbonded pairs updated: ', i18) 

!    --- cell index and bookkeeping
     ELSE

!---    Calculation cell length and number of cells

!       - cell count
        cellcount(1:3) = INT(box(1:3) / rrcutbook_hfewr)

        if (cellcount(1) > max_cellx) then
           write(6,*) 'Error: number of cells (x)', &
                &     ' exceeds max num. of cells'
           stop
        end if

        if (cellcount(2) > max_celly) then
           write(6,*) 'Error: number of cells (y)', &
                &     ' exceeds max num. of cells'
           stop
        end if

        if (cellcount(3) > max_cellz) then
           write(6,*) 'Error: number of cells (z)', &
                &     ' exceeds max num. of cells'
           stop
        end if

        cellindex_last = cellcount(1) - 1 &
             &         + cellcount(1) &
             &         * ((cellcount(2) - 1) &
             &          + cellcount(2) * (cellcount(3) - 1))

!       - cell size
        cellsize(1:3) = box(1:3) / DBLE(cellcount(1:3))

!----   Initialize cell-index variables  ----
        cell_natom(0:cellindex_last) = 0

!----   Cell-index registration

        do i = 1, natom

           atom_to_cell(1:3,i) = MOD(INT((atmcor(1:3,i) + box(1:3)) &
                &                    / cellsize(1:3)), cellcount(1:3))

           cellindex = atom_to_cell(1,i) &
                &    + cellcount(1) &
                &    * (atom_to_cell(2,i) &
                &     + cellcount(2) * atom_to_cell(3,i))

           cell_natom(cellindex) = cell_natom(cellindex) + 1
           if (cell_natom(cellindex) > max_atomcell) then
              write(6,*) 'Error: number of atoms in cell exceeds max', &
                   &     ' (cstmnb)'
              write(6,*) '       cell index= ',atom_to_cell(1,i),',', &
                   &     atom_to_cell(2,i),',',atom_to_cell(3,i)
              stop
           end if

           cell_atomid(cell_natom(cellindex),cellindex) = i

        end do

#if defined(_CELL_DEBUG)
        write(6,*) 'DEBUG: cellindex for custom NB interaction'
        write(6,*) 'Number of cells= ',cellcount(1),'x',cellcount(2),'x', &
             &                         cellcount(3)
        write(6,*) 'Number of cellindex= ',cellindex_last

        do i = 0, cellindex_last
           write(6,*) i, cell_natom(i)
        end do
#endif

!      ---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!     - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!        DO i = 1, natom

        icut = 0
        nlistcstmnball_lcl = 0

        looplast = natom - 1

        DO i = loopinit, looplast, loopstep ! loop over ncstmnblist(cstmnb1)
           itype = atmindex(i)

           itmp = 0               ! temporal number of neighboring atoms
                                  ! of atom(i)

!           DO cstmnb2 = cstmnb1+1, ncstmnb     ! loop over ncstmnblist(cstmnb2)
!              j = ncstmnblist(cstmnb2)
!              jtype = atmindex(j)

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

! !----                  if not custom NB interaction type, skip the loop

!                        if (inter_inttyp(itype,jtype) /= INTTYPE_CSTMNB) cycle

!                      - CHECK IF ATOM(i) is close to ATOM(j) -

                       rij(1:3) = atmcor(1:3,i) - atmcor(1:3,j)

!                      - periodic boundary -

                       rij(1:3) = rij(1:3) - box(1:3) &
                            &   * anint(rij(1:3)*box_inv(1:3))

!                      * quite large box has been set in case of
!                        NON-periodicity to avoid imaging

                       dr = rij(1)**2 + rij(2)**2 + rij(3)**2

                       if (dr < sqcut) then

                          itmp = itmp + 1
                          list_tmp_i(itmp) = j

                       end if

!                    END DO             ! end of loop over ncstmnblist(cstmnb2)
                    end do

                 end do
              end do
           end do

!          - delete excluded atoms from list_tmp_i -

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

!          - transfer data of list_tmp_i to cstmnb_list

           cstmnb_indexall(i) = icut + 1

           IF (iitmp == 0) then

!             - atom(i) has no neighbor but still 0 is put to cstmnb_list -

              cstmnb_listall(icut+1) = 0
              icut = icut + 1
              nlistcstmnball_lcl = nlistcstmnball_lcl + 1

           ELSE

!             - atom(i) has at least one neighbor -

              jtmp = 0
              do j = 1, itmp
                 if (list_tmp_i(j) > 0) then
                    jtmp = jtmp + 1
                    cstmnb_listall(icut+jtmp) = list_tmp_i(j)
                 end if
              end do
              icut = icut + iitmp
              nlistcstmnball_lcl = nlistcstmnball_lcl + iitmp

           END IF

        END DO                    ! end of loop over ncstmnblist(cstmnb1)


!       - end of loop -

        nlistcstmnball = icut
!!! non-existence last loop index cstmnb1
        cstmnb_indexall(i) = nlistcstmnball + 1

#if defined(MPI)
        call mpi_reduce(nlistcstmnball_lcl,nlistcstmnball,1,MPI_INTEGER, &
             &          MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

        if (irank == 0) then
           write(6,999) nlistcstmnball
        end if

        if (maxcstmnblist*maxnatom / nproc < nlistcstmnball_lcl) then
           write(6,*) 'Error: nlistcstmnball(local) exceed number', &
                &     ' of custom NB list at process',irank
           stop
        end if

999     format(/3x,'Custom nonbonded pairs updated: ', i18) 

     END IF

  END IF

!      +      +      +      +      +      +

end subroutine mklist2a_cstmnb
