!*****************************
!*  mklist_cstmnb.f Ver.1.1  *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*      Coded by T.Nakano    *
!*****************************
! Time-stamp: <2015-03-12 11:24:37 gota>

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
  real*8:: rij(3), rik(3), rjk(3)
                                ! cooridnate of cent of res i & j 
                                ! difference, etc.
  real*8:: atmcor_sig(3,natom)  ! non-dimenstioinal atom coodinate
  real*8:: drij, drik, drjk

  integer:: list_tmp_i2(maxnatom)
  integer,allocatable:: list_tmp_i3(:,:)
                                !  temporal storage for list for atom(i)

  integer:: icut2, icut3           ! current size of lists

  real(8):: box(3,ncstmnbtyp)      ! BOX size
  real(8):: box_inv(3,ncstmnbtyp)  ! inverse of BOX size
  real(8):: box_min(3)             ! max BOX size

!--- for cell-index method
  integer:: cellcount(3)    ! number of cell for each direction
  real(8):: cellsize(3)      ! cell length for each direction
  integer:: atom_to_cell(3,natom) ! 3-dimensional cell-index
  integer:: cellindex,cellindex1,cellindex2       ! cell-index

  integer:: cell_natom(0:max_cellindex) ! number of atoms in each cell
  integer:: cell_atomid(max_atomcell,0:max_cellindex)
                                ! atom ID for each cell
  integer:: cellindex_last

  integer:: cell_exist1(3), cell_exist2(3)
  integer:: x1,y1,z1, x1d,y1d,z1d
  integer:: x2,y2,z2, x2d,y2d,z2d
  integer:: posyz1, posyz2

!     -- misc --

  integer:: i, j, k         ! atoms index
  integer:: jj, kk          ! do loop index for atoms
  integer:: j1, j2          ! do loop index for atoms
  integer:: k1, k2          ! do loop index for atoms
  integer:: j_old
  integer:: itmp2, jtmp2,iitmp2
  integer:: itmp3(maxcstmnblist2), iitmp3(maxcstmnblist2)
  integer:: cstmnb1, cstmnb2, cstmnb3
  integer:: itype, jtype, ktype
  integer:: itype_inttyp, jtype_inttyp, ktype_inttyp
  integer:: i2btype, i3btype

! MPI
  integer:: nlistcstmnball_lcl2 ! local variable nlistcstmnball
  integer:: nlistcstmnball_lcl3 ! local variable nlistcstmnball

!     +     +     +     +     +     +     +

!---- dynamic memory allocation for neighbor lists ----

  allocate(list_tmp_i3(maxcstmnblist2,maxcstmnblist))

!---- Generate NB list for SW interaction

  IF ((mod(current_step,nstep_bookSW) == 0) .or. &
       & (nstep_bookSW == 1)) THEN

     sqcut = rcut_bookSW * rcut_bookSW

     if (.not. ifbookcstmnb) then          ! if bookkeeping is not used
        sqcut =  (xcel/2.0d0)*(xcel/2.0d0) &
             & + (ycel/2.0d0)*(ycel/2.0d0) &
             & + (zcel/2.0d0)*(zcel/2.0d0)
     end if

!---- Initialize nb_list or other ----

     cstmnb_indexall2(1:maxnatom+maxnproc) = 0

     cstmnb_listall2(1:maxcstmnblist2) = 0
     cstmnb_indexall3(1:maxcstmnblist2) = 0
     cstmnb_listall3(1:maxcstmnblist2 * maxcstmnblist) = 0

     DO cstmnb1 = 1, ncstmnb
        i = ncstmnblist(cstmnb1)
        itype = atmindex_ncstmnb(i)
        atmcor_sig(1:3,i) = atmcor(1:3,i) / sigma2_SW(itype)

        box(1,itype) = xcel / sigma2_SW(itype)
        box(2,itype) = ycel / sigma2_SW(itype)
        box(3,itype) = zcel / sigma2_SW(itype)
        box_inv(1:3,itype) = 1.0d0/box(1:3,itype)
     END DO

!    --- ordinary bookkeeping
     IF (.not.ifcellindex_cstmnb) THEN

!---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!       - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, natom

        icut2 = 0
        icut3 = 0
        nlistcstmnball_lcl2 = 0
        nlistcstmnball_lcl3 = 0

        looplast = ncstmnb - 1

        DO cstmnb1 = loopinit, looplast, loopstep
           i = ncstmnblist(cstmnb1)
           itype_inttyp = atmindex(i)
           itype = atmindex_ncstmnb(i)

           itmp2 = 0   ! temporal number of neighboring atoms of atom(i)
           itmp3(1:maxcstmnblist2) = 0  ! temporal number of neighboring atoms of atom(j)
           iitmp3(1:maxcstmnblist2) = 0  ! temporal number of neighboring atoms of atom(j)
           j_old = 0

           DO cstmnb2 = cstmnb1+1, ncstmnb     ! loop over ncstmnblist(cstmnb2)
              j = ncstmnblist(cstmnb2)
              jtype_inttyp = atmindex(j)
              jtype = atmindex_ncstmnb(j)

!----         if not custom NB interaction type, skip the loop

!             - CHECK IF ATOM(i) is close to ATOM(j) -

              rij(1:3) = atmcor_sig(1:3,i) - atmcor_sig(1:3,j)

!             - periodic boundary -

              i2btype = cstmnbtypeindex2(itype,jtype)

              rij(1:3) = rij(1:3) - box(1:3,i2btype) * &
                   &                anint(rij(1:3)*box_inv(1:3,i2btype))

!                 * quite large box has been set in case of
!                   NON-periodicity to avoid imaging

              drij = rij(1)**2 + rij(2)**2 + rij(3)**2

              if (drij < sqcut) then

                 itmp2 = itmp2 + 1
                 list_tmp_i2(itmp2) = j

              end if

              DO cstmnb3 = cstmnb2+1, ncstmnb  ! loop over ncstmnblist(cstmnb2)
                 k = ncstmnblist(cstmnb3)
                 ktype_inttyp = atmindex(k)
                 ktype = atmindex_ncstmnb(k)

!----            if not custom NB interaction type, skip the loop

!                 if (inter_inttyp(jtype_inttyp,ktype_inttyp) /= INTTYPE_CSTMNB) cycle
!                 if (inter_inttyp(itype_inttyp,ktype_inttyp) /= INTTYPE_CSTMNB) cycle

!                - CHECK IF two of ATOM(i), ATOM(j) and ATOM(k) are close -

                 i3btype = cstmnbtypeindex3(itype,jtype,ktype)

                 rik(1:3) = atmcor_sig(1:3,i) - atmcor_sig(1:3,k)
                 rik(1:3) = rik(1:3) - box(1:3,i3btype) &
                      &              * anint(rik(1:3)*box_inv(1:3,i3btype))
                 drik = rik(1)**2 + rik(2)**2 + rik(3)**2

                 rjk(1:3) = atmcor_sig(1:3,j) - atmcor_sig(1:3,k)
                 rjk(1:3) = rjk(1:3) - box(1:3,i3btype) &
                      &              * anint(rjk(1:3)*box_inv(1:3,i3btype))
                 drjk = rjk(1)**2 + rjk(2)**2 + rjk(3)**2

                 if (drij < sqcut .and. (drik < sqcut .or. drjk < sqcut)) then
                    itmp3(itmp2) = itmp3(itmp2) + 1
                    list_tmp_i3(itmp2,itmp3(itmp2)) = k
                    cycle
                 end if

                 if (drik < sqcut .and. drjk < sqcut) then

                    ! not to count up 'itmp2' if only this condition matches more than twice for the same j
                    if (j/=j_old) itmp2 = itmp2 + 1
                    j_old = j

                    list_tmp_i2(itmp2) = j
                    itmp3(itmp2) = itmp3(itmp2) + 1
                    list_tmp_i3(itmp2,itmp3(itmp2)) = k
                    cycle
                 end if

              END DO                 ! end of loop over ncstmnblist()

              if (itmp2>0) iitmp3(itmp2) = itmp3(itmp2)           ! = itmp3 - excluded atoms

           END DO                 ! end of loop over ncstmnblist(2)

           iitmp2 = itmp2           ! = itmp2 - excluded atoms

!          - transfer data of list_tmp_i to cstmnb_list

           cstmnb_indexall2(cstmnb1) = icut2 + 1
           IF (iitmp2 == 0) then

!             - atom(i) has no neighbor but still 0 is put to cstmnb_list -

              cstmnb_indexall3(icut2+1) = icut3 + 1

              icut2 = icut2 + 1
              cstmnb_listall2(icut2) = 0
              nlistcstmnball_lcl2 = nlistcstmnball_lcl2 + 1

              icut3 = icut3 + 1
              cstmnb_listall3(icut3) = 0
              nlistcstmnball_lcl3 = nlistcstmnball_lcl3 + 1

           ELSE

!             - atom(i) has at least one neighbor -

              do j = 1, itmp2
                 if (list_tmp_i2(j) > 0) then
                    cstmnb_listall2(icut2+j)  = list_tmp_i2(j)
                 end if

                 cstmnb_indexall3(icut2+j) = icut3 + 1

                 IF (iitmp3(j) == 0) then

                    icut3 = icut3 + 1

                    cstmnb_listall3(icut3) = 0
                    nlistcstmnball_lcl3 = nlistcstmnball_lcl3 + 1

                 ELSE

                    do k = 1, iitmp3(j)
                       if (list_tmp_i3(j,k) > 0) then
                          cstmnb_listall3(icut3+k) = list_tmp_i3(j,k)
                       end if
                    end do
                    icut3 = icut3 + iitmp3(j)
                    nlistcstmnball_lcl3 = nlistcstmnball_lcl3 + iitmp3(j)

                 END IF

              end do

              icut2 = icut2 + iitmp2
              nlistcstmnball_lcl2 = nlistcstmnball_lcl2 + iitmp2
           END IF

        END DO                    ! end of loop over ncstmnblist(cstmnb1)

!       - end of loop -

        nlistcstmnball2 = icut2
!!! non-existent last loop index cstmnb1
        cstmnb_indexall2(cstmnb1) = nlistcstmnball2 + 1

#if defined(MPI)
        call mpi_reduce(nlistcstmnball_lcl2,nlistcstmnball2,1,MPI_INTEGER, &
             &          MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

!        if (irank == 0) then
!           write(6,998) nlistcstmnball2
!        end if

        if (maxcstmnblist*maxnatom / nproc < nlistcstmnball_lcl2) then
           write(6,*) 'Error: nlistcstmnball2(local) exceed number', &
                &     ' of custom NB list at process ',irank
           stop
        end if

        nlistcstmnball3 = icut3
!!! non-existent last loop index cstmnb1
        cstmnb_indexall3(icut2+1) = nlistcstmnball3 + 1

#if defined(MPI)
        call mpi_reduce(nlistcstmnball_lcl3,nlistcstmnball3,1,MPI_INTEGER, &
             &          MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

        if (irank == 0) then
           write(6,998) nlistcstmnball2,nlistcstmnball3
           if (maxcstmnblist*maxnatom/2 < nlistcstmnball3) then
              write(6,*) 'Error: nlistcstmnball3 exceed number', &
                   &     ' of custom NB list'
              stop
           end if
        end if

998     format(/3x,'Custom nonbonded pairs updated (i-j,i-j-k): ', 2i18)


!    --- cell index and bookkeeping
     ELSE

!---    Calculation cell length and number of cells

!       - cell count

        box_min(1:3) = 1.0d5
        do cstmnb1 = 1, ncstmnb
           i = ncstmnblist(cstmnb1)
           itype = atmindex_ncstmnb(i)
           box_min(1:3) = min(box_min(1:3),box(1:3,itype))
        end do
        cellcount(1:3) = INT(box_min(1:3) / rcut_bookSW)

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

        cellindex_last = cellcount(1) * cellcount(2) * cellcount(3) - 1
!        cellindex_last = cellcount(1) - 1 &
!             &         + cellcount(1) &
!             &         * ((cellcount(2) - 1) &
!             &          + cellcount(2) * (cellcount(3) - 1))

!       - cell size
        cellsize(1:3) = box_min(1:3) / dble(cellcount(1:3))

!----   Initialize cell-index variables  ----
        cell_natom(0:cellindex_last) = 0

!----   Cell-index registration

        do cstmnb1 = 1, ncstmnb
           i = ncstmnblist(cstmnb1)

           atom_to_cell(1:3,i) = mod(int((atmcor_sig(1:3,i) + box_min(1:3)) &
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

!      ---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM DISTANCE ---

!     - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!        DO i = 1, natom

        icut2 = 0
        icut3 = 0
        nlistcstmnball_lcl2 = 0
        nlistcstmnball_lcl3 = 0

        looplast = ncstmnb - 1

        DO cstmnb1 = loopinit, looplast, loopstep ! loop over ncstmnblist(cstmnb1)
           i = ncstmnblist(cstmnb1)
           itype_inttyp = atmindex(i)
           itype = atmindex_ncstmnb(i)

           itmp2 = 0   ! temporal number of neighboring atoms of atom(i)
           itmp3(1:maxcstmnblist2) = 0  ! temporal number of neighboring atoms of atom(j)
           iitmp3(1:maxcstmnblist2) = 0  ! temporal number of neighboring atoms of atom(j)
           j_old = 0

!           DO cstmnb2 = cstmnb1+1, ncstmnb     ! loop over ncstmnblist(cstmnb2)
!              j = ncstmnblist(cstmnb2)
!              jtype_inttyp = atmindex(j)
!              jtype = atmindex_ncstmnb(j)

           cell_exist1(1:3) = atom_to_cell(1:3,i)

           do z1d = cell_exist1(3)-2, cell_exist1(3)+2
              z1 = mod(z1d + cellcount(3), cellcount(3))
              do y1d = cell_exist1(2)-2, cell_exist1(2)+2
                 y1 = mod(y1d + cellcount(2), cellcount(2))
                 posyz1 = y1 * cellcount(1) + z1 * cellcount(1)*cellcount(2)
                 do x1d = cell_exist1(1)-2, cell_exist1(1)+2
                    x1 = mod(x1d + cellcount(1), cellcount(1))
                    cellindex1 = x1 + posyz1

                    do jj = 1, cell_natom(cellindex1)

                       j = cell_atomid(jj,cellindex1) ! atom(j)
                       if (i >= j) cycle
                       jtype_inttyp= atmindex(j)
                       jtype = atmindex_ncstmnb(j)

!----                  if not custom NB interaction type, skip the loop

!                       if (inter_inttyp(itype_inttyp,jtype_inttyp) /= INTTYPE_CSTMNB) cycle

!                      - CHECK IF ATOM(i) is close to ATOM(j) -

                       rij(1:3) = atmcor_sig(1:3,i) - atmcor_sig(1:3,j)

!                      - periodic boundary -

                       i2btype = cstmnbtypeindex2(itype,jtype)

                       rij(1:3) = rij(1:3) - box(1:3,i2btype) &
                            &   * anint(rij(1:3)*box_inv(1:3,i2btype))

!                      * quite large box has been set in case of
!                        NON-periodicity to avoid imaging

                       drij = rij(1)**2 + rij(2)**2 + rij(3)**2

                       if (drij < sqcut) then

                          itmp2 = itmp2 + 1
                          list_tmp_i2(itmp2) = j

                       end if

                       cell_exist2(1:3) = atom_to_cell(1:3,j)

                       do z2d = cell_exist2(3)-2, cell_exist2(3)+2
                          z2 = mod(z2d + cellcount(3), cellcount(3))
                          do y2d = cell_exist2(2)-2, cell_exist2(2)+2
                             y2 = mod(y2d + cellcount(2), cellcount(2))
                             posyz2 = y2 * cellcount(1) &
                                  & + z2 * cellcount(1)*cellcount(2)
                             do x2d = cell_exist2(1)-2, cell_exist2(1)+2
                                x2 = mod(x2d + cellcount(1), cellcount(1))
                                cellindex2 = x2 + posyz2

                                do kk = 1, cell_natom(cellindex2)

                                   k = cell_atomid(kk,cellindex2) ! atom(k)
                                   if (j >= k) cycle
                                   ktype_inttyp = atmindex(k)
                                   ktype = atmindex_ncstmnb(k)

!----                              if not custom NB interaction type, skip the loop

!                                   if (inter_inttyp(jtype_inttyp,ktype_inttyp) /= INTTYPE_CSTMNB) cycle
!                                   if (inter_inttyp(itype_inttyp,ktype_inttyp) /= INTTYPE_CSTMNB) cycle

!                                  - CHECK IF two of ATOM(i), ATOM(j) and ATOM(k) are close -

                                   i3btype = cstmnbtypeindex3(itype,jtype,ktype)

                                   rik(1:3) = atmcor_sig(1:3,i) - atmcor_sig(1:3,k)
                                   rik(1:3) = rik(1:3) - box(1:3,i3btype) &
                                        &              * anint(rik(1:3)*box_inv(1:3,i3btype))
                                   drik = rik(1)**2 + rik(2)**2 + rik(3)**2

                                   rjk(1:3) = atmcor_sig(1:3,j) - atmcor_sig(1:3,k)
                                   rjk(1:3) = rjk(1:3) - box(1:3,i3btype) &
                                        &              * anint(rjk(1:3)*box_inv(1:3,i3btype))
                                   drjk = rjk(1)**2 + rjk(2)**2 + rjk(3)**2

                                   if (drij < sqcut .and. (drik < sqcut .or. drjk < sqcut)) then
                                      itmp3(itmp2) = itmp3(itmp2) + 1
                                      list_tmp_i3(itmp2,itmp3(itmp2)) = k
                                      cycle
                                   end if

                                   if (drik < sqcut .and. drjk < sqcut) then

                                       ! not to count up 'itmp2' if only this condition matches more than twice for the same j
                                      if (j/=j_old) itmp2 = itmp2 + 1
                                      j_old = j

                                      list_tmp_i2(itmp2) = j
                                      itmp3(itmp2) = itmp3(itmp2) + 1
                                      list_tmp_i3(itmp2,itmp3(itmp2)) = k
                                      cycle
                                   end if
                                end do

                             end do
                          end do
                       end do

                       if (itmp2>0) iitmp3(itmp2) = itmp3(itmp2)           ! = itmp3 - excluded atoms

!                    END DO             ! end of loop over ncstmnblist(cstmnb2)
                    end do

                 end do
              end do
           end do

!           if (itmp2>0) iitmp3(itmp2) = itmp3(itmp2)           ! = itmp3 - excluded atoms


           iitmp2 = itmp2           ! = itmp2 - excluded atoms

!          - transfer data of list_tmp_i to cstmnb_list

           cstmnb_indexall2(cstmnb1) = icut2 + 1
           IF (iitmp2 == 0) then

!             - atom(i) has no neighbor but still 0 is put to cstmnb_list -

              cstmnb_indexall3(icut2+1) = icut3 + 1

              icut2 = icut2 + 1
              cstmnb_listall2(icut2) = 0
              nlistcstmnball_lcl2 = nlistcstmnball_lcl2 + 1

              icut3 = icut3 + 1
              cstmnb_listall3(icut3) = 0
              nlistcstmnball_lcl3 = nlistcstmnball_lcl3 + 1

           ELSE

!             - atom(i) has at least one neighbor -

              do j = 1, itmp2
                 if (list_tmp_i2(j) > 0) then
                    cstmnb_listall2(icut2+j)  = list_tmp_i2(j)
                 end if

                 cstmnb_indexall3(icut2+j) = icut3 + 1

                 IF (iitmp3(j) == 0) then

                    icut3 = icut3 + 1

                    cstmnb_listall3(icut3) = 0
                    nlistcstmnball_lcl3 = nlistcstmnball_lcl3 + 1

                 ELSE

                    do k = 1, iitmp3(j)
                       if (list_tmp_i3(j,k) > 0) then
                          cstmnb_listall3(icut3+k) = list_tmp_i3(j,k)
                       end if
                    end do
                    icut3 = icut3 + iitmp3(j)
                    nlistcstmnball_lcl3 = nlistcstmnball_lcl3 + iitmp3(j)

                 END IF

              end do

              icut2 = icut2 + iitmp2
              nlistcstmnball_lcl2 = nlistcstmnball_lcl2 + iitmp2
           END IF

        END DO                    ! end of loop over ncstmnblist(cstmnb1)


!       - end of loop -

        nlistcstmnball2 = icut2
!!! non-existent last loop index cstmnb1
        cstmnb_indexall2(cstmnb1) = nlistcstmnball2 + 1

#if defined(MPI)
        call mpi_reduce(nlistcstmnball_lcl2,nlistcstmnball2,1,MPI_INTEGER, &
             &          MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

!        if (irank == 0) then
!           write(6,998) nlistcstmnball2
!        end if

        if (maxcstmnblist*maxnatom / nproc < nlistcstmnball_lcl2) then
           write(6,*) 'Error: nlistcstmnball2(local) exceed number', &
                &     ' of custom NB list at process ',irank
           stop
        end if

        nlistcstmnball3 = icut3
!!! non-existent last loop index cstmnb1
        cstmnb_indexall3(icut2+1) = nlistcstmnball3 + 1

#if defined(MPI)
        call mpi_reduce(nlistcstmnball_lcl3,nlistcstmnball3,1,MPI_INTEGER, &
             &          MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

        if (irank == 0) then
           write(6,999) nlistcstmnball2,nlistcstmnball3
           if (maxcstmnblist*maxnatom/2 < nlistcstmnball3) then
              write(6,*) 'Error: nlistcstmnball3 exceed number', &
                   &     ' of custom NB list'
              stop
           end if
        end if
999     format(/3x,'Custom nonbonded pairs updated (i-j,i-j-k): ', 2i18)

     END IF

  END IF

!---- dynamic memory release ----
  deallocate(list_tmp_i3)

!      +      +      +      +      +      +

end subroutine mklist2a_cstmnb
