!***************************************
!*  mklist_cnpvw.f90 Ver.1.3           *
!*      for peachgk_md.f               *
!*            by N.Yamamoto            *
!*            modified by G.Kikugawa   *
!***************************************
! Time-stamp: <2015-01-22 00:53:28 gota>

subroutine mklist2a_rpvw( rcut_bookrpvw, zcel, ifbookrpvw)

  use md_global
  use mpi_global

  implicit none

!    subroutine to generate ATOM-BASED RP-VW interaction ATOM LIST
!
! ARGUMENT:
!     INPUT
  real(8),intent(in):: rcut_bookrpvw     ! cut off radius of bookkeeping of RP-VW

!  real(8),intent(in):: xcel             ! x cell length[non-d]
!  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: ifbookrpvw      ! flag for bookkeeping of RP-VW interaction

!        ATTENTION: list_excl and index_ecxl should have been 
!                   reconstructed so that only i < j are contained.
        
!               i.e.  list for atom(i) begins at nb_indexall(i) of nb_listall
!        If atom(i) does not have any neighbor, then 0 is input in nb_listall.
!
! MPI version by G.Kikugawa

! LOCAL:

!---- for general use --

  real(8):: sqcut            ! rcut2**2
  real(8):: zij_2
                                ! cooridnate of cent of res i & j 
                                ! difference, etc.

  integer:: list_tmp_i(maxnatom)
                                !  temporal storage for list for atom(i)

  integer:: icut            ! current size of lists  

!  real*8:: box(3)           ! BOX size
!  real*8:: box_inv(3)       ! inverse of BOX size

!     -- misc --

  integer:: i, j            ! do loop index for atoms
  integer:: itmp, jtmp,iitmp
  integer:: i_rpvw, j_rpvw
  integer:: itype, jtype

!      logical:: ifmorsepair
!      logical:: ifmorse_former, ifmorse_latter

! MPI
  integer:: nlistrpvwall_lcl ! local variable nlistrpvwall

!     +     +     +     +     +     +     +

!  box(1) = xcel
!  box(2) = ycel
!  box(3) = zcel
!  box_inv(1) = 1.0d0/xcel
!  box_inv(2) = 1.0d0/ycel
!  box_inv(3) = 1.0d0/zcel

  sqcut = rcut_bookrpvw * rcut_bookrpvw

  if (.not. ifbookrpvw) then          ! if bookkeeping is not used
     sqcut = (zcel/2.0d0)*(zcel/2.0d0)
  end if

!---- Initialize nb_list or other ----

  rpvw_indexall(1:maxnatom+maxnproc) = 0

  rpvw_listall(1:maxnatom*maxrpvwlist / nproc) = 0

!---- MAKE ATOM-ATOM NONBONDED LIST BASED ON ATOM & GP DISTANCE ---

!     - MAKE NONBONDED LIST FOR ATOM(i) -

! MPI  each process calculates only a part of contribution.
!     tasks are distributed in a cyclic fashion to each process.
!
!      DO i = 1, natom

  icut = 0
  nlistrpvwall_lcl = 0

  looplast = nrpvw - 1

  DO i_rpvw = loopinit, looplast, loopstep ! loop over nrpvwselist(i_rpvw)
     i = nrpvwlist(i_rpvw)
     itype = atmindex(i)

     itmp = 0               ! temporal number of neighboring atoms
                                ! of atom(i)
            
     DO j_rpvw = i_rpvw+1, nrpvw     ! loop over nrpvwselist(j_rpvw)
        j = nrpvwlist(j_rpvw)
        jtype = atmindex(j)

!----       if not RP-VW interaction type, skip the loop

        if (inter_inttyp(itype,jtype) /= INTTYPE_RPVW) cycle

!       - CHECK IF ATOM(i) is close to ATOM(j) -

        zij_2 = ( atmcor(3,i) - atmcor(3,j) ) **2

!           - periodic boundary -
            
!        rij(1:3) = rij(1:3) - box(1:3) * &
!             &                dnint(rij(1:3)*box_inv(1:3))

!                 * quite large box has been set in case of
!                   NON-periodicity to avoid imaging

!        dr = rij(1)**2 + rij(2)**2 + rij(3)**2

        if (zij_2 < sqcut) then
 
           itmp = itmp + 1
           list_tmp_i(itmp) = j

        end if

     END DO                 ! end of loop over nrpvwselist(j_rpvw)
                  

     iitmp = itmp           ! = itmp - excluded atoms

!        - transfer data of list_tmp_i to rpvw_list

     rpvw_indexall(i_rpvw) = icut + 1

     IF (iitmp == 0) then

!        - atom(i) has no neighbor but still 0 is put to rpvw_list -

        rpvw_listall(icut+1) = 0
        icut = icut + 1
        nlistrpvwall_lcl = nlistrpvwall_lcl + 1

     ELSE

!           - atom(i) has at least one neighbor -

        jtmp = 0
        do j = 1, itmp
           if (list_tmp_i(j) > 0) then
              jtmp = jtmp + 1
              rpvw_listall(icut+jtmp) = list_tmp_i(j)
           end if
        end do
        icut = icut + iitmp
        nlistrpvwall_lcl = nlistrpvwall_lcl + iitmp

     END IF

  END DO                    ! end of loop over nrpvwselist(i_rpvw)


!     - end of loop -

  nlistrpvwall = icut
!!! non-existence last loop index i_rpvw
  rpvw_indexall(i_rpvw) = nlistrpvwall + 1

#if defined(MPI)
  call mpi_reduce(nlistrpvwall_lcl,nlistrpvwall,1,MPI_INTEGER, &
       &                MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif

  if (irank == 0) then
     write(6,999) nlistrpvwall
  end if

  if (maxrpvwlist*maxnatom / nproc < nlistrpvwall_lcl) then
     write(6,*) 'Error: nlistrpvwall(local) exceed number', &
          &     ' of RP-VW list at process ',irank
     stop
  end if

999 format(/3x,'RP-VW interaction pairs updated: ', i18) 

!      +      +      +      +      +      +

end subroutine mklist2a_rpvw
