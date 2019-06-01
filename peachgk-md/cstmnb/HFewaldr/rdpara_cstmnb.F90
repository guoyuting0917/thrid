!*******************************
!*  rdpara_cstmnb.f90 Ver.1.3  *
!*      for peachgk_md.f       *
!*            by G.Kikugawa    *
!*******************************
! Time-stamp: <2015-01-22 01:04:34 gota>

subroutine rdpara_cstmnb(iuparacstmnb, &
     &                   ifcellindex_cstmnb,ifbookcstmnb, &
     &                   xref,eref,mref,qref, &
     &                   vref,timeref,tempref,pref,fref,eps0ref, &
     &                   xcel,ycel,zcel)

  use interface_tools

  use md_global
  use mpi_global
  use cstmnb

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: iuparacstmnb    ! input custom NB parameter file unit

  logical,intent(inout):: ifcellindex_cstmnb ! flag for cell index (custom NB)
  logical,intent(in):: ifbookcstmnb       
                               ! flag for bookkeeping of custom NB interaction

  real(8),intent(in):: xref            ! distanse base value [m]
  real(8),intent(in):: eref            ! energy base value [J]
  real(8),intent(in):: mref            ! mass base value [kg]
  real(8),intent(in):: qref            ! charge base value [C]
  real(8),intent(in):: vref            ! velocity base value [m/s]
  real(8),intent(in):: timeref         ! time base value [sec]
  real(8),intent(in):: tempref         ! temperature base value [K]
  real(8),intent(in):: pref            ! pressure base value [Pa]
  real(8),intent(in):: fref            ! force base value [N]

  real(8),intent(in):: eps0ref         ! dielectric constant base value [c^2/Jm]

  real(8),intent(in):: xcel            ! x cell length [non-d]
  real(8),intent(in):: ycel            ! y cell length [non-d]
  real(8),intent(in):: zcel            ! z cell length [non-d]

! LOCAL:
  character(80):: fredat(maxnword)

  real(8):: an=6.0221367d23  ! Avogadro's number
  integer:: i,j,k           ! do loop index
  integer:: nword

  ! logical:: ifcstmnb_former
  ! logical:: ifcstmnb_latter

  real(8):: hcelmin
  integer:: cellcount(3)

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- dynamic memory allocation for neighbor lists ----

  allocate(cstmnb_listall(maxnatom*maxcstmnblist / nproc))

  write(6,*) 'Dynamically allocate memory for cstmnb lists at process: ',irank

!-------- Read parameter of custom NB interaction --------

  ! ncstmnbtyp = 0
  READCSTMNB:DO
     call rdfree(iuparacstmnb, maxnword, fredat)

     if (fredat(1) == '<END>') then
        exit
     else if (fredat(1) == ' ') then
        cycle READCSTMNB

     else if (fredat(1) == '<CUSTOM_NB>') then

        do
           call rdfree_w(iuparacstmnb, maxnword, fredat, nword)
           if (fredat(1) == ' ') cycle READCSTMNB
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

!!! Describe your column definition here
!!!  (variables should be declared in cstmnb.F90)
           read(fredat(2),*) rrcut_hfewr
           read(fredat(3),*) rrcutbook_hfewr
           read(fredat(4),*) nstep_bookhfewr
           read(fredat(5),*) alpha_hfewr

        end do

     end if

  END DO READCSTMNB

!---- non-dimensionalize
  rrcut_hfewr = rrcut_hfewr / xref
  rrcutbook_hfewr = rrcutbook_hfewr / xref
  alpha_hfewr = alpha_hfewr * xref

!---- checking cutoff parameters

  hcelmin = min(xcel*0.5d0,ycel*0.5d0,zcel*0.5d0)

  if (abs(rrcut_hfewr) < 1.0d-16) then
     rrcut_hfewr = hcelmin
  end if

  if (rrcut_hfewr > hcelmin) then
     write(6,*) 'Error: find incorrect rrcut_hfewr = ',rrcut_hfewr
     stop
  end if

  if ((ifbookcstmnb) .and. (rrcutbook_hfewr > hcelmin)) then
     write(6,*) 'Error: find incorrect rrcutbook_hfewr= ', &
          &     rrcutbook_hfewr
     stop
  end if

!---- check cellindex method
  if (.not.ifbookcstmnb .and. ifcellindex_cstmnb) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       choose ifbookcstmnb'
     stop
  end if

  if (ifcellindex_cstmnb) then
     cellcount(1) = INT(xcel/rrcutbook_hfewr)
     cellcount(2) = INT(ycel/rrcutbook_hfewr)
     cellcount(3) = INT(zcel/rrcutbook_hfewr)

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (custom NB)'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_cstmnb = .false.
     end if
  end if

#if defined(_DO_NOT_USE_THIS)
!---- list of atoms of custom NB interaction
!** usually not needed to change here
  ncstmnb = 0

  do i = 1, natom
     ifcstmnb(i) = .false.
     do j = 1, ncstmnbtyp

        if (atmtyp(i)(1:2) == para_cstmnbtyp(j)(1:2)) then
           ncstmnb = ncstmnb + 1
           ncstmnblist(ncstmnb) = i
           ifcstmnb(i) = .true.
           exit
        else if (atmtyp(i)(1:2) == para_cstmnbtyp(j)(4:5)) then
           ncstmnb = ncstmnb + 1
           ncstmnblist(ncstmnb) = i
           ifcstmnb(i) = .true.
           exit
        end if

     end do
  end do
#endif

#if defined(_DO_NOT_USE_THIS)
!---- register intermolecular interaction type for particle pair
!** usually not needed to change here

!     - detect custom NB pair
  do i = 1, natmtyp
     do j = i, natmtyp

        if (inter_inttyp(i,j) == INTTYPE_VDW) then

           inter_inttyp(i,j) = INTTYPE_CSTMNB
           inter_inttyp(j,i) = inter_inttyp(i,j)

        end if

     end do
  end do
#endif

#if defined(_INTTYP_DEBUG)
  do i = 1, natmtyp
     do j = i, natmtyp
        write(6,*) i,j,inter_inttyp(i,j)
     end do
  end do
#endif

  close(iuparacstmnb)

!     +     +     +     +     +     +     +

end subroutine rdpara_cstmnb
