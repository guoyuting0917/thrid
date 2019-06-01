!*******************************
!*  rdpara_cstmnb.f90 Ver.1.3  *
!*      for peachgk_md.f       *
!*            by G.Kikugawa    *
!*******************************
! Time-stamp: <2015-01-22 01:11:31 gota>

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

  logical:: ifcstmnb_i
  logical:: ifcstmnb_j
  integer:: cstmnb_i
  integer:: cstmnb_j

  integer:: nline

  integer:: cellcount(3)

  real(8):: hcelmin

  real(8):: A_MATSUI(maxncstmnbtyp)   ! A parameter for Matsui potential
  real(8):: B_MATSUI(maxncstmnbtyp)   ! B parameter for Matsui potential
  real(8):: C_MATSUI(maxncstmnbtyp)   ! CiCi parameter for Matsui potential

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- dynamic memory allocation for neighbor lists ----

  allocate(cstmnb_listall(maxnatom*maxcstmnblist / nproc))

  write(6,*) 'Dynamically allocate memory for cstmnb lists at process: ',irank

!-------- Read parameter of custom NB interaction --------

  ncstmnbtyp = 0
  READCSTMNB:DO
     call rdfree(iuparacstmnb, maxnword, fredat)

     if (fredat(1) == '<END>') then
        exit
     else if (fredat(1) == ' ') then
        cycle READCSTMNB

     else if (fredat(1) == '<CUSTOM_NB>') then
        nline = 0
        do
           call rdfree_w(iuparacstmnb, maxnword, fredat, nword)
           if (fredat(1) == ' ') cycle READCSTMNB
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

!!! Describe your column definition here
!!!  (variables should be declared in cstmnb.F90)
           nline = nline + 1

           ncstmnbtyp = ncstmnbtyp + 1
           para_cstmnbtyp(ncstmnbtyp) = fredat(2)(1:2)
           if (nline == 1) then 
              read(fredat(3),*) rcut_MATSUI
              read(fredat(4),*) rcut_bookMATSUI
              read(fredat(5),*) nstep_bookMATSUI
           endif
           read(fredat(6),*) A_MATSUI(ncstmnbtyp)
           read(fredat(7),*) B_MATSUI(ncstmnbtyp)
           read(fredat(8),*) C_MATSUI(ncstmnbtyp)
           if (nline == 1) then
              read(fredat(9),*) D_MATSUI
           endif

        end do

     end if

  END DO READCSTMNB

!---- non-dimensionalize
  rcut_MATSUI = rcut_MATSUI / xref
  rcut_bookMATSUI = rcut_bookMATSUI / xref
  A_MATSUI(1:ncstmnbtyp) = A_MATSUI(1:ncstmnbtyp) / xref
  B_MATSUI(1:ncstmnbtyp) = B_MATSUI(1:ncstmnbtyp) / xref
  C_MATSUI(1:ncstmnbtyp) = C_MATSUI(1:ncstmnbtyp) / (xref**6 * eref)
  D_MATSUI = D_MATSUI / fref

!---- cell index checking
  if (.not.ifbookcstmnb .and. ifcellindex_cstmnb) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       enable ifbookcstmnb'
     stop
  end if

  if (ifcellindex_cstmnb) then

     cellcount(1) = INT(xcel/rcut_bookMATSUI)
     cellcount(2) = INT(ycel/rcut_bookMATSUI)
     cellcount(3) = INT(zcel/rcut_bookMATSUI)

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (Matsui)'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_cstmnb = .false.
     end if

  endif

!---- checking cutoff parameters

  hcelmin = min(xcel*0.5d0,ycel*0.5d0,zcel*0.5d0)

  if (abs(rcut_MATSUI) < 1.0d-16) then
     rcut_MATSUI = hcelmin
  end if

  if (rcut_MATSUI > hcelmin) then
     write(6,*) 'Error: find incorrect rcut_MATSUI = ',rcut_MATSUI
     stop
  end if

  if ((ifbookcstmnb) .and. (rcut_bookMATSUI > hcelmin)) then
     write(6,*) 'Error: find incorrect rcut_bookMATSUI= ', &
          &     rcut_bookMATSUI
     stop
  end if

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
        end if

     end do
  end do

#if defined(_DO_NOT_USE_THIS)
!---- order and link custom NB parameter
  ifcalMATSUI(1:natmtyp,1:natmtyp) = .false.

  do i = 1, natmtyp
     do j = 1, natmtyp

        ifcstmnb_i = .false.
        ifcstmnb_j = .false.

        do k = 1, ncstmnbtyp

           !-- scan i
           if (para_atmtyp(i)(1:2) == para_cstmnbtyp(k)(1:2)) then
              ifcstmnb_i = .true.
           endif

           !-- scan j
           if (para_atmtyp(j)(1:2) == para_cstmnbtyp(k)(1:2)) then
              ifcstmnb_j = .true.
           endif

        end do

        if (ifcstmnb_i .and. ifcstmnb_j) then
           ifcalMATSUI(i,j) = .true.   ! this pair calculated in cstmnb
        endif

     end do
  end do
#endif

!---- register intermolecular interaction type for particle pair
!** usually not needed to change here

!     - detect custom NB pair
  do i = 1, natmtyp
     do j = i, natmtyp

        ifcstmnb_i = .false.
        ifcstmnb_j = .false.

        do k = 1, ncstmnbtyp
           !-- scan i
           if (para_atmtyp(i)(1:2) == para_cstmnbtyp(k)(1:2)) then
              ifcstmnb_i = .true.
              cstmnb_i = k
           end if

           !-- scan j
           if (para_atmtyp(j)(1:2) == para_cstmnbtyp(k)(1:2)) then
              ifcstmnb_j = .true.
              cstmnb_j = k
           end if

        end do

        if (ifcstmnb_i .and. ifcstmnb_j) then
           inter_inttyp(i,j) = INTTYPE_CSTMNB
           inter_inttyp(j,i) = inter_inttyp(i,j)

           vdw_welij(i,j) = 0.0d0
           vdw_welij(j,i) = 0.0d0
           vdw_radij(i,j) = 0.0d0
           vdw_radij(j,i) = 0.0d0


           !-- Exception: omit short range interaction between AL-HI and HI-HI
           if ((para_atmtyp(i)(1:2) == 'AL') .and. (para_atmtyp(j)(1:2) == 'HI') &
        & .or. (para_atmtyp(i)(1:2) == 'HI') .and. (para_atmtyp(j)(1:2) == 'AL') &
        & .or. (para_atmtyp(i)(1:2) == 'HI') .and. (para_atmtyp(j)(1:2) == 'HI')) then
              AA_MATSUI(i,j) = 0.0d0
              AA_MATSUI(j,i) = AA_MATSUI(i,j)
              BB_MATSUI(i,j) = 0.0d0
              BB_MATSUI(j,i) = BB_MATSUI(i,j)
              CC_MATSUI(i,j) = 0.0d0
              CC_MATSUI(j,i) = CC_MATSUI(i,j)
           else
              AA_MATSUI(i,j) = A_MATSUI(cstmnb_i) + A_MATSUI(cstmnb_j)
              AA_MATSUI(j,i) = AA_MATSUI(i,j)
              BB_MATSUI(i,j) = B_MATSUI(cstmnb_i) + B_MATSUI(cstmnb_j)
              BB_MATSUI(j,i) = BB_MATSUI(i,j)
              CC_MATSUI(i,j) = sqrt(C_MATSUI(cstmnb_i) * C_MATSUI(cstmnb_j))
              CC_MATSUI(j,i) = CC_MATSUI(i,j)
           end if

        end if

     end do
  end do

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
