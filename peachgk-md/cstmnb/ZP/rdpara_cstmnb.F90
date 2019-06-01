!*******************************
!*  rdpara_cstmnb.f90 Ver.1.2  *
!*      for peachgk_md.f       *
!*            by G.Kikugawa    *
!*******************************
! Time-stamp: <2015-01-22 00:32:28 gota>

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

  real(8):: para_alphaZP_aniso_tmp(maxncstmnbtyp)
  real(8):: para_cZP_iso_tmp(maxncstmnbtyp)
  real(8):: para_welZP_tmp(maxncstmnbtyp)
  real(8):: para_radZP_tmp(maxncstmnbtyp)

  real(8):: an=6.0221367d23  ! Avogadro's number
  integer:: i,j,k           ! do loop index
  integer:: nword

  logical:: ifcstmnb_M

  logical:: ifcstmnb_former
  logical:: ifcstmnb_latter

  integer:: cellcount(3)

  real(8):: hcelmin

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- dynamic memory allocation for neighbor lists ----

  allocate(cstmnb_listall(maxnatom*maxcstmnblist / nproc))
  allocate(ZPcond1_listall(maxnatom*maxcstmnblist / nproc))
  allocate(ZPcond2_listall(maxnatom*maxcstmnblist / nproc))

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
        do
           call rdfree_w(iuparacstmnb, maxnword, fredat, nword)
           if (fredat(1) == ' ') cycle READCSTMNB
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           ncstmnbtyp = ncstmnbtyp + 1

!!! Describe your column definition here
!!!  (variables should be declared in cstmnb.F90)
           para_cstmnbtyp(ncstmnbtyp) = fredat(2)(1:5)
           read(fredat(3),*) rcut_ZPcond(ncstmnbtyp)
           read(fredat(4),*) rcut_bookZPcond(ncstmnbtyp)
           read(fredat(5),*) nstep_bookZPcond(ncstmnbtyp)
           read(fredat(6),*) para_ZPint_p1(ncstmnbtyp)
           read(fredat(7),*) para_ZPint_p2(ncstmnbtyp)
           read(fredat(8),*) rcut_ZPaniso(ncstmnbtyp)
           read(fredat(9),*) rcut_bookZPaniso(ncstmnbtyp)
           read(fredat(10),*) nstep_bookZPaniso(ncstmnbtyp)
           read(fredat(11),*) para_alphaZP_aniso_tmp(ncstmnbtyp)
           read(fredat(12),*) rcut_ZPiso(ncstmnbtyp)
           read(fredat(13),*) rcut_bookZPiso(ncstmnbtyp)
           read(fredat(14),*) nstep_bookZPiso(ncstmnbtyp)
           read(fredat(15),*) para_cZP_iso_tmp(ncstmnbtyp)
           read(fredat(16),*) para_welZP_tmp(ncstmnbtyp)
           read(fredat(17),*) para_radZP_tmp(ncstmnbtyp)
           read(fredat(18),*) ifrenew_msurf(ncstmnbtyp)

        end do

     end if

  END DO READCSTMNB

!---- non-dimensionalize
  rcut_ZPcond(1:ncstmnbtyp) = rcut_ZPcond(1:ncstmnbtyp) / xref
  rcut_bookZPcond(1:ncstmnbtyp) = rcut_bookZPcond(1:ncstmnbtyp) / xref
  para_ZPint_p1(1:ncstmnbtyp) = para_ZPint_p1(1:ncstmnbtyp) / xref
  para_ZPint_p2(1:ncstmnbtyp) = para_ZPint_p2(1:ncstmnbtyp) / xref

  rcut_ZPaniso(1:ncstmnbtyp) = rcut_ZPaniso(1:ncstmnbtyp) / xref
  rcut_bookZPaniso(1:ncstmnbtyp) = rcut_bookZPaniso(1:ncstmnbtyp) / xref

  rcut_ZPiso(1:ncstmnbtyp) = rcut_ZPiso(1:ncstmnbtyp) / xref
  rcut_bookZPiso(1:ncstmnbtyp) = rcut_bookZPiso(1:ncstmnbtyp) / xref
  para_welZP_tmp(1:ncstmnbtyp) = para_welZP_tmp(1:ncstmnbtyp) / eref
  para_radZP_tmp(1:ncstmnbtyp) = para_radZP_tmp(1:ncstmnbtyp) / xref

!---- cell index checking
  if (.not.ifbookcstmnb .and. ifcellindex_cstmnb) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       enable ifbookcstmnb'
     stop
  end if

  ifcellindex_ZPaniso = ifcellindex_cstmnb
!  ifcellindex_ZPiso   = ifcellindex_cstmnb

  if (ifcellindex_cstmnb) then
     !!! The cutoff parameter in the first line is used for all species.

!    -- for ZP anisotropic term
     cellcount(1) = INT(xcel/rcut_bookZPaniso(1))
     cellcount(2) = INT(ycel/rcut_bookZPaniso(1))
     cellcount(3) = INT(zcel/rcut_bookZPaniso(1))

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (ZP aniso)'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_ZPaniso = .false.
        ifcellindex_cstmnb = ifcellindex_ZPaniso
     end if

#if defined(_DO_NOT_USE_THIS)
!    -- for ZP isotropic term
     cellcount(1) = INT(xcel/rcut_bookZPiso(1))
     cellcount(2) = INT(ycel/rcut_bookZPiso(1))
     cellcount(3) = INT(zcel/rcut_bookZPiso(1))

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (ZP iso)'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_ZPiso = .false.
     end if

#endif
  end if

!---- checking cutoff parameters
  !!! The cutoff parameter in the first line is used for all species.

  hcelmin = min(xcel*0.5d0,ycel*0.5d0,zcel*0.5d0)

! -- for ZP cond
  if (abs(rcut_ZPcond(1)) < 1.0d-16) then
     rcut_ZPcond(1) = zcel
  end if

!  if (rcut_ZPcond(1) > hcelmin) then
!     write(6,*) 'Error: find incorrect rcut_ZPcond= ',rcut_ZPcond(1)
!     stop
!  end if

!  if ((ifbookcstmnb) .and. (rcut_bookZPcond(1) > hcelmin)) then
!     write(6,*) 'Error: find incorrect rcut_bookZPcond= ', &
!          &     rcut_bookZPcond(1)
!     stop
!  end if

! -- for ZP anisotropic term
  if (abs(rcut_ZPaniso(1)) < 1.0d-16) then
     rcut_ZPaniso(1) = hcelmin
  end if

  if (rcut_ZPaniso(1) > hcelmin) then
     write(6,*) 'Error: find incorrect rcut_ZPaniso= ',rcut_ZPaniso(1)
     stop
  end if

  if ((ifbookcstmnb) .and. (rcut_bookZPaniso(1) > hcelmin)) then
     write(6,*) 'Error: find incorrect rcut_bookZPaniso= ', &
          &     rcut_bookZPaniso(1)
     stop
  end if

#if defined(_DO_NOT_USE_THIS)
! -- for ZP isotropic term
  if (abs(rcut_ZPiso(1)) < 1.0d-16) then
     rcut_ZPiso(1) = hcelmin
  end if

  if (rcut_ZPiso(1) > hcelmin) then
     write(6,*) 'Error: find incorrect rcut_ZPiso= ',rcut_ZPiso(1)
     stop
  end if

  if ((ifbookcstmnb) .and. (rcut_bookZPiso(1) > hcelmin)) then
     write(6,*) 'Error: find incorrect rcut_bookZPiso= ', &
          &     rcut_bookZPiso(1)
     stop
  end if
#endif

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

!---- list of atoms for ZP cond interaction
  nZPcond = 0

  do i = 1, natom
     ifcstmnb_M = .false.
     do j = 1, ncstmnbtyp

        if (atmtyp(i)(1:2) == para_cstmnbtyp(j)(1:2)) then  ! this means metal
           ifcstmnb_M = .true.
           exit
        end if

     end do

     if (ifcstmnb_M .or. abs(atmchrg(i)) < 1.0d-16) cycle
                                ! this atom is metal or does not have a charge

     nZPcond = nZPcond + 1
     nZPcondlist(nZPcond) = i

  end do

!---- order and link custom NB parameter
  do i = 1, natmtyp
     do j = 1, natmtyp

        do k = 1, ncstmnbtyp
           ifcstmnb_former = .false.
           ifcstmnb_latter = .false.

           !-- scan i
           if (para_atmtyp(i)(1:2) == para_cstmnbtyp(k)(1:2)) then
              ifcstmnb_former = .true.
           else if (para_atmtyp(i)(1:2) == para_cstmnbtyp(k)(4:5)) then
              ifcstmnb_latter = .true.
           end if

           !-- scan j
           if (ifcstmnb_former) then

              if (para_atmtyp(j)(1:2) == para_cstmnbtyp(k)(4:5)) then

                 para_alphaZP_aniso(i,j) = para_alphaZP_aniso_tmp(k)
                 para_cZP_iso(i,j) = para_cZP_iso_tmp(k)
                 para_welZP(i,j) = para_welZP_tmp(k)
                 para_radZP(i,j) = para_radZP_tmp(k)
                 exit

              end if

           else if (ifcstmnb_latter) then

              if (para_atmtyp(j)(1:2) == para_cstmnbtyp(k)(1:2)) then

                 para_alphaZP_aniso(i,j) = para_alphaZP_aniso_tmp(k)
                 para_cZP_iso(i,j) = para_cZP_iso_tmp(k)
                 para_welZP(i,j) = para_welZP_tmp(k)
                 para_radZP(i,j) = para_radZP_tmp(k)
                 exit

              end if

           end if

           para_alphaZP_aniso(i,j) = 1.0d0
           para_cZP_iso(i,j) = 0.0d0
           para_welZP(i,j) = 0.0d0
           para_radZP(i,j) = 0.0d0

        end do

     end do
  end do

!---- register intermolecular interaction type for particle pair
!** usually not needed to change here

!     - detect custom NB pair
  do i = 1, natmtyp
     do j = i, natmtyp

        do k = 1, ncstmnbtyp
           ifcstmnb_former = .false.
           ifcstmnb_latter = .false.

           !-- scan i
           if (para_atmtyp(i)(1:2) == para_cstmnbtyp(k)(1:2)) then
              ifcstmnb_former = .true.
           else if (para_atmtyp(i)(1:2) == para_cstmnbtyp(k)(4:5)) then
              ifcstmnb_latter = .true.
           else
              cycle
           end if

           !-- scan j
           if (ifcstmnb_former) then

              if (para_atmtyp(j)(1:2) == para_cstmnbtyp(k)(4:5)) then
                 inter_inttyp(i,j) = INTTYPE_CSTMNB
                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           else if (ifcstmnb_latter) then

              if (para_atmtyp(j)(1:2) == para_cstmnbtyp(k)(1:2)) then

                 inter_inttyp(i,j) = INTTYPE_CSTMNB
                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           end if

        end do

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
