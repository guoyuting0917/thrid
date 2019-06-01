!*******************************
!*  rdpara_cstmnb.f90 Ver.1.3  *
!*      for peachgk_md.f       *
!*            by G.Kikugawa    *
!*******************************
! Time-stamp: <2015-01-22 01:00:59 gota>

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

  logical:: ifcstmnb_former
  logical:: ifcstmnb_latter

  integer:: nline

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- dynamic memory allocation for neighbor lists ----

! no execution
#if defined(_DO_NOT_USE_THIS)
  allocate(cstmnb_listall(maxnatom*maxcstmnblist / nproc))

  write(6,*) 'Dynamically allocate memory for cstmnb lists at process: ',irank
#endif

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

           if (nline == 1) then   ! read rcut and other parameters
              read(fredat(2),*) rrcut_spLJ
              read(fredat(3),*) rcut_spLJ
              read(fredat(4),*) alpha_spLJ

           else
              ncstmnbtyp = ncstmnbtyp + 1
              para_cstmnbtyp(ncstmnbtyp) = fredat(2)(1:5)

           end if

        end do

     end if

  END DO READCSTMNB

!---- non-dimensionalize
  rcut_spLJ = rcut_spLJ / xref
  rrcut_spLJ = rrcut_spLJ / xref
  alpha_spLJ = alpha_spLJ * xref

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

!---- order and link custom NB parameter
  ifcalspLJ(1:natmtyp,1:natmtyp) = .false.

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

                 ifcalspLJ(i,j) = .true.   ! this pair calculated in cstmnb
                 exit

              end if

           else if (ifcstmnb_latter) then

              if (para_atmtyp(j)(1:2) == para_cstmnbtyp(k)(1:2)) then

                 ifcalspLJ(i,j) = .true.   ! this pair calculated in cstmnb
                 exit

              end if

           end if

        end do

     end do
  end do

#if defined(_DO_NOT_USE_THIS)
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

#endif

  close(iuparacstmnb)

!     +     +     +     +     +     +     +

end subroutine rdpara_cstmnb
