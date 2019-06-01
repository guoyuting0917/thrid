!*******************************
!*  rdpara_cstmnb.f90 Ver.1.2  *
!*      for peachgk_md.f       *
!*            by G.Kikugawa    *
!*******************************
! Time-stamp: <2015-01-28 14:15:37 gota>

!!! Virtual wall (continuous media) interaction by LJ potential
!!!   originally coded by Jo Suzuki, IFS, Tohoku University

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

  integer:: i,j           ! do loop index
  integer:: nword

  integer:: nline

  real(8):: hcelmin

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

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
              read(fredat(3),*) zcut_WALL
              read(fredat(4),*) zcut_bookWALL
              read(fredat(5),*) nstep_bookWALL
           endif
           read(fredat(6),*) E_WALL(ncstmnbtyp)
           read(fredat(7),*) S_WALL(ncstmnbtyp)
           if (nline == 1) then
              read(fredat(8),*) R_WALL
              read(fredat(9),*) dstnc
           endif

        end do

     end if

  END DO READCSTMNB

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

! set cstmnb type indexes on all atoms
  DO i = 1, natom
     do j = 1, ncstmnbtyp
        if(atmtyp(i)(1:2) == para_cstmnbtyp(j)) then
           index_cstmnbtyp(i) = j
        end if
     end do
  END DO
    
!---- non-dimensionalize
  zcut_WALL = zcut_WALL / xref
  zcut_bookWALL = zcut_bookWALL / xref
  E_WALL(1:ncstmnbtyp) = E_WALL(1:ncstmnbtyp) / eref
  S_WALL(1:ncstmnbtyp) = S_WALL(1:ncstmnbtyp) / xref
  R_WALL = R_WALL * xref**3
  dstnc = dstnc / xref

!---- checking cutoff parameters

  hcelmin = zcel*0.5d0

  if (abs(zcut_WALL) < 1.0d-16) then
     zcut_WALL = hcelmin
  end if

!---- register intermolecular interaction type for particle pair
!** usually not needed to change here

!     - detect custom NB pair

  close(iuparacstmnb)

!     +     +     +     +     +     +     +

end subroutine rdpara_cstmnb
