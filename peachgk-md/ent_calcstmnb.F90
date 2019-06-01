!*******************************
!*  ent_calcstmnb.f90 Ver.1.0  *
!*     for peachgk_md.f        *
!*            by G.Kikugawa    *
!*******************************
! Time-stamp: <2015-07-17 19:02:45 gota>

subroutine ent_calcstmnb(mts_flag,mts_cstmnb, &
     &                   xcel,ycel,zcel, &
     &                   npoly,npolytyp,npoly_mole,npoly_atom, &
     &                   nwater,nmatom,nmatyp,nmatomtyp, &
     &                   force,pot_cstmnb, &
     &                   for_viri_cstmnb,pot_viri_cstmnb, &
     &                   pot_virit_cstmnb, &
     &                   ncstmnbex, &
     &                   pot_cstmnbex, &
     &                   for_viri_cstmnbex,pot_viri_cstmnbex, &
     &                   pot_virit_cstmnbex, &
     &                   ifcalpreatom,ifcalpremole, &
     &                   current_step, &
     &                   nstep_short,istep_short)

  use md_global

  use cstmnb, only: calcstmnb, calcstmnbp

  implicit none

! ARGUMENTS:
!     INPUT
  integer,intent(in):: mts_flag        ! flag for MTS integration
                                       ! 1 long-range force mode
                                       ! 2 medium-range force mode
                                       ! 3 short-range force mode
  integer,intent(in):: mts_cstmnb      ! MTS flag for custom NB interaction

  integer,intent(in):: npoly           ! all number of poly
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly
  
  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatom          ! number of monatomic molecules
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  real(8),intent(in):: xcel            ! x cell length[non-d]
  real(8),intent(in):: ycel            ! y cell length[non-d]
  real(8),intent(in):: zcel            ! z cell length[non-d]

  logical,intent(in):: ifcalpremole    ! pressure calculation of molecule
  logical,intent(in):: ifcalpreatom    ! pressure calculation of atom

  integer,intent(in):: current_step    ! current time step
  integer,intent(in):: nstep_short     ! number of step for short force
  integer,intent(in):: istep_short

  integer,intent(in):: ncstmnbex       ! number of extra custom NB output

!     OUTPUT
  real(8),intent(inout):: force(:,:)   ! force calculated here

  real(8),intent(inout):: pot_cstmnb   ! Custom NB potential

  real(8),intent(inout):: for_viri_cstmnb(:,:)
                                         ! virial(custom NB force) of each atom
  real(8),intent(inout):: pot_viri_cstmnb
                                   ! virial(custom NB potential) of each atom
  real(8),intent(inout):: pot_virit_cstmnb(:,:) ! virial tensor (custom NB)

  real(8),intent(inout):: pot_cstmnbex(:) ! Custom NB extra potential

  real(8),intent(inout):: for_viri_cstmnbex(:,:,:)
                                   ! extra virial(custom NB force) of each atom
  real(8),intent(inout):: pot_viri_cstmnbex(:)
                               ! extra virial(custom NB potential) of each atom
  real(8),intent(inout):: pot_virit_cstmnbex(:,:,:) 
                                          ! extra virial tensor (custom NB)

! LOCAL:

!     +     +     +     +     +     +     +

!---- CSTMNB function Ver. 1
#if defined(_CSTMNB_V1)
  if (mts_cstmnb == mts_flag) THEN

     if (ifcalpremole .or. ifcalpreatom) then

        call calcstmnbp(xcel,ycel,zcel, &
             &          npoly,npolytyp,npoly_mole,npoly_atom, &
             &          nwater,nmatom,nmatyp,nmatomtyp, &
             &          force,pot_cstmnb, &
             &          for_viri_cstmnb,pot_viri_cstmnb, &
             &          pot_virit_cstmnb)

     else

        call calcstmnb(xcel,ycel,zcel, &
             &         npoly,npolytyp,npoly_mole,npoly_atom, &
             &         nwater,nmatom,nmatyp,nmatomtyp, &
             &         force,pot_cstmnb)

     end if

  end if

!---- CSTMNB function Ver. 2
#elif defined(_CSTMNB_V2)
  if (ifcalpremole .or. ifcalpreatom) then

     call calcstmnbp(mts_flag,mts_cstmnb, &
          &          xcel,ycel,zcel, &
          &          npoly,npolytyp,npoly_mole,npoly_atom, &
          &          nwater,nmatom,nmatyp,nmatomtyp, &
          &          force,pot_cstmnb, &
          &          for_viri_cstmnb,pot_viri_cstmnb, &
          &          pot_virit_cstmnb, &
          &          ncstmnbex, &
          &          pot_cstmnbex, &
          &          for_viri_cstmnbex,pot_viri_cstmnbex, &
          &          pot_virit_cstmnbex, &
          &          current_step, &
          &          nstep_short,istep_short)

  else

     call calcstmnb(mts_flag,mts_cstmnb, &
          &         xcel,ycel,zcel, &
          &         npoly,npolytyp,npoly_mole,npoly_atom, &
          &         nwater,nmatom,nmatyp,nmatomtyp, &
          &         force,pot_cstmnb, &
          &         ncstmnbex, &
          &         pot_cstmnbex, &
          &         current_step, &
          &         nstep_short,istep_short)

  end if

#endif

!     +     +     +     +     +     +     +     +

end subroutine ent_calcstmnb
