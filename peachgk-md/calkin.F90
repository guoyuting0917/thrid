!*****************************
!*  calkin.f90 Ver.2.1       *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-02-17 18:49:57 gota>

subroutine calkin(npoly,npolytyp,npoly_mole, &
     &            nwater, &
     &            nmatom,nmatyp,nmatomtyp, &
     &            degfree_poly,degfree_water,degfree_ma, &
     &            ene_kin_poly,ene_kin_water,ene_kin_ma, &
     &            temp_poly,temp_water,temp_ma, &
     &            degfree_all,ene_kin_all,temp_all, &
     &            mchain,text_c, &
     &            ene_kin_th,ene_pot_th, &
     &            pext_c, &
     &            ene_kin_ba,ene_pot_ba,extra_pot_ba, &
     &            pcont_axis)

  use md_global

  implicit none

!    subroutine to calculate kinetic energy
!
!    temperature T = 2<K>/(Df.R)
!    where <K> is the total kinetic energy,
!           Df  is the degree of freedom (usually 3N)

! ARGUMENTS:
!     INPUT
  integer,intent(in):: npoly           ! number of polymer1
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatom          ! number of monatomic molecules
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  integer,intent(in):: degfree_poly(:) ! degree of freedom of polymer1
  integer,intent(in):: degfree_water   ! degree of freedom of H2O
  integer,intent(in):: degfree_ma(:)   ! degree of freedom of monatomic mole.
  integer,intent(in):: degfree_all     ! degree of freedom of all atoms

  integer,intent(in):: mchain          ! Nose-Hoover chain number (>1 in NVT)
  real(8),intent(in):: text_c           ! external temp. [K] (Nose-Hoover chain)

  real(8),intent(in):: pext_c           ! external pressure [non-d]

  character(5),intent(in),optional:: pcont_axis  ! axis for pressure control (iso, aniso, etc.)

!     OUTPUT

  real(8),intent(out):: ene_kin_poly(:)  ! kinetic energy of polymer1
  real(8),intent(out):: ene_kin_water    ! kinetic energy of water
  real(8),intent(out):: ene_kin_ma(:)    ! kinetic energy of monatomic mole.
  real(8),intent(out):: ene_kin_all      ! kinetic energy
  real(8),intent(out):: temp_poly(:)     ! temperature of polymer1
  real(8),intent(out):: temp_water       ! temperature of H2O
  real(8),intent(out):: temp_ma(:)       ! temperature of MA
  real(8),intent(out):: temp_all         ! temperature of all atoms

  real(8),intent(out):: ene_kin_th       ! kinetic energy of NHC thermostat
  real(8),intent(out):: ene_pot_th       ! potential energy of NHC thermostat

  real(8),intent(out):: ene_kin_ba       ! kinetic energy of Andersen barostat
  real(8),intent(out):: ene_pot_ba       ! potential energy of Andersen barostat
  real(8),intent(out):: extra_pot_ba     ! extra potentail (=kTxi(1))

! LOCAL:

  real(8):: factor           ! factor to multipy to Ekin
  real(8):: veloc2           ! |veloc|**2
  integer:: i,j              ! do loop index
  integer:: j1,j2
  integer:: jj
  integer:: ipolytyp
  integer:: imatyp
  integer:: imole1

  real(8):: vol              ! volume

  character(5):: pcont_axis_here
  real(8):: dfba             ! DoF of barostat

  real(8):: atmvel_tmp(3,natom) ! temporary atom velocity without streaming vel.

!     +     +     +     +     +     +     +

!---- CALCULATE KINETIC ENERGY ----

! - some preparation -

  atmvel_tmp(1:3,1:natom) = atmvel(1:3,1:natom)
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) - atmvel_strm(1:3,1:natom)

  ene_kin_poly(1:npolytyp)   = 0.0d0
  ene_kin_water  = 0.0d0
  ene_kin_ma(1:nmatyp)     = 0.0d0
  ene_kin_all    = 0.0d0

  ene_kin_th     = 0.0d0
  ene_pot_th     = 0.0d0

  ene_kin_ba     = 0.0d0
  ene_pot_ba     = 0.0d0

! NOT in MTK simulation
  pcont_axis_here = 'iso'
  if (present(pcont_axis)) then   ! if optional argument 'pcont_axis' is present
      pcont_axis_here = pcont_axis
  end if

  if (pcont_axis_here == 'iso') then   ! isotropic NPT
      dfba = 1.0d0
  else if (pcont_axis_here == 'aniso') then   ! full aniisotropic NPT
      dfba = 3.0d0
  else   ! uniaxial NPT (pcont_axis = 'x', 'y', 'z')
      dfba = 1.0d0
  end if

! --- kinetic energy of polymer1 ---

  imole1 = 0

  DO ipolytyp = 1, npolytyp

     DO i = 1, npoly_mole(ipolytyp)
        imole1 = imole1 + 1
        j1 = molept_index(imole1)
        j2 = molept_index(imole1+1) - 1

        do j = j1, j2
           jj = molept_list(j)
           veloc2 = atmvel(1,jj)**2 + atmvel(2,jj)**2 &
                & + atmvel(3,jj)**2
           ene_kin_poly(ipolytyp) = ene_kin_poly(ipolytyp) &
                &                 + 0.5d0 * atmmass(jj) * veloc2
        end do

     END DO

     ene_kin_all = ene_kin_all + ene_kin_poly(ipolytyp)

  END DO

! --- kinetic energy of water ---
  DO i = npoly+1,npoly+nwater
     j1 = molept_index(i)
     j2 = molept_index(i+1) - 1

     do j = j1, j2
        jj = molept_list(j)
        veloc2 = atmvel(1,jj)**2 + atmvel(2,jj)**2 + atmvel(3,jj)**2
        ene_kin_water = ene_kin_water + 0.5d0 * atmmass(jj) * veloc2
     end do

  END DO

  ene_kin_all = ene_kin_all + ene_kin_water

! --- kinetic energy of MA ---
  imole1 = 0

  DO imatyp = 1, nmatyp

     DO i = 1,nmatomtyp(imatyp)
        imole1 = imole1 + 1
        j1 = molept_index(npoly+nwater+imole1)
        j2 = molept_index(npoly+nwater+imole1+1) - 1

        do j = j1, j2
           jj = molept_list(j)
           veloc2 = atmvel(1,jj)**2 + atmvel(2,jj)**2  &
                & + atmvel(3,jj)**2
           ene_kin_ma(imatyp) = ene_kin_ma(imatyp)  &
                &             + 0.5d0 * atmmass(jj) * veloc2
        end do

     END DO

     ene_kin_all = ene_kin_all + ene_kin_ma(imatyp)

  END DO

! --- kinetic energy of thermostat ---
  DO i = 1, mchain
     veloc2 = vlogs(i)*vlogs(i)
     ene_kin_th = ene_kin_th + 0.5d0*qmass(i)*veloc2
  END DO

! --- potential energy of thermostat ---
! first thermostat
  ene_pot_th = ene_pot_th + degfree_all*text_c*xlogs(1)
! 2 to M thermostat
  DO i = 2,mchain
     ene_pot_th = ene_pot_th + text_c*xlogs(i)
  END DO

! --- kinetic energy of barostat ---
  if (pcont_axis_here == 'iso') then   ! isotropic NPT
      veloc2 = vlogv*vlogv
  else   ! anisotropic NPT
      veloc2 = vboxg(1)*vboxg(1) + vboxg(2)*vboxg(2) + vboxg(3)*vboxg(3)
  end if
  ene_kin_ba = ene_kin_ba + 0.5d0*vmass*veloc2

! --- potential energy of barostat ---
  if (pcont_axis_here == 'iso') then   ! isotropic NPT
      vol = exp(3.0d0*xlogv)
      ene_pot_ba = ene_pot_ba + pext_c*vol
  else   ! anisotropic NPT
      vol = xboxh(1)*xboxh(2)*xboxh(3)
      ene_pot_ba = ene_pot_ba + pext_c*vol
  end if
  extra_pot_ba = dfba*text_c*xlogs(1)   != dfba * kB.T * xi(1)

!---- CALCULATE TEMPERATURE ----

!     --- temperature of polymer1 ---
  do ipolytyp = 1, npolytyp
     if (degfree_poly(ipolytyp) /= 0) then
        factor = 2.0d0 / dble(degfree_poly(ipolytyp))
     else
        factor = 0.0d0
     end if
     temp_poly(ipolytyp) = factor*ene_kin_poly(ipolytyp)
  end do

!     --- temperature of H2O ---
  if (degfree_water /= 0) then
     factor = 2.0d0 / dble(degfree_water)
  else
     factor = 0.0d0
  end if
  temp_water = factor*ene_kin_water

!     --- temperature of MA ---
  do imatyp = 1, nmatyp
     if (degfree_ma(imatyp) .ne. 0) then
        factor = 2.0d0 / dble(degfree_ma(imatyp))
     else
        factor = 0.0d0
     end if
     temp_ma(imatyp) = factor*ene_kin_ma(imatyp)
  end do

!     --- temperature of all atoms ---
  if (degfree_all /= 0) then
     factor = 2.0d0 / dble(degfree_all)
  else
     factor = 0.0d0
  end if
  temp_all = factor*ene_kin_all

!--- restore velocity
  atmvel(1:3,1:natom) = atmvel_tmp(1:3,1:natom)

!     +     +     +     +     +     +     +

end subroutine calkin
