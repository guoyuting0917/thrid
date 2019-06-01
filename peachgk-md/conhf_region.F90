!******************************
!*  conhf_region.f90 Ver.1.6  *
!*      for peachgk_md.f      *
!*            by G.Kikugawa   *
!******************************
! Time-stamp: <2015-02-17 19:24:44 gota>

subroutine conhf_region(npoly,npolytyp, &
     &                  npoly_mole,npoly_atom, &
     &                  nwater, &
     &                  nmatom,nmatyp,nmatomtyp, &
     &                  xcel,ycel,zcel, &
     &                  iftcratom, &
     &                  nhfcregion, &
     &                  hfczpos1,hfczpos2, &
     &                  r_hfcont)

  use md_global

  implicit none

!     subroutine to control heat flux by means proposed by Jund and Jullien
!       (Jund, P. and Jullien, R., Phys. Rev. B, 59 (1999), pp. 13707-13711.)
!       in specific regions (z-slab only) described in "hfcont.ini"
!
!--------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!! Caution !!!!!!!!!!!!!!!!!!!!!!!!!
!  In this version, do not use atom-based control
!    if some molecules have contraint bonding.
!    (excluding water molecules)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ARGUMENTS:
!     INPUT
  integer,intent(in):: npoly           ! number of polymer1
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatom          ! number of monatomic molecules
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: iftcratom    ! region temp. control based on atom or mole.

  integer,intent(in):: nhfcregion      ! number of region to control heat flux
  real(8),intent(in):: hfczpos1(:),hfczpos2(:)
                                ! z-position of heat flux control region
  real(8),intent(in):: r_hfcont(:)      ! magnitude of heat flux in each region
                                ! (converted to the input energy)

! LOCAL:       
  real(8):: veloc2           ! |veloc|**2                 
  real(8):: sfact            ! factor to multiply velocity   

  integer:: i,j             ! do loop index
  integer:: j1,j2           ! do loop index
  integer:: ipoly
  integer:: iwater
  integer:: imatyp
  integer:: atm_index

  real(8):: molemass
  real(8):: molegrav(3)

  integer:: ihfcr

  integer:: natom_hfc(nhfcregion) ! number of atoms in region
  integer:: atmindex_hfc(maxnatom,nhfcregion) ! index of atoms in region
  real(8):: mass_region(nhfcregion) ! total mass in region
  real(8):: vel_com(3,nhfcregion) ! velocity of COM in region
  real(8):: ene_kin(nhfcregion) ! kinetic energy in region
  real(8):: com_kin(nhfcregion) ! COM kinetic energy in region

  real(8):: atmvel_tmp(3,natom) ! temporary atom velocity without streaming vel.

!     +     +     +     +     +     +     +

!---- some preparation

  atmvel_tmp(1:3,1:natom) = atmvel(1:3,1:natom)
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) - atmvel_strm(1:3,1:natom)

  natom_hfc(1:nhfcregion) = 0
  mass_region(1:nhfcregion) = 0.0d0
  vel_com(1:3,1:nhfcregion) = 0.0d0
  ene_kin(1:nhfcregion) = 0.0d0

  atm_index = 1

!-------- region-based temperature control --------

!     - poly type

!---- atom-based control
  IF (iftcratom) then

     DO ipoly = 1, npoly
        j1 = molept_index(ipoly)
        j2 = molept_index(ipoly+1) - 1
        do j = j1, j2
           atm_index = molept_list(j)

           molegrav(1:3) = atmcor(1:3,atm_index)

!              - P.B.C.
           if (molegrav(1) < 0.0d0) then
              molegrav(1) = molegrav(1) + xcel
           else if (molegrav(1) >= xcel) then
              molegrav(1) = molegrav(1) - xcel
           end if

           if (molegrav(2) < 0.0d0) then
              molegrav(2) = molegrav(2) + ycel
           else if (molegrav(2) >= ycel) then
              molegrav(2) = molegrav(2) - ycel
           end if

           if (molegrav(3) < 0.0d0) then
              molegrav(3) = molegrav(3) + zcel
           else if (molegrav(3) >= zcel) then
              molegrav(3) = molegrav(3) - zcel
           end if

!              - if this molecule within h.f. control region
           do ihfcr = 1, nhfcregion

              if (((hfczpos1(ihfcr) <= molegrav(3)) .and. &
                   & (molegrav(3) < hfczpos2(ihfcr))) .or. &
                   & ((hfczpos1(ihfcr) <= molegrav(3) - zcel) .and. &
                   & (molegrav(3) - zcel < hfczpos2(ihfcr)))) then

                 natom_hfc(ihfcr) = natom_hfc(ihfcr) + 1
                 atmindex_hfc(natom_hfc(ihfcr),ihfcr) = atm_index

                 mass_region(ihfcr) = mass_region(ihfcr) &
                      &             + atmmass(atm_index)
                 vel_com(1:3,ihfcr) = vel_com(1:3,ihfcr) &
                      &             + atmmass(atm_index) &
                      &             * atmvel(1:3,atm_index)
                 veloc2 = atmvel(1,atm_index)**2 &
                      & + atmvel(2,atm_index)**2 &
                      & + atmvel(3,atm_index)**2
                 ene_kin(ihfcr) = ene_kin(ihfcr) &
                      &         + atmmass(atm_index) * veloc2

                 exit
              end if

           end do

        end do

     END DO

!---- molecule-based control
  ELSE

     DO ipoly = 1, npoly
        j1 = molept_index(ipoly)
        j2 = molept_index(ipoly+1) - 1

        molemass = 0.0d0
        molegrav(1) = 0.0d0
        molegrav(2) = 0.0d0
        molegrav(3) = 0.0d0

!           - calculate center of mass
        do j = j1, j2
           atm_index = molept_list(j)

           molemass = molemass + atmmass(atm_index)
           molegrav(1:3) = molegrav(1:3) &
                &        + atmmass(atm_index)*atmcor(1:3,atm_index)
        end do

        molegrav(1:3) = molegrav(1:3) / molemass

!           - P.B.C.
        if (molegrav(1) < 0.0d0) then
           molegrav(1) = molegrav(1) + xcel
        else if (molegrav(1) >= xcel) then
           molegrav(1) = molegrav(1) - xcel
        end if

        if (molegrav(2) < 0.0d0) then
           molegrav(2) = molegrav(2) + ycel
        else if (molegrav(2) >= ycel) then
           molegrav(2) = molegrav(2) - ycel
        end if

        if (molegrav(3) < 0.0d0) then
           molegrav(3) = molegrav(3) + zcel
        else if (molegrav(3) >= zcel) then
           molegrav(3) = molegrav(3) - zcel
        end if

!           - if this molecule within h.f. control region
        do ihfcr = 1, nhfcregion

           if (((hfczpos1(ihfcr) <= molegrav(3)) .and. &
                & (molegrav(3) < hfczpos2(ihfcr))) .or. &
                & ((hfczpos1(ihfcr) <= molegrav(3) - zcel) .and. &
                & (molegrav(3) - zcel < hfczpos2(ihfcr)))) then

              do j = j1, j2
                 atm_index = molept_list(j)

                 natom_hfc(ihfcr) = natom_hfc(ihfcr) + 1
                 atmindex_hfc(natom_hfc(ihfcr),ihfcr) = atm_index

                 mass_region(ihfcr) = mass_region(ihfcr) &
                      &             + atmmass(atm_index)
                 vel_com(1:3,ihfcr) = vel_com(1:3,ihfcr) &
                      &             + atmmass(atm_index) &
                      &             * atmvel(1:3,atm_index)
                 veloc2 = atmvel(1,atm_index)**2 &
                      & + atmvel(2,atm_index)**2 &
                      & + atmvel(3,atm_index)**2
                 ene_kin(ihfcr) = ene_kin(ihfcr) &
                      &         + atmmass(atm_index) * veloc2

              end do

              exit
           end if

        end do

     END DO

  END IF

!     - water type
!      do iwater = 1, nwater
  DO iwater = npoly+1, npoly+nwater
     j1 = molept_index(iwater)
     j2 = molept_index(iwater+1) - 1

     molemass = 0.0d0
     molegrav(1:3) = 0.0d0

     do j = j1, j2
        atm_index = molept_list(j)

!           - calculate center of mass
        molemass = molemass + atmmass(atm_index)
        molegrav(1:3) = molegrav(1:3) &
             &        + atmmass(atm_index)*atmcor(1:3,atm_index)
     end do

     molegrav(1:3) = molegrav(1:3) / molemass

!        - P.B.C.
     if (molegrav(1) < 0.0d0) then
        molegrav(1) = molegrav(1) + xcel
     else if (molegrav(1) >= xcel) then
        molegrav(1) = molegrav(1) - xcel
     end if

     if (molegrav(2) < 0.0d0) then
        molegrav(2) = molegrav(2) + ycel
     else if (molegrav(2) >= ycel) then
        molegrav(2) = molegrav(2) - ycel
     end if

     if (molegrav(3) < 0.0d0) then
        molegrav(3) = molegrav(3) + zcel
     else if (molegrav(3) >= zcel) then
        molegrav(3) = molegrav(3) - zcel
     end if

!        - if this molecule within h.f. control region
     do ihfcr = 1, nhfcregion

        if (((hfczpos1(ihfcr) <= molegrav(3)) .and. &
             & (molegrav(3) < hfczpos2(ihfcr))) .or. &
             & ((hfczpos1(ihfcr) <= molegrav(3) - zcel) .and. &
             & (molegrav(3) - zcel < hfczpos2(ihfcr)))) then

           do j = j1, j2
              atm_index = molept_list(j)

              natom_hfc(ihfcr) = natom_hfc(ihfcr) + 1
              atmindex_hfc(natom_hfc(ihfcr),ihfcr) = atm_index

              mass_region(ihfcr) = mass_region(ihfcr) &
                   &             + atmmass(atm_index)
              vel_com(1:3,ihfcr) = vel_com(1:3,ihfcr) &
                   &             + atmmass(atm_index) &
                   &             * atmvel(1:3,atm_index)
              veloc2 = atmvel(1,atm_index)**2 &
                   & + atmvel(2,atm_index)**2 &
                   & + atmvel(3,atm_index)**2
              ene_kin(ihfcr) = ene_kin(ihfcr) &
                   &         + atmmass(atm_index) * veloc2

           end do

           exit
        end if

     end do

  END DO

!     - matom type
!      do imatyp = 1, nmatyp
  DO imatyp = npoly+nwater+1, npoly+nwater+nmatom
     j1 = molept_index(imatyp)
     j2 = molept_index(imatyp+1) - 1

     do j = j1, j2
        atm_index = molept_list(j)

        molegrav(1:3) = atmcor(1:3,atm_index)

!           - P.B.C.
        if (molegrav(1) < 0.0d0) then
           molegrav(1) = molegrav(1) + xcel
        else if (molegrav(1) >= xcel) then
           molegrav(1) = molegrav(1) - xcel
        end if

        if (molegrav(2) < 0.0d0) then
           molegrav(2) = molegrav(2) + ycel
        else if (molegrav(2) >= ycel) then
           molegrav(2) = molegrav(2) - ycel
        end if

        if (molegrav(3) < 0.0d0) then
           molegrav(3) = molegrav(3) + zcel
        else if (molegrav(3) >= zcel) then
           molegrav(3) = molegrav(3) - zcel
        end if

!           - if this molecule within temp. control region
        do ihfcr = 1, nhfcregion

           if (((hfczpos1(ihfcr) <= molegrav(3)) .and. &
                & (molegrav(3) < hfczpos2(ihfcr))) .or. &
                & ((hfczpos1(ihfcr) <= molegrav(3) - zcel) .and. &
                & (molegrav(3) - zcel < hfczpos2(ihfcr)))) then

              natom_hfc(ihfcr) = natom_hfc(ihfcr) + 1
              atmindex_hfc(natom_hfc(ihfcr),ihfcr) = atm_index

              mass_region(ihfcr) = mass_region(ihfcr) &
                   &             + atmmass(atm_index)
              vel_com(1:3,ihfcr) = vel_com(1:3,ihfcr) &
                   &             + atmmass(atm_index) &
                   &             * atmvel(1:3,atm_index)
              veloc2 = atmvel(1,atm_index)**2 &
                   & + atmvel(2,atm_index)**2 &
                   & + atmvel(3,atm_index)**2
              ene_kin(ihfcr) = ene_kin(ihfcr) &
                   &         + atmmass(atm_index) * veloc2

              exit
           end if

        end do

     end do

  END DO

#if defined(_HFCONT_DEBUG)
  write(6,*) '***** h.f. control debug info *****'
#endif

!---- Calculate kinetic energy and reset velocity
  do ihfcr = 1, nhfcregion

!        - calculate velocity of COM in eachregion

#if defined(_HFCONT_NOMOM)
     ! no momentum conservation
     vel_com(1:3,ihfcr) = 0.0d0
     com_kin(ihfcr) = 0.0d0
#else
     if (natom_hfc(ihfcr) /= 0) then
        vel_com(1:3,ihfcr) = vel_com(1:3,ihfcr) / mass_region(ihfcr)
        com_kin(ihfcr) = mass_region(ihfcr) &
             &         * (vel_com(1,ihfcr)*vel_com(1,ihfcr) &
             &          + vel_com(2,ihfcr)*vel_com(2,ihfcr) &
             &          + vel_com(3,ihfcr)*vel_com(3,ihfcr))
     else
        vel_com(1:3,ihfcr) = 0.0d0
        com_kin(ihfcr) = 0.0d0
     end if
#endif

!        - reset velocity
     if (natom_hfc(ihfcr) /= 0) then

!           - sqrt(1 + input energy / relative kinetic energy)
        sfact = sqrt(1.0d0 &
             &     + r_hfcont(ihfcr) &
             &     / (0.5d0 * (ene_kin(ihfcr) - com_kin(ihfcr))))
                                                   ! ene_kin = 2.0 * K
#if defined(_HFCONT_DEBUG)
        write(6,*) 'region No. ',ihfcr, &
             &     ', curr energy: ',0.5d0*ene_kin(ihfcr), &
             &     ', scaling fact: ',sfact
#endif

     else

        sfact = 0.0d0
#if defined(_HFCONT_DEBUG)
        write(6,*) 'region No. ',ihfcr, &
             &     ', curr energy: ',0.5d0*ene_kin(ihfcr), &
             &     ', scaling fact: ',sfact
#endif

     end if

     do i = 1, natom_hfc(ihfcr)
        atm_index = atmindex_hfc(i,ihfcr)

!           - vi(new) = vg + sfact*(vi - vg)
        atmvel(1:3,atm_index) = vel_com(1:3,ihfcr) &
             &                + sfact * (atmvel(1:3,atm_index) &
             &                - vel_com(1:3,ihfcr))

     end do

  end do

!--- reset velocity
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) + atmvel_strm(1:3,1:natom)

!     +     +     +     +     +     +     +

end subroutine conhf_region
