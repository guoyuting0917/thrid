!******************************
!*  contemp_region.f Ver.1.6  *
!*      for peachgk_md.f      *
!*            by G.Kikugawa   *
!******************************
! Time-stamp: <2015-02-17 19:02:58 gota>

subroutine contemp_region(npoly,npolytyp, &
     &                    npoly_mole,npoly_atom, &
     &                    nwater, &
     &                    nmatom,nmatyp,nmatomtyp, &
     &                    degfree_poly,degfree_water, &
     &                    degfree_ma, &
     &                    xcel,ycel,zcel, &
     &                    iftcratom, &
     &                    ntcregion, &
     &                    tcxpos1,tcxpos2, &
     &                    tcypos1,tcypos2, &
     &                    tczpos1,tczpos2, &
     &                    r_tcont, &
     &                    ifoutthc, &
     &                    det_ene_kin)

  use md_global

  implicit none

!     subroutine to control temperature by WOODCOCK's method
!       in specific regions described in "tempcont.ini"
!
!     temperature T = 2<K>/(Df.R)
!     where <K> is the total kinetic energy, 
!           Df  is the degree of freedom
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

  integer,intent(in):: degfree_poly(:) ! degree of freedom of poly
  integer,intent(in):: degfree_water   ! degree of freedom of H2O
  integer,intent(in):: degfree_ma(:)   ! degree of freedom of monatomic mole.

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  logical,intent(in):: iftcratom    ! region temp. control based on atom or mole.
  integer,intent(in):: ntcregion       ! number of region to control temp.
  real(8),intent(in):: tcxpos1(:),tcxpos2(:)
                                ! x-position of temp. control region
  real(8),intent(in):: tcypos1(:),tcypos2(:)
                                ! y-position of temp. control region
  real(8),intent(in):: tczpos1(:),tczpos2(:)
                                ! z-position of temp. control region
  real(8),intent(in):: r_tcont(:)       ! control temp. in each region

  logical,intent(in):: ifoutthc     ! flag for outputting thermal control file

!     OUTPUT
  real(8),intent(out):: det_ene_kin(:)          ! for outputting thermal control data

! LOCAL:

  real(8):: factor           ! factor to multipy to Ekin 
  real(8):: veloc2           ! |veloc|**2
  real(8):: real_temp        ! real temperature
  real(8):: sfact            ! factor to multiply velocity

  integer:: i,j             ! do loop index
  integer:: j1,j2           ! do loop index
  integer:: ipolytyp
  integer:: ipoly
  integer:: iwater
  integer:: imatyp
  integer:: atm_index
  integer:: atm_index_tmp

  real(8):: molemass
  real(8):: molegrav(3)

  integer:: itcr

  integer:: natom_tc(ntcregion) ! number of atoms in region
  integer:: degfree_tc(ntcregion) ! degree of freedom in region
  integer:: atmindex_tc(maxnatom,ntcregion) ! index of atoms in region
  real(8):: ene_kin(ntcregion) ! kinetic energy in region

  real(8):: ene_kin_new(ntcregion) ! kinetic energy in region (updated)

  real(8):: atmvel_tmp(3,natom) ! temporary atom velocity without streaming vel.

!     +     +     +     +     +     +     +

!---- some preparation

  atmvel_tmp(1:3,1:natom) = atmvel(1:3,1:natom)
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) - atmvel_strm(1:3,1:natom)

  do itcr = 1, ntcregion
     natom_tc(itcr) = 0
     degfree_tc(itcr) = 0
     ene_kin(itcr) = 0.0d0
  end do

  atm_index = 0

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
           if (atmcor(1,atm_index) < 0.0d0) then
              molegrav(1) = atmcor(1,atm_index) + xcel
           else if (atmcor(1,atm_index) >= xcel) then
              molegrav(1) = atmcor(1,atm_index) - xcel
           end if

           if (atmcor(2,atm_index) < 0.0d0) then
              molegrav(2) = atmcor(2,atm_index) + ycel
           else if (atmcor(2,atm_index) >= ycel) then
              molegrav(2) = atmcor(2,atm_index) - ycel
           end if

           if (atmcor(3,atm_index) < 0.0d0) then
              molegrav(3) = atmcor(3,atm_index) + zcel
           else if (atmcor(3,atm_index) >= zcel) then
              molegrav(3) = atmcor(3,atm_index) - zcel
           end if

!              - if this molecule within temp. control region
           do itcr = 1, ntcregion

              if (((tcxpos1(itcr) <= molegrav(1)) .and. &
                   & (molegrav(1) < tcxpos2(itcr))) .or. &
                   & ((tcxpos1(itcr) <= molegrav(1) - xcel) .and. &
                   & (molegrav(1) - xcel < tcxpos2(itcr)))) then ! xcel
              if (((tcypos1(itcr) <= molegrav(2)) .and. &
                   & (molegrav(2) < tcypos2(itcr))) .or. &
                   & ((tcypos1(itcr) <= molegrav(2) - ycel) .and. &
                   & (molegrav(2) - ycel < tcypos2(itcr)))) then ! ycel
              if (((tczpos1(itcr) <= molegrav(3)) .and. &
                   & (molegrav(3) < tczpos2(itcr))) .or. &
                   & ((tczpos1(itcr) <= molegrav(3) - zcel) .and. &
                   & (molegrav(3) - zcel < tczpos2(itcr)))) then ! zcel

                 degfree_tc(itcr) = degfree_tc(itcr) + 3   ! atom-based control

                 natom_tc(itcr) = natom_tc(itcr) + 1
                 atmindex_tc(natom_tc(itcr),itcr) = atm_index

                 veloc2 = atmvel(1,atm_index)**2 &
                      & + atmvel(2,atm_index)**2 &
                      & + atmvel(3,atm_index)**2
                 ene_kin(itcr) = ene_kin(itcr) &
                      &        + atmmass(atm_index) * veloc2

                 exit

              end if
              end if
              end if

           end do

        end do

     end do

!---- molecule-based control
  ELSE

     do ipolytyp = 1, npolytyp

        do i = 1, npoly_mole(ipolytyp)
           molemass = 0.0d0
           molegrav(1:3) = 0.0d0

!              - calculate center of mass
           do j = 1, npoly_atom(ipolytyp)
              atm_index = atm_index + 1
              molemass = molemass + atmmass(atm_index)
              molegrav(1:3) = molegrav(1:3) &
                   &        + atmmass(atm_index)*atmcor(1:3,atm_index)
           end do

           molegrav(1:3) = molegrav(1:3) / molemass

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

!              - if this molecule within temp. control region
           do itcr = 1, ntcregion

              if (((tcxpos1(itcr) <= molegrav(1)) .and. &
                   & (molegrav(1) < tcxpos2(itcr))) .or. &
                   & ((tcxpos1(itcr) <= molegrav(1) - xcel) .and. &
                   & (molegrav(1) - xcel < tcxpos2(itcr)))) then ! xcel
              if (((tcypos1(itcr) <= molegrav(2)) .and. &
                   & (molegrav(2) < tcypos2(itcr))) .or. &
                   & ((tcypos1(itcr) <= molegrav(2) - ycel) .and. &
                   & (molegrav(2) - ycel < tcypos2(itcr)))) then ! ycel
              if (((tczpos1(itcr) <= molegrav(3)) .and. &
                   & (molegrav(3) < tczpos2(itcr))) .or. &
                   & ((tczpos1(itcr) <= molegrav(3) - zcel) .and. &
                   & (molegrav(3) - zcel < tczpos2(itcr)))) then ! zcel

                 degfree_tc(itcr) = degfree_tc(itcr) &
                      &           + degfree_poly(ipolytyp) &
                      &           / npoly_mole(ipolytyp) ! per molecule

                 atm_index_tmp = atm_index - npoly_atom(ipolytyp)
                 do j = 1, npoly_atom(ipolytyp)
                    atm_index_tmp = atm_index_tmp + 1
                    natom_tc(itcr) = natom_tc(itcr) + 1
                    atmindex_tc(natom_tc(itcr),itcr) = atm_index_tmp

                    veloc2 = atmvel(1,atm_index_tmp)**2 &
                         & + atmvel(2,atm_index_tmp)**2 &
                         & + atmvel(3,atm_index_tmp)**2
                    ene_kin(itcr) = ene_kin(itcr) &
                         &        + atmmass(atm_index_tmp) * veloc2
                 end do

                 exit

              end if
              end if
              end if

           end do

        end do

     end do

  END IF

!     - water type
!      do iwater = 1, nwater
  do iwater = npoly+1, npoly+nwater
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

!        - if this molecule within temp. control region
     do itcr = 1, ntcregion

        if (((tcxpos1(itcr) <= molegrav(1)) .and. &
             & (molegrav(1) < tcxpos2(itcr))) .or. &
             & ((tcxpos1(itcr) <= molegrav(1) - xcel) .and. &
             & (molegrav(1) - xcel < tcxpos2(itcr)))) then ! xcel
        if (((tcypos1(itcr) <= molegrav(2)) .and. &
             & (molegrav(2) < tcypos2(itcr))) .or. &
             & ((tcypos1(itcr) <= molegrav(2) - ycel) .and. &
             & (molegrav(2) - ycel < tcypos2(itcr)))) then ! ycel
        if (((tczpos1(itcr) <= molegrav(3)) .and. &
             & (molegrav(3) < tczpos2(itcr))) .or. &
             & ((tczpos1(itcr) <= molegrav(3) - zcel) .and. &
             & (molegrav(3) - zcel < tczpos2(itcr)))) then ! zcel

           degfree_tc(itcr) = degfree_tc(itcr) &
                &           + degfree_water / nwater ! per molecule

           do j = j1, j2
              atm_index = molept_list(j)

              natom_tc(itcr) = natom_tc(itcr) + 1
              atmindex_tc(natom_tc(itcr),itcr) = atm_index
              veloc2 = atmvel(1,atm_index)**2 &
                   & + atmvel(2,atm_index)**2 &
                   & + atmvel(3,atm_index)**2
              ene_kin(itcr) = ene_kin(itcr) &
                   &        + atmmass(atm_index) * veloc2
           end do

           exit

        end if
        end if
        end if

     end do

  end do

!     - matom type
!      do imatyp = 1, nmatyp
  do imatyp = npoly+nwater+1, npoly+nwater+nmatom
     j1 = molept_index(imatyp)
     j2 = molept_index(imatyp+1) - 1

     do j = j1, j2
        atm_index = molept_list(j)

        molegrav(1:3) = atmcor(1:3,atm_index)

!           - P.B.C.
        if (atmcor(1,atm_index) < 0.0d0) then
           molegrav(1) = atmcor(1,atm_index) + xcel
        else if (atmcor(1,atm_index) >= xcel) then
           molegrav(1) = atmcor(1,atm_index) - xcel
        end if

        if (atmcor(2,atm_index) < 0.0d0) then
           molegrav(2) = atmcor(2,atm_index) + ycel
        else if (atmcor(2,atm_index) >= ycel) then
           molegrav(2) = atmcor(2,atm_index) - ycel
        end if

        if (atmcor(3,atm_index) < 0.0d0) then
           molegrav(3) = atmcor(3,atm_index) + zcel
        else if (atmcor(3,atm_index) >= zcel) then
           molegrav(3) = atmcor(3,atm_index) - zcel
        end if

!           - if this molecule within temp. control region
        do itcr = 1, ntcregion

           if (((tcxpos1(itcr) <= molegrav(1)) .and. &
                & (molegrav(1) < tcxpos2(itcr))) .or. &
                & ((tcxpos1(itcr) <= molegrav(1) - xcel) .and. &
                & (molegrav(1) - xcel < tcxpos2(itcr)))) then ! xcel
           if (((tcypos1(itcr) <= molegrav(2)) .and. &
                & (molegrav(2) < tcypos2(itcr))) .or. &
                & ((tcypos1(itcr) <= molegrav(2) - ycel) .and. &
                & (molegrav(2) - ycel < tcypos2(itcr)))) then ! ycel
           if (((tczpos1(itcr) <= molegrav(3)) .and. &
                & (molegrav(3) < tczpos2(itcr))) .or. &
                & ((tczpos1(itcr) <= molegrav(3) - zcel) .and. &
                & (molegrav(3) - zcel < tczpos2(itcr)))) then ! zcel

              natom_tc(itcr) = natom_tc(itcr) + 1
              atmindex_tc(natom_tc(itcr),itcr) = atm_index
              degfree_tc(itcr) = degfree_tc(itcr) + 3   ! matom type

              veloc2 = atmvel(1,atm_index)**2 &
                   & + atmvel(2,atm_index)**2 &
                   & + atmvel(3,atm_index)**2
              ene_kin(itcr) = ene_kin(itcr) &
                   &        + atmmass(atm_index) * veloc2

              exit

           end if
           end if
           end if

        end do

     end do

  end do

!---- Calculate temperature and reset velocity
!     - calculate temperature
  do itcr = 1, ntcregion

     if (degfree_tc(itcr) /= 0) then
        factor = 1.0d0 / dble(degfree_tc(itcr)) ! ene_kin = 2.0 * K
     else
        factor = 0.0d0
     end if
     real_temp = factor*ene_kin(itcr)

!        - reset velocity
     if (real_temp > 1.0d-5) then
        sfact = sqrt(r_tcont(itcr) / real_temp)
     else
        sfact = 1.0d0
     end if

     do i = 1, natom_tc(itcr)
        atm_index = atmindex_tc(i,itcr)
        atmvel(1:3,atm_index) = 1.0d0 * (sfact-1.0d0)  &
             &                * atmvel(1:3,atm_index) &
             &                + atmvel(1:3,atm_index)
     end do

  end do

!---- Calculate imposed (or extracted) kinetic energy 
  if (ifoutthc) then
     ene_kin_new(1:ntcregion) = 0.0d0

     do itcr = 1, ntcregion
        do i = 1, natom_tc(itcr)
           atm_index = atmindex_tc(i,itcr)
           ene_kin_new(itcr) = ene_kin_new(itcr) &
                &        + atmmass(atm_index) &
                &        * (atmvel(1,atm_index)**2 &
                &         + atmvel(2,atm_index)**2 &
                &         + atmvel(3,atm_index)**2)
        end do

        det_ene_kin(itcr) = det_ene_kin(itcr) &
             &            + 0.5d0 * (ene_kin_new(itcr) - ene_kin(itcr))
     end do

  end if

!--- reset velocity
  atmvel(1:3,1:natom) = atmvel(1:3,1:natom) + atmvel_strm(1:3,1:natom)

!     +     +     +     +     +     +     +

end subroutine contemp_region
