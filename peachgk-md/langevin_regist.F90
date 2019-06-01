!*********************************
!*  langevin_regist.f90 Ver.1.0  *
!*      for peachgk_md.f         *
!*            by G.Kikugawa      *
!*********************************
! Time-stamp: <2015-03-01 13:09:38 gota>

subroutine langevin_regist(npoly,npolytyp, &
     &                     npoly_mole,npoly_atom, &
     &                     nwater, &
     &                     nmatom,nmatyp,nmatomtyp, &
     &                     degfree_poly,degfree_water, &
     &                     degfree_ma, &
     &                     xcel,ycel,zcel, &
     &                     iftcratom, &
     &                     nlangeregion, &
     &                     ltxpos1,ltxpos2, &
     &                     ltypos1,ltypos2, &
     &                     ltzpos1,ltzpos2)

  use md_global

  implicit none

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

  integer,intent(in):: nlangeregion    ! number of region for Langevin thermo.
  real(8),intent(in):: ltxpos1(:),ltxpos2(:) 
                                       ! x-position of temp. control region
  real(8),intent(in):: ltypos1(:),ltypos2(:)
                                       ! y-position of temp. control region
  real(8),intent(in):: ltzpos1(:),ltzpos2(:)
                                       ! z-position of temp. control region

! LOCAL:
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

  integer:: iltr

!     +     +     +     +     +     +     +

!---- some preparation

  atm_index = 0

  natom_lt(1:nlangeregion) = 0

!-------- region-based Langevin thermostat --------

!---- registration of atoms and calculate kinetic energy in the region

! - poly type

  !---- atom-based control
  IF (iftcratom) then

     DO ipoly = 1, npoly
        j1 = molept_index(ipoly)
        j2 = molept_index(ipoly+1) - 1
        do j = j1, j2
           atm_index = molept_list(j)

           molegrav(1:3) = atmcor(1:3,atm_index)

!          - P.B.C.
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

!          - if this molecule within temp. control region
           do iltr = 1, nlangeregion

              if (((ltxpos1(iltr) <= molegrav(1)) .and. &
                   & (molegrav(1) < ltxpos2(iltr))) .or. &
                   & ((ltxpos1(iltr) <= molegrav(1) - xcel) .and. &
                   & (molegrav(1) - xcel < ltxpos2(iltr)))) then ! xcel
              if (((ltypos1(iltr) <= molegrav(2)) .and. &
                   & (molegrav(2) < ltypos2(iltr))) .or. &
                   & ((ltypos1(iltr) <= molegrav(2) - ycel) .and. &
                   & (molegrav(2) - ycel < ltypos2(iltr)))) then ! ycel
              if (((ltzpos1(iltr) <= molegrav(3)) .and. &
                   & (molegrav(3) < ltzpos2(iltr))) .or. &
                   & ((ltzpos1(iltr) <= molegrav(3) - zcel) .and. &
                   & (molegrav(3) - zcel < ltzpos2(iltr)))) then ! zcel

!                 degfree_lt(iltr) = degfree_lt(iltr) + 3   ! atom-based control
                 natom_lt(iltr) = natom_lt(iltr) + 1
                 atmindex_lt(natom_lt(iltr),iltr) = atm_index

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

!          - calculate center of mass
           do j = 1, npoly_atom(ipolytyp)
              atm_index = atm_index + 1
              molemass = molemass + atmmass(atm_index)
              molegrav(1:3) = molegrav(1:3) &
                   &        + atmmass(atm_index)*atmcor(1:3,atm_index)
           end do

           molegrav(1:3) = molegrav(1:3) / molemass

!          - P.B.C.
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

!          - if this molecule within temp. control region
           do iltr = 1, nlangeregion

              if (((ltxpos1(iltr) <= molegrav(1)) .and. &
                   & (molegrav(1) < ltxpos2(iltr))) .or. &
                   & ((ltxpos1(iltr) <= molegrav(1) - xcel) .and. &
                   & (molegrav(1) - xcel < ltxpos2(iltr)))) then ! xcel
              if (((ltypos1(iltr) <= molegrav(2)) .and. &
                   & (molegrav(2) < ltypos2(iltr))) .or. &
                   & ((ltypos1(iltr) <= molegrav(2) - ycel) .and. &
                   & (molegrav(2) - ycel < ltypos2(iltr)))) then ! ycel
              if (((ltzpos1(iltr) <= molegrav(3)) .and. &
                   & (molegrav(3) < ltzpos2(iltr))) .or. &
                   & ((ltzpos1(iltr) <= molegrav(3) - zcel) .and. &
                   & (molegrav(3) - zcel < ltzpos2(iltr)))) then ! zcel

                 ! degfree_lt(iltr) = degfree_lt(iltr) &
                 !      &           + degfree_poly(ipolytyp) &
                 !      &           / npoly_mole(ipolytyp) ! per molecule

                 atm_index_tmp = atm_index - npoly_atom(ipolytyp)
                 do j = 1, npoly_atom(ipolytyp)
                    atm_index_tmp = atm_index_tmp + 1
                    natom_lt(iltr) = natom_lt(iltr) + 1
                    atmindex_lt(natom_lt(iltr),iltr) = atm_index_tmp
                 end do

                 exit

              end if
              end if
              end if

           end do

        end do

     end do

  END IF

! - water type
  do iwater = npoly+1, npoly+nwater
     j1 = molept_index(iwater)
     j2 = molept_index(iwater+1) - 1

     molemass = 0.0d0
     molegrav(1:3) = 0.0d0

     do j = j1, j2
        atm_index = molept_list(j)

!       - calculate center of mass
        molemass = molemass + atmmass(atm_index)
        molegrav(1:3) = molegrav(1:3) &
             &        + atmmass(atm_index)*atmcor(1:3,atm_index)
     end do

     molegrav(1:3) = molegrav(1:3) / molemass

!    - P.B.C.
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

!    - if this molecule within temp. control region
     do iltr = 1, nlangeregion

        if (((ltxpos1(iltr) <= molegrav(1)) .and. &
             & (molegrav(1) < ltxpos2(iltr))) .or. &
             & ((ltxpos1(iltr) <= molegrav(1) - xcel) .and. &
             & (molegrav(1) - xcel < ltxpos2(iltr)))) then ! xcel
        if (((ltypos1(iltr) <= molegrav(2)) .and. &
             & (molegrav(2) < ltypos2(iltr))) .or. &
             & ((ltypos1(iltr) <= molegrav(2) - ycel) .and. &
             & (molegrav(2) - ycel < ltypos2(iltr)))) then ! ycel
        if (((ltzpos1(iltr) <= molegrav(3)) .and. &
             & (molegrav(3) < ltzpos2(iltr))) .or. &
             & ((ltzpos1(iltr) <= molegrav(3) - zcel) .and. &
             & (molegrav(3) - zcel < ltzpos2(iltr)))) then ! zcel

           ! degfree_lt(iltr) = degfree_lt(iltr) &
           !      &           + degfree_water / nwater ! per molecule

           do j = j1, j2
              atm_index = molept_list(j)

              natom_lt(iltr) = natom_lt(iltr) + 1
              atmindex_lt(natom_lt(iltr),iltr) = atm_index
           end do

           exit

        end if
        end if
        end if

     end do

  end do

! - matom type
!  do imatyp = 1, nmatyp
  do imatyp = npoly+nwater+1, npoly+nwater+nmatom
     j1 = molept_index(imatyp)
     j2 = molept_index(imatyp+1) - 1

     do j = j1, j2
        atm_index = molept_list(j)

        molegrav(1:3) = atmcor(1:3,atm_index)

!       - P.B.C.
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

!       - if this molecule within temp. control region
        do iltr = 1, nlangeregion

           if (((ltxpos1(iltr) <= molegrav(1)) .and. &
                & (molegrav(1) < ltxpos2(iltr))) .or. &
                & ((ltxpos1(iltr) <= molegrav(1) - xcel) .and. &
                & (molegrav(1) - xcel < ltxpos2(iltr)))) then ! xcel
           if (((ltypos1(iltr) <= molegrav(2)) .and. &
                & (molegrav(2) < ltypos2(iltr))) .or. &
                & ((ltypos1(iltr) <= molegrav(2) - ycel) .and. &
                & (molegrav(2) - ycel < ltypos2(iltr)))) then ! ycel
           if (((ltzpos1(iltr) <= molegrav(3)) .and. &
                & (molegrav(3) < ltzpos2(iltr))) .or. &
                & ((ltzpos1(iltr) <= molegrav(3) - zcel) .and. &
                & (molegrav(3) - zcel < ltzpos2(iltr)))) then ! zcel

              natom_lt(iltr) = natom_lt(iltr) + 1
              atmindex_lt(natom_lt(iltr),iltr) = atm_index
              ! degfree_lt(iltr) = degfree_lt(iltr) + 3   ! matom type

              exit

           end if
           end if
           end if

        end do

     end do

  end do

!     +     +     +     +     +     +     +

end subroutine langevin_regist
