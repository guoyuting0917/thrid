!*****************************
!*  ass_strmvel.f90 Ver.1.0  *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-02-17 18:31:01 gota>

subroutine ass_strmvel(npoly,npolytyp, &
     &                 npoly_mole,npoly_atom, &
     &                 nwater, &
     &                 nmatom,nmatyp,nmatomtyp, &
     &                 xcel,ycel,zcel, &
     &                 iftcratom, &
     &                 ifstrmvel)

  use md_global

  implicit none

!     subroutine to assign streaming velocity to each atom

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

  logical,intent(in):: ifstrmvel     ! flag to input and use streaming velocity

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

  integer:: ireg

!     +     +     +     +     +     +     +

!---- initialization of array

  atmvel_strm(1:3,1:natom) = 0.0d0

!---- ifstrmvel is not used, just return
  if (.not. ifstrmvel) return

!---- assign streaming velocity to each atom

  atm_index = 0

! - poly type

!---- atom-based control
  if (iftcratom) then

     do ipoly = 1, npoly
        j1 = molept_index(ipoly)
        j2 = molept_index(ipoly+1) - 1
        do j = j1, j2
           atm_index = molept_list(j)

           molegrav(3) = atmcor(3,atm_index)

!              - P.B.C.
           ! if (atmcor(1,atm_index) < 0.0d0) then
           !    molegrav(1) = atmcor(1,atm_index) + xcel
           ! else if (atmcor(1,atm_index) >= xcel) then
           !    molegrav(1) = atmcor(1,atm_index) - xcel
           ! end if

           ! if (atmcor(2,atm_index) < 0.0d0) then
           !    molegrav(2) = atmcor(2,atm_index) + ycel
           ! else if (atmcor(2,atm_index) >= ycel) then
           !    molegrav(2) = atmcor(2,atm_index) - ycel
           ! end if

           if (atmcor(3,atm_index) < 0.0d0) then
              molegrav(3) = atmcor(3,atm_index) + zcel
           else if (atmcor(3,atm_index) >= zcel) then
              molegrav(3) = atmcor(3,atm_index) - zcel
           end if

!          - if this molecule within defined region
           do ireg = 1, nvelregion

              if ((strmzpos1(ireg) <= molegrav(3)) .and. &
                   & (molegrav(3) < strmzpos2(ireg))) then

                 atmvel_strm(1:3,atm_index) = strmvel(1:3,ireg)

                 exit

              end if

           end do

        end do

     end do

!---- molecule-based control
  else

     do ipolytyp = 1, npolytyp

        do i = 1, npoly_mole(ipolytyp)
           molemass = 0.0d0
           molegrav(3) = 0.0d0

!          - calculate center of mass
           do j = 1, npoly_atom(ipolytyp)
              atm_index = atm_index + 1
              molemass = molemass + atmmass(atm_index)
              molegrav(3) = molegrav(3) &
                   &        + atmmass(atm_index)*atmcor(3,atm_index)
           end do

           molegrav(3) = molegrav(3) / molemass

!              - P.B.C.
           ! if (molegrav(1) < 0.0d0) then
           !    molegrav(1) = molegrav(1) + xcel
           ! else if (molegrav(1) >= xcel) then
           !    molegrav(1) = molegrav(1) - xcel
           ! end if

           ! if (molegrav(2) < 0.0d0) then
           !    molegrav(2) = molegrav(2) + ycel
           ! else if (molegrav(2) >= ycel) then
           !    molegrav(2) = molegrav(2) - ycel
           ! end if

           if (molegrav(3) < 0.0d0) then
              molegrav(3) = molegrav(3) + zcel
           else if (molegrav(3) >= zcel) then
              molegrav(3) = molegrav(3) - zcel
           end if

!          - if this molecule within defined region
           do ireg = 1, nvelregion

              if ((strmzpos1(ireg) <= molegrav(3)) .and. &
                   & (molegrav(3) < strmzpos2(ireg))) then

                 atm_index_tmp = atm_index - npoly_atom(ipolytyp)
                 do j = 1, npoly_atom(ipolytyp)
                    atm_index_tmp = atm_index_tmp + 1

                    atmvel_strm(1:3,atm_index_tmp) = strmvel(1:3,ireg)

                 end do

                 exit

              end if

           end do

        end do

     end do

  end if

!     - water type
!      do iwater = 1, nwater
  do iwater = npoly+1, npoly+nwater
     j1 = molept_index(iwater)
     j2 = molept_index(iwater+1) - 1

     molemass = 0.0d0
     molegrav(3) = 0.0d0

     do j = j1, j2
        atm_index = molept_list(j)

!       - calculate center of mass
        molemass = molemass + atmmass(atm_index)
        molegrav(3) = molegrav(3) &
             &        + atmmass(atm_index)*atmcor(3,atm_index)
     end do

     molegrav(3) = molegrav(3) / molemass

!        - P.B.C.
     ! if (molegrav(1) < 0.0d0) then
     !    molegrav(1) = molegrav(1) + xcel
     ! else if (molegrav(1) >= xcel) then
     !    molegrav(1) = molegrav(1) - xcel
     ! end if

     ! if (molegrav(2) < 0.0d0) then
     !    molegrav(2) = molegrav(2) + ycel
     ! else if (molegrav(2) >= ycel) then
     !    molegrav(2) = molegrav(2) - ycel
     ! end if

     if (molegrav(3) < 0.0d0) then
        molegrav(3) = molegrav(3) + zcel
     else if (molegrav(3) >= zcel) then
        molegrav(3) = molegrav(3) - zcel
     end if

!    - if this molecule within defined region
     do ireg = 1, nvelregion

        if ((strmzpos1(ireg) <= molegrav(3)) .and. &
             & (molegrav(3) < strmzpos2(ireg))) then

           do j = j1, j2
              atm_index = molept_list(j)

              atmvel_strm(1:3,atm_index) = strmvel(1:3,ireg)

           end do

           exit

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

        molegrav(3) = atmcor(3,atm_index)

!           - P.B.C.
        ! if (atmcor(1,atm_index) < 0.0d0) then
        !    molegrav(1) = atmcor(1,atm_index) + xcel
        ! else if (atmcor(1,atm_index) >= xcel) then
        !    molegrav(1) = atmcor(1,atm_index) - xcel
        ! end if

        ! if (atmcor(2,atm_index) < 0.0d0) then
        !    molegrav(2) = atmcor(2,atm_index) + ycel
        ! else if (atmcor(2,atm_index) >= ycel) then
        !    molegrav(2) = atmcor(2,atm_index) - ycel
        ! end if

        if (atmcor(3,atm_index) < 0.0d0) then
           molegrav(3) = atmcor(3,atm_index) + zcel
        else if (atmcor(3,atm_index) >= zcel) then
           molegrav(3) = atmcor(3,atm_index) - zcel
        end if

!           - if this molecule within temp. control region
        do ireg = 1, nvelregion

           if ((strmzpos1(ireg) <= molegrav(3)) .and. &
                & (molegrav(3) < strmzpos2(ireg))) then

              atmvel_strm(1:3,atm_index) = strmvel(1:3,ireg)

              exit

           end if

        end do

     end do

  end do

!     +     +     +     +     +     +     +

end subroutine ass_strmvel
