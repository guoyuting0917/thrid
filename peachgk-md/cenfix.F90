!*****************************
!*  cenfix.f90 Ver.2.2       *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>
subroutine cenfix(ifcenterfix_all, &
     &            ifcenterfix_poly, &
     &            ifcenterfix_water, &
     &            ifcenterfix_ma, &
     &            cenfix_free, &
     &            ifcenterfix_polytyp, &
     &            ifcenterfix_watertyp, &
     &            ifcenterfix_matyp, &
     &            npoly,npolytyp,npoly_mole,npoly_atom, &
     &            nwater, &
     &            nmatom,nmatyp,nmatomtyp, &
     &            iflocalfix, &
     &            nlfix,index_nlfix,   &
     &            iflocalfixz,iflocalfixzg)

  use md_global

  implicit none

!     This subroutine fixes the motion of centroid to 0

! ARGUMENT:
!     INPUT
  logical,intent(in):: ifcenterfix_all   ! center fix for all
  logical,intent(in):: ifcenterfix_poly  ! center fix for polymer
  logical,intent(in):: ifcenterfix_water ! center fix for water
  logical,intent(in):: ifcenterfix_ma    ! center fix for monatomic mole.
  character(4),intent(in):: cenfix_free  ! COM not fixed in this direction

  logical,intent(in):: ifcenterfix_polytyp(:) ! center fix for each polymer
  logical,intent(in):: ifcenterfix_watertyp   ! center fix for each water
  logical,intent(in):: ifcenterfix_matyp(:)
                                          ! center fix for each monatomic mole.

  integer,intent(in):: npoly           ! all number of poly
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatom          ! number of monatomic molecules
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  logical,intent(in):: iflocalfix      ! fix atoms flag

  integer,intent(in):: nlfix           ! number of fix atoms
  integer,intent(in):: index_nlfix(:)  ! index of fix atoms

  logical,intent(in):: iflocalfixz     ! flag for fixing z coordinate of atoms

  logical,intent(in):: iflocalfixzg    ! flag for fixing z coordinate of COM of molecules

! LOCAL:
  integer:: i,j             ! do loop index
  integer:: i1,i2,jj        ! temporary valiables
!      integer:: wc              ! temporary counter
  real(8):: all_mass(3)         ! summation of all atom mass
  real(8):: vcent(3)         ! current vel of centroid

!  real(8):: all_mass_z       ! summation of all atom mass (for z coordinate)

  integer:: ilfix

  integer:: ipoly, imatom

  integer:: index_mole, index_mole2

!     +     +     +     +     +     +     +

!---- fix velocity of COM of whole system
  IF (ifcenterfix_all) THEN
     all_mass(1:3) = 0.0d0
!     all_mass_z = 0.0d0
     vcent(1:3) = 0.0d0

     do i=1,npoly+nwater+nmatom

        i1 = molept_index(i)
        i2 = molept_index(i+1)-1
        do j=i1,i2
           jj = molept_list(j)
           if( atmtyp(jj) == 'VW' ) then   ! for virtual wall degree (VW)
                                           ! only z has to be treated
              vcent(3) = vcent(3) + atmmass(jj)*atmvel(3,jj)
              all_mass(3) = all_mass(3) + atmmass(jj)
           else   ! normal case
              vcent(1:3) = vcent(1:3) + atmmass(jj)*atmvel(1:3,jj)
              all_mass(1:3) = all_mass(1:3) + atmmass(jj)
           end if
        end do

     end do

!    - local fix in z-coordinate -
     if (iflocalfixz) then
        do i=1, nlfixz
           jj = index_nlfixz(i)
           all_mass(3) = all_mass(3) - atmmass(jj)
        end do
     end if

!    - local fix in z-coordinate of COM -
     if (iflocalfixzg) then

        do i=1, nlfixzg
           i1 = molept_index(index_nlfixzg(i))
           i2 = molept_index(index_nlfixzg(i)+1)-1
           do j = i1, i2
              jj = molept_list(j)
              all_mass(3) = all_mass(3) - atmmass(jj)
           end do
        end do

     end if

!    - local fix -
     if (iflocalfix) then
        do ilfix = 1, nlfix
           i = index_nlfix(ilfix)
           all_mass(1:3) = all_mass(1:3) - atmmass(i)
        end do
     end if

     do i = 1,3
        if (all_mass(i) < 1.0d-32) then   ! if all_mass is minute
           vcent(i) = 0.0d0
        else
           vcent(i) = vcent(i) / all_mass(i)
        end if
     end do

!    - correct velocity of all atoms
     do i=1,npoly+nwater+nmatom

        i1 = molept_index(i)
        i2 = molept_index(i+1)-1
        do j=i1,i2
           jj = molept_list(j)
           if (atmtyp(jj) == 'VW') then   ! for virtual wall degree (VW)
                                           ! only z has to be treated
              atmvel(3,jj) = atmvel(3,jj) - vcent(3)
           else   ! normal case
              if (cenfix_free == 'none') then
                 atmvel(1:3,jj) = atmvel(1:3,jj) - vcent(1:3)
              else if (cenfix_free == 'x') then
                 atmvel(2:3,jj) = atmvel(2:3,jj) - vcent(2:3)
              else if (cenfix_free == 'y') then
                 atmvel(1,jj) = atmvel(1,jj) - vcent(1)
                 atmvel(3,jj) = atmvel(3,jj) - vcent(3)
              else if (cenfix_free == 'z') then
                 atmvel(1:2,jj) = atmvel(1:2,jj) - vcent(1:2)
             end if

           end if
        end do

     end do

  ELSE
!---- fix velocity of COM of each component
     ! polymer type
     if (ifcenterfix_poly .and. (npoly /= 0)) then

        index_mole = 0
        index_mole2 = index_mole

        do ipoly=1, npolytyp

           if (.not. ifcenterfix_polytyp(ipoly)) then
              index_mole = index_mole + npoly_mole(ipoly)
              index_mole2 = index_mole2 + npoly_mole(ipoly)
              cycle
           end if

           all_mass(1:3) = 0.0d0
           vcent(1:3) = 0.0d0

           do i=1, npoly_mole(ipoly)
              index_mole = index_mole + 1

              i1 = molept_index(index_mole)
              i2 = molept_index(index_mole+1)-1
              do j=i1,i2
                 jj = molept_list(j)

                 vcent(1:3) = vcent(1:3) + atmmass(jj)*atmvel(1:3,jj)
                 all_mass(1:3) = all_mass(1:3) + atmmass(jj)

!!! '11.02.08    This is very slow process, and should be improved...
                 !    - local fix in z-coordinate -
                 if (iflocalfixz) then
                    do ilfix=1, nlfixz
                       if (jj == index_nlfixz(ilfix)) then
                          all_mass(3) = all_mass(3) - atmmass(jj)
                          exit
                       end if
                    end do
                 end if

                 !    - local fix in z-coordinate of COM -
                 if (iflocalfixzg) then
                    do ilfix=1, nlfixzg
                       if (index_mole == index_nlfixzg(ilfix)) then
                          all_mass(3) = all_mass(3) - atmmass(jj)
                          exit
                       end if
                    end do
                 end if

                 !    - local fix -
                 if (iflocalfix) then
                    do ilfix=1, nlfix
                       if (jj == index_nlfix(ilfix)) then
                          all_mass(1:3) = all_mass(1:3) - atmmass(jj)
                          exit
                       end if
                    end do
                 end if
!!! '11.02.08    till here

              end do

           end do

           do i = 1,3
              if (all_mass(i) < 1.0d-32) then   ! if all_mass is minute
                 vcent(i) = 0.0d0
              else
                 vcent(i) = vcent(i) / all_mass(i)
              end if
           end do

!          - correct velocity of polymeric type
           do i=1, npoly_mole(ipoly)
              index_mole2 = index_mole2 + 1

              i1 = molept_index(index_mole2)
              i2 = molept_index(index_mole2+1)-1
              do j=i1,i2
                 jj = molept_list(j)

                 atmvel(1:3,jj) = atmvel(1:3,jj) - vcent(1:3)
              end do

           end do

        end do

     end if

     ! water type
     if (ifcenterfix_water .and. (nwater /= 0)   &
          & .and. ifcenterfix_watertyp) then

        index_mole = npoly
        index_mole2 = index_mole

        all_mass(1:3) = 0.0d0
        vcent(1:3) = 0.0d0

        do i=1, nwater
           index_mole = index_mole + 1

           i1 = molept_index(index_mole)
           i2 = molept_index(index_mole+1)-1
           do j=i1,i2
              jj = molept_list(j)

              vcent(1:3) = vcent(1:3) + atmmass(jj)*atmvel(1:3,jj)
              all_mass(1:3) = all_mass(1:3) + atmmass(jj)

!!! '11.02.08 This is very slow process, and should be improved...
              !    - local fix in z-coordinate -
              if (iflocalfixz) then
                 do ilfix=1, nlfixz
                    if (jj == index_nlfixz(ilfix)) then
                       all_mass(3) = all_mass(3) - atmmass(jj)
                       exit
                    end if
                 end do
              end if

              !    - local fix in z-coordinate of COM -
              if (iflocalfixzg) then
                 do ilfix=1, nlfixzg
                    if (index_mole == index_nlfixzg(ilfix)) then
                       all_mass(3) = all_mass(3) - atmmass(jj)
                       exit
                    end if
                 end do
              end if

              !    - local fix -
              if (iflocalfix) then
                 do ilfix=1, nlfix
                    if (jj == index_nlfix(ilfix)) then
                       all_mass(1:3) = all_mass(1:3) - atmmass(jj)
                       exit
                    end if
                 end do
              end if
!!! '11.02.08 till here

           end do

        end do

        do i = 1,3
           if (all_mass(i) < 1.0d-32) then   ! if all_mass is minute
              vcent(i) = 0.0d0
           else
              vcent(i) = vcent(i) / all_mass(i)
           end if
        end do

!       - correct velocity of water type
        do i=1, nwater
           index_mole2 = index_mole2 + 1

           i1 = molept_index(index_mole2)
           i2 = molept_index(index_mole2+1)-1
           do j=i1,i2
              jj = molept_list(j)

              atmvel(1:3,jj) = atmvel(1:3,jj) - vcent(1:3)
           end do

        end do

     end if

     ! matom type
     if (ifcenterfix_ma .and. (nmatom /= 0)) then

        index_mole = npoly + nwater
        index_mole2 = index_mole

        do imatom = 1, nmatyp

           if (.not. ifcenterfix_matyp(imatom)) then
              index_mole = index_mole + nmatomtyp(imatom)
              index_mole2 = index_mole2 + nmatomtyp(imatom)
              cycle
           end if

           all_mass(1:3) = 0.0d0
           vcent(1:3) = 0.0d0

           do i=1, nmatomtyp(imatom)
              index_mole = index_mole + 1

              i1 = molept_index(index_mole)
              i2 = molept_index(index_mole+1)-1
              do j=i1,i2
                 jj = molept_list(j)
                 if( atmtyp(jj) == 'VW' ) then   ! for virtual wall degree (VW)
                                                 ! only z has to be treated
                    vcent(3) = vcent(3) + atmmass(jj)*atmvel(3,jj)
                    all_mass(3) = all_mass(3) + atmmass(jj)
                 else   ! normal case
                    vcent(1:3) = vcent(1:3) + atmmass(jj)*atmvel(1:3,jj)
                    all_mass(1:3) = all_mass(1:3) + atmmass(jj)
                 end if

!!! '11.02.08    This is very slow process, and should be improved...
                 !    - local fix in z-coordinate -
                 if (iflocalfixz) then
                    do ilfix=1, nlfixz
                       if (jj == index_nlfixz(ilfix)) then
                          all_mass(3) = all_mass(3) - atmmass(jj)
                          exit
                       end if
                    end do
                 end if

                 !    - local fix in z-coordinate of COM -
                 if (iflocalfixzg) then
                    do ilfix=1, nlfixzg
                       if (index_mole == index_nlfixzg(ilfix)) then
                          all_mass(3) = all_mass(3) - atmmass(jj)
                          exit
                       end if
                    end do
                 end if

                 !    - local fix -
                 if (iflocalfix) then
                    do ilfix=1, nlfix
                       if (jj == index_nlfix(ilfix)) then
                          all_mass(1:3) = all_mass(1:3) - atmmass(jj)
                          exit
                       end if
                    end do
                 end if
!!! '11.02.08    till here

              end do

           end do

           do i = 1,3
              if (all_mass(i) < 1.0d-32) then   ! if all_mass is minute
                 vcent(i) = 0.0d0
              else
                 vcent(i) = vcent(i) / all_mass(i)
              end if
           end do

!          - correct velocity of monoatomic type
           do i=1, nmatomtyp(imatom)
              index_mole2 = index_mole2 + 1

              i1 = molept_index(index_mole2)
              i2 = molept_index(index_mole2+1)-1
              do j=i1,i2
                 jj = molept_list(j)

                 if( atmtyp(jj) == 'VW' ) then   ! for virtual wall degree (VW)
                                                 ! only z has to be treated
                    atmvel(3,jj) = atmvel(3,jj) - vcent(3)
                 else   ! normal case
                    atmvel(1:3,jj) = atmvel(1:3,jj) - vcent(1:3)
                 end if
              end do

           end do

        end do

     end if

  END IF

! - local fix in z-coordinate -
  if (iflocalfixz) then
     do i=1, nlfixz
        jj = index_nlfixz(i)
        atmvel(3,jj) = 0.0d0
     end do
  end if

!    - local fix in z-coordinate of COM -
  if (iflocalfixzg) then
     all_mass(3) = 0.0d0
     vcent(3) = 0.0d0
     do i=1, nlfixzg
        i1 = molept_index(index_nlfixzg(i))
        i2 = molept_index(index_nlfixzg(i)+1)-1

        do j = i1, i2
           jj = molept_list(j)
           vcent(3) = vcent(3) + atmmass(jj)*atmvel(3,jj)
           all_mass(3) = all_mass(3) + atmmass(jj)
        end do
        vcent(3) = vcent(3) / all_mass(3)

        do j = i1, i2
           jj = molept_list(j)
           atmvel(3,jj) = atmvel(3,jj) - vcent(3)
        end do

     end do
  end if

! - local fix
  if (iflocalfix) then
     do i = 1, nlfix
        jj = index_nlfix(i)
        atmvel(1:3,jj) = 0.0d0
     end do
  end if

!     +     +     +     +     +     +     +

end subroutine cenfix
