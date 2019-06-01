!*************************************
!*  linkbond.f90 Ver.2.3 '13.11.05   *
!*      for peachgk_md.f             *
!*            by G.Kikugawa          *
!*************************************
subroutine linkbond(npolytyp,npoly_mole,npoly_atom, &
     &              polytyp_free, &
     &              atmchrg_tmp, &
     &              ifsetatmmass,ifsetatmchrg)

  use md_global

  implicit none

! ARGUMENTS:
! INPUT
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

  character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
  real(8),intent(in):: atmchrg_tmp(:)  ! temporary atom charge read from cor file

  logical,intent(in):: ifsetatmmass(:) ! flag for setting mass by add_top
  logical,intent(in):: ifsetatmchrg(:) ! flag for setting charge by add_top

! LOCAL
  integer:: i,j             ! do loop index

  logical:: ifpara_exist    ! parameter check

  integer:: ipoly
  integer:: iword

  integer:: index_atom

!     +     +     +     +     +     +

!---- link atom
  do i = 1, natom

     ifpara_exist = .false.
     do j = 1, natmtyp
            
        if (atmtyp(i) == para_atmtyp(j)) then
           atmindex(i) = j
           if (.not. ifsetatmmass(i)) atmmass(i) = para_atmmass(j)
           if (.not. ifsetatmchrg(i)) atmchrg(i) = para_atmchrg(j)
           ifpara_exist = .true.
        end if

     end do
!        - error check -
     if (.not. ifpara_exist) then
        write(6,*) 'Error: atmtyp ',atmtyp(i),' does not exist'
        stop
     end if

  end do

!     --- set partial charge from cor file
  do ipoly = 1, npolytyp
     do iword = 1, maxnword

        if (polytyp_free(ipoly,iword) == 'setcharge') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           do i = 1, npoly_mole(ipoly)
              do j = 1, npoly_atom(ipoly)
                 index_atom = index_atom + 1
                 atmchrg(index_atom) = atmchrg_tmp(index_atom)
              end do
           end do

           exit
        end if

     end do
  end do

#if defined(_LINKBOND_DEBUG)
  write(6,*) '***** linkbond DEBUG information *****'
  write(6,*) '*** atom info ***'
  do i = 1,natom
     write(6,'(I8,A3,2E16.8)') i,atmtyp(i),atmmass(i),atmchrg(i)
  end do
  write(6,*) '**************************************'

#endif

!---- link bond
  do i = 1, nbond

     ifpara_exist = .false.
     do j=1,nbondtyp

        if (bondtyp(i) == para_bondtyp(j)) then
           indexbondtyp(i) = j
           ifpara_exist = .true.
        end if

     end do
!        - error check -
     if  (.not. ifpara_exist) then
        write(6,*) 'Error: bondtyp ',bondtyp(i),' does not exist'
        stop
     end if

  end do

!---- link angle
  do i = 1, nangl

     ifpara_exist = .false.
     do j=1,nangltyp

        if (angltyp(i) == para_angltyp(j)) then
           indexangltyp(i) = j
           ifpara_exist = .true.
        end if

     end do
!        - error check -
     if  (.not. ifpara_exist) then
        write(6,*) 'Error: angltyp ',angltyp(i),' does not exist'
        stop
     end if

  end do

!---- link Urey-Bradley angle
  do i = 1, nanglub

     ifpara_exist = .false.
     do j=1,nanglubtyp

        if (anglubtyp(i) == para_anglubtyp(j)) then
           indexanglubtyp(i) = j
           ifpara_exist = .true.
        end if

     end do
!        - error check -
     if  (.not. ifpara_exist) then
        write(6,*) 'Error: anglubtyp ',anglubtyp(i),' does not exist'
        stop
     end if

  end do

!---- link torsion
  do i = 1, ntors

     ifpara_exist = .false.
     do j = 1, ntorstyp

        if (torstyp(i) == para_torstyp(j)) then
           indextorstyp(i) = j
           ifpara_exist = .true.
        end if

     end do
!        - error check -
     if  (.not. ifpara_exist) then
        write(6,*) 'Error: torstyp ',torstyp(i),' does not exist'
        stop
     end if

  end do

!---- link torsionrb
  do i = 1, ntorsrb

     ifpara_exist = .false.
     do j = 1, ntorsrbtyp

        if (torsrbtyp(i) == para_torsrbtyp(j)) then
           indextorsrbtyp(i) = j
           ifpara_exist = .true.
        end if

     end do
!        - error check -
     if  (.not. ifpara_exist) then
        write(6,*) 'Error: torsrbtyp ',torsrbtyp(i), ' does not exist'
        stop
     end if

  end do

!---- link torsionim
  do i = 1, ntorsim

     ifpara_exist = .false.
     do j = 1, ntorsimtyp

        if (torsimtyp(i) == para_torsimtyp(j)) then
           indextorsimtyp(i) = j
           ifpara_exist = .true.
        end if

     end do
!        - error check -
     if  (.not. ifpara_exist) then
        write(6,*) 'Error: torsimtyp ',torsimtyp(i), ' does not exist'
        stop
     end if

  end do

!---- list of atoms of morse interaction
  nmorse = 0

  do i = 1, natom
     ifmorse(i) = .false.
     do j = 1, nmortyp

        if (atmtyp(i)(1:2) == para_mortyp(j)(1:2)) then
           nmorse = nmorse + 1
           nmorselist(nmorse) = i
           ifmorse(i) = .true.
           exit
        else if (atmtyp(i)(1:2) == para_mortyp(j)(4:5)) then
           nmorse = nmorse + 1
           nmorselist(nmorse) = i
           ifmorse(i) = .true.
           exit
        end if

     end do
  end do

!---- list of atoms of SH interaction
  nsho = 0
  nshh = 0

  do i = 1, natom
     do j = 1, nshtyp 

        if (atmtyp(i)(1:2) == para_shtyp(j)(1:2)) then
           if (j == 1) then               ! if j = 1, O atoms
              nsho = nsho + 1
              nsholist(nsho) = i
              exit
           else if (j == 2) then          ! if j = 2, H atoms
              nshh = nshh + 1
              nshhlist(nshh) = i
              exit
           end if
        else if (atmtyp(i)(1:2) == para_shtyp(j)(4:5)) then ! PT atoms
           nsho = nsho + 1
           nshh = nshh + 1
           nsholist(nsho) = i
           nshhlist(nshh) = i
           exit
        end if

     end do
  end do

!---- list of atoms of RFH interaction
  nrfh = 0

  do i = 1, natom
     do j = 1, nrfhtyp

        if (atmtyp(i)(1:2) == para_rfhtyp(j)(1:2)) then
           nrfh = nrfh + 1
           nrfhlist(nrfh) = i
           exit
        else if (atmtyp(i)(1:2) == para_rfhtyp(j)(4:5)) then
           nrfh = nrfh + 1
           nrfhlist(nrfh) = i
           exit
        end if

     end do
  end do

!---- list of atoms of DOU interaction
  ndou = 0

  do i = 1, natom
     do j = 1, ndoutyp 

        if (atmtyp(i)(1:2) == para_doutyp(j)(1:2)) then ! O or H atoms
           ndou = ndou + 1
           ndoulist(ndou) = i
           exit
        else if (atmtyp(i)(1:2) == para_doutyp(j)(4:5)) then ! Au atoms
           ndou = ndou + 1
           ndoulist(ndou) = i
           exit
        end if

     end do
  end do

!---- list of atoms of RP-VW interaction
  nrpvw = 0

  do i = 1, natom
     do j = 1, nrpvwtyp

        if (atmtyp(i)(1:2) == para_rpvwtyp(j)(1:2)) then
           nrpvw = nrpvw + 1
           nrpvwlist(nrpvw) = i
           exit
        else if (atmtyp(i)(1:2) == para_rpvwtyp(j)(4:5)) then
           nrpvw = nrpvw + 1
           nrpvwlist(nrpvw) = i
           exit
        end if

     end do
  end do

!---- list of atoms of VW

  j = 0

  do i = 1, natom

     if (atmtyp(i)(1:2) == 'VW') then
        j = j + 1
        nvwlist(j) = i
     end if

  end do

!     +     +     +     +     +     +

end subroutine linkbond
