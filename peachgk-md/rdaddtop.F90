!*****************************
!*  rdaddtop.f90 Ver.1.2     *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2014-06-13 16:37:23 gota>

subroutine rdaddtop(iuaddtop, &
     &              mref, &
     &              ifsetatmmass,ifsetatmchrg)

  use interface_tools

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: iuaddtop        ! input additional topology file unit

  real(8),intent(in):: mref            ! mass base value [kg] (C atom)

!     OUTPUT
  logical,intent(out):: ifsetatmmass(:) ! flag for setting mass by add_top
  logical,intent(out):: ifsetatmchrg(:) ! flag for setting charge by add_top

! LOCAL:
  character(80):: fredat(20)
  integer:: nword

  integer:: atom_index
  real(8):: atom_mass
  real(8):: atom_chrg
  integer:: bond_index
  integer:: angl_index
  integer:: anglub_index
  integer:: tors_index
  integer:: torsrb_index
  integer:: torsim_index

  integer:: nbond_old
  integer:: nangl_old
  integer:: nanglub_old
  integer:: ntors_old
  integer:: ntorsrb_old
  integer:: ntorsim_old

  logical:: ifoverwrite
  integer:: ibond_tmp,jbond_tmp
  integer:: iangl_tmp,jangl_tmp,kangl_tmp
  integer:: ianglub_tmp,janglub_tmp,kanglub_tmp
  integer:: itors_tmp,jtors_tmp,ktors_tmp,ltors_tmp
  integer:: itorsrb_tmp,jtorsrb_tmp,ktorsrb_tmp,ltorsrb_tmp
  integer:: itorsim_tmp,jtorsim_tmp,ktorsim_tmp,ltorsim_tmp

  real(8):: an = 6.0221367d+23  ! Avogadro's number

  integer:: i

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- initialization

  nbond_old = nbond
  nangl_old = nangl
  nanglub_old = nanglub
  ntors_old = ntors
  ntorsrb_old = ntorsrb
  ntorsim_old = ntorsim

!-------- Poly topology --------

  DOREAD:DO
     call rdfree(iuaddtop,20,fredat)

     if (fredat(1) == '<END>') then
        exit
     else if (fredat(1) == ' ') then
        cycle DOREAD

     else if (fredat(1) == '<ATOM>') then
        do
           call rdfree_w(iuaddtop,20,fredat,nword)
           if (fredat(1) == ' ') cycle DOREAD
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           !!! Current version does not support adding new atoms,
           !!! but only overwriting.
           read(fredat(1),*) atom_index
           atmtyp(atom_index) = fredat(2)(1:2)

           !---- read overwrite mass data
           if (nword >= 3) then
              read(fredat(3),*) atom_mass
              ifsetatmmass(atom_index) = .true.
              atmmass(atom_index) = atom_mass * 1.0d-3/an / mref
           endif

           !---- read overwrite charge data
           if (nword >= 4) then
              read(fredat(4),*) atom_chrg
              ifsetatmchrg(atom_index) = .true.
              atmchrg(atom_index) = atom_chrg
           endif

        end do

     else if (fredat(1) == '<BOND>') then
        do
           call rdfree(iuaddtop,20,fredat)
           if (fredat(1) == ' ') cycle DOREAD
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           read(fredat(1),*) bond_index
           read(fredat(2),*) ibond_tmp
           read(fredat(3),*) jbond_tmp

           ifoverwrite = .false.
           if (bond_index <= nbond_old) then ! force to overwrite
              ifoverwrite = .true.
           else ! search existing bond pair
              do i = 1, nbond_old
                 if (((ibond_tmp == ibond(i)) .and. (jbond_tmp == jbond(i))) &
              & .or. ((ibond_tmp == jbond(i)) .and. (jbond_tmp == ibond(i)))) then
                    ifoverwrite = .true.
                    bond_index = i
                 end if
              end do
           end if

           if (.not. ifoverwrite) then ! when adding new bond pair
              nbond = nbond + 1
              ibond(nbond) = ibond_tmp
              jbond(nbond) = jbond_tmp
              bondtyp(nbond) = fredat(4)(1:5)
           else ! when overwriting existing bond pair
              ibond(bond_index) = ibond_tmp
              jbond(bond_index) = jbond_tmp
              bondtyp(bond_index) = fredat(4)(1:5)
           end if

        end do

     else if (fredat(1) == '<ANGLE>') then
        do
           call rdfree(iuaddtop,20,fredat)
           if (fredat(1) == ' ') cycle DOREAD
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           read(fredat(1),*) angl_index
           read(fredat(2),*) iangl_tmp
           read(fredat(3),*) jangl_tmp
           read(fredat(4),*) kangl_tmp

           ifoverwrite = .false.
           if (angl_index <= nangl_old) then ! force to overwrite
              ifoverwrite = .true.
           else ! search existing angl pair
              do i = 1, nangl_old
                 if (((iangl_tmp == iangl(i)) .and. (jangl_tmp == jangl(i)) &
              & .and. (kangl_tmp == kangl(i))) &
              & .or. ((iangl_tmp == kangl(i)) .and. (jangl_tmp == jangl(i)) &
              & .and. (kangl_tmp == iangl(i)))) then
                    ifoverwrite = .true.
                    angl_index = i
                 end if
              end do
           end if

           if (.not. ifoverwrite) then ! when adding new angl pair
              nangl = nangl + 1
              iangl(nangl) = iangl_tmp
              jangl(nangl) = jangl_tmp
              kangl(nangl) = kangl_tmp
              angltyp(nangl) = fredat(5)(1:8)
           else ! when overwriting existing angl pair
              iangl(angl_index) = iangl_tmp
              jangl(angl_index) = jangl_tmp
              kangl(angl_index) = kangl_tmp
              angltyp(angl_index) = fredat(5)(1:8)
           end if

        end do

     else if (fredat(1) == '<ANGLE_UB>') then
        do
           call rdfree(iuaddtop,20,fredat)
           if (fredat(1) == ' ') cycle DOREAD
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           read(fredat(1),*) anglub_index
           read(fredat(2),*) ianglub_tmp
           read(fredat(3),*) janglub_tmp
           read(fredat(4),*) kanglub_tmp

           ifoverwrite = .false.
           if (anglub_index <= nanglub_old) then ! force to overwrite
              ifoverwrite = .true.
           else ! search existing anglub pair
              do i = 1, nanglub_old
                 if (((ianglub_tmp == ianglub(i)) &
              & .and. (janglub_tmp == janglub(i)) &
              & .and. (kanglub_tmp == kanglub(i))) &
              & .or. ((ianglub_tmp == kanglub(i)) &
              & .and. (janglub_tmp == janglub(i)) &
              & .and. (kanglub_tmp == ianglub(i)))) then
                    ifoverwrite = .true.
                    anglub_index = i
                 end if
              end do
           end if

           if (.not. ifoverwrite) then ! when adding new anglub pair
              nanglub = nanglub + 1
              ianglub(nanglub) = ianglub_tmp
              janglub(nanglub) = janglub_tmp
              kanglub(nanglub) = kanglub_tmp
              anglubtyp(nanglub) = fredat(5)(1:8)
           else ! when overwriting existing anglub pair
              ianglub(anglub_index) = ianglub_tmp
              janglub(anglub_index) = janglub_tmp
              kanglub(anglub_index) = kanglub_tmp
              anglubtyp(anglub_index) = fredat(5)(1:8)
           end if

        end do

     else if (fredat(1) == '<TORSION>') then
        do
           call rdfree(iuaddtop,20,fredat)
           if (fredat(1) == ' ') cycle DOREAD
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           read(fredat(1),*) tors_index
           read(fredat(2),*) itors_tmp
           read(fredat(3),*) jtors_tmp
           read(fredat(4),*) ktors_tmp
           read(fredat(5),*) ltors_tmp

           ifoverwrite = .false.
           if (tors_index <= ntors_old) then ! force to overwrite
              ifoverwrite = .true.
           else ! search existing tors pair
              do i = 1, ntors_old
                 if (((itors_tmp == itors(i)) .and. (jtors_tmp == jtors(i)) &
              & .and. (ktors_tmp == ktors(i)) .and. (ltors_tmp == ltors(i))) &
              & .or. ((itors_tmp == ltors(i)) .and. (jtors_tmp == ktors(i)) &
              & .and. (ktors_tmp == jtors(i)) .and. (ltors_tmp == itors(i)))) then
                    ifoverwrite = .true.
                    tors_index = i
                 end if
              end do
           end if

           if (.not. ifoverwrite) then ! when adding new tors pair
              ntors = ntors + 1
              itors(ntors) = itors_tmp
              jtors(ntors) = jtors_tmp
              ktors(ntors) = ktors_tmp
              ltors(ntors) = ltors_tmp
              torstyp(ntors) = fredat(6)(1:11)
           else ! when overwriting existing tors pair
              itors(tors_index) = itors_tmp
              jtors(tors_index) = jtors_tmp
              ktors(tors_index) = ktors_tmp
              ltors(tors_index) = ltors_tmp
              torstyp(tors_index) = fredat(6)(1:11)
           end if

        end do

     else if (fredat(1) == '<TORSION_RB>') then
        do
           call rdfree(iuaddtop,20,fredat)
           if (fredat(1) == ' ') cycle DOREAD
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           read(fredat(1),*) torsrb_index
           read(fredat(2),*) itorsrb_tmp
           read(fredat(3),*) jtorsrb_tmp
           read(fredat(4),*) ktorsrb_tmp
           read(fredat(5),*) ltorsrb_tmp

           ifoverwrite = .false.
           if (torsrb_index <= ntorsrb_old) then ! force to overwrite
              ifoverwrite = .true.
           else ! search existing tors pair
              do i = 1, ntorsrb_old
                 if (((itorsrb_tmp == itorsrb(i)) .and. (jtorsrb_tmp == jtorsrb(i)) &
              & .and. (ktorsrb_tmp == ktorsrb(i)) .and. (ltorsrb_tmp == ltorsrb(i))) &
              & .or. ((itorsrb_tmp == ltorsrb(i)) .and. (jtorsrb_tmp == ktorsrb(i)) &
              & .and. (ktorsrb_tmp == jtorsrb(i)) .and. (ltorsrb_tmp == itorsrb(i)))) then
                    ifoverwrite = .true.
                    torsrb_index = i
                 end if
              end do
           end if

           if (.not. ifoverwrite) then ! when adding new torsrb pair
              ntorsrb = ntorsrb + 1
              itorsrb(ntorsrb) = itorsrb_tmp
              jtorsrb(ntorsrb) = jtorsrb_tmp
              ktorsrb(ntorsrb) = ktorsrb_tmp
              ltorsrb(ntorsrb) = ltorsrb_tmp
              torsrbtyp(ntorsrb) = fredat(6)(1:11)
           else ! when overwriting existing torsrb pair
              itorsrb(torsrb_index) = itorsrb_tmp
              jtorsrb(torsrb_index) = jtorsrb_tmp
              ktorsrb(torsrb_index) = ktorsrb_tmp
              ltorsrb(torsrb_index) = ltorsrb_tmp
              torsrbtyp(torsrb_index) = fredat(6)(1:11)
           end if

        end do

     else if (fredat(1) == '<TORSION_IM>') then
        do
           call rdfree(iuaddtop,20,fredat)
           if (fredat(1) == ' ') cycle DOREAD
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           read(fredat(1),*) torsim_index
           read(fredat(2),*) itorsim_tmp
           read(fredat(3),*) jtorsim_tmp
           read(fredat(4),*) ktorsim_tmp
           read(fredat(5),*) ltorsim_tmp

           ifoverwrite = .false.
           if (torsim_index <= ntorsim_old) then ! force to overwrite
              ifoverwrite = .true.
           else ! search existing tors pair
              do i = 1, ntorsim_old
                 if (((itorsim_tmp == itorsim(i)) .and. (jtorsim_tmp == jtorsim(i)) &
              & .and. (ktorsim_tmp == ktorsim(i)) .and. (ltorsim_tmp == ltorsim(i))) &
              & .or. ((itorsim_tmp == ltorsim(i)) .and. (jtorsim_tmp == ktorsim(i)) &
              & .and. (ktorsim_tmp == jtorsim(i)) .and. (ltorsim_tmp == itorsim(i)))) then
                    ifoverwrite = .true.
                    torsim_index = i
                 end if
              end do
           end if

           if (.not. ifoverwrite) then ! when adding new torsim pair
              ntorsim = ntorsrb + 1
              itorsim(ntorsim) = itorsim_tmp
              jtorsim(ntorsim) = jtorsim_tmp
              ktorsim(ntorsim) = ktorsim_tmp
              ltorsim(ntorsim) = ltorsim_tmp
              torsimtyp(ntorsim) = fredat(6)(1:11)
           else ! when overwriting existing torsrb pair
              itorsim(torsim_index) = itorsim_tmp
              jtorsim(torsim_index) = jtorsim_tmp
              ktorsim(torsim_index) = ktorsim_tmp
              ltorsim(torsim_index) = ltorsim_tmp
              torsimtyp(torsim_index) = fredat(6)(1:11)
           end if

        end do

     end if

  END DO DOREAD

#if defined(_RDADDTOP_DEBUG)
  write(6,*) '***** rdaddtop DEBUG information *****'
  write(6,*) 'total number of bonds:', nbond
  write(6,*)
  do i = 1, nbond
     write(6,'(3I4,A6)') i,ibond(i),jbond(i),bondtyp(i)
  end do
  write(6,*) '**************************************'
#endif

  close(iuaddtop)

!     +     +     +     +     +     +     +

end subroutine rdaddtop
