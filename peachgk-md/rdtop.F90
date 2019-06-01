!*****************************
!*  rdtop.f90 Ver.2.1        *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>

subroutine rdtop(iutop,iuwtop, &
     &           npolytyp,npoly_mole,npoly_atom, &
     &           nwater)

  use interface_tools

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: iutop(:)        ! input poly topology file unit
  integer,intent(in):: iuwtop          ! input water topology file unit
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)
                          ! number of atoms belonging to poly
  integer,intent(in):: nwater          ! number of H2O molecules

! LOCAL:
  character(80):: fredat(20)
  integer:: i,j,m,n             ! do loop index

  integer:: nbondindex, nbondend
  integer:: nanglindex, nanglend
  integer:: nanglubindex, nanglubend
  integer:: ntorsindex, ntorsend
  integer:: ntorsrbindex, ntorsrbend
  integer:: ntorsimindex, ntorsimend

  integer:: natomindex

  integer:: ipoly
  integer:: iatom, latom

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- initialization
  nbond = 0
  nangl = 0
  nanglub = 0
  ntors = 0
  ntorsrb = 0
  ntorsim = 0

  natomindex = 0

!-------- Poly topology --------

  DO ipoly=1,npolytyp

     nbondindex = nbond
     nanglindex = nangl
     nanglubindex = nanglub
     ntorsindex = ntors
     ntorsrbindex = ntorsrb
     ntorsimindex = ntorsim

     DOREAD:DO
        call rdfree( iutop(ipoly), 20, fredat )

        if (fredat(1) == '<END>') then
           exit
        else if (fredat(1) == ' ') then
           cycle DOREAD

        else if (fredat(1) == '<BOND>') then
           do
              call rdfree(iutop(ipoly),20,fredat)
              if (fredat(1) == ' ') cycle DOREAD
              if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                               ! comment line

              nbond = nbond + 1
              if (nbond > maxnbond) then
                 write(6,*) 'Error : nbond exceeds maxnbond= ',maxnbond
                 stop
              end if

              read(fredat(2),*) ibond(nbond)
              read(fredat(3),*) jbond(nbond)
              bondtyp(nbond) = fredat(4)(1:5)
              ibond(nbond) = ibond(nbond) + natomindex
              jbond(nbond) = jbond(nbond) + natomindex

           end do

        else if (fredat(1) == '<ANGLE>') then
           do
              call rdfree(iutop(ipoly),20,fredat)
              if (fredat(1) == ' ') cycle DOREAD
              if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                               ! comment line

              nangl = nangl + 1
              if (nangl > maxnangl) then
                 write(6,*) 'Error : nangl exceeds maxnangl= ',maxnangl
                 stop
              end if

              read(fredat(2),*) iangl(nangl)
              read(fredat(3),*) jangl(nangl)
              read(fredat(4),*) kangl(nangl)
              angltyp(nangl) = fredat(5)(1:8)
              iangl(nangl) = iangl(nangl) + natomindex
              jangl(nangl) = jangl(nangl) + natomindex
              kangl(nangl) = kangl(nangl) + natomindex

           end do

        else if (fredat(1) == '<ANGLE_UB>') then
           do
              call rdfree(iutop(ipoly),20,fredat)
              if (fredat(1) == ' ') cycle DOREAD
              if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                               ! comment line

              nanglub = nanglub + 1
              if (nanglub > maxnangl) then
                 write(6,*) 'Error : nangl exceeds maxnangl= ',maxnangl
                 stop
              end if

              read(fredat(2),*) ianglub(nanglub)
              read(fredat(3),*) janglub(nanglub)
              read(fredat(4),*) kanglub(nanglub)
              anglubtyp(nanglub) = fredat(5)(1:8)
              ianglub(nanglub) = ianglub(nanglub) + natomindex
              janglub(nanglub) = janglub(nanglub) + natomindex
              kanglub(nanglub) = kanglub(nanglub) + natomindex

           end do

#if !defined(_CHARMM_NONB14)
        else if (fredat(1) == '<TORSION>') then
           do
              call rdfree(iutop(ipoly),20,fredat)
              if (fredat(1) == ' ') cycle DOREAD
              if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                               ! comment line

              ntors = ntors + 1
              if (ntors > maxntors) then
                 write(6,*) 'Error : ntors exceeds maxntors= ',maxntors
                 stop
              end if

              read(fredat(2),*) itors(ntors)
              read(fredat(3),*) jtors(ntors)
              read(fredat(4),*) ktors(ntors)
              read(fredat(5),*) ltors(ntors)
              torstyp(ntors) = fredat(6)(1:11)
              itors(ntors) = itors(ntors) + natomindex
              jtors(ntors) = jtors(ntors) + natomindex
              ktors(ntors) = ktors(ntors) + natomindex
              ltors(ntors) = ltors(ntors) + natomindex

           end do

#else
        else if (fredat(1) == '<TORSION_CH>') then
           do
              call rdfree(iutop(ipoly),20,fredat)
              if (fredat(1) == ' ') cycle DOREAD
              if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                               ! comment line

              ntors = ntors + 1
              if (ntors > maxntors) then
                 write(6,*) 'Error : ntors exceeds maxntors= ',maxntors
                 stop
              end if

              read(fredat(2),*) itors(ntors)
              read(fredat(3),*) jtors(ntors)
              read(fredat(4),*) ktors(ntors)
              read(fredat(5),*) ltors(ntors)
              torstyp(ntors) = fredat(6)(1:11)
              itors(ntors) = itors(ntors) + natomindex
              jtors(ntors) = jtors(ntors) + natomindex
              ktors(ntors) = ktors(ntors) + natomindex
              ltors(ntors) = ltors(ntors) + natomindex

           end do
#endif

        else if (fredat(1) == '<TORSION_RB>') then
           do
              call rdfree(iutop(ipoly),20,fredat)
              if (fredat(1) == ' ') cycle DOREAD
              if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                               ! comment line

              ntorsrb = ntorsrb + 1
              if (ntorsrb > maxntors) then
                 write(6,*) 'Error : ntorsrb exceeds maxntors= ',maxntors
                 stop
              end if

              read(fredat(2),*) itorsrb(ntorsrb)
              read(fredat(3),*) jtorsrb(ntorsrb)
              read(fredat(4),*) ktorsrb(ntorsrb)
              read(fredat(5),*) ltorsrb(ntorsrb)
              torsrbtyp(ntorsrb) = fredat(6)(1:11)
              itorsrb(ntorsrb) = itorsrb(ntorsrb) + natomindex
              jtorsrb(ntorsrb) = jtorsrb(ntorsrb) + natomindex
              ktorsrb(ntorsrb) = ktorsrb(ntorsrb) + natomindex
              ltorsrb(ntorsrb) = ltorsrb(ntorsrb) + natomindex

           end do

        else if (fredat(1) == '<TORSION_IM>') then
           do
              call rdfree(iutop(ipoly),20,fredat)
              if (fredat(1) == ' ') cycle DOREAD
              if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                               ! comment line

              ntorsim = ntorsim + 1
              if (ntorsim > maxntors) then
                 write(6,*) 'Error : ntorsim exceeds maxntors= ',maxntors
                 stop
              end if

              read(fredat(2),*) itorsim(ntorsim)
              read(fredat(3),*) jtorsim(ntorsim)
              read(fredat(4),*) ktorsim(ntorsim)
              read(fredat(5),*) ltorsim(ntorsim)
              torsimtyp(ntorsim) = fredat(6)(1:11)
              itorsim(ntorsim) = itorsim(ntorsim) + natomindex
              jtorsim(ntorsim) = jtorsim(ntorsim) + natomindex
              ktorsim(ntorsim) = ktorsim(ntorsim) + natomindex
              ltorsim(ntorsim) = ltorsim(ntorsim) + natomindex

           end do

        end if

     END DO DOREAD

!-------- duplicate topology to all poly --------

     nbondend = nbond - nbondindex
     nanglend = nangl - nanglindex
     nanglubend = nanglub - nanglubindex
     ntorsend = ntors - ntorsindex
     ntorsrbend = ntorsrb - ntorsrbindex
     ntorsimend = ntorsim - ntorsimindex

     do n=2,npoly_mole(ipoly)

        do m=1,nbondend
           nbond = nbond + 1
           if (nbond > maxnbond) then
              write(6,*) 'Error : nbond exceeds maxnbond= ',maxnbond
              stop
           end if

           ibond(nbond) = ibond(nbondindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           jbond(nbond) = jbond(nbondindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           bondtyp(nbond) = bondtyp(nbondindex+m)
        end do

        do m=1,nanglend
           nangl = nangl + 1
           if (nangl > maxnangl) then
              write(6,*) 'Error : nangl exceeds maxnangl= ',maxnangl
              stop
           end if

           iangl(nangl) = iangl(nanglindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           jangl(nangl) = jangl(nanglindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           kangl(nangl) = kangl(nanglindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           angltyp(nangl) = angltyp(nanglindex+m)
        end do

        do m=1,nanglubend
           nanglub = nanglub + 1
           if (nanglub > maxnangl) then
              write(6,*) 'Error : nanglub exceeds maxnangl= ',maxnangl
              stop
           end if

           ianglub(nanglub) = ianglub(nanglubindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           janglub(nanglub) = janglub(nanglubindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           kanglub(nanglub) = kanglub(nanglubindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           anglubtyp(nanglub) = anglubtyp(nanglubindex+m)
        end do

        do m=1,ntorsend
           ntors = ntors + 1
           if (ntors > maxntors) then
              write(6,*) 'Error : ntors exceeds maxntors= ',maxntors
              stop
           end if

           itors(ntors) = itors(ntorsindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           jtors(ntors) = jtors(ntorsindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           ktors(ntors) = ktors(ntorsindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           ltors(ntors) = ltors(ntorsindex+m) + &
                &         npoly_atom(ipoly)*(n-1)
           torstyp(ntors) = torstyp(ntorsindex+m)
        end do

        do m=1,ntorsrbend
           ntorsrb = ntorsrb + 1
           if (ntorsrb > maxntors) then
              write(6,*) 'Error : ntorsrb exceeds maxntors= ',maxntors
              stop
           end if

           itorsrb(ntorsrb) = itorsrb(ntorsrbindex+m) + &
                &             npoly_atom(ipoly)*(n-1)
           jtorsrb(ntorsrb) = jtorsrb(ntorsrbindex+m) + &
                &             npoly_atom(ipoly)*(n-1)
           ktorsrb(ntorsrb) = ktorsrb(ntorsrbindex+m) + &
                &             npoly_atom(ipoly)*(n-1)
           ltorsrb(ntorsrb) = ltorsrb(ntorsrbindex+m) + &
                &             npoly_atom(ipoly)*(n-1)
           torsrbtyp(ntorsrb) = torsrbtyp(ntorsrbindex+m)
        end do

        do m=1,ntorsimend
           ntorsim = ntorsim + 1
           if (ntorsim > maxntors) then
              write(6,*) 'Error : ntorsim exceeds maxntors= ',maxntors
              stop
           end if

           itorsim(ntorsim) = itorsim(ntorsimindex+m) + &
                &             npoly_atom(ipoly)*(n-1)
           jtorsim(ntorsim) = jtorsim(ntorsimindex+m) + &
                &             npoly_atom(ipoly)*(n-1)
           ktorsim(ntorsim) = ktorsim(ntorsimindex+m) + &
                &             npoly_atom(ipoly)*(n-1)
           ltorsim(ntorsim) = ltorsim(ntorsimindex+m) + &
                &             npoly_atom(ipoly)*(n-1)
           torsimtyp(ntorsim) = torsimtyp(ntorsimindex+m)
        end do

     end do

     natomindex = natomindex + npoly_mole(ipoly)*npoly_atom(ipoly)

     close(iutop(ipoly))

  END DO

! ---- detect ring structure
  ! 0.0 = within 4- or 5- membered ring
  ! 0.5 = within 6-membered ring
  ! 1.0 = others
  ! initialization to others
  ring14(1:ntors) = 1.0d0
  ring14_rb(1:ntors) = 1.0d0

  ! check periodic type torsion
  do i = 1, ntors
      iatom = itors(i)
      latom = ltors(i)
      ! find 4-membered ring
      do j = 1, nbond
          if ((iatom == ibond(j) .and. latom == jbond(j)) .or. &
          &   (iatom == jbond(j) .and. latom == ibond(j))) then
              ring14(i) = 0.0d0
          endif
      end do

      ! find 5-membered ring
      do j = 1, nangl
          if ((iatom == iangl(j) .and. latom == kangl(j)) .or. &
          &   (iatom == kangl(j) .and. latom == iangl(j))) then
              ring14(i) = 0.0d0
          endif
      end do

      ! find 6-membered ring
      do j = 1, ntors
          if (i == j) cycle   ! skip the same torsion
          if ((iatom == itors(j) .and. latom == ltors(j)) .or. &
          &   (iatom == ltors(j) .and. latom == itors(j))) then
              ring14(i) = 0.5d0
          endif
      end do

  end do

  ! check RB type torsion
  do i = 1, ntorsrb
      iatom = itorsrb(i)
      latom = ltorsrb(i)
      ! find 4-membered ring
      do j = 1, nbond
          if ((iatom == ibond(j) .and. latom == jbond(j)) .or. &
          &   (iatom == jbond(j) .and. latom == ibond(j))) then
              ring14_rb(i) = 0.0d0
          endif
      end do

      ! find 5-membered ring
      do j = 1, nangl
          if ((iatom == iangl(j) .and. latom == kangl(j)) .or. &
          &   (iatom == kangl(j) .and. latom == iangl(j))) then
              ring14_rb(i) = 0.0d0
          endif
      end do

      ! find 6-membered ring
      do j = 1, ntors
          if (i == j) cycle   ! skip the same torsion
          if ((iatom == itors(j) .and. latom == ltors(j)) .or. &
          &   (iatom == ltors(j) .and. latom == itors(j))) then
              ring14_rb(i) = 0.5d0
          endif
      end do

  end do

!-------- H2O topology --------

  nbondindex = nbond
  nanglindex = nangl

  if (nwater == 0) return
  DOREADW:DO
     call rdfree( iuwtop, 20, fredat )

     if (fredat(1) == '<END>') then
        exit
     else if (fredat(1) == ' ') then
        cycle DOREADW

     else if (fredat(1) == '<BOND>') then
        do
           call rdfree( iuwtop, 20, fredat )
           if (fredat(1) == ' ') cycle DOREADW
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nbond = nbond + 1
           if (nbond > maxnbond) then
              write(6,*) 'Error : nbond exceeds maxnbond= ',maxnbond
              stop
           end if

           read(fredat(2),*) ibond(nbond)
           read(fredat(3),*) jbond(nbond)
           bondtyp(nbond) = fredat(4)(1:5)

           ibond(nbond) = ibond(nbond) + natomindex
           jbond(nbond) = jbond(nbond) + natomindex

        end do

     else if (fredat(1) == '<ANGLE>') then
        do
           call rdfree( iuwtop, 20, fredat)
           if (fredat(1) == ' ') cycle DOREADW
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nangl = nangl + 1
           if (nangl > maxnangl) then
              write(6,*) 'Error : nangl exceeds maxnangl= ',maxnangl
              stop
           end if

           read(fredat(2),*) iangl(nangl)
           read(fredat(3),*) jangl(nangl)
           read(fredat(4),*) kangl(nangl)
           angltyp(nangl) = fredat(5)(1:8)

           iangl(nangl) = iangl(nangl) + natomindex
           jangl(nangl) = jangl(nangl) + natomindex
           kangl(nangl) = kangl(nangl) + natomindex

        end do

     end if

  END DO DOREADW

!-------- duplicate topology to all water --------

  nbondend = nbond - nbondindex
  nanglend = nangl - nanglindex

  do n=2,nwater

     do m=1, nbondend
        nbond = nbond + 1
        if (nbond > maxnbond) then
           write(6,*) 'Error : nbond exceeds maxnbond= ',maxnbond
           stop
        end if

        ibond(nbond) = ibond(nbondindex+m) + 3*(n-1)
        jbond(nbond) = jbond(nbondindex+m) + 3*(n-1)
        bondtyp(nbond) = bondtyp(nbondindex+m)
     end do

     do m=1, nanglend
        nangl = nangl + 1
        if (nangl > maxnangl) then
           write(6,*) 'Error : nangl exceeds maxnangl= ',maxnangl
           stop
        end if

        iangl(nangl) = iangl(nanglindex+m) + 3*(n-1)
        jangl(nangl) = jangl(nanglindex+m) + 3*(n-1)
        kangl(nangl) = kangl(nanglindex+m) + 3*(n-1)
        angltyp(nangl) = angltyp(nanglindex+m)
     end do

  end do

  close(iuwtop)

!     +     +     +     +     +     +     +

  return
end subroutine rdtop
