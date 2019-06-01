!*******************************
!*  rdpara_cstmnb.f90 Ver.1.1  *
!*      for peachgk_md.f       *
!*            by G.Kikugawa    *
!*     Coded by T.Nakano       *
!*******************************
! Time-stamp: <Aug 31 2017>

subroutine rdpara_cstmnb(iuparacstmnb, &
     &                   ifcellindex_cstmnb,ifbookcstmnb, &
     &                   xref,eref,mref,qref, &
     &                   vref,timeref,tempref,pref,fref,eps0ref, &
     &                   xcel,ycel,zcel)

  use interface_tools

  use md_global
  use mpi_global
  use cstmnb

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: iuparacstmnb    ! input custom NB parameter file unit

  logical,intent(inout):: ifcellindex_cstmnb ! flag for cell index (custom NB)
  logical,intent(in):: ifbookcstmnb       
                               ! flag for bookkeeping of custom NB interaction

  real(8),intent(in):: xref            ! distanse base value [m]
  real(8),intent(in):: eref            ! energy base value [J]
  real(8),intent(in):: mref            ! mass base value [kg]
  real(8),intent(in):: qref            ! charge base value [C]
  real(8),intent(in):: vref            ! velocity base value [m/s]
  real(8),intent(in):: timeref         ! time base value [sec]
  real(8),intent(in):: tempref         ! temperature base value [K]
  real(8),intent(in):: pref            ! pressure base value [Pa]
  real(8),intent(in):: fref            ! force base value [N]

  real(8),intent(in):: eps0ref         ! dielectric constant base value [c^2/Jm]

  real(8),intent(in):: xcel            ! x cell length [non-d]
  real(8),intent(in):: ycel            ! y cell length [non-d]
  real(8),intent(in):: zcel            ! z cell length [non-d]

! LOCAL:
  character(80):: fredat(maxnword)

  real(8):: an=6.0221367d23  ! Avogadro's number
  real(8):: pi = dacos(-1.0d0)
  integer:: i,j,k,l          ! do loop index
  integer:: nword

  logical:: ifcstmnb_i, ifcstmnb3_i
  logical:: ifcstmnb_j, ifcstmnb3_j
  logical:: ifcstmnb_k, ifcstmnb3_k
  integer:: cstmnb_i, cstmnb_j, cstmnb_k
  integer:: itype, jtype

  integer:: nline

  integer:: cellcount(3)

  real(8):: hcelmin
  real(8):: rcut_maxSW, rcut_book

  integer:: ncstmnbtyp2, ncstmnbtyp3
  character(5):: para_cstmnbtyp2(maxncstmnbtyp)
  character(8):: para_cstmnbtyp3(maxncstmnbtyp)
  character(8):: char_cstmnbtyp

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!-------- Read parameter of custom NB interaction --------

  ncstmnbtyp = 0
  READCSTMNB:DO
     call rdfree(iuparacstmnb, maxnword, fredat)

     if (fredat(1) == '<END>') then
        exit
     else if (fredat(1) == ' ') then
        cycle READCSTMNB

     else if (fredat(1) == '<CUSTOM_NB>') then
!        nline = 0
        ncstmnbtyp = 0

        do
           call rdfree_w(iuparacstmnb, maxnword, fredat, nword)
           if (fredat(1) == ' ') cycle READCSTMNB
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

!!! Describe your column definition here
!!!  (variables should be declared in cstmnb.F90)
!           nline = nline + 1
           ncstmnbtyp = ncstmnbtyp + 1
           read(fredat(2),*) para_cstmnbtyp(ncstmnbtyp)
        end do

     else if (fredat(1) == '<COMMON_2&3>') then
        do
           call rdfree_w(iuparacstmnb, maxnword, fredat, nword)
           if (fredat(1) == ' ') cycle READCSTMNB
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle

           read(fredat(2),*) rcut_bookSW
           read(fredat(3),*) nstep_bookSW
        end do

     else if (fredat(1) == '<2-BODY>') then

        ncstmnbtyp2 = 0
        do
           call rdfree_w(iuparacstmnb, maxnword, fredat, nword)
           if (fredat(1) == ' ') cycle READCSTMNB
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle

           ncstmnbtyp2 = ncstmnbtyp2 + 1
           para_cstmnbtyp2(ncstmnbtyp2) = fredat(2)(1:5)

           read(fredat(3),*) rcut_SW(ncstmnbtyp2)
           read(fredat(4),*) A_SW(ncstmnbtyp2)
           read(fredat(5),*) B_SW(ncstmnbtyp2)
           read(fredat(6),*) p_SW(ncstmnbtyp2)
           read(fredat(7),*) q_SW(ncstmnbtyp2)
           read(fredat(8),*) epsilon_SW(ncstmnbtyp2)
           read(fredat(9),*) sigma2_SW(ncstmnbtyp2)
        end do

     else if (fredat(1) == '<3-BODY>') then

        ncstmnbtyp3 = 0

        do
           call rdfree_w(iuparacstmnb, maxnword, fredat, nword)
           if (fredat(1) == ' ') cycle READCSTMNB
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle

           ncstmnbtyp3 = ncstmnbtyp3 + 1
           para_cstmnbtyp3(ncstmnbtyp3) = fredat(2)(1:8)

           read(fredat(3),*) lambda_SW(ncstmnbtyp3)
           read(fredat(4),*) gamma_SW(ncstmnbtyp3)
           read(fredat(5),*) theta0_SW(ncstmnbtyp3)
           read(fredat(6),*) sigma3_SW(ncstmnbtyp3)
        end do

     end if

  END DO READCSTMNB

!---- non-dimensionalize
  epsilon_SW(1:ncstmnbtyp2) = epsilon_SW(1:ncstmnbtyp2) / eref
  sigma2_SW(1:ncstmnbtyp2) = sigma2_SW(1:ncstmnbtyp2) / xref
  A_SW(1:ncstmnbtyp2) = A_SW(1:ncstmnbtyp2) * epsilon_SW(1:ncstmnbtyp2)

  theta0_SW(1:ncstmnbtyp3) = theta0_SW(1:ncstmnbtyp3) * pi / 180.d0
  sigma3_SW(1:ncstmnbtyp3) = sigma3_SW(1:ncstmnbtyp3) / xref

  if (ncstmnbtyp3 /= ncstmnbtyp2) then
     write(6,*)'Error: mismatch between num. of lambda and epsilon'
     write(6,*)'       in 3-body and 2-body, respectively'
     write(6,*)'       Add epsilon term in <3-BODY>'
     stop
  end if
  lambda_SW(1:ncstmnbtyp3) = lambda_SW(1:ncstmnbtyp3) &
       &                   * epsilon_SW(1:ncstmnbtyp3)
!  gamma_SW(1:ncstmnbtyp3) = gamma_SW(1:ncstmnbtyp3) * sigma3_SW(1:ncstmnbtyp3)

!---- cell index checking
  if (.not.ifbookcstmnb .and. ifcellindex_cstmnb) then
     write(6,*) 'Error: if you want to use cellindex method,'
     write(6,*) '       enable ifbookcstmnb'
     stop
  end if

  if (ifcellindex_cstmnb) then

     cellcount(1) = INT(xcel/rcut_bookSW)
     cellcount(2) = INT(ycel/rcut_bookSW)
     cellcount(3) = INT(zcel/rcut_bookSW)

     if (cellcount(1) < 3 .or. cellcount(2) < 3 .or. &
          & cellcount(3) < 3 .or. &
          & cellcount(1)*cellcount(2)*cellcount(3) <= 27) then
        if (irank == 0) then
           write(6,*) 'Warning: cell length is too short to execute'
           write(6,*) '         cell-index method (SW)'
           write(6,*) '         force to use normal book-keeping'
        end if
        ifcellindex_cstmnb = .false.
     end if

  endif

!---- checking cutoff parameters

  hcelmin = min(xcel*0.5d0,ycel*0.5d0,zcel*0.5d0)

  rcut_maxSW = 0.0d0
  do i = 1, ncstmnbtyp2
    rcut_maxSW = max(rcut_maxSW,rcut_SW(i))
  end do

  if (abs(rcut_maxSW) < 1.0d-16) then
     rcut_maxSW = hcelmin
  end if

  if (rcut_maxSW > hcelmin) then
     write(6,*) 'Error: find incorrect rcut_SW = ',rcut_maxSW
     stop
  end if

  if ((ifbookcstmnb) .and. (rcut_bookSW > hcelmin)) then
     write(6,*) 'Error: find incorrect rcut_bookSW= ', &
          &     rcut_bookSW
     stop
  end if

!---- list of atoms of custom NB interaction
!** usually not needed to change here
  ncstmnb = 0

  do i = 1, natom
     ifcstmnb(i) = .false.
     atmindex_ncstmnb(i) = 0
     do j = 1, ncstmnbtyp

        if (atmtyp(i)(1:2) == para_cstmnbtyp(j)(1:2)) then
           ncstmnb = ncstmnb + 1
           ncstmnblist(ncstmnb) = i
           atmindex_ncstmnb(i) = j
           ifcstmnb(i) = .true.
           exit
        end if

     end do
  end do

!---- order and link custom NB parameter
  ifcalSW2(1:ncstmnbtyp,1:ncstmnbtyp)              = .false.
  ifcalSW3(1:ncstmnbtyp,1:ncstmnbtyp,1:ncstmnbtyp) = .false.

  do i = 1, ncstmnbtyp
     do j = 1, ncstmnbtyp


        ifcstmnb_i = .false.
        ifcstmnb_j = .false.

        do l = 1, ncstmnbtyp2

           !-- scan i, j
           if (para_cstmnbtyp(i)(1:2) == para_cstmnbtyp2(l)(1:2)) then
              ifcstmnb_i = .true.
           end if
           if (para_cstmnbtyp(j)(1:2) == para_cstmnbtyp2(l)(4:5)) then
              ifcstmnb_j = .true.
           end if

           if (ifcstmnb_i .and. ifcstmnb_j) then
              char_cstmnbtyp(1:5) = &
                   &   para_cstmnbtyp(i)(1:2) // '-' // para_cstmnbtyp(j)(1:2)

              if (para_cstmnbtyp2(l)(1:5) == char_cstmnbtyp(1:5)) then
                 ifcalSW2(i,j) = .true.   ! this pair calculated in cstmnb
                 cstmnbtypeindex2(i,j) = l
                 exit
              endif
           endif

        end do

        do k = 1, ncstmnbtyp

           ifcstmnb3_i = .false.
           ifcstmnb3_j = .false.
           ifcstmnb3_k = .false.

           do l = 1, ncstmnbtyp3

              !-- scan i, j, k
              if (para_cstmnbtyp(i)(1:2) == para_cstmnbtyp3(l)(1:2)) then
                 ifcstmnb3_i = .true.
              end if
              if (para_cstmnbtyp(j)(1:2) == para_cstmnbtyp3(l)(4:5)) then
                 ifcstmnb3_j = .true.
              end if
              if (para_cstmnbtyp(k)(1:2) == para_cstmnbtyp3(l)(7:8)) then
                 ifcstmnb3_k = .true.
              end if

              if (ifcstmnb3_i .and. ifcstmnb3_j .and. ifcstmnb3_k) then
                 char_cstmnbtyp(1:8) = para_cstmnbtyp(i)(1:2) &
                      &              // '-' // para_cstmnbtyp(j)(1:2) &
                      &              // '-' // para_cstmnbtyp(k)(1:2)
                 if (para_cstmnbtyp3(l) == char_cstmnbtyp) then
                    ifcalSW3(i,j,k) = .true.
                    cstmnbtypeindex3(i,j,k) = l
                    exit
                 endif
              endif

           end do
        end do

     end do
  end do

!---- register intermolecular interaction type for particle pair
!** usually not needed to change here

!     - detect custom NB pair
  do i = 1, natmtyp

     itype = 0
     do l = 1, ncstmnbtyp
        if (para_atmtyp(i)(1:2) == para_cstmnbtyp(l)(1:2)) then
           itype = l
           exit
        endif
     end do

     do j = i, natmtyp

        jtype = 0
        do l = 1, ncstmnbtyp
           if (para_atmtyp(j)(1:2) == para_cstmnbtyp(l)(1:2)) then
              jtype = l
              exit
           endif
        end do

        if (ifcalSW2(itype,jtype)) then
           inter_inttyp(i,j) = INTTYPE_CSTMNB
           inter_inttyp(j,i) = INTTYPE_CSTMNB

           vdw_welij(i,j) = 0.0d0
           vdw_radij(i,j) = 0.0d0
           vdw_welij(j,i) = 0.0d0
           vdw_radij(j,i) = 0.0d0
        end if

     end do
  end do

#if defined(_INTTYP_DEBUG)
  do i = 1, natmtyp
     do j = i, natmtyp
        write(6,*) i,j,inter_inttyp(i,j)
     end do
  end do
#endif

  close(iuparacstmnb)

!---- dynamic memory allocation for neighbor lists ----

!  maxcstmnblist2 = maxnatom*maxcstmnblist / nproc
  maxcstmnblist2 = ncstmnb*maxcstmnblist / nproc
  maxcstmnblist3 = maxcstmnblist2 * maxcstmnblist
                       ! memory allocation size for cstmnblistall for 1 process

  allocate(cstmnb_listall2(maxcstmnblist2))
  allocate(cstmnb_indexall3(maxcstmnblist2))
  allocate(cstmnb_listall3(maxcstmnblist3))

  write(6,*) 'Dynamically allocate memory for cstmnb lists at process: ',irank

!     +     +     +     +     +     +     +

end subroutine rdpara_cstmnb
