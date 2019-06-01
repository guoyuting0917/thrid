!*****************************
!*  prepmd.f90 Ver.3.9       *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-02-12 14:57:29 gota>

subroutine prepmd(npoly, npolytyp, npoly_mole, npoly_atom, &
     &            nwater, nmatom,nmatyp, nmatomtyp, &
     &            ifcenterfix_all, &
     &            ifcenterfix_poly, ifcenterfix_water, ifcenterfix_ma, &
     &            cenfix_free, &
     &            ifcenterfix_polytyp, &
     &            ifcenterfix_watertyp, &
     &            ifcenterfix_matyp, &
     &            degfree_poly, degfree_water, &
     &            degfree_ma, degfree_all, &
     &            iflocalfix, iflocalfixz, iflocalfixzg, &
     &            nlfix_deg_poly, &
     &            nlfix_deg_water, &
     &            nlfix_deg_ma, &
     &            nlfixz_deg_poly, nlfixz_deg_water, nlfixz_deg_ma, &
     &            nlfixzg_deg_poly, nlfixzg_deg_water, nlfixzg_deg_ma, &
     &            iflocalvel, &
     &            nlvel_deg_poly,nlvel_deg_water,nlvel_deg_ma, &
     &            dt_short_cal, dt_med_cal, dt_long_cal, &
     &            nstep_short, nstep_med, &
     &            xref, timeref, eps0, &
     &            rcut, ifbook, rcut_book, nstep_book, &
     &            rcutmor, ifbookmor, rcut_bookmor, nstep_bookmor, &
     &            rcutsh, ifbooksh, rcut_booksh, nstep_booksh, &
     &            rcutrfhfo, ifbookrfhfo, rcut_bookrfhfo, &
     &            nstep_bookrfhfo, &
     &            rcutrfhoo, ifbookrfhoo, rcut_bookrfhoo, &
     &            nstep_bookrfhoo, &
     &            rcutrfhoh, ifbookrfhoh, rcut_bookrfhoh, &
     &            nstep_bookrfhoh, &
     &            rcutdouo, ifbookdouo, rcut_bookdouo, &
     &            nstep_bookdouo, &
     &            rcutdouh, ifbookdouh, rcut_bookdouh, &
     &            nstep_bookdouh, &
     &            rcutrpvw, ifbookrpvw, rcut_bookrpvw, nstep_bookrpvw, &
     &            tcont_poly, tcont_water, tcont_ma, &
     &            xcel, ycel, zcel, &
     &            ifrattle, &
     &            mchain, tfreq, text, &
     &            vfreq, &
     &            pcont_axis, &
     &            nyosh, &
     &            solvetyp, solveindex, &
     &            netchrgsq)

  use md_global

  implicit none

! ARGUMENTS:
!     INPUT
  integer,intent(in):: npoly           ! all number of poly
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatom          ! number of monatomic molecules
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.


  logical,intent(in):: ifcenterfix_all   ! center fix for all
  logical,intent(in):: ifcenterfix_poly  ! center fix for polymer
  logical,intent(in):: ifcenterfix_water ! center fix for water
  logical,intent(in):: ifcenterfix_ma    ! center fix for monatomic mole.
  character(4),intent(in):: cenfix_free  ! COM not fixed in this direction

  logical,intent(in):: ifcenterfix_polytyp(:) ! center fix for each polymer
  logical,intent(in):: ifcenterfix_watertyp   ! center fix for each water
  logical,intent(in):: ifcenterfix_matyp(:)
                                          ! center fix for each monatomic mole.

  integer,intent(in):: nstep_med       ! number of step for medium force
  integer,intent(in):: nstep_short     ! number of step for short force

  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: timeref          ! time base value [sec]
  real(8),intent(in):: eps0             ! dielectric constant of vacuum [non-d]

  real(8),intent(in):: tcont_poly(:)    ! poly Temp. [non-d] in NVT
  real(8),intent(in):: tcont_water      ! H2O Temp. [non-d] in NVT
  real(8),intent(in):: tcont_ma(:)      ! monatomic mole. Temp. [non-d] in NVT

  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  real(8),intent(inout):: rcut             ! vdw cutoff length [non-d]
  logical,intent(in):: ifbook          ! flag for bookkeeping
  real(8),intent(inout):: rcut_book    ! cut off radius of bookkeeping [m->non-d]
  integer,intent(in):: nstep_book      ! bookkeeping interval

  real(8),intent(inout):: rcutmor          ! Morse cutoff length [non-d]
  logical,intent(in):: ifbookmor      ! flag for bookkeeping of Morse interaction
  real(8),intent(inout):: rcut_bookmor
                            ! cut off radius of bookkeeping[m] of Morse
  integer,intent(in):: nstep_bookmor  ! bookkeeping interval of Morse interaction

  real(8),intent(inout):: rcutsh           ! SH cutoff length [non-d]
  logical,intent(in):: ifbooksh        ! flag for bookkeeping of SH interaction
  real(8),intent(inout):: rcut_booksh   ! cut off radius of bookkeeping[m] of SH
  integer,intent(in):: nstep_booksh    ! bookkeeping interval of SH interaction

  real(8),intent(inout):: rcutrfhfo        ! RFH(FO) cutoff length [non-d]
  logical,intent(in):: ifbookrfhfo  ! flag for bookkeeping of RFH(FO) interaction
  real(8),intent(inout):: rcut_bookrfhfo
                         ! cut off radius of bookkeeping[m] of RFH(FO)
  integer,intent(in):: nstep_bookrfhfo
                         ! bookkeeping interval of RFH(FO) interaction

  real(8),intent(inout):: rcutrfhoo        ! RFH(OO) cutoff length [non-d]
  logical,intent(in):: ifbookrfhoo  ! flag for bookkeeping of RFH(OO) interaction
  real(8),intent(inout):: rcut_bookrfhoo
                         ! cut off radius of bookkeeping[m] of RFH(OO)
  integer,intent(in):: nstep_bookrfhoo
                         ! bookkeeping interval of RFH(OO) interaction

  real(8),intent(inout):: rcutrfhoh        ! RFH(OH) cutoff length [non-d]
  logical,intent(in):: ifbookrfhoh  ! flag for bookkeeping of RFH(OH) interaction
  real(8),intent(inout):: rcut_bookrfhoh
                         ! cut off radius of bookkeeping[m] of RFH(OH)
  integer,intent(in):: nstep_bookrfhoh
                         ! bookkeeping interval of RFH(OH) interaction

  real(8),intent(inout):: rcutdouo         ! DOU cutoff length [non-d] for O-Au
  logical,intent(in):: ifbookdouo
                         ! flag for bookkeeping of DOU interaction (O-Au)
  real(8),intent(inout):: rcut_bookdouo
                         ! cut off radius of bookkeep[non-d] of DOU (O-Au)
  integer,intent(in):: nstep_bookdouo
                         ! bookkeeping interval of DOU interaction (O-Au)

  real(8),intent(inout):: rcutdouh         ! DOU cutoff length [non-d] for H-Au
  logical,intent(in):: ifbookdouh
                         ! flag for bookkeeping of DOU interaction (H-Au)
  real(8),intent(inout):: rcut_bookdouh
                         ! cut off radius of bookkeep[non-d] of DOU (H-Au)
  integer,intent(in):: nstep_bookdouh
                         ! bookkeeping interval of DOU interaction (H-Au)

  real(8),intent(inout):: rcutrpvw         ! RP-VW cutoff length [non-d]
  logical,intent(in):: ifbookrpvw   ! flag for bookkeeping of RP-VW interaction
  real(8),intent(inout):: rcut_bookrpvw
                         ! cut off radius of bookkeeping[m] of RP-VW
  integer,intent(in):: nstep_bookrpvw
                         ! bookkeeping interval of RP-VW interaction

  logical,intent(in):: ifrattle        ! rattle flag

  integer,intent(in):: mchain          ! Nose-Hoover chain number
  real(8),intent(in):: tfreq            ! temperature frequency [non-d]
  real(8),intent(in):: text          ! external temp. [non-d] (Nose-Hoover chain)

  real(8),intent(in):: vfreq            ! volume change frequency [non-d]

  character(5),intent(in):: pcont_axis
                                 ! axis for pressure control (iso, aniso, etc.)

  integer,intent(in):: nyosh           ! expansion order of Yoshida-Suzuki method

  character(2),intent(in):: solvetyp ! solvent molecule

  logical,intent(in):: iflocalfix      ! fix atoms flag
  integer,intent(in):: nlfix_deg_poly(:) ! fixed degree of freedom of poly
  integer,intent(in):: nlfix_deg_water ! fixed degree of freedom of H2O
  integer,intent(in):: nlfix_deg_ma(:) ! fixed degree of freedom of matom

  logical,intent(in):: iflocalfixz     ! flag for fixing z coordinate of atoms
  integer,intent(in):: nlfixz_deg_poly(:) ! fixed degree of freedom of poly
  integer,intent(in):: nlfixz_deg_water   ! fixed degree of freedom of H2O
  integer,intent(in):: nlfixz_deg_ma(:)   ! fixed degree of freedom of matom

  logical,intent(in):: iflocalfixzg    ! flag for fixing z coordinate of COM of molecules
  integer,intent(in):: nlfixzg_deg_poly(:) ! fixed degree of freedom of poly
  integer,intent(in):: nlfixzg_deg_water   ! fixed degree of freedom of H2O
  integer,intent(in):: nlfixzg_deg_ma(:)   ! fixed degree of freedom of matom

  logical,intent(in):: iflocalvel          ! flag to force local atomic velocity

  integer,intent(in):: nlvel_deg_poly(:)   ! fixed degree of freedom of poly
  integer,intent(in):: nlvel_deg_water     ! fixed degree of freedom of H2O
  integer,intent(in):: nlvel_deg_ma(:)     ! fixed degree of freedom of matom

!     INPUT&OUTPUT
  real(8),intent(inout):: dt_long_cal      ! time step of long force [sec->non-d]

!     OUTPUT
  integer,intent(out):: degfree_poly(:) ! degree of freedom of polymer1
  integer,intent(out):: degfree_water   ! degree of freedom of H2O
  integer,intent(out):: degfree_ma(:)   ! degree of freedom of monatomic mole.
  integer,intent(out):: degfree_all     ! degree of freedom of all molecules

  real(8),intent(out):: dt_med_cal       ! time step of medium force [non-d]
  real(8),intent(out):: dt_short_cal     ! time step of short force [non-d]

!      real*8:: vdw_welij_solve  ! well depth of vdw parameter of solvent
!      real*8:: vdw_radij_solve  ! radius of vdw parameter of solvent
  integer,intent(out):: solveindex      ! atmindex of solvent atom

  real(8),intent(out):: netchrgsq           ! = (sum(qi))**2

! LOCAL
  integer:: i               ! do loop index
  real(8):: pi               ! = 3.14159...
  real(8):: tempmax          ! maximum temp. among control temp.
  real(8):: hcelmin          ! minimum half-cell length

!     +     +     +     +     +     +

!---- some preparation
  pi = acos(-1.0d0)

!-------- Convert equilibrium angle from degree to radian --------
  do i=1,nangltyp
     para_eqangl(i) = para_eqangl(i) * pi/180.0d0
  end do

  do i=1,nanglubtyp
     para_eqanglub(i) = para_eqanglub(i) * pi/180.0d0
  end do

  do i=1,ntorstyp
     para_phsang(i) = para_phsang(i) * pi/180.0d0
  end do

  do i=1,ntorsimtyp
     para_phsangim(i) = para_phsangim(i) * pi/180.0d0
  end do

!-------- Convert atomic charge * e/sqrt(4*pi*eps0) --------
  do i=1,natom
     atmchrg(i) = atmchrg(i) * 1.0d0/dsqrt(4.0d0*pi*eps0)
  end do

!-------- calculate (sum(qi))**2 --------
  netchrgsq = 0.0d0
  do i=1,natom
     netchrgsq = netchrgsq + atmchrg(i)
  end do
  netchrgsq = netchrgsq * netchrgsq

!-------- Calculate degree of freedom --------
!     flex molecule

  degfree_all = 0

  do i = 1, npolytyp
     degfree_poly(i) = npoly_mole(i)*npoly_atom(i) * 3
     if (iflocalfix) then
        degfree_poly(i) = degfree_poly(i) - nlfix_deg_poly(i)
     end if
     if (iflocalfixz) then
        degfree_poly(i) = degfree_poly(i) - nlfixz_deg_poly(i)
     end if
     if (iflocalfixzg) then
        degfree_poly(i) = degfree_poly(i) - nlfixzg_deg_poly(i)
     end if
     if (iflocalvel) then
        degfree_poly(i) = degfree_poly(i) - nlvel_deg_poly(i)
     end if
     degfree_all = degfree_all + degfree_poly(i)
  end do

  degfree_water = nwater*3 * 3
  if (iflocalfix) then
     degfree_water = degfree_water - nlfix_deg_water
  end if
  if (iflocalfixz) then
     degfree_water = degfree_water - nlfixz_deg_water
  end if
  if (iflocalfixzg) then
     degfree_water = degfree_water - nlfixzg_deg_water
  end if
  if (iflocalvel) then
     degfree_water = degfree_water - nlvel_deg_water
  end if
  degfree_all = degfree_all + degfree_water

  do i = 1, nmatyp
     degfree_ma(i) = nmatomtyp(i) * 3
     if (iflocalfix) then
        degfree_ma(i) = degfree_ma(i) - nlfix_deg_ma(i)
     end if
     if (iflocalfixz) then
        degfree_ma(i) = degfree_ma(i) - nlfixz_deg_ma(i)
     end if
     if (iflocalfixzg) then
        degfree_ma(i) = degfree_ma(i) - nlfixzg_deg_ma(i)
     end if
     if (iflocalvel) then
        degfree_ma(i) = degfree_ma(i) - nlvel_deg_ma(i)
     end if
     degfree_all = degfree_all + degfree_ma(i)
  end do

!      degfree_all = degfree_poly + degfree_water + degfree_ma

  if (ifcenterfix_all) then
     degfree_all = degfree_all - 3

     if (cenfix_free /= 'none') then
         degfree_all = degfree_all + 1
     end if
  else
     if (ifcenterfix_poly) then
        do i = 1, npolytyp
           if (.not. ifcenterfix_polytyp(i)) cycle
           degfree_all = degfree_all - 3
           degfree_poly(i) = degfree_poly(i) - 3
        end do
     endif

     if (ifcenterfix_water .and. ifcenterfix_watertyp) then
        degfree_all = degfree_all - 3
        degfree_water = degfree_water - 3
     endif

     if (ifcenterfix_ma) then
        do i = 1, nmatyp
           if (.not. ifcenterfix_matyp(i)) cycle
           degfree_all = degfree_all - 3
           degfree_ma(i) = degfree_ma(i) - 3
        end do
     endif
  end if

  if (degfree_all < 0) degfree_all = 0
!     rigid molecule (for rattle)
  if (ifrattle) then
     degfree_all = degfree_all - nconst
     if (degfree_all < 0) degfree_all = 0
!        !!! Caution! This routine is only for RATTLE! !!!
     do i = 1, npolytyp
        degfree_poly(i) = degfree_poly(i) - nconst_poly(i)
        if (degfree_poly(i) < 0) degfree_poly(i) = 0
     end do
     degfree_water = degfree_water - nconst_water
     if (degfree_water < 0) degfree_water = 0
  end if

!-------- Calculate time step --------
  dt_med_cal = dt_long_cal / dble(nstep_med)
  dt_short_cal = dt_med_cal / dble(nstep_short)

!---- non-dimensonalize
  dt_long_cal = dt_long_cal / timeref
  dt_med_cal = dt_med_cal / timeref
  dt_short_cal = dt_short_cal / timeref

!---- Calculate radius of bookkeeping

  tempmax = max(tcont_poly(1),tcont_water,tcont_ma(1))
!      rcut_book =  rcut
!     &           + dble(nstep_book)*dsqrt(3.0d0*tempmax)*dt_long_cal
!      rcut_book = 0.49999d0 * dmin1(xcel,ycel,zcel)
!      rcut_book = 1.5d-9/xref
  rcut_book = rcut_book/xref

  hcelmin = min(xcel*0.5d0,ycel*0.5d0,zcel*0.5d0)

  if (abs(rcut) < 1.0d-16) then
     rcut = hcelmin
  end if

  if (rcut > hcelmin) then
     write(6,*) 'Warning: calculate incorrect rcut= ',rcut
     stop
  end if

  if ((ifbook) .and. (rcut_book > hcelmin)) then
     write(6,*) 'Warning: calculate incorrect rcut_book= ',rcut_book
     stop
  end if


  rcut_bookmor = rcut_bookmor / xref

  if (abs(rcutmor) < 1.0d-16) then
     rcutmor = hcelmin
  end if

  if (rcutmor > hcelmin) then
     write(6,*) 'Warning: calculate incorrect rcutmor= ',rcutmor
     stop
  end if

  if ((ifbookmor) .and. (rcut_bookmor > hcelmin)) then
     write(6,*) 'Warning: calculate incorrect rcut_bookmor= ', rcut_bookmor
     stop
  end if

  rcut_booksh = rcut_booksh / xref

  if (abs(rcutsh) < 1.0d-16) then
     rcutsh = hcelmin
  end if

  if (rcutsh > hcelmin) then
     write(6,*) 'Warning: calculate incorrect rcutsh= ',rcutsh
     stop
  end if

  if ((ifbooksh) .and. (rcut_booksh > hcelmin)) then
       write(6,*) 'Warning: calculate incorrect rcut_booksh= ', rcut_booksh
     stop
  end if

  rcut_bookrfhfo = rcut_bookrfhfo / xref

  if (abs(rcutrfhfo) < 1.0d-16) then
     rcutrfhfo = hcelmin
  end if

  if (rcutrfhfo > hcelmin) then
     write(6,*) 'Warning: calculate incorrect rcutrfhfo= ',rcutrfhfo
     stop
  end if

  if ((ifbookrfhfo) .and. (rcut_bookrfhfo > hcelmin)) then
     write(6,*) 'Warning: calculate incorrect rcut_bookrfhfo= ', rcut_bookrfhfo
     stop
  end if

  rcut_bookrfhoo = rcut_bookrfhoo / xref

  if (abs(rcutrfhoo) < 1.0d-16) then
     rcutrfhoo = hcelmin
  end if

  if (rcutrfhoo > hcelmin) then
     write(6,*) 'Warning: calculate incorrect rcutrfhoo= ',rcutrfhoo
     stop
  end if

  if ((ifbookrfhoo) .and. (rcut_bookrfhoo > hcelmin)) then
     write(6,*) 'Warning: calculate incorrect rcut_bookrfhoo= ', rcut_bookrfhoo
     stop
  end if

  rcut_bookrfhoh = rcut_bookrfhoh / xref

  if (abs(rcutrfhoh) < 1.0d-16) then
     rcutrfhoh = hcelmin
  end if

  if (rcutrfhoh > hcelmin) then
     write(6,*) 'Warning: calculate incorrect rcutrfhoh= ',rcutrfhoh
     stop
  end if

  if ((ifbookrfhoh) .and. (rcut_bookrfhoh > hcelmin)) then
     write(6,*) 'Warning: calculate incorrect rcut_bookrfhoh= ', rcut_bookrfhoh
     stop
  end if

  rcut_bookdouo = rcut_bookdouo / xref

  if (abs(rcutdouo) < 1.0d-16) then
     rcutdouo = hcelmin
  end if

  if (rcutdouo > hcelmin) then
     write(6,*) 'Warning: calculate incorrect rcutdouo= ', rcutdouo
     stop
  end if

  if ((ifbookdouo) .and. (rcut_bookdouo > hcelmin)) then
     write(6,*) 'Warning: calculate incorrect rcut_bookdouo= ', rcut_bookdouo
     stop
  end if

  rcut_bookdouh = rcut_bookdouh / xref

  if (abs(rcutdouh) < 1.0d-16) then
     rcutdouh = hcelmin
  end if

  if (rcutdouh > hcelmin) then
     write(6,*) 'Warning: calculate incorrect rcutdouh= ',rcutdouh
     stop
  end if

  if ((ifbookdouh) .and. (rcut_bookdouh > hcelmin)) then
     write(6,*) 'Warning: calculate incorrect rcut_bookdouh= ', rcut_bookdouh
     stop
  end if

  rcut_bookrpvw = rcut_bookrpvw / xref

  if (abs(rcutrpvw) < 1.0d-16) then
     rcutrpvw = hcelmin
  end if

  if (rcutrpvw > hcelmin) then
     write(6,*) 'Warning: calculate incorrect rcutrpvw= ',rcutrpvw
     stop
  end if

  if ((ifbookrpvw) .and. (rcut_bookrpvw > hcelmin)) then
     write(6,*) 'Warning: calculate incorrect rcut_bookrpvw= ', rcut_bookrpvw
     stop
  end if

!---- Calculate masses of Nose-Hoover chain thermostat

  qmass(1) = dble(degfree_all)*text/tfreq/tfreq ! = NfkT/wp**2
  do i=2,mchain
     qmass(i) = text/tfreq/tfreq          ! = kT/wp**2
  end do

!---- Calculate Andersen (Hoover type) barostat

  if (pcont_axis == 'iso') then
      vmass = dble(degfree_all+3)*text/(vfreq*vfreq) ! = (Nf+d)kT/wb**2
  else
      vmass = dble(degfree_all+3)*text/(3.0d0*vfreq*vfreq) ! = (Nf+d)kT/d*wb**2
  end if
  xlogv = log(xcel*ycel*zcel) / 3.0d0     ! = ln(Lx*Ly*Lz)/d
  xboxh(1) = xcel   ! = Lx
  xboxh(2) = ycel   ! = Ly
  xboxh(3) = zcel   ! = Lz

!---- Input coefficient of Yoshida-Suzuki method

  if (nyosh == 1) then
     yscoeff(1) = 1.0d0
  end if
  if (nyosh == 3) then
     yscoeff(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0)) ! = 1/(2-2**(1/3))
     yscoeff(2) = 1.0d0 - 2.0d0*yscoeff(1)             ! = 1-2*w(1)
     yscoeff(3) = yscoeff(1)                           ! = 1/(2-2**(1/3))
  end if
  if (nyosh == 5) then
     yscoeff(1) = 1.0d0 / (4.0d0-4.0d0**(1.0d0/3.0d0)) ! = 1/(4-4**(1/3))
     yscoeff(2) = yscoeff(1)                           ! = 1/(4-4**(1/3))
     yscoeff(3) = 1.0d0 - 4.0d0*yscoeff(1)             ! = 1-4*w(1)
     yscoeff(4) = yscoeff(1)                           ! = 1/(4-4**(1/3))
     yscoeff(5) = yscoeff(1)                           ! = 1/(4-4**(1/3))
  end if

!---- Input coefficient of Long-range correction

  do i = 1, natmtyp

     if (solvetyp(1:2) == para_atmtyp(i)) then
        solveindex = i
     end if

  end do

!     +     +     +     +     +     +

end subroutine prepmd
