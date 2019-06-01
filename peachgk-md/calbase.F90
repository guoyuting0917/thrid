!**************************************
!*  calbase.f90 Ver.2.7               *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
! Time-stamp: <2015-01-06 14:26:44 gota>

subroutine calbase(xref,eref,mref,qref, &
     &             vref,timeref,tempref,pref,fref, &
     &             eps0ref, &
     &             npolytyp,nmatyp, &
     &             eps0, &
     &             alpha,rrcut, &
     &             xcel,ycel,zcel, &
     &             rcut, &
     &             tcont_poly,tcont_water,tcont_ma, &
     &             tcont_poly_ini,tcont_water_ini,tcont_ma_ini, &
     &             text,tfreq, &
     &             vfreq,pext, &
     &             rcutmor, &
     &             rcutsh, &
     &             rcutrfhfo,rcutrfhoo,rcutrfhoh, &
     &             rcutdouo,rcutindouo,rcutdouh, &
     &             rcutrpvw, &
     &             d_rini,d_rmax,d_econv,d_rmsf, &
     &             limitdist)

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: eref             ! energy base value [J]
  real(8),intent(in):: mref             ! mass base value [kg]
  real(8),intent(in):: qref             ! charge base value [C]

  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: nmatyp          ! number of species of monatomic mole.

!     OUTPUT
  real(8),intent(out):: vref             ! velocity base value [m/s]
  real(8),intent(out):: timeref          ! time base value [sec]
  real(8),intent(out):: tempref          ! temperature base value [K]
  real(8),intent(out):: pref             ! pressure base value [Pa]
  real(8),intent(out):: fref             ! force base value [N]

  real(8),intent(out):: eps0ref         ! dielectric constant base value [c^2/Jm]

!     INPUT&OUTPUT
  real(8),intent(inout):: alpha            ! parameter alpha [1/m->non-d]
  real(8),intent(inout):: rrcut
                            ! ewald real space cutoff length [m->non-d]
  real(8),intent(inout):: eps0             ! parameter eps0 [c^2/Jm->non-d]
  real(8),intent(inout):: xcel             ! x cell length [m->non-d]
  real(8),intent(inout):: ycel             ! y cell length [m->non-d]
  real(8),intent(inout):: zcel             ! z cell length [m->non-d]
  real(8),intent(inout):: rcut             ! vdw cutoff length [m->non-d]
  real(8),intent(inout):: tcont_poly(:)    ! poly Temp. [K->non-d]
  real(8),intent(inout):: tcont_water      ! H2O Temp. [K->non-d]
  real(8),intent(inout):: tcont_ma(:)      ! monatomic mole. Temp. [K->non-d]
  real(8),intent(inout):: tcont_poly_ini   ! poly Temp. [K->non-d]
  real(8),intent(inout):: tcont_water_ini  ! H2O Temp. [K->non-d]
  real(8),intent(inout):: tcont_ma_ini     ! MA Temp. [K->non-d]
  real(8),intent(inout):: text             ! external temperature [K->non-d]

  real(8),intent(inout):: tfreq            ! temperature frequency [1/s->non-d]

  real(8),intent(inout):: vfreq            ! volume change frequency [1/s->non-d]
  real(8),intent(inout):: pext             ! external pressure [Pa->non-d]      

  real(8),intent(inout):: rcutmor          ! Morse cutoff length [m->non-d]

  real(8),intent(inout):: rcutsh           ! SH cutoff length [m->non-d]

  real(8),intent(inout):: rcutrfhfo        ! RFH(FO) cutoff length [m->non-d]
  real(8),intent(inout):: rcutrfhoo        ! RFH(OO) cutoff length [m->non-d]
  real(8),intent(inout):: rcutrfhoh        ! RFH(OH) cutoff length [m->non-d]

  real(8),intent(inout):: rcutdouo  
                            ! DOU cutoff length [m->non-d] for O-Au
  real(8),intent(inout):: rcutindouo       ! DOU cutin length [m->non-d] for O-Au
  real(8),intent(inout):: rcutdouh
                            ! DOU cutoff length [m->non-d] for H-Au

  real(8),intent(inout):: rcutrpvw ! cutoff length for RP-VW interaction

  real(8),intent(inout):: d_rini   ! initial displacement dr for EM [m->non-d]
  real(8),intent(inout):: d_rmax   ! maximum displacement dr for EM [m->non-d]

  real(8),intent(inout):: d_econv  ! convergence condition for energy
                                   !    in EM [J->non-d]
  real(8),intent(inout):: d_rmsf   ! convergence condition for 
                                   !    root mean square force in EM [N->non-d]

  real(8),intent(inout):: limitdist   ! maximum atomic displacement [m->non-d]

! LOCAL:
  real(8):: kb = 1.38066d-23 ! Boltzmann constant [J/K]
  integer:: i

!     +     +     +     +     +     +

!---- calculate base value of physical quantity
  vref = dsqrt(eref/mref)
  timeref = xref / vref
  tempref = eref / kb
  pref = eref / xref**3
  fref = eref / xref

!---- calculate base value of coefficient
  eps0ref = qref**2 / eref / xref

!---- non-dimensionalize coefficients & variables
  alpha = alpha * xref
  rrcut = rrcut / xref
  eps0 = eps0 / eps0ref
  xcel = xcel / xref
  ycel = ycel / xref
  zcel = zcel / xref
  rcut = rcut / xref
  tcont_poly(1) = tcont_poly(1) / tempref
  do i = 2, npolytyp
     tcont_poly(i) = tcont_poly(i) / tempref
  end do
  tcont_water = tcont_water / tempref
  tcont_ma(1) = tcont_ma(1) / tempref
  do i = 2, nmatyp
     tcont_ma(i) = tcont_ma(i) / tempref
  end do
  tcont_poly_ini = tcont_poly_ini / tempref
  tcont_water_ini = tcont_water_ini / tempref
  tcont_ma_ini = tcont_ma_ini / tempref
  text = text /tempref
  tfreq = tfreq * timeref
  vfreq = vfreq * timeref
  pext = pext / pref

  rcutmor = rcutmor / xref

  rcutsh = rcutsh / xref

  rcutrfhfo = rcutrfhfo / xref
  rcutrfhoo = rcutrfhoo / xref
  rcutrfhoh = rcutrfhoh / xref

  rcutdouo = rcutdouo / xref
  rcutindouo = rcutindouo / xref
  rcutdouh = rcutdouh / xref

  rcutrpvw = rcutrpvw / xref

  d_rini = d_rini / xref
  d_rmax = d_rmax / xref
  d_econv = d_econv / eref
  d_rmsf = d_rmsf / fref

  limitdist = limitdist / xref

!     +     +     +     +     +     +

end subroutine calbase
