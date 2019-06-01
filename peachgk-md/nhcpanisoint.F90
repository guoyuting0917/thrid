!*****************************
!*  nhcpisoint.f90 Ver.1.0   *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>

subroutine nhcpanisoint(degfree_all, &
     &                  dt_long_cal, &
     &                  mchain,text_c, &
     &                  next,nyosh, &
     &                  pintt,pext_c,pcont_dir)

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: degfree_all     ! degree of freedom of all atoms

  real(8),intent(in):: dt_long_cal      ! time step of long force [non-d]

  integer,intent(in):: mchain          ! Nose-Hoover chain number (>1 in NVT)
  real(8),intent(in):: text_c           ! external temp. [K] (Nose-Hoover chain)

  integer,intent(in):: next            ! iteration number of extended system
  integer,intent(in):: nyosh           ! expansion order of Yoshida-Suzuki method

  real(8),intent(in):: pintt(:,:)       ! internal pressure tensor in NPT (MTK)
  real(8),intent(in):: pext_c           ! external pressure in NPT(MTK)
  integer,intent(in):: pcont_dir        ! direction for uniaxial NPT control

! LOCAL:
  real(8):: kint             ! kinetic energy of physical space (=mv**2)

  real(8):: wdti2(maxnyosh)  ! yscoeff()*dt/next/2
  real(8):: wdti4(maxnyosh)  ! yscoeff()*dt/next/4
  real(8):: wdti8(maxnyosh)  ! yscoeff()*dt/next/8

  real(8):: glogs(maxnhc)    ! Glogs() in thermostat velocity update
  real(8):: gboxv(3)         ! Gboxv for barostat velocity update

  real(8):: gndfkt           ! (Nf+df)kT (df=1 for uniaxial, 3 for aniso)
  real(8):: gkt              ! kT

  integer:: i               ! do loop index
  integer:: inc             ! do loop index of iteration number (next)
  integer:: iys             ! do loop index of Yoshida-Suzuki method (nyosh)

  real(8):: expwdti8         ! = exp(wdti8*vxi(m))

  real(8):: dnf              ! = 1/Nf

  real(8):: vol              ! = Lx*Ly*Lz

  real(8):: pintv(3)         ! = (Pxx, Pyy, Pzz)

!     +     +     +     +     +     +     +     +

! --- some preparation ---

! substitute wdtis
  wdti2(1:nyosh) = yscoeff(1:nyosh) * dt_long_cal / dble(next) / 2.0d0
  wdti4(1:nyosh) = wdti2(1:nyosh) / 2.0d0
  wdti8(1:nyosh) = wdti4(1:nyosh) / 2.0d0

! calculate (Nf+df)kT and kT
  if (pcont_dir == 0) then   ! full anisotropic NPT
      gndfkt = text_c * dble(degfree_all+3)
  else   ! uniaxial NPT
      gndfkt = text_c * dble(degfree_all+1)
  end if
  gkt  = text_c

! calculate dnf
  dnf = 1.0d0/dble(degfree_all)

! calculate cell volume
  vol = xboxh(1) * xboxh(2) * xboxh(3)

! pressure diagonal component
  pintv(1) = pintt(1,1)
  pintv(2) = pintt(2,2)
  pintv(3) = pintt(3,3)

!     ----- Main loop of Nose-Hoover chain - Andersen (Hoover type) -----

!     --- calculate kinetic energy of physical system ---
  kint = 0.0d0
  do i=1,natom
     kint = kint + atmmass(i)*(  atmvel(1,i)*atmvel(1,i) &
          &                    + atmvel(2,i)*atmvel(2,i) &
          &                    + atmvel(3,i)*atmvel(3,i))
  end do

! --- update glogs(1) ---
  glogs(1) = (kint + &
           &  vmass*(vboxg(1)*vboxg(1)+vboxg(2)*vboxg(2)+vboxg(3)*vboxg(3)) &
           & - gndfkt) / qmass(1)

! --- update gboxv ---
  gboxv(1:3) = (dnf*kint + vol*(pintv(1:3)- pext_c)) / vmass

! --- Start MTS in NHCPaniso procedure ---

  DO inc = 1, next
     DO iys = 1, nyosh

!       --- update thermostat velocities ---
!       last thermostat
        glogs(mchain) = (qmass(mchain-1) &
             &           *vlogs(mchain-1)*vlogs(mchain-1) &
             &           - gkt) / qmass(mchain)
                                ! G(M)=(Q(M-1)*vxi(M-1)**2-kT)/Q(M)
        vlogs(mchain) = vlogs(mchain) + glogs(mchain)*wdti4(iys)

!       M-1 to 2 thermostat
        do i=mchain-1,2,-1
           glogs(i) = (qmass(i-1) &
                &      *vlogs(i-1)*vlogs(i-1) &
                &      - gkt) / qmass(i)
           expwdti8 = exp(-1.0d0*wdti8(iys)*vlogs(i+1))
           vlogs(i) = vlogs(i) * expwdti8*expwdti8 &
                  & + wdti4(iys)*glogs(i)*expwdti8
        end do

!       first thermostat (to physical system)
        expwdti8 = exp(-1.0d0*wdti8(iys)*vlogs(2))
        vlogs(1) = vlogs(1) * expwdti8*expwdti8 &
               & + wdti4(iys)*glogs(1)*expwdti8

!       --- update barostat velocity ---
        expwdti8 = exp(-1.0d0*wdti8(iys)*vlogs(1))
        if (pcont_dir == 0) then   ! full anisotropic NPT
            vboxg(1:3) = vboxg(1:3) * expwdti8*expwdti8 &
                     & + wdti4(iys)*gboxv(1:3)*expwdti8
        else   ! uniaxial NPT
            vboxg(pcont_dir) = vboxg(pcont_dir) * expwdti8*expwdti8 &
                     & + wdti4(iys)*gboxv(pcont_dir)*expwdti8
        end if

!       --- update thermostat positions ---
        xlogs(1:mchain) = xlogs(1:mchain) + vlogs(1:mchain)*wdti2(iys)

!       --- update barostat velocity ---
!        expwdti8 = exp(-1.0d0*wdti8(iys)*vlogs(1))
        if (pcont_dir == 0) then   ! full anisotropic NPT
            vboxg(1:3) = vboxg(1:3) * expwdti8*expwdti8 &
                     & + wdti4(iys)*gboxv(1:3)*expwdti8
        else   ! uniaxial NPT
            vboxg(pcont_dir) = vboxg(pcont_dir) * expwdti8*expwdti8 &
                     & + wdti4(iys)*gboxv(pcont_dir)*expwdti8
        end if

!       --- update thermostat velocities ---
!       --- update glogs(1) ---
        glogs(1) = (kint + &
            &  vmass*(vboxg(1)*vboxg(1)+vboxg(2)*vboxg(2)+vboxg(3)*vboxg(3)) &
            & - gndfkt) / qmass(1)

!       first thermostat (to physical system)
        expwdti8 = exp(-1.0d0*wdti8(iys)*vlogs(2))
        vlogs(1) = vlogs(1) * expwdti8*expwdti8 &
               & + wdti4(iys)*glogs(1)*expwdti8

!       2 to M-1 thermostat
        do i=2,mchain-1
           glogs(i) = (qmass(i-1) &
                &      *vlogs(i-1)*vlogs(i-1) &
                &      - gkt) / qmass(i)
           expwdti8 = exp(-1.0d0*wdti8(iys)*vlogs(i+1))
           vlogs(i) = vlogs(i) * expwdti8*expwdti8 &
                  & + wdti4(iys)*glogs(i)*expwdti8
        end do

!       last thermostat
        glogs(mchain) = (qmass(mchain-1) &
             &           *vlogs(mchain-1)*vlogs(mchain-1) &
             &           - gkt) / qmass(mchain)
                                ! G(M)=(Q(M-1)*vxi(M-1)**2-kT)/Q(M)
        vlogs(mchain) = vlogs(mchain) + glogs(mchain)*wdti4(iys)

     END DO
  END DO

!     +     +     +     +     +     +     +     +

end subroutine nhcpanisoint
