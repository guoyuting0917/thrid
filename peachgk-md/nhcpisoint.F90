!**************************************
!*  nhcpisoint_.f Ver.1.2 '10.04.09   *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine nhcpisoint( degfree_all, &
     &                 dt_long_cal, &
     &                 mchain,text_c, &
     &                 next,nyosh, &
     &                 pint,pext_c )

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

  real(8),intent(in):: pint             ! internal pressure in NPT(MTK)
  real(8),intent(in):: pext_c           ! external pressure in NPT(MTK)

! LOCAL:
  real(8):: kint             ! kinetic energy of physical space (=mv**2)

  real(8):: wdti2(maxnyosh)  ! yscoeff()*dt/next/2
  real(8):: wdti4(maxnyosh)  ! yscoeff()*dt/next/4
  real(8):: wdti8(maxnyosh)  ! yscoeff()*dt/next/8
      
  real(8):: glogs(maxnhc)    ! Glogs() in thermostat velocity update
  real(8):: glogv            ! Glogv in barostat velocity update

  real(8):: gn1kt            ! (Nf+1)kT
  real(8):: gkt              ! kT

  integer:: i               ! do loop index
  integer:: inc             ! do loop index of iteration number (next)
  integer:: iys             ! do loop index of Yoshida-Suzuki method (nyosh)

  real(8):: expwdti8         ! = exp(wdti8*vxi(m))

  real(8):: dnf              ! = d/Nf
      
  real(8):: vol              ! = exp(d*xlogv)

!     +     +     +     +     +     +     +     +

!     --- some preparation ---

!     substitute wdtis
  wdti2(1:nyosh) = yscoeff(1:nyosh) * dt_long_cal / dble(next) / 2.0d0
  wdti4(1:nyosh) = wdti2(1:nyosh) / 2.0d0
  wdti8(1:nyosh) = wdti4(1:nyosh) / 2.0d0

!     calculate NfkT and kT
  gn1kt = text_c * dble(degfree_all+1)
  gkt  = text_c

!     calculate dnf
  dnf = 3.0d0/dble(degfree_all)

!     calculate vol
  vol = dexp(3.0d0*xlogv)

!     ----- Main loop of Nose-Hoover chain - Andersen (Hoover type) -----

!     --- calculate kinetic energy of physical system ---
  kint = 0.0d0
  do i=1,natom
     kint = kint + atmmass(i)*(  atmvel(1,i)*atmvel(1,i) &
          &                    + atmvel(2,i)*atmvel(2,i) &
          &                    + atmvel(3,i)*atmvel(3,i))
  end do

!     --- update glogs(1) ---
  glogs(1) = (kint + vmass*vlogv*vlogv - gn1kt) / qmass(1)

!     --- update glogv ---
  glogv = (dnf*kint + 3.0d0*vol*(pint- pext_c)) / vmass

!     --- Start MTS in NHCPiso procedure ---

  DO inc = 1, next
     DO iys = 1, nyosh

!           --- update thermostat velocities ---
!           last thermostat
        glogs(mchain) = (qmass(mchain-1) &
             &           *vlogs(mchain-1)*vlogs(mchain-1) &
             &           - gkt) / qmass(mchain)
                                ! G(M)=(Q(M-1)*vxi(M-1)**2-kT)/Q(M)
        vlogs(mchain) = vlogs(mchain) + glogs(mchain)*wdti4(iys)

!           M-1 to 2 thermostat
        do i=mchain-1,2,-1
           glogs(i) = (qmass(i-1) &
                &      *vlogs(i-1)*vlogs(i-1) &
                &      - gkt) / qmass(i)
           expwdti8 = dexp(-1.0d0*wdti8(iys)*vlogs(i+1))
           vlogs(i) =  vlogs(i) * expwdti8*expwdti8 &
                &    + wdti4(iys)*glogs(i)*expwdti8
        end do

!           first thermostat (to physical system)
        expwdti8 = dexp(-1.0d0*wdti8(iys)*vlogs(2))
        vlogs(1) =  vlogs(1) * expwdti8*expwdti8 &
             &    + wdti4(iys)*glogs(1)*expwdti8

!           --- update barostat velocity ---
        expwdti8 = dexp(-1.0d0*wdti8(iys)*vlogs(1))
        vlogv =  vlogv * expwdti8*expwdti8 &
             & + wdti4(iys)*glogv*expwdti8

!           --- update thermostat positions ---
        xlogs(1:mchain) = xlogs(1:mchain) + vlogs(1:mchain)*wdti2(iys)
           
!           --- update barostat velocity ---
!            expwdti8 = dexp(-1.0d0*wdti8(iys)*vlogs(1))
        vlogv =  vlogv * expwdti8*expwdti8 &
             & + wdti4(iys)*glogv*expwdti8

!           --- update thermostat velocities ---
!           --- update glogs(1) ---
        glogs(1) = (kint + vmass*vlogv*vlogv - gn1kt) / qmass(1)

!           first thermostat (to physical system)
        expwdti8 = exp(-1.0d0*wdti8(iys)*vlogs(2))
        vlogs(1) = vlogs(1) * expwdti8*expwdti8 &
             &   + wdti4(iys)*glogs(1)*expwdti8

!           2 to M-1 thermostat
        do i=2,mchain-1
           glogs(i) = (qmass(i-1) &
                &      *vlogs(i-1)*vlogs(i-1) &
                &      - gkt) / qmass(i)
           expwdti8 = exp(-1.0d0*wdti8(iys)*vlogs(i+1))
           vlogs(i) =  vlogs(i) * expwdti8*expwdti8 &
                &    + wdti4(iys)*glogs(i)*expwdti8
        end do

!           last thermostat
        glogs(mchain) = (qmass(mchain-1) &
             &           *vlogs(mchain-1)*vlogs(mchain-1) &
             &           - gkt) / qmass(mchain)
                                ! G(M)=(Q(M-1)*vxi(M-1)**2-kT)/Q(M)
        vlogs(mchain) = vlogs(mchain) + glogs(mchain)*wdti4(iys)

     END DO
  END DO

!     +     +     +     +     +     +     +     +

  return
end subroutine nhcpisoint
