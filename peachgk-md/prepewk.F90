!******************************************
!*  prepewk.f90 Ver.1.6                   *
!*      for peachgk_md.f                  *
!*            by G.Kikugawa               *
!*  Orinially coded by Y.Komeiji (PEACH)  *
!******************************************
! Time-stamp: <2015-06-02 17:18:08 gota>

subroutine prepewk(pot_ewc, alpha, kmax, &
     &             xcel, ycel, zcel, &
     &             yratio, zratio)

  use md_global

  implicit none

!    This subroutine prepares Ewald wave number vectors
!    Revision 2.01
!     modified to expand to the rectanglar box
! 
!    Ewald wave number space 
!
!          Potential = 1/2pi/V * sigma(k) A(k) *
!                    [ (sigma(i) qi cos(2pi/L*k.ri))^2
!                    + (sigma(i) qi sin(2pi/L*k.ri))^2 ]
!          Force(i)  = 2qi sigma(K) A(K) *
!                    [ sin(2pi/L k.ri) sigma(j) qj cos(2pi/L*k.rj)
!                    - cos(2pi/L k.ri) sigma(j) qj sin(2pi/L*k.rj) ] K
!
!       where L    = box size
!             k    = (kx, ky, kz); 0, +-1, +-2, +-3,.....
!             A(k) = exp(- (pi*ehta* k/L)^2)/(k^2/L^2)
!
!       For General purpose computers, 
!             K    = (kx/Lx, ky/Ly, kz/Lz) will be transferred to kwave
!             
!             In both cases, A(k) = A(K) = exp(- (pi*ehta*K)^2)/K^2)
!             will be used.
!             
! ARGUMENT:
!     INPUT
  real(8),intent(in):: alpha            ! parameter alpha [non-d]
  integer,intent(in):: kmax           ! parameter kmax

  real(8),intent(in):: xcel             ! x cell length [non-d]
  real(8),intent(in):: ycel             ! y cell length [non-d]
  real(8),intent(in):: zcel             ! z cell length [non-d]

  real(8),intent(in):: yratio           ! y cell ratio of y to x
  real(8),intent(in):: zratio           ! z cell ratio of z to x

!    OUTPUT
  real(8),intent(out):: pot_ewc          ! potential of ewald self-energy [non-d]
      
! LOCAL:
  integer:: i               ! do loop index
  integer:: nwave_app       ! approximate number of wave vectors
  integer:: kmax2           ! kmax**2 
  integer:: kx,ky,kz        ! x,y,z of the wave vectors
  integer:: ykmax,zkmax     ! ykmax,zkmax
  integer:: dnyratio,dnzratio

  real(8)::  kwave_abs2      ! |kwave|**2
  real(8)::  pi              ! 3,1415...  
  real(8)::  coeff           ! coefficient
  real(8)::  rtmp            ! for temporal use 

      
!   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!---- some preparation

  dnyratio = anint(yratio)
  dnzratio = anint(zratio)
  ykmax = kmax * dnyratio
  zkmax = kmax * dnzratio

!   --- CALCULATE CONSTANT PART OF EWALD POTENTIAL ---
!          -1/sqrt(pi)*alpha * sigmj(qj^2)          
! 

  pi = acos(-1.0d0)

  pot_ewc = 0.0d0
        
  do i = 1, natom
     pot_ewc = pot_ewc + atmchrg(i)**2
  end do

!---- use at standard ewald
  pot_ewc = -1.0d0/sqrt(pi)*alpha * pot_ewc 
!---- use at ewald-chart
!      pot_ewc = -1.0d0/dsqrt(pi)*alpha/xcel * pot_ewc 


!   --- PREPARE THE WAVE NUMBER VECTORS ---

!     - find approximate number of vectors and allocate -

  nwave_app = 5 * kmax**3

  kmax2 = kmax * kmax

!     - make wave vectors kwave (ktmp) -
      
  nvec(1:3,1) = 0.0d0       ! nwave = 1 is reserved for (0,0,0)
  nwave = 1                 ! count of kwave (ktmp)


  DO kx = -kmax, kmax
     DO ky = -ykmax, ykmax
        DO kz = -zkmax, zkmax

           IF (.not.(kx == 0 .and. ky == 0 .and. kz == 0)) THEN

              kwave_abs2 =  dble(kx*kx) &
                   &     + dble(ky*ky)/dble(dnyratio*dnyratio) &
                   &     + dble(kz*kz)/dble(dnzratio*dnzratio)
              if (kwave_abs2 <= dble(kmax2)) then
                 nwave = nwave + 1
                 nvec(1,nwave) = dble(kx)
                 nvec(2,nwave) = dble(ky)
                 nvec(3,nwave) = dble(kz)
              end if

           END IF

        END DO
     END DO
  END DO

!   --- PREPARE THE COEFFICIENT A(k) ---
!             A(k) = exp(- (pi/alpha* k/L)^2)/(k/L)^2
!       or    A(K) = exp(- (pi/alpha* K  )^2)/    K^2

  coeff = (pi/alpha) ** 2
  awave(1)   = 0.0d0
 
  do i = 2, nwave
     kwave_abs2 =  (nvec(1,i)/xcel)**2 + (nvec(2,i)/ycel)**2 &
          &     + (nvec(3,i)/zcel)**2 
     rtmp       = dexp(-coeff * kwave_abs2)
     awave(i)   = rtmp/kwave_abs2
  end do


!   --- K = (kx/Lx, ky/Ly, kz/Lz) is set to kwave
!       instead of             k = (kx,    ky,    kz) ---

  do i = 1, nwave
     kwave(1,i) = nvec(1,i)/xcel
     kwave(2,i) = nvec(2,i)/ycel
     kwave(3,i) = nvec(3,i)/zcel
  end do

!----------------------------------------------------------------------

end subroutine prepewk
