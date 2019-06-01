!**************************************
!*  upewk_.f Ver.1.3 '11.01.31        *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine upewk( alpha, &
     &            xcel, ycel, zcel )

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

  real(8),intent(in):: xcel             ! x cell length [non-d]
  real(8),intent(in):: ycel             ! y cell length [non-d]
  real(8),intent(in):: zcel             ! z cell length [non-d]

! LOCAL:
  integer:: i                ! do loop index
  real(8):: kwave_abs2       ! |kwave|**2
  real(8):: pi               ! 3,1415...  
  real(8):: coeff            ! coefficient
  real(8):: rtmp             ! for temporal use 
  real(8):: box_inv(3)       ! inverse of cell length (rectangular)

!   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!---- some preparation
  pi = dacos(-1.0d0)

  box_inv(1) = 1.0d0 / xcel
  box_inv(2) = 1.0d0 / ycel
  box_inv(3) = 1.0d0 / zcel

!   --- CALCULATE THE COEFFICIENT A(k) ---
!             A(k) = exp(- (pi/alpha* k/L)^2)/(k/L)^2
!       or    A(K) = exp(- (pi/alpha* K  )^2)/    K^2

  coeff = (pi/alpha) ** 2
  awave(1)   = 0.0d0
 
  do i = 2, nwave
     kwave_abs2 =  (nvec(1,i)/xcel)**2 + (nvec(2,i)/ycel)**2 &
          &      + (nvec(3,i)/zcel)**2 
     rtmp       = dexp(-coeff * kwave_abs2)
     awave(i)   = rtmp/kwave_abs2
  end do


!   --- K = (kx/Lx, ky/Ly, kz/Lz) is set to kwave
!       instead of             k = (kx,    ky,    kz) ---

  do i = 1, nwave
     kwave(1:3,i) = nvec(1:3,i)*box_inv(1:3)
  end do

!----------------------------------------------------------------------

  return
end subroutine upewk
