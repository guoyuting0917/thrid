!**********************************
!*  fft_pme.f Ver.1.5 '10.06.28   *
!*      for peachgk_md.f          *
!*            by G.Kikugawa       *
!**********************************

!---- original code Copyright
!*-*
!*-* Copyright (C) 2000 Massimo Marchi and Piero Procacci
!*-* Full copyright notice at http://www.chim.unifi.it/orac/copyright4.0.html
!*-* Contact for info M. Marchi, CEA,  Gif Sur Yvette 91191 (FRANCE)
!*-* Email:marchi@villon.saclay.cea.fr
!*-*

subroutine fft_pme(xcel,ycel,zcel,   &
     &             alpha,   &
     &             nfft1,nfft2,nfft3,pme_order,   &
     &             ifnetqcorrp,   &
     &             netchrgsq,   &
     &             force,pot_ewk)

  use interface_interact, only: get_scaled_fractionals, get_bspline_coeffs,   &
       &                        fill_charge_grid, fft_back, fft_forward,   &
       &                        scalar_sum, grad_sum

  use md_global
  use spme_global

  implicit none

! ARGUMENT:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  real(8),intent(in):: alpha            ! parameter alpha [non-d]

  integer,intent(in):: nfft1, nfft2, nfft3 ! grid points in SPME
  integer,intent(in):: pme_order        ! B-spline order

  logical,intent(in):: ifnetqcorrp      ! net charge correction for pressure
  real(8),intent(in):: netchrgsq        ! = (sum(qi))**2

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! total force

!   OUTPUT
  real(8),intent(out):: pot_ewk         ! electrostatic potential

! LOCAL:
  real(8):: theta1(mth),dtheta1(mth),theta2(mth),dtheta2(mth),   &
       &    theta3(mth),dtheta3(mth)
  real(8):: Q(maxt)
  real(8):: fr1(max_atm),fr2(max_atm),fr3(max_atm)
                                        ! scaled fractional coordinate

  real(8):: vir(3,3)                    ! virial tensor
  real(8):: vir_sca                     ! virial scalar for charge correction

  real(8):: box(3)
  real(8):: recip(3)
  real(8):: vol

  integer:: nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw

  real(8):: vir_for(3,natom)

!     +     +     +     +     +     +     +

! --- SOME PREPARATIONS ---

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel

  vol = box(1) * box(2) * box(3)

  vir(1:3,1:3) = 0.0d0

  vir_for(1:3,1:natom) = 0.0d0

!-- call stand-alone Darden routine

  recip(1:3) = 1.0d0/box(1:3)

!---- get some integer array dimensions
  call get_fftdims(nfft1,nfft2,nfft3,   &
       &           nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
                                        ! in fft_pme_init.f

  call get_scaled_fractionals(nfft1,nfft2,nfft3,   &
       &                      recip,   &
       &                      fr1,fr2,fr3)

  call get_bspline_coeffs(pme_order,   &
       &                  fr1,fr2,fr3,   &
       &                  theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)

  call fill_charge_grid(theta1,theta2,theta3,fr1,fr2,fr3,   &
       &                pme_order,nfft1,nfft2,nfft3,   &
       &                nfftdim1,nfftdim2,nfftdim3,Q)

#if defined(_FFTW3)
  call fft_back(Q,   &
       &        nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3)
#else
  call fft_back(Q,   &
       &        nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,   &
       &        nfftable,nffwork)
#endif

  call scalar_sum(alpha,pot_ewk,   &
       &          Q,vol,recip,   &
       &          nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,vir,   &
       &          vir_sca,   &
       &          ifnetqcorrp,   &
       &          netchrgsq)

#if defined(_FFTW3)
  call fft_forward(Q,   &
       &           nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3)
#else
  call fft_forward(Q,   &
       &           nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,   &
       &           nfftable,nffwork)
#endif

  call grad_sum(recip,   &
       &        theta1,theta2,theta3,   &
       &        dtheta1,dtheta2,dtheta3,   &
       &        force,   &
       &        fr1,fr2,fr3,   &
       &        pme_order,nfft1,nfft2,nfft3,   &
       &        nfftdim1,nfftdim2,nfftdim3,   &
       &        Q,   &
       &        vir_for)

!     +     +     +     +     +     +     +

end subroutine fft_pme

!----------------------------------------------------------------------
subroutine fft_pmep(xcel,ycel,zcel,   &
     &              alpha,   &
     &              nfft1,nfft2,nfft3,pme_order,   &
     &              force,pot_ewk,   &
     &              ifnetqcorrp,   &
     &              netchrgsq,   &
     &              for_viri_coul,pot_viri_coul,pot_virit_coul)

  use interface_interact, only: get_scaled_fractionals, get_bspline_coeffs,   &
       &                        fill_charge_grid, fft_back, fft_forward,   &
       &                        scalar_sum, grad_sum

  use md_global
  use spme_global

  implicit none

! ARGUMENT:
!   INPUT
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  real(8),intent(in):: alpha            ! parameter alpha [non-d]

  integer,intent(in):: nfft1, nfft2, nfft3 ! grid points in SPME
  integer,intent(in):: pme_order        ! B-spline order

  logical,intent(in):: ifnetqcorrp      ! net charge correction for pressure
  real(8),intent(in):: netchrgsq        ! = (sum(qi))**2

!   INPUT & OUTPUT
  real(8),intent(inout):: force(:,:)    ! total force

!   OUTPUT
  real(8),intent(out):: pot_ewk         ! electrostatic potential

  real(8),intent(out):: for_viri_coul(:,:) ! virial(coulomb force) of each atom
  real(8),intent(out):: pot_viri_coul   ! virial(coulomb potential) of each atom
  real(8),intent(out):: pot_virit_coul(:,:) ! virial tensor (coulomb)

! LOCAL:
  real(8):: theta1(mth),dtheta1(mth),theta2(mth),dtheta2(mth),   &
       &    theta3(mth),dtheta3(mth)
  real(8):: Q(maxt)
  real(8):: fr1(max_atm),fr2(max_atm),fr3(max_atm)
                                        ! scaled fractional coordinate

  real(8):: box(3)
  real(8):: recip(3)
  real(8):: vol

  integer:: nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw

  real(8):: vir_sca                     ! virial scalar for charge correction

!     +     +     +     +     +     +     +

! --- SOME PREPARATIONS ---

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel

  vol = box(1) * box(2) * box(3)

!--   call stand-alone Darden routine

  recip(1:3) = 1.0d0/box(1:3)

!---- get some integer array dimensions
  call get_fftdims(nfft1,nfft2,nfft3,   &
       &       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
                                        ! in fft_pme_init.f

  call get_scaled_fractionals(nfft1,nfft2,nfft3,   &
       &                      recip,   &
       &                      fr1,fr2,fr3)

  call get_bspline_coeffs(pme_order,   &
       &                  fr1,fr2,fr3,   &
       &                  theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)

  call fill_charge_grid(theta1,theta2,theta3,fr1,fr2,fr3,   &
       &                pme_order,nfft1,nfft2,nfft3,   &
       &                nfftdim1,nfftdim2,nfftdim3,Q)

#if defined(_FFTW3)
  call fft_back(Q,   &
       &        nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3)
#else
  call fft_back(Q,   &
       &        nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,   &
       &        nfftable,nffwork)
#endif

  call scalar_sum(alpha,pot_ewk,   &
       &          Q,vol,recip,   &
       &          nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,pot_virit_coul, &
       &          vir_sca,   &
       &          ifnetqcorrp,   &
       &          netchrgsq)

#if defined(_FFTW3)
  call fft_forward(Q,   &
       &           nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3)
#else
  call fft_forward(Q,   &
       &           nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,   &
       &           nfftable,nffwork)
#endif

  call grad_sum(recip,   &
       &        theta1,theta2,theta3,   &
       &        dtheta1,dtheta2,dtheta3,   &
       &        force,   &
       &        fr1,fr2,fr3,   &
       &        pme_order,nfft1,nfft2,nfft3,   &
       &        nfftdim1,nfftdim2,nfftdim3,   &
       &        Q,   &
       &        for_viri_coul)

  pot_viri_coul = pot_viri_coul + pot_ewk + vir_sca

!     +     +     +     +     +     +     +

end subroutine fft_pmep

!----------------------------------------------------------------------
subroutine get_scaled_fractionals(nfft1,nfft2,nfft3,   &
     &                            recip,   &
     &                            fr1,fr2,fr3)

  use md_global

  implicit none

! ARGUMENT:
!   INPUT
  integer,intent(in):: nfft1, nfft2, nfft3 ! number of grid points in SPME
  real(8),intent(in):: recip(:)

!   OUTPUT
  real(8),intent(out):: fr1(:),fr2(:),fr3(:)

! LOCAL:
  integer:: n
  real(8):: w

!     +     +     +     +     +     +     +

  do n = 1, natom

     w = atmcor(1,n)*recip(1)
     w = w - floor(w)
     fr1(n) = nfft1*w

     w = atmcor(2,n)*recip(2)
     w = w - floor(w)
     fr2(n) = nfft2*w

     w = atmcor(3,n)*recip(3)
     w = w - floor(w)
     fr3(n) = nfft3*w

  end do

!     +     +     +     +     +     +     +

end subroutine get_scaled_fractionals

!---------------------------------------------------------------------
subroutine get_bspline_coeffs(pme_order,   &
     &                        fr1,fr2,fr3,   &
     &                        theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
!---------------------------------------------------------------------
! INPUT:
!      fr1,fr2,fr3 the scaled and shifted fractional coords
! OUTPUT
!      theta1,theta2,theta3: the spline coeff arrays
!      dtheta1,dtheta2,dtheta3: the 1st deriv of spline coeff arrays
!---------------------------------------------------------------------

  use md_global

  implicit none

  integer,intent(in):: pme_order       ! B-spline order
  real(8),intent(in):: fr1(:),fr2(:),fr3(:)

  real(8),intent(out):: theta1(pme_order,natom),theta2(pme_order,natom),   &
       &                theta3(pme_order,natom)
  real(8),intent(out):: dtheta1(pme_order,natom),   &
       &                dtheta2(pme_order,natom),dtheta3(pme_order,natom)

  real(8):: w
  integer:: n

!     +     +     +     +     +     +     +     +     +

  do n = 1,natom
     w = fr1(n)-dint(fr1(n))
     call fill_bspline(w,pme_order,theta1(1,n),dtheta1(1,n))
     w = fr2(n)-dint(fr2(n))
     call fill_bspline(w,pme_order,theta2(1,n),dtheta2(1,n))
     w = fr3(n)-dint(fr3(n))
     call fill_bspline(w,pme_order,theta3(1,n),dtheta3(1,n))
  end do
!     subroutine fill_bspline is in fft_pme_init.f

!     +     +     +     +     +     +     +     +     +

end subroutine get_bspline_coeffs

!------------------------------------------------------------------------
subroutine fill_charge_grid(theta1,theta2,theta3,fr1,fr2,fr3,   &
     &                      pme_order,nfft1,nfft2,nfft3,   &
     &                      nfftdim1,nfftdim2,nfftdim3,Q)
!---------------------------------------------------------------------
! INPUT:
!      theta1,theta2,theta3: the spline coeff arrays
!      fr1,fr2,fr3 the scaled and shifted fractional coords
!      nfft1,nfft2,nfft3: the charge grid dimensions
!      nfftdim1,nfftdim2,nfftdim3: physical charge grid dims
!      pme_order: the order of spline interpolation
! OUTPUT:
!      Q the charge grid
!---------------------------------------------------------------------

  use md_global

  implicit none

  integer,intent(in):: pme_order,nfft1,nfft2,nfft3
  integer,intent(in):: nfftdim1,nfftdim2,nfftdim3
  real(8),intent(in):: fr1(:),fr2(:),fr3(:)
  real(8),intent(out):: Q(2,nfftdim1,nfftdim2,nfftdim3)

  real(8),intent(in):: theta1(pme_order,natom),theta2(pme_order,natom),   &
       &               theta3(pme_order,natom)

  integer:: n,ith1,ith2,ith3,i0,j0,k0,i,j,k
  real(8):: prod, prod1

!     +     +     +     +     +     +     +     +     +

  Q(1:2,1:nfftdim1,1:nfftdim2,1:nfftdim3) = 0.0d0

  do n = 1,natom
     k0 = int(fr3(n)) - pme_order

     do ith3 = 1,pme_order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        prod1 = theta3(ith3,n) * atmchrg(n)

        j0 = int(fr2(n)) - pme_order

        do ith2 = 1,pme_order
           j0 = j0 + 1
           j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
           prod = theta2(ith2,n) * prod1
           i0 = int(fr1(n)) - pme_order

           do ith1 = 1,pme_order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
              Q(1,i,j,k) = Q(1,i,j,k) + theta1(ith1,n) * prod
           end do

        end do

     end do

  end do

!     +     +     +     +     +     +     +     +     +

end subroutine fill_charge_grid

!-----------------------------------------------------------
#if defined(_FFTW3)
subroutine fft_back(array,   &
     &              nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3)

  use spme_global

  implicit none

  include 'fftw3.f'

  real(8),intent(out):: array(2,nfftdim1,nfftdim2,nfftdim3)

  integer,intent(in):: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

  integer:: dim1, dim2, dim3

!     +     +     +     +     +     +     +     +     +

!
! 1) X-direction FFT
!
  do dim3=1,nfft3
     do dim2=1,nfft2
        do dim1=1,nfft1

           fftWork1(dim1) = CMPLX(array(1,dim1,dim2,dim3),   &
                &                 array(2,dim1,dim2,dim3),kind=8)

        enddo

        call DFFTW_EXECUTE(planB(1))

        do dim1=1,nfft1
           array(1,dim1,dim2,dim3) = DBLE(fftWork1(dim1))
           array(2,dim1,dim2,dim3) = AIMAG(fftWork1(dim1))
        enddo

     enddo
  enddo

!
! 2) Y-direction FFT
!
  do dim1=1,nfft1
     do dim3=1,nfft3
        do dim2=1,nfft2

           fftWork1(dim2) = CMPLX(array(1,dim1,dim2,dim3),   &
                &                 array(2,dim1,dim2,dim3),kind=8)

        enddo

        call DFFTW_EXECUTE(planB(2))

        do dim2=1,nfft2
           array(1,dim1,dim2,dim3) = DBLE(fftWork1(dim2))
           array(2,dim1,dim2,dim3) = AIMAG(fftWork1(dim2))
        enddo

     enddo
  enddo

!
! 3) Z-direction FFT
!
  do dim2=1,nfft2
     do dim1=1,nfft1
        do dim3=1,nfft3

           fftWork1(dim3) = CMPLX(array(1,dim1,dim2,dim3),   &
                &                 array(2,dim1,dim2,dim3),kind=8)

        enddo

        call DFFTW_EXECUTE(planB(3))

        do dim3=1,nfft3
           array(1,dim1,dim2,dim3) = DBLE(fftWork1(dim3))
           array(2,dim1,dim2,dim3) = AIMAG(fftWork1(dim3))
        enddo

     enddo
  enddo

!     +     +     +     +     +     +     +     +     +

end subroutine fft_back

#else
!-----------------------------------------------------------
subroutine fft_back(array,   &
     &              nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,   &
     &              nfftable,nffwork)

  use spme_global

  implicit none

  real(8),intent(out):: array(:)
  integer,intent(in):: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

  integer,intent(out):: nfftable,nffwork

  integer:: isign

!     +     +     +     +     +     +     +     +     +

! implement forward fft3d
  isign = -1

  call pubz3d(isign,nfft1,nfft2,nfft3,array,   &
       &      nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)

!     +     +     +     +     +     +     +     +     +

end subroutine fft_back
#endif

!----------------------------------------------------
subroutine scalar_sum(alpha,pot_ewk,   &
     &     Q,vol,recip,   &
     &     nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,vir,   &
     &     vir_sca,   &
     &     ifnetqcorrp,   &
     &     netchrgsq)

  use spme_global

  implicit none

  real(8),intent(in):: alpha           ! parameter alpha [non-d]

  integer,intent(in):: nfft1, nfft2, nfft3 ! grid points in SPME
  integer,intent(in):: nfftdim1,nfftdim2,nfftdim3

  real(8),intent(out):: pot_ewk        ! electrostatic potential

  real(8),intent(inout):: Q(2,nfftdim1,nfftdim2,nfftdim3)

  real(8),intent(in):: vol
  real(8),intent(in):: recip(:)

  real(8),intent(out):: vir(:,:)       ! virial tensor
  real(8),intent(out):: vir_sca        ! virial scalar for charge correction

  logical,intent(in):: ifnetqcorrp     ! net charge correction for pressure
  real(8),intent(in):: netchrgsq       ! = (sum(qi))**2

  real(8):: pi
  real(8):: fac                        ! = pi^2/alpha^2
  real(8):: denom                      ! = pi*vol*msq/B
  real(8):: eterm                      ! = B/(pi*vol)*exp(-fac*msq)/msq
  real(8):: eterm2                     ! = 1/(pi*vol)*netchrgsq*exp(-fac*msq)/msq
  real(8):: eterm3                     ! = eterm2*(-1 + 2*pi^2/alpha^2*msq)
  real(8):: vterm                      ! = 2/msq*(pi^2/alpha^2*msq+1)
  real(8):: energy                     ! = 1/(pi*vol)*exp(-fac*msq)/msq*|S(m)|^2

  integer:: k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
  integer:: nf1,nf2,nf3
  real(8):: mhat1,mhat2,mhat3,msq,struc2,rkcut2
  real(8):: vir_tmp(3,3)

!     +     +     +     +     +     +     +     +     +

  vir_tmp(1:3,1:3) = 0.0d0
  vir_sca = 0.0d0

  indtop = nfft1*nfft2*nfft3
  pi = 3.14159265358979323846d0
  rkcut2=rkcut*rkcut
  fac = pi**2/alpha**2
  nff = nfft1*nfft2
  nf1 = nfft1/2
  if ( 2*nf1 < nfft1 ) nf1 = nf1+1
  nf2 = nfft2/2
  if ( 2*nf2 < nfft2 ) nf2 = nf2+1
  nf3 = nfft3/2
  if ( 2*nf3 < nfft3 ) nf3 = nf3+1
  energy = 0.d0

  do ind = 1,indtop-1

! get k1,k2,k3 from the relationship
!          ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1

     k3 = ind/nff + 1
     jnd = ind - (k3-1)*nff
     k2 = jnd/nfft1 + 1
     k1 = jnd - (k2-1)*nfft1 +1
     m1 = k1 - 1
     if ( k1 > nf1 )m1 = k1 - 1 - nfft1
     m2 = k2 - 1
     if ( k2 > nf2 )m2 = k2 - 1 - nfft2
     m3 = k3 - 1
     if ( k3 > nf3 )m3 = k3 - 1 - nfft3
     mhat1 = recip(1)*m1
     mhat2 = recip(2)*m2
     mhat3 = recip(3)*m3
     msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3

     IF(msq > rkcut2) THEN
        eterm = 0.d0
     ELSE

        if (.not. ifnetqcorrp) then    ! no charge correction

           denom = pi*vol*bsp_mod1(k1)*bsp_mod2(k2)*bsp_mod3(k3)*msq
           eterm = dexp(-fac*msq)/denom
           vterm = 2.d0*(fac*msq + 1.d0)/msq
           struc2 = Q(1,k1,k2,k3)**2 + Q(2,k1,k2,k3)**2
           energy = energy + eterm * struc2
           vir_tmp(1,1) = vir_tmp(1,1) + eterm * struc2   &
                &       * (vterm*mhat1*mhat1 - 1.d0)
           vir_tmp(1,2) = vir_tmp(1,2) + eterm * struc2   &
                &       * (vterm*mhat1*mhat2)
           vir_tmp(1,3) = vir_tmp(1,3) + eterm * struc2   &
                &       * (vterm*mhat1*mhat3)
           vir_tmp(2,2) = vir_tmp(2,2) + eterm * struc2   &
                &       * (vterm*mhat2*mhat2 - 1.d0)
           vir_tmp(2,3) = vir_tmp(2,3) + eterm * struc2   &
                &       * (vterm*mhat2*mhat3)
           vir_tmp(3,3) = vir_tmp(3,3) + eterm * struc2   &
                &       * (vterm*mhat3*mhat3 - 1.d0)

        else                           ! charge correction

           denom = pi*vol*bsp_mod1(k1)*bsp_mod2(k2)*bsp_mod3(k3)*msq
           eterm2 = dexp(-fac*msq)
           eterm = eterm2 / denom
           eterm2 = eterm2 * netchrgsq / (pi*vol*msq)
           eterm3 = eterm2 * (-1.0d0 + 2.0d0*fac*msq)
           vterm = 2.d0*(fac*msq + 1.d0)/msq
           struc2 = Q(1,k1,k2,k3)**2 + Q(2,k1,k2,k3)**2
           energy = energy + eterm * struc2

           vir_sca = vir_sca + eterm3
           vir_tmp(1,1) = vir_tmp(1,1) + (eterm * struc2 - eterm2)   &
                &       * (vterm*mhat1*mhat1 - 1.d0)
           vir_tmp(1,2) = vir_tmp(1,2) + (eterm * struc2 - eterm2)   &
                &       * (vterm*mhat1*mhat2)
           vir_tmp(1,3) = vir_tmp(1,3) + (eterm * struc2 - eterm2)   &
                &       * (vterm*mhat1*mhat3)
           vir_tmp(2,2) = vir_tmp(2,2) + (eterm * struc2 - eterm2)   &
                &       * (vterm*mhat2*mhat2 - 1.d0)
           vir_tmp(2,3) = vir_tmp(2,3) + (eterm * struc2 - eterm2)   &
                &       * (vterm*mhat2*mhat3)
           vir_tmp(3,3) = vir_tmp(3,3) + (eterm * struc2 - eterm2)   &
                &       * (vterm*mhat3*mhat3 - 1.d0)
        end if

     END IF
     Q(1:2,k1,k2,k3) = eterm * Q(1:2,k1,k2,k3)

  end do
  pot_ewk = 0.5d0 * energy

  vir_sca = 0.5d0 * vir_sca
  vir_tmp(2,1)=vir_tmp(1,2)
  vir_tmp(3,1)=vir_tmp(1,3)
  vir_tmp(3,2)=vir_tmp(2,3)

  vir(1:3,1:3) = vir(1:3,1:3) - 0.5d0*vir_tmp(1:3,1:3)

!     +     +     +     +     +     +     +     +     +

end subroutine scalar_sum

!-----------------------------------------------------------
#if defined(_FFTW3)
subroutine fft_forward(array,   &
     &                 nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3)

  use spme_global

  implicit none

  include 'fftw3.f'

  real(8),intent(out):: array(2,nfftdim1,nfftdim2,nfftdim3)
  integer,intent(in):: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

  integer:: dim1, dim2, dim3

!     +     +     +     +     +     +     +     +     +

!
! 1) X-direction FFT
!
  do dim3=1,nfft3
     do dim2=1,nfft2
        do dim1=1,nfft1

           fftWork1(dim1) = CMPLX(array(1,dim1,dim2,dim3),   &
                &                 array(2,dim1,dim2,dim3),kind=8)

        enddo

        call DFFTW_EXECUTE(planF(1))

        do dim1=1,nfft1
           array(1,dim1,dim2,dim3) = DBLE(fftWork1(dim1))
           array(2,dim1,dim2,dim3) = AIMAG(fftWork1(dim1))
        enddo

     enddo
  enddo

!
! 2) Y-direction FFT
!
  do dim1=1,nfft1
     do dim3=1,nfft3
        do dim2=1,nfft2

           fftWork1(dim2) = CMPLX(array(1,dim1,dim2,dim3),   &
                &                 array(2,dim1,dim2,dim3),kind=8)

        enddo

        call DFFTW_EXECUTE(planF(2))

        do dim2=1,nfft2
           array(1,dim1,dim2,dim3) = DBLE(fftWork1(dim2))
           array(2,dim1,dim2,dim3) = AIMAG(fftWork1(dim2))
        enddo

     enddo
  enddo

!
! 3) Z-direction FFT
!
  do dim2=1,nfft2
     do dim1=1,nfft1
        do dim3=1,nfft3

           fftWork1(dim3) = CMPLX(array(1,dim1,dim2,dim3),   &
                &                 array(2,dim1,dim2,dim3),kind=8)

        enddo

        call DFFTW_EXECUTE(planF(3))

        do dim3=1,nfft3
           array(1,dim1,dim2,dim3) = DBLE(fftWork1(dim3))
           array(2,dim1,dim2,dim3) = AIMAG(fftWork1(dim3))
        enddo

     enddo
  enddo

!     +     +     +     +     +     +     +     +     +

end subroutine fft_forward

#else
!-----------------------------------------------------------
subroutine fft_forward(array,   &
     &                 nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,   &
     &                 nfftable,nffwork)

  use spme_global

  implicit none


  real(8),intent(out):: array(:)
  integer,intent(in):: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
  integer,intent(out):: nfftable,nffwork

  integer:: isign

!     +     +     +     +     +     +     +     +     +

! implement backward fft3d
  isign = 1

  call pubz3d(isign,nfft1,nfft2,nfft3,array,   &
       &      nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)

end subroutine fft_forward
#endif

!----------------------------------------------------
subroutine pubz3d(isign,n1,n2,n3,w,ld1,ld2,table,ntable,   &
     &            work,nwork)

  implicit none

  integer,intent(in):: n1,n2,n3,ld1,ld2,isign,ntable,nwork
  complex(8),intent(inout):: w(ld1,ld2,n3)
  complex(8),intent(inout):: work(nwork)
  real(8),intent(in):: table(ntable,3)

  integer:: i,j,k
! ntable should be 4*max(n1,n2,n3) +15
! nwork should be max(n1,n2,n3)

!
! transform along X  first ...
!
  do k = 1, n3
     do j = 1, n2
        do i = 1,n1
           work(i) = w(i,j,k)
        end do
        if ( isign == -1) call cfftf(n1,work,table(1,1))
        if ( isign == 1) call cfftb(n1,work,table(1,1))
        do i = 1,n1
           w(i,j,k) = work(i)
        end do
     end do
  end do
!
! transform along Y then ...
!
  do k = 1,n3
     do i = 1,n1
        do j = 1,n2
           work(j) = w(i,j,k)
        end do
        if ( isign == -1) call cfftf(n2,work,table(1,2))
        if ( isign == 1) call cfftb(n2,work,table(1,2))
        do j = 1,n2
           w(i,j,k) = work(j)
        end do
     end do
  end do
!
! transform along Z finally ...
!
  do i = 1, n1
     do j = 1, n2
        do k = 1,n3
           work(k) = w(i,j,k)
        end do
        if ( isign == -1) call cfftf(n3,work,table(1,3))
        if ( isign == 1) call cfftb(n3,work,table(1,3))
        do k = 1,n3
           w(i,j,k) = work(k)
        end do
     end do
  end do

end subroutine pubz3d

!-----------------------------------------------------------
subroutine grad_sum(recip,   &
     &              theta1,theta2,theta3,   &
     &              dtheta1,dtheta2,dtheta3,   &
     &              force,   &
     &              fr1,fr2,fr3,   &
     &              pme_order,nfft1,nfft2,nfft3,   &
     &              nfftdim1,nfftdim2,nfftdim3,   &
     &              Q,   &
     &              vir_for)

  use md_global

  implicit none

  integer,intent(in):: pme_order,nfft1,nfft2,nfft3
  integer,intent(in):: nfftdim1,nfftdim2,nfftdim3
  real(8),intent(in):: recip(:)
  real(8),intent(in):: fr1(:),fr2(:),fr3(:)

  real(8),intent(inout):: force(:,:)

  real(8),intent(in):: theta1(pme_order,natom),theta2(pme_order,natom),   &
       &               theta3(pme_order,natom)
  real(8),intent(in):: dtheta1(pme_order,natom),dtheta2(pme_order,natom),   &
       &               dtheta3(pme_order,natom)
  real(8),intent(in):: Q(2,nfftdim1,nfftdim2,nfftdim3)

  real(8),intent(inout):: vir_for(:,:)

  integer:: n,ith1,ith2,ith3,i0,j0,k0,i,j,k
  real(8):: f_tmp(3),term

!     +     +     +     +     +     +     +     +     +

  do n = 1,natom
     f_tmp(1:3) = 0.0d0
     k0 = int(fr3(n)) - pme_order
     do ith3 = 1,pme_order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        j0 = int(fr2(n)) - pme_order
        do ith2 = 1,pme_order
           j0 = j0 + 1
           j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
           i0 = int(fr1(n)) - pme_order
           do ith1 = 1,pme_order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
              term = atmchrg(n)*Q(1,i,j,k)
              ! force is negative of grad
              f_tmp(1) = f_tmp(1) - nfft1 * term * dtheta1(ith1,n)   &
                     & * theta2(ith2,n) * theta3(ith3,n)
              f_tmp(2) = f_tmp(2) - nfft2 * term * theta1(ith1,n)   &
                     & * dtheta2(ith2,n) * theta3(ith3,n)
              f_tmp(3) = f_tmp(3) - nfft3 * term * theta1(ith1,n)   &
                     & * theta2(ith2,n) * dtheta3(ith3,n)
           end do
        end do
     end do
     force(1:3,n) = force(1:3,n) + recip(1:3)*f_tmp(1:3)
     vir_for(1:3,n) = vir_for(1:3,n) + recip(1:3)*f_tmp(1:3)
  end do

!     +     +     +     +     +     +     +     +     +

end subroutine grad_sum
