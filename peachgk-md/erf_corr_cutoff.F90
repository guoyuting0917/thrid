!****************************************
!*  erf_corr_cutoff.f Ver.1.0 '09.12.01 *
!*      for peachgk_md.f                *
!*            by G.Kikugawa             *
!****************************************

!---- original code Copyright
!*-*
!*-* Copyright (C) 2000 Massimo Marchi and Piero Procacci
!*-* Full copyright notice at http://www.chim.unifi.it/orac/copyright4.0.html
!*-* Contact for info M. Marchi, CEA,  Gif Sur Yvette 91191 (FRANCE) 
!*-* Email:marchi@villon.saclay.cea.fr
!*-* 
!************************************************************************
!*   Time-stamp: <98/02/10 12:14:48 marchi>                             *
!*                                                                      *
!*                                                                      *
!*                                                                      *
!*======================================================================*
!*                                                                      *
!*              Author:  Massimo Marchi                                 *
!*              CEA/Centre d'Etudes Saclay, FRANCE                      *
!*                                                                      *
!*              - Tue Feb 10 1998 -                                     *
!*                                                                      *
!************************************************************************

subroutine erf_corr_cutoff( xcel, ycel, zcel, &
     &                      nfft1, nfft2, nfft3 )

  use spme_global

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: xcel             ! initial x0 cell length[m]
  real(8),intent(in):: ycel             ! initial y0 cell length[m]
  real(8),intent(in):: zcel             ! initial z0 cell length[m]
  integer,intent(in):: nfft1,nfft2,nfft3

! LOCAL:
  real(8):: aux,mhat1,mhat2,mhat3

!     +     +     +     +     +     +     +     +     +

!--find out kcut-off in PME
!--To find out Kmax must divide by 2 nfftx
  aux = 0.5d0
  mhat1 = aux*(nfft1/xcel)
  rkcut = mhat1
  mhat2 = aux*(nfft2/ycel)
  if(mhat2 < mhat1) rkcut = mhat2
  mhat3 = aux*(nfft3/zcel)
  if(mhat3 < mhat2) rkcut = mhat3
!--decrease cut-off to avoid error at the boundary 
  rkcut = 0.9d0*rkcut 

!     +     +     +     +     +     +     +     +     +

  return
end subroutine erf_corr_cutoff
