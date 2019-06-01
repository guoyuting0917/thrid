**********************************************
*  createcor_.f Ver.2.5magnetite  '08.01.30  *
*      for peachgk_md.f                      *
*            by G.Kikugawa                   *
**********************************************
      subroutine createcor(npoly,npolytyp,npoly_mole,npoly_atom,
     &                     nwater,nmatom,nmatyp,nmatomtyp,
     &                     xmaxpo,ymaxpo,zmaxpo,
     &                     xmaxw,ymaxw,zmaxw,
     &                     xmaxma,ymaxma,zmaxma,
     &                     inicorpo,inicorw,inicorma,
     &                     ncrecorpo,index_crecorpo,
     &                     ncrecorw,index_crecorw,
     &                     ncrecorma,index_crecorma,
     &                     ifsetcor,
     &                     xcel,ycel,zcel,
     &                     xref,
     &                     mt,mti,
     &                     compfact,
     &                     mchain)

C#################################################################
C# Usage of the inicorpo(), inicorw(), inicorma()                #
C#   inicorpo(), inicorw(), inicorma() are the real*8 type.      #
C#   Therefore, you have to cast the type which you want to use. #
C#   For example, int(inicorpo()) to cast the int type           #
C#################################################################

      implicit none

      include 'md_common.h'
#if defined(MPI)
      include 'mpi_common.h'
#endif

c ARGUMENT:
c     INPUT
      integer:: npoly           ! all number of poly
      integer:: npolytyp        ! number of poly type
      integer:: npoly_mole(*)   ! number of molecules of each poly
      integer:: npoly_atom(*)   ! number of atoms belonging to poly

      integer:: nwater          ! number of H2O molecules

      integer:: nmatom          ! number of monatomic molecules
      integer:: nmatyp          ! number of species of monatomic mole.
      integer:: nmatomtyp(*)    ! each number of monatomic mole.

      integer:: xmaxpo(*)       ! use for positioning of polymer1
      integer:: ymaxpo(*)       ! use for positioning of polymer1
      integer:: zmaxpo(*)       ! use for positioning of polymer1
      integer:: xmaxw           ! use for positioning of water
      integer:: ymaxw           ! use for positioning of water
      integer:: zmaxw           ! use for positioning of water
      integer:: xmaxma(*)       ! use for positioning of monatomic mole.
      integer:: ymaxma(*)       ! use for positioning of monatomic mole.
      integer:: zmaxma(*)       ! use for positioning of monatomic mole.
      real*8:: inicorpo(maxnpolytyp,*) ! use for positioning of poly
      real*8:: inicorw(*)       ! use for positioning of water
      real*8:: inicorma(maxnmatyp,*) ! use for positioning of ma
      integer:: ncrecorpo       ! max number of poly type for createcor
      integer:: index_crecorpo(*) ! index of polymer for createcor
      integer:: ncrecorw        ! water for createcor
      integer:: index_crecorw   ! index of water type for createcor
      integer:: ncrecorma       ! max number of matom type for createcor
      integer:: index_crecorma(*) ! index of matom type for createcor

      logical:: ifsetcor(*)     ! frag for checking if coordinate has set

      real*8:: xref             ! distanse base value [m]

      real*8:: xcel             ! x cell length
      real*8:: ycel             ! y cell length
      real*8:: zcel             ! z cell length
      
!      real*8:: randseed         ! random seed for createcor etc.
      integer:: nperiod         ! period parameter
      parameter(nperiod = 624)
      integer:: mt(0:nperiod-1)
      integer:: mti

      real*8:: compfact         ! compact factor using at poly arrange(<1.0)

      integer:: mchain          ! Nose-Hoover chain number

c LOCAL:
      real*8:: rmin = 2.0d-10   ! minimum of atom distance 
      integer:: maxicount = 1000000 ! max count for arrange water
      real*8:: dmin = 0.05d-10  ! difference of rmin

      integer:: i
      integer:: j,jj            ! do loop index
      integer:: m,n             ! do loop index
      integer:: n1,n2,nn
      integer:: j1,j2,k,kk
      integer:: ix,iy,iz
      integer:: imole
      integer:: icount

      real*8:: rot_m(3,3)       ! rotational matrix
      real*8:: psi,phi,theta    ! Euler angle
      real*8:: rxp,ryp,rzp      ! coordinate of poly1
      real*8:: rxw,ryw,rzw      ! coordinate of water
!      real*8:: rxma,ryma,rzma   ! coordinate of monatomic mole.

!      real*8:: ixx              ! random seed

      real*8:: pi               ! = 3.14159...

      real*8:: atmcor_t(3)      ! atmcor tmp

      real*8:: patmcor(3,maxnatom) ! tmp atmcor for poly1
      real*8:: watmcor(3,3)     ! tmp atmcor for water
      real*8:: maatmcor(3)      ! tmp atmcor for MA

      real*8:: r_ij(3)          ! xi - xj etc.
      real*8:: rij              ! distance of r_ij
      real*8:: r_min            ! distance of r_min

      integer:: ipoly
      integer:: imatom

      integer:: icrepo, icrewa, icrema

c     parameters for spinel crystal lattice of magnetite
      real*8:: A0, A1
      real*8:: snx,sny,snz

      real*8:: lower_origin     ! origin for lower magnetite crystal layer
      real*8:: upper_origin     ! origin for upper magnetite crystal layer

      real*8:: lower_water      ! lower z-origin for water layer
      real*8:: water_layer      ! width of water layer

      real*8:: base_latt_O(3,64)   ! Oxide ion
      real*8:: base_latt_FeO(3,32) ! Ferric ion 2.5+ (Octohedral site)
      real*8:: base_latt_FeT(3,16)  ! Ferric ion 3+ (Tetrahedral site)

      real*8:: baseloop

      integer:: index_magx,index_magy,index_magz

c FUNCTIONS:
!      real*8:: ranq            ! function creating random number
!      external ranq             ! random number (0,1)
      real*8:: grnd             ! function creating random number
      
c     +     +     +     +     +     +     +

c---- notation ----
#ifdef MPI
      if (irank == 0) then
#endif

      write(6,*) 'Create atmcor from createcor'

#ifdef MPI
      end if
#endif

      pi = dacos(-1.0d0)

      A0 = 8.6600d-10 / xref    ! lattice constant
!      A0 = 8.3925d-10 / xref    ! lattice constant at 299K
                                ! in Okudera, Kihara, Matsumoto,
                                !     Acta Cryst., B52 (1996), pp. 450-457.
      A1 = dsqrt(2.0d0) * A0

      lower_origin = 0.0d0
      lower_water = lower_origin + A0*DBLE(zmaxma(1)/4 - 1)
      water_layer = 50.0d-10 / xref
      upper_origin = lower_water + water_layer

c---- create base lattice for magnetite (full layer)
!     1st layer
      base_latt_O(1, 1) = 1.0d0/8.0d0 * A1
      base_latt_O(2, 1) = 1.0d0/8.0d0 * A1
      base_latt_O(3, 1) = 0.0d0       * A0
      base_latt_O(1, 2) = 3.0d0/8.0d0 * A1
      base_latt_O(2, 2) = 1.0d0/8.0d0 * A1
      base_latt_O(3, 2) = 0.0d0       * A0
      base_latt_O(1, 3) = 5.0d0/8.0d0 * A1
      base_latt_O(2, 3) = 1.0d0/8.0d0 * A1
      base_latt_O(3, 3) = 0.0d0       * A0
      base_latt_O(1, 4) = 7.0d0/8.0d0 * A1
      base_latt_O(2, 4) = 1.0d0/8.0d0 * A1
      base_latt_O(3, 4) = 0.0d0       * A0
      base_latt_O(1, 5) = 1.0d0/8.0d0 * A1
      base_latt_O(2, 5) = 3.0d0/8.0d0 * A1
      base_latt_O(3, 5) = 0.0d0       * A0
      base_latt_O(1, 6) = 3.0d0/8.0d0 * A1
      base_latt_O(2, 6) = 3.0d0/8.0d0 * A1
      base_latt_O(3, 6) = 0.0d0       * A0
      base_latt_O(1, 7) = 5.0d0/8.0d0 * A1
      base_latt_O(2, 7) = 3.0d0/8.0d0 * A1
      base_latt_O(3, 7) = 0.0d0       * A0
      base_latt_O(1, 8) = 7.0d0/8.0d0 * A1
      base_latt_O(2, 8) = 3.0d0/8.0d0 * A1
      base_latt_O(3, 8) = 0.0d0       * A0
      base_latt_O(1, 9) = 1.0d0/8.0d0 * A1
      base_latt_O(2, 9) = 5.0d0/8.0d0 * A1
      base_latt_O(3, 9) = 0.0d0       * A0
      base_latt_O(1,10) = 3.0d0/8.0d0 * A1
      base_latt_O(2,10) = 5.0d0/8.0d0 * A1
      base_latt_O(3,10) = 0.0d0       * A0
      base_latt_O(1,11) = 5.0d0/8.0d0 * A1
      base_latt_O(2,11) = 5.0d0/8.0d0 * A1
      base_latt_O(3,11) = 0.0d0       * A0
      base_latt_O(1,12) = 7.0d0/8.0d0 * A1
      base_latt_O(2,12) = 5.0d0/8.0d0 * A1
      base_latt_O(3,12) = 0.0d0       * A0
      base_latt_O(1,13) = 1.0d0/8.0d0 * A1
      base_latt_O(2,13) = 7.0d0/8.0d0 * A1
      base_latt_O(3,13) = 0.0d0       * A0
      base_latt_O(1,14) = 3.0d0/8.0d0 * A1
      base_latt_O(2,14) = 7.0d0/8.0d0 * A1
      base_latt_O(3,14) = 0.0d0       * A0
      base_latt_O(1,15) = 5.0d0/8.0d0 * A1
      base_latt_O(2,15) = 7.0d0/8.0d0 * A1
      base_latt_O(3,15) = 0.0d0       * A0
      base_latt_O(1,16) = 7.0d0/8.0d0 * A1
      base_latt_O(2,16) = 7.0d0/8.0d0 * A1
      base_latt_O(3,16) = 0.0d0       * A0

      base_latt_FeO(1,1) = 2.0d0/8.0d0 * A1
      base_latt_FeO(2,1) = 0.0d0       * A1
      base_latt_FeO(3,1) = 0.0d0       * A0
      base_latt_FeO(1,2) = 6.0d0/8.0d0 * A1
      base_latt_FeO(2,2) = 0.0d0       * A1
      base_latt_FeO(3,2) = 0.0d0       * A0
      base_latt_FeO(1,3) = 2.0d0/8.0d0 * A1
      base_latt_FeO(2,3) = 2.0d0/8.0d0 * A1
      base_latt_FeO(3,3) = 0.0d0       * A0
      base_latt_FeO(1,4) = 6.0d0/8.0d0 * A1
      base_latt_FeO(2,4) = 2.0d0/8.0d0 * A1
      base_latt_FeO(3,4) = 0.0d0       * A0
      base_latt_FeO(1,5) = 2.0d0/8.0d0 * A1
      base_latt_FeO(2,5) = 4.0d0/8.0d0 * A1
      base_latt_FeO(3,5) = 0.0d0       * A0
      base_latt_FeO(1,6) = 6.0d0/8.0d0 * A1
      base_latt_FeO(2,6) = 4.0d0/8.0d0 * A1
      base_latt_FeO(3,6) = 0.0d0       * A0
      base_latt_FeO(1,7) = 2.0d0/8.0d0 * A1
      base_latt_FeO(2,7) = 6.0d0/8.0d0 * A1
      base_latt_FeO(3,7) = 0.0d0       * A0
      base_latt_FeO(1,8) = 6.0d0/8.0d0 * A1
      base_latt_FeO(2,8) = 6.0d0/8.0d0 * A1
      base_latt_FeO(3,8) = 0.0d0       * A0

!     2nd layer
      base_latt_FeT(1,1) = 0.0d0       * A1
      base_latt_FeT(2,1) = 1.0d0/8.0d0 * A1
      base_latt_FeT(3,1) = 1.0d0/8.0d0 * A0
      base_latt_FeT(1,2) = 4.0d0/8.0d0 * A1
      base_latt_FeT(2,2) = 1.0d0/8.0d0 * A1
      base_latt_FeT(3,2) = 1.0d0/8.0d0 * A0
      base_latt_FeT(1,3) = 0.0d0       * A1
      base_latt_FeT(2,3) = 5.0d0/8.0d0 * A1
      base_latt_FeT(3,3) = 1.0d0/8.0d0 * A0
      base_latt_FeT(1,4) = 4.0d0/8.0d0 * A1
      base_latt_FeT(2,4) = 5.0d0/8.0d0 * A1
      base_latt_FeT(3,4) = 1.0d0/8.0d0 * A0

!     3rd layer
      base_latt_O(1,17) = 0.0d0       * A1
      base_latt_O(2,17) = 0.0d0       * A1
      base_latt_O(3,17) = 2.0d0/8.0d0 * A0
      base_latt_O(1,18) = 2.0d0/8.0d0 * A1
      base_latt_O(2,18) = 0.0d0       * A1
      base_latt_O(3,18) = 2.0d0/8.0d0 * A0
      base_latt_O(1,19) = 4.0d0/8.0d0 * A1
      base_latt_O(2,19) = 0.0d0       * A1
      base_latt_O(3,19) = 2.0d0/8.0d0 * A0
      base_latt_O(1,20) = 6.0d0/8.0d0 * A1
      base_latt_O(2,20) = 0.0d0       * A1
      base_latt_O(3,20) = 2.0d0/8.0d0 * A0
      base_latt_O(1,21) = 0.0d0       * A1
      base_latt_O(2,21) = 2.0d0/8.0d0 * A1
      base_latt_O(3,21) = 2.0d0/8.0d0 * A0
      base_latt_O(1,22) = 2.0d0/8.0d0 * A1
      base_latt_O(2,22) = 2.0d0/8.0d0 * A1
      base_latt_O(3,22) = 2.0d0/8.0d0 * A0
      base_latt_O(1,23) = 4.0d0/8.0d0 * A1
      base_latt_O(2,23) = 2.0d0/8.0d0 * A1
      base_latt_O(3,23) = 2.0d0/8.0d0 * A0
      base_latt_O(1,24) = 6.0d0/8.0d0 * A1
      base_latt_O(2,24) = 2.0d0/8.0d0 * A1
      base_latt_O(3,24) = 2.0d0/8.0d0 * A0
      base_latt_O(1,25) = 0.0d0       * A1
      base_latt_O(2,25) = 4.0d0/8.0d0 * A1
      base_latt_O(3,25) = 2.0d0/8.0d0 * A0
      base_latt_O(1,26) = 2.0d0/8.0d0 * A1
      base_latt_O(2,26) = 4.0d0/8.0d0 * A1
      base_latt_O(3,26) = 2.0d0/8.0d0 * A0
      base_latt_O(1,27) = 4.0d0/8.0d0 * A1
      base_latt_O(2,27) = 4.0d0/8.0d0 * A1
      base_latt_O(3,27) = 2.0d0/8.0d0 * A0
      base_latt_O(1,28) = 6.0d0/8.0d0 * A1
      base_latt_O(2,28) = 4.0d0/8.0d0 * A1
      base_latt_O(3,28) = 2.0d0/8.0d0 * A0
      base_latt_O(1,29) = 0.0d0       * A1
      base_latt_O(2,29) = 6.0d0/8.0d0 * A1
      base_latt_O(3,29) = 2.0d0/8.0d0 * A0
      base_latt_O(1,30) = 2.0d0/8.0d0 * A1
      base_latt_O(2,30) = 6.0d0/8.0d0 * A1
      base_latt_O(3,30) = 2.0d0/8.0d0 * A0
      base_latt_O(1,31) = 4.0d0/8.0d0 * A1
      base_latt_O(2,31) = 6.0d0/8.0d0 * A1
      base_latt_O(3,31) = 2.0d0/8.0d0 * A0
      base_latt_O(1,32) = 6.0d0/8.0d0 * A1
      base_latt_O(2,32) = 6.0d0/8.0d0 * A1
      base_latt_O(3,32) = 2.0d0/8.0d0 * A0

      base_latt_FeO(1, 9) = 1.0d0/8.0d0 * A1
      base_latt_FeO(2, 9) = 3.0d0/8.0d0 * A1
      base_latt_FeO(3, 9) = 2.0d0/8.0d0 * A0
      base_latt_FeO(1,10) = 3.0d0/8.0d0 * A1
      base_latt_FeO(2,10) = 3.0d0/8.0d0 * A1
      base_latt_FeO(3,10) = 2.0d0/8.0d0 * A0
      base_latt_FeO(1,11) = 5.0d0/8.0d0 * A1
      base_latt_FeO(2,11) = 3.0d0/8.0d0 * A1
      base_latt_FeO(3,11) = 2.0d0/8.0d0 * A0
      base_latt_FeO(1,12) = 7.0d0/8.0d0 * A1
      base_latt_FeO(2,12) = 3.0d0/8.0d0 * A1
      base_latt_FeO(3,12) = 2.0d0/8.0d0 * A0
      base_latt_FeO(1,13) = 1.0d0/8.0d0 * A1
      base_latt_FeO(2,13) = 7.0d0/8.0d0 * A1
      base_latt_FeO(3,13) = 2.0d0/8.0d0 * A0
      base_latt_FeO(1,14) = 3.0d0/8.0d0 * A1
      base_latt_FeO(2,14) = 7.0d0/8.0d0 * A1
      base_latt_FeO(3,14) = 2.0d0/8.0d0 * A0
      base_latt_FeO(1,15) = 5.0d0/8.0d0 * A1
      base_latt_FeO(2,15) = 7.0d0/8.0d0 * A1
      base_latt_FeO(3,15) = 2.0d0/8.0d0 * A0
      base_latt_FeO(1,16) = 7.0d0/8.0d0 * A1
      base_latt_FeO(2,16) = 7.0d0/8.0d0 * A1
      base_latt_FeO(3,16) = 2.0d0/8.0d0 * A0

!     4th layer
      base_latt_FeT(1,5) = 2.0d0/8.0d0 * A1
      base_latt_FeT(2,5) = 1.0d0/8.0d0 * A1
      base_latt_FeT(3,5) = 3.0d0/8.0d0 * A0
      base_latt_FeT(1,6) = 6.0d0/8.0d0 * A1
      base_latt_FeT(2,6) = 1.0d0/8.0d0 * A1
      base_latt_FeT(3,6) = 3.0d0/8.0d0 * A0
      base_latt_FeT(1,7) = 2.0d0/8.0d0 * A1
      base_latt_FeT(2,7) = 5.0d0/8.0d0 * A1
      base_latt_FeT(3,7) = 3.0d0/8.0d0 * A0
      base_latt_FeT(1,8) = 6.0d0/8.0d0 * A1
      base_latt_FeT(2,8) = 5.0d0/8.0d0 * A1
      base_latt_FeT(3,8) = 3.0d0/8.0d0 * A0

!     5th layer
      base_latt_O(1,33) = 1.0d0/8.0d0 * A1
      base_latt_O(2,33) = 1.0d0/8.0d0 * A1
      base_latt_O(3,33) = 4.0d0/8.0d0 * A0
      base_latt_O(1,34) = 3.0d0/8.0d0 * A1
      base_latt_O(2,34) = 1.0d0/8.0d0 * A1
      base_latt_O(3,34) = 4.0d0/8.0d0 * A0
      base_latt_O(1,35) = 5.0d0/8.0d0 * A1
      base_latt_O(2,35) = 1.0d0/8.0d0 * A1
      base_latt_O(3,35) = 4.0d0/8.0d0 * A0
      base_latt_O(1,36) = 7.0d0/8.0d0 * A1
      base_latt_O(2,36) = 1.0d0/8.0d0 * A1
      base_latt_O(3,36) = 4.0d0/8.0d0 * A0
      base_latt_O(1,37) = 1.0d0/8.0d0 * A1
      base_latt_O(2,37) = 3.0d0/8.0d0 * A1
      base_latt_O(3,37) = 4.0d0/8.0d0 * A0
      base_latt_O(1,38) = 3.0d0/8.0d0 * A1
      base_latt_O(2,38) = 3.0d0/8.0d0 * A1
      base_latt_O(3,38) = 4.0d0/8.0d0 * A0
      base_latt_O(1,39) = 5.0d0/8.0d0 * A1
      base_latt_O(2,39) = 3.0d0/8.0d0 * A1
      base_latt_O(3,39) = 4.0d0/8.0d0 * A0
      base_latt_O(1,40) = 7.0d0/8.0d0 * A1
      base_latt_O(2,40) = 3.0d0/8.0d0 * A1
      base_latt_O(3,40) = 4.0d0/8.0d0 * A0
      base_latt_O(1,41) = 1.0d0/8.0d0 * A1
      base_latt_O(2,41) = 5.0d0/8.0d0 * A1
      base_latt_O(3,41) = 4.0d0/8.0d0 * A0
      base_latt_O(1,42) = 3.0d0/8.0d0 * A1
      base_latt_O(2,42) = 5.0d0/8.0d0 * A1
      base_latt_O(3,42) = 4.0d0/8.0d0 * A0
      base_latt_O(1,43) = 5.0d0/8.0d0 * A1
      base_latt_O(2,43) = 5.0d0/8.0d0 * A1
      base_latt_O(3,43) = 4.0d0/8.0d0 * A0
      base_latt_O(1,44) = 7.0d0/8.0d0 * A1
      base_latt_O(2,44) = 5.0d0/8.0d0 * A1
      base_latt_O(3,44) = 4.0d0/8.0d0 * A0
      base_latt_O(1,45) = 1.0d0/8.0d0 * A1
      base_latt_O(2,45) = 7.0d0/8.0d0 * A1
      base_latt_O(3,45) = 4.0d0/8.0d0 * A0
      base_latt_O(1,46) = 3.0d0/8.0d0 * A1
      base_latt_O(2,46) = 7.0d0/8.0d0 * A1
      base_latt_O(3,46) = 4.0d0/8.0d0 * A0
      base_latt_O(1,47) = 5.0d0/8.0d0 * A1
      base_latt_O(2,47) = 7.0d0/8.0d0 * A1
      base_latt_O(3,47) = 4.0d0/8.0d0 * A0
      base_latt_O(1,48) = 7.0d0/8.0d0 * A1
      base_latt_O(2,48) = 7.0d0/8.0d0 * A1
      base_latt_O(3,48) = 4.0d0/8.0d0 * A0

      base_latt_FeO(1,17) = 0.0d0       * A1
      base_latt_FeO(2,17) = 0.0d0       * A1
      base_latt_FeO(3,17) = 4.0d0/8.0d0 * A0
      base_latt_FeO(1,18) = 4.0d0/8.0d0 * A1
      base_latt_FeO(2,18) = 0.0d0       * A1
      base_latt_FeO(3,18) = 4.0d0/8.0d0 * A0
      base_latt_FeO(1,19) = 0.0d0       * A1
      base_latt_FeO(2,19) = 2.0d0/8.0d0 * A1
      base_latt_FeO(3,19) = 4.0d0/8.0d0 * A0
      base_latt_FeO(1,20) = 4.0d0/8.0d0 * A1
      base_latt_FeO(2,20) = 2.0d0/8.0d0 * A1
      base_latt_FeO(3,20) = 4.0d0/8.0d0 * A0
      base_latt_FeO(1,21) = 0.0d0       * A1
      base_latt_FeO(2,21) = 4.0d0/8.0d0 * A1
      base_latt_FeO(3,21) = 4.0d0/8.0d0 * A0
      base_latt_FeO(1,22) = 4.0d0/8.0d0 * A1
      base_latt_FeO(2,22) = 4.0d0/8.0d0 * A1
      base_latt_FeO(3,22) = 4.0d0/8.0d0 * A0
      base_latt_FeO(1,23) = 0.0d0       * A1
      base_latt_FeO(2,23) = 6.0d0/8.0d0 * A1
      base_latt_FeO(3,23) = 4.0d0/8.0d0 * A0
      base_latt_FeO(1,24) = 4.0d0/8.0d0 * A1
      base_latt_FeO(2,24) = 6.0d0/8.0d0 * A1
      base_latt_FeO(3,24) = 4.0d0/8.0d0 * A0

!     6th layer
      base_latt_FeT(1, 9) = 2.0d0/8.0d0 * A1
      base_latt_FeT(2, 9) = 3.0d0/8.0d0 * A1
      base_latt_FeT(3, 9) = 5.0d0/8.0d0 * A0
      base_latt_FeT(1,10) = 6.0d0/8.0d0 * A1
      base_latt_FeT(2,10) = 3.0d0/8.0d0 * A1
      base_latt_FeT(3,10) = 5.0d0/8.0d0 * A0
      base_latt_FeT(1,11) = 2.0d0/8.0d0 * A1
      base_latt_FeT(2,11) = 7.0d0/8.0d0 * A1
      base_latt_FeT(3,11) = 5.0d0/8.0d0 * A0
      base_latt_FeT(1,12) = 6.0d0/8.0d0 * A1
      base_latt_FeT(2,12) = 7.0d0/8.0d0 * A1
      base_latt_FeT(3,12) = 5.0d0/8.0d0 * A0

!     7th layer
      base_latt_O(1,49) = 0.0d0       * A1
      base_latt_O(2,49) = 0.0d0       * A1
      base_latt_O(3,49) = 6.0d0/8.0d0 * A0
      base_latt_O(1,50) = 2.0d0/8.0d0 * A1
      base_latt_O(2,50) = 0.0d0       * A1
      base_latt_O(3,50) = 6.0d0/8.0d0 * A0
      base_latt_O(1,51) = 4.0d0/8.0d0 * A1
      base_latt_O(2,51) = 0.0d0       * A1
      base_latt_O(3,51) = 6.0d0/8.0d0 * A0
      base_latt_O(1,52) = 6.0d0/8.0d0 * A1
      base_latt_O(2,52) = 0.0d0       * A1
      base_latt_O(3,52) = 6.0d0/8.0d0 * A0
      base_latt_O(1,53) = 0.0d0       * A1
      base_latt_O(2,53) = 2.0d0/8.0d0 * A1
      base_latt_O(3,53) = 6.0d0/8.0d0 * A0
      base_latt_O(1,54) = 2.0d0/8.0d0 * A1
      base_latt_O(2,54) = 2.0d0/8.0d0 * A1
      base_latt_O(3,54) = 6.0d0/8.0d0 * A0
      base_latt_O(1,55) = 4.0d0/8.0d0 * A1
      base_latt_O(2,55) = 2.0d0/8.0d0 * A1
      base_latt_O(3,55) = 6.0d0/8.0d0 * A0
      base_latt_O(1,56) = 6.0d0/8.0d0 * A1
      base_latt_O(2,56) = 2.0d0/8.0d0 * A1
      base_latt_O(3,56) = 6.0d0/8.0d0 * A0
      base_latt_O(1,57) = 0.0d0       * A1
      base_latt_O(2,57) = 4.0d0/8.0d0 * A1
      base_latt_O(3,57) = 6.0d0/8.0d0 * A0
      base_latt_O(1,58) = 2.0d0/8.0d0 * A1
      base_latt_O(2,58) = 4.0d0/8.0d0 * A1
      base_latt_O(3,58) = 6.0d0/8.0d0 * A0
      base_latt_O(1,59) = 4.0d0/8.0d0 * A1
      base_latt_O(2,59) = 4.0d0/8.0d0 * A1
      base_latt_O(3,59) = 6.0d0/8.0d0 * A0
      base_latt_O(1,60) = 6.0d0/8.0d0 * A1
      base_latt_O(2,60) = 4.0d0/8.0d0 * A1
      base_latt_O(3,60) = 6.0d0/8.0d0 * A0
      base_latt_O(1,61) = 0.0d0       * A1
      base_latt_O(2,61) = 6.0d0/8.0d0 * A1
      base_latt_O(3,61) = 6.0d0/8.0d0 * A0
      base_latt_O(1,62) = 2.0d0/8.0d0 * A1
      base_latt_O(2,62) = 6.0d0/8.0d0 * A1
      base_latt_O(3,62) = 6.0d0/8.0d0 * A0
      base_latt_O(1,63) = 4.0d0/8.0d0 * A1
      base_latt_O(2,63) = 6.0d0/8.0d0 * A1
      base_latt_O(3,63) = 6.0d0/8.0d0 * A0
      base_latt_O(1,64) = 6.0d0/8.0d0 * A1
      base_latt_O(2,64) = 6.0d0/8.0d0 * A1
      base_latt_O(3,64) = 6.0d0/8.0d0 * A0

      base_latt_FeO(1,25) = 1.0d0/8.0d0 * A1
      base_latt_FeO(2,25) = 1.0d0/8.0d0 * A1
      base_latt_FeO(3,25) = 6.0d0/8.0d0 * A0
      base_latt_FeO(1,26) = 3.0d0/8.0d0 * A1
      base_latt_FeO(2,26) = 1.0d0/8.0d0 * A1
      base_latt_FeO(3,26) = 6.0d0/8.0d0 * A0
      base_latt_FeO(1,27) = 5.0d0/8.0d0 * A1
      base_latt_FeO(2,27) = 1.0d0/8.0d0 * A1
      base_latt_FeO(3,27) = 6.0d0/8.0d0 * A0
      base_latt_FeO(1,28) = 7.0d0/8.0d0 * A1
      base_latt_FeO(2,28) = 1.0d0/8.0d0 * A1
      base_latt_FeO(3,28) = 6.0d0/8.0d0 * A0
      base_latt_FeO(1,29) = 1.0d0/8.0d0 * A1
      base_latt_FeO(2,29) = 5.0d0/8.0d0 * A1
      base_latt_FeO(3,29) = 6.0d0/8.0d0 * A0
      base_latt_FeO(1,30) = 3.0d0/8.0d0 * A1
      base_latt_FeO(2,30) = 5.0d0/8.0d0 * A1
      base_latt_FeO(3,30) = 6.0d0/8.0d0 * A0
      base_latt_FeO(1,31) = 5.0d0/8.0d0 * A1
      base_latt_FeO(2,31) = 5.0d0/8.0d0 * A1
      base_latt_FeO(3,31) = 6.0d0/8.0d0 * A0
      base_latt_FeO(1,32) = 7.0d0/8.0d0 * A1
      base_latt_FeO(2,32) = 5.0d0/8.0d0 * A1
      base_latt_FeO(3,32) = 6.0d0/8.0d0 * A0

!     8th layer
      base_latt_FeT(1,13) = 0.0d0       * A1
      base_latt_FeT(2,13) = 3.0d0/8.0d0 * A1
      base_latt_FeT(3,13) = 7.0d0/8.0d0 * A0
      base_latt_FeT(1,14) = 4.0d0/8.0d0 * A1
      base_latt_FeT(2,14) = 3.0d0/8.0d0 * A1
      base_latt_FeT(3,14) = 7.0d0/8.0d0 * A0
      base_latt_FeT(1,15) = 0.0d0       * A1
      base_latt_FeT(2,15) = 7.0d0/8.0d0 * A1
      base_latt_FeT(3,15) = 7.0d0/8.0d0 * A0
      base_latt_FeT(1,16) = 4.0d0/8.0d0 * A1
      base_latt_FeT(2,16) = 7.0d0/8.0d0 * A1
      base_latt_FeT(3,16) = 7.0d0/8.0d0 * A0

c-------- Locate molecules of monatomic mole. --------

      do icrema = 1, ncrecorma
         imatom = index_crecorma(icrema)
         
         imole = npoly + nwater
         do i = 1, imatom-1
            imole = imole + nmatomtyp(i)
         end do

         r_min = rmin / xref
#ifdef MPI
         if (irank == 0) then
#endif
         write(6,*) 'r_min= ',r_min*xref
#ifdef MPI
         end if
#endif

c----    storage of MA p-coordinate
         j1 = molept_index(imole+1) ! first monatomic mole.
         j2 = molept_index(imole+1+1) - 1
         if ((nmatom .gt. 0) .and. (j1 .ne. j2)) then
            write(6,*) 'Failure: molept_index is wrong to arrange MA'
            stop
         end if
         k = molept_list(j1)
         maatmcor(1) = atmcor(1,k)
         maatmcor(2) = atmcor(2,k)
         maatmcor(3) = atmcor(3,k)

c----    arrange MA molecules (for magnetite(001) surface)
         do iz=1,zmaxma(imatom)
            do iy=1,ymaxma(imatom)
               do ix=1,xmaxma(imatom)

                  icount = 0
                  imole = imole + 1
                  ifsetcor(imole) = .true.
                     
                  j1 = molept_index(imole)
                  j2 = molept_index(imole+1) - 1

 90               continue

                  if (icrema == 1) then ! Oxide ion (OM)

                     index_magz = mod(iz-1,4) ! 0,1,2,3,0,...
                     index_magy = mod(iy-1,4) ! 0,1,2,3,0,...
                     index_magx = mod(ix-1,4) ! 0,1,2,3,0,...
                     
                     snx = DBLE((ix-1) / 4) * A1 ! shift of unit lattice
     &                   + base_latt_O(1,1 + index_magx
     &                                     + 4*index_magy
     &                                     + 4*4*index_magz)
                     sny = DBLE((iy-1) / 4) * A1 ! shift of unit lattice
     &                   + base_latt_O(2,1 + index_magx
     &                                     + 4*index_magy
     &                                     + 4*4*index_magz)
                     snz = DBLE((iz-1) / 4) * A0 ! shift of unit lattice
     &                   + base_latt_O(3,1 + index_magx
     &                                     + 4*index_magy
     &                                     + 4*4*index_magz)

                  else if (icrema == 2) then ! Ferric ion 2.5+ (FO)

                     index_magz = mod(iz-1,4) ! 0,1,2,3,0,...
                     index_magy = mod(iy-1,4) ! 0,1,2,3,0,...
                     index_magx = mod(ix-1,2) ! 0,1,0,...

                     snx = DBLE((ix-1) / 2) * A1 ! shift of unit lattice
     &                   + base_latt_FeO(1,1 + index_magx
     &                                       + 2*index_magy
     &                                       + 2*4*index_magz)
                     sny = DBLE((iy-1) / 4) * A1 ! shift of unit lattice
     &                   + base_latt_FeO(2,1 + index_magx
     &                                       + 2*index_magy
     &                                       + 2*4*index_magz)
                     snz = DBLE((iz-1) / 4) * A0 ! shift of unit lattice
     &                   + base_latt_FeO(3,1 + index_magx
     &                                       + 2*index_magy
     &                                       + 2*4*index_magz)

                  else if (icrema == 3) then ! Ferric ion 3+ (FT)
                     
                     index_magz = mod(iz-1,4) ! 0,1,2,3,0,...
                     index_magy = mod(iy-1,2) ! 0,1,0,...
                     index_magx = mod(ix-1,2) ! 0,1,0,...

                     snx = DBLE((ix-1) / 2) * A1 ! shift of unit lattice
     &                   + base_latt_FeT(1,1 + index_magx
     &                                       + 2*index_magy
     &                                       + 2*2*index_magz)
                     sny = DBLE((iy-1) / 2) * A1 ! shift of unit lattice
     &                   + base_latt_FeT(2,1 + index_magx
     &                                       + 2*index_magy
     &                                       + 2*2*index_magz)
                     snz = DBLE((iz-1) / 4) * A0 ! shift of unit lattice
     &                   + base_latt_FeT(3,1 + index_magx
     &                                       + 2*index_magy
     &                                       + 2*2*index_magz)

                  end if

                  do j= j1, j2
                     k = molept_list(j)

                     atmcor_t(1) = atmcor(1,k)
                     atmcor_t(2) = atmcor(2,k)
                     atmcor_t(3) = atmcor(3,k)

                     atmcor(1,k) = snx
                     atmcor(2,k) = sny
                     atmcor(3,k) = snz

c---- check the atom is very close to the other atoms
                     do m = 1, nmole
                        if (.not. ifsetcor(m) .or. imole .eq. m) cycle

                        n1 = molept_index(m)
                        n2 = molept_index(m+1) - 1
                        do n = n1, n2
                           nn = molept_list(n)

                           r_ij(1) = atmcor(1,k) - atmcor(1,nn)
                           r_ij(2) = atmcor(2,k) - atmcor(2,nn)
                           r_ij(3) = atmcor(3,k) - atmcor(3,nn)
                           r_ij(1) = r_ij(1) - xcel
     &                             * dnint(r_ij(1)/xcel)
                           r_ij(2) = r_ij(2) - ycel
     &                             * dnint(r_ij(2)/ycel)
                           r_ij(3) = r_ij(3) - zcel
     &                             * dnint(r_ij(3)/zcel)

                           rij = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2
                           rij = dsqrt(rij)

                           if (rij .lt. r_min) then
                              kk = molept_list(j1)
                              atmcor(1,kk) = maatmcor(1)
                              atmcor(2,kk) = maatmcor(2)
                              atmcor(3,kk) = maatmcor(3)
                              icount = icount + 1
                              if (icount .gt. maxicount) then
                                 r_min = (r_min*xref - dmin)/xref
#ifdef MPI
                                 if (irank == 0) then
#endif
                                 write(6,*) 'r_min= ',r_min*xref
#ifdef MPI
                                 end if
#endif
                                 icount = 0
                              end if
                              goto 90
                           end if
                        end do
                     end do

#ifdef MPI
                     if (irank == 0) then
#endif
                     write(6,*) 'Set atmcor of mole No. ',imole,',',
     &                          icount
#ifdef MPI
                     end if
#endif
                  end do

               end do
            end do
         end do

      end do

c-------- Locate molecules of poly (for arranging thiol) --------

      do icrepo = 1, ncrecorpo
         ipoly = index_crecorpo(icrepo)

         imole = 0
         do i = 1, ipoly-1
            imole = imole + npoly_mole(i)
         end do

         r_min = rmin / xref
#if defined(MPI)
         if (irank == 0) then
#endif
         write(6,*) 'r_min= ',r_min*xref
#if defined(MPI)
         end if
#endif

         do iz=1,zmaxpo(ipoly)
            do iy=1,ymaxpo(ipoly)
               do ix=1,xmaxpo(ipoly)

                  do i = 1, baseloop ! loop for base lattice
                     icount = 0
                     imole = imole + 1
                     ifsetcor(imole) = .true.

                     j1 = molept_index(imole)
                     j2 = molept_index(imole+1) - 1

 110                 do j = j1, j2
                        k = molept_list(j)
                        patmcor(1,k) = atmcor(1,k)
                        patmcor(2,k) = atmcor(2,k)
                        patmcor(3,k) = atmcor(3,k)
                     end do

c                    --- make gravity coordinate ---

                     rxp = compfact*grnd(mt,mti) * xcel
                     ryp = compfact*grnd(mt,mti) * ycel
                     rzp = compfact*grnd(mt,mti) * zcel

c                    --- make rotational matrix ---
                     psi =   (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi
                     phi =   (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi
                     theta = (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi

                     rot_m(1,1) =   dcos(psi)*dcos(phi) 
     &                          -   dcos(theta)*dsin(phi)*dsin(psi)
                     rot_m(1,2) =   dcos(psi)*dsin(phi)
     &                          +   dcos(theta)*dcos(phi)*dsin(psi)
                     rot_m(1,3) =   dsin(psi)*dsin(theta)
                     rot_m(2,1) = - dsin(psi)*dcos(phi)
     &                          -   dcos(theta)*dsin(phi)*dcos(psi)
                     rot_m(2,2) = - dsin(psi)*dsin(phi)
     &                          +   dcos(theta)*dcos(phi)*dcos(psi)
                     rot_m(2,3) =   dcos(psi)*dsin(theta)
                     rot_m(3,1) =   dsin(theta)*dsin(phi)
                     rot_m(3,2) = - dsin(theta)*dcos(phi)
                     rot_m(3,3) =   dcos(theta)

c                    --- make coordinate ---
                     do j=j1,j2
                        k = molept_list(j)

                        atmcor_t(1) = atmcor(1,k)
                        atmcor_t(2) = atmcor(2,k)
                        atmcor_t(3) = atmcor(3,k)
                     
                        atmcor(1,k) =  rot_m(1,1)*atmcor_t(1)
     &                               + rot_m(2,1)*atmcor_t(2)
     &                               + rot_m(3,1)*atmcor_t(3)
     &                               + rxp
                        atmcor(2,k) =  rot_m(1,2)*atmcor_t(1)
     &                               + rot_m(2,2)*atmcor_t(2)
     &                               + rot_m(3,2)*atmcor_t(3)
     &                               + ryp
                        atmcor(3,k) =  rot_m(1,3)*atmcor_t(1)
     &                               + rot_m(2,3)*atmcor_t(2)
     &                               + rot_m(3,3)*atmcor_t(3)
     &                               + rzp
                     
c---- check the atom is very close to the other atoms
                        do m = 1, nmole
                          if (.not. ifsetcor(m) .or. imole .eq. m) cycle

                           n1 = molept_index(m)
                           n2 = molept_index(m+1) - 1
                           do n = n1, n2
                              nn = molept_list(n)

                              r_ij(1) = atmcor(1,k) - atmcor(1,nn)
                              r_ij(2) = atmcor(2,k) - atmcor(2,nn)
                              r_ij(3) = atmcor(3,k) - atmcor(3,nn)
                              r_ij(1) = r_ij(1) - xcel
     &                                * dnint(r_ij(1)/xcel)
                              r_ij(2) = r_ij(2) - ycel
     &                                * dnint(r_ij(2)/ycel)
                              r_ij(3) = r_ij(3) - zcel
     &                                * dnint(r_ij(3)/zcel)

                              rij = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2
                              rij = dsqrt(rij)

                              if (rij .lt. r_min) then
                                 do jj = j1,j2
                                    kk = molept_list(jj)
                                    atmcor(1,kk) = patmcor(1,kk)
                                    atmcor(2,kk) = patmcor(2,kk)
                                    atmcor(3,kk) = patmcor(3,kk)
                                 end do
                                 icount = icount + 1
                                 if (icount .gt. maxicount) then
                                    r_min = (r_min*xref - dmin)/xref
#ifdef MPI
                                    if (irank == 0) then
#endif
                                    write(6,*) 'r_min= ',r_min*xref
#ifdef MPI
                                    end if
#endif
                                    icount = 0
                                 end if
                                 goto 110
                              end if
                           end do
                        end do

                     end do
#ifdef MPI
                     if (irank == 0) then
#endif
                     write(6,*) 'Set atmcor of mole No. ',imole,',',
     &                          icount
#ifdef MPI
                     end if
#endif
                  end do

               end do
            end do
         end do

      end do

c-------- Locate molecules of H2O --------

      do icrewa = 1, ncrecorw

c---- some preparation
         r_min = rmin / xref
#ifdef MPI
         if (irank == 0) then
#endif
         write(6,*) 'r_min= ',r_min*xref
#ifdef MPI
         end if
#endif

         imole = npoly

c---- storage of water p-coordinate
         j1 = molept_index(imole+1) ! first water molecule
         j2 = molept_index(imole+1+1) - 1
         if ((nwater .gt. 0) .and. (j1+2 .ne. j2)) then
            write(6,*) 'Failure: molept_index is wrong to arrange H2O'
            stop
         end if
         k = molept_list(j1)
         watmcor(1,1) = atmcor(1,k)
         watmcor(2,1) = atmcor(2,k)
         watmcor(3,1) = atmcor(3,k)
         k = molept_list(j1+1)
         watmcor(1,2) = atmcor(1,k)
         watmcor(2,2) = atmcor(2,k)
         watmcor(3,2) = atmcor(3,k)
         k = molept_list(j1+2)
         watmcor(1,3) = atmcor(1,k)
         watmcor(2,3) = atmcor(2,k)
         watmcor(3,3) = atmcor(3,k)

c---- arrange water molecules
         do ix=0,xmaxw-1
            do iy=0,ymaxw-1
               do iz=0,zmaxw-1

                  icount = 0
                  imole = imole + 1
                  ifsetcor(imole) = .true.

                  j1 = molept_index(imole)
                  j2 = molept_index(imole+1) - 1

c              --- make rotational matrix ---           
 100              psi =   (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi
                  phi =   (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi
                  theta = (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi

                  rxw = grnd(mt,mti)*xcel
                  ryw = grnd(mt,mti)*ycel
                  rzw = grnd(mt,mti)*water_layer + lower_water

                  rot_m(1,1) =   dcos(psi)*dcos(phi) 
     &                         - dcos(theta)*dsin(phi)*dsin(psi)
                  rot_m(1,2) =   dcos(psi)*dsin(phi)
     &                         + dcos(theta)*dcos(phi)*dsin(psi)
                  rot_m(1,3) =   dsin(psi)*dsin(theta)
                  rot_m(2,1) = - dsin(psi)*dcos(phi)
     &                         - dcos(theta)*dsin(phi)*dcos(psi)
                  rot_m(2,2) = - dsin(psi)*dsin(phi)
     &                         + dcos(theta)*dcos(phi)*dcos(psi)
                  rot_m(2,3) =   dcos(psi)*dsin(theta)
                  rot_m(3,1) =   dsin(theta)*dsin(phi)
                  rot_m(3,2) = - dsin(theta)*dcos(phi)
                  rot_m(3,3) =   dcos(theta)

                  do j= j1, j2
                     k = molept_list(j)

                     atmcor_t(1) = atmcor(1,k)
                     atmcor_t(2) = atmcor(2,k)
                     atmcor_t(3) = atmcor(3,k)

                     atmcor(1,k) =  rot_m(1,1)*atmcor_t(1)
     &                         + rot_m(2,1)*atmcor_t(2)
     &                         + rot_m(3,1)*atmcor_t(3)
     &                         + rxw
                     atmcor(2,k) =  rot_m(1,2)*atmcor_t(1)
     &                         + rot_m(2,2)*atmcor_t(2)
     &                         + rot_m(3,2)*atmcor_t(3)
     &                         + ryw
                     atmcor(3,k) =  rot_m(1,3)*atmcor_t(1)
     &                         + rot_m(2,3)*atmcor_t(2)
     &                         + rot_m(3,3)*atmcor_t(3)
     &                         + rzw

c---- check the atom is very close to the other atoms
                     do m = 1, nmole
                        if (.not. ifsetcor(m) .or. imole .eq. m) cycle

                        n1 = molept_index(m)
                        n2 = molept_index(m+1) - 1
                        do n = n1, n2
                           nn = molept_list(n)

                           r_ij(1) = atmcor(1,k) - atmcor(1,nn)
                           r_ij(2) = atmcor(2,k) - atmcor(2,nn)
                           r_ij(3) = atmcor(3,k) - atmcor(3,nn)
                           r_ij(1) = r_ij(1) 
     &                             - xcel * dnint(r_ij(1)/xcel)
                           r_ij(2) = r_ij(2)
     &                             - ycel * dnint(r_ij(2)/ycel)
                           r_ij(3) = r_ij(3)
     &                             - zcel * dnint(r_ij(3)/zcel)

                           rij = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2
                           rij = dsqrt(rij)

                           if (rij .lt. r_min) then
                              kk = molept_list(j1)
                              atmcor(1,kk) = watmcor(1,1)
                              atmcor(2,kk) = watmcor(2,1)
                              atmcor(3,kk) = watmcor(3,1)
                              kk = molept_list(j1+1)
                              atmcor(1,kk) = watmcor(1,2)
                              atmcor(2,kk) = watmcor(2,2)
                              atmcor(3,kk) = watmcor(3,2)
                              kk = molept_list(j1+2)
                              atmcor(1,kk) = watmcor(1,3)
                              atmcor(2,kk) = watmcor(2,3)
                              atmcor(3,kk) = watmcor(3,3)
                              icount = icount + 1
                              if (icount .gt. maxicount) then
                                 r_min = (r_min*xref - dmin)/xref
#if defined(MPI)
                                 if (irank == 0) then
#endif
                                 write(6,*) 'r_min= ',r_min*xref
#if defined(MPI)
                                 end if
#endif
                                 icount = 0
                              end if
                              goto 100
                           end if
                        end do
                     end do

                  end do
#if defined(MPI)
                  if (irank == 0) then
#endif
                  write(6,*) 'Set atmcor of mole No. ',imole,',',icount
#if defined(MPI)
                  end if
#endif
               end do
            end do
         end do

      end do

c     ----- Preset the Nose-Hoover chain thermostat -----

#ifdef MPI
      if (irank == 0) then
#endif
      write(6,*) 'Set Nose-Hoover chain thermostat'
      write(6,*) 'chain number= ',mchain
#ifdef MPI
      end if
#endif
      do m=1,mchain
         xlogs(m) = 0.0d0
      end do

c     +     +     +     +     +     +     +

      return
      end
