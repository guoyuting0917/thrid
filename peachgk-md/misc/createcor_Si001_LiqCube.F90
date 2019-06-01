!*************************************
!*  createcor.f Ver.2.7 '15.01.15    *
!*      for peachgk_md.f             *
!*            by G.Kikugawa          *
!*      modified by T.Nakano         *
!*      for Diamond lattice (100)    *
!*************************************
subroutine createcor(npoly,npolytyp,npoly_mole,npoly_atom,   &
     &               nwater,nmatom,nmatyp,nmatomtyp,   &
     &               xmaxpo,ymaxpo,zmaxpo,   &
     &               xmaxw,ymaxw,zmaxw,   &
     &               xmaxma,ymaxma,zmaxma,   &
     &               inicorpo,inicorw,inicorma,   &
     &               ncrecorpo,index_crecorpo,   &
     &               ncrecorw,index_crecorw,   &
     &               ncrecorma,index_crecorma,   &
     &               ifsetcor,   &
     &               xcel,ycel,zcel,   &
     &               xref,   &
     &               compfact,   &
     &               mchain)

!#################################################################
!# Usage of the inicorpo(), inicorw(), inicorma()                #
!#   inicorpo(), inicorw(), inicorma() are the real*8 type.      #
!#   Therefore, you have to cast the type which you want to use. #
!#   For example, int(inicorpo()) to cast the int type           #
!#################################################################

  use md_global
#if defined(MPI)
  use mpi_global
#endif

  implicit none

! ARGUMENT:
!   INPUT
  integer,intent(in):: npoly            ! all number of poly
  integer,intent(in):: npolytyp         ! number of poly type
  integer,intent(in):: npoly_mole(:)    ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)    ! number of atoms belonging to poly

  integer,intent(in):: nwater           ! number of H2O molecules

  integer,intent(in):: nmatom           ! number of monatomic molecules
  integer,intent(in):: nmatyp           ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)     ! each number of monatomic mole.

  integer,intent(in):: xmaxpo(:)        ! use for positioning of polymer1
  integer,intent(in):: ymaxpo(:)        ! use for positioning of polymer1
  integer,intent(in):: zmaxpo(:)        ! use for positioning of polymer1
  integer,intent(in):: xmaxw            ! use for positioning of water
  integer,intent(in):: ymaxw            ! use for positioning of water
  integer,intent(in):: zmaxw            ! use for positioning of water
  integer,intent(in):: xmaxma(:)        ! use for positioning of monatomic mole.
  integer,intent(in):: ymaxma(:)        ! use for positioning of monatomic mole.
  integer,intent(in):: zmaxma(:)        ! use for positioning of monatomic mole.
  real(8),intent(in):: inicorpo(:,:) ! use for positioning of poly
  real(8),intent(in):: inicorw(:)       ! use for positioning of water
  real(8),intent(in):: inicorma(:,:) ! use for positioning of ma
  integer,intent(in):: ncrecorpo        ! max number of poly type for createcor
  integer,intent(in):: index_crecorpo(:) ! index of polymer for createcor
  integer,intent(in):: ncrecorw         ! water for createcor
  integer,intent(in):: index_crecorw    ! index of water type for createcor
  integer,intent(in):: ncrecorma        ! max number of matom type for createcor
  integer,intent(in):: index_crecorma(:) ! index of matom type for createcor

  logical,intent(out):: ifsetcor(:)     ! frag for checking if coordinate has set

  real(8),intent(in):: xref             ! distanse base value [m]

  real(8),intent(in):: xcel             ! x cell length
  real(8),intent(in):: ycel             ! y cell length
  real(8),intent(in):: zcel             ! z cell length

!!! old MT random generator
!      real*8:: randseed         ! random seed for createcor etc.
!      integer:: nperiod         ! period parameter
!      parameter(nperiod = 624)
!      integer:: mt(0:nperiod-1)
!      integer:: mti

  real(8),intent(in):: compfact     ! compact factor using at poly arrange(<1.0)

  integer,intent(in):: mchain           ! Nose-Hoover chain number

! LOCAL:
  real(8):: rmin = 2.0d-10              ! minimum of atom distance 
  integer:: maxicount = 1000000         ! max count for arrange water
  real(8):: dmin = 0.05d-10             ! difference of rmin

  integer:: i
  integer:: j,jj                        ! do loop index
  integer:: m,n                         ! do loop index
  integer:: n1,n2,nn
  integer:: j1,j2,k,kk
  integer:: ix,iy,iz
  integer:: imole
  integer:: icount

  real(8):: rot_m(3,3)                  ! rotational matrix
  real(8):: psi,phi,theta               ! Euler angle
  real(8):: rxp,ryp,rzp                 ! coordinate of poly1
  real(8):: rxw,ryw,rzw                 ! coordinate of water
  real(8):: rxma,ryma,rzma              ! coordinate of monatomic mole.

  real(8):: pi                          ! = 3.14159...

  real(8):: atmcor_t(3)                 ! atmcor tmp

  real(8):: patmcor(3,maxnatom)         ! tmp atmcor for poly1
  real(8):: watmcor(3,3)                ! tmp atmcor for water
  real(8):: maatmcor(3)                 ! tmp atmcor for MA

  real(8):: r_ij(3)                     ! xi - xj etc.
  real(8):: rij                         ! distance of r_ij
  real(8):: r_min                       ! distance of r_min

  integer:: ipoly
  integer:: imatom

  integer:: icrepo, icrewa, icrema

! parameters for diamond lattice
  real(8):: Lc                          ! lattice constant
  real(8):: A0(3)
  real(8):: snx,sny,snz

  integer:: nxlatt, nylatt, nzlatt      ! number of lattices in each direction
  integer,parameter:: natomlatt = 8     ! number of mole. in a basic lattice

  real(8):: base_latt(3,natomlatt)
  real(8):: shift_latt(3)

  real(8):: expand_coeff

! FUNCTIONS:
!  real(8):: grnd                        ! function creating random number
  real(8):: genrand_res53                ! SFMT random generator
      
!     +     +     +     +     +     +     +

!---- notation ----
#if defined(MPI)
  if (irank == 0) then
#endif

     write(6,*) 'Create atmcor from createcor'

#if defined(MPI)
  end if
#endif

  pi = acos(-1.0d0)

  Lc = 5.43d-10 / xref    ! lattice constant of diamond structure of silicon
!  Lc = 3.56d-10 / xref    ! lattice constant of diamond structure of carbon
!  Lc = 5.65d-10 / xref    ! lattice constant of diamond structure of germanium

  A0(1:3) = Lc * 0.25d0

  shift_latt(1:3) = Lc

! #       x       y       z
! 1   0.000   0.000   0.000
! 2   2.715   2.715   0.000
! 3   1.358   1.358   1.358
! 4   4.073   4.073   1.358
! 5   2.715   0.000   2.715
! 6   0.000   2.715   2.715
! 7   1.358   4.073   4.073
! 8   4.073   1.358   4.073

!---- create base lattice for Si
  base_latt(1,1) = 0.0d0 * A0(1)
  base_latt(2,1) = 0.0d0 * A0(2)
  base_latt(3,1) = 0.0d0 * A0(3)

  base_latt(1,2) = 2.0d0 * A0(1)
  base_latt(2,2) = 2.0d0 * A0(2)
  base_latt(3,2) = 0.0d0 * A0(3)

  base_latt(1,3) = 1.0d0 * A0(1)
  base_latt(2,3) = 1.0d0 * A0(2)
  base_latt(3,3) = 1.0d0 * A0(3)

  base_latt(1,4) = 3.0d0 * A0(1)
  base_latt(2,4) = 3.0d0 * A0(2)
  base_latt(3,4) = 1.0d0 * A0(3)

  base_latt(1,5) = 2.0d0 * A0(1)
  base_latt(2,5) = 0.0d0 * A0(2)
  base_latt(3,5) = 2.0d0 * A0(3)

  base_latt(1,6) = 0.0d0 * A0(1)
  base_latt(2,6) = 2.0d0 * A0(2)
  base_latt(3,6) = 2.0d0 * A0(3)

  base_latt(1,7) = 1.0d0 * A0(1)
  base_latt(2,7) = 3.0d0 * A0(2)
  base_latt(3,7) = 3.0d0 * A0(3)

  base_latt(1,8) = 3.0d0 * A0(1)
  base_latt(2,8) = 1.0d0 * A0(2)
  base_latt(3,8) = 3.0d0 * A0(3)

!-------- Locate molecules of monatomic Si mole. --------

     icrema = 1
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

!----    storage of MA p-coordinate
     j1 = molept_index(imole+1) ! first monatomic mole.
     j2 = molept_index(imole+1+1) - 1
     if ((nmatom .gt. 0) .and. (j1 .ne. j2)) then
        write(6,*) 'Failure: molept_index is wrong to arrange MA'
        stop
     end if
     k = molept_list(j1)
     maatmcor(1:3) = atmcor(1:3,k)

     nxlatt = xmaxma(imatom) / 2
     nylatt = ymaxma(imatom) / 2
     nzlatt = zmaxma(imatom) / 2

!----    arrange MA molecules (for (100) surface)
     do iz=1,nzlatt
        do iy=1,nylatt
           do ix=1,nxlatt

              do i = 1, natomlatt  ! loop for basic lattice
                 icount = 0
                 imole = imole + 1
                 ifsetcor(imole) = .true.
                     
                 j1 = molept_index(imole)
                 j2 = molept_index(imole+1) - 1

                 snx = base_latt(1,i) + shift_latt(1) * DBLE(ix-1)
                 sny = base_latt(2,i) + shift_latt(2) * DBLE(iy-1)
                 snz = base_latt(3,i) + shift_latt(3) * DBLE(iz-1)

                 do j= j1, j2
                    k = molept_list(j)

!                    atmcor_t(1:3) = atmcor(1:3,k)

                    atmcor(1,k) = snx
                    atmcor(2,k) = sny
                    atmcor(3,k) = snz - 0.5d0*DBLE(nzlatt)*shift_latt(3)

                 end do
#ifdef MPI
                 if (irank == 0) then
#endif
                    write(6,*) 'Set atmcor of mole No. ',imole,',', icount
#ifdef MPI
                 end if
#endif
              end do

           end do
        end do
     end do


!-------- Locate molecules of poly --------

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

     do ix=0,xmaxpo(ipoly)-1
        do iy=0,ymaxpo(ipoly)-1
           do iz=0,zmaxpo(ipoly)-1

              icount = 0
              imole = imole + 1
              ifsetcor(imole) = .true.

              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

              do j = j1, j2
                 k = molept_list(j)
                 patmcor(1:3,k) = atmcor(1:3,k)
              end do

!             --- make gravity coordinate ---
!110           rxp = compfact*(genrand_res53()-0.5d0)*xcel + 0.5d0*xcel
!              ryp = compfact*(genrand_res53()-0.5d0)*ycel + 0.5d0*ycel
!              rzp = compfact*(genrand_res53()-0.5d0)*zcel + 0.5d0*zcel

110           rxp = compfact*(genrand_res53()-0.5d0)*xcel * 0.75d0 + xcel*0.5d0
              ryp = compfact*(genrand_res53()-0.5d0)*ycel + 0.5d0*ycel
!              rzp = compfact* genrand_res53()*zcel / 3.0d0 + 0.5d0*DBLE(nzlatt+1)*shift_latt(3)
              rzp = compfact* genrand_res53()*30.0d0 + 0.5d0*DBLE(nzlatt+1)*shift_latt(3)

!             --- make rotational matrix ---
              psi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
              phi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
              theta = (genrand_res53() - 0.5d0) * 2.0d0 * pi
!              psi   = 0.5d0 * pi
!              phi   = 0.0d0
!              theta = 0.5d0 * pi

              rot_m(1,1) =   cos(psi)*cos(phi)   &
                   &       - cos(theta)*sin(phi)*sin(psi)
              rot_m(1,2) =   cos(psi)*sin(phi)   &
                   &       + cos(theta)*cos(phi)*sin(psi)
              rot_m(1,3) =   sin(psi)*sin(theta)
              rot_m(2,1) = - sin(psi)*cos(phi)   &
                   &       - cos(theta)*sin(phi)*cos(psi)
              rot_m(2,2) = - sin(psi)*sin(phi)   &
                   &       + cos(theta)*cos(phi)*cos(psi)
              rot_m(2,3) =   cos(psi)*sin(theta)
              rot_m(3,1) =   sin(theta)*sin(phi)
              rot_m(3,2) = - sin(theta)*cos(phi)
              rot_m(3,3) =   cos(theta)

              !--- make coordinate ---
              do j=j1,j2

                 k = molept_list(j)

                 atmcor_t(1:3) = atmcor(1:3,k)

                 atmcor(1,k) = rot_m(1,1)*atmcor_t(1)   &
                      &      + rot_m(2,1)*atmcor_t(2)   &
                      &      + rot_m(3,1)*atmcor_t(3)   &
                      &      + rxp
                 atmcor(2,k) = rot_m(1,2)*atmcor_t(1)   &
                      &      + rot_m(2,2)*atmcor_t(2)   &
                      &      + rot_m(3,2)*atmcor_t(3)   &
                      &      + ryp
                 atmcor(3,k) = rot_m(1,3)*atmcor_t(1)   &
                      &      + rot_m(2,3)*atmcor_t(2)   &
                      &      + rot_m(3,3)*atmcor_t(3)   &
                      &      + rzp

                 !---- check the atom is very close to the other atoms
                 do m = 1, nmole
                    if (.not. ifsetcor(m) .or. (imole == m)) cycle

                    n1 = molept_index(m)
                    n2 = molept_index(m+1) - 1
                    do n = n1, n2
                       nn = molept_list(n)

                       r_ij(1:3) = atmcor(1:3,k) - atmcor(1:3,nn)
                       r_ij(1) = r_ij(1) - xcel * anint(r_ij(1)/xcel)
                       r_ij(2) = r_ij(2) - ycel * anint(r_ij(2)/ycel)
                       r_ij(3) = r_ij(3) - zcel * anint(r_ij(3)/zcel)

                       rij = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2
                       rij = sqrt(rij)

                       if (rij < r_min) then
                          do jj = j1,j2
                             kk = molept_list(jj)
                             atmcor(1:3,kk) = patmcor(1:3,kk)
                          end do
                          icount = icount + 1
                          if (icount > maxicount) then
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
                          goto 110
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

!-------- Locate molecules of H2O --------

  do icrewa = 1, ncrecorw

     !---- some preparation
     r_min = rmin / xref
#if defined(MPI)
     if (irank == 0) then
#endif
        write(6,*) 'r_min= ',r_min*xref
#if defined(MPI)
     end if
#endif

     imole = npoly

     !---- storage of water p-coordinate
     j1 = molept_index(imole+1) ! first water molecule
     j2 = molept_index(imole+1+1) - 1
     if ((nwater > 0) .and. (j1+2 /= j2)) then
        write(6,*) 'Failure: molept_index is wrong to arrange H2O'
        stop
     end if
     k = molept_list(j1)
     watmcor(1:3,1) = atmcor(1:3,k)
     k = molept_list(j1+1)
     watmcor(1:3,2) = atmcor(1:3,k)
     k = molept_list(j1+2)
     watmcor(1:3,3) = atmcor(1:3,k)

     !---- arrange water molecules
     do ix=0,xmaxw-1
        do iy=0,ymaxw-1
           do iz=0,zmaxw-1

              icount = 0
              imole = imole + 1
              ifsetcor(imole) = .true.

              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

              !--- make rotational matrix ---           
100           psi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
              phi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
              theta = (genrand_res53() - 0.5d0) * 2.0d0 * pi

!              rxw = genrand_res53()*xcel * 0.5d0 + xcel*0.25d0
              rxw =(genrand_res53()-0.5d0)*xcel * 0.75d0 + xcel*0.5d0
              ryw = genrand_res53()*ycel
              rzw = genrand_res53()*zcel / 2.0d0 + 0.5d0*DBLE(nzlatt+1)*shift_latt(3)
!              rzw = genrand_res53()*zcel + (1.0d0-genrand_res53())*DBLE(nzlatt+1)*shift_latt(3)

              rot_m(1,1) =   cos(psi)*cos(phi)   &
                   &       - cos(theta)*sin(phi)*sin(psi)
              rot_m(1,2) =   cos(psi)*sin(phi)   &
              &            + cos(theta)*cos(phi)*sin(psi)
              rot_m(1,3) =   sin(psi)*sin(theta)
              rot_m(2,1) = - sin(psi)*cos(phi)   &
                   &       - cos(theta)*sin(phi)*cos(psi)
              rot_m(2,2) = - sin(psi)*sin(phi)   &
                   &       + cos(theta)*cos(phi)*cos(psi)
              rot_m(2,3) =   cos(psi)*sin(theta)
              rot_m(3,1) =   sin(theta)*sin(phi)
              rot_m(3,2) = - sin(theta)*cos(phi)
              rot_m(3,3) =   cos(theta)

              do j= j1, j2
                 k = molept_list(j)
                    
                 atmcor_t(1:3) = atmcor(1:3,k)

                 atmcor(1,k) = rot_m(1,1)*atmcor_t(1)   &
                      &      + rot_m(2,1)*atmcor_t(2)   &
                      &      + rot_m(3,1)*atmcor_t(3)   &
                      &      + rxw
                 atmcor(2,k) = rot_m(1,2)*atmcor_t(1)   &
                      &      + rot_m(2,2)*atmcor_t(2)   &
                      &      + rot_m(3,2)*atmcor_t(3)   &
                      &      + ryw
                 atmcor(3,k) = rot_m(1,3)*atmcor_t(1)   &
                      &      + rot_m(2,3)*atmcor_t(2)   &
                      &      + rot_m(3,3)*atmcor_t(3)   &
                      &      + rzw

                 !---- check the atom is very close to the other atoms
                 do m = 1, nmole
                    if (.not. ifsetcor(m) .or. (imole == m)) cycle

                    n1 = molept_index(m)
                    n2 = molept_index(m+1) - 1
                    do n = n1, n2
                       nn = molept_list(n)

                       r_ij(1:3) = atmcor(1:3,k) - atmcor(1:3,nn)
                       r_ij(1) = r_ij(1) - xcel * anint(r_ij(1)/xcel)
                       r_ij(2) = r_ij(2) - ycel * anint(r_ij(2)/ycel)
                       r_ij(3) = r_ij(3) - zcel * anint(r_ij(3)/zcel)

                       rij = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2
                       rij = sqrt(rij)

                       if (rij < r_min) then
                          kk = molept_list(j1)
                          atmcor(1:3,kk) = watmcor(1:3,1)
                          kk = molept_list(j1+1)
                          atmcor(1:3,kk) = watmcor(1:3,2)
                          kk = molept_list(j1+2)
                          atmcor(1:3,kk) = watmcor(1:3,3)
                          icount = icount + 1
                          if (icount > maxicount) then
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

!-------- Locate molecules of monatomic mole. --------

  do icrema = 2, ncrecorma
     imatom = index_crecorma(icrema)

     imole = npoly + nwater
     do i = 1, imatom-1
        imole = imole + nmatomtyp(i)
     end do

     r_min = rmin / xref
#if defined(MPI)
     if (irank == 0) then
#endif
        write(6,*) 'r_min= ',r_min*xref
#if defined(MPI)
     end if
#endif

     !---- storage of MA p-coordinate
     j1 = molept_index(imole+1)         ! first monatomic mole.
     j2 = molept_index(imole+1+1) - 1
     if ((nmatom > 0) .and. (j1 /= j2)) then
        write(6,*) 'Failure: molept_index is wrong to arrange MA'
        stop
     end if
     k = molept_list(j1)
     maatmcor(1:3) = atmcor(1:3,k)

     !---- arrange MA molecules
     do ix=0,xmaxma(imatom)-1
        do iy=0,ymaxma(imatom)-1
           do iz=0,zmaxma(imatom)-1

              icount = 0
              imole = imole + 1
              ifsetcor(imole) = .true.

              j1 = molept_index(imole)
              j2 = molept_index(imole+1) - 1

90            continue

              rxma = genrand_res53()*xcel
              ryma = genrand_res53()*ycel
              rzma = genrand_res53()*zcel

              do j= j1, j2
                 k = molept_list(j)

                 atmcor_t(1:3) = atmcor(1:3,k)

                 atmcor(1,k) = rxma
                 atmcor(2,k) = ryma
                 atmcor(3,k) = rzma

                 !---- check the atom is very close to the other atoms
                 do m = 1, nmole
                    if (.not. ifsetcor(m) .or. (imole == m)) cycle

                    n1 = molept_index(m)
                    n2 = molept_index(m+1) - 1
                    do n = n1, n2
                       nn = molept_list(n)

                       r_ij(1:3) = atmcor(1:3,k) - atmcor(1:3,nn)
                       r_ij(1) = r_ij(1) - xcel * anint(r_ij(1)/xcel)
                       r_ij(2) = r_ij(2) - ycel * anint(r_ij(2)/ycel)
                       r_ij(3) = r_ij(3) - zcel * anint(r_ij(3)/zcel)

                       rij = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2
                       rij = sqrt(rij)

                       if (rij < r_min) then
                          kk = molept_list(j1)
                          atmcor(1:3,kk) = maatmcor(1:3)
                          icount = icount + 1
                          if (icount > maxicount) then
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
                          goto 90
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

!----- Preset the Nose-Hoover chain thermostat -----

#if defined(MPI)
  if (irank == 0) then
#endif
     write(6,*) 'Set Nose-Hoover chain thermostat'
     write(6,*) 'chain number= ',mchain
#if defined(MPI)
  end if
#endif
  xlogs(1:mchain) = 0.0d0

!     +     +     +     +     +     +     +

end subroutine createcor
