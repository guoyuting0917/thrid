!*************************************************
!*  createcor.f Ver.2.5cu-thi-water-cu '13.12.25 *
!*      for peachgk_md.f                         *
!*            by G.Kikugawa                      *
!**************************************************
subroutine createcor(npoly,npolytyp,npoly_mole,npoly_atom, &
     &               nwater,nmatom,nmatyp,nmatomtyp, &
     &               xmaxpo,ymaxpo,zmaxpo, &
     &               xmaxw,ymaxw,zmaxw, &
     &               xmaxma,ymaxma,zmaxma, &
     &               inicorpo,inicorw,inicorma, &
     &               ncrecorpo,index_crecorpo, &
     &               ncrecorw,index_crecorw, &
     &               ncrecorma,index_crecorma, &
     &               ifsetcor, &
     &               xcel,ycel,zcel, &
     &               xref, &
     &               compfact, &
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

!   INPUT
  integer,intent(in):: npoly            ! all number of poly
  integer,intent(in):: npolytyp         ! number of poly type
  integer,intent(in):: npoly_mole(:)    ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)    ! number of atoms belonging to poly

  integer,intent(in):: nwater           ! number of H2O molecules

  integer,intent(in):: nmatom           ! number of monatomic molecules
  integer,intent(in):: nmatyp           ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)     ! each number of monatomic mole.

  integer,intent(inout):: xmaxpo(:)        ! use for positioning of polymer1
  integer,intent(inout):: ymaxpo(:)        ! use for positioning of polymer1
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
!  real(8):: rxma,ryma,rzma              ! coordinate of monatomic mole.

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

!     parameters for fcc lattice
  real(8):: A0,A1,A2,A3
  real(8):: snx,sny,snz

  real(8):: vacuum_layer
  real(8):: water_layer
  real(8):: thiol_layer
  real(8):: lower_S
  real(8):: upper_S 
  real(8):: lower_water_z
  real(8):: water_layer_z
  real(8):: upper_AU_z

  real(8),allocatable:: base_latt(:,:)
  real(8):: shift_latt(2,0:2)

  real(8),allocatable:: base_latt_S1(:,:)
  real(8),allocatable:: base_latt_S2(:,:)
  real(8):: x_shift,y_shift

  real(8):: baseloop

  integer:: nunit_x,nunit_y
  integer:: nunit_x_org,nunit_y_org
  real(8):: s_scale

! FUNCTIONS:
!  real(8):: grnd                        ! function creating random number
  real(8):: genrand_res53                ! SFMT random generator

!     +     +     +     +     +     +     +

!---- notation ----
#ifdef MPI
  if (irank == 0) then
#endif

     write(6,*) 'Create atmcor from createcor'

#ifdef MPI
  end if
#endif

  pi = dacos(-1.0d0)

  A0=2.556d-10 / xref
  A1=A0*sqrt(3.0)*0.5
  A2=sqrt(6.0)/3.0*A0

  A3=1.9d-10 / xref

  vacuum_layer = 0.0d-10 / xref
  water_layer = 90.0d-10 / xref
  thiol_layer = 14.08d-10 / xref

  lower_S = vacuum_layer + DBLE(zmaxma(1)-1) * A2 + A3
  upper_S = lower_S + thiol_layer*2.0d0 + water_layer

  lower_water_z = lower_S + thiol_layer + 5.0d0
  water_layer_z = water_layer - 2.0d0*5.0d0

  upper_AU_z = upper_S + A3

  ! for definition of surface unit
  nunit_y = 30
  nunit_y_org = 32
  nunit_x = nunit_y / 2           ! x unit length is half of y
  nunit_x_org = nunit_y_org / 2

  s_scale = dble(nunit_x)/dble(nunit_x_org)

  ! memory allocation
  allocate(base_latt(2,nunit_x*nunit_y))
  allocate(base_latt_S1(2,4*nunit_x_org/4*nunit_y_org/4))
  allocate(base_latt_S2(2,4*nunit_x_org/4*nunit_y_org/4))

!---- create base lattice for Cu (0 layer)
  do ix =1,nunit_x
     base_latt(1,ix) = dble(ix-1)*A0
     base_latt(2,ix) = 0.0d0*A1
  end do

  icount = nunit_x
  do iy=2,nunit_y
     do ix=1,nunit_x
        icount = icount + 1
        base_latt(1,icount) = dble(mod(iy-1,2))*0.5d0*A0 + base_latt(1,ix)
        base_latt(2,icount) = dble(iy-1)*A1 + base_latt(2,ix)
     end do
  end do

!---- shift amount of base lattice for Cu
  shift_latt(1,0) = 0.0d0*A0
  shift_latt(2,0) = 0.0d0*A1

  shift_latt(1,1) = 0.5d0*A0
  shift_latt(2,1) = 1.0d0/3.0d0*A1

  shift_latt(1,2) = 0.0d0*A0
  shift_latt(2,2) = 2.0d0/3.0d0*A1

!---- create base lattice for thiol S (2 layer over 1 layer)
  base_latt_S1(1,1) = 0.5d0 * A0
  base_latt_S1(2,1) = (0.0d0+1.0d0/3.0d0) * A1
  base_latt_S1(1,2) = 2.5d0 * A0
  base_latt_S1(2,2) = (0.0d0+1.0d0/3.0d0) * A1
  base_latt_S1(1,3) = 1.5d0 * A0
  base_latt_S1(2,3) = (2.0d0+1.0d0/3.0d0) * A1
  base_latt_S1(1,4) = 3.5d0 * A0
  base_latt_S1(2,4) = (2.0d0+1.0d0/3.0d0) * A1

  ! copy to full unit range
  icount = 0
  do iy = 1, nunit_y_org/4
     do ix = 1, nunit_x_org/4
        do i = 1, 4
           icount = icount + 1
           base_latt_S1(1,icount) = base_latt_S1(1,i) + dble(ix-1)*4.0*A0
           base_latt_S1(2,icount) = base_latt_S1(2,i) + dble(iy-1)*4.0*A1
        end do
     end do
  end do

  baseloop = icount

  ! copy to base_latt_S2
  base_latt_S2(1:2,1:baseloop) = base_latt_S1(1:2,1:baseloop)

  ! x,y shift
  x_shift = 0.0d0 * A0
  y_shift = 2.0d0/3.0d0 * A1

  base_latt_S1(1,1:baseloop) = base_latt_S1(1,1:baseloop) + x_shift
  base_latt_S1(2,1:baseloop) = base_latt_S1(2,1:baseloop) + y_shift

!---- create base lattice for thiol S (0 layer over 1 layer)
  x_shift = -0.5d0 * A0
  y_shift = 1.0d0/3.0d0 * A1

  base_latt_S2(1,1:baseloop) = base_latt_S2(1,1:baseloop) + x_shift
  base_latt_S2(2,1:baseloop) = base_latt_S2(2,1:baseloop) + y_shift

!---- scale by scaling factor to shrink surface unit

  base_latt_S1(1:2,1:baseloop) = base_latt_S1(1:2,1:baseloop) * s_scale

  base_latt_S2(1:2,1:baseloop) = base_latt_S2(1:2,1:baseloop) * s_scale

!-------- Locate molecules of monatomic mole. --------

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

!----    storage of MA p-coordinate
     j1 = molept_index(imole+1) ! first monatomic mole.
     j2 = molept_index(imole+1+1) - 1
     if ((nmatom .gt. 0) .and. (j1 .ne. j2)) then
        write(6,*) 'Failure: molept_index is wrong to arrange MA'
        stop
     end if
     k = molept_list(j1)
     maatmcor(1:3) = atmcor(1:3,k)

!----    arrange MA molecules (for (111) surface)
     do iz=1,zmaxma(imatom)
        do iy=1,ymaxma(imatom) / nunit_y
           do ix=1,xmaxma(imatom) / nunit_x

              do i = 1, nunit_x*nunit_y  ! loop for basic lattice
                 icount = 0
                 imole = imole + 1
                 ifsetcor(imole) = .true.
                     
                 j1 = molept_index(imole)
                 j2 = molept_index(imole+1) - 1

90               continue

                 if (icrema == 1) then ! lower layer
                  
                    snx = DBLE(ix-1) * A0*dble(nunit_x) &
                         & + base_latt(1,i) &
                         & + shift_latt(1,mod(iz-1,3))
                    sny = DBLE(iy-1) * A1*dble(nunit_y) &
                         & + base_latt(2,i) &
                         & + shift_latt(2,mod(iz-1,3))
                    snz = DBLE(iz-1) * A2 + vacuum_layer

                 else if (icrema == 2) then ! upper layer
                        
                    snx = DBLE(ix-1) * A0*dble(nunit_x) &
                         & + base_latt(1,i) &
                         & + shift_latt(1,mod(iz-1,3))
                    sny = DBLE(iy-1) * A1*dble(nunit_y) &
                         & + base_latt(2,i) &
                         & + shift_latt(2,mod(iz-1,3))
                    snz = DBLE(iz-1) * A2 &
                         & + upper_AU_z

                 end if

                 do j= j1, j2
                    k = molept_list(j)

                    atmcor_t(1:3) = atmcor(1:3,k)

                    atmcor(1,k) = snx
                    atmcor(2,k) = sny
                    atmcor(3,k) = snz

!---- check the atom is very close to the other atoms
                    do m = 1, nmole
                       if (.not. ifsetcor(m) .or. imole == m) cycle

                       n1 = molept_index(m)
                       n2 = molept_index(m+1) - 1
                       do n = n1, n2
                          nn = molept_list(n)

                          r_ij(1) = atmcor(1,k) - atmcor(1,nn)
                          r_ij(2) = atmcor(2,k) - atmcor(2,nn)
                          r_ij(3) = atmcor(3,k) - atmcor(3,nn)
                          r_ij(1) = r_ij(1) - xcel*dnint(r_ij(1)/xcel)
                          r_ij(2) = r_ij(2) - ycel*dnint(r_ij(2)/ycel)
                          r_ij(3) = r_ij(3) - zcel*dnint(r_ij(3)/zcel)

                          rij = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2
                          rij = dsqrt(rij)

                          if (rij < r_min) then
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

  end do

!-------- Locate molecules of poly (for arranging thiol) --------

  do icrepo = 1, ncrecorpo
     ipoly = index_crecorpo(icrepo)

     imole = 0
     do i = 1, ipoly-1
        imole = imole + npoly_mole(i)
     end do

     r_min = rmin / xref
#ifdef MPI
     if (irank == 0) then
#endif
        write(6,*) 'r_min= ',r_min*xref
#ifdef MPI
     end if
#endif
     if ((ipoly == 1) .or. (ipoly == 2)) then
        ymaxpo(ipoly) = ymaxpo(ipoly) / (nunit_y_org/2)
        xmaxpo(ipoly) = xmaxpo(ipoly) / (nunit_x_org/2)
!        baseloop = baseloop
     else
        baseloop = 1
     end if

     do iz=1,zmaxpo(ipoly)
        do iy=1,ymaxpo(ipoly)
           do ix=1,xmaxpo(ipoly)

              do i = 1, baseloop ! loop for base lattice
                 icount = 0
                 imole = imole + 1
                 ifsetcor(imole) = .true.

                 j1 = molept_index(imole)
                 j2 = molept_index(imole+1) - 1

110              do j = j1, j2
                    k = molept_list(j)
                    patmcor(1,k) = atmcor(1,k)
                    patmcor(2,k) = atmcor(2,k)
                    patmcor(3,k) = atmcor(3,k)
                 end do

!                    --- make gravity coordinate ---

                 if (ipoly == 1)  then ! lower thiol
                    rxp = DBLE(ix-1) * A0*dble(nunit_x) &
                         & + base_latt_S1(1,i)
                    ryp = DBLE(iy-1) * A1*dble(nunit_y) &
                         & + base_latt_S1(2,i)
                    rzp = lower_S

!                       --- make rotational matrix ---
                    psi =   0.0d0
                    phi =   0.0d0
                    theta = 0.0d0

                 else if (ipoly == 2)  then ! upper thiol

                    rxp = DBLE(ix-1) * A0*dble(nunit_x) &
                         & + base_latt_S2(1,i)
                    ryp = DBLE(iy-1) * A1*dble(nunit_y) &
                         & + base_latt_S2(2,i)
                    rzp = upper_S

!                       --- make rotational matrix ---
                    psi =   0.0d0
                    phi =   0.0d0
                    theta = pi

                 else if (ipoly == 3)  then ! toluene

                    rxp = genrand_res53() * xcel
                    ryp = genrand_res53() * ycel
                    rzp = genrand_res53() * water_layer_z + lower_water_z

!                       --- make rotational matrix ---
                    psi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
                    phi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
                    theta = (genrand_res53() - 0.5d0) * 2.0d0 * pi

                 end if

                 rot_m(1,1) =   dcos(psi)*dcos(phi) &
                      &     -   dcos(theta)*dsin(phi)*dsin(psi)
                 rot_m(1,2) =   dcos(psi)*dsin(phi) &
                      &     +   dcos(theta)*dcos(phi)*dsin(psi)
                 rot_m(1,3) =   dsin(psi)*dsin(theta)
                 rot_m(2,1) = - dsin(psi)*dcos(phi) &
                      &     -   dcos(theta)*dsin(phi)*dcos(psi)
                 rot_m(2,2) = - dsin(psi)*dsin(phi) &
                      &     +   dcos(theta)*dcos(phi)*dcos(psi)
                 rot_m(2,3) =   dcos(psi)*dsin(theta)
                 rot_m(3,1) =   dsin(theta)*dsin(phi)
                 rot_m(3,2) = - dsin(theta)*dcos(phi)
                 rot_m(3,3) =   dcos(theta)

!                    --- make coordinate ---
                 do j=j1,j2
                    k = molept_list(j)

                    atmcor_t(1:3) = atmcor(1:3,k)
                     
                    atmcor(1,k) =  rot_m(1,1)*atmcor_t(1) &
                         &       + rot_m(2,1)*atmcor_t(2) &
                         &       + rot_m(3,1)*atmcor_t(3) &
                         &       + rxp
                    atmcor(2,k) =  rot_m(1,2)*atmcor_t(1) &
                         &       + rot_m(2,2)*atmcor_t(2) &
                         &       + rot_m(3,2)*atmcor_t(3) &
                         &       + ryp
                    atmcor(3,k) =  rot_m(1,3)*atmcor_t(1) &
                         &       + rot_m(2,3)*atmcor_t(2) &
                         &       + rot_m(3,3)*atmcor_t(3) &
                         &       + rzp
                     
!---- check the atom is very close to the other atoms

                    if (ipoly > 2) then
                    do m = 1, nmole
                       if (.not. ifsetcor(m) .or. imole == m) cycle

                       n1 = molept_index(m)
                       n2 = molept_index(m+1) - 1
                       do n = n1, n2
                          nn = molept_list(n)

                          r_ij(1:3) = atmcor(1:3,k) - atmcor(1:3,nn)
                          r_ij(1) = r_ij(1) - xcel*dnint(r_ij(1)/xcel)
                          r_ij(2) = r_ij(2) - ycel*dnint(r_ij(2)/ycel)
                          r_ij(3) = r_ij(3) - zcel*dnint(r_ij(3)/zcel)

                          rij = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2
                          rij = dsqrt(rij)

                          if (rij < r_min) then
                             do jj = j1,j2
                                kk = molept_list(jj)
                                atmcor(1:3,kk) = patmcor(1:3,kk)
                             end do
                             icount = icount + 1
                             if (icount > maxicount) then
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
                    end if

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

  end do

!-------- Locate molecules of H2O --------

!---- some preparation
  r_min = rmin / xref
#ifdef MPI
  if (irank == 0) then
#endif
     write(6,*) 'r_min= ',r_min*xref
#ifdef MPI
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

!---- arrange water molecules
  do ix=0,xmaxw-1
     do iy=0,ymaxw-1
        do iz=0,zmaxw-1

           icount = 0
           imole = imole + 1
           ifsetcor(imole) = .true.

           j1 = molept_index(imole)
           j2 = molept_index(imole+1) - 1

!              --- make rotational matrix ---           
100        psi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
           phi =   (genrand_res53() - 0.5d0) * 2.0d0 * pi
           theta = (genrand_res53() - 0.5d0) * 2.0d0 * pi

           rxw = genrand_res53()*xcel
           ryw = genrand_res53()*ycel
           rzw = genrand_res53() * water_layer_z + lower_water_z


           rot_m(1,1) =   dcos(psi)*dcos(phi) &
                &       - dcos(theta)*dsin(phi)*dsin(psi)
           rot_m(1,2) =   dcos(psi)*dsin(phi) &
                &       + dcos(theta)*dcos(phi)*dsin(psi)
           rot_m(1,3) =   dsin(psi)*dsin(theta)
           rot_m(2,1) = - dsin(psi)*dcos(phi) &
                &       - dcos(theta)*dsin(phi)*dcos(psi)
           rot_m(2,2) = - dsin(psi)*dsin(phi) &
                &       + dcos(theta)*dcos(phi)*dcos(psi)
           rot_m(2,3) =   dcos(psi)*dsin(theta)
           rot_m(3,1) =   dsin(theta)*dsin(phi)
           rot_m(3,2) = - dsin(theta)*dcos(phi)
           rot_m(3,3) =   dcos(theta)

           do j= j1, j2
              k = molept_list(j)

              atmcor_t(1) = atmcor(1,k)
              atmcor_t(2) = atmcor(2,k)
              atmcor_t(3) = atmcor(3,k)

              atmcor(1,k) =  rot_m(1,1)*atmcor_t(1) &
                   &       + rot_m(2,1)*atmcor_t(2) &
                   &       + rot_m(3,1)*atmcor_t(3) &
                   &       + rxw
              atmcor(2,k) =  rot_m(1,2)*atmcor_t(1) &
                   &       + rot_m(2,2)*atmcor_t(2) &
                   &       + rot_m(3,2)*atmcor_t(3) &
                   &       + ryw
              atmcor(3,k) =  rot_m(1,3)*atmcor_t(1) &
                   &       + rot_m(2,3)*atmcor_t(2) &
                   &       + rot_m(3,3)*atmcor_t(3) &
                   &       + rzw

!---- check the atom is very close to the other atoms
              do m = 1, nmole
                 if (.not. ifsetcor(m) .or. imole == m) cycle

                 n1 = molept_index(m)
                 n2 = molept_index(m+1) - 1
                 do n = n1, n2
                    nn = molept_list(n)

                    r_ij(1:3) = atmcor(1:3,k) - atmcor(1:3,nn)
                    r_ij(1) = r_ij(1) - xcel * dnint(r_ij(1)/xcel)
                    r_ij(2) = r_ij(2) - ycel * dnint(r_ij(2)/ycel)
                    r_ij(3) = r_ij(3) - zcel * dnint(r_ij(3)/zcel)

                    rij = r_ij(1)**2 + r_ij(2)**2 + r_ij(3)**2
                    rij = dsqrt(rij)

                    if (rij < r_min) then
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
                       if (icount > maxicount) then
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
                       goto 100
                    end if
                 end do
              end do

           end do
#ifdef MPI
           if (irank == 0) then
#endif
              write(6,*) 'Set atmcor of mole No. ',imole,',',icount
#ifdef MPI
           end if
#endif
        end do
     end do
  end do

!     ----- Preset the Nose-Hoover chain thermostat -----

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

!----- memory release
  deallocate(base_latt)
  deallocate(base_latt_S1)
  deallocate(base_latt_S2)

!     +     +     +     +     +     +     +

end subroutine createcor
