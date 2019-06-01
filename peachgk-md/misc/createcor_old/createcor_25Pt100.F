************************************************
*  createcor_25Pt100.f Ver.2.5 '07.02.11       *
*      for peachgk_md.f                        *
*            by G.Kikugawa                     *
************************************************
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
#ifdef MPI
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
      real*8:: rxma,ryma,rzma   ! coordinate of monatomic mole.

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

      integer:: icrepo, icrema

c     parameters for fcc lattice
      real*8:: A0               ! lattice constant

      real*8:: a1(3),a2(3),a3(3) ! lattice vector for fcc (100) surface
      
      real*8:: u(3,4)           ! vector to the particle position in unit lattice

c     parameters for coordinate arrangement
      real*8:: water_layer      ! z length of water layer
      real*8:: upper_Pt_z       ! lower bound of upper Pt layer
      

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

c     set lattice vector
      A0 = 3.9160d-10 / xref

      a1(1) = 1.0d0 * A0
      a1(2) = 0.0d0 * A0
      a1(3) = 0.0d0 * A0

      a2(1) = 0.0d0 * A0
      a2(2) = 1.0d0 * A0
      a2(3) = 0.0d0 * A0

      a3(1) = 0.0d0 * A0
      a3(2) = 0.0d0 * A0
      a3(3) = 1.0d0 * A0

c     vectors in unit lattice
      u(1,1) = 0.0d0
      u(2,1) = 0.0d0
      u(3,1) = 0.0d0

      u(1,2) = 0.5d0*a1(1) + 0.5d0*a2(1)
      u(2,2) = 0.5d0*a1(2) + 0.5d0*a2(2)
      u(3,2) = 0.5d0*a1(3) + 0.5d0*a2(3)

      u(1,3) = 0.5d0*a1(1) + 0.5d0*a3(1)
      u(2,3) = 0.5d0*a1(2) + 0.5d0*a3(2)
      u(3,3) = 0.5d0*a1(3) + 0.5d0*a3(3)

      u(1,4) = 0.5d0*a2(1) + 0.5d0*a3(1)
      u(2,4) = 0.5d0*a2(2) + 0.5d0*a3(2)
      u(3,4) = 0.5d0*a2(3) + 0.5d0*a3(3)

c     define water layer
      water_layer = 4.0d-9 / xref
      upper_Pt_z = (4.0d0 - 0.5d0) * A0 + water_layer

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

c----    arrange MA molecules
         do iz=1,zmaxma(imatom)

            do iy=1,ymaxma(imatom) / 2

               do ix=1,xmaxma(imatom) / 2

                  do i = 1, 4  ! per unit lattice

                     icount = 0
                     imole = imole + 1
                     ifsetcor(imole) = .true.

                     j1 = molept_index(imole)
                     j2 = molept_index(imole+1) - 1

 90                  continue

                     rxma = 0.0d0
                     ryma = 0.0d0
                     rzma = 0.0d0
                     
                     if (imatom == 1) then ! lower Pt layer
                  
                        rxma = u(1,i) + DBLE(ix-1)*a1(1)
     &                                + DBLE(iy-1)*a2(1)
     &                                + DBLE(iz-1)*a3(1)
                        ryma = u(2,i) + DBLE(ix-1)*a1(2)
     &                                + DBLE(iy-1)*a2(2)
     &                                + DBLE(iz-1)*a3(2)
                        rzma = u(3,i) + DBLE(ix-1)*a1(3)
     &                                + DBLE(iy-1)*a2(3)
     &                                + DBLE(iz-1)*a3(3)
                  
                     else if (imatom == 2) then ! upper Pt layer

                        rxma = u(1,i) + DBLE(ix-1)*a1(1)
     &                                + DBLE(iy-1)*a2(1)
     &                                + DBLE(iz-1)*a3(1)
                        ryma = u(2,i) + DBLE(ix-1)*a1(2)
     &                                + DBLE(iy-1)*a2(2)
     &                                + DBLE(iz-1)*a3(2)
                        rzma = u(3,i) + DBLE(ix-1)*a1(3)
     &                                + DBLE(iy-1)*a2(3)
     &                                + DBLE(iz-1)*a3(3)
     &                                + upper_Pt_z

                     end if

                     do j= j1, j2
                        k = molept_list(j)

                        atmcor_t(1) = atmcor(1,k)
                        atmcor_t(2) = atmcor(2,k)
                        atmcor_t(3) = atmcor(3,k)

                        atmcor(1,k) = rxma
                        atmcor(2,k) = ryma
                        atmcor(3,k) = rzma

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
                              r_ij(1) =  r_ij(1) - xcel
     &                                * dnint(r_ij(1)/xcel)
                              r_ij(2) =  r_ij(2) - ycel
     &                                * dnint(r_ij(2)/ycel)
                              r_ij(3) =  r_ij(3) - zcel
     &                                * dnint(r_ij(3)/zcel)

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

                     end do
#ifdef MPI
                     if (irank == 0) then
#endif
                   write(6,*) 'Set atmcor of mole No. ',imole,',',icount
#ifdef MPI
                     end if
#endif
                  end do    ! loop i
      
               end do
            end do
         end do

      end do

c-------- Locate molecules of poly --------

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

         do ix=0,xmaxpo(ipoly)-1
            do iy=0,ymaxpo(ipoly)-1
               do iz=0,zmaxpo(ipoly)-1

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

c                 --- make gravity coordinate ---
                  rxp = compfact*(grnd(mt,mti)-0.5d0)*xcel 
     &                 + 0.5d0*xcel
                  ryp = compfact*(grnd(mt,mti)-0.5d0)*ycel 
     &                 + 0.5d0*ycel
                  rzp = compfact*(grnd(mt,mti)-0.5d0)*zcel 
     &                 + 0.5d0*zcel

c                 --- make rotational matrix ---
                  psi =   (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi
                  phi =   (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi
                  theta = (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi
!               psi   = 0.5d0 * pi
!               phi   = 0.0d0
!               theta = 0.5d0 * pi

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

c                 --- make coordinate ---
                  do j=j1,j2
                     k = molept_list(j)

                     atmcor_t(1) = atmcor(1,k)
                     atmcor_t(2) = atmcor(2,k)
                     atmcor_t(3) = atmcor(3,k)
                     
                     atmcor(1,k) =  rot_m(1,1)*atmcor_t(1)
     &                            + rot_m(2,1)*atmcor_t(2)
     &                            + rot_m(3,1)*atmcor_t(3)
     &                            + rxp
                     atmcor(2,k) =  rot_m(1,2)*atmcor_t(1)
     &                            + rot_m(2,2)*atmcor_t(2)
     &                            + rot_m(3,2)*atmcor_t(3)
     &                            + ryp
                     atmcor(3,k) =  rot_m(1,3)*atmcor_t(1)
     &                            + rot_m(2,3)*atmcor_t(2)
     &                            + rot_m(3,3)*atmcor_t(3)
     &                            + rzp
                     
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
                           r_ij(1) =  r_ij(1) - xcel 
     &                              * dnint(r_ij(1)/xcel)
                           r_ij(2) =  r_ij(2) - ycel
     &                              * dnint(r_ij(2)/ycel)
                           r_ij(3) =  r_ij(3) - zcel
     &                              * dnint(r_ij(3)/zcel)

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
                  write(6,*) 'Set atmcor of mole No. ',imole,',',icount
#ifdef MPI
                  end if
#endif
               end do
            end do
         end do

      end do

c-------- Locate molecules of H2O --------

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
 100           psi =   (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi
               phi =   (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi
               theta = (grnd(mt,mti) - 0.5d0) * 2.0d0 * pi

               rxw = grnd(mt,mti)*xcel
               ryw = grnd(mt,mti)*ycel
!               rzw = grnd(mt,mti)*zcel
               rzw = grnd(mt,mti)*(water_layer - 6.0d0)
     &             + A0*(4.0d0-0.5d0) + 3.0d0

               rot_m(1,1) =   dcos(psi)*dcos(phi) 
     &                      - dcos(theta)*dsin(phi)*dsin(psi)
               rot_m(1,2) =   dcos(psi)*dsin(phi)
     &                      + dcos(theta)*dcos(phi)*dsin(psi)
               rot_m(1,3) =   dsin(psi)*dsin(theta)
               rot_m(2,1) = - dsin(psi)*dcos(phi)
     &                      - dcos(theta)*dsin(phi)*dcos(psi)
               rot_m(2,2) = - dsin(psi)*dsin(phi)
     &                      + dcos(theta)*dcos(phi)*dcos(psi)
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
                        r_ij(1) = r_ij(1) - xcel * dnint(r_ij(1)/xcel)
                        r_ij(2) = r_ij(2) - ycel * dnint(r_ij(2)/ycel)
                        r_ij(3) = r_ij(3) - zcel * dnint(r_ij(3)/zcel)

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
