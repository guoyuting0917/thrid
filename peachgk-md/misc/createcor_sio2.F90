!*************************************
!*  createcor.f Ver.2.8 '10.12.11    *
!*      for peachgk_md.f             *
!*            by T.Nakano            *
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

  real(8):: box(3)                     ! cell length
  real(8):: atmcor_maxz                 ! 
  real(8):: atmcor_minz                 ! 

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

!-------- Locate atoms of SiO2 solid --------

  atmcor_maxz = 0.0d0
  atmcor_minz = 0.0d0

  ipoly = 1
  j1 = 1
  j2 = npoly_atom(ipoly)*npoly_mole(ipoly)

  do j = j1, j2
     atmcor_maxz = max(atmcor_maxz,atmcor(3,j))
     atmcor_minz = min(atmcor_minz,atmcor(3,j))
     if (atmcor(3,j) < 0.d0) then
       atmcor(3,j) = atmcor(3,j) + zcel
     end if
  end do

!-------- Locate molecules of poly --------

  box(1) = xcel
  box(2) = ycel
  box(3) = zcel - (atmcor_maxz + abs(atmcor_minz) + 3.0d0)

  do icrepo = 2, ncrecorpo
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
110           rxp = compfact*(genrand_res53()-0.5d0)*box(1) + 0.5d0*box(1)
              ryp = compfact*(genrand_res53()-0.5d0)*box(2) + 0.5d0*box(2)
              rzp = compfact*(genrand_res53()-0.5d0)*box(3) + 0.5d0*box(3)

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
                       r_ij(1) = r_ij(1) - box(1) * anint(r_ij(1)/box(1))
                       r_ij(2) = r_ij(2) - box(2) * anint(r_ij(2)/box(2))
                       r_ij(3) = r_ij(3) - box(3) * anint(r_ij(3)/box(3))

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

              rxw = genrand_res53()*box(1)
              ryw = genrand_res53()*box(2)
              rzw = genrand_res53()*box(3)

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
                       r_ij(1) = r_ij(1) - box(1) * anint(r_ij(1)/box(1))
                       r_ij(2) = r_ij(2) - box(2) * anint(r_ij(2)/box(2))
                       r_ij(3) = r_ij(3) - box(3) * anint(r_ij(3)/box(3))

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

  do icrema = 1, ncrecorma
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

              rxma = genrand_res53()*box(1)
              ryma = genrand_res53()*box(2)
              rzma = genrand_res53()*box(3)

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
                       r_ij(1) = r_ij(1) - box(1) * anint(r_ij(1)/box(1))
                       r_ij(2) = r_ij(2) - box(2) * anint(r_ij(2)/box(2))
                       r_ij(3) = r_ij(3) - box(3) * anint(r_ij(3)/box(3))

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

  ipoly = 1
  j1 = npoly_atom(ipoly)*npoly_mole(ipoly) + 1
  j2 = natom
  atmcor(3,j1:j2) = atmcor(3,j1:j2) + (atmcor_maxz + 1.5d0)

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
