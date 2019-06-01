!**************************************
!*  cell_expand.F90 Ver.1.4 '11.01.20 *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine cell_expand( xref,current_step, &
     &                  xcel,ycel,zcel, &
     &                  yratio,zratio, &
     &                  rcut_book,ifbook,rcut, &
     &                  npoly,nwater,nmatom, &
     &                  r_expand )

  use md_global

#if defined(MPI)
  use mpi_global
#endif

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: xref             ! distanse base value [m]
  integer,intent(in):: current_step    ! current time step

  real(8),intent(in):: yratio           ! y cell ratio of y to x
  real(8),intent(in):: zratio           ! z cell ratio of z to x      

  real(8),intent(in):: rcut_book        ! cut off radius of bookkeeping
  logical,intent(in):: ifbook          ! flag for bookkeeping

  real(8),intent(in):: rcut             ! vdw cutoff length [non-d]

  integer,intent(in):: npoly           ! number of polymer1
  integer,intent(in):: nwater          ! number of H2O molecules
  integer,intent(in):: nmatom          ! number of monatomic molecules

  real(8),intent(in):: r_expand         ! expansion ratio of cell (volume)

!     INPUT&OUTPUT
  real(8),intent(inout):: xcel             ! x cell length[non-d]
  real(8),intent(inout):: ycel             ! y cell length[non-d]
  real(8),intent(inout):: zcel             ! z cell length[non-d]

! LOCAL:
  integer:: i,j,k           ! do loop index
  integer:: j1,j2

  real(8):: vol              ! cell volume

  real(8):: hcelmin          ! minimum among half of (Lx,Ly,Lz)

  real(8):: coeff_expand     ! = r_expand^(1/3)

  real(8):: molemass         ! mass of molecule
  real(8):: molgcor(3)       ! coordinate of molecular center of mass

  integer:: nmolatm
  real(8):: r_atmcor(3,1000) ! relative coordinate of atom to molecular center

!     +     +     +     +     +     +     +     +

!---- change xcel,ycel,zcel

  coeff_expand = r_expand**(1.0d0/3.0d0)

  xcel = xcel*coeff_expand ! xcel=xcel*r_expand^(1/3)
  ycel = yratio*xcel
  zcel = zratio*xcel

  vol = xcel*ycel*zcel
  xlogv = 1.0d0/3.0d0*dlog(vol)

#if defined(MPI)
  if (irank == 0) then
#endif
     write(6,*) 'Change cell volume'
     write(6,'(3E16.8)') xcel*xref,ycel*xref,zcel*xref
#if defined(MPI)
  end if
#endif

!---- if rcut or rcut_book is over half cell length

  hcelmin = dmin1(xcel*0.5d0,ycel*0.5d0,zcel*0.5d0)

  if (rcut > hcelmin) then
     write(6,*) 'Error: rcut is over half cell length'
     write(6,*) 'Timestep= ',current_step,' ','rcut= ',rcut*xref
     stop
  end if

  if ((ifbook) .and. (rcut_book > hcelmin)) then
     write(6,*) 'Error: rcut_book is over half cell length'
     write(6,*) 'Timestep= ',current_step,' ','rcut_book= ',rcut_book*xref
     stop
  end if

!---- change molecular center of mass with change of cell volume

!---- molecule polymer1
!     - calculate center of mass -
  DO i = 1, npoly

     molgcor(1:3) = 0.0d0
     molemass     = 0.0d0

!        - find the first (j1) & last (j2) atoms  of mol(i) -
     j1 = molept_index(i)       ! the first atom of mol(i)  
     j2 = molept_index(i+1) - 1 ! the first atom of mol(i+1)

     do j=j1,j2
        k = molept_list(j)

        molgcor(1:3) = molgcor(1:3) + atmcor(1:3,k)*atmmass(k)
        molemass = molemass + atmmass(k)

     end do
         
     molgcor(1:3) = molgcor(1:3) / molemass

     nmolatm = 0
     do j=j1,j2
        k = molept_list(j)

        nmolatm = nmolatm + 1
        r_atmcor(1:3,nmolatm) = atmcor(1:3,k) - molgcor(1:3)

     end do

!        - translate molecule -
     molgcor(1:3) = molgcor(1:3) * coeff_expand

     nmolatm = 0
     do j=j1,j2
        k = molept_list(j)

        nmolatm = nmolatm + 1
        atmcor(1:3,k) = r_atmcor(1:3,nmolatm) + molgcor(1:3)

     end do

  END DO

!---- molecule water
!     - calculate center of mass -
  DO i = npoly+1, npoly+nwater

     molgcor(1:3) = 0.0d0
     molemass     = 0.0d0

!        - find the first (j1) & last (j2) atoms  of mol(i) -
     j1 = molept_index(i)       ! the first atom of mol(i)  
     j2 = molept_index(i+1) - 1 ! the first atom of mol(i+1)

     do j=j1,j2
        k = molept_list(j)

        molgcor(1:3) = molgcor(1:3) + atmcor(1:3,k)*atmmass(k)
        molemass = molemass + atmmass(k)

     end do
         
     molgcor(1:3) = molgcor(1:3) / molemass

     nmolatm = 0
     do j=j1,j2
        k = molept_list(j)

        nmolatm = nmolatm + 1
        r_atmcor(1:3,nmolatm) = atmcor(1:3,k) - molgcor(1:3)

     end do

!        - translate molecule -
     molgcor(1:3) = molgcor(1:3) * coeff_expand

     nmolatm = 0
     do j=j1,j2
        k = molept_list(j)

        nmolatm = nmolatm + 1
        atmcor(1:3,k) = r_atmcor(1:3,nmolatm) + molgcor(1:3)

     end do

  END DO

!---- monatomic molecule
  DO i = npoly+nwater+1, npoly+nwater+nmatom

!        - find the first (j1) & last (j2) atoms  of mol(i) -
     j1 = molept_index(i)       ! the first atom of mol(i)  
     j2 = molept_index(i+1) - 1 ! the first atom of mol(i+1)

     do j=j1,j2
        k = molept_list(j)

        atmcor(1:3,k) = atmcor(1:3,k) * coeff_expand

     end do
         
  END DO

!     +     +     +     +     +     +     +     +

  return
end subroutine cell_expand
