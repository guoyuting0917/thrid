!**************************************
!*  calvol_.f Ver.1.1 '05.07.06       *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine calvol( xref,current_step, &
     &             xcel,ycel,zcel, &
     &             yratio,zratio, &
     &             xlogv, &
     &             rcut_book,ifbook,rcut, &
     &             rcutmor,ifbookmor,rcut_bookmor)

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: xref             ! distanse base value [m]
  integer,intent(in):: current_step    ! current time step

  real(8),intent(in):: yratio           ! y cell ratio of y to x
  real(8),intent(in):: zratio           ! z cell ratio of z to x      

  real(8),intent(in):: xlogv            ! epsilon of the barostat coordinate

  real(8),intent(in):: rcut_book        ! cut off radius of bookkeeping
  logical,intent(in):: ifbook          ! flag for bookkeeping

  real(8),intent(in):: rcut             ! vdw cutoff length [non-d]

  real(8),intent(in):: rcutmor          ! Morse cutoff length [m]

  logical,intent(in):: ifbookmor      ! flag for bookkeeping of Morse interaction
  real(8),intent(in):: rcut_bookmor   ! cut off radius of bookkeeping[m] of Morse

!     INPUT&OUTPUT
  real(8),intent(inout):: xcel             ! x cell length[non-d]
  real(8),intent(inout):: ycel             ! y cell length[non-d]
  real(8),intent(inout):: zcel             ! z cell length[non-d]

! LOCAL:
  real(8):: vol              ! cell volume

  real(8):: hcelmin          ! minimum among half of (Lx,Ly,Lz)

!     +     +     +     +     +     +     +     +

!---- Calculate xcel,ycel,zcel
  vol = dexp(3.0d0*xlogv)
  xcel = (vol/yratio/zratio)**(1.0d0/3.0d0) ! xcel=(vol/(yratio*zratio))^1/3
  ycel = yratio*xcel
  zcel = zratio*xcel

!---- if rcut or rcut_book is over half cell length

  hcelmin = dmin1(xcel*0.5d0,ycel*0.5d0,zcel*0.5d0)

  if (rcut > hcelmin) then
     write(6,*) 'Error: rcut is over half cell length'
     write(6,*) 'Timestep= ',current_step,' ','rcut= ',rcut*xref
     stop
  end if

  if ((ifbook) .and. (rcut_book > hcelmin)) then
     write(6,*) 'Error: rcut_book is over half cell length'
     write(6,*) 'Timestep= ',current_step,' ' &
          &     ,'rcut_book= ',rcut_book*xref
     stop
  end if

!---- if rcutmor or rcut_bookmor is over half cell length
  if (rcutmor > hcelmin) then
     write(6,*) 'Error: rcutmor is over half cell length'
     write(6,*) 'Timestep= ',current_step,' ', &
          &     'rcutmor= ',rcutmor*xref
     stop
  end if

  if ((ifbookmor) .and. (rcut_bookmor > hcelmin)) then
     write(6,*) 'Error: rcut_bookmor is over half cell length'
     write(6,*) 'Timestep= ',current_step,' ' &
          &     ,'rcut_bookmor= ',rcut_bookmor*xref
     stop
  end if

!     +     +     +     +     +     +     +     +

  return
end subroutine calvol
