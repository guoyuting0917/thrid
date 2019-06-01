!**************************************
!*  output_umb.f90 Ver.1.1 '10.07.09  *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine outumb( ouumb, &
     &             xref,eref )

  use md_global

  implicit none

!     subroutine for outputting bias potential data

! ARGUMENTS:
!     INPUT
  integer,intent(in):: ouumb         ! output unit for outumb bias potential data

!---- base value for non-dimensionalize
  real(8),intent(in):: xref             ! distanse base value [m]

  real(8),intent(in):: eref             ! energy base value [J]

! LOCAL:
  integer:: i

  real(8):: d_coord          ! Delta coordinate
  real(8):: coord            ! position
  real(8):: pot_bias         ! bias potential

!     +     +     +     +     +     +     +

!---- write header

  write(ouumb,'(A8,A16,A16)') '#index','coord','pot_bias'

!---- write data
      
  if (potbias_typ == 'HARM') then

     d_coord = (para_potbias(4) - para_potbias(3)) &
          &  / DBLE(nint_pbias)

     do i = 1, nint_pbias
        coord = para_potbias(3) + (DBLE(i-1) + 0.5d0) * d_coord
        pot_bias = para_potbias(2) * (coord - para_potbias(1))**2

        coord = coord * xref
        pot_bias = pot_bias * eref
            
        write(ouumb,'(I8,E16.8,E16.8)') i, coord, pot_bias
     end do

  else if (potbias_typ == 'SWAL') then

     d_coord = (para_potbias(2) - para_potbias(1)) &
          &  / DBLE(nint_pbias)

     do i = 1, nint_pbias
        coord = para_potbias(1) + (DBLE(i-1) + 0.5d0) * d_coord
        pot_bias = 0.0d0

        coord = coord * xref
        pot_bias = pot_bias * eref
            
        write(ouumb,'(I8,E16.8,E16.8)') i, coord, pot_bias
     end do

  else if (potbias_typ == 'HM2D') then
! No execution
  end if

!     +     +     +     +     +     +     +

  return
end subroutine outumb
