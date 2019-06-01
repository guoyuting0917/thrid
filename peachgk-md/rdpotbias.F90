!**************************************
!*  rdpotbias.f Ver.1.3 '11.02.11     *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine rdpotbias( ouumbname, &
     &                xref, eref )

  use interface_tools
  
  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: eref             ! energy base value [J]

!     OUTPUT
  character(80),intent(out):: ouumbname ! output file name

! LOCAL:
  character(80):: fredat(20)

  integer:: nword

  integer:: i               ! do loop index

  integer:: iuini
  integer:: ios

  character(6):: head_potbias  ! header of potential type <....>

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!---- intialize variable --------

  head_potbias = ' '

!---- read parameter from MD script (peachgk.ini)
      
  iuini = 99
  open(iuini,file='potbias.ini',status='old',iostat=ios)
  if (ios /= 0) then
     write(6,*) 'Failure in opening file: potbias.ini'
     stop
  end if

  READPBIAS:DO
     call rdfree_w( iuini,20, fredat, nword )

     if (fredat(1)(1:3) == 'END') then
        exit

     else if (fredat(1) == ' ') then
        continue

!        output file name
     else if (fredat(1) == 'ouumbname') then
        ouumbname = fredat(2)

!        specify potential type
     else if (fredat(1) == 'potbias_typ') then
        potbias_typ = fredat(2)
        head_potbias = '<' // potbias_typ // '>'

!        sambling number according to the reaction coordinate
     else if (fredat(1) == 'nint_pbias') then
        read(fredat(2),*) nint_pbias

!        potential parameters for selected potential type
     else if (fredat(1) == head_potbias) then
        do
           call rdfree_w( iuini, 20, fredat, nword )
           if (fredat(1) .eq. ' ') cycle READPBIAS

           do i=2,nword
              read(fredat(i),*) para_potbias(i-1)
           end do

        end do

     end if

  END DO READPBIAS

  close(iuini)

!---- non-dimensionalize parameters

  ! harmonic potential
  if (potbias_typ == 'HARM') then 
     para_potbias(1) = para_potbias(1) / xref
     para_potbias(2) = para_potbias(2) * xref*xref / eref

  ! harmonic potential based on COM
  else if (potbias_typ == 'HMCM') then 
     para_potbias(1) = para_potbias(1) / xref
     para_potbias(2) = para_potbias(2) * xref*xref / eref

  ! soft wall potential
  else if (potbias_typ == 'SWAL') then
     para_potbias(1) = para_potbias(1) / xref
     para_potbias(2) = para_potbias(2) / xref
     para_potbias(3) = para_potbias(3) / eref * xref**3

  ! 2D harmonic potential
  else if (potbias_typ == 'HM2D') then
     para_potbias(1) = para_potbias(1) / xref
     para_potbias(2) = para_potbias(2) / xref
     para_potbias(3) = para_potbias(3) * xref*xref / eref

  ! 2D harmonic potential
  else if (potbias_typ == 'HC2D') then
     para_potbias(1) = para_potbias(1) / xref
     para_potbias(2) = para_potbias(2) / xref
     para_potbias(3) = para_potbias(3) * xref*xref / eref

  end if

#if defined(_POTBIAS_DEBUG)
  write(6,*) '***** potbias debug info *****'
  write(6,*) 'Bias potential type= ', potbias_typ
  write(6,*) 'Sampling number= ',nint_pbias
#endif

!     +     +     +     +     +     +     +

  return

end subroutine rdpotbias
