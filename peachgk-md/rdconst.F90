!*****************************
!*  rdconst.f Ver.1.4        *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2014-06-13 15:30:25 gota>

subroutine rdconst( iuparaconst )

  use interface_tools

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: iuparaconst     ! input const parameter file unit

! LOCAL:
  character(80):: fredat(20)

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!-------- Read parameter of intramolecular constraint --------

  nconsttyp = 0
  nconsttyp_p = 0
  nconsttyp_w = 0

  READCONST:DO
     call rdfree( iuparaconst, 20, fredat )

     if (fredat(1) == '<END>') then
        exit
     else if (fredat(1) == ' ') then
        cycle READCONST

     else if (fredat(1) == '<CONST_POLY>') then
        do
           call rdfree( iuparaconst, 20, fredat )
           if (fredat(1) == ' ') cycle READCONST
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nconsttyp = nconsttyp + 1
           nconsttyp_p = nconsttyp_p + 1
           para_consttyp(nconsttyp) = fredat(2)(1:5)

        end do

     else if (fredat(1) == '<CONST_WATER>') then
        do
           call rdfree(iuparaconst,20,fredat)
           if (fredat(1) == ' ') cycle READCONST
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nconsttyp = nconsttyp + 1
           nconsttyp_w = nconsttyp_w + 1
           para_consttyp(nconsttyp) = fredat(2)(1:5)

        end do

     end if

  END DO READCONST

  close(iuparaconst)

!     +     +     +     +     +     +     +

  return
end subroutine rdconst
