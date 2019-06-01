!**************************************
!*  prepdbout.f Ver.1.2 '10.06.30     *
!*      for peachgk_md.f              *
!*            by G.Kikugawa           *
!**************************************
subroutine prepdbout( npolytyp, polytyp_free, &
     &                watertyp_free, &
     &                nmatyp, matomtyp_free, monoatmtyp, &
     &                resname_poly_pdb, &
     &                resname_water_pdb, &
     &                resname_matom_pdb)

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: npolytyp        ! number of poly type
  character(80),intent(in):: polytyp_free(:,:)

  character(80),intent(in):: watertyp_free(:) ! use for water type control

  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  character(80),intent(in):: matomtyp_free(:,:)

  character(2),intent(in):: monoatmtyp(:) ! monatomic mole. type

!     OUTPUT
  character(4),intent(out):: resname_poly_pdb(:)
  character(4),intent(out):: resname_water_pdb ! residue name for poly
  character(4),intent(out):: resname_matom_pdb(:) ! residue name for matom

! LOCAL:
  integer:: ipoly
  integer:: iword
  integer:: imatom
      
!     +     +     +     +     +     +

!---- set resname for poly type

  do ipoly = 1, npolytyp
 
!        input default strings
     write(resname_poly_pdb(ipoly),'(A2,I2.2)') 'PO',ipoly

!        input user defined resname
     do iword = 1, maxnword
        if (polytyp_free(ipoly,iword) == 'pdbresname') then
           resname_poly_pdb(ipoly) = polytyp_free(ipoly,iword+1)(1:4)
           exit
        end if
     end do

  end do

!---- set resname for water type

!     input default strings
  resname_water_pdb = 'WAT '

!        input user defined resname
  do iword = 1, maxnword
     if (watertyp_free(iword) == 'pdbresname') then
        resname_water_pdb = watertyp_free(iword+1)(1:4)
        exit
     end if
  end do

!---- set resname for matom type

  do imatom = 1, nmatyp
 
!        input default strings
     resname_matom_pdb(imatom) = monoatmtyp(imatom) // '  '

!        input user defined resname
     do iword = 1, maxnword
        if (matomtyp_free(imatom,iword) == 'pdbresname') then
           resname_matom_pdb(imatom) = matomtyp_free(imatom,iword+1)(1:4)
           exit
        end if
     end do

  end do


!     +     +     +     +     +     +     +

  return
end subroutine prepdbout
