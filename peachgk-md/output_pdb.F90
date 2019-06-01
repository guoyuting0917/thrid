!*****************************
!*  output_pdb.f Ver.1.4     *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2014-06-24 14:13:16 gota>

subroutine outpdb(oupdb,current_step,   &
     &            npolytyp,npoly_mole,npoly_atom,   &
     &            nwater,   &
     &            nmatyp,nmatomtyp,   &
     &            resname_poly_pdb,   &
     &            resname_water_pdb,   &
     &            resname_matom_pdb)

  use md_global

  implicit none

!     subroutine for outputting energy data

! ARGUMENTS:
!     INPUT
  integer,intent(in):: oupdb           ! output unit for outpdb PDB data

  integer,intent(in):: current_step    ! current time step

  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  character(4),intent(in):: resname_poly_pdb(:) ! residue name for poly
  character(4),intent(in):: resname_water_pdb ! residue name for water
  character(4),intent(in):: resname_matom_pdb(:) ! residue name for matom

! LOCAL:

  integer:: i

  integer:: ipoly,ipmole,ipatom
  integer:: iwater,iwatom
  integer:: ima,imatom

  integer:: index_res
  integer:: index_atom

  character(80):: str_remark
  integer:: len_remark

  character(4):: atmtyp_tmp

#if defined(_PDB_CONECT)
  integer:: nbpair(maxnatom)  ! number of atoms to be connected to atom
  integer:: bpair(12,maxnatom)
                      ! index of atoms to be connected to atom (max. 12 atoms)

  integer:: npairloop         ! number of CONECT lines for an atom

  integer:: k,j
#endif

!     +     +     +     +     +     +     +

!---- Initialization
  index_res = 0
  index_atom = 0

!---- Write REMARK entry
  write(str_remark,'(A6,1X,I8)') 'NSTEP=',current_step
  len_remark = LEN_TRIM(str_remark)

  write(oupdb,'(A6,1X,A)') 'REMARK',str_remark(1:len_remark)

!---- Write polytyp configuration
  do ipoly = 1, npolytyp

     do ipmole = 1, npoly_mole(ipoly)
        index_res = index_res + 1
        do ipatom = 1, npoly_atom(ipoly)
           index_atom = index_atom + 1

           atmtyp_tmp = ' ' // atmtyp(index_atom) // ' '

           ! index_atom < 100000 and index_res < 100000
           if ((index_atom < 100000) .and. (index_res < 100000)) then
              write(oupdb,201) index_atom,atmtyp_tmp, &
                   &           resname_poly_pdb(ipoly),index_res, &
                   &           (atmcor(i,index_atom),i=1,3)
           ! index_atom >= 100000 and index_res < 100000
           else if ((index_atom >= 100000) .and. (index_res < 100000)) then
              write(oupdb,202) mod(index_atom,100000),atmtyp_tmp, &
                   &           resname_poly_pdb(ipoly),index_res, &
                   &           (atmcor(i,index_atom),i=1,3)
           ! index_atom >= 100000 and index_res >= 100000
           else
              write(oupdb,203) mod(index_atom,100000),atmtyp_tmp, &
                   &           resname_poly_pdb(ipoly),mod(index_res,100000), &
                   &           (atmcor(i,index_atom),i=1,3)
           end if

        end do
     end do

  end do

!---- Write water configuration
  do iwater = 1, nwater
     index_res = index_res + 1
     do iwatom = 1, 3
        index_atom = index_atom + 1

        atmtyp_tmp = ' ' // atmtyp(index_atom) // ' '

        ! index_atom < 100000 and index_res < 100000
        if ((index_atom < 100000) .and. (index_res < 100000)) then
           write(oupdb,201) index_atom,atmtyp_tmp, &
                &           resname_water_pdb,index_res, &
                &           (atmcor(i,index_atom),i=1,3)
        ! index_atom >= 100000 and index_res < 100000
        else if ((index_atom >= 100000) .and. (index_res < 100000)) then
           write(oupdb,202) mod(index_atom,100000),atmtyp_tmp, &
                &           resname_water_pdb,index_res, &
                &           (atmcor(i,index_atom),i=1,3)
        ! index_atom >= 100000 and index_res >= 100000
        else
           write(oupdb,203) mod(index_atom,100000),atmtyp_tmp, &
                &           resname_water_pdb,mod(index_res,100000), &
                &           (atmcor(i,index_atom),i=1,3)
        end if

     end do
  end do

!---- Write monatomic configuration
  do ima = 1, nmatyp
     do imatom = 1, nmatomtyp(ima)
        index_res = index_res + 1
        index_atom = index_atom + 1

        atmtyp_tmp = ' ' // atmtyp(index_atom) // ' '

        ! index_atom < 100000 and index_res < 100000
        if ((index_atom < 100000) .and. (index_res < 100000)) then
           write(oupdb,201) index_atom,atmtyp_tmp, &
                &           resname_matom_pdb(ima),index_res, &
                &           (atmcor(i,index_atom),i=1,3)
        ! index_atom >= 100000 and index_res < 100000
        else if ((index_atom >= 100000) .and. (index_res < 100000)) then
           write(oupdb,202) mod(index_atom,100000),atmtyp_tmp, &
                &           resname_matom_pdb(ima),index_res, &
                &           (atmcor(i,index_atom),i=1,3)
        ! index_atom >= 100000 and index_res >= 100000
        else
           write(oupdb,203) mod(index_atom,100000),atmtyp_tmp, &
                &           resname_matom_pdb(ima),mod(index_res,100000), &
                &           (atmcor(i,index_atom),i=1,3)
        end if

     end do
  end do

!---- PDB format for HETATM entry
201 format('HETATM',I5,X,A4,1X,A4,I5,4X,3F8.3)
202 format('HETATM',I5.5,X,A4,1X,A4,I5,4X,3F8.3)
203 format('HETATM',I5.5,X,A4,1X,A4,I5.5,4X,3F8.3)

#if defined(_PDB_CONECT)
!---- construct connect information
  if (natom >= 100000) then
     write(6,*) 'Warning: When natom exceeds 100000, do not output CONECT entry.'
     return
  end if

  nbpair(1:natom) = 0
  bpair(1:12,1:natom) = 0

  do i=1, nbond
     nbpair(ibond(i)) = nbpair(ibond(i)) + 1
     nbpair(jbond(i)) = nbpair(jbond(i)) + 1
     bpair(nbpair(ibond(i)),ibond(i)) = jbond(i)
     bpair(nbpair(jbond(i)),jbond(i)) = ibond(i)
  end do

!---- write CONECT entry
  do i=1, natom
     if (nbpair(i) == 0) cycle   ! no connection

     npairloop = int((nbpair(i)-1)/4) + 1   ! number of CONECT lines for an atom
     do j=1, npairloop
        write(oupdb,fmt='(A6,I5)',advance='no') 'CONECT',i
        do k=1,4   ! in each line, only 4 serials are allowed
           if (bpair(k+(j-1)*4,i) /= 0) then
              write(oupdb,fmt='(I5)',advance='no') bpair(k+(j-1)*4,i)
           end if
        end do
        write(oupdb,*)
     end do
  end do
#endif

!---- write END entry
  write(oupdb,'(3A)') 'END'

!     +     +     +     +     +     +     +

end subroutine outpdb
