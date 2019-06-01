!*************************************
!*  mkconst.f Ver.1.6 '14.03.28      *
!*      for peachgk_md.f             *
!*            by G.Kikugawa          *
!*************************************
subroutine mkconst(npolytyp, npoly_mole)

  use interface_tools

  use md_global
#if defined(MPI)
  use mpi_global
#endif

  implicit none

! ARGUMENTS:
!     INPUT
  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly

! LOCAL
  integer:: i,j             ! do loop index

  integer:: i1,i2,ii
  integer:: m
  integer:: iatm
  integer:: nmconst_tmp
  logical:: ifconstmol
  integer:: list_mconst_tmp(maxnatom)
  integer:: ipolytyp
  integer:: imole
  integer:: index_mole, index_mole2

!     +     +     +     +     +     +

!---- link const

  nconst = 0
!      nconst_poly = 0
  do i=1,npolytyp
     nconst_poly(i) = 0
  end do
  nconst_water = 0

  do i = 1, nbond
     do j = 1, nconsttyp_p

        if (bondtyp(i) == para_consttyp(j)) then
           nconst = nconst + 1
           if (nconst > maxnconst) then
              write(6,*) 'Error : nconst exceeds maxnconst= ',maxnconst
              stop
           end if

!               nconst_poly = nconst_poly + 1
           iconst(nconst) = ibond(i)
           jconst(nconst) = jbond(i)
           dconst(nconst) = para_eqbond(indexbondtyp(i))

!              count nconst_poly of each polytyp
           index_mole = irmolept_list(ibond(i))
           index_mole2 = 0

           do ipolytyp=1,npolytyp
              do imole=1,npoly_mole(ipolytyp)
                 index_mole2 = index_mole2 + 1
                 if (index_mole == index_mole2) then
                    nconst_poly(ipolytyp) = nconst_poly(ipolytyp)+1
                 end if
              end do
           end do

           if (para_cbond(indexbondtyp(i)) /= 0.0d0) then
              write(6,*) 'Error: para_cbond of constraint pair is', &
                   &     ' not zero.'
              stop
           end if

        end if

     end do
  end do

  do i = 1, nbond
     do j = 1, nconsttyp_w

        if (bondtyp(i) == para_consttyp(j+nconsttyp_p)) then
           nconst = nconst + 1
           if (nconst > maxnconst) then
              write(6,*) 'Error : nconst exceeds maxnconst= ',maxnconst
              stop
           end if

           nconst_water = nconst_water + 1
           iconst(nconst) = ibond(i)
           jconst(nconst) = jbond(i)
           dconst(nconst) = para_eqbond(indexbondtyp(i))

           if (para_cbond(indexbondtyp(i)) /= 0.0d0) then
              write(6,*) 'Error: para_cbond of constraint pair is', &
                   &     ' not zero.'
              stop
           end if

        end if

     end do
  end do

!---- make const-list based on molecules
!     This list has the nconst of each molecules and used for rattle_mpi
!     --- initialize list ---
  do i = 1, maxnatom+1
     index_mconst(i) = 0
  end do
  do i = 1, maxnbond
     list_mconst(i) = 0
  end do

!     --- make list_mconst ---
  nmconstindex = 0
  nmconst = 0
  do i = 1, nmoleptindex - 1

     i1 = molept_index(i)
     i2 = molept_index(i+1) - 1

     ifconstmol = .false.
     nmconst_tmp = 0
     do ii = i1, i2
        iatm = molept_list(ii)
        do j = 1, nconst
           if ((iconst(j) == iatm) .or. (jconst(j) == iatm)) then
              ifconstmol = .true.
              nmconst_tmp = nmconst_tmp + 1
              if (nmconst_tmp > maxnatom) then
                 write(6,*) 'Error: nmconst_tmp exceed maxnatom', &
                      &     ' in mkconst'
                 stop
              end if
              list_mconst_tmp(nmconst_tmp) = j
           end if
        end do
     end do
         
     if (ifconstmol) then
        call isort( list_mconst_tmp, nmconst_tmp )
        call iskip( list_mconst_tmp, nmconst_tmp )

        nmconstindex = nmconstindex + 1
        index_mconst(nmconstindex) = nmconst + 1

        do m = 1, nmconst_tmp
           list_mconst(nmconst+m) = list_mconst_tmp(m)
        end do

        nmconst = nmconst + nmconst_tmp

     end if

  end do
      
  index_mconst(nmconstindex+1) = nmconst + 1

!---- notation of nconst

#if defined(MPI)
  if (irank == 0) then
#endif
     write(6,999) nconst
#if defined(MPI)
  end if
#endif

999 format(/3x,'#const: ', i8) 

!     +     +     +     +     +     +

  return
end subroutine mkconst
