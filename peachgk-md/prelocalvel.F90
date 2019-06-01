!******************************
!*  prelocalvel.f90 Ver.1.0   *
!*      for peachgk_md.f      *
!*            by G.Kikugawa   *
!******************************
! Time-stamp: <2015-02-12 14:54:42 gota>

!!! Originally coded by Abdul Rafeq Bin Saleman, IFS, Tohoku Univ.
!!! Modified by Gota Kikugawa, IFS, Tohoku Univ.

subroutine prelocalvel(vref, &
     &                 npolytyp,npoly_mole,npoly_atom, &
     &                 nwater, &
     &                 nmatyp,nmatomtyp, &
     &                 polytyp_free,watertyp_free,matomtyp_free, &
     &                 nlvel,index_nlvel,v_nlvel, &
     &                 nlvel_deg_poly, &
     &                 nlvel_deg_water, &
     &                 nlvel_deg_ma)

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  real(8),intent(in):: vref            ! velocity base value [m/s]

  integer,intent(in):: npolytyp        ! number of poly type
  integer,intent(in):: npoly_mole(:)   ! number of molecules of each poly
  integer,intent(in):: npoly_atom(:)   ! number of atoms belonging to poly

  integer,intent(in):: nwater          ! number of H2O molecules

  integer,intent(in):: nmatyp          ! number of species of monatomic mole.
  integer,intent(in):: nmatomtyp(:)    ! each number of monatomic mole.

  character(80),intent(in):: polytyp_free(:,:)
                                ! use for poly type control
  character(80),intent(in):: watertyp_free(:) ! use for water type control
  character(80),intent(in):: matomtyp_free(:,:)
                                ! use for matom type control

!     OUTPUT
  integer,intent(out):: nlvel           ! number of atoms for velocity fix
  integer,intent(out):: index_nlvel(:)  ! index of vel-fix atoms
  real(8),intent(out):: v_nlvel(:,:)    ! local velocity values

  integer,intent(out):: nlvel_deg_poly(:) ! fixed degree of freedom of poly
  integer,intent(out):: nlvel_deg_water ! fixed degree of freedom of H2O
  integer,intent(out):: nlvel_deg_ma(:) ! fixed degree of freedom of matom

! LOCAL:
  integer:: ipoly
  integer:: imatom
  integer:: iword
      
  integer:: index_atom

  integer:: i,j

  logical:: ifsetlocalvel
  integer:: nlvel_lcl,nlvel_reg

!     +     +     +     +     +     +

!---- initialization
  nlvel_deg_poly(1:npolytyp) = 0
  nlvel_deg_water = 0
  nlvel_deg_ma(1:nmatyp) = 0

  v_nlvel(1:3,1:natom) = 0.0d0

!---- make list for fix atoms of poly type

  nlvel = 0
  do ipoly = 1, npolytyp
     ifsetlocalvel = .false.

     do iword = 1, maxnword

        if (polytyp_free(ipoly,iword)(1:8) == 'localvel') then
           index_atom = 0
           do i = 1, ipoly-1
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           if (.not. ifsetlocalvel) then    ! first time for this poly type
              ifsetlocalvel = .true.

              nlvel_lcl = nlvel
              do i = 1, npoly_mole(ipoly)
                 do j = 1, npoly_atom(ipoly)
                    nlvel = nlvel + 1
                    index_atom = index_atom + 1
                    index_nlvel(nlvel) = index_atom
                    nlvel_deg_poly(ipoly) = nlvel_deg_poly(ipoly) + 3
                 end do
              end do

           end if

           nlvel_reg = nlvel_lcl
           do i = 1, npoly_mole(ipoly)
              do j = 1, npoly_atom(ipoly)

                 nlvel_reg = nlvel_reg + 1
                 if (polytyp_free(ipoly,iword)(9:9) == 'x') then
                    read(polytyp_free(ipoly,iword+1),*) v_nlvel(1,nlvel_reg)
                    v_nlvel(1,nlvel_reg) = v_nlvel(1,nlvel_reg) / vref
                 else if (polytyp_free(ipoly,iword)(9:9) == 'y') then
                    read(polytyp_free(ipoly,iword+1),*) v_nlvel(2,nlvel_reg)
                    v_nlvel(2,nlvel_reg) = v_nlvel(2,nlvel_reg) / vref
                 else if (polytyp_free(ipoly,iword)(9:9) == 'z') then
                    read(polytyp_free(ipoly,iword+1),*) v_nlvel(3,nlvel_reg)
                    v_nlvel(3,nlvel_reg) = v_nlvel(3,nlvel_reg) / vref
                 end if

              end do
           end do

        end if

     end do

  end do

!---- make list for fix atoms of water type

  ifsetlocalvel = .false.

  do iword = 1, maxnword

     if (watertyp_free(iword)(1:8) == 'localvel') then

        index_atom = 0
        do i = 1, npolytyp
           index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
        end do

        if (.not. ifsetlocalvel) then    ! first time for this poly type
           ifsetlocalvel = .true.

           nlvel_lcl = nlvel
           do i = 1, nwater
              do j = 1, 3
                 nlvel = nlvel + 1
                 index_atom = index_atom + 1
                 index_nlvel(nlvel) = index_atom
                 nlvel_deg_water = nlvel_deg_water + 3
              end do
           end do

        end if

        nlvel_reg = nlvel_lcl
        do i = 1, nwater
           do j = 1, 3

              nlvel_reg = nlvel_reg + 1
              if (watertyp_free(iword)(9:9) == 'x') then
                 read(watertyp_free(iword+1),*) v_nlvel(1,nlvel_reg)
                 v_nlvel(1,nlvel_reg) = v_nlvel(1,nlvel_reg) / vref
              else if (watertyp_free(iword)(9:9) == 'y') then
                 read(watertyp_free(iword+1),*) v_nlvel(2,nlvel_reg)
                 v_nlvel(2,nlvel_reg) = v_nlvel(2,nlvel_reg) / vref
              else if (watertyp_free(iword)(9:9) == 'z') then
                 read(watertyp_free(iword+1),*) v_nlvel(3,nlvel_reg)
                 v_nlvel(3,nlvel_reg) = v_nlvel(3,nlvel_reg) / vref
              end if

           end do
        end do

     end if

  end do

!---- make list for fix atoms of matom type

  do imatom = 1, nmatyp
     ifsetlocalvel = .false.

     do iword = 1, maxnword

        if (matomtyp_free(imatom,iword)(1:8) == 'localvel') then
           index_atom = 0
           do i = 1, npolytyp
              index_atom = index_atom + npoly_mole(i)*npoly_atom(i)
           end do

           index_atom = index_atom + nwater*3

           do i = 1, imatom - 1
              index_atom = index_atom + nmatomtyp(i)
           end do

           if (.not. ifsetlocalvel) then    ! first time for this poly type
              ifsetlocalvel = .true.

              nlvel_lcl = nlvel
              do i = 1, nmatomtyp(imatom)
                 nlvel = nlvel + 1
                 index_atom = index_atom + 1
                 index_nlvel(nlvel) = index_atom
                 nlvel_deg_ma(imatom) = nlvel_deg_ma(imatom) + 3
              end do

           end if

           nlvel_reg = nlvel_lcl
           do i = 1, nmatomtyp(imatom)

              nlvel_reg = nlvel_reg + 1
              if (matomtyp_free(imatom,iword)(9:9) == 'x') then
                 read(matomtyp_free(imatom,iword+1),*) v_nlvel(1,nlvel_reg)
                 v_nlvel(1,nlvel_reg) = v_nlvel(1,nlvel_reg) / vref
              else if (matomtyp_free(imatom,iword)(9:9) == 'y') then
                 read(matomtyp_free(imatom,iword+1),*) v_nlvel(2,nlvel_reg)
                 v_nlvel(2,nlvel_reg) = v_nlvel(2,nlvel_reg) / vref
              else if (matomtyp_free(imatom,iword)(9:9) == 'z') then
                 read(matomtyp_free(imatom,iword+1),*) v_nlvel(3,nlvel_reg)
                 v_nlvel(3,nlvel_reg) = v_nlvel(3,nlvel_reg) / vref
              end if

           end do

        end if

     end do

  end do

!     +     +     +     +     +     +     +

end subroutine prelocalvel
