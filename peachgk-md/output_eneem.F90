!******************************
!*  output_eneem.f90 Ver.1.4  *
!*      for peachgk_md.f      *
!*            by G.Kikugawa   *
!******************************
! Time-stamp: <2015-07-17 16:14:51 gota>

subroutine outene_em(ouene,current_step,xref,eref,tempref,fref, &
     &               npolytyp,nmatyp, &
     &               pot_tot,pot_nonbon,pot_vdw, &
     &               pot_elc,pot_ewk,pot_ewc,pot_ewnq, &
     &               pot_vdw14,pot_elc14, &
     &               pot_bond,pot_angl,pot_anglub, &
     &               pot_tors,pot_torsrb,pot_torsim, &
     &               pot_mor, &
     &               pot_sh, &
     &               pot_rfh, &
     &               pot_dou, &
     &               pot_rpvw,pot_vw, &
     &               pot_cstmnb, &
     &               ncstmnbex,pot_cstmnbex, &
     &               pot_posres, &
     &               pot_pbias, &
     &               ene_kin_poly,ene_kin_water, &
     &               ene_kin_ma,ene_kin_all, &
     &               ene_tot, &
     &               temp_poly,temp_water,temp_ma,temp_all, &
     &               ene_kin_th,ene_pot_th, &
     &               ene_kin_ba,ene_pot_ba, &
     &               ene_conserve, &
     &               d_r,d_ene,drms)

  use interface_tools

  use md_global

  implicit none

!     subroutine for outputting energy data

! ARGUMENTS:
!     INPUT
  integer,intent(in):: ouene           ! output unit for output energy data

  integer,intent(in):: current_step    ! current time step

!---- base value for non-dimensionalize
  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: eref             ! energy base value [J]

  real(8),intent(in):: tempref          ! temperature base value [K]

  real(8),intent(in):: fref             ! force base value [N]

!---- number of type of poly and matom
  integer,intent(in):: npolytyp        ! number of poly type

  integer,intent(in):: nmatyp          ! number of species of monatomic mole.

!---- potential valiables
  real(8),intent(inout):: pot_tot          ! all potential
  real(8),intent(inout):: pot_nonbon       ! all non-bonded potential
  real(8),intent(inout):: pot_vdw          ! vdw potential
  real(8),intent(inout):: pot_elc          ! coulomb potential(ewald real)
  real(8),intent(inout):: pot_ewk          ! coulomb potential(ewald wave)
  real(8),intent(in):: pot_ewc             ! coulomb potential(ewald self)
  real(8),intent(inout):: pot_ewnq         ! coulomb potential(ewald netq)
  real(8),intent(inout):: pot_vdw14        ! 1-4vdw potential
  real(8),intent(inout):: pot_elc14        ! 1-4elc potential
  real(8),intent(inout):: pot_bond         ! bond potential
  real(8),intent(inout):: pot_angl         ! angle potential
  real(8),intent(inout):: pot_anglub       ! Urey-Bradley angle potential
  real(8),intent(inout):: pot_tors         ! torsion potential
  real(8),intent(inout):: pot_torsrb       ! RBtorsion potential
  real(8),intent(inout):: pot_torsim       ! improper torsion potential
  real(8),intent(inout):: pot_mor          ! Morse potential
  real(8),intent(inout):: pot_sh           ! SH potential
  real(8),intent(inout):: pot_rfh          ! RFH potential
  real(8),intent(inout):: pot_dou          ! DOU potential
  real(8),intent(inout):: pot_rpvw         ! RP-VW interaction
  real(8),intent(inout):: pot_vw           ! potential of constant force
  real(8),intent(inout):: pot_cstmnb       ! Custom NB potential
  real(8),intent(inout):: pot_cstmnbex(:)  ! extra Custom NB potential
  real(8),intent(inout):: pot_posres       ! position restraint potential
  real(8),intent(inout):: pot_pbias        ! bias potential

  integer,intent(in):: ncstmnbex           ! number of extra custom NB output

!---- kinematic energy
  real(8),intent(inout):: ene_kin_poly(:)  ! kinetic energy
  real(8),intent(inout):: ene_kin_water    ! kinetic energy
  real(8),intent(inout):: ene_kin_ma(:)    ! kinetic energy
  real(8),intent(inout):: ene_kin_all      ! kinetic energy

  real(8),intent(inout):: ene_kin_th       ! kinetic energy of NHC thermostat
  real(8),intent(inout):: ene_pot_th       ! potential energy of NHC thermostat

  real(8),intent(inout):: ene_kin_ba       ! kinetic energy of Andersen barostat
  real(8),intent(inout):: ene_pot_ba      ! potential energy of Andersen barostat

!---- all energy

  real(8),intent(inout):: ene_tot          ! total Hamiltonian
  real(8),intent(inout):: ene_conserve     ! total conserved quantity

!---- temperature

  real(8),intent(inout):: temp_poly(:)     ! temperature of polymer1
  real(8),intent(inout):: temp_water       ! temperature of H2O
  real(8),intent(inout):: temp_ma(:)       ! temperature of monatomic mole.
  real(8),intent(inout):: temp_all         ! temperature of all atoms

!---- variables for EM
  real(8),intent(in):: d_r                 ! delta r
  real(8),intent(in):: d_ene               ! delta pot_tot
  real(8),intent(in):: drms                ! root mean square force / natom

! LOCAL:

  real(8):: pot_bon          ! total bond energy
  real(8):: pot_coul         ! total coulomb energy
  real(8):: ene_kin          ! total kinetic energy
  real(8):: pot_ewco

  character(80):: sub_tmp
  integer:: sub_len
  character(80):: sub_ene_kin_poly(maxnpolytyp)
  character(80):: sub_ene_kin_ma(maxnpolytyp)
  character(80):: sub_temp_poly(maxnpolytyp)
  character(80):: sub_temp_ma(maxnpolytyp)

  character(80):: sub_pot_cstmnbex(1:ncstmnbex)

  character(25):: fmt
  integer:: nenefield

  integer:: i

!     +     +     +     +     +     +     +

!---- some preparation

  !- count number of output fields used below
  nenefield = 36               ! fixed number of fields
  nenefield = nenefield + ncstmnbex
  nenefield = nenefield + npolytyp
  nenefield = nenefield + nmatyp
  nenefield = nenefield + npolytyp
  nenefield = nenefield + nmatyp

!---- write header

  if (current_step == 0) then

     do i = 1, npolytyp
        sub_ene_kin_poly(i) = ' '
        sub_temp_poly(i) = ' '
        write(sub_tmp,*) 'ene_kin_poly_',i
        call excl_sp(sub_tmp,sub_len,sub_ene_kin_poly(i))
        write(sub_tmp,*) 'temp_poly_',i
        call excl_sp(sub_tmp,sub_len,sub_temp_poly(i))
     end do

     do i = 1, nmatyp
        sub_ene_kin_ma(i) = ' '
        sub_temp_ma(i) = ' '
        write(sub_tmp,*) 'ene_kin_ma_',i
        call excl_sp(sub_tmp,sub_len,sub_ene_kin_ma(i))
        write(sub_tmp,*) 'temp_ma_',i
        call excl_sp(sub_tmp,sub_len,sub_temp_ma(i))
     end do

     !- creating pot_cstmnbex subject
     do i = 1, ncstmnbex
        sub_pot_cstmnbex(i) = ' '
        write(sub_tmp,*) 'pot_cstmnbex_',i
        call excl_sp(sub_tmp,sub_len,sub_pot_cstmnbex(i))
     end do

     !- creating output format
#if defined(_OUTENE_HIGH)
     write(fmt,'(''(A9,1X,''I8''(1X,A23))'')') nenefield
     ! 100  format(A9,1X,<n>(1X,A23))
#else
     write(fmt,'(''(A9,1X,''I8''(1X,A15))'')') nenefield
     ! 100  format(A9,1X,<n>(1X,A15))
#endif

     write(ouene,fmt) 'nstep','ene_conserve','ene_tot','pot_tot', &
          &           'pot_nonbon','pot_bon','pot_coul','pot_vdw', &
          &           'pot_elc','pot_ewk','pot_ewc','pot_ewnq', &
          &           'pot_vdw14','pot_elc14', &
          &           'pot_bond','pot_angl','pot_anglub', &
          &           'pot_tors','pot_torsrb','pot_torsim', &
          &           'pot_mor', &
          &           'pot_sh', &
          &           'pot_rfh', &
          &           'pot_dou', &
          &           'pot_rpvw','pot_vw', &
          &           'pot_cstmnb', &
          &           (sub_pot_cstmnbex(i),i=1,ncstmnbex), &
          &           'pot_posres', &
          &           'pot_pbias', &
          &           'ene_kin', &
          &           (sub_ene_kin_poly(i),i=1,npolytyp), &
          &           'ene_kin_water', &
          &           (sub_ene_kin_ma(i),i=1,nmatyp), &
          &           'ene_kin_th','ene_pot_th', &
          &           'ene_kin_ba','ene_pot_ba', &
          &           (sub_temp_poly(i),i=1,npolytyp), &
          &           'temp_water', &
          &           (sub_temp_ma(i),i=1,nmatyp), &
          &           'temp_all'


  end if

!---- dimensionalize and so on
  ene_tot = ene_tot * eref
  pot_tot = pot_tot * eref
  pot_nonbon = pot_nonbon * eref
  pot_bon =  (pot_bond + pot_angl + pot_anglub &
       &    + pot_tors + pot_torsrb + pot_torsim + pot_mor) * eref
  pot_coul = (pot_elc + pot_ewk + pot_ewc) * eref
  pot_vdw = pot_vdw * eref
  pot_elc = pot_elc * eref
  pot_ewk = pot_ewk * eref
  pot_ewco = pot_ewc * eref
  pot_ewnq = pot_ewnq * eref
  pot_vdw14 = pot_vdw14 * eref
  pot_elc14 = pot_elc14 * eref
  pot_bond = pot_bond * eref
  pot_angl = pot_angl * eref
  pot_anglub = pot_anglub * eref
  pot_tors = pot_tors * eref
  pot_torsrb = pot_torsrb * eref
  pot_torsim = pot_torsim * eref
  pot_mor = pot_mor * eref
  pot_sh = pot_sh * eref
  pot_rfh = pot_rfh * eref
  pot_dou = pot_dou * eref
  pot_rpvw = pot_rpvw * eref
  pot_vw = pot_vw * eref
  pot_cstmnb = pot_cstmnb * eref
  pot_cstmnbex(1:ncstmnbex) = pot_cstmnbex(1:ncstmnbex) * eref
  pot_posres = pot_posres * eref
  pot_pbias = pot_pbias * eref

!      ene_kin = (ene_kin_poly + ene_kin_water + ene_kin_ma) * eref
  ene_kin = 0.0d0
  do i = 1, npolytyp
     ene_kin_poly(i) = ene_kin_poly(i) * eref
     ene_kin = ene_kin + ene_kin_poly(i)
  end do
  ene_kin_water = ene_kin_water * eref
  ene_kin = ene_kin + ene_kin_water
  do i = 1, nmatyp
     ene_kin_ma(i) = ene_kin_ma(i) * eref
     ene_kin = ene_kin + ene_kin_ma(i)
  end do
!      ene_kin = ene_kin * eref
  ene_kin_all = ene_kin_all * eref
  temp_poly(1:npolytyp) = temp_poly(1:npolytyp) * tempref
  temp_water = temp_water * tempref
  temp_ma(1:nmatyp) = temp_ma(1:nmatyp) * tempref
  temp_all = temp_all * tempref

  ene_kin_th = ene_kin_th * eref
  ene_pot_th = ene_pot_th * eref

  ene_kin_ba = ene_kin_ba * eref
  ene_pot_ba = ene_pot_ba * eref

  ene_conserve = ene_conserve * eref

!---- write data

  !- creating output format
#if defined(_OUTENE_HIGH)
  write(fmt,'(''(I10,''I8''(1X,E23.16))'')') nenefield
  ! 101 format(I10,<n>(1X,E23.16))
#else
  write(fmt,'(''(I10,''I8''(1X,E15.8))'')') nenefield
  ! 101 format(I10,<n>(1X,E15.8))
#endif
      
  write(ouene,fmt) current_step,ene_conserve,ene_tot,pot_tot, &
       &           pot_nonbon,pot_bon,pot_coul,pot_vdw, &
       &           pot_elc,pot_ewk,pot_ewco,pot_ewnq, &
       &           pot_vdw14,pot_elc14, &
       &           pot_bond,pot_angl,pot_anglub, &
       &           pot_tors,pot_torsrb,pot_torsim, &
       &           pot_mor, &
       &           pot_sh, &
       &           pot_rfh, &
       &           pot_dou, &
       &           pot_rpvw,pot_vw, &
       &           pot_cstmnb, &
       &           (pot_cstmnbex(i),i=1,ncstmnbex), &
       &           pot_posres, &
       &           pot_pbias, &
       &           ene_kin, &
       &           (ene_kin_poly(i),i=1,npolytyp), &
       &           ene_kin_water, &
       &           (ene_kin_ma(i),i=1,nmatyp), &
       &           ene_kin_th,ene_pot_th, &
       &           ene_kin_ba,ene_pot_ba, &
       &           (temp_poly(i),i=1,npolytyp),temp_water, &
       &           (temp_ma(i),i=1,nmatyp),temp_all

!     output energy
  write(6,'(I10,8E13.5)') current_step,ene_conserve,ene_tot,pot_tot,ene_kin, &
       &                  ene_kin_th,ene_pot_th, &
       &                  ene_kin_ba,ene_pot_ba

!     output temperature

  !- creating output format
  nenefield = 2
  nenefield = nenefield + npolytyp
  nenefield = nenefield + nmatyp
  write(fmt,'(''(10X,''I5''F9.2)'')') nenefield
  ! 103 format(10X,100F9.2)

  write(6,fmt) (temp_poly(i),i=1,npolytyp),temp_water, &
       &       (temp_ma(i),i=1,nmatyp),temp_all

!     output EM information
  write(6,*) '        *** EM information'
  write(6,'(10X,4E13.5)') d_r*xref,pot_tot,d_ene*eref,drms*fref
  write(6,*)

!     +     +     +     +     +     +     +

end subroutine outene_em
