!*****************************
!*  output_data.f90 Ver.4.4  *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2015-07-17 16:13:24 gota>

!----------------------------------------------------------
subroutine outene(ouene,current_step,eref,tempref, &
     &            npolytyp,nmatyp, &
     &            pot_tot,pot_nonbon,pot_vdw, &
     &            pot_elc,pot_ewk,pot_ewc,pot_ewnq, &
     &            pot_vdw14,pot_elc14, &
     &            pot_bond,pot_angl,pot_anglub, &
     &            pot_tors,pot_torsrb,pot_torsim, &
     &            pot_mor, &
     &            pot_sh, &
     &            pot_rfh, &
     &            pot_dou, &
     &            pot_rpvw,pot_vw, &
     &            pot_cstmnb, &
     &            ncstmnbex,pot_cstmnbex, &
     &            pot_posres, &
     &            pot_pbias, &
     &            ene_kin_poly,ene_kin_water, &
     &            ene_kin_ma,ene_kin_all, &
     &            ene_tot, &
     &            temp_poly,temp_water,temp_ma,temp_all, &
     &            ene_kin_th,ene_pot_th, &
     &            ene_kin_ba,ene_pot_ba, &
     &            ene_conserve)

  use interface_tools

  use md_global

  implicit none

!     subroutine for outputting energy data

! ARGUMENTS:
!     INPUT
  integer,intent(in):: ouene           ! output unit for output energy data

  integer,intent(in):: current_step    ! current time step

!---- base value for non-dimensionalize
  real(8),intent(in):: eref             ! energy base value [J]

  real(8),intent(in):: tempref          ! temperature base value [K]

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
  real(8),intent(inout):: ene_pot_ba       ! potential energy of Andersen barostat

!---- all energy

  real(8),intent(inout):: ene_tot          ! total Hamiltonian
  real(8),intent(inout):: ene_conserve     ! total conserved quantity

!---- temperature

  real(8),intent(inout):: temp_poly(:)     ! temperature of polymer1
  real(8),intent(inout):: temp_water       ! temperature of H2O
  real(8),intent(inout):: temp_ma(:)       ! temperature of monatomic mole.
  real(8),intent(inout):: temp_all         ! temperature of all atoms

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
  nenefield = 36+5               ! fixed number of fields
  nenefield = nenefield + ncstmnbex
  nenefield = nenefield + npolytyp
  nenefield = nenefield + nmatyp
  nenefield = nenefield + npolytyp
  nenefield = nenefield + nmatyp

!---- write header

  if (current_step == 1) then

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
        write(sub_tmp,*) 'pot_cstmex_',i
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
          &           'temp_all',&             
          &           'wal_int_1_a1', &
          &           'wal_int_1_a2', &
          &           'wal_int_2_a1', &
          &           'wal_int_2_a2', &
          &           'wal_int_n'             



  end if

!---- dimensionalize and so on
  ene_tot = ene_tot * eref
  pot_tot = pot_tot * eref
  pot_nonbon = pot_nonbon * eref
  pot_bon =  (pot_bond + pot_angl + pot_anglub &
  &         + pot_tors + pot_torsrb + pot_torsim + pot_mor) * eref
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


  wall_pot_interface_1_ar1 =  wall_pot_interface_1_ar1 * eref                                    !variable for specific solid-liquid interface-by guo
  wall_pot_interface_1_ar2 =  wall_pot_interface_1_ar2 * eref 
  wall_pot_interface_2_ar1 =  wall_pot_interface_2_ar1 * eref 
  wall_pot_interface_2_ar2 =  wall_pot_interface_2_ar2 * eref 
  wall_pot_interface_no = wall_pot_interface_no * eref

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
       &           (temp_ma(i),i=1,nmatyp),temp_all, &
       &           wall_pot_interface_1_ar1, &
       &           wall_pot_interface_1_ar2, &
       &           wall_pot_interface_2_ar1, &
       &           wall_pot_interface_2_ar2, &
       &           wall_pot_interface_no            

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


!     +     +     +     +     +     +     +

end subroutine outene

!----------------------------------------------------------
subroutine outpos(oupos,current_step,xref)

  use md_global

  implicit none

!     subroutine for outputting trajectory data

! ARGUMENTS:
!     INPUT
  integer,intent(in):: oupos           ! output unit for output position data

  integer,intent(in):: current_step    ! current time step

!---- base value for non-dimensionalize
  real(8),intent(in):: xref             ! distanse base value [m]

! LOCAL:
  integer:: i,m             ! do loop index

#if defined(_DOUBLE_OUTPOS)
  real(8),allocatable:: atmcor_out(:,:)   ! output type for double precision
#else
  real,allocatable:: atmcor_out(:,:)   ! single precision is default
#endif

!     +     +     +     +     +     +     +     +

  !---- allocalte memory
  allocate(atmcor_out(3,natom))

  !---- coordinates are copied to arrays for output

  atmcor_out(1:3,1:natom) = atmcor(1:3,1:natom)*xref

  !---- output coordinate
  write(oupos) current_step

  do i=1,natom
     write(oupos) i,(atmcor_out(m,i),m=1,3),atmtyp(i)
  end do

  !---- release memory
  deallocate(atmcor_out)

! 103  format(I5,3E16.8,A3)

!     +     +     +     +     +     +     +     +

end subroutine outpos

!----------------------------------------------------------
subroutine outvel(ouvel,current_step,vref)

  use md_global

  implicit none

!     subroutine for outputting velocity data

!ARGUMENTS:
!     INPUT
  integer,intent(in):: ouvel           ! output unit for output velocity data

  integer,intent(in):: current_step    ! current time step

!---- base value for non-dimensionalize
  real(8),intent(in):: vref             ! velocity base value [m/s]

! LOCAL:
  integer:: i,m             ! do loop index

#if defined(_DOUBLE_OUTVEL)
  real(8),allocatable:: atmvel_out(:,:)   ! output type for double precision
#else
  real,allocatable:: atmvel_out(:,:)   ! single precision is default
#endif

!     +     +     +     +     +     +     +     +

  !---- allocalte memory
  allocate(atmvel_out(3,natom))

  !---- velocities are copied to arrays for output

  atmvel_out(1:3,1:natom) = atmvel(1:3,1:natom)*vref

  !---- output velocity
  write(ouvel) current_step

  do i=1,natom
     write(ouvel) i,(atmvel_out(m,i),m=1,3),atmtyp(i)
  end do

  !---- release memory
  deallocate(atmvel_out)

! 104  format(I5,3E16.8,A3)

!     +     +     +     +     +     +     +     +

end subroutine outvel

!----------------------------------------------------------
subroutine outthe(outhe,current_step,timeref, &
     &            mchain,xlogs,vlogs)

  implicit none

!     subroutine for outputting thermostat data

! ARGUMENTS:
!     INPUT
  integer,intent(in):: outhe           ! output unit for output thermostat data

  integer,intent(in):: current_step    ! current time step

  real(8),intent(in):: timeref          ! time base value [sec]

!---- thermostat data
  integer,intent(in):: mchain          ! Nose-Hoover chain number (>1 in NVT)
  real(8),intent(in):: xlogs(:)         ! xi of the thermostat coordinate
  real(8),intent(in):: vlogs(:)         ! vxi of the thermostat velocity

! LOCAL:
  integer:: i               ! do loop index

!     +     +     +     +     +     +     +     +

  write(outhe,107) current_step,(xlogs(i),i=1,mchain), &
       &                        (vlogs(i)/timeref,i=1,mchain)

107 format(I10,30E16.8)

!     +     +     +     +     +     +     +     +

end subroutine outthe

!----------------------------------------------------------
subroutine outbar(oubar,current_step,xref,timeref, &
     &            xcel,ycel,zcel,xlogv,vlogv,xboxh,vboxg, &
     &            pcont_axis)

  implicit none

!     subroutine for outputting barostat data

! ARGUMENTS:
!     INPUT
  integer,intent(in):: oubar           ! output unit for output barostat data

  integer,intent(in):: current_step    ! current time step

  real(8),intent(in):: xref             ! distanse base value [m]
  real(8),intent(in):: timeref          ! time base value [sec]

!---- barostat data
  real(8),intent(in):: xcel             ! x cell length[non-d]
  real(8),intent(in):: ycel             ! y cell length[non-d]
  real(8),intent(in):: zcel             ! z cell length[non-d]

  real(8),intent(in):: xlogv            ! epsilon of the barostat coordinate
  real(8),intent(in):: vlogv            ! vepsilon of the barostat velocity

  real(8),intent(in):: xboxh(:)         ! cell vector h of the barostat coordinate
  real(8),intent(in):: vboxg(:)         ! vg vector of the barostat coordinate

  character(5),intent(in):: pcont_axis  ! axis for pressure control (iso, aniso, etc.)

! LOCAL:
  real(8):: vol              ! volume

!     +     +     +     +     +     +     +     +

!---- some preparation
  vol = xcel * ycel * zcel *xref*xref*xref

!---- output barostat data
  if (current_step == 1) then
     if (pcont_axis == 'iso') then   ! isotropic NPT
        write(oubar,'(A8,6A16)') 'nstep','volume','xcel','ycel','zcel' &
        &                  ,'xlogv','vlogv'
    else   ! anisotropic NPT
        write(oubar,'(A8,7A16)') 'nstep','volume','xcel','ycel','zcel', &
        &                        'vboxg_x','vboxg_y','vboxg_z'
    end if
  end if

  if (pcont_axis == 'iso') then   ! isotropic NPT
      write(oubar,108) current_step,vol,xcel*xref,ycel*xref,zcel*xref, &
      &                xlogv,vlogv/timeref
  else   ! anisotropic NPT
      write(oubar,109) current_step,vol,xcel*xref,ycel*xref,zcel*xref, &
      &                vboxg(1)/timeref,vboxg(2)/timeref,vboxg(3)/timeref
  end if

108 format(I10,6E16.8)
109 format(I10,7E16.8)

!     +     +     +     +     +     +     +     +

end subroutine outbar

!----------------------------------------------------------
subroutine outpre(oupre,current_step,pref,eref, &
     &            pressmol_ktot,pressmolt_ktot, &
     &            pressmol_vtot, pressmolt_vtot, &
     &            pressmol_tot,pressmolt_tot, &
     &            pressatm_ktot,pressatmt_ktot, &
     &            pressatm_vtot,pressatmt_vtot, &
     &            pressatm_tot,pressatmt_tot, &
     &            pot_viric_all,pot_virict_all, &
     &            pot_virilj_all,pot_viriljt_all, &
     &            pot_virimor_all,pot_virimort_all, &
     &            pot_virish_all,pot_virisht_all, &
     &            pot_virirfh_all,pot_virirfht_all, &
     &            pot_viridou_all,pot_viridout_all, &
     &            pot_viricstmnb_all,pot_viricstmnbt_all, &
     &            ncstmnbex, &
     &            pot_viricstmnbex_all,pot_viricstmnbext_all, &
     &            pot_viri_all,pot_virit_all, &
     &            viri_fdotd,virit_fdotd, &
     &            atm_viribo_all,atm_viribot_all, &
     &            atm_virian_all,atm_viriant_all, &
     &            atm_virito_all,atm_viritot_all, &
     &            atm_viri14_all,atm_viri14t_all, &
     &            atm_viri_const,atm_virit_const, &
     &            atm_viri_corr,atm_virit_corr)

  use interface_tools

  implicit none

!     subroutine for outputting pressure data

! ARGUMENTS:
!     INPUT
  integer,intent(in):: oupre           ! output unit for outpre velocity data

  integer,intent(in):: current_step    ! current time step

!---- base value for non-dimensionalize
  real(8),intent(in):: pref             ! pressure base value [Pa]
  real(8),intent(in):: eref             ! energy base value [J]

!---- valiables for pressure
  real(8),intent(in):: pressmol_ktot       ! pressure kinetic part of molecule
  real(8),intent(in):: pressmolt_ktot(:,:)
                    ! pressure tensor kinetic part of molecule
  real(8),intent(in):: pressmol_vtot       ! pressure virial part of molecule
  real(8),intent(in):: pressmolt_vtot(:,:)
                    ! pressure tensor virial part of molecule
  real(8),intent(in):: pressmol_tot        ! pressure total of molecule
  real(8),intent(in):: pressmolt_tot(:,:)  ! pressure tensor total of molecule
  real(8),intent(in):: pressatm_ktot       ! pressure kinetic part of atm
  real(8),intent(in):: pressatmt_ktot(:,:) ! pressure tensor kinetic part of atm
  real(8),intent(in):: pressatm_vtot       ! pressure virial part of atm
  real(8),intent(in):: pressatmt_vtot(:,:) ! pressure tensor virial part of atm
  real(8),intent(in):: pressatm_tot        ! pressure total of atm
  real(8),intent(in):: pressatmt_tot(:,:)  ! pressure tensor total of atm

  real(8),intent(in):: pot_viric_all        ! virial(coulomb potential)
  real(8),intent(in):: pot_virict_all(:,:)  ! virial tensor (coulomb potential)
  real(8),intent(in):: pot_virilj_all       ! virial(L-J potential)
  real(8),intent(in):: pot_viriljt_all(:,:) ! virial tensor (L-J potential)
  real(8),intent(in):: pot_virimor_all       ! virial(Morse potential)
  real(8),intent(in):: pot_virimort_all(:,:) ! virial tensor (Morse potential)
  real(8),intent(in):: pot_virish_all       ! virial(SH potential)
  real(8),intent(in):: pot_virisht_all(:,:) ! virial tensor (SH potential)
  real(8),intent(in):: pot_virirfh_all       ! virial(RFH potential)
  real(8),intent(in):: pot_virirfht_all(:,:) ! virial tensor (RFH potential)
  real(8),intent(in):: pot_viridou_all       ! virial(DOU potential)
  real(8),intent(in):: pot_viridout_all(:,:) ! virial tensor (DOU potential)
  real(8),intent(in):: pot_viricstmnb_all    ! virial(custom NB potential)
  real(8),intent(in):: pot_viricstmnbt_all(:,:)
                                          ! virial tensor (custom NB potential)
  real(8),intent(in):: pot_viricstmnbex_all(:)
                                     ! virial(extra custom NB potential)
  real(8),intent(in):: pot_viricstmnbext_all(:,:,:)
                                     ! virial tensor (extra custom NB potential)

  real(8),intent(in):: pot_viri_all         ! all virial of potential term
  real(8),intent(in):: pot_virit_all(:,:)   ! all virial tensor of potential term
  real(8),intent(in):: viri_fdotd           ! virial of correction term F.d
  real(8),intent(in):: virit_fdotd(:,:)    ! virial tensor of correction term F.d

  real(8),intent(in):: atm_viribo_all       ! virial(bond potential) of each atom
  real(8),intent(in):: atm_viribot_all(:,:) ! virial(bond potential) of each atom
  real(8),intent(in):: atm_virian_all      ! virial(angle potential) of each atom
  real(8),intent(in):: atm_viriant_all(:,:)
                    ! virial(angle potential) of each atom
  real(8),intent(in):: atm_virito_all    ! virial(torsion potential) of each atom
  real(8),intent(in):: atm_viritot_all(:,:)
                    ! virial(torsion potential) of each atom
  real(8),intent(in):: atm_viri14_all  ! virial(1-4 force potential) of each atom
  real(8),intent(in):: atm_viri14t_all(:,:)
                    ! virial(1-4 force potential) of each atom

  real(8),intent(in):: atm_viri_const     ! virial(constraint force) of each atom
  real(8),intent(in):: atm_virit_const(:,:)
                    ! virial tensor (constraint) of each atom

  real(8),intent(in):: atm_viri_corr          ! virial(L-J correction)
  real(8),intent(in):: atm_virit_corr(:,:)    ! virial tensor (L-J correction)

  integer,intent(in):: ncstmnbex           ! number of extra custom NB output

! LOCAL:
  real(8):: pressmol_trPdiv3     ! tr(Patm)/3
  real(8):: pressmol_trK         ! tr(Katm)
  real(8):: pressmol_trV         ! tr(Vatm)

  real(8):: pressatm_trPdiv3     ! tr(Patm)/3
  real(8):: pressatm_trK         ! tr(Katm)
  real(8):: pressatm_trV         ! tr(Vatm)

  character(80):: sub_tmp
  integer:: sub_len
  character(80):: sub_pre_cstmnbex(1:ncstmnbex)

  character(25):: fmt
  integer:: nprefield

  integer:: i

!     +     +     +     +     +     +     +

!---- some preparation

  !- count number of output fields used below
  nprefield = 81               ! fixed number of fields
  nprefield = nprefield + ncstmnbex

!     --- Calculate isotropic pressure from pressure tensor ---

  pressmol_trPdiv3 =  1.0d0/3.0d0 * (  pressmolt_tot(1,1) &
       &                             + pressmolt_tot(2,2) &
       &                             + pressmolt_tot(3,3))
  pressmol_trK     =  pressmolt_ktot(1,1) &
       &            + pressmolt_ktot(2,2) &
       &            + pressmolt_ktot(3,3)
  pressmol_trV     =  pressmolt_vtot(1,1) &
       &            + pressmolt_vtot(2,2) &
       &            + pressmolt_vtot(3,3)

  pressatm_trPdiv3 = 1.0d0/3.0d0 * (  pressatmt_tot(1,1) &
       &                            + pressatmt_tot(2,2) &
       &                            + pressatmt_tot(3,3))
  pressatm_trK     =  pressatmt_ktot(1,1) &
       &            + pressatmt_ktot(2,2) &
       &            + pressatmt_ktot(3,3)
  pressatm_trV     =  pressatmt_vtot(1,1) &
       &            + pressatmt_vtot(2,2) &
       &            + pressatmt_vtot(3,3)

!     --- Output pressure data ---

  if (current_step == 1) then

     !- creating pot_cstmnbex subject
     do i = 1, ncstmnbex
        sub_pre_cstmnbex(i) = ' '
        write(sub_tmp,*) 'pot_viricstmnbex_',i
        call excl_sp(sub_tmp,sub_len,sub_pre_cstmnbex(i))
     end do

     !- creating output format
#if defined(_OUTENE_HIGH)
     write(fmt,'(''(A9,1X,''I8''(1X,A23))'')') nprefield
     ! 105  format(A9,1X,<n>(1X,A23))
#else
     write(fmt,'(''(A9,1X,''I8''(1X,A23))'')') nprefield
     ! 105  format(A9,1X,<n>(1X,A23))
#endif

     write(oupre,fmt) 'nstep', &
          &        'pressmol_tot','pressmol_ktot','pressmol_vtot', &
          &        'pressatm_tot','pressatm_ktot','pressatm_vtot', &
          &        'pot_viric_all','pot_virilj_all', &
          &        'pot_virimor_all','pot_virish_all', &
          &        'pot_virirfh_all','pot_viridou_all', &
          &        'pot_viricstmnb', &
          &        (sub_pre_cstmnbex(i),i=1,ncstmnbex), &
          &        'pot_viri_all','viri_fdotd', &
          &        'atm_viribo_all','atm_virian_all', &
          &        'atm_virito_all','atm_viri14_all', &
          &        'atm_viri_const','atm_viri_corr', &
          &        'premolt_tot11','premolt_tot12','premolt_tot13', &
          &        'premolt_tot21','premolt_tot22','premolt_tot23', &
          &        'premolt_tot31','premolt_tot32','premolt_tot33', &
          &        'premolt_ktot11','premolt_ktot12','premolt_ktot13', &
          &        'premolt_ktot21','premolt_ktot22','premolt_ktot23', &
          &        'premolt_ktot31','premolt_ktot32','premolt_ktot33', &
          &        'premolt_vtot11','premolt_vtot12','premolt_vtot13', &
          &        'premolt_vtot21','premolt_vtot22','premolt_vtot23', &
          &        'premolt_vtot31','premolt_vtot32','premolt_vtot33', &
          &        'preatmt_tot11','preatmt_tot12','preatmt_tot13', &
          &        'preatmt_tot21','preatmt_tot22','preatmt_tot23', &
          &        'preatmt_tot31','preatmt_tot32','preatmt_tot33', &
          &        'preatmt_ktot11','preatmt_ktot12','preatmt_ktot13', &
          &        'preatmt_ktot21','preatmt_ktot22','preatmt_ktot23', &
          &        'preatmt_ktot31','preatmt_ktot32','preatmt_ktot33', &
          &        'preatmt_vtot11','preatmt_vtot12','preatmt_vtot13', &
          &        'preatmt_vtot21','preatmt_vtot22','preatmt_vtot23', &
          &        'preatmt_vtot31','preatmt_vtot32','preatmt_vtot33', &
          &        'premol_trPdiv3','premol_trK','premol_trV', &
          &        'preatm_trPdiv3','preatm_trK','preatm_trV'

  end if


  !- creating output format
#if defined(_OUTENE_HIGH)
  write(fmt,'(''(I10,''I8''(1X,E23.16))'')') nprefield
  ! 106  format(I10,<n>(1X,E23.16))
#else
  write(fmt,'(''(I10,''I8''(1X,E23.8))'')') nprefield
  ! 106  format(I10,<n>(1X,E23.8))
#endif

  write(oupre,fmt) current_step, &
       &     pressmol_tot*pref,pressmol_ktot*eref,pressmol_vtot*eref, &
       &     pressatm_tot*pref,pressatm_ktot*eref,pressatm_vtot*eref, &
       &     pot_viric_all*eref,pot_virilj_all*eref, &
       &     pot_virimor_all*eref,pot_virish_all*eref, &
       &     pot_virirfh_all*eref,pot_viridou_all*eref, &
       &     pot_viricstmnb_all*eref, &
       &     (pot_viricstmnbex_all(i)*eref,i=1,ncstmnbex), &
       &     pot_viri_all*eref,viri_fdotd*eref, &
       &     atm_viribo_all*eref,atm_virian_all*eref, &
       &     atm_virito_all*eref,atm_viri14_all*eref, &
       &     atm_viri_const*eref,atm_viri_corr*eref, &
       &     pressmolt_tot(1,1)*pref,pressmolt_tot(1,2)*pref, &
       &     pressmolt_tot(1,3)*pref, &
       &     pressmolt_tot(2,1)*pref,pressmolt_tot(2,2)*pref, &
       &     pressmolt_tot(2,3)*pref, &
       &     pressmolt_tot(3,1)*pref,pressmolt_tot(3,2)*pref, &
       &     pressmolt_tot(3,3)*pref, &
       &     pressmolt_ktot(1,1)*eref,pressmolt_ktot(1,2)*eref, &
       &     pressmolt_ktot(1,3)*eref, &
       &     pressmolt_ktot(2,1)*eref,pressmolt_ktot(2,2)*eref, &
       &     pressmolt_ktot(2,3)*eref, &
       &     pressmolt_ktot(3,1)*eref,pressmolt_ktot(3,2)*eref, &
       &     pressmolt_ktot(3,3)*eref, &
       &     pressmolt_vtot(1,1)*eref,pressmolt_vtot(1,2)*eref, &
       &     pressmolt_vtot(1,3)*eref, &
       &     pressmolt_vtot(2,1)*eref,pressmolt_vtot(2,2)*eref, &
       &     pressmolt_vtot(2,3)*eref, &
       &     pressmolt_vtot(3,1)*eref,pressmolt_vtot(3,2)*eref, &
       &     pressmolt_vtot(3,3)*eref, &
       &     pressatmt_tot(1,1)*pref,pressatmt_tot(1,2)*pref, &
       &     pressatmt_tot(1,3)*pref, &
       &     pressatmt_tot(2,1)*pref,pressatmt_tot(2,2)*pref, &
       &     pressatmt_tot(2,3)*pref, &
       &     pressatmt_tot(3,1)*pref,pressatmt_tot(3,2)*pref, &
       &     pressatmt_tot(3,3)*pref, &
       &     pressatmt_ktot(1,1)*eref,pressatmt_ktot(1,2)*eref, &
       &     pressatmt_ktot(1,3)*eref, &
       &     pressatmt_ktot(2,1)*eref,pressatmt_ktot(2,2)*eref, &
       &     pressatmt_ktot(2,3)*eref, &
       &     pressatmt_ktot(3,1)*eref,pressatmt_ktot(3,2)*eref, &
       &     pressatmt_ktot(3,3)*eref, &
       &     pressatmt_vtot(1,1)*eref,pressatmt_vtot(1,2)*eref, &
       &     pressatmt_vtot(1,3)*eref, &
       &     pressatmt_vtot(2,1)*eref,pressatmt_vtot(2,2)*eref, &
       &     pressatmt_vtot(2,3)*eref, &
       &     pressatmt_vtot(3,1)*eref,pressatmt_vtot(3,2)*eref, &
       &     pressatmt_vtot(3,3)*eref, &
       &     pressmol_trPdiv3*pref,pressmol_trK*eref,pressmol_trV*eref, &
       &     pressatm_trPdiv3*pref,pressatm_trK*eref,pressatm_trV*eref

!     +     +     +     +     +     +     +     +

end subroutine outpre
