!*****************************
!*  rdpara.f90 Ver.3.2       *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <2014-06-13 18:58:14 gota>

subroutine rdpara(iuparavdw,iuparabond, &
     &            xref,eref,mref, &
     &            xcel,ycel,zcel, &
     &            ifljari,ifljgeo)

  use interface_tools

  use md_global

  implicit none

! ARGUMENT:
!     INPUT
  integer,intent(in):: iuparavdw       ! input vdw parameter file unit
  integer,intent(in):: iuparabond      ! input bond parameter file unit

  real(8),intent(in):: xref            ! distanse base value [m]
  real(8),intent(in):: eref            ! energy base value [J]
  real(8),intent(in):: mref            ! mass base value [kg]

  real(8),intent(in):: xcel            ! x cell length [non-d]
  real(8),intent(in):: ycel            ! y cell length [non-d]
  real(8),intent(in):: zcel            ! z cell length [non-d]

  logical,intent(in):: ifljari         ! arithmetic mean for LJ cross parameter
  logical,intent(in):: ifljgeo         ! geometric mean for LJ cross parameter

! LOCAL:
  character(80):: fredat(maxnword)

  real(8):: para_vdw_welij(maxnvdwtyp)
  real(8):: para_vdw_radij(maxnvdwtyp)
  real(8):: tvdw_welij1 = 0.0d0, tvdw_radij1 = 0.0d0
  real(8):: tvdw_welij2 = 0.0d0, tvdw_radij2 = 0.0d0

  integer:: nvdwsptyp
  character(len=5):: para_vdwsptyp(maxnvdwtyp)
  real(8):: para_vdw_welsp(maxnvdwtyp)
  real(8):: para_vdw_radsp(maxnvdwtyp)
  logical:: ifvdwsp_former
  logical:: ifvdwsp_latter

  integer:: nvdwspmlttyp
  character(len=5):: para_vdwspmlttyp(maxnvdwtyp)
  real(8):: para_vdw_welspmlt(maxnvdwtyp)
  real(8):: para_vdw_radspmlt(maxnvdwtyp)

  real(8):: para_mor_wel(maxnmortyp)
  real(8):: para_mor_rad(maxnmortyp)
  real(8):: para_mor_alpha(maxnmortyp)

  real(8):: para_rpvw_rho(maxnrpvwtyp)
  real(8):: para_rpvw_wel(maxnrpvwtyp)
  real(8):: para_rpvw_rad(maxnrpvwtyp)
  real(8):: para_rpvw_alpha(maxnrpvwtyp)

  real(8):: an=6.0221367d23  ! Avogadro's number
  integer:: i,j,k           ! do loop index
  integer:: iflag


  logical:: ifmorse_former
  logical:: ifmorse_latter

  logical:: ifsh_former
  logical:: ifsh_latter

  logical:: ifrfh_former
  logical:: ifrfh_latter

  logical:: ifdou_former
  logical:: ifdou_latter

  logical:: ifrpvw_former
  logical:: ifrpvw_latter

! FUNCTION:

!     rdfree(iu,ndata,fredat)

!     +     +     +     +     +     +     +

!-------- Read parameter of vdw interaction --------

  nvdwtyp = 0
  READVDW:DO
     call rdfree( iuparavdw, maxnword, fredat )

     if (fredat(1)(1:3) == 'END') then
        exit

     else if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) then
        cycle                        ! comment line

     else if (fredat(1) /= ' ') then

        nvdwtyp = nvdwtyp + 1
        vdw_typ(nvdwtyp) = fredat(2)(1:2)
        read(fredat(3),*) para_vdw_welij(nvdwtyp)
        read(fredat(4),*) para_vdw_radij(nvdwtyp)
        vdw_comb(nvdwtyp) = fredat(5)(1:1)

     end if

  END DO READVDW

!---- non-dimensionalize
  do i=1,nvdwtyp
     para_vdw_welij(i) = para_vdw_welij(i) / eref
     para_vdw_radij(i) = para_vdw_radij(i) / xref
  end do

!-------- Read parameter of intramolecular interaction --------

  natmtyp = 0
  nvdwsptyp = 0
  nvdwspmlttyp = 0
  nbondtyp = 0
  nangltyp = 0
  nanglubtyp = 0
  ntorstyp = 0
  ntorsrbtyp = 0
  ntorsimtyp = 0
  nmortyp = 0
  nshtyp = 0
  nrfhtyp = 0
  ndoutyp = 0
  nrpvwtyp = 0

  READBOND:DO
     call rdfree(iuparabond, maxnword, fredat)

     if (fredat(1) == '<END>') then
        exit
     else if (fredat(1) == ' ') then
        cycle READBOND

     else if (fredat(1) == '<ATOM>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           natmtyp = natmtyp + 1
           para_atmtyp(natmtyp) = fredat(2)(1:2)
           read(fredat(3),*) para_atmmass(natmtyp)
           read(fredat(4),*) para_atmchrg(natmtyp)

        end do

     else if (fredat(1) == '<VDW_SP>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nvdwsptyp = nvdwsptyp + 1
           para_vdwsptyp(nvdwsptyp) = fredat(2)(1:5)
           read(fredat(3),*) para_vdw_welsp(nvdwsptyp)
           read(fredat(4),*) para_vdw_radsp(nvdwsptyp)

        end do

     else if (fredat(1) == '<VDW_SP_MLT>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nvdwspmlttyp = nvdwspmlttyp + 1
           para_vdwspmlttyp(nvdwspmlttyp) = fredat(2)(1:5)
           read(fredat(3),*) para_vdw_welspmlt(nvdwspmlttyp)
           read(fredat(4),*) para_vdw_radspmlt(nvdwspmlttyp)

        end do

     else if (fredat(1) == '<BOND>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nbondtyp = nbondtyp + 1
           para_bondtyp(nbondtyp) = fredat(2)(1:5)
           read(fredat(3),*) para_cbond(nbondtyp)
           read(fredat(4),*) para_eqbond(nbondtyp)

        end do

     else if (fredat(1) == '<ANGLE>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nangltyp = nangltyp + 1
           para_angltyp(nangltyp) = fredat(2)(1:8)
           read(fredat(3),*) para_cangl(nangltyp)
           read(fredat(4),*) para_eqangl(nangltyp)

        end do

     else if (fredat(1) == '<ANGLE_UB>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nanglubtyp = nanglubtyp + 1
           para_anglubtyp(nanglubtyp) = fredat(2)(1:8)
           read(fredat(3),*) para_canglub(nanglubtyp)
           read(fredat(4),*) para_eqanglub(nanglubtyp)
           read(fredat(5),*) para_cbondub(nanglubtyp)
           read(fredat(6),*) para_eqbondub(nanglubtyp)

        end do

#if !defined(_CHARMM_NONB14)
     else if (fredat(1) == '<TORSION>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           ntorstyp = ntorstyp + 1
           para_torstyp(ntorstyp) = fredat(2)(1:11)
           read(fredat(3),*) para_divfac(ntorstyp)
           read(fredat(4),*) para_barhig(ntorstyp)
           read(fredat(5),*) para_phsang(ntorstyp)
           read(fredat(6),*) para_period(ntorstyp)
           read(fredat(7),*) para_divfac_vdw(ntorstyp)
           read(fredat(8),*) para_divfac_elc(ntorstyp)

        end do

#else
     else if (fredat(1) == '<TORSION_CH>') then !!! by KAWAGOE
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           ntorstyp = ntorstyp + 1
           para_torstyp(ntorstyp) = fredat(2)(1:11)
           read(fredat(3),*) para_divfac(ntorstyp)
           read(fredat(4),*) para_barhig(ntorstyp)
           read(fredat(5),*) para_phsang(ntorstyp)
           read(fredat(6),*) para_period(ntorstyp)
           read(fredat(7),*) para_divfac_vdw(ntorstyp)
           read(fredat(8),*) para_divfac_elc(ntorstyp)
           read(fredat(9),*) para_wfac(ntorstyp)

        end do
#endif

     else if (fredat(1) == '<TORSION_RB>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           ntorsrbtyp = ntorsrbtyp + 1
           para_torsrbtyp(ntorsrbtyp) = fredat(2)(1:11)
           read(fredat(3),*) para_barhigrb(1,ntorsrbtyp)
           read(fredat(4),*) para_barhigrb(2,ntorsrbtyp)
           read(fredat(5),*) para_barhigrb(3,ntorsrbtyp)
           read(fredat(6),*) para_barhigrb(4,ntorsrbtyp)
           read(fredat(7),*) para_barhigrb(5,ntorsrbtyp)
           read(fredat(8),*) para_barhigrb(6,ntorsrbtyp)
           read(fredat(9),*) para_divfac_vdwrb(ntorsrbtyp)
           read(fredat(10),*) para_divfac_elcrb(ntorsrbtyp)

        end do

     else if (fredat(1) == '<TORSION_IM>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           ntorsimtyp = ntorsimtyp + 1
           para_torsimtyp(ntorsimtyp) = fredat(2)(1:11)
           read(fredat(3),*) para_barhigim(ntorsimtyp)
           read(fredat(4),*) para_phsangim(ntorsimtyp)

        end do

     else if (fredat(1) == '<MORSE>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nmortyp = nmortyp + 1
           para_mortyp(nmortyp) = fredat(2)(1:5)
           read(fredat(3),*) para_mor_wel(nmortyp)
           read(fredat(4),*) para_mor_rad(nmortyp)
           read(fredat(5),*) para_mor_alpha(nmortyp)

        end do

     else if (fredat(1) == '<SH>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nshtyp = nshtyp + 1
           if (fredat(2)(1:1) == 'O') then
              para_shtyp(1) = fredat(2)(1:5)
              read(fredat(3),*) para_sh_a(1)
              read(fredat(4),*) para_sh_a(2)
              read(fredat(5),*) para_sh_a(3)
              read(fredat(6),*) para_sh_b(1)
              read(fredat(7),*) para_sh_b(2)
              read(fredat(8),*) para_sh_b(3)
              read(fredat(9),*) para_sh_c

           else if (fredat(2)(1:1) == 'H') then
              para_shtyp(2) = fredat(2)(1:5)
              read(fredat(3),*) para_sh_a(4)
              read(fredat(4),*) para_sh_b(4)

           else
              write(6,*) 'Error: something is wrong with SH entry'
              stop
           end if

        end do

     else if (fredat(1) == '<RFH>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nrfhtyp = nrfhtyp + 1
           para_rfhtyp(nrfhtyp) = fredat(2)(1:5)
           if (fredat(2)(1:1) == 'F') then   ! Fe-O interaction
              read(fredat(3),*) para_rfhfo_a
              read(fredat(4),*) para_rfhfo_b
              read(fredat(5),*) para_rfhfo_c
              read(fredat(6),*) para_rfhfo_d

           else if ((fredat(2)(1:1) == 'O') .and. &
                &   (fredat(2)(4:4) == 'O')) then   ! O-O interaction
              read(fredat(3),*) para_rfhoo_a
              read(fredat(4),*) para_rfhoo_b
              read(fredat(5),*) para_rfhoo_c
              read(fredat(6),*) para_rfhoo_d
              read(fredat(7),*) para_rfhoo_lja
              read(fredat(8),*) para_rfhoo_ljb

           else if ((fredat(2)(1:1) == 'O') .and. &
                &   (fredat(2)(4:4) == 'H')) then   ! O-H interaction
              read(fredat(3),*) para_rfhoh_a
              read(fredat(4),*) para_rfhoh_b
              read(fredat(5),*) para_rfhoh_c
              read(fredat(6),*) para_rfhoh_d
              read(fredat(7),*) para_rfhoh_e
              read(fredat(8),*) para_rfhoh_req

           else
              write(6,*) 'Error: something is wrong with RFH entry'
              stop
           end if

        end do

     else if (fredat(1) == '<DOU>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           ndoutyp = ndoutyp + 1

           para_doutyp(ndoutyp) = fredat(2)(1:5)

           if (fredat(2)(1:1) == 'O') then
              read(fredat(3),*) para_dou_wel_o
              read(fredat(4),*) para_dou_rad_o
              read(fredat(5),*) para_dou_beta_o

           else if (fredat(2)(1:1) == 'H') then
              read(fredat(3),*) para_dou_wel_h
              read(fredat(4),*) para_dou_rad_h
              read(fredat(5),*) para_dou_beta_h
              read(fredat(6),*) para_dou_gamma

           else
              write(6,*) 'Error: something is wrong with DOU entry'
              stop
           end if

        end do

     else if (fredat(1) == '<POSRES>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           read(fredat(2),*) para_posres

        end do

     else if (fredat(1) == '<CNP_VW>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           read(fredat(2),*) alpha_wall

        end do

     else if (fredat(1) == '<RP_VW>') then
        do
           call rdfree(iuparabond,maxnword,fredat)
           if (fredat(1) == ' ') cycle READBOND
           if ((fredat(1)(1:1) == '#') .or. (fredat(1)(1:1) == ';')) cycle
                                                            ! comment line

           nrpvwtyp = nrpvwtyp + 1
           para_rpvwtyp(nrpvwtyp) = fredat(2)(1:5)
           read(fredat(3),*) para_rpvw_rho(nrpvwtyp)
           read(fredat(4),*) para_rpvw_wel(nrpvwtyp)
           read(fredat(5),*) para_rpvw_rad(nrpvwtyp)
           read(fredat(6),*) para_rpvw_alpha(nrpvwtyp)

        end do

     end if

  END DO READBOND

!---- non-dimensionalize
  para_atmmass(1:natmtyp) = para_atmmass(1:natmtyp) * 1.0d-3/an / mref

  para_vdw_welsp(1:nvdwsptyp) = para_vdw_welsp(1:nvdwsptyp) / eref
  para_vdw_radsp(1:nvdwsptyp) = para_vdw_radsp(1:nvdwsptyp) / xref

  para_cbond(1:nbondtyp) = para_cbond(1:nbondtyp) / (eref/xref**2)
  para_eqbond(1:nbondtyp) = para_eqbond(1:nbondtyp) / xref

  para_cangl(1:nangltyp) = para_cangl(1:nangltyp) / eref

  para_canglub(1:nanglubtyp) = para_canglub(1:nanglubtyp) / eref
  para_cbondub(1:nanglubtyp) = para_cbondub(1:nanglubtyp) / (eref/xref**2)
  para_eqbondub(1:nanglubtyp) = para_eqbondub(1:nanglubtyp) / xref

  para_barhig(1:ntorstyp) = para_barhig(1:ntorstyp) / eref
#if defined(_CHARMM_NONB14)
! For CHARMM F.F., para_divfac_vdw is regarded as 1-4 LJ epsilon and
!                  para_divfac_elc is regarded as 1-4 LJ radius
  para_divfac_vdw(1:ntorstyp) = para_divfac_vdw(1:ntorstyp) / eref
  para_divfac_elc(1:ntorstyp) = para_divfac_elc(1:ntorstyp) / xref
#endif

  para_barhigrb(1:6,1:ntorsrbtyp) = para_barhigrb(1:6,1:ntorsrbtyp) / eref

  para_barhigim(1:ntorsimtyp) = para_barhigim(1:ntorsimtyp) / eref

  para_mor_wel(1:nmortyp) = para_mor_wel(1:nmortyp) / eref
  para_mor_rad(1:nmortyp) = para_mor_rad(1:nmortyp) / xref
  para_mor_alpha(1:nmortyp) = para_mor_alpha(1:nmortyp) * xref

  para_sh_a(1:4) = para_sh_a(1:4) / eref
  para_sh_b(1:4) = para_sh_b(1:4) / xref
  para_sh_c    = para_sh_c    / xref

  para_rfhfo_a = para_rfhfo_a / eref
  para_rfhfo_b = para_rfhfo_b * xref
  para_rfhfo_c = para_rfhfo_c / eref ! [J \AA^6]
  para_rfhfo_d = para_rfhfo_d / eref ! [J \AA^12]
  para_rfhoo_a = para_rfhoo_a / eref
  para_rfhoo_b = para_rfhoo_b * xref
!      para_rfhoo_c = para_rfhoo_c
  para_rfhoo_d = para_rfhoo_d / xref
  para_rfhoo_lja = para_rfhoo_lja / eref ! [J \AA^12]
  para_rfhoo_ljb = para_rfhoo_ljb / eref ! [J \AA^6]
  para_rfhoh_a = para_rfhoh_a / eref / xref
  para_rfhoh_b = para_rfhoh_b * xref
  para_rfhoh_c = para_rfhoh_c / eref * xref**2
  para_rfhoh_d = para_rfhoh_d / eref * xref
  para_rfhoh_e = para_rfhoh_e * xref**2
  para_rfhoh_req = para_rfhoh_req / xref

  para_dou_wel_o = para_dou_wel_o / eref
  para_dou_wel_h = para_dou_wel_h / eref
  para_dou_rad_o = para_dou_rad_o / xref
  para_dou_rad_h = para_dou_rad_h / xref
  para_dou_beta_o = para_dou_beta_o * xref
  para_dou_beta_h = para_dou_beta_h * xref

  alpha_wall = (alpha_wall * xref**3 / eref) * xcel * ycel
                                          ! convert pressure to force

  para_rpvw_rho(1:nrpvwtyp)   = para_rpvw_rho * (xref**3)
  para_rpvw_wel(1:nrpvwtyp)   = para_rpvw_wel / eref
  para_rpvw_rad(1:nrpvwtyp)   = para_rpvw_rad / xref
  para_rpvw_alpha(1:nrpvwtyp) = para_rpvw_alpha * xref

!---- order and link vdw parameter
  do i = 1, natmtyp
     do j = 1, natmtyp
        iflag = 0
        do k = 1, nvdwtyp

           if (para_atmtyp(i)(1:2) == vdw_typ(k)(1:2)) then
              tvdw_welij1 = para_vdw_welij(k)
              tvdw_radij1 = para_vdw_radij(k)
              iflag = iflag + 1
           end if
           if (para_atmtyp(j)(1:2) == vdw_typ(k)(1:2)) then
              tvdw_welij2 = para_vdw_welij(k)
              tvdw_radij2 = para_vdw_radij(k)
              iflag = iflag + 1
           end if

        end do
        if (iflag /= 2) then
           write(6,*) 'Warning: Failure in vdw parameter'
           stop
        end if

        vdw_welij(i,j) = sqrt(tvdw_welij1*tvdw_welij2)

        if ((vdw_comb(i) == 'A') .and. (vdw_comb(j) == 'A')) then

           vdw_radij(i,j) = 0.5d0 * (tvdw_radij1 + tvdw_radij2)

        else if ((vdw_comb(i) == 'G') .and. (vdw_comb(j) == 'G')) then

           vdw_radij(i,j) = sqrt(tvdw_radij1*tvdw_radij2)

        else

           if (ifljari) then
              vdw_radij(i,j) = 0.5d0 * (tvdw_radij1 + tvdw_radij2)
           else if (ifljgeo) then
              vdw_radij(i,j) = sqrt(tvdw_radij1*tvdw_radij2)
           end if

        end if

     end do
  end do

!---- specific vdW parameter
  do i = 1, natmtyp
     do j = 1, natmtyp

        do k = 1, nvdwsptyp
           ifvdwsp_former = .false.
           ifvdwsp_latter = .false.

!----      scan i
           if (para_atmtyp(i)(1:2) == para_vdwsptyp(k)(1:2)) then
              ifvdwsp_former = .true.
           else if (para_atmtyp(i)(1:2) == para_vdwsptyp(k)(4:5)) then
              ifvdwsp_latter = .true.
           end if

!----      scan j
           if (ifvdwsp_former) then

              if (para_atmtyp(j)(1:2) == para_vdwsptyp(k)(4:5)) then
                 vdw_welij(i,j) = para_vdw_welsp(k)
                 vdw_radij(i,j) = para_vdw_radsp(k)
                 exit

              end if

           else if (ifvdwsp_latter) then

              if (para_atmtyp(j)(1:2) == para_vdwsptyp(k)(1:2)) then

                 vdw_welij(i,j) = para_vdw_welsp(k)
                 vdw_radij(i,j) = para_vdw_radsp(k)
                 exit

              end if

           end if

        end do

!----   vdW sp multiplier (VDW_SP_MLT)
        do k = 1, nvdwspmlttyp
           ifvdwsp_former = .false.
           ifvdwsp_latter = .false.

!----      scan i
           if (para_atmtyp(i)(1:2) == para_vdwspmlttyp(k)(1:2)) then
              ifvdwsp_former = .true.
           else if (para_atmtyp(i)(1:2) == para_vdwspmlttyp(k)(4:5)) then
              ifvdwsp_latter = .true.
           end if

!----      scan j
           if (ifvdwsp_former) then

              if (para_atmtyp(j)(1:2) == para_vdwspmlttyp(k)(4:5)) then

                 vdw_welij(i,j) = vdw_welij(i,j) * para_vdw_welspmlt(k)
                 vdw_radij(i,j) = vdw_radij(i,j) * para_vdw_radspmlt(k)
                 exit

              end if

           else if (ifvdwsp_latter) then

              if (para_atmtyp(j)(1:2) == para_vdwspmlttyp(k)(1:2)) then

                 vdw_welij(i,j) = vdw_welij(i,j) * para_vdw_welspmlt(k)
                 vdw_radij(i,j) = vdw_radij(i,j) * para_vdw_radspmlt(k)
                 exit

              end if

           end if

        end do

     end do
  end do

!---- order and link Morse parameter
  do i = 1, natmtyp
     do j = 1, natmtyp

        do k = 1, nmortyp
           ifmorse_former = .false.
           ifmorse_latter = .false.

!----          scan i
           if (para_atmtyp(i)(1:2) == para_mortyp(k)(1:2)) then
              ifmorse_former = .true.
           else if (para_atmtyp(i)(1:2) == para_mortyp(k)(4:5)) then
              ifmorse_latter = .true.
           end if

!----          scan j
           if (ifmorse_former) then

              if (para_atmtyp(j)(1:2) == para_mortyp(k)(4:5)) then
                 para_welmor(i,j) = para_mor_wel(k)
                 para_radmor(i,j) = para_mor_rad(k)
                 para_alphamor(i,j) = para_mor_alpha(k)
                 exit

              end if

           else if (ifmorse_latter) then

              if (para_atmtyp(j)(1:2) == para_mortyp(k)(1:2)) then

                 para_welmor(i,j) = para_mor_wel(k)
                 para_radmor(i,j) = para_mor_rad(k)
                 para_alphamor(i,j) = para_mor_alpha(k)
                 exit

              end if

           end if

           para_welmor(i,j) = 0.0d0
           para_radmor(i,j) = 0.0d0
           para_alphamor(i,j) = 0.0d0

        end do

     end do
  end do

!---- order and link RealParticle-VirtualWall parameter
  do i = 1, natmtyp
     do j = 1, natmtyp

        do k = 1, nrpvwtyp
           ifrpvw_former = .false.
           ifrpvw_latter = .false.

!----          scan i
           if (para_atmtyp(i)(1:2) == para_rpvwtyp(k)(1:2)) then
              ifrpvw_former = .true.
           else if (para_atmtyp(i)(1:2) == para_rpvwtyp(k)(4:5)) then
              ifrpvw_latter = .true.
           end if

!----          scan j
           if (ifrpvw_former) then

              if (para_atmtyp(j)(1:2) == para_rpvwtyp(k)(4:5)) then

                 para_rhorpvw(i,j) = para_rpvw_rho(k)
                 para_welrpvw(i,j) = para_rpvw_wel(k)
                 para_radrpvw(i,j) = para_rpvw_rad(k)
                 para_alpharpvw(i,j) = para_rpvw_alpha(k)

                 exit

              end if

           else if (ifrpvw_latter) then

              if (para_atmtyp(j)(1:2) == para_rpvwtyp(k)(1:2)) then

                 para_rhorpvw(i,j) = para_rpvw_rho(k)
                 para_welrpvw(i,j) = para_rpvw_wel(k)
                 para_radrpvw(i,j) = para_rpvw_rad(k)
                 para_alpharpvw(i,j) = para_rpvw_alpha(k)

                 exit

              end if

           end if

           para_rhorpvw(i,j) = 0.0d0
           para_welrpvw(i,j) = 0.0d0
           para_radrpvw(i,j) = 0.0d0
           para_alpharpvw(i,j) = 1.0d0

        end do

     end do
  end do

!---- register intermolecular interaction type for particle pair

!     - clear all parameter to vdw type
  do i = 1, natmtyp
     do j = i, natmtyp
        inter_inttyp(i,j) = INTTYPE_VDW
        inter_inttyp(j,i) = inter_inttyp(i,j)
     end do
  end do

!     - detect morse pair
  do i = 1, natmtyp
     do j = i, natmtyp

        do k = 1, nmortyp
           ifmorse_former = .false.
           ifmorse_latter = .false.

!----          scan i
           if (para_atmtyp(i)(1:2) == para_mortyp(k)(1:2)) then
              ifmorse_former = .true.
           else if (para_atmtyp(i)(1:2) == para_mortyp(k)(4:5)) then
              ifmorse_latter = .true.
           else
              cycle
           end if

!----          scan j
           if (ifmorse_former) then

              if (para_atmtyp(j)(1:2) == para_mortyp(k)(4:5)) then
                 inter_inttyp(i,j) = INTTYPE_MOR
                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           else if (ifmorse_latter) then

              if (para_atmtyp(j)(1:2) == para_mortyp(k)(1:2)) then

                 inter_inttyp(i,j) = INTTYPE_MOR
                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           end if

        end do

     end do
  end do

!     - detect SH pair
  do i = 1, natmtyp
     do j = i, natmtyp

        do k = 1, nshtyp
           ifsh_former = .false.
           ifsh_latter = .false.

!----          scan i
           if (para_atmtyp(i)(1:2) == para_shtyp(k)(1:2)) then
              ifsh_former = .true.
           else if (para_atmtyp(i)(1:2) == para_shtyp(k)(4:5)) then
              ifsh_latter = .true.
           else
              cycle
           end if

!----          scan j
           if (ifsh_former) then

              if (para_atmtyp(j)(1:2) == para_shtyp(k)(4:5)) then
                 inter_inttyp(i,j) = INTTYPE_SH
                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           else if (ifsh_latter) then

              if (para_atmtyp(j)(1:2) == para_shtyp(k)(1:2)) then

                 inter_inttyp(i,j) = INTTYPE_SH
                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           end if

        end do

     end do
  end do

!     - detect RFH pair
  do i = 1, natmtyp
     do j = i, natmtyp

        do k = 1, nrfhtyp
           ifrfh_former = .false.
           ifrfh_latter = .false.

!----          scan i
           if (para_atmtyp(i)(1:2) == para_rfhtyp(k)(1:2)) then
              ifrfh_former = .true.
           else if (para_atmtyp(i)(1:2) == para_rfhtyp(k)(4:5)) then
              ifrfh_latter = .true.
           else
              cycle
           end if

!----          scan j
           if (ifrfh_former) then

              if (para_atmtyp(j)(1:2) == para_rfhtyp(k)(4:5)) then

                 if (para_atmtyp(i)(1:1) == 'F') then ! Fe-O
                    inter_inttyp(i,j) = INTTYPE_RFHFO

                 else if ((para_atmtyp(i)(1:1) == 'O') .and. &
                      &   (para_atmtyp(j)(1:1) == 'O')) then ! O-O
                    inter_inttyp(i,j) = INTTYPE_RFHOO

                 else if ((para_atmtyp(i)(1:1) == 'O') .and. &
                      &   (para_atmtyp(j)(1:1) == 'H')) then ! O-H
                    inter_inttyp(i,j) = INTTYPE_RFHOH

                 end if

                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           else if (ifrfh_latter) then

              if (para_atmtyp(j)(1:2) == para_rfhtyp(k)(1:2)) then

                 if (para_atmtyp(j)(1:1) == 'F') then ! Fe-O
                    inter_inttyp(i,j) = INTTYPE_RFHFO

                 else if ((para_atmtyp(i)(1:1) == 'O') .and. &
                      &   (para_atmtyp(j)(1:1) == 'O')) then ! O-O
                    inter_inttyp(i,j) = INTTYPE_RFHOO

                 else if ((para_atmtyp(j)(1:1) == 'O') .and. &
                      &   (para_atmtyp(i)(1:1) == 'H')) then ! O-H
                    inter_inttyp(i,j) = INTTYPE_RFHOH

                 end if

                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           end if

        end do

     end do
  end do

!     - detect DOU pair
  do i = 1, natmtyp
     do j = i, natmtyp

        do k = 1, ndoutyp
           ifdou_former = .false.
           ifdou_latter = .false.

!----          scan i
           if (para_atmtyp(i)(1:2) == para_doutyp(k)(1:2)) then
              ifdou_former = .true.
           else if (para_atmtyp(i)(1:2) == para_doutyp(k)(4:5)) then
              ifdou_latter = .true.
           else
              cycle
           end if

!----          scan j
           if (ifdou_former) then

              if (para_atmtyp(j)(1:2) == para_doutyp(k)(4:5)) then

                 if (para_atmtyp(i)(1:1) == 'O') then ! O atoms
                    inter_inttyp(i,j) = INTTYPE_DOUO
                 else if (para_atmtyp(i)(1:1) == 'H') then ! H atoms
                    inter_inttyp(i,j) = INTTYPE_DOUH
                 end if

                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           else if (ifdou_latter) then

              if (para_atmtyp(j)(1:2) == para_doutyp(k)(1:2)) then

                 if (para_atmtyp(j)(1:1) == 'O') then ! O atoms
                    inter_inttyp(i,j) = INTTYPE_DOUO
                 else if (para_atmtyp(j)(1:1) == 'H') then ! H atoms
                    inter_inttyp(i,j) = INTTYPE_DOUH
                 end if

                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           end if

        end do

     end do
  end do

!     - detect RealParticle-VirtualWall pair
  do i = 1, natmtyp
     do j = i, natmtyp

        do k = 1, nrpvwtyp
           ifrpvw_former = .false.
           ifrpvw_latter = .false.

!----          scan i
           if (para_atmtyp(i)(1:2) == para_rpvwtyp(k)(1:2)) then
              ifrpvw_former = .true.
           else if (para_atmtyp(i)(1:2) == para_rpvwtyp(k)(4:5)) then
              ifrpvw_latter = .true.
           else
              cycle
           end if

!----          scan j
           if (ifrpvw_former) then

              if (para_atmtyp(j)(1:2) == para_rpvwtyp(k)(4:5)) then
                 inter_inttyp(i,j) = INTTYPE_RPVW
                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           else if (ifrpvw_latter) then

              if (para_atmtyp(j)(1:2) == para_rpvwtyp(k)(1:2)) then
                 inter_inttyp(i,j) = INTTYPE_RPVW
                 inter_inttyp(j,i) = inter_inttyp(i,j)
                 vdw_welij(i,j) = 0.0d0
                 vdw_welij(j,i) = 0.0d0
                 vdw_radij(i,j) = 0.0d0
                 vdw_radij(j,i) = 0.0d0
                 exit

              end if

           end if

        end do

     end do
  end do

#if defined(_INTTYP_DEBUG)
  do i = 1, natmtyp
     do j = i, natmtyp
        write(6,*) i,j,inter_inttyp(i,j)
     end do
  end do
#endif

!---- dividing factor of 1-4 nonb potential for torsion (periodic type, CHARMM)

  do i = 1, ntorstyp
     if ((abs(para_divfac_vdw(i)) < 1.0d-8) .and. &
          & (abs(para_divfac_elc(i)) < 1.0d-8)) then
        ifcalnonb14(i) = .false.
     else
        ifcalnonb14(i) = .true.
     end if
  end do


!---- dividing factor of 1-4 nonb potential for torsion RB (for OPLS)

  do i = 1, ntorsrbtyp
     if ((abs(para_divfac_vdwrb(i)) < 1.0d-8) .and. &
          & (abs(para_divfac_elcrb(i)) < 1.0d-8)) then
        ifcalnonb14_rb(i) = .false.
     else
        ifcalnonb14_rb(i) = .true.
     end if
  end do

  close(iuparavdw)
  close(iuparabond)

!     +     +     +     +     +     +     +

end subroutine rdpara
