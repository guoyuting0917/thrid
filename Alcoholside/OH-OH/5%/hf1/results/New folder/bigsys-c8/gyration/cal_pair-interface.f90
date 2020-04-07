!#####################################################################
!     Calculation of configuration property
!     using output files via LAMMPS
!     (dcd file and data file)
!     * Radius of gyration
!     * Orientation order parameter
!     * Charge distribution
! 
!     Version 3.4
!     Final update 2019/12/13 (Y. Kawagoe)
!#####################################################################
!    <PROBLEMS>
!    Treatment of min in x, y, z dire
!    dcd file do not have xmin, ymin, zmin info (just box size),
!    so slab calculation is not corrrect when box size is changing and 
!    have atoms in the vicinity of z boundary.
!    In unchanged system or system having vaccum regs in the both ends, slab is available.
!    PBC treatment is correct because it use box size not min/max value.
!
!    Temporary, min value was set using data file
!
!    CV in x-dire is applied only for Rg, not for charge,
!    becasue charge out values are sumation, not average.
!
!#####################################################################
program calc_config
  implicit none

  character(30):: datafile
  character(30):: dcdfile

! I/O unit
  integer,parameter:: iu_data=11
  integer,parameter:: iu_dcd=12
  integer,parameter:: iu_in=13
  integer,parameter:: iu_tmp=14
  
  integer,parameter:: iu_out=21
  
  integer:: ios
  
  integer,parameter:: max_ndata=20 !max ndata for fredat
  character(80):: fredat(max_ndata)
  integer:: end_flag

  integer:: n_atom
  integer:: n_atom_type
  integer:: n_bond
  integer:: n_mol
  integer,dimension(:),allocatable:: atom_type  ! atom type per atom
  integer,dimension(:),allocatable:: mol_id  ! molecule id per atom
  real*8,dimension(:),allocatable:: charge  ! charge per atom
  integer,dimension(:,:),allocatable:: cnct  !connection info per atom
  real*8,dimension(:),allocatable:: mass_type  ! mass type per atomtype
  integer,dimension(:,:),allocatable:: bond

  character(4):: hdr 
  integer:: nset,istart,iend,nsavc,icntrl(5),namnf,jcntrl(9),natoms 
  integer:: ntitle
  real*8:: delta, xtable(6), dummy
  real,dimension(:,:,:),allocatable :: pos ! (xyz,atom_id,step)
  character(80):: title(10) 

  real*8:: x_min,x_max
  real*8:: y_min,y_max
  real*8:: z_min,z_max
  real*8:: lx,ly,lz
  real*8:: x_half,y_half,z_half

  ! for slab
  real*8:: zmin_slab,zmax_slab
  real*8:: lz_slab
  real*8:: pos_z

  real*8,dimension(:,:),allocatable:: mpos ! mass center of molecule
  real*8,dimension(:),allocatable:: mass_mol ! molecular mass
  real*8:: pos0(3)
  integer,dimension(:,:),allocatable:: mol_set ! list of atoms in molecule
  integer,dimension(:),allocatable:: min_mol_set ! minval of mol_set in each molecule
  integer,dimension(:),allocatable:: mol_count
  integer,dimension(:),allocatable:: atom_flag

  integer:: s_step,e_step
  real*8:: delta_t
  integer:: i,j,k,it,it0,i2,j2,k2,l,i3
  integer:: ityp,imol
  integer:: dummy_i
  integer:: count,count2

  real*8::loop_length
  integer::os_count,h2_count,m,n,o1,h1,o11,h11,loop,loop_start,loop_length1,index,in_flag
  integer::o_number,h_number,os_type,hs_type,in_maxmol
  integer,dimension(:),allocatable::os,h2
  integer,dimension(:,:),allocatable::hydron_symbol
  integer,dimension(:,:),allocatable::hydron_symbol_start
  integer,dimension(:,:),allocatable::hydron_acf
  integer,dimension(:),allocatable::hydron_acf_sum
  integer,dimension(:),allocatable::hydron_acf_list
  real*8::pos_o1(3),pos_h1(3),pos_o11(3),pos_h11(3)
  real*8::dis_ao_ah,dis_so_ah,dis_so_sh,dis_ao_sh,dis_oh11
  real*8::vec_ao_ah(3),vec_so_ah(3),vec_so_sh(3),vec_ao_sh(3),vec_oh11(3),cos_oh1,cos_oh2

  integer,dimension(:),allocatable::count_hyd_bond

  real*8:: dx,dy,dz

! Radius of gyration
  real*8,dimension(:),allocatable:: Rg
  real*8,dimension(:,:),allocatable:: Rg_xyz
  real*8,dimension(:),allocatable:: ave_Rg,in_Rg
  real*8,dimension(:,:),allocatable:: ave_Rg_xyz,in_Rg_xyz
  integer,dimension(:),allocatable:: count_Rg,in_count_Rg

  integer:: n_f
  real*8,dimension(:,:),allocatable:: f_Rg,in_f_Rg
  real*8,dimension(:,:,:),allocatable:: f_Rg_xyz,in_f_Rg_xyz
  real*8:: df_Rg
  real*8:: max_Rg

  integer:: n_moltype
  integer,dimension(:),allocatable:: st_molid
  integer,dimension(:),allocatable:: ed_molid
  character(10),dimension(:),allocatable:: molname
  integer:: nslab_Rg
  real*8:: divz_Rg
  real*8,dimension(:,:,:),allocatable:: slab_Rg !(xyz,slabs,moltype)
  real*8,dimension(:,:,:),allocatable:: slab_Rg2 !(xyz,slabs,moltype) for st err
  integer,dimension(:,:),allocatable:: count_smp_slab !(slabs,moltype)

! Orientation oder parameter
  real*8,parameter:: pi=4.0d0*atan(1.0d0)
  integer,parameter:: n_theta=36
  real*8,parameter:: dtheta=pi/dble(n_theta)*0.50d0
  integer,parameter:: n_cosine=100
  real*8,parameter:: dcosine=1.0d0/dble(n_cosine)

  integer:: oopflag
  character(50):: dummy_a
  character(100),dimension(:),allocatable:: vecfile
  integer,dimension(:,:),allocatable:: veclist_i,veclist_j
  integer:: max_nlist
  integer,dimension(:),allocatable:: n_veclist
  integer,dimension(:),allocatable::f_cosine1_list,f_cosine2_list
  
  real*8:: theta
  real*8:: cosine,cosine_oh,cosine_oh_sum,f_cosine_oh_sum
  real,dimension(:),allocatable::f_cosine1_sum,cosine1_sum,f_cosine2_sum,cosine2_sum
  real*8:: oop_vec(4)

  real*8,dimension(:,:),allocatable:: f_theta
  real*8,dimension(:,:,:),allocatable:: f_theta1
  real*8,dimension(:,:),allocatable:: f_cosine
  real*8,dimension(:,:,:),allocatable:: f_cosine_decom
  real*8,dimension(:),allocatable::f_cosine_oh
  real*8,dimension(:),allocatable:: oop
  integer,dimension(:),allocatable:: count_oop
  integer:: count_oh
  real*8,dimension(:,:),allocatable:: slab_oop !(slabs,moltype)
  real*8,dimension(:,:),allocatable:: slab_oop2 !(slabs,moltype)
  real*8,dimension(:,:,:),allocatable:: slab_tensor !(9,slabs,moltype) for tensor order parameter
  integer:: nslab_oop
  integer,dimension(:,:),allocatable:: count_slab_oop !(slabs,moltype)
  real*8:: divz_oop

  integer,parameter:: n_meas=2
  real*8:: ee_vec(4)
  real*8:: mchain_vec(4)



  real*8:: f_ee1(n_theta)
  real*8:: f_ee2(n_theta)
  real*8:: S_ee1
  real*8:: S_ee2

  real*8:: f_bv1(n_theta)
  real*8:: f_bv2(n_theta)
  real*8:: S_bv1
  real*8:: S_bv2

  real*8:: f_mc1(n_theta)
  real*8:: f_mc2(n_theta)
  real*8:: S_mc1
  real*8:: S_mc2

  real*8:: f_mcskip1(n_theta)
  real*8:: f_mcskip2(n_theta)
  real*8:: S_mcskip1
  real*8:: S_mcskip2

  real*8:: z_meas1(2)
  real*8:: z_meas2(2)
  real*8:: x_meas1(2)
  real*8:: x_meas2(2)
  integer:: smp_flag
  !temp
  integer:: count_ee1,count_ee2
  integer:: count_bv1,count_bv2
  integer:: count_mc1,count_mc2
  integer:: count_mcbond
  integer:: count_mcskip1,count_mcskip2

  integer:: sid,eid,nid,idir
  integer:: s_atom,e_atom

  integer:: iHdel
  integer:: icheck

  integer:: n_mc
  integer:: n_mc_bond
  integer:: n_mcskip_bond
  integer,dimension(:,:),allocatable:: mc_ij
  integer,dimension(:),allocatable:: mc_count_atom
  integer,dimension(:,:),allocatable:: mc_bond
  integer,dimension(:,:),allocatable:: mcskip_bond
  integer,dimension(:),allocatable:: end_atom
  integer:: mc_flag

  ! for check
  integer,parameter:: n_check=5000
  real*8:: check_vec(4)
  real*8:: cos_theta
  real*8:: sin_theta
  real*8:: phi
  real*8,parameter:: r_check=2.0d0/pi
  real*8:: rnd
  real*8:: f_check(n_theta)
  real*8:: S_check
  integer:: count_check

! Charge distribution
  integer:: nslab_ch
  integer:: count_smp
  real*8:: divz_ch
  real*8,dimension(:,:),allocatable:: slab_charge
  real*8,dimension(:,:),allocatable:: slab_charge2
  real*8,dimension(:,:),allocatable:: slab_charge2_tmp

! molecule id sort
  integer:: n_molmax
  integer:: flag
  integer,dimension(:),allocatable:: read2seri
  integer,dimension(:),allocatable:: seri2read


  oopflag=0
  in_flag=0

!### READ INPUT FILE ###
  print*,'### Read input file: calc_config.ini ###'
  open(iu_in,file='calc_config.ini',status='old',iostat=ios)
  if(ios.ne.0)then
     print*,'ERROR: no ini file'
     stop
  endif
  do
     call rdfree(iu_in,max_ndata,fredat,end_flag)

     if((fredat(1)(1:1).eq.'#').or.(fredat(1)(1:1).eq.';')) cycle !comment line

     if(fredat(1)(1:9).eq.'data_file')read(fredat(2),'(a)')datafile
     if(fredat(1)(1:8).eq.'dcd_file')read(fredat(2),'(a)')dcdfile
     if(fredat(1)(1:9).eq.'start_set')read(fredat(2),*)s_step
     if(fredat(1)(1:7).eq.'end_set')read(fredat(2),*)e_step     
     if(fredat(1)(1:7).eq.'delta_t')read(fredat(2),*)delta_t

     !Rg
     if(fredat(1)(1:6).eq.'max_Rg')read(fredat(2),*)max_Rg
     if(fredat(1)(1:5).eq.'n_fRg')read(fredat(2),*)n_f
     if(fredat(1)(1:9).eq.'n_moltype')then
        read(fredat(2),*)n_moltype
        allocate(st_molid(n_moltype))
        allocate(ed_molid(n_moltype))
        allocate(molname(n_moltype))
        if(n_moltype.gt.1)then
           do i=1,n_moltype
              call rdfree(iu_in,max_ndata,fredat,end_flag)
              read(fredat(2),*)molname(i)
              read(fredat(3),*)st_molid(i)
              read(fredat(5),*)ed_molid(i)
           end do
        end if
     end if
     if(fredat(1)(1:12).eq.'nslab_for_Rg')read(fredat(2),*)nslab_Rg

     !orientation
     if(fredat(1)(1:7).eq.'oopflag')read(fredat(2),*)oopflag
     if(oopflag.eq.1)then
        if(fredat(1)(1:9).eq.'direction')read(fredat(2),*)idir
        if(fredat(1)(1:13).eq.'nslab_for_oop')read(fredat(2),*)nslab_oop
        if(fredat(1)(1:11).eq.'oop_vecfile')then
           allocate(vecfile(n_moltype))
           do i=1,n_moltype
              call rdfree(iu_in,max_ndata,fredat,end_flag)
              read(fredat(1),*)dummy_a

              if(n_moltype.ne.1)then
                 if(trim(dummy_a).ne.trim(molname(i)))then
                    print*,"ERROR: miss matching of molname"
                    stop
                 end if
              end if

              read(fredat(2),'(a)')vecfile(i)
           end do
        end if
     end if

     !CV
     if(fredat(1)(1:7).eq.'z_meas1')read(fredat(2),*)z_meas1(1)
     if(fredat(1)(1:7).eq.'z_meas1')read(fredat(4),*)z_meas1(2)
     if(fredat(1)(1:7).eq.'z_meas2')read(fredat(2),*)z_meas2(1)
     if(fredat(1)(1:7).eq.'z_meas2')read(fredat(4),*)z_meas2(2)

     if(fredat(1)(1:7).eq.'x_meas1')read(fredat(2),*)x_meas1(1)
     if(fredat(1)(1:7).eq.'x_meas1')read(fredat(4),*)x_meas1(2)
     if(fredat(1)(1:7).eq.'x_meas2')read(fredat(2),*)x_meas2(1)
     if(fredat(1)(1:7).eq.'x_meas2')read(fredat(4),*)x_meas2(2)

     !charge
     if(fredat(1)(1:16).eq.'nslab_for_charge')read(fredat(2),*)nslab_ch
     if(end_flag.eq.1)exit

     !hydrogen bond
    if(fredat(1)(1:8).eq.'o_number')read(fredat(2),*)o_number
     if(fredat(1)(1:8).eq.'h_number')read(fredat(2),*)h_number
     if(fredat(1)(1:7).eq.'os_type')read(fredat(2),*)os_type
     if(fredat(1)(1:7).eq.'hs_type')read(fredat(2),*)hs_type
     if(fredat(1)(1:11).eq.'loop_length')read(fredat(2),*)loop_length
     if(fredat(1)(1:11).eq.'loop_length')read(fredat(2),*)loop_length1
     if(fredat(1)(1:14).eq.'interface_flag')read(fredat(2),*)in_flag
     if(fredat(1)(1:16).eq.'interface_maxmol')read(fredat(2),*)in_maxmol
 !    print *,'----------------------------------------------',o_number,h_number,os_type,hs_type
  end do
  
  print'("  ->start set:",i10)',s_step
  print'("  ->end set:  ",i10)',e_step
  print'("  ->delta t(fs):",f8.2)',delta_t
  print*,'# Rg'
  print'("  ->range of Rg:",f8.2," -",f7.2," (",i4,")")',0.0,max_Rg,n_f
  print'("  ->Nmoltype: ",i10)',n_moltype
  if(n_moltype.gt.1)then
     do i=1,n_moltype
        print'("      molid(",a4,"): ",i5," -",i5)',trim(molname(i)),st_molid(i),ed_molid(i)
     end do
  end if
  print'("  ->Nslab (Rg): ",i8)',nslab_Rg
  if(oopflag.eq.1)then
     print*,'# Orientation oder parameter'
     print'("  ->direction:",i10)',idir
     print'("  ->Nslab (oop):",i8)',nslab_oop
     print*," ->vecfile:"
     do i=1,n_moltype
        print*,"     "//trim(vecfile(i))
     end do
  end if
  print*,'# Control Volume'
  print'("  ->zCV1:      ",f10.2," -",f7.2)',z_meas1(1:2)
  print'("  ->zCV2:      ",f10.2," -",f7.2)',z_meas2(1:2)
  print'("  ->xCV1:      ",f10.2," -",f7.2)',x_meas1(1:2)
  print'("  ->xCV2:      ",f10.2," -",f7.2)',x_meas2(1:2)
  print*,'# Charge distribution'
  print'("  ->Nslab (charge):",i5)',nslab_ch

  allocate(f_Rg(n_f,n_moltype))
  allocate(f_Rg_xyz(3,n_f,n_moltype)) 
  allocate(in_f_Rg(n_f,n_moltype))
  allocate(in_f_Rg_xyz(3,n_f,n_moltype))
  allocate(hydron_acf_list(loop_length1))
  
  df_Rg=max_Rg/dble(n_f)

  print*,''
!### READ DCD FILE ###
  open(iu_dcd,file=trim(dcdfile),form='unformatted',status='old',iostat=ios)
  if(ios.ne.0)then
     print*,'ERROR: no file of ',trim(dcdfile)
     stop
  endif
  print*,'### Read dcd file: '//trim(dcdfile)//' ###'
! read in the header part (only once) 
  read(iu_dcd) hdr,nset,istart,nsavc,(icntrl(i),i=1,5),namnf,delta, & 
       (jcntrl(i),i=1,9) 

  read(iu_dcd) ntitle,(title(i),i=1,ntitle) 
  read(iu_dcd) natoms 

  print'(x"Number of set: ",i8)',nset
  print'(x"Step info: ",i10," - ",i10," (evey",i6," steps)")',istart,icntrl(1),nsavc    
  print'(x"Number of atoms: ",i10)',natoms
  print'(" Time average: ",i10," -",i10,", during ",f10.5,"[ns]")',&
       &s_step*nsavc+istart,e_step*nsavc+istart,dble(e_step-s_step+1)*delta_t*1.0d-6*nsavc
  ! allocate (pos(3,natoms,nset))
  allocate (pos(3,natoms,1))
  allocate(count_hyd_bond(nset))

  ! read initial (t=0) frame in the dcd 
  it=1
  read(iu_dcd) (xtable(i),i=1,6) 
  read(iu_dcd) (pos(1,i,it),i=1,natoms) 
  read(iu_dcd) (pos(2,i,it),i=1,natoms) 
  read(iu_dcd) (pos(3,i,it),i=1,natoms) 



  print*,''
  print*,'# Initial box size'
  if(xtable(1).ge.xtable(2))then
     x_max=xtable(1)
     x_min=xtable(2)
  else
     x_max=xtable(2)
     x_min=xtable(1)
  end if

  if(xtable(3).ge.xtable(4))then
     y_max=xtable(3)
     y_min=xtable(4)
  else
     y_max=xtable(4)
     y_min=xtable(3)
  end if

  if(xtable(5).ge.xtable(6))then
     z_max=xtable(5)
     z_min=xtable(6)
  else
     z_max=xtable(6)
     z_min=xtable(5)
  end if
  ! print'(x"x direction:",f10.5," - ",f10.5)',x_min,x_max
  ! print'(x"y direction:",f10.5," - ",f10.5)',y_min,y_max
  ! print'(x"z direction:",f10.5," - ",f10.5)',z_min,z_max

  ! This is initial length at step=0
  ! (Note: DCD file donot have min info, just length, z_min=0)
  lx=(x_max-x_min)
  ly=(y_max-y_min)
  lz=(z_max-z_min)

  print'(x"x direction:",f10.5)',lx
  print'(x"y direction:",f10.5)',ly
  print'(x"z direction:",f10.5)',lz

  x_half=(x_max-x_min)*0.50d0
  y_half=(y_max-y_min)*0.50d0
  z_half=(z_max-z_min)*0.50d0

!### READ BOND INFOMATION ###
  print*,''
  print*,'### Read dada file: '//trim(datafile)//' ###'
  open(iu_data,file=trim(datafile),status='old',iostat=ios)
  if(ios.ne.0)then
     print*,'ERROR: no file of ',trim(datafile)
     stop
  endif
  do 
     call rdfree(iu_data,max_ndata,fredat,end_flag)

     if((fredat(1)(1:1).eq.'#').or.(fredat(1)(1:1).eq.';')) cycle !comment line
     
     if(fredat(3)(1:3).eq.'zlo')then
        ! This is correct min/max value at this time
        ! (z_min,z_max is not correct)
        read(fredat(1),*)zmin_slab
        read(fredat(2),*)zmax_slab
     end if

     if(fredat(2)(1:5).eq.'atoms')then
        read(fredat(1),*)n_atom
        if(n_atom.ne.natoms)then
           print*,"ERROR: miss matching of n_atom & natoms"
           stop
        end if
        print'(x,i10,x,"atoms")',n_atom
        allocate(atom_type(n_atom))
        allocate(mol_id(n_atom))
        allocate(charge(n_atom))
        allocate(cnct(4,n_atom))
        allocate(mc_count_atom(n_atom))
        cnct=-100
     end if

     if(fredat(2)(1:4).eq.'atom'.and.fredat(3)(1:5).eq.'types')then
        read(fredat(1),*)n_atom_type
        print'(x,i10,x,"atom types")',n_atom_type
        allocate(mass_type(n_atom_type))
     end if
     
     if(fredat(2)(1:5).eq.'bonds')then
        read(fredat(1),*)n_bond
        print'(x,i10,x,"bonds")',n_bond
        allocate(bond(3,n_bond))
        allocate(mc_bond(2,n_bond))
        allocate(mcskip_bond(2,n_bond))
     end if

     if(fredat(1)(1:5).eq.'Atoms')then
        read(iu_data,*)
        do i=1,n_atom
           call rdfree(iu_data,max_ndata,fredat,end_flag)
           read(fredat(1),*)j
           read(fredat(2),*)mol_id(j)
           read(fredat(3),*)atom_type(j)
           read(fredat(4),*)charge(j)
        end do
     end if

     if(fredat(1)(1:6).eq.'Masses')then
        read(iu_data,*)
        do i=1,n_atom_type
           call rdfree(iu_data,max_ndata,fredat,end_flag)           
           read(fredat(2),*)mass_type(i)
        end do
     end if
     
     if(fredat(1)(1:5).eq.'Bonds')then
        read(iu_data,*)
        do i=1,n_bond
           call rdfree(iu_data,max_ndata,fredat,end_flag)
           read(fredat(2),*)bond(1,i)
           read(fredat(3),*)j
           read(fredat(4),*)k

           bond(2,i)=j
           bond(3,i)=k

           if(cnct(1,j).eq.-100)then
              cnct(1,j)=k
           else if(cnct(2,j).eq.-100)then
              cnct(2,j)=k
           else if(cnct(3,j).eq.-100)then
              cnct(3,j)=k
           else if(cnct(4,j).eq.-100)then
              cnct(4,j)=k
           end if

           if(cnct(1,k).eq.-100)then
              cnct(1,k)=j
           else if(cnct(2,k).eq.-100)then
              cnct(2,k)=j
           else if(cnct(3,k).eq.-100)then
              cnct(3,k)=j
           else if(cnct(4,k).eq.-100)then
              cnct(4,k)=j
           end if
           
        end do
     end if

     if(end_flag.eq.1)exit
  end do

  close(iu_data)

  n_molmax=maxval(mol_id)+10
  allocate(read2seri(0:n_molmax)) ! read mol_id -> serial mol_id
  allocate(seri2read(0:n_molmax)) ! serial mol_id -> read mol_id

  ! for initial molecule
  seri2read(1)=mol_id(1)
  read2seri(mol_id(1))=1
  count=1
  
  do i=2,n_atom
     flag=0
     do j=1,count
        if(mol_id(i).eq.seri2read(j))then
           flag=1
        end if
     end do
     
     if(flag.eq.0)then
        count=count+1
        seri2read(count)=mol_id(i)  !seri2read(1029) molecules number
        read2seri(mol_id(i))=count  !read2seri(38286) atoms number  

     end if
  end do
  
  n_mol=count
  

  print'(x,i10,x,"molecules")',n_mol
  allocate(mpos(3,n_mol))
  allocate(mass_mol(n_mol))
  allocate(mol_count(n_mol))
  allocate(Rg(n_mol))
  allocate(Rg_xyz(3,n_mol))
  allocate(ave_Rg(n_moltype))
  allocate(ave_Rg_xyz(3,n_moltype))
  allocate(count_Rg(n_moltype))
  
  allocate(in_Rg(n_moltype))
  allocate(in_Rg_xyz(3,n_moltype))
  allocate(in_count_Rg(n_moltype))

  mol_count=0
  do i=1,n_atom
        mol_id(i)=read2seri(mol_id(i))
        j=mol_id(i)
        mol_count(j)=mol_count(j)+1
     end do
     
  dummy_i=maxval(mol_count)
  allocate(mol_set(dummy_i,n_mol))
  allocate(min_mol_set(n_mol))
  allocate(atom_flag(n_atom))
  print'(x,i10,x,"max_natoms/mol")',dummy_i
  allocate(end_atom(n_mol))

  if(n_moltype.eq.1)then
     st_molid(1)=1
     ed_molid(1)=n_mol
  end if

  lz_slab=zmax_slab-zmin_slab
  print'(x,"zmin for slab:",x,f10.5)',zmin_slab
  print'(x,"zmax for slab:",x,f10.5)',zmax_slab
  print'(x,"lz   for slab:",x,f10.5)',lz_slab

  ! slab definition
  divz_Rg=lz_slab/dble(nslab_Rg)
  allocate(slab_Rg(3,nslab_Rg,n_moltype))
  allocate(slab_Rg2(3,nslab_Rg,n_moltype))
  allocate(count_smp_slab(nslab_Rg,n_moltype))

  divz_ch=lz_slab/dble(nslab_ch)
  allocate(slab_charge(nslab_ch,0:n_moltype)) 
  allocate(slab_charge2(nslab_ch,0:n_moltype)) 
  allocate(slab_charge2_tmp(nslab_ch,0:n_moltype)) 

!### READING VECFILE (OOP) ###
  if(oopflag.eq.1)then
     print*,""
     print*,"### Read vecfiles ###"
     allocate(n_veclist(n_moltype))
     max_nlist=0

     ! read n_veclist
     do ityp=1,n_moltype
        open(iu_tmp,file=trim(vecfile(ityp)),status='old',iostat=ios)
        if(ios.ne.0)then
           print*,'ERROR: no file of ',trim(vecfile(ityp))
           stop
        endif
        read(iu_tmp,*) ! skip 1st line
        read(iu_tmp,*)n_veclist(ityp)
        close(iu_tmp)
     end do


     max_nlist=maxval(n_veclist)
     allocate(veclist_i(max_nlist,n_moltype))
     allocate(veclist_j(max_nlist,n_moltype))
 

     ! read pair list
     do ityp=1,n_moltype
        open(iu_tmp,file=trim(vecfile(ityp)),status='old',iostat=ios)
        read(iu_tmp,*) ! skip 1st line
        read(iu_tmp,*) ! skip n_veclist line
        do i=1,n_veclist(ityp)
           read(iu_tmp,*)veclist_i(i,ityp),veclist_j(i,ityp)
        end do
        close(iu_tmp)
     end do

     do ityp=1,n_moltype
        print*,trim(molname(ityp))//": "//trim(vecfile(ityp))
        print*,"    nlist         i         j"
        do i=1,n_veclist(ityp)
           print'(3i10)',i,veclist_i(i,ityp),veclist_j(i,ityp)
        end do
     end do

     allocate(f_theta(n_theta,n_moltype))
 !!!!------------------------------------------------------------------    
     allocate(f_theta1(n_theta,n_moltype,max_nlist))
     allocate(f_cosine_decom(n_cosine,n_moltype,max_nlist))
     
     allocate(f_cosine1_list(n_veclist(1)))
     allocate(f_cosine2_list(n_veclist(2)))
     allocate(f_cosine1_sum(n_veclist(1)))
     allocate(f_cosine2_sum(n_veclist(2)))
     allocate(cosine1_sum(n_veclist(1)))
     allocate(cosine2_sum(n_veclist(2)))

     do i=1,n_veclist(1)     
        f_cosine1_list(i)=i
     end do
     
     do i=1,n_veclist(2)     
        f_cosine2_list(i)=i
     end do
     
 !!!!---------------------------------------------------    
     allocate(f_cosine(n_cosine,n_moltype))
     allocate(f_cosine_oh(n_cosine))
     allocate(oop(n_moltype))
     allocate(count_oop(n_moltype))
     allocate(slab_oop(nslab_oop,n_moltype))
     allocate(slab_oop2(nslab_oop,n_moltype))
     allocate(slab_tensor(9,nslab_oop,n_moltype))
     allocate(count_slab_oop(nslab_oop,n_moltype))
     divz_oop=lz_slab/dble(nslab_oop)

  else
     
     allocate(f_theta(1,1))
     allocate(f_cosine(1,1))
     allocate(oop(1))
     allocate(count_oop(1))
     allocate(slab_oop(1,1))
     allocate(slab_oop2(1,1))
     allocate(slab_tensor(1,1,1))
     allocate(count_slab_oop(1,1))
     divz_oop=lz_slab

  end if
  
   

!### READING POS IN EACH STEP ###
  mass_mol=0.0d0
  count_Rg=0
  ave_Rg=0.0d0
  ave_Rg_xyz=0.0d0
  f_Rg=0.0d0
  f_Rg_xyz=0.0d0

  in_count_Rg=0
  in_Rg=0.0d0
  in_Rg_xyz=0.0d0
  in_f_Rg=0.0d0
  in_f_Rg_xyz=0.0d0
  
  slab_Rg=0.0d0
  slab_Rg2=0.0d0
  count_smp_slab=0
  slab_charge=0.0d0
  slab_charge2=0.0d0
  count_smp=0
  oop=0.0d0
  count_oop=0
  count_oh=0
  slab_oop=0.0d0
  slab_oop2=0.0d0
  slab_tensor=0.0d0
  count_slab_oop=0
  f_theta=0.0d0
  f_theta1=0.0d0
  f_cosine=0.0d0
  f_cosine_decom=0.0d0
  os_count=0.0d0
  h2_count=0.0d0
  f_cosine_oh=0.0d0 
  cosine_oh_sum=0.0d0
  f_cosine_oh_sum=0.0d0 
  cosine1_sum=0.0d0
  f_cosine1_sum=0.0d0
  cosine2_sum=0.0d0
  f_cosine2_sum=0.0d0
  
  count_hyd_bond=0.0d0
  hydron_acf_list(1:loop_length1)=0
 
  
  ! make mol_count: number of atoms in each molecules
  ! make mol_set: list of atom-ID in each molecules
  ! min_make mol_set: minimum ID in  each moleset
  mol_count=0
  mol_set=-100
  do i=1,n_atom
     j=mol_id(i)
     !----------------------------------------------------------------------select type
     k=atom_type(i)
     if (k==1)then
        h2_count=h2_count+1
     end if
        
     if (k==3)then
        os_count=os_count+1
     end if
     
     
     mol_count(j)=mol_count(j)+1
     mol_set(mol_count(j),j)=i

  end do

  !n_typ2=(ed_molid(2)-st_molid(2))
  allocate(os(os_count))
  allocate(h2(h2_count))
  allocate(hydron_symbol(in_maxmol,in_maxmol))
  allocate(hydron_symbol_start(in_maxmol,in_maxmol))
  allocate(hydron_acf(in_maxmol,in_maxmol))
  allocate(hydron_acf_sum(nset))
 ! print *,'--------------------------------os_count',os_count,h2_count

  m=0
  n=0
  
  do i=1,n_atom
     k=atom_type(i)
     if (k==hs_type)then
        m=m+1
        h2(m)=i
     end if
     
     if (k==os_type)then
        n=n+1
        os(n)=i               
     end if
  end do

!  do i=1,os_count
 !    print *,i,h2(i),os(i)
 ! end do
  

  !-----------------------------------------------------------------------select type
  min_mol_set=minval(mol_set, dim = 1, mask = mol_set.gt.0)


  print*,''
  print*,'### Read each step ###',nset

  do it0=1,nset-1

     it=1
     read(iu_dcd) (xtable(i),i=1,6)
     read(iu_dcd) (pos(1,i,it),i=1,natoms) 
     read(iu_dcd) (pos(2,i,it),i=1,natoms) 
     read(iu_dcd) (pos(3,i,it),i=1,natoms)

     ! Update of box size
     if(xtable(1).ge.xtable(2))then
        x_max=xtable(1)
        x_min=xtable(2)
     else
        x_max=xtable(2)
        x_min=xtable(1)
     end if

     if(xtable(3).ge.xtable(4))then
        y_max=xtable(3)
        y_min=xtable(4)
     else
        y_max=xtable(4)
        y_min=xtable(3)
     end if

     if(xtable(5).ge.xtable(6))then
        z_max=xtable(5)
        z_min=xtable(6)
     else
        z_max=xtable(6)
        z_min=xtable(5)
     end if

     lx=(x_max-x_min)
     ly=(y_max-y_min)
     lz=(z_max-z_min)

     x_half=(x_max-x_min)*0.50d0
     y_half=(y_max-y_min)*0.50d0
     z_half=(z_max-z_min)*0.50d0



     ! time average: s_step - e_step
     if(it0.lt.s_step)cycle
     if(it0.gt.e_step)exit
     
     if(it0.eq.s_step)print*,'s_step',s_step*nsavc+istart,lz
     if(it0.eq.e_step)print*,'e_step',e_step*nsavc+istart,lz

     if(abs(lz-lz_slab).gt.min(divz_ch,divz_Rg,divz_oop))&
          & print'(x,a,x,f9.5,"  >",f9.5)',"!!! miss matching of lz for slab !!!",&
          & abs(lz-lz_slab),min(divz_ch,divz_Rg,divz_oop)

     !--------------------------------------count hydrogen bond-------------------by guo

     if (in_flag.eq.1) then
     
        loop=ceiling(it0/loop_length)
        loop_start=(loop-1)*loop_length1+1    
        count_hyd_bond(it0)=0 !bug


        hydron_symbol(1:in_maxmol,1:in_maxmol)=0.0d0
        if (it0==loop_start)then
           hydron_symbol_start(1:in_maxmol,1:in_maxmol)=0.0d0
        end if
        hydron_acf(1:in_maxmol,1:in_maxmol)=0.0d0     
        hydron_acf_sum(it0)=0.0d0

        !     print *,nset,loop,loop_start
        !     goto 30

        do i=st_molid(2),ed_molid(2)
              !print *,min_mol_set(i)
              o1=min_mol_set(i)+o_number-1 !start from 1,->> -1 for O 
              h1=min_mol_set(i)+h_number-1 !H

              !print *,atom_type(o1),atom_type(h1)
              pos_o1(1:3)=pos(1:3,o1,it) !alcohol-o
              pos_h1(1:3)=pos(1:3,h1,it) !alcohol-h
              
              

              do j=1,os_count
                 o11=os(j) !si-o
                 h11=h2(j) !si-h
                 !print *,atom_type(o11),atom_type(h11)
                 pos_o11(1:3)=pos(1:3,o11,it)
                 pos_h11(1:3)=pos(1:3,h11,it)

                 dis_ao_ah=sqrt((pos_o1(1)-pos_h1(1))**2+(pos_o1(2)-pos_h1(2))**2+(pos_o1(3)-pos_h1(3))**2) !distance of O(alcohol)-H(alcohol)
                 dis_so_ah=sqrt((pos_o11(1)-pos_h1(1))**2+(pos_o11(2)-pos_h1(2))**2+(pos_o11(3)-pos_h1(3))**2) !distance of O(si)-H(alcohol)
                 dis_so_sh=sqrt((pos_o11(1)-pos_h11(1))**2+(pos_o11(2)-pos_h11(2))**2+(pos_o11(3)-pos_h11(3))**2) !distance of O(si)-H(si)
                 dis_ao_sh=sqrt((pos_o1(1)-pos_h11(1))**2+(pos_o1(2)-pos_h11(2))**2+(pos_o1(3)-pos_h11(3))**2) !distance of O(alcohol)-H(si)

                 vec_ao_ah(1:3)=pos_o1(1:3)-pos_h1(1:3)
                 vec_so_ah(1:3)=pos_o11(1:3)-pos_h1(1:3)
                 vec_so_sh(1:3)=pos_o11(1:3)-pos_h11(1:3)
                 vec_ao_sh(1:3)=pos_o1(1:3)-pos_h11(1:3)
                 

                 !cosa=(OH*OO)/(|OH|*|OO|)
                 cos_oh1=(vec_ao_ah(1)*vec_so_ah(1)+vec_ao_ah(2)*vec_so_ah(2)+vec_ao_ah(3)*vec_so_ah(3))/(dis_so_ah*dis_ao_ah)
                 cos_oh2=(vec_so_sh(1)*vec_ao_sh(1)+vec_so_sh(2)*vec_ao_sh(2)+vec_so_sh(3)*vec_ao_sh(3))/(dis_ao_sh*dis_so_sh)


                 !if(dis_oo<3.5d0.and.dis_hh<2.45d0.and.(abs(cos_oh1)<(pi/6.0d0).or.abs(cos_oh2)<(pi/6.0d0)))then
                ! if(dis_ao_sh<2.75d0 .and.cos_oh2 <-0.642d0)then
                 if(pos_o11(3)>100 .and. dis_ao_sh<2.75d0 .and. cos_oh2 < -0.642d0)then   !
                   ! print*,'cosoh1=', cos_oh1

                    count_hyd_bond(it0)=count_hyd_bond(it0)+1
                    hydron_symbol(i,j)=1
                    if (it0==loop_start)then
                       hydron_symbol_start(i,j)=hydron_symbol(i,j)
                    end if
                    hydron_acf(i,j)=hydron_symbol(i,j)*hydron_symbol_start(i,j)! 1*n (current * start)
                    hydron_acf_sum(it0)=hydron_acf_sum(it0)+hydron_acf(i,j)
                    index=it0-(loop-1)*loop_length1
                    hydron_acf_list(index)=hydron_acf_list(index)+hydron_acf(i,j)
                     print *,it0,i,j,count_hyd_bond(it0),hydron_acf(i,j),hydron_acf_list(index)
                    !  print *,cos_oh1,'hyd-bond',it0,i,dis_oo,dis_hh,count_hyd_bond(it0)
                    !  exit
                 end if
              end do

        end do
    !    print *,'hyd-bond',it0,count_hyd_bond(it0),loop,index,hydron_acf_sum(it0)
        
     end if !inflag
     
    
  !  30 continue 
!--------------------------------------count justice hydrogen bond-------------------by guo

     ! counter of number of sampling
     count_smp=count_smp+1

     slab_charge2_tmp=0.0d0
     ! sampling total charge in slab

     
     do i=1,n_atom
        pos_z=pos(3,i,it)
        if(pos_z.lt.zmin_slab)pos_z=pos_z+lz_slab
        if(pos_z.gt.zmax_slab)pos_z=pos_z-lz_slab
        k=int((pos_z-zmin_slab)/divz_ch)+1
        if(nslab_ch.gt.1.and.k.lt.1)print*,'!!! under slab(ch) !!!',i
        if(nslab_ch.gt.1.and.k.gt.nslab_ch)print*,'!!! over slab(ch) !!!',i
        if(k.lt.1)cycle
        if(k.gt.nslab_ch)cycle
        slab_charge(k,0)=slab_charge(k,0)+charge(i)
        slab_charge2_tmp(k,0)=slab_charge2_tmp(k,0)+charge(i)
     end do

     ! sampling charge of each moltype in slab
     do ityp=1,n_moltype
        do j=st_molid(ityp),ed_molid(ityp)
           do j2=1,mol_count(j)
              i=mol_set(j2,j)

              pos_z=pos(3,i,it)
              if(pos_z.lt.zmin_slab)pos_z=pos_z+lz_slab
              if(pos_z.gt.zmax_slab)pos_z=pos_z-lz_slab
              k=int((pos_z-zmin_slab)/divz_ch)+1
              if(k.lt.1)cycle
              if(k.gt.nslab_ch)cycle
              slab_charge(k,ityp)=slab_charge(k,ityp)+charge(i)
              slab_charge2_tmp(k,ityp)=slab_charge2_tmp(k,ityp)+charge(i)
           end do
        end do
     end do
     atom_flag=0

     slab_charge2=slab_charge2+slab_charge2_tmp*slab_charge2_tmp

     ! open(88,file='test0.xyz')
     ! write(88,*)mol_count(1)
     ! write(88,*)"test1"
     ! do j2=1,mol_count(1)
     !    i=mol_set(j2,1)
     !    write(88,*)'A',pos(1:3,i,it)
     ! end do
     ! close(88)


     ! tratment of pbc
     do j=1,n_mol  
        count=1

        ! for non-bonded mol (ex. wall atoms)
        if(cnct(1,mol_set(1,j)).eq.-100)then
           
           do j2=1,mol_count(j)
              i=mol_set(j2,j)

              if(pos(1,i,it).gt.lx)pos(1,i,it)=pos(1,i,it)-lx
              if(pos(1,i,it).lt.0.0d0)pos(1,i,it)=pos(1,i,it)+lx
              if(pos(2,i,it).gt.ly)pos(2,i,it)=pos(2,i,it)-ly
              if(pos(2,i,it).lt.0.0d0)pos(2,i,it)=pos(2,i,it)+ly
              if(pos(3,i,it).gt.lz)pos(3,i,it)=pos(3,i,it)-lz
              if(pos(3,i,it).lt.0.0d0)pos(3,i,it)=pos(3,i,it)+lz
              
           end do
           
           cycle

        end if

! for bonded poly-type

        do

           do j2=1,mol_count(j)
              i=mol_set(j2,j)
              if(j2.eq.1)atom_flag(i)=1
              if(atom_flag(i).eq.0)cycle ! no-sifted atom
              pos0(1)=pos(1,i,it)
              pos0(2)=pos(2,i,it)
              pos0(3)=pos(3,i,it)

              do k=1,4
                 if(cnct(k,i).eq.-100)cycle

                 i2=cnct(k,i)
                 if(atom_flag(i2).eq.1)cycle

                 dx=pos(1,i2,it)-pos0(1)
                 dy=pos(2,i2,it)-pos0(2)
                 dz=pos(3,i2,it)-pos0(3)

                 if(dx.ge.x_half)then
                    pos(1,i2,it)=pos(1,i2,it)-lx
                 else if(dx.le.-x_half)then
                    pos(1,i2,it)=pos(1,i2,it)+lx
                 end if

                 if(dy.ge.y_half)then
                    pos(2,i2,it)=pos(2,i2,it)-ly
                 else if(dy.le.-y_half)then
                    pos(2,i2,it)=pos(2,i2,it)+ly
                 end if

                 if(dz.ge.z_half)then
                    pos(3,i2,it)=pos(3,i2,it)-lz
                 else if(dz.le.-z_half)then
                    pos(3,i2,it)=pos(3,i2,it)+lz
                 end if

                 atom_flag(i2)=1
                 count=count+1

              end do



           end do


           if(count.eq.mol_count(j)) exit

        end do ! infdo

     end do ! j=1,n_mol

     ! open(88,file='test.xyz')
     ! write(88,*)mol_count(1)
     ! write(88,*)"test1"
     ! do j2=1,mol_count(1)
     !    i=mol_set(j2,1)
     !    write(88,*)'A',pos(1:3,i,it)
     ! end do
     ! close(88)

!### CALCULATION OF RADIUS OF GYRATION ###
     mpos=0.0d0
     do i=1,n_atom
        if(cnct(1,i).ne.-100)then
           j=mol_id(i)
           mpos(1:3,j)=mpos(1:3,j)+mass_type(atom_type(i))*pos(1:3,i,it)
           if(it0.eq.s_step)mass_mol(j)=mass_mol(j)+mass_type(atom_type(i))

        end if
     end do
     
     ! center-of-mass
     do j=1,n_mol
        if(mass_mol(j).le.1.0d-8)then
           mpos(1:3,j)=0.0d0
        else
           mpos(1:3,j)=mpos(1:3,j)/mass_mol(j)
        end if
     end do
     

     Rg=0.0d0
     Rg_xyz=0.0d0
     do i=1,n_atom
        if(cnct(1,i).ne.-100)then
           j=mol_id(i)
     
           Rg(j)=Rg(j)+mass_type(atom_type(i))*&
                &((pos(1,i,it)-mpos(1,j))*(pos(1,i,it)-mpos(1,j))&
                &+(pos(2,i,it)-mpos(2,j))*(pos(2,i,it)-mpos(2,j))&
                &+(pos(3,i,it)-mpos(3,j))*(pos(3,i,it)-mpos(3,j)))
           
           Rg_xyz(1:3,j)=Rg_xyz(1:3,j)+mass_type(atom_type(i))*&
                &((pos(1:3,i,it)-mpos(1:3,j))*(pos(1:3,i,it)-mpos(1:3,j)))

        end if
     end do
     
     do ityp=1,n_moltype
        do j=st_molid(ityp),ed_molid(ityp)

           if(mass_mol(j).le.1.0d-8)then
              !print*,"nonbonded mol",j,st_molid(ityp),ed_molid(ityp)
              cycle
           end if
           Rg(j)=sqrt(Rg(j)/mass_mol(j))
           Rg_xyz(1:3,j)=sqrt(Rg_xyz(1:3,j)/mass_mol(j))
        
           smp_flag=0
           if(mpos(1,j).ge.x_meas1(1).and.mpos(1,j).le.x_meas1(2)) smp_flag=1
           if(mpos(1,j).ge.x_meas2(1).and.mpos(1,j).le.x_meas2(2)) smp_flag=1
           if(smp_flag.eq.0)cycle

           ! sampling Rg in slab
           pos_z=mpos(3,j)
           if(pos_z.lt.zmin_slab)pos_z=pos_z+lz_slab
           if(pos_z.gt.zmax_slab)pos_z=pos_z-lz_slab
           k=int((pos_z-zmin_slab)/divz_Rg)+1
           if(nslab_Rg.gt.1.and.k.lt.1)print*,'!!! under slab(Rg) !!!',j
           if(nslab_Rg.gt.1.and.k.gt.nslab_Rg)print*,'!!! over slab(Rg) !!!',j
           if(k.lt.1)cycle
           if(k.gt.nslab_Rg)cycle
           slab_Rg(1:3,k,ityp)=slab_Rg(1:3,k,ityp)+Rg_xyz(1:3,j)
           slab_Rg2(1:3,k,ityp)=slab_Rg2(1:3,k,ityp)+Rg_xyz(1:3,j)*Rg_xyz(1:3,j)
           count_smp_slab(k,ityp)=count_smp_slab(k,ityp)+1

           !interface----------------------------------------------------------------------------#
           if (ityp.eq.2 .and. in_flag.eq.1) then
                 o1=min_mol_set(j)+o_number-1 !start from 1,->> -1 for O 
                 h1=min_mol_set(j)+h_number-1 !H
                ! CV1 
                 if(pos(3,o1,it).ge.z_meas1(1).and.pos(3,o1,it).le.z_meas1(2))then
                    in_Rg(ityp)=in_Rg(ityp)+Rg(j)
                    in_Rg_xyz(1:3,ityp)=in_Rg_xyz(1:3,ityp)+Rg_xyz(1:3,j)
                    in_count_Rg(ityp)=in_count_Rg(ityp)+1

                    !iso
                    k=int(Rg(j)/df_Rg)+1
                    in_f_Rg(k,ityp)=in_f_Rg(k,ityp)+1.0d0
                    !x
                    k=int(Rg_xyz(1,j)/df_Rg)+1
                    in_f_Rg_xyz(1,k,ityp)=in_f_Rg_xyz(1,k,ityp)+1.0d0
                    !y
                    k=int(Rg_xyz(2,j)/df_Rg)+1
                    in_f_Rg_xyz(2,k,ityp)=in_f_Rg_xyz(2,k,ityp)+1.0d0
                    !z
                    k=int(Rg_xyz(3,j)/df_Rg)+1
                    in_f_Rg_xyz(3,k,ityp)=in_f_Rg_xyz(3,k,ityp)+1.0d0
                 end if
                 
                 ! CV2
                  if(pos(3,o1,it).ge.z_meas2(1).and.pos(3,o1,it).le.z_meas2(2))then
                    in_Rg(ityp)=in_Rg(ityp)+Rg(j)
                    in_Rg_xyz(1:3,ityp)=in_Rg_xyz(1:3,ityp)+Rg_xyz(1:3,j)
                    in_count_Rg(ityp)=in_count_Rg(ityp)+1

                    !iso
                    k=int(Rg(j)/df_Rg)+1
                    in_f_Rg(k,ityp)=in_f_Rg(k,ityp)+1.0d0
                    !x
                    k=int(Rg_xyz(1,j)/df_Rg)+1
                    in_f_Rg_xyz(1,k,ityp)=in_f_Rg_xyz(1,k,ityp)+1.0d0
                    !y
                    k=int(Rg_xyz(2,j)/df_Rg)+1
                    in_f_Rg_xyz(2,k,ityp)=in_f_Rg_xyz(2,k,ityp)+1.0d0
                    !z
                    k=int(Rg_xyz(3,j)/df_Rg)+1
                    in_f_Rg_xyz(3,k,ityp)=in_f_Rg_xyz(3,k,ityp)+1.0d0
                 end if
     
           end if

           ! CV1
           if(mpos(3,j).ge.z_meas1(1).and.mpos(3,j).le.z_meas1(2))then
              ave_Rg(ityp)=ave_Rg(ityp)+Rg(j)
              ave_Rg_xyz(1:3,ityp)=ave_Rg_xyz(1:3,ityp)+Rg_xyz(1:3,j)
              count_Rg(ityp)=count_Rg(ityp)+1

              !iso
              k=int(Rg(j)/df_Rg)+1
              f_Rg(k,ityp)=f_Rg(k,ityp)+1.0d0
              !x
              k=int(Rg_xyz(1,j)/df_Rg)+1
              f_Rg_xyz(1,k,ityp)=f_Rg_xyz(1,k,ityp)+1.0d0
              !y
              k=int(Rg_xyz(2,j)/df_Rg)+1
              f_Rg_xyz(2,k,ityp)=f_Rg_xyz(2,k,ityp)+1.0d0
              !z
              k=int(Rg_xyz(3,j)/df_Rg)+1
              f_Rg_xyz(3,k,ityp)=f_Rg_xyz(3,k,ityp)+1.0d0
           end if

           ! CV2
           if(mpos(3,j).ge.z_meas2(1).and.mpos(3,j).le.z_meas2(2))then
              ave_Rg(ityp)=ave_Rg(ityp)+Rg(j)
              ave_Rg_xyz(1:3,ityp)=ave_Rg_xyz(1:3,ityp)+Rg_xyz(1:3,j)
              count_Rg(ityp)=count_Rg(ityp)+1

              !iso
              k=int(Rg(j)/df_Rg)+1
              f_Rg(k,ityp)=f_Rg(k,ityp)+1.0d0
              !x
              k=int(Rg_xyz(1,j)/df_Rg)+1
              f_Rg_xyz(1,k,ityp)=f_Rg_xyz(1,k,ityp)+1.0d0
              !y
              k=int(Rg_xyz(2,j)/df_Rg)+1
              f_Rg_xyz(2,k,ityp)=f_Rg_xyz(2,k,ityp)+1.0d0
              !z
              k=int(Rg_xyz(3,j)/df_Rg)+1
              f_Rg_xyz(3,k,ityp)=f_Rg_xyz(3,k,ityp)+1.0d0
           end if

        end do
     end do

     

     !### CALCULATION OF OOP ###
     if(oopflag.eq.1)then
        
!#----------------------------------------Si-OH-cos-----------------------------------------#
        do j=1,os_count
           o11=os(j) !si-o
           h11=h2(j) !si-h          
           pos_o11(1:3)=pos(1:3,o11,it)
           pos_h11(1:3)=pos(1:3,h11,it)
                    
           vec_oh11(1:3)=pos_o11(1:3)-pos_h11(1:3) !#vector of OH(si)
           dis_oh11=sqrt((pos_o11(1)-pos_h11(1))**2+(pos_o11(2)-pos_h11(2))**2+(pos_o11(3)-pos_h11(3))**2) !distance of OH(si)
           cosine_oh=abs(vec_oh11(idir)/dis_oh11)
 
           i3=int(cosine_oh/dcosine)+1
           if(i3.gt.n_cosine)i3=n_cosine
           if(i3.lt.1)i3=1
         !  print *,'-----------------------si-oh: ',atom_type(o11),atom_type(h11),pos_o11(3),pos_h11(3)
           f_cosine_oh(i3)=f_cosine_oh(i3)+1.0d0
           count_oh=count_oh+1
        !    print*,cosine_oh,i3,f_cosine_oh(i3)                         
        end do
        ! print*,'oh',f_cosine_oh(10)/count_oh
!#----------------------------------------Si-OH-cos-----------------------------------------#
        
        do ityp=1,n_moltype
           do imol=st_molid(ityp),ed_molid(ityp)

              if(mass_mol(imol).le.1.0d-8)cycle ! skip if non-bonded mol
              
              do k=1,n_veclist(ityp)
                ! print *,"------------------------------------------min-mol-set",min_mol_set(imol)
                !print *,"------------------------------------------veclist",veclist_i(k,ityp)
                 i=min_mol_set(imol)+veclist_i(k,ityp)-1
                 j=min_mol_set(imol)+veclist_j(k,ityp)-1

                 oop_vec(1:3)=pos(1:3,i,it)-pos(1:3,j,it)

                 !### for cheack (isotoropic)
                 ! call random_number(rnd)
                 ! cos_theta=1.0d0-2.0d0*rnd
                 ! sin_theta=sqrt(1.0d0-cos_theta*cos_theta)
                 ! call random_number(rnd)
                 ! phi=2.0d0*pi*rnd
                 ! oop_vec(1)=sin_theta*cos(phi)
                 ! oop_vec(2)=sin_theta*sin(phi)
                 ! oop_vec(3)=cos_theta
                 ! write(123,*)oop_vec(1:3)

                 !### for check (uniaxial)
                 ! oop_vec(1)=1
                 ! oop_vec(2)=0
                 ! oop_vec(3)=0

                 oop_vec(4)=sqrt(oop_vec(1)*oop_vec(1)+oop_vec(2)*oop_vec(2)+oop_vec(3)*oop_vec(3))

                 cosine=abs(oop_vec(idir))/oop_vec(4)
                 theta=acos(cosine)

                 pos_z=(pos(3,i,it)+pos(3,j,it))*0.5d0
                 if(pos_z.lt.zmin_slab)pos_z=pos_z+lz_slab
                 if(pos_z.gt.zmax_slab)pos_z=pos_z-lz_slab

                 k2=int((pos_z-zmin_slab)/divz_oop)+1
                 if(nslab_oop.gt.1.and.k2.lt.1)print*,'!!! under slab(oop) !!!',j
                 if(nslab_oop.gt.1.and.k2.gt.nslab_oop)print*,'!!! over slab(oop) !!!',j
                 if(k2.lt.1)cycle
                 if(k2.gt.nslab_oop)cycle
                 slab_oop(k2,ityp)=slab_oop(k2,ityp)+3.0d0*cosine*cosine-1.0d0
                 slab_oop2(k2,ityp)=slab_oop2(k2,ityp)+(3.0d0*cosine*cosine-1.0d0)*(3.0d0*cosine*cosine-1.0d0)
                 count_slab_oop(k2,ityp)=count_slab_oop(k2,ityp)+1

                 ! Tensor order parameter
                 ! dyadic (uu)ij, u: unit vector of oop_vec
                 slab_tensor(1,k2,ityp)=slab_tensor(1,k2,ityp)+3.0d0*oop_vec(1)*oop_vec(1)/(oop_vec(4)*oop_vec(4))-1.0d0
                 slab_tensor(2,k2,ityp)=slab_tensor(2,k2,ityp)+3.0d0*oop_vec(1)*oop_vec(2)/(oop_vec(4)*oop_vec(4))
                 slab_tensor(3,k2,ityp)=slab_tensor(3,k2,ityp)+3.0d0*oop_vec(1)*oop_vec(3)/(oop_vec(4)*oop_vec(4))
                 slab_tensor(4,k2,ityp)=slab_tensor(4,k2,ityp)+3.0d0*oop_vec(2)*oop_vec(1)/(oop_vec(4)*oop_vec(4))
                 slab_tensor(5,k2,ityp)=slab_tensor(5,k2,ityp)+3.0d0*oop_vec(2)*oop_vec(2)/(oop_vec(4)*oop_vec(4))-1.0d0
                 slab_tensor(6,k2,ityp)=slab_tensor(6,k2,ityp)+3.0d0*oop_vec(2)*oop_vec(3)/(oop_vec(4)*oop_vec(4))
                 slab_tensor(7,k2,ityp)=slab_tensor(7,k2,ityp)+3.0d0*oop_vec(3)*oop_vec(1)/(oop_vec(4)*oop_vec(4))
                 slab_tensor(8,k2,ityp)=slab_tensor(8,k2,ityp)+3.0d0*oop_vec(3)*oop_vec(2)/(oop_vec(4)*oop_vec(4))
                 slab_tensor(9,k2,ityp)=slab_tensor(9,k2,ityp)+3.0d0*oop_vec(3)*oop_vec(3)/(oop_vec(4)*oop_vec(4))-1.0d0

                 ! CV1
                 if(pos(3,i,it).ge.z_meas1(1).and.&
                      & pos(3,i,it).le.z_meas1(2).and.&
                      & pos(3,j,it).ge.z_meas1(1).and.&
                      & pos(3,j,it).le.z_meas1(2))then

                    i2=int(theta/dtheta)+1
                    if(i2.gt.n_theta)i2=n_theta
                    if(i2.lt.1)i2=1
                    f_theta(i2,ityp)=f_theta(i2,ityp)+1.0d0

                    i3=int(cosine/dcosine)+1
                    if(i3.gt.n_cosine)i3=n_cosine
                    if(i3.lt.1)i3=1
                    f_cosine(i3,ityp)=f_cosine(i3,ityp)+1.0d0
  !!!!-------------------------------------------------------------------------                  
                    f_cosine_decom(i3,ityp,k)=f_cosine_decom(i3,ityp,k)+1.0d0

                    
                    oop(ityp)=oop(ityp)+3.0d0*cosine*cosine-1.0d0
                    count_oop(ityp)=count_oop(ityp)+1
                    
                 end if
                 

                 ! CV2
                 if(pos(3,i,it).ge.z_meas2(1).and.&
                      & pos(3,i,it).le.z_meas2(2).and.&
                      & pos(3,j,it).ge.z_meas2(1).and.&
                      & pos(3,j,it).le.z_meas2(2))then

                    i2=int(theta/dtheta)+1
                    if(i2.gt.n_theta)i2=n_theta
                    if(i2.lt.1)i2=1
                    f_theta(i2,ityp)=f_theta(i2,ityp)+1.0d0
                   

                    i3=int(cosine/dcosine)+1
                    if(i3.gt.n_cosine)i3=n_cosine
                    if(i3.lt.1)i3=1
                    f_cosine(i3,ityp)=f_cosine(i3,ityp)+1.0d0
                    
                    f_cosine_decom(i3,ityp,k)=f_cosine_decom(i3,ityp,k)+1.0d0
  !!!!-----------------------------------------------------------------------
                    oop(ityp)=oop(ityp)+3.0d0*cosine*cosine-1.0d0
                    count_oop(ityp)=count_oop(ityp)+1
                 end if

              end do !k

           end do !imol

        end do !ityp
   !     print*,'--------------------------',it0,f_cosine(10,2)/count_oop(2)

        ! ! calc check vector (isotropic)
        ! ityp=1
        ! do i=1,n_check
        !    call random_number(rnd)
        !    cos_theta=1.0d0-2.0d0*rnd
        !    sin_theta=sqrt(1.0d0-cos_theta*cos_theta)
        !    call random_number(rnd)
        !    phi=2.0d0*pi*rnd


        !    oop_vec(1)=r_check*sin_theta*cos(phi)
        !    oop_vec(2)=r_check*sin_theta*sin(phi)
        !    oop_vec(3)=r_check*cos_theta

        !    oop_vec(4)=sqrt(oop_vec(1)*oop_vec(1)+oop_vec(2)*oop_vec(2)+oop_vec(3)*oop_vec(3))
        !    cosine=abs(oop_vec(idir))/oop_vec(4)
        !    theta=acos(cosine)

        !    i2=int(theta/dtheta)+1
        !    if(i2.gt.n_theta)i2=n_theta
        !    if(i2.lt.1)i2=1
        !    f_theta(i2,ityp)=f_theta(i2,ityp)+1.0d0

        !    i3=int(cosine/dcosine)+1
        !    if(i3.gt.n_cosine)i3=n_cosine
        !    if(i3.lt.1)i3=1
        !    f_cosine(i3,ityp)=f_cosine(i3,ityp)+1.0d0

        !    oop(ityp)=oop(ityp)+3.0d0*cosine*cosine-1.0d0
        !    count_oop(ityp)=count_oop(ityp)+1

        ! end do

     end if ! oopflag

  end do ! it0=1,nset

!### TIME AVERAGE ###
  print*,''
  print*,'### Time average ###'
  ! average of Rg
  do ityp=1,n_moltype

     print*,'!!',ave_Rg(ityp), ave_Rg_xyz(1,ityp)+ave_Rg_xyz(2,ityp)+ave_Rg_xyz(3,ityp)
     
     ave_Rg(ityp)=ave_Rg(ityp)/dble(count_Rg(ityp))
     ave_Rg_xyz(1:3,ityp)=ave_Rg_xyz(1:3,ityp)/dble(count_Rg(ityp))
     ! distribution of Rg
     f_Rg(1:n_f,ityp)=f_Rg(1:n_f,ityp)/(df_Rg*dble(count_Rg(ityp)))
     f_Rg_xyz(1,1:n_f,ityp)=f_Rg_xyz(1,1:n_f,ityp)/(df_Rg*dble(count_Rg(ityp)))
     f_Rg_xyz(2,1:n_f,ityp)=f_Rg_xyz(2,1:n_f,ityp)/(df_Rg*dble(count_Rg(ityp)))
     f_Rg_xyz(3,1:n_f,ityp)=f_Rg_xyz(3,1:n_f,ityp)/(df_Rg*dble(count_Rg(ityp)))

     in_Rg(ityp)=in_Rg(ityp)/dble(in_count_Rg(ityp))
     in_Rg_xyz(1:3,ityp)=in_Rg_xyz(1:3,ityp)/dble(in_count_Rg(ityp))
     ! distribution of Rg
     in_f_Rg(1:n_f,ityp)=in_f_Rg(1:n_f,ityp)/(df_Rg*dble(in_count_Rg(ityp)))
     in_f_Rg_xyz(1,1:n_f,ityp)=in_f_Rg_xyz(1,1:n_f,ityp)/(df_Rg*dble(in_count_Rg(ityp)))
     in_f_Rg_xyz(2,1:n_f,ityp)=in_f_Rg_xyz(2,1:n_f,ityp)/(df_Rg*dble(in_count_Rg(ityp)))
     in_f_Rg_xyz(3,1:n_f,ityp)=in_f_Rg_xyz(3,1:n_f,ityp)/(df_Rg*dble(in_count_Rg(ityp)))

     dummy=0
     do i=1,n_f
        dummy=dummy+f_Rg(i,ityp)*df_Rg
     end do
     print*,'average Rg: ',ave_Rg(ityp),dummy
     
     dummy=0
     do i=1,n_f
        dummy=dummy+f_Rg_xyz(1,i,ityp)*df_Rg
     end do
     print*,'average Rg_x: ',ave_Rg_xyz(1,ityp),dummy
     
     dummy=0
     do i=1,n_f
        dummy=dummy+f_Rg_xyz(2,i,ityp)*df_Rg
     end do
     print*,'average Rg_y: ',ave_Rg_xyz(2,ityp),dummy
     
     dummy=0
     do i=1,n_f
        dummy=dummy+f_Rg_xyz(3,i,ityp)*df_Rg
     end do
     print*,'average Rg_z: ',ave_Rg_xyz(3,ityp),dummy
     print*,'ave_Rg vs sqrt(ave_Rg_xyz): ',ave_Rg(ityp)**2,ave_Rg_xyz(1,ityp)**2+ave_Rg_xyz(2,ityp)**2+ave_Rg_xyz(3,ityp)**2

     print*,''
  end do

  open(iu_out,file='f_Rg.dat')
  write(iu_out,*)"# Frequency of Rg"
  write(iu_out,'(a,3e24.16)')' # average Rg[A]:',ave_Rg(1:n_moltype)
  write(iu_out,'(a,9e24.16)')' # average Rg_xyz[A]:',(ave_Rg_xyz(1:3,ityp),ityp=1,n_moltype)
  write(iu_out,*)'# nsample:   ',(count_Rg(ityp),ityp=1,n_moltype)
  if(n_moltype.eq.1)then
     write(iu_out,*)'Rg[A] f(Rg) f(Rg_x) f(Rg_y) f(Rg_z)'
  else
     write(iu_out,'(a,3(2x,a,"(r)",2x,a,"(x)",2x,a,"(y)",2x,a,"(z)"))')" Rg[A]",&
          &(trim(molname(ityp)),trim(molname(ityp)),trim(molname(ityp)),trim(molname(ityp)),ityp=1,n_moltype)
  end if
  do k=1,n_f
     write(iu_out,'(100e24.16)')df_Rg*(k-1)+df_Rg*0.50d0,(f_Rg(k,ityp),f_Rg_xyz(1:3,k,ityp),ityp=1,n_moltype)
  end do
  close(iu_out)
  print*,'  --> output file: f_Rg.dat'
  
!-------------------------------------------------INTERFACE-------------------------------------------------#
  open(iu_out,file='in_f_Rg.dat')
  write(iu_out,*)"# Frequency of Rg"
  write(iu_out,'(a,3e24.16)')' # average Rg[A]:',in_Rg(1:n_moltype)
  write(iu_out,'(a,9e24.16)')' # average Rg_xyz[A]:',(in_Rg_xyz(1:3,ityp),ityp=1,n_moltype)
  write(iu_out,*)'# nsample:   ',(in_count_Rg(ityp),ityp=1,n_moltype)
  if(n_moltype.eq.1)then
     write(iu_out,*)'Rg[A] f(Rg) f(Rg_x) f(Rg_y) f(Rg_z)'
  else
     write(iu_out,'(a,3(2x,a,"(r)",2x,a,"(x)",2x,a,"(y)",2x,a,"(z)"))')" Rg[A]",&
          &(trim(molname(ityp)),trim(molname(ityp)),trim(molname(ityp)),trim(molname(ityp)),ityp=1,n_moltype)
  end if
  do k=1,n_f
     write(iu_out,'(100e24.16)')df_Rg*(k-1)+df_Rg*0.50d0,(in_f_Rg(k,ityp),in_f_Rg_xyz(1:3,k,ityp),ityp=1,n_moltype)
  end do
  close(iu_out)
  print*,'  --> output file: in_f_Rg.dat'

!----------------------------------------------------------------------------------------------------------#
  open(iu_out,file='slab_Rg.dat')
  write(iu_out,*)"# Local Rg in each slab"
  if(n_moltype.eq.1)then
     write(iu_out,*)"z Rg(x) Rg(y) Rg(z) sterr(x) sterr(y) sterr(z) count"
  else
     write(iu_out,'(a,2(2x,a,"(x)",2x,a,"(y)",2x,a,"(z)"),&
          &2(2x,a,"_sterr(x)",2x,a,"_sterr(y)",2x,a,"_sterr(z)"),2(2x,a,"_count"))')" z",&
          &(trim(molname(ityp)),trim(molname(ityp)),trim(molname(ityp)),ityp=1,n_moltype),&
          &(trim(molname(ityp)),trim(molname(ityp)),trim(molname(ityp)),ityp=1,n_moltype),&
          &(trim(molname(ityp)),ityp=1,n_moltype)
  end if
  do k=1,nslab_Rg
     do ityp=1,n_moltype
        slab_Rg(1:3,k,ityp)=slab_Rg(1:3,k,ityp)/dble(count_smp_slab(k,ityp))
        slab_Rg2(1:3,k,ityp)=slab_Rg2(1:3,k,ityp)/dble(count_smp_slab(k,ityp)-1)
        ! st dev^2
        slab_Rg2(1:3,k,ityp)=slab_Rg2(1:3,k,ityp)&
             & -slab_Rg(1:3,k,ityp)*slab_Rg(1:3,k,ityp)*dble(count_smp_slab(k,ityp))/dble(count_smp_slab(k,ityp)-1)
        ! st err
        slab_Rg2(1:3,k,ityp)=sqrt(slab_Rg2(1:3,k,ityp)/dble(count_smp_slab(k,ityp)))
        if(count_smp_slab(k,ityp).le.10)slab_Rg(1:3,k,ityp)=0.0d0
        if(count_smp_slab(k,ityp).le.10)slab_Rg2(1:3,k,ityp)=0.0d0
     end do
     write(iu_out,'(100e24.16)')divz_Rg*(k-1)+divz_Rg*0.50d0+zmin_slab,&
          & (slab_Rg(1:3,k,ityp),ityp=1,n_moltype),&
          & (slab_Rg2(1:3,k,ityp),ityp=1,n_moltype),&
          & (dble(count_smp_slab(k,ityp)),ityp=1,n_moltype)
  end do
  close(iu_out)
  print*,'  --> output file: slab_Rg.dat'


  open(iu_out,file='slab_charge.dat')
  if(n_moltype.eq.1)then
     write(iu_out,*)'z total mol'
  else
     write(iu_out,'(a,2(2x,a),2(2x,a,"_sterr"))')" z total",&
          & (trim(molname(ityp)),ityp=1,n_moltype),&
          & (trim(molname(ityp)),ityp=1,n_moltype)
  end if
  do k=1,nslab_ch
     slab_charge(k,0:n_moltype)=slab_charge(k,0:n_moltype)/dble(count_smp)
     slab_charge2(k,0:n_moltype)=slab_charge2(k,0:n_moltype)/dble(count_smp-1)
     ! st dev^2
     slab_charge2(k,0:n_moltype)=slab_charge2(k,0:n_moltype)&
          & -slab_charge(k,0:n_moltype)*slab_charge(k,0:n_moltype)*dble(count_smp)/dble(count_smp-1)
     ! st err
     slab_charge2(k,0:n_moltype)=sqrt(slab_charge2(k,0:n_moltype)/dble(count_smp))
     write(iu_out,'(100e24.16)')divz_ch*(k-1)+divz_ch*0.50d0+zmin_slab,&
          & slab_charge(k,0:n_moltype),slab_charge2(k,0:n_moltype)
  end do
  close(iu_out)
  print*,'  --> output file: slab_charge.dat'


! ### OOP ###
  if(oopflag.eq.1)then

     do ityp=1,n_moltype
        ! average of oop
        oop(ityp)=oop(ityp)*0.50d0/dble(count_oop(ityp))

        ! normalization

        f_theta(1:n_theta,ityp)=f_theta(1:n_theta,ityp)/(pi*sin(dtheta*(i-1)+dtheta*0.50d0)*dtheta)*2.0d0/count_oop(ityp)
        f_cosine(1:n_cosine,ityp)=f_cosine(1:n_cosine,ityp)/(dcosine*dble(count_oop(ityp)))
        print*,'--------------------------',f_cosine(10,2),count_oop
        !---------------------------------------------
        f_cosine_decom(1:n_cosine,ityp,1:max_nlist)=f_cosine_decom(1:n_cosine,ityp,1:max_nlist)/(dcosine*dble(count_oop(ityp)))       
        f_cosine_oh(1:n_cosine)=f_cosine_oh(1:n_cosine)/(dcosine*dble(count_oh))

        dummy=0
        do i=1,n_theta
           dummy=dummy+f_theta(i,ityp)*dtheta
        end do
        print*,'integral of f_theta',dummy

        dummy=0
        do i=1,n_cosine
           dummy=dummy+f_cosine(i,ityp)*dcosine
        end do
        print*,'integral of f_cosine',dummy

     end do
     
     do i=1,n_cosine
           dummy=dummy+f_cosine_oh(i)*dcosine
        end do
        print*,'integral of f_cosine_oh',dummy

     open(iu_out,file='f_theta.dat')
     write(iu_out,*)'# distribution of theta'
     write(iu_out,'(a,3e24.16)')' # average oop:',oop(1:n_moltype)
     write(iu_out,*)'# nsample: ',count_oop(1:n_moltype)
     if(n_moltype.eq.1)then
        write(iu_out,*)'theta  f(theta)'
     else
        write(iu_out,'(a,3(2x,a))')" theta",(trim(molname(ityp)),ityp=1,n_moltype)
     end if
     do k=1,n_theta
        write(iu_out,'(100e24.16)')dtheta*(k-1)+dtheta*0.50d0,(f_theta(k,ityp),ityp=1,n_moltype)
     end do
     close(iu_out)
     print*,'  --> output file: f_tehta.dat'

     open(iu_out,file='f_cosine.dat')
     write(iu_out,*)'# distribution of cosine'
     write(iu_out,'(a,3e24.16)')' # average oop:',oop(1:n_moltype)
     write(iu_out,*)'# nsample: ',count_oop(1:n_moltype)
     if(n_moltype.eq.1)then
        write(iu_out,*)'cosine  f(cosine)'
     else
        write(iu_out,'(a,3(2x,a))')" cosine",(trim(molname(ityp)),ityp=1,n_moltype)
     end if
     do k=1,n_cosine
        write(iu_out,'(100e24.16)')dcosine*(k-1)+dcosine*0.50d0,(f_cosine(k,ityp),ityp=1,n_moltype)
     end do
     close(iu_out)
     print*,'  --> output file: f_cosine.dat'
     !---------------------------------------------------------------
!#cos1
     open(iu_out,file='f_cosine1.dat')
     write(iu_out,*)'# distribution of separate cosine'
     write(iu_out,'(a,3e24.16)')' # average oop:',oop(1:n_moltype)
     write(iu_out,*)'# nsample: ',count_oop(1:n_moltype)
     if(n_moltype.eq.1)then
        write(iu_out,*)'cosine  f(cosine)'
     else
        
        write(iu_out,*)" cosine_",trim(molname(1)),(f_cosine1_list(i),i=1,n_veclist(1))
  !      write(iu_out,'(a,3(2x,a))')" cosine",i=1,max_nlist
     end if
     do k=1,n_cosine
        write(iu_out,'(100e24.16)')dcosine*(k-1)+dcosine*0.50d0,(f_cosine_decom(k,1,i),i=1,n_veclist(1))
        do i=1,n_veclist(1)
           cosine1_sum(i)=cosine1_sum(i)+(dcosine*(k-1)+dcosine*0.50d0)*f_cosine_decom(k,1,i)
           f_cosine1_sum(i)=f_cosine1_sum(i)+f_cosine_decom(k,1,i)
        end do
     end do
     write(iu_out,*)'# cosine_ave:',(cosine1_sum(i)/f_cosine1_sum(i),i=1,n_veclist(1))
     
     close(iu_out)
     print*,'  --> output file: f_cosine1.dat'
     
!#cos2
     open(iu_out,file='f_cosine2.dat')
     write(iu_out,*)'# distribution of separate cosine'
     write(iu_out,'(a,3e24.16)')' # average oop:',oop(1:n_moltype)
     write(iu_out,*)'# nsample: ',count_oop(1:n_moltype)
     if(n_moltype.eq.1)then
        write(iu_out,*)'cosine  f(cosine)'
     else
        
        write(iu_out,*)" cosine_",trim(molname(2)),(f_cosine2_list(i),i=1,n_veclist(2))
  !      write(iu_out,'(a,3(2x,a))')" cosine",i=1,max_nlist
     end if

     do k=1,n_cosine
         write(iu_out,'(100e24.16)')dcosine*(k-1)+dcosine*0.50d0,(f_cosine_decom(k,2,i),i=1,n_veclist(2)) 
        do i=1,n_veclist(2)
           cosine2_sum(i)=cosine2_sum(i)+(dcosine*(k-1)+dcosine*0.50d0)*f_cosine_decom(k,2,i)
           f_cosine2_sum(i)=f_cosine2_sum(i)+f_cosine_decom(k,2,i)
        end do    
     end do
      write(iu_out,*)'# cosine_ave:',(cosine2_sum(i)/f_cosine2_sum(i),i=1,n_veclist(2))
     
     close(iu_out)
     print*,'  --> output file: f_cosine2.dat'
     
!#cos-oh(si)
     open(iu_out,file='f_cosine-oh-si.dat')
     write(iu_out,*)'# distribution of separate cosine'
     write(iu_out,*)'# nsample: ',count_oh

     write(iu_out,*)'cosine  f(cosine)'

 
     do k=1,n_cosine
        write(iu_out,'(100e24.16)')dcosine*(k-1)+dcosine*0.50d0,(f_cosine_oh(k))
        cosine_oh_sum=cosine_oh_sum+(dcosine*(k-1)+dcosine*0.50d0)*f_cosine_oh(k)
        f_cosine_oh_sum=f_cosine_oh_sum+f_cosine_oh(k)
     end do
     
     write(iu_out,*)'# cosine_ave:',cosine_oh_sum/f_cosine_oh_sum
     print *,'# cosine_ave:',cosine_oh_sum/f_cosine_oh_sum    
        
     close(iu_out)
     print*,'  --> output file: f_cosine-oh-si.dat'

!----------------------------------------------------------------------------------

     open(iu_out,file='slab_oop.dat')
     write(iu_out,*)"# Local orientation order parameter in each slab"
     if(n_moltype.eq.1)then
        write(iu_out,*)"z oop sterr"
     else
        write(iu_out,'(a,2(2x,a),2(2x,a,"_sterr"))')" z",&
             & (trim(molname(ityp)),ityp=1,n_moltype),&
             & (trim(molname(ityp)),ityp=1,n_moltype)
     end if
     do k=1,nslab_oop
        do ityp=1,n_moltype
           slab_oop(k,ityp)=slab_oop(k,ityp)/dble(count_slab_oop(k,ityp))*0.5d0
           slab_oop2(k,ityp)=slab_oop2(k,ityp)/dble(count_slab_oop(k,ityp)-1)*0.5d0*0.5d0
           ! st dev^2
           slab_oop2(k,ityp)=slab_oop2(k,ityp)&
                & -slab_oop(k,ityp)*slab_oop(k,ityp)*dble(count_slab_oop(k,ityp))/dble(count_slab_oop(k,ityp)-1)
           ! st err
           slab_oop2(k,ityp)=sqrt(slab_oop2(k,ityp)/dble(count_slab_oop(k,ityp)))
           if(count_slab_oop(k,ityp).le.10)slab_oop(k,ityp)=0.0d0
           if(count_slab_oop(k,ityp).le.10)slab_oop2(k,ityp)=0.0d0
        end do
        write(iu_out,'(100e24.16)')divz_oop*(k-1)+divz_oop*0.50d0+zmin_slab,&
             & (slab_oop(k,ityp),ityp=1,n_moltype),&
             & (slab_oop2(k,ityp),ityp=1,n_moltype)
     end do
     close(iu_out)
     print*,'  --> output file: slab_oop.dat'


     open(iu_out,file='slab_tensor.dat')
     write(iu_out,*)"# Local tensor order parameter for director"
     if(n_moltype.eq.1)then
        write(iu_out,*)"z S11 S12 S13 S21 S22 S23 S31 S32 S33 sterroop"
     else
        write(iu_out,'(a,2(2x,a,"11",2x,a,"12",2x,a,"13",2x,a,"21",2x,a,"22",2x,a,"23",2x,&
             &a,"31",2x,a,"32",2x,a,"33"),2(2x,a,"_sterroop"))')&
             & " z",(trim(molname(ityp)),trim(molname(ityp)),trim(molname(ityp)),&
             & trim(molname(ityp)),trim(molname(ityp)),trim(molname(ityp)),&
             & trim(molname(ityp)),trim(molname(ityp)),trim(molname(ityp)),ityp=1,n_moltype),&
             & (trim(molname(ityp)),ityp=1,n_moltype)
     end if
     do k=1,nslab_oop
        do ityp=1,n_moltype
           slab_tensor(1:9,k,ityp)=slab_tensor(1:9,k,ityp)/dble(count_slab_oop(k,ityp))*0.5d0
           if(count_slab_oop(k,ityp).le.10)slab_tensor(1:9,k,ityp)=0.0d0
        end do
        write(iu_out,'(100e24.16)')divz_oop*(k-1)+divz_oop*0.50d0+zmin_slab,&
             & (slab_tensor(1:9,k,ityp),ityp=1,n_moltype),&
             & (slab_oop2(k,ityp),ityp=1,n_moltype)
     end do
     close(iu_out)
     print*,'  --> output file: slab_tensor.dat'

     
     open(iu_out,file='hyd_count.dat')
     write(iu_out,*)"# the number of hydrogen bonds at interface"
     write(iu_out,'(a,2x,a)')"step", "hyd_count"
     do k=1,e_step
        write(iu_out,'(i5,2x,i5)')k,count_hyd_bond(k)
     end do
     print*,'  --> output file: hyd_count.dat'

     open(iu_out,file='hyd_acf.dat')
     write(iu_out,*)"# the loop length = ",loop_length,"loop time = ",loop_length*delta_t,'fs'
     write(iu_out,'(a,2x,a)')"timestep", "hydron_bond_acf"
     do k=1,loop_length1
        write(iu_out,'(i5,2x,e24.16)')k,real(hydron_acf_list(k))/(real(hydron_acf_list(1)))
     end do
     print*,'  --> output file: hyd_acf.dat'

  end if
















100     format(i6,3(x,f8.4))
      end program calc_config

!-----------------------------------------------------------------------
subroutine rdfree(iuin, ndata, fredat,flag)
!     
!     subroutine to read in a free format 
!
!     This subroutine reads in a free format a line containing various
!     types of data separated by a space ' ' or a comma ','.
!     For instance,
!         'ATOM 500 RES 101  50.2 30.3 1.99 10e-10'  
!     DATA ARE RETURNED AS CHARACTERS
!     If ndata > actual number of data, then spaces are input in the
!     rest of fredat.

  implicit none

!     ARGUMENT

!        INPUT: 
  integer,intent(in)::     iuin          ! unit to read
  integer,intent(in)::     ndata         ! number of data to be read.

!        OUTPUT:
  character(80),intent(inout):: fredat(ndata) ! data read as characters
!                   length of each datum is assumed =<80
!                   each datum begins at the left end of data(i)
!

!     LOCAL

!        * LINE 
  integer,parameter:: lnleng =2000  ! line length
  character(len=lnleng):: line ! line is read as a character
 
!        * COUNTS
  integer::     idata         ! temporal count of data
  integer::     ldata         ! length of the datum
 
!        * FLAGS
  logical::     spcom0, spcom1 ! true  if i-1 (or i) is space
!                                                  or comma

!        * OHTERS
  character(len=1):: space,  comma,  tab
                                ! separators
  integer::       i           ! do loop index

  integer:: flag

!          +    +    +    +    +


!     --- INITIALIZE THE DATA ---

  space=' '
  comma=','
  tab='	'
  flag=1
  do i=1,ndata
     fredat(i) = space
  end do

!     --- READ THE LINE FROM IUIN AS A LONG CHARACTER ---
  
  read(iuin,'(a)',end=999) line
!  print*,trim(line)
  flag = 0

!     --- FIND DATA SEPARATED BY SPACE OR COMMA ---

  idata = 0       ! temporal number of data
  ldata = 0
  spcom0 = .true.

  do i = 1, lnleng 
     spcom1 = (line(i:i)==space.or.line(i:i)==comma &
          &    .or.line(i:i)==tab) 
     if (spcom0) then
!           * the last character was space or comma
        if (.not.spcom1) then
!              * i is the begining of a data
           idata = idata + 1       ! temporal number of data
           if (idata>ndata) return
           ldata = 1               ! length of the datum
           fredat(idata)(ldata:ldata) = line(i:i)
        end if
     else 
!           * the last character was not space or comma
        if (.not.spcom1) then
           ldata = ldata + 1                               
           fredat(idata)(ldata:ldata) = line(i:i)
        end if
     end if
     spcom0 = spcom1
  end do

999 continue
!          +    +    +    +    +

  return
end subroutine rdfree
!-----------------------------------------------------------------------
