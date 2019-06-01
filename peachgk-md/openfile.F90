!*****************************
!*  openfile.f90 Ver.3.4     *
!*      for peachgk_md.f     *
!*            by G.Kikugawa  *
!*****************************
! Time-stamp: <>

subroutine openfile(iucor,iucorname,iutop,iutopname, &
     &              iuwtop,iuwtopname, &
     &              iuparavdw,iuparavdwname,iuparabond,iuparabondname, &
     &              iuparaconst,iuparaconstname, &
     &              iuparacstmnb,iuparacstmnbname,ifcstmnb, &
     &              iuaddtop,iuaddtopname,ifrdaddtop, &
     &              iustrmvel,iustrmvelname,ifstrmvel, &
     &              iuposres,iuposresname,ifposres, &
     &              iostarec,iostarecname,ifstarec, &
     &              npolytyp, &
     &              ousum,ousumname, &
     &              ouene,ouenename,oupos,ouposname,ouvel,ouvelname, &
     &              oufor,ouforname, &
     &              outhe,outhename,oubar,oubarname, &
     &              oupre,ouprename,outhc,outhcname,ifoutthc, &
     &              oupdb,oupdbname,ifoutpdb, &
     &              ouhtf,ouhtfname, &
     &              oumtf,oumtfname, &
     &              ouumb,ouumbname,ifpotbias, &
     &              ifoutene,ifoutpos,ifoutvel,ifoutfor,ifoutthe, &
     &              ifoutbar,ifoutpre)

#if defined(MPI)
  use mpi_global    ! include mpi parameter
#endif

  implicit none

! ARGUMENT:
!     INPUT
  character(80),intent(in):: iucorname(:)     ! input file name
  character(80),intent(in):: iutopname(:)     !       "
  character(80),intent(in):: iuwtopname       !       "
  character(80),intent(in):: iuparavdwname    !       "
  character(80),intent(in):: iuparabondname   !       "
  character(80),intent(in):: iuparaconstname  !       "
  character(80),intent(in):: iuparacstmnbname !       "
  character(80),intent(in):: iuaddtopname     !       "
  character(80),intent(in):: iustrmvelname    !       "
  character(80),intent(in):: iuposresname     !       "

  character(80),intent(in):: iostarecname ! state record file name
  logical,intent(in):: ifstarec        ! read old state

  character(80),intent(in):: ousumname ! output file name
  character(80),intent(in):: ouenename ! output file name
  character(80),intent(in):: ouposname ! output file name
  character(80),intent(in):: ouvelname ! output file name
  character(80),intent(in):: ouforname ! output file name
  character(80),intent(in):: outhename ! output file name
  character(80),intent(in):: oubarname ! output file name
  character(80),intent(in):: ouprename ! output file name
  character(80),intent(in):: outhcname ! output file name
  character(80),intent(in):: oupdbname ! output file name
  character(80),intent(in):: ouhtfname ! output file name
  character(80),intent(in):: oumtfname ! output file name
  character(80),intent(in):: ouumbname ! output file name

  integer,intent(in):: npolytyp       ! number of poly type

  logical,intent(in):: ifcstmnb       ! flag if using custom NB interaction

  logical,intent(in):: ifrdaddtop     ! input additional topology information

  logical,intent(in):: ifoutthc       ! flag for outputting thermal control file

  logical,intent(in):: ifoutpdb       ! flag for outputting PDB format file

  logical,intent(in):: ifposres       ! position restraint flag

  logical,intent(in):: ifpotbias      ! bias potential flag
  logical,intent(in):: ifstrmvel      ! flag to input and use streaming velocity

  logical,intent(in):: ifoutene       ! if ouput energy file
  logical,intent(in):: ifoutpos       ! if ouput position file
  logical,intent(in):: ifoutvel       ! if ouput velocity file
  logical,intent(in):: ifoutfor       ! if ouput force file
  logical,intent(in):: ifoutthe       ! if ouput NVT file
  logical,intent(in):: ifoutbar       ! if ouput NPT file
  logical,intent(in):: ifoutpre       ! if ouput pressure file

!     OUTPUT
  integer,intent(out):: iucor(:)      ! input poly coordinate file unit
  integer,intent(out):: iutop(:)      ! input poly topology file unit
  integer,intent(out):: iuwtop        ! input poly topology file unit (water)
  integer,intent(out):: iuparavdw     ! input vdw parameter file unit
  integer,intent(out):: iuparabond    ! input bond parameter file unit
  integer,intent(out):: iuparaconst   ! input const parameter file unit
  integer,intent(out):: iuparacstmnb  ! input custom NB parameter file unit
  integer,intent(out):: iuaddtop      ! input additional topology file unit
  integer,intent(out):: iustrmvel     ! input streaming velocity file unit
  integer,intent(out):: iuposres      ! input position restraint ref. file unit

  integer,intent(out):: iostarec      ! state record file unit

  integer,intent(out):: ousum         ! output parameter summarization file unit
  integer,intent(out):: ouene         ! output unit for output energy data
  integer,intent(out):: oupos         ! output unit for output position data
  integer,intent(out):: ouvel         ! output unit for output velocity data
  integer,intent(out):: oufor         ! output unit for output force data
  integer,intent(out):: outhe         ! output unit for output thermostat data
  integer,intent(out):: oubar         ! output unit for output barostat data
  integer,intent(out):: oupre         ! output unit for output pressure data
  integer,intent(out):: outhc         ! output unit for outthc thermal control data
  integer,intent(out):: oupdb         ! output unit for outpdb PDB data
  integer,intent(out):: ouhtf         ! output unit for outhtf heat flux data
  integer,intent(out):: oumtf         ! output unit for outmtf momentum flux data
  integer,intent(out):: ouumb         ! output unit for outumb bias potential data

! LOCAL:
  integer:: ios             ! i/o status
  integer:: i               ! loop index

!     +     +     +     +     +     +     +

!---- define i/o unit

  if (npolytyp > 20) then
     write(6,*) 'Failure in openfile: too much poly type!'
     stop
  end if

!**** unit number 51-70 is reserved for iucor ****
  do i = 1, npolytyp
     iucor(i) = 50 + i
  end do
!**** unit number 71-90 is reserved for iutop ****
  do i= 1, npolytyp
     iutop(i) = 70 + i
  end do

  iuwtop = 18
  iuparavdw = 12
  iuparabond = 13
  iuparaconst = 22
  iuparacstmnb = 31
  iuaddtop = 30
  iustrmvel = 33
  iuposres = 27

  iostarec = 20

  ousum = 14
  ouene = 15
  oupos = 16
  ouvel = 17
  oufor = 34
  outhe = 23
  oubar = 24
  oupre = 21
  outhc = 32
  oupdb = 26
  ouhtf = 25
  oumtf = 29
  ouumb = 28

!-------- Opening input file --------

  do i = 1, npolytyp
     open(iucor(i),file=iucorname(i),status='old',iostat=ios)
     if (ios /= 0) then
        write(6,*) 'Failure in opening file: ',iucorname(i)
        stop
     end if
  end do

  do i = 1, npolytyp
     open(iutop(i),file=iutopname(i),status='old',iostat=ios)
     if (ios /= 0) then
        write(6,*) 'Failure in opening file: ',iutopname(i)
        stop
     end if
  end do

  open(iuwtop,file=iuwtopname,status='old',iostat=ios)
  if (ios /= 0) then
     write(6,*) 'Failure in opening file: ',iuwtopname
     stop
  end if

  open(iuparavdw,file=iuparavdwname,status='old',iostat=ios)
  if (ios /= 0) then
     write(6,*) 'Failure in opening file: ',iuparavdwname
     stop
  end if

  open(iuparabond,file=iuparabondname,status='old',iostat=ios)
  if (ios /= 0) then
     write(6,*) 'Failure in opening file: ',iuparabondname
     stop
  end if

  open(iuparaconst,file=iuparaconstname,status='old',iostat=ios)
  if (ios /= 0) then
     write(6,*) 'Failure in opening file: ',iuparaconstname
     stop
  end if

  if (ifcstmnb) then
     open(iuparacstmnb,file=iuparacstmnbname,status='old',iostat=ios)
     if (ios /= 0) then
        write(6,*) 'Failure in opening file: ',iuparacstmnbname
        stop
     end if
  end if

  if (ifrdaddtop) then
     open(iuaddtop,file=iuaddtopname,status='old',iostat=ios)
     if (ios /= 0) then
        write(6,*) 'Failure in opening file: ',iuaddtopname
        stop
     end if
  end if

  if (ifstrmvel) then
     open(iustrmvel,file=iustrmvelname,status='old',iostat=ios)
     if (ios /= 0) then
        write(6,*) 'Failure in opening file: ',iustrmvelname
        stop
     end if
  end if

  if (ifposres) then
     open(iuposres,file=iuposresname,status='old',iostat=ios)
     if (ios /= 0) then
        write(6,*) 'Failure in opening file: ',iuposresname
        stop
     end if
  end if

!-------- Opening state file --------

  open(iostarec,file=iostarecname,status='unknown',iostat=ios)
  if (ios /= 0) then
     write(6,*) 'Failure in opening file: ',iostarecname
     stop
  end if

!-------- Opening output file --------

#if defined(MPI)
  if (irank == 0) then
#endif

     open(ousum,file=ousumname,status='new',iostat=ios)
     if (ios /= 0) then
        write(6,*) 'Failure in opening file: ',ousumname
        stop
     end if

     if (ifoutene) then
        open(ouene,file=ouenename,status='new',iostat=ios)
        if (ios /= 0) then
           write(6,*) 'Failure in opening file: ',ouenename
           stop
        end if
     end if

     if (ifoutpos) then
        open(oupos,file=ouposname,form='unformatted', &
             & status='new',iostat=ios)
        if (ios /= 0) then
           write(6,*) 'Failure in opening file: ',ouposname
           stop
        end if
     end if

     if (ifoutvel) then
        open(ouvel,file=ouvelname,form='unformatted', &
             & status='new',iostat=ios)
        if (ios /= 0) then
           write(6,*) 'Failure in opening file: ',ouvelname
           stop
        end if
     end if

     if (ifoutfor) then
        open(oufor,file=ouforname,form='unformatted', &
             & status='new',iostat=ios)
        if (ios /= 0) then
           write(6,*) 'Failure in opening file: ',ouforname
           stop
        end if
     end if

     if (ifoutthe) then
        open(outhe,file=outhename,status='new',iostat=ios)
        if (ios /= 0) then
           write(6,*) 'Failure in opening file: ',outhename
           stop
        end if
     end if

     if (ifoutbar) then
        open(oubar,file=oubarname,status='new',iostat=ios)
        if (ios /= 0) then
           write(6,*) 'Failure in opening file: ',oubarname
           stop
        end if
     end if

     if (ifoutpre) then
        open(oupre,file=ouprename,status='new',iostat=ios)
        if (ios /= 0) then
           write(6,*) 'Failure in opening file: ',ouprename
           stop
        end if
     end if

     if (ifoutthc) then
        open(outhc,file=outhcname,status='new',iostat=ios)
        if (ios /= 0) then
           write(6,*) 'Failure in opening file: ',outhcname
           stop
        end if
     end if

     if (ifoutpdb) then
        open(oupdb,file=oupdbname,status='new',iostat=ios)
        if (ios /= 0) then
           write(6,*) 'Failure in opening file: ',oupdbname
           stop
        end if
     end if

#if defined(HF)
     open(ouhtf,file=ouhtfname,status='new',iostat=ios)
     if (ios /= 0) then
        write(6,*) 'Failure in opening file: ',ouhtfname
        stop
     end if

     open(oumtf,file=oumtfname,status='new',iostat=ios)
     if (ios /= 0) then
        write(6,*) 'Failure in opening file: ',oumtfname
        stop
     end if

#endif

     if (ifpotbias) then
        open(ouumb,file=ouumbname,status='new',iostat=ios)
        if (ios /= 0) then
           write(6,*) 'Failure in opening file: ',ouumbname
           stop
        end if
     end if

#if defined(MPI)
  end if
#endif

!     +     +     +     +     +     +     +

  return
end subroutine openfile
