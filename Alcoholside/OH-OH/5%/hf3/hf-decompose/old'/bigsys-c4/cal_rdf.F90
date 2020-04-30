program main
  implicit none

  integer::ios
  integer::ios1
  integer::ios2
  integer::ios3
  integer::nchunk
  real(8)::atom_number

  integer::realstep

  integer::chunk_number
  real(8),allocatable ::rdf(:)

  real(8),allocatable ::chunk_all(:)  !!!!!!!!!!!!!!!!!!!
  real(8),allocatable ::ave_rdf(:,:)
  
  integer::step=0
  integer::i,j,z,k
  integer::startstep



  character(len=30)::filename
  character(len=30)::new_filename
  character(1000)::long_line

  filename='rdf1.dat'
  new_filename='ave_rdf.dat'
  z=0
 ! open(unit=19,file='../density_step.dat',status='old',action='read')
 ! read(19,*) string0,startstep

  
  
  !open(unit=11,file='../data.txt',status='old',action='read',iostat=ios3)
  !read(11,*) filename,new_filename



  open(unit=10,file=filename,status='old',action='read')
  open(unit=102,file=new_filename,form='formatted',status='replace')


  ! skip headers
  do i=1,2
     read(10,*) 
  end do
  ! read the variable name
  read(10,'(a)') long_line
  do i=1,len_trim(long_line)
   if(long_line(i:i+1) == 'c_') then
      z=z+1
   end if
  end do

  print *,'z=',z

  write(102,*) trim(adjustl(long_line))

  !read chunk number


  read(10,*) realstep,nchunk
  write(*,*) 'read the chunk number:',nchunk

  allocate(ave_rdf(nchunk,z))
  allocate(rdf(z))
  allocate(chunk_all(nchunk))
  
  rdf=0.0d0
  ave_rdf=0.0d0


  do while(.true.)
     write(*,*) 'realstep=', realstep !,'startstep=',startstep

     do i=1,nchunk

        read(10,*,iostat=ios) chunk_all(i),(rdf(k),k=1,z)
        if (ios/=0) exit
        ave_rdf(i,1:z) = ave_rdf(i,1:z)+rdf(1:z)
     !   print*,'i=',i,'chunk=',rdf(1),rdf(2),ave_rdf(i,1)

     end do


     read(10,*,iostat=ios1) realstep,nchunk
     if (ios1/=0) exit
     step=step+1
  end do

  print*,step

  !write(*,*) 'density=',ave_rdf(1:nchunk)!,ave_ncount(1:nchunk)

  ave_rdf = ave_rdf /step

  do i=1,nchunk

     write(102,*) chunk_all(i),(ave_rdf(i,k),k=1,z)

  end do

  stop




end program main
