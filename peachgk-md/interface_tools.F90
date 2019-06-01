!********************************************
!*  interface_tools.f90 Ver.1.0  '09.11.16  *
!*      for peachgk_md.f                    *
!*            by G.Kikugawa                 *
!********************************************

!***** This module is interface module for routines of computational tools *****

module interface_tools

  interface
     
     subroutine rdfree(iuin, ndata, fredat)

       ! ARGUMENT:
       !        INPUT: 
       integer,intent(in)::     iuin          ! unit to read
       integer,intent(in)::     ndata         ! number of data to be read.

       !        OUTPUT:
       character(80),intent(inout):: fredat(ndata) ! data read as characters

     end subroutine rdfree

     subroutine rdfree_w(iuin, ndata, fredat, nword)

       ! ARGUMENT:
       !        INPUT: 
       integer,intent(in)::     iuin          ! unit to read
       integer,intent(in)::     ndata         ! number of data to be read.

       !        OUTPUT:
       character(80),intent(inout):: fredat(ndata) ! data read as characters
       !                   length of each datum is assumed =<80
       !                   each datum begins at the left end of data(i)
       !
       integer,intent(inout)::     nword         ! actual number of data

     end subroutine rdfree_w

     subroutine isort(idata, ndata)

       ! ARGUMENT:
       !     INPUT
       integer,intent(in)::   ndata         ! number of data 

       !     INPUT & OUTPUT 
       integer,intent(inout)::   idata(ndata)  ! data array to sort

     end subroutine isort

     subroutine iskip(idata, ndata)

       ! ARGUMENT:
       !     INPUT & OUTPUT 
       integer,intent(inout)::   ndata         ! number of data 
       integer,intent(inout)::   idata(ndata)  ! data array to sort

     end subroutine iskip

     subroutine excl_sp(in_word,out_len,out_word)

       ! ARGUMENT:
       !     INPUT
       integer,parameter:: in_len = 80

       character(in_len),intent(in):: in_word ! input character

       !     OUTPUT
       integer,intent(out):: out_len         ! output word length
       character(in_len),intent(out):: out_word ! output character

     end subroutine excl_sp

  end interface

end module interface_tools
