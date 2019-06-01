!**************************************
!*  excl_sp.f Ver.1.0 '03.05.19       *
!*            by G.Kikugawa           *
!**************************************
subroutine excl_sp(in_word,out_len,out_word)

  implicit none

! This subroutine is for excluding space from characters.
!
! usage
! input word length = 80
! in_word: input character
! out_len: output word length (spaces are excluded)
! out_word: output character (spaces are excluded)
!
! ARGUMENT:
!     INPUT
  integer,parameter:: in_len = 80

  character(in_len),intent(in):: in_word ! input character

!     OUTPUT
  integer,intent(out):: out_len         ! output word length
  character(in_len),intent(out):: out_word ! output character

! LOCAL:
  integer:: i
  logical:: ifsp

!     +     +     +     +     +     +     +     +

  out_len = 0
  do i=1,80
     ifsp = (in_word(i:i)==' ')
     if (.not.ifsp) then
        out_len = out_len + 1
        out_word(out_len:out_len) = in_word(i:i)
     end if
  end do

  return
end subroutine excl_sp
