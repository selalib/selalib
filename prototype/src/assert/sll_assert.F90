module sll_assertion
  implicit none

contains

subroutine sll_assert(msg, file, line)
  intrinsic :: trim
  character(len=*), intent(in) :: msg
  character(len=*), intent(in) :: file
  integer, intent(in)          :: line
  character(len=8)            :: local_line
  write(local_line, '(i8)') line ! hoping that I could trim this later, but no..
  write (*,'(a, a, a, a, a)') msg, ': Assertion error triggered in file ', file, ' in line ', trim(local_line)
  stop 'ASSERTION FAILURE'
end subroutine sll_assert

end module sll_assertion
