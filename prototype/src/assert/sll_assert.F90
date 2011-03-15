module sll_assert
  implicit none

contains
subroutine sll_assert(msg)
  character(*), intent(in) :: msg
  write (*,'(a)') msg
  stop 'ASSERTION FAILURE'
end subroutine sll_assert

end module sll_assert
