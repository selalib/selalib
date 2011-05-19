module sll_memory
  implicit none
  intrinsic :: merge
! Useful for 'alternative implementation'. Ignore for now.
#if 0  
 interface sll_allocate
     module procedure sll_allocate_int_1d, sll_allocate_real_1d, &
          sll_allocate_complex_1d
  end interface
#endif

contains 

  subroutine test_error_code(err_code, descriptor, file_name, line_number)
    integer                      :: err_code
    character(len=*), intent(in) :: descriptor
    character(len=*), intent(in) :: file_name
    integer, intent(in)          :: line_number
    if (err_code .ne. 0) then
       write(*, '(a, a, i8)') descriptor, ' Triggered in FILE '//file_name// &
            ', in LINE: ', line_number
       stop 'ERROR: test_error_code(): exiting program'
    end if
  end subroutine test_error_code


end module sll_memory
