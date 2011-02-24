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
       write(*, '(a, a, i8)') descriptor, 'Triggered in FILE '//file_name// &
            ', in LINE: ', line_number
       stop 'ERROR: test_error_code(): exiting program'
    end if
  end subroutine test_error_code

#if 0
  ! Notice how the following macro could work for arrays of any dimension
  ! if the use of arbitrary array indexing were not necessary. Thus, if
  ! it were desired to use only indices that start from 1, this would be
  ! a way to ensure that all arrays that are declared with this allocator
  ! would have that property.

! I'll need a macro for every dimension... sucks!!
#define TEST_MACRO(arry, lo, hi) \
  sll_allocate_int_1d(arry, lo, hi,"Memory allocation failure in file ".conc.__FILE__.conc.", line: ".conc. XSTRNG(__LINE__))

#define DECLARE_ALLOCATION_FUNC(f_name, arry_type, dim)               \
  subroutine f_name(a_ptr, x1_lo, x1_hi, file_name, line_number);     \
    integer, parameter           :: long_int = selected_int_kind(12); \
    integer, intent(in)                     :: line_number;           \
    character(len=*), intent(in)            :: file_name;             \
    integer                                 :: err;                   \
    arry_type, dimension dim, pointer       :: a_ptr;                 \
    integer, intent(in)                     :: x1_lo;                 \
    integer(kind=long_int), intent(in)      :: x1_hi;                 \
    allocate(a_ptr(x1_lo:x1_hi),stat=err);                            \
    call test_error_code(err, file_name, line_number);                \
  end subroutine f_name

DECLARE_ALLOCATION_FUNC(sll_allocate_int_1d,     integer, (:))
DECLARE_ALLOCATION_FUNC(sll_allocate_real_1d,       real, (:))
DECLARE_ALLOCATION_FUNC(sll_allocate_complex_1d, complex, (:))

#endif

end module sll_memory
