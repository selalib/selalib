module sll_working_precision
  implicit none
  intrinsic :: kind, selected_real_kind
  ! From the Fortran Standard (2.4.1.1): "The kind type parameter indicates 
  ! the decimal exponent range for the integer type (4.4.1), the decimal 
  ! precision and exponent range for the real and complex types 
  ! (4.4.2, 4.4.3), and the representation methods for the character and 
  ! logical types (4.4.4, 4.4.5)."

  ! The intent is that i32 will hold values up to 2**32-1
  integer, parameter :: sll_i32 = kind(0) 
  integer, parameter :: sll_i64 = kind(2_8**32)  ! 1.0d0 should be enough...
  integer, parameter :: sll_f32 = selected_real_kind(1,37)
  integer, parameter :: sll_f64 = selected_real_kind(1,99)

end module sll_working_precision
