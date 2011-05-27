module sll_misc_utils
  implicit none
  intrinsic :: selected_int_kind

  ! Tentative implementation of a standard-compliant way to get the
  ! memory footprint of a variable. This is our yardstick...
  !
  ! selected_int_kind(r) returns the default integer scalar that is the kind 
  ! type parameter value for an integer data type able to represent all 
  ! integer values n in the range -10^(-r) < n < 10^r, where r is a scalar 
  ! integer. If more than one is available, a kind with least decimal exponent 
  ! range is chosen (and least kind value if several have least decimal 
  ! exponent range). If no corresponding kind is availalble, the result is -1. 
  ! (Metcalf & Reid. F90/95 2nd ed. p. 176).
  !
  ! We are, maybe dangerously, relying on the common practice of many compilers
  ! of using the kind values to indicate the number of bytes of storage
  ! occupied by a value. But this is not mandated by the standard. For 
  ! instance a compiler that only has available 4-byte integers may still
  ! support kind values of 1, 2 and 4 to 'ease portability' from other 
  ! platforms. The size of k1 will become our basic 'yardstick' to measure
  ! the size of a memory footprint of a variable. When we ask for the size
  ! of 'var', the answer will be given in terms of how many 'yardsticks'
  ! are needed to represent 'var'. The implications of this assumption
  ! need to be checked further.

  integer, parameter :: k1 = selected_int_kind(0)

end module sll_misc_utils
