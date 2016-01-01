program test_crash

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  character( len=* ), parameter :: my_name = "test_crash"

  SLL_WARNING( my_name, "Let's test a warning." )
  SLL_ERROR  ( my_name, "Let's test an error."  )

end program test_crash
