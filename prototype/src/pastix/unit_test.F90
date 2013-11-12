program test_pastix

use sll_collective
use sll_pastix 
use sll_murge

implicit none

call sll_boot_collective()

call test_pastix_fortran()

!call test_pastix_murge()

call sll_halt_collective()

contains

subroutine test_pastix_fortran()

   call initialize_pastix()
   call delete_pastix()

end subroutine test_pastix_fortran

subroutine test_pastix_murge()

   call initialize_murge()
   call delete_murge()

end subroutine test_pastix_murge

end program test_pastix
