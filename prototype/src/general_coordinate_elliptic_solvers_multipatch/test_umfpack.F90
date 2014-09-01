program test_umfpack

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use mod_umfpack
  implicit none

  integer(umf_void) :: symbolic
  integer(umf_void) :: numeric

  sll_real64, dimension(umfpack_control) :: control
  sll_real64, dimension(umfpack_info) :: info
  sll_int32, parameter :: sys = 1

  integer, parameter :: n=4, nz=6
  sll_real64, dimension(nz) :: Ax
  sll_int32, dimension(nz) :: Ai
  sll_int32, dimension(n+1) :: Ap
  sll_real64, dimension(n) :: x, b

  call umf4def(control)  ! get the default configuration
  control(umfpack_prl) = real( 5 , umf_dp ) ! change verbosity
  call umf4pcon(control) ! update the umfpack configuration

  print*, umf_void
  ! Assemble matrix
  Ap = (/ 0,2,3,4,6/)
  Ai = (/0,1,2,3,0,1/)
  Ax = (/1.0,-3.5,-1.0,-2.0,1.0,7.0/)

  call umf4sym (n,n, Ap, Ai, Ax, symbolic, control, info)

  call umf4num (Ap, Ai, Ax, symbolic, numeric, control, info)

 call umf4pinf(control,info) ! print info 

  b=(/1.0,1.0,1.0,1.0/)
  call umf4sol (sys, x, b, numeric, control, info)
  print*, 'sol', x
end program test_umfpack
