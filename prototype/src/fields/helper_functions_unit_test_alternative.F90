module helper_functions
  implicit none
  
  ! This module was made necessary because the assumed-shape dummy argument 
  ! params requires an explicit interface for the functions, so the usual
  ! method of putting these functions at the end of the unit test does not
  ! work.

contains 

 ! ------------->FUNCTION PERIODIC- PERIODIC
function test_function_perper( eta1, eta2, params ) result(res)
   use sll_constants
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  intrinsic :: cos
  real(8), dimension(:), intent(in) :: params
  res = cos(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_perper

function test_function_perper_der1( eta1, eta2, params) result(res)
  use sll_constants
  intrinsic :: cos,sin
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  res = -2*sll_pi*sin(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_perper_der1

function test_function_perper_der2( eta1, eta2, params ) result(res)
  use sll_constants
  intrinsic :: sin
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  res = -2*sll_pi*cos(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_perper_der2


 ! ------------->FUNCTION PERIODIC- DIRICHLET
function test_function_perdir( eta1, eta2, params) result(res)
  use sll_constants
  intrinsic :: cos,sin
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  res = cos(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_perdir

function test_function_perdir_der1( eta1, eta2, params ) result(res)
  use sll_constants
  intrinsic :: sin
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  res = -2.0*sll_pi*sin(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_perdir_der1

function test_function_perdir_der2( eta1, eta2, params ) result(res)
  use sll_constants
  intrinsic :: cos,sin
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  res = 2.0*sll_pi*cos(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_perdir_der2

! ------------->FUNCTION DIRICHLET- PERIODIC
function test_function_dirper( eta1, eta2, params) result(res)
  use sll_constants
  intrinsic :: cos,sin
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  res = sin(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_dirper

function test_function_dirper_der1( eta1, eta2, params ) result(res)
  use sll_constants
  intrinsic :: cos
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  res = 2.0*sll_pi*cos(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_dirper_der1

function test_function_dirper_der2( eta1, eta2, params ) result(res)
  use sll_constants
  intrinsic :: sin
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  res = -2.0*sll_pi*sin(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_dirper_der2


! ------------->FUNCTION DIRICHLET-DIRICHLET 


function test_function_dirdir( eta1, eta2, params ) result(res)
  use sll_constants
  intrinsic :: sin
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  res = sin(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_dirdir

function test_function_dirdir_der1( eta1, eta2, params ) result(res)
  use sll_constants
  intrinsic :: cos,sin
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  res = 2.0*sll_pi*cos(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_dirdir_der1

function test_function_dirdir_der2( eta1, eta2, params ) result(res)
  use sll_constants
  intrinsic :: cos,sin
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  res = 2.0*sll_pi*sin(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_dirdir_der2


end module helper_functions
