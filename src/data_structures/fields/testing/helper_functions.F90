module helper_functions

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_constants, only: &
    sll_pi

  implicit none

  public :: &
    test_function_dirdir, &
    test_function_dirdir_der1, &
    test_function_dirdir_der2, &
    test_function_dirper, &
    test_function_dirper_der1, &
    test_function_dirper_der2, &
    test_function_perdir, &
    test_function_perdir_der1, &
    test_function_perdir_der2, &
    test_function_perper, &
    test_function_perper_der1, &
    test_function_perper_der2

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ! This module was made necessary because the assumed-shape dummy argument 
  ! params requires an explicit interface for the functions, so the usual
  ! method of putting these functions at the end of the unit test does not
  ! work.

contains 

! -------------> FUNCTION PERIODIC-PERIODIC
function test_function_perper( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = cos(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_perper

function test_function_perper_der1( eta1, eta2, params) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = -2*sll_pi*sin(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_perper_der1

function test_function_perper_der2( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = -2*sll_pi*cos(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_perper_der2


! -------------> FUNCTION PERIODIC-DIRICHLET
function test_function_perdir( eta1, eta2, params) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = cos(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_perdir

function test_function_perdir_der1( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = -2.0*sll_pi*sin(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_perdir_der1

function test_function_perdir_der2( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = 2.0*sll_pi*cos(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_perdir_der2


!-------------> FUNCTION DIRICHLET-PERIODIC
function test_function_dirper( eta1, eta2, params) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = sin(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_dirper

function test_function_dirper_der1( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = 2.0*sll_pi*cos(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_dirper_der1

function test_function_dirper_der2( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = -2.0*sll_pi*sin(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_dirper_der2


!-------------> FUNCTION DIRICHLET-DIRICHLET 
function test_function_dirdir( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = sin(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_dirdir

function test_function_dirdir_der1( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = 2.0*sll_pi*cos(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
end function test_function_dirdir_der1

function test_function_dirdir_der2( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), intent(in) :: params(:)
  real(8) :: res

#ifdef DEBUG
  real(8) :: dummy
  dummy = params(1)
#endif
  res = 2.0*sll_pi*sin(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
end function test_function_dirdir_der2


end module helper_functions
