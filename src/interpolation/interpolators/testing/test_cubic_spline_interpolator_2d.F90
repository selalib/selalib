program unit_test_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_hermite, &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_cubic_spline_interpolator_2d, only: &
    sll_f_new_cubic_spline_interpolator_2d, &
    sll_t_cubic_spline_interpolator_2d

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define NPTS1 65 
#define NPTS2 65 

  !type(sll_t_cubic_spline_interpolator_2d) :: cs2d
  class(sll_c_interpolator_2d), pointer :: cs2d
  sll_real64, dimension(:,:), allocatable    :: x1
  sll_real64, dimension(:), allocatable      :: x1_eta1_min
  sll_real64, dimension(:), allocatable      :: x1_eta1_max
  sll_real64, dimension(2) :: params ! for the transformation
  sll_int32  :: i, j
  sll_real64 :: eta1, eta2, h1, h2, acc, acc1, acc2, node_val, ref, deriv1_val 
  sll_real64 :: deriv2_val

#define RMIN 0.1_f64
#define RMAX 1.0_f64

  params(:) = (/RMIN,RMAX/)

  print *,  'filling out discrete arrays for x1 '
  h1 = 1.0_f64/real(NPTS1-1,f64)
  h2 = 1.0_f64/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  allocate(x1(NPTS1,NPTS2))
  allocate(x1_eta1_min(NPTS2))
  allocate(x1_eta1_max(NPTS2))

  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1          = real(i,f64)*h1
        eta2          = real(j,f64)*h2
        x1(i+1,j+1)   = x1_polar_f(eta1,eta2,params) 
     end do
  end do
  print *, 'eta1, eta2 = ', real(NPTS1-1,f64)*h1, real(NPTS2-1,f64)*h2
  print *, 'x1_polar_f(eta1=1, eta2=1) = ', x1_polar_f(1.0_f64,1.0_f64,params)
  ! Fill out the transformation's slopes at the borders
  do j=0,NPTS2-1
     eta1           = 0.0_f64
     eta2           = real(j,f64)*h2
     x1_eta1_min(j+1) = deriv_x1_polar_f_eta1(eta1,eta2,params)
     eta1           = 1.0_f64
     x1_eta1_max(j+1) = deriv_x1_polar_f_eta1(eta1,eta2,params)
  end do

  ! Test the 2D transformation:
  !
  ! X1 = (r1 + (r2-r1)*eta1)*cos(2*pi*eta2)
  
  cs2d =>sll_f_new_cubic_spline_interpolator_2d(&
       NPTS1, &
       NPTS2, &
       0.0_f64, &
       1.0_f64, &
       0.0_f64, &
       1.0_f64, &
       sll_p_hermite, &
       sll_p_periodic, &
       eta1_min_slopes=x1_eta1_min, &
       eta1_max_slopes=x1_eta1_max )
  
  
!  call cs2d%initialize( &
!       NPTS1, &
!       NPTS2, &
!       0.0_f64, &
!       1.0_f64, &
!       0.0_f64, &
!       1.0_f64, &
!       sll_p_hermite, &
!       sll_p_periodic, &
!       eta1_min_slopes=x1_eta1_min, &
!       eta1_max_slopes=x1_eta1_max )

  call cs2d%compute_interpolants(x1)
  print *, 'Compare the values of the transformation at the nodes: '
  acc  = 0.0_f64
  acc1 = 0.0_f64
  acc2 = 0.0_f64
  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1       = real(i,f64)*h1
        eta2       = real(j,f64)*h2
        node_val   = cs2d%interpolate_from_interpolant_value(eta1,eta2)
        ref        = x1_polar_f(eta1,eta2,params)
        acc        = acc + abs(node_val-ref)
        deriv1_val = cs2d%interpolate_from_interpolant_derivative_eta1(eta1,eta2)
        ref        = deriv_x1_polar_f_eta1(eta1,eta2,params)
        acc1       = acc1 + abs(deriv1_val-ref)
        deriv2_val = cs2d%interpolate_from_interpolant_derivative_eta2(eta1,eta2)
        ref        = deriv_x1_polar_f_eta2(eta1,eta2,params)
        acc2       = acc2 + abs(deriv2_val-ref)
     end do
  end do
  print *, 'Average error in nodes, x1 transformation = ', acc/real(NPTS1*NPTS2,f64)
  print *, 'Average error, x1 deriv eta1 = ', acc1/real(NPTS1*NPTS2,f64)
  print *, 'Average error, x1 deriv eta2 = ', acc2/real(NPTS1*NPTS2,f64)
  print *, 'PASSED'

  call test_interpolator_2d()

contains

subroutine test_interpolator_2d()
  class(sll_c_interpolator_2d),    pointer   :: interp
  type(sll_t_cubic_spline_interpolator_2d), target    :: spline
  sll_real64, dimension(NPTS1,NPTS2) :: xx1
  sll_real64, dimension(NPTS1,NPTS2) :: xx2
  sll_real64, dimension(NPTS1,NPTS2) :: data_in
  sll_real64, dimension(NPTS1,NPTS2) :: data_out

  call spline%initialize(NPTS1,NPTS2, &
                         0.0_f64,2.0*sll_p_pi,0.0_f64,2.*sll_p_pi, &
                         sll_p_periodic, sll_p_periodic )
  interp =>  spline
  do j = 1, NPTS2
  do i = 1, NPTS1
     xx1(i,j) = 2.*sll_p_pi*real(i-1,f64)/real(NPTS1-1,f64)
     xx2(i,j) = 2.*sll_p_pi*real(j-1,f64)/real(NPTS2-1,f64)
  end do
  end do
  data_in = cos(xx1)*sin(xx2)

  do j = 1, NPTS2
  do i = 1, NPTS1
     xx1(i,j) = 2.*sll_p_pi*real(i-1,f64)/real(NPTS1,f64)
     xx2(i,j) = 2.*sll_p_pi*real(j-1,f64)/real(NPTS2,f64)
  end do
  end do
  call interp%interpolate_array(NPTS1, NPTS2, data_in, xx1, xx2, data_out)

  print*, " error = ", maxval(abs(data_out-cos(xx1)*sin(xx2)))
end subroutine test_interpolator_2d

  function deriv_x1_polar_f_eta1( eta1, eta2, params )
    sll_real64 :: deriv_x1_polar_f_eta1
    sll_real64, intent(in) :: eta1, eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    r1 = params(1)
    r2 = params(2)
    deriv_x1_polar_f_eta1 = eta1
    deriv_x1_polar_f_eta1 = (r2-r1)*cos(2.0_f64*sll_p_pi*eta2)
  end function deriv_x1_polar_f_eta1

  function deriv_x1_polar_f_eta2( eta1, eta2, params )
    sll_real64 :: deriv_x1_polar_f_eta2
    sll_real64, intent(in) :: eta1, eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: k
    sll_real64 :: r1
    sll_real64 :: r2
    deriv_x1_polar_f_eta2  = eta1+eta2

    r1 = params(1)
    r2 = params(2)
    k = 2.0_f64*sll_p_pi
    deriv_x1_polar_f_eta2 = -(r1+(r2-r1)*eta1)*sin(k*eta2)*k
  end function deriv_x1_polar_f_eta2

  function x1_polar_f( eta1, eta2, params )
    sll_real64 :: x1_polar_f
    sll_real64, intent(in) :: eta1, eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    r1 = params(1)
    r2 = params(2)
    x1_polar_f = (r1 + (r2-r1)*eta1)*cos(2.0_f64*sll_p_pi*eta2)
  end function x1_polar_f

end program unit_test_2d
