!> @ingroup operator_splitting
!> @brief Unit test for operator splitting. Constant coefficient advection.
!> 
program test_const_coef_adv_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_m_const_coef_advection_2d
  use sll_m_cubic_spline_interpolator_1d
  use sll_m_interpolators_1d_base
  use sll_m_operator_splitting
  use sll_m_hdf5_io_serial
  implicit none
#define N1 50
#define N2 60
#define XMIN (-1.0_f64)
#define XMAX 1.0_f64
  class(const_coef_advection_2d), pointer :: split
  sll_real64, dimension(:,:), pointer :: data
  sll_real64 :: x_1, x_2
  sll_int32 :: i, j
  sll_real64 :: dt
  sll_int32 :: ierr, file_id
  character(len=20) :: filename

  type(sll_cubic_spline_interpolator_1d), target  :: interp_eta1
  type(sll_cubic_spline_interpolator_1d), target  :: interp_eta2
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr

  ! initialize interpolator
  call interp_eta1%initialize( N1, XMIN, XMAX, SLL_PERIODIC )
  call interp_eta2%initialize( N2, XMIN, XMAX, SLL_PERIODIC )
  interp_eta1_ptr => interp_eta1
  interp_eta2_ptr => interp_eta2

  ! initialize data
  SLL_ALLOCATE(data(N1,N2), ierr)
  do j=1, N2
     x_2 = XMIN + (j-1)*(XMAX-XMIN)/(N2-1)
     do i=1, N1
        x_1 =  XMIN + (i-1)*(XMAX-XMIN)/(N1-1)
        data(i,j) = exp(-10*(x_1*x_1+x_2*x_2))
     end do
  end do

  ! initialize time splitting method
  split => new_const_coef_advection_2d( data, N1, N2, 0.1_f64, 0.2_f64, &
       interp_eta1_ptr, interp_eta2_ptr, SLL_STRANG_TVT)

  ! do some steps of lie_splitting
  dt = 0.5
  call do_split_steps(split, dt, 4)

  ! save results
  filename = "data.h5"
  call sll_hdf5_file_create(filename, file_id, ierr)
  call sll_hdf5_write_array_2d(file_id, data, "data", ierr)
  call sll_hdf5_file_close(file_id, ierr)
  
end program test_const_coef_adv_2d
