!> @ingroup sll_t_operator_splitting
!> @brief Unit test for operator splitting. Constant coefficient advection.
!> 
program test_const_coef_adv_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_const_coef_advection_2d, only: &
    sll_t_const_coef_advection_2d, &
    sll_f_new_const_coef_advection_2d

  use sll_m_cubic_spline_interpolator_1d, only: &
    sll_t_cubic_spline_interpolator_1d

  use sll_m_hdf5_io_serial, only: &
    sll_t_hdf5_ser_handle, &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close, &
    sll_o_hdf5_ser_write_array

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_operator_splitting, only: &
    sll_s_do_split_steps, &
    sll_p_strang_tvt

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define N1 50
#define N2 60
#define XMIN (-1.0_f64)
#define XMAX 1.0_f64
  class(sll_t_const_coef_advection_2d), pointer :: split
  sll_real64, dimension(:,:), pointer :: data
  sll_real64 :: x_1, x_2
  sll_int32 :: i, j
  sll_real64 :: dt
  sll_int32 :: ierr!, file_id
  character(len=20) :: filename

  type(sll_t_cubic_spline_interpolator_1d), target  :: interp_eta1
  type(sll_t_cubic_spline_interpolator_1d), target  :: interp_eta2
  class(sll_c_interpolator_1d), pointer :: interp_eta1_ptr
  class(sll_c_interpolator_1d), pointer :: interp_eta2_ptr

  type(sll_t_hdf5_ser_handle) :: hfile_id

  ! initialize interpolator
  call interp_eta1%initialize( N1, XMIN, XMAX, sll_p_periodic )
  call interp_eta2%initialize( N2, XMIN, XMAX, sll_p_periodic )
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
  split => sll_f_new_const_coef_advection_2d( data, N1, N2, 0.1_f64, 0.2_f64, &
       interp_eta1_ptr, interp_eta2_ptr, sll_p_strang_tvt)

  ! do some steps of lie_splitting
  dt = 0.5_f64
  call sll_s_do_split_steps(split, dt, 4)

  ! save results
  filename = "data.h5"
  call sll_s_hdf5_ser_file_create( filename, hfile_id, ierr )
  call sll_o_hdf5_ser_write_array( hfile_id, data, "data", ierr )
  call sll_s_hdf5_ser_file_close( hfile_id, ierr )
  
end program test_const_coef_adv_2d
