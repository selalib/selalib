program test_time_splitting
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_field_2d.h"
  use sll_const_coef_advection_2d
  use sll_cubic_spline_interpolator_1d
#ifndef STDF95
  use sll_module_interpolators_1d_base
  use sll_time_splitting
#endif
  use sll_hdf5_io
  implicit none
#define N1 50
#define N2 60
#define XMIN (-1.0_f64)
#define XMAX 1.0_f64
  type(const_coef_advection_2d), target :: const_adv
#ifdef STDF95
  type(const_coef_advection_2d), pointer :: time_split
#else
  class(time_splitting), pointer :: time_split
#endif
  sll_real64, dimension(N1,N2) :: data
  sll_real64 :: x_1, x_2
  sll_int32 :: i, j
  sll_real64 :: dt
  sll_int32 :: ierr, file_id
  character(len=20) :: filename

  type(cubic_spline_1d_interpolator), target  :: interp_eta1
  type(cubic_spline_1d_interpolator), target  :: interp_eta2
#ifdef STDF95
  type(cubic_spline_1d_interpolator), pointer :: interp_eta1_ptr
  type(cubic_spline_1d_interpolator), pointer :: interp_eta2_ptr
#else
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr
#endif

  ! initialize interpolator
#ifdef STDF95
  call cubic_spline_1d_interpolator_initialize(interp_eta1, N1, XMIN, XMAX, SLL_PERIODIC )
  call cubic_spline_1d_interpolator_initialize(interp_eta2, N2, XMIN, XMAX, SLL_PERIODIC )
#else
  call interp_eta1%initialize( N1, XMIN, XMAX, SLL_PERIODIC )
  call interp_eta2%initialize( N2, XMIN, XMAX, SLL_PERIODIC )
#endif
  interp_eta1_ptr => interp_eta1
  interp_eta2_ptr => interp_eta2

  ! initialize data
  do j=1, N2
     x_2 = XMIN + (j-1)*(XMAX-XMIN)/(N2-1)
     do i=1, N1
        x_1 =  XMIN + (i-1)*(XMAX-XMIN)/(N1-1)
        data(i,j) = exp(-10*(x_1*x_1+x_2*x_2))
     end do
  end do

  ! initialize time splitting method
  call const_coef_advection_2d_initialize(const_adv, data, N1, N2, 0.1_f64, 0.2_f64, &
       interp_eta1_ptr, interp_eta2_ptr)
  time_split => const_adv

  ! do some steps of lie_splitting
  dt = 0.5
#ifdef STDF95
  call lie_splitting(time_split, dt, 4)
#else
  call time_split%lie_splitting(dt, 4)
#endif

#ifndef NOHDF5
  ! save results
  filename = "data.h5"
  call sll_hdf5_file_create(filename, file_id, ierr)
  call sll_hdf5_write_array_2d(file_id, data, "data", ierr)
  call sll_hdf5_file_close(file_id, ierr)
#endif
  
end program test_time_splitting
