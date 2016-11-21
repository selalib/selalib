!> @ingroup sll_t_operator_splitting
!> @brief Unit test for operator splitting. Constant coefficient advection.
!> 
program test_adv_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d

  use sll_m_characteristics_1d_base, only: &
    sll_f_process_outside_point_periodic, &
    sll_i_signature_process_outside_point_1d, &
    sll_c_characteristics_1d_base

  use sll_m_characteristics_1d_explicit_euler, only: &
    sll_f_new_explicit_euler_1d_charac

  use sll_m_cubic_spline_interpolator_1d, only: &
    sll_f_new_cubic_spline_interpolator_1d

  use hdf5, only: hid_t
  use sll_m_hdf5_io_serial, only: &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close, &
    sll_o_hdf5_ser_write_array

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_operator_splitting, only: &
    sll_s_do_split_steps, &
    sll_p_strang_tvt

  use sll_m_split_advection_2d, only: &
    sll_f_new_split_advection_2d, &
    sll_p_advective, &
    sll_t_split_advection_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  class(sll_t_split_advection_2d), pointer :: split
  class(sll_c_interpolator_1d), pointer :: interp1
  class(sll_c_interpolator_1d), pointer :: interp2
  class(sll_c_characteristics_1d_base), pointer :: charac1
  class(sll_c_characteristics_1d_base), pointer :: charac2
  sll_real64, dimension(:,:), pointer :: f
  sll_real64,dimension(:,:), pointer :: A1
  sll_real64,dimension(:,:), pointer :: A2
  sll_real64 :: x_1
  sll_real64 :: x_2
  sll_int32 :: i
  sll_int32 :: j
  sll_real64 :: dt
  !sll_int32 :: file_id
  integer(hid_t) :: hfile_id
  character(len=20) :: filename
  sll_real64 :: x1_min
  sll_real64 :: x1_max
  sll_real64 :: x2_min
  sll_real64 :: x2_max
  sll_int32 :: num_cells_x1
  sll_int32 :: num_cells_x2
  sll_real64, dimension(:), pointer :: x1_mesh
  sll_real64, dimension(:), pointer :: x2_mesh
  type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
  sll_real64 :: err
  sll_int32 :: ierr
  sll_real64 :: delta_x1
  sll_real64 :: delta_x2
  procedure(sll_i_signature_process_outside_point_1d), pointer :: process_outside_point1
  procedure(sll_i_signature_process_outside_point_1d), pointer :: process_outside_point2
  
  x1_min = 0._f64
  x1_max = 1._f64
  x2_min = 0._f64
  x2_max = 1._f64
  num_cells_x1 = 32
  num_cells_x2 = 32
  dt = 0.1_f64
  
  delta_x1 = (x1_max-x1_min)/real(num_cells_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(num_cells_x2,f64)
  SLL_ALLOCATE(x1_mesh(num_cells_x1+1),ierr)
  SLL_ALLOCATE(x2_mesh(num_cells_x2+1),ierr)
  SLL_ALLOCATE(f(num_cells_x1+1,num_cells_x2+1),ierr)
  SLL_ALLOCATE(A1(num_cells_x1+1,num_cells_x2+1),ierr)
  SLL_ALLOCATE(A2(num_cells_x1+1,num_cells_x2+1),ierr)

  do i=1,num_cells_x1+1
    x1_mesh(i) = x1_min+real(i-1,f64)*delta_x1
  enddo

  do i=1,num_cells_x2+1
    x2_mesh(i) = x2_min+real(i-1,f64)*delta_x2
  enddo

  mesh_2d => sll_f_new_cartesian_mesh_2d( &
    num_cells_x1, &
    num_cells_x2, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max)
  
  process_outside_point1 =>  sll_f_process_outside_point_periodic
  process_outside_point2 =>  sll_f_process_outside_point_periodic


  do j=1,num_cells_x2+1
    x_2 = x2_mesh(j)
    do i=1,num_cells_x1+1
      x_1 = x1_mesh(i)
      f(i,j) = exp(-10._f64*(x_1*x_1+x_2*x_2))
    end do
  end do

  
  A1 = 1._f64
  A2 = 1._f64

  err=0._f64



  interp1 => sll_f_new_cubic_spline_interpolator_1d( &
    num_cells_x1+1, &
    x1_min, &
    x1_max, &
    sll_p_periodic)


  charac1 => sll_f_new_explicit_euler_1d_charac(&
    num_cells_x1+1, &
    sll_p_periodic)
  

  interp2 => sll_f_new_cubic_spline_interpolator_1d( &
    num_cells_x2+1, &
    x2_min, &
    x2_max, &
    sll_p_periodic)


  charac2 => sll_f_new_explicit_euler_1d_charac(&
    num_cells_x2+1, &
    sll_p_periodic)
  


  ! initialize time splitting method
  split => sll_f_new_split_advection_2d( &
      f, &
      A1, &
      A2, &
      interp1, &
      charac1, &
      process_outside_point1, &
      interp2, &
      charac2, &
      process_outside_point2, &
      mesh_2d, &
      sll_p_advective, &
      sll_p_strang_tvt) 
  
  
  
  
!  split => sll_f_new_split_advection_2d( &
!    f, &
!    num_cells_x1+1, &
!    num_cells_x2+1, &
!    num_cells_x1+1, &
!    num_cells_x2+1, &
!    A1, &
!    A2, &
!    adv_x1, &
!    adv_x2, &
!    sll_p_advective, &
!    sll_p_strang_tvt)

  ! do some steps of lie_splitting
  dt = 0.5_f64
  call sll_s_do_split_steps(split, dt, 4)

  ! save results
  filename = "data.h5"
  call sll_s_hdf5_ser_file_create( filename, hfile_id, ierr )
  call sll_o_hdf5_ser_write_array( hfile_id, f, "data", ierr )
  call sll_s_hdf5_ser_file_close( hfile_id, ierr )


  
    
  
  !err=maxval(abs(input-output))
  !to complete...
  print *,'#err=',err
  if(err<1.e-15_f64)then  
    print *,'#PASSED' 
  endif



  
end program test_adv_2d
