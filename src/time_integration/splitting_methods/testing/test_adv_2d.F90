!> @ingroup operator_splitting
!> @brief Unit test for operator splitting. Constant coefficient advection.
!> 
program test_adv_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_m_split_advection_2d
  use sll_m_advection_1d_BSL
  use sll_m_characteristics_1d_explicit_euler
  use sll_m_cubic_spline_interpolator_1d
  use sll_m_operator_splitting
  use sll_m_hdf5_io_serial
  implicit none
  class(split_advection_2d), pointer :: split
  class(sll_c_interpolator_1d), pointer :: interp1
  class(sll_c_interpolator_1d), pointer :: interp2
  class(sll_characteristics_1d_base), pointer :: charac1
  class(sll_characteristics_1d_base), pointer :: charac2
  sll_real64, dimension(:,:), pointer :: f
  sll_real64,dimension(:,:), pointer :: A1
  sll_real64,dimension(:,:), pointer :: A2
  sll_real64 :: x_1
  sll_real64 :: x_2
  sll_int32 :: i
  sll_int32 :: j
  sll_real64 :: dt
  sll_int32 :: file_id
  character(len=20) :: filename
  sll_real64 :: x1_min
  sll_real64 :: x1_max
  sll_real64 :: x2_min
  sll_real64 :: x2_max
  sll_int32 :: num_cells_x1
  sll_int32 :: num_cells_x2
  sll_real64, dimension(:), pointer :: x1_mesh
  sll_real64, dimension(:), pointer :: x2_mesh
  type(sll_cartesian_mesh_2d), pointer :: mesh_2d
  sll_real64 :: err
  sll_int32 :: ierr
  sll_real64 :: delta_x1
  sll_real64 :: delta_x2
  procedure(signature_process_outside_point_1d), pointer :: process_outside_point1
  procedure(signature_process_outside_point_1d), pointer :: process_outside_point2
  
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

  mesh_2d => new_cartesian_mesh_2d( &
    num_cells_x1, &
    num_cells_x2, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max)
  
  process_outside_point1 =>  process_outside_point_periodic
  process_outside_point2 =>  process_outside_point_periodic


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



  interp1 => new_cubic_spline_interpolator_1d( &
    num_cells_x1+1, &
    x1_min, &
    x1_max, &
    SLL_PERIODIC)


  charac1 => new_explicit_euler_1d_charac(&
    num_cells_x1+1, &
    SLL_PERIODIC)
  

  interp2 => new_cubic_spline_interpolator_1d( &
    num_cells_x2+1, &
    x2_min, &
    x2_max, &
    SLL_PERIODIC)


  charac2 => new_explicit_euler_1d_charac(&
    num_cells_x2+1, &
    SLL_PERIODIC)
  


  ! initialize time splitting method
  split => new_split_advection_2d( &
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
      SLL_ADVECTIVE, &
      SLL_STRANG_TVT) 
  
  
  
  
!  split => new_split_advection_2d( &
!    f, &
!    num_cells_x1+1, &
!    num_cells_x2+1, &
!    num_cells_x1+1, &
!    num_cells_x2+1, &
!    A1, &
!    A2, &
!    adv_x1, &
!    adv_x2, &
!    SLL_ADVECTIVE, &
!    SLL_STRANG_TVT)

  ! do some steps of lie_splitting
  dt = 0.5_f64
  call do_split_steps(split, dt, 4)

  ! save results
  filename = "data.h5"
  call sll_hdf5_file_create(filename, file_id, ierr)
  call sll_hdf5_write_array_2d(file_id, f, "data", ierr)
  call sll_hdf5_file_close(file_id, ierr)


  
    
  
  !err=maxval(abs(input-output))
  !to complete...
  print *,'#err=',err
  if(err<1.e-15_f64)then  
    print *,'#PASSED' 
  endif



  
end program test_adv_2d
