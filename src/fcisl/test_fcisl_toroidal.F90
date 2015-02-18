program fcisl_toroidal_test
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_file_io.h"
  use sll_fcisl_toroidal_module
  use sll_constants

implicit none

  sll_real64 :: R0
  sll_real64 :: dt
  sll_int32 :: num_time_points
  sll_real64 :: F0
  sll_real64 :: smallr
  sll_real64 :: psipr
  sll_real64, dimension(:), allocatable :: time_points
  sll_real64, dimension(:), allocatable :: theta_euler
  sll_real64, dimension(:), allocatable :: phi_euler
  sll_real64, dimension(:), allocatable :: theta_analytic
  sll_real64, dimension(:), allocatable :: phi_analytic
  sll_real64, dimension(:), allocatable :: theta
  sll_real64, dimension(:), allocatable :: phi
  sll_real64 :: theta0
  sll_real64 :: phi0
  sll_int32 :: ierr
  sll_int32 :: file_id
  sll_real64 :: err
  
  R0 = 10._f64
  smallr = 2._f64
  psipr = 4._f64
  dt = 0.01
  num_time_points = 15000
  theta0 = 0._f64
  phi0 = 0._f64

  F0=-psipr*R0/smallr+1
  SLL_ALLOCATE(time_points(num_time_points),ierr)
  SLL_ALLOCATE(theta_euler(num_time_points),ierr)
  SLL_ALLOCATE(phi_euler(num_time_points),ierr)
  SLL_ALLOCATE(theta_analytic(num_time_points),ierr)
  SLL_ALLOCATE(phi_analytic(num_time_points),ierr)
  SLL_ALLOCATE(theta(num_time_points),ierr)
  SLL_ALLOCATE(phi(num_time_points),ierr)
  
  call compute_time_points(dt,num_time_points,time_points)
  
  call compute_euler_field( &
    R0, &
    time_points, &
    num_time_points, &
    psipr, &
    F0, &
    smallr, &
    theta_euler, &
    phi_euler, &
    theta0, &
    phi0)   
  call compute_modulo_vect(phi_euler,phi,num_time_points,2._f64*sll_pi)
  call compute_modulo_vect(theta_euler,theta,num_time_points,2._f64*sll_pi)
  call sll_ascii_file_create("euler_field.dat", file_id, ierr)    
  call sll_ascii_write_array_1d(file_id,time_points,ierr,num_time_points,phi,theta)  
  call sll_ascii_file_close(file_id,ierr)
  !for visualization: gnuplot -e 'plot "euler_field.dat" u 2:3 w p ps 1 lt 5'  


  call compute_analytic_field( &
    R0, &
    time_points, &
    num_time_points, &
    psipr, &
    F0, &
    smallr, &
    theta_analytic, &
    phi_analytic, &
    theta0, &
    phi0)   
  call compute_modulo_vect(phi_analytic,phi,num_time_points,2._f64*sll_pi)
  call compute_modulo_vect(theta_analytic,theta,num_time_points,2._f64*sll_pi)
  call sll_ascii_file_create("analytic_field.dat", file_id, ierr)    
  call sll_ascii_write_array_1d(file_id,time_points,ierr,num_time_points,phi,theta) 
  call sll_ascii_file_close(file_id,ierr)
  !for visualization: gnuplot -e 'plot "analytic_field.dat" u 2:3 w p ps 1 lt 5'  

  
  err=max(maxval(abs(phi_euler-phi_analytic)),maxval(abs(theta_euler-theta_analytic)))
  
  print *,'#err between euler and analytic=',err
  
  print*,'#PASSED'


end program