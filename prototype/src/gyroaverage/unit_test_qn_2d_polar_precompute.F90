!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

program unit_test_qn_2d_polar_precompute
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!#include "sll_field_2d.h"
use sll_module_qn_2d_polar_precompute

implicit none

  sll_int32 :: ierr
  sll_real64,dimension(:), allocatable :: Ti
  sll_real64 :: r_min
  sll_real64 :: r_max
  sll_int32 :: num_cells_r
  sll_int32 :: num_cells_theta
  sll_real64, dimension(:), allocatable :: mu_points
  sll_real64, dimension(:), allocatable :: mu_weights
  sll_int32 :: N_mu
  sll_int32 :: N_points
  sll_comp64, dimension(:,:,:), allocatable :: mat
  sll_comp64, dimension(:,:,:), allocatable :: mat_loc
  sll_real64, dimension(:), allocatable :: rho  
  sll_real64, dimension(:,:), allocatable :: points
  sll_real64, dimension(:), allocatable :: lambda
  sll_int32 :: i
  
  r_min = 0.2_f64
  r_max = 0.8_f64
  num_cells_r = 32
  num_cells_theta = 32
  N_mu = 2
  N_points = 8
  
  SLL_ALLOCATE(rho(num_cells_r+1),ierr)
  
  
  rho = 0.1_f64


  
  SLL_ALLOCATE(Ti(num_cells_r+1),ierr)
  SLL_ALLOCATE(mu_points(N_mu),ierr)
  SLL_ALLOCATE(mu_weights(N_mu),ierr)
  SLL_ALLOCATE(lambda(num_cells_r+1),ierr)
  
  lambda = 1._f64
  
  !SLL_ALLOCATE(mat(num_cells_r+1,num_cells_r+1,num_cells_theta),ierr)
  
  
  SLL_ALLOCATE(points(3,N_points),ierr)
  
  call compute_shape_circle(points,N_points)
      
  SLL_ALLOCATE(mat_loc(0:num_cells_theta-1,0:num_cells_r,0:num_cells_r), ierr)
  SLL_ALLOCATE(mat(0:num_cells_theta-1,0:num_cells_r,0:num_cells_r), ierr)
  
  mat = (0._f64,0._f64)
  
  
  do i=1,N_mu
    rho(1:num_cells_r+1) = sqrt(2._f64*mu_points(i))*Ti(1:num_cells_r+1)
    call compute_qns_matrix_polar_splines( &
      r_min, &
      r_max, &
      num_cells_r, &
      num_cells_theta, &
      rho, &
      points, &
      N_points, &
      mat_loc)  
    mat = mat+mat_loc*mu_weights(i)
  enddo

  call compute_qns_inverse_polar_splines( &
    mat, &
    lambda, &
    num_cells_r, &
    num_cells_theta)
  
  
!  call precompute_matrix( &
!    Ti, &
!    r_min, &
!    r_max, &
!    num_cells_r, &
!    num_cells_theta, &
!    mu_points, &
!    mu_weights, &
!    N_mu, &
!    N_points, &
!    mat)
  
  print *,'#PASSED'
  

end program