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
use sll_timer
use sll_ascii_io
use gauss_lobatto_integration
use gauss_legendre_integration
use sll_gnuplot

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
  sll_real64 :: mu_max
  character(len=256) :: filename
  sll_int32             :: IO_stat
  sll_int32, parameter  :: input_file = 99
  type(sll_time_mark)  :: t0
  sll_real64 :: time
  sll_int32 :: j
  sll_int32 :: m
  sll_int32 :: mat_id
  sll_int32 :: nml_id
  sll_real64 :: lambda_coeff
  character(len=256) :: quadrature_case
  sll_int32 :: quadrature_points_per_cell
  sll_int32 :: num_cells
  sll_int32 :: s
  sll_int32 :: mode_r
  sll_int32 :: mode_theta
  sll_real64 :: eta_min(2)
  sll_real64 :: eta_max(2)
  sll_int32 :: mode(2)
  sll_int32 :: Nc(2)
  sll_real64,dimension(:,:),allocatable :: phi,phi_restr,phi_init,phi_qn
  sll_real64,dimension(:,:),allocatable :: phi_quotient
  sll_real64 :: gamma0
  sll_real64 :: gamma0_num
    
  namelist /params/ &
    r_min, &
    r_max, &
    num_cells_r, &
    num_cells_theta, &
    N_mu, &
    N_points, &
    mu_max, &
    lambda_coeff, &
    quadrature_case, &
    quadrature_points_per_cell, &
    mode_r, &
    mode_theta

  call sll_set_time_mark(t0)


  !default parameters
  r_min = 0.2_f64
  r_max = 0.8_f64
  num_cells_r = 128
  num_cells_theta = 256
  N_mu = 2
  N_points = 32
  mu_max = 3._f64
  lambda_coeff = 1._f64
  quadrature_case = "SLL_RECTANGLE"
  quadrature_points_per_cell = 10
  mode_r = 1
  mode_theta = 1
  
  call get_command_argument(1, filename)
  if (len_trim(filename) == 0)then
    print *,'#initialization with default parameters'
  else  
    open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
      if( IO_stat /= 0 ) then
        print *, '#program unit_test_qn_2d_polar_precompute failed to open file ', &
          trim(filename)//'.nml'
        print *,'#at line',__LINE__
        print *,'#in file',__FILE__  
        STOP
      end if
    print *,'#initialization with filename:'
    print *,'#',trim(filename)//'.nml'
    read(input_file, params) 
    close(input_file)
  endif


  
  SLL_ALLOCATE(rho(num_cells_r+1),ierr)
  SLL_ALLOCATE(Ti(num_cells_r+1),ierr)
  SLL_ALLOCATE(mu_points(N_mu),ierr)
  SLL_ALLOCATE(mu_weights(N_mu),ierr)
  SLL_ALLOCATE(lambda(num_cells_r+1),ierr)
  SLL_ALLOCATE(points(3,N_points),ierr)
  SLL_ALLOCATE(mat_loc(num_cells_theta,num_cells_r+1,num_cells_r+1), ierr)
  SLL_ALLOCATE(mat(num_cells_theta,num_cells_r+1,num_cells_r+1), ierr)


  
  select case(quadrature_case)
    case ("SLL_RECTANGLE")
      do i=1,N_mu
        mu_points(i) = real(i-1,f64)*mu_max/real(N_mu,f64)
        mu_weights(i) = mu_max/real(N_mu,f64)*exp(-mu_points(i))
      enddo
    case ("SLL_GAUSS_LOBATTO")
      num_cells = N_mu/quadrature_points_per_cell
      mu_points(1:N_mu) = mu_max
      mu_weights(1:N_mu) = 0._f64
      s=1
      do i=1,num_cells
        mu_points(s:s+quadrature_points_per_cell-1) = &
          gauss_lobatto_points( &
            quadrature_points_per_cell, &
            real(i-1,f64)/real(num_cells,f64)*mu_max, &
            real(i,f64)/real(num_cells,f64)*mu_max )
        mu_weights(s:s+quadrature_points_per_cell-1) = &
          gauss_lobatto_weights( &
            quadrature_points_per_cell, &
            real(i-1,f64)/real(num_cells,f64)*mu_max, &
            real(i,f64)/real(num_cells,f64)*mu_max )
        s=s+quadrature_points_per_cell        
      enddo
      !mu_points(1:N_mu) = gauss_lobatto_points( N_mu, 0._f64, mu_max )
      !mu_weights(1:N_mu) = gauss_lobatto_weights( N_mu, 0._f64, mu_max )
      do i=1,N_mu
         mu_weights(i) = mu_weights(i)*exp(-mu_points(i))
      enddo       
    case ("SLL_GAUSS_LEGENDRE")
      num_cells = N_mu/quadrature_points_per_cell
      mu_points(1:N_mu) = mu_max
      mu_weights(1:N_mu) = 0._f64
      s=1
      do i=1,num_cells
        mu_points(s:s+quadrature_points_per_cell-1) = &
          gauss_legendre_points( &
            quadrature_points_per_cell, &
            real(i-1,f64)/real(num_cells,f64)*mu_max, &
            real(i,f64)/real(num_cells,f64)*mu_max )
        mu_weights(s:s+quadrature_points_per_cell-1) = &
          gauss_legendre_weights( &
            quadrature_points_per_cell, &
            real(i-1,f64)/real(num_cells,f64)*mu_max, &
            real(i,f64)/real(num_cells,f64)*mu_max )
        s=s+quadrature_points_per_cell        
      enddo
      !mu_points(1:N_mu) = gauss_lobatto_points( N_mu, 0._f64, mu_max )
      !mu_weights(1:N_mu) = gauss_lobatto_weights( N_mu, 0._f64, mu_max )
      do i=1,N_mu
         mu_weights(i) = mu_weights(i)*exp(-mu_points(i))
      enddo       
    case default
      print *,'#bad quadrature_case',trim(quadrature_case)
      print *,'#not implemented'
      print *,'#in initialize_analytic_field_2d_curvilinear'
      print*,'#at line and file:',__LINE__,__FILE__
      stop
  end select
  

  do i=1,num_cells_r+1
    Ti(i) = 1._f64
  enddo


  
  lambda = lambda_coeff*1._f64
  
  !SLL_ALLOCATE(mat(num_cells_r+1,num_cells_r+1,num_cells_theta),ierr)
  
  
  
  call compute_shape_circle(points,N_points)
      
  
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

  call sll_ascii_file_create("mat.dat", mat_id, ierr)

   do m = 1,num_cells_theta
     do i= 1,num_cells_r+1
       do j= 1,num_cells_r+1
         write(mat_id,*) mat(m,i,j)
       enddo
     enddo
   enddo

  call sll_ascii_file_close(mat_id,ierr)

  call compute_qns_inverse_polar_splines( &
    mat, &
    lambda, &
    num_cells_r, &
    num_cells_theta)


  Nc = (/num_cells_r,num_cells_theta/)
  eta_min = (/r_min,0._f64/)
  eta_max = (/r_max,2._f64*sll_pi/)
  mode = (/mode_r,mode_theta/)

  SLL_ALLOCATE(phi(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_init(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_qn(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_restr(1:Nc(1)+1,1:Nc(2)),ierr)
  SLL_ALLOCATE(phi_quotient(1:Nc(1)+1,1:Nc(2)+1),ierr)


  
  call compute_init_f_polar(phi,mode,Nc,eta_min,eta_max)
  phi_init = phi
  phi_restr(1:Nc(1)+1,1:Nc(2)) = phi(1:Nc(1)+1,1:Nc(2))



  call solve_qns_polar_splines( &
    mat, &
    phi_restr, &
    num_cells_r, &
    num_cells_theta)

  phi_qn(1:Nc(1)+1,1:Nc(2)) = phi_restr(1:Nc(1)+1,1:Nc(2))
  phi_qn(1:Nc(1)+1,Nc(2)+1) = phi_restr(1:Nc(1)+1,1)

  do j=1,Nc(2)+1
    do i=1,Nc(1)+1
      if(abs(phi_init(i,j))>0.2_f64)then
        phi_quotient(i,j) = abs(phi_qn(i,j))/(abs(phi_init(i,j))+1e-6_f64)
      else
        phi_quotient(i,j) = 0._f64
      endif  
    enddo
  enddo  
  !do i=1,N_mu
  !  mu_weights(i) = mu_weights(i)*exp(mu_points(i))
  !enddo       

  call sll_gnuplot_corect_2d( &
    r_min, &
    r_max, &
    num_cells_r, &
    0._f64, &
    2._f64*sll_pi, &
    num_cells_theta, &
    phi_init, &
    "phi_init", &
    1, &
    ierr)

  call sll_gnuplot_corect_2d( &
    r_min, &
    r_max, &
    num_cells_r, &
    0._f64, &
    2._f64*sll_pi, &
    num_cells_theta, &
    phi_qn, &
    "phi_qn", &
    1, &
    ierr)

  call sll_gnuplot_corect_2d( &
    r_min, &
    r_max, &
    num_cells_r, &
    0._f64, &
    2._f64*sll_pi, &
    num_cells_theta, &
    phi_quotient, &
    "phi_quotient", &
    1, &
    ierr)


  gamma0_num = compute_gamma0_quadrature( &
    Nc, &
    eta_min, &
    eta_max, &
    mu_points, &
    mu_weights, &
    N_mu, &
    mode )
  
  print *,'gamma0 num=',gamma0_num
  
  call compute_gamma0(mode,eta_min,eta_max,gamma0)

  print *,'gamma0=',gamma0,gamma0-gamma0_num
  
  
  call test_solve_qn_polar_new(Nc,eta_min,eta_max, &
    mu_points(1:N_mu),mu_weights(1:N_mu),N_mu,mode,phi_init,phi_qn)

  
  
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

    print *, '#reached end of vp2d test'
    time = sll_time_elapsed_since(t0)
    print *, '#time elapsed since t0 : ',time

  
  print *,'#PASSED'
  

end program