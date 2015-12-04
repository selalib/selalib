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

program test_qn_2d_polar_precompute
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_m_qn_2d_polar_precompute
use sll_m_timer
use sll_m_ascii_io
use sll_m_gauss_lobatto_integration
use sll_m_gauss_legendre_integration
use sll_m_gnuplot
use sll_m_qn_2d_polar
use sll_m_gyroaverage_2d_polar_hermite_solver

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
  sll_real64 :: mu_min
  sll_real64 :: mu_max
  character(len=256) :: filename
  sll_int32             :: IO_stat
  sll_int32, parameter  :: input_file = 99
  type(sll_time_mark)  :: t0
  sll_real64 :: time
  sll_int32 :: j
  sll_int32 :: m
  sll_int32 :: mat_id
  !sll_int32 :: nml_id
  sll_real64 :: lambda_coeff
  character(len=256) :: quadrature_case
  sll_int32 :: quadrature_points_per_cell
  !sll_int32 :: num_cells
  !sll_int32 :: s
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
  sll_real64 :: error(3)
  sll_real64,dimension(:,:),allocatable :: val
  sll_real64 :: r
  sll_real64 :: theta
  sll_real64 :: r_minus
  sll_real64 :: r_plus
  sll_real64 :: eps
  sll_real64 :: num_buffer_cell
  sll_int32  :: hermite_case
  sll_int32  :: interp_degree(2)
  !class(sll_gyroaverage_2d_base), pointer :: gyroaverage
  sll_int32 :: interp_x1
  sll_int32 :: interp_x2
  !sll_real64 :: mass 
  !sll_real64 :: vol
    
  namelist /params/ &
    r_min, &
    r_max, &
    num_cells_r, &
    num_cells_theta, &
    N_mu, &
    N_points, &
    mu_min, &
    mu_max, &
    lambda_coeff, &
    quadrature_case, &
    quadrature_points_per_cell, &
    mode_r, &
    mode_theta, &
    r_minus, &
    r_plus, &
    eps, &
    num_buffer_cell, &
    interp_x1, &
    interp_x2

  call sll_set_time_mark(t0)


  !default parameters
  r_min = 2._f64
  r_max = 18._f64
  r_minus = 0.4_f64
  r_plus = 0.5_f64
  num_cells_r = 32
  num_cells_theta = 64
  N_mu = 2
  N_points = 32
  mu_min = 0._f64
  mu_max = 3._f64
  lambda_coeff = 1._f64
  quadrature_case = "SLL_RECTANGLE"
  quadrature_points_per_cell = 10
  mode_r = 1
  mode_theta = 1
  eps = 1.e-3_f64
  num_buffer_cell = 1._f64
  interp_x1 = 3
  interp_x2 = 3

  
  hermite_case = 2

  
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

  interp_degree(1) = interp_x1
  interp_degree(2) = interp_x2


  
  SLL_ALLOCATE(rho(num_cells_r+1),ierr)
  SLL_ALLOCATE(Ti(num_cells_r+1),ierr)
  SLL_ALLOCATE(mu_points(N_mu),ierr)
  SLL_ALLOCATE(mu_weights(N_mu),ierr)
  SLL_ALLOCATE(lambda(num_cells_r+1),ierr)
  SLL_ALLOCATE(points(3,N_points),ierr)
  SLL_ALLOCATE(mat_loc(num_cells_theta,num_cells_r+1,num_cells_r+1), ierr)
  SLL_ALLOCATE(mat(num_cells_theta,num_cells_r+1,num_cells_r+1), ierr)
  SLL_ALLOCATE(val(2,num_cells_r+1), ierr)

  
  call compute_mu( &
    quadrature_case, &
    mu_points, &
    mu_weights, &
    N_mu, &
    mu_min, &
    mu_max, &
    quadrature_points_per_cell)
  

  do i=1,num_cells_r+1
    Ti(i) = 1._f64
  enddo


  
  lambda = lambda_coeff*1._f64
  
  !SLL_ALLOCATE(mat(num_cells_r+1,num_cells_r+1,num_cells_theta),ierr)
  
  
  
  call compute_shape_circle(points,N_points)


  Nc = (/num_cells_r,num_cells_theta/)
  eta_min = (/r_min,0._f64/)
  eta_max = (/r_max,2._f64*sll_pi/)
  mode = (/mode_r,mode_theta/)


  gamma0_num = compute_gamma0_quadrature( &
    Nc, &
    eta_min, &
    eta_max, &
    mu_points, &
    mu_weights, &
    N_mu, &
    mode )
  
  print *,'gamma0 num=',gamma0_num
  
  gamma0 = 1._f64
  
  call compute_gamma0(mode,eta_min,eta_max,gamma0)

  print *,'gamma0=',gamma0,gamma0-gamma0_num,gamma0-(1._f64-gamma0_num)
  
  if(gamma0==1)then
    print *,'#value is not ok'
    print *,'#to test you should have zeros_bessel.txt'
    print *,'#in build directory'
    print *,'#and rmin=2 and rmax=18'
    print *,'#PASSED'
    stop
  endif
  
  !stop
      
  
  mat = (0._f64,0._f64)
  
  
  !print *,mu_points,mu_weights
  do i=1,N_mu
    rho(1:num_cells_r+1) = sqrt(2._f64*mu_points(i))*Ti(1:num_cells_r+1)
!    do j=1,num_cells_r+1
!      r = r_min+real(j-1,f64)*(r_max-r_min)/real(num_cells_r,f64)
!      rho(j) = minval((/rho(j),abs(r-r_min)/num_buffer_cell,abs(r-r_max)/num_buffer_cell/))
!    enddo
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

  mat_loc = mat

  call compute_qns_inverse_polar_splines( &
    mat, &
    lambda, &
    num_cells_r, &
    num_cells_theta)

!  call compute_qns_inverse_polar_splines( &
!    mat, &
!    lambda, &
!    num_cells_r, &
!    num_cells_theta)
!
!  print *,'#error for matrix fdgrs=',maxval(abs(mat-mat_loc))



  SLL_ALLOCATE(phi(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_init(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_qn(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_restr(1:Nc(1)+1,1:Nc(2)),ierr)
  SLL_ALLOCATE(phi_quotient(1:Nc(1)+1,1:Nc(2)+1),ierr)


  
  call compute_init_f_polar(phi,mode,Nc,eta_min,eta_max)
  
  print *,'#mode=',mode
  
  do i=1,Nc(1)+1
    r = r_min+real(i-1,f64)*(r_max-r_min)/real(Nc(1),f64)
    do j=1,Nc(2)+1
      theta = real(j-1,f64)*2._f64*sll_pi/real(Nc(2),f64)
      phi(i,j) = 1._f64
!if((r>=r_minus).and.(r<=r_plus))then
!  phi(i,j) = (1.0_f64+eps*cos(mode_theta*theta))
!else
!  phi(i,j) = 0._f64  
!endif 
phi(i,j) = (1.0_f64+eps*cos(mode_theta*theta)) *exp(-(r-10._f64)**2/4._f64) !*sin(2._f64*sll_pi*(r-r_min)/(r_max-r_min))

      
    enddo
  enddo

!  vol = 0._f64
!  mass = 0._f64
!  do i=1,Nc(1)+1
!    r = r_min+real(i-1,f64)*(r_max-r_min)/real(Nc(1),f64)
!    do j=1,Nc(2)
!      mass = mass+r*phi(i,j)
!      vol = vol+r
!    enddo
!  enddo
!  
!  !mass = mass*(r_max-r_min)/real(Nc(1),f64)*(2._f64*sll_pi)/real(Nc(2),f64)
!  
!  mass = mass/vol
!  
!  phi = phi-mass !/(sll_pi*(r_max**2-r_min**2))
  
  phi_init = phi
  phi_restr(1:Nc(1)+1,1:Nc(2)) = phi(1:Nc(1)+1,1:Nc(2))

! *** old version before prototype/src move
!  call sll_gnuplot_corect_2d( &
!    r_min, &
!    r_max, &
!    num_cells_r+1, &
!    0._f64, &
!    2._f64*sll_pi, &
!    num_cells_theta+1, &
!    phi_init, &
!    "phi_init", &
!    1, &
!    ierr)



!  gyroaverage => new_gyroaverage_2d_polar_hermite_solver( &
!    eta_min, &
!    eta_max, &
!    Nc, &
!    N_points, &
!    interp_degree, &
!    rho(Nc(1)/2), &
!    hermite_case)
!
!  call gyroaverage%compute_gyroaverage( rho(Nc(1)/2), phi)
!  phi_restr(1:Nc(1)+1,1:Nc(2)) = phi_restr(1:Nc(1)+1,1:Nc(2))-phi(1:Nc(1)+1,1:Nc(2))


  call solve_qns_polar_splines( &
    mat, &
    phi_restr, &
    num_cells_r, &
    num_cells_theta)


!  vol = 0._f64
!  mass = 0._f64
!  do i=1,Nc(1)+1
!    r = r_min+real(i-1,f64)*(r_max-r_min)/real(Nc(1),f64)
!    do j=1,Nc(2)
!      mass = mass+r*phi_restr(i,j)
!      vol = vol+r
!    enddo
!  enddo
!  
!  !mass = mass*(r_max-r_min)/real(Nc(1),f64)*(2._f64*sll_pi)/real(Nc(2),f64)
!  
!  mass = mass/vol
!  
!  phi_restr = phi_restr-mass !/(sll_pi*(r_max**2-r_min**2))
!


!  call solve_qns_polar_splines( &
!    mat_loc, &
!    phi_restr, &
!    num_cells_r, &
!    num_cells_theta)

  do i=1,Nc(2)
    
    do j=1,Nc(1)+1
      val(1,j) = sum(abs(mat(i,j,:)))
      val(2,j) = sum(abs(mat_loc(i,j,:)))
    enddo
    print *,maxval(val(1,:))*maxval(val(2,:))
  enddo

!  print *,'#error with phi',maxval(abs(phi_restr-phi_init(1:Nc(1)+1,1:Nc(2))))


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


! *** old version before prototype/src move
!  call sll_gnuplot_corect_2d( &
!    r_min, &
!    r_max, &
!    num_cells_r+1, &
!    0._f64, &
!    2._f64*sll_pi, &
!    num_cells_theta+1, &
!    phi_qn, &
!    "phi_qn", &
!    1, &
!    ierr)
!
!  call sll_gnuplot_corect_2d( &
!    r_min, &
!    r_max, &
!    num_cells_r+1, &
!    0._f64, &
!    2._f64*sll_pi, &
!    num_cells_theta+1, &
!    phi_quotient, &
!    "phi_quotient", &
!    1, &
!    ierr)




    print *,'#test double gyroaverage'
    call compute_error(phi_qn,phi_init,1-gamma0_num,error,(/1,1/),Nc)
    print *,'#error whole domain=',error

!stop  
  
!  call test_solve_qn_polar_new(Nc,eta_min,eta_max, &
!    mu_points(1:N_mu),mu_weights(1:N_mu),N_mu,mode,phi_init,phi_qn)

  
  
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
  

end program test_qn_2d_polar_precompute
