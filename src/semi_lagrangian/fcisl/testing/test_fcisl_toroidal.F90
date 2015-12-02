program test_fcisl_toroidal
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_m_fcisl_toroidal
  use sll_m_constants
  use sll_m_interpolators_2d_base
  use sll_m_cubic_spline_interpolator_2d
  use sll_m_ascii_io

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
  sll_real64, dimension(:), allocatable :: theta_rk4
  sll_real64, dimension(:), allocatable :: phi_rk4
  sll_real64, dimension(:), allocatable :: theta_analytic
  sll_real64, dimension(:), allocatable :: phi_analytic
  sll_real64, dimension(:), allocatable :: theta
  sll_real64, dimension(:), allocatable :: phi
  sll_real64 :: theta0
  sll_real64 :: phi0
  sll_int32 :: ierr
  sll_int32 :: file_id
  sll_real64 :: err
  sll_real64, dimension(:), allocatable :: x_array
  sll_real64 :: xmin
  sll_real64 :: xmax
  sll_int32 :: Npts_x
  sll_int32 :: i
  sll_int32 :: j
  sll_real64, dimension(:), allocatable :: phi_array
  sll_int32 :: Npts_phi
  sll_real64, dimension(:), allocatable :: theta_array
  sll_int32 :: Npts_theta
  sll_real64 :: params(4)
  sll_real64, dimension(:,:), allocatable :: charac_theta_euler
  sll_real64, dimension(:,:), allocatable :: charac_phi_euler
  sll_real64, dimension(:,:), allocatable :: charac_theta_rk4
  sll_real64, dimension(:,:), allocatable :: charac_phi_rk4
  sll_real64, dimension(:,:), allocatable :: charac_theta_analytic
  sll_real64, dimension(:,:), allocatable :: charac_phi_analytic
  sll_real64, dimension(:,:), allocatable :: charac_theta
  sll_real64, dimension(:,:), allocatable :: charac_phi
  sll_real64, dimension(:,:), allocatable :: f_init
  sll_real64, dimension(:,:), allocatable :: f_classic
  sll_real64, dimension(:,:), allocatable :: f_aligned
  sll_real64, dimension(:,:), allocatable :: f
  sll_real64, dimension(:,:), allocatable :: f_exact
  sll_real64, dimension(:), allocatable :: params_aligned
  
  class(sll_c_interpolator_2d), pointer :: interp_classic
  sll_int32 :: mode_m
  sll_int32 :: mode_n
  sll_int32 :: nb_iter
  sll_real64 :: th
  sll_real64 :: ph
  sll_int32 :: iter
  sll_int32 :: hermite_p
  sll_int32 :: hermite_r_left
  sll_int32 :: hermite_r_right
  sll_int32 :: hermite_s_left
  sll_int32 :: hermite_s_right
  sll_int32 :: lag_p
  sll_int32 :: lag_r
  sll_int32 :: lag_s
  sll_real64 :: iota  
  sll_real64, dimension(:,:,:), allocatable :: hermite_w_aligned
  sll_int32, dimension(:,:), allocatable :: hermite_w_cell_aligned
  sll_real64, dimension(:,:), allocatable :: theta_pos_left
  sll_real64, dimension(:,:,:), allocatable :: buf
  
  
  nb_iter = 10
  iota = 0.5_f64
  
  mode_m = 20
  mode_n = 20
  
  hermite_p = 6
  lag_p = 9
  lag_r = -lag_p/2
  lag_s = (lag_p+1)/2
  
  R0 = 10._f64
  smallr = 2._f64
  psipr = 4._f64
  dt = 0.01_f64
  num_time_points = 15000
  theta0 = 0._f64
  phi0 = 0._f64
  
  xmax=10._f64
  xmin=-xmax
  Npts_x=1000000
  Npts_phi=64
  Npts_theta=512
  

  F0=-psipr*R0/(smallr*iota)


psipr=4._f64;
F0=-2.4_f64*psipr;
R0=10._f64;
smallr=4._f64; 

iota=-psipr*R0/(smallr*F0)

R0=10._f64;
F0=-psipr*R0/(smallr*iota)


print*,"iota=",iota  
  
  SLL_ALLOCATE(time_points(num_time_points),ierr)
  SLL_ALLOCATE(theta_euler(num_time_points),ierr)
  SLL_ALLOCATE(phi_euler(num_time_points),ierr)
  SLL_ALLOCATE(theta_rk4(num_time_points),ierr)
  SLL_ALLOCATE(phi_rk4(num_time_points),ierr)
  SLL_ALLOCATE(theta_analytic(num_time_points),ierr)
  SLL_ALLOCATE(phi_analytic(num_time_points),ierr)
  SLL_ALLOCATE(theta(num_time_points),ierr)
  SLL_ALLOCATE(phi(num_time_points),ierr)
  SLL_ALLOCATE(x_array(Npts_x),ierr)
  SLL_ALLOCATE(phi_array(Npts_phi),ierr)
  SLL_ALLOCATE(theta_array(Npts_theta),ierr)
  SLL_ALLOCATE(charac_theta_euler(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(charac_phi_euler(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(charac_theta_rk4(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(charac_phi_rk4(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(charac_theta_analytic(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(charac_phi_analytic(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(charac_theta(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(charac_phi(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(f_init(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(f_exact(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(f_classic(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(f_aligned(Npts_theta,Npts_phi),ierr)
  SLL_ALLOCATE(f(Npts_theta,Npts_phi),ierr)
  
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

  call compute_rk4_field( &
    R0, &
    time_points, &
    num_time_points, &
    psipr, &
    F0, &
    smallr, &
    theta_rk4, &
    phi_rk4, &
    theta0, &
    phi0)   
  call compute_modulo_vect(phi_rk4,phi,num_time_points,2._f64*sll_pi)
  call compute_modulo_vect(theta_rk4,theta,num_time_points,2._f64*sll_pi)
  call sll_ascii_file_create("rk4_field.dat", file_id, ierr)    
  call sll_ascii_write_array_1d(file_id,time_points,ierr,num_time_points,phi,theta)  
  call sll_ascii_file_close(file_id,ierr)
  !for visualization: gnuplot -e 'plot "rk4_field.dat" u 2:3 w p ps 1 lt 5'  


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
  err=max(maxval(abs(phi_rk4-phi_analytic)),maxval(abs(theta_rk4-theta_analytic)))
  
  print *,'#err between rk4 and analytic=',err
  
  call compute_linspace(x_array,xmin,xmax,Npts_x)
  err = 0._f64
  do i=1,Npts_x
    err= max(err,abs(compute_inverse_invR_integral(R0,smallr,compute_invR_integral(R0,smallr,x_array(i)))-x_array(i)))
    err= max(err,abs(compute_invR_integral(R0,smallr,compute_inverse_invR_integral(R0,smallr,x_array(i)))-x_array(i)))    
  enddo

  print *,'#err for invR integral and inverse=',err

  call compute_linspace(phi_array,0._f64,2._f64*sll_pi,Npts_phi)
  call compute_linspace(theta_array,0._f64,2._f64*sll_pi,Npts_theta)
  params(1) = R0
  params(2) = psipr 
  params(3) = F0
  params(4) = smallr
  
  call compute_feet_euler( &
    theta_array, &
    phi_array, &
    Npts_theta, &
    Npts_phi, &
    dt, &
    params, &
    charac_theta_euler, &
    charac_phi_euler) 

  call compute_feet_rk4( &
    theta_array, &
    phi_array, &
    Npts_theta, &
    Npts_phi, &
    dt, &
    params, &
    charac_theta_rk4, &
    charac_phi_rk4) 


  call compute_feet_analytic( &
    theta_array, &
    phi_array, &
    Npts_theta, &
    Npts_phi, &
    dt, &
    params, &
    charac_theta_analytic, &
    charac_phi_analytic) 
    
  err=max(maxval(abs(charac_theta_euler-charac_theta_analytic)), &
    maxval(abs(charac_phi_euler-charac_phi_analytic)))
  
  print *,'#err between euler and analytic for charac=',err

  err=max(maxval(abs(charac_theta_rk4-charac_theta_analytic)), &
    maxval(abs(charac_phi_rk4-charac_phi_analytic)))
  
  print *,'#err between rk4 and analytic for charac=',err
  
  
  !set initial function
  
  do j=1,Npts_phi
    do i=1,Npts_theta
      th = theta_array(i)
      ph = phi_array(j)
      f_init(i,j) = cos(real(mode_m,f64)*th-real(mode_n,f64)*ph)
    enddo
  enddo 


  call compute_feet_analytic( &
    theta_array, &
    phi_array, &
    Npts_theta, &
    Npts_phi, &
    dt, &
    params, &
    charac_theta, &
    charac_phi) 
  call compute_modulo_vect2d_inplace(charac_theta,Npts_theta,Npts_phi,2._f64*sll_pi)
  call compute_modulo_vect2d_inplace(charac_phi,Npts_theta,Npts_phi,2._f64*sll_pi)


  interp_classic => new_cubic_spline_interpolator_2d( &
    Npts_theta, &
    Npts_phi, &
    0._f64, &
    2._f64*sll_pi, &
    0._f64, &
    2._f64*sll_pi, &
    SLL_PERIODIC, &
    SLL_PERIODIC)
  f_classic = f_init
  err = 0._f64
  do iter=1,nb_iter
    call interp_classic%interpolate_array( &
      Npts_theta, &
      Npts_phi, &
      f_classic, &
      charac_theta, &
      charac_phi, &
      f)
    f_classic = f
    
    call compute_feet_analytic( &
      theta_array, &
      phi_array, &
      Npts_theta, &
      Npts_phi, &
      real(iter,f64)*dt, &
      params, &
      charac_theta_analytic, &
      charac_phi_analytic) 
    do j=1,Npts_phi
      do i=1,Npts_theta
        th = charac_theta_analytic(i,j)
        ph = charac_phi_analytic(i,j)
        f_exact(i,j) = cos(real(mode_m,f64)*th-real(mode_n,f64)*ph)
      enddo
    enddo 
    
    err = max(err,maxval(abs(f_exact-f)))
      
  enddo
  
  
  
  print *,'#err for classical method=',err
  
  SLL_ALLOCATE(params_aligned(11),ierr)
  params_aligned(1) = R0
  params_aligned(2) = psipr
  params_aligned(3) = F0 
  params_aligned(4) = smallr
  params_aligned(5) = real(hermite_p,f64) 
  params_aligned(6) = real(lag_r,f64)
  params_aligned(7) = real(lag_s,f64)
  params_aligned(8) = 0._f64
  params_aligned(9) = 2._f64*sll_pi
  params_aligned(10) = 0._f64
  params_aligned(11) = 2._f64*sll_pi

  f_aligned = f_init
  err = 0._f64
  do iter=1,nb_iter
    f = interpolate2d_toroidal( &
      Npts_theta, &
      Npts_phi, &
      f_aligned, &
      charac_theta, &
      charac_phi, &
      params_aligned)
    f_aligned = f
    
    call compute_feet_analytic( &
      theta_array, &
      phi_array, &
      Npts_theta, &
      Npts_phi, &
      real(iter,f64)*dt, &
      params, &
      charac_theta_analytic, &
      charac_phi_analytic) 
    do j=1,Npts_phi
      do i=1,Npts_theta
        th = charac_theta_analytic(i,j)
        ph = charac_phi_analytic(i,j)
        f_exact(i,j) = cos(real(mode_m,f64)*th-real(mode_n,f64)*ph)
      enddo
    enddo 
    
    err = max(err,maxval(abs(f_exact-f)))
      
  enddo

  print *,'#err for aligned method=',err

  
!  hermite_r_left = compute_hermite_r_left(hermite_p)
!  hermite_s_left = compute_hermite_s_left(hermite_p)
!  hermite_r_right = compute_hermite_r_right(hermite_p)
!  hermite_s_right = compute_hermite_s_right(hermite_p)
!  
!  SLL_ALLOCATE(hermite_w_aligned(4,hermite_s_left-hermite_r_left,Npts_theta),ierr)  
!  SLL_ALLOCATE(hermite_w_cell_aligned(hermite_s_left-hermite_r_left,Npts_theta),ierr)  
!  SLL_ALLOCATE(theta_pos_left(hermite_s_left-hermite_r_left,Npts_theta),ierr)  
!
!  call compute_w_hermite_aligned( &
!    hermite_w_aligned, &
!    hermite_w_cell_aligned, &
!    Npts_theta, &
!    hermite_r_left, &
!    hermite_s_left, &
!    theta_pos_left, &
!    0._f64, &
!    2._f64*sll_pi )
!
!
!  SLL_ALLOCATE(buf(9,Npts_theta,Npts_phi),ierr)  
!
!  call compute_hermite_derivatives_aligned( &
!    f, &
!    Npts_theta, &
!    Npts_phi, &
!    hermite_p, &
!    hermite_w_aligned, &
!    hermite_w_cell_aligned, &
!    hermite_w_aligned, &
!    hermite_w_cell_aligned, &
!    buf)
  
  print*,'#PASSED'


end program test_fcisl_toroidal
