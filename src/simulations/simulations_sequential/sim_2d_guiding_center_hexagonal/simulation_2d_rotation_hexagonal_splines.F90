program rotation_2d_hexagonal_splines

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_constants
  use sll_hex_meshes
  use sll_box_splines
  implicit none

  type(sll_hex_mesh_2d),   pointer        :: mesh
  type(sll_box_spline_2d), pointer        :: spline
  sll_int32 :: num_cells
  sll_real64 :: center_mesh_x1
  sll_real64 :: center_mesh_x2
  sll_real64 :: radius
  sll_real64 :: t_init
  sll_real64 :: t_end
  sll_int32 :: spline_degree
  sll_real64 :: dt
  sll_int32 :: number_iterations
  sll_real64 :: gauss_x1
  sll_real64 :: gauss_x2
  sll_real64 :: gauss_sig
  sll_real64 :: gauss_amp
  
  
  sll_int32 :: n_points
  sll_int32 :: ierr
  sll_real64, dimension(:), allocatable :: rho_tn
  sll_real64, dimension(:), allocatable :: rho_tn1
  sll_real64, dimension(:), allocatable :: rho_exact
  sll_int32 :: step
  sll_int32 :: i
  sll_real64 :: x
  sll_real64 :: y
  sll_real64 :: xx
  sll_real64 :: yy
  sll_real64 :: r11
  sll_real64 :: r12
  sll_real64 :: r21
  sll_real64 :: r22
  sll_real64 :: det
  logical :: inside
  sll_real64 :: h1
  sll_real64 :: h2
  sll_real64 :: t
  sll_real64 :: linf_err

  call cpu_time(t_init)
  linf_err = 0._f64
  
  center_mesh_x1 = 0._F64
  center_mesh_x2 = 0._F64
  radius = 8._f64
  num_cells = 80
  spline_degree = 2
  number_iterations = 300
  dt = 0.1_f64

  gauss_x1  = 2._f64
  gauss_x2  = 2._f64
  gauss_sig = 1._f64/( 2._f64 * sqrt(2._f64)) 
  gauss_amp = 1.0_f64
  
  mesh => new_hex_mesh_2d( &
    num_cells, &
    center_mesh_x1, &
    center_mesh_x2,&
    radius=radius, &
    EXTRA_TABLES = 0 )
  !call display_hex_mesh_2d(mesh)
  n_points = mesh%num_pts_tot
  
  
  spline => new_box_spline_2d(mesh, SLL_DIRICHLET)

  SLL_ALLOCATE(rho_tn(n_points),ierr)
  SLL_ALLOCATE(rho_tn1(n_points),ierr)
  SLL_ALLOCATE(rho_exact(n_points),ierr)

  !initialization
  call gaussian_at_t( &
    rho_tn, &
    0._f64, &
    mesh, &
    gauss_x1, &
    gauss_x2, &
    gauss_sig, &
    gauss_amp)

  det = (mesh%r1_x1*mesh%r2_x2-mesh%r1_x2*mesh%r2_x1)/mesh%delta
  r11 = + mesh%r2_x2/det
  r12 = - mesh%r2_x1/det
  r21 = - mesh%r1_x2/det
  r22 = + mesh%r1_x1/det
  
  
  
  
  do step=1,number_iterations
  
    call compute_coeff_box_spline_2d( rho_tn, spline_degree, spline )
    do i=1, n_points
      x = mesh%cartesian_coord(1,i)
      y = mesh%cartesian_coord(2,i)
      xx = x*cos(dt) - y*sin(dt)
      yy = x*sin(dt) + y*cos(dt)


      inside = .true.
      h1 =  xx*r11 + yy*r12
      h2 =  xx*r21 + yy*r22

      if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
      if ( abs(xx) > (radius-mesh%delta)*sll_sqrt3*0.5_f64) inside = .false.

      if ( inside ) then
        rho_tn1(i) = hex_interpolate_value(mesh, xx, yy, spline, spline_degree)
      else
        rho_tn1(i) = 0._f64 ! dirichlet boundary condition
      endif

  
    enddo
    rho_tn = rho_tn1
    
    t = real(step,f64)*dt
    call gaussian_at_t( &
      rho_exact, &
      t, &
      mesh, &
      gauss_x1, &
      gauss_x2, &
      gauss_sig, &
      gauss_amp)
    linf_err = max(linf_err,maxval(abs(rho_exact-rho_tn)))
      
  enddo
  call cpu_time(t_end)

  
  print *,'#time/linf_error',t_end-t_init,linf_err
  print *,"#PASSED"

contains

  subroutine gaussian_at_t(f_tn, t, mesh, center_x1, center_x2, sigma, amplitude)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:)       :: f_tn
    sll_real64, intent(in) :: t
    sll_real64, intent(in) :: center_x1
    sll_real64, intent(in) :: center_x2
    sll_real64, intent(in) :: sigma
    sll_real64, intent(in) :: amplitude
    sll_real64 :: x
    sll_real64 :: y
    sll_int32  :: i
    sll_real64 :: xx
    sll_real64 :: yy

    do i = 1,mesh%num_pts_tot
       x = mesh%cartesian_coord(1,i)
       y = mesh%cartesian_coord(2,i)
       xx = x*cos(t) - y*sin(t)
       yy = x*sin(t) + y*cos(t)       
       f_tn(i) = amplitude * exp(-0.5_f64* &
            ((xx-center_x1)**2 + (yy-center_x2)**2) / sigma**2 )
    enddo
  end subroutine gaussian_at_t


end program rotation_2d_hexagonal_splines