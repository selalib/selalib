!>$L_x$ domain dimensions and M is an integer.
!>$
!>$
!>$
!>B_z(x,y,t) =   \cos(\frac{2 M \pi}{L_x} x)  \cos(\frac{2 M \pi}{L_x} t)
!>$
!>$
!>E_y(x,y,t) = \sin(\frac{2 M \pi}{L_x} x)  \sin(\frac{2 M \pi}{L_x} t)
!>$
!
!  Contact : Eric Sonnendrucker, Katharina Kormann
!
program test_maxwell_clamped_1d_fem_sm
  !------------------------------------------------------------------------
  !  test 1D Maxwell_Clamped spline finite element solver on a periodic grid
  !------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_low_level_bsplines, only: &
       sll_s_eval_uniform_periodic_spline_curve

  use sll_m_splines_pp

  use sll_m_constants, only: &
    sll_p_pi, sll_p_twopi

  use sll_m_maxwell_1d_base, only: &
    sll_s_plot_two_fields_1d

  use sll_m_maxwell_clamped_1d_fem_sm, only: &
    sll_t_maxwell_clamped_1d_fem_sm

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !  type(sll_t_arbitrary_degree_spline_1d), pointer :: aspl
  !  sll_real64, dimension(:), allocatable :: knots

  sll_real64 :: eta1_max, eta1_min
  sll_real64 :: delta_eta1

  sll_int32  :: nc_eta1
  sll_int32  :: error
  sll_real64 :: tol

  type(sll_t_maxwell_clamped_1d_fem_sm)  :: maxwell_clamped_1d

  sll_real64, dimension(:), allocatable :: ex
  sll_real64, dimension(:), allocatable :: ey
  sll_real64, dimension(:), allocatable :: bz
  sll_real64, dimension(:), allocatable :: bz_exact
  sll_real64, dimension(:), allocatable :: ex_exact
  sll_real64, dimension(:), allocatable :: ey_exact
  sll_real64, dimension(:), allocatable :: rho
  sll_real64, dimension(:), allocatable :: current
  sll_real64, dimension(:), allocatable :: sval
  sll_real64, dimension(:), allocatable :: x
  sll_real64, allocatable :: efield(:), bfield(:)

  sll_int32                               :: i!, j
  sll_real64                              :: time
  sll_int32                               :: istep, nstep
  sll_real64                              :: err_ex
  sll_real64                              :: err_ex2
  sll_real64                              :: err_ey
  sll_real64                              :: err_bz
  sll_real64                              :: err_l2norm
  sll_real64                              :: err_inner
  
  sll_real64                              :: dt
  sll_real64                              :: Lx
  sll_real64, dimension(2)                :: domain
  sll_int32                               :: deg, boundary
  sll_int32,  parameter                   :: mode = 1
  sll_real64                              :: l2norm
  sll_real64 :: alpha = 0.01_f64
  sll_real64 :: k = 0.8_f64

  ! Define computational domain
  eta1_min = .0_f64; eta1_max = sll_p_twopi/k
  nc_eta1 = 128
  Lx = eta1_max-eta1_min
  delta_eta1 = Lx/real(nc_eta1,f64)
  domain = [eta1_min, eta1_max]
  ! Set spline degree of 0-forms
  deg = 3
  if( deg == 2 ) then
     boundary =  sll_p_boundary_clamped_square
  else if( deg == 3 ) then
     boundary = sll_p_boundary_clamped_cubic
  else
     boundary = sll_p_boundary_clamped
  end if
  ! Initialise maxwell_clamped FEM object
  call maxwell_clamped_1d%init(domain, nc_eta1, deg, boundary, mass_tolerance=1d-12, poisson_tolerance=1d-12)

  allocate( x(1:nc_eta1+1) )

  do i = 1, nc_eta1+1
    x(i) = eta1_min + real(i-1,f64)*delta_eta1
  end do

  ! Allocate arrays
  SLL_CLEAR_ALLOCATE(ex(1:nc_eta1+deg-1),error)
  SLL_CLEAR_ALLOCATE(ey(1:nc_eta1+deg),error)
  SLL_CLEAR_ALLOCATE(bz(1:nc_eta1+deg-1),error)
  SLL_CLEAR_ALLOCATE(rho(1:nc_eta1+deg),error)
  SLL_CLEAR_ALLOCATE(current(1:nc_eta1+deg-1),error)
  
  SLL_CLEAR_ALLOCATE(bz_exact(1:nc_eta1+1),error)
  SLL_CLEAR_ALLOCATE(ex_exact(1:nc_eta1+1),error)
  SLL_CLEAR_ALLOCATE(ey_exact(1:nc_eta1+1),error)
  
  SLL_CLEAR_ALLOCATE(sval(1:nc_eta1+1),error)

  SLL_CLEAR_ALLOCATE(efield(1:nc_eta1+deg),error)
  SLL_CLEAR_ALLOCATE(bfield(1:nc_eta1+deg-1),error)

  
  ex = 0._f64
  ey = 0._f64
  rho = 0._f64
!!$  call maxwell_clamped_1d%L2projection( polynom, deg-1, ex )
!!$  call maxwell_clamped_1d%compute_rho_from_e(ex, ey)
!!$  call maxwell_clamped_1d%linear_solver_mass0%solve(ey, rho)
!!$  call sll_s_spline_evaluate_basis_b_form_1d_clamped( maxwell_clamped_1d%spline0_pp, nc_eta1,  rho, sval)   
!!$  ! Set exact solution
!!$  do i = 1, nc_eta1+1
!!$     ex_exact(i) = constant(x(i))
!!$  end do
!!$  err_ex = maxval(abs(sval-ex_exact))
!!$  print*, 'Error L2projectionE1',  err_ex
!!$  call sll_s_plot_two_fields_1d('rho',nc_eta1+1,sval,ex_exact,0,0.0_f64)
!!$
!!$  ey = 0._f64
!!$  ex = 0._f64
!!$  call maxwell_clamped_1d%L2projection( polynom, deg, ey )
!!$  call maxwell_clamped_1d%multiply_g( ey, ex )
!!$  call sll_s_spline_evaluate_basis_b_form_1d_clamped( maxwell_clamped_1d%spline1_pp, nc_eta1,  ex, sval) 
!!$  ! Set exact solution
!!$  do i = 1, nc_eta1+1
!!$     ey_exact(i) = constant(x(i))
!!$  end do
!!$  err_ey = maxval(abs(sval-ey_exact))
!!$  print*, 'Error L2projectionE2',  err_ey
!!$  call sll_s_plot_two_fields_1d('g',nc_eta1+1,sval,ey_exact,0,0.0_f64)
!!$ 
!!$  bz = 0._f64
!!$  ey = 0._f64
!!$  call maxwell_clamped_1d%L2projection( polynom, deg-1, bz )
!!$  call sll_s_spline_evaluate_basis_b_form_1d_clamped( maxwell_clamped_1d%spline1_pp, nc_eta1,  bz, sval) 
!!$  ! Set exact solution
!!$  do i = 1, nc_eta1+1
!!$     bz_exact(i) = polynom(x(i))
!!$  end do
!!$  err_bz = maxval(abs(sval-bz_exact))
!!$  print*, 'Error L2projectionB3',  err_bz
!!$  call sll_s_plot_two_fields_1d('bzcos',nc_eta1+1,sval,bz_exact,0,0.0_f64)


  ex = 0._f64
  rho = 0._f64
  ! Test Poisson
  !-------------
  ! Set exact solution
  do i = 1, nc_eta1+1
     ex_exact(i) =  E_k(x(i))!polynomo(x(i))!
  end do

  call maxwell_clamped_1d%compute_rhs_from_function( rho_k, deg, rho)!polynom, deg, rho)!
  !print*, 'rho', rho
  !rho(1)=rho(1)+E_k(domain(1))!+polynomo(domain(1))!
  !rho(nc_eta1+deg)=rho(nc_eta1+deg)-E_k(domain(2))!-polynomo(domain(2))!
  call maxwell_clamped_1d%compute_e_from_rho( rho, ex )
  ! Evaluate spline curve at grid points and compute error
  ! Ex is a 1-form, i.e. one spline degree lower
  call sll_s_spline_evaluate_basis_b_form_1d_clamped( maxwell_clamped_1d%spline1_pp, nc_eta1,  ex, sval) 
  err_ex = maxval(abs(sval-ex_exact))
  print*, 'error Poisson',  err_ex

  !call sll_s_plot_two_fields_1d('poisson',nc_eta1+1,sval,ex_exact,0,0.0_f64)
 
  ! Test Ampere
  !-------------
  ! Set time step
  dt = .5_f64 * delta_eta1
  ! Set exact solution
  do i = 1, nc_eta1+1
     ex_exact(i) = -sin_k(x(i))*dt
  end do

  call maxwell_clamped_1d%compute_rhs_from_function(sin_k, deg-1, current)
  ex = 0.0_f64
  call maxwell_clamped_1d%compute_E_from_j(current*dt, 1, ex )

  ! Evaluate spline curve at grid points and compute error
  ! Ex is a 1-form, i.e. one spline degree lower
  call sll_s_spline_evaluate_basis_b_form_1d_clamped( maxwell_clamped_1d%spline1_pp, nc_eta1,  ex, sval)
  err_ex2 = maxval(abs(sval-ex_exact))
  print*, 'error Ampere',  err_ex2
  !call sll_s_plot_two_fields_1d('current',nc_eta1+1,sval,ex_exact,0,0.0_f64)

  l2norm =  maxwell_clamped_1d%l2norm_squared(ex,deg-1)
  err_l2norm = abs(l2norm - dt**2*sll_p_pi)
  print*, 'error l2 norm', err_l2norm


!!$  ! Test mixed mass
!!$  ex = 0.0_f64
!!$  ey = 0.0_f64
!!$  call maxwell_clamped_1d%L2projection( sin_k, deg-1, ex )
!!$  call maxwell_clamped_1d%L2projection( sin_k, deg, ey )
!!$  err_inner =  abs(maxwell_clamped_1d%inner_product( ex, ey, deg-1, deg ) - sll_p_pi)
!!$  print*, 'error inner product', err_inner
  err_inner = 0._f64

  ! Test Maxwell_Clamped on By and Ez 
  !--------------------------
  ! Set time stepping parameters
  time  = 0.0_f64
  dt = .5_f64 * delta_eta1
  nstep = 10

  ! Compute initial fields 
  ex = 0.0_f64
  call maxwell_clamped_1d%L2projection(cos_k, deg-1, bz) ! 1-form -> splines of degree deg-1
  ey = 0.0_f64  ! 0-form -> splines of degree deg
  
  ! Time loop. Second order Strang splitting
  do istep = 1, nstep
     efield = ey
     bfield = bz

!!$     call maxwell_clamped_1d%compute_b_from_e( 0.5_f64*dt, ey, bz)
!!$     call maxwell_clamped_1d%compute_e_from_b(         dt, bz, ey)
!!$     call maxwell_clamped_1d%compute_b_from_e( 0.5_f64*dt, ey, bz)
     call maxwell_clamped_1d%compute_curl_part( dt, ey, bz, 1._f64 )

  
     time = time + dt

     do i = 1, nc_eta1+1
        ey_exact(i) =   sin(mode*sll_p_twopi*x(i)/Lx) * sin(mode*sll_p_twopi*time/Lx)
        bz_exact(i) =   cos(mode*sll_p_twopi*x(i)/Lx) * cos(mode*sll_p_twopi*time/Lx)
     end do
     call sll_s_spline_evaluate_basis_b_form_1d_clamped( maxwell_clamped_1d%spline0_pp, nc_eta1,  ey, sval)
     err_ey = maxval(abs(sval-ey_exact))
     !call sll_s_plot_two_fields_1d('ey',nc_eta1+1,sval,ey_exact,istep,time)
     
     call sll_s_spline_evaluate_basis_b_form_1d_clamped( maxwell_clamped_1d%spline1_pp, nc_eta1,  bz, sval)
     err_bz = maxval(abs(sval-bz_exact))
     !call sll_s_plot_two_fields_1d('bz',nc_eta1+1,sval,bz_exact,istep,time)

     write(*,"(10x,' istep = ',I6)",advance="no") istep
     write(*,"(' time = ',g12.3,' sec')",advance="no") time
     write(*,"(' erreur L2 = ',2g15.5)") err_ey, err_bz
  end do

  tol = 1.0d-3

  if ((err_bz < tol) .and. (err_ey < tol) .and. (err_ex < tol) .and. (err_ex2 < 4.0d-8) .and. (err_l2norm < 1.0d-3).and. (err_inner < 1.0d-10)) then
     print*,'PASSED'
  endif

  call maxwell_clamped_1d%free()

  DEALLOCATE(ex)
  DEALLOCATE(ey)
  DEALLOCATE(bz)
  DEALLOCATE(bz_exact)
  DEALLOCATE(ex_exact)
  DEALLOCATE(ey_exact)
  DEALLOCATE(rho) 

contains
  function rho_k(x)
    sll_real64             :: rho_k
    sll_real64, intent(in) :: x

    rho_k = - alpha*sin(k*x) 
  end function rho_k
  function phi_k(x)
    sll_real64             :: phi_k
    sll_real64, intent(in) :: x

    phi_k = - alpha/k**2*sin(k*x) 
  end function phi_k
  function E_k(x)
    sll_real64             :: E_k
    sll_real64, intent(in) :: x
    
    E_k = alpha/k*cos(k*x) 
  end function E_k

  function polynomo(x)
    sll_real64             :: polynomo
    sll_real64, intent(in) :: x

    polynomo = 2._f64*x**2
  end function polynomo

  function polynom(x)
    sll_real64             :: polynom
    sll_real64, intent(in) :: x

    polynom = 4._f64*x
  end function polynom

 function constant(x)
    sll_real64             :: constant
    sll_real64, intent(in) :: x

    constant = 4._f64
  end function constant
  
  
  function cos_k(x)
    sll_real64             :: cos_k
    sll_real64, intent(in) :: x

    cos_k = cos(mode*sll_p_twopi*x/Lx) 
  end function cos_k
  function sin_k(x)
    sll_real64             :: sin_k
    sll_real64, intent(in) :: x

    sin_k = sin(mode*sll_p_twopi*x/Lx) 
  end function sin_k


  
end program test_maxwell_clamped_1d_fem_sm
