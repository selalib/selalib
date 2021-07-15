!>$ physical domain is [0,2pi]  and M is an integer.
!>$ logical domian is [0,1], linear trafo from [0,1] to [0,2pi]
!>$
!>$
!>B_z(x,y,t) =   \cos(\frac{2 M \pi}{L_x} x)  \cos(\frac{2 M \pi}{L_x} t)
!>$
!>$
!>E_y(x,y,t) = \sin(\frac{2 M \pi}{L_x} x)  \sin(\frac{2 M \pi}{L_x} t)
!>$
!
!  Contact : Benedikt Perse
!
program test_maxwell_clamped_1d_trafo
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

  use sll_m_maxwell_clamped_1d_trafo, only: &
       sll_t_maxwell_clamped_1d_trafo

  use sll_m_3d_coordinate_transformations

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_real64 :: eta1_max, eta1_min
  sll_real64 :: delta_eta1

  sll_int32  :: nc_eta1
  sll_int32  :: error
  sll_real64 :: tol

  type(sll_t_maxwell_clamped_1d_trafo)    :: maxwell_clamped_1d
  type(sll_t_mapping_3d), pointer :: map  


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
  sll_int32                               :: deg, boundary
  sll_real64                              :: l2norm
  sll_real64 :: alpha = 1._f64!0.01_f64
  sll_real64 :: k = 1.0_f64

  sll_real64 :: params(12)
  sll_real64 :: factor(3,3)
  !Initialize mapping
  params=0._f64
  params(1)=sll_p_twopi/k
  params(2)=1._f64
  params(3)=1._f64

  params(4)=0.1_f64


  allocate(map)
  call map%init( params,&
       sll_f_orthogonal_x1,&
       sll_f_orthogonal_x2,&
       sll_f_orthogonal_x3,&
       sll_f_orthogonal_jac11,&
       sll_f_orthogonal_jac12,&
       sll_f_orthogonal_jac13,&
       sll_f_orthogonal_jac21,&
       sll_f_orthogonal_jac22,&
       sll_f_orthogonal_jac23,&
       sll_f_orthogonal_jac31,&
       sll_f_orthogonal_jac32,&
       sll_f_orthogonal_jac33,&
       sll_f_orthogonal_jacobian)
  ! Define computational domain in physical coordinates
  eta1_min = 0.0_f64; eta1_max = params(1)
  !eta1_min=map%get_x1(eta1_min,0.0_f64, 0.0_f64)
  !eta1_max=map%get_x1(eta1_max, 0.0_f64, 0.0_f64)

  nc_eta1 = 128
  Lx = eta1_max-eta1_min
  delta_eta1 = 1._f64/real(nc_eta1,f64)
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
  call maxwell_clamped_1d%init( [eta1_min, eta1_min+Lx],nc_eta1, deg, boundary, map, 1d-12, 1d-12)

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

  ex = 0._f64
  rho = 0._f64
  ! Test Poisson
  !-------------
  call maxwell_clamped_1d%compute_rhs_from_function( rho_k, deg, rho)
  !rho(1)=rho(1)-(cos_l(0._f64))
  !rho(nc_eta1+deg)=rho(nc_eta1+deg)+(cos_l(1._f64))
  call maxwell_clamped_1d%compute_e_from_rho( rho, ex )
  ! Evaluate spline curve at grid points and compute error
  ! Ex is a 1-form, i.e. one spline degree lower
  call sll_s_spline_evaluate_basis_b_form_1d_clamped( maxwell_clamped_1d%spline1_pp, nc_eta1,  ex, sval)
  ! Set exact solution
  do i = 1, nc_eta1+1
     factor=map%jacobian_matrix_inverse_transposed( [x(i), 0._f64, 0._f64] )
     sval(i)=sval(i)*factor(1,1)
     ex_exact(i) = E_l(x(i))
  end do
  err_ex = maxval(abs(sval-ex_exact))
  print*, 'Error Poisson',  err_ex

  !call sll_s_plot_two_fields_1d('poisson',nc_eta1+1,sval,ex_exact,0,0.0_f64)

  ! Test Ampere
  !-------------
  ! Set time step
  dt = .5 * delta_eta1
  rho = 0._f64
  call maxwell_clamped_1d%compute_rhs_from_function(sin_k, deg-1, current)
  ex = 0.0_f64
  call maxwell_clamped_1d%compute_E_from_j(dt*current, 1, ex )
  ! Evaluate spline curve at grid points and compute error
  ! Ex is a 1-form, i.e. one spline degree lower
  call sll_s_spline_evaluate_basis_b_form_1d_clamped( maxwell_clamped_1d%spline1_pp, nc_eta1,  ex, sval)
  ! Set exact solution
  do i = 1, nc_eta1+1
     factor=map%jacobian_matrix_inverse_transposed( [x(i), 0._f64, 0._f64] )
     sval(i)=sval(i)*factor(1,1)
     ex_exact(i) = -dt*sin_l(x(i))
  end do
  err_ex2 = maxval(abs(sval-ex_exact))
  print*, 'error Ampere',  err_ex2
  !call sll_s_plot_two_fields_1d('current',nc_eta1+1,sval,ex_exact,0,0.0_f64)

  l2norm =  maxwell_clamped_1d%l2norm_squared(ex,deg-1)
  err_l2norm = abs(l2norm - dt**2*sll_p_pi)
  print*, 'error l2 norm', err_l2norm

!!$  ! Test mixed mass
!!$  ex = 0.0_f64
!!$  ey = 0.0_f64
!!$  call maxwell_clamped_1d%L2projection( cos_k, deg-1, ex )
!!$  call maxwell_clamped_1d%L2projection( cos_k, deg, ey )
!!$  err_inner =  abs(maxwell_clamped_1d%inner_product( ex, ey, deg-1, deg ) - sll_p_pi)
!!$  print*, 'error inner product', err_inner
  err_inner = 0._f64

  ! Test Maxwell_Clamped on By and Ez 
  !--------------------------
  ! Set time stepping parameters
  time  = 0.0_f64
  dt = .5 * delta_eta1
  nstep = 10

  ! Compute initial fields 
  ex = 0.0_f64 ! 0-form -> splines of degree deg
  call maxwell_clamped_1d%L2projection(cos_k, deg-1, bz) 
  ey = 0.0_f64 

  !d/dt E_y=-D_x B_z, d/dt B_z= -D_x E_y
  ! Time loop. Second order Strang splitting
  do istep = 1, nstep
!!$    call maxwell_clamped_1d%compute_b_from_e( 0.5_f64*dt, ey, bz)
!!$    call maxwell_clamped_1d%compute_e_from_b(         dt, bz, ey)
!!$    call maxwell_clamped_1d%compute_b_from_e( 0.5_f64*dt, ey, bz)
     call maxwell_clamped_1d%compute_curl_part( dt, ey, bz )


     call sll_s_spline_evaluate_basis_b_form_1d_clamped( maxwell_clamped_1d%spline0_pp, nc_eta1,  ey, sval)
     time = real(istep,f64)*dt
     do i = 1, nc_eta1+1
        factor = map%jacobian_matrix_inverse_transposed([x(i),0._f64, 0._f64])
        sval(i) = factor(2,2)*sval(i)
        ey_exact(i) =   sin_l(x(i)) * sin(k*time)
     end do
     err_ey = maxval(abs(sval-ey_exact))
     !call sll_s_plot_two_fields_1d('ey',nc_eta1+1,sval,ey_exact,istep,time)
     call sll_s_spline_evaluate_basis_b_form_1d_clamped( maxwell_clamped_1d%spline1_pp, nc_eta1,  bz, sval)

     do i = 1, nc_eta1+1
        factor = map%jacobian_matrix([x(i),0._f64, 0._f64])/map%jacobian( [x(i), 0._f64, 0._f64])
        sval(i) = sval(i)*factor(3,3)
        bz_exact(i) =   cos_l(x(i)) *  cos(k*time)
     end do
     err_bz = maxval(abs(sval-bz_exact))
     !call sll_s_plot_two_fields_1d('bz',nc_eta1+1,sval,bz_exact,istep,time)
     write(*,"(10x,' istep = ',I6)",advance="no") istep
     write(*,"(' time = ',g12.3,' sec')",advance="no") time
     write(*,"(' error L2 = ',2g15.5)") err_ey, err_bz

  end do

  tol = 1.0d-3

  if ((err_bz < tol) .and. (err_ey < tol) .and. (err_ex < tol) .and. (err_ex2 < 1.0d-6) .and. (err_l2norm < 1.0d-5).and. (err_inner < tol)) then
     print*,'PASSED'
  endif

  call maxwell_clamped_1d%free()
  map=>null()

  DEALLOCATE(ex)
  DEALLOCATE(ey)
  DEALLOCATE(bz)
  DEALLOCATE(bz_exact)
  DEALLOCATE(ex_exact)
  DEALLOCATE(ey_exact)
  DEALLOCATE(rho)
  deallocate(x)
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

  function E_l(eta)
    sll_real64             :: E_l
    sll_real64, intent(in) :: eta
    sll_real64             :: x

    x=map%get_x1( [eta, 0._f64, 0._f64] )
    E_l = alpha/k*cos(k*x) 
  end function E_l


  function cos_k(x1)
    sll_real64             :: cos_k
    sll_real64, intent(in) :: x1

    cos_k = cos(k*x1)
  end function cos_k

  function sin_k(x1)
    sll_real64             :: sin_k
    sll_real64, intent(in) :: x1

    sin_k = sin(k*x1)
  end function sin_k


  function cos_l(eta)
    sll_real64             :: cos_l
    sll_real64, intent(in) :: eta
    !local variable
    sll_real64             :: x1

    x1=map%get_x1( [eta, 0._f64, 0._f64] )
    cos_l = cos(k*x1)
  end function cos_l


  function sin_l(eta)
    sll_real64             :: sin_l
    sll_real64, intent(in) :: eta
    !local variable
    sll_real64             :: x1

    x1=map%get_x1( [eta, 0._f64, 0._f64] )
    sin_l = sin(k*x1)
  end function sin_l

  function constant(x)
    sll_real64             :: constant
    sll_real64, intent(in) :: x

    constant=1._f64
  end function constant

  function poly1_k(x)
    sll_real64             :: poly1_k
    sll_real64, intent(in) :: x

    poly1_k= x
  end function poly1_k

  function poly1_l(eta)
    sll_real64             :: poly1_l
    sll_real64, intent(in) :: eta
    !local variable
    sll_real64             :: x1

    x1=map%get_x1( [eta, 0._f64, 0._f64] )

    poly1_l= x1
  end function poly1_l

  function poly2_k(x)
    sll_real64             :: poly2_k
    sll_real64, intent(in) :: x

    poly2_k= x**2
  end function poly2_k

  function poly2_l(eta)
    sll_real64             :: poly2_l
    sll_real64, intent(in) :: eta
    !local variable
    sll_real64             :: x1

    x1=map%get_x1( [eta, 0._f64, 0._f64] )

    poly2_l= x1**2
  end function poly2_l

  function poly3_k(x)
    sll_real64             :: poly3_k
    sll_real64, intent(in) :: x

    poly3_k= x**3
  end function poly3_k

  function poly3_l(eta)
    sll_real64             :: poly3_l
    sll_real64, intent(in) :: eta
    !local variable
    sll_real64             :: x1

    x1=map%get_x1( [eta, 0._f64, 0._f64] )

    poly3_l= x1**3
  end function poly3_l

  function poly4_k(x)
    sll_real64             :: poly4_k
    sll_real64, intent(in) :: x

    poly4_k= x**4
  end function poly4_k

  function poly4_l(eta)
    sll_real64             :: poly4_l
    sll_real64, intent(in) :: eta
    !local variable
    sll_real64             :: x1

    x1=map%get_x1( [eta, 0._f64, 0._f64] )

    poly4_l= x1**4
  end function poly4_l


  function zero(x)
    sll_real64             :: zero
    sll_real64, intent(in) :: x

    zero=0._f64
  end function zero

end program test_maxwell_clamped_1d_trafo
