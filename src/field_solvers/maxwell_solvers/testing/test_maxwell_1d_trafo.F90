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
program test_maxwell_1d_trafo
  !------------------------------------------------------------------------
  !  test 1D Maxwell spline finite element solver on a periodic grid
  !------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_low_level_bsplines, only: &
    sll_s_eval_uniform_periodic_spline_curve

  use sll_m_constants, only: &
    sll_p_pi, sll_p_twopi

  use sll_m_maxwell_1d_base, only: &
    sll_s_plot_two_fields_1d

  use sll_m_maxwell_1d_trafo, only: &
       sll_t_maxwell_1d_trafo

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

  type(sll_t_maxwell_1d_trafo)    :: maxwell_1d
  type(sll_t_mapping_3d), pointer :: map  


  sll_real64, dimension(:), allocatable :: ex
  sll_real64, dimension(:), allocatable :: ey
  sll_real64, dimension(:), allocatable :: bz
  sll_real64, dimension(:), allocatable :: bz_exact
  sll_real64, dimension(:), allocatable :: ex_exact
  sll_real64, dimension(:), allocatable :: ey_exact
  sll_real64, dimension(:), allocatable :: rho
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
  sll_real64                              :: factor(3,3)
  sll_real64                              :: dt
  sll_real64                              :: Lx
  sll_int32                               :: deg
  sll_real64                              :: l2norm
 
  sll_real64 :: params(6)
  !Initialize mapping
  params = 0._f64
  params(1)=sll_p_twopi
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
  eta1_min = 0.0_f64; eta1_max = 1.0_f64
  !eta1_min=map%get_x1(eta1_min,0.0_f64, 0.0_f64)
  !eta1_max=map%get_x1(eta1_max, 0.0_f64, 0.0_f64)

  nc_eta1 = 256
  Lx = eta1_max-eta1_min
  delta_eta1 = 1._f64/real(nc_eta1,f64)
  ! Set spline degree of 0-forms
  deg = 3
  ! Initialise maxwell FEM object
  call maxwell_1d%init( [eta1_min, eta1_min+params(1)],nc_eta1, deg, map, 1d-12, 1d-12)

  allocate( x(1:nc_eta1+1) )

  do i = 1, nc_eta1+1
    x(i) = eta1_min + real(i-1,f64)*delta_eta1
  end do

 
  ! Allocate arrays
  SLL_CLEAR_ALLOCATE(ex(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(ey(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(bz(1:nc_eta1),error)

  SLL_CLEAR_ALLOCATE(bz_exact(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(ex_exact(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(ey_exact(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(rho(1:nc_eta1),error)

  SLL_CLEAR_ALLOCATE(sval(1:nc_eta1),error)

  call maxwell_1d%L2projection( sin_k, deg-1, ex )
  call sll_s_eval_uniform_periodic_spline_curve(deg-1, ex, sval)
  ! Set exact solution
  do i = 1, nc_eta1
    factor=map%jacobian_matrix_inverse_transposed( [x(i), 0._f64, 0._f64] )
    sval(i)=sval(i)*factor(1,1)
    ex_exact(i) = sin_l(x(i))
  end do
  err_ex = maxval(abs(sval-ex_exact))
  print*, 'Error L2projectionE1',  err_ex

  ! call sll_s_plot_two_fields_1d('ex',nc_eta1,sval,ex_exact,0,0.0_f64)

  call maxwell_1d%L2projection( sin_k, deg, ey )
  call sll_s_eval_uniform_periodic_spline_curve(deg, ey, sval)
  ! Set exact solution
  do i = 1, nc_eta1
    factor=map%jacobian_matrix_inverse_transposed( [x(i), 0._f64, 0._f64] )
    sval(i)=sval(i)*factor(2,2)
    ey_exact(i) = sin_l(x(i))
  end do
  err_ey = maxval(abs(sval-ey_exact))
  print*, 'Error L2projectionE2',  err_ey

  ! call sll_s_plot_two_fields_1d('ey',nc_eta1,sval,ey_exact,0,0.0_f64)

  call maxwell_1d%L2projection( sin_k, deg-1, bz )
  call sll_s_eval_uniform_periodic_spline_curve(deg-1, bz, sval)
  ! Set exact solution
  do i = 1, nc_eta1
    factor=map%jacobian_matrix( [x(i), 0._f64, 0._f64] )/map%jacobian( [x(i), 0._f64, 0._f64])
    sval(i)=sval(i)*factor(3,3)
    bz_exact(i) = sin_l(x(i))
  end do
  err_bz = maxval(abs(sval-bz_exact))
  print*, 'Error L2projectionB3',  err_bz

  !call sll_s_plot_two_fields_1d('bz',nc_eta1,sval,bz_exact,0,0.0_f64)

  ! Test Poisson
  !-------------
   call maxwell_1d%compute_rhs_from_function( cos_k, deg, rho)
   call maxwell_1d%compute_E_from_rho( rho, ex )
   ! Evaluate spline curve at grid points and compute error
  ! Ex is a 1-form, i.e. one spline degree lower
   call sll_s_eval_uniform_periodic_spline_curve(deg-1, ex, sval)
   ! Set exact solution
   do i = 1, nc_eta1
     factor=map%jacobian_matrix_inverse_transposed( [x(i), 0._f64, 0._f64] )
     sval(i)=sval(i)*factor(1,1)
     ex_exact(i) = sin_l(x(i))
   end do
   err_ex = maxval(abs(sval-ex_exact))
  print*, 'Error Poisson',  err_ex

  !call sll_s_plot_two_fields_1d('poisson',nc_eta1,sval,ex_exact,0,0.0_f64)
 
  ! Test Ampere
  !-------------
  ! Set time step
  dt = .5 * delta_eta1

  call maxwell_1d%compute_rhs_from_function(cos_k, deg-1, rho)
  ex = 0.0_f64
  call maxwell_1d%compute_E_from_j(dt*rho, 1, ex )
  ! Evaluate spline curve at grid points and compute error
  ! Ex is a 1-form, i.e. one spline degree lower
  call sll_s_eval_uniform_periodic_spline_curve(deg-1, ex, sval)
   ! Set exact solution
   do i = 1, nc_eta1
     factor=map%jacobian_matrix_inverse_transposed( [x(i), 0._f64, 0._f64] )
     sval(i)=sval(i)*factor(1,1)
     ex_exact(i) = -dt*cos_l(x(i))
   end do
  err_ex2 = maxval(abs(sval-ex_exact))
  print*, 'error Ampere',  err_ex2
  !call sll_plot_two_fields_1d('ex',nc_eta1,sval,ex_exact,0,0.0_f64)

  l2norm =  maxwell_1d%l2norm_squared(ex,deg-1)
  err_l2norm = abs(l2norm - dt**2*sll_p_pi)
  print*, 'error l2 norm', err_l2norm

  ! Test mixed mass
  ex = 0.0_f64
  ey = 0.0_f64
  call maxwell_1d%L2projection( cos_k, deg-1, ex )
  call maxwell_1d%L2projection( cos_k, deg, ey )
  err_inner =  abs(maxwell_1d%inner_product( ex, ey, deg-1, deg ) - sll_p_pi)
  print*, 'error inner product', err_inner


  ! Test Maxwell on By and Ez 
  !--------------------------
  ! Set time stepping parameters
  time  = 0.0_f64
  dt = .5 * delta_eta1
  nstep = 10

  ! Compute initial fields 
  ex = 0.0_f64 ! 0-form -> splines of degree deg
  call maxwell_1d%L2projection(cos_k, deg-1, bz) 
  ey = 0.0_f64 

  !d/dt E_y=-D_x B_z, d/dt B_z= -D_x E_y
  ! Time loop. Second order Strang splitting
  do istep = 1, nstep
!!$    call maxwell_1d%compute_b_from_e( 0.5_f64*dt, ey, bz)
!!$    call maxwell_1d%compute_e_from_b(         dt, bz, ey)
!!$    call maxwell_1d%compute_b_from_e( 0.5_f64*dt, ey, bz)
    call maxwell_1d%compute_curl_part( dt, ey, bz )
   
  call sll_s_eval_uniform_periodic_spline_curve(deg, ey, sval)
   time = real(istep,f64)*dt
   do i = 1, nc_eta1
     factor = map%jacobian_matrix_inverse_transposed( [x(i), 0._f64, 0._f64] )
     sval(i) = factor(2,2)*sval(i)
     ey_exact(i) =   sin_l(x(i)) * sin(time)

   end do
   err_ey = maxval(abs(sval-ey_exact))

   
  call sll_s_eval_uniform_periodic_spline_curve(deg-1, bz, sval)
      
   do i = 1, nc_eta1
     factor = map%jacobian_matrix( [x(i), 0._f64, 0._f64])/map%jacobian( [x(i), 0._f64, 0._f64])
     sval(i) = sval(i)*factor(3,3)
     bz_exact(i) =   cos_l(x(i)) *  cos(time)
   end do
   err_bz = maxval(abs(sval-bz_exact))
   
   !call sll_s_plot_two_fields_1d('bz',nc_eta1,sval,bz_exact,istep,time)


     write(*,"(10x,' istep = ',I6)",advance="no") istep
     write(*,"(' time = ',g12.3,' sec')",advance="no") time
     write(*,"(' error L2 = ',2g15.5)") err_ey, err_bz

    

   end do

  tol = 1.0d-3

  if ((err_bz < tol) .and. (err_ey < tol) .and. (err_ex < 1.0d-4) .and. (err_ex2 < 1.0d-6) .and. (err_l2norm < 1.0d-9).and. (err_inner < 1.0d-3)) then
     print*,'PASSED'
  endif

  call maxwell_1d%free()
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
    function cos_k(x1)
    sll_real64             :: cos_k
    sll_real64, intent(in) :: x1

    cos_k = cos(x1)
  end function cos_k
  
  function sin_k(x1)
    sll_real64             :: sin_k
    sll_real64, intent(in) :: x1
  
    sin_k = sin(x1)
  end function sin_k


  function cos_l(eta)
    sll_real64             :: cos_l
    sll_real64, intent(in) :: eta
    !local variable
    sll_real64             :: x1

    x1=map%get_x1( [eta, 0._f64, 0._f64] )
    cos_l = cos(x1)
  end function cos_l


  function sin_l(eta)
    sll_real64             :: sin_l
    sll_real64, intent(in) :: eta
    !local variable
    sll_real64             :: x1

    x1=map%get_x1( [eta, 0._f64, 0._f64] )
    sin_l = sin(x1)
  end function sin_l


  function zero(x)
    sll_real64             :: zero
    sll_real64, intent(in) :: x

    zero=0._f64
  end function zero

end program test_maxwell_1d_trafo
