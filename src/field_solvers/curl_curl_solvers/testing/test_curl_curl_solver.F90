!
!  Contact : Benedikt Perse
!
program test_curl_curl_solver
  !------------------------------------------------------------------------
  !  test 3D Curl_curl finite element solver on a periodic grid
  !------------------------------------------------------------------------
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_low_level_bsplines, only: &
       sll_s_eval_uniform_periodic_spline_curve_with_zero

  use sll_m_constants, only: &
       sll_p_pi, sll_p_twopi

  use sll_m_maxwell_1d_base, only: &
       sll_s_plot_two_fields_1d

  use sll_m_constants, only: sll_p_pi, sll_p_twopi

  use sll_m_splines_pp

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark
  
  use sll_m_preconditioner_curl_solver_fft, only : &
       sll_t_preconditioner_curl_solver_fft
  
  use sll_m_preconditioner_jacobi, only : &
       sll_t_preconditioner_jacobi

  use sll_m_maxwell_3d_fem, only : &
       sll_t_maxwell_3d_fem
  
  use sll_m_linear_operator_curl_3d

  use sll_m_linear_operator_GTM

  use sll_m_linear_operator_MG

  use sll_m_linear_operator_penalized, only : &
       sll_t_linear_operator_penalized

  use sll_m_linear_solver_cg, only : &
       sll_t_linear_solver_cg

  use sll_m_spline_fem_utilities, only : &
       sll_s_spline_fem_mass_line, &
       sll_s_spline_fem_mixedmass_line, &
       sll_s_multiply_g, &
       sll_s_multiply_gt, &
       sll_s_multiply_c, &
       sll_s_multiply_ct
  
  use sll_m_uzawa_iterator

  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_real64 :: x_max, x_min
  sll_real64 :: delta_x(3)
  type(sll_t_maxwell_3d_fem)  :: maxwell_3d
  sll_int32  :: n_dofs(3), n_total, n_total0, n_total1, j

  type(sll_t_linear_operator_curl_3d) :: curl_matrix  !< curl matrix
  type(sll_t_linear_operator_penalized)  :: curl_operator !< curl matrix with constraint on constant vector
  type( sll_t_preconditioner_curl_solver_fft ) :: preconditioner_curl_fft
  type( sll_t_preconditioner_jacobi ) :: preconditioner_curl_jacobi
  type(sll_t_linear_solver_cg)  :: curl_solver !< CG solver to invert curl matrix

  type(sll_t_linear_operator_MG) :: MG_operator

  type(sll_t_linear_operator_GTM) :: GTM_operator

  type(sll_t_uzawa_iterator) :: uzawa_iterator
  
  sll_real64, allocatable :: afield(:)
  sll_real64, allocatable :: afield_ref(:)
  sll_real64, allocatable :: afield_val(:), scratch(:)
  sll_real64, allocatable :: rho(:), rho_ref(:), current(:), noise(:)
  sll_real64, allocatable :: nullspace(:,:)

  sll_real64                              :: eps
  sll_real64, dimension(3,2)              :: domain
  sll_real64                              :: error(4)
  sll_int32                               :: deg(3)

  type(sll_t_time_mark) :: start, end

  call sll_s_set_time_mark( start )

  ! Define computational domain
  x_min = .0_f64; x_max = 2.0_f64*sll_p_pi
  n_dofs = 16![16, 16, 8]
  n_total = product(n_dofs)
  domain(1,:) = [x_min, x_max]
  domain(2,:) = [x_min, x_max]
  domain(3,:) = [x_min, x_max]
  ! Set spline degree of 0-forms
  deg = 3![3,3,2]

  delta_x = (domain(:,2)-domain(:,1))/real(n_dofs,f64)

  ! Initialise maxwell FEM object
  call maxwell_3d%init(domain, n_dofs, deg)
  
  n_total0 = maxwell_3d%n_total0
  n_total1 = maxwell_3d%n_total1

  allocate(nullspace(1,1:3*n_total0))
    nullspace(1,:) = 1.0_f64
  call preconditioner_curl_fft%create( n_dofs, delta_x, deg )

  call curl_matrix%create( maxwell_3d%mass1_operator, maxwell_3d%mass2_operator, n_dofs, delta_x   )
  call curl_operator%create( linear_operator=curl_matrix, vecs=nullspace, n_dim_nullspace=1 )
  call preconditioner_curl_jacobi%create( curl_matrix )

  call curl_solver%create( curl_operator, preconditioner_curl_fft )
  curl_solver%null_space = .true.
  curl_solver%atol = 1.0d-12
  !curl_solver%verbose = .true.
  curl_solver%n_maxiter=1000

  call MG_operator%create( maxwell_3d%mass1_operator, n_dofs, delta_x )
  call GTM_operator%create( maxwell_3d%mass1_operator, n_dofs, delta_x )
  call uzawa_iterator%create( curl_solver, MG_operator, GTM_operator )
  !uzawa_iterator%verbose = .true.
  uzawa_iterator%atol = 1.0d-10
  uzawa_iterator%n_maxiter=100


  allocate( afield(1:n_total1+n_total0*2) )
  allocate( afield_ref(1:n_total1+n_total0*2) )
  allocate( afield_val(1:n_total1+n_total0*2) )
  allocate( scratch(1:n_total1+n_total0*2) )

  allocate( rho( n_total0 ) )
  allocate( rho_ref( n_total0 ) )
  allocate( current(1:n_total1+n_total0*2) )
  allocate( noise(1:n_total1+n_total0*2) )


  current = 0._f64
!!$  call maxwell_3d%compute_rhs_from_function( 1, 1, current(1:n_total1), j1x )
!!$  call maxwell_3d%compute_rhs_from_function( 1, 2, current(n_total1+1:n_total1+n_total0), j1y )
!!$  call maxwell_3d%compute_rhs_from_function( 1, 3, current(n_total1+n_total0+1:n_total1+n_total0*2), j1z )
!!$  
  call maxwell_3d%compute_rhs_from_function( 1, 1, current(1:n_total1), j2x )
  call maxwell_3d%compute_rhs_from_function( 1, 2, current(n_total1+1:n_total1+n_total0), j2y )
  call maxwell_3d%compute_rhs_from_function( 1, 3, current(n_total1+n_total0+1:n_total1+n_total0*2), j2z )
!!$  call random_number( afield )
!!$   call sll_s_multiply_ct(n_dofs, delta_x, afield, noise )
  call random_number( noise )
  eps = 1d-8 !0._f64!
  current = (1._f64-eps) * current +  eps* noise
!!$  rho = 0._f64
  call sll_s_multiply_gt(n_dofs, delta_x, current, rho )
  print*, 'rhs divergence', maxval(abs(rho))

  !curl_matrix%epsilon = 2.0_f64!1d-0
  call uzawa_iterator%solve( current, afield )

  afield_ref = 0._f64
!!$  call maxwell_3d%L2projection( 1, 1, afield_ref(1:n_total1), a1x )
!!$  call maxwell_3d%L2projection( 1, 2, afield_ref(n_total1+1:n_total1+n_total0), a1y )
!!$  call maxwell_3d%L2projection( 1, 3, afield_ref(n_total1+n_total0+1:n_total1+n_total0*2), a1z )

  call maxwell_3d%L2projection( 1, 1, afield_ref(1:n_total1), a2x )
  call maxwell_3d%L2projection( 1, 2, afield_ref(n_total1+1:n_total1+n_total0), a2y )
  call maxwell_3d%L2projection( 1, 3, afield_ref(n_total1+n_total0+1:n_total1+n_total0*2), a2z )
!!$  
  rho = 0._f64
  call maxwell_3d%compute_rho_from_E( afield_ref, rho )
!!$  print*, 'lhs divergence', maxval(abs(rho))
!!$
  rho_ref = 0._f64
  call maxwell_3d%L2projection( 0, 1, rho_ref, cos_k )
  curl_matrix%epsilon = 0._f64
  call curl_matrix%dot(afield, afield_val )
  call MG_operator%dot(uzawa_iterator%x_0, scratch )

  afield_val = afield_val + scratch
  
  print*, 'error operator', maxval(abs(afield_val(1:n_total1) - current(1:n_total1) ) ), maxval(abs(afield_val(n_total1+1:n_total1+n_total0) - current(n_total1+1:n_total1+n_total0) ) ), maxval(abs(afield_val(n_total1+n_total0+1:n_total0+n_total*2) - current(n_total1+n_total0+1:n_total0+n_total*2) ) )



  error(1) = maxval(abs(afield(1:n_total1) -afield_ref(1:n_total1) ) )
  error(2) = maxval(abs(afield(n_total1+1:n_total1+n_total0) -afield_ref(n_total1+1:n_total1+n_total0) ) )
  error(3) = maxval(abs(afield(n_total1+n_total0+1:n_total0+n_total*2) -afield_ref(n_total1+n_total0+1:n_total0+n_total*2) ) )
  error(4) = maxval(abs(rho_ref- uzawa_iterator%x_0))


  print*, 'error solver', error

!!$  call sll_s_plot_two_fields_1d('currentx',n_total,afield_val(1:n_total),current(1:n_total),0._f64,0._f64)
!!$  call sll_s_plot_two_fields_1d('currenty',n_total,afield_val(n_total+1:n_total*2),current(n_total+1:n_total*2),0._f64,0._f64)
!!$  call sll_s_plot_two_fields_1d('currentz',n_total,afield_val(n_total*2+1:n_total*3),current(n_total*2+1:n_total*3),0._f64,0._f64)
  
  ! Clean up
  call maxwell_3d%free()
  deallocate( afield )
  deallocate( afield_ref )
  deallocate( afield_val )
  deallocate( scratch )
  deallocate( rho )
  deallocate( rho_ref )
  deallocate( current )
  deallocate( noise )

  if ( error(1) < 3.d-4 .AND. error(2) < 8.d-4 .AND. error(3) < 1.d-3 .AND. error(4) < 4.d-4 ) then
     print*, 'PASSED.'
  else
     print*, 'FAILED.'
  end if


  call sll_s_set_time_mark( end )
  write(*, "(A, F10.3)") "Main part run time [s] = ", sll_f_time_elapsed_between( start, end)



contains
  function cos_k(x)
    sll_real64             :: cos_k
    sll_real64, intent(in) :: x(3)

    cos_k = cos((x(1)+x(2)+x(3))) 
  end function cos_k


  function sin_k(x)
    sll_real64             :: sin_k
    sll_real64, intent(in) :: x(3)

    sin_k = sin((x(1)+x(2)+x(3))) 
  end function sin_k

    function j1x(x)
    sll_real64             :: j1x
    sll_real64, intent(in) :: x(3)

    j1x = 3._f64*cos(x(1)+x(2)+x(3)) -sin(x(1)+x(2)+x(3))
  end function j1x

  function j1y(x)
    sll_real64             :: j1y
    sll_real64, intent(in) :: x(3)

    j1y = -sin(x(1)+x(2)+x(3))
  end function j1y

  function j1z(x)
    sll_real64             :: j1z
    sll_real64, intent(in) :: x(3)

    j1z = -3._f64*cos(x(1)+x(2)+x(3))  -sin(x(1)+x(2)+x(3))
  end function j1z

  function j2x(x)
    sll_real64             :: j2x
    sll_real64, intent(in) :: x(3)

    j2x = 3._f64*( sin(x(1)+x(2)+x(3)) + cos(x(1)+x(2)+x(3)) ) -sin(x(1)+x(2)+x(3))
  end function j2x

  function j2y(x)
    sll_real64             :: j2y
    sll_real64, intent(in) :: x(3)

    j2y = -3._f64*( cos(x(1)+x(2)+x(3)) - 4._f64* cos(2._f64*(x(1)+x(2)+x(3)) ) )  -sin(x(1)+x(2)+x(3))
  end function j2y

  function j2z(x)
    sll_real64             :: j2z
    sll_real64, intent(in) :: x(3)

    j2z = -3._f64*( sin(x(1)+x(2)+x(3)) + 4._f64* cos(2._f64*(x(1)+x(2)+x(3)) ) )  -sin(x(1)+x(2)+x(3))
  end function j2z

  function a1x(x)
    sll_real64             :: a1x
    sll_real64, intent(in) :: x(3)

    a1x = cos(x(1)+x(2)+x(3)) 
  end function a1x

  function a1y(x)
    sll_real64             :: a1y
    sll_real64, intent(in) :: x(3)

    a1y = 0._f64
  end function a1y

  function a1z(x)
    sll_real64             :: a1z
    sll_real64, intent(in) :: x(3)

    a1z = -cos(x(1)+x(2)+x(3)) 
  end function a1z

   function a2x(x)
    sll_real64             :: a2x
    sll_real64, intent(in) :: x(3)

    a2x =  sin(x(1)+x(2)+x(3)) + cos(x(1)+x(2)+x(3)) 
  end function a2x

  function a2y(x)
    sll_real64             :: a2y
    sll_real64, intent(in) :: x(3)

    a2y = -sin(x(1)+x(2)+x(3))**2 + cos(x(1)+x(2)+x(3))**2 -cos(x(1)+x(2)+x(3)) 
  end function a2y

  function a2z(x)
    sll_real64             :: a2z
    sll_real64, intent(in) :: x(3)

    a2z = ( sin(x(1)+x(2)+x(3)) -1._f64) * sin(x(1)+x(2)+x(3)) - cos(x(1)+x(2)+x(3))**2 
  end function a2z


end program test_curl_curl_solver
