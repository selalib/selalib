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
program test_curl_curl_part
  !------------------------------------------------------------------------
  !  test 3D Maxwell spline finite element solver on a periodic grid
  !------------------------------------------------------------------------
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_low_level_bsplines, only: &
       sll_s_eval_uniform_periodic_spline_curve_with_zero

  use sll_m_constants, only: &
       sll_p_pi, sll_p_twopi

  use sll_m_maxwell_3d_fem, only: &
       sll_t_maxwell_3d_fem

  use sll_m_maxwell_clamped_3d_fem

  use sll_m_maxwell_1d_base, only: &
       sll_s_plot_two_fields_1d

  use sll_m_constants, only: sll_p_pi, sll_p_twopi

  use sll_m_splines_pp

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark

  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !  type(sll_t_arbitrary_degree_spline_1d), pointer :: aspl
  !  sll_real64, dimension(:), allocatable :: knots

  sll_real64 :: eta1_max, eta1_min
  sll_real64 :: delta_eta(3)

  sll_int32  :: nc_eta(3), nc_total, nc_total0, nc_total1

  type(sll_t_maxwell_3d_fem)  :: maxwell_3d
  !type(sll_t_maxwell_clamped_3d_fem)  :: maxwell_3d
  sll_real64, allocatable :: afield(:)
  sll_real64, allocatable :: afield_ref(:)
  sll_real64, allocatable :: afield_val(:), scratch(:)
  sll_real64, allocatable :: rho(:), rho_ref(:), current(:), noise(:)

  sll_real64                              :: eps
  sll_real64, dimension(3,2)              :: domain
  sll_int32                               :: deg(3), boundary(3)

  type(sll_t_time_mark) :: start, end

  call sll_s_set_time_mark( start )

  ! Define computational domain
  eta1_min = .0_f64; eta1_max = 2.0_f64*sll_p_pi
  nc_eta = 16![8, 16, 32]
  nc_total = product(nc_eta)
  domain(1,:) = [eta1_min, eta1_max]
  domain(2,:) = [eta1_min, eta1_max]
  domain(3,:) = [eta1_min, eta1_max]
  ! Set spline degree of 0-forms
  deg = 3![2,2,3]


  if( deg(1) == 2 ) then
     boundary = [ sll_p_boundary_clamped_square, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else if( deg(1) == 3 ) then
     boundary = [ sll_p_boundary_clamped_cubic, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else
     boundary = [ sll_p_boundary_clamped, sll_p_boundary_periodic, sll_p_boundary_periodic]
  end if

  ! Initialise maxwell FEM object
  call maxwell_3d%init(domain, nc_eta, deg)
  !call maxwell_3d%init(domain, nc_eta, deg, boundary)

  nc_total0 = maxwell_3d%n_total0
  nc_total1 = maxwell_3d%n_total1

  allocate( afield(1:nc_total1+nc_total0*2) )
  allocate( afield_ref(1:nc_total1+nc_total0*2) )
  allocate( afield_val(1:nc_total1+nc_total0*2) )
  allocate( scratch(1:nc_total1+nc_total0*2) )

  allocate( rho( nc_total0 ) )
  allocate( rho_ref( nc_total0 ) )
  allocate( current(1:nc_total1+nc_total0*2) )
  allocate( noise(1:nc_total1+nc_total0*2) )


  current = 0._f64
!!$  call maxwell_3d%compute_rhs_from_function( 1, 1, current(1:nc_total1), j1x )
!!$  call maxwell_3d%compute_rhs_from_function( 1, 2, current(nc_total1+1:nc_total1+nc_total0), j1y )
!!$  call maxwell_3d%compute_rhs_from_function( 1, 3, current(nc_total1+nc_total0+1:nc_total1+nc_total0*2), j1z )
!!$  
  call maxwell_3d%compute_rhs_from_function( 1, 1, current(1:nc_total1), j2x )
  call maxwell_3d%compute_rhs_from_function( 1, 2, current(nc_total1+1:nc_total1+nc_total0), j2y )
  call maxwell_3d%compute_rhs_from_function( 1, 3, current(nc_total1+nc_total0+1:nc_total1+nc_total0*2), j2z )
  call random_number( afield )
  call maxwell_3d%multiply_ct( afield, noise )
!!$  call random_number( noise )
  eps = 1d-8 !0._f64!
  current = (1._f64-eps) * current +  eps* noise
!!$  rho = 0._f64
  call maxwell_3d%multiply_gt( current, rho )
  print*, 'rhs divergence', maxval(abs(rho))

 ! maxwell_3d%curl_matrix%epsilon = 6.0_f64!1d-0
 ! call maxwell_3d%curl_solver%solve( current, afield )
  call maxwell_3d%preconditioner_curl_fft%solve( current, afield )
  !call maxwell_3d%uzawa_iterator%solve( current, afield )

  !print*, 'afield', afield

  !print*, 'p',maxwell_3d%uzawa_iterator%x_0

  afield_ref = 0._f64
!!$  call maxwell_3d%L2projection( 1, 1, afield_ref(1:nc_total1), a1x )
!!$  call maxwell_3d%L2projection( 1, 2, afield_ref(nc_total1+1:nc_total1+nc_total0), a1y )
!!$  call maxwell_3d%L2projection( 1, 3, afield_ref(nc_total1+nc_total0+1:nc_total1+nc_total0*2), a1z )

  call maxwell_3d%L2projection( 1, 1, afield_ref(1:nc_total1), a2x )
  call maxwell_3d%L2projection( 1, 2, afield_ref(nc_total1+1:nc_total1+nc_total0), a2y )
  call maxwell_3d%L2projection( 1, 3, afield_ref(nc_total1+nc_total0+1:nc_total1+nc_total0*2), a2z )
!!$  
  rho = 0._f64
  call maxwell_3d%compute_rho_from_E( afield_ref, rho )
!!$  print*, 'lhs divergence', maxval(abs(rho))
!!$
  rho_ref = 0._f64
  !call maxwell_3d%L2projection( 0, 1, rho_ref, cos_k )
  maxwell_3d%curl_matrix%epsilon = 0._f64
  call maxwell_3d%curl_matrix%dot(afield, afield_val )
  call maxwell_3d%MG_operator%dot(maxwell_3d%uzawa_iterator%x_0, scratch )

  afield_val = afield_val + scratch
  
  print*, 'error operator', maxval(abs(afield_val(1:nc_total1) - current(1:nc_total1) ) ), maxval(abs(afield_val(nc_total1+1:nc_total1+nc_total0) - current(nc_total1+1:nc_total1+nc_total0) ) ), maxval(abs(afield_val(nc_total1+nc_total0+1:nc_total0+nc_total*2) - current(nc_total1+nc_total0+1:nc_total0+nc_total*2) ) )



  print*, 'error solver', maxval(abs(afield(1:nc_total1) -afield_ref(1:nc_total1) ) ), maxval(abs(afield(nc_total1+1:nc_total1+nc_total0) -afield_ref(nc_total1+1:nc_total1+nc_total0) ) ), maxval(abs(afield(nc_total1+nc_total0+1:nc_total0+nc_total*2) -afield_ref(nc_total1+nc_total0+1:nc_total0+nc_total*2) ) ), maxval(abs(rho_ref- maxwell_3d%uzawa_iterator%x_0))


!!$  call sll_s_plot_two_fields_1d('currentx',nc_total,afield_val(1:nc_total),current(1:nc_total),0._f64,0._f64)
!!$  call sll_s_plot_two_fields_1d('currenty',nc_total,afield_val(nc_total+1:nc_total*2),current(nc_total+1:nc_total*2),0._f64,0._f64)
!!$  call sll_s_plot_two_fields_1d('currentz',nc_total,afield_val(nc_total*2+1:nc_total*3),current(nc_total*2+1:nc_total*3),0._f64,0._f64)
  
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

    j1x = 3._f64*cos(x(1)+x(2)+x(3))  !-sin(x(1)+x(2)+x(3))
  end function j1x

  function j1y(x)
    sll_real64             :: j1y
    sll_real64, intent(in) :: x(3)

    j1y = 0._f64!-sin(x(1)+x(2)+x(3))
  end function j1y

  function j1z(x)
    sll_real64             :: j1z
    sll_real64, intent(in) :: x(3)

    j1z = -3._f64*cos(x(1)+x(2)+x(3))  !-sin(x(1)+x(2)+x(3))
  end function j1z

  function j2x(x)
    sll_real64             :: j2x
    sll_real64, intent(in) :: x(3)

    j2x = 3._f64*( sin(x(1)+x(2)+x(3)) + cos(x(1)+x(2)+x(3)) ) !-sin(x(1)+x(2)+x(3))
  end function j2x

  function j2y(x)
    sll_real64             :: j2y
    sll_real64, intent(in) :: x(3)

    j2y = -3._f64*( cos(x(1)+x(2)+x(3)) - 4._f64* cos(2._f64*(x(1)+x(2)+x(3)) ) )  !-sin(x(1)+x(2)+x(3))
  end function j2y

  function j2z(x)
    sll_real64             :: j2z
    sll_real64, intent(in) :: x(3)

    j2z = -3._f64*( sin(x(1)+x(2)+x(3)) + 4._f64* cos(2._f64*(x(1)+x(2)+x(3)) ) )  !-sin(x(1)+x(2)+x(3))
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


end program test_curl_curl_part
