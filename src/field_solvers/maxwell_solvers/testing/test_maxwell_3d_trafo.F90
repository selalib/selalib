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
!  Contact : Benedikt Perse
!
program test_maxwell_3d_trafo
  !------------------------------------------------------------------------
  !  test 3D Maxwell spline finite element solver on a periodic grid
  !------------------------------------------------------------------------
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_deriv, &
       sll_s_eval_uniform_periodic_spline_curve, &
       sll_s_eval_uniform_periodic_spline_curve_with_zero

  use sll_m_constants, only: &
       sll_p_pi, sll_p_twopi, sll_p_fourpi

  use sll_m_maxwell_3d_trafo, only: &
       sll_t_maxwell_3d_trafo

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_3d_coordinate_transformations

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark

  use sll_m_maxwell_1d_base, only: &
       sll_s_plot_two_fields_1d

  !use omp_lib


  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !  type(sll_t_arbitrary_degree_spline_1d), pointer :: aspl
  !  sll_real64, dimension(:), allocatable :: knots

  sll_real64                      :: eta1_max, eta1_min
  sll_real64                      :: delta_eta(3)
  sll_int32                       :: nc_eta(3), nc_total
  type(sll_t_maxwell_3d_trafo)    :: maxwell_3d
  type(sll_t_mapping_3d), pointer :: map
  sll_real64, allocatable         :: x(:,:,:), y(:,:,:), z(:,:,:)
  sll_real64, allocatable         :: efield(:), bfield(:)
  sll_real64, allocatable         :: efield_ref(:), bfield_ref(:)
  sll_real64, allocatable         :: efield_val(:), bfield_val(:)
  sll_real64, allocatable         :: rho(:), rho_ref(:), current(:)
  sll_int32                       :: i, j, k, istep, nsteps, ind
  sll_real64                      :: w1, w2
  sll_real64                      :: time
  sll_real64                      :: dt
  sll_int32                       :: deg(3)
  sll_real64                      :: params(6),jmatrix(3,3)
  sll_real64                      :: error(6)
  sll_real64, dimension(3,2)      :: domain
  sll_real64                      :: xvec(3)
  type(sll_t_time_mark) :: start, end

  params=0._f64
  params(1)= sll_p_twopi
  params(2)= sll_p_twopi
  params(3)= sll_p_twopi
  
  allocate(map)
  call map%init(params,&
       sll_f_scaling_x1,&
       sll_f_scaling_x2,&
       sll_f_scaling_x3,&
       sll_f_scaling_jac11,&
       sll_f_scaling_jac12,&
       sll_f_scaling_jac13,&
       sll_f_scaling_jac21,&
       sll_f_scaling_jac22,&
       sll_f_scaling_jac23,&
       sll_f_scaling_jac31,&
       sll_f_scaling_jac32,&
       sll_f_scaling_jac33,&
       sll_f_scaling_jacobian, Lx=params(1:3))

  print*,'Initialisation mapping finished'

  ! Define computational domain in physical coordinates
  eta1_min = 0._f64
  eta1_max = 1._f64
  nc_eta = [16, 16, 2 ]
  nc_total = product(nc_eta)
  delta_eta = 1._f64/real(nc_eta,f64)
  domain(1,:) = [eta1_min, params(1)]
  domain(2,:) = [eta1_min, params(2)]
  domain(3,:) = [eta1_min, params(3)]
  ! Set spline degree of 0-forms
  deg = [3, 3, 1]
  ! Time loop
  dt = 0.01_f64
  nsteps = 2
  w1 = sqrt( 3._f64 )
  w2 = + sqrt(3.0_f64)
  time = 0.0_f64
  
  call sll_s_set_time_mark( start )
  ! Initialise maxwell FEM object
  call maxwell_3d%init( domain, nc_eta, deg, map)
  print*,'Initialisation finished'
  call sll_s_set_time_mark( end )
  write(*, "(A, F10.3)") "Maxwell init run time [s] = ", sll_f_time_elapsed_between( start, end)
   
  call sll_s_set_time_mark( start )
  allocate(x(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(y(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(z(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(efield(1:nc_total*3))
  allocate(bfield(1:nc_total*3))
  allocate(efield_ref(1:nc_total*3))
  allocate(bfield_ref(1:nc_total*3))
  allocate(efield_val(1:nc_total*3))
  allocate(bfield_val(1:nc_total*3))
  allocate( rho( nc_total) )
  allocate( rho_ref( nc_total) )
  allocate( current(1:nc_total*3))

  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           x(i,j,k) = eta1_min + real(i-1,f64)*delta_eta(1)
           y(i,j,k) = eta1_min + real(j-1,f64)*delta_eta(2)
           z(i,j,k) = eta1_min + real(k-1,f64)*delta_eta(3)
        end do
     end do
  end do

  !test of L2projection 1form
  call maxwell_3d%L2projection( 1, 0, efield, sin_k, cos_k, constant )
  call evaluate_spline_3d( nc_eta, [deg(1)-1,deg(2),deg(3)], efield(1:nc_total), efield_val(1:nc_total) )
  call evaluate_spline_3d( nc_eta, [deg(1),deg(2)-1,deg(3)], efield(nc_total+1:2*nc_total), efield_val(nc_total+1:2*nc_total) )
  call evaluate_spline_3d( nc_eta, [deg(1),deg(2),deg(3)-1], efield(2*nc_total+1:3*nc_total), efield_val(2*nc_total+1:3*nc_total) )

  !reference solution
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           jmatrix=map%jacobian_matrix_inverse_transposed([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield(ind)           = efield_val(ind)*jmatrix(1,1) + efield_val(nc_total+ind)*jmatrix(1,2) + efield_val(2*nc_total+ind)*jmatrix(1,3)
           efield(nc_total+ind)  = efield_val(ind)*jmatrix(2,1) + efield_val(nc_total+ind)*jmatrix(2,2) + efield_val(2*nc_total+ind)*jmatrix(2,3)
           efield(2*nc_total+ind)= efield_val(ind)*jmatrix(3,1) + efield_val(nc_total+ind)*jmatrix(3,2) + efield_val(2*nc_total+ind)*jmatrix(3,3)

           xvec = map%get_x([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind)= sin_k(xvec)
           efield_ref(nc_total+ind)= cos_k(xvec)
           efield_ref(2*nc_total+ind)= constant([x(i,j,k), y(i,j,k), z(i,j,k)])
           ind = ind+1
        end do
     end do
  end do

  error(1) = maxval(abs(efield(1:nc_total)-efield_ref(1:nc_total)))
  print*, 'Error L2 projection 1form:', maxval(abs(efield(1:nc_total)-efield_ref(1:nc_total))), maxval(abs(efield(nc_total+1:2*nc_total)-efield_ref(nc_total+1:2*nc_total))), maxval(abs(efield(2*nc_total+1:3*nc_total)-efield_ref(2*nc_total+1:3*nc_total)))

!!$  call sll_s_plot_two_fields_1d('Ex3dcos',nc_total, efield(1:nc_total),efield_ref(1:nc_total),nc_eta(1),0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Ey3dsin',nc_total, efield(nc_total+1:2*nc_total),efield_ref(nc_total+1:2*nc_total),nc_eta(1),0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Ez3dconstant',nc_total, efield(2*nc_total+1:3*nc_total),efield_ref(2*nc_total+1:3*nc_total),nc_eta(1),0.0_f64)


  !test of L2projection 2form
  call maxwell_3d%L2projection( 2, 0, bfield, sin_k, cos_k, constant )
  call evaluate_spline_3d( nc_eta, [deg(1),deg(2)-1,deg(3)-1], bfield(1:nc_total), bfield_val(1:nc_total) )
  call evaluate_spline_3d( nc_eta, [deg(1)-1,deg(2),deg(3)-1], bfield(nc_total+1:2*nc_total), bfield_val(nc_total+1:2*nc_total) )
  call evaluate_spline_3d( nc_eta, [deg(1)-1,deg(2)-1,deg(3)], bfield(2*nc_total+1:3*nc_total), bfield_val(2*nc_total+1:3*nc_total) )

  !reference solution
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           jmatrix=map%jacobian_matrix([x(i,j,k), y(i,j,k), z(i,j,k)])/map%jacobian([x(i,j,k), y(i,j,k), z(i,j,k)])
           bfield(ind)           = bfield_val(ind)*jmatrix(1,1) + bfield_val(nc_total+ind)*jmatrix(1,2)+ bfield_val(2*nc_total+ind)*jmatrix(1,3)
           bfield(nc_total+ind)  = bfield_val(ind)*jmatrix(2,1) + bfield_val(nc_total+ind)*jmatrix(2,2)+ bfield_val(2*nc_total+ind)*jmatrix(2,3)
           bfield(2*nc_total+ind)= bfield_val(ind)*jmatrix(3,1) + bfield_val(nc_total+ind)*jmatrix(3,2)+ bfield_val(2*nc_total+ind)*jmatrix(3,3)

           xvec = map%get_x([x(i,j,k), y(i,j,k), z(i,j,k)])
           bfield_ref(ind)= sin_k(xvec)
           bfield_ref(nc_total+ind)= cos_k(xvec)
           bfield_ref(2*nc_total+ind)= constant([x(i,j,k), y(i,j,k), z(i,j,k)])
           ind = ind+1
        end do
     end do
  end do

  error(1) = maxval(abs(bfield(1:nc_total)-bfield_ref(1:nc_total)))
  print*, 'Error L2 projection 2form:', maxval(abs(bfield(1:nc_total)-bfield_ref(1:nc_total))), maxval(abs(bfield(nc_total+1:2*nc_total)-bfield_ref(nc_total+1:2*nc_total))), maxval(abs(bfield(2*nc_total+1:3*nc_total)-bfield_ref(2*nc_total+1:3*nc_total)))
!!$
!!$  call sll_s_plot_two_fields_1d('Bx3dcos',nc_total, bfield(1:nc_total),bfield_ref(1:nc_total),nc_eta(1),0.0_f64)
!!$  call sll_s_plot_two_fields_1d('By3dsin',nc_total, bfield(nc_total+1:2*nc_total),bfield_ref(nc_total+1:2*nc_total),nc_eta(1),0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Bz3dconstant',nc_total, bfield(2*nc_total+1:3*nc_total),bfield_ref(2*nc_total+1:3*nc_total),nc_eta(1),0.0_f64)


  ! Poisson problem
  call maxwell_3d%compute_rhs_from_function( 0, 1, rho, cos_k )
  call maxwell_3d%compute_E_from_rho( rho, efield )
  call evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)], efield(1:nc_total),  &
       efield_val(1:nc_total))
  call evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)], efield(1+nc_total:nc_total*2),  &
       efield_val(1+nc_total:nc_total*2))
  call evaluate_spline_3d ( nc_eta, [deg(1),deg(2),deg(3)-1], efield(1+nc_total*2:nc_total*3),  &
       efield_val(1+nc_total*2:nc_total*3))
  ! Reference solution
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           jmatrix=map%jacobian_matrix_inverse_transposed([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield(ind)           = efield_val(ind)*jmatrix(1,1) + efield_val(nc_total+ind)*jmatrix(1,2) + efield_val(2*nc_total+ind)*jmatrix(1,3)
           efield(nc_total+ind)  = efield_val(ind)*jmatrix(2,1) + efield_val(nc_total+ind)*jmatrix(2,2) + efield_val(2*nc_total+ind)*jmatrix(2,3)
           efield(2*nc_total+ind)= efield_val(ind)*jmatrix(3,1) + efield_val(nc_total+ind)*jmatrix(3,2) + efield_val(2*nc_total+ind)*jmatrix(3,3)

           xvec = map%get_x([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind) = sin_3k(xvec)
           efield_ref(ind+nc_total) = efield_ref(ind)
           efield_ref(ind+nc_total*2) = 0._f64!efield_ref(ind)
           ind = ind+1
        end do
     end do
  end do
  error(1) = maxval(abs(efield-efield_ref))
  print*, 'Error Poisson:', error(1)

  !L2 norm
  call maxwell_3d%L2projection( 1, 0, efield, cos_k, cos_k, cos_k )
  error(2) = abs(maxwell_3d%inner_product( efield, efield, 1, 1 )- 0.5_f64*params(1)*params(2)*params(3)) 
  print*, 'Error in L2 norm squared:', error(2)


  ! Solve the two equations
  ! \partial_t E = \curl B
  ! \partial_t B = -\curl E
  ! on a staggered time grid
  ! We use the solution
  ! E(x,t) =  \begin{pmatrix} \cos(x_1+x_2+x_3 - \sqrt{3} t) \\ -2\cos(x_1+x_2+x_3) \\ \cos(x_1+x_2+x_3 - \sqrt{3} t) \end{pmatrix}
  ! B(x,t) = \begin{pmatrix} \sqrt{3} \cos(x_1+x_2+x_3 - \sqrt{3} t) \\ 0 \\ -\sqrt{3} \cos(x_1+x_2+x_3 - \sqrt{3} t) \end{pmatrix}

  ! Now assemble initial efield and bfield
  ! B starts at time -0.5\Delta t
  time=-0.5_f64*dt
  call maxwell_3d%L2projection( 2, 0, bfield, cos_k, zero, cos_k )
  bfield(1:nc_total) =  w2 *  bfield(1:nc_total)
  bfield(nc_total*2+1:nc_total*3) = - w2 *  bfield(nc_total*2+1:nc_total*3)


  ! E starts at time 0
  time = 0.0_f64

  call maxwell_3d%L2projection( 1, 0, efield, cos_k, cos_k, cos_k )
  efield(nc_total+1:nc_total*2) = - 2.0_f64 *  efield(nc_total+1:nc_total*2)




  !curl part of maxwell
  ! Time stepping
  do istep = 1, nsteps
!!$     call maxwell_3d%compute_b_from_e( dt, efield, bfield )
!!$     call maxwell_3d%compute_e_from_b( dt, bfield, efield )
     call maxwell_3d%compute_curl_part( dt, efield, bfield )
  end do


  ! Evaluate E and B at the grid points
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)-1], bfield(1:nc_total), bfield_val(1:nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)-1], bfield(1+nc_total:2*nc_total), bfield_val(1+nc_total:2*nc_total)  )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2)-1,deg(3)], bfield(1+nc_total*2:3*nc_total), bfield_val(1+nc_total*2:3*nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)], efield(1:nc_total), efield_val(1:nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)], efield(1+nc_total:2*nc_total), efield_val(1+nc_total:2*nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2),deg(3)-1], efield(1+nc_total*2:3*nc_total), efield_val(1+nc_total*2:3*nc_total) )

  ! Reference solutions
  time = real(nsteps,f64)*dt
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           jmatrix=map%jacobian_matrix_inverse_transposed([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield(ind)           = efield_val(ind)*jmatrix(1,1) + efield_val(nc_total+ind)*jmatrix(1,2) + efield_val(2*nc_total+ind)*jmatrix(1,3)
           efield(nc_total+ind)  = efield_val(ind)*jmatrix(2,1) + efield_val(nc_total+ind)*jmatrix(2,2) + efield_val(2*nc_total+ind)*jmatrix(2,3)
           efield(2*nc_total+ind)= efield_val(ind)*jmatrix(3,1) + efield_val(nc_total+ind)*jmatrix(3,2) + efield_val(2*nc_total+ind)*jmatrix(3,3)

           xvec = map%get_x([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind) = cos_k(xvec)
           efield_ref(ind+nc_total) = -2.0_f64*efield_ref(ind)
           efield_ref(ind+nc_total*2) = efield_ref(ind)
           ind = ind+1
        end do
     end do
  end do

  time = (real(nsteps,f64)-0.5_f64)*dt
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           jmatrix=map%jacobian_matrix([x(i,j,k), y(i,j,k), z(i,j,k)])/map%jacobian([x(i,j,k), y(i,j,k), z(i,j,k)])
           bfield(ind)           = bfield_val(ind)*jmatrix(1,1) + bfield_val(nc_total+ind)*jmatrix(1,2)+ bfield_val(2*nc_total+ind)*jmatrix(1,3)
           bfield(nc_total+ind)  = bfield_val(ind)*jmatrix(2,1) + bfield_val(nc_total+ind)*jmatrix(2,2)+ bfield_val(2*nc_total+ind)*jmatrix(2,3)
           bfield(2*nc_total+ind)= bfield_val(ind)*jmatrix(3,1) + bfield_val(nc_total+ind)*jmatrix(3,2)+ bfield_val(2*nc_total+ind)*jmatrix(3,3)

           xvec = map%get_x([x(i,j,k), y(i,j,k), z(i,j,k)])
           bfield_ref(ind) = w2*cos_k(xvec)
           bfield_ref(ind+nc_total) = 0.0_f64
           bfield_ref(ind+nc_total*2) = -bfield_ref(ind)
           ind = ind+1
        end do
     end do
  end do
  error(3) = maxval(abs(efield-efield_ref))
  error(4) = maxval(abs(bfield-bfield_ref))
  print*, 'Error efield:', error(3)
  print*, 'Error bfield:', error(4)


  ! Test compute_e_from_j
  call maxwell_3d%L2projection( 1, 0, efield, cos_k, zero, zero )
  call maxwell_3d%compute_rhs_from_function( 1, 1, current(1:nc_total), sin_k, zero, zero )
  call maxwell_3d%compute_rhs_from_function( 1, 2, current(1+nc_total:2*nc_total), sin_k, zero, zero )
  call maxwell_3d%compute_rhs_from_function( 1, 3, current(1+2*nc_total:3*nc_total), sin_k, zero, zero )
  current=dt*current
  call maxwell_3d%compute_E_from_j( current, efield )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)], efield(1:nc_total),  &
       efield_val(1:nc_total))
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)], efield(1+nc_total:2*nc_total), efield_val(1+nc_total:2*nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2),deg(3)-1], efield(1+nc_total*2:3*nc_total), efield_val(1+nc_total*2:3*nc_total) )

  ! Reference solution
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           jmatrix=map%jacobian_matrix_inverse_transposed([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield(ind)           = efield_val(ind)*jmatrix(1,1) + efield_val(nc_total+ind)*jmatrix(1,2) + efield_val(2*nc_total+ind)*jmatrix(1,3)
           efield(nc_total+ind)  = efield_val(ind)*jmatrix(2,1) + efield_val(nc_total+ind)*jmatrix(2,2) + efield_val(2*nc_total+ind)*jmatrix(2,3)
           efield(2*nc_total+ind)= efield_val(ind)*jmatrix(3,1) + efield_val(nc_total+ind)*jmatrix(3,2) + efield_val(2*nc_total+ind)*jmatrix(3,3)

           xvec = map%get_x([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind) =(cos_k(xvec)-&
                dt*sin_k(xvec))
           ind = ind+1
        end do
     end do
  end do
  error(5) = maxval(abs(efield(1:nc_total)-efield_ref(1:nc_total)))
  print*, 'Error compute_e_from_j:', error(5)

  ! Test compute_rho_from_e
  call maxwell_3d%compute_rhs_from_function( 0, 1, rho_ref, cos_k )
  call maxwell_3d%L2projection( 1, 0,  efield, sin_3k, sin_3k, sin_3k )

  call maxwell_3d%compute_rho_from_e( efield, rho )

  error(6) =  maxval( abs( rho - rho_ref ) )
  print*, 'Error compute_rho_from_e:', error(6)


  ! Clean up
  call maxwell_3d%free()
  call map%free()
  map => null()
  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(efield)
  deallocate(bfield)
  deallocate(efield_ref)
  deallocate(bfield_ref)
  deallocate(efield_val)
  deallocate(bfield_val)
  deallocate( rho )
  deallocate( rho_ref )
  deallocate( current )
 
  if ( error(1) < 5.d-5 .AND. error(2) < 3.d-5 .AND. error(3) < 5.0d-2 .AND. error(4) < 5.d-2 .AND. error(5)<8.d-5 .AND. error(6)<1.d-8) then
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

    cos_k = cos((x(1)+x(2))-w1*time) 
  end function cos_k

  function sin_k(x)
    sll_real64             :: sin_k
    sll_real64, intent(in) :: x(3)

    sin_k = sin((x(1)+x(2))-w1*time) 
  end function sin_k

   function sin_3k(x)
    sll_real64             :: sin_3k
    sll_real64, intent(in) :: x(3)

    sin_3k = sin(x(1)+x(2)-w1*time)/2._f64 
  end function sin_3k
 
  function cos_3k(x)
    sll_real64             :: cos_3k
    sll_real64, intent(in) :: x(3)

    cos_3k = cos(x(1)+x(2)-w1*time)/2._f64
  end function cos_3k

  function constant(x)
    sll_real64             :: constant
    sll_real64, intent(in) :: x(3)

    constant=1._f64
  end function constant


  function zero(x)
    sll_real64             :: zero
    sll_real64, intent(in) :: x(3)

    zero=0._f64
  end function zero



  subroutine evaluate_spline_3d ( ndofs, deg, dofs, vals )
    sll_int32, intent( in ) :: ndofs(3)
    sll_int32, intent( in ) :: deg(3)
    sll_real64, intent( in ) :: dofs(:)
    sll_real64, intent( out ) :: vals(:)

    sll_int32 :: i,j,k,istart,iend
    sll_real64 :: a_in(ndofs(2)), a_out(ndofs(2)),b_in(ndofs(3)), b_out(ndofs(3))


    istart = 1
    iend = ndofs(1)
    do k=1,ndofs(3)
       do j=1,ndofs(2)
          call sll_s_eval_uniform_periodic_spline_curve_with_zero(deg(1), dofs(istart:iend), vals(istart:iend))
          istart = iend+1
          iend = iend + ndofs(1)
       end do
    end do

    do k=1,ndofs(3)
       do i=1,ndofs(1)
          istart = (k-1)*ndofs(2)*ndofs(1)+i
          do j=1,ndofs(2)
             a_in(j) = vals(istart+(j-1)*ndofs(1))
          end do
          call sll_s_eval_uniform_periodic_spline_curve_with_zero(deg(2), a_in, a_out)
          do j=1,ndofs(2)
             vals(istart+(j-1)*ndofs(1)) = a_out(j)
          end do
       end do
    end do

    do j=1,ndofs(2)
       do i=1,ndofs(1)
          istart = (j-1)*ndofs(1)+i
          do k=1,ndofs(3)
             b_in(k) = vals(istart+(k-1)*ndofs(1)*ndofs(2))
          end do
          call sll_s_eval_uniform_periodic_spline_curve_with_zero(deg(3), b_in, b_out)
          do k=1,ndofs(3)
             vals(istart+(k-1)*ndofs(1)*ndofs(2)) = b_out(k)
          end do
       end do
    end do


  end subroutine evaluate_spline_3d


!!$  function sin_k(x)
!!$    sll_real64             :: sin_k
!!$    sll_real64, intent(in) :: x
!!$
!!$    sin_k = sin(mode*2*sll_p_pi*x/Lx) 
!!$  end function sin_k
end program test_maxwell_3d_trafo
