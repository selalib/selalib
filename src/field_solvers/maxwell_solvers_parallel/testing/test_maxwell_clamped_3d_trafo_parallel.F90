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
program test_maxwell_clamped_3d_trafo_parallel
  !------------------------------------------------------------------------
  !  test 3D Maxwell_Clamped spline finite element solver on a periodic grid
  !------------------------------------------------------------------------
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_low_level_bsplines, only: &
       sll_s_eval_uniform_periodic_spline_curve

  use sll_m_splines_pp

  use sll_m_constants, only: &
       sll_p_pi, sll_p_twopi, sll_p_fourpi

  use sll_m_maxwell_clamped_3d_trafo_parallel

  use sll_m_maxwell_clamped_3d_trafo

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_3d_coordinate_transformations

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark

  use sll_m_maxwell_1d_base, only: &
       sll_s_plot_two_fields_1d

  use sll_m_collective, only: &
       sll_s_boot_collective, &
       sll_s_halt_collective

  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !  type(sll_t_arbitrary_degree_spline_1d), pointer :: aspl
  !  sll_real64, dimension(:), allocatable :: knots

  sll_real64                      :: eta1_max, eta1_min
  sll_real64                      :: delta_eta(3)
  sll_int32                       :: nc_eta(3), nc_total, nc_total0, nc_total1
  type(sll_t_maxwell_clamped_3d_trafo_parallel)    :: maxwell_solver
  type(sll_t_mapping_3d), pointer :: map
  sll_real64, allocatable         :: x(:,:,:), y(:,:,:), z(:,:,:)
  sll_real64, allocatable         :: efield(:), bfield(:)
  sll_real64, allocatable         :: efield_ref(:), bfield_ref(:)
  sll_real64, allocatable         :: efield_val(:), bfield_val(:)
  sll_real64, allocatable         :: efield_val1(:), bfield_val1(:)
  sll_real64, allocatable         :: rho(:), rho_ref(:), current(:)
  sll_real64, allocatable         :: rho_val(:), rho_val1(:), phi(:)
  sll_int32                       :: i, j, k, istep, nsteps, ind
  sll_real64                      :: xvec(3)
  sll_real64                      :: time
  sll_real64                      :: dt
  sll_int32                       :: deg(3), boundary(3)
  sll_real64                      :: params(12),jmatrix(3,3)
  sll_real64                      :: error(6) = 0._f64
  sll_real64, dimension(3,2)      :: domain
  sll_real64 :: Lx(3)
  type(sll_t_time_mark) :: start, end

  type(sll_t_spline_pp_3d) :: spline_pp_0
  type(sll_t_spline_pp_3d) :: spline_pp_11
  type(sll_t_spline_pp_3d) :: spline_pp_12
  type(sll_t_spline_pp_3d) :: spline_pp_13
  type(sll_t_spline_pp_3d) :: spline_pp_21
  type(sll_t_spline_pp_3d) :: spline_pp_22
  type(sll_t_spline_pp_3d) :: spline_pp_23

  call sll_s_boot_collective()

  params=0._f64
  params(1)= 0.1_f64
  params(2)= sll_p_twopi+params(1)
  params(3)= 1.0_f64!sll_p_twopi

  Lx(1:2) = 2._f64*(params(2)-params(1))
  Lx(3) = params(3)

  allocate(map)
  call map%init(params,&
       sll_f_cylindrical_x1, sll_f_cylindrical_x2, sll_f_cylindrical_x3,&
       sll_f_cylindrical_jac11, sll_f_cylindrical_jac12, sll_f_cylindrical_jac13,&
       sll_f_cylindrical_jac21, sll_f_cylindrical_jac22, sll_f_cylindrical_jac23,&
       sll_f_cylindrical_jac31, sll_f_cylindrical_jac32, sll_f_cylindrical_jac33,&
       sll_f_cylindrical_jacobian, Lx=Lx)


  ! Define computational domain in logical coordinates
  eta1_min = 0._f64
  eta1_max = 1._f64
  nc_eta = [32,32,2]
  nc_total = product(nc_eta+1)
  delta_eta = 1._f64/real(nc_eta,f64)
  domain(1,:) = [eta1_min, params(1)]
  domain(2,:) = [eta1_min, params(2)]
  domain(3,:) = [eta1_min, params(3)]
  ! Set spline degree of 0-forms
  deg = [3,3,1]
  nc_total0 = (nc_eta(1)+deg(1))*nc_eta(2)*nc_eta(3)
  nc_total1 = (nc_eta(1)+deg(1)-1)*nc_eta(2)*nc_eta(3)
  ! Time loop
  dt = 0.1_f64
  nsteps = 1
  time = 0.0_f64
  
  if( deg(1) == 2 ) then
     boundary = [ sll_p_boundary_clamped_square, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else if( deg(1) == 3 ) then
     boundary = [ sll_p_boundary_clamped_cubic, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else
     boundary = [ sll_p_boundary_clamped, sll_p_boundary_periodic, sll_p_boundary_periodic]
  end if

  call sll_s_set_time_mark( start )
  ! Initialise maxwell_clamped FEM object
  call maxwell_solver%init( domain, nc_eta, deg, boundary, map)

  call sll_s_spline_pp_init_3d(spline_pp_0, deg, nc_eta)
  call sll_s_spline_pp_init_3d(spline_pp_11, [deg(1)-1,deg(2),deg(3)], nc_eta, boundary)
  call sll_s_spline_pp_init_3d(spline_pp_12, [deg(1),deg(2)-1,deg(3)], nc_eta, boundary)
  call sll_s_spline_pp_init_3d(spline_pp_13, [deg(1),deg(2),deg(3)-1], nc_eta, boundary)
  call sll_s_spline_pp_init_3d(spline_pp_21, [deg(1),deg(2)-1,deg(3)-1], nc_eta, boundary)
  call sll_s_spline_pp_init_3d(spline_pp_22, [deg(1)-1,deg(2),deg(3)-1], nc_eta, boundary)
  call sll_s_spline_pp_init_3d(spline_pp_23, [deg(1)-1,deg(2)-1,deg(3)], nc_eta, boundary)

  print*,'Initialisation finished'
  call sll_s_set_time_mark( end )
  write(*, "(A, F10.3)") "Maxwell_Clamped init run time [s] = ", sll_f_time_elapsed_between( start, end)
  call sll_s_set_time_mark( start )
  allocate(x(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(y(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(z(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(efield(1:nc_total1+nc_total0*2))
  allocate(bfield(1:nc_total0+nc_total1*2))
  allocate(efield_ref(1:nc_total*3))
  allocate(bfield_ref(1:nc_total*3))
  allocate(efield_val(1:nc_total*3))
  allocate(bfield_val(1:nc_total*3))
  allocate(efield_val1(1:nc_total*3))
  allocate(bfield_val1(1:nc_total*3))

  allocate( rho( nc_total0) )
  allocate( rho_ref( nc_total0 ) )
  allocate( phi( nc_total0) )
  allocate( current(1:nc_total1+nc_total0*2) )
  allocate( rho_val( nc_total ) )
  allocate( rho_val1( nc_total ) )

  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           x(i,j,k) = eta1_min + real(i-1,f64)*delta_eta(1)
           y(i,j,k) = eta1_min + real(j-1,f64)*delta_eta(2)
           z(i,j,k) = eta1_min + real(k-1,f64)*delta_eta(3)
        end do
     end do
  end do

!!$
!!$  call maxwell_solver%L2projection( 1, 0, efield, E_1k, E_2k, zero)
!!$
!!$  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_11, nc_eta, efield(1:nc_total1),  &
!!$       efield_val(1:nc_total))
!!$  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_12, nc_eta, efield(1+nc_total1:nc_total1+nc_total0),  &
!!$       efield_val(1+nc_total:2*nc_total))
!!$  call sll_s_spline_evaluate_basis_b_form_3d_clamped ( spline_pp_13, nc_eta, efield(1+nc_total1+nc_total0:nc_total1+nc_total0*2),  &
!!$       efield_val(1+2*nc_total:3*nc_total))
!!$ 
!!$  ! Reference solution
!!$  ind = 1
!!$  do k = 1, nc_eta(3)+1
!!$     do j = 1, nc_eta(2)+1
!!$        do i = 1, nc_eta(1)+1
!!$           jmatrix=map%jacobian_matrix_inverse_transposed([x(i,j,k), y(i,j,k), z(i,j,k)])
!!$           efield_val1(ind)           = efield_val(ind)*jmatrix(1,1) + efield_val(nc_total+ind)*jmatrix(1,2) + efield_val(2*nc_total+ind)*jmatrix(1,3)
!!$           efield_val1(nc_total+ind)  = efield_val(ind)*jmatrix(2,1) + efield_val(nc_total+ind)*jmatrix(2,2) + efield_val(2*nc_total+ind)*jmatrix(2,3)
!!$           efield_val1(2*nc_total+ind)= efield_val(ind)*jmatrix(3,1) + efield_val(nc_total+ind)*jmatrix(3,2) + efield_val(2*nc_total+ind)*jmatrix(3,3)
!!$
!!$           xvec = map%get_x([x(i,j,k), y(i,j,k), z(i,j,k)])
!!$           
!!$           efield_ref(ind) = E_1k(xvec)
!!$           efield_ref(ind+nc_total) = E_2k(xvec)
!!$           efield_ref(ind+nc_total*2) = 0._f64
!!$           ind = ind+1
!!$        end do
!!$     end do
!!$  end do
!!$  print*, maxval(abs(efield_val1(1:nc_total)-efield_ref(1:nc_total))), maxval(abs(efield_val1(1+nc_total:2*nc_total)-efield_ref(1+nc_total:2*nc_total))), maxval(abs(efield_val1(1+2*nc_total:3*nc_total)-efield_ref(1+2*nc_total:3*nc_total)))


  
  ! Poisson problem
  call maxwell_solver%compute_rhs_from_function( 0, 1, rho, rho_k )
  call maxwell_solver%compute_E_from_rho( rho, efield )


  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_11, nc_eta, efield(1:nc_total1),  &
       efield_val(1:nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_12, nc_eta, efield(1+nc_total1:nc_total1+nc_total0),  &
       efield_val(1+nc_total:2*nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped ( spline_pp_13, nc_eta, efield(1+nc_total1+nc_total0:nc_total1+nc_total0*2),  &
       efield_val(1+2*nc_total:3*nc_total))

  ! Reference solution
  ind = 1
  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           jmatrix=map%jacobian_matrix_inverse_transposed([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_val1(ind)           = efield_val(ind)*jmatrix(1,1) + efield_val(nc_total+ind)*jmatrix(1,2) + efield_val(2*nc_total+ind)*jmatrix(1,3)
           efield_val1(nc_total+ind)  = efield_val(ind)*jmatrix(2,1) + efield_val(nc_total+ind)*jmatrix(2,2) + efield_val(2*nc_total+ind)*jmatrix(2,3)
           efield_val1(2*nc_total+ind)= efield_val(ind)*jmatrix(3,1) + efield_val(nc_total+ind)*jmatrix(3,2) + efield_val(2*nc_total+ind)*jmatrix(3,3)

           xvec = map%get_x([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind) = E_1k(xvec)
           efield_ref(ind+nc_total) = E_2k(xvec)
           efield_ref(ind+nc_total*2) = 0._f64
           ind = ind+1
        end do
     end do
  end do
  error(1) = maxval(abs(efield_val1-efield_ref))
  print*, 'Error Poisson:', maxval(abs(efield_val1(1:nc_total)-efield_ref(1:nc_total))), maxval(abs(efield_val1(nc_total+1:2*nc_total)-efield_ref(nc_total+1:2*nc_total))), maxval(abs(efield_val1(2*nc_total+1:3*nc_total)-efield_ref(2*nc_total+1:3*nc_total)))

!!$  call sll_s_plot_two_fields_1d('Poissonx',nc_total, efield_val1(1:nc_total),efield_ref(1:nc_total),nc_eta(1),0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Poissony',nc_total, efield_val1(nc_total+1:2*nc_total),efield_ref(nc_total+1:2*nc_total),nc_eta(1),0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Poissonz',nc_total, efield_val1(2*nc_total+1:3*nc_total),efield_ref(2*nc_total+1:3*nc_total),nc_eta(1),0.0_f64)


  ! Solve the two equations
  ! \partial_t E = \curl B
  ! \partial_t B = -\curl E
  ! Now assemble initial efield and bfield
  !L2 error
  efield = 0._f64
  call maxwell_solver%L2projection( 1, 0, efield, poly_1k, poly_2k, zero )
  error(2) = abs(maxwell_solver%inner_product( efield, efield, 1, 1 )-sll_p_twopi*0.25_f64*(params(2)**4-params(1)**4))
  print*, 'Error in L2 norm squared:', error(2)



  efield = 0._f64
  call maxwell_solver%L2projection( 1, 0, efield, sin_k, zero, zero )
  !curl part of maxwell_clamped
  ! Time stepping
  do istep = 1, nsteps
!!$     call maxwell_solver%compute_curl_part( dt, efield, bfield )
     call maxwell_solver%compute_b_from_e( dt, efield, bfield )
     call maxwell_solver%compute_e_from_b( dt, bfield, efield )
  end do

  ! Evaluate E and B at the grid points
   call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_11, nc_eta, efield(1:nc_total1),  &
       efield_val(1:nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_12, nc_eta, efield(1+nc_total1:nc_total1+nc_total0),  &
       efield_val(1+nc_total:2*nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped ( spline_pp_13, nc_eta , efield(1+nc_total1+nc_total0:nc_total1+nc_total0*2),  &
       efield_val(1+2*nc_total:3*nc_total))
  
  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_21, nc_eta, bfield(1:nc_total0),  &
       bfield_val(1:nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_22, nc_eta, bfield(1+nc_total0:nc_total0+nc_total1),  &
       bfield_val(1+nc_total:2*nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped ( spline_pp_23, nc_eta , bfield(1+nc_total0+nc_total1:nc_total0+nc_total1*2),  &
       bfield_val(1+2*nc_total:3*nc_total))

  ! Reference solutions
  !time = real(nsteps,f64)*dt
  ind = 1
  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           jmatrix=map%jacobian_matrix_inverse_transposed([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_val1(ind)           = efield_val(ind)*jmatrix(1,1) + efield_val(nc_total+ind)*jmatrix(1,2) + efield_val(2*nc_total+ind)*jmatrix(1,3)
           efield_val1(nc_total+ind)  = efield_val(ind)*jmatrix(2,1) + efield_val(nc_total+ind)*jmatrix(2,2) + efield_val(2*nc_total+ind)*jmatrix(2,3)
           efield_val1(2*nc_total+ind)= efield_val(ind)*jmatrix(3,1) + efield_val(nc_total+ind)*jmatrix(3,2) + efield_val(2*nc_total+ind)*jmatrix(3,3)

           xvec = map%get_x([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind) = ex_k(xvec)
           efield_ref(ind+nc_total) = ey_k(xvec)
           efield_ref(ind+nc_total*2) = 0._f64
           ind = ind+1
        end do
     end do
  end do

  ind = 1
  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           jmatrix=map%jacobian_matrix([x(i,j,k), y(i,j,k), z(i,j,k)])/map%jacobian([x(i,j,k), y(i,j,k), z(i,j,k)])
           bfield_val1(ind)           = bfield_val(ind)*jmatrix(1,1) + bfield_val(nc_total+ind)*jmatrix(1,2)+ bfield_val(2*nc_total+ind)*jmatrix(1,3)
           bfield_val1(nc_total+ind)  = bfield_val(ind)*jmatrix(2,1) + bfield_val(nc_total+ind)*jmatrix(2,2)+ bfield_val(2*nc_total+ind)*jmatrix(2,3)
           bfield_val1(2*nc_total+ind)= bfield_val(ind)*jmatrix(3,1) + bfield_val(nc_total+ind)*jmatrix(3,2)+ bfield_val(2*nc_total+ind)*jmatrix(3,3)

           xvec = map%get_x([x(i,j,k), y(i,j,k), z(i,j,k)])
           bfield_ref(ind) = 0._f64
           bfield_ref(ind+nc_total) = 0.0_f64
           bfield_ref(ind+nc_total*2) = bz_k(xvec)
           ind = ind+1
        end do
     end do
  end do
  error(3) = maxval(abs(efield_val1-efield_ref))
  error(4) = maxval(abs(bfield_val1-bfield_ref))
  print*, 'Error efield:', error(3)
  print*, 'Error bfield:', error(4)

!!$  call sll_s_plot_two_fields_1d('Ex3d',nc_total, efield_val1(1:nc_total),efield_ref(1:nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Ey3d',nc_total, efield_val1(nc_total+1:2*nc_total),efield_ref(nc_total+1:2*nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Ez3d',nc_total, efield_val1(2*nc_total+1:3*nc_total),efield_ref(2*nc_total+1:3*nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Bx3d',nc_total, bfield_val1(1:nc_total),bfield_ref(1:nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('By3d',nc_total, bfield_val1(nc_total+1:2*nc_total),bfield_ref(nc_total+1:2*nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Bz3d',nc_total, bfield_val1(2*nc_total+1:3*nc_total),bfield_ref(2*nc_total+1:3*nc_total),0,0.0_f64)



  ! Test compute_e_from_j
  current = 0._f64
  efield = 0._f64
  call maxwell_solver%L2projection( 1, 0, efield, sin_k, zero, zero )
  call maxwell_solver%compute_rhs_from_function( 1, 1, current(1:nc_total1), sin_k, zero, zero )
  call maxwell_solver%compute_rhs_from_function( 1, 2, current(1+nc_total1:nc_total1+nc_total0), sin_k, zero, zero )
  call maxwell_solver%compute_rhs_from_function( 1, 3, current(1+nc_total1+nc_total0:nc_total0*2+nc_total1), sin_k, zero, zero )

  call maxwell_solver%compute_E_from_j( dt*current, efield )
  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_11, nc_eta, efield(1:nc_total1),  &
       efield_val(1:nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_12, nc_eta, efield(1+nc_total1:nc_total1+nc_total0),  &
       efield_val(1+nc_total:2*nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped ( spline_pp_13, nc_eta , efield(1+nc_total1+nc_total0:nc_total1+nc_total0*2),  &
       efield_val(1+2*nc_total:3*nc_total))

  ! Reference solution
  ind = 1
  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           jmatrix=map%jacobian_matrix_inverse_transposed([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_val1(ind)           = efield_val(ind)*jmatrix(1,1) + efield_val(nc_total+ind)*jmatrix(1,2) + efield_val(2*nc_total+ind)*jmatrix(1,3)
           efield_val1(nc_total+ind)  = efield_val(ind)*jmatrix(2,1) + efield_val(nc_total+ind)*jmatrix(2,2) + efield_val(2*nc_total+ind)*jmatrix(2,3)
           efield_val1(2*nc_total+ind)= efield_val(ind)*jmatrix(3,1) + efield_val(nc_total+ind)*jmatrix(3,2) + efield_val(2*nc_total+ind)*jmatrix(3,3)

           xvec = map%get_x([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind) = sin_k(xvec) - dt*sin_k(xvec) 
           ind = ind+1
        end do
     end do
  end do
  error(5) = maxval(abs(efield_val1(1:nc_total)-efield_ref(1:nc_total)))
  print*, 'Error compute_e_from_j:', error(5)

  !call sll_s_plot_two_fields_1d('Current',nc_total, efield_val1(1:nc_total),efield_ref(1:nc_total),0,0.0_f64)

  ! Test compute_rho_from_e
  call maxwell_solver%compute_rhs_from_function( 0, 1, rho_ref, rho_k )
  call maxwell_solver%L2projection( 1, 0,  efield, E_1k, E_2k, zero )

  call maxwell_solver%compute_rho_from_e( efield, rho )

  error(6) =  maxval( abs( rho - rho_ref ) )
  print*, 'Error compute_rho_from_e:', error(6)

  !call sll_s_plot_two_fields_1d('Rho',nc_total0, rho,rho_ref,nc_eta(1),0.0_f64)


  ! Clean up
  call maxwell_solver%free()
  deallocate(x )
  deallocate(y )
  deallocate(z )
  deallocate(efield )
  deallocate(bfield )
  deallocate(efield_ref )
  deallocate(bfield_ref )
  deallocate(efield_val )
  deallocate(bfield_val )
  deallocate(efield_val1 )
  deallocate(bfield_val1 )

  deallocate( rho )
  deallocate( rho_ref )
  deallocate( current )
  deallocate( rho_val )
  deallocate( rho_val1 )
  deallocate( map )

  if ( error(1) < 2.1d-10 .AND. error(2) < 6.0d-12 .AND. error(3) < 9.d-3 .AND. error(4) < 3.d-4 .AND. error(5)<7.d-5 .AND. error(6)<4.d-12) then
     print*, 'PASSED.'
  else
     print*, 'FAILED.'
  end if
  call sll_s_set_time_mark( end )
  write(*, "(A, F10.3)") "Main part run time [s] = ", sll_f_time_elapsed_between( start, end)

  call sll_s_halt_collective()
contains
  
  function sin_k(x)
    sll_real64             :: sin_k
    sll_real64, intent(in) :: x(3)

    sin_k = sin(sqrt(x(1)**2+x(2)**2)-params(1)) 
  end function sin_k
  

  function cos_k(x)
    sll_real64             :: cos_k
    sll_real64, intent(in) :: x(3)

    cos_k = cos(sqrt(x(1)**2+x(2)**2)-params(1))
  end function cos_k

  function rho_k(x)
    sll_real64             :: rho_k
    sll_real64, intent(in) :: x(3)

    rho_k = -4._f64 + (params(1)+params(2))/sqrt(x(1)**2+x(2)**2)!sin(sqrt(x(1)**2+x(2)**2)-params(1)) -cos(sqrt(x(1)**2+x(2)**2)-params(1))/sqrt(x(1)**2+x(2)**2)
    !
    
  end function rho_k

  function phi_k(x)
    sll_real64             :: phi_k
    sll_real64, intent(in) :: x(3)
    
    phi_k = (sqrt(x(1)**2+x(2)**2)-params(1))*(sqrt(x(1)**2+x(2)**2)-params(2))!sin(sqrt(x(1)**2+x(2)**2)-params(1))
    !
    !x(1)**2+x(2)**2-sqrt(x(1)**2+x(2)**2)*(params(1)+params(2))+params(1)*params(2)
  end function phi_k

  function E_1k(x)
    sll_real64             :: E_1k
    sll_real64, intent(in) :: x(3)
   
    E_1k = -2._f64*x(1)+(x(1)/sqrt(x(1)**2+x(2)**2))*(params(1)+params(2))!-cos(sqrt(x(1)**2+x(2)**2)-params(1))*x(1)/sqrt(x(1)**2+x(2)**2)
    !
    
  end function E_1k

  function E_2k(x)
    sll_real64             :: E_2k
    sll_real64, intent(in) :: x(3)
  
    E_2k = -2._f64*x(2)+(x(2)/sqrt(x(1)**2+x(2)**2))*(params(1)+params(2))!-cos(sqrt(x(1)**2+x(2)**2)-params(1))*x(2)/sqrt(x(1)**2+x(2)**2)
    !
    
  end function E_2k

  function bz_k(x)
    sll_real64             :: bz_k
    sll_real64, intent(in) :: x(3)
  
    bz_k =  dt*cos_k(x)*x(2)/sqrt(x(1)**2+x(2)**2)
    
  end function bz_k

  function ex_k(x)
    sll_real64             :: ex_k
    sll_real64, intent(in) :: x(3)
  
    ex_k = sin_k(x)*(1-dt**2*x(2)**2/(x(1)**2+x(2)**2)) + dt**2 * cos_k(x)/sqrt(x(1)**2+x(2)**2) *(1-x(2)**2/(x(1)**2+x(2)**2)) 
    
  end function ex_k

  function ey_k(x)
    sll_real64             :: ey_k
    sll_real64, intent(in) :: x(3)
  
    ey_k = dt**2*x(1)*x(2)/(x(1)**2+x(2)**2) * ( sin_k(x) + cos_k(x)/sqrt(x(1)**2+x(2)**2) )
    
  end function ey_k


  function poly_1k(x)
    sll_real64             :: poly_1k
    sll_real64, intent(in) :: x(3)

    poly_1k = x(1)
  end function poly_1k

   function poly_2k(x)
    sll_real64             :: poly_2k
    sll_real64, intent(in) :: x(3)

    poly_2k = x(2)
  end function poly_2k

  function constant(x)
    sll_real64             :: constant
    sll_real64, intent(in) :: x(3)

    constant=2._f64
  end function constant


  function zero(x)
    sll_real64             :: zero
    sll_real64, intent(in) :: x(3)

    zero=0._f64
  end function zero


!!$  function phi_ex(x)
!!$    sll_real64             :: phi_ex
!!$    sll_real64, intent(in) :: x(3)
!!$
!!$    phi_ex = (1._f64-x(1)**2-x(2)**2)*cos(sll_p_twopi*x(1))*sin(sll_p_twopi*x(2))
!!$  end function phi_ex
!!$
!!$  function Ex_ex(x)
!!$    sll_real64             :: Ex_ex
!!$    sll_real64, intent(in) :: x(3)
!!$
!!$    Ex_ex = 2._f64*sin(sll_p_twopi*x(2))*(-sll_p_pi*(x(1)**2+x(2)**2-1._f64)*sin(sll_p_twopi*x(1))+x(1)*cos(sll_p_twopi*x(1)))
!!$  end function Ex_ex
!!$
!!$  function Ey_ex(x)
!!$    sll_real64             :: Ey_ex
!!$    sll_real64, intent(in) :: x(3)
!!$
!!$    Ey_ex = 2._f64*cos(sll_p_twopi*x(1))*(sll_p_pi*(x(1)**2+x(2)**2-1._f64)*cos(sll_p_twopi*x(2))+x(2)*sin(sll_p_twopi*x(2)))
!!$  end function Ey_ex
!!$
!!$
!!$  function rho_ex(x)
!!$    sll_real64             :: rho_ex
!!$    sll_real64, intent(in) :: x(3)
!!$
!!$    rho_ex = cos(sll_p_twopi*x(1))*(4._f64*(1._f64-2._f64*sll_p_pi**2*(x(1)**2+x(2)**2-1._f64))*sin(sll_p_twopi*x(2))+4._f64*sll_p_twopi*x(2)*cos(sll_p_twopi*x(2))) - 4._f64*sll_p_twopi*x(1)*sin(sll_p_twopi*x(1))*sin(sll_p_twopi*x(2))
!!$  end function rho_ex


end program test_maxwell_clamped_3d_trafo_parallel
