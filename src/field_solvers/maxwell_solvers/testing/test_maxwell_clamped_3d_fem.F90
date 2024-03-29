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
program test_maxwell_clamped_3d_fem
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

  use sll_m_maxwell_1d_base, only: &
       sll_s_plot_two_fields_1d

  use sll_m_constants, only: &
       sll_p_pi, sll_p_twopi

  use sll_m_maxwell_clamped_3d_fem

  use sll_m_constants, only: sll_p_pi, sll_p_twopi

   use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark

  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_real64 :: eta1_max, eta1_min
  sll_real64 :: delta_eta(3)

  sll_int32  :: nc_eta(3), nc_total, nc_total0, nc_total1

  type(sll_t_maxwell_clamped_3d_fem)  :: maxwell_solver
  sll_real64, allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
  sll_real64, allocatable :: efield(:), bfield(:), efield1(:), bfield1(:) 
  sll_real64, allocatable :: efield_ref(:), bfield_ref(:)
  sll_real64, allocatable :: efield_val(:), bfield_val(:)
  sll_real64, allocatable :: rho(:), rho_ref(:), current(:)
  sll_real64, allocatable :: rho_val(:), rho_val1(:)

  sll_int32                               :: i, j, k, istep, nsteps, ind
  sll_real64                                :: time
  sll_real64                                :: delta_t

  sll_real64                              :: Lx(3)
  sll_real64, dimension(3,2)              :: domain
  sll_int32                               :: deg(3), boundary(3)

  sll_real64                              :: error(6), error1(3), energy(2)
  sll_real64                              :: alpha, kx
  type(sll_t_time_mark) :: start, end
  type(sll_t_spline_pp_3d) :: spline_pp_0
  type(sll_t_spline_pp_3d) :: spline_pp_11
  type(sll_t_spline_pp_3d) :: spline_pp_12
  type(sll_t_spline_pp_3d) :: spline_pp_13
  type(sll_t_spline_pp_3d) :: spline_pp_21
  type(sll_t_spline_pp_3d) :: spline_pp_22
  type(sll_t_spline_pp_3d) :: spline_pp_23

  call sll_s_set_time_mark( start )

  alpha = 0.01_f64
  kx = 1._f64
  
  ! Define computational domain
  eta1_min = .0_f64; eta1_max = sll_p_twopi/kx
  nc_eta = [16,8,8]
  nc_total = product(nc_eta+1)
  Lx(1) = eta1_max-eta1_min
  Lx(2) = Lx(1); Lx(3) = Lx(2)
  delta_eta = Lx/real(nc_eta,f64)
  domain(1,:) = [eta1_min, eta1_max]
  domain(2,:) = [eta1_min, eta1_max]
  domain(3,:) = [eta1_min, eta1_max]
  ! Set spline degree of 0-forms
  deg = [3,3,3]
  nc_total0 = (nc_eta(1)+deg(1))*nc_eta(2)*nc_eta(3)
  nc_total1 = (nc_eta(1)+deg(1)-1)*nc_eta(2)*nc_eta(3)
  ! Time loop
  delta_t = 0.01_f64
  nsteps = 12
  time = 0.0_f64

  if( deg(1) == 2 ) then
     boundary = [ sll_p_boundary_clamped_square, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else if( deg(1) == 3 ) then
     boundary = [ sll_p_boundary_clamped_cubic, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else
     boundary = [ sll_p_boundary_clamped, sll_p_boundary_periodic, sll_p_boundary_periodic]
  end if
 
  ! Initialise maxwell_clamped FEM object
  call maxwell_solver%init(domain, nc_eta, deg,boundary )!,  boundary_func1= cos_3k, boundary_func2=zero, boundary_func3=cos_3k )

  call sll_s_spline_pp_init_3d(spline_pp_0, deg, nc_eta, boundary)
  call sll_s_spline_pp_init_3d(spline_pp_11, [deg(1)-1,deg(2),deg(3)], nc_eta, boundary)
  call sll_s_spline_pp_init_3d(spline_pp_12, [deg(1),deg(2)-1,deg(3)], nc_eta, boundary)
  call sll_s_spline_pp_init_3d(spline_pp_13, [deg(1),deg(2),deg(3)-1], nc_eta, boundary)
  call sll_s_spline_pp_init_3d(spline_pp_21, [deg(1),deg(2)-1,deg(3)-1], nc_eta, boundary)
  call sll_s_spline_pp_init_3d(spline_pp_22, [deg(1)-1,deg(2),deg(3)-1], nc_eta, boundary)
  call sll_s_spline_pp_init_3d(spline_pp_23, [deg(1)-1,deg(2)-1,deg(3)], nc_eta, boundary)

  allocate( x(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate( y(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate( z(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate( efield(1:nc_total1+nc_total0*2) )
  allocate( bfield(1:nc_total0+nc_total1*2) )
  allocate( efield1(1:nc_total1+nc_total0*2) )
  allocate( bfield1(1:nc_total0+nc_total1*2) )
  allocate( efield_ref(1:nc_total*3) )
  allocate( bfield_ref(1:nc_total*3) )
  allocate( efield_val(1:nc_total*3) )
  allocate( bfield_val(1:nc_total*3) )

  allocate( rho( nc_total0 ) )
  allocate( rho_ref( nc_total0 ) )
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


  ! Poisson problem
  call maxwell_solver%compute_rhs_from_function( 0, 1, rho, rho_k )
  call maxwell_solver%compute_e_from_rho( rho, efield )

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
           efield_ref(ind) = E_1k([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind+nc_total) = E_2k([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind+nc_total*2) = E_3k([x(i,j,k), y(i,j,k), z(i,j,k)])
           ind = ind+1
        end do
     end do
  end do
  error(1) = maxval(abs(efield_val-efield_ref))
  print*, 'Error Poisson:', error(1)

!!$  call sll_s_plot_two_fields_1d('Poissonx',nc_total, efield_val(1:nc_total),efield_ref(1:nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Poissony',nc_total, efield_val(nc_total+1:2*nc_total),efield_ref(nc_total+1:2*nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Poissonz',nc_total, efield_val(2*nc_total+1:3*nc_total),efield_ref(2*nc_total+1:3*nc_total),0,0.0_f64)



  time = 0.0_f64
  call maxwell_solver%L2projection( 1, 0, efield, cos_k, zero, zero )

  error(2) = abs(maxwell_solver%inner_product( efield(1:nc_total1), efield(1:nc_total1), 1, 1 ) - 4._f64*sll_p_pi**3)
  print*, 'Error in L2 norm squared:', error(2)


  ! Solve the two equations
  ! \partial_t E = \curl B
  ! \partial_t B = -\curl E
  ! on a staggered time grid
  ! We use the solution
  ! E(x,t) =  \begin{pmatrix} \cos(x_1) \sin(x_2) \sin(x_3) \sin( \sqrt{3} t) \\ -2\sin(x_1) \cos(x_2) \sin(x_3)  \sin( \sqrt{3} t) \\ \sin(x_1) \sin(x_2)\cos(x_3)  \sin( \sqrt{3} t) \end{pmatrix}
  ! B(x,t) = \begin{pmatrix} \sqrt{3} \sin(x_1) \cos(x_2) \cos(x_3) \cos( \sqrt{3} t) \\ 0 \\ -\sqrt{3} \cos(x_1) \cos(x_2) \sin(x_3) \cos( \sqrt{3} t) \end{pmatrix}


  ! Now assemble initial efield and bfield
  time = 0._f64
  !call maxwell_solver%L2projection( 1, 0, efield, E_1k, E_2k, E_3k )
  !efield(nc_total1+1:nc_total1+nc_total0) = 2.0_f64 *  efield(nc_total1+1:nc_total1+nc_total0)
  efield = 0._f64
  call maxwell_solver%L2projection( 2, 0, bfield, B_1k, zero, B_3k )


  ! Time stepping
  do istep = 1, nsteps
     call maxwell_solver%compute_field_energy( efield, bfield, energy(1) )
     efield1 = efield
     bfield1 = bfield

!!$     call maxwell_solver%compute_b_from_e( 0.5_f64*delta_t, efield, bfield )
!!$     call maxwell_solver%compute_e_from_b(         delta_t, bfield, efield )
!!$     call maxwell_solver%compute_b_from_e( 0.5_f64*delta_t, efield, bfield )
     call maxwell_solver%compute_curl_part( delta_t, efield, bfield)


!!$
     call maxwell_solver%compute_field_energy( efield, bfield, energy(2) )
     bfield1 = 0.5_f64*( bfield1 + bfield )
     efield1 = 0.5_f64*( efield1 + efield )
!!$
!!$     print*, 'error total energy', abs(energy(2)-energy(1))
  end do

  ! Evaluate E and B at the grid points
  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_21, nc_eta, bfield(1:nc_total0),  &
       bfield_val(1:nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_22, nc_eta, bfield(1+nc_total0:nc_total0+nc_total1),  &
       bfield_val(1+nc_total:2*nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped ( spline_pp_23, nc_eta , bfield(1+nc_total0+nc_total1:nc_total0+nc_total1*2),  &
       bfield_val(1+2*nc_total:3*nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_11, nc_eta, efield(1:nc_total1),  &
       efield_val(1:nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_12, nc_eta, efield(1+nc_total1:nc_total1+nc_total0),  &
       efield_val(1+nc_total:2*nc_total))
  call sll_s_spline_evaluate_basis_b_form_3d_clamped ( spline_pp_13, nc_eta , efield(1+nc_total1+nc_total0:nc_total1+nc_total0*2),  &
       efield_val(1+2*nc_total:3*nc_total))

  error1 = 0._f64
  ! Reference solutions
  time = real(nsteps,f64)*delta_t
  ind = 1
  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           efield_ref(ind) = E_1k([x(i,j,k), y(i,j,k), z(i,j,k)]) *sin(sqrt(3._f64)*time) 
           efield_ref(ind+nc_total) = -2.0_f64*E_2k([x(i,j,k), y(i,j,k), z(i,j,k)]) *sin(sqrt(3._f64)*time) 
           efield_ref(ind+nc_total*2) = E_3k([x(i,j,k), y(i,j,k), z(i,j,k)])*sin(sqrt(3._f64)*time)

           error1(1)=error1(1)+(efield_val(ind)-efield_ref(ind))**2
           error1(2)=error1(2)+(efield_val(nc_total+ind)-efield_ref(nc_total+ind))**2
           error1(3)=error1(3)+(efield_val(2*nc_total+ind)-efield_ref(2*nc_total+ind))**2
           ind = ind+1
        end do
     end do
  end do
  error1(1:3) = sqrt(error1(1:3)/real(nc_total,f64))
  print*, 'L2 Error efield', error1(1:3)
  
  error1 = 0._f64
  ind = 1
  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           bfield_ref(ind) = B_1k([x(i,j,k), y(i,j,k), z(i,j,k)])*cos(sqrt(3._f64)*time) 
           bfield_ref(ind+nc_total) = 0.0_f64
           bfield_ref(ind+nc_total*2) = B_3k([x(i,j,k), y(i,j,k), z(i,j,k)])*cos(sqrt(3._f64)*time)
           error1(1)=error1(1)+(bfield_val(ind)-bfield_ref(ind))**2
           error1(2)=error1(2)+(bfield_val(nc_total+ind)-bfield_ref(nc_total+ind))**2
           error1(3)=error1(3)+(bfield_val(2*nc_total+ind)-bfield_ref(2*nc_total+ind))**2
           ind = ind+1
        end do
     end do
  end do
  error1(1:3) = sqrt(error1(1:3)/real(nc_total,f64))
  print*, 'L2 Error bfield', error1(1:3)

  error(3) = maxval(abs(efield_val-efield_ref))
  error(4) = maxval(abs(bfield_val-bfield_ref))
  print*, 'Error efield:', error(3)
  print*, 'Error bfield:', error(4)

!!$  call sll_s_plot_two_fields_1d('Ex3d',nc_total, efield_val(1:nc_total),efield_ref(1:nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Ey3d',nc_total, efield_val(nc_total+1:2*nc_total),efield_ref(nc_total+1:2*nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Ez3d',nc_total, efield_val(2*nc_total+1:3*nc_total),efield_ref(2*nc_total+1:3*nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Bx3d',nc_total, bfield_val(1:nc_total),bfield_ref(1:nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('By3d',nc_total, bfield_val(nc_total+1:2*nc_total),bfield_ref(nc_total+1:2*nc_total),0,0.0_f64)
!!$  call sll_s_plot_two_fields_1d('Bz3d',nc_total, bfield_val(2*nc_total+1:3*nc_total),bfield_ref(2*nc_total+1:3*nc_total),0,0.0_f64)

  time = 0._f64
  ! Test compute_e_from_j
  call maxwell_solver%L2projection( 1, 1, efield(1:nc_total1), cos_k )


  current = 0._f64
  call maxwell_solver%compute_rhs_from_function( 1, 1, current(1:nc_total1), sin_k )

  call maxwell_solver%compute_E_from_j( current(1:nc_total1), efield(1:nc_total1), 1 )
  call sll_s_spline_evaluate_basis_b_form_3d_clamped( spline_pp_11, nc_eta, efield(1:nc_total1),  &
       efield_val(1:nc_total))
  ! Reference solution
  ind = 1
  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           efield_ref(ind) = cos_k([x(i,j,k), y(i,j,k), z(i,j,k)]) - sin_k([x(i,j,k), y(i,j,k), z(i,j,k)])
           ind = ind+1
        end do
     end do
  end do
  error(5) = maxval(abs(efield_val(1:nc_total)-efield_ref(1:nc_total)))
  print*, 'Error compute_e_from_j:', error(5)

!!$  call sll_s_plot_two_fields_1d('E_from_j',nc_total, efield_val(1:nc_total),efield_ref(1:nc_total),0,0.0_f64)

  ! Test compute_rho_from_e
  call maxwell_solver%compute_rhs_from_function( 0, 1, rho_ref, rho_k )
  call maxwell_solver%L2projection( 1, 0, efield, E_1k, E_2k, E_3k )
  call maxwell_solver%compute_rho_from_E( efield, rho )
  error(6) =  maxval( abs( rho - rho_ref ) )
  print*, 'Error compute_rho_from_e:', error(6)

  !call sll_s_plot_two_fields_1d('rho_from_E',nc_total0, rho, rho_ref, 0,0.0_f64)


  ! Clean up
  call maxwell_solver%free()
  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(efield)
  deallocate(bfield)
  deallocate(efield_ref)
  deallocate(bfield_ref)
  deallocate(efield_val)
  deallocate(bfield_val)

  
  if ( error(1) < 5.0d-3 .AND. error(2) < 8.0d-5 .AND. error(3) < 6.d-4 .AND. error(4) < 7.d-3 .AND. error(5)<4.0d-3 .AND. error(6)<3.0d-14) then
     print*, 'PASSED.'
  else
     print*, 'FAILED.'
  end if
  call sll_s_set_time_mark( end )
  write(*, "(A, F10.3)") "Main part run time [s] = ", sll_f_time_elapsed_between( start, end)
contains
  function rho_k(x)
    sll_real64             :: rho_k
    sll_real64, intent(in) :: x(3)

    rho_k = -3._f64*sin(x(1))*sin(x(2))*sin(x(3))
  end function rho_k
  function phi_k(x)
    sll_real64             :: phi_k
    sll_real64, intent(in) :: x(3)

    phi_k = - sin(x(1))*sin(x(2))*sin(x(3))
  end function phi_k

  function E_1k(x)
    sll_real64             :: E_1k
    sll_real64, intent(in) :: x(3)

    E_1k = cos(x(1))*sin(x(2))*sin(x(3)) 
  end function E_1k

  function E_2k(x)
    sll_real64             :: E_2k
    sll_real64, intent(in) :: x(3)

    E_2k = sin(x(1))*cos(x(2))*sin(x(3))
  end function E_2k

  function E_3k(x)
    sll_real64             :: E_3k
    sll_real64, intent(in) :: x(3)

    E_3k = sin(x(1))*sin(x(2))*cos(x(3))
  end function E_3k

  function B_1k(x)
    sll_real64             :: B_1k
    sll_real64, intent(in) :: x(3)

    B_1k = sqrt(3._f64)*sin(x(1))*cos(x(2))*cos(x(3))
  end function B_1k

  function B_3k(x)
    sll_real64             :: B_3k
    sll_real64, intent(in) :: x(3)

    B_3k = -sqrt(3._f64)*cos(x(1))*cos(x(2))*sin(x(3))
  end function B_3k


  function E_k(x)
    sll_real64             :: E_k
    sll_real64, intent(in) :: x(3)

    E_k = alpha/kx*cos(kx*x(1)) 
  end function E_k

  function cos_k(x)
    sll_real64             :: cos_k
    sll_real64, intent(in) :: x(3)

    cos_k = cos(x(1)+x(2)+x(3))
  end function cos_k


  function sin_k(x)
    sll_real64             :: sin_k
    sll_real64, intent(in) :: x(3)

    sin_k = sin(x(1)+x(2)+x(3))
  end function sin_k



  function zero(x)
    sll_real64             :: zero
    sll_real64, intent(in) :: x(3)

    zero=0._f64
  end function zero

  function constant(x)
    sll_real64             :: constant
    sll_real64, intent(in) :: x(3)

    constant=5._f64
  end function constant


  function cos_l(x)
    sll_real64             :: cos_l
    sll_real64, intent(in) :: x(3)

    cos_l = cos(kx*x(1)) 
  end function cos_l
  function sin_l(x)
    sll_real64             :: sin_l
    sll_real64, intent(in) :: x(3)

    sin_l = sin(kx*x(1)) 
  end function sin_l


end program test_maxwell_clamped_3d_fem
