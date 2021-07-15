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
program test_maxwell_3d_fem
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
    sll_p_pi

  use sll_m_maxwell_3d_fem, only: &
       sll_t_maxwell_3d_fem
  
  use sll_m_constants, only: sll_p_pi, sll_p_twopi

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !  type(sll_t_arbitrary_degree_spline_1d), pointer :: aspl
  !  sll_real64, dimension(:), allocatable :: knots

  sll_real64 :: eta1_max, eta1_min
  sll_real64 :: delta_eta(3)

  sll_int32  :: nc_eta(3), nc_total

  type(sll_t_maxwell_3d_fem)  :: maxwell_3d
  sll_real64, allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
  sll_real64, allocatable :: efield(:), bfield(:)
  sll_real64, allocatable :: efield_ref(:), bfield_ref(:)
  sll_real64, allocatable :: efield_val(:), bfield_val(:)
  sll_real64, allocatable :: rho(:), rho_ref(:)
  
  sll_int32                               :: i, j, k, istep, nsteps, ind
  sll_real64                                :: w1, w2
  sll_real64                                :: time
  sll_real64                                :: delta_t
  
  sll_real64                              :: Lx(3)
  sll_real64, dimension(3,2)              :: domain
  sll_int32                               :: deg(3)

  sll_real64                              :: error(6)
  sll_real64                              :: energy(2)

  ! Define computational domain
  eta1_min = .0_f64; eta1_max = 2.0_f64*sll_p_pi
  nc_eta = [16, 8, 32]
  nc_total = product(nc_eta)
  Lx(1) = eta1_max-eta1_min
  Lx(2) = Lx(1); Lx(3) = Lx(2)
  delta_eta = Lx/real(nc_eta,f64)
  domain(1,:) = [eta1_min, eta1_max]
  domain(2,:) = [eta1_min, eta1_max]
  domain(3,:) = [eta1_min, eta1_max]
  ! Set spline degree of 0-forms
  deg = [3,3,3]
  ! Time loop
  delta_t = 0.01_f64
  nsteps = 3
  w1 = sqrt(3.0_f64)
  w2 = + sqrt(3.0_f64)
  time = 0.0_f64
  ! Initialise maxwell FEM object
  call maxwell_3d%init(domain, nc_eta, deg)


  allocate(x(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(y(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(z(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(efield(1:nc_total*3))
  allocate(bfield(1:nc_total*3))
  allocate(efield_ref(1:nc_total*3))
  allocate(bfield_ref(1:nc_total*3))
  allocate(efield_val(1:nc_total*3))
  allocate(bfield_val(1:nc_total*3))
  allocate( rho(1:nc_total) )
  allocate( rho_ref( 1:nc_total) )
   
  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           x(i,j,k) = eta1_min + real(i-1,f64)*delta_eta(1)
           y(i,j,k) = eta1_min + real(j-1,f64)*delta_eta(2)
           z(i,j,k) = eta1_min + real(k-1,f64)*delta_eta(3)
        end do
     end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now assemble initial efield and bfield
  ! B starts at time -0.5\Delta t
  time = -0.5_f64*delta_t
   
  call maxwell_3d%L2projection( 2, 1, bfield(1:nc_total), cos_k )
  bfield(1:nc_total) =  w2 *  bfield(1:nc_total)
  bfield(nc_total+1:nc_total*2) = 0.0_f64
  call maxwell_3d%L2projection( 2, 3, bfield(nc_total*2+1:nc_total*3), cos_k )
  bfield(nc_total*2+1:nc_total*3) = - w2 *  bfield(nc_total*2+1:nc_total*3)

  ! E starts at time 0
  time = 0.0_f64

  call maxwell_3d%L2projection( 1, 1, efield(1:nc_total), cos_k )
  call maxwell_3d%L2projection( 1, 2, efield(nc_total+1:nc_total*2), cos_k )
  call maxwell_3d%L2projection( 1, 3, efield(nc_total*2+1:nc_total*3), cos_k )
  efield(nc_total+1:nc_total*2) = - 2.0_f64 *  efield(nc_total+1:nc_total*2)

  call maxwell_3d%compute_field_energy( efield, bfield, energy(1) )
  print*, 'total energy1', energy(1)
  
  ! Time stepping
  do istep = 1, nsteps
     call maxwell_3d%compute_curl_part( delta_t, efield, bfield )
!!$     call maxwell_3d%compute_B_from_E( delta_t, efield, bfield )
!!$     call maxwell_3d%compute_E_from_B( delta_t, bfield, efield )
  end do

  ! Evaluate E and B at the grid points
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)-1], bfield(1:nc_total), bfield_val(1:nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)-1], bfield(1+nc_total:2*nc_total), bfield_val(1+nc_total:2*nc_total)  )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2)-1,deg(3)], bfield(1+nc_total*2:3*nc_total), bfield_val(1+nc_total*2:3*nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)], efield(1:nc_total), efield_val(1:nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)], efield(1+nc_total:2*nc_total), efield_val(1+nc_total:2*nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2),deg(3)-1], efield(1+nc_total*2:3*nc_total), efield_val(1+nc_total*2:3*nc_total) )
  

  ! Reference solutions
  time = real(nsteps,f64)*delta_t
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           efield_ref(ind) = cos_k([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind+nc_total) = -2.0_f64*efield_ref(ind)
           efield_ref(ind+nc_total*2) = efield_ref(ind)
           ind = ind+1
        end do
     end do
  end do
  time = (real(nsteps,f64)-0.5_f64)*delta_t
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           bfield_ref(ind) = w2*cos_k([x(i,j,k), y(i,j,k), z(i,j,k)])
           bfield_ref(ind+nc_total) = 0.0_f64
           bfield_ref(ind+nc_total*2) = -bfield_ref(ind)
           ind = ind+1
        end do
     end do
  end do
  error(3) = maxval(abs(efield_val-efield_ref))
  error(4) = maxval(abs(bfield_val-bfield_ref))
  print*, 'Error efield:', error(3)
  print*, 'Error bfield:', error(4)

  call maxwell_3d%compute_field_energy( efield, bfield, energy(2) )
  print*, 'total energy2', energy(2)
  print*, 'error total energy', abs(energy(1)-energy(2))
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Poisson problem
  call maxwell_3d%compute_rhs_from_function( 0, 1, rho, cos_k )
  call maxwell_3d%compute_e_from_rho( rho, efield )
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
           efield_ref(ind) = sin_k([x(i,j,k), y(i,j,k), z(i,j,k)])/3.0_f64
           efield_ref(ind+nc_total) = efield_ref(ind)
           efield_ref(ind+nc_total*2) = efield_ref(ind)
           ind = ind+1
        end do
     end do
  end do
  error(1) = maxval(abs(efield_val-efield_ref))
  print*, 'Error Poisson:', error(1)

  ! Solve the two equations
  ! \partial_t E = \curl B
  ! \partial_t B = -\curl E
  ! on a staggered time grid
  ! We use the solution
  ! E(x,t) =  \begin{pmatrix} \cos(x_1+x_2+x_3 - \sqrt{3} t) \\ -2\cos(x_1+x_2+x_3) \\ \cos(x_1+x_2+x_3 - \sqrt{3} t) \end{pmatrix}
  ! B(x,t) = \begin{pmatrix} \sqrt{3} \cos(x_1+x_2+x_3 - \sqrt{3} t) \\ 0 \\ -\sqrt{3} \cos(x_1+x_2+x_3 - \sqrt{3} t) \end{pmatrix}

  
  ! Now assemble initial efield and bfield
  ! B starts at time -0.5\Delta t
  time = -0.5_f64*delta_t
   
  call maxwell_3d%L2projection( 2, 1, bfield(1:nc_total), cos_k )
  bfield(1:nc_total) =  w2 *  bfield(1:nc_total)
  bfield(nc_total+1:nc_total*2) = 0.0_f64
  call maxwell_3d%L2projection( 2, 3, bfield(nc_total*2+1:nc_total*3), cos_k )
  bfield(nc_total*2+1:nc_total*3) = - w2 *  bfield(nc_total*2+1:nc_total*3)

  ! E starts at time 0
  time = 0.0_f64

  call maxwell_3d%L2projection( 1, 1, efield(1:nc_total), cos_k )
  call maxwell_3d%L2projection( 1, 2, efield(nc_total+1:nc_total*2), cos_k )
  call maxwell_3d%L2projection( 1, 3, efield(nc_total*2+1:nc_total*3), cos_k )
  efield(nc_total+1:nc_total*2) = - 2.0_f64 *  efield(nc_total+1:nc_total*2)
  error(2) = abs(maxwell_3d%inner_product( efield(1:nc_total), efield(1:nc_total), 1, 1 ) - 4.0_f64*sll_p_pi**3)
  print*, 'Error in L2 norm squared:', error(2)

  ! Time stepping
  do istep = 1, nsteps
     call maxwell_3d%compute_B_from_E( delta_t, efield, bfield )
     call maxwell_3d%compute_E_from_B( delta_t, bfield, efield )
  end do

  ! Evaluate E and B at the grid points
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)-1], bfield(1:nc_total), bfield_val(1:nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)-1], bfield(1+nc_total:2*nc_total), bfield_val(1+nc_total:2*nc_total)  )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2)-1,deg(3)], bfield(1+nc_total*2:3*nc_total), bfield_val(1+nc_total*2:3*nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)], efield(1:nc_total), efield_val(1:nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)], efield(1+nc_total:2*nc_total), efield_val(1+nc_total:2*nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2),deg(3)-1], efield(1+nc_total*2:3*nc_total), efield_val(1+nc_total*2:3*nc_total) )
  

  ! Reference solutions
  time = real(nsteps,f64)*delta_t
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           efield_ref(ind) = cos_k([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind+nc_total) = -2.0_f64*efield_ref(ind)
           efield_ref(ind+nc_total*2) = efield_ref(ind)
           ind = ind+1
        end do
     end do
  end do
  time = (real(nsteps,f64)-0.5_f64)*delta_t
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           bfield_ref(ind) = w2*cos_k([x(i,j,k), y(i,j,k), z(i,j,k)])
           bfield_ref(ind+nc_total) = 0.0_f64
           bfield_ref(ind+nc_total*2) = -bfield_ref(ind)
           ind = ind+1
        end do
     end do
  end do
  error(3) = maxval(abs(efield_val-efield_ref))
  error(4) = maxval(abs(bfield_val-bfield_ref))
  print*, 'Error efield:', error(3)
  print*, 'Error bfield:', error(4)


  ! Test compute_e_from_j
  call maxwell_3d%L2projection( 1, 1, efield(1:nc_total), cos_k )
  call maxwell_3d%compute_rhs_from_function( 1, 1, rho, sin_k )

  call maxwell_3d%compute_E_from_j( rho, efield(1:nc_total), 1 )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)], efield(1:nc_total),  &
       efield_val(1:nc_total))
  ! Reference solution
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           efield_ref(ind) = cos_k([x(i,j,k), y(i,j,k), z(i,j,k)])-&
                sin_k([x(i,j,k), y(i,j,k), z(i,j,k)])
           ind = ind+1
        end do
     end do
  end do
  error(5) = maxval(abs(efield_val(1:nc_total)-efield_ref(1:nc_total)))
  print*, 'Error compute_e_from_j:', error(5)

  ! Test compute_rho_from_e
  call maxwell_3d%compute_rhs_from_function( 0, 1, rho_ref, cos_k )
  rho_ref = 3.0_f64*rho_ref
  call maxwell_3d%L2projection( 1, 1, efield(1:nc_total), sin_k )
  call maxwell_3d%L2projection( 1, 2, efield(nc_total+1:nc_total*2), sin_k )
  call maxwell_3d%L2projection( 1, 3, efield(nc_total*2+1:nc_total*3), sin_k )

  call maxwell_3d%compute_rho_from_E( efield, rho )
  
  error(6) =  maxval( abs( rho - rho_ref ) )
  print*, 'Error compute_rho_from_e:', error(6)
 
  

  ! Clean up
  call maxwell_3d%free()
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


  if ( error(1) < 5.5d-4 .AND. error(2) < 3.6d-5 .AND. error(3) < 3.3d-3 .AND. error(4) < 3.1d-3 .AND. error(5)<7.3d-4 .AND. error(6)<6.6d-8) then
     print*, 'PASSED.'
  else
     print*, 'FAILED.'
  end if
  !second part with different spline degrees
  nc_eta = [16, 64, 32]
  nc_total = product(nc_eta)
  delta_eta = Lx/real(nc_eta,f64)
  domain(1,:) = [eta1_min, eta1_max]
  domain(2,:) = [eta1_min, eta1_max]
  domain(3,:) = [eta1_min, eta1_max]
  ! Set spline degree of 0-forms
  deg = [3,2,3]
  ! Time loop
  delta_t = 0.01_f64
  time = 0.0_f64
  ! Initialise maxwell FEM object
  call maxwell_3d%init(domain, nc_eta, deg)


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
  call maxwell_3d%compute_rhs_from_function( 0, 1, rho, cos_k )
  call maxwell_3d%compute_e_from_rho( rho, efield )
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
           efield_ref(ind) = sin_k([x(i,j,k), y(i,j,k), z(i,j,k)])/3.0_f64
           efield_ref(ind+nc_total) = efield_ref(ind)
           efield_ref(ind+nc_total*2) = efield_ref(ind)
           ind = ind+1
        end do
     end do
  end do
  error(1) = maxval(abs(efield_val-efield_ref))
  print*, 'Error Poisson:', error(1)

  ! Solve the two equations
  ! \partial_t E = \curl B
  ! \partial_t B = -\curl E
  ! on a staggered time grid
  ! We use the solution
  ! E(x,t) =  \begin{pmatrix} \cos(x_1+x_2+x_3 - \sqrt{3} t) \\ -2\cos(x_1+x_2+x_3) \\ \cos(x_1+x_2+x_3 - \sqrt{3} t) \end{pmatrix}
  ! B(x,t) = \begin{pmatrix} \sqrt{3} \cos(x_1+x_2+x_3 - \sqrt{3} t) \\ 0 \\ -\sqrt{3} \cos(x_1+x_2+x_3 - \sqrt{3} t) \end{pmatrix}

  
  ! Now assemble initial efield and bfield
  ! B starts at time -0.5\Delta t
  time = -0.5_f64*delta_t
   
  call maxwell_3d%L2projection( 2, 1, bfield(1:nc_total), cos_k )
  bfield(1:nc_total) =  w2 *  bfield(1:nc_total)
  bfield(nc_total+1:nc_total*2) = 0.0_f64
  call maxwell_3d%L2projection( 2, 3, bfield(nc_total*2+1:nc_total*3), cos_k )
  bfield(nc_total*2+1:nc_total*3) = - w2 *  bfield(nc_total*2+1:nc_total*3)

  ! E starts at time 0
  time = 0.0_f64

  call maxwell_3d%L2projection( 1, 1, efield(1:nc_total), cos_k )
  call maxwell_3d%L2projection( 1, 2, efield(nc_total+1:nc_total*2), cos_k )
  call maxwell_3d%L2projection( 1, 3, efield(nc_total*2+1:nc_total*3), cos_k )
  efield(nc_total+1:nc_total*2) = - 2.0_f64 *  efield(nc_total+1:nc_total*2)
  error(2) = abs(maxwell_3d%inner_product( efield(1:nc_total), efield(1:nc_total), 1, 1 ) - 4.0_f64*sll_p_pi**3)
  print*, 'Error in L2 norm squared:', error(2)

  ! Time stepping
  do istep = 1, nsteps
     call maxwell_3d%compute_B_from_E( delta_t, efield, bfield )
     call maxwell_3d%compute_E_from_B( delta_t, bfield, efield )
  end do

  ! Evaluate E and B at the grid points
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)-1], bfield(1:nc_total), bfield_val(1:nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)-1], bfield(1+nc_total:2*nc_total), bfield_val(1+nc_total:2*nc_total)  )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2)-1,deg(3)], bfield(1+nc_total*2:3*nc_total), bfield_val(1+nc_total*2:3*nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)], efield(1:nc_total), efield_val(1:nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)], efield(1+nc_total:2*nc_total), efield_val(1+nc_total:2*nc_total) )
  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2),deg(3)-1], efield(1+nc_total*2:3*nc_total), efield_val(1+nc_total*2:3*nc_total) )
  

  ! Reference solutions
  time = real(nsteps,f64)*delta_t
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           efield_ref(ind) = cos_k([x(i,j,k), y(i,j,k), z(i,j,k)])
           efield_ref(ind+nc_total) = -2.0_f64*efield_ref(ind)
           efield_ref(ind+nc_total*2) = efield_ref(ind)
           ind = ind+1
        end do
     end do
  end do
  time = (real(nsteps,f64)-0.5_f64)*delta_t
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           bfield_ref(ind) = w2*cos_k([x(i,j,k), y(i,j,k), z(i,j,k)])
           bfield_ref(ind+nc_total) = 0.0_f64
           bfield_ref(ind+nc_total*2) = -bfield_ref(ind)
           ind = ind+1
        end do
     end do
  end do
  error(3) = maxval(abs(efield_val-efield_ref))
  error(4) = maxval(abs(bfield_val-bfield_ref))
  print*, 'Error efield:', error(3)
  print*, 'Error bfield:', error(4)


  ! Test compute_e_from_j
  call maxwell_3d%L2projection( 1, 1, efield(1:nc_total), cos_k )
  call maxwell_3d%compute_rhs_from_function( 1, 1, rho, sin_k )

  call maxwell_3d%compute_E_from_j( rho, efield(1:nc_total), 1 )
  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)], efield(1:nc_total),  &
       efield_val(1:nc_total))
  ! Reference solution
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           efield_ref(ind) = cos_k([x(i,j,k), y(i,j,k), z(i,j,k)])-&
                sin_k([x(i,j,k), y(i,j,k), z(i,j,k)])
           ind = ind+1
        end do
     end do
  end do
  error(5) = maxval(abs(efield_val(1:nc_total)-efield_ref(1:nc_total)))
  print*, 'Error compute_e_from_j:', error(5)

  ! Test compute_rho_from_e
  call maxwell_3d%compute_rhs_from_function( 0, 1, rho_ref, cos_k )
  rho_ref = 3.0_f64*rho_ref
  call maxwell_3d%L2projection( 1, 1, efield(1:nc_total), sin_k )
  call maxwell_3d%L2projection( 1, 2, efield(nc_total+1:nc_total*2), sin_k )
  call maxwell_3d%L2projection( 1, 3, efield(nc_total*2+1:nc_total*3), sin_k )

  call maxwell_3d%compute_rho_from_E( efield, rho )
  
  error(6) =  maxval( abs( rho - rho_ref ) )
  print*, 'Error compute_rho_from_e:', error(6)
 
  

  ! Clean up
  call maxwell_3d%free()
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

  if ( error(1) < 5.5d-4 .AND. error(2) < 3.6d-5 .AND. error(3) < 3.3d-3 .AND. error(4) < 3.1d-3 .AND. error(5)<7.3d-4 .AND. error(6)<6.6d-8) then
     print*, 'PASSED.'
  else
     print*, 'FAILED.'
  end if

  
contains
  function cos_k(x)
    sll_real64             :: cos_k
    sll_real64, intent(in) :: x(3)

    cos_k = cos((x(1)+x(2)+x(3))-w1*time) 
  end function cos_k


  function sin_k(x)
    sll_real64             :: sin_k
    sll_real64, intent(in) :: x(3)

    sin_k = sin((x(1)+x(2)+x(3))-w1*time) 
  end function sin_k

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
end program test_maxwell_3d_fem
