!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations with coordinate transformation in 3D
!> @details
!> 
!> @author
!> Benedikt Perse

module sll_m_maxwell_3d_trafo
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
       sll_f_get_collective_rank, &
       sll_v_world_collective

  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_linear_operator_block, only : &
       sll_t_linear_operator_block

  use sll_m_linear_operator_penalized, only : &
       sll_t_linear_operator_penalized

  use sll_m_linear_operator_poisson_3d, only : &
       sll_t_linear_operator_poisson_3d

  use sll_m_linear_operator_schur_eb_3d, only : &
       sll_t_linear_operator_schur_eb_3d

  use sll_m_linear_solver_cg, only : &
       sll_t_linear_solver_cg

  use sll_m_linear_solver_mgmres, only : &
       sll_t_linear_solver_mgmres

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_maxwell_3d_base, only: &
       sll_i_function_3d_real64, &
       sll_c_maxwell_3d_base

  use sll_m_preconditioner_fft, only : &
       sll_t_preconditioner_fft

  use sll_m_profile_functions, only: &
       sll_t_profile_functions

  use sll_m_spline_fem_utilities, only : &
       sll_s_multiply_g, &
       sll_s_multiply_gt, &
       sll_s_multiply_c, &
       sll_s_multiply_ct

  use sll_m_spline_fem_utilities_3d, only: &
       sll_s_spline_fem_mass3d, &
       sll_s_spline_fem_mixedmass3d

  implicit none

  public :: &
       sll_t_maxwell_3d_trafo

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_c_maxwell_3d_base) :: sll_t_maxwell_3d_trafo

     sll_real64, allocatable           :: work0(:)         !< scratch data
     sll_real64, allocatable           :: work(:)          !< scratch data
     sll_real64, allocatable           :: work2(:)         !< scratch data
     type(sll_t_matrix_csr)            :: mass0            !< mass matrix
     type(sll_t_matrix_csr)            :: mass1(3,3)       !< mass matrix
     type(sll_t_matrix_csr)            :: mass2(3,3)       !< mass matrix
     type(sll_t_linear_operator_block) :: mass1_operator   !< block mass matrix
     type(sll_t_linear_operator_block) :: mass2_operator   !< block mass matrix
     type(sll_t_linear_solver_cg) :: mass0_solver     !< mass matrix solver
     type(sll_t_linear_solver_cg) :: mass1_solver     !< mass matrix solver
     type(sll_t_linear_solver_cg) :: mass2_solver     !< mass matrix solver
     type(sll_t_preconditioner_fft) :: preconditioner_fft !< preconditioner for mass matrices

     type(sll_t_linear_operator_poisson_3d) :: poisson_matrix !< Poisson matrix
     type(sll_t_linear_operator_penalized)  :: poisson_operator !< Poisson matrix with constraint on constant vector
     type(sll_t_linear_solver_cg)  :: poisson_solver        !< CG solver to invert Poisson matrix
     type( sll_t_linear_operator_schur_eb_3d ) :: linear_op_schur_eb !< Schur complement operator for advect_eb
     type( sll_t_linear_solver_mgmres )        :: linear_solver_schur_eb !< Schur complement solver for advect_eb
     type(sll_t_mapping_3d), pointer    :: map              !< coordinate transformation

     logical :: adiabatic_electrons = .false. !< Set true if solver with adiabatic electrions

   contains

     procedure :: &
          compute_E_from_B => sll_s_compute_E_from_B_3d_trafo!< Solve E and B part of Amperes law with B constant in time
     procedure :: &
          compute_B_from_E => sll_s_compute_B_from_E_3d_trafo!< Solve Faraday equation with E constant in time
     procedure :: &
          compute_curl_part => sll_s_compute_curl_part_3d_trafo!< Solve curl part of Maxwell's equations
     procedure :: &
          compute_E_from_rho => sll_s_compute_E_from_rho_3d_trafo!< Solve E from rho using Gauss/Poisson
     procedure :: &
          compute_rho_from_E => sll_s_compute_rho_from_E_3d_trafo !< Compute rho from E 
     procedure :: &
          compute_E_from_j => sll_s_compute_E_from_j_3d_trafo !< Compute E from the current j
     procedure :: &
          compute_phi_from_rho => sll_s_compute_phi_from_rho_3d_trafo !< Compute phi from rho (by solving the quasi neutrality equation)
     procedure :: &
          compute_phi_from_j => sll_s_compute_phi_from_j_3d_trafo !< Compute phi from j (dynamic from of quasineutrality equation for adiabatic electrons)
     procedure :: &
          compute_rhs_from_function => sll_s_compute_rhs_trafo !< Compute integral over given function tested by the basis
     procedure :: &
          L2projection => L2projection_3d_trafo !< Compute L_2 projection of a given function
     procedure :: &
          L2norm_squared => L2norm_squared_3d_trafo !< Compute the square of the L2 norm of a given vector
     procedure :: &
          inner_product => inner_product_3d_trafo  !< Inner product of two dof-vectors with mass matrix
     procedure :: &
          init => init_3d_trafo !< Initialize the Maxwell class
     procedure :: &
          init_from_file => init_from_file_3d_trafo !< Initialize the Maxwell class with parameters read from nml-file
     procedure :: &
          free => free_3d_trafo !< Free Maxwell class
     procedure :: &
          multiply_g !< Multiplication with gradient matrix 
     procedure :: &
          multiply_gt !< Multiplication with transposed gradient matrix  
     procedure :: &
          multiply_c !< Multiplication with curl matrix
     procedure :: &
          multiply_ct !< Multiplication with transposed curl matrix
     procedure :: &
          multiply_mass => multiply_mass_trafo !< Product with the mass matrix
     procedure :: &
          multiply_mass_inverse => multiply_mass_inverse_trafo !< Invert mass matrix

  end type sll_t_maxwell_3d_trafo

contains


  !> compute E from B using weak Ampere formulation 
  subroutine sll_s_compute_e_from_b_3d_trafo( self, delta_t, field_in, field_out )
    class(sll_t_maxwell_3d_trafo) :: self         !< Maxwell solver class
    sll_real64, intent( in    )   :: delta_t      !< Time step
    sll_real64, intent( in    )   :: field_in(:)  !< B
    sll_real64, intent( inout )   :: field_out(:) !< E

    call multiply_mass_trafo( self, [2], field_in, self%work )

    call multiply_ct( self, self%work, self%work2 )

    call self%mass1_solver%solve( self%work2, self%work )
    ! Update e from self value
    field_out = field_out + delta_t*self%work

  end subroutine sll_s_compute_e_from_b_3d_trafo


  !> Compute B from E using strong 3D Faraday equation for spline coefficients
  subroutine sll_s_compute_b_from_e_3d_trafo( self, delta_t, field_in, field_out )
    class(sll_t_maxwell_3d_trafo) :: self         !< Maxwell solver class
    sll_real64, intent( in    )   :: delta_t      !< Time step
    sll_real64, intent( in    )   :: field_in(:)  !< E
    sll_real64, intent( inout )   :: field_out(:) !< B

    call multiply_c( self, field_in, self%work )
    ! Update b from self value
    field_out = field_out - delta_t * self%work

  end subroutine sll_s_compute_b_from_e_3d_trafo


  !> Solve curl part of Maxwell's equations
  subroutine sll_s_compute_curl_part_3d_trafo( self, delta_t, efield, bfield, betar )
    class(sll_t_maxwell_3d_trafo) :: self      !< Maxwell solver class
    sll_real64, intent( in    )   :: delta_t   !< Time step
    sll_real64, intent( inout )   :: efield(:) !< E
    sll_real64, intent( inout )   :: bfield(:) !< B
    sll_real64, optional          :: betar     !< 1/beta
    !local variables
    sll_real64 :: factor

    if( present(betar) ) then
       factor = betar
    else
       factor = 1._f64
    end if

    ! Compute C^T M2 b
    call self%mass2_operator%dot( bfield, self%work )
    call self%multiply_ct( self%work, self%work2 ) 

    self%linear_op_schur_eb%sign = -delta_t**2*factor*0.25_f64
    call self%linear_op_schur_eb%dot( efield, self%work )
    self%work = self%work + delta_t*factor*self%work2

    ! Save efield dofs from previous time step for B field update
    self%work2 = efield

    ! Invert Schur complement matrix
    self%linear_op_schur_eb%sign = delta_t**2*factor*0.25_f64
    call self%linear_solver_schur_eb%set_guess( efield )
    call self%linear_solver_schur_eb%solve( self%work, efield)

    ! Update B field
    self%work2 = self%work2 + efield
    call self%compute_b_from_e( delta_t*0.5_f64, self%work2, bfield)

  end subroutine sll_s_compute_curl_part_3d_trafo


  !> Compute E_i from rho_i integrated over the time interval using weak Poisson's equation
  subroutine sll_s_compute_e_from_rho_3d_trafo( self, field_in, field_out )  
    class(sll_t_maxwell_3d_trafo) :: self         !< Maxwell solver class
    sll_real64, intent( in    )   :: field_in(:)  !< rho
    sll_real64, intent(   out )   :: field_out(:) !< E

    call self%poisson_solver%solve( field_in, self%work0 )
    call multiply_g(self, self%work0, field_out)
    field_out = -field_out

  end subroutine sll_s_compute_e_from_rho_3d_trafo


  !> compute rho from e using weak Gauss law ( rho = G^T M_1 e ) 
  subroutine sll_s_compute_rho_from_e_3d_trafo( self, field_in, field_out ) 
    class(sll_t_maxwell_3d_trafo) :: self         !< Maxwell solver class
    sll_real64, intent( in    )   :: field_in(:)  !< E
    sll_real64, intent(   out )   :: field_out(:) !< rho

    call multiply_mass_trafo( self, [1], field_in, self%work )

    call multiply_gt( self, self%work, field_out )

    field_out = - field_out

  end subroutine sll_s_compute_rho_from_e_3d_trafo


  !> Compute E_i from j_i integrated over the time interval using weak Ampere formulation, delta_t is already included 
  subroutine sll_s_compute_e_from_j_3d_trafo( self, current, E, component )
    class(sll_t_maxwell_3d_trafo)         :: self       !< Maxwell solver class
    sll_real64, intent( in    )           :: current(:) !< Current integrated over time interval
    sll_real64, intent( inout )           :: E(:)       !< Updated electric field
    sll_int32,  intent( in    ), optional :: component  !< component

    call self%mass1_solver%solve( current, self%work )
    E = E - self%work

  end subroutine sll_s_compute_e_from_j_3d_trafo


  !> Compute phi from rho_i integrated over the time interval
  subroutine sll_s_compute_phi_from_rho_3d_trafo( self, field_in, field_out, efield_dofs )  
    class(sll_t_maxwell_3d_trafo) :: self         !< Maxwell solver class
    sll_real64, intent( in    )           :: field_in(:)  !< rho
    sll_real64, intent( inout )           :: field_out(:) !< phi
    sll_real64, intent(   out )           :: efield_dofs(:) !< E

    call self%mass0_solver%solve( field_in, field_out )
    call self%multiply_g( field_out, efield_dofs )
    efield_dofs = -efield_dofs

  end subroutine sll_s_compute_phi_from_rho_3d_trafo


  !> Compute phi from j_i integrated over the time interval, delta_t is already included 
  subroutine sll_s_compute_phi_from_j_3d_trafo( self, field_in, field_out, efield_dofs )
    class(sll_t_maxwell_3d_trafo) :: self       !< Maxwell solver class
    sll_real64, intent( in    )           :: field_in(:) !< Current integrated over time interval
    sll_real64, intent( inout )           :: field_out(:) !< phi
    sll_real64, intent(   out )           :: efield_dofs(:) !< E

    call self%multiply_gt( field_in, self%work0 ) 
    call self%mass0_solver%solve( self%work0, self%work2(1:self%n_total) )
    field_out = field_out + self%work2(1:self%n_total)

    call self%multiply_g( field_out, efield_dofs )
    efield_dofs = -efield_dofs

  end subroutine sll_s_compute_phi_from_j_3d_trafo


  !> Compute the FEM right-hand-side for a given function f and periodic splines of given degree
  !> Its components are $\int f(F(\xi)) N_i(F(\xi)) J_F d\xi$ where $N_i$ is the B-spline starting at $\xi_i$ 
  subroutine sll_s_compute_rhs_trafo( self, form,  component, coefs_dofs, func1, func2, func3 )
    class(sll_t_maxwell_3d_trafo)                  :: self          !< Maxwell solver class
    sll_int32,  intent( in    )                    :: form          !< Specify if the function is a0,1,2 or 3-form
    sll_int32,  intent( in    )                    :: component     !< Specify the component of the function
    sll_real64, intent(   out )                    :: coefs_dofs(:) !< Finite Element right-hand-side
    procedure(sll_i_function_3d_real64)            :: func1         !< Function first component
    procedure(sll_i_function_3d_real64), optional  :: func2         !< Function second component
    procedure(sll_i_function_3d_real64), optional  :: func3         !< Function third component
    ! local variables
    sll_int32 :: degree(3)

    ! Define the spline degree in the 3 dimensions, depending on form and component of the form
    if ( form == 0 ) then
       call rhs_zeroform( self, self%s_deg_0, func1, coefs_dofs )
    elseif (form == 1 ) then
       degree = self%s_deg_0
       degree(component) = self%s_deg_1(component)
       call rhs_oneform( self, degree, func1, func2, func3, component, coefs_dofs )
    elseif( form == 2) then
       degree =  self%s_deg_1
       degree(component) = self%s_deg_0(component)
       call rhs_twoform( self, degree, func1, func2, func3, component, coefs_dofs )
    elseif( form == 3) then
       call rhs_threeform( self, self%s_deg_1, func1, coefs_dofs )
    else 
       print*, 'Wrong form.'
    end if

  end subroutine sll_s_compute_rhs_trafo


  !> Compute $\int f(F(\xi)) N_i(\xi) J_F(\xi) d\xi$ for scalar function f, where $N_i$ is the B-spline starting at $\xi_i$, replace modulo(a,b) by a-floor(a/b)*b
  subroutine rhs_zeroform( self, deg, func, coefs_dofs )
    class(sll_t_maxwell_3d_trafo)       :: self          !< Maxwell solver class
    sll_int32,  intent( in    )         :: deg(3)           !< spline degree
    procedure(sll_i_function_3d_real64) :: func          !< function
    sll_real64, intent(   out )         :: coefs_dofs(:) !< Finite Element right-hand-side
    ! local variables
    sll_int32               :: i1, i2, i3, j1, j2, j3, k1, k2, k3, q(3), counter
    sll_real64              :: coef, jacobian, c(3)
    sll_real64, allocatable :: xw_gauss_d1(:,:), xw_gauss_d2(:,:), xw_gauss_d3(:,:)
    sll_real64, allocatable :: bspl_d1(:,:), bspl_d2(:,:), bspl_d3(:,:)
    sll_real64 :: scratch(maxval(deg+1))

    ! take enough Gauss points so that projection is exact for splines of degree deg
    q = deg+1
    ! rescale on [0,1] for compatibility with B-splines
    allocate(xw_gauss_d1(2, q(1)))
    allocate(bspl_d1(q(1), deg(1)+1 ))
    xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights(q(1), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k1 = 1, q(1)
       call sll_s_uniform_bsplines_eval_basis(deg(1),xw_gauss_d1(1,k1), scratch(1:deg(1)+1))
       bspl_d1(k1,:) = scratch(1:deg(1)+1)
    end do

    allocate(xw_gauss_d2(2, q(2)))
    allocate(bspl_d2(q(2), deg(2)+1))
    xw_gauss_d2 = sll_f_gauss_legendre_points_and_weights(q(2), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k2 = 1, q(2)
       call sll_s_uniform_bsplines_eval_basis(deg(2),xw_gauss_d2(1,k2), scratch(1:deg(2)+1))
       bspl_d2(k2,:) = scratch(1:deg(2)+1)
    end do

    allocate(xw_gauss_d3(2,q(3)))
    allocate(bspl_d3(q(3), deg(3)+1 ))
    xw_gauss_d3 = sll_f_gauss_legendre_points_and_weights(q(3), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k3 = 1, q(3)
       call sll_s_uniform_bsplines_eval_basis(deg(3),xw_gauss_d3(1,k3), scratch(1:deg(3)+1) )
       bspl_d3(k3,:) = scratch(1:deg(3)+1)
    end do


    counter = 1
    ! Compute coefs_dofs = int f(\xi)N_i(\xi) J_F(\xi) d\xi 
    do i3 = 1, self%n_dofs(3)
       do i2 = 1, self%n_dofs(2)
          do i1 = 1, self%n_dofs(1)
             coef = 0.0_f64
             ! loop over support of B spline
             do j3 = 1, deg(3)+1
                do j2 = 1, deg(2)+1
                   do j1 = 1, deg(1)+1
                      ! loop over Gauss points
                      do k3 = 1, q(3)
                         c(3) = self%delta_x(3) * (xw_gauss_d3(1,k3) + real(i3 + j3 - 2-((i3 + j3 - 2)/self%n_dofs(3)) * self%n_dofs(3),f64))
                         do k2 = 1, q(2)
                            c(2) = self%delta_x(2) * (xw_gauss_d2(1,k2) + real(i2 + j2 - 2-((i2 + j2 - 2)/self%n_dofs(2)) * self%n_dofs(2),f64))
                            do k1 = 1, q(1)
                               c(1) = self%delta_x(1) * (xw_gauss_d1(1,k1) + real(i1 + j1 - 2-((i1 + j1 - 2)/self%n_dofs(1)) * self%n_dofs(1),f64))
                               jacobian = self%map%jacobian( c ) 
                               coef = coef + xw_gauss_d1(2,k1) * xw_gauss_d2(2,k2) * xw_gauss_d3(2,k3) *&
                                    func(self%map%get_x(c) ) * &
                                    bspl_d1(k1,deg(1)+2-j1) *&
                                    bspl_d2(k2,deg(2)+2-j2) *&
                                    bspl_d3(k3,deg(3)+2-j3) * abs(jacobian)
                            enddo
                         enddo
                      end do
                   end do
                end do
             end do
             ! rescale by cell size
             coefs_dofs(counter) = coef*self%volume
             counter = counter+1
          enddo
       end do
    end do

  end subroutine rhs_zeroform


  !> Compute $\int f(F(\xi))\cdot DF^{-T}(\xi) N_i(\xi) J_F(\xi) d\xi$ for vector function f, where $N_i$ is the B-spline starting at $\xi_i$ 
  subroutine rhs_oneform( self, deg, func1, func2, func3, component, coefs_dofs )
    class(sll_t_maxwell_3d_trafo)       :: self          !< Maxwell solver class
    sll_int32,  intent( in    )         :: deg(3)        !< spline degree
    procedure(sll_i_function_3d_real64) :: func1         !< Function first component
    procedure(sll_i_function_3d_real64) :: func2         !< Function second component
    procedure(sll_i_function_3d_real64) :: func3         !< Function third component
    sll_int32,  intent( in    )         :: component     !< Specify the component of the function
    sll_real64, intent(   out )         :: coefs_dofs(:) !< Finite Element right-hand-side
    ! local variables
    sll_int32               :: i1, i2, i3,j1, j2, j3, k1, k2, k3, q(3), counter
    sll_real64              :: c(3)
    sll_real64              :: coef, jacobian, jmatrix(3,3)
    sll_real64, allocatable :: xw_gauss_d1(:,:), xw_gauss_d2(:,:), xw_gauss_d3(:,:)
    sll_real64, allocatable :: bspl_d1(:,:), bspl_d2(:,:), bspl_d3(:,:)
    sll_real64 :: scratch(maxval(deg+1))

    ! take enough Gauss points so that projection is exact for splines of degree deg
    q = deg+1
    ! rescale on [0,1] for compatibility with B-splines
    allocate(xw_gauss_d1(2, q(1)))
    allocate(bspl_d1(q(1), deg(1)+1 ))
    xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights(q(1), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k1 = 1, q(1)
       call sll_s_uniform_bsplines_eval_basis(deg(1),xw_gauss_d1(1,k1), scratch(1:deg(1)+1))
       bspl_d1(k1,:) = scratch(1:deg(1)+1)
    end do

    allocate(xw_gauss_d2(2, q(2)))
    allocate(bspl_d2(q(2), deg(2)+1))
    xw_gauss_d2 = sll_f_gauss_legendre_points_and_weights(q(2), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k2 = 1, q(2)
       call sll_s_uniform_bsplines_eval_basis(deg(2),xw_gauss_d2(1,k2), scratch(1:deg(2)+1))
       bspl_d2(k2,:) = scratch(1:deg(2)+1)
    end do

    allocate(xw_gauss_d3(2,q(3)))
    allocate(bspl_d3(q(3), deg(3)+1 ))
    xw_gauss_d3 = sll_f_gauss_legendre_points_and_weights(q(3), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k3 = 1, q(3)
       call sll_s_uniform_bsplines_eval_basis(deg(3),xw_gauss_d3(1,k3), scratch(1:deg(3)+1) )
       bspl_d3(k3,:) = scratch(1:deg(3)+1)
    end do


    counter = 1
    ! Compute coefs_dofs = int f(\xi)N_i(\xi) J_F(\xi) d\xi 
    do i3 = 1, self%n_dofs(3)
       do i2 = 1, self%n_dofs(2)
          do i1 = 1, self%n_dofs(1)
             coef = 0.0_f64
             ! loop over support of B spline
             do j3 = 1, deg(3)+1
                do j2 = 1, deg(2)+1
                   do j1 = 1, deg(1)+1
                      ! loop over Gauss points
                      do k3 = 1, q(3)
                         c(3) = self%delta_x(3) * (xw_gauss_d3(1,k3) + real(i3 + j3 - 2-((i3 + j3 - 2)/self%n_dofs(3)) * self%n_dofs(3),f64))
                         do k2 = 1, q(2)
                            c(2) = self%delta_x(2) * (xw_gauss_d2(1,k2) + real(i2 + j2 - 2-((i2 + j2 - 2)/self%n_dofs(2)) * self%n_dofs(2),f64))
                            do k1 = 1, q(1)
                               c(1) = self%delta_x(1) * (xw_gauss_d1(1,k1) + real(i1 + j1 - 2-((i1 + j1 - 2)/self%n_dofs(1)) * self%n_dofs(1),f64))
                               jacobian = self%map%jacobian( c ) 
                               jmatrix = self%map%jacobian_matrix_inverse_transposed( c )
                               coef = coef + xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* xw_gauss_d3(2,k3) *&
                                    ( func1( self%map%get_x(c) )*jmatrix(1,component) + &
                                    func2( self%map%get_x(c) )*jmatrix(2,component) + &
                                    func3( self%map%get_x(c) )*jmatrix(3,component)) * &
                                    (bspl_d1(k1,deg(1)+2-j1)*&
                                    bspl_d2(k2,deg(2)+2-j2)*&
                                    bspl_d3(k3,deg(3)+2-j3)) * abs(jacobian)
                            enddo
                         enddo
                      end do
                   end do
                end do
             end do
             ! rescale by cell size
             coefs_dofs(counter) = coef*self%volume
             counter = counter+1
          enddo
       end do
    end do

  end subroutine rhs_oneform


  !> Compute $\int f(F(\xi))\cdot DF(\xi) N_i(\xi) d\xi$ for vector function f, where $N_i$ is the B-spline starting at $\xi_i$ 
  subroutine rhs_twoform( self, deg, func1, func2, func3, component, coefs_dofs )
    class(sll_t_maxwell_3d_trafo)       :: self          !< Maxwell solver class
    sll_int32,  intent( in    )         :: deg(3)        !< spline degree
    procedure(sll_i_function_3d_real64) :: func1         !< Function first component
    procedure(sll_i_function_3d_real64) :: func2         !< Function second component
    procedure(sll_i_function_3d_real64) :: func3         !< Function third component
    sll_int32,  intent( in    )         :: component     !< Specify the component of the function
    sll_real64, intent(   out )         :: coefs_dofs(:) !< Finite Element right-hand-side
    ! local variables
    sll_int32 :: i1, i2, i3,j1, j2, j3, k1, k2, k3, q(3), counter
    sll_real64 :: c(3)
    sll_real64 :: coef, jmatrix(3,3)
    sll_real64, allocatable :: xw_gauss_d1(:,:), xw_gauss_d2(:,:), xw_gauss_d3(:,:)
    sll_real64, allocatable :: bspl_d1(:,:), bspl_d2(:,:), bspl_d3(:,:)
    sll_real64 :: scratch(maxval(deg+1))

    ! take enough Gauss points so that projection is exact for splines of degree deg
    q = deg+1
    ! rescale on [0,1] for compatibility with B-splines
    allocate(xw_gauss_d1(2, q(1)))
    allocate(bspl_d1(q(1), deg(1)+1 ))
    xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights(q(1), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k1 = 1, q(1)
       call sll_s_uniform_bsplines_eval_basis(deg(1),xw_gauss_d1(1,k1), scratch(1:deg(1)+1))
       bspl_d1(k1,:) = scratch(1:deg(1)+1)
    end do

    allocate(xw_gauss_d2(2, q(2)))
    allocate(bspl_d2(q(2), deg(2)+1))
    xw_gauss_d2 = sll_f_gauss_legendre_points_and_weights(q(2), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k2 = 1, q(2)
       call sll_s_uniform_bsplines_eval_basis(deg(2),xw_gauss_d2(1,k2), scratch(1:deg(2)+1))
       bspl_d2(k2,:) = scratch(1:deg(2)+1)
    end do

    allocate(xw_gauss_d3(2,q(3)))
    allocate(bspl_d3(q(3), deg(3)+1 ))
    xw_gauss_d3 = sll_f_gauss_legendre_points_and_weights(q(3), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k3 = 1, q(3)
       call sll_s_uniform_bsplines_eval_basis(deg(3),xw_gauss_d3(1,k3), scratch(1:deg(3)+1) )
       bspl_d3(k3,:) = scratch(1:deg(3)+1)
    end do


    counter = 1
    ! Compute coefs_dofs = int f(\xi)N_i(\xi) J_F(\xi) d\xi 
    do i3 = 1, self%n_dofs(3)
       do i2 = 1, self%n_dofs(2)
          do i1 = 1, self%n_dofs(1)
             coef = 0.0_f64
             ! loop over support of B spline
             do j3 = 1, deg(3)+1
                do j2 = 1, deg(2)+1
                   do j1 = 1, deg(1)+1
                      ! loop over Gauss points
                      do k3 = 1, q(3)
                         c(3) = self%delta_x(3) * (xw_gauss_d3(1,k3) + real(i3 + j3 - 2-((i3 + j3 - 2)/self%n_dofs(3)) * self%n_dofs(3),f64))
                         do k2 = 1, q(2)
                            c(2) = self%delta_x(2) * (xw_gauss_d2(1,k2) + real(i2 + j2 - 2-((i2 + j2 - 2)/self%n_dofs(2)) * self%n_dofs(2),f64))
                            do k1 = 1, q(1)
                               c(1) = self%delta_x(1) * (xw_gauss_d1(1,k1) + real(i1 + j1 - 2-((i1 + j1 - 2)/self%n_dofs(1)) * self%n_dofs(1),f64))

                               jmatrix = self%map%jacobian_matrix( c ) * sign( 1._f64, self%map%jacobian( c ) )
                               coef = coef + xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* xw_gauss_d3(2,k3) *&
                                    (func1( self%map%get_x(c) ) * jmatrix(1,component) + &
                                    func2( self%map%get_x(c) ) * jmatrix(2,component) + &
                                    func3( self%map%get_x(c) ) * jmatrix(3,component)) *&
                                    (bspl_d1(k1,deg(1)+2-j1)*&
                                    bspl_d2(k2,deg(2)+2-j2)*&
                                    bspl_d3(k3,deg(3)+2-j3))
                            enddo
                         enddo
                      end do
                   end do
                end do
             end do
             ! rescale by cell size
             coefs_dofs(counter) = coef*self%volume
             counter = counter+1
          enddo
       end do
    end do

  end subroutine rhs_twoform


  !> Compute $\int f(F(\xi)) N_i(\xi) d\xi$ for scalar function f, where $N_i$ is the B-spline starting at $\xi_i$ 
  subroutine rhs_threeform( self, deg, func, coefs_dofs )
    class(sll_t_maxwell_3d_trafo)       :: self          !< Maxwell solver class
    sll_int32,  intent( in    )         :: deg(3)           !< spline degree  
    procedure(sll_i_function_3d_real64) :: func          !< Function
    sll_real64, intent(   out )         :: coefs_dofs(:) !< Finite Element right-hand-side

    ! local variables
    sll_int32               :: i1, i2, i3,j1, j2, j3, k1, k2, k3, q(3), counter
    sll_real64              :: coef, c(3)
    sll_real64, allocatable :: xw_gauss_d1(:,:), xw_gauss_d2(:,:), xw_gauss_d3(:,:)
    sll_real64, allocatable :: bspl_d1(:,:), bspl_d2(:,:), bspl_d3(:,:)

    ! take enough Gauss points so that projection is exact for splines of degree deg
    q = deg+1
    ! rescale on [0,1] for compatibility with B-splines
    allocate(xw_gauss_d1(2, q(1)))
    allocate(bspl_d1(q(1), deg(1)+1 ))
    xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights(q(1), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k1 = 1, q(1)
       call sll_s_uniform_bsplines_eval_basis(deg(1),xw_gauss_d1(1,k1), bspl_d1(k1,:))
    end do

    allocate(xw_gauss_d2(2, q(2)))
    allocate(bspl_d2(q(2), deg(2)+1))
    xw_gauss_d2 = sll_f_gauss_legendre_points_and_weights(q(2), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k2 = 1, q(2)
       call sll_s_uniform_bsplines_eval_basis(deg(2),xw_gauss_d2(1,k2), bspl_d2(k2,:))
    end do

    allocate(xw_gauss_d3(2,q(3)))
    allocate(bspl_d3(q(3), deg(3)+1 ))
    xw_gauss_d3 = sll_f_gauss_legendre_points_and_weights(q(3), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k3 = 1, q(3)
       call sll_s_uniform_bsplines_eval_basis(deg(3),xw_gauss_d3(1,k3), bspl_d3(k3,:))
    end do


    counter = 1
    ! Compute coefs_dofs = int f(\xi)N_i(\xi) J_F(\xi) d\xi 
    do i3 = 1, self%n_dofs(3)
       do i2 = 1, self%n_dofs(2)
          do i1 = 1, self%n_dofs(1)
             coef = 0.0_f64
             ! loop over support of B spline
             do j3 = 1, deg(3)+1
                do j2 = 1, deg(2)+1
                   do j1 = 1, deg(1)+1
                      ! loop over Gauss points
                      do k3 = 1, q(3)
                         c(3) = self%delta_x(3) * (xw_gauss_d3(1,k3) + real(i3 + j3 - 2-((i3 + j3 - 2)/self%n_dofs(3)) * self%n_dofs(3),f64))
                         do k2 = 1, q(2)
                            c(2) = self%delta_x(2) * (xw_gauss_d2(1,k2) + real(i2 + j2 - 2-((i2 + j2 - 2)/self%n_dofs(2)) * self%n_dofs(2),f64))
                            do k1 = 1, q(1)
                               c(1) = self%delta_x(1) * (xw_gauss_d1(1,k1) + real(i1 + j1 - 2-((i1 + j1 - 2)/self%n_dofs(1)) * self%n_dofs(1),f64))
                               coef = coef + xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* xw_gauss_d3(2,k3) *&
                                    func(self%map%get_x(c) ) * &
                                    bspl_d1(k1,deg(1)+2-j1)*&
                                    bspl_d2(k2,deg(2)+2-j2)*&
                                    bspl_d3(k3,deg(3)+2-j3)
                            enddo
                         enddo
                      end do
                   end do
                end do
             end do
             ! rescale by cell size
             coefs_dofs(counter) = coef*self%volume
             counter = counter+1
          enddo
       end do
    end do

  end subroutine rhs_threeform


  !> Compute the L2 projection of a given function f on periodic splines of given degree
  subroutine L2projection_3d_trafo( self, form,  component, coefs_dofs, func1, func2, func3 )
    class(sll_t_maxwell_3d_trafo)                  :: self          !< Maxwell solver class
    sll_int32,  intent( in    )                    :: form          !< Specify if the function is a0,1,2 or 3-form
    sll_int32,  intent( in    )                    :: component     !< Specify the component of the function
    sll_real64, intent(   out )                    :: coefs_dofs(:) !< Finite Element right-hand-side
    procedure(sll_i_function_3d_real64)            :: func1         !< Function first component
    procedure(sll_i_function_3d_real64), optional  :: func2         !< Function second component
    procedure(sll_i_function_3d_real64), optional  :: func3         !< Function third component

    if( present(func2) .and. present(func3) ) then
       ! Compute right-hand-side
       call sll_s_compute_rhs_trafo( self, form, 1, self%work(1:self%n_total), func1, func2, func3 )
       call sll_s_compute_rhs_trafo( self, form, 2, self%work(1+self%n_total:2*self%n_total), func1, func2, func3 )
       call sll_s_compute_rhs_trafo( self, form, 3, self%work(1+2*self%n_total:3*self%n_total), func1, func2, func3 )

       select case( form )
       case( 1 )
          call self%mass1_solver%solve( self%work, coefs_dofs )
       case( 2 )
          call self%mass2_solver%solve( self%work, coefs_dofs )
       case  default
          print*, 'L2projection for', form, '-form not implemented.'
       end select
    else
       select case( form )
       case(0)
          call sll_s_compute_rhs_trafo( self, form, 0, self%work0, func1 )
          call self%mass0_solver%solve(self%work0, coefs_dofs)
       case default
          print*,'l2 projection not for single function or component implemented'
       end select
    end if


  end subroutine L2projection_3d_trafo


  !> Compute square of the L2norm 
  function L2norm_squared_3d_trafo( self, coefs, form, component ) result ( r )
    class(sll_t_maxwell_3d_trafo) :: self      !< Maxwell solver class
    sll_real64                    :: coefs(:)  !< Coefficient for each DoF
    sll_int32                     :: form      !< Specify 0,1,2 or 3-form
    sll_int32                     :: component !< Specify the component of the form
    sll_real64                    :: r         !< Result: squared L2 norm

    r = inner_product_3d_trafo( self, coefs, coefs, form, component  )

  end function L2norm_squared_3d_trafo


  !> Compute inner product
  function inner_product_3d_trafo( self, coefs1, coefs2, form, component ) result ( r )
    class(sll_t_maxwell_3d_trafo) :: self      !< Maxwell solver class
    sll_real64                    :: coefs1(:) !< Coefficient for each DoF
    sll_real64                    :: coefs2(:) !< Coefficient for each DoF
    sll_int32                     :: form      !< Specify 0,1,2 or 3-form
    sll_int32, optional           :: component !< Specify the component of the form
    sll_real64                    :: r         !< Result: squared L2 norm

    !local variables
    sll_int32  :: istart, iend

    if ( form == 0 ) then
       call self%multiply_mass( [form], coefs2, self%work0)
       r = sum( coefs1*self%work0)
    else
       if (present(component)) then
          istart = 1+(component-1)*self%n_total
          iend =  component*self%n_total
       else
          istart=1
          iend=3*self%n_total
       end if
       call multiply_mass_trafo( self, [form], coefs2, self%work(1:3*self%n_total) )
       r = sum(coefs1(istart:iend)*self%work(istart:iend))
    end if

  end function inner_product_3d_trafo


  !> Initialization
  subroutine init_3d_trafo( self, domain, n_dofs, s_deg_0, map, mass_tolerance, poisson_tolerance, solver_tolerance, adiabatic_electrons, profile )
    class(sll_t_maxwell_3d_trafo),   intent( inout ) :: self      !< Maxwell solver class
    sll_real64, intent(in) :: domain(3,2) 
    sll_int32,                       intent( in    ) :: n_dofs(3) !< number of degrees of freedom (here number of cells and grid points)
    sll_int32,                       intent( in    ) :: s_deg_0(3)   !< highest spline degree
    type(sll_t_mapping_3d), target, intent( inout    ) :: map       !< coordinate transformation
    sll_real64, intent(in), optional :: mass_tolerance !< tolerance for mass solver
    sll_real64, intent(in), optional :: poisson_tolerance !< tolerance for Poisson solver
    sll_real64, intent(in), optional :: solver_tolerance !< tolerance for Schur complement solver
    logical, intent(in), optional :: adiabatic_electrons !< flag if adiabatic electrons are used
    type(sll_t_profile_functions), intent(in), optional :: profile !< temperature and density profiles
    ! local variables
    sll_int32  :: j
    sll_int32  :: deg1(3), deg2(3)
    sll_real64, allocatable :: nullspace(:,:)!, ones(:)

    self%Lx = map%Lx

    if (present( mass_tolerance) ) then
       self%mass_solver_tolerance = mass_tolerance
    else
       self%mass_solver_tolerance = 1d-12
    end if
    if (present( poisson_tolerance) ) then
       self%poisson_solver_tolerance = poisson_tolerance
    else
       self%poisson_solver_tolerance = 1d-12
    end if
    if (present( solver_tolerance) ) then
       self%solver_tolerance = solver_tolerance
    else
       self%solver_tolerance = 1d-12
    end if

    if( present( adiabatic_electrons ) ) then
       self%adiabatic_electrons = adiabatic_electrons
    end if

    if( present( profile ) ) then
       self%profile = profile
    end if

    self%map=>map
    self%n_cells = n_dofs
    self%n_dofs = n_dofs
    self%n_total = product(n_dofs)
    self%n_total0 = self%n_total
    self%n_total1 = self%n_total

    self%delta_x = 1._f64 / real(n_dofs,f64)
    self%s_deg_0 = s_deg_0
    self%s_deg_1 = s_deg_0 - 1
    self%volume = product(self%delta_x)

    ! Allocate scratch data
    allocate( self%work0(1:self%n_total) )
    allocate( self%work(1:self%n_total*3) )
    allocate( self%work2(1:self%n_total*3) )

    ! Assemble the sparse mass matrices
    !0-form
    deg1=self%s_deg_0
    if(self%adiabatic_electrons) then

       call sll_s_spline_fem_mass3d( n_dofs, deg1, -1, self%mass0, profile_m0 )
    else

       call sll_s_spline_fem_mass3d( n_dofs, deg1, 0, self%mass0, profile_0 )
    end if
    do j=1, 3
       deg1=self%s_deg_0
       deg1(j)=self%s_deg_1(j)
       call sll_s_spline_fem_mass3d( n_dofs, deg1, j, self%mass1(j,j), profile_1 )

       deg1 =  self%s_deg_1
       deg1(j) = self%s_deg_0(j)
       call sll_s_spline_fem_mass3d( n_dofs, deg1, j, self%mass2(j,j), profile_2 )

    end do
    if(self%map%flag2d)then
       deg1=self%s_deg_0
       deg1(1)=self%s_deg_1(1)
       deg2=self%s_deg_0
       deg2(2) = self%s_deg_1(2)
       call sll_s_spline_fem_mixedmass3d( n_dofs, deg1, deg2, [1,2], self%mass1(1,2), profile_m1 )
       call sll_s_spline_fem_mixedmass3d( n_dofs, deg2, deg1, [2,1], self%mass1(2,1), profile_m1 )

       deg1 =  self%s_deg_1
       deg1(1) = self%s_deg_0(1)
       deg2 = self%s_deg_1
       deg2(2) = self%s_deg_0(2)
       call sll_s_spline_fem_mixedmass3d( n_dofs, deg1, deg2, [1,2], self%mass2(1,2), profile_m2 )
       call sll_s_spline_fem_mixedmass3d( n_dofs, deg2, deg1, [2,1], self%mass2(2,1), profile_m2 )

       if(self%map%flag3d)then
          do j = 2, 3
             deg1=self%s_deg_0
             deg1(j)=self%s_deg_1(j)
             deg2=self%s_deg_0
             deg2(modulo(j,3)+1) = self%s_deg_1(modulo(j,3)+1)
             call sll_s_spline_fem_mixedmass3d( n_dofs, deg1, deg2, [j,modulo(j,3)+1], self%mass1(j,modulo(j,3)+1), profile_m1 )
             call sll_s_spline_fem_mixedmass3d( n_dofs, deg2, deg1, [modulo(j,3)+1,j], self%mass1(modulo(j,3)+1,j), profile_m1 )


             deg1 =  self%s_deg_1
             deg1(j) = self%s_deg_0(j)
             deg2 = self%s_deg_1
             deg2(modulo(j,3)+1) = self%s_deg_0(modulo(j,3)+1)
             call sll_s_spline_fem_mixedmass3d( n_dofs, deg1, deg2, [j,modulo(j,3)+1], self%mass2(j,modulo(j,3)+1), profile_m2 )
             call sll_s_spline_fem_mixedmass3d( n_dofs, deg2, deg1, [modulo(j,3)+1,j], self%mass2(modulo(j,3)+1,j), profile_m2 )
          end do
       end if
    end if

    call self%mass1_operator%create( 3, 3 )
    call self%mass2_operator%create( 3, 3 )

    do j=1,3
       call self%mass1_operator%set( j, j, self%mass1(j,j) )
       call self%mass2_operator%set( j, j, self%mass2(j,j) )
    end do
    if(self%map%flag2d)then
       call self%mass1_operator%set( 1, 2, self%mass1(1, 2) )
       call self%mass1_operator%set( 2, 1, self%mass1(2, 1) )
       call self%mass2_operator%set( 1, 2, self%mass2(1, 2) )
       call self%mass2_operator%set( 2, 1, self%mass2(2, 1) )
       if(self%map%flag3d)then
          do j = 2, 3
             call self%mass1_operator%set( j, modulo(j,3)+1, self%mass1(j, modulo(j,3)+1) )
             call self%mass1_operator%set( modulo(j,3)+1, j, self%mass1(modulo(j,3)+1, j) )
             call self%mass2_operator%set( j, modulo(j,3)+1, self%mass2(j, modulo(j,3)+1) )
             call self%mass2_operator%set( modulo(j,3)+1, j, self%mass2(modulo(j,3)+1, j) )
          end do
       end if
    end if

    call self%preconditioner_fft%init( self%Lx, n_dofs, s_deg_0 )
    call self%mass0_solver%create( self%mass0 )
    call self%mass1_solver%create( self%mass1_operator, self%preconditioner_fft%inverse_mass1_3d)
    call self%mass2_solver%create( self%mass2_operator, self%preconditioner_fft%inverse_mass2_3d)
    self%mass1_solver%atol = self%mass_solver_tolerance/maxval(self%Lx)
    self%mass2_solver%atol = self%mass_solver_tolerance/maxval(self%Lx)**2
    !self%mass1_solver%verbose = .true.
    !self%mass2_solver%verbose = .true.

    call self%poisson_matrix%create( self%mass1_operator, self%n_dofs, self%delta_x )
    ! Penalized Poisson operator
    allocate(nullspace(1,1:self%n_total))
    nullspace(1,:) = 1.0_f64
    call self%poisson_operator%create( self%poisson_matrix, nullspace, 1 )
    ! Poisson solver
    call self%poisson_solver%create( self%poisson_operator )
    self%poisson_solver%null_space = .true.
    self%poisson_solver%atol = self%poisson_solver_tolerance
    !self%poisson_solver%verbose = .true.
    self%poisson_solver%n_maxiter=40000

    ! Only for Schur complement eb solver
    call self%linear_op_schur_eb%create( self%mass1_operator, self%mass2_operator, self%n_total, self%n_dofs, self%delta_x   )
    call self%linear_solver_schur_eb%create( self%linear_op_schur_eb, self%preconditioner_fft%inverse_mass1_3d)
    self%linear_solver_schur_eb%atol = self%solver_tolerance/maxval(self%Lx)
    self%linear_solver_schur_eb%rtol = self%solver_tolerance/maxval(self%Lx)
    !self%linear_solver_schur_eb%verbose = .true.
    self%linear_solver_schur_eb%n_maxiter = 2000


  contains
    function profile_m0( x, component)
      sll_real64 :: profile_m0
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_m0 = self%map%jacobian( x ) * self%profile%rho_0( x(1) )/ self%profile%T_e( x(1) )

    end function profile_m0

    function profile_0(x, component)
      sll_real64 :: profile_0
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_0 = self%map%jacobian( x ) 

    end function profile_0

    function profile_1(x, component)
      sll_real64 :: profile_1
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_1 = self%map%metric_inverse_single_jacobian( x, component(1), component(1) )

    end function profile_1

    function profile_2(x, component)
      sll_real64 :: profile_2
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_2 = self%map%metric_single_jacobian( x, component(1), component(1) )

    end function profile_2

    function profile_m1(x, component)
      sll_real64 :: profile_m1
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_m1 = self%map%metric_inverse_single_jacobian( x, component(1), component(2) )
    end function profile_m1

    function profile_m2(x, component)
      sll_real64 :: profile_m2
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_m2 = self%map%metric_single_jacobian( x, component(1), component(2) )

    end function profile_m2

  end subroutine init_3d_trafo


  !> Initialization from nml file
  subroutine init_from_file_3d_trafo( self, domain, n_dofs, s_deg_0, map, nml_file, adiabatic_electrons, profile )
    class(sll_t_maxwell_3d_trafo),   intent( inout ) :: self     !< Maxwell solver class
    sll_real64, intent(in) :: domain(3,2) 
    sll_int32,                       intent( in    ) :: n_dofs(3) !< number of degrees of freedom (here number of cells and grid points)
    sll_int32,                       intent( in    ) :: s_deg_0(3)   !< highest spline degree
    type(sll_t_mapping_3d), target, intent( inout    ) :: map       !< coordinate transformation
    character(len=*), intent(in) :: nml_file !< nml-file
    logical, intent(in), optional :: adiabatic_electrons !< flag if adiabatic electrons are used
    type(sll_t_profile_functions), intent(in), optional :: profile !< temperature and density profiles
    ! local variables
    character(len=256) :: file_prefix
    sll_int32 :: input_file
    sll_int32 :: io_stat, io_stat0, io_stat1, rank, file_id
    sll_real64 :: mass_tolerance
    sll_real64 :: poisson_tolerance
    sll_real64 :: maxwell_tolerance

    namelist /output/ file_prefix
    namelist /maxwell_solver/ mass_tolerance, poisson_tolerance
    namelist /time_solver/ maxwell_tolerance 

    rank = sll_f_get_collective_rank(sll_v_world_collective)

    if( present( adiabatic_electrons ) ) then
       self%adiabatic_electrons = adiabatic_electrons
    end if

    if( present( profile ) ) then
       self%profile = profile
    end if

    ! Read in solver tolerance
    open(newunit = input_file, file=nml_file, status='old',IOStat=io_stat)
    if (io_stat /= 0) then
       if (rank == 0 ) then
          print*, 'sll_m_maxwell_3d_trafo: Input file does not exist. Set default tolerance.'
          open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
          write(file_id, *) 'mass solver tolerance:', 1d-12
          write(file_id, *) 'poisson solver tolerance:', 1d-12
          close(file_id)
       end if
       call self%init( domain, n_dofs, s_deg_0, map )
    else
       read(input_file, output,IOStat=io_stat0)
       read(input_file, maxwell_solver,IOStat=io_stat)
       read(input_file, time_solver,IOStat=io_stat1)
       if (io_stat /= 0 .and. io_stat1 /= 0) then
          if (rank == 0 ) then
             print*, 'sll_m_maxwell_3d_trafo: Input parameter does not exist. Set default tolerance.'
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'mass solver tolerance:', 1d-12
             write(file_id, *) 'poisson solver tolerance:', 1d-12
             close(file_id)
          end if
          call self%init( domain, n_dofs, s_deg_0, map )
       else if (io_stat /= 0 .and. io_stat1 == 0) then
          call self%init( domain, n_dofs, s_deg_0, map, solver_tolerance=maxwell_tolerance)
       else if (io_stat == 0 .and. io_stat1 /= 0) then
          if (rank == 0 ) then
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'mass solver tolerance:', mass_tolerance
             write(file_id, *) 'poisson solver tolerance:', poisson_tolerance
             close(file_id)
          end if
          call self%init( domain, n_dofs, s_deg_0, map, mass_tolerance, poisson_tolerance )
       else
          if (rank == 0 ) then
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'mass solver tolerance:', mass_tolerance
             write(file_id, *) 'poisson solver tolerance:', poisson_tolerance
             close(file_id)
          end if
          call self%init( domain, n_dofs, s_deg_0, map, mass_tolerance, poisson_tolerance, maxwell_tolerance )  
       end if
       close(input_file)
    end if

  end subroutine init_from_file_3d_trafo


  !> Finalization
  subroutine free_3d_trafo( self )
    class(sll_t_maxwell_3d_trafo) :: self !< Maxwell solver class
    self%map => null()
    call self%preconditioner_fft%free()
    call self%mass0_solver%free()
    call self%mass1_solver%free()
    call self%mass2_solver%free()
    call self%mass1_operator%free()
    call self%mass2_operator%free()
    call self%poisson_matrix%free()
    call self%poisson_solver%free()
    call self%poisson_operator%free()
    call self%linear_solver_schur_eb%free()
    call self%linear_op_schur_eb%free()
    deallocate(self%work0)
    deallocate(self%work)
    deallocate(self%work2)

  end subroutine free_3d_trafo


  !> Multiply by dicrete gradient matrix
  subroutine multiply_g( self, field_in, field_out )
    class(sll_t_maxwell_3d_trafo) :: self !< Maxwell solver class
    sll_real64, intent( in    )   :: field_in(:) !< field_in  
    sll_real64, intent(   out )   :: field_out(:) !< G*field_in

    call sll_s_multiply_g(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_g


  !> Multiply by transpose of dicrete gradient matrix
  subroutine multiply_gt(self, field_in, field_out)
    class(sll_t_maxwell_3d_trafo)  :: self !< Maxwell solver class
    sll_real64, intent( in    )  :: field_in(:) !< field_in  
    sll_real64, intent(   out )  :: field_out(:) !< G^T*field_in

    call sll_s_multiply_gt(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_gt


  !> Multiply by discrete curl matrix
  subroutine multiply_c(self, field_in, field_out)
    class(sll_t_maxwell_3d_trafo)  :: self !< Maxwell solver class
    sll_real64, intent( in    )  :: field_in(:) !< field_in   
    sll_real64, intent(   out )  :: field_out(:) !< C*field_in  

    call sll_s_multiply_c(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_c


  !> Multiply by transpose of discrete curl matrix
  subroutine multiply_ct(self, field_in, field_out)
    class(sll_t_maxwell_3d_trafo)  :: self !< Maxwell solver class
    sll_real64, intent( in    )  :: field_in(:) !< field_in
    sll_real64, intent(   out )  :: field_out(:) !< C^T*field_in

    call sll_s_multiply_ct(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_ct


  !> Multiply by the mass matrix 
  subroutine multiply_mass_trafo( self, deg, coefs_in, coefs_out )
    class(sll_t_maxwell_3d_trafo) :: self         !< Maxwell solver class
    sll_int32,  intent( in    )   :: deg(:)       !< \a deg/form specifies if we multiply the mass to a  1- or 2-form or a mix of both
    sll_real64, intent( in    )   :: coefs_in(:)  !< coefs_in
    sll_real64, intent(   out )   :: coefs_out(:) !< M*coefs_in


    SLL_ASSERT(deg(1)>=0 .and. deg(1)<=2)
    select case(deg(1))
    case(0)
       call self%mass0%dot( coefs_in, coefs_out ) 
    case(1)
       call self%mass1_operator%dot( coefs_in, coefs_out )
    case(2)
       call self%mass2_operator%dot( coefs_in, coefs_out )
    case default
       print*, 'multiply mass for other form not yet implemented'
       stop
    end select
  end subroutine multiply_mass_trafo

  
  !> Multiply by the inverse mass matrix 
  subroutine multiply_mass_inverse_trafo( self, form, coefs_in, coefs_out )
    class(sll_t_maxwell_3d_trafo) :: self         !< Maxwell solver class
    sll_int32,  intent( in    )   :: form       !< \a form specifies the form (Note: 0 for 0-form, 1 for 1-form, 2 for 2-form, 3 for 0-1-form mix)
    sll_real64, intent( in    )   :: coefs_in(:)  !< coefs_in
    sll_real64, intent(   out )   :: coefs_out(:) !< M^{-1}*coefs_in


    SLL_ASSERT(form>=0 .and. form<=2)
    select case(form)
    case(0)
       call self%mass0_solver%solve( coefs_in, coefs_out )
    case(1)
       call self%mass1_solver%solve( coefs_in, coefs_out )
    case(2)
       call self%mass2_solver%solve( coefs_in, coefs_out )
    case default
       print*, 'multiply inverse mass for other form not yet implemented'
       stop
    end select
  end subroutine multiply_mass_inverse_trafo


end module sll_m_maxwell_3d_trafo
