!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations in 3D
!> The linear systems are solved using iterative linear solvers
!> @details
!> 
!> @author
!> Katharina Kormann

module sll_m_maxwell_3d_fem
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

  use sll_m_linear_operator_kron, only : &
       sll_t_linear_operator_kron

  use sll_m_linear_operator_penalized, only : &
       sll_t_linear_operator_penalized

  use sll_m_linear_operator_poisson_3d, only : &
       sll_t_linear_operator_poisson_3d

  use sll_m_linear_operator_schur_eb_3d, only : &
       sll_t_linear_operator_schur_eb_3d

  use sll_m_linear_solver_cg, only : &
       sll_t_linear_solver_cg

  use sll_m_linear_solver_kron, only : &
       sll_t_linear_solver_kron

  use sll_m_linear_solver_mgmres, only : &
       sll_t_linear_solver_mgmres

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_maxwell_3d_base, only: &
       sll_i_function_3d_real64, &
       sll_c_maxwell_3d_base

  use sll_m_preconditioner_curl_solver_fft, only : &
       sll_t_preconditioner_curl_solver_fft

  use sll_m_preconditioner_fft, only : &
       sll_t_preconditioner_fft

  use sll_m_profile_functions, only: &
       sll_t_profile_functions

  use sll_m_spline_fem_utilities, only : &
       sll_s_spline_fem_mass_line, &
       sll_s_spline_fem_mixedmass_line, &
       sll_s_multiply_g, &
       sll_s_multiply_gt, &
       sll_s_multiply_c, &
       sll_s_multiply_ct

  use sll_m_spline_fem_utilities_3d, only: &
       sll_s_spline_fem_mass3d, &
       sll_s_spline_fem_mixedmass3d

  use sll_m_spline_fem_utilities_sparse, only : &
       sll_s_spline_fem_mass1d, &
       sll_s_spline_fem_mixedmass1d

  use sll_m_linear_operator_curl_3d

  use sll_m_linear_operator_GTM

  use sll_m_linear_operator_MG

  use sll_m_uzawa_iterator


  implicit none

  public :: &
       sll_t_maxwell_3d_fem

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_c_maxwell_3d_base) :: sll_t_maxwell_3d_fem

     sll_real64, allocatable :: work0(:)  !< scratch data
     sll_real64, allocatable :: work(:)  !< scratch data
     sll_real64, allocatable :: work3d(:,:,:)  !< scratch data
     sll_real64, allocatable :: work2(:)  !< scratch data
     sll_real64, allocatable :: work_d1(:) !< scratch data
     sll_real64, allocatable :: work_d2_in(:) !< scratch data
     sll_real64, allocatable :: work_d2_out(:) !< scratch data
     sll_real64, allocatable :: work_d3_in(:) !< scratch data
     sll_real64, allocatable :: work_d3_out(:) !< scratch data

     type(sll_t_matrix_csr)  :: mass0 !< 0-form mass matrix
     type(sll_t_matrix_csr)  :: mass1d(3,3) !< 1D mass matrices
     type(sll_t_linear_solver_cg) :: mass1d_solver(2,3) !< 1D mass matrix solvers
     type(sll_t_linear_solver_kron) :: mass_1_solver(3) !< Tensorproduct solver for 3D 1-form mass matrices
     type(sll_t_linear_solver_kron) :: mass_2_solver(3) !< Tensorproduct solver for 3D 2-form mass matrices

     type(sll_t_linear_operator_kron)  :: mass1(3) !< Tensorproduct 1-form mass matrix
     type(sll_t_linear_operator_kron)  :: mass2(3) !< Tensorproduct 2-form mass matrix
     type(sll_t_linear_operator_block) :: mass1_operator   !< block mass matrix
     type(sll_t_linear_operator_block) :: mass2_operator   !< block mass matrix
     type(sll_t_linear_solver_cg) :: mass0_solver     !< mass matrix solver
     type(sll_t_linear_solver_cg) :: mass1_solver     !< mass matrix solver
     type(sll_t_preconditioner_fft) :: preconditioner_fft !< preconditioner for mass matrices
     type(sll_t_linear_operator_poisson_3d) :: poisson_matrix  !< Poisson matrix
     type(sll_t_linear_operator_penalized)  :: poisson_operator !< Poisson matrix with constraint on constant vector
     type(sll_t_linear_solver_cg)  :: poisson_solver !< CG solver to invert Poisson matrix
     type( sll_t_linear_operator_schur_eb_3d ) :: linear_op_schur_eb !< Schur complement operator for advect_eb
     type( sll_t_linear_solver_mgmres )        :: linear_solver_schur_eb !< Schur complement solver for advect_eb

     type(sll_t_linear_operator_curl_3d) :: curl_matrix  !< curl matrix
     type(sll_t_linear_operator_penalized)  :: curl_operator !< curl matrix with constraint on constant vector
     type( sll_t_preconditioner_curl_solver_fft ) :: preconditioner_curl_fft
     type(sll_t_linear_solver_cg)  :: curl_solver !< CG solver to invert curl matrix
     !type(sll_t_linear_solver_mgmres)  :: curl_solver !< CG solver to invert curl matrix
     
     type(sll_t_linear_operator_MG) :: MG_operator
     
     type(sll_t_linear_operator_GTM) :: GTM_operator
     
     type(sll_t_uzawa_iterator) :: uzawa_iterator

     logical :: adiabatic_electrons = .false. !< flag if adiabatic electrons are used
     
   contains
     procedure :: &
          compute_E_from_B => sll_s_compute_e_from_b_3d_fem !< Solve E and B part of Amperes law with B constant in time
     procedure :: &
          compute_B_from_E => sll_s_compute_b_from_e_3d_fem !< Solve Faraday equation with E constant in time
     procedure :: &
          compute_curl_part => sll_s_compute_curl_part_3d_fem !< Solve source-free Maxwell's equations
     procedure :: &
          compute_E_from_rho => sll_s_compute_e_from_rho_3d_fem !< Solve E from rho using Poisson
     procedure :: &
          compute_rho_from_e => sll_s_compute_rho_from_e_3d_fem !< Compute rho from E (by Poisson matrix multiply)
     procedure :: &
          compute_E_from_j => sll_s_compute_E_from_j_3d_fem !< Compute E from the current j
     procedure :: &
          compute_phi_from_rho => sll_s_compute_phi_from_rho_3d_fem !< Compute phi from rho (by solving the quasi-neutrality equation)
     procedure :: &
          compute_phi_from_j => sll_s_compute_phi_from_j_3d_fem !< Compute phi from j (dynamic of quasi-neutrality equation for adiabatic electrons)
     procedure :: &
          compute_rhs_from_function => sll_s_compute_rhs_fem !< Compute integral over given function tested by the basis 
     procedure :: &
          L2projection => L2projection_3d_fem !< Compute L_2 projection of a given function
     procedure :: &
          L2norm_squared => L2norm_squared_3d_fem !< Compute the square of the L2 norm of a given vector
     procedure :: &
          inner_product => inner_product_3d_fem !< Inner product of two dof-vectors with mass matrix
     procedure :: &
          init => init_3d_fem !< Initialize the Maxwell class
     procedure :: &
          init_from_file => init_from_file_3d_fem !< Initialize the Maxwell class with parameters read from nml-file
     procedure :: &
          free => free_3d_fem !< Free Maxwell class
     procedure :: &
          multiply_g  !< Multiplication with gradient matrix 
     procedure :: &
          multiply_gt !< Multiplication with transposed gradient matrix 
     procedure :: &
          multiply_c  !< Multiplication with curl matrix
     procedure :: &
          multiply_ct !< Multiplication with transposed curl matrix
     procedure :: &
          multiply_mass => multiply_mass_3d_fem  !< Product with the mass matrix
     procedure :: &
          multiply_mass_inverse => multiply_mass_inverse_3dkron !< Invert mass matrix
     procedure :: &
          compute_field_energy !< Compute field energy

  end type sll_t_maxwell_3d_fem

contains


  !> compute Ey from Bz using weak Ampere formulation 
  subroutine sll_s_compute_e_from_b_3d_fem(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_3d_fem) :: self        !< Maxwell solver class
    sll_real64, intent(in)     :: delta_t      !< Time step
    sll_real64, intent(in)     :: field_in(:)  !< B
    sll_real64, intent(inout)  :: field_out(:) !< E

    call self%mass2_operator%dot( field_in, self%work )
    call multiply_ct(self, self%work, self%work2)

    call self%mass1_solver%solve( self%work2, self%work )
    ! Update b from self value
    field_out = field_out + delta_t*self%work

  end subroutine sll_s_compute_e_from_b_3d_fem


  !> Compute Bz from Ey using strong 1D Faraday equation for spline coefficients
  !> $B_z^{new}(x_j) = B_z^{old}(x_j) - \frac{\Delta t}{\Delta x} (E_y(x_j) - E_y(x_{j-1})  $
  subroutine sll_s_compute_b_from_e_3d_fem(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_3d_fem)  :: self       !< Maxwell solver class
    sll_real64, intent(in)     :: delta_t      !< Time step
    sll_real64, intent(in)     :: field_in(:)  ! E
    sll_real64, intent(inout)  :: field_out(:) ! B 

    call multiply_c(self, field_in, self%work)
    ! Update b from self value
    field_out = field_out - delta_t * self%work

  end subroutine sll_s_compute_b_from_e_3d_fem


  !> Solve curl part of Maxwell's equations
  subroutine sll_s_compute_curl_part_3d_fem( self, delta_t, efield, bfield, betar )
    class(sll_t_maxwell_3d_fem) :: self         !< Maxwell solver class
    sll_real64, intent( in    )   :: delta_t    !< Time step
    sll_real64, intent( inout )   :: efield(:)  !< E
    sll_real64, intent( inout )   :: bfield(:)  !< B
    sll_real64, optional          :: betar      !< 1/beta
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

  end subroutine sll_s_compute_curl_part_3d_fem


  !> Compute E_i from rho_i integrated over the time interval using weak Poisson's equation
  subroutine sll_s_compute_E_from_rho_3d_fem( self, field_in, field_out )  
    class(sll_t_maxwell_3d_fem) :: self         !< Maxwell solver class
    sll_real64, intent( in    ) :: field_in(:)  !< rho
    sll_real64, intent(   out ) :: field_out(:) !< E

    ! Version with iterative solver
    call self%poisson_solver%solve( field_in, self%work(1:self%n_total) )
    call multiply_g(self, self%work(1:self%n_total), field_out)
    field_out = -field_out

  end subroutine sll_s_compute_e_from_rho_3d_fem


  !> compute rho from e using weak Gauss law ( rho = G^T M_1 e ) 
  subroutine sll_s_compute_rho_from_E_3d_fem( self, field_in, field_out ) 
    class(sll_t_maxwell_3d_fem) :: self         !< Maxwell solver class
    sll_real64, intent( in    ) :: field_in(:)  !< E
    sll_real64, intent(   out ) :: field_out(:) !< rho

    call self%mass1_operator%dot( field_in, self%work )

    call multiply_gt( self, self%work, field_out )
    field_out = - field_out

  end subroutine sll_s_compute_rho_from_e_3d_fem


  !> Compute E_i from j_i integrated over the time interval using weak Ampere formulation
  subroutine sll_s_compute_E_from_j_3d_fem( self, current, E, component )
    class(sll_t_maxwell_3d_fem)           :: self       !< Maxwell solver class
    sll_real64, intent( in    )           :: current(:) !< Current integrated over time interval
    sll_real64, intent( inout )           :: E(:)       !< Updated electric field
    sll_int32,  intent( in    ), optional :: component  !< component of the Efield to be computed

    if(present(component)) then
       call self%mass_1_solver(component)%solve( current, self%work(1:self%n_total) )
       E = E - self%work(1:self%n_total)
    else
       call self%mass1_solver%solve( current, self%work )
       E = E - self%work
    end if

  end subroutine sll_s_compute_E_from_j_3d_fem


  !> Compute phi from rho_i integrated over the time interval
  subroutine sll_s_compute_phi_from_rho_3d_fem( self, field_in, field_out, efield_dofs )  
    class(sll_t_maxwell_3d_fem)           :: self           !< Maxwell solver class
    sll_real64, intent( in    )           :: field_in(:)    !< rho
    sll_real64, intent( inout )           :: field_out(:)   !< phi
    sll_real64, intent(   out )           :: efield_dofs(:) !< E

    call self%mass0_solver%solve( field_in, field_out )
    call self%multiply_g( field_out, efield_dofs )
    efield_dofs = -efield_dofs

  end subroutine sll_s_compute_phi_from_rho_3d_fem


  !> Compute phi from j_i integrated over the time interval, delta_t is already included 
  subroutine sll_s_compute_phi_from_j_3d_fem( self, field_in, field_out, efield_dofs )
    class(sll_t_maxwell_3d_fem)           :: self           !< Maxwell solver class
    sll_real64, intent( in    )           :: field_in(:)    !< Current integrated over time interval
    sll_real64, intent( inout )           :: field_out(:)   !< phi
    sll_real64, intent(   out )           :: efield_dofs(:) !< E

    self%work0 = 0._f64
    self%work = 0._f64
    call self%multiply_gt( field_in, self%work0 ) 
    call self%mass0_solver%solve( self%work0, self%work(1:self%n_total) )
    field_out = field_out + self%work(1:self%n_total)

    call self%multiply_g( field_out, efield_dofs )
    efield_dofs = -efield_dofs

  end subroutine sll_s_compute_phi_from_j_3d_fem


  !> Compute the FEM right-hand-side for a given function f and periodic splines of given degree
  !> Its components are $\int f N_i dx$ where $N_i$ is the B-spline starting at $x_i$ 
  subroutine sll_s_compute_rhs_fem(self, form, component, coefs_dofs, func1, func2, func3)
    class(sll_t_maxwell_3d_fem)                    :: self          !< Maxwell solver class
    sll_int32,  intent( in    )                    :: form          !< Specify if the function is a0,1,2 or 3-form
    sll_int32,  intent( in    )                    :: component     !< Specify the component of the function
    sll_real64, intent(   out )                    :: coefs_dofs(:) !< Finite Element right-hand-side
    procedure(sll_i_function_3d_real64)            :: func1         !< Function first component
    procedure(sll_i_function_3d_real64), optional  :: func2         !< Function second component
    procedure(sll_i_function_3d_real64), optional  :: func3         !< Function third component
    ! local variables
    sll_int32 :: i1, i2, i3,j1, j2, j3, k1, k2, k3, q(3), counter
    sll_int32 :: degree(3)
    sll_real64 :: c(3)
    sll_real64 :: coef
    sll_real64, allocatable :: xw_gauss_d1(:,:), xw_gauss_d2(:,:), xw_gauss_d3(:,:)
    sll_real64, allocatable :: bspl_d1(:,:), bspl_d2(:,:), bspl_d3(:,:)
    sll_real64 :: scratch(maxval(self%s_deg_0+1))

    ! Define the spline degree in the 3 dimensions, depending on form and component of the form
    if ( form == 0 ) then
       degree = self%s_deg_0
    elseif (form == 1 ) then
       degree = self%s_deg_0
       degree(component) = self%s_deg_1(component)
    elseif( form == 2) then
       degree =  self%s_deg_1
       degree(component) = self%s_deg_0(component)
    elseif( form == 3) then
       degree =  self%s_deg_1
    else 
       print*, 'Wrong form.'
    end if

    ! take enough Gauss points so that projection is exact for splines of degree deg
    q = degree+1
    ! rescale on [0,1] for compatibility with B-splines
    allocate(xw_gauss_d1(2,q(1)))
    allocate(bspl_d1(q(1), degree(1)+1))
    xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights(q(1), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k1 = 1, q(1)
       call sll_s_uniform_bsplines_eval_basis(degree(1),xw_gauss_d1(1,k1), scratch(1:degree(1)+1))
       bspl_d1(k1,:) = scratch(1:degree(1)+1)
    end do

    allocate(xw_gauss_d2(2,q(2)))
    allocate(bspl_d2(q(2), degree(2)+1))
    xw_gauss_d2 = sll_f_gauss_legendre_points_and_weights(q(2), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k2 = 1, q(2)
       call sll_s_uniform_bsplines_eval_basis(degree(2),xw_gauss_d2(1,k2), scratch(1:degree(2)+1))
       bspl_d2(k2,:) = scratch(1:degree(2)+1)
    end do

    allocate(xw_gauss_d3(2,q(3)))
    allocate(bspl_d3(q(3), degree(3)+1))
    xw_gauss_d3 = sll_f_gauss_legendre_points_and_weights(q(3), 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k3 = 1, q(3)
       call sll_s_uniform_bsplines_eval_basis(degree(3),xw_gauss_d3(1,k3), scratch(1:degree(3)+1))
       bspl_d3(k3,:) = scratch(1:degree(3)+1)
    end do

    counter = 1
    ! Compute coefs_dofs = int f(x)N_i(x) 
    do i3 = 1, self%n_dofs(3)
       do i2 = 1, self%n_dofs(2)
          do i1 = 1, self%n_dofs(1)
             coef=0.0_f64
             ! loop over support of B spline
             do j3 = 1, degree(3)+1
                do j2 = 1, degree(2)+1
                   do j1 = 1, degree(1)+1
                      ! loop over Gauss points
                      do k3 = 1, q(3)
                         c(3) = self%delta_x(3)*(xw_gauss_d3(1,k3) + real(i3 + j3 - 2,f64))
                         do k2 = 1, q(2)
                            c(2) = self%delta_x(2)*(xw_gauss_d2(1,k2) + real(i2 + j2 - 2,f64))
                            do k1 = 1, q(1)
                               c(1) = self%delta_x(1)*(xw_gauss_d1(1,k1) + real(i1 + j1 - 2,f64))
                               coef = coef + xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* &
                                    xw_gauss_d3(2,k3) *&
                                    func1( c ) * &
                                    bspl_d1(k1,degree(1)+2-j1)*&
                                    bspl_d2(k2,degree(2)+2-j2)*&
                                    bspl_d3(k3,degree(3)+2-j3)

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

  end subroutine sll_s_compute_rhs_fem


  !> Compute the L2 projection of a given function f on periodic splines of given degree
  subroutine L2projection_3d_fem(self, form,  component, coefs_dofs, func1, func2, func3 )
    class(sll_t_maxwell_3d_fem)                    :: self          !< Maxwell solver class
    sll_int32,  intent( in    )                    :: form          !< Specify if the function is a0,1,2 or 3-form
    sll_int32,  intent( in    )                    :: component     !< Specify the component of the function
    sll_real64, intent(   out )                    :: coefs_dofs(:) !< Finite Element right-hand-side
    procedure(sll_i_function_3d_real64)            :: func1         !< Function first component
    procedure(sll_i_function_3d_real64), optional  :: func2         !< Function second component
    procedure(sll_i_function_3d_real64), optional  :: func3         !< Function third component

    ! Compute right-hand-side
    call sll_s_compute_rhs_fem( self, form,  component, self%work(1:self%n_total), func1 )

    select case( form )
    case( 0 )
       call self%mass0_solver%solve( self%work(1:self%n_total), coefs_dofs )
    case( 1 )
       call self%mass_1_solver(component)%solve( self%work(1:self%n_total), coefs_dofs )
    case( 2 )
       call self%mass_2_solver(component)%solve( self%work(1:self%n_total), coefs_dofs )
    case  default
       print*, 'L2projection for', form, '-form not implemented.'
    end select

  end subroutine L2projection_3d_fem

  subroutine compute_dofs(self, form,  component, coefs_dofs, func1, func2, func3 )
    class(sll_t_maxwell_3d_fem)                    :: self          !< Maxwell solver class
    sll_int32,  intent( in    )                    :: form          !< Specify if the function is a0,1,2 or 3-form
    sll_int32,  intent( in    )                    :: component     !< Specify the component of the function
    sll_real64, intent(   out )                    :: coefs_dofs(:) !< Finite Element right-hand-side
    procedure(sll_i_function_3d_real64)            :: func1         !< Function first component
    procedure(sll_i_function_3d_real64), optional  :: func2         !< Function second component
    procedure(sll_i_function_3d_real64), optional  :: func3         !< Function third component

    select case( form )
    case( 0 )
       !call sll_compute_point_integral( self, func1, coefs_dofs )
    case( 1 )
       call sll_compute_edge_integral( self, func1, func2, func3, coefs_dofs )
    case( 2 )
       !call sll_compute_face_integral( self, func1, func2, func3, coefs_dofs )
    case  default
       print*, 'projection for', form, '-form not implemented.'
    end select
    
  end subroutine compute_dofs
  
  subroutine sll_compute_edge_integral( self, func1, func2, func3, coefs_dofs )
    class(sll_t_maxwell_3d_fem)                    :: self          !< Maxwell solver class
    procedure(sll_i_function_3d_real64)            :: func1         !< Function first component
    procedure(sll_i_function_3d_real64)            :: func2         !< Function second component
    procedure(sll_i_function_3d_real64)            :: func3         !< Function third component
    sll_real64, intent(   out )                    :: coefs_dofs(:) !< Finite Element right-hand-side
    !local variables
    sll_int32 :: i1, i2, i3, j, q, ind
    sll_real64 :: c(3), e(3)
    sll_real64, allocatable :: xw_gauss(:,:)

    q = 4
    ! rescale on [0,1] for compatibility with B-splines
    allocate(xw_gauss(2,q))
    xw_gauss = sll_f_gauss_legendre_points_and_weights(q, 0.0_f64, 1.0_f64)

    ind = 0
    do i3 = 1, self%n_dofs(3)
       c(3) = self%delta_x(3)* real(i3-1,f64)
       do i2 = 1, self%n_dofs(2)
          c(2) = self%delta_x(2)* real(i2-1,f64)
          do i1 = 1, self%n_dofs(1)
             c(1) = self%delta_x(1)* real(i1-1,f64)
             ind = ind +1 !i + (j-1)*self%n_dofs(1) + (k-1)*self%n_dofs(1)*self%n_dofs(2)
             do j = 1, q
                e(1) = self%delta_x(1)* (xw_gauss(1,q) + real(i1-1,f64) )
                e(2) = self%delta_x(2)* (xw_gauss(1,q) + real(i2-1,f64) )
                e(3) = self%delta_x(3)* (xw_gauss(1,q) + real(i3-1,f64) )
                coefs_dofs(ind) = coefs_dofs(ind) + self%delta_x(1)*xw_gauss(2,j)*c(1)*func1([e(1),c(2),c(3)])
                coefs_dofs(ind+self%n_total) = coefs_dofs(ind+self%n_total) + self%delta_x(2)*xw_gauss(2,j)*c(2)*func2([c(1),e(2),c(3)])
                coefs_dofs(ind+2*self%n_total) = coefs_dofs(ind+2*self%n_total) + self%delta_x(3)*xw_gauss(2,j)*c(3)*func3([c(1),c(2),e(3)])
             end do
          end do
       end do
    end do

  end subroutine sll_compute_edge_integral


  !> Compute square of the L2norm 
  function L2norm_squared_3d_fem( self, coefs, form, component) result (r)
    class(sll_t_maxwell_3d_fem) :: self !< Maxwell solver class
    sll_real64 :: coefs(:) !< Coefficient for each DoF
    sll_int32  :: form !< Specify 0,1,2 or 3 form
    sll_int32  :: component !< Specify the component
    sll_real64 :: r !< Result: squared L2 norm


    r = inner_product_3d_fem(self, coefs, coefs, form, component)


  end function L2norm_squared_3d_fem


  !> Compute inner product
  function inner_product_3d_fem( self, coefs1, coefs2, form, component ) result ( r )
    class(sll_t_maxwell_3d_fem)   :: self      !< Maxwell solver class
    sll_real64                    :: coefs1(:) !< Coefficient for each DoF
    sll_real64                    :: coefs2(:) !< Coefficient for each DoF
    sll_int32                     :: form      !< Specify 0,1,2 or 3-form
    sll_int32, optional           :: component !< Specify the component of the form
    sll_real64                    :: r         !< Result: squared L2 norm
    !local variables
    sll_int32 :: deg(3)

    if ( form == 0 ) then
       !deg = 0
       call self%multiply_mass( [0], coefs2, self%work0 )
    elseif (form == 1 ) then
       deg = 1
       deg(component) = 2
       call multiply_mass_3dkron( self, deg, coefs2, self%work0 )
    elseif( form == 2) then
       deg = 2
       deg(component) = 1
       call multiply_mass_3dkron( self, deg, coefs2, self%work0 )
    elseif( form == 3) then
       deg = 2
       call multiply_mass_3dkron( self, deg, coefs2, self%work0 )
    else
       print*, 'Wrong form.'
    end if


    r = sum(coefs1*self%work0)

  end function inner_product_3d_fem


  !> Initialization
  subroutine init_3d_fem( self, domain, n_dofs, s_deg_0, mass_tolerance, poisson_tolerance, solver_tolerance, adiabatic_electrons, profile  )
    class(sll_t_maxwell_3d_fem), intent(inout) :: self !< Maxwell solver class
    sll_real64, intent(in) :: domain(3,2)     !< xmin, xmax
    sll_int32, intent(in) :: n_dofs(3)  !< number of degrees of freedom (here number of cells and grid points)
    sll_int32, intent(in) :: s_deg_0(3) !< highest spline degree
    sll_real64, intent(in), optional :: mass_tolerance !< tolerance for mass solver
    sll_real64, intent(in), optional :: poisson_tolerance !< tolerance for Poisson solver
    sll_real64, intent(in), optional :: solver_tolerance !< tolerance for Schur complement solver
    logical, intent(in), optional :: adiabatic_electrons !< flag if adiabatic electrons are used
    type(sll_t_profile_functions), intent(in), optional :: profile !< temperature and density profiles
    ! local variables
    sll_int32 :: j,k
    sll_real64, allocatable ::  mass_line_0(:), mass_line_1(:), mass_line_mixed(:)
    sll_real64, allocatable :: nullspace(:,:)

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

    self%n_cells = n_dofs
    self%n_dofs = n_dofs
    self%n_total = product(n_dofs)
    self%n_total0 = self%n_total
    self%n_total1 = self%n_total

    self%Lx = domain(:,2) - domain(:,1)
    self%delta_x = self%Lx / real(n_dofs, f64)
    self%s_deg_0 = s_deg_0
    self%s_deg_1 = s_deg_0 - 1

    self%volume = product(self%delta_x)

    ! Allocate scratch data
    allocate( self%work3d(n_dofs(1), n_dofs(2), n_dofs(3)) )
    allocate( self%work0(self%n_total) )
    allocate( self%work(self%n_total*3) )
    allocate( self%work2(self%n_total*3) )
    allocate( self%work_d1( n_dofs(1) ) ) 
    allocate( self%work_d2_in( n_dofs(2) ) ) 
    allocate( self%work_d2_out( n_dofs(2) ) ) 
    allocate( self%work_d3_in( n_dofs(3) ) ) 
    allocate( self%work_d3_out( n_dofs(3) ) )

    ! Sparse matrices
    ! Assemble the mass matrices
    ! First assemble a mass line for both degrees
    ! Next put together the 1d parts of the 3d Kronecker product
    do j=1, 3
       allocate( mass_line_0(s_deg_0(j)+1) ) 
       allocate( mass_line_1(s_deg_0(j)) )
       allocate( mass_line_mixed(s_deg_0(j)*2) )
       call sll_s_spline_fem_mass_line ( self%s_deg_0(j), mass_line_0 )
       call sll_s_spline_fem_mass_line ( self%s_deg_1(j), mass_line_1 )

       call sll_s_spline_fem_mixedmass_line ( self%s_deg_0(j), mass_line_mixed )


       call sll_s_spline_fem_mass1d( self%n_dofs(j), self%s_deg_0(j), mass_line_0, self%mass1d(1,j) )
       call sll_s_spline_fem_mass1d( self%n_dofs(j), self%s_deg_1(j), mass_line_1, self%mass1d(2,j) )
       self%mass1d(1,j)%arr_a = self%mass1d(1,j)%arr_a*self%delta_x(j)
       self%mass1d(2,j)%arr_a = self%mass1d(2,j)%arr_a*self%delta_x(j)
       call self%mass1d_solver(1,j)%create( self%mass1d(1,j) )
       call self%mass1d_solver(2,j)%create( self%mass1d(2,j) )

       call sll_s_spline_fem_mixedmass1d( self%n_dofs(j), self%s_deg_0(j), mass_line_mixed*self%delta_x(j), &
            self%mass1d(3,j) )

       deallocate( mass_line_0 )
       deallocate( mass_line_1 )
       deallocate( mass_line_mixed )
    end do

    do j=1,3
       do k=1,2
          self%mass1d_solver(k,j)%atol = self%mass_solver_tolerance
          !self%mass1d_solver(k,j)%verbose = .true.
       end do
    end do

    ! Put together the componentwise matrices as Kronecker products
    call self%mass_1_solver(1)%create( linear_solver_a=self%mass1d_solver(2,1), &
         linear_solver_b=self%mass1d_solver(1,2), &
         linear_solver_c=self%mass1d_solver(1,3) )
    call self%mass_1_solver(2)%create( linear_solver_a=self%mass1d_solver(1,1), &
         linear_solver_b=self%mass1d_solver(2,2), &
         linear_solver_c=self%mass1d_solver(1,3) )
    call self%mass_1_solver(3)%create( linear_solver_a=self%mass1d_solver(1,1), &
         linear_solver_b=self%mass1d_solver(1,2), &
         linear_solver_c=self%mass1d_solver(2,3))
    call self%mass_2_solver(1)%create( linear_solver_a=self%mass1d_solver(1,1), &
         linear_solver_b=self%mass1d_solver(2,2), &
         linear_solver_c=self%mass1d_solver(2,3) )
    call self%mass_2_solver(2)%create( linear_solver_a=self%mass1d_solver(2,1), &
         linear_solver_b=self%mass1d_solver(1,2), &
         linear_solver_c=self%mass1d_solver(2,3) )
    call self%mass_2_solver(3)%create( linear_solver_a=self%mass1d_solver(2,1), &
         linear_solver_b=self%mass1d_solver(2,2), &
         linear_solver_c=self%mass1d_solver(1,3) )


    if(self%adiabatic_electrons) then
       call sll_s_spline_fem_mass3d( self%n_dofs, s_deg_0, -1, self%mass0, profile_m0 )
    else
       call sll_s_spline_fem_mass3d( self%n_dofs, s_deg_0, 0, self%mass0, profile_0 )
    end if
    call self%mass0_solver%create( self%mass0 )
    self%mass0_solver%atol = self%mass_solver_tolerance
    !self%mass0_solver%verbose = .true.

    call self%mass1(1)%create( linop_a=self%mass1d(1,3), &
         linop_b=self%mass1d(1,2), &
         linop_c=self%mass1d(2,1) )
    call self%mass1(2)%create( linop_a=self%mass1d(1,3), &
         linop_b=self%mass1d(2,2), &
         linop_c=self%mass1d(1,1) )
    call self%mass1(3)%create( linop_a=self%mass1d(2,3), &
         linop_b=self%mass1d(1,2), &
         linop_c=self%mass1d(1,1))

    call self%mass2(1)%create( linop_a=self%mass1d(2,3), &
         linop_b=self%mass1d(2,2), &
         linop_c=self%mass1d(1,1) )
    call self%mass2(2)%create( linop_a=self%mass1d(2,3), &
         linop_b=self%mass1d(1,2), &
         linop_c=self%mass1d(2,1) )
    call self%mass2(3)%create( linop_a=self%mass1d(1,3), &
         linop_b=self%mass1d(2,2), &
         linop_c=self%mass1d(2,1))

    call self%mass1_operator%create( 3, 3 )
    call self%mass2_operator%create( 3, 3 )
    do j= 1, 3
       call self%mass1_operator%set( j, j, self%mass1(j) )
       call self%mass2_operator%set( j, j, self%mass2(j) )
    end do
    call self%preconditioner_fft%init( self%Lx, n_dofs, s_deg_0 )
    call self%mass1_solver%create( self%mass1_operator, self%preconditioner_fft%inverse_mass1_3d)
    self%mass1_solver%atol = self%mass_solver_tolerance
    !self%mass1_solver%verbose = .true.

    call self%poisson_matrix%create( self%mass1_operator, self%n_dofs, self%delta_x )
    ! Penalized Poisson operator
    allocate(nullspace(1,1:3*self%n_total))
    nullspace(1,:) = 1.0_f64
    call self%poisson_operator%create( linear_operator=self%poisson_matrix, vecs=nullspace(:,1:self%n_total), n_dim_nullspace=1 )
    ! Poisson solver
    call self%poisson_solver%create( self%poisson_operator )
    self%poisson_solver%null_space = .true.
    self%poisson_solver%atol = self%poisson_solver_tolerance
    !self%poisson_solver%verbose = .true.
    self%poisson_solver%n_maxiter=40000


    ! Only for Schur complement eb solver
    call self%linear_op_schur_eb%create( self%mass1_operator, self%mass2_operator, self%n_total, self%n_dofs, self%delta_x )
    call self%linear_solver_schur_eb%create( self%linear_op_schur_eb, self%preconditioner_fft%inverse_mass1_3d )
    self%linear_solver_schur_eb%atol = self%solver_tolerance
    self%linear_solver_schur_eb%rtol = self%solver_tolerance
    !self%linear_solver_schur_eb%verbose = .true.

    !call self%preconditioner_curl_fft%create( self%n_dofs, self%delta_x, self%s_deg_0 )

    call self%curl_matrix%create( self%mass1_operator, self%mass2_operator, self%n_dofs, self%delta_x   )
    call self%curl_operator%create( linear_operator=self%curl_matrix, vecs=nullspace, n_dim_nullspace=1 )
    call self%curl_solver%create( self%curl_operator )!, self%preconditioner_curl_fft )
    self%curl_solver%null_space = .true.
    self%curl_solver%atol = self%solver_tolerance
    self%curl_solver%verbose = .true.
    !self%curl_solver%n_maxiter=1000

    call self%MG_operator%create( self%mass1_operator, self%n_dofs, self%delta_x )
    call self%GTM_operator%create( self%mass1_operator, self%n_dofs, self%delta_x )
    call self%uzawa_iterator%create( self%curl_solver, self%MG_operator, self%GTM_operator )
    self%uzawa_iterator%verbose = .true.
    self%uzawa_iterator%atol = 1.0d-10
    !self%uzawa_iterator%n_maxiter=5000
  contains
    function profile_m0( x, component)
      sll_real64 :: profile_m0
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_m0 = product(self%Lx) * self%profile%rho_0( x(1) )/self%profile%T_e( x(1) )

    end function profile_m0

    function profile_0( x, component)
      sll_real64 :: profile_0
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_0 = product(self%Lx) 

    end function profile_0


  end subroutine init_3d_fem


  !> Initialization from nml file
  subroutine init_from_file_3d_fem( self, domain, n_dofs, s_deg_0, nml_file, adiabatic_electrons, profile )
    class(sll_t_maxwell_3d_fem), intent(inout) :: self !< Maxwell solver class
    sll_real64, intent(in) :: domain(3,2)     !< xmin, xmax
    sll_int32, intent(in) :: n_dofs(3)  !< number of degrees of freedom (here number of cells and grid points)
    sll_int32, intent(in) :: s_deg_0(3) !< highest spline degree
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
          print*, 'sll_m_maxwell_3d_fem: Input file does not exist. Set default tolerance.'
          open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
          write(file_id, *) 'mass solver tolerance:', 1d-12
          write(file_id, *) 'poisson solver tolerance:', 1d-12
          close(file_id)
       end if
       call self%init( domain, n_dofs, s_deg_0 )
    else
       read(input_file, output,IOStat=io_stat0)
       read(input_file, maxwell_solver,IOStat=io_stat)
       read(input_file, time_solver,IOStat=io_stat1)
       if (io_stat /= 0) then
          if (rank == 0 ) then
             print*, 'sll_m_maxwell_3d_fem: Input parameter does not exist. Set default tolerance.'
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'mass solver tolerance:', 1d-12
             write(file_id, *) 'poisson solver tolerance:', 1d-12
             close(file_id)
          end if
          call self%init( domain, n_dofs, s_deg_0 )
       else
          if (rank == 0 ) then
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'mass solver tolerance:', mass_tolerance
             write(file_id, *) 'poisson solver tolerance:', poisson_tolerance
             close(file_id)
          end if
          call self%init( domain, n_dofs, s_deg_0, mass_tolerance, poisson_tolerance, maxwell_tolerance )
       end if
       close(input_file)
    end if

  end subroutine init_from_file_3d_fem


  !> Finalization
  subroutine free_3d_fem(self)
    class(sll_t_maxwell_3d_fem) :: self !< Maxwell solver class
    !local variable
    sll_int32 :: j

    !call self%poisson_fft%free()
    call self%poisson_solver%free()
    call self%poisson_operator%free()
    call self%poisson_matrix%free
    do j=1, 3
       call self%mass_1_solver(j)%free()
       call self%mass_2_solver(j)%free()
       call self%mass1d_solver(1,j)%free()
       call self%mass1d_solver(2,j)%free()
       call self%mass1d(1,j)%free()
       call self%mass1d(2,j)%free()
       call self%mass1d(3,j)%free()
    end do
    call self%mass1_solver%free()
    call self%mass1_operator%free()
    call self%mass2_operator%free()
    do j=1, 3
       call self%mass1(j)%free()
       call self%mass2(j)%free()
    end do

    call self%linear_solver_schur_eb%free()
    call self%linear_op_schur_eb%free()

    deallocate( self%work3d )
    deallocate( self%work0 )
    deallocate( self%work )
    deallocate( self%work2 )
    deallocate( self%work_d1 ) 
    deallocate( self%work_d2_in ) 
    deallocate( self%work_d2_out ) 
    deallocate( self%work_d3_in ) 
    deallocate( self%work_d3_out )

  end subroutine free_3d_fem


  !> Multiply by dicrete gradient matrix
  subroutine multiply_g( self, field_in, field_out )
    class(sll_t_maxwell_3d_fem) :: self  !< Maxwell solver class
    sll_real64, intent( in    )   :: field_in(:) !< field_in
    sll_real64, intent(   out )   :: field_out(:)! G*field_in

    call sll_s_multiply_g(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_g


  !> Multiply by transpose of dicrete gradient matrix
  subroutine multiply_gt(self, field_in, field_out)
    class(sll_t_maxwell_3d_fem)  :: self !< Maxwell solver class
    sll_real64, intent( in    )  :: field_in(:)  !< field_in
    sll_real64, intent(   out )  :: field_out(:) !< G^T*field_in

    call sll_s_multiply_gt(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_gt


  !> Multiply by discrete curl matrix
  subroutine multiply_c(self, field_in, field_out)
    class(sll_t_maxwell_3d_fem)  :: self         !< Maxwell solver class
    sll_real64, intent( in    )  :: field_in(:)  !< field_in
    sll_real64, intent(   out )  :: field_out(:) !< C*field_in

    call sll_s_multiply_c(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_c


  !> Multiply by transpose of discrete curl matrix
  subroutine multiply_ct(self, field_in, field_out)
    class(sll_t_maxwell_3d_fem)  :: self !< Maxwell solver class
    sll_real64, intent( in    )  :: field_in(:)  !< field_in
    sll_real64, intent(   out )  :: field_out(:) !< C^T*field_in

    call sll_s_multiply_ct(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_ct


  !> Multiply by the mass matrix 
  subroutine multiply_mass_3d_fem( self, deg, coefs_in, coefs_out )
    class(sll_t_maxwell_3d_fem) :: self         !< Maxwell solver class
    sll_int32,  intent( in    )   :: deg(:)       !< \a deg/form specifies if we multiply the mass to a  1- or 2-form or a mix of both
    sll_real64, intent( in    )   :: coefs_in(:)  !< Coefficient for each DoF
    sll_real64, intent(   out )   :: coefs_out(:) !< Coefficient for each DoF


    if( size(deg) ==1 )then
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
    else if( size(deg) == 3 ) then
       call multiply_mass_3dkron(  self, deg, coefs_in, coefs_out )
    end if

  end subroutine multiply_mass_3d_fem


  !> Multiply by the mass matrix 
  subroutine multiply_mass_3dkron(  self, deg, coefs_in, coefs_out )
    class(sll_t_maxwell_3d_fem) :: self !< Maxwell solver class
    sll_int32,  intent( in   )  :: deg(:) !< \a deg(i) specifies the degree of the 1d mass matrix in dimension \a i (Note: 1 for 0-form, 2 for 1-form, 3 for 0-1-form mix)
    sll_real64, intent( in   )  :: coefs_in(:)
    sll_real64, intent(  out )  :: coefs_out(:)  
    ! Local variables
    sll_int32 :: i,j,k,istart,iend

    if( deg(1) == 0 ) then
       call self%mass0%dot( coefs_in, coefs_out )
    else
       istart = 1
       iend = self%n_dofs(1)
       do k=1,self%n_dofs(3)
          do j=1,self%n_dofs(2)

             call self%mass1d(deg(1),1)%dot( coefs_in(istart:iend), self%work_d1 )
             self%work3d(:,j,k) = self%work_d1
             istart = iend+1
             iend = iend + self%n_dofs(1)
          end do
       end do

       do k=1,self%n_dofs(3)
          do i =1,self%n_dofs(1)
             self%work_d2_in = self%work3d(i,:,k)
             call self%mass1d(deg(2),2)%dot( self%work_d2_in, self%work_d2_out )
             self%work3d(i,:,k) = self%work_d2_out
          end do
       end do

       istart = 1
       do j=1,self%n_dofs(2)
          do i =1,self%n_dofs(1)
             self%work_d3_in = self%work3d(i,j,:)
             call self%mass1d(deg(3),3)%dot( self%work_d3_in, self%work_d3_out )
             do k=1,self%n_dofs(3)
                coefs_out(istart+(k-1)*self%n_dofs(1)*self%n_dofs(2)) = self%work_d3_out(k)
             end do
             istart = istart +1
          end do
       end do
    end if


  end subroutine multiply_mass_3dkron


  !> Multiply by the inverse mass matrix 
  subroutine multiply_mass_inverse_3dkron(  self, form, coefs_in, coefs_out )
    class(sll_t_maxwell_3d_fem) :: self !< Maxwell solver class
    sll_int32,  intent( in   )  :: form !< \a form specifies the form (Note: 0 for 0-form, 1 for 1-form, 2 for 2-form, 3 for 0-1-form mix)
    sll_real64, intent( in   )  :: coefs_in(:) !< Coefficient for each DoF
    sll_real64, intent(  out )  :: coefs_out(:) !< Coefficient for each DoF
    ! Local variables
    sll_int32 :: comp, istart, iend

    select case(form)
    case(1)
       do comp=1,3
          istart = 1+(comp-1)*self%n_total
          iend =  comp*self%n_total
          call self%mass_1_solver(comp)%solve( coefs_in(istart:iend), coefs_out(istart:iend) )
       end do
    case(2)
       do comp=1,3
          istart = 1+(comp-1)*self%n_total
          iend =  comp*self%n_total
          call self%mass_2_solver(comp)%solve( coefs_in(istart:iend), coefs_out(istart:iend) )
       end do
    case default
       print*, 'multiply inverse mass for other form not yet implemented'
       stop
    end select


  end subroutine multiply_mass_inverse_3dkron


  !> Compute field energy
  subroutine compute_field_energy( self, efield_dofs, bfield_dofs, energy)
    class(sll_t_maxwell_3d_fem)  :: self !< Maxwell solver class
    sll_real64, intent( in    )          :: efield_dofs(:) !< E
    sll_real64, intent( in    )          :: bfield_dofs(:) !< B
    sll_real64, intent(   out )          :: energy !< field energy
    !local variables
    sll_real64 :: field_energy(6)


    field_energy(1) = self%l2norm_squared &
         ( efield_dofs(1:self%n_total), 1, 1 )
    field_energy(2) = self%l2norm_squared &
         ( efield_dofs(self%n_total+1:2*self%n_total), 1, 2 )
    field_energy(3) = self%l2norm_squared &
         ( efield_dofs(2*self%n_total+1:3*self%n_total), 1, 3 )
    field_energy(4) = self%l2norm_squared &
         ( bfield_dofs(1:self%n_total), 2, 1 )
    field_energy(5) = self%l2norm_squared &
         ( bfield_dofs(self%n_total+1:2*self%n_total), 2, 2 )
    field_energy(6) =self%l2norm_squared &
         ( bfield_dofs(2*self%n_total+1:3*self%n_total), 2, 3 )


    energy = 0.5_f64*sum(field_energy) 

  end subroutine compute_field_energy


end module sll_m_maxwell_3d_fem
