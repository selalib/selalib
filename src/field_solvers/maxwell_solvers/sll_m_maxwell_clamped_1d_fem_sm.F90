!> @ingroup maxwell_solvers
!> @brief
!> Solve Maxwell's equations with boundary conditions in 1D based on spline FEM, version based on sparse matrices
!> @author
!> Benedikt Perse

module sll_m_maxwell_clamped_1d_fem_sm
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
       sll_f_get_collective_rank, &
       sll_v_world_collective

  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_linear_operator_poisson_clamped_1d

  use sll_m_linear_operator_schur_eb_cl_1d, only : &
       sll_t_linear_operator_schur_eb_cl_1d

  use sll_m_linear_solver_cg, only : &
       sll_t_linear_solver_cg

  use sll_m_linear_solver_mgmres, only : &
       sll_t_linear_solver_mgmres

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_maxwell_1d_base, only: &
       sll_i_function_1d_real64, &
       sll_c_maxwell_1d_base

  use sll_m_spline_fem_utilities, only : &
       sll_s_spline_fem_mass_line, &
       sll_s_spline_fem_mass_line_boundary

  use sll_m_spline_fem_utilities_sparse, only : &
       sll_s_spline_fem_mass1d, &
       sll_s_spline_fem_mass1d_clamped

  use sll_m_spline_fem_utilities_3d_clamped, only: &
       sll_s_multiply_g_clamped_1d, &
       sll_s_multiply_gt_clamped_1d

  use sll_m_splines_pp, only: &
       sll_t_spline_pp_1d, &
       sll_s_spline_pp_init_1d, &
       sll_s_spline_pp_free_1d, &
       sll_f_spline_pp_horner_1d

  implicit none

  public :: &
       sll_t_maxwell_clamped_1d_fem_sm

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_c_maxwell_1d_base) :: sll_t_maxwell_clamped_1d_fem_sm

     type(sll_t_spline_pp_1d) :: spline0_pp !< spline for 0-form
     type(sll_t_spline_pp_1d) :: spline1_pp !< spline for 1-from
     sll_real64, allocatable :: work1(:)  !< scratch data
     sll_real64, allocatable :: work0(:)  !< scratch data
     sll_real64, allocatable :: work01(:)  !< scratch data
     type(sll_t_matrix_csr)  :: mass0 !< spline mass matrix for 0-form
     type(sll_t_matrix_csr)  :: mass1 !< spline mass matrix for 1-form
     type(sll_t_linear_solver_cg) :: linear_solver_mass0 !< linear solver to invert 0-form mass matrix
     type(sll_t_linear_solver_cg) :: linear_solver_mass1 !< linear solver to invert 1-form mass matrix

     type(sll_t_linear_operator_poisson_clamped_1d) :: poisson_matrix !< Poisson matrix 
     type(sll_t_linear_solver_cg) :: poisson_solver !< CG solver to invert Poisson matrix

     type( sll_t_linear_operator_schur_eb_cl_1d ) :: linear_op_schur_eb !< Schur complement operator for advect_eb
     type( sll_t_linear_solver_mgmres )        :: linear_solver_schur_eb !< Schur complement solver for advect_eb

   contains
     procedure :: &
          compute_E_from_B => sll_s_compute_e_from_b_1d_fem_sm !< Solve E and B part of Amperes law with B constant in time
     procedure :: &
          compute_B_from_E => sll_s_compute_b_from_e_1d_fem_sm !< Solve Faraday equation with E constant in time
     procedure :: &
          compute_curl_part => sll_s_compute_curl_part_1d_fem_sm !< Solve source-free Maxwell's equations
     procedure :: &
          compute_E_from_rho => sll_s_compute_e_from_rho_1d_fem_sm !< Solve E from rho using Poisson
     procedure :: &
          compute_rho_from_e => compute_rho_from_e_1d_fem_sm !< Compute rho from E 
     procedure :: &
          compute_E_from_j => compute_E_from_j_1d_fem_sm !< Compute E from the current j
     procedure :: &
          compute_phi_from_rho => compute_phi_from_rho_1d_fem_sm !< Compute phi from rho (by solving the quasi-neutrality equation)
     procedure :: &
          compute_phi_from_j => compute_phi_from_j_1d_fem_sm !< Compute phi from j (dynamic of quasi-neutrality equation for adiabatic electrons)
     procedure :: &
          compute_rhs_from_function => sll_s_compute_rhs_fem_sm !< Compute integral over given function tested by the basis 
     procedure :: &
          L2projection => L2projection_1d_fem_sm !< Compute L_2 projection of a given function
     procedure :: &
          L2norm_squared => L2norm_squared_1d_fem_sm !< Compute the square of the L2 norm of a given vector
     procedure :: &
          inner_product => inner_product_1d_fem_sm !< Inner product of two dof-vectors with mass matrix
     procedure :: &
          init => init_1d_fem_sm !< Initialize the Maxwell class
     procedure :: &
          init_from_file => init_from_file_1d_fem_sm  !< Initialize the Maxwell class with parameters read from nml-file
     procedure :: &
          free => free_1d_fem_sm !< Free Maxwell class
     procedure :: &
          multiply_g !< Multiplication with gradient matrix 
     procedure :: &
          multiply_gt !< Multiplication with divergence matrix
     procedure :: &
          multiply_mass => multiply_mass_1d_fem_sm !< Product with the mass matrix
     procedure :: &
          invert_mass => invert_mass_1d_fem_sm !< Invert mass matrix
     procedure :: &
          compute_field_energy !< Compute field energy


  end type sll_t_maxwell_clamped_1d_fem_sm

contains


  !> compute Ey from Bz using weak Ampere formulation 
  subroutine sll_s_compute_e_from_b_1d_fem_sm(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_clamped_1d_fem_sm) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)     :: delta_t   !< Time step
    sll_real64, intent(in)     :: field_in(:)  !< Bz
    sll_real64, intent(inout)  :: field_out(:)  !< Ey

    ! Compute potential weak curl of Bz
    call self%multiply_mass( field_in, self%work1, self%s_deg_1 )
    call self%multiply_gt( self%work1, self%work0 )
    call self%invert_mass( self%work0, self%work01, self%s_deg_0 )

    ! Update Ey from self value
    field_out = field_out + delta_t*self%work01

  end subroutine sll_s_compute_e_from_b_1d_fem_sm


  !> Compute Bz from Ey using strong 1D Faraday equation for spline coefficients
  !> $B_z^{new}(x_j) = B_z^{old}(x_j) - \frac{\Delta t}{\Delta x} (E_y(x_j) - E_y(x_{j-1})  $
  subroutine sll_s_compute_b_from_e_1d_fem_sm(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_clamped_1d_fem_sm)  :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)     :: delta_t !< Time step
    sll_real64, intent(in)     :: field_in(:)  ! Ey
    sll_real64, intent(inout)  :: field_out(:) ! Bz 

    call self%multiply_g(field_in, self%work1)
    ! Update Bz from self value
    field_out = field_out - delta_t * self%work1 

  end subroutine sll_s_compute_b_from_e_1d_fem_sm


  !> Solve curl part of Maxwell's equations
  subroutine sll_s_compute_curl_part_1d_fem_sm( self, delta_t, efield, bfield, betar )
    class(sll_t_maxwell_clamped_1d_fem_sm) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)     :: delta_t   !< Time step
    sll_real64, intent(inout)  :: efield(:)  !< Ey
    sll_real64, intent(inout)  :: bfield(:)  !< Bz
    sll_real64, optional       :: betar      !< 1/beta
    !local variables
    sll_real64 :: factor

    if( present(betar) ) then
       factor = betar
    else
       factor = 1._f64
    end if

    self%work1 = 0._f64
    self%work0 = 0._f64
    self%work01 = 0._f64

    ! Compute D^T M2 b
    call self%multiply_mass( bfield, self%work1, self%s_deg_1 )
    call self%multiply_gt( self%work1, self%work0 )

    self%linear_op_schur_eb%sign = -delta_t**2*0.25_f64*factor
    call self%linear_op_schur_eb%dot( efield, self%work01 )
    self%work01 = self%work01 + delta_t*factor*self%work0

    ! Save efield dofs from previous time step for B field update
    self%work0 = efield

    ! Invert Schur complement matrix
    self%linear_op_schur_eb%sign = delta_t**2*0.25_f64*factor
    call self%linear_solver_schur_eb%set_guess( efield )
    call self%linear_solver_schur_eb%solve( self%work01, efield )

    ! Update B field
    self%work0 = self%work0 + efield
    call self%compute_B_from_E( delta_t*0.5_f64, self%work0, bfield )

  end subroutine sll_s_compute_curl_part_1d_fem_sm


  !> compute e from rho using weak Poisson's equation ( rho = G^T M_1 G \phi, e = G \phi ) 
  subroutine sll_s_compute_e_from_rho_1d_fem_sm(self, field_in, field_out )       
    class(sll_t_maxwell_clamped_1d_fem_sm) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)  :: field_in(:) !< rho
    sll_real64, intent(out) :: field_out(:) !< E

    call self%poisson_solver%solve( field_in, self%work0 )
    call multiply_g( self, self%work0, field_out )
    field_out = -field_out

  end subroutine sll_s_compute_e_from_rho_1d_fem_sm


  !> Compute rho from Gauss law for given efield
  subroutine compute_rho_from_e_1d_fem_sm(self, field_in, field_out)
    class(sll_t_maxwell_clamped_1d_fem_sm)  :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)     :: field_in(:)  ! E
    sll_real64, intent(out)  :: field_out(:) ! rho

    call self%multiply_mass( field_in, self%work1, self%s_deg_1 )
    call multiply_gt( self, self%work1, field_out )
    field_out = - field_out 

  end subroutine compute_rho_from_e_1d_fem_sm


  !> Compute E_i from j_i integrated over the time interval using weak Ampere formulation
  subroutine compute_e_from_j_1d_fem_sm(self, current, component, E)
    class(sll_t_maxwell_clamped_1d_fem_sm)             :: self !< Maxwell_Clamped solver class
    sll_real64,dimension(:),intent(in)    :: current !< Component \a component of the current integrated over time interval
    sll_int32, intent(in)                 :: component !< Component of the Efield to be computed
    sll_real64,dimension(:),intent(inout) :: E !< Updated electric field

    ! Multiply by inverse mass matrix  
    if (component == 1) then
       call self%invert_mass( current, self%work1, self%s_deg_1  )
       ! Update the electric field and scale
       E = E -  self%work1
    elseif (component == 2) then
       call self%invert_mass( current, self%work0, self%s_deg_0  )
       E = E -  self%work0
    else
       print*, 'Component ', component, 'not implemented in compute_E_from_j_1d_fem_sm.'
    end if

  end subroutine compute_e_from_j_1d_fem_sm


  !> For model with adiabatic electrons
  subroutine compute_phi_from_rho_1d_fem_sm( self, in, phi, efield )
    class(sll_t_maxwell_clamped_1d_fem_sm) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)                     :: in(:) !< rho 
    sll_real64, intent(out)                    :: phi(:) !< phi 
    sll_real64, intent(out)                    :: efield(:) !< E

    ! Compute phi by inverting the mass matrix
    call self%invert_mass( in, phi, self%s_deg_0 )

    ! Compute the degrees of freedom of the electric field as -G phi
    call multiply_g( self, phi, efield )
    efield = -efield

  end subroutine compute_phi_from_rho_1d_fem_sm


  !> For model with adiabatic electrons
  subroutine compute_phi_from_j_1d_fem_sm( self, in, phi, efield )
    class(sll_t_maxwell_clamped_1d_fem_sm) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)                     :: in(:) !< Current integrated over time interval
    sll_real64, intent(out)                    :: phi(:) !< phi 
    sll_real64, intent(out)                    :: efield(:) !< E

    ! Compute divergence of the current (G^T current) (assuming a 1v model)
    call multiply_gt( self, in, self%work01 )

    ! Compute phi by inverting the mass matrix
    call self%invert_mass( self%work01, self%work0, self%s_deg_0 )

    phi = phi + self%work0

    ! Compute the degrees of freedom of the electric field as -G phi
    call multiply_g( self, phi, efield )
    efield = -efield

  end subroutine compute_phi_from_j_1d_fem_sm


  !> Compute the FEM right-hand-side for a given function f and clamped splines of given degree
  !> Its components are $\int f N_i dx$ where $N_i$ is the B-spline starting at $x_i$ 
  subroutine sll_s_compute_rhs_fem_sm(self, func, degree, coefs_dofs)
    class(sll_t_maxwell_clamped_1d_fem_sm) :: self !< Maxwell_Clamped solver class
    procedure(sll_i_function_1d_real64) :: func !< Function
    sll_int32, intent(in) :: degree !< Specify the degree of the basis functions
    sll_real64, intent(out) :: coefs_dofs(:)  !< Finite Element right-hand-side
    ! local variables
    sll_int32 :: i,j,quad,q
    sll_real64, allocatable :: xw_gauss(:,:)
    sll_real64, allocatable :: bspl(:,:)

    q = 2*degree+1
    allocate( xw_gauss(2,q) )
    allocate( bspl(degree+1,q) )
    ! rescale on [0,1] for compatibility with B-splines
    xw_gauss = sll_f_gauss_legendre_points_and_weights(q, 0.0_f64, 1.0_f64)

    if( degree == self%s_deg_0 ) then
       coefs_dofs = 0._f64
       !first degree-1 cells are different
       !loop over cells
       do i = 1, degree-1
          !loop over splines
          do j = 1, degree+1
             ! loop over Gauss points
             do quad = 1, q
                coefs_dofs(i-1+j) = coefs_dofs(i-1+j) + &
                     self%delta_x * xw_gauss(2,quad) * func( self%delta_x * (xw_gauss(1,quad)+ real(i - 1,f64)) )   *&
                     sll_f_spline_pp_horner_1d( degree, self%spline0_pp%poly_coeffs_boundary_left(:,:,i), xw_gauss(1,quad), j)
             enddo
          enddo
       enddo
       bspl=0._f64
       ! Compute bsplines at gauss_points
       do quad=1,q
          do j=1, degree+1
             bspl(j,quad) = sll_f_spline_pp_horner_1d( degree, self%spline0_pp%poly_coeffs, xw_gauss(1,quad), j)
          end do
       end do

       ! Compute coefs_dofs = int f(x)N_i(x)
       !loop over cells
       do i = degree, self%n_cells-degree+1
          ! loop over splines
          do j = 1, degree+1
             ! loop over Gauss points
             do quad = 1, q
                coefs_dofs(i-1+j) = coefs_dofs(i-1+j) + &
                     self%delta_x * xw_gauss(2,quad) * func(self%delta_x * (xw_gauss(1,quad) + real(i - 1,f64)) ) * bspl(j,quad)
             enddo
          enddo
       enddo

       !last degree-1 cells are different
       !loop over cells
       do i = self%n_cells-degree+2, self%n_cells
          !loop over splines
          do j = 1, degree+1
             ! loop over Gauss points
             do quad = 1, q
                coefs_dofs(i-1+j) = coefs_dofs(i-1+j) + &
                     self%delta_x * xw_gauss(2,quad) * func( self%delta_x * (xw_gauss(1,quad)+real(i-1,f64)) ) *&
                     sll_f_spline_pp_horner_1d( degree, self%spline0_pp%poly_coeffs_boundary_right(:,:,i-self%n_cells+degree-1), xw_gauss(1,quad), j)
             enddo
          enddo
       end do

    else if (degree == self%s_deg_1)then

       coefs_dofs = 0._f64
       !first degree-1 cells are different
       !loop over cells
       do i = 1, degree-1
          !loop over splines
          do j = 1, degree+1
             ! loop over Gauss points
             do quad = 1, q
                coefs_dofs(i-1+j) = coefs_dofs(i-1+j) + &
                     self%delta_x * xw_gauss(2,quad)*func( self%delta_x * (xw_gauss(1,quad)+ real(i - 1,f64)) )   *&
                     sll_f_spline_pp_horner_1d( degree, self%spline1_pp%poly_coeffs_boundary_left(:,:,i), xw_gauss(1,quad), j)
             enddo
          enddo
       enddo

       bspl=0._f64
       ! Compute bsplines at gauss_points
       do quad=1, q
          do j=1, degree+1
             bspl(j,quad) = sll_f_spline_pp_horner_1d( degree, self%spline1_pp%poly_coeffs, xw_gauss(1,quad), j)
          end do
       end do

       ! Compute coefs_dofs = int f(x)N_i(x)
       !loop over cells
       do i = max(1,degree), min(self%n_cells,self%n_cells-degree+1)
          !loop over splines
          do j = 1, degree+1
             ! loop over Gauss points
             do quad = 1, q
                coefs_dofs(i-1+j) = coefs_dofs(i-1+j) + &
                     self%delta_x * xw_gauss(2,quad)*func(self%delta_x * (xw_gauss(1,quad) + real(i - 1,f64)) ) * bspl(j,quad)
             enddo
          enddo
       enddo

       !last degree-1 cells are different
       !loop over cells
       do i = self%n_cells-degree+2, self%n_cells
          !loop over splines
          do j = 1, degree+1
             ! loop over Gauss points
             do quad = 1, q
                coefs_dofs(i-1+j) = coefs_dofs(i-1+j) + &
                     self%delta_x * xw_gauss(2,quad)*func( self%delta_x * (xw_gauss(1,quad)+real(i-1,f64)) ) *&
                     sll_f_spline_pp_horner_1d( degree, self%spline1_pp%poly_coeffs_boundary_right(:,:,i-self%n_cells+degree-1), xw_gauss(1,quad), j)
             enddo
          enddo
       end do
    end if

  end subroutine sll_s_compute_rhs_fem_sm


  !> Compute the L2 projection of a given function f on periodic splines of given degree
  subroutine L2projection_1d_fem_sm(self, func, degree, coefs_dofs)
    class(sll_t_maxwell_clamped_1d_fem_sm) :: self !< Maxwell_Clamped solver class
    procedure(sll_i_function_1d_real64) :: func !< Function
    sll_int32, intent(in) :: degree !< Specify the degree of the basis functions
    sll_real64, intent(out) :: coefs_dofs(:)  ! spline coefficients of projection

    ! Multiply by inverse mass matrix 
    if (degree == self%s_deg_0) then
       ! Compute right-hand-side
       call sll_s_compute_rhs_fem_sm(self, func, degree, self%work0)
       call self%invert_mass( self%work0, coefs_dofs, degree )

    elseif  (degree == self%s_deg_1) then
       ! Compute right-hand-side
       call sll_s_compute_rhs_fem_sm(self, func, degree, self%work1)
       call self%invert_mass( self%work1, coefs_dofs, degree )
    else
       print*, 'degree ', degree, 'not availlable in maxwell_clamped_1d_fem_sm object' 
    endif

  end subroutine L2projection_1d_fem_sm


  !> Compute square of the L2norm 
  function L2norm_squared_1d_fem_sm(self, coefs_dofs, degree) result (r)
    class(sll_t_maxwell_clamped_1d_fem_sm) :: self !< Maxwell_Clamped solver class
    sll_real64 :: coefs_dofs(:) !< Coefficient for each DoF
    sll_int32  :: degree !< Specify the degree of the basis functions
    sll_real64 :: r !< Result: squared L2 norm

    ! Multiply coefficients by mass matrix 
    if (degree == self%s_deg_0 ) then
       call self%multiply_mass( coefs_dofs, self%work0, self%s_deg_0 )
       ! Multiply by the coefficients from the left (inner product)
       r = sum(coefs_dofs*self%work0)
    elseif (degree == self%s_deg_1) then
       call self%multiply_mass( coefs_dofs, self%work1, self%s_deg_1 )
       ! Multiply by the coefficients from the left (inner product)
       r = sum(coefs_dofs*self%work1)
    end if

  end function L2norm_squared_1d_fem_sm


  !> Compute inner product
  function inner_product_1d_fem_sm(self, coefs1_dofs, coefs2_dofs, degree, degree2) result (r)
    class(sll_t_maxwell_clamped_1d_fem_sm) :: self !< Maxwell_Clamped solver class
    sll_real64 :: coefs1_dofs(:) !< Coefficient for each DoF
    sll_real64 :: coefs2_dofs(:) !< Coefficient for each DoF
    sll_int32  :: degree !< Specify the degree of the basis functions
    sll_int32, optional  :: degree2 !< Specify the degree of the basis functions
    sll_real64 :: r !< Result: squared L2 norm

    ! Multiply coefficients by mass matrix (use diagonalization FFT and mass matrix eigenvalues)
    if (present(degree2)) then
       if (degree == degree2) then
          if (degree == self%s_deg_0 ) then
             call self%multiply_mass( coefs2_dofs, self%work0, self%s_deg_0 )
             ! Multiply by the coefficients from the left (inner product)
             r = sum(coefs1_dofs*self%work0)
          elseif (degree == self%s_deg_1) then
             call self%multiply_mass( coefs2_dofs, self%work1, self%s_deg_1 )
             ! Multiply by the coefficients from the left (inner product)
             r = sum(coefs1_dofs*self%work1)
          end if
       else
          if (degree == self%s_deg_0) then
             call self%multiply_mass( coefs2_dofs, self%work01, 10 )
             ! Multiply by the coefficients from the left (inner product)
             r = sum(coefs1_dofs*self%work01)
          else
             call self%multiply_mass( coefs1_dofs, self%work01, 10 )
             ! Multiply by the coefficients from the left (inner product)
             r = sum(coefs2_dofs*self%work01)
          end if
       end if
    else
       if (degree == self%s_deg_0 ) then
          call self%multiply_mass( coefs2_dofs, self%work0, self%s_deg_0 )
          ! Multiply by the coefficients from the left (inner product)
          r = sum(coefs1_dofs*self%work0)
       elseif (degree == self%s_deg_1) then
          call self%multiply_mass( coefs2_dofs, self%work1, self%s_deg_1 )
          ! Multiply by the coefficients from the left (inner product)
          r = sum(coefs1_dofs*self%work1)
       end if
    end if

  end function inner_product_1d_fem_sm


  !> Initialization
  subroutine init_1d_fem_sm( self, domain, n_cells, s_deg_0, boundary, mass_tolerance, poisson_tolerance, solver_tolerance )
    class(sll_t_maxwell_clamped_1d_fem_sm), intent(out) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in) :: domain(2)     !< xmin, xmax
    sll_int32, intent(in) :: n_cells  !< number of degrees of freedom (here number of cells and grid points)
    sll_int32, intent(in) :: s_deg_0 !< highest spline degree
    sll_int32, intent(in) :: boundary !< field boundary conditions
    sll_real64, intent(in), optional :: mass_tolerance !< tolerance for mass solver
    sll_real64, intent(in), optional :: poisson_tolerance !< tolerance for Poisson solver
    sll_real64, intent(in), optional :: solver_tolerance !< tolerance for Schur complement solver
    ! local variables
    sll_int32 :: ierr
    sll_real64 :: mass_line_0(s_deg_0+1), mass_line_1(s_deg_0)
    sll_real64 :: mass_line_b0((s_deg_0+1)*s_deg_0), mass_line_b1(s_deg_0*(s_deg_0-1))

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

    self%n_cells = n_cells
    self%n_dofs0 = n_cells+s_deg_0
    self%n_dofs1 = n_cells+s_deg_0-1
    self%Lx = domain(2) - domain(1)
    self%delta_x = self%Lx /real(n_cells,f64)
    self%s_deg_0 = s_deg_0
    self%s_deg_1 = s_deg_0 - 1

    call sll_s_spline_pp_init_1d( self%spline0_pp, self%s_deg_0, n_cells, boundary)
    call sll_s_spline_pp_init_1d( self%spline1_pp, self%s_deg_1, n_cells, boundary)

    SLL_ALLOCATE(self%work1(self%n_dofs1),ierr)
    SLL_ALLOCATE(self%work0(self%n_dofs0),ierr)
    SLL_ALLOCATE(self%work01(self%n_dofs0),ierr)


    ! Sparse matrices
    call sll_s_spline_fem_mass_line ( self%s_deg_0, mass_line_0 )
    call sll_s_spline_fem_mass_line_boundary ( self%s_deg_0, self%spline0_pp, mass_line_b0 )
    !scale with delta_x=L/n_cells
    mass_line_0= mass_line_0*self%delta_x
    mass_line_b0= mass_line_b0*self%delta_x
    call sll_s_spline_fem_mass1d_clamped( n_cells, self%s_deg_0, mass_line_0, mass_line_b0, self%mass0 )

    if(self%s_deg_1 > 0)then
       call sll_s_spline_fem_mass_line ( self%s_deg_1, mass_line_1 )
       call sll_s_spline_fem_mass_line_boundary ( self%s_deg_1, self%spline1_pp, mass_line_b1 )
       !scale with delta_x=L/n_cells
       mass_line_1= mass_line_1*self%delta_x
       mass_line_b1= mass_line_b1*self%delta_x
       call sll_s_spline_fem_mass1d_clamped( n_cells, self%s_deg_1, mass_line_1, mass_line_b1, self%mass1 )
    else
       call sll_s_spline_fem_mass_line ( self%s_deg_1, mass_line_1 )
       mass_line_1 = mass_line_1*self%delta_x
       call sll_s_spline_fem_mass1d( n_cells, self%s_deg_1, mass_line_1, self%mass1 )
    end if
    call self%linear_solver_mass0%create( self%mass0 )
    call self%linear_solver_mass1%create( self%mass1 )
    self%linear_solver_mass0%atol = self%mass_solver_tolerance
    self%linear_solver_mass1%atol = self%mass_solver_tolerance
    !self%linear_solver_mass0%verbose = .true.
    !self%linear_solver_mass1%verbose = .true.

    ! Poisson solver
    call self%poisson_matrix%create( self%mass1, self%s_deg_0, self%n_dofs0, self%delta_x)
    call self%poisson_solver%create( self%poisson_matrix )
    self%poisson_solver%atol = self%poisson_solver_tolerance
    !self%poisson_solver%verbose = .true.
    !self%poisson_solver%n_maxiter=40000

    ! Only for Schur complement eb solver
    call self%linear_op_schur_eb%create( self%mass0, self%mass1, self%n_dofs0, self%delta_x, self%s_deg_0 )
    call self%linear_solver_schur_eb%create( self%linear_op_schur_eb )
    self%linear_solver_schur_eb%atol = self%solver_tolerance
    self%linear_solver_schur_eb%rtol = self%solver_tolerance
    !self%linear_solver_schur_eb%verbose = .true.

  contains
    function profile_0( x)
      sll_real64 :: profile_0
      sll_real64, intent(in) :: x

      profile_0 = self%Lx

    end function profile_0
  end subroutine init_1d_fem_sm


  !> Initialization from nml file
  subroutine init_from_file_1d_fem_sm( self, domain, n_cells, s_deg_0, boundary, nml_file)
    class(sll_t_maxwell_clamped_1d_fem_sm), intent(out) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in) :: domain(2)     !< xmin, xmax
    sll_int32, intent(in) :: n_cells  !< number of degrees of freedom (here number of cells and grid points)
    !sll_real64 :: delta_x ! cell size
    sll_int32, intent(in) :: s_deg_0 !< highest spline degree
    sll_int32, intent(in) :: boundary !< field boundary conditions
    character(len=*) :: nml_file !< nml-file
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

    ! Read in solver tolerance
    open(newunit = input_file, file=nml_file, status='old',IOStat=io_stat)
    if (io_stat /= 0) then
       if (rank == 0 ) then
          print*, 'sll_m_maxwell_clamped_1d_fem_sm: Input file does not exist. Set default tolerance.'
          open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
          write(file_id, *) 'mass solver tolerance:', 1d-12
          write(file_id, *) 'poisson solver tolerance:', 1d-12
          close(file_id)
       end if
       call self%init( domain, n_cells, s_deg_0, boundary )
    else
       read(input_file, output,IOStat=io_stat0)
       read(input_file, maxwell_solver,IOStat=io_stat)
       read(input_file, time_solver,IOStat=io_stat1)
       if (io_stat /= 0) then
          if (rank == 0 ) then
             print*, 'sll_m_maxwell_clamped_1d_fem_sm: Input parameter does not exist. Set default tolerance.'
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'mass solver tolerance:', 1d-12
             write(file_id, *) 'poisson solver tolerance:', 1d-12
             close(file_id)
          end if
          call self%init( domain, n_cells, s_deg_0, boundary )
       else
          if (rank == 0 ) then
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'mass solver tolerance:', mass_tolerance
             write(file_id, *) 'poisson solver tolerance:', poisson_tolerance
             close(file_id)
          end if
          call self%init( domain, n_cells, s_deg_0, boundary, mass_tolerance, poisson_tolerance, maxwell_tolerance )
       end if
       close(input_file)
    end if

  end subroutine init_from_file_1d_fem_sm


  !> Finalization
  subroutine free_1d_fem_sm(self)
    class(sll_t_maxwell_clamped_1d_fem_sm) :: self !< Maxwell_Clamped solver class

    deallocate(self%work1)
    deallocate(self%work0)
    deallocate(self%work01)
    call self%linear_solver_mass0%free()
    call self%linear_solver_mass1%free()
    call self%poisson_solver%free()
    call self%poisson_matrix%free()
    call self%linear_solver_schur_eb%free()
    call self%linear_op_schur_eb%free()

    call sll_s_spline_pp_free_1d( self%spline0_pp )
    call sll_s_spline_pp_free_1d( self%spline1_pp )
    
  end subroutine free_1d_fem_sm


  !> Multiply by dicrete gradient matrix
  subroutine multiply_g( self,  in, out)
    class(sll_t_maxwell_clamped_1d_fem_sm), intent(in) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)  :: in(:) !< field_in 
    sll_real64, intent(out) :: out(:) !< G*field_in

    call sll_s_multiply_g_clamped_1d( self%n_dofs0, self%delta_x, self%s_deg_0, in, out)

  end subroutine multiply_g


  !> Multiply by transpose of dicrete gradient matrix
  subroutine multiply_gt( self,  in, out)
    class(sll_t_maxwell_clamped_1d_fem_sm), intent(in) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)  :: in(:) !< field_in 
    sll_real64, intent(out) :: out(:) !< G^T*field_in

    call sll_s_multiply_gt_clamped_1d( self%n_dofs0, self%delta_x, self%s_deg_0, in, out)

  end subroutine multiply_gt


  !> Multiply by the mass matrix 
  subroutine multiply_mass_1d_fem_sm( self,  in, out, degree)
    class(sll_t_maxwell_clamped_1d_fem_sm), intent(inout) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)  :: in(:) !< Coefficient for each DoF
    sll_real64, intent(out) :: out(:) !< Coefficient for each DoF
    sll_int32, intent(in)  :: degree !< Specify the degree of the basis functions

    ! Multiply coefficients by mass matrix 
    if (degree == self%s_deg_0 ) then

       call self%mass0%dot ( in, out)

    elseif (degree == self%s_deg_1) then

       call self%mass1%dot ( in, out)
    else
       print*, 'multiply mass for other form not yet implemented'
       stop
    end if

  end subroutine multiply_mass_1d_fem_sm


  !> Multiply by the inverse mass matrix 
  subroutine invert_mass_1d_fem_sm( self,  in, out, degree)
    class(sll_t_maxwell_clamped_1d_fem_sm), intent(inout) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)  :: in(:) !< Coefficient for each DoF
    sll_real64, intent(out) :: out(:) !< Coefficient for each DoF
    sll_int32, intent(in)  :: degree !< Specify the degree of the basis functions

    ! Multiply coefficients by mass matrix 
    if (degree == self%s_deg_0 ) then

       call self%linear_solver_mass0%solve ( in, out)
       !perfect conductor boundary
       out(1) = 0._f64
       out(self%n_dofs0) = 0._f64
    elseif (degree == self%s_deg_1) then
       call self%linear_solver_mass1%solve ( in, out)
    else
       print*, 'Invert mass for other form not yet implemented'
       stop
    end if

  end subroutine invert_mass_1d_fem_sm

  !> Compute the field energy
  subroutine compute_field_energy( self, efield_dofs1, efield_dofs2, bfield_dofs, energy)
    class(sll_t_maxwell_clamped_1d_fem_sm)  :: self !< Maxwell_Clamped solver class
    sll_real64, intent( in    )          :: efield_dofs1(:) !< Ex
    sll_real64, intent( in    )          :: efield_dofs2(:) !< Ey
    sll_real64, intent( in    )          :: bfield_dofs(:)  !< Bz
    sll_real64, intent(   out )          :: energy          !< field energy
    !local variables
    sll_real64 :: field_energy(3)


    field_energy(1) = self%l2norm_squared &
         ( efield_dofs1(1:self%n_dofs1), self%s_deg_1 )
    field_energy(2) = self%l2norm_squared &
         ( efield_dofs2(1:self%n_dofs0), self%s_deg_0 )
    field_energy(3) =self%l2norm_squared &
         ( bfield_dofs(1:self%n_dofs1), self%s_deg_1 )

    energy = 0.5_f64*sum(field_energy) 

  end subroutine compute_field_energy

  
end module sll_m_maxwell_clamped_1d_fem_sm
