!> @ingroup maxwell_solvers
!> @brief
!> Solve Maxwell's equations in curvilinear coordinates in 1D based on spline FEM, version based on sparse matrices
!> @author
!> Benedikt Perse

module sll_m_maxwell_1d_trafo
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
       sll_f_get_collective_rank, &
       sll_v_world_collective

  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_linear_operator_penalized, only : &
       sll_t_linear_operator_penalized

  use sll_m_linear_operator_poisson_1d, only : &
       sll_t_linear_operator_poisson_1d

  use sll_m_linear_operator_schur_eb_1d, only : &
       sll_t_linear_operator_schur_eb_1d

  use sll_m_linear_solver_cg, only : &
       sll_t_linear_solver_cg

  use sll_m_linear_solver_mgmres, only : &
       sll_t_linear_solver_mgmres

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_maxwell_1d_base, only: &
       sll_i_function_1d_real64, &
       sll_c_maxwell_1d_base

  use sll_m_spline_fem_utilities, only : &
       sll_s_multiply_g_1d, &
       sll_s_multiply_gt_1d

  use sll_m_spline_fem_utilities_sparse, only : &
       sll_s_spline_fem_mass1d_full, &
       sll_s_spline_fem_mixedmass1d_full

  implicit none

  public :: &
       sll_t_maxwell_1d_trafo

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_c_maxwell_1d_base) :: sll_t_maxwell_1d_trafo

     sll_real64, allocatable :: work1(:)   !< scratch data
     sll_real64, allocatable :: work2(:)   !< scratch data
     type(sll_t_matrix_csr)  :: mass0      !< spline mass matrix for 0-form
     type(sll_t_matrix_csr)  :: mass1      !< spline mass matrix for 1-form
     type(sll_t_matrix_csr)  :: mixed_mass !< spline mass matrix for mix 0-1-form
     type(sll_t_linear_solver_cg) :: linear_solver_mass0 !< linear solver to invert 0-form mass matrix
     type(sll_t_linear_solver_cg) :: linear_solver_mass1 !< linear solver to invert 1-form mass matrix

     type(sll_t_linear_operator_poisson_1d) :: poisson_matrix  !< Poisson matrix 
     type(sll_t_linear_operator_penalized) :: poisson_operator !< Poisson matrix with constraint on constant vector
     type(sll_t_linear_solver_cg) :: poisson_solver !< CG solver to invert Poisson matrix
     type( sll_t_linear_operator_schur_eb_1d ) :: linear_op_schur_eb !< Schur complement operator for advect_eb
     type( sll_t_linear_solver_mgmres )        :: linear_solver_schur_eb !< Schur complement solver for advect_eb
     type(sll_t_mapping_3d), pointer :: map !< coordinate transformation
     sll_real64 :: force_sign = 1._f64 !< sign of particle force, + for attraction, - for repulsion

   contains

     procedure :: &
          compute_E_from_B => sll_s_compute_e_from_b_1d_trafo !< Solve E and B part of Amperes law with B constant in time
     procedure :: &
          compute_B_from_E => sll_s_compute_b_from_e_1d_trafo !< Solve Faraday equation with E constant in time
     procedure :: &
          compute_curl_part => sll_s_compute_curl_part_1d_trafo !< Solve source-free Maxwell's equations
     procedure :: &
          compute_E_from_rho => sll_s_compute_e_from_rho_1d_trafo !< Solve E from rho using Poisson
     procedure :: &
          compute_rho_from_e => compute_rho_from_e_1d_trafo !< Compute rho from E
     procedure :: &
          compute_E_from_j => compute_E_from_j_1d_trafo !< Compute E from the current j 
     procedure :: &
          compute_phi_from_rho => compute_phi_from_rho_1d_trafo !< Compute phi from rho (by solving the quasi neutrality equation)
     procedure :: &
          compute_phi_from_j => compute_phi_from_j_1d_trafo !< Compute phi from j (dynamic from of quasineutrality equation for adiabatic electrons)

     procedure :: &
          compute_rhs_from_function => sll_s_compute_rhs_trafo !< Compute integral over given function tested by the basis
     procedure :: &
          L2projection => L2projection_1d_trafo !< Compute L_2 projection of a given function
     procedure :: &
          L2norm_squared => L2norm_squared_1d_trafo !< Compute the square of the L2 norm of a given vector
     procedure :: &
          inner_product => inner_product_1d_trafo !< Inner product of two dof-vectors with mass matrix
     procedure :: &
          init => init_1d_trafo !< Initialize the Maxwell class
     procedure :: &
          init_from_file => init_from_file_1d_trafo !< Initialize the Maxwell class with parameters read from nml-file
     procedure :: &
          free => free_1d_trafo !< Free Maxwell class
     procedure :: &
          multiply_g !< Multiplication with gradient matrix 
     procedure :: &
          multiply_gt !< Multiplication with divergence matrix
     procedure :: &
          multiply_mass => multiply_mass_1d_trafo !< Product with the mass matrix
     procedure :: &
          invert_mass => invert_mass_1d_trafo !< Invert mass matrix

  end type sll_t_maxwell_1d_trafo

contains

  !> compute Ey from Bz using weak Ampere formulation 
  subroutine sll_s_compute_E_from_B_1d_trafo(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_1d_trafo) :: self         !< Maxwell solver class
    sll_real64, intent( in    )   :: delta_t      !< Time step
    sll_real64, intent( in    )   :: field_in(:)  !< Bz
    sll_real64, intent( inout )   :: field_out(:) !< Ey

    call self%multiply_mass( field_in, self%work2, self%s_deg_1 )
    call self%multiply_gt( self%work2, self%work1 )
    self%work2=0._f64
    call self%linear_solver_mass0%solve ( self%work1, self%work2 )
    ! Update ey from self value
    field_out =  field_out + delta_t*self%work2

  end subroutine sll_s_compute_E_from_B_1d_trafo


  !> Compute Bz from Ey using strong 1D Faraday equation for spline coefficients
  !> $B_z^{new}(x_j) = B_z^{old}(x_j) - \frac{\Delta t}{\Delta x} (E_y(x_j) - E_y(x_{j-1})  $
  subroutine sll_s_compute_B_from_E_1d_trafo(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_1d_trafo) :: self         !< Maxwell solver class
    sll_real64, intent( in    )   :: delta_t      !< Time step
    sll_real64, intent( in    )   :: field_in(:)  ! Ey
    sll_real64, intent( inout )   :: field_out(:) ! Bz

    call self%multiply_g(field_in, self%work1)
    ! Update Bz from self value
    field_out = field_out - delta_t * self%work1

  end subroutine sll_s_compute_B_from_E_1d_trafo


  !> Solve curl part of Maxwell's equations
  subroutine sll_s_compute_curl_part_1d_trafo( self, delta_t, efield, bfield, betar )
    class(sll_t_maxwell_1d_trafo) :: self   !< Maxwell solver class
    sll_real64, intent(in)     :: delta_t   !< Time step
    sll_real64, intent(inout)  :: efield(:) !< E
    sll_real64, intent(inout)  :: bfield(:) !< B
    sll_real64, optional       :: betar     !< 1/beta
    !local variables
    sll_real64 :: factor

    if( present(betar) ) then
       factor = betar
    else
       factor = 1._f64
    end if

    self%work1 = 0._f64
    self%work2 = 0._f64

    ! Compute D^T M2 b
    call self%multiply_mass( bfield, self%work1, self%s_deg_1 )
    call self%multiply_gt( self%work1, self%work2 )

    self%linear_op_schur_eb%sign = -delta_t**2*0.25_f64*factor
    call self%linear_op_schur_eb%dot( efield, self%work1 )
    self%work1 = self%work1 + delta_t*factor*self%work2

    ! Save efield dofs from previous time step for B field update
    self%work2 = efield

    ! Invert Schur complement matrix
    self%linear_op_schur_eb%sign = delta_t**2*0.25_f64*factor
    call self%linear_solver_schur_eb%set_guess( efield )
    call self%linear_solver_schur_eb%solve( self%work1, efield )

    ! Update B field
    self%work2 = self%work2 + efield
    call self%compute_B_from_E( delta_t*0.5_f64, self%work2, bfield )

  end subroutine sll_s_compute_curl_part_1d_trafo


  !> compute e from rho using weak Poisson's equation ( rho = G^T M_1 G \phi, e = -G \phi ) 
  subroutine sll_s_compute_e_from_rho_1d_trafo(self, field_in, field_out )       
    class(sll_t_maxwell_1d_trafo) :: self         !< Maxwell solver class
    sll_real64, intent( in    )   :: field_in(:)  !< rho
    sll_real64, intent(   out )   :: field_out(:) !< E

    call self%poisson_solver%solve( field_in, self%work1 )
    call multiply_g( self, self%work1, field_out )
    field_out = -field_out

  end subroutine sll_s_compute_e_from_rho_1d_trafo


  !> Compute rho from Gauss law for given efield
  subroutine compute_rho_from_e_1d_trafo(self, field_in, field_out)
    class(sll_t_maxwell_1d_trafo) :: self         !< Maxwell solver class
    sll_real64, intent( in    )   :: field_in(:)  !< E
    sll_real64, intent(   out )   :: field_out(:) !< rho

    call self%multiply_mass( field_in, self%work1, self%s_deg_1 )
    call multiply_gt( self, self%work1, field_out )
    field_out = - self%force_sign * field_out

  end subroutine compute_rho_from_E_1d_trafo


  !> Compute E_i from j_i integrated over the time interval using weak Ampere formulation
  subroutine compute_E_from_j_1d_trafo(self, current, component, E)
    class(sll_t_maxwell_1d_trafo)         :: self !< Maxwell solver class
    sll_real64,dimension(:),intent(in)    :: current !< Component \a component of the current integrated over time interval
    sll_int32, intent(in)                 :: component !< Component of the Efield to be computed
    sll_real64,dimension(:),intent(inout) :: E !< Updated electric field

    ! Multiply by inverse mass matrix  using the eigenvalues of the circulant inverse matrix
    if (component == 1) then        
       call self%linear_solver_mass1%solve( current, self%work1  )

    elseif (component == 2) then
       call self%linear_solver_mass0%solve( current, self%work1  )
    else
       print*, 'Component ', component, 'not implemented in compute_E_from_j_1d_trafo.'
    end if

    ! Update the electric field and scale
    E = E - self%force_sign * self%work1

  end subroutine compute_E_from_j_1d_trafo


  !> For model with adiabatic electrons
  subroutine compute_phi_from_rho_1d_trafo( self, in, phi, efield )
    class(sll_t_maxwell_1d_trafo)              :: self      !< Maxwell solver class
    sll_real64, intent(in)                     :: in(:)     !< rho
    sll_real64, intent(out)                    :: phi(:)    !< phi 
    sll_real64, intent(out)                    :: efield(:) !< E

    ! Compute phi by inverting the mass matrix
    call self%invert_mass( in, phi, self%s_deg_0 )

    ! Compute the degrees of freedom of the electric field as G phi
    call multiply_g( self, phi, efield )
    efield = -efield

  end subroutine compute_phi_from_rho_1d_trafo


  !> For model with adiabatic electrons
  subroutine compute_phi_from_j_1d_trafo( self, in, phi, efield )
    class(sll_t_maxwell_1d_trafo)              :: self      !< Maxwell solver class 
    sll_real64, intent(in)                     :: in(:)     !< Current integrated over time interval
    sll_real64, intent(out)                    :: phi(:)    !< phi 
    sll_real64, intent(out)                    :: efield(:) !< E

    ! Compute divergence of the current (G^T current) (assuming a 1v model)
    call multiply_gt( self, in, self%work1 )

    ! Compute phi by inverting the mass matrix
    call self%invert_mass( self%work1, self%work2, self%s_deg_0 )
    phi = phi + self%work2

    ! Compute the degrees of freedom of the electric field as -G phi
    call multiply_g( self, phi, efield )
    efield = -efield

  end subroutine compute_phi_from_j_1d_trafo


  !> Compute the FEM right-hand-side for a given function f and periodic splines of given degree
  !> Its components are $\int f N_i dx$ where $N_i$ is the B-spline starting at $x_i$
  subroutine sll_s_compute_rhs_trafo(self, func, degree, coefs_dofs)
    class(sll_t_maxwell_1d_trafo)       :: self  !< Maxwell solver class
    procedure(sll_i_function_1d_real64) :: func !< Function
    sll_int32,  intent( in  ) :: degree !< Specify the degree of the basis functions
    sll_real64, intent( out ) :: coefs_dofs(:)  !< Finite Element right-hand-side
    ! local variables
    sll_int32 :: i,j,k
    sll_real64 :: coef, c
    sll_real64, dimension(2,degree+1) :: xw_gauss
    sll_real64, dimension(degree+1,degree+1) :: bspl
    sll_real64 :: jacobian

    ! take enough Gauss points so that projection is exact for splines of degree deg
    ! rescale on [0,1] for compatibility with B-splines
    xw_gauss = sll_f_gauss_legendre_points_and_weights(degree+1, 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k=1,degree+1
       call sll_s_uniform_bsplines_eval_basis(degree,xw_gauss(1,k), bspl(k,:))
    end do

    if( degree == self%s_deg_0 ) then
       ! Compute coefs_dofs = int f(x)N_i(x) 
       do i = 1, self%n_cells
          coef=0.0_f64
          ! loop over support of B spline
          do j = 1, degree+1
             ! loop over Gauss points
             do k=1, degree+1
                c = self%delta_x * (xw_gauss(1,k) + real(i + j - 2 - ((i + j - 2)/self%n_cells) * self%n_cells,f64))
                jacobian= self%map%jacobian( [[c, 0.0_f64, 0._f64]])
                coef = coef + xw_gauss(2,k)*func(self%map%get_x1( [c, 0._f64, 0._f64])) * bspl(k,degree+2-j)*abs(jacobian)
             enddo
          enddo
          ! rescale by cell size
          coefs_dofs(i) = coef*self%delta_x
       enddo
    else
       ! Compute coefs_dofs = int f(x)N_i(x) 
       do i = 1, self%n_cells
          coef=0.0_f64
          ! loop over support of B spline
          do j = 1, degree+1
             ! loop over Gauss points
             do k=1, degree+1
                c = self%delta_x * (xw_gauss(1,k) + real(i + j - 2 - ((i + j - 2)/self%n_cells) * self%n_cells,f64))
                coef = coef + xw_gauss(2,k)*func(self%map%get_x1( [c, 0._f64, 0._f64])) * bspl(k,degree+2-j)* sign( 1._f64, self%map%jacobian( [c, 0._f64, 0._f64] ) )
             enddo
          enddo
          ! rescale by cell size
          coefs_dofs(i) = coef*self%delta_x
       enddo
    endif

  end subroutine sll_s_compute_rhs_trafo


  !> Compute the L2 projection of a given function f on periodic splines of given degree
  subroutine L2projection_1d_trafo(self, func, degree, coefs_dofs)
    class(sll_t_maxwell_1d_trafo) :: self !< Maxwell solver class
    procedure(sll_i_function_1d_real64) :: func !< Function
    sll_int32,  intent( in  ) :: degree !< Specify the degree of the basis functions
    sll_real64, intent( out ) :: coefs_dofs(:)  !< spline coefficients of projection

    ! Compute right-hand-side
    call self%compute_rhs_from_function( func, degree, self%work1)

    ! Multiply by inverse mass matrix 
    if (degree == self%s_deg_0) then
       call self%linear_solver_mass0%solve ( self%work1, coefs_dofs )

    elseif  (degree == self%s_deg_1) then
       call self%linear_solver_mass1%solve ( self%work1, coefs_dofs )
    else
       print*, 'degree ', degree, 'not availlable in maxwell_1d_trafo object' 
    endif
  end subroutine L2projection_1d_trafo


  !> Compute square of the L2norm 
  function L2norm_squared_1d_trafo(self, coefs_dofs, degree) result (r)
    class(sll_t_maxwell_1d_trafo) :: self !< Maxwell solver class
    sll_real64 :: coefs_dofs(:) !< Coefficient for each DoF
    sll_int32  :: degree !< Specify the degree of the basis functions
    sll_real64 :: r !< Result: squared L2 norm

    ! Multiply coefficients by mass matrix (use diagonalization FFT and mass matrix eigenvalues)
    if (degree == self%s_deg_0 ) then

       call self%multiply_mass( coefs_dofs, self%work1, self%s_deg_0 )

    elseif (degree == self%s_deg_1) then

       call self%multiply_mass( coefs_dofs, self%work1, self%s_deg_1 )

    end if
    ! Multiply by the coefficients from the left (inner product)
    r = sum(coefs_dofs*self%work1)

  end function L2norm_squared_1d_trafo


  !> Compute inner product
  function inner_product_1d_trafo(self, coefs1_dofs, coefs2_dofs, degree, degree2) result (r)
    class(sll_t_maxwell_1d_trafo) :: self !< Maxwell solver class
    sll_real64 :: coefs1_dofs(:) !< Coefficient for each DoF
    sll_real64 :: coefs2_dofs(:) !< Coefficient for each DoF
    sll_int32  :: degree !< Specify the degree of the basis functions
    sll_int32, optional  :: degree2 !< Specify the degree of the basis functions
    sll_real64 :: r !< Result: squared L2 norm

    ! Multiply coefficients by mass matrix 
    if (present(degree2)) then
       if (degree == degree2) then
          if (degree == self%s_deg_0 ) then

             call self%multiply_mass( coefs2_dofs, self%work1, self%s_deg_0 )

          elseif (degree == self%s_deg_1) then

             call self%multiply_mass( coefs2_dofs, self%work1, self%s_deg_1 )

          end if
          ! Multiply by the coefficients from the left (inner product)
          r = sum(coefs1_dofs*self%work1)
       else
          if (degree == self%s_deg_0) then

             call self%multiply_mass( coefs2_dofs, self%work1, 10 )
             ! Multiply by the coefficients from the left (inner product)
             r = sum(coefs1_dofs*self%work1)
          else
             call self%multiply_mass( coefs1_dofs, self%work1, 10 )
             ! Multiply by the coefficients from the left (inner product)
             r = sum(coefs2_dofs*self%work1)
          end if
       end if
    else
       if (degree == self%s_deg_0 ) then

          call self%multiply_mass( coefs2_dofs, self%work1, self%s_deg_0 )

       elseif (degree == self%s_deg_1) then

          call self%multiply_mass( coefs2_dofs, self%work1, self%s_deg_1 )

       end if
       ! Multiply by the coefficients from the left (inner product)
       r = sum(coefs1_dofs*self%work1)
    end if
  end function inner_product_1d_trafo


  !> Initialization
  subroutine init_1d_trafo( self, domain, n_dofs, s_deg_0, map, mass_tolerance, poisson_tolerance, solver_tolerance, force_sign )
    class(sll_t_maxwell_1d_trafo), intent(out) :: self !< solver class
    sll_real64, intent(in) :: domain(2)     !< xmin, xmax
    sll_int32, intent(in) :: n_dofs  !< number of degrees of freedom (here number of cells and grid points)
    sll_int32, intent(in) :: s_deg_0 !< highest spline degree
    type(sll_t_mapping_3d), target, intent(inout) :: map !< coordinate transformation
    sll_real64, intent(in), optional :: mass_tolerance !< tolerance for mass solver
    sll_real64, intent(in), optional :: poisson_tolerance !< tolerance for Poisson solver
    sll_real64, intent(in), optional :: solver_tolerance !< tolerance for Schur complement solver
    sll_real64, intent(in), optional :: force_sign !< sign of particle force
    ! local variables
    sll_real64, allocatable :: nullspace(:,:)
    sll_int32 :: ierr

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

    if (present( force_sign) ) then
       self%force_sign = force_sign
    end if

    self%n_dofs0 = n_dofs
    self%n_dofs1 = n_dofs
    self%n_cells = n_dofs
    self%Lx =  map%params(1)
    self%delta_x = 1._f64 / real(n_dofs,f64)
    self%s_deg_0 = s_deg_0
    self%s_deg_1 = s_deg_0 - 1
    self%map=>map

    SLL_ALLOCATE(self%work1(n_dofs),ierr)
    SLL_ALLOCATE(self%work2(n_dofs),ierr)

    ! Sparse matrices
    call sll_s_spline_fem_mass1d_full ( n_dofs, self%s_deg_0, self%mass0, profile_0 )
    call sll_s_spline_fem_mass1d_full ( n_dofs, self%s_deg_1, self%mass1, profile_1 )

    call self%linear_solver_mass0%create( self%mass0)
    self%linear_solver_mass0%atol = self%mass_solver_tolerance

    call self%linear_solver_mass1%create( self%mass1)
    self%linear_solver_mass1%atol = self%mass_solver_tolerance /self%Lx
    !self%linear_solver_mass0%verbose = .true.
    !self%linear_solver_mass1%verbose = .true.

    ! Penalized Poisson operator
    allocate(nullspace(1,1:self%n_cells))
    nullspace(1,:) = 1.0_f64
    call self%poisson_matrix%create( self%mass1, self%n_cells, self%delta_x )
    call self%poisson_operator%create( self%poisson_matrix, vecs=nullspace, &
         n_dim_nullspace=1 )

    call self%poisson_solver%create( self%poisson_operator )
    self%poisson_solver%null_space = .true.
    self%poisson_solver%atol = self%poisson_solver_tolerance
    !self%poisson_solver%verbose = .true.

    call sll_s_spline_fem_mixedmass1d_full( n_dofs, self%s_deg_0, self%s_deg_1, self%mixed_mass, profile_m0 )

    ! Only for Schur complement eb solver
    call self%linear_op_schur_eb%create( self%mass0, self%mass1, self%n_cells, self%delta_x )
    call self%linear_solver_schur_eb%create( self%linear_op_schur_eb )
    self%linear_solver_schur_eb%atol = self%solver_tolerance
    self%linear_solver_schur_eb%rtol = self%solver_tolerance
    !self%linear_solver_schur_eb%verbose = .true.

  contains
    function profile_0(x)
      sll_real64 :: profile_0
      sll_real64, intent(in) :: x

      profile_0 = map%jacobian( [x, 0._f64, 0._f64] )

    end function profile_0

    function profile_1(x)
      sll_real64 :: profile_1
      sll_real64, intent(in) :: x

      profile_1 = 1._f64/map%jacobian( [x, 0._f64, 0._f64] )

    end function profile_1

    function profile_m0(x)
      sll_real64 :: profile_m0
      sll_real64, intent(in) :: x

      profile_m0 = 1._f64
    end function profile_m0

  end subroutine init_1d_trafo


  !> Initialization from nml file
  subroutine init_from_file_1d_trafo( self, domain, n_dofs, s_deg_0, map, nml_file, force_sign )
    class(sll_t_maxwell_1d_trafo), intent(out) :: self !< solver class
    sll_real64, intent(in) :: domain(2)     !< xmin, xmax
    sll_int32, intent(in) :: n_dofs  !< number of degrees of freedom (here number of cells and grid points)
    sll_int32, intent(in) :: s_deg_0 !< highest spline degree
    type(sll_t_mapping_3d), target, intent(inout) :: map !< coordinate transformation
    character(len=*), intent(in) :: nml_file !< nml-file
    sll_real64, intent(in), optional :: force_sign !< sign of particle force
    ! local variables
    sll_int32 :: input_file
    sll_int32 :: io_stat, io_stat1, rank, file_id
    sll_real64 :: mass_tolerance
    sll_real64 :: poisson_tolerance
    sll_real64 :: maxwell_tolerance

    namelist /maxwell_solver/ mass_tolerance, poisson_tolerance

    namelist /time_solver/ maxwell_tolerance 

    rank = sll_f_get_collective_rank(sll_v_world_collective)

    if (present( force_sign) ) then
       self%force_sign = force_sign
    end if

    ! Read in solver tolerance
    open(newunit = input_file, file=nml_file, status='old',IOStat=io_stat)
    if (io_stat /= 0) then
       if (rank == 0 ) then
          print*, 'sll_m_maxwell_1d_trafo: Input file does not exist. Set default tolerance.'
          open(newunit=file_id, file=trim(nml_file)//'_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
          write(file_id, *) 'mass solver tolerance:', 1d-12
          write(file_id, *) 'poisson solver tolerance:', 1d-12
          close(file_id)
       end if
       call self%init( domain, n_dofs, s_deg_0, map )
    else       
       read(input_file, maxwell_solver,IOStat=io_stat)
       read(input_file, time_solver,IOStat=io_stat1)
       if (io_stat /= 0 .and. io_stat1 /= 0) then
          if (rank == 0 ) then
             print*, 'sll_m_maxwell_1d_trafo: Input parameter does not exist. Set default tolerance.'
             open(newunit=file_id, file=trim(nml_file)//'_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'mass solver tolerance:', 1d-12
             write(file_id, *) 'poisson solver tolerance:', 1d-12
             close(file_id)
          end if
          call self%init( domain, n_dofs, s_deg_0, map )
       else if (io_stat /= 0 .and. io_stat1 == 0) then
          call self%init( domain, n_dofs, s_deg_0, map, solver_tolerance=maxwell_tolerance )
       else if (io_stat == 0 .and. io_stat1 /= 0) then
          if (rank == 0 ) then
             open(newunit=file_id, file=trim(nml_file)//'_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'mass solver tolerance:', mass_tolerance
             write(file_id, *) 'poisson solver tolerance:', poisson_tolerance
             close(file_id)
          end if
          call self%init( domain, n_dofs, s_deg_0, map, mass_tolerance, poisson_tolerance )
       else
          if (rank == 0 ) then
             open(newunit=file_id, file=trim(nml_file)//'_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'mass solver tolerance:', mass_tolerance
             write(file_id, *) 'poisson solver tolerance:', poisson_tolerance
             close(file_id)
          end if
          call self%init( domain, n_dofs, s_deg_0, map, mass_tolerance, poisson_tolerance, solver_tolerance=maxwell_tolerance )
       end if
       close(input_file)
    end if


  end subroutine init_from_file_1d_trafo


  !> Finalization
  subroutine free_1d_trafo(self)
    class(sll_t_maxwell_1d_trafo) :: self !< Maxwell solver class

    deallocate(self%work1)
    deallocate(self%work2)
    self%map=> null()
    call self%linear_solver_mass0%free()
    call self%linear_solver_mass1%free()
    call self%poisson_solver%free()
    call self%poisson_operator%free()
    call self%poisson_matrix%free()
    call self%linear_solver_schur_eb%free()
    call self%linear_op_schur_eb%free()

  end subroutine free_1d_trafo


  !> Multiply by dicrete gradient matrix
  subroutine multiply_g( self,  in, out)
    class(sll_t_maxwell_1d_trafo), intent(in) :: self !< Maxwell solver class
    sll_real64, intent(in)  :: in(:) !< field_in  
    sll_real64, intent(out) :: out(:) !< G*field_in

    call sll_s_multiply_g_1d( self%n_cells, self%delta_x, in, out )

  end subroutine multiply_g


  !> Multiply by transpose of dicrete gradient matrix
  subroutine multiply_gt( self,  in, out)
    class(sll_t_maxwell_1d_trafo), intent(in) :: self !< Maxwell solver class
    sll_real64, intent(in)  :: in(:) !< field_in  
    sll_real64, intent(out) :: out(:) !< G^T*field_in

    call sll_s_multiply_gt_1d( self%n_cells, self%delta_x, in, out )

  end subroutine multiply_gt


  !> Multiply by the mass matrix 
  subroutine multiply_mass_1d_trafo( self,  in, out, degree)
    class(sll_t_maxwell_1d_trafo), intent(inout) :: self !< Maxwell solver class
    sll_real64, intent(in)  :: in(:) !< Coefficient for each DoF
    sll_real64, intent(out) :: out(:) !< Coefficient for each DoF
    sll_int32, intent(in)  :: degree !< Specify the degree of the basis functions

    if (degree == self%s_deg_0 ) then 
       call self%mass0%dot( in, out )
    elseif (degree == self%s_deg_1) then
       call self%mass1%dot( in, out )
    elseif(degree == 10) then
       call self%mixed_mass%dot( in, out )
    else
       print*, 'multiply mass for other form not yet implemented'
       stop
    end if

  end subroutine multiply_mass_1d_trafo


  !> Multiply by the inverse mass matrix 
  subroutine invert_mass_1d_trafo( self,  in, out, degree)
    class(sll_t_maxwell_1d_trafo), intent(inout) :: self !< Maxwell solver class
    sll_real64, intent(in)  :: in(:) !< Coefficient for each DoF
    sll_real64, intent(out) :: out(:) !< Coefficient for each DoF
    sll_int32, intent(in)  :: degree !< Specify the degree of the basis functions

    if (degree == self%s_deg_0 ) then 
       call self%linear_solver_mass0%solve( in, out )
    elseif (degree == self%s_deg_1) then
       call self%linear_solver_mass1%solve( in, out )
    else
       print*, 'Invert mass for other form not yet implemented'
       stop
    end if

  end subroutine invert_mass_1d_trafo


end module sll_m_maxwell_1d_trafo
