!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations in 3D
!> The linear systems are solved based on FFT diagnoalization
!> @details
!> 
!> @author
!> Katharina Kormann


! TODO: Write FFT-based mass solver: There is such a solver already defined as linear_solver_mass1 in particle_methods. Reuse? Can we put the parallelization in this solver?
! Remove all parts that belong the PLF
! Add also solver for combined e and b (first step for AVF-based algorithm)


module sll_m_maxwell_3d_fem_fft
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"


  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_linear_operator_kron, only : &
       sll_t_linear_operator_kron

  use sll_m_linear_operator_maxwell_eb_schur, only : &
       sll_t_linear_operator_maxwell_eb_schur

  use sll_m_linear_solver_cg, only : &
       sll_t_linear_solver_cg

  use sll_m_linear_solver_mgmres, only : &
       sll_t_linear_solver_mgmres

  use sll_m_linear_solver_spline_mass_fft, only : &
       sll_t_linear_solver_spline_mass_fft

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_maxwell_3d_base, only: &
       sll_i_function_3d_real64, &
       sll_c_maxwell_3d_base

  use sll_m_poisson_3d_fem_fft, only : &
       sll_t_poisson_3d_fem_fft

  use sll_m_profile_functions, only: &
       sll_t_profile_functions

  use sll_m_spline_fem_utilities, only : &
       sll_s_spline_fem_mass_line, &
       sll_s_spline_fem_mixedmass_line, &
       sll_s_spline_fem_multiply_mass, &
       sll_s_spline_fem_multiply_massmixed, &
       sll_s_spline_fem_compute_mass_eig, &
       sll_s_multiply_g, &
       sll_s_multiply_gt, &
       sll_s_multiply_c, &
       sll_s_multiply_ct

  use sll_m_spline_fem_utilities_3d, only: &
       sll_s_spline_fem_mass3d

  use sll_m_spline_fem_utilities_sparse, only : &
       sll_s_spline_fem_mass1d

  implicit none

  public :: &
       sll_t_maxwell_3d_fem_fft

  private

  type, extends(sll_c_maxwell_3d_base) :: sll_t_maxwell_3d_fem_fft

     sll_real64, allocatable :: work(:)  !< scratch data
     sll_real64, allocatable :: work2(:)  !< scratch data
     sll_real64, allocatable :: work3d(:,:,:)  !< scratch data
     sll_real64, allocatable :: work_d1(:) !< scratch data
     sll_real64, allocatable :: work_d2_in(:) !< scratch data
     sll_real64, allocatable :: work_d2_out(:) !< scratch data
     sll_real64, allocatable :: work_d3_in(:) !< scratch data
     sll_real64, allocatable :: work_d3_out(:) !< scratch data

     sll_real64, allocatable :: mass_line0_1(:) !< massline
     sll_real64, allocatable :: mass_line0_2(:) !< massline
     sll_real64, allocatable :: mass_line0_3(:) !< massline
     sll_real64, allocatable :: mass_line1_1(:) !< massline
     sll_real64, allocatable :: mass_line1_2(:) !< massline
     sll_real64, allocatable :: mass_line1_3(:) !< massline
     sll_real64, allocatable :: mass_line_mixed_1(:) !< mixed massline
     sll_real64, allocatable :: mass_line_mixed_2(:) !< mixed massline
     sll_real64, allocatable :: mass_line_mixed_3(:) !< mixed massline
     type(sll_t_linear_solver_spline_mass_fft) :: inverse_mass_0 !< Fourier solver for 0-form mass matrix
     type(sll_t_linear_solver_spline_mass_fft) :: inverse_mass_1(3) !< Fourier solver for 1-form mass matrix
     type(sll_t_linear_solver_spline_mass_fft) :: inverse_mass_2(3) !< Fourier solver for 2-form mass matrix
     type(sll_t_poisson_3d_fem_fft ) :: poisson_fft !< Fourier solver for Poisson matrix

     logical :: adiabatic_electrons = .false. !< Set true if solver with adiabatic electrions

     type(sll_t_matrix_csr)  :: mass0 !< 0-form mass matrix
     type(sll_t_matrix_csr)  :: mass1d(3,3) !< 1D mass matrices 
     type(sll_t_linear_operator_kron)  :: mass1(3) !< Tensorproduct 1-form mass matrix
     type(sll_t_linear_operator_kron)  :: mass2(3) !< Tensorproduct 2-form mass matrix
     type(sll_t_linear_solver_cg) :: mass0_solver     !< mass matrix solver
     type( sll_t_linear_operator_maxwell_eb_schur ) :: linear_op_schur_eb !< Schur complement operator for advect_eb
    
   contains
     procedure :: &
          compute_E_from_B => sll_s_compute_e_from_b_3d_fem_fft !< Solve E and B part of Amperes law with B constant in time
     procedure :: &
          compute_B_from_E => sll_s_compute_b_from_e_3d_fem_fft !< Solve Faraday equation with E constant in time
     procedure :: &
          compute_curl_part => sll_s_compute_curl_part_3d_fem_fft !< Solve source-free Maxwell's equations
     procedure :: &
          compute_E_from_rho => sll_s_compute_e_from_rho_3d_fem_fft !< Solve E from rho using Poisson
     procedure :: &
          compute_rho_from_e => sll_s_compute_rho_from_e_3d_fem_fft !< Compute rho from E 
     procedure :: &
          compute_E_from_j => sll_s_compute_E_from_j_3d_fem_fft !< Compute E from the current j
     procedure :: &
          compute_phi_from_rho => sll_s_compute_phi_from_rho_3d_fem_fft !< Compute phi from rho (by solving the quasi-neutrality equation)
     procedure :: &
          compute_phi_from_j => sll_s_compute_phi_from_j_3d_fem_fft !< Compute phi from j (dynamic of quasi-neutrality equation for adiabatic electrons)
     procedure :: &
          compute_rhs_from_function => sll_s_compute_rhs_fem_fft !< Compute integral over given function tested by the basis
     procedure :: &
          L2projection => L2projection_3d_fem_fft  !< Compute L_2 projection of a given function
     procedure :: &
          L2norm_squared => L2norm_squared_3d_fem_fft  !< Compute the square of the L2 norm of a given vector
     procedure :: &
          inner_product => inner_product_3d_fem_fft !< Inner product of two dof-vectors with mass matrix
     procedure :: &
          init => init_3d_fem_fft !< Initialize the Maxwell class
     procedure :: &
          free => free_3d_fem_fft !< Free Maxwell class
     procedure :: &
          multiply_g !< Multiplication with gradient matrix 
     procedure :: &
          multiply_gt !< Multiplication with transposed gradient matrix  
     procedure :: &
          multiply_c !< Multiplication with curl matrix
     procedure :: &
          multiply_ct !< Multiplication with transposed curl matrix
     procedure :: &
          multiply_mass => multiply_mass_3d_fem_fft !< Product with the mass matrix
     procedure :: &
          multiply_mass_inverse => multiply_mass_inverse_all !< Invert mass matrix

  end type sll_t_maxwell_3d_fem_fft

contains


  !> compute Ey from Bz using weak Ampere formulation 
  subroutine sll_s_compute_e_from_b_3d_fem_fft(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_3d_fem_fft) :: self !< Maxwell solver class
    sll_real64, intent(in)     :: delta_t   !< Time step
    sll_real64, intent(in)     :: field_in(:)  !< B
    sll_real64, intent(inout)  :: field_out(:)  !< E
    ! local variables
    sll_real64 :: coef
    sll_int32  :: comp, istart, iend

    call multiply_mass_2form( self, field_in, self%work )

    call multiply_ct(self, self%work, self%work2)

    do comp=1,3
       istart = 1+(comp-1)*self%n_total
       iend =  comp*self%n_total
       call self%inverse_mass_1(comp)%solve( self%work2(istart:iend), self%work(istart:iend) ) 
    end do
    ! Update b from self value
    coef = delta_t
    field_out = field_out + coef*self%work
  end subroutine sll_s_compute_e_from_b_3d_fem_fft


  !> Compute Bz from Ey using strong 1D Faraday equation for spline coefficients
  !> $B_z^{new}(x_j) = B_z^{old}(x_j) - \frac{\Delta t}{\Delta x} (E_y(x_j) - E_y(x_{j-1})  $
  subroutine sll_s_compute_b_from_e_3d_fem_fft(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_3d_fem_fft) :: self         !< Maxwell solver class
    sll_real64, intent( in    )     :: delta_t      !< time step
    sll_real64, intent( in    )     :: field_in(:)  !< E
    sll_real64, intent( inout )     :: field_out(:) !< B

    call multiply_c(self, field_in, self%work)

    field_out = field_out - delta_t * self%work

  end subroutine sll_s_compute_b_from_e_3d_fem_fft


  !> Solve curl part of Maxwell's equations
  subroutine sll_s_compute_curl_part_3d_fem_fft( self, delta_t, efield, bfield, betar )
    class(sll_t_maxwell_3d_fem_fft) :: self       !< Maxwell solver class
    sll_real64, intent( in    )     :: delta_t    !< Time step
    sll_real64, intent( inout )     :: efield(:)  !< E
    sll_real64, intent( inout )     :: bfield(:)  !< B
    sll_real64, optional            :: betar      !< 1/beta
    !local variables
    sll_real64 :: factor

    if( present(betar) ) then
       factor = betar
    else
       factor = 1._f64
    end if

    ! Compute C^T M2 b
    call multiply_mass_2form( self, bfield, self%work )
    call self%multiply_ct( self%work, self%work2 ) 

    self%linear_op_schur_eb%factor = -delta_t**2*factor*0.25_f64
    call self%linear_op_schur_eb%dot( efield, self%work )
    self%work = self%work + delta_t*factor*self%work2

    ! Save efield dofs from previous time step for B field update
    self%work2 = efield

    ! Invert Schur complement matrix
    self%linear_op_schur_eb%factor = delta_t**2*factor*0.25_f64
    call self%linear_op_schur_eb%dot_inverse( self%work, efield )

    ! Update B field
    self%work2 = self%work2 + efield
    call self%compute_b_from_e( delta_t*0.5_f64, self%work2, bfield )

  end subroutine sll_s_compute_curl_part_3d_fem_fft


  !> Compute E_i from rho_i integrated over the time interval using weak Poisson's equation
  subroutine sll_s_compute_e_from_rho_3d_fem_fft(self, field_in, field_out )  
    class(sll_t_maxwell_3d_fem_fft) :: self         !< Maxwell solver class
    sll_real64, intent( in    )     :: field_in(:)  !< rho
    sll_real64, intent(   out )     :: field_out(:) !< E


    call self%poisson_fft%compute_e_from_rho( field_in, field_out )

  end subroutine sll_s_compute_e_from_rho_3d_fem_fft


  !> compute rho from e using weak Gauss law ( rho = G^T M_1 e ) 
  subroutine sll_s_compute_rho_from_e_3d_fem_fft(self, field_in, field_out ) 
    class(sll_t_maxwell_3d_fem_fft) :: self         !< Maxwell solver class
    sll_real64, intent( in    )     :: field_in(:)  !< E
    sll_real64, intent(   out )     :: field_out(:) !< rho

    call multiply_mass_1form( self, field_in, self%work )

    call multiply_gt( self, self%work,  field_out )
    field_out = - field_out
  end subroutine sll_s_compute_rho_from_e_3d_fem_fft


  !> Compute E_i from j_i integrated over the time interval using weak Ampere formulation
  subroutine sll_s_compute_E_from_j_3d_fem_fft(self, current, E, component)
    class(sll_t_maxwell_3d_fem_fft)       :: self       !< Maxwell solver class
    sll_real64, intent( in    )           :: current(:) !< Current integrated over time interval
    sll_real64, intent( inout )           :: E(:)       !< Updated electric field
    sll_int32,  intent( in    ), optional :: component  !< component of the Efield to be computed

    if(present(component) ) then
       call self%inverse_mass_1(component)%solve( current, self%work(1:self%n_total) )
       E = E - self%work(1:self%n_total)
    else
       call self%inverse_mass_1(1)%solve( current(1:self%n_total), self%work(1:self%n_total) )
       call self%inverse_mass_1(2)%solve( current(self%n_total+1:2*self%n_total), self%work(self%n_total+1:2*self%n_total) )
       call self%inverse_mass_1(3)%solve( current(2*self%n_total+1:3*self%n_total), self%work(2*self%n_total+1:3*self%n_total) )
       E = E - self%work
    end if

  end subroutine sll_s_compute_E_from_j_3d_fem_fft


  !> Compute phi from rho_i integrated over the time interval
  subroutine sll_s_compute_phi_from_rho_3d_fem_fft( self, field_in, field_out, efield_dofs )  
    class(sll_t_maxwell_3d_fem_fft) :: self         !< Maxwell solver class
    sll_real64, intent( in    )           :: field_in(:)  !< rho
    sll_real64, intent( inout )           :: field_out(:) !< phi
    sll_real64, intent(   out )           :: efield_dofs(:) !< E

    call self%mass0_solver%solve( field_in, field_out )
    call self%multiply_g( field_out, efield_dofs )
    efield_dofs = -efield_dofs

  end subroutine sll_s_compute_phi_from_rho_3d_fem_fft


  !> Compute phi from j_i integrated over the time interval, delta_t is already included 
  subroutine sll_s_compute_phi_from_j_3d_fem_fft( self, field_in, field_out, efield_dofs )
    class(sll_t_maxwell_3d_fem_fft) :: self       !< Maxwell solver class
    sll_real64, intent( in    )           :: field_in(:) !< Current integrated over time interval
    sll_real64, intent( inout )           :: field_out(:) !< phi
    sll_real64, intent(   out )           :: efield_dofs(:) !< E

    call self%multiply_gt( field_in, self%work(1:self%n_total) ) 
    call self%mass0_solver%solve( self%work(1:self%n_total), self%work2(1:self%n_total) )
    field_out = field_out + self%work2(1:self%n_total)

    call self%multiply_g( field_out, efield_dofs )
    efield_dofs = -efield_dofs

  end subroutine sll_s_compute_phi_from_j_3d_fem_fft


  !> Compute the FEM right-hand-side for a given function f and periodic splines of given degree
  !> Its components are $\int f N_i dx$ where $N_i$ is the B-spline starting at $x_i$ 
  subroutine sll_s_compute_rhs_fem_fft(self, form,  component, coefs_dofs, func1, func2, func3 )
    class(sll_t_maxwell_3d_fem_fft)                :: self          !< Maxwell solver class
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

  end subroutine sll_s_compute_rhs_fem_fft


  !> Compute the L2 projection of a given function f on periodic splines of given degree
  subroutine L2projection_3d_fem_fft(self, form,  component, coefs_dofs, func1, func2, func3 )
    class(sll_t_maxwell_3d_fem_fft)                :: self          !< Maxwell solver class
    sll_int32,  intent( in    )                    :: form          !< Specify if the function is a0,1,2 or 3-form
    sll_int32,  intent( in    )                    :: component     !< Specify the component of the function
    sll_real64, intent(   out )                    :: coefs_dofs(:) !< Finite Element right-hand-side
    procedure(sll_i_function_3d_real64)            :: func1         !< Function first component
    procedure(sll_i_function_3d_real64), optional  :: func2         !< Function second component
    procedure(sll_i_function_3d_real64), optional  :: func3         !< Function third component

    ! Compute right-hand-side
    call sll_s_compute_rhs_fem_fft( self, form,  component, self%work(1:self%n_total), func1 )

    select case( form )
    case( 0 )
       call self%inverse_mass_0%solve( self%work(1:self%n_total), coefs_dofs )
    case( 1 )
       call self%inverse_mass_1(component)%solve( self%work(1:self%n_total), coefs_dofs )
    case( 2 )
       call self%inverse_mass_2(component)%solve( self%work(1:self%n_total), coefs_dofs )
    case  default
       print*, 'L2projection for', form, '-form not implemented.'
    end select

  end subroutine L2projection_3d_fem_fft


  !> Compute square of the L2norm 
  function L2norm_squared_3d_fem_fft(self, coefs, form, component ) result ( r )
    class(sll_t_maxwell_3d_fem_fft) :: self      !< Maxwell solver class
    sll_real64                      :: coefs(:)  !< Coefficient for each DoF
    sll_int32                       :: form      !< Specify 0,1,2 or 3-form
    sll_int32                       :: component !< Specify the component of the form
    sll_real64                      :: r         !< Result: squared L2 norm

    r = inner_product_3d_fem_fft(self, coefs, coefs, form, component)

  end function L2norm_squared_3d_fem_fft


  !> Compute inner product
  function inner_product_3d_fem_fft(self, coefs1, coefs2, form, component) result (r)
    class(sll_t_maxwell_3d_fem_fft) :: self !< Maxwell solver class
    sll_real64 :: coefs1(:) !< Coefficient for each DoF
    sll_real64 :: coefs2(:) !< Coefficient for each DoF
    sll_int32  :: form !< Specify 0,1,2 or 3 form
    sll_int32, optional  :: component !< Specify the component
    sll_real64 :: r !< Result: squared L2 norm

    if ( form == 0 ) then
       call multiply_mass_3dkron( self, self%mass_line0_1, &
            self%mass_line0_2, self%mass_line0_3, &
            coefs2, self%work(1:self%n_total) )
    elseif (form == 1 ) then
       select case(component)
       case (1)
          call multiply_mass_3dkron( self, self%mass_line1_1, &
               self%mass_line0_2, self%mass_line0_3, &
               coefs2, self%work(1:self%n_total) )
       case(2)
          call multiply_mass_3dkron( self, self%mass_line0_1, &
               self%mass_line1_2, self%mass_line0_3, &
               coefs2, self%work(1:self%n_total) )
       case(3)
          call multiply_mass_3dkron( self, self%mass_line0_1, &
               self%mass_line0_2, self%mass_line1_3, &
               coefs2, self%work(1:self%n_total) )
       case default
          print*, 'wrong component.'
       end select
    elseif( form == 2) then
       select case(component)
       case (1)
          call multiply_mass_3dkron( self, self%mass_line0_1, &
               self%mass_line1_2, self%mass_line1_3, &
               coefs2, self%work(1:self%n_total) )
       case(2)
          call multiply_mass_3dkron( self, self%mass_line1_1, &
               self%mass_line0_2, self%mass_line1_3, &
               coefs2, self%work(1:self%n_total) )
       case(3)
          call multiply_mass_3dkron( self, self%mass_line1_1, &
               self%mass_line1_2, self%mass_line0_3, &
               coefs2, self%work(1:self%n_total) )
       case default
          print*, 'wrong component.'
       end select
    elseif( form == 3) then
       call multiply_mass_3dkron( self, self%mass_line1_1, &
            self%mass_line1_2, self%mass_line1_3, &
            coefs2, self%work(1:self%n_total) )
    else
       print*, 'Wrong form.'
    end if

    r = sum(coefs1*self%work(1:self%n_total))
    !r = r*self%volume


  end function inner_product_3d_fem_fft


  !> Initialization
  subroutine init_3d_fem_fft( self, domain, n_dofs, s_deg_0, adiabatic_electrons, profile  )
    class(sll_t_maxwell_3d_fem_fft), intent(out) :: self !< Maxwell solver class
    sll_real64, intent(in) :: domain(3,2)     !< xmin, xmax
    sll_int32, intent(in) :: n_dofs(3)  !< number of degrees of freedom (here number of cells and grid points)
    sll_int32, intent(in) :: s_deg_0(3) !< highest spline degree
    logical, intent(in), optional :: adiabatic_electrons !< flag if adiabatic electrons are used
    type(sll_t_profile_functions), intent(in), optional :: profile !< temperature and density profiles
    ! local variables
    sll_int32 :: j
    sll_real64, allocatable :: eig_values_mass_0_1(:)
    sll_real64, allocatable :: eig_values_mass_0_2(:)
    sll_real64, allocatable :: eig_values_mass_0_3(:)
    sll_real64, allocatable :: eig_values_mass_1_1(:)
    sll_real64, allocatable :: eig_values_mass_1_2(:)
    sll_real64, allocatable :: eig_values_mass_1_3(:)
    sll_real64, allocatable :: inv_eig_values_mass_0_1(:)
    sll_real64, allocatable :: inv_eig_values_mass_0_2(:)
    sll_real64, allocatable :: inv_eig_values_mass_0_3(:)
    sll_real64, allocatable :: inv_eig_values_mass_1_1(:)
    sll_real64, allocatable :: inv_eig_values_mass_1_2(:)
    sll_real64, allocatable :: inv_eig_values_mass_1_3(:)

    self%n_cells = n_dofs
    self%n_dofs = n_dofs
    self%n_total = product(n_dofs)
    self%n_total0 = self%n_total
    self%n_total1 = self%n_total

    self%Lx = domain(:,2) - domain(:,1)
    self%delta_x = self%Lx / real(n_dofs,f64)
    self%s_deg_0 = s_deg_0
    self%s_deg_1 = s_deg_0 - 1

    self%volume = product(self%delta_x)

    if( present( adiabatic_electrons ) ) then
       self%adiabatic_electrons = adiabatic_electrons
    end if

    if( present( profile ) ) then
       self%profile = profile
    end if

    ! Allocate scratch data
    allocate( self%work3d(n_dofs(1), n_dofs(2), n_dofs(3)) )
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
    allocate( self%mass_line0_1(s_deg_0(1)+1) )
    allocate( self%mass_line1_1(s_deg_0(1)) )
    allocate( self%mass_line0_2(s_deg_0(2)+1) )
    allocate( self%mass_line1_2(s_deg_0(2)) )
    allocate( self%mass_line0_3(s_deg_0(3)+1) )
    allocate( self%mass_line1_3(s_deg_0(3)) )
    allocate( self%mass_line_mixed_1(2*s_deg_0(1)) )
    allocate( self%mass_line_mixed_2(2*s_deg_0(2)) )
    allocate( self%mass_line_mixed_3(2*s_deg_0(3)) )


    call sll_s_spline_fem_mass_line ( self%s_deg_0(1), self%mass_line0_1 )
    call sll_s_spline_fem_mass_line ( self%s_deg_1(1), self%mass_line1_1 )
    call sll_s_spline_fem_mixedmass_line ( self%s_deg_0(1), self%mass_line_mixed_1 )

    call sll_s_spline_fem_mass_line ( self%s_deg_0(2), self%mass_line0_2 )
    call sll_s_spline_fem_mass_line ( self%s_deg_1(2), self%mass_line1_2 )
    call sll_s_spline_fem_mixedmass_line ( self%s_deg_0(2), self%mass_line_mixed_2 )

    call sll_s_spline_fem_mass_line ( self%s_deg_0(3), self%mass_line0_3 )
    call sll_s_spline_fem_mass_line ( self%s_deg_1(3), self%mass_line1_3 )
    call sll_s_spline_fem_mixedmass_line ( self%s_deg_0(3), self%mass_line_mixed_3 )

    self%mass_line0_1 = self%delta_x(1) * self%mass_line0_1
    self%mass_line1_1 = self%delta_x(1) *  self%mass_line1_1 
    self%mass_line_mixed_1 = self%delta_x(1) * self%mass_line_mixed_1

    self%mass_line0_2 = self%delta_x(2) * self%mass_line0_2
    self%mass_line1_2 = self%delta_x(2) * self%mass_line1_2
    self%mass_line_mixed_2 = self%delta_x(2) * self%mass_line_mixed_2

    self%mass_line0_3 = self%delta_x(3) * self%mass_line0_3
    self%mass_line1_3 = self%delta_x(3) * self%mass_line1_3
    self%mass_line_mixed_3 = self%delta_x(3) * self%mass_line_mixed_3

    allocate( eig_values_mass_0_1( n_dofs(1) ) )
    allocate( eig_values_mass_0_2( n_dofs(2) ) )
    allocate( eig_values_mass_0_3( n_dofs(3) ) )
    allocate( eig_values_mass_1_1( n_dofs(1) ) )
    allocate( eig_values_mass_1_2( n_dofs(2) ) )
    allocate( eig_values_mass_1_3( n_dofs(3) ) )
    allocate( inv_eig_values_mass_0_1( n_dofs(1) ) )
    allocate( inv_eig_values_mass_0_2( n_dofs(2) ) )
    allocate( inv_eig_values_mass_0_3( n_dofs(3) ) )
    allocate( inv_eig_values_mass_1_1( n_dofs(1) ) )
    allocate( inv_eig_values_mass_1_2( n_dofs(2) ) )
    allocate( inv_eig_values_mass_1_3( n_dofs(3) ) )

    call sll_s_spline_fem_compute_mass_eig( n_dofs(1), self%s_deg_0(1), self%mass_line0_1, &
         eig_values_mass_0_1 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(2), self%s_deg_0(2), self%mass_line0_2, &
         eig_values_mass_0_2 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(3), self%s_deg_0(3), self%mass_line0_3, &
         eig_values_mass_0_3 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(1), self%s_deg_1(1), self%mass_line1_1, &
         eig_values_mass_1_1 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(2), self%s_deg_1(2), self%mass_line1_2, &
         eig_values_mass_1_2 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(3), self%s_deg_1(3), self%mass_line1_3, &
         eig_values_mass_1_3 )

    inv_eig_values_mass_0_1(1) = 1._f64/eig_values_mass_0_1(1)
    inv_eig_values_mass_0_1(n_dofs(1)/2+1) = 1._f64/eig_values_mass_0_1(n_dofs(1)/2+1)
    inv_eig_values_mass_1_1(1) = 1._f64/eig_values_mass_1_1(1)
    inv_eig_values_mass_1_1(n_dofs(1)/2+1) = 1._f64/eig_values_mass_1_1(n_dofs(1)/2+1)
    do j = 2, n_dofs(1)/2
       inv_eig_values_mass_0_1(j) = 1._f64/eig_values_mass_0_1(j)
       inv_eig_values_mass_0_1(n_dofs(1)+2-j) = 1._f64/eig_values_mass_0_1(n_dofs(1)+2-j)
       inv_eig_values_mass_1_1(j) = 1._f64/eig_values_mass_1_1(j)
       inv_eig_values_mass_1_1(n_dofs(1)+2-j) = 1._f64/eig_values_mass_1_1(n_dofs(1)+2-j)
    end do

    inv_eig_values_mass_0_2(1) = 1._f64/eig_values_mass_0_2(1)
    inv_eig_values_mass_0_2(n_dofs(2)/2+1) = 1._f64/eig_values_mass_0_2(n_dofs(2)/2+1)
    inv_eig_values_mass_1_2(1) = 1._f64/eig_values_mass_1_2(1)
    inv_eig_values_mass_1_2(n_dofs(2)/2+1) = 1._f64/eig_values_mass_1_2(n_dofs(2)/2+1)
    do j = 2, n_dofs(2)/2
       inv_eig_values_mass_0_2(j) = 1._f64/eig_values_mass_0_2(j)
       inv_eig_values_mass_0_2(n_dofs(2)+2-j) = 1._f64/eig_values_mass_0_2(n_dofs(2)+2-j)
       inv_eig_values_mass_1_2(j) = 1._f64/eig_values_mass_1_2(j)
       inv_eig_values_mass_1_2(n_dofs(2)+2-j) = 1._f64/eig_values_mass_1_2(n_dofs(2)+2-j)
    end do

    inv_eig_values_mass_0_3(1) = 1._f64/eig_values_mass_0_3(1)
    inv_eig_values_mass_0_3(n_dofs(3)/2+1) = 1._f64/eig_values_mass_0_3(n_dofs(3)/2+1)
    inv_eig_values_mass_1_3(1) = 1._f64/eig_values_mass_1_3(1)
    inv_eig_values_mass_1_3(n_dofs(3)/2+1) = 1._f64/eig_values_mass_1_3(n_dofs(3)/2+1)
    do j = 2, n_dofs(3)/2
       inv_eig_values_mass_0_3(j) = 1._f64/eig_values_mass_0_3(j)
       inv_eig_values_mass_0_3(n_dofs(3)+2-j) = 1._f64/eig_values_mass_0_3(n_dofs(3)+2-j)
       inv_eig_values_mass_1_3(j) = 1._f64/eig_values_mass_1_3(j)
       inv_eig_values_mass_1_3(n_dofs(3)+2-j) = 1._f64/eig_values_mass_1_3(n_dofs(3)+2-j)
    end do


    call self%inverse_mass_0%create( n_dofs, inv_eig_values_mass_0_1, inv_eig_values_mass_0_2, inv_eig_values_mass_0_3 )

    call self%inverse_mass_1(1)%create( n_dofs, inv_eig_values_mass_1_1, inv_eig_values_mass_0_2, inv_eig_values_mass_0_3 )
    call self%inverse_mass_1(2)%create( n_dofs, inv_eig_values_mass_0_1, inv_eig_values_mass_1_2, inv_eig_values_mass_0_3 )
    call self%inverse_mass_1(3)%create( n_dofs, inv_eig_values_mass_0_1, inv_eig_values_mass_0_2, inv_eig_values_mass_1_3 )

    call self%inverse_mass_2(1)%create( n_dofs, inv_eig_values_mass_0_1, inv_eig_values_mass_1_2, inv_eig_values_mass_1_3 )
    call self%inverse_mass_2(2)%create( n_dofs, inv_eig_values_mass_1_1, inv_eig_values_mass_0_2, inv_eig_values_mass_1_3 )
    call self%inverse_mass_2(3)%create( n_dofs, inv_eig_values_mass_1_1, inv_eig_values_mass_1_2, inv_eig_values_mass_0_3 )


    ! Poisson solver based on fft inversion
    call self%poisson_fft%init( self%n_dofs, self%s_deg_0, self%delta_x )

    ! Next put together the 1d parts of the 3d Kronecker product
    call sll_s_spline_fem_mass1d( self%n_dofs(1), self%s_deg_0(1), self%mass_line0_1, self%mass1d(1,1) )
    call sll_s_spline_fem_mass1d( self%n_dofs(2), self%s_deg_0(2), self%mass_line0_2, self%mass1d(1,2) )
    call sll_s_spline_fem_mass1d( self%n_dofs(3), self%s_deg_0(3), self%mass_line0_3, self%mass1d(1,3) )

    call sll_s_spline_fem_mass1d( self%n_dofs(1), self%s_deg_1(1), self%mass_line1_1, self%mass1d(2,1) )
    call sll_s_spline_fem_mass1d( self%n_dofs(2), self%s_deg_1(2), self%mass_line1_2, self%mass1d(2,2) )
    call sll_s_spline_fem_mass1d( self%n_dofs(3), self%s_deg_1(3), self%mass_line1_3, self%mass1d(2,3) )

    ! Only for Schur complement eb solver
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

 
    call self%linear_op_schur_eb%create( eig_values_mass_0_1, eig_values_mass_0_2, eig_values_mass_0_3, eig_values_mass_1_1, eig_values_mass_1_2, eig_values_mass_1_3, self%n_dofs, self%delta_x )
   
    if(self%adiabatic_electrons) then
       call sll_s_spline_fem_mass3d( self%n_dofs, s_deg_0, -1, self%mass0, profile_m0 )
       call self%mass0_solver%create( self%mass0 )
       self%mass0_solver%atol = self%mass_solver_tolerance
       !self%mass0_solver%verbose = .true.
    end if

  contains
    function profile_m0( x, component)
      sll_real64 :: profile_m0
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_m0 = product(self%Lx) * self%profile%rho_0( x(1) )/ self%profile%T_e( x(1) )

    end function profile_m0

    function profile_0( x, component)
      sll_real64 :: profile_0
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_0 = product(self%Lx) 

    end function profile_0


  end subroutine init_3d_fem_fft


  !> Finalization
  subroutine free_3d_fem_fft(self)
    class(sll_t_maxwell_3d_fem_fft) :: self !< Maxwell solver class
    !local variable
    sll_int32 :: j

    call self%inverse_mass_0%free()
    do j=1, 3
       call self%inverse_mass_1(j)%free()
       call self%inverse_mass_2(j)%free()
    end do
    call self%poisson_fft%free()
    deallocate( self%work )
    deallocate( self%work3d )
    deallocate( self%work2 )
    deallocate( self%work_d1 ) 
    deallocate( self%work_d2_in ) 
    deallocate( self%work_d2_out ) 
    deallocate( self%work_d3_in ) 
    deallocate( self%work_d3_out )

    deallocate( self%mass_line0_1 )
    deallocate( self%mass_line1_1 )
    deallocate( self%mass_line0_2 )
    deallocate( self%mass_line1_2 )
    deallocate( self%mass_line0_3 )
    deallocate( self%mass_line1_3 )
    deallocate( self%mass_line_mixed_1 )
    deallocate( self%mass_line_mixed_2 )
    deallocate( self%mass_line_mixed_3 )

  end subroutine free_3d_fem_fft


  !> Multiply by dicrete gradient matrix
  subroutine multiply_g( self, field_in, field_out )
    class(sll_t_maxwell_3d_fem_fft) :: self !< Maxwell solver class
    sll_real64, intent( in    )   :: field_in(:)  !< field_in
    sll_real64, intent(   out )   :: field_out(:) !< G*field_in

    call sll_s_multiply_g(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_g


  !> Multiply by transpose of dicrete gradient matrix
  subroutine multiply_gt(self, field_in, field_out)
    class(sll_t_maxwell_3d_fem_fft)  :: self !< Maxwell solver class
    sll_real64, intent( in    )  :: field_in(:) !< field_in
    sll_real64, intent(   out )  :: field_out(:) !< G^T*field_in

    call sll_s_multiply_gt(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_gt


  !> Multiply by discrete curl matrix
  subroutine multiply_c(self, field_in, field_out)
    class(sll_t_maxwell_3d_fem_fft)  :: self     !< Maxwell solver class
    sll_real64, intent( in    )  :: field_in(:)  !< field_in    
    sll_real64, intent(   out )  :: field_out(:) !< C*field_in

    call sll_s_multiply_c(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_c


  !> Multiply by transpose of discrete curl matrix
  subroutine multiply_ct(self, field_in, field_out)
    class(sll_t_maxwell_3d_fem_fft)  :: self     !< Maxwell solver class
    sll_real64, intent( in    )  :: field_in(:)  !< field_in
    sll_real64, intent(   out )  :: field_out(:) !< C^T*field_in

    call sll_s_multiply_ct(self%n_dofs, self%delta_x, field_in, field_out)

  end subroutine multiply_ct


  !> Multiply by the mass matrix 
  subroutine multiply_mass_3d_fem_fft( self, deg, coefs_in, coefs_out )
    class(sll_t_maxwell_3d_fem_fft) :: self         !< Maxwell solver class
    sll_int32,  intent( in    )   :: deg(:)       !< \a deg/form specifies if we multiply the mass to a  1- or 2-form or a mix of both
    sll_real64, intent( in    )   :: coefs_in(:)  !< Coefficient for each DoF
    sll_real64, intent(   out )   :: coefs_out(:) !< Coefficient for each DoF

    if( size(deg) ==1 )then
       select case(deg(1))
       case(0)
          call multiply_mass_3dkron(  self, self%mass_line0_1, &
               self%mass_line0_2, self%mass_line0_3, &
               coefs_in, coefs_out )
       case(1)
          call multiply_mass_1form( self, coefs_in, coefs_out )
       case(2)
          call multiply_mass_2form( self, coefs_in, coefs_out )
       case default
          SLL_ERROR('maxwell_3d_fem_fft','multiply mass for other form not yet implemented')
       end select
    else if( size(deg) == 3 ) then
       call multiply_mass_all( self, deg, coefs_in, coefs_out )
    end if

  end subroutine multiply_mass_3d_fem_fft


  !> Helper function for multiply_mass
  subroutine multiply_mass_1form( self, coefs_in, coefs_out )
    class(sll_t_maxwell_3d_fem_fft), intent( inout ) :: self !< Maxwell solver class
    sll_real64, intent( in    )                      :: coefs_in(:) !< Coefficient for each DoF
    sll_real64, intent(   out )                      :: coefs_out(:) !< Coefficient for each DoF
    !local variables
    sll_int32:: iend, istart

    istart = 1
    iend = self%n_total
    call multiply_mass_3dkron(  self, self%mass_line1_1, &
         self%mass_line0_2, self%mass_line0_3, &
         coefs_in(istart:iend), coefs_out(istart:iend) )
    istart = iend+1
    iend = iend + self%n_total
    call multiply_mass_3dkron(  self, self%mass_line0_1, &
         self%mass_line1_2, self%mass_line0_3, &
         coefs_in(istart:iend), coefs_out(istart:iend) )
    istart = iend+1
    iend = iend + self%n_total
    call multiply_mass_3dkron(  self, self%mass_line0_1, &
         self%mass_line0_2, self%mass_line1_3, &
         coefs_in(istart:iend), coefs_out(istart:iend) )

  end subroutine multiply_mass_1form


  !> Helper function for multiply_mass
  subroutine multiply_mass_2form( self, coefs_in, coefs_out )
    class(sll_t_maxwell_3d_fem_fft), intent( inout )  :: self !< Maxwell solver class
    sll_real64, intent(in)     :: coefs_in(:) !< Coefficient for each DoF
    sll_real64, intent(out)  :: coefs_out(:) !< Coefficient for each DoF
    !local variables
    sll_int32:: iend, istart

    istart = 1
    iend = self%n_total
    call multiply_mass_3dkron(  self, self%mass_line0_1, &
         self%mass_line1_2, self%mass_line1_3, &
         coefs_in(istart:iend), coefs_out(istart:iend) )
    istart = iend+1
    iend = iend + self%n_total
    call multiply_mass_3dkron(  self, self%mass_line1_1, &
         self%mass_line0_2, self%mass_line1_3, &
         coefs_in(istart:iend), coefs_out(istart:iend) )
    istart = iend+1
    iend = iend + self%n_total
    call multiply_mass_3dkron(  self, self%mass_line1_1, &
         self%mass_line1_2, self%mass_line0_3, &
         coefs_in(istart:iend), coefs_out(istart:iend) )

  end subroutine multiply_mass_2form


  !> Multiply by the mass matrix 
  subroutine multiply_mass_3dkron(  self, mass_line_1, mass_line_2, mass_line_3, coefs_in, coefs_out )
    class(sll_t_maxwell_3d_fem_fft), intent( inout )  :: self !< Maxwell solver class
    !sll_int32, intent( in ) :: deg(3) !< \a deg(i) specifies the degree of the 1d mass matrix in dimension \a i (Note: 1 for 0-form, 2 for 1-form, 3 for 0-1-form mix)
    sll_real64, intent(in)    :: mass_line_1(:) !< massline
    sll_real64, intent(in)    :: mass_line_2(:) !< massline
    sll_real64, intent(in)    :: mass_line_3(:) !< massline
    sll_real64, intent(in)     :: coefs_in(:) !< Coefficient for each DoF
    sll_real64, intent(out)  :: coefs_out(:)  !< Coefficient for each DoF 
    ! Local variables
    sll_int32 :: i,j,k,istart,iend
    sll_int32 :: deg(3)

    deg(1) = size(mass_line_1)-1
    deg(2) = size(mass_line_2)-1
    deg(3) = size(mass_line_3)-1

    istart = 1
    iend = self%n_dofs(1)
    do k=1,self%n_dofs(3)
       do j=1,self%n_dofs(2)
          !print*, coefs_in(istart:iend)
          call sll_s_spline_fem_multiply_mass ( self%n_dofs(1), deg(1), &
               mass_line_1, coefs_in(istart:iend), self%work_d1 )           
          !print*, self%mass_line_1(:,1)
          !print*, self%work_d1
          !stop
          self%work3d(:,j,k) = self%work_d1
          istart = iend+1
          iend = iend + self%n_dofs(1)
       end do
    end do

    do k=1,self%n_dofs(3)
       do i =1,self%n_dofs(1)
          self%work_d2_in = self%work3d(i,:,k)
          call sll_s_spline_fem_multiply_mass ( self%n_dofs(2), deg(2), &
               mass_line_2, self%work_d2_in, self%work_d2_out )  
          self%work3d(i,:,k) = self%work_d2_out
       end do
    end do

    istart = 1
    do j=1,self%n_dofs(2)
       do i =1,self%n_dofs(1)
          self%work_d3_in = self%work3d(i,j,:)
          call sll_s_spline_fem_multiply_mass ( self%n_dofs(3), deg(3), &
               mass_line_3, self%work_d3_in, self%work_d3_out ) 

          do k=1,self%n_dofs(3)
             coefs_out(istart+(k-1)*self%n_dofs(1)*self%n_dofs(2)) = self%work_d3_out(k)
          end do
          istart = istart +1
       end do
    end do

  end subroutine multiply_mass_3dkron


  !> Multiply by the mass matrix 
  subroutine multiply_mass_all(  self, deg, coefs_in, coefs_out )
    class(sll_t_maxwell_3d_fem_fft) :: self !< Maxwell solver class
    sll_int32,  intent( in    )     :: deg(:) !< \a deg(i) specifies the degree of the 1d mass matrix in dimension \a i (Note: 1 for 0-form, 2 for 1-form, 3 for 0-1-form mix)
    sll_real64, intent( in    )    :: coefs_in(:) !< Coefficient for each DoF
    sll_real64, intent(   out )    :: coefs_out(:) !< Coefficient for each DoF 
    ! Local variables
    sll_int32 :: i,j,k,istart,iend

    istart = 1
    iend = self%n_dofs(1)
    do k=1,self%n_dofs(3)
       do j=1,self%n_dofs(2)
          select case ( deg(1) )
          case( 1 )
             call sll_s_spline_fem_multiply_mass ( self%n_dofs(1), self%s_deg_0(1), &
                  self%mass_line0_1, coefs_in(istart:iend), self%work_d1 )
          case( 2 )
             call sll_s_spline_fem_multiply_mass ( self%n_dofs(1), self%s_deg_1(1), &
                  self%mass_line1_1, coefs_in(istart:iend), self%work_d1 )
          case ( 3 )
             call sll_s_spline_fem_multiply_massmixed ( self%n_dofs(1), self%s_deg_0(1), &
                  self%mass_line_mixed_1, coefs_in(istart:iend), self%work_d1 )
          case   default
             print*, 'multiply_mass is not implemented for that spline degree.'
          end select
          self%work3d(:,j,k) = self%work_d1
          istart = iend+1
          iend = iend + self%n_dofs(1)
       end do
    end do

    do k=1,self%n_dofs(3)
       do i =1,self%n_dofs(1)
          self%work_d2_in = self%work3d(i,:,k)
          select case ( deg(2) )
          case( 1 )
             call sll_s_spline_fem_multiply_mass ( self%n_dofs(2), self%s_deg_0(2), &
                  self%mass_line0_2, self%work_d2_in, self%work_d2_out )
          case( 2 )
             call sll_s_spline_fem_multiply_mass ( self%n_dofs(2), self%s_deg_1(2), &
                  self%mass_line1_2, self%work_d2_in, self%work_d2_out )
          case ( 3 )
             call sll_s_spline_fem_multiply_massmixed ( self%n_dofs(2), self%s_deg_0(2), &
                  self%mass_line_mixed_2, self%work_d2_in, self%work_d2_out )
          case   default
             print*, 'multiply_mass is not implemented for that spline degree.'
          end select
          self%work3d(i,:,k) = self%work_d2_out
       end do
    end do

    istart = 1
    do j=1,self%n_dofs(2)
       do i =1,self%n_dofs(1)
          self%work_d3_in = self%work3d(i,j,:)
          select case ( deg(3) )
          case( 1 )
             call sll_s_spline_fem_multiply_mass ( self%n_dofs(3), self%s_deg_0(3), &
                  self%mass_line0_3, self%work_d3_in, self%work_d3_out )
          case( 2 )
             call sll_s_spline_fem_multiply_mass ( self%n_dofs(3), self%s_deg_1(3), &
                  self%mass_line1_3, self%work_d3_in, self%work_d3_out )
          case ( 3 )
             call sll_s_spline_fem_multiply_massmixed ( self%n_dofs(3), self%s_deg_0(3), &
                  self%mass_line_mixed_3, self%work_d3_in, self%work_d3_out )
          case   default
             print*, 'multiply_mass is not implemented for that spline degree.'
          end select

          do k=1,self%n_dofs(3)
             coefs_out(istart+(k-1)*self%n_dofs(1)*self%n_dofs(2)) = self%work_d3_out(k)
          end do
          istart = istart +1
       end do
    end do


  end subroutine multiply_mass_all


  !> Multiply by the inverse mass matrix 
  subroutine multiply_mass_inverse_all(  self, form, coefs_in, coefs_out )
    class(sll_t_maxwell_3d_fem_fft) :: self !< Maxwell solver class
    sll_int32,  intent( in    )     :: form !< \a form specifies the form (Note: 0 for 0-form, 1 for 1-form, 2 for 2-form, 3 for 0-1-form mix)
    sll_real64, intent( in    )    :: coefs_in(:) !< Coefficient for each DoF
    sll_real64, intent(   out )    :: coefs_out(:) !< Coefficient for each DoF

    select case( form )
    case( 0 )
       call self%inverse_mass_0%solve( coefs_in, coefs_out )
    case( 1 )
       call self%inverse_mass_1(1)%solve( coefs_in(1:self%n_total), coefs_out(1:self%n_total) )
       call self%inverse_mass_1(2)%solve( coefs_in(self%n_total+1:2*self%n_total), coefs_out(self%n_total+1:2*self%n_total) )
       call self%inverse_mass_1(3)%solve( coefs_in(2*self%n_total+1:3*self%n_total), coefs_out(2*self%n_total+1:3*self%n_total) )
    case( 2 )
       call self%inverse_mass_2(1)%solve( coefs_in(1:self%n_total), coefs_out(1:self%n_total) )
       call self%inverse_mass_2(2)%solve( coefs_in(self%n_total+1:2*self%n_total), coefs_out(self%n_total+1:2*self%n_total) )
       call self%inverse_mass_2(3)%solve( coefs_in(2*self%n_total+1:3*self%n_total), coefs_out(2*self%n_total+1:3*self%n_total) )
    case  default
       print*, 'inverse_mass for', form, '-form not implemented.'
    end select



  end subroutine multiply_mass_inverse_all


end module sll_m_maxwell_3d_fem_fft
