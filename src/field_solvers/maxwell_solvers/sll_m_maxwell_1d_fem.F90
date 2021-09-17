!> @ingroup maxwell_solvers
!> @brief
!> Solve Maxwell's equations in 1D 
!> @details
!> This version corresponds to conforming spline finite elements 
!> @authors
!> Katharina Kormann
!> Eric Sonnendrücker
!>
module sll_m_maxwell_1d_fem
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
       sll_p_pi, &
       sll_p_twopi

  use sll_m_fft, only: &
       sll_t_fft, &
       sll_s_fft_init_r2r_1d, &
       sll_s_fft_exec_r2r_1d, &
       sll_s_fft_free, &
       sll_p_fft_forward, &
       sll_p_fft_backward

  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_linear_operator_schur_eb_1d, only : &
       sll_t_linear_operator_schur_eb_1d

  use sll_m_linear_solver_mgmres, only : &
       sll_t_linear_solver_mgmres

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_maxwell_1d_base, only: &
       sll_i_function_1d_real64, &
       sll_c_maxwell_1d_base

  use sll_m_spline_fem_utilities, only: &
       sll_s_spline_fem_mass_line, &
       sll_s_spline_fem_mixedmass_line, &
       sll_s_spline_fem_interpolation_eigenvalues, &
       sll_s_multiply_g_1d, &
       sll_s_multiply_gt_1d

  use sll_m_spline_fem_utilities_sparse, only: &
       sll_s_spline_fem_mass1d, &
       sll_s_spline_fem_mixedmass1d


  implicit none

  public :: &
       sll_t_maxwell_1d_fem

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Maxwell solver class
  type, extends(sll_c_maxwell_1d_base) :: sll_t_maxwell_1d_fem

     type( sll_t_linear_operator_schur_eb_1d ) :: linear_op_schur_eb !< Schur complement operator for advect_eb
     type( sll_t_linear_solver_mgmres )        :: linear_solver_schur_eb !< Schur complement solver for advect_eb

     type(sll_t_matrix_csr)  :: mass0 !< 0-form mass matrix
     type(sll_t_matrix_csr)  :: mass1 !< 1-form mass matrix
     type(sll_t_matrix_csr)  :: mixed_mass !< mixed mass matrix

     sll_real64, allocatable :: mass_0(:)      !< coefficients of 0-form mass matrix
     sll_real64, allocatable :: mass_1(:)      !< coefficients of 1-form mass matrix
     sll_real64, allocatable :: eig_mass0(:)   !< eigenvalues of circulant 0-form mass matrix
     sll_real64, allocatable :: eig_mass1(:)   !< eigenvalues of circulant 1-form mass matrix
     sll_real64, allocatable :: eig_weak_ampere(:)  !< eigenvalues of circulant update matrix for Ampere
     sll_real64, allocatable :: eig_weak_poisson(:) !< eigenvalues of circulant update matrix for Poisson
     sll_real64, allocatable :: eig_strong_poisson(:) !< eigenvalues of circulant update matrix for Poisson
     sll_real64, allocatable :: eig_splines0(:) !< eigenvalues of interpolation matrix
     sll_real64, allocatable :: eig_splines1(:) !< eigenvalues of interpolation matrix
     sll_real64, allocatable :: eig_isplines0(:) !< eigenvalues of inverse interpolation matrix
     sll_real64, allocatable :: eig_itsplines0(:) !< eigenvalues of inverse interpolation matrix
     sll_real64, allocatable :: eig_kit_mass_0(:) !< eigenvalues of inverse transposed interpolation matrix times inverse mass matrix
     sll_real64, allocatable :: eig_kit_mass_1(:) !< eigenvalues of inverse transposed interpolation matrix times inverse mass matrix

     sll_real64, allocatable :: eig_eb(:)      !< eigenvalues for the Schur complement in a E-B solver
     sll_real64, allocatable :: eig_mass0_inv(:) !< eigenvalues of the inverse 0-form mass matrix
     sll_real64, allocatable :: eig_mass1_inv(:) !< eigenvalues of the inverse 1-form mass matrix

     type(sll_t_fft) :: plan_fw !< fft plan (forward) 
     type(sll_t_fft) :: plan_bw !< fft plan (backward)
     sll_real64, allocatable :: wsave(:) !< scratch data
     sll_real64, allocatable :: work(:)  !< scratch data

     sll_real64 :: force_sign = 1._f64  !< coefficient for the charge (multiplies rho and j)
     logical ::  adiabatic_electrons  !< Set true if solver with adiabatic electrions


   contains
     procedure :: &
          compute_E_from_B => sll_s_compute_e_from_b_1d_fem !< Solve E and B part of Amperes law with B constant in time
     procedure :: &
          compute_B_from_E => sll_s_compute_b_from_e_1d_fem !< Solve Faraday equation with E constant in time
     procedure :: &
          compute_curl_part => sll_s_compute_curl_part_1d_fem  !< Solve source-free Maxwell's equations 
     procedure :: &
          compute_E_from_rho => sll_s_compute_e_from_rho_1d_fem !< Solve E from rho using Poisson
     procedure :: &
          compute_rho_from_e => compute_rho_from_e_1d_fem !< Compute rho from E
     procedure :: &
          compute_E_from_j => choose_interpolation !< choose compute_E_from_j_1d_fem or compute_E_from_j_1d_fem_shape
     procedure :: &
          compute_phi_from_rho => compute_phi_from_rho_1d_fem !< Compute phi from rho (by solving the quasi-neutrality equation)
     procedure :: &
          compute_phi_from_j => compute_phi_from_j_1d_fem !< Compute phi from j (dynamic of quasi-neutrality equation for adiabatic electrons)
     procedure :: &
          compute_rhs_from_function => sll_s_compute_fem_rhs !< Compute integral over given function tested by the basis
     procedure :: &
          L2projection => L2projection_1d_fem !< Compute L_2 projection of a given function
     procedure :: &
          L2norm_squared => L2norm_squared_1d_fem !< Compute the square of the L2 norm of a given vector
     procedure :: &
          inner_product => inner_product_1d_fem !< Inner product of two dof-vectors with mass matrix
     procedure :: &
          init => init_1d_fem !< Initialize the Maxwell class
     procedure :: &
          free => free_1d_fem !< Free Maxwell class
     procedure :: &
          multiply_g !< Multiplication with gradient matrix
     procedure :: &
          multiply_gt !< Multiplication with divergence matrix 
     procedure :: &
          invert_mass => invert_mass_1d_fem !< Invert mass matrix
     procedure :: &
          multiply_mass => multiply_mass_1d_fem !< Product with the mass matrix
     procedure :: &
          transform_dofs => transform_dofs_1d_fem
     procedure :: &
          multiply_interpolation_inverse_transpose
     procedure :: &
          solve_e_b => solve_e_b_1d_fem


  end type sll_t_maxwell_1d_fem

contains

  !> compute Ey from Bz using weak Ampere formulation 
  subroutine sll_s_compute_e_from_b_1d_fem(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_1d_fem) :: self  !< Maxwell solver class
    sll_real64, intent(in)      :: delta_t   !< Time step
    sll_real64, intent(in)      :: field_in(:)  !< Bz
    sll_real64, intent(inout)   :: field_out(:)  !< Ey


    if(self%strong_ampere .eqv. .true.) then 
       call strong_curl(self, delta_t, field_in, field_out)
    else
       call weak_curl(self, delta_t, field_in, field_out)
    end if

  end subroutine sll_s_compute_e_from_b_1d_fem


  subroutine sll_s_compute_b_from_e_1d_fem(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_1d_fem) :: self  !< Maxwell solver class
    sll_real64, intent(in)      :: delta_t !< Time step
    sll_real64, intent(in)      :: field_in(:)  !< Ey
    sll_real64, intent(inout)   :: field_out(:) !< Bz

    if(self%strong_ampere .eqv. .true.) then 
       call weak_curl(self, delta_t, field_in, field_out)
    else
       call strong_curl(self, delta_t, field_in, field_out)
    end if

  end subroutine sll_s_compute_b_from_e_1d_fem


  subroutine sll_s_compute_curl_part_1d_fem( self, delta_t, efield, bfield, betar )
    class(sll_t_maxwell_1d_fem) :: self  !< Maxwell solver class
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

    self%work = 0._f64
    self%wsave = 0._f64

    ! Compute D^T M2 b
    call self%multiply_mass( bfield, self%work, self%s_deg_1 )
    call self%multiply_gt( self%work, self%wsave )

    self%linear_op_schur_eb%sign = -delta_t**2*0.25_f64*factor
    call self%linear_op_schur_eb%dot( efield, self%work )
    self%work = self%work + delta_t*factor*self%wsave

    ! Save efield dofs from previous time step for B field update
    self%wsave = efield

    ! Invert Schur complement matrix
    self%linear_op_schur_eb%sign = delta_t**2*0.25_f64*factor
    call self%linear_solver_schur_eb%set_guess( efield )
    call self%linear_solver_schur_eb%solve( self%work, efield )

    ! Update B field
    self%wsave = self%wsave + efield
    call self%compute_B_from_E( delta_t*0.5_f64, self%wsave, bfield )

  end subroutine sll_s_compute_curl_part_1d_fem


  !> compute Ey from Bz using weak Ampere formulation 
  subroutine weak_curl(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_1d_fem) :: self !< Maxwell solver class
    sll_real64, intent(in)     :: delta_t   !< Time step
    sll_real64, intent(in)     :: field_in(:)  !< input field
    sll_real64, intent(inout)  :: field_out(:)  !< output_field
    ! local variables
    sll_real64 :: coef

    ! Compute potential weak curl of Bz using eigenvalue of circulant inverse matrix
    call solve_circulant(self, self%eig_weak_ampere, field_in, self%work)
    ! Update bz from self value
    coef = delta_t/self%delta_x
    field_out =  field_out + coef*self%work

  end subroutine weak_curl


  !> Compute Bz from Ey using strong 1D Faraday equation for spline coefficients
  !> $B_z^{new}(x_j) = B_z^{old}(x_j) - \frac{\Delta t}{\Delta x} (E_y(x_j) - E_y(x_{j-1})  $
  !subroutine sll_s_compute_b_from_e_1d_fem
  subroutine strong_curl(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_1d_fem)  :: self !< Maxwell solver class
    sll_real64, intent(in)     :: delta_t !< Time step
    sll_real64, intent(in)     :: field_in(:)  !< Ey
    sll_real64, intent(inout)  :: field_out(:) !< Bz 
    ! local variables
    sll_real64 :: coef
    sll_int32 :: i

    coef = delta_t/self%delta_x
    ! relation betwen spline coefficients for strong Ampere
    do i=2,self%n_cells
       field_out(i) = field_out(i) + coef * ( field_in(i-1) - field_in(i) )
    end do
    ! treat Periodic point
    field_out(1) = field_out(1) + coef * ( field_in(self%n_cells) - field_in(1) )
  end subroutine strong_curl


  subroutine sll_s_compute_e_from_rho_1d_fem(self, field_in, field_out )       
    class(sll_t_maxwell_1d_fem) :: self !< Maxwell solver class
    sll_real64, intent(in)  :: field_in(:)  !< rho
    sll_real64, intent(out) :: field_out(:) !< E
    ! local variables
    sll_int32 :: i

    self%wsave = self%force_sign * field_in

    if ( self%strong_ampere .eqv. .false. ) then
       ! Compute potential phi from rho, using eigenvalue of circulant inverse matrix
       call solve_circulant(self, self%eig_weak_poisson, self%wsave, self%work)
       ! Compute spline coefficients of Ex from those of phi
       do i=2,self%n_cells
          field_out(i) =  (self%work(i-1) -  self%work(i)) !* (self%delta_x)
       end do
       ! treat Periodic point
       field_out(1) = (self%work(self%n_cells) - self%work(1)) !* (self%delta_x)
    else
       call solve_circulant( self, self%eig_isplines0, self%wsave, self%work )
       call solve_circulant( self, self%eig_strong_poisson, self%work, field_out )

    end if


  end subroutine sll_s_compute_e_from_rho_1d_fem


  !> Compute rho from Gauss law for given efield
  subroutine compute_rho_from_e_1d_fem(self, field_in, field_out)
    class(sll_t_maxwell_1d_fem)  :: self !< Maxwell solver class
    sll_real64, intent(in)     :: field_in(:)  !< field_in
    sll_real64, intent(out)  :: field_out(:) !< rho

    sll_int32 :: i

    if ( self%strong_ampere .eqv. .true. ) then

       ! The differentiation matrix
       do i=2,self%n_cells
          self%work(i) = self%force_sign * ( field_in(i) - field_in(i-1) )/ self%delta_x
       end do
       ! treat Periodic point
       self%work(1) =  self%force_sign * (  field_in(1) - field_in(self%n_cells) )/ self%delta_x

       call solve_circulant( self, self%eig_splines0 , self%work, field_out )
    else
       call solve_circulant(self, self%eig_mass1, field_in, self%work)

       ! relation betwen spline coefficients for strong Ampere
       do i=1,self%n_cells-1
          field_out(i) = self%force_sign * ( self%work(i+1) - self%work(i) )
       end do
       ! treat Periodic point
       field_out(self%n_cells) =  self%force_sign * ( self%work(1) - self%work(self%n_cells) )
    end if

  end subroutine compute_rho_from_e_1d_fem


  !> Choose between compute_E_from_j_1d_fem and compute_E_from_j_1d_fem_shape
  subroutine choose_interpolation(self, current, component, E)
    class(sll_t_maxwell_1d_fem)             :: self !< Maxwell solver class
    sll_real64,dimension(:),intent(in)    :: current !< Component \a component of the current integrated over time interval
    sll_int32, intent(in)                 :: component !< Component of the Efield to be computed
    sll_real64,dimension(:),intent(inout) :: E !< Updated electric field

    if(self%strong_ampere .eqv. .true.) then 
       call compute_E_from_j_1d_fem_shape(self, current, component, E)
    else
       call compute_E_from_j_1d_fem(self, current, component, E)
    end if

  end subroutine choose_interpolation


  !> Compute E_i from j_i integrated over the time interval using weak Ampere formulation  
  subroutine compute_E_from_j_1d_fem(self, current, component, E)
    class(sll_t_maxwell_1d_fem)           :: self !< Maxwell solver class
    sll_real64,dimension(:),intent(in)    :: current !< Component \a component of the current integrated over time interval
    sll_int32, intent(in)                 :: component !< Component of the Efield to be computed
    sll_real64,dimension(:),intent(inout) :: E !< Updated electric field

    ! Multiply by inverse mass matrix  using the eigenvalues of the circulant inverse matrix
    if (component == 1) then
       call solve_circulant(self, self%eig_mass1_inv, current, self%work)
    elseif (component == 2) then
       call solve_circulant(self, self%eig_mass0_inv, current, self%work)
    else
       print*, 'Component ', component, 'not implemented in compute_E_from_j_1d_fem.'
    end if


    ! Update the electric field and scale
    E = E - self%force_sign * self%work/self%delta_x

  end subroutine compute_E_from_j_1d_fem

  subroutine compute_E_from_j_1d_fem_shape(self, current, component, E)
    class(sll_t_maxwell_1d_fem)             :: self !< Maxwell solver class
    sll_real64,dimension(:),intent(in)    :: current !< Component \a component of the current integrated over time interval
    sll_int32, intent(in)                 :: component !< Component of the Efield to be computed
    sll_real64,dimension(:),intent(inout) :: E !< Updated electric field

    ! Multiply by inverse mass matrix  using the eigenvalues of the circulant inverse matrix
    call solve_circulant(self, self%eig_isplines0, current, self%work)

    ! Update the electric field and scale
    E = E - self%work!/self%delta_x

  end subroutine compute_E_from_j_1d_fem_shape


  subroutine multiply_interpolation_inverse_transpose(self, in, out)
    class(sll_t_maxwell_1d_fem)             :: self !< Maxwell solver class
    sll_real64,dimension(:),intent(in)    :: in !< input
    sll_real64,dimension(:),intent(inout) :: out !< output

    ! Multiply by inverse mass matrix  using the eigenvalues of the circulant inverse matrix
    call solve_circulant(self, self%eig_itsplines0, in, out)

  end subroutine multiply_interpolation_inverse_transpose


  !> Solves for the efield part in a implicit curl part solve
  subroutine solve_e_b_1d_fem(self, delta_t, bfield, efield)
    class(sll_t_maxwell_1d_fem)             :: self !< Maxwell solver class
    sll_real64, intent(in)   :: delta_t !< Time step
    sll_real64,intent(in)    :: bfield(:) !< Component \a component of the current integrated over time interval
    sll_real64,intent(inout) :: efield(:) !< Updated electric field


    call self%compute_e_from_b( delta_t, bfield, efield )

    call solve_circulant( self, self%eig_eb, efield, self%work )
    efield = self%work

  end subroutine solve_e_b_1d_fem


  subroutine solve_circulant(self, eigvals, rhs, res)
    class(sll_t_maxwell_1d_fem) :: self !< Maxwell solver class
    sll_real64, intent(in) :: eigvals(:)    !< eigenvalues of circulant matrix
    sll_real64, intent(in) :: rhs(:) !< rhs
    sll_real64, intent(out) :: res(:) !< result
    ! local variables
    sll_int32 :: k
    sll_real64 :: re, im 

    ! Compute res from rhs, using eigenvalue of circulant  matrix
    res = rhs
    ! Forward FFT
    call sll_s_fft_exec_r2r_1d ( self%plan_fw, res, self%wsave )
    self%wsave(1) = self%wsave(1) * eigvals(1)
    do k=2, (self%n_cells+1)/2
       re = self%wsave(k) * eigvals(k) - &
            self%wsave(self%n_cells-k+2) * eigvals(self%n_cells-k+2)
       im = self%wsave(k) * eigvals(self%n_cells-k+2) + &
            self%wsave(self%n_cells-k+2) * eigvals(k)
       self%wsave(k) = re
       self%wsave(self%n_cells-k+2) = im
    end do
    if (modulo(self%n_cells,2) == 0 ) then
       self%wsave(self%n_cells/2+1) = self%wsave(self%n_cells/2+1)*eigvals(self%n_cells/2+1)
    end if
    ! Backward FFT 
    call sll_s_fft_exec_r2r_1d( self%plan_bw, self%wsave, res )
    ! normalize
    res = res / real(self%n_cells, f64)
  end subroutine solve_circulant


  !> For model with adiabatic electrons
  subroutine compute_phi_from_rho_1d_fem( self, in, phi, efield )
    class(sll_t_maxwell_1d_fem)                :: self  !< Maxwell solver class
    sll_real64, intent(in)                     :: in(:) !< rho
    sll_real64, intent(out)                    :: phi(:) !< phi 
    sll_real64, intent(out)                    :: efield(:) !< E

    ! Compute phi by inverting the mass matrix
    call self%invert_mass( in, phi, self%s_deg_0 )

    ! Compute the degrees of freedom of the electric field as -G phi
    call multiply_g( self, phi, efield )
    efield = -efield

  end subroutine compute_phi_from_rho_1d_fem


  !> For model with adiabatic electrons
  subroutine compute_phi_from_j_1d_fem( self, in, phi, efield )
    class(sll_t_maxwell_1d_fem)                :: self  !< Maxwell solver class
    sll_real64, intent(in)                     :: in(:) !< Current integrated over time interval
    sll_real64, intent(out)                    :: phi(:) !< phi 
    sll_real64, intent(out)                    :: efield(:) !< E

    ! Compute divergence of the current (G^T current) (assuming a 1v model)
    call multiply_gt( self, in, self%work )

    ! Compute phi by inverting the mass matrix
    call self%invert_mass( self%work, self%wsave, self%s_deg_0 )
    phi = phi + self%wsave

    ! Compute the degrees of freedom of the electric field as -G phi
    call multiply_g( self, phi, efield )
    efield = -efield

  end subroutine compute_phi_from_j_1d_fem


  !> Compute the FEM right-hand-side for a given function f and periodic splines of given degree
  !> Its components are $\int f N_i dx$ where $N_i$ is the B-spline starting at $x_i$ 
  subroutine sll_s_compute_fem_rhs(self, func, degree, coefs_dofs)
    class(sll_t_maxwell_1d_fem)             :: self !< Maxwell solver class
    procedure(sll_i_function_1d_real64) :: func !< function
    sll_int32, intent(in) :: degree !< Specify the degree of the basis functions
    sll_real64, intent(out) :: coefs_dofs(:)  ! Finite Element right-hand-side
    ! local variables
    sll_int32 :: i,j,k
    sll_real64 :: coef
    sll_real64, dimension(2,degree+1) :: xw_gauss
    sll_real64, dimension(degree+1,degree+1) :: bspl

    ! take enough Gauss points so that projection is exact for splines of degree deg
    ! rescale on [0,1] for compatibility with B-splines
    xw_gauss = sll_f_gauss_legendre_points_and_weights(degree+1, 0.0_f64, 1.0_f64)
    ! Compute bsplines at gauss_points
    do k=1,degree+1
       call sll_s_uniform_bsplines_eval_basis(degree,xw_gauss(1,k), bspl(:,k))
       !print*, 'bs', bspl(k,:)
    end do

    ! Compute coefs_dofs = int f(x)N_i(x) 
    do i = 1, self%n_cells
       coef=0.0_f64
       ! loop over support of B spline
       do j = 1, degree+1
          ! loop over Gauss points
          do k=1, degree+1
             coef = coef + xw_gauss(2,k)*func(self%delta_x * (xw_gauss(1,k) + real(i + j - 2 - ((i + j - 2)/self%n_cells) * self%n_cells,f64))) * bspl(degree+2-j,k)
             !print*, i,j,k, xw_gauss(2,k), xw_gauss(1,k),f(self%delta_x*(xw_gauss(1,k) + i + j - 2)) 
          enddo
       enddo
       ! rescale by cell size
       coefs_dofs(i) = coef*self%delta_x
    enddo

  end subroutine sll_s_compute_fem_rhs


  !> Compute square of the L2norm 
  function L2norm_squared_1d_fem(self, coefs_dofs, degree) result (r)
    class(sll_t_maxwell_1d_fem) :: self !< Maxwell solver class
    sll_real64 :: coefs_dofs(:) !< Coefficient for each DoF
    sll_int32  :: degree !< Specify the degree of the basis functions
    sll_real64 :: r !< Result: squared L2 norm

    ! Multiply coefficients by mass matrix (use diagonalization FFT and mass matrix eigenvalues)
    if (degree == self%s_deg_0 ) then

       call solve_circulant(self, self%eig_mass0, coefs_dofs, self%work)

    elseif (degree == self%s_deg_1) then

       call solve_circulant(self, self%eig_mass1, coefs_dofs, self%work)

    end if
    ! Multiply by the coefficients from the left (inner product)
    r = sum(coefs_dofs*self%work)
    ! Scale by delt_x
    r = r*self%delta_x

  end function L2norm_squared_1d_fem



  function inner_product_1d_fem( self, coefs1_dofs, coefs2_dofs, degree, degree2 ) result (r)
    class(sll_t_maxwell_1d_fem) :: self !< Maxwell solver class
    sll_real64 :: coefs1_dofs(:) !< Coefficient for each DoF
    sll_real64 :: coefs2_dofs(:) !< Coefficient for each DoF
    sll_int32  :: degree !< Specify the degree of the basis functions
    sll_int32, optional  :: degree2 !< Specify the degree of the basis functions
    sll_real64 :: r !< Result: squared L2 norm

    ! Multiply coefficients by mass matrix (use diagonalization FFT and mass matrix eigenvalues)
    if (present(degree2)) then
       if (degree == degree2) then
          if (degree == self%s_deg_0 ) then
             call solve_circulant(self, self%eig_mass0, coefs2_dofs, self%work)
          elseif (degree == self%s_deg_1) then
             call solve_circulant(self, self%eig_mass1, coefs2_dofs, self%work)
          end if
          ! Multiply by the coefficients from the left (inner product)
          r = sum(coefs1_dofs*self%work)
          ! Scale by delt_x
          r = r*self%delta_x
       else
          if (degree == self%s_deg_0) then
             call self%mixed_mass%dot( coefs2_dofs, self%work )

             ! Multiply by the coefficients from the left (inner product)
             r = sum(coefs1_dofs*self%work)
          else
             call self%mixed_mass%dot( coefs1_dofs, self%work )

             ! Multiply by the coefficients from the left (inner product)
             r = sum(coefs2_dofs*self%work)
          end if
       end if
    else
       if (degree == self%s_deg_0 ) then
          call solve_circulant(self, self%eig_mass0, coefs2_dofs, self%work)
       elseif (degree == self%s_deg_1) then
          call solve_circulant(self, self%eig_mass1, coefs2_dofs, self%work)
       end if
       ! Multiply by the coefficients from the left (inner product)
       r = sum(coefs1_dofs*self%work)
       ! Scale by delt_x
       r = r*self%delta_x
    end if

  end function inner_product_1d_fem


  !> Compute the L2 projection of a given function f on periodic splines of given degree
  subroutine L2projection_1d_fem(self, func, degree, coefs_dofs)
    class(sll_t_maxwell_1d_fem) :: self !< Maxwell solver class
    procedure(sll_i_function_1d_real64) :: func !< function
    sll_int32, intent(in) :: degree !< Specify the degree of the basis functions
    sll_real64, intent(out) :: coefs_dofs(:)  ! spline coefficients of projection

    ! Compute right-hand-side
    call sll_s_compute_fem_rhs(self, func, degree, self%work)

    ! Multiply by inverse mass matrix (! complex numbers stored in real array with fftpack ordering)

    if (degree == self%s_deg_0) then
       call solve_circulant(self, self%eig_mass0_inv, self%work, coefs_dofs)
    elseif  (degree == self%s_deg_0-1) then
       call solve_circulant( self, self%eig_mass1_inv, self%work, coefs_dofs)
    else
       print*, 'degree ', degree, 'not availlable in maxwell_1d_fem object' 
    endif

    ! Account for scaling in the mass matrix by dx
    coefs_dofs = coefs_dofs/self%delta_x

  end subroutine L2projection_1d_fem


  subroutine init_1d_fem( self, domain, n_dofs, s_deg_0, delta_t, strong_ampere, solver_tolerance, force_sign, adiabatic_electrons )
    class(sll_t_maxwell_1d_fem), intent(out) :: self  !< Maxwell solver class
    sll_real64, intent(in) :: domain(2)     !< xmin, xmax
    sll_int32, intent(in) :: n_dofs  !< number of degrees of freedom (here number of cells and grid points)
    sll_int32, intent(in) :: s_deg_0 !< highest spline degree
    sll_real64, intent(in) :: delta_t !< Time step
    logical, intent(in),optional :: strong_ampere !< flag to switch between strong and weak Ampere formulation
    sll_real64, intent(in), optional :: solver_tolerance !< tolerance for Schur complement solver
    sll_real64, intent(in), optional :: force_sign !< sign of particle force
    logical, intent(in), optional :: adiabatic_electrons !< flag if adiabatic electrons are used
    ! local variables
    sll_int32 :: ierr
    sll_real64 :: mass_line_0(s_deg_0+1), mass_line_1(s_deg_0), mass_line_mixed(s_deg_0*2)
    sll_int32 :: j, k ! loop variables
    sll_real64 :: coef0, coef1, sin_mode, cos_mode 
    sll_real64 :: ni !1/n_dofs in double precision
    sll_real64 :: dt
    sll_real64 :: T_e = 1._f64!10000._f64

    ni=1.0_f64/real(n_dofs,f64)

    self%n_dofs0 = n_dofs
    self%n_dofs1 = n_dofs
    self%n_cells = n_dofs
    self%Lx = domain(2) - domain(1)
    self%delta_x = self%Lx *ni
    !print*, 'dx', self%delta_x
    self%s_deg_0 = s_deg_0
    self%s_deg_1 = s_deg_0 - 1

    if (present(strong_ampere)) then
       self%strong_ampere=strong_ampere
    else
       self%strong_ampere=.false.
    endif

    dt = delta_t

    if (present( solver_tolerance) ) then
       self%solver_tolerance = solver_tolerance
    else
       self%solver_tolerance = 1d-12
    end if

    if (present( force_sign) ) then
       self%force_sign = force_sign
    end if

    if (present( adiabatic_electrons) ) then
       self%adiabatic_electrons = adiabatic_electrons
    end if

    SLL_ALLOCATE(self%mass_0(s_deg_0+1), ierr)
    SLL_ALLOCATE(self%mass_1(s_deg_0), ierr)

    select case(s_deg_0)
    case(1) ! linear and constant splines
       ! Upper diagonal coeficients  of linear spline mass matrix (from Eulerian numbers)
       self%mass_0(1) = 4.0_f64/6.0_f64 
       self%mass_0(2) = 1.0_f64/6.0_f64
       ! Upper diagonal coeficients  of constant spline mass matrix
       self%mass_1(1) = 1.0_f64 
    case(2) ! quadratic and linear splines
       ! Upper diagonal coeficients  of quadratic spline mass matrix (from Eulerian numbers)
       self%mass_0(1) = 66.0_f64/120.0_f64 
       self%mass_0(2) = 26.0_f64/120.0_f64
       self%mass_0(3) = 1.0_f64/120.0_f64
       ! Upper diagonal coeficients  of linear spline mass matrix (from Eulerian numbers)
       self%mass_1(1) = 4.0_f64/6.0_f64 
       self%mass_1(2) = 1.0_f64/6.0_f64
    case(3)
       ! Upper diagonal coeficients  of cubic spline mass matrix (from Eulerian numbers)
       self%mass_0(1) = 2416.0_f64/5040.0_f64 
       self%mass_0(2) = 1191.0_f64/5040.0_f64
       self%mass_0(3) = 120.0_f64/5040.0_f64
       self%mass_0(4) = 1.0_f64/5040.0_f64
       ! Upper diagonal coeficients  of quadratic spline mass matrix (from Eulerian numbers)
       self%mass_1(1) = 66.0_f64/120.0_f64 
       self%mass_1(2) = 26.0_f64/120.0_f64
       self%mass_1(3) = 1.0_f64/120.0_f64

    case default
       call sll_s_spline_fem_mass_line ( self%s_deg_0, self%mass_0 )
       call sll_s_spline_fem_mass_line ( self%s_deg_1, self%mass_1 )
    end select

    if(self%adiabatic_electrons) then
       mass_line_0=self%mass_0*self%delta_x/T_e
    else
       mass_line_0= self%mass_0*self%delta_x
    end if
    mass_line_1= self%mass_1*self%delta_x
    call sll_s_spline_fem_mass1d( n_dofs, self%s_deg_0, mass_line_0, self%mass0 )
    call sll_s_spline_fem_mass1d( n_dofs, self%s_deg_1, mass_line_1, self%mass1 )

    if(self%adiabatic_electrons) then
       self%mass_0= self%mass_0/T_e
    end if

    SLL_ALLOCATE(self%eig_mass0(n_dofs), ierr)
    SLL_ALLOCATE(self%eig_mass1(n_dofs), ierr)
    SLL_ALLOCATE(self%eig_mass0_inv(n_dofs), ierr)
    SLL_ALLOCATE(self%eig_mass1_inv(n_dofs), ierr)
    SLL_ALLOCATE(self%eig_weak_ampere(n_dofs), ierr)
    SLL_ALLOCATE(self%eig_weak_poisson(n_dofs), ierr)
    SLL_ALLOCATE(self%eig_strong_poisson(n_dofs), ierr)
    ! for interpolation matrix
    SLL_ALLOCATE(self%eig_splines0(n_dofs), ierr)
    SLL_ALLOCATE(self%eig_splines1(n_dofs), ierr)
    SLL_ALLOCATE(self%eig_isplines0(n_dofs), ierr)
    SLL_ALLOCATE(self%eig_itsplines0(n_dofs), ierr)
    SLL_ALLOCATE(self%eig_kit_mass_0(n_dofs), ierr)
    SLL_ALLOCATE(self%eig_kit_mass_1(n_dofs), ierr)

    self%eig_mass0 = 0.0_f64
    self%eig_mass1 = 0.0_f64
    self%eig_mass0_inv = 0.0_f64
    self%eig_mass1_inv = 0.0_f64
    self%eig_weak_ampere = 0.0_f64 
    self%eig_weak_poisson = 0.0_f64
    ! for the e, b solver
    SLL_ALLOCATE(self%eig_eb(n_dofs), ierr )
    self%eig_eb = 0.0_f64

    ! Initialise FFT     
    SLL_ALLOCATE(self%work(n_dofs),ierr)
    SLL_ALLOCATE(self%wsave(n_dofs),ierr)
    call sll_s_fft_init_r2r_1d( self%plan_fw, n_dofs, self%work, self%wsave, sll_p_fft_forward, normalized=.false. )
    call sll_s_fft_init_r2r_1d( self%plan_bw, n_dofs, self%work, self%work, sll_p_fft_backward, normalized=.false. )

    ! Compute eigenvalues of circulant Ampere update matrix M_0^{-1} D^T M_1
    ! and circulant Poisson Matrix (D^T M_1 D)^{-1}
    ! zero mode vanishes due to derivative matrix D^T
    self%eig_weak_ampere(1) = 0.0_f64 
    self%eig_weak_poisson(1) = 0.0_f64  ! Matrix is not invertible: 0-mode is set to 0
    self%eig_strong_poisson(1) = 0.0_f64  ! Matrix is not invertible: 0-mode is set to 0
    self%eig_mass0(1) = 1.0_f64  ! sum of coefficents is one
    self%eig_mass1(1) = 1.0_f64  ! sum of coefficents is one

    call sll_s_spline_fem_interpolation_eigenvalues( self%s_deg_0, self%n_cells, self%eig_splines0 )
    self%eig_isplines0(1) = 1.0_f64/ self%eig_splines0(1)
    self%eig_eb(1) = 1.0_f64

    do k=1, (n_dofs+1)/2 - 1
       coef0 =  self%mass_0(1)
       coef1 =  self%mass_1(1)
       do j=1,s_deg_0 - 1
          cos_mode = cos(sll_p_twopi*ni*real(j*k, f64))
          coef0 = coef0 + 2* self%mass_0(j+1)*cos_mode
          coef1 = coef1 + 2* self%mass_1(j+1)*cos_mode
       enddo
       ! add last term for larger matrix
       j = s_deg_0

       coef0 = coef0 + 2* self%mass_0(j+1)*cos(sll_p_twopi*ni*real(j*k, f64))

       ! compute eigenvalues
       self%eig_mass0(k+1) = coef0 ! real part
       self%eig_mass0(n_dofs-k+1) = 0.0_f64 ! imaginary part
       self%eig_mass1(k+1) = coef1 ! real part
       self%eig_mass1(n_dofs-k+1) = 0.0_f64 ! imaginary part

       cos_mode = cos(sll_p_twopi*ni*real(k, f64))
       sin_mode = sin(sll_p_twopi*ni*real(k, f64))

       self%eig_weak_ampere(k+1) =  (coef1 / coef0) * (1-cos_mode) ! real part
       self%eig_weak_ampere(n_dofs-k+1) =  -(coef1 / coef0) * sin_mode   ! imaginary part
       self%eig_weak_poisson(k+1) = 1.0_f64 / (coef1 * ((1.0_f64-cos_mode)**2 + &
            sin_mode**2))  ! real part
       self%eig_weak_poisson(n_dofs-k+1) = 0.0_f64  ! imaginary part

       coef0 = (1.0_f64 - cos_mode) * self%eig_splines0(k+1) - &
            sin_mode * self%eig_splines0(n_dofs-k+1)
       coef1 = (1.0_f64 - cos_mode) * self%eig_splines0(n_dofs-k+1) + &
            sin_mode * self%eig_splines0(k+1)

       self%eig_strong_poisson(k+1) = (1.0_f64 - cos_mode)/((1.0_f64-cos_mode)**2 + sin_mode**2)!coef0 /(coef0**2 + coef1**2)
       self%eig_strong_poisson(n_dofs-k+1) = -sin_mode/((1.0_f64-cos_mode)**2 + sin_mode**2)!-coef1 /(coef0**2 + coef1**2)

       self%eig_isplines0(k+1) = self%eig_splines0(k+1)/(self%eig_splines0(k+1)**2 + self%eig_splines0(n_dofs-k+1)**2)
       self%eig_isplines0(n_dofs-k+1) = -self%eig_splines0(n_dofs-k+1)/(self%eig_splines0(k+1)**2 + self%eig_splines0(n_dofs-k+1)**2)


       self%eig_eb(k+1) = 1.0_f64/ ( 1.0_f64 + (dt/self%delta_x)**2*self%eig_mass1(k+1)/self%eig_mass0(k+1)*2.0_f64 *(1.0_f64 -cos_mode)  )
       self%eig_eb(n_dofs-k+1) = 0.0_f64

    enddo

    if ( modulo(n_dofs,2) == 0 ) then
       ! N/2 mode
       coef0 =  self%mass_0(1)
       coef1 =  self%mass_1(1)
       do j=1, s_deg_0 - 1
          coef0 = coef0 + 2 * self%mass_0(j+1)*cos(sll_p_pi*real(j,f64))
          coef1 = coef1 + 2 * self%mass_1(j+1)*cos(sll_p_pi*real(j,f64))
       enddo

       ! add last term for larger matrix
       j = s_deg_0
       coef0 = coef0 + 2 * self%mass_0(j+1)*cos(sll_p_pi*real(j, f64))

       ! compute eigenvalues
       self%eig_mass0(n_dofs/2+1) = coef0
       self%eig_mass1(n_dofs/2+1) = coef1
       self%eig_weak_ampere(n_dofs/2+1) = 2.0_f64 * (coef1 / coef0)
       self%eig_weak_poisson(n_dofs/2+1) = 1.0_f64 / (coef1 *4.0_f64)
       self%eig_strong_poisson(n_dofs/2+1) = 0.5_f64!0.25_f64!self%eig_splines1(n_dofs/2+1)/ 4.0_f64
       self%eig_isplines0(n_dofs/2+1) = 1.0_f64 / self%eig_splines0(n_dofs/2+1)
       self%eig_eb(n_dofs/2+1) = 1.0_f64/ ( 1.0_f64 + (dt/self%delta_x)**2*self%eig_mass1(n_dofs/2+1)/self%eig_mass0(n_dofs/2+1)*4.0_f64  )
    end if

    self%eig_strong_poisson = self%eig_strong_poisson* self%delta_x

    !print*, n_dofs/2+1, size(self%eig_kit_mass_0), size(
    do k=1,n_dofs/2+1
       self%eig_kit_mass_0(k) = self%eig_isplines0(k) * self%eig_mass0(k) * self%delta_x
       self%eig_kit_mass_1(k) = self%eig_isplines0(k) * self%eig_mass1(k) * self%delta_x
       self%eig_itsplines0(k) =  self%eig_isplines0(k) 
    end do
    do k=1, (n_dofs-1)/2 
       self%eig_kit_mass_0(n_dofs-k+1) = -self%eig_isplines0(n_dofs-k+1) * self%eig_mass0(k+1)*self%delta_x
       self%eig_kit_mass_1(n_dofs-k+1) = -self%eig_isplines0(n_dofs-k+1) * self%eig_mass1(k+1)*self%delta_x
       self%eig_itsplines0(n_dofs-k+1) =  -self%eig_isplines0(n_dofs-k+1)
    end do

    do j=1,(self%n_cells+1)/2!+1 !2
       self%eig_mass0_inv(j) = 1.0_f64 / self%eig_mass0(j)
       self%eig_mass1_inv(j) = 1.0_f64 / self%eig_mass1(j)
    end do
    if ( modulo(n_dofs,2) == 0 ) then
       self%eig_mass0_inv(self%n_cells/2+1) = 1.0_f64 / self%eig_mass0(self%n_cells/2+1)
       self%eig_mass1_inv(self%n_cells/2+1) = 1.0_f64 / self%eig_mass1(self%n_cells/2+1)
    end if

    call sll_s_spline_fem_mixedmass_line( self%s_deg_0, mass_line_mixed)
    mass_line_mixed=mass_line_mixed*self%delta_x
    call sll_s_spline_fem_mixedmass1d( n_dofs, self%s_deg_0, mass_line_mixed, self%mixed_mass )

    ! Only for Schur complement eb solver
    call self%linear_op_schur_eb%create( self%mass0, self%mass1, self%n_cells, self%delta_x )
    call self%linear_solver_schur_eb%create( self%linear_op_schur_eb )
    self%linear_solver_schur_eb%atol = self%solver_tolerance
    self%linear_solver_schur_eb%rtol = self%solver_tolerance
    !self%linear_solver_schur_eb%verbose = .true.

  end subroutine init_1d_fem


  subroutine free_1d_fem(self)
    class(sll_t_maxwell_1d_fem) :: self !< Maxwell solver class

    call sll_s_fft_free( self%plan_fw )
    call sll_s_fft_free( self%plan_bw )
    deallocate(self%mass_0)
    deallocate(self%mass_1)
    deallocate(self%eig_mass0)
    deallocate(self%eig_mass1)
    deallocate(self%eig_mass0_inv)
    deallocate(self%eig_mass1_inv)
    deallocate(self%eig_weak_ampere)
    deallocate(self%eig_weak_poisson)
    deallocate(self%wsave)
    deallocate(self%work)
    deallocate(self%eig_splines0)
  end subroutine free_1d_fem

  subroutine multiply_g( self,  in, out)
    class(sll_t_maxwell_1d_fem), intent(in) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)  :: in(:)  !< field_in
    sll_real64, intent(out) :: out(:) !< G*field_in

    call sll_s_multiply_g_1d( self%n_cells, self%delta_x, in, out )

  end subroutine multiply_g


  subroutine multiply_gt( self,  in, out)
    class(sll_t_maxwell_1d_fem), intent(in) :: self !< Maxwell_Clamped solver class
    sll_real64, intent(in)  :: in(:)  !< field_in
    sll_real64, intent(out) :: out(:) !< G^T*field_in

    call sll_s_multiply_gt_1d( self%n_cells, self%delta_x, in, out )

  end subroutine multiply_gt


  subroutine multiply_mass_1d_fem( self,  in, out, degree)
    class(sll_t_maxwell_1d_fem), intent(inout) :: self !< Maxwell solver class
    sll_real64, intent(in)  :: in(:) !< Coefficient for each DoF
    sll_real64, intent(out) :: out(:) !< Coefficient for each DoF
    sll_int32,  intent(in)  :: degree !< Specify the degree of the basis functions

    ! Multiply coefficients by mass matrix (use diagonalization FFT and mass matrix eigenvalues)
    if (degree == self%s_deg_0 ) then
       call solve_circulant(self, self%eig_mass0, in, out)
       out = out*self%delta_x
    elseif (degree == self%s_deg_1) then
       call solve_circulant(self, self%eig_mass1, in, out)
       out = out*self%delta_x
    elseif(degree == 10) then
       call self%mixed_mass%dot( in, out )
    else
       print*, 'maxwell_solver_1d_fem: multiply mass for other form not yet implemented'
       stop
    end if

   


  end subroutine multiply_mass_1d_fem


  !> Invert the mass matrix
  subroutine invert_mass_1d_fem(self, in, out, degree)
    class(sll_t_maxwell_1d_fem), intent(inout) :: self !< Maxwell solver class
    sll_int32, intent(in) :: degree !< Specify the degree of the basis functions
    sll_real64, intent(in) :: in(:)  !< spline coefficients of projection
    sll_real64, intent(out) :: out(:)  !< spline coefficients of projection

    ! Multiply by inverse mass matrix (! complex numbers stored in real array with fftpack ordering)
    if (degree == self%s_deg_0) then
       call solve_circulant( self, self%eig_mass0_inv, in, out )
    elseif  (degree == self%s_deg_0-1) then
       call solve_circulant( self, self%eig_mass1_inv, in, out )
    else
       print*, 'degree ', degree, 'not availlable in maxwell_1d_fem object' 
    endif

    ! Account for scaling in the mass matrix by dx
    out = out/self%delta_x

  end subroutine invert_mass_1d_fem


  !> Invert the mass matrix
  subroutine transform_dofs_1d_fem(self, in, out, degree)
    class(sll_t_maxwell_1d_fem), intent(inout) :: self !< Maxwell solver class
    sll_int32, intent(in) :: degree !< Specify the degree of the basis functions ! this is 0 for 0 form and 1 for 1 form
    sll_real64, intent(in) :: in(:)  !< spline coefficients of projection
    sll_real64, intent(out) :: out(:)  !< spline coefficients of projection
    ! local variables
    sll_int32 :: i

    if(self%strong_ampere .eqv. .true.) then
       if (degree == 0) then
          call solve_circulant(self, self%eig_kit_mass_0, in, out)
       elseif  (degree == 1) then
          call solve_circulant(self, self%eig_kit_mass_1, in, out)
       elseif (degree == 2 ) then
          call self%multiply_mass( in , self%work, self%s_deg_0 )
          do i=1,self%n_cells-1
             self%wsave(i+1) = 0.5_f64 * ( self%work(i) + self%work(i+1) )
          end do
          self%wsave(1) = 0.5_f64 * ( self%work(1) + self%work(self%n_cells) )
          call self%multiply_interpolation_inverse_transpose( self%wsave, out )
       elseif (degree == 3) then
          call self%multiply_mass( in , self%work, self%s_deg_0-1 )
          do i=1,self%n_cells-1
             self%wsave(i+1) = 0.5_f64 * ( self%work(i) + self%work(i+1) )
          end do
          self%wsave(1) = 0.5_f64 * ( self%work(1) + self%work(self%n_cells) )
          call self%multiply_interpolation_inverse_transpose( self%wsave, out )
       else
          print*, 'degree ', degree, 'not availlable in maxwell_1d_fem object' 
       endif
    else
       out = in
    end if

  end subroutine transform_dofs_1d_fem


end module sll_m_maxwell_1d_fem
