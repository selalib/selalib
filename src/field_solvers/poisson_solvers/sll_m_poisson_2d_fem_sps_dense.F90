module sll_m_poisson_2d_fem_sps_dense
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

   use sll_m_working_precision, only: f64

   use sll_m_bsplines, only: sll_c_bsplines

   use sll_m_spline_2d, only: sll_t_spline_2d

   use sll_m_polar_bsplines_2d, only: sll_t_polar_bsplines_2d

   use sll_m_singular_mapping_discrete, only: sll_t_singular_mapping_discrete

   use sll_m_poisson_2d_fem_sps_weak_form, only: sll_t_poisson_2d_fem_sps_weak_form

   use sll_m_poisson_2d_fem_sps_dense_assembler, only: sll_t_poisson_2d_fem_sps_dense_assembler

   use sll_m_poisson_2d_fem_sps_dense_projector, only: sll_t_poisson_2d_fem_sps_dense_projector

   use sll_m_vector_space_real_array_1d, only: sll_t_vector_space_real_array_1d

   use sll_m_linear_operator_matrix_dense_to_dense, only: sll_t_linear_operator_matrix_dense_to_dense

   use sll_m_conjugate_gradient, only: sll_t_conjugate_gradient

   use sll_m_gauss_legendre_integration, only: &
      sll_f_gauss_legendre_points, &
      sll_f_gauss_legendre_weights

   implicit none

   public :: sll_t_poisson_2d_fem_sps_dense

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Working precision
   integer, parameter :: wp = f64

   ! Abstract interface for right hand side
   abstract interface
      SLL_PURE function i_fun_rhs(x) result(rhs)
         import wp
         real(wp), intent(in) :: x(2)
         real(wp) :: rhs
      end function i_fun_rhs
   end interface

   type :: sll_t_poisson_2d_fem_sps_dense

      ! To initialize B-splines (p1,p2 degrees)
      integer :: mm, n1, n2, ncells1, ncells2

      ! Analytical and discrete mappings
      type(sll_t_singular_mapping_discrete), pointer :: mapping

      ! Number of finite elements
      integer :: Nk1, Nk2

      ! Number of Gauss-Legendre quadrature points
      integer :: Nq1, Nq2

      ! Quadrature points
      real(wp), allocatable :: quad_points_eta1(:, :)
      real(wp), allocatable :: quad_points_eta2(:, :)

      ! 1D data (B-splines values and derivatives)
      real(wp), allocatable :: data_1d_eta1(:, :, :, :)
      real(wp), allocatable :: data_1d_eta2(:, :, :, :)

      ! 2D data (inverse metric tensor, integral volume, right hand side)
      real(wp), allocatable :: int_volume(:, :, :, :)
      real(wp), allocatable :: inv_metric(:, :, :, :, :, :)
      real(wp), allocatable :: data_2d_rhs(:, :, :, :)

      ! Stiffness and mass matrices and C1 projections
      real(wp), allocatable :: A(:, :)
      real(wp), allocatable :: M(:, :)
      real(wp), allocatable :: Ap(:, :)
      real(wp), allocatable :: Mp(:, :)

      ! Matrix for barycentric coordinates
      real(wp), allocatable :: L(:, :)

      ! Right hand side vector and C1 projection
      real(wp), allocatable :: b(:)
      real(wp), allocatable :: bp(:)

      ! Solution and C1 projection
      real(wp), allocatable :: x(:)
      real(wp), allocatable :: xp(:)

      ! Weak form
      type(sll_t_poisson_2d_fem_sps_weak_form) :: weak_form

      ! Assembler
      type(sll_t_poisson_2d_fem_sps_dense_assembler) :: assembler

      ! C1 projector
      type(sll_t_poisson_2d_fem_sps_dense_projector) :: projector

      ! Linear solver
      type(sll_t_vector_space_real_array_1d)            :: bp_vecsp
      type(sll_t_vector_space_real_array_1d)            :: xp_vecsp
      type(sll_t_linear_operator_matrix_dense_to_dense) :: Ap_linop
      type(sll_t_conjugate_gradient)                    :: cjsolver
      real(wp) :: tol = 1.0e-14_wp  ! default value, can be overwritten from init method
      logical  :: verbose = .false. ! default value, can be overwritten from init method

   contains

      procedure :: init => s_poisson_2d_fem_sps_dense__init
      procedure :: free => s_poisson_2d_fem_sps_dense__free

      ! Generic procedure with multiple implementations
      procedure :: solve_1 => s_poisson_2d_fem_sps_dense__solve_1
      procedure :: solve_2 => s_poisson_2d_fem_sps_dense__solve_2
      generic   :: solve => solve_1, solve_2

   end type sll_t_poisson_2d_fem_sps_dense

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Initializer
   subroutine s_poisson_2d_fem_sps_dense__init( &
      self, &
      bsplines_eta1, &
      bsplines_eta2, &
      breaks_eta1, &
      breaks_eta2, &
      mapping, &
      tol, &
      verbose)
      class(sll_t_poisson_2d_fem_sps_dense), intent(inout) :: self
      class(sll_c_bsplines), target, intent(in) :: bsplines_eta1
      class(sll_c_bsplines), target, intent(in) :: bsplines_eta2
      real(wp), allocatable, intent(in) :: breaks_eta1(:)
      real(wp), allocatable, intent(in) :: breaks_eta2(:)
      type(sll_t_singular_mapping_discrete), target, intent(in) :: mapping
      real(wp), optional, intent(in) :: tol
      logical, optional, intent(in) :: verbose

      ! Polar B-splines
      type(sll_t_polar_bsplines_2d) :: polar_bsplines

      ! Quadrature weights
      real(wp), allocatable :: quad_weights_eta1(:, :)
      real(wp), allocatable :: quad_weights_eta2(:, :)

      ! Auxiliary variables
      integer  :: i, j, k, idx, k1, k2, q1, q2, p1, p2
      real(wp) :: cx, cy, jdet, eta(2), jmat(2, 2)

      if (present(tol)) self%tol = tol
      if (present(verbose)) self%verbose = verbose

      ! Smooth spline mapping
      self%mapping => mapping

      ! Polar B-splines
      call polar_bsplines%init(bsplines_eta1, bsplines_eta2, mapping)

      ! Number of degrees of freedom and spline degrees
      self%n1 = bsplines_eta1%nbasis
      self%n2 = bsplines_eta2%nbasis

      !---------------------------------------------------------------------------
      ! Allocations
      !---------------------------------------------------------------------------

      ! Number of finite elements
      self%Nk1 = size(breaks_eta1) - 1
      self%Nk2 = size(breaks_eta2) - 1

      ! Spline degrees
      p1 = bsplines_eta1%degree
      p2 = bsplines_eta2%degree

      ! Number of Gauss-Legendre quadrature points
      self%Nq1 = 1 + p1
      self%Nq2 = 1 + p2

      associate (n1 => self%n1, &
                 n2 => self%n2, &
                 nn => 3 + (self%n1 - 2)*self%n2, &
                 Nk1 => self%Nk1, &
                 Nk2 => self%Nk2, &
                 Nq1 => self%Nq1, &
                 Nq2 => self%Nq2)

         ! Quadrature points
         allocate (self%quad_points_eta1(Nq1, Nk1))
         allocate (self%quad_points_eta2(Nq2, Nk2))

         ! Quadrature weights
         allocate (quad_weights_eta1(Nq1, Nk1))
         allocate (quad_weights_eta2(Nq2, Nk2))

         ! 1D data
         allocate (self%data_1d_eta1(Nq1, 2, 1 + p1, Nk1))
         allocate (self%data_1d_eta2(Nq2, 2, 1 + p2, Nk2))

         ! 2D data
         allocate (self%int_volume(Nq1, Nq2, Nk1, Nk2))
         allocate (self%inv_metric(Nq1, Nq2, Nk1, Nk2, 2, 2))
         allocate (self%data_2d_rhs(Nq1, Nq2, Nk1, Nk2))

         ! Barycentric coordinates
         allocate (self%L(2*n2, 3))

         ! Stiffness and mass matrices
         allocate (self%A(n1*n2, n1*n2))
         allocate (self%M(n1*n2, n1*n2))

         ! Right hand side
         allocate (self%b(n1*n2))

         ! Solution
         allocate (self%x(n1*n2))

         ! C1 projections of stiffness and mass matrices
         allocate (self%Ap(nn, nn))
         allocate (self%Mp(nn, nn))

         ! C1 projection of right hand side
         allocate (self%bp(nn))

         ! C1 projection of solution
         allocate (self%xp(nn))

         !-------------------------------------------------------------------------
         ! Initialize quadrature points and weights
         !-------------------------------------------------------------------------

         ! Quadrature points and weights along s
         do k1 = 1, Nk1
            self%quad_points_eta1(:, k1) = sll_f_gauss_legendre_points(Nq1, breaks_eta1(k1), breaks_eta1(k1 + 1))
            quad_weights_eta1(:, k1) = sll_f_gauss_legendre_weights(Nq1, breaks_eta1(k1), breaks_eta1(k1 + 1))
         end do

         ! Quadrature points and weights along theta
         do k2 = 1, Nk2
            self%quad_points_eta2(:, k2) = sll_f_gauss_legendre_points(Nq2, breaks_eta2(k2), breaks_eta2(k2 + 1))
            quad_weights_eta2(:, k2) = sll_f_gauss_legendre_weights(Nq2, breaks_eta2(k2), breaks_eta2(k2 + 1))
         end do

         !-------------------------------------------------------------------------
         ! Fill in 1D data
         !-------------------------------------------------------------------------

         ! 1D data along s
         do k1 = 1, Nk1
            do q1 = 1, Nq1
               call bsplines_eta1%eval_basis_and_n_derivs( &
                  x=self%quad_points_eta1(q1, k1), &
                  n=1, &
                  derivs=self%data_1d_eta1(q1, :, :, k1), &
                  jmin=i)
            end do
         end do

         ! 1D data along theta
         do k2 = 1, Nk2
            do q2 = 1, Nq2
               call bsplines_eta2%eval_basis_and_n_derivs( &
                  x=self%quad_points_eta2(q2, k2), &
                  n=1, &
                  derivs=self%data_1d_eta2(q2, :, :, k2), &
                  jmin=i)
            end do
         end do

         !-------------------------------------------------------------------------
         ! Fill in 2D data
         !-------------------------------------------------------------------------

         ! 2D data
         do k2 = 1, Nk2
            do k1 = 1, Nk1
               do q2 = 1, Nq2
                  do q1 = 1, Nq1

                     eta = [self%quad_points_eta1(q1, k1), self%quad_points_eta2(q2, k2)]

                     jdet = mapping%jdet(eta)
                     jmat = mapping%jmat(eta)

                     ! Area associated to each quadrature point
                     self%int_volume(q1, q2, k1, k2) = abs(jdet)*quad_weights_eta1(q1, k1)*quad_weights_eta2(q2, k2)

                     ! Inverse metric tensor
                     self%inv_metric(q1, q2, k1, k2, 1, 1) = jmat(1, 2)**2 + jmat(2, 2)**2
                     self%inv_metric(q1, q2, k1, k2, 1, 2) = -jmat(1, 1)*jmat(1, 2) - jmat(2, 1)*jmat(2, 2)
                     self%inv_metric(q1, q2, k1, k2, 2, 1) = self%inv_metric(q1, q2, k1, k2, 1, 2) ! symmetry
                     self%inv_metric(q1, q2, k1, k2, 2, 2) = jmat(1, 1)**2 + jmat(2, 1)**2
                     self%inv_metric(q1, q2, k1, k2, :, :) = self%inv_metric(q1, q2, k1, k2, :, :)/jdet**2

                  end do
               end do
            end do
         end do

         !-------------------------------------------------------------------------
         ! Matrix of barycentric coordinates
         !-------------------------------------------------------------------------

         do i = 1, 2*n2

            ! Indices to access control points
            j = (i - 1)/n2 + 1
            k = modulo((i - 1), n2) + 1

            cx = mapping%spline_2d_x1%bcoef(j, k)
            cy = mapping%spline_2d_x2%bcoef(j, k)

            self%L(i, 1) = polar_bsplines%eval_l0(cx, cy)
            self%L(i, 2) = polar_bsplines%eval_l1(cx, cy)
            self%L(i, 3) = polar_bsplines%eval_l2(cx, cy)

         end do

         !-------------------------------------------------------------------------
         ! Initialize assembler and projector
         !-------------------------------------------------------------------------

         call self%assembler%init(n1, n2, self%weak_form)
         call self%projector%init(n1, n2, self%L)

         !-------------------------------------------------------------------------
         ! Fill in stiffness and mass matrices
         !-------------------------------------------------------------------------

         self%A = 0.0_wp
         self%M = 0.0_wp

         ! Cycle over finite elements
         do k2 = 1, Nk2
            do k1 = 1, Nk1
               call self%assembler%add_element_mat( &
                  k1, &
                  k2, &
                  self%data_1d_eta1, &
                  self%data_1d_eta2, &
                  self%int_volume, &
                  self%inv_metric, &
                  self%A, &
                  self%M)
            end do
         end do

         !-------------------------------------------------------------------------
         ! Compute C1 projections of stiffness and mass matrices
         !-------------------------------------------------------------------------

         call self%projector%change_basis_matrix(self%A, self%Ap)
         call self%projector%change_basis_matrix(self%M, self%Mp)

         !-------------------------------------------------------------------------
         ! Initialize linear system (homogeneous Dirichlet boundary conditions)
         !-------------------------------------------------------------------------

         idx = 3 + (n1 - 3)*n2

         ! Construct linear operator from matrix Ap
         call self%Ap_linop%init(idx, idx)
         self%Ap_linop%A = self%Ap(:idx, :idx)

         ! Construct vector space for solution
         self%xp = 0.0_wp
         allocate (self%xp_vecsp%array(size(self%xp(:idx))), source=self%xp(:idx))

         ! Initialize conjugate gradient solver
         call self%cjsolver%init( &
            tol=self%tol, &
            verbose=self%verbose, &
            template_vector=self%xp_vecsp)

      end associate

      call polar_bsplines%free()

   end subroutine s_poisson_2d_fem_sps_dense__init

   ! Solver: signature #1
   subroutine s_poisson_2d_fem_sps_dense__solve_1(self, rhs, sol)
      class(sll_t_poisson_2d_fem_sps_dense), intent(inout) :: self
      procedure(i_fun_rhs)                                 :: rhs
      type(sll_t_spline_2d), intent(inout) :: sol

      ! Auxiliary variables
      integer  :: idx, k1, k2, q1, q2
      real(wp) :: eta(2), x(2)

      associate (n1 => self%n1, &
                 n2 => self%n2, &
                 nn => 3 + (self%n1 - 2)*self%n2, &
                 Nk1 => self%Nk1, &
                 Nk2 => self%Nk2, &
                 Nq1 => self%Nq1, &
                 Nq2 => self%Nq2)

         ! 2D discrete RHS
         do k2 = 1, Nk2
            do k1 = 1, Nk1
               do q2 = 1, Nq2
                  do q1 = 1, Nq1

                     eta = [self%quad_points_eta1(q1, k1), self%quad_points_eta2(q2, k2)]

                     x = self%mapping%eval(eta)

                     self%data_2d_rhs(q1, q2, k1, k2) = rhs(x)

                  end do
               end do
            end do
         end do

         !-------------------------------------------------------------------------
         ! Fill in right hand side
         !-------------------------------------------------------------------------

         self%b = 0.0_wp

         ! Cycle over finite elements
         do k2 = 1, Nk2
            do k1 = 1, Nk1
               call self%assembler%add_element_rhs( &
                  k1, &
                  k2, &
                  self%data_1d_eta1, &
                  self%data_1d_eta2, &
                  self%data_2d_rhs, &
                  self%int_volume, &
                  self%b)
            end do
         end do

         !-------------------------------------------------------------------------
         ! Compute C1 projection of right hand side
         !-------------------------------------------------------------------------

         call self%projector%change_basis_vector(self%b, self%bp)

         !-------------------------------------------------------------------------
         ! Solve linear system (homogeneous Dirichlet boundary conditions)
         !-------------------------------------------------------------------------

         idx = 3 + (n1 - 3)*n2

         ! Construct vector space from vector bp
         allocate (self%bp_vecsp%array(size(self%bp(:idx))), source=self%bp(:idx))

         ! Solve linear system Ap*xp=bp
         call self%cjsolver%solve( &
            A=self%Ap_linop, &
            b=self%bp_vecsp, &
            x=self%xp_vecsp)

         ! Copy solution into array xp
         self%xp(:idx) = self%xp_vecsp%array

         !-------------------------------------------------------------------------
         ! Compute solution in tensor-product space
         !-------------------------------------------------------------------------

         call self%projector%change_basis_vector_inv(self%xp, self%x)

         sol%bcoef = reshape(self%x, (/n2, n1/))

      end associate

   end subroutine s_poisson_2d_fem_sps_dense__solve_1

   ! Solver: signature #2
   subroutine s_poisson_2d_fem_sps_dense__solve_2(self, rhs, sol)
      class(sll_t_poisson_2d_fem_sps_dense), intent(inout) :: self
      type(sll_t_spline_2d), intent(in) :: rhs
      type(sll_t_spline_2d), intent(inout) :: sol

      integer :: idx, i, i1, i2

      associate (n1 => self%n1, n2 => self%n2)

         ! Store temporarily spline coefficients of right hand side into self % b
         do i2 = 1, n2
            do i1 = 1, n1
               i = (i1 - 1)*n2 + i2
               self%b(i) = rhs%bcoef(i1, i2)
            end do
         end do

         ! Compute right hand side
         self%b = matmul(self%M, self%b)

         ! Compute C1 projection of right hand side
         call self%projector%change_basis_vector(self%b, self%bp)

         !-------------------------------------------------------------------------
         ! Solve linear system (homogeneous Dirichlet boundary conditions)
         !-------------------------------------------------------------------------

         idx = 3 + (n1 - 3)*n2

         ! Construct vector space from vector bp
         allocate (self%bp_vecsp%array(size(self%bp(:idx))), source=self%bp(:idx))

         ! Solve linear system Ap*xp=bp
         call self%cjsolver%solve( &
            A=self%Ap_linop, &
            b=self%bp_vecsp, &
            x=self%xp_vecsp)

         ! Copy solution into array xp
         self%xp(:idx) = self%xp_vecsp%array

         !-------------------------------------------------------------------------
         ! Compute solution in tensor-product space
         !-------------------------------------------------------------------------

         call self%projector%change_basis_vector_inv(self%xp, self%x)

         sol%bcoef = reshape(self%x, (/n2, n1/))

      end associate

   end subroutine s_poisson_2d_fem_sps_dense__solve_2

   ! Free
   subroutine s_poisson_2d_fem_sps_dense__free(self)
      class(sll_t_poisson_2d_fem_sps_dense), intent(inout) :: self

      deallocate (self%quad_points_eta1)
      deallocate (self%quad_points_eta2)
      deallocate (self%data_1d_eta1)
      deallocate (self%data_1d_eta2)
      deallocate (self%int_volume)
      deallocate (self%inv_metric)
      deallocate (self%A)
      deallocate (self%M)
      deallocate (self%Ap)
      deallocate (self%Mp)
      deallocate (self%b)
      deallocate (self%bp)
      deallocate (self%L)

      call self%projector%free()

      call self%cjsolver%free()
      call self%Ap_linop%free()
      deallocate (self%bp_vecsp%array)
      deallocate (self%xp_vecsp%array)

   end subroutine s_poisson_2d_fem_sps_dense__free

end module sll_m_poisson_2d_fem_sps_dense
