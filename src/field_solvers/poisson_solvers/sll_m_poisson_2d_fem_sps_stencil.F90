module sll_m_poisson_2d_fem_sps_stencil
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_bsplines, only: sll_c_bsplines

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_polar_bsplines_2d, only: sll_t_polar_bsplines_2d

  use sll_m_singular_mapping_discrete, only: sll_t_singular_mapping_discrete

  use sll_m_poisson_2d_fem_sps_weak_form, only: sll_t_poisson_2d_fem_sps_weak_form

  use sll_m_poisson_2d_fem_sps_stencil_assembler, only: sll_t_poisson_2d_fem_sps_stencil_assembler

  use sll_m_poisson_2d_fem_sps_stencil_projector, only: sll_t_poisson_2d_fem_sps_stencil_projector

  use sll_m_vector_space_real_array_2d, only: sll_t_vector_space_real_array_2d

  use sll_m_vector_space_c1_block, only: sll_t_vector_space_c1_block

  use sll_m_linear_operator_matrix_stencil_to_stencil, only: sll_t_linear_operator_matrix_stencil_to_stencil

  use sll_m_linear_operator_matrix_c1_block, only: sll_t_linear_operator_matrix_c1_block

  use sll_m_conjugate_gradient, only: sll_t_conjugate_gradient

  use sll_m_gauss_legendre_integration, only: &
    sll_f_gauss_legendre_points, &
    sll_f_gauss_legendre_weights

  use sll_m_boundary_condition_descriptors, only: sll_p_dirichlet

  implicit none

  public :: sll_t_poisson_2d_fem_sps_stencil

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  ! Abstract interface for right hand side
  abstract interface
    SLL_PURE function i_fun_rhs( x ) result( rhs )
      import wp
      real(wp), intent(in) :: x(2)
      real(wp) :: rhs
    end function i_fun_rhs
  end interface

  type :: sll_t_poisson_2d_fem_sps_stencil

    ! To initialize B-splines
    integer :: mm, n1, n2, ncells1, ncells2, p1, p2

    ! Analytical and discrete mappings
    type(sll_t_singular_mapping_discrete), pointer :: mapping

    ! Number of finite elements
    integer :: Nk1, Nk2

    ! Number of Gauss-Legendre quadrature points
    integer :: Nq1, Nq2

    ! Quadrature points
    real(wp), allocatable :: quad_points_eta1(:,:)
    real(wp), allocatable :: quad_points_eta2(:,:)

    ! 1D data (B-splines values and derivatives)
    real(wp), allocatable :: data_1d_eta1(:,:,:,:)
    real(wp), allocatable :: data_1d_eta2(:,:,:,:)

    ! 2D data (inverse metric tensor, integral volume, right hand side)
    real(wp), allocatable :: int_volume(:,:,:,:)
    real(wp), allocatable :: inv_metric(:,:,:,:,:,:)
    real(wp), allocatable :: data_2d_rhs(:,:,:,:)

    ! Stiffness and mass matrices and C1 projections
    type(sll_t_linear_operator_matrix_stencil_to_stencil) :: A_linop_stencil
    type(sll_t_linear_operator_matrix_stencil_to_stencil) :: M_linop_stencil

    ! Matrix for barycentric coordinates
    real(wp), allocatable :: L(:,:)

    ! Right hand side vector and C1 projection
    real(wp), allocatable :: b(:)
    type(sll_t_vector_space_real_array_2d) :: bs_tmp1
    type(sll_t_vector_space_real_array_2d) :: bs_tmp2

    ! Solution and C1 projection
    real(wp), allocatable :: x(:)

    ! Weak form
    type(sll_t_poisson_2d_fem_sps_weak_form) :: weak_form

    ! Assembler
    type(sll_t_poisson_2d_fem_sps_stencil_assembler) :: assembler

    ! C1 projector
    type(sll_t_poisson_2d_fem_sps_stencil_projector) :: projector

    ! Linear solver
    type(sll_t_conjugate_gradient) :: cjsolver
    real(wp) :: tol = 1.0e-14_wp  ! default value, can be overwritten from init method
    logical  :: verbose = .false. ! default value, can be overwritten from init method

    ! New C1 block format
    type(sll_t_linear_operator_matrix_c1_block) :: Ap_linop_c1_block
    type(sll_t_linear_operator_matrix_c1_block) :: Mp_linop_c1_block
    type(sll_t_vector_space_c1_block) :: xp_vecsp_c1_block
    type(sll_t_vector_space_c1_block) :: bp_vecsp_c1_block

  contains

    procedure :: init                    => s_poisson_2d_fem_sps_stencil__init
    procedure :: set_boundary_conditions => s_poisson_2d_fem_sps_stencil__set_boundary_conditions
    procedure :: free                    => s_poisson_2d_fem_sps_stencil__free

    ! Generic procedure with multiple implementations
    procedure :: solve_1 => s_poisson_2d_fem_sps_stencil__solve_1
    procedure :: solve_2 => s_poisson_2d_fem_sps_stencil__solve_2
    generic   :: solve   => solve_1, solve_2

  end type sll_t_poisson_2d_fem_sps_stencil

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initializer
  subroutine s_poisson_2d_fem_sps_stencil__init( &
    self         , &
    bsplines_eta1, &
    bsplines_eta2, &
    breaks_eta1  , &
    breaks_eta2  , &
    mapping      , &
    tol          , &
    verbose )
    class(sll_t_poisson_2d_fem_sps_stencil)      , intent(inout) :: self
    class(sll_c_bsplines)                , target, intent(in   ) :: bsplines_eta1
    class(sll_c_bsplines)                , target, intent(in   ) :: bsplines_eta2
    real(wp), allocatable                        , intent(in   ) :: breaks_eta1(:)
    real(wp), allocatable                        , intent(in   ) :: breaks_eta2(:)
    type(sll_t_singular_mapping_discrete), target, intent(in   ) :: mapping
    real(wp), optional                           , intent(in   ) :: tol
    logical , optional                           , intent(in   ) :: verbose

    ! Polar B-splines
    type(sll_t_polar_bsplines_2d) :: polar_bsplines

    ! Quadrature weights
    real(wp), allocatable :: quad_weights_eta1(:,:)
    real(wp), allocatable :: quad_weights_eta2(:,:)

    ! Auxiliary variables
    integer  :: i, j, k, k1, k2, q1, q2
    real(wp) :: cx, cy, jdet, eta(2), jmat(2,2)

    ! New auxiliary variables
    integer :: nb, m1(4), m2(4)

    if ( present( tol     ) ) self % tol     = tol
    if ( present( verbose ) ) self % verbose = verbose

    ! Smooth spline mapping
    self % mapping => mapping

    ! Polar B-splines
    call polar_bsplines % init( bsplines_eta1, bsplines_eta2, mapping )

    ! Number of degrees of freedom and spline degrees
    self % n1 = bsplines_eta1 % nbasis
    self % n2 = bsplines_eta2 % nbasis

    !---------------------------------------------------------------------------
    ! Allocations
    !---------------------------------------------------------------------------

    ! Number of finite elements
    self % Nk1 = size( breaks_eta1 )-1
    self % Nk2 = size( breaks_eta2 )-1

    ! Spline degrees
    self % p1 = bsplines_eta1 % degree
    self % p2 = bsplines_eta2 % degree

    ! Number of Gauss-Legendre quadrature points
    self % Nq1 = 1 + self % p1
    self % Nq2 = 1 + self % p2

    associate( n1  => self % n1            , &
               n2  => self % n2            , &
               p1  => self % p1            , &
               p2  => self % p2            , &
               nn  => 3+(self%n1-2)*self%n2, &
               Nk1 => self % Nk1           , &
               Nk2 => self % Nk2           , &
               Nq1 => self % Nq1           , &
               Nq2 => self % Nq2 )

      ! Quadrature points
      allocate( self % quad_points_eta1( Nq1, Nk1 ) )
      allocate( self % quad_points_eta2( Nq2, Nk2 ) )

      ! Quadrature weights
      allocate( quad_weights_eta1( Nq1, Nk1 ) )
      allocate( quad_weights_eta2( Nq2, Nk2 ) )

      ! 1D data
      allocate( self % data_1d_eta1( Nq1, 2, 1+p1, Nk1 ) )
      allocate( self % data_1d_eta2( Nq2, 2, 1+p2, Nk2 ) )

      ! 2D data
      allocate( self % int_volume ( Nq1, Nq2, Nk1, Nk2 ) )
      allocate( self % inv_metric ( Nq1, Nq2, Nk1, Nk2, 2, 2 ) )
      allocate( self % data_2d_rhs( Nq1, Nq2, Nk1, Nk2 ) )

      ! Barycentric coordinates
      allocate( self % L( 2*n2, 3 ) )

      ! Right hand side
      allocate( self % b( n1*n2 ) )
      allocate( self % bs_tmp1 % array( 1-p1:n1+p1, 1-p2:n2+p2 ) )
      allocate( self % bs_tmp2 % array( 1-p1:n1+p1, 1-p2:n2+p2 ) )

      ! Solution
      allocate( self % x( n1*n2 ) )

      !-------------------------------------------------------------------------
      ! Initialize quadrature points and weights
      !-------------------------------------------------------------------------

      ! Quadrature points and weights along s
      do k1 = 1, Nk1
        self % quad_points_eta1 (:,k1) = sll_f_gauss_legendre_points ( Nq1, breaks_eta1(k1), breaks_eta1(k1+1) )
        quad_weights_eta1(:,k1) = sll_f_gauss_legendre_weights( Nq1, breaks_eta1(k1), breaks_eta1(k1+1) )
      end do

      ! Quadrature points and weights along theta
      do k2 = 1, Nk2
        self % quad_points_eta2 (:,k2) = sll_f_gauss_legendre_points ( Nq2, breaks_eta2(k2), breaks_eta2(k2+1) )
        quad_weights_eta2(:,k2) = sll_f_gauss_legendre_weights( Nq2, breaks_eta2(k2), breaks_eta2(k2+1) )
      end do

      !-------------------------------------------------------------------------
      ! Fill in 1D data
      !-------------------------------------------------------------------------

      ! 1D data along s
      do k1 = 1, Nk1
        do q1 = 1, Nq1
          call bsplines_eta1 % eval_basis_and_n_derivs( &
            x      = self % quad_points_eta1(q1,k1), &
            n      = 1, &
            derivs = self % data_1d_eta1(q1,:,:,k1), &
            jmin   = i )
        end do
      end do

      ! 1D data along theta
      do k2 = 1, Nk2
        do q2 = 1, Nq2
          call bsplines_eta2 % eval_basis_and_n_derivs( &
            x      = self % quad_points_eta2(q2,k2), &
            n      = 1, &
            derivs = self % data_1d_eta2(q2,:,:,k2), &
            jmin   = i )
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

              eta  = [ self % quad_points_eta1(q1,k1), self % quad_points_eta2(q2,k2) ]

              jdet = mapping % jdet( eta )
              jmat = mapping % jmat( eta )

              ! Area associated to each quadrature point
              self % int_volume(q1,q2,k1,k2) = abs( jdet ) * quad_weights_eta1(q1,k1) * quad_weights_eta2(q2,k2)

              ! Inverse metric tensor
              self % inv_metric(q1,q2,k1,k2,1,1) = jmat(1,2)**2 + jmat(2,2)**2
              self % inv_metric(q1,q2,k1,k2,1,2) = - jmat(1,1)*jmat(1,2) - jmat(2,1)*jmat(2,2)
              self % inv_metric(q1,q2,k1,k2,2,1) = self % inv_metric(q1,q2,k1,k2,1,2) ! symmetry
              self % inv_metric(q1,q2,k1,k2,2,2) = jmat(1,1)**2 + jmat(2,1)**2
              self % inv_metric(q1,q2,k1,k2,:,:) = self % inv_metric(q1,q2,k1,k2,:,:) / jdet**2

            end do
          end do
        end do
      end do

      !-------------------------------------------------------------------------
      ! Matrix of barycentric coordinates
      !-------------------------------------------------------------------------

      do i = 1, 2*n2

        ! Indices to access control points
        j = (i-1) / n2 + 1
        k = modulo( (i-1), n2 ) + 1

        cx = mapping % spline_2d_x1 % bcoef(j,k)
        cy = mapping % spline_2d_x2 % bcoef(j,k)

        self % L(i,1) = polar_bsplines % eval_l0( cx, cy )
        self % L(i,2) = polar_bsplines % eval_l1( cx, cy )
        self % L(i,3) = polar_bsplines % eval_l2( cx, cy )

      end do

      !-------------------------------------------------------------------------
      ! Initialize assembler and projector
      !-------------------------------------------------------------------------

      call self % assembler % init( n1, n2, p1, p2, self % weak_form )
      call self % projector % init( n1, n2, p1, p2, self % L )

      !-------------------------------------------------------------------------
      ! Fill in stiffness and mass matrices
      !-------------------------------------------------------------------------

      ! Construct stencil linear operators
      call self % A_linop_stencil % init( n1, n2, p1, p2 )
      call self % M_linop_stencil % init( n1, n2, p1, p2 )

      ! Cycle over finite elements
      do k2 = 1, Nk2
        do k1 = 1, Nk1
          call self % assembler % add_element_mat( &
            k1                        , &
            k2                        , &
            self % data_1d_eta1       , &
            self % data_1d_eta2       , &
            self % int_volume         , &
            self % inv_metric         , &
            self % A_linop_stencil % A, &
            self % M_linop_stencil % A)
        end do
      end do

      !-------------------------------------------------------------------------
      ! Initialize linear system (homogeneous Dirichlet boundary conditions)
      !-------------------------------------------------------------------------

      ! Construct C1 block linear operators

      nb = (n1-2)*n2
      m1 = (/ 3, 3 , nb, (n1-2) /)
      m2 = (/ 3, nb, 3 , n2 /)

      call self % Ap_linop_c1_block % init( m1, m2, p1, p2 )
      call self % Mp_linop_c1_block % init( m1, m2, p1, p2 )

      !-------------------------------------------------------------------------
      ! Compute C1 projections of stiffness and mass matrices
      !-------------------------------------------------------------------------

      call self % projector % change_basis_matrix( self % A_linop_stencil, self % Ap_linop_c1_block )
      call self % projector % change_basis_matrix( self % M_linop_stencil, self % Mp_linop_c1_block )

      ! Construct C1 block vector space from vector xp

      call self % xp_vecsp_c1_block % init( n1-2, n2, p1, p2 )

      allocate( self % xp_vecsp_c1_block % vd % array( 1:3 ) )
      allocate( self % xp_vecsp_c1_block % vs % array( 1-p1:(n1-2)+p1, 1-p2:n2+p2 ) )
      self % xp_vecsp_c1_block % vd % array(:  ) = 0.0_wp
      self % xp_vecsp_c1_block % vs % array(:,:) = 0.0_wp

      ! Initialize conjugate gradient solver
      call self % cjsolver % init( &
        tol             = self % tol    , &
        verbose         = self % verbose, &
        template_vector = self % xp_vecsp_c1_block )

    end associate

    call polar_bsplines % free()

  end subroutine s_poisson_2d_fem_sps_stencil__init

  ! Set boundary conditions
  subroutine s_poisson_2d_fem_sps_stencil__set_boundary_conditions( self, bc )
    class(sll_t_poisson_2d_fem_sps_stencil), intent(inout) :: self
    integer                                , intent(in   ) :: bc

    integer :: i, j, i1, i2, j1, j2, k1, k2, nb

    character(len=*), parameter :: this_sub_name = "sll_t_poisson_2d_fem_sps_stencil % set_boundary_conditions"

    if ( bc == sll_p_dirichlet ) then

      associate( n1  => self % n1, &
                 n2  => self % n2, &
                 p1  => self % p1, &
                 p2  => self % p2 )

        nb = (n1-2) * n2

        ! Homogeneous Dirichlet boundary conditions
        self % Ap_linop_c1_block % block2 % A(:,nb-n2+1:nb) = 0.0_wp
        self % Ap_linop_c1_block % block3 % A(nb-n2+1:nb,:) = 0.0_wp
        do i2 = 1, n2
          do i1 = 1, n1-2
            do k2 = -p2, p2
              do k1 = -p1, p1
                j1 = modulo( i1 - 1 + k1, n1-2 ) + 1
                j2 = modulo( i2 - 1 + k2, n2   ) + 1
                i  = (i1-1) * n2 + i2
                j  = (j1-1) * n2 + j2
                ! Check range of indices i and j
                if ( i > nb - n2 .or. j > nb - n2 ) &
                  self % Ap_linop_c1_block % block4 % A(k1,k2,i1,i2) = 0.0_wp
              end do
            end do
          end do
        end do


      end associate

    else
      SLL_ERROR( this_sub_name, "Boundary conditions not implemented" )

    end if

  end subroutine s_poisson_2d_fem_sps_stencil__set_boundary_conditions

  ! Solver: signature #1
  subroutine s_poisson_2d_fem_sps_stencil__solve_1( self, rhs, sol )
    class(sll_t_poisson_2d_fem_sps_stencil), intent(inout) :: self
    procedure(i_fun_rhs)                                   :: rhs
    type(sll_t_spline_2d)                  , intent(inout) :: sol

    ! Auxiliary variables
    integer  :: k1, k2, q1, q2, j2
    real(wp) :: eta(2), x(2)

    associate( n1  => self % n1 , &
               n2  => self % n2 , &
               p1  => self % p1 , &
               p2  => self % p2 , &
               Nk1 => self % Nk1, &
               Nk2 => self % Nk2, &
               Nq1 => self % Nq1, &
               Nq2 => self % Nq2 )

      ! 2D discrete RHS
      do k2 = 1, Nk2
        do k1 = 1, Nk1
          do q2 = 1, Nq2
            do q1 = 1, Nq1

              eta = [ self % quad_points_eta1(q1,k1), self % quad_points_eta2(q2,k2) ]

              x   = self % mapping % eval( eta )

              self % data_2d_rhs(q1,q2,k1,k2) = rhs( x )

            end do
          end do
        end do
      end do

      !-------------------------------------------------------------------------
      ! Fill in right hand side
      !-------------------------------------------------------------------------

      self % b = 0.0_wp

      ! Cycle over finite elements
      do k2 = 1, Nk2
        do k1 = 1, Nk1
          call self % assembler % add_element_rhs( &
            k1                 , &
            k2                 , &
            self % data_1d_eta1, &
            self % data_1d_eta2, &
            self % data_2d_rhs , &
            self % int_volume  , &
            self % b )
        end do
      end do

      !-------------------------------------------------------------------------
      ! Solve linear system (homogeneous Dirichlet boundary conditions)
      !-------------------------------------------------------------------------

      ! Construct C1 block vector space
      call self % bp_vecsp_c1_block % init( n1-2, n2, p1, p2 )
      allocate( self % bp_vecsp_c1_block % vd % array( 1:3 ) )
      allocate( self % bp_vecsp_c1_block % vs % array( 1-p1:(n1-2)+p1, 1-p2:n2+p2 ) )

      ! Compute C1 projection of right hand side
      call self % projector % change_basis_vector( self % b, self % bp_vecsp_c1_block )

      ! Homogeneous Dirichlet boundary conditions
      do j2 = 1, n2
        self % bp_vecsp_c1_block % vs % array(n1-2,j2) = 0.0_wp
      end do

      ! Update buffer regions
      self % bp_vecsp_c1_block % vs % array(1-p1:0            ,:) = 0.0_wp ! no periodicity
      self % bp_vecsp_c1_block % vs % array((n1-2)+1:(n1-2)+p1,:) = 0.0_wp ! no periodicity
      self % bp_vecsp_c1_block % vs % array(:,1-p2:0    ) = self % bp_vecsp_c1_block % vs % array(:,n2-p2+1:n2)
      self % bp_vecsp_c1_block % vs % array(:,n2+1:n2+p2) = self % bp_vecsp_c1_block % vs % array(:,1:p2      )

      ! Solve linear system Ap*xp=bp
      call self % cjsolver % solve( &
        A = self % Ap_linop_c1_block, &
        b = self % bp_vecsp_c1_block, &
        x = self % xp_vecsp_c1_block )

      ! Compute solution in tensor-product space
      call self % projector % change_basis_vector_inv( self % xp_vecsp_c1_block, self % x )

      sol % bcoef = reshape( self % x, (/ n2, n1 /) )

    end associate

  end subroutine s_poisson_2d_fem_sps_stencil__solve_1

  ! Solver: signature #2
  subroutine s_poisson_2d_fem_sps_stencil__solve_2( self, rhs, sol )
    class(sll_t_poisson_2d_fem_sps_stencil), intent(inout) :: self
    type(sll_t_spline_2d)                  , intent(in   ) :: rhs
    type(sll_t_spline_2d)                  , intent(inout) :: sol

    integer :: i, i1, i2, j2

    associate( n1 => self % n1, &
               n2 => self % n2, &
               p1 => self % p1, &
               p2 => self % p2 )

      ! Store temporarily spline coefficients of right hand side into self % b
      self % bs_tmp1 % array(:,:) = 0.0_wp
      self % bs_tmp1 % array(1:n1,1:n2) = rhs % bcoef(1:n1,1:n2)

      call self % M_linop_stencil % dot( self % bs_tmp1, self % bs_tmp2 )

      do i2 = 1, n2
        do i1 = 1, n1
          i = (i1-1) * n2 + i2
          self % b(i) = self % bs_tmp2 % array(i1,i2)
        end do
      end do

      !-------------------------------------------------------------------------
      ! Solve linear system (homogeneous Dirichlet boundary conditions)
      !-------------------------------------------------------------------------

      ! Construct C1 block vector space
      call self % bp_vecsp_c1_block % init( n1-2, n2, p1, p2 )
      allocate( self % bp_vecsp_c1_block % vd % array( 1:3 ) )
      allocate( self % bp_vecsp_c1_block % vs % array( 1-p1:(n1-2)+p1, 1-p2:n2+p2 ) )

      ! Compute C1 projection of right hand side
      call self % projector % change_basis_vector( self % b, self % bp_vecsp_c1_block )

      ! Homogeneous Dirichlet boundary conditions
      do j2 = 1, n2
        self % bp_vecsp_c1_block % vs % array(n1-2,j2) = 0.0_wp
      end do

      ! Update buffer regions
      self % bp_vecsp_c1_block % vs % array(1-p1:0            ,:) = 0.0_wp ! no periodicity
      self % bp_vecsp_c1_block % vs % array((n1-2)+1:(n1-2)+p1,:) = 0.0_wp ! no periodicity
      self % bp_vecsp_c1_block % vs % array(:,1-p2:0    ) = self % bp_vecsp_c1_block % vs % array(:,n2-p2+1:n2)
      self % bp_vecsp_c1_block % vs % array(:,n2+1:n2+p2) = self % bp_vecsp_c1_block % vs % array(:,1:p2      )

      ! Solve linear system Ap*xp=bp
      call self % cjsolver % solve( &
        A = self % Ap_linop_c1_block, &
        b = self % bp_vecsp_c1_block, &
        x = self % xp_vecsp_c1_block )

      ! Compute solution in tensor-product space
      call self % projector % change_basis_vector_inv( self % xp_vecsp_c1_block, self % x )

      sol % bcoef = reshape( self % x, (/ n2, n1 /) )

    end associate

  end subroutine s_poisson_2d_fem_sps_stencil__solve_2

  ! Free
  subroutine s_poisson_2d_fem_sps_stencil__free( self )
    class(sll_t_poisson_2d_fem_sps_stencil), intent(inout) :: self

    deallocate( self % quad_points_eta1 )
    deallocate( self % quad_points_eta2 )
    deallocate( self % data_1d_eta1 )
    deallocate( self % data_1d_eta2 )
    deallocate( self % int_volume )
    deallocate( self % inv_metric )
    deallocate( self % L  )
    deallocate( self % x  )
    deallocate( self % b  )
    deallocate( self % bs_tmp1 % array )
    deallocate( self % bs_tmp2 % array )

    call self % Ap_linop_c1_block % free()
    call self % Mp_linop_c1_block % free()
    deallocate( self % bp_vecsp_c1_block % vd % array )
    deallocate( self % bp_vecsp_c1_block % vs % array )
    deallocate( self % xp_vecsp_c1_block % vd % array )
    deallocate( self % xp_vecsp_c1_block % vs % array )

    call self % projector % free()

    call self % cjsolver % free()

  end subroutine s_poisson_2d_fem_sps_stencil__free

end module sll_m_poisson_2d_fem_sps_stencil
