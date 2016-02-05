!**************************************************************
!  Author: Jakob Ameres, jakob.ameres@tum.de
!**************************************************************

!> Module to solve Poisson equation on one dimensional mesh using Finite Elements
module sll_m_poisson_1d_fem
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

! use F77_lapack, only: &
!   dgbsv

  use sll_m_arbitrary_degree_splines, only: &
    sll_t_arbitrary_degree_spline_1d, &
    sll_f_spline_derivatives_at_x, &
    sll_f_splines_at_x, &
    sll_f_new_arbitrary_degree_spline_1d, &
    sll_p_periodic_arbitrary_deg_spline

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_periodic

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_1d, &
    sll_o_cell, &
    sll_o_cell_margin, &
    sll_o_mesh_area

  use sll_m_fft, only: &
    sll_s_fft_exec_c2r_1d, &
    sll_s_fft_exec_r2c_1d, &
    sll_s_fft_free, &
    sll_s_fft_init_c2r_1d, &
    sll_s_fft_init_r2c_1d, &
    sll_t_fft

  use sll_m_gauss_legendre_integration, only: &
    sll_f_gauss_legendre_points_and_weights

  use sll_m_utilities, only: &
    sll_f_is_power_of_two

  implicit none

  public :: &
    sll_f_new_poisson_1d_fem, &
    sll_t_poisson_1d_fem, &
    sll_i_poisson_1d_fem_rhs_function

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  private
    !  public :: initialize, new, solve

    !  !> Solver data structure
    !  type :: poisson_1d_periodic
    !     sll_int32                         :: nc_eta1 !< number of cells
    !     sll_real64                        :: eta1_min !< left corner
    !     sll_real64                        :: eta1_max !< right corner
    !     sll_real64, dimension(:), pointer :: wsave !< array used by fftpack
    !     sll_real64, dimension(:), pointer :: work  !< array used by fftpack
    !  end type poisson_1d_periodic

    !> Structure to solve Poisson equation on 1d domain. Mesh is cartesian and
    !> could be irregular. Numerical method is using finite elements.
    type :: sll_t_poisson_1d_fem
        class(sll_t_cartesian_mesh_1d), private, pointer  :: cartesian_mesh

        !> Spline basis
        sll_int32 , private :: spline_degree !<Degree of the Bsplines
        type(sll_t_arbitrary_degree_spline_1d), pointer, private :: bspline  !<Bspline object

        sll_real64, dimension(:), allocatable :: fem_solution
        sll_real64, dimension(:), allocatable, private :: stiffn_matrix_first_line
        !For the L2-Norm
        sll_real64, dimension(:), allocatable, private :: mass_matrix_first_line

        !For Precomputation do one FFT of stiffn_matrix_first_line
        sll_comp64, dimension(:), allocatable, private :: stiffn_matrix_first_line_fourier
        sll_comp64, dimension(:), allocatable, private :: mass_matrix_first_line_fourier
        !sll_comp64 , dimension(:),allocatable, private :: circulant_first_line_fourier
        sll_real64, private :: seminorm

        !FFT Plans
        type(sll_t_fft), private :: forward_fftplan 
        type(sll_t_fft), private :: backward_fftplan
        sll_int32,private ::num_cells

        sll_int32,private :: boundarycondition

        sll_real64, private :: tr_stiffinv_massm
        sll_int32 :: num_sample     !number of samples for the rhs, for statistic reasons
    contains

        procedure,  pass(self) :: initialize =>sll_initialize_poisson_1d_fem
        !procedure,  pass(poisson) :: new =>initialize_poisson_1d_fem
        !procedure,  pass(poisson) :: initialize =>initialize_poisson_1d_fem
        procedure,  pass(self) :: delete =>sll_delete_poisson_1d_fem
        !procedure,  pass(self) :: solve =>solve_poisson_1d_fem
        procedure,  pass(self) :: eval_solution=>poisson_1d_fem_eval_solution
        procedure,  pass(self) :: eval_solution_derivative=>poisson_1d_fem_eval_solution_derivative
        procedure, pass(self) :: H1seminorm_solution=>poisson_1d_fem_H1seminorm_solution
        procedure, pass(self) :: calc_H1seminorm_solution=>poisson_1d_fem_calc_H1seminorm_solution

        procedure, pass(self) :: L2norm_solution=>poisson_1d_fem_L2norm_solution
        procedure, pass(self) :: set_solution=>poisson_1d_fem_set_solution

        procedure, pass(self) :: H1seminorm_solution_stoch_approx=>&
                                            poisson_1d_fem_H1seminorm_solution_stoch_approx
        !        module interface solve
        !        procedure,  pass(self) :: solve_poisson_1d_fem
        !        end interface

        generic ,  public:: solve => solve_poisson_1d_fem_rhs,&
            solve_poisson_1d_fem_functionrhs,&
            solve_poisson_1d_fem_rhs_and_get
        procedure,  pass(self), private :: solve_poisson_1d_fem_rhs
        procedure,pass(self),  public ::  solve_poisson_1d_fem_functionrhs
        procedure,pass(self),  private ::  solve_poisson_1d_fem_rhs_and_get

        procedure, pass(self), private :: set_tr_stiffinv_massm=>poisson_1d_fem_tr_stiffinv_massm

        procedure, pass(self), public :: get_rhs_from_klimontovich_density_weighted=>&
            poisson_1d_fem_get_rhs_from_klimontovich_density_weighted
        procedure, pass(self), public :: get_rhs_from_klimontovich_density=>&
            poisson_1d_fem_get_rhs_from_klimontovich_density
        procedure, pass(self), public :: get_rhs_from_klimontovich_density_moment=>&
            poisson_1d_fem_get_rhs_from_klim_dens_weight_moment
        procedure,pass(self), public :: get_rhs_from_function=>&
            sll_poisson_1d_fem_get_rhs_from_function
    end type sll_t_poisson_1d_fem

    !Interface for one dimensional function for right hand side
    abstract interface
        function sll_i_poisson_1d_fem_rhs_function(x) result(y)
            use sll_m_working_precision
            sll_real64, dimension(:),intent(in) :: x
            sll_real64, dimension(size(x)) :: y
        endfunction
    endinterface

    !     interface solve_poisson_1d_fem
    !                  module procedure solve_poisson_1d_fem_rhs
    !                  module procedure solve_poisson_1d_fem_functionrhs
    !                  module procedure solve_poisson_1d_fem_rhs_and_get
    !    end interface solve_poisson_1d_fem
    interface delete
        module procedure sll_delete_poisson_1d_fem
    endinterface
!    interface new
!        module procedure sll_f_new_poisson_1d_fem
!    endinterface
contains

    !>Destructor
    subroutine sll_delete_poisson_1d_fem(self,ierr)
        class(sll_t_poisson_1d_fem),intent(inout) :: self     !< Solver data structure
        sll_int32, intent(out)                :: ierr    !< error code
        SLL_DEALLOCATE_ARRAY(self%stiffn_matrix_first_line,ierr)
        SLL_DEALLOCATE_ARRAY(self%mass_matrix_first_line,ierr)
        SLL_DEALLOCATE_ARRAY(self%fem_solution,ierr)
        SLL_DEALLOCATE_ARRAY(self%mass_matrix_first_line_fourier,ierr)
        SLL_DEALLOCATE_ARRAY(self%stiffn_matrix_first_line_fourier,ierr)
        call sll_s_fft_free(self%backward_fftplan)
        call sll_s_fft_free(self%forward_fftplan)
    endsubroutine


    function  sll_f_new_poisson_1d_fem(cartesian_mesh_1d,spline_degree,boundarycondition, ierr) &
            result(solver)
        type(sll_t_poisson_1d_fem), pointer :: solver     !< Solver data structure
        class(sll_t_cartesian_mesh_1d), intent(in),pointer  :: cartesian_mesh_1d !< Logical mesh
        sll_int32, intent(out)                :: ierr    !< error code
        sll_int32, intent(in)                :: spline_degree !<Degree of the bsplines
        sll_int32, intent(in) :: boundarycondition !<Boundary Condition Descriptor

        SLL_ALLOCATE(solver,ierr)
        call sll_initialize_poisson_1d_fem(solver,  cartesian_mesh_1d,spline_degree,boundarycondition,ierr)
    endfunction


    subroutine sll_initialize_poisson_1d_fem(self, cartesian_mesh_1d,spline_degree,boundarycondition, ierr)
        class(sll_t_poisson_1d_fem),intent(inout) :: self     !< Solver data structure
        class(sll_t_cartesian_mesh_1d), intent(in),pointer  :: cartesian_mesh_1d !< Logical mesh
        sll_int32, intent(out)                :: ierr    !< error code
        sll_int32, intent(in)                :: spline_degree !<Degree of the bsplines
        sll_int32, intent(in) :: boundarycondition !<Boundary Condition Descriptor
        ierr=0
        self%boundarycondition=boundarycondition !sll_p_periodic
        !self%boundarycondition=sll_p_dirichlet !sll_p_periodic

        !scale_matrix_equation=1.0_f64
        self%cartesian_mesh=>cartesian_mesh_1d
        self%num_cells=self%cartesian_mesh%num_cells
        self%spline_degree=spline_degree
        !Allocate spline
        selectcase(self%boundarycondition)
            case(sll_p_periodic)
                self%bspline=>sll_f_new_arbitrary_degree_spline_1d(spline_degree,self%cartesian_mesh%eta1_nodes(), self%cartesian_mesh%num_nodes(), &
                    sll_p_periodic_arbitrary_deg_spline)
            case(sll_p_dirichlet)
                self%bspline=>sll_f_new_arbitrary_degree_spline_1d( spline_degree, &
                    self%cartesian_mesh%eta1_nodes(), self%cartesian_mesh%num_nodes(), &
                    sll_p_periodic_arbitrary_deg_spline)
            case default
                print *, "This boundary condition is not supported!"
                stop
        endselect

        SLL_ALLOCATE(self%fem_solution(self%num_cells),ierr)

        !FFT--------------------------------------------------------------------------------------
        !To get the same output as in the MATLAB example use
        SLL_ALLOCATE(self%stiffn_matrix_first_line_fourier(1:self%num_cells/2+1), ierr)
        self%stiffn_matrix_first_line_fourier = (0._f64,0._f64)
        call sll_s_fft_init_r2c_1d(self%forward_fftplan, self%num_cells, &
             self%fem_solution,self%stiffn_matrix_first_line_fourier)
        call sll_s_fft_init_c2r_1d(self%backward_fftplan, &
             self%num_cells, &
             self%stiffn_matrix_first_line_fourier,self%fem_solution)
        SLL_DEALLOCATE_ARRAY(self%stiffn_matrix_first_line_fourier, ierr)

        !------------------------------------------------------------------------------------------
        !Assemble stiffnes matrix


        call sll_poisson_1d_fem_assemble_stiffn_matrix(self,ierr)
        !prepare fourier transformation for stiffness matrix

        SLL_ALLOCATE(self%stiffn_matrix_first_line_fourier(1:self%num_cells/2+1), ierr)
        self%stiffn_matrix_first_line_fourier=(0._f64,0._f64)
        call sll_poisson_1d_fft_precomputation(self,self%stiffn_matrix_first_line, &
            self%stiffn_matrix_first_line_fourier, ierr)
        !Assemble mass matrix for L2-Norm
        call sll_poisson_1d_fem_assemble_mass_matrix(self, ierr)

        SLL_ALLOCATE(self%mass_matrix_first_line_fourier(1:self%num_cells/2+1), ierr)
        self%mass_matrix_first_line_fourier=(0._f64,0._f64)
        call sll_poisson_1d_fft_precomputation(self,self%mass_matrix_first_line, &
            self%mass_matrix_first_line_fourier, ierr)

        call self%set_tr_stiffinv_massm()
        
    endsubroutine

    subroutine poisson_1d_fem_set_solution(self, solution_vector)
        class(sll_t_poisson_1d_fem),intent(inout) :: self     !< Solver data structure
        sll_real64, dimension(:) :: solution_vector
        SLL_ASSERT(size(solution_vector)==self%num_cells)

        self%fem_solution=solution_vector
        !Since the solution has been set we have to reset the seminorm
        self%seminorm=-1._f64
    endsubroutine

    !<Calculates the inhomogenity b={<f, \phi_i>, i =1, ... N} with given function f
    !<by Gauss Legendre integration
    function sll_poisson_1d_fem_get_rhs_from_function(self, eval_function,  n_quadrature_points_user) &
            result( rhs )
        implicit none
        class(sll_t_poisson_1d_fem),intent(in) :: self     !< Solver data structure
        procedure (sll_i_poisson_1d_fem_rhs_function) :: eval_function
        sll_real64, dimension(self%num_cells ) :: rhs !<Right hand side
        sll_real64, dimension(:,:), allocatable   :: quadrature_points_weights
        sll_real64, dimension(:,:), allocatable   :: bspline_qpoint_values
        integer, intent(in), optional      :: n_quadrature_points_user  !<Number of quadrature points for GL Integration
        sll_int32                             :: n_quadrature_points
        sll_int32 :: cell
        sll_int32 :: i,j
        sll_real64,dimension(2) :: cellmargin
        sll_int32 :: ierr=0
        !Set right hand side to zero
        rhs=0._f64
        !Get Gauss Legendre points and weights to be exact for the selected spline degree
        !Note a Bspline is a piecewise polynom
        if ( .not. present(n_quadrature_points_user)) then
            !Gauss Legendre, 2m-1
            n_quadrature_points=ceiling(0.5_f64*real(2*self%spline_degree+1,8))
            !            !Gauss Lobatto, 2m-3
            !            n_quadrature_points=ceiling(0.5_f64*real(2*spline_degree+3))
        else
            n_quadrature_points=n_quadrature_points_user
        endif


        SLL_ALLOCATE(quadrature_points_weights(2,1:n_quadrature_points),ierr)
        SLL_CLEAR_ALLOCATE(bspline_qpoint_values(1:self%spline_degree+1,1:n_quadrature_points), ierr)

        !Quadrature for each support cell
        do cell=1, self%num_cells
            cellmargin=sll_o_cell_margin(self%cartesian_mesh, cell)
            quadrature_points_weights=sll_f_gauss_legendre_points_and_weights(n_quadrature_points,&
                cellmargin(1),cellmargin(2))
            !quadrature_points_weights=gauss_lobatto_points_and_weights(n_quadrature_points, knots_mesh(cell), knots_mesh(cell+1))
            do i=1,n_quadrature_points
                bspline_qpoint_values(:,i)=sll_f_splines_at_x(self%bspline, cell, quadrature_points_weights(1,i ))
                !reorientate to the right side
                !bspline_qpoint_values(:,i) = bspline_qpoint_values((spline_degree+1):1:-1,i)
                !call sll_o_display(bspline_qpoint_values(:,i), "(F8.4)")
            enddo

            !Loop over all splines with support in self cell
            !            do j=0,spline_degree
            !                i=mod(cell +j, n_cells)+1
            !                fem_inhomogenity_steady(i)= fem_inhomogenity_steady(i) + &
                !                    dot_product(bspline_qpoint_values(1+j, :)*eval_function(quadrature_points_weights(1,:) ), quadrature_points_weights(2,:))
            !            enddo
            do j=0,self%spline_degree
                !i=cell - j
                i=cell - self%spline_degree + j
                if (i > self%num_cells) then
                    i= i - self%num_cells
                elseif (i < 1) then
                    i = i + self%num_cells
                endif
                rhs(i)= rhs(i) + &
                    dot_product(bspline_qpoint_values(1+j, :)*eval_function(quadrature_points_weights(1,:) ), quadrature_points_weights(2,:))
            enddo
            !enddo

        enddo
        !call sll_o_display(knots(1:degree+2), "(F10.5)" )
        !call sll_o_display(quadrature_points_weights(1,:), "(F10.5)" )
        !call sll_o_display(quadrature_points_weights(2,:), "(F10.5)" )

        !call sll_o_display(bspline_qpoint_values, "(F10.5)" )
        !write(*,*)
        !call sll_o_display(fem_matrix_period, "(F10.5)" )

        !!<SLL_DEALLOCATE_ARRAY( bspline_qpoint_values,ierr)
        SLL_DEALLOCATE_ARRAY( quadrature_points_weights,ierr)

    endfunction




    !<Does precomputation for the FFT solver
    subroutine sll_poisson_1d_fft_precomputation(self,circulant_matrix_first_line, &
            circulant_matrix_first_line_fourier,ierr)
        class(sll_t_poisson_1d_fem),intent(inout) :: self     !< Solver data structure
        sll_int32, intent(out) :: ierr
        sll_real64, dimension(:),allocatable, intent(in)   :: circulant_matrix_first_line
        sll_comp64, dimension(self%num_cells/2 +1),intent(out)  ::circulant_matrix_first_line_fourier
        sll_real64, dimension(:) :: circulantvector(size(circulant_matrix_first_line))
        integer :: N
        ierr=0
        !Determine dimension of problem
        N=size(circulant_matrix_first_line)
        SLL_ASSERT(sll_f_is_power_of_two(int(N,i64)))
        !Generate Circulant Seed

        !Circulant seed c to generate stiffness matrix as a circular matrix
        !Remember c is not the first line of the circulant matrix
        !self is the first line of the matrix
        ! (c_1  0 0  0  .... c_4   c_3  c_2)
        !Note self only works because the period is symmetric
        circulantvector=0._f64
        circulantvector(1)=circulant_matrix_first_line(1)
        circulantvector(2:N)=circulant_matrix_first_line(N:2:-1)
        circulantvector=circulantvector

        !        if (allocated(circulant_matrix_first_line_fourier)) then
        !            SLL_DEALLOCATE_ARRAY(circulant_matrix_first_line_fourier, ierr)
        !        endif
        !        SLL_CLEAR_ALLOCATE(circulant_matrix_first_line_fourier(1:N/2+1), ierr)

        !Fourier transform the circulant seed
        call sll_s_fft_exec_r2c_1d(self%forward_fftplan,circulantvector,circulant_matrix_first_line_fourier)

        !circulant_matrix_first_line_fourier(1)=1.0_f64
        !        sll_real64:: matrix_condition
        !        matrix_condition=maxval(abs(circulant_matrix_first_line_fourier))/minval(abs(circulant_matrix_first_line_fourier))
        !        print *, "Corrected Stiffnes Matrix Condition: ", matrix_condition
    end subroutine


    subroutine sll_poisson_1d_fem_assemble_mass_matrix(self,ierr)
        class(sll_t_poisson_1d_fem),intent(inout) :: self     !< Solver data structure
        sll_int32, intent(out)                :: ierr    !< error code

        sll_real64, dimension(self%spline_degree+1) :: mass_matrix_period
        mass_matrix_period=sll_gen_fem_bspline_matrix_period ( self%bspline , &
            self%cartesian_mesh%eta1_nodes(), sll_f_splines_at_x , ierr)
        SLL_CLEAR_ALLOCATE(self%mass_matrix_first_line(1:self%num_cells),ierr)
        !First line of stiffness matrix
        self%mass_matrix_first_line=0.0_f64
        self%mass_matrix_first_line(self%num_cells- (size(mass_matrix_period) -1) +1: self%num_cells)=&
            mass_matrix_period(size(mass_matrix_period):2:-1)
        self%mass_matrix_first_line(1:size(mass_matrix_period))=mass_matrix_period
    endsubroutine

    !<Generates the FEM matrix for arbitrary spline evaluation functions
    function sll_gen_fem_bspline_matrix_period (bspline_arbitrary_degree, knots, spline_evaluation,ierr, quadrature_points  ) &
            result( fem_matrix_period)
        type(sll_t_arbitrary_degree_spline_1d), pointer, intent(in) :: bspline_arbitrary_degree
        sll_int32, intent(out)                :: ierr    !< error code
        sll_real64, dimension(:), intent(in)     :: knots
        sll_real64, dimension(:), allocatable    :: fem_matrix_period
        procedure (sll_f_splines_at_x) :: spline_evaluation
        integer                   :: degree
        sll_real64, dimension(:,:), allocatable   :: quadrature_points_weights
        sll_real64, dimension(:,:), allocatable   :: bspline_qpoint_values
        integer                             :: quadrature_npoints
        integer, intent(in), optional      :: quadrature_points
        integer :: cell
        integer :: i,j
        degree=bspline_arbitrary_degree%degree
        !Get Gauss Legendre points and weights to be exact for the selected spline degree
        !Note a Bspline is a piecewise polynom
        if ( .not. present(quadrature_points)) then
            quadrature_npoints=ceiling(0.5_f64*real(2*degree+1,8))
            !!!JUST FOR DEBUGGING
            !quadrature_npoints=10
        else
            quadrature_npoints=quadrature_points
        endif

        SLL_ALLOCATE(quadrature_points_weights(2,1:quadrature_npoints),ierr)
        !Get spline values at
        SLL_ALLOCATE(bspline_qpoint_values(1:degree+1,1:quadrature_npoints), ierr)
        SLL_ASSERT(size(knots) >degree+1)

        SLL_CLEAR_ALLOCATE(fem_matrix_period(1:degree+1),ierr)

        !fem_matrix_period=0_f64
        !Quadrature for each support cell
        do cell=1,degree+1
            quadrature_points_weights=sll_f_gauss_legendre_points_and_weights(quadrature_npoints, knots(cell), knots(cell+1))
            do i=1,quadrature_npoints
                bspline_qpoint_values(:,i)=spline_evaluation(bspline_arbitrary_degree, cell, quadrature_points_weights(1,i ))
                !reorientate to the right side
                !bspline_qpoint_values(:,i) = bspline_qpoint_values((degree+1):1:-1,i)
                !shift in order to account for intersecting and limited support
                bspline_qpoint_values(:,i) = eoshift(bspline_qpoint_values(:,i), (cell-1))
            enddo
            do j=1,degree+1
                !forall (j=1:degree+1)
                fem_matrix_period(j)= fem_matrix_period(j) + &
                    dot_product(bspline_qpoint_values(1, :)*bspline_qpoint_values(j,:), quadrature_points_weights(2,:)  )
                !fem_matrix_period(j)= fem_matrix_period(j) + &
                    !   dot_product(bspline_qpoint_values(j, :), quadrature_points_weights(2,:)  )
            enddo
            !endforall
        enddo
        !!This depends on the bspline library
        !!Essentially we scale so that each spline has volume 1
        !!fem_matrix_period=fem_matrix_period*n_cells*n_cells*n_cells

        SLL_DEALLOCATE_ARRAY( quadrature_points_weights,ierr)
    endfunction


    subroutine sll_poisson_1d_fem_assemble_stiffn_matrix(self,ierr)
        class(sll_t_poisson_1d_fem),intent(inout) :: self     !< Solver data structure
        sll_int32, intent(out)                :: ierr    !< error code
        !sll_real64 :: fem_matrix_period_sum

        sll_real64, dimension(self%spline_degree+1) :: stiffn_matrix_period
        stiffn_matrix_period=sll_gen_fem_bspline_matrix_period ( self%bspline , &
            self%cartesian_mesh%eta1_nodes(), sll_f_spline_derivatives_at_x , ierr)
        SLL_CLEAR_ALLOCATE(self%stiffn_matrix_first_line(1:self%num_cells),ierr)
        !First line of stiffness matrix
        self%stiffn_matrix_first_line=0.0_f64

        !        self%stiffn_matrix_first_line(self%num_cells-1 -1- (size(stiffn_matrix_period) -1) +1: self%num_cells-1-1)=&
            !            stiffn_matrix_period(size(stiffn_matrix_period):2:-1)
        !
        self%stiffn_matrix_first_line(self%num_cells:self%num_cells-(self%spline_degree -1):-1)=&
            stiffn_matrix_period(2:self%spline_degree +1)
        self%stiffn_matrix_first_line(1:size(stiffn_matrix_period))=stiffn_matrix_period

        ! !First line of stiffness matrix
        !        stiffn_matrix_first_line=0.0_f64
        !        stiffn_matrix_first_line(n_knots-1- (size(stiffn_matrix_period) -1) +1: n_knots-1)=&
            !            stiffn_matrix_period(size(stiffn_matrix_period):2:-1)
        !        stiffn_matrix_first_line(1:size(stiffn_matrix_period))=stiffn_matrix_period
        !
        !        fem_matrix_period_sum=2*sum(stiffn_matrix_period(2:size(stiffn_matrix_period)))+ stiffn_matrix_period(1)

        !        fem_matrix_period_sum=2*sum(stiffn_matrix_period(2:size(stiffn_matrix_period)))+ stiffn_matrix_period(1)
        !
        !                if  ((fem_matrix_period_sum > FEM_MASS_MATRIX_TOL) &
            !                            .OR. (fem_matrix_period_sum - sum(stiffn_matrix_first_line) > FEM_MASS_MATRIX_TOL)  ) then
        !                    print *, "FEM matrix assembly failed: "
        !                    print *, fem_matrix_period_sum    , sum(stiffn_matrix_first_line)
        !!                    call sll_o_display(stiffn_matrix_period, "(F12.8)" )
        !                    stop
        !                endif
    endsubroutine

    subroutine solve_poisson_1d_fem_rhs_and_get(self, field, rhs)
        class(sll_t_poisson_1d_fem),intent(inout) :: self
        sll_real64, dimension(:), intent(in)      :: rhs
        sll_real64, dimension(self%num_cells), intent(out)     :: field

        call  solve_poisson_1d_fem_rhs(self, rhs)
        field=self%fem_solution

    end subroutine

    subroutine solve_poisson_1d_fem_rhs(self, rhs)
        class(sll_t_poisson_1d_fem),intent(inout) :: self
        sll_real64, dimension(:), intent(in)      :: rhs
      !  sll_real64, dimension(size(rhs)) :: rhs_bc

       ! rhs_bc=rhs

        if (  self%boundarycondition==sll_p_dirichlet ) then
            !rhs_bc(1:self%bspline%degree*2)=0
            !rhs_bc(self%num_cells-(self%bspline%degree*2)-2:self%num_cells)=0
            !self%fem_solution(1:self%bspline%degree+3)=0
            !self%fem_solution(self%num_cells-(self%bspline%degree)-2:self%num_cells )=0
        call solve_poisson_1d_fem_rhs_dirichlet(self, rhs)

        else
        SLL_ASSERT(size(rhs)==self%num_cells)
            self%fem_solution=poisson_1d_fem_solve_circulant_matrix_equation(self, &
            self%stiffn_matrix_first_line_fourier ,rhs )
        endif

        !Set the seminorm

        !STOCHASTIC APPROX
        call self%H1seminorm_solution_stoch_approx(rhs, self%num_sample)
        !self%seminorm=dot_product(self%fem_solution, rhs)

!        if (  self%boundarycondition==sll_p_dirichlet ) then
!            !         rhs_bc(1:self%bspline%degree+3)=0
!            !         rhs_bc(self%num_cells-(self%bspline%degree)-3:self%num_cells)=0
!            !self%fem_solution(1:self%bspline%degree+3)=0
!            !self%fem_solution(self%num_cells-(self%bspline%degree)-3:self%num_cells )=0
!        endif

    endsubroutine

    !Not efficient yet, in the future store LU decomposition
    subroutine solve_poisson_1d_fem_rhs_dirichlet(self, rhs)
        class(sll_t_poisson_1d_fem),intent(inout) :: self
        sll_real64, dimension(:), intent(in)      :: rhs
        !sll_real64, dimension(size(rhs)) :: rhs_bc
        sll_real64, dimension(self%spline_degree+1) :: femperiod
        sll_int32 :: N ,KD
        sll_real64, dimension(self%spline_degree*3 +1, self%num_cells):: AB
        sll_int32 :: idx,ierr
        sll_int32, dimension( self%num_cells ) :: IPIV

       femperiod=sll_gen_fem_bspline_matrix_period ( self%bspline , &
            self%cartesian_mesh%eta1_nodes(), sll_f_spline_derivatives_at_x , ierr)
        N=self%num_cells
        KD=self%spline_degree
        self%fem_solution =rhs

        AB=0._f64
        !Set up the diaognal for dirichlet
        AB(2*KD+1,1)=1.0_f64
         self%fem_solution(1)=0._f64
        AB(2*KD+1,N)=1.0_f64
         self%fem_solution(N)=0._f64
        AB(2*KD+1,2:N-1)=femperiod(1)


!        !Fill up the lower diaognals
!        do idx=2,KD+1
!            AB(3*KD+1-(idx-1), idx:N)=femperiod(idx)
!        enddo
!


        !Lower diagonals
        do idx=1,KD
            AB(3*KD+1-(KD-idx), 1:N-idx -1)=femperiod(idx+1)
        enddo
        !Upper diagonals
        do idx=1,KD
            AB(2*KD-(idx-1), 2+idx:N)=femperiod(idx+1)
        enddo


        call DGBSV( N, KD, KD, 1, AB, self%spline_degree*3 +1, IPIV,  self%fem_solution, N, ierr )

        !call DPBSV( 'U' , N , KD, 1, AB, KD+1,  self%fem_solution, N, ierr )

           if (ierr/=0) then
            print *, ierr

            endif

    endsubroutine

    !> @brief Solves the poisson equation for a given right hand side function
    !> @param self pointer to a sll_t_poisson_1d_fem object.
    !> @param rhs right hand side function of type sll_i_poisson_1d_fem_rhs_function
    subroutine solve_poisson_1d_fem_functionrhs(self, rhs_fun)
        implicit none
        class(sll_t_poisson_1d_fem),intent(inout) :: self
        procedure (sll_i_poisson_1d_fem_rhs_function) :: rhs_fun
        sll_real64, dimension(self%num_cells) :: rhs

        rhs=sll_poisson_1d_fem_get_rhs_from_function(self, rhs_fun,10)
        call solve_poisson_1d_fem_rhs(self,rhs  )
    endsubroutine




    function poisson_1d_fem_solve_circulant_matrix_equation(self, matrix_fl_fourier,rightside ) &
            result(solution)
        class(sll_t_poisson_1d_fem),intent(inout) :: self
        sll_real64, dimension(:), intent(in)   :: rightside
        sll_comp64, dimension(:), intent(in)   :: matrix_fl_fourier
        sll_real64, dimension(:) :: constant_factor(size(rightside))
        sll_real64 , dimension(:) :: solution(size(rightside))
        integer :: N,ierr

        !Since the input data is real, the data in fourier space satisfies
        !X_{N-k}=X_k^* (complex conjugate) which can be stored more efficiently,
        !so for the complex data there will be only allocated space for N/2 +1 values
        sll_comp64 , dimension(:),allocatable ::constant_factor_fourier
        sll_comp64 , dimension(:),allocatable :: data_complex

        N=size(rightside)
        SLL_ASSERT(N==self%num_cells)
        SLL_ASSERT(sll_f_is_power_of_two(int(N,i64)))

        SLL_ASSERT(N/2+1==size(matrix_fl_fourier))

        !Determine dimension of problem

        SLL_ALLOCATE(constant_factor_fourier(1:N/2+1), ierr)
        constant_factor_fourier = (0._f64,0._f64)
        SLL_ALLOCATE(data_complex(1:N/2+1),ierr)
        data_complex = (0._f64,0._f64)
        constant_factor=0._f64
        constant_factor=rightside

        call sll_s_fft_exec_r2c_1d(self%forward_fftplan,constant_factor,constant_factor_fourier)
        !constant_factor_fourier(1)=0
        data_complex=constant_factor_fourier/(matrix_fl_fourier)
        data_complex(1)=(0._f64,0._f64)
        SLL_DEALLOCATE_ARRAY(constant_factor_fourier,ierr)

        call sll_s_fft_exec_c2r_1d(self%backward_fftplan,data_complex,solution)
        SLL_DEALLOCATE_ARRAY(data_complex,ierr)
        !Somehow the normalization does not do what it should do:
        solution=solution/N
    endfunction


    !< knots are the interpolant points and have nothing to do with the mesh
    function poisson_1d_fem_bspline_basis_to_realvals ( self , bspline_vector, knots_eval_in, interpolfun_user ) &
            result(realvals)
        class(sll_t_poisson_1d_fem),intent(in) :: self
        procedure (sll_f_splines_at_x),pointer :: interpolfun
        procedure (sll_f_splines_at_x), optional :: interpolfun_user
        sll_real64, dimension(:), intent(in)     :: bspline_vector
        sll_real64, dimension(:), intent(in)     :: knots_eval_in
        sll_real64, dimension(size(knots_eval_in))    :: knots_eval
        sll_real64 ::  realvals(size(knots_eval_in))
        sll_real64 :: b_contribution(self%bspline%degree+1,size(knots_eval_in))
        sll_int32,  dimension(size(knots_eval_in))  :: b_idx
        sll_int64 :: eval_idx
        sll_int :: b_contrib_idx
        sll_int, dimension(size(knots_eval_in)) :: cell

        !Check Boundary Conditions for evaluation
        !knots_eval=ensure_boundaryc(knots_eval_in)
        knots_eval=knots_eval_in

        if (present(interpolfun_user)) then
            interpolfun=>interpolfun_user
        else
            print *, "Warning no Interpolation function given!"
            interpolfun=>sll_f_splines_at_x
        endif
        SLL_ASSERT( size(knots_eval)==size(realvals))
        SLL_ASSERT( self%num_cells==size(bspline_vector))

        realvals=0._f64

        cell=sll_o_cell(self%cartesian_mesh, knots_eval)
        !cell= bspline_fem_solver_1d_cell_number(knots_mesh, knots_eval(eval_idx))
        !Loop over all points to evaluate
        !This should be vectorizzed
        do eval_idx=1,int(size(knots_eval),i64)
            !Get the values for the spline at the eval point
            b_contribution(:,eval_idx)=interpolfun(self%bspline,cell(eval_idx), knots_eval(eval_idx))
        enddo

        do b_contrib_idx=1,self%bspline%degree+1
            !Determine which value belongs to which spline
            !b_idx=cell(eval_idx)  - (b_contrib_idx-1)
            b_idx=cell -(self%bspline%degree+1)  +(b_contrib_idx)

            !Periodicity
            where (b_idx > self%num_cells)
                b_idx = b_idx - self%num_cells
            elsewhere (b_idx < 1)
                b_idx = b_idx + self%num_cells
            endwhere

            realvals=realvals  + bspline_vector(b_idx)*b_contribution(b_contrib_idx,:)
        enddo

        !
        !        do eval_idx=1,size(knots_eval)
        !
        !            !Get the values for the spline at the eval point
        !            b_contribution=interpolfun(self%bspline,cell(eval_idx), knots_eval(eval_idx))
        !
        !            do b_contrib_idx=1,self%bspline%degree+1
        !                !Determine which value belongs to which spline
        !                !b_idx=cell(eval_idx)  - (b_contrib_idx-1)
        !                b_idx=cell(eval_idx) - (self%bspline%degree+1)  +(b_contrib_idx)
        !                !Periodicity
        !                if (b_idx > self%num_cells) then
        !                    b_idx = b_idx - self%num_cells
        !                elseif (b_idx < 1) then
        !                    b_idx = b_idx + self%num_cells
        !                endif
        !                realvals(eval_idx)=realvals(eval_idx)  + bspline_vector(b_idx)*b_contribution(b_contrib_idx)
        !            enddo

        !enddo

    endfunction

    !< Evaluates the first derivative of the solution at the given points knots_eval
    !< The Result is written into eval_solution
    subroutine poisson_1d_fem_eval_solution(self, knots_eval, eval_solution)
        class(sll_t_poisson_1d_fem),intent(in) :: self
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        SLL_ASSERT(size(knots_eval)==size(eval_solution))

        eval_solution=poisson_1d_fem_bspline_basis_to_realvals(self, self%fem_solution,&
            knots_eval, sll_f_splines_at_x)
    endsubroutine

    !< Evaluates the solution at the given points knots_eval
    !< The Result is written into eval_solution
    subroutine poisson_1d_fem_eval_solution_derivative(self, knots_eval, eval_solution)
        class(sll_t_poisson_1d_fem),intent(in) :: self
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        SLL_ASSERT(size(knots_eval)==size(eval_solution))
        eval_solution=poisson_1d_fem_bspline_basis_to_realvals(self, self%fem_solution,&
            knots_eval, sll_f_spline_derivatives_at_x)
    endsubroutine


        !<Returns the squared H1-seminorm of the solution $\Phi$: $|\nabla \Phi|^2$
        function  poisson_1d_fem_H1seminorm_solution(self) result(seminorm)
            class(sll_t_poisson_1d_fem),intent(inout) :: self
            sll_real64 :: seminorm

            seminorm=self%seminorm
        endfunction

    !<Sets the squared H1-seminorm of the solution $\Phi$: $|\nabla \Phi|^2$
    subroutine poisson_1d_fem_calc_H1seminorm_solution(self)
        class(sll_t_poisson_1d_fem),intent(inout) :: self
        sll_real64 :: seminorm
        sll_int32 :: N
        sll_real64 , dimension(:) :: rhs(self%num_cells)
        !Since the input data is real, the data in fourier space satisfies
        !X_{N-k}=X_k^* (complex conjugate) which can be stored more efficiently,
        !so for the complex data there will be only allocated space for N/2 +1 values
        sll_comp64 , dimension(:) :: data_complex(self%num_cells/2+1)

        rhs=self%fem_solution

        !Determine dimension of problem
        N=size(self%fem_solution)
        SLL_ASSERT(N==self%num_cells)
        SLL_ASSERT(sll_f_is_power_of_two(int(N,i64)))

        call sll_s_fft_exec_r2c_1d(self%forward_fftplan,rhs,data_complex)
        data_complex=data_complex*(self%stiffn_matrix_first_line_fourier)
        data_complex(1)=(0._f64,0._f64)
        call sll_s_fft_exec_c2r_1d(self%backward_fftplan,data_complex,rhs)

        !!!
        !Somehow the normalization does not do what it should do:
        rhs=rhs/(N)

        seminorm=dot_product(self%fem_solution, rhs )
        !seminorm=sqrt(dot_product(fem_solution, fem_solution ))
        !        matrix_product=bspline_fem_solver_1d_circulant_matrix_vector_product&
            !            (stiffn_matrix_first_line,fem_solution)
            self%seminorm=seminorm
    endsubroutine

!    !<Sets the squared H1-seminorm of the solution $\Phi$: $|\nabla \Phi|^2$
!    !<This function uses stochastic approximates for the covariances
!    !<matrix of the rhs as
!    subroutine poisson_1d_fem_calc_H1seminorm_solution_stochastic(self)
!        class(sll_t_poisson_1d_fem),intent(inout) :: self
!        sll_real64 :: seminorm
!        sll_int32 :: N
!        sll_int32 :: Nmark !Number of markers
!        sll_real64 , dimension(:) :: rhs(self%num_cells)
!        !Since the input data is real, the data in fourier space satisfies
!        !X_{N-k}=X_k^* (complex conjugate) which can be stored more efficiently,
!        !so for the complex data there will be only allocated space for N/2 +1 values
!        sll_comp64 , dimension(:) :: data_complex(self%num_cells/2+1)
!
!        !Reconstruct the rhs
!        rhs=self%fem_solution
!
!        !Determine dimension of problem
!        N=size(self%fem_solution)
!        SLL_ASSERT(N==self%num_cells)
!        SLL_ASSERT(sll_f_is_power_of_two(int(N,i64)))
!
!        !Ge
!        call fft_apply_plan(self%forward_fftplan,rhs,data_complex)
!        data_complex=data_complex*(self%stiffn_matrix_first_line_fourier)
!        data_complex(1)=0.0_f64
!        call fft_apply_plan(self%backward_fftplan,data_complex,rhs)
!
!        !Somehow the normalization does not do what it should do:
!        rhs=rhs/(N)
!
!        seminorm=dot_product(self%fem_solution, rhs )
!        !seminorm=sqrt(dot_product(fem_solution, fem_solution ))
!        !        matrix_product=bspline_fem_solver_1d_circulant_matrix_vector_product&
!            !            (stiffn_matrix_first_line,fem_solution)
!            self%seminorm=seminorm
!    endsubroutine

    !<Sets the squared H1-seminorm of the solution $\Phi$: $|\nabla \Phi|^2$
    !<This function uses stochastic approximates for the covariances
    !<matrix of the rhs as
    subroutine poisson_1d_fem_H1seminorm_solution_stoch_approx &
                    (self, rhs, num_sample)
              class(sll_t_poisson_1d_fem),intent(inout) :: self
              sll_real64, dimension(:) :: rhs
              sll_int32,  intent(in) :: num_sample
                sll_int32 :: N
            sll_real64 :: bias
        !Determine dimension of problem
        N=size(self%fem_solution)
        SLL_ASSERT(N==size(rhs))
        SLL_ASSERT(sll_f_is_power_of_two(int(N,i64)))

        bias= self%tr_stiffinv_massm&
                 - dot_product(self%fem_solution, rhs )
        !Covariance matrix under the CLT
        !Independent from weights
        bias=bias/num_sample
                ! print *,self%tr_stiffinv_massm

        self%seminorm=dot_product(self%fem_solution, rhs )
                   !bias=0
!  write (*, "(F8.4,A,F8.4)")  bias/ self%seminorm, "   ",  self%tr_stiffinv_massm/ self%seminorm
        self%seminorm=self%seminorm-bias
    endsubroutine

    !Calculates  trace(K^{-1} M)
    subroutine poisson_1d_fem_tr_stiffinv_massm(self)
              class(sll_t_poisson_1d_fem),intent(inout) :: self
        sll_comp64 , dimension(:) :: data_complex(self%num_cells/2+1)
         !sll_int32 :: ierr


!        call sll_poisson_1d_fft_precomputation&
!                    (self,self%mass_matrix_first_line/, &
!                                    data_complex, ierr)
             !print *,self%mass_matrix_first_line_fourier
             !print *, self%mass_matrix_first_line, "#"
!             print *, self%stiffn_matrix_first_line
!
!             stop
        data_complex=self%mass_matrix_first_line_fourier/self%stiffn_matrix_first_line_fourier
        !Remove the offset
        data_complex(1)=(0._f64,0._f64)


        !The factor two comes from the real to complex fourier transform and
        !the symmmetry we introduce
        !We have double Eigenvalues therefore we must multiply with two
        !We also know that sum(image(data_complex))==0
        self%tr_stiffinv_massm=2.0_f64*real(sum(data_complex))

        !Scale
        self%tr_stiffinv_massm=self%tr_stiffinv_massm/sll_o_mesh_area(self%cartesian_mesh)
    endsubroutine


    !<Gives the squared L2-norm of the solution $\Phi$: $|\Phi|^2$
    function  poisson_1d_fem_L2norm_solution(self) result(l2norm)
        class(sll_t_poisson_1d_fem),intent(inout) :: self
        sll_real64 :: l2norm
        sll_int32 :: N
        sll_real64 , dimension(:) :: solution(self%num_cells)
        sll_comp64 , dimension(:) :: data_complex(self%num_cells/2+1)
        solution=self%fem_solution
        !Determine dimension of problem
        N=size(self%fem_solution)
        SLL_ASSERT(N==self%num_cells)
        SLL_ASSERT(sll_f_is_power_of_two(int(N,i64)))

        call sll_s_fft_exec_r2c_1d(self%forward_fftplan,solution,data_complex)
        data_complex=data_complex*(self%mass_matrix_first_line_fourier)
        data_complex(1)=(0._f64,0._f64)
        call sll_s_fft_exec_c2r_1d(self%backward_fftplan,data_complex,solution)
        solution=solution/(N)
        l2norm=dot_product(self%fem_solution, solution )
    endfunction


    function poisson_1d_fem_get_rhs_from_klimontovich_density(self, &
            ppos)  result(rhs)
        class(sll_t_poisson_1d_fem),intent(inout) :: self
        sll_real64, dimension(:), intent(in) ::ppos
        !procedure (sll_f_splines_at_x), optional:: interpolfun_user
        sll_real64, dimension(self%num_cells) :: rhs

        sll_real64, dimension(1) :: pweight
        pweight=1.0_f64

        !        if (present(interpolfun_user)) then
        !            rhs=poisson_1d_fem_get_rhs_from_klimontovich_density_weighted(self,&
            !                ppos,pweight, interpolfun_user)
        !        else
        rhs=poisson_1d_fem_get_rhs_from_klimontovich_density_weighted(self,&
            ppos,pweight)
        !        endif
    endfunction


    function poisson_1d_fem_get_rhs_from_klimontovich_density_weighted( self,&
            ppos, pweight) result(rhs)
        class(sll_t_poisson_1d_fem),intent(inout) :: self
        sll_real64, dimension(:), intent(in) ::ppos
        sll_real64, dimension(:), intent(in) ::pweight
        sll_real64, dimension(self%num_cells) :: rhs

        rhs=poisson_1d_fem_get_rhs_from_klim_dens_weight_moment( self,&
            ppos, pweight,1)
    endfunction



    function poisson_1d_fem_get_rhs_from_klim_dens_weight_moment( self,&
            ppos, pweight, moment) result(rhs)
        class(sll_t_poisson_1d_fem),intent(inout) :: self
        sll_real64, dimension(:), intent(in) ::ppos
        sll_real64, dimension(:), intent(in) ::pweight
        sll_real64, dimension(self%num_cells) :: rhs
        sll_int32 :: moment
        !procedure (sll_f_splines_at_x), optional:: interpolfun_user
        procedure (sll_f_splines_at_x), pointer :: interpolfun
        sll_int32 :: idx
        sll_int32, dimension( size(ppos) ) :: cell
        sll_real64, dimension(self%bspline%degree+1) :: b_contribution
        sll_int :: b_idx
        sll_int :: b_contrib_idx
        rhs=0._f64

        SLL_ASSERT( (size(pweight)==size(ppos) .OR. size(pweight)==1))
        !        if (present(interpolfun_user)) then
        !            interpolfun=>interpolfun_user
        !        else
        !            interpolfun=>sll_f_splines_at_x
        !        endif

        interpolfun=>sll_f_splines_at_x

        cell= sll_o_cell(self%cartesian_mesh, ppos)

        do idx=1,size(ppos)


            b_contribution=interpolfun(self%bspline,cell(idx), ppos(idx))**moment
            !b_contribution=interpolfun(bspline,5, particleposition(idx) - (cell-1)*(knots(2)-knots(1)))

            !if  (self%boundarycondition==sll_p_dirichlet) then
            !                if (cell(idx)>self%num_cells-(self%bspline%degree+1) ) then
            !                b_contribution((self%bspline%degree+1):(self%bspline%degree+1)-(self%num_cells-cell(idx)):-1)=0.0_f64
            !                endif
            !endif


            if (size(pweight)==1) then
                b_contribution=b_contribution*pweight(1)
            else
                b_contribution=b_contribution*pweight(idx)
            endif

            !Amount of splines equals the amount of
            !b_contribution(1:bspline%degree+1)=b_period
            !b_contribution(bspline%degree+2: bspline%degree*2 +1)=b_period(size(b_period)-1:1:-1)
            !First spline has first support in first cell of interval
            do b_contrib_idx=1, self%bspline%degree+1
                !b_idx=cell + (b_contrib_idx- (bspline%degree+1))
                b_idx=cell(idx) - b_contrib_idx+1


                if (b_idx > self%num_cells) then
                    b_idx = b_idx - self%num_cells
                elseif (b_idx < 1) then
                    b_idx = b_idx + self%num_cells
                endif

                SLL_ASSERT(b_idx <= self%num_cells)
                SLL_ASSERT(b_idx >= 1 )
                rhs(b_idx)=rhs(b_idx)+b_contribution(b_contrib_idx)
            enddo

        enddo

        if (  self%boundarycondition==sll_p_dirichlet ) then
            rhs(1:self%bspline%degree+1)=0._f64
            rhs(self%num_cells-(self%bspline%degree):self%num_cells )=0._f64
        endif
    endfunction

    !    function poisson_1d_fem_calculate_residual(self) result(residuum)
    !                class(sll_t_poisson_1d_fem),intent(inout) :: self
    !
    !        sll_real64 :: residuum
    !        residuum= sqrt(sum(( &
        !            poisson_1d_fem_circulant_matrix_vector_product&
        !            ( self%stiffn_matrix_first_line, self%fem_solution)&
        !            - self%fem_inhomogenity &
        !            )**2)) !/sqrt(sum(fem_inhomogenity**2 ))
    !    endfunction
    !
    !
    !
    !
    !    !<can be optimized, with a sparse first line
    !    function  poisson_1d_fem_circulant_matrix_vector_product( circulant_matrix_first_line, vector ) &
        !            result(solution)
    !        sll_real64, dimension(:), intent(in)   :: circulant_matrix_first_line
    !        sll_real64, dimension(:), intent(in)   :: vector
    !        sll_real64, dimension(size(vector)) ::solution
    !
    !        integer idx, N
    !        N=size(circulant_matrix_first_line)
    !        SLL_ASSERT(size(circulant_matrix_first_line)==size(vector))
    !        solution=0
    !        do idx=1, N
    !            solution(idx)=dot_product(cshift(circulant_matrix_first_line, -(idx-1)), vector)
    !            !print *,cshift(circulant_matrix_first_line, -(idx-1)),"#"
    !
    !        enddo
    !    endfunction

    !Interpolate right side on splines
    !In parallel the following steps will be performed:
    !Each core is responsible for a number of cells
    !OR
    !Each core is responsible for his particles
    !interpolfun is optional and sll_f_splines_at_x by default
    !This is the delta function
    !    function interpolate_particles_bsplines( bspline, knots, particleposition_in, interpolfun_user, particleweight) &
        !            result(b)
    !
    !    endfunction

end module
