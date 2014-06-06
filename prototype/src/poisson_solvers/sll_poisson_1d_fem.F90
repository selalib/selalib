
!> Module to solve Poisson equation on one dimensional mesh using Finite Elements
module sll_poisson_1d_fem
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

    use sll_constants
    use sll_logical_meshes
    use sll_arbitrary_degree_splines
    use sll_boundary_condition_descriptors
    use gauss_legendre_integration
    use sll_utilities
    use sll_fft
    implicit none
    !  private
    !  public :: initialize, new, solve

    !  !> Solver data structure
    !  type, public :: poisson_1d_periodic
    !     sll_int32                         :: nc_eta1 !< number of cells
    !     sll_real64                        :: eta1_min !< left corner
    !     sll_real64                        :: eta1_max !< right corner
    !     sll_real64, dimension(:), pointer :: wsave !< array used by fftpack
    !     sll_real64, dimension(:), pointer :: work  !< array used by fftpack
    !  end type poisson_1d_periodic

    !> Structure to solve Poisson equation on 1d domain. Mesh is cartesian and
    !> could be irregular. Numerical method is using finite elements.
    type :: poisson_1d_fem
        type(sll_logical_mesh_1d), private, pointer  :: logical_mesh

        !> Spline basis
        sll_int32 , private :: spline_degree !<Degree of the Bsplines
        type(arbitrary_degree_spline_1d), pointer, private :: bspline  !<Bspline object

        sll_real64, dimension(:), allocatable :: fem_solution
        sll_real64, dimension(:), allocatable, private :: stiffn_matrix_first_line
        !For the L2-Norm
        sll_real64, dimension(:), allocatable, private :: mass_matrix_first_line

        !For Precomputation do one FFT of stiffn_matrix_first_line
        sll_comp64, dimension(:), allocatable, private :: stiffn_matrix_first_line_fourier
        sll_comp64, dimension(:), allocatable, private :: mass_matrix_first_line_fourier
        !sll_comp64 , dimension(:),allocatable, private :: circulant_first_line_fourier


        !FFT Plans
        type(sll_fft_plan), pointer, private :: forward_fftplan => null()
        type(sll_fft_plan), pointer, private :: backward_fftplan => null()
        sll_int32,private ::num_cells

        sll_int32,private :: boundarycondition
    contains
        procedure,  pass(this) :: initialize =>sll_initialize_poisson_1d_fem
        !procedure,  pass(poisson) :: new =>initialize_poisson_1d_fem
        !procedure,  pass(poisson) :: initialize =>initialize_poisson_1d_fem
        procedure,  pass(this) :: delete =>sll_delete_poisson_1d_fem
        !procedure,  pass(this) :: solve =>solve_poisson_1d_fem
        procedure,  pass(this) :: eval_solution=>poisson_1d_fem_eval_solution
        procedure,  pass(this) :: eval_solution_derivative=>poisson_1d_fem_eval_solution_derivative
        procedure, pass(this) :: H1seminorm_solution=>poisson_1d_fem_H1seminorm_solution
        procedure, pass(this) :: L2norm_solution=>poisson_1d_fem_L2norm_solution
        procedure, pass(this) :: set_solution=>poisson_1d_fem_set_solution

        !        module interface solve
        !        procedure,  pass(this) :: solve_poisson_1d_fem
        !        end interface

        generic ,  public:: solve => solve_poisson_1d_fem_rhs,&
            solve_poisson_1d_fem_functionrhs,&
            solve_poisson_1d_fem_rhs_and_get
        procedure,  pass(this), private :: solve_poisson_1d_fem_rhs
        procedure,pass(this),  public ::  solve_poisson_1d_fem_functionrhs
        procedure,pass(this),  private ::  solve_poisson_1d_fem_rhs_and_get


        procedure, pass(this), public :: get_rhs_from_klimontovich_density_weighted=>&
                                    poisson_1d_fem_get_rhs_from_klimontovich_density_weighted
        procedure, pass(this), public :: get_rhs_from_klimontovich_density=>&
                                poisson_1d_fem_get_rhs_from_klimontovich_density
        procedure,pass(this), public :: get_rhs_from_function=>&
                                    sll_poisson_1d_fem_get_rhs_from_function

    end type poisson_1d_fem

    !Interface for one dimensional function for right hand side
    abstract interface
        function poisson_1d_fem_rhs_function(x) result(y)
            use sll_working_precision
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

contains

    !>Destructor
    subroutine sll_delete_poisson_1d_fem(this,ierr)
        class(poisson_1d_fem),intent(inout) :: this     !< Solver data structure
        sll_int32, intent(out)                :: ierr    !< error code
        SLL_DEALLOCATE_ARRAY(this%stiffn_matrix_first_line,ierr)
        SLL_DEALLOCATE_ARRAY(this%mass_matrix_first_line,ierr)
        SLL_DEALLOCATE_ARRAY(this%fem_solution,ierr)
        SLL_DEALLOCATE_ARRAY(this%mass_matrix_first_line_fourier,ierr)
        SLL_DEALLOCATE_ARRAY(this%stiffn_matrix_first_line_fourier,ierr)
        call fft_delete_plan(this%backward_fftplan)
        call fft_delete_plan(this%forward_fftplan)
    endsubroutine


    function  new_poisson_1d_fem(logical_mesh_1d,spline_degree, ierr) &
            result(solver)
        type(poisson_1d_fem), pointer :: solver     !< Solver data structure
        type(sll_logical_mesh_1d), intent(in),pointer  :: logical_mesh_1d !< Logical mesh
        sll_int32, intent(out)                :: ierr    !< error code
        sll_int32, intent(in)                :: spline_degree !<Degree of the bsplines

        SLL_ALLOCATE(solver,ierr)
        call sll_initialize_poisson_1d_fem(solver,  logical_mesh_1d,spline_degree,ierr)
    endfunction


    subroutine sll_initialize_poisson_1d_fem(this, logical_mesh_1d,spline_degree, ierr)
        class(poisson_1d_fem),intent(inout) :: this     !< Solver data structure
        type(sll_logical_mesh_1d), intent(in),pointer  :: logical_mesh_1d !< Logical mesh
        sll_int32, intent(out)                :: ierr    !< error code
        sll_int32, intent(in)                :: spline_degree !<Degree of the bsplines
        ierr=0
        this%boundarycondition=SLL_PERIODIC !SLL_PERIODIC

        !scale_matrix_equation=1.0_f64
        this%logical_mesh=>logical_mesh_1d
        this%num_cells=this%logical_mesh%num_cells
        this%spline_degree=spline_degree
        !Allocate spline
        selectcase(this%boundarycondition)
        case(SLL_PERIODIC)
            this%bspline=>new_arbitrary_degree_spline_1d( spline_degree, &
            sll_mesh_nodes(this%logical_mesh), sll_mesh_num_nodes(this%logical_mesh), &
            PERIODIC_ARBITRARY_DEG_SPLINE)
        case(SLL_DIRICHLET)
            this%bspline=>new_arbitrary_degree_spline_1d( spline_degree, &
            sll_mesh_nodes(this%logical_mesh), sll_mesh_num_nodes(this%logical_mesh), &
            PERIODIC_ARBITRARY_DEG_SPLINE)
        endselect

        SLL_ALLOCATE(this%fem_solution(this%num_cells),ierr)

        !FFT--------------------------------------------------------------------------------------
        !To get the same output as in the MATLAB example use
        SLL_CLEAR_ALLOCATE(this%stiffn_matrix_first_line_fourier(1:this%num_cells/2+1), ierr)
        this%forward_fftplan => fft_new_plan(this%num_cells,this%fem_solution,this%stiffn_matrix_first_line_fourier, &
            FFT_FORWARD+FFT_NORMALIZE)
        this%backward_fftplan=>fft_new_plan(this%num_cells,this%stiffn_matrix_first_line_fourier,this%fem_solution, &
            FFT_INVERSE  + FFT_NORMALIZE_INVERSE)
        SLL_DEALLOCATE_ARRAY(this%stiffn_matrix_first_line_fourier, ierr)

        !------------------------------------------------------------------------------------------
        !Assemble stiffnes matrix

        call sll_poisson_1d_fem_assemble_stiffn_matrix(this,ierr)
        !prepare fourier transformation for stiffness matrix

        SLL_CLEAR_ALLOCATE(this%stiffn_matrix_first_line_fourier(1:this%num_cells/2+1), ierr)
        call sll_poisson_1d_fft_precomputation(this,this%stiffn_matrix_first_line, &
            this%stiffn_matrix_first_line_fourier, ierr)

        !Assemble mass matrix for L2-Norm
        call sll_poisson_1d_fem_assemble_mass_matrix(this, ierr)

        SLL_CLEAR_ALLOCATE(this%mass_matrix_first_line_fourier(1:this%num_cells/2+1), ierr)
        call sll_poisson_1d_fft_precomputation(this,this%mass_matrix_first_line, &
            this%mass_matrix_first_line_fourier, ierr)


    endsubroutine

    subroutine poisson_1d_fem_set_solution(this, solution_vector)
        class(poisson_1d_fem),intent(inout) :: this     !< Solver data structure
        sll_real64, dimension(:) :: solution_vector
        SLL_ASSERT(size(solution_vector)==this%num_cells)

        this%fem_solution=solution_vector
    endsubroutine

    !<Calculates the inhomogenity b={<f, \phi_i>, i =1, ... N} with given function f
    !<by Gauss Legendre integration
    function sll_poisson_1d_fem_get_rhs_from_function(this, eval_function,  n_quadrature_points_user) &
            result( rhs )
        implicit none
        class(poisson_1d_fem),intent(in) :: this     !< Solver data structure
        procedure (poisson_1d_fem_rhs_function) :: eval_function
        sll_real64, dimension(this%num_cells ) :: rhs !<Right hand side
        sll_real64, dimension(:,:), allocatable   :: quadrature_points_weights
        sll_real64, dimension(:,:), allocatable   :: bspline_qpoint_values
        integer, intent(in), optional      :: n_quadrature_points_user  !<Number of quadrature points for GL Integration
        sll_int32                             :: n_quadrature_points
        sll_int32 :: cell
        sll_int32 :: i,j
        sll_real64,dimension(2) :: cellmargin
        sll_int32 :: ierr=0
        !Set right hand side to zero
        rhs=0
        !Get Gauss Legendre points and weights to be exact for the selected spline degree
        !Note a Bspline is a piecewise polynom
        if ( .not. present(n_quadrature_points_user)) then
            !Gauss Legendre, 2m-1
            n_quadrature_points=ceiling(0.5_f64*real(2*this%spline_degree+1))
            !            !Gauss Lobatto, 2m-3
            !            n_quadrature_points=ceiling(0.5_f64*real(2*spline_degree+3))
        else
            n_quadrature_points=n_quadrature_points_user
        endif


        SLL_ALLOCATE(quadrature_points_weights(2,1:n_quadrature_points),ierr)
        SLL_CLEAR_ALLOCATE(bspline_qpoint_values(1:this%spline_degree+1,1:n_quadrature_points), ierr)

        !Quadrature for each support cell
        do cell=1, this%num_cells
            cellmargin=sll_cell_margin(this%logical_mesh, cell)
            quadrature_points_weights=gauss_legendre_points_and_weights(n_quadrature_points,&
                cellmargin(1),cellmargin(2))
            !quadrature_points_weights=gauss_lobatto_points_and_weights(n_quadrature_points, knots_mesh(cell), knots_mesh(cell+1))
            do i=1,n_quadrature_points
                bspline_qpoint_values(:,i)=b_splines_at_x(this%bspline, cell, quadrature_points_weights(1,i ))
                !reorientate to the right side
                !bspline_qpoint_values(:,i) = bspline_qpoint_values((spline_degree+1):1:-1,i)
                !call sll_display(bspline_qpoint_values(:,i), "(F8.4)")
            enddo

            !Loop over all splines with support in this cell
            !            do j=0,spline_degree
            !                i=mod(cell +j, n_cells)+1
            !                fem_inhomogenity_steady(i)= fem_inhomogenity_steady(i) + &
                !                    dot_product(bspline_qpoint_values(1+j, :)*eval_function(quadrature_points_weights(1,:) ), quadrature_points_weights(2,:))
            !            enddo
            do j=0,this%spline_degree
                !i=cell - j
                i=cell - this%spline_degree + j
                if (i > this%num_cells) then
                    i= i - this%num_cells
                elseif (i < 1) then
                    i = i + this%num_cells
                endif
                rhs(i)= rhs(i) + &
                    dot_product(bspline_qpoint_values(1+j, :)*eval_function(quadrature_points_weights(1,:) ), quadrature_points_weights(2,:))
            enddo
            !enddo

        enddo
        !call sll_display(knots(1:degree+2), "(F10.5)" )
        !call sll_display(quadrature_points_weights(1,:), "(F10.5)" )
        !call sll_display(quadrature_points_weights(2,:), "(F10.5)" )

        !call sll_display(bspline_qpoint_values, "(F10.5)" )
        !write(*,*)
        !call sll_display(fem_matrix_period, "(F10.5)" )

        !!<SLL_DEALLOCATE_ARRAY( bspline_qpoint_values,ierr)
        SLL_DEALLOCATE_ARRAY( quadrature_points_weights,ierr)

    endfunction




    !<Does precomputation for the FFT solver
    subroutine sll_poisson_1d_fft_precomputation(this,circulant_matrix_first_line, &
            circulant_matrix_first_line_fourier,ierr)
        class(poisson_1d_fem),intent(inout) :: this     !< Solver data structure
        sll_int32, intent(out) :: ierr
        sll_real64, dimension(:),allocatable, intent(in)   :: circulant_matrix_first_line
        sll_comp64, dimension(this%num_cells/2 +1),intent(out)  ::circulant_matrix_first_line_fourier
        sll_real64, dimension(:) :: circulantvector(size(circulant_matrix_first_line))
        integer :: N
        ierr=0
        !Determine dimension of problem
        N=size(circulant_matrix_first_line)
        SLL_ASSERT(is_power_of_two(int(N,i64)))
        !Generate Circulant Seed

        !Circulant seed c to generate stiffness matrix as a circular matrix
        !Remember c is not the first line of the circulant matrix
        !this is the first line of the matrix
        ! (c_1  0 0  0  .... c_4   c_3  c_2)
        !Note this only works because the period is symmetric
        circulantvector=0
        circulantvector(1)=circulant_matrix_first_line(1)
        circulantvector(2:N)=circulant_matrix_first_line(N:2:-1)
        circulantvector=circulantvector

        !        if (allocated(circulant_matrix_first_line_fourier)) then
        !            SLL_DEALLOCATE_ARRAY(circulant_matrix_first_line_fourier, ierr)
        !        endif
        !        SLL_CLEAR_ALLOCATE(circulant_matrix_first_line_fourier(1:N/2+1), ierr)

        !Fourier transform the circulant seed
        call fft_apply_plan(this%forward_fftplan,circulantvector,circulant_matrix_first_line_fourier)

        !circulant_matrix_first_line_fourier(1)=1.0_f64
        !        sll_real64:: matrix_condition
        !        matrix_condition=maxval(abs(circulant_matrix_first_line_fourier))/minval(abs(circulant_matrix_first_line_fourier))
        !        print *, "Corrected Stiffnes Matrix Condition: ", matrix_condition
    end subroutine


    subroutine sll_poisson_1d_fem_assemble_mass_matrix(this,ierr)
        class(poisson_1d_fem),intent(inout) :: this     !< Solver data structure
        sll_int32, intent(out)                :: ierr    !< error code

        sll_real64, dimension(this%spline_degree+1) :: mass_matrix_period
        mass_matrix_period=sll_gen_fem_bspline_matrix_period ( this%bspline , &
            sll_mesh_nodes(this%logical_mesh), b_splines_at_x , ierr)
        SLL_CLEAR_ALLOCATE(this%mass_matrix_first_line(1:this%num_cells),ierr)
        !First line of stiffness matrix
        this%mass_matrix_first_line=0.0_f64
        this%mass_matrix_first_line(this%num_cells-1 -1- (size(mass_matrix_period) -1) +1: this%num_cells-1-1)=&
            mass_matrix_period(size(mass_matrix_period):2:-1)
        this%mass_matrix_first_line(1:size(mass_matrix_period))=mass_matrix_period
    endsubroutine

    !<Generates the FEM matrix for arbitrary spline evaluation functions
    function sll_gen_fem_bspline_matrix_period (bspline_arbitrary_degree, knots, spline_evaluation,ierr, quadrature_points  ) &
            result( fem_matrix_period)
        type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline_arbitrary_degree
        sll_int32, intent(out)                :: ierr    !< error code
        sll_real64, dimension(:), intent(in)     :: knots
        sll_real64, dimension(:), allocatable    :: fem_matrix_period
        procedure (b_splines_at_x) :: spline_evaluation
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
            quadrature_npoints=ceiling(0.5_f64*real(2*degree+1))
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
            quadrature_points_weights=gauss_legendre_points_and_weights(quadrature_npoints, knots(cell), knots(cell+1))
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


    subroutine sll_poisson_1d_fem_assemble_stiffn_matrix(this,ierr)
        class(poisson_1d_fem),intent(inout) :: this     !< Solver data structure
        sll_int32, intent(out)                :: ierr    !< error code
        !sll_real64 :: fem_matrix_period_sum

        sll_real64, dimension(this%spline_degree+1) :: stiffn_matrix_period
        stiffn_matrix_period=sll_gen_fem_bspline_matrix_period ( this%bspline , &
            sll_mesh_nodes(this%logical_mesh), b_spline_derivatives_at_x , ierr)
        SLL_CLEAR_ALLOCATE(this%stiffn_matrix_first_line(1:this%num_cells),ierr)
        !First line of stiffness matrix
        this%stiffn_matrix_first_line=0.0_f64

        !        this%stiffn_matrix_first_line(this%num_cells-1 -1- (size(stiffn_matrix_period) -1) +1: this%num_cells-1-1)=&
            !            stiffn_matrix_period(size(stiffn_matrix_period):2:-1)
        !
        this%stiffn_matrix_first_line(this%num_cells:this%num_cells-(this%spline_degree -1):-1)=&
            stiffn_matrix_period(2:this%spline_degree +1)
        this%stiffn_matrix_first_line(1:size(stiffn_matrix_period))=stiffn_matrix_period

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
        !!                    call sll_display(stiffn_matrix_period, "(F12.8)" )
        !                    stop
        !                endif
    endsubroutine

    subroutine solve_poisson_1d_fem_rhs_and_get(this, field, rhs)
        class(poisson_1d_fem),intent(inout) :: this
        sll_real64, dimension(:), intent(in)      :: rhs
        sll_real64, dimension(this%num_cells), intent(out)     :: field

        call  solve_poisson_1d_fem_rhs(this, rhs)
        field=this%fem_solution

    end subroutine

    subroutine solve_poisson_1d_fem_rhs(this, rhs)
        class(poisson_1d_fem),intent(inout) :: this
        sll_real64, dimension(:), intent(in)      :: rhs

        SLL_ASSERT(size(rhs)==this%num_cells)
        this%fem_solution=poisson_1d_fem_solve_circulant_matrix_equation(this, &
            this%stiffn_matrix_first_line_fourier ,rhs )

        if (  this%boundarycondition==SLL_DIRICHLET ) then
        this%fem_solution(1:this%bspline%degree+2)=0
        this%fem_solution(this%num_cells-(this%bspline%degree)-1:this%num_cells )=0
        endif
    endsubroutine

    !> @brief Solves the poisson equation for a given right hand side function
    !> @param this pointer to a poisson_1d_fem object.
    !> @param rhs right hand side function of type poisson_1d_fem_rhs_function
    subroutine solve_poisson_1d_fem_functionrhs(this, rhs_fun)
        implicit none
        class(poisson_1d_fem),intent(inout) :: this
        procedure (poisson_1d_fem_rhs_function) :: rhs_fun
        sll_real64, dimension(this%num_cells) :: rhs

        rhs=sll_poisson_1d_fem_get_rhs_from_function(this, rhs_fun,10)
        call solve_poisson_1d_fem_rhs(this,rhs  )
    endsubroutine




    function poisson_1d_fem_solve_circulant_matrix_equation(this, matrix_fl_fourier,rightside ) &
            result(solution)
        class(poisson_1d_fem),intent(inout) :: this
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
        SLL_ASSERT(N==this%num_cells)
        SLL_ASSERT(is_power_of_two(int(N,i64)))

        SLL_ASSERT(N/2+1==size(matrix_fl_fourier))

        !Determine dimension of problem

        SLL_CLEAR_ALLOCATE(constant_factor_fourier(1:N/2+1), ierr)
        SLL_CLEAR_ALLOCATE(data_complex(1:N/2+1),ierr)
        constant_factor=0
        constant_factor=rightside

        call fft_apply_plan(this%forward_fftplan,constant_factor,constant_factor_fourier)
        !constant_factor_fourier(1)=0
        data_complex=constant_factor_fourier/(matrix_fl_fourier)
        data_complex(1)=0
        SLL_DEALLOCATE_ARRAY(constant_factor_fourier,ierr)

        call fft_apply_plan(this%backward_fftplan,data_complex,solution)
        SLL_DEALLOCATE_ARRAY(data_complex,ierr)
        !Somehow the normalization does not do what it should do:
        solution=solution/N
    endfunction


    !< knots are the interpolant points and have nothing to do with the mesh
    function poisson_1d_fem_bspline_basis_to_realvals ( this , bspline_vector, knots_eval_in, interpolfun_user ) &
            result(realvals)
        class(poisson_1d_fem),intent(in) :: this
        procedure (b_splines_at_x),pointer :: interpolfun
        procedure (b_splines_at_x), optional :: interpolfun_user
        sll_real64, dimension(:), intent(in)     :: bspline_vector
        sll_real64, dimension(:), intent(in)     :: knots_eval_in
        sll_real64, dimension(size(knots_eval_in))    :: knots_eval
        sll_real64 ::  realvals(size(knots_eval_in))
        sll_real64 :: b_contribution(this%bspline%degree+1,size(knots_eval_in))
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
            interpolfun=>b_splines_at_x
        endif
        SLL_ASSERT( size(knots_eval)==size(realvals))
        SLL_ASSERT( this%num_cells==size(bspline_vector))

        realvals=0

        cell=sll_cell(this%logical_mesh, knots_eval)
        !cell= bspline_fem_solver_1d_cell_number(knots_mesh, knots_eval(eval_idx))
        !Loop over all points to evaluate
        !This should be vectorizzed
        do eval_idx=1,size(knots_eval)
            !Get the values for the spline at the eval point
            b_contribution(:,eval_idx)=interpolfun(this%bspline,cell(eval_idx), knots_eval(eval_idx))
        enddo

            do b_contrib_idx=1,this%bspline%degree+1
                !Determine which value belongs to which spline
                !b_idx=cell(eval_idx)  - (b_contrib_idx-1)
                b_idx=cell -(this%bspline%degree+1)  +(b_contrib_idx)

                !Periodicity
                where (b_idx > this%num_cells)
                    b_idx = b_idx - this%num_cells
                elsewhere (b_idx < 1)
                    b_idx = b_idx + this%num_cells
                endwhere

                realvals=realvals  + bspline_vector(b_idx)*b_contribution(b_contrib_idx,:)
            enddo

!
!        do eval_idx=1,size(knots_eval)
!
!            !Get the values for the spline at the eval point
!            b_contribution=interpolfun(this%bspline,cell(eval_idx), knots_eval(eval_idx))
!
!            do b_contrib_idx=1,this%bspline%degree+1
!                !Determine which value belongs to which spline
!                !b_idx=cell(eval_idx)  - (b_contrib_idx-1)
!                b_idx=cell(eval_idx) - (this%bspline%degree+1)  +(b_contrib_idx)
!                !Periodicity
!                if (b_idx > this%num_cells) then
!                    b_idx = b_idx - this%num_cells
!                elseif (b_idx < 1) then
!                    b_idx = b_idx + this%num_cells
!                endif
!                realvals(eval_idx)=realvals(eval_idx)  + bspline_vector(b_idx)*b_contribution(b_contrib_idx)
!            enddo

        !enddo

    endfunction

    !< Evaluates the first derivative of the solution at the given points knots_eval
    !< The Result is written into eval_solution
    subroutine poisson_1d_fem_eval_solution(this, knots_eval, eval_solution)
        class(poisson_1d_fem),intent(in) :: this
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        SLL_ASSERT(size(knots_eval)==size(eval_solution))

        eval_solution=poisson_1d_fem_bspline_basis_to_realvals(this, this%fem_solution,&
            knots_eval, b_splines_at_x)
    endsubroutine

    !< Evaluates the solution at the given points knots_eval
    !< The Result is written into eval_solution
    subroutine poisson_1d_fem_eval_solution_derivative(this, knots_eval, eval_solution)
        class(poisson_1d_fem),intent(in) :: this
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        SLL_ASSERT(size(knots_eval)==size(eval_solution))
        eval_solution=poisson_1d_fem_bspline_basis_to_realvals(this, this%fem_solution,&
            knots_eval, b_spline_derivatives_at_x)
    endsubroutine

    !<Gives the squared H1-seminorm of the solution $\Phi$: $|\nabla \Phi|^2$
    function  poisson_1d_fem_H1seminorm_solution(this) result(seminorm)
        class(poisson_1d_fem),intent(inout) :: this
        sll_real64 :: seminorm
        sll_int32 :: N
        sll_real64 , dimension(:) :: solution(this%num_cells)
        !Since the input data is real, the data in fourier space satisfies
        !X_{N-k}=X_k^* (complex conjugate) which can be stored more efficiently,
        !so for the complex data there will be only allocated space for N/2 +1 values
        sll_comp64 , dimension(:) :: data_complex(this%num_cells/2+1)

        !!type(sll_fft_plan), pointer :: forward_fftplan => null()
        !!type(sll_fft_plan), pointer :: backward_fftplan => null()
        solution=this%fem_solution

        !Determine dimension of problem
        N=size(this%fem_solution)
        SLL_ASSERT(N==this%num_cells)
        SLL_ASSERT(is_power_of_two(int(N,i64)))

        call fft_apply_plan(this%forward_fftplan,solution,data_complex)
        data_complex=data_complex*(this%stiffn_matrix_first_line_fourier)
        data_complex(1)=0.0_f64
        call fft_apply_plan(this%backward_fftplan,data_complex,solution)

        !!!
        !Somehow the normalization does not do what it should do:
        solution=solution/(N)

        seminorm=dot_product(this%fem_solution, solution )
        !seminorm=sqrt(dot_product(fem_solution, fem_solution ))
        !        matrix_product=bspline_fem_solver_1d_circulant_matrix_vector_product&
            !            (stiffn_matrix_first_line,fem_solution)
    endfunction

    !<Gives the squared L2-norm of the solution $\Phi$: $|\Phi|^2$
    function  poisson_1d_fem_L2norm_solution(this) result(l2norm)
        class(poisson_1d_fem),intent(inout) :: this
        sll_real64 :: l2norm
        sll_int32 :: N
        sll_real64 , dimension(:) :: solution(this%num_cells)
        sll_comp64 , dimension(:) :: data_complex(this%num_cells/2+1)
        solution=this%fem_solution
        !Determine dimension of problem
        N=size(this%fem_solution)
        SLL_ASSERT(N==this%num_cells)
        SLL_ASSERT(is_power_of_two(int(N,i64)))

        call fft_apply_plan(this%forward_fftplan,solution,data_complex)
        data_complex=data_complex*(this%mass_matrix_first_line_fourier)
        data_complex(1)=0.0_f64
        call fft_apply_plan(this%backward_fftplan,data_complex,solution)
        solution=solution/(N)
        l2norm=dot_product(this%fem_solution, solution )
    endfunction


    function poisson_1d_fem_get_rhs_from_klimontovich_density(this, &
            ppos)  result(rhs)
        class(poisson_1d_fem),intent(inout) :: this
        sll_real64, dimension(:), intent(in) ::ppos
        !procedure (b_splines_at_x), optional:: interpolfun_user
        sll_real64, dimension(this%num_cells) :: rhs

        sll_real64, dimension(1) :: pweight
        pweight=1.0_f64

!        if (present(interpolfun_user)) then
!            rhs=poisson_1d_fem_get_rhs_from_klimontovich_density_weighted(this,&
!                ppos,pweight, interpolfun_user)
!        else
            rhs=poisson_1d_fem_get_rhs_from_klimontovich_density_weighted(this,&
                ppos,pweight)
!        endif
    endfunction


    function poisson_1d_fem_get_rhs_from_klimontovich_density_weighted( this,&
            ppos, pweight) result(rhs)
        class(poisson_1d_fem),intent(inout) :: this
        sll_real64, dimension(:), intent(in) ::ppos
        sll_real64, dimension(:), intent(in) ::pweight
        sll_real64, dimension(this%num_cells) :: rhs
        !procedure (b_splines_at_x), optional:: interpolfun_user
        procedure (b_splines_at_x), pointer :: interpolfun
        sll_int32 :: idx
        sll_int32, dimension( size(ppos) ) :: cell
        sll_real64, dimension(this%bspline%degree+1) :: b_contribution
        sll_int :: b_idx
        sll_int :: b_contrib_idx
        rhs=0

        SLL_ASSERT( (size(pweight)==size(ppos) .OR. size(pweight)==1))
!        if (present(interpolfun_user)) then
!            interpolfun=>interpolfun_user
!        else
!            interpolfun=>b_splines_at_x
!        endif

            interpolfun=>b_splines_at_x

        cell= sll_cell(this%logical_mesh, ppos)

        do idx=1,size(ppos)


            b_contribution=interpolfun(this%bspline,cell(idx), ppos(idx))
            !b_contribution=interpolfun(bspline,5, particleposition(idx) - (cell-1)*(knots(2)-knots(1)))

            !if  (this%boundarycondition==SLL_DIRICHLET) then
!                if (cell(idx)>this%num_cells-(this%bspline%degree+1) ) then
!                b_contribution((this%bspline%degree+1):(this%bspline%degree+1)-(this%num_cells-cell(idx)):-1)=0.0_f64
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
            do b_contrib_idx=1, this%bspline%degree+1
                !b_idx=cell + (b_contrib_idx- (bspline%degree+1))
                b_idx=cell(idx) - b_contrib_idx+1


                if (b_idx > this%num_cells) then
                    b_idx = b_idx - this%num_cells
                elseif (b_idx < 1) then
                    b_idx = b_idx + this%num_cells
                endif

                SLL_ASSERT(b_idx <= this%num_cells)
                SLL_ASSERT(b_idx >= 1 )
                rhs(b_idx)=rhs(b_idx)+b_contribution(b_contrib_idx)
            enddo

        enddo

        if (  this%boundarycondition==SLL_DIRICHLET ) then
        rhs(1:this%bspline%degree+1)=0
        rhs(this%num_cells-(this%bspline%degree):this%num_cells )=0
        endif
    endfunction

!    function poisson_1d_fem_calculate_residual(this) result(residuum)
!                class(poisson_1d_fem),intent(inout) :: this
!
!        sll_real64 :: residuum
!        residuum= sqrt(sum(( &
!            poisson_1d_fem_circulant_matrix_vector_product&
!            ( this%stiffn_matrix_first_line, this%fem_solution)&
!            - this%fem_inhomogenity &
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
    !interpolfun is optional and b_splines_at_x by default
    !This is the delta function
    !    function interpolate_particles_bsplines( bspline, knots, particleposition_in, interpolfun_user, particleweight) &
        !            result(b)
    !
    !    endfunction

end module
