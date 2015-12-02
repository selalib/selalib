
!> Module to solve Poisson equation on one dimensional mesh using Finite Elements
module sll_m_poisson_1d_fd
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

    use sll_m_constants
    use sll_m_cartesian_meshes
    use sll_m_arbitrary_degree_splines
    use sll_m_boundary_condition_descriptors
    use sll_m_arbitrary_degree_spline_interpolator_1d
    use sll_m_fft
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
    type :: poisson_1d_fd
        class(sll_cartesian_mesh_1d), private, pointer  :: cartesian_mesh

        !> Spline basis
        !sll_int32 , private :: fd_degree !<Degree of the Bsplines
        sll_int32 , private :: fd_degree !<Order of the finite differences
        sll_int32 , private :: spline_degree !<Degree of the B-splines
        
        
        type(arbitrary_degree_spline_1d), pointer, private :: bspline  !<Bspline object
        class(sll_arbitrary_degree_spline_interpolator_1d), pointer, private :: interpolator


        !Finite difference stencils with weights
        sll_real64, dimension(:), allocatable :: stenw_dx    !gradient 
        sll_real64, dimension(:), allocatable :: stenw_dxdx  !Laplace
        
        sll_real64, dimension(:), allocatable :: fd_solution
        sll_real64, dimension(:), allocatable,private :: fd_solution_deriv

        !For Precomputation do one FFT of fd_matrix_first_line
        sll_real64, dimension(:), allocatable, private :: fd_matrix_first_line
        sll_comp64, dimension(:), allocatable, private :: fd_matrix_first_line_fourier

        
        
        !FFT Plans
        type(sll_fft_plan), pointer, private :: forward_fftplan => null()
        type(sll_fft_plan), pointer, private :: backward_fftplan => null()
        sll_int32,private ::num_cells
        sll_int32,private ::problem_size !Problemsize

        sll_int32,private :: boundarycondition
    contains
        procedure,  pass(this) :: initialize =>sll_initialize_poisson_1d_fd
        !procedure,  pass(poisson) :: new =>initialize_poisson_1d_fd
        !procedure,  pass(poisson) :: initialize =>initialize_poisson_1d_fd
        procedure,  pass(this) :: delete =>sll_delete_poisson_1d_fd
        procedure,  pass(this) :: solve =>solve_poisson_1d_fd_rhs
        procedure,  pass(this) :: eval_solution=>poisson_1d_fd_eval_solution
        procedure,  pass(this) :: eval_solution_derivative=>poisson_1d_fd_eval_solution_derivative
        procedure, pass(this) :: H1seminorm_solution=>poisson_1d_fd_H1seminorm_solution
        procedure, pass(this) :: L2norm_solution=>poisson_1d_fd_L2norm_solution
        procedure, pass(this) :: interp=>poisson_1d_fd_interp
        procedure, pass(this) :: set_solution=>poisson_1d_fd_set_solution
        procedure, pass(this) :: get_rhs_klimontovich=>poisson_1d_fd_get_rhs_from_klimontovich_density_weighted
        procedure, pass(this),private :: set_stencil=>poisson_1d_fd_set_stencil
        procedure,pass(this), public :: get_rhs_fun=>poisson_1d_fd_get_rhs_fun

    end type poisson_1d_fd

    !Interface for one dimensional function for right hand side
    abstract interface
        function poisson_1d_fd_rhs_function(x) result(y)
            use sll_m_working_precision
            sll_real64, dimension(:),intent(in) :: x
            sll_real64, dimension(size(x)) :: y
        endfunction
    endinterface



    !Interface for one dimensional function for right hand side
    abstract interface
        function poisson_1d_fd_particle_shape_function(ppos, nodes) result(nodevals)
            use sll_m_working_precision
            sll_real64, dimension(:),intent(in) :: nodes
            sll_real64, intent(in) :: ppos
            sll_real64, dimension(size(nodes)) :: nodevals
        endfunction
    endinterface



    !     interface solve_poisson_1d_fd
    !                  module procedure solve_poisson_1d_fd_rhs
    !                  module procedure solve_poisson_1d_fd_functionrhs
    !                  module procedure solve_poisson_1d_fd_rhs_and_get
    !    end interface solve_poisson_1d_fd
    interface delete
        module procedure sll_delete_poisson_1d_fd
    endinterface

    interface new
         module procedure new_poisson_1d_fd
    endinterface

contains

    !>Destructor
    subroutine sll_delete_poisson_1d_fd(this,ierr)
        class(poisson_1d_fd),intent(inout) :: this     !< Solver data structure
        sll_int32, intent(out)                :: ierr    !< error code
        SLL_DEALLOCATE_ARRAY(this%fd_matrix_first_line,ierr)
        SLL_DEALLOCATE_ARRAY(this%fd_solution,ierr)
        SLL_DEALLOCATE_ARRAY(this%fd_matrix_first_line_fourier,ierr)
        call fft_delete_plan(this%backward_fftplan)
        call fft_delete_plan(this%forward_fftplan)
    endsubroutine


    function  new_poisson_1d_fd(cartesian_mesh_1d,fd_degree, sp_degree,bc_type,ierr) &
            result(solver)
        type(poisson_1d_fd), pointer :: solver     !< Solver data structure
        class(sll_cartesian_mesh_1d), intent(in),pointer  :: cartesian_mesh_1d !< Logical mesh
        sll_int32, intent(out)                :: ierr    !< error code
        sll_int32, intent(in)                :: fd_degree !<Degree of the finite differences approximation
        sll_int32, intent(in)                :: sp_degree !<Degree of the bsplines

        sll_int32, intent(in)               :: bc_type !< type of boundary connditions

        SLL_ALLOCATE(solver,ierr)
        call sll_initialize_poisson_1d_fd(solver,  cartesian_mesh_1d,fd_degree,sp_degree,bc_type,ierr)
    endfunction


    subroutine sll_initialize_poisson_1d_fd(this, cartesian_mesh_1d,fd_degree,sp_degree, bc_type, ierr)
        class(poisson_1d_fd),intent(inout) :: this     !< Solver data structure
        class(sll_cartesian_mesh_1d), intent(in),pointer  :: cartesian_mesh_1d !< Logical mesh
        sll_int32, intent(out)                :: ierr    !< error code
        sll_int32, intent(in)                :: sp_degree !<Degree of the bsplines
        sll_int32, intent(in)                :: fd_degree !<Degree of the finite differences

        sll_int32, intent(in)               :: bc_type !< type of boundary connditions
        ierr=0
        this%boundarycondition=bc_type !SLL_PERIODIC

        !scale_matrix_equation=1.0_f64
        this%cartesian_mesh=>cartesian_mesh_1d
        this%num_cells=this%cartesian_mesh%num_cells
        this%spline_degree=sp_degree
        this%fd_degree=fd_degree
        
        !Allocate spline
        selectcase(this%boundarycondition)
            case(SLL_PERIODIC)
                this%bspline=>new_arbitrary_degree_spline_1d(this%spline_degree,this%cartesian_mesh%eta1_nodes(), this%cartesian_mesh%num_nodes(), &
                    PERIODIC_ARBITRARY_DEG_SPLINE)
                this%problem_size=this%num_cells
                this%interpolator=>new_arbitrary_degree_1d_interpolator(this%num_cells, &
                    this%cartesian_mesh%eta_min, &
                    this%cartesian_mesh%eta_max, &
                    SLL_PERIODIC, &
                    SLL_PERIODIC, &
                    1)

            case(SLL_DIRICHLET)
                this%problem_size=this%num_cells
                this%interpolator=>new_arbitrary_degree_1d_interpolator( this%num_cells, &
                    this%cartesian_mesh%eta_min, &
                    this%cartesian_mesh%eta_max, &
                    SLL_DIRICHLET, &
                    SLL_DIRICHLET, &
                    1)
        endselect

        SLL_ALLOCATE(this%fd_solution(this%problem_size),ierr)
        SLL_ALLOCATE(this%fd_solution_deriv(this%problem_size),ierr)

        !FFT--------------------------------------------------------------------------------------
        !To get the same output as in the MATLAB example use
        SLL_ALLOCATE(this%fd_matrix_first_line_fourier(1:this%num_cells/2+1), ierr)
        this%fd_matrix_first_line_fourier = (0._f64,0._f64)
        this%forward_fftplan => fft_new_plan(this%num_cells,this%fd_solution,this%fd_matrix_first_line_fourier, &
            FFT_FORWARD+FFT_NORMALIZE)
        this%backward_fftplan=>fft_new_plan(this%num_cells,this%fd_matrix_first_line_fourier,this%fd_solution, &
            FFT_INVERSE  + FFT_NORMALIZE_INVERSE)
        SLL_DEALLOCATE_ARRAY(this%fd_matrix_first_line_fourier, ierr)

        !------------------------------------------------------------------------------------------
        !prepare fourier transformation for stiffness matrix

        SLL_ALLOCATE(this%fd_matrix_first_line_fourier(1:this%num_cells/2+1), ierr)
        this%fd_matrix_first_line_fourier = (0._f64,0._f64)

        call sll_poisson_1d_fd_assemble_fd_matrix(this,ierr)

        call sll_poisson_1d_fd_fft_precomputation(this,this%fd_matrix_first_line, &
            this%fd_matrix_first_line_fourier, ierr)
    endsubroutine

    subroutine poisson_1d_fd_set_solution(this, solution_vector)
        class(poisson_1d_fd),intent(inout) :: this     !< Solver data structure
        sll_real64, dimension(:) :: solution_vector
        SLL_ASSERT(size(solution_vector)==this%num_cells)

        this%fd_solution=solution_vector
    endsubroutine

    !<Calculates the inhomogenity rhs with given function f
    !<by Gauss Legendre integration
    function sll_poisson_1d_fd_get_rhs_from_function(this, eval_function) &
            result( rhs )
        implicit none
        class(poisson_1d_fd),intent(in) :: this     !< Solver data structure
        procedure (poisson_1d_fd_rhs_function) :: eval_function
        sll_real64, dimension(this%num_cells ) :: rhs !<Right hand side
        !sll_int32 :: ierr=0
        sll_real64, dimension(this%num_cells+1) :: evalpoints
        evalpoints=eval_function(this%cartesian_mesh%eta1_nodes())
        !Map to right hand side according to boundary condition

        rhs=evalpoints(1:this%num_cells)
        selectcase(this%boundarycondition)
            case(SLL_PERIODIC)
            case(SLL_DIRICHLET)
                rhs(1)=0._f64
        endselect

    endfunction



    !<Does precomputation for the FFT solver
    subroutine sll_poisson_1d_fd_fft_precomputation(this,circulant_matrix_first_line, &
            circulant_matrix_first_line_fourier,ierr)
        class(poisson_1d_fd),intent(inout) :: this     !< Solver data structure
        sll_int32, intent(out) :: ierr
        sll_real64, dimension(:),allocatable, intent(in)   :: circulant_matrix_first_line
        sll_comp64, dimension(this%num_cells/2 +1),intent(out)  ::circulant_matrix_first_line_fourier
        sll_real64, dimension(:) :: circulantvector(size(circulant_matrix_first_line))
        integer :: N
        ierr=0
        !Determine dimension of problem
        N=size(circulant_matrix_first_line)
        !!!SLL_ASSERT(is_power_of_two(int(N,i64)))
        !Generate Circulant Seed

        !Circulant seed c to generate stiffness matrix as a circular matrix
        !Remember c is not the first line of the circulant matrix
        !this is the first line of the matrix
        ! (c_1  0 0  0  .... c_4   c_3  c_2)
        !Note this only works because the period is symmetric
        circulantvector=0._f64
        !circulantvector(1)=circulant_matrix_first_line(1)
        !circulantvector(2:N)=circulant_matrix_first_line(N:2:-1)
!        circulantvector=circulantvector
        circulantvector=circulant_matrix_first_line

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


    subroutine sll_poisson_1d_fd_assemble_fd_matrix(this,ierr)
        class(poisson_1d_fd),intent(inout) :: this     !< Solver data structure
        sll_int32, intent(out)                :: ierr    !< error code

        !sll_real64, dimension(this%fd_degree*3) :: period
        !sll_int32 :: periodlen

        SLL_CLEAR_ALLOCATE(this%fd_matrix_first_line(1:this%num_cells),ierr)
        this%fd_matrix_first_line=0.0_f64
        !selectcase(this%fd_degree)
            !case(1)
                !                  period=(/ 1.0_f64, 2.0_f64, 1.0_f64 /)
                !                 periodlen=3
                this%fd_matrix_first_line(1:2:1)=(/ -2.0_f64, 1.0_f64 /)
                this%fd_matrix_first_line(this%num_cells)=1.0_f64
         !   case default
          !      print *, "Higher order Finite Differences are not implemented!"
           !     stop
        !endselect

    endsubroutine

    subroutine solve_poisson_1d_fd_rhs_and_get(this, field, rhs)
        class(poisson_1d_fd),intent(inout) :: this
        sll_real64, dimension(:), intent(in)      :: rhs
        sll_real64, dimension(this%num_cells), intent(out)     :: field

        call  solve_poisson_1d_fd_rhs(this, rhs)
        field=this%fd_solution

    end subroutine

    subroutine solve_poisson_1d_fd_rhs(this, rhs)
        class(poisson_1d_fd),intent(inout) :: this
        sll_real64, dimension(:), intent(in)      :: rhs
        sll_int32 :: N
        SLL_ASSERT(size(rhs)==this%num_cells)
            N=size(rhs)
        
        this%fd_solution=poisson_1d_fd_solve_circulant_matrix_equation(this, &
            this%fd_matrix_first_line_fourier ,-rhs )


        if (  this%boundarycondition==SLL_DIRICHLET ) then
            this%fd_solution(1)=0._f64
        endif

            this%fd_solution=this%fd_solution*this%cartesian_mesh%delta_eta*2
        
        
        !Calculate the Electric field from the Potential
        !Central Difference
        this%fd_solution_deriv=(cshift(this%fd_solution,-1)- &
                                        cshift(this%fd_solution,1))&
                                     /(-2*(this%cartesian_mesh%delta_eta))

        !Set up Interpolation
        !call this%interpolator%compute_interpolants(this%fd_solution)

    endsubroutine


    !> @brief Solves the poisson equation for a given right hand side function
    !> @param this pointer to a poisson_1d_fd object.
    !> @param rhs right hand side function of type poisson_1d_fd_rhs_function
    subroutine solve_poisson_1d_fd_functionrhs(this, rhs_fun)
        implicit none
        class(poisson_1d_fd),intent(inout) :: this
        procedure (poisson_1d_fd_rhs_function) :: rhs_fun
        sll_real64, dimension(this%num_cells) :: rhs

        rhs=sll_poisson_1d_fd_get_rhs_from_function(this, rhs_fun)
        call solve_poisson_1d_fd_rhs(this,rhs  )
    endsubroutine


    function poisson_1d_fd_solve_circulant_matrix_equation(this, matrix_fl_fourier,rightside ) &
            result(solution)
        class(poisson_1d_fd),intent(inout) :: this
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
        !!!SLL_ASSERT(is_power_of_two(int(N,i64)))

        SLL_ASSERT(N/2+1==size(matrix_fl_fourier))

        !Determine dimension of problem

        SLL_ALLOCATE(constant_factor_fourier(1:N/2+1), ierr)
        constant_factor_fourier = (0._f64,0._f64)
        SLL_ALLOCATE(data_complex(1:N/2+1),ierr)
        data_complex = (0._f64,0._f64)
        constant_factor=0._f64
        constant_factor=rightside

        call fft_apply_plan(this%forward_fftplan,constant_factor,constant_factor_fourier)
        !constant_factor_fourier(1)=0
        data_complex=constant_factor_fourier/(matrix_fl_fourier)
        data_complex(1)=(0._f64,0._f64)
        SLL_DEALLOCATE_ARRAY(constant_factor_fourier,ierr)

        call fft_apply_plan(this%backward_fftplan,data_complex,solution)
        SLL_DEALLOCATE_ARRAY(data_complex,ierr)
        !Somehow the normalization does not do what it should do:
        solution=solution/N
    endfunction


    !< knots are the interpolant points and have nothing to do with the mesh
    function poisson_1d_fd_bspline_basis_to_realvals ( this , bspline_vector, knots_eval_in, interpolfun_user ) &
            result(realvals)
        class(poisson_1d_fd),intent(in) :: this
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

        realvals=0._f64

        cell=sll_cell(this%cartesian_mesh, knots_eval)
        !cell= bspline_fem_solver_1d_cell_number(knots_mesh, knots_eval(eval_idx))
        !Loop over all points to evaluate
        !This should be vectorizzed
        do eval_idx=1,int(size(knots_eval),i64)
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
    subroutine poisson_1d_fd_eval_solution(this, knots_eval, eval_solution)
        class(poisson_1d_fd),intent(inout) :: this
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        !sll_real64, dimension(this%bspline%degree+1) :: b_contribution
        !sll_real64 :: celldx
        !sll_int :: b_idx
        !sll_int :: idx
        !sll_int :: b_contrib_idx
      
        SLL_ASSERT(size(knots_eval)==size(eval_solution))

        call this%interp( knots_eval, eval_solution, this%fd_solution)
    endsubroutine

    !< Evaluates the solution at the given points knots_eval
    !< The Result is written into eval_solution
    subroutine poisson_1d_fd_eval_solution_derivative(this, knots_eval, eval_solution)
        class(poisson_1d_fd),intent(inout) :: this
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        SLL_ASSERT(size(knots_eval)==size(eval_solution))

        call this%interp( knots_eval, eval_solution, this%fd_solution_deriv)
       
    endsubroutine

    !    !<Gives the squared H1-seminorm of the solution $\Phi$: $|\nabla \Phi|^2$
        function  poisson_1d_fd_H1seminorm_solution(this) result(seminorm)
            class(poisson_1d_fd),intent(inout) :: this
            sll_real64 :: seminorm
            sll_int32 :: N
            N=size(this%fd_solution)

            !This is a really bad way to do it
            seminorm=sum(this%fd_solution_deriv**2)/real(N**2,8) !&
                            !/(this%cartesian_mesh%eta_max -this%cartesian_mesh%eta_min))/N
                            !+(this%fd_solution(1:N-1)-this%fd_solution(2:N))
        endfunction
!
        !<Gives the squared L2-norm of the solution $\Phi$: $|\Phi|^2$
        function  poisson_1d_fd_L2norm_solution(this) result(l2norm)
            class(poisson_1d_fd),intent(inout) :: this
            sll_real64 :: l2norm

            l2norm=sum(this%fd_solution**2)
        endfunction


    function poisson_1d_fd_get_rhs_from_klimontovich_density(this, &
            ppos)  result(rhs)
        class(poisson_1d_fd),intent(inout) :: this
        sll_real64, dimension(:), intent(in) ::ppos
        !procedure (b_splines_at_x), optional:: interpolfun_user
        sll_real64, dimension(this%num_cells) :: rhs

        sll_real64, dimension(1) :: pweight
        pweight=1.0_f64

        !        if (present(interpolfun_user)) then
        !            rhs=poisson_1d_fd_get_rhs_from_klimontovich_density_weighted(this,&
            !                ppos,pweight, interpolfun_user)
        !        else
        rhs=poisson_1d_fd_get_rhs_from_klimontovich_density_weighted(this,&
            ppos,pweight)
        !        endif
    endfunction


    function poisson_1d_fd_get_rhs_particle_shape_function( this, shf, ppos, pweight) result(rhs)
        class(poisson_1d_fd),intent(in) :: this
        procedure (poisson_1d_fd_particle_shape_function) :: shf
        sll_real64, dimension(:), intent(in) ::ppos
        sll_real64, dimension(:), intent(in) ::pweight
        sll_real64, dimension(this%problem_size) :: rhs
        sll_real64, dimension(this%num_cells+1) :: knotsvals
        sll_real64, dimension(this%num_cells+1) :: knots
        sll_int32 :: idx,npart

        npart=size(ppos)
        SLL_ASSERT(size(ppos)==size(pweight))

        knots=this%cartesian_mesh%eta1_nodes()
        knotsvals=0._f64
        do idx=1,npart
            knotsvals=shf(ppos(idx), knots)*pweight(idx)
        enddo

        rhs=0._f64
        selectcase(this%boundarycondition)
            case(SLL_PERIODIC)
                !Last knot is identical with the first knot
                SLL_ASSERT(this%problem_size==this%num_cells)
                rhs=knotsvals(1:this%num_cells)
                rhs(1)=rhs(1)+knotsvals(this%num_cells+1)

            case(SLL_DIRICHLET)
                SLL_ASSERT(this%problem_size==this%num_cells)
                rhs=knotsvals(1:this%num_cells)
                rhs(1)=0.0_f64 !homogenous dirichlet
                !works only because we have first order scheme
        endselect


    endfunction

subroutine poisson_1d_fd_interp(this, knots_eval, eval_solution , solution)
    class(poisson_1d_fd),intent(inout) :: this
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        sll_real64, dimension(:), intent(in) :: solution
        sll_real64, dimension(this%bspline%degree+1) :: b_contribution
        sll_real64 :: celldx
        sll_int :: b_idx, idx!, jdx
        sll_int :: b_contrib_idx
      
        SLL_ASSERT(size(knots_eval)==size(eval_solution))
         print *, b_splines_at_x(this%bspline,1, this%cartesian_mesh%delta_eta*0)

         eval_solution=0._f64
 
        do idx=1,size(knots_eval)

        
           celldx=(knots_eval(idx)+this%cartesian_mesh%delta_eta*(this%bspline%degree+1)/2.0_f64)
           celldx=this%cartesian_mesh%period(celldx) !periodicity
           
           b_idx=floor(celldx/this%cartesian_mesh%delta_eta)
           celldx=mod(celldx,this%cartesian_mesh%delta_eta)
           
           
           b_contribution=b_splines_at_x(this%bspline,1, celldx)
           b_contribution=b_contribution(this%bspline%degree+1:1:-1)

            !First spline has first support in first cell of interval
            do b_contrib_idx=1, this%bspline%degree+1

                
                if (b_idx > this%num_cells) then
                    b_idx = b_idx - this%num_cells
                elseif (b_idx < 1) then
                    b_idx = b_idx + this%num_cells
                endif
                
                eval_solution(idx)=eval_solution(idx)+  b_contribution(b_contrib_idx)*solution(b_idx)
                b_idx=b_idx+1
            enddo
        enddo
    
    !eval_solution=-eval_solution
    endsubroutine

    function poisson_1d_fd_get_rhs_from_klimontovich_density_weighted( this,&
            ppos, pweight) result(rhs)
        class(poisson_1d_fd),intent(inout) :: this
        sll_real64, dimension(:), intent(in) ::ppos
        sll_real64, dimension(:), intent(in) ::pweight
        sll_real64, dimension(this%num_cells) :: rhs
         sll_int32 :: idx
        !sll_int32, dimension( size(ppos) ) :: cell
        sll_real64, dimension(this%bspline%degree+1) :: b_contribution
        sll_real64 :: celldx
        sll_int :: b_idx
        sll_int :: b_contrib_idx
        rhs=0._f64

        SLL_ASSERT( (size(pweight)==size(ppos) .OR. size(pweight)==1))
 

        do idx=1,size(ppos)
           celldx=(ppos(idx)+this%cartesian_mesh%delta_eta*(this%bspline%degree+1)/2.0_f64)
           celldx=this%cartesian_mesh%period(celldx) !periodicity
           
           b_idx=floor(celldx/this%cartesian_mesh%delta_eta)
           celldx=mod(celldx,this%cartesian_mesh%delta_eta)
           
           
           b_contribution=b_splines_at_x(this%bspline,1, celldx)
           b_contribution=b_contribution(this%bspline%degree+1:1:-1)
           

            if (size(pweight)==1) then
                b_contribution=b_contribution*pweight(1)
            else
                b_contribution=b_contribution*pweight(idx)
            endif

            !First spline has first support in first cell of interval
            do b_contrib_idx=1, this%bspline%degree+1

                
                if (b_idx > this%num_cells) then
                    b_idx = b_idx - this%num_cells
                elseif (b_idx < 1) then
                    b_idx = b_idx + this%num_cells
                endif
                
                rhs(b_idx)=rhs(b_idx)+b_contribution(b_contrib_idx)
                b_idx=b_idx+1
            enddo
        enddo
        
        !normalize
        rhs=rhs/this%cartesian_mesh%delta_eta  
       
    endfunction


    subroutine poisson_1d_fd_linear_interpol(this, solution, knots_eval, eval_solution )
                                           !periodic
        class(poisson_1d_fd),intent(inout) :: this
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        sll_real64, dimension(:), intent(in)     :: solution
        sll_int32 :: idx, idx_r, idx_l
        sll_real64, dimension(size(knots_eval)) :: pw
        sll_real64 :: weight
        SLL_ASSERT(size(knots_eval)==size(eval_solution))

        pw=(knots_eval-this%cartesian_mesh%eta_min)/this%cartesian_mesh%delta_eta

        SLL_ASSERT(size(solution)==this%problem_size)
            do idx=1,size(eval_solution)
                    idx_l=floor(pw(idx)) +1
                    idx_r=mod(idx_l,this%problem_size)+1
                    weight=pw(idx)-idx_l
                eval_solution(idx)= eval_solution(idx) &
                                        +(1.0_f64- weight)*solution(idx_l ) &
                                        +(weight)*solution(idx_r )
            enddo

    endsubroutine
 
 
    function poisson_1d_fd_get_rhs_fun(this, eval_function)      result( rhs )
        class(poisson_1d_fd),intent(in) :: this     !< Solver data structure
        procedure (poisson_1d_fd_rhs_function) :: eval_function
        sll_real64, dimension(this%num_cells ) :: rhs !<Right hand side
        !sll_int32 :: idx 
        sll_real64, dimension(this%num_cells +1)  :: knots
        
        
        knots=this%cartesian_mesh%eta1_nodes()
        rhs=eval_function(knots(1:this%num_cells))
    endfunction


    subroutine poisson_1d_fd_set_stencil(this, order)
            class(poisson_1d_fd),intent(inout) :: this     !< Solver data structure
        sll_int32, intent(in) :: order
        sll_int32:: ierr
        SLL_ALLOCATE(this%stenw_dx(1:order+1),ierr)
        SLL_ALLOCATE(this%stenw_dxdx(1:order+1),ierr)
        
 
 SELECT CASE (order)
 CASE (2)
      ! Order = 2 
this%stenw_dx=(/-1._f64/2._f64, 0._f64, 1._f64/2._f64/)
this%stenw_dxdx=(/1._f64, -2._f64, 1._f64/)

CASE (4)
  ! Order = 4
this%stenw_dx=(/1._f64/12._f64, -2._f64/3._f64, 0._f64, 2._f64/3._f64, -1._f64/12._f64/)
this%stenw_dxdx=(/-1._f64/12._f64, 4._f64/3._f64, -5._f64/2._f64, 4._f64/3._f64, -1._f64/12._f64/)
    
   CASE (6)
    ! Order = 6
this%stenw_dx=(/-1._f64/60._f64, 3._f64/20._f64, -3._f64/4._f64, 0._f64, 3._f64/4._f64, -3._f64/20._f64, 1._f64/60._f64/)
this%stenw_dxdx=(/1._f64/90._f64, -3._f64/20._f64, 3._f64/2._f64, -49._f64/18._f64, 3._f64/2._f64, -3._f64/20._f64, 1._f64/90._f64/)

CASE (8)
! Order = 8
this%stenw_dx=(/1._f64/280._f64, -4._f64/105._f64, 1._f64/5._f64, -4._f64/5._f64, 0._f64, 4._f64/5._f64, -1._f64/5._f64, 4._f64/105._f64, -1._f64/280._f64/)
this%stenw_dxdx=(/-1._f64/560._f64, 8._f64/315._f64, -1._f64/5._f64, 8._f64/5._f64, -205._f64/72._f64, 8._f64/5._f64, -1._f64/5._f64, 8._f64/315._f64, -1._f64/560._f64/)

CASE(10)
! Order = 10
this%stenw_dx=(/-1._f64/1260._f64, 5._f64/504._f64, -5._f64/84._f64, 5._f64/21._f64, -5._f64/6._f64, 0._f64, 5._f64/6._f64, -5._f64/21._f64, 5._f64/84._f64, -5._f64/504._f64, 1._f64/1260._f64/)
this%stenw_dxdx=(/1._f64/3150._f64, -5._f64/1008._f64, 5._f64/126._f64, -5._f64/21._f64, 5._f64/3._f64, -5269._f64/1800._f64, 5._f64/3._f64, -5._f64/21._f64, 5._f64/126._f64, -5._f64/1008._f64, 1._f64/3150._f64/)

CASE(12)
! Order = 12
this%stenw_dx=(/1._f64/5544._f64, -1._f64/385._f64, 1._f64/56._f64, -5._f64/63._f64, 15._f64/56._f64, -6._f64/7._f64, 0._f64, 6._f64/7._f64, -15._f64/56._f64, 5._f64/63._f64, -1._f64/56._f64, 1._f64/385._f64, -1._f64/5544._f64/)
this%stenw_dxdx=(/-1._f64/16632._f64, 2._f64/1925._f64, -1._f64/112._f64, 10._f64/189._f64, -15._f64/56._f64, 12._f64/7._f64, -5369._f64/1800._f64, 12._f64/7._f64, -15._f64/56._f64, 10._f64/189._f64, -1._f64/112._f64, 2._f64/1925._f64, -1._f64/16632._f64/)
   CASE DEFAULT
      
END SELECT
    end subroutine
    
    
    !    function poisson_1d_fd_calculate_residual(this) result(residuum)
    !                class(poisson_1d_fd),intent(inout) :: this
    !
    !        sll_real64 :: residuum
    !        residuum= sqrt(sum(( &
        !            poisson_1d_fd_circulant_matrix_vector_product&
        !            ( this%fd_matrix_first_line, this%fd_solution)&
        !            - this%fem_inhomogenity &
        !            )**2)) !/sqrt(sum(fem_inhomogenity**2 ))
    !    endfunction
    !
    !
    !
    !
    !    !<can be optimized, with a sparse first line
    !    function  poisson_1d_fd_circulant_matrix_vector_product( circulant_matrix_first_line, vector ) &
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

end module sll_m_poisson_1d_fd
