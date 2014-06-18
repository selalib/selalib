
!> Module to solve Poisson equation on one dimensional mesh using Finite Elements
module sll_poisson_1d_fourier
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

    use sll_constants
    use sll_logical_meshes
    use sll_arbitrary_degree_splines
    use sll_boundary_condition_descriptors
    use sll_arbitrary_degree_spline_interpolator_1d_module
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
    type :: poisson_1d_fourier
        type(sll_logical_mesh_1d), private, pointer  :: logical_mesh

        !> Fourier modes
        sll_int32 , private :: num_modes !<Number of fourier modes
        sll_real64, private :: Ilength
        sll_real64, private :: eta_min

        sll_comp64,  dimension(:), allocatable :: fourier_fmode


    contains
        procedure,  pass(this) :: initialize =>sll_initialize_poisson_1d_fourier
        !procedure,  pass(poisson) :: new =>initialize_poisson_1d_fourier
        !procedure,  pass(poisson) :: initialize =>initialize_poisson_1d_fourier
        procedure,  pass(this) :: delete =>sll_delete_poisson_1d_fourier
        procedure,  pass(this) :: solve =>solve_poisson_1d_fourier_rhs
        procedure,  pass(this) :: eval_solution=>poisson_1d_fourier_eval_solution
        procedure,  pass(this) :: eval_solution_derivative=>poisson_1d_fourier_eval_solution_derivative

        procedure, pass(this) :: H1seminorm_solution=>poisson_1d_fourier_H1seminorm_solution
        procedure, pass(this) :: L2norm_solution=>poisson_1d_fourier_L2norm_solution
        procedure, pass(this) ::  get_rhs_from_klimontovich_density_weighted=>&
                                            poisson_1d_fourier_get_rhs_from_klimontovich_density_weighted
    end type poisson_1d_fourier

    !Interface for one dimensional function for right hand side
    abstract interface
        function poisson_1d_fourier_rhs_function(x) result(y)
            use sll_working_precision
            sll_real64, dimension(:),intent(in) :: x
            sll_real64, dimension(size(x)) :: y
        endfunction
    endinterface


    !     interface solve_poisson_1d_fourier
    !                  module procedure solve_poisson_1d_fourier_rhs
    !                  module procedure solve_poisson_1d_fourier_functionrhs
    !                  module procedure solve_poisson_1d_fourier_rhs_and_get
    !    end interface solve_poisson_1d_fourier
    interface delete
        module procedure sll_delete_poisson_1d_fourier
    endinterface

    !    interface new
    !         module procedure new_poisson_1d_fourier
    !    endinterface

contains

    !>Destructor
    subroutine sll_delete_poisson_1d_fourier(this,ierr)
        class(poisson_1d_fourier),intent(inout) :: this     !< Solver data structure
        sll_int32, intent(out)                :: ierr    !< error code
        SLL_DEALLOCATE_ARRAY(this%fourier_fmode,ierr)

    endsubroutine


    function  new_poisson_1d_fourier(logical_mesh_1d,num_modes, bc_type,ierr) &
            result(solver)
        type(poisson_1d_fourier), pointer :: solver     !< Solver data structure
        type(sll_logical_mesh_1d), intent(in),pointer  :: logical_mesh_1d !< Logical mesh
        sll_int32, intent(out)                :: ierr    !< error code
        sll_int32, intent(in)                :: num_modes !<Degree of the finite differences approximation
        sll_int32, intent(in)               :: bc_type !< type of boundary connditions

        SLL_ALLOCATE(solver,ierr)
        call sll_initialize_poisson_1d_fourier(solver,  logical_mesh_1d,num_modes,bc_type,ierr)
    endfunction


    subroutine sll_initialize_poisson_1d_fourier(this, logical_mesh_1d,num_modes,bc_type, ierr)
        class(poisson_1d_fourier),intent(inout) :: this     !< Solver data structure
        type(sll_logical_mesh_1d), intent(in),pointer  :: logical_mesh_1d !< Logical mesh
        sll_int32, intent(out)                :: ierr    !< error code
        sll_int32, intent(in)                :: num_modes !<Degree of the bsplines
        sll_int32, intent(in)               :: bc_type !< type of boundary connditions
        ierr=0
        !this%boundarycondition=bc_type !SLL_PERIODIC

        this%logical_mesh=>logical_mesh_1d
        this%Ilength=this%logical_mesh%eta_max-this%logical_mesh%eta_min
        this%eta_min=this%logical_mesh%eta_min

        this%num_modes=num_modes

        !        selectcase(this%boundarycondition)
        !            case(SLL_PERIODIC)
        !
        !
        !            case(SLL_DIRICHLET)
        !
        !        endselect
        SLL_ALLOCATE(this%fourier_fmode(this%num_modes),ierr)


    endsubroutine

    !    subroutine poisson_1d_fourier_set_solution(this, solution_vector)
    !        class(poisson_1d_fourier),intent(inout) :: this     !< Solver data structure
    !        sll_real64, dimension(:) :: solution_vector
    !        SLL_ASSERT(size(solution_vector)==this%num_cells)
    !
    !        this%fd_solution=solution_vector
    !    endsubroutine

    !    !<Calculates the inhomogenity rhs with given function f
    !    !<by Gauss Legendre integration
    !    function sll_poisson_1d_fourier_get_rhs_from_function(this, eval_function) &
        !            result( rhs )
    !        implicit none
    !        class(poisson_1d_fourier),intent(in) :: this     !< Solver data structure
    !        procedure (poisson_1d_fourier_rhs_function) :: eval_function
    !        sll_real64, dimension(this%num_cells ) :: rhs !<Right hand side
    !        sll_int32 :: ierr=0
    !        sll_real64, dimension(this%num_cells+1) :: evalpoints
    !        evalpoints=eval_function(nodes_logical_mesh_1d(this%logical_mesh))
    !        !Map to right hand side according to boundary condition
    !
    !        rhs=evalpoints(1:this%num_cells)
    !        selectcase(this%boundarycondition)
    !            case(SLL_PERIODIC)
    !            case(SLL_DIRICHLET)
    !                rhs(1)=0
    !        endselect
    !
    !    endfunction






    subroutine solve_poisson_1d_fourier_rhs_and_get(this, field, rhs)
        class(poisson_1d_fourier),intent(inout) :: this
        sll_comp64, dimension(:), intent(in)      :: rhs
        sll_comp64, dimension(this%num_modes), intent(out)  :: field
        sll_real64 :: coeff
        coeff=2.0_f64*sll_pi/this%Ilength

        call  solve_poisson_1d_fourier_rhs(this, rhs)

        field=-this%fourier_fmode/coeff/solve_poisson_1d_fourier_get_modes(this)/sll_i1

    end subroutine


    function solve_poisson_1d_fourier_get_modes(this)  result(modes)
        class(poisson_1d_fourier),intent(inout) :: this
        sll_int32, dimension(this%num_modes) :: modes
        sll_int32 :: idx

        do idx=1,this%num_modes
            modes(idx)=idx
        enddo
    endfunction

    subroutine solve_poisson_1d_fourier_rhs(this, rhs)
        class(poisson_1d_fourier),intent(inout) :: this
        sll_comp64, dimension(:), intent(in)      :: rhs
        SLL_ASSERT(size(rhs)==this%num_modes)
        SLL_ASSERT(size(rhs)==size(this%fourier_fmode))

        !Nothing to be done here
        !rhs=this%fourier_fmode
        this%fourier_fmode=rhs
    endsubroutine

    !
    !    !> @brief Solves the poisson equation for a given right hand side function
    !    !> @param this pointer to a poisson_1d_fourier object.
    !    !> @param rhs right hand side function of type poisson_1d_fourier_rhs_function
    !    subroutine solve_poisson_1d_fourier_functionrhs(this, rhs_fun)
    !        implicit none
    !        class(poisson_1d_fourier),intent(inout) :: this
    !        procedure (poisson_1d_fourier_rhs_function) :: rhs_fun
    !        sll_real64, dimension(this%num_cells) :: rhs
    !
    !        rhs=sll_poisson_1d_fourier_get_rhs_from_function(this, rhs_fun)
    !        call solve_poisson_1d_fourier_rhs(this,rhs  )
    !    endsubroutine




    !< Evaluates the first derivative of the solution at the given points knots_eval
    !< The Result is written into eval_solution
    subroutine poisson_1d_fourier_eval_solution(this, knots_eval, eval_solution)
        class(poisson_1d_fourier),intent(inout) :: this
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        sll_int32, dimension(this%num_modes) :: fmode
        sll_int32 :: idx
        sll_real64 :: coeff
        coeff=2.0_f64*sll_pi/this%Ilength

        SLL_ASSERT(size(knots_eval)==size(eval_solution))
        do idx=1,this%num_modes
            fmode(idx)=idx
        enddo

        do idx=1, size(knots_eval)
                eval_solution(idx)= -2.0_f64* sum( &
                    real(  this%fourier_fmode/(coeff*fmode)**2 )*cos(knots_eval(idx)*fmode*coeff) &
                    - imag( this%fourier_fmode/(coeff*fmode)**2 )*sin(knots_eval(idx)*fmode*coeff))

        enddo

    endsubroutine

    !< Evaluates the solution at the given points knots_eval
    !< The Result is written into eval_solution
    subroutine poisson_1d_fourier_eval_solution_derivative(this, knots_eval, eval_solution)
        class(poisson_1d_fourier),intent(inout) :: this
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        sll_int32, dimension(this%num_modes) :: fmode
        sll_int32 :: idx
        sll_real64 :: coeff
        coeff=2.0_f64*sll_pi/this%Ilength

        SLL_ASSERT(size(knots_eval)==size(eval_solution))

        do idx=1,this%num_modes
            fmode(idx)=idx
        enddo
        do idx=1, size(knots_eval)
                eval_solution(idx)= -2.0_f64* sum( &
                    real(   this%fourier_fmode/(coeff*fmode*sll_i1) )*cos(knots_eval(idx)*fmode*coeff) &
                    - imag( this%fourier_fmode/(coeff*fmode*sll_i1) )*sin(knots_eval(idx)*fmode*coeff))

        enddo

    endsubroutine

    !    !<Gives the squared H1-seminorm of the solution $\Phi$: $|\nabla \Phi|^2$
    function  poisson_1d_fourier_H1seminorm_solution(this) result(seminorm)
        class(poisson_1d_fourier),intent(inout) :: this
        sll_real64 :: seminorm
        sll_int32 :: fmode_a, fmode_b
        sll_real64 :: coeff
        coeff=2.0_f64*sll_pi/this%Ilength

        seminorm=0
        do fmode_a=1,this%num_modes
            do fmode_b=1,this%num_modes
            seminorm=seminorm+ 2.0_f64*abs( (this%fourier_fmode(fmode_a)/coeff/fmode_a/sll_i1)&
                                *(this%fourier_fmode(fmode_b)/coeff/fmode_b/sll_i1))

            !seminorm=seminorm+ 2.0_f64*abs((this%fourier_fmode(fmode_a)/coeff/fmode_a/sll_i1 )**2)
            enddo
        enddo
    endfunction
    !
    !<Gives the squared L2-norm of the solution $\Phi$: $|\Phi|^2$
    function  poisson_1d_fourier_L2norm_solution(this) result(l2norm)
        class(poisson_1d_fourier),intent(inout) :: this
        sll_real64 :: l2norm
        sll_int32 :: fmode
        sll_real64 :: coeff
        coeff=2.0_f64*sll_pi/this%Ilength
        do fmode=1,this%num_modes
            l2norm=l2norm+ real((this%fourier_fmode(fmode)/(coeff*fmode*sll_i1)**2 )**2*2.0_f64)
        enddo
    endfunction


    function poisson_1d_fourier_get_rhs_from_klimontovich_density(this, &
            ppos)  result(rhs)
        class(poisson_1d_fourier),intent(inout) :: this
        sll_real64, dimension(:), intent(in) ::ppos
        sll_comp64, dimension(this%num_modes) :: rhs
        sll_int32 ::fmode
        do fmode=1,this%num_modes
            rhs(fmode)=sum(exp(-sll_i1 * fmode*ppos*2.0_f64*sll_pi/this%Ilength));
        enddo
    endfunction



    function poisson_1d_fourier_get_rhs_from_klimontovich_density_weighted( this,&
            ppos, pweight) result(rhs)
        class(poisson_1d_fourier),intent(inout) :: this
        sll_real64, dimension(:), intent(in) ::ppos
        sll_real64, dimension(:), intent(in) ::pweight
        sll_comp64, dimension(this%num_modes) :: rhs
        sll_int32 :: fmode

        SLL_ASSERT(size(ppos)==size(pweight))
        do fmode=1,this%num_modes
            rhs(fmode)=dot_product(exp(-sll_i1*  fmode*ppos*2.0_f64*sll_pi/this%Ilength), pweight );

        enddo
    endfunction



    !    function poisson_1d_fourier_calculate_residual(this) result(residuum)
    !                class(poisson_1d_fourier),intent(inout) :: this
    !
    !        sll_real64 :: residuum
    !        residuum= sqrt(sum(( &
        !            poisson_1d_fourier_circulant_matrix_vector_product&
        !            ( this%fd_matrix_first_line, this%fd_solution)&
        !            - this%fem_inhomogenity &
        !            )**2)) !/sqrt(sum(fem_inhomogenity**2 ))
    !    endfunction
    !
    !
    !
    !
    !    !<can be optimized, with a sparse first line
    !    function  poisson_1d_fourier_circulant_matrix_vector_product( circulant_matrix_first_line, vector ) &
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

end module sll_poisson_1d_fourier
