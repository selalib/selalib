#ifndef FEM_MASS_MATRIX_TOL
#define  FEM_MASS_MATRIX_TOL 0.0000000001_f64
#endif



module sll_bspline_fem_solver_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"
    !  use sll_arbitrary_degree_spline_interpolator_1d_module
    use sll_arbitrary_degree_splines
    use sll_boundary_condition_descriptors
    use gauss_legendre_integration
    !!use gauss_lobatto_integration
    use sll_fft
    use sll_constants
    use sll_collective
    use sll_module_poisson_1d_periodic_solver
    implicit none
    !initialize
    !solve
    !interpolate_solution
    !destroy

    sll_int32, parameter :: SLL_SOLVER_FEM = 1
    sll_int32, parameter :: SLL_SOLVER_FD = 2

    sll_int32, private :: sll_pic1d_poisson_solver=SLL_SOLVER_FEM
    sll_real64, DIMENSION(:), allocatable, private :: knots_mesh

    integer, private :: ierr !!!FIX THIS , DECIDE THIS

    sll_real64, private :: test_mode
    sll_int, private :: n_knots
    sll_int, private :: n_cells
    sll_int, protected :: n_steadyparticles
    type(arbitrary_degree_spline_1d), pointer, private :: bspline_arbitrary_degree
    integer, private :: spline_degree
    sll_real64, private :: scalar_inhomogenity

    sll_int, protected :: n_particles !number of particles
    sll_real64, dimension(:), allocatable, protected :: fem_inhomogenity_steady
    sll_real64, dimension(:), allocatable, private :: fem_inhomogenity
    !sll_real64, dimension(:), allocatable, private :: particleweights

    sll_real64, dimension(:), allocatable, private :: fem_solution
    sll_real64, dimension(:), allocatable, private :: stiffn_matrix_first_line

    !For the L2-Norm
    sll_real64, dimension(:), allocatable, private :: mass_matrix_first_line

    !For Precomputation do one FFT of stiffn_matrix_first_line

    sll_comp64 , dimension(:),allocatable, private :: circulant_first_line_fourier
    sll_real64, dimension(:), allocatable, private :: weights
    sll_real64, private :: scale_matrix_equation    !<Scale for the stiffnes matrix and the inhomogenity

    private :: gen_fem_matrix_period
    !    private :: sll_bspline_fem_solve_poisson_1d
    private :: interpolate_particles_bsplines
    private :: sll_bspline_fem_solver_1d_fft_precomputation
    private :: bspline_fem_solver_1d_solve_matrix_equation
    private :: ensure_boundaryc

    !Boundary Description handeled by the solver
    sll_int32, private :: sll_bspline_fem_solver_boundary_type = SLL_PERIODIC


    !Poisson Finite Differences Solver
    type (poisson_1d_periodic),pointer,private           :: poissonsolverFD => null()



    !Parallelization
    integer, private :: collective_rank, collective_size
    type(sll_collective_t), pointer , private:: collective

    !FFT Plans
    type(sll_fft_plan), pointer, private :: forward_fftplan => null()
    type(sll_fft_plan), pointer, private :: backward_fftplan => null()

    interface bspline_fem_solver_1d_cell_number
        module procedure         bspline_fem_solver_1d_cell_number_single
        module procedure     bspline_fem_solver_1d_cell_number_array
    endinterface


    !  integer, private :: startsbefore=0

contains
    subroutine sll_bspline_fem_solver_1d_destroy()
        SLL_DEALLOCATE_ARRAY(fem_inhomogenity_steady, ierr)
        SLL_DEALLOCATE_ARRAY(stiffn_matrix_first_line,ierr)
        SLL_DEALLOCATE_ARRAY(mass_matrix_first_line,ierr)
        SLL_DEALLOCATE_ARRAY(fem_solution,ierr)
        SLL_DEALLOCATE_ARRAY(circulant_first_line_fourier,ierr)
        SLL_DEALLOCATE_ARRAY(fem_inhomogenity, ierr)
        SLL_DEALLOCATE_ARRAY(knots_mesh,ierr)
        call fft_delete_plan(backward_fftplan)
        call fft_delete_plan(forward_fftplan)

        !parallel
        !        call sll_halt_collective( )

    endsubroutine

    function  sll_bspline_fem_solver_1d_get_solution() result(solution)
        sll_real64, dimension(size(fem_solution)) :: solution

        solution=fem_solution
    endfunction

    function  sll_bspline_fem_solver_1d_get_inhomogenity() result(inhomogenity)
        sll_real64, dimension(size(fem_solution)) :: inhomogenity

        inhomogenity=fem_inhomogenity
    endfunction

    function sll_bspline_fem_solver_1d_calculate_residual() result(residuum)
        sll_real64 :: residuum
        residuum= sqrt(sum(( &
            bspline_fem_solver_1d_circulant_matrix_vector_product&
            ( stiffn_matrix_first_line, fem_solution)&
            - fem_inhomogenity &
            )**2)) !/sqrt(sum(fem_inhomogenity**2 ))
    endfunction


    subroutine sll_bspline_fem_solver_1d_set_weights(weights_user )
        sll_real64, dimension(:), allocatable, intent(in):: weights_user
        if (.NOT. allocated(weights)) then
            SLL_CLEAR_ALLOCATE(weights(1:size(weights_user)),ierr)
        endif
        weights=weights_user
    endsubroutine

    subroutine sll_bspline_fem_solver_1d_get_weights(weights_user )
        sll_real64, dimension(:), allocatable, intent(out):: weights_user
        weights_user=weights
    endsubroutine

    !<The weights have to be set later
    !<No weighting for the steady background heavy ions
    subroutine sll_bspline_fem_solver_1d_initialize(knots, splinedegree, &
            steadyparticleposition,scalar_inhomogenity_user, collective_user)
        sll_real64, dimension(:), intent(in)  :: knots
        sll_real64, dimension(:), optional, intent(in)  :: steadyparticleposition
        sll_real64, optional, intent(in)   :: scalar_inhomogenity_user
        type(sll_collective_t), pointer,optional :: collective_user
        integer, intent(in) :: splinedegree
        scale_matrix_equation=1.0_f64

        if (present(collective_user)) then
            collective=>collective_user
            collective_rank = sll_get_collective_rank( collective )
            collective_size = sll_get_collective_size( collective )
        else
            nullify(collective)
        endif



        SLL_CLEAR_ALLOCATE(knots_mesh(size(knots)),ierr)
        knots_mesh=knots
        n_knots=size(knots)
        n_cells=n_knots-1
        spline_degree=splinedegree
        bspline_arbitrary_degree=>new_arbitrary_degree_spline_1d( spline_degree, &
            knots_mesh, size(knots), PERIODIC_ARBITRARY_DEG_SPLINE)


        selectcase(sll_pic1d_poisson_solver)
            case(SLL_SOLVER_FEM)
                !Do precalculations, especially for the FEM solver
                SLL_CLEAR_ALLOCATE(fem_inhomogenity_steady(1:n_cells),ierr)

                if (present(steadyparticleposition )) then
                    n_steadyparticles=size(steadyparticleposition)
                    fem_inhomogenity_steady=interpolate_particles_bsplines( bspline_arbitrary_degree, &
                        knots_mesh, steadyparticleposition, b_splines_at_x)

                    !Normalize background
                    fem_inhomogenity_steady=fem_inhomogenity_steady/n_steadyparticles
                    !Check if particle loading has to be in parallel
                    if (present(collective_user)) then
                        !we want the result to be on all nodes
                        fem_inhomogenity_steady=fem_inhomogenity_steady/collective_size
                        call sll_collective_globalsum_array_real64(collective,fem_inhomogenity_steady)
                        !                call sll_collective_globalsum_array_real64(collective,fem_inhomogenity_steady,0)
                    endif

                endif

                if (present(scalar_inhomogenity_user)) then
                    scalar_inhomogenity=scalar_inhomogenity_user
                else
                    scalar_inhomogenity = 1.0_f64
                endif


                SLL_ALLOCATE(fem_solution(n_cells),ierr)
                SLL_CLEAR_ALLOCATE(fem_inhomogenity(1:n_cells),ierr)

                !FFT--------------------------------------------------------------------------------------
                !To get the same output as in the MATLAB example use
                SLL_CLEAR_ALLOCATE(circulant_first_line_fourier(1:n_cells/2+1), ierr)

                forward_fftplan => fft_new_plan(n_cells,fem_solution,circulant_first_line_fourier, &
                    FFT_FORWARD+FFT_NORMALIZE)
                backward_fftplan=>fft_new_plan(n_cells,circulant_first_line_fourier,fem_solution, &
                    FFT_INVERSE  + FFT_NORMALIZE_INVERSE)
                SLL_DEALLOCATE_ARRAY(circulant_first_line_fourier, ierr)

                !------------------------------------------------------------------------------------------
                !Assemble stiffnes matrix
                call sll_bspline_fem_solver_1d_assemble_stiffnes_matrix
                !prepare fourier transformation for stiffness matrix
                call sll_bspline_fem_solver_1d_fft_precomputation(stiffn_matrix_first_line)

                !Assemble mass matrix for L2-Norm
                call sll_bspline_fem_solver_1d_assemble_mass_matrix

            case(SLL_SOLVER_FD)
                poissonsolverFD=>new(knots(1),knots(n_knots),n_knots,ierr)
                SLL_ALLOCATE(fem_solution(n_knots),ierr)

        endselect

    endsubroutine



    !<Does precomputation for the FFT solver
    subroutine sll_bspline_fem_solver_1d_fft_precomputation(circulant_matrix_first_line)
        sll_real64, dimension(:), intent(in)   :: circulant_matrix_first_line

        sll_real64, dimension(:) :: circulantvector(size(circulant_matrix_first_line))
        integer :: N
        sll_real64:: matrix_condition


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
        circulantvector=circulantvector*scale_matrix_equation



        if (allocated(circulant_first_line_fourier)) then
            SLL_DEALLOCATE_ARRAY(circulant_first_line_fourier, ierr)
        endif
        SLL_CLEAR_ALLOCATE(circulant_first_line_fourier(1:N/2+1), ierr)


        !Fourier transform the circulant seed
        call fft_apply_plan(forward_fftplan,circulantvector,circulant_first_line_fourier)


        !Check on Problems
        !        if (minval(abs(circulant_first_line_fourier))==0.0_f64) then
        !            print *, "WARNING: Calculating Eigenvalues of Circulant Matrix seed via fft: "
        !            print *, "Eigenvalues near machine precision, can't continue"
        !            print *, "Machine precision: ", epsilon(REALPART(circulant_first_line_fourier(1)))
        !            print *, "Min Eigenvalue: ", minval(REALPART(circulant_first_line_fourier)) ," + " , minval(IMAGPART(circulant_first_line_fourier)) ,"i"
        !            !Install preconditioner
        !            !print *, REALPART(circulant_first_line_fourier)
        !            !call sll_display(REALPART(circulant_first_line_fourier), "(F14.14)")
        !            !call sll_display(IMAGPART(circulant_first_line_fourier), "(F16.14)")
        !            !call sll_display(circulant_matrix_first_line, "(F12.8)")
        !            circulant_first_line_fourier(1)=1.0_f64
        !            !return
        !        endif

        !matrix_condition=maxval(abs(circulant_first_line_fourier))/minval(abs(circulant_first_line_fourier))
        !print *, "Stiffnes Matrix Condition: ", matrix_condition
        circulant_first_line_fourier(1)=1.0_f64

        matrix_condition=maxval(abs(circulant_first_line_fourier))/minval(abs(circulant_first_line_fourier))
        print *, "Corrected Stiffnes Matrix Condition: ", matrix_condition


    end subroutine

    subroutine sll_bspline_fem_solver_1d_assemble_mass_matrix
        sll_real64, dimension(spline_degree+1) :: mass_matrix_period

        SLL_ASSERT(n_cells==n_knots-1)
        SLL_CLEAR_ALLOCATE(mass_matrix_first_line(1:n_cells),ierr)

        mass_matrix_period=gen_fem_matrix_period (bspline_arbitrary_degree, knots_mesh, b_splines_at_x )
        !First line of stiffness matrix
        mass_matrix_first_line=0.0_f64
        mass_matrix_first_line(n_knots-1- (size(mass_matrix_period) -1) +1: n_knots-1)=&
            mass_matrix_period(size(mass_matrix_period):2:-1)
        mass_matrix_first_line(1:size(mass_matrix_period))=mass_matrix_period

    endsubroutine


    subroutine sll_bspline_fem_solver_1d_assemble_stiffnes_matrix
        sll_real64, dimension(spline_degree+1) :: stiffn_matrix_period
        sll_real64 :: fem_matrix_period_sum

        SLL_ASSERT(n_cells==n_knots-1)
        SLL_CLEAR_ALLOCATE(stiffn_matrix_first_line(1:n_cells),ierr)
        stiffn_matrix_period=gen_fem_matrix_period (bspline_arbitrary_degree, knots_mesh, b_spline_derivatives_at_x )


        !First line of stiffness matrix
        stiffn_matrix_first_line=0.0_f64
        stiffn_matrix_first_line(n_knots-1- (size(stiffn_matrix_period) -1) +1: n_knots-1)=&
            stiffn_matrix_period(size(stiffn_matrix_period):2:-1)
        stiffn_matrix_first_line(1:size(stiffn_matrix_period))=stiffn_matrix_period

        fem_matrix_period_sum=2*sum(stiffn_matrix_period(2:size(stiffn_matrix_period)))+ stiffn_matrix_period(1)

        if  ((fem_matrix_period_sum > FEM_MASS_MATRIX_TOL) &
                .OR. (fem_matrix_period_sum - sum(stiffn_matrix_first_line) > FEM_MASS_MATRIX_TOL)  ) then
            print *, "FEM matrix assembly failed: "
            print *, fem_matrix_period_sum    , sum(stiffn_matrix_first_line)
            call sll_display(stiffn_matrix_period, "(F12.8)" )
            stop
        endif

    endsubroutine


    !< Evaluates the first derivative of the solution at the given points knots_eval
    !< The Result is written into eval_solution
    subroutine sll_bspline_fem_solver_1d_eval_solution_derivative(knots_eval, eval_solution)
        implicit none
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution



            selectcase(sll_pic1d_poisson_solver)
            case(SLL_SOLVER_FEM)
                  eval_solution=bspline_basis_to_realvals(bspline_arbitrary_degree, knots_mesh,&
            fem_solution, knots_eval, b_spline_derivatives_at_x)
            case(SLL_SOLVER_FD)
               !Interpolate solution

eval_solution=0
               !eval_solution=sll_cloudincell_1d(knots_eval, (/1.0_f64/))
        endselect


    endsubroutine


    !< Evaluates the solution at the given points knots_eval
    !< The Result is written into eval_solution
    subroutine sll_bspline_fem_solver_1d_eval_solution(knots_eval, eval_solution)
        implicit none
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        SLL_ASSERT(size(knots_eval)==size(eval_solution))

        selectcase(sll_pic1d_poisson_solver)
            case(SLL_SOLVER_FEM)
                eval_solution=bspline_basis_to_realvals(bspline_arbitrary_degree, knots_mesh,&
                    fem_solution, knots_eval, b_splines_at_x)
            case(SLL_SOLVER_FD)
                print *, "EVALUATION OF SOLUTION NOT IMPLEMENTED"
                eval_solution=0
        endselect

    endsubroutine


    !<Solves the given FEM-Matrix, with right side as sum of deltafunctions at particleposition
    !<with sign (-1), so particleposition are the positions of the negative charges
    !<If none given solves the
    !< solves  - \laplace u = f See test case or more help
    !<
    subroutine sll_bspline_fem_solver_1d_solve( particleposition)
        sll_real64, DIMENSION(:), intent(in),optional :: particleposition
        !sll_real64, DIMENSION(:), intent(in),optional   :: particleweight_user
        !SLL_ALLOCATE(fem_inhomogenity(n_cells),ierr)



        selectcase(sll_pic1d_poisson_solver)
            case(SLL_SOLVER_FEM)
                if (associated(collective)) then
                    fem_inhomogenity= interpolate_particles_bsplines( bspline_arbitrary_degree, knots_mesh,&
                        particleposition, b_splines_at_x,weights )

                    !Here one would also have the option to solve on each core and
                    !omit the broadcast
                    call sll_collective_globalsum(collective, fem_inhomogenity, 0)

                    if (collective_rank==0) then
                        !fem_inhomogenity=scalar_inhomogenity *( fem_inhomogenity_steady-fem_inhomogenity)
                        !print *, "Inhom:",  sum(fem_inhomogenity)

                        !fem_inhomogenity=scalar_inhomogenity*(1.0_f64/n_cells-(fem_inhomogenity/fem_inhomogenity_steady))
                        fem_inhomogenity=scalar_inhomogenity*(fem_inhomogenity_steady-fem_inhomogenity)
                        !   print *,  sum(fem_inhomogenity)
                        !fem_inhomogenity=scalar_inhomogenity*(1.0_f64/n_cells-fem_inhomogenity&
                            !           *(knots_mesh(n_knots)-knots_mesh(1)) )
                        !Delta f, where fem_inhomogenity_steady is delta f
                        !                fem_inhomogenity=scalar_inhomogenity*(1.0_f64/n_cells - &
                            !                                (fem_inhomogenity-fem_inhomogenity_steady)*n_cells)
                        !                fem_inhomogenity=scalar_inhomogenity*(fem_inhomogenity*n_cells -1.0_f64/n_cells)

                        fem_solution=bspline_fem_solver_1d_solve_matrix_equation(fem_inhomogenity)
                    endif
                    ! call sll_collective_barrier(sll_world_collective)
                    !Broadcast result so it lies on every core for interpolation
                    call sll_collective_bcast_real64( collective, fem_solution, n_cells, 0 )
                else

                    if (present(particleposition)) then
                        !Calculate Electric Potential (gather)
                        !> \frac{q}{\epsilon_0}(n_i - n_e)
                        if (size(particleposition)==0) then
                            fem_inhomogenity=0
                        else
                            fem_inhomogenity=scalar_inhomogenity   &
                                *( fem_inhomogenity_steady - &
                                interpolate_particles_bsplines( bspline_arbitrary_degree, knots_mesh,&
                                particleposition, b_splines_at_x,weights ))
                        endif
                    else
                        fem_inhomogenity=scalar_inhomogenity*fem_inhomogenity_steady
                    endif
                    !Solve
                    fem_solution=bspline_fem_solver_1d_solve_matrix_equation(fem_inhomogenity)

                endif
            case(SLL_SOLVER_FD)
                if (present(particleposition)) then
                    fem_inhomogenity=sll_cloudincell_1d(particleposition,weights)
                    if (associated(collective)) call sll_collective_globalsum(collective, fem_inhomogenity, 0)
                endif

                if ( ( .NOT. associated(collective) ) .OR. collective_rank==0) then
                    fem_inhomogenity=scalar_inhomogenity*(fem_inhomogenity_steady-fem_inhomogenity)

                    call solve(poissonsolverFD,fem_solution, fem_inhomogenity)
                endif

                if (associated(collective)) call sll_collective_bcast_real64( collective, fem_solution, n_knots, 0 )
        endselect



        !fem_solution=fem_solution-sum(fem_solution)/n_cells
        !print *, sum(fem_solution), sum (fem_inhomogenity)
        !SLL_DEALLOCATE_ARRAY(fem_inhomogenity, ierr)
    endsubroutine



    !<Linear interpolator for particles, which corresponds to a cloud in cell scheme in 1d
    function sll_cloudincell_1d(ppos, pweight) result(rho)
        sll_real64, dimension(:), intent(in) :: ppos
        sll_real64, dimension(:), intent(in) :: pweight
        sll_real64, dimension(n_knots)  :: rho
        sll_int32, dimension(size(ppos)) :: cell
        !!SLL_ASSERT(size(pweight)==size(ppos))
        rho=0
        cell=bspline_fem_solver_1d_cell_number(knots_mesh,ppos)
        rho(cell)=(ppos-knots_mesh(cell+1)  )/(knots_mesh(cell+1)-knots_mesh(cell))*pweight
        rho(cell+1)=(knots_mesh(cell)-ppos)/(knots_mesh(cell+1)-knots_mesh(cell))*pweight


        !periodic
        !rho=rho(1:end-1)
    endfunction

    !<Only solve the fem equation on one core
    subroutine sll_bspline_fem_solver_1d_solve_collective( particleposition_collective,weights_collective)
        sll_real64, DIMENSION(:), intent(in) :: particleposition_collective
        sll_real64, DIMENSION(:), intent(in) :: weights_collective
        !sll_real64, dimension(:), allocatable :: recieveinhomogenity
        !sll_int32 :: collective_rank, collective_size
        !Stuff to have on every core
        !n_cells
        !bspline_arbitrary_degree
        !
        !allocate fem_inhomogenity

        fem_inhomogenity=interpolate_particles_bsplines( bspline_arbitrary_degree, knots_mesh,&
            particleposition_collective, b_splines_at_x,weights_collective )

        !Calculate the sum of fem_inhomogenity over the collective and save it in node 0
        !#######################################################################
        !        if (collective_rank ==0) then
        !                SLL_ALLOCATE(recieveinhomogenity(1:n_cells),ierr)
        !        endif
        !        call sll_collective_reduce_real64( collective, fem_inhomogenity, n_cells, MPI_SUM, 0, &
            !        recieveinhomogenity )

        call sll_collective_globalsum(collective, fem_inhomogenity, 0)
        if (collective_rank ==0) then
            fem_inhomogenity=scalar_inhomogenity*(fem_inhomogenity_steady- fem_inhomogenity )
            !SLL_DEALLOCATE_ARRAY(recieveinhomogenity,ierr)
        endif
        !#######################################################################

        !Get the actual FEM solution only on the local node
        if (collective_rank ==0) then
            fem_solution=bspline_fem_solver_1d_solve_matrix_equation(fem_inhomogenity)
        endif


    endsubroutine




    function bspline_particle_density(particleposition) result(density)
        sll_real64, dimension(n_cells) :: density
        sll_real64, DIMENSION(:), intent(in):: particleposition

        density=interpolate_particles_bsplines( bspline_arbitrary_degree, knots_mesh,&
            particleposition, b_splines_at_x,weights )

    endfunction

    function test_poisson_solver_testfunction(x) result(y)
        sll_real64, dimension(:), intent(in) :: x
        !sll_real:: test_mode
        sll_real64, dimension(size(x)) :: y
        !test_mode=4
        y= ((test_mode*sll_kx)**2 )*sin(test_mode*sll_kx*x)
    endfunction



    subroutine set_test_mode( test_mode_user)
        sll_real64, intent(in):: test_mode_user
        test_mode=test_mode_user
    endsubroutine

    function test_poisson_solver_testfunction2(x) result(y)
        sll_real, dimension(:), intent(in) :: x
        sll_real:: test_mode
        sll_real, dimension(size(x)) :: y
        y= (x+1)*(x-0.5)*(x-1)
    endfunction

    !<Simple Wrapper
    !<Calculates the inhomogenity b={<f, \phi_i>, i =1, ... N} with given function f
    !<by Gauss Legendre integration
    subroutine sll_bspline_fem_solver_1d_calc_inhomogenity_function(eval_function, n_quadrature_points_user)
        integer, intent(in), optional      :: n_quadrature_points_user
        procedure (test_poisson_solver_testfunction) :: eval_function
        fem_inhomogenity_steady=0
        call sll_bspline_fem_solver_1d_add_inhomogenity_function(eval_function, n_quadrature_points_user)
    endsubroutine

    !<Calculates the inhomogenity b={<f, \phi_i>, i =1, ... N} with given function f
    !<by Gauss Legendre integration
    subroutine sll_bspline_fem_solver_1d_add_inhomogenity_function(eval_function, n_quadrature_points_user)
        implicit none
        procedure (test_poisson_solver_testfunction) :: eval_function
        !procedure (function1d) :: eval_function
        sll_real64, dimension(:,:), allocatable   :: quadrature_point_vals
        sll_real64, dimension(:,:), allocatable   :: quadrature_points_weights
        sll_real64, dimension(:,:), allocatable   :: bspline_qpoint_values
        integer, intent(in), optional      :: n_quadrature_points_user
        sll_real64, dimension(:), allocatable    :: fem_matrix_period
        integer                             :: n_quadrature_points
        integer :: cell
        integer :: i,j

        !Get Gauss Legendre points and weights to be exact for the selected spline degree
        !Note a Bspline is a piecewise polynom
        if ( .not. present(n_quadrature_points_user)) then
            !Gauss Legendre, 2m-1
            n_quadrature_points=ceiling(0.5_f64*real(2*spline_degree+1))
            !            !Gauss Lobatto, 2m-3
            !            n_quadrature_points=ceiling(0.5_f64*real(2*spline_degree+3))
        else
            n_quadrature_points=n_quadrature_points_user
        endif


        SLL_ALLOCATE(quadrature_points_weights(2,1:n_quadrature_points),ierr)
        SLL_CLEAR_ALLOCATE(bspline_qpoint_values(1:spline_degree+1,1:n_quadrature_points), ierr)

        SLL_ASSERT(size(knots_mesh) >spline_degree+1)

        !Quadrature for each support cell
        do cell=1, n_cells
            quadrature_points_weights=gauss_legendre_points_and_weights(n_quadrature_points, knots_mesh(cell), knots_mesh(cell+1))
            !quadrature_points_weights=gauss_lobatto_points_and_weights(n_quadrature_points, knots_mesh(cell), knots_mesh(cell+1))
            do i=1,n_quadrature_points
                bspline_qpoint_values(:,i)=b_splines_at_x(bspline_arbitrary_degree, cell, quadrature_points_weights(1,i ))
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
            do j=0,spline_degree
                !i=cell - j
                i=cell - spline_degree + j
                if (i > n_cells) then
                    i= i - n_cells
                elseif (i < 1) then
                    i = i + n_cells
                endif
                fem_inhomogenity_steady(i)= fem_inhomogenity_steady(i) + &
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


    endsubroutine


    !<Mostly used to set the heavy ion contribution to a constant level
    subroutine sll_bspline_fem_solver_1d_set_inhomogenity_constant( average )
        sll_real64, intent(in) :: average
        fem_inhomogenity_steady=average
    endsubroutine
    !
    !
    !
    !    function sll_bspline_fem_solve_poisson_1d(bspline_arbitrary_degree, knots,b)  result(solution)
    !        type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline_arbitrary_degree
    !        sll_real64, dimension(:), intent(in)  :: knots
    !        sll_real64, dimension(:), intent(in)   :: b
    !        sll_real64, dimension(:) , allocatable:: stiffn_matrix_period
    !        sll_real64 ::stiffn_matrix_first_line(size(knots)-1)
    !        sll_real64  :: solution(size(b))
    !        sll_real64 :: fem_matrix_period_sum
    !        stiffn_matrix_period=gen_stiffn_matrix_period (bspline_arbitrary_degree, knots)
    !
    !        !First line of stiffness matrix
    !        stiffn_matrix_first_line=0.0_f64
    !        stiffn_matrix_first_line(size(knots)-1- (size(stiffn_matrix_period) -1) +1: size(knots)-1)=&
        !            stiffn_matrix_period(size(stiffn_matrix_period):2:-1)
    !        stiffn_matrix_first_line(1:size(stiffn_matrix_period))=stiffn_matrix_period
    !
    !    fem_matrix_period_sum=2*sum(stiffn_matrix_period(2:size(stiffn_matrix_period)))+ stiffn_matrix_period(1)
    !
    !    if  ((fem_matrix_period_sum > FEM_MASS_MATRIX_TOL) &
        !    .OR. (fem_matrix_period_sum - sum(stiffn_matrix_first_line) > FEM_MASS_MATRIX_TOL)  ) then
    !           print *, "FEM matrix assembly failed: "
    !           print *,fem_matrix_period_sum    , sum(stiffn_matrix_first_line)
    !           call sll_display(stiffn_matrix_period, "(F12.8)" )
    !        stop
    !    endif
    !
    !
    !    SLL_DEALLOCATE_ARRAY(stiffn_matrix_period,ierr)
    !        !call sll_display(stiffn_matrix_circular_seed,  "(F8.4)" )
    !
    !        !call sll_display(stiffn_matrix_first_line,  "(F8.4)" )
    !        solution=bspline_fem_solver_1d_solve_circulant_matrix_equation( stiffn_matrix_first_line, b)
    !    endfunction

    !Interpolate right side on splines
    !In parallel the following steps will be performed:
    !Each core is responsible for a number of cells
    !OR
    !Each core is responsible for his particles
    !interpolfun is optional and b_splines_at_x by default
    !This is the delta function
    function interpolate_particles_bsplines( bspline, knots, particleposition_in, interpolfun_user, particleweight) &
            result(b)
        sll_real64, DIMENSION(:), intent(in):: particleposition_in
        sll_real64, DIMENSION(size(particleposition_in)) :: particleposition
        type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline
        sll_real64, dimension(:), intent(in)     :: knots
        sll_real64, dimension(:), intent(in), optional :: particleweight
        !sll_int, dimension(:), intent(in)     :: cells
        sll_real64, dimension(size(knots)-1)  :: b
        procedure (b_splines_at_x), optional:: interpolfun_user
        procedure (b_splines_at_x), pointer :: interpolfun
        integer :: idx
        sll_int :: cell
        sll_real64, dimension(bspline%degree+1) :: b_contribution
        sll_int :: b_idx
        sll_int :: b_contrib_idx


        if (present(interpolfun_user)) then
            interpolfun=>interpolfun_user
        else
            interpolfun=>b_splines_at_x
        endif

        b(:)=0

        !Check boundary Conditions
        particleposition=ensure_boundaryc(particleposition_in)

        do idx=1,size(particleposition)

            cell= bspline_fem_solver_1d_cell_number(knots, particleposition(idx))

            b_contribution=interpolfun(bspline,cell, particleposition(idx))
            !b_contribution=interpolfun(bspline,5, particleposition(idx) - (cell-1)*(knots(2)-knots(1)))


            if (present(particleweight))  then
                b_contribution=b_contribution*particleweight(idx)
            endif


            !Amount of splines equals the amount of
            !b_contribution(1:bspline%degree+1)=b_period
            !b_contribution(bspline%degree+2: bspline%degree*2 +1)=b_period(size(b_period)-1:1:-1)
            !First spline has first support in first cell of interval
            do b_contrib_idx=1, bspline%degree+1
                !b_idx=cell + (b_contrib_idx- (bspline%degree+1))
                b_idx=cell - b_contrib_idx+1
                if (b_idx > n_cells) then
                    b_idx = b_idx - n_cells
                elseif (b_idx < 1) then
                    b_idx = b_idx + n_cells
                endif

                SLL_ASSERT(b_idx <= n_cells)
                SLL_ASSERT(b_idx >= 1 )
                !                       print *, b_idx, particleposition(idx) , cell, b_contribution(b_contrib_idx)
                b(b_idx)=b(b_idx)+b_contribution(b_contrib_idx)
            enddo
            !b(cell)=b(cell)+1
            !
            !            do b_contrib_idx=1, bspline%degree+1
            !                b_idx=cell + (b_contrib_idx) !- (bspline%degree+1))
            !
            !                b_idx=mod(b_idx -1 + 2*size(b), size(b))+1
            !                       print *, b_idx, particleposition(idx) , cell, b_contribution(b_contrib_idx)
            !                b(b_idx)=b(b_idx)+b_contribution(b_contrib_idx)
            !            enddo
            !MOD(cell:
            !b(cell:cell+ bspline_arbitrary_degree%degree+1:1)
        enddo

    endfunction

    !< knots are the interpolant points and have nothing to do with the mesh
    function bspline_basis_to_realvals (  bspline, knots_mesh, bspline_vector, knots_eval_in, interpolfun_user ) &
            result(realvals)
        type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline
        procedure (b_splines_at_x),pointer :: interpolfun
        procedure (b_splines_at_x), optional :: interpolfun_user
        sll_real64, dimension(:), intent(in)     :: knots_mesh
        sll_real64, dimension(:), intent(in)     :: bspline_vector
        sll_real64, dimension(:), intent(in)     :: knots_eval_in
        sll_real64, dimension(size(knots_eval_in))    :: knots_eval
        sll_real64 ::  realvals(size(knots_eval_in))
        sll_real64 :: b_contribution(bspline%degree+1)
        sll_int32   :: b_idx
        sll_int64 :: eval_idx
        sll_int :: b_contrib_idx
        sll_int :: cell

        !Check Boundary Conditions for evaluation
        knots_eval=ensure_boundaryc(knots_eval_in)

        if (present(interpolfun_user)) then
            interpolfun=>interpolfun_user
        else
            print *, "Warning no Interpolation function given!"
            interpolfun=>b_splines_at_x
        endif
        SLL_ASSERT( size(knots_eval)==size(realvals))
        SLL_ASSERT( size(knots_mesh)==size(bspline_vector)+1)

        realvals=0



        !Loop over all points to evaluate
        do eval_idx=1,size(knots_eval)
            cell= bspline_fem_solver_1d_cell_number(knots_mesh, knots_eval(eval_idx))
            !Get the values for the spline at the eval point
            b_contribution=interpolfun(bspline,cell, knots_eval(eval_idx))

            realvals(eval_idx)=0
            do b_contrib_idx=1,bspline%degree+1
                !Determine which value belongs to which spline
                !b_idx=cell + (b_contrib_idx- (bspline%degree+1))
                b_idx=cell - (bspline%degree+1)  +(b_contrib_idx)
                !Periodicity
                if (b_idx > n_cells) then
                    b_idx = b_idx - n_cells
                elseif (b_idx < 1) then
                    b_idx = b_idx + n_cells
                endif

                realvals(eval_idx)=realvals(eval_idx)  + bspline_vector(b_idx)*b_contribution(b_contrib_idx)

                !                !forall (b_contrib_idx=1:bspline%degree+1)
                !                b_idx(b_contrib_idx)=cell  + (b_contrib_idx)!- (bspline%degree+1))
                !                !periodicity
                !                b_idx(b_contrib_idx)=1 + mod(b_idx(b_contrib_idx) -1 + 2*size( bspline_vector ), &
                    !                    size( bspline_vector))
                !                    print *, eval_idx, b_contrib_idx,  knots_eval(eval_idx), b_idx(b_contrib_idx)
                !                realvals(eval_idx)=realvals(eval_idx) &
                    !                    +b_contribution(b_contrib_idx)*bspline_vector( b_idx(b_contrib_idx)  )
                !endforall
            enddo

            !  do eval_idx=1,size(knots_eval)
            !            cell= bspline_fem_solver_1d_cell_number(knots_mesh, knots_eval(eval_idx))
            !            b_contribution=interpolfun(bspline,cell, knots_eval(eval_idx))
            !
            !            do b_contrib_idx=1,bspline%degree+1
            !
            !
            !
            !                !forall (b_contrib_idx=1:bspline%degree+1)
            !                b_idx(b_contrib_idx)=cell  + (b_contrib_idx)!- (bspline%degree+1))
            !                !periodicity
            !                b_idx(b_contrib_idx)=1 + mod(b_idx(b_contrib_idx) -1 + 2*size( bspline_vector ), &
                !                    size( bspline_vector))
            !                    print *, eval_idx, b_contrib_idx,  knots_eval(eval_idx), b_idx(b_contrib_idx)
            !                realvals(eval_idx)=realvals(eval_idx) &
                !                    +b_contribution(b_contrib_idx)*bspline_vector( b_idx(b_contrib_idx)  )
            !                !endforall
            !            enddo

        enddo


    endfunction


    !<First Cell is in analogy to fortran arrays Cell Number 1
    function  bspline_fem_solver_1d_cell_number_array(knots_mesh, particle_position) result(cell)
        sll_real64, dimension(:), intent(in)     :: knots_mesh
        sll_real64, dimension(:),  intent(in) ::particle_position
        sll_int32, dimension(size(particle_position)) :: cell

        cell=floor((particle_position - knots_mesh(1))/(knots_mesh(2)-knots_mesh(1)))+1
    endfunction


    function  bspline_fem_solver_1d_cell_number_single(knots_mesh, particle_position) result(cell)
        sll_real64, dimension(:), intent(in)     :: knots_mesh
        sll_real64, intent(in) ::particle_position
        sll_int32 :: cell
        SLL_ASSERT(size(knots_mesh)>=3)
        SLL_ASSERT(particle_position>=knots_mesh(1))
        SLL_ASSERT(particle_position<=knots_mesh(size(knots_mesh)))

        cell=floor((particle_position - knots_mesh(1))/(knots_mesh(2)-knots_mesh(1)))+1




        !Periodicity Exception: last knot belongs to last cell
        if (cell==size(knots_mesh))  then
            cell=cell-1
            !cell=1
            print *, "Particle sitting on Last knot of Mesh", particle_position
            !stop
        else
            do while (particle_position>=knots_mesh(cell+1))
                cell=cell +1
            enddo
            do while (particle_position<knots_mesh(cell))
                cell=cell -1
            enddo
            SLL_ASSERT( particle_position < knots_mesh(cell+1))
        endif


        SLL_ASSERT(cell>=1)
        SLL_ASSERT(cell<=n_cells)
        SLL_ASSERT( particle_position >= knots_mesh(cell))

    endfunction




    !function interpolate_particles_bsplines( particleposition, knots , spline_degree ) &
        !                                                        result(fem_vector)
    !sll_real, DIMENSION(:), allocatable, intent(in):: particleposition
    !sll_real64, dimension(:), intent(in)     :: knots
    !sll_int, intent(in) :: spline_degree
    !sll_real64, dimension(1:size(knots)) :: fem_vector
    !class(arb_deg_1d_interpolator),pointer :: interpolator

    !type(arb_deg_1d_interpolator) :: interpolator

    !!SLL_ASSERT(knots(1)<=min(particleposition,1)   )
    !!SLL_ASSERT(knots(size(knots))>=max(particleposition,1)   )


    !call initialize_ad1d_interpolator(interpolator,&
        !       size(knots), &
        !       knots(1), &
        !       knots(size(knots)), &
        !       SLL_PERIODIC, &
        !       SLL_PERIODIC, &
        !       spline_degree)
    !interpolator&


    !delete()
    !endfunction


    !<can be optimized, with a sparse first line
    function bspline_fem_solver_1d_circulant_matrix_vector_product( circulant_matrix_first_line, vector ) &
            result(solution)
        sll_real64, dimension(:), intent(in)   :: circulant_matrix_first_line
        sll_real64, dimension(:), intent(in)   :: vector
        sll_real64, dimension(size(vector)) ::solution

        integer idx, N
        N=size(circulant_matrix_first_line)
        SLL_ASSERT(size(circulant_matrix_first_line)==size(vector))
        solution=0
        do idx=1, N
            solution(idx)=dot_product(cshift(circulant_matrix_first_line, -(idx-1)), vector)
            !print *,cshift(circulant_matrix_first_line, -(idx-1)),"#"

        enddo
    endfunction

    !<Gives the squared H1-seminorm of the solution $\Phi$: $|\nabla \Phi|^2$
    function bspline_fem_solver_1d_H1seminorm_solution() result(seminorm)
        sll_real64 :: seminorm
        integer :: idx, N
        sll_real64 , dimension(:) :: solution(size(fem_solution))
        sll_real64, dimension(n_cells) :: matrix_product

        !Since the input data is real, the data in fourier space satisfies
        !X_{N-k}=X_k^* (complex conjugate) which can be stored more efficiently,
        !so for the complex data there will be only allocated space for N/2 +1 values
        sll_comp64 , dimension(:) :: data_complex(size(fem_solution)/2+1)

        !!type(sll_fft_plan), pointer :: forward_fftplan => null()
        !!type(sll_fft_plan), pointer :: backward_fftplan => null()
        solution=fem_solution

        !Determine dimension of problem
        N=size(fem_solution)
        SLL_ASSERT(N==n_cells)
        SLL_ASSERT(is_power_of_two(int(N,i64)))



        call fft_apply_plan(forward_fftplan,solution,data_complex)

        data_complex=data_complex*(circulant_first_line_fourier)
        data_complex(1)=0.0_f64

        call fft_apply_plan(backward_fftplan,data_complex,solution)


        !!!
        !Somehow the normalization does not do what it should do:
        solution=solution/(N)

        seminorm=dot_product(fem_solution, solution )
        !seminorm=sqrt(dot_product(fem_solution, fem_solution ))
        !        matrix_product=bspline_fem_solver_1d_circulant_matrix_vector_product&
            !            (stiffn_matrix_first_line,fem_solution)
        !        seminorm=0
        !        seminorm=dot_product(fem_solution, matrix_product )

    endfunction



    function bspline_fem_solver_1d_L2norm_solution() result(l2norm)
        sll_real64 :: l2norm
        integer :: idx
        sll_real64, dimension(n_cells) :: matrix_product
        matrix_product=bspline_fem_solver_1d_circulant_matrix_vector_product&
            (mass_matrix_first_line,fem_solution)
        l2norm=dot_product(fem_solution, matrix_product )
        !seminorm=sqrt(sum(fem_solution**2))
    endfunction

    !    Circular seed to generate stiffness matrix as a circular matrix
    !    Note this only works because the period is symmetric
    !    <Solve Cx=b
    function bspline_fem_solver_1d_solve_circulant_matrix_equation( circulant_matrix_first_line, &
            rightside  )  result(solution)

        sll_real64, dimension(:), intent(in)   :: circulant_matrix_first_line
        sll_real64, dimension(:), intent(in)   :: rightside
        sll_real64, dimension(:) :: constant_factor(size(rightside))
        sll_real64 , dimension(:) :: solution(size(rightside))
        integer :: N
        type(sll_fft_plan), pointer :: forward_fftplan => null()
        type(sll_fft_plan), pointer :: backward_fftplan => null()
        !Since the input data is real, the data in fourier space satisfies
        !X_{N-k}=X_k^* (complex conjugate) which can be stored more efficiently,
        !so for the complex data there will be only allocated space for N/2 +1 values
        !sll_real64 , dimension(:),allocatable :: circulantv
        sll_comp64 , dimension(:),allocatable ::circulantv_fourier

        sll_comp64 , dimension(:),allocatable ::constant_factor_fourier
        sll_comp64 , dimension(:),allocatable :: data_complex

        sll_real64, dimension(:) :: circulantv(size(circulant_matrix_first_line))
        !Determine dimension of problem
        N=size(rightside)
        SLL_ASSERT(is_power_of_two(int(N,i64)))
        SLL_ASSERT(size(circulant_matrix_first_line)==N)
        stop
        !        !Generate Circulant Seed
        !        !!SLL_CLEAR_ALLOCATE(circulantv(N),ierr)
        !        !!SLL_ASSERT(size(c_seed)<=N)
        !        !circulantv(1:size(c_seed))=c_seed
        !        !Circulant seed c to generate stiffness matrix as a circular matrix
        !        !Remember c is not the first line of the circulant matrix
        !        !this is the first line of the matrix
        !        ! (c_1  0 0  0  .... c_4   c_3  c_2)
        !        !Note this only works because the period is symmetric
        !        circulantv=0
        !        circulantv(1)=circulant_matrix_first_line(1)
        !        circulantv(2:N)=circulant_matrix_first_line(N:2:-1)
        !
        !        SLL_CLEAR_ALLOCATE(circulantv_fourier(1:N/2+1), ierr)
        !        !!SLL_ALLOCATE(circulantv_fourier(N), ierr)
        !        SLL_ALLOCATE(constant_factor_fourier(N/2+1), ierr)
        !        SLL_ALLOCATE(data_complex(N/2+1),ierr)
        !
        !
        !
        !        constant_factor=rightside
        !
        !        !forward_fftplan => fft_new_plan(N,circulantv,circulantv_fourier,FFT_FORWARD)
        !        !To get the same output as in the MATLAB example use
        !        forward_fftplan => fft_new_plan(N,circulantv,circulantv_fourier,FFT_FORWARD + FFT_NORMALIZE)
        !        !forward_fftplan => fft_new_plan(N,circulantv,circulantv_fourier,FFT_FORWARD)
        !        backward_fftplan =>fft_new_plan(N,data_complex,solution, FFT_INVERSE+ FFT_NORMALIZE_INVERSE )
        !
        !        !Go into Fourier space and do division there, circulantv will be changed after this
        !        call fft_apply_plan(forward_fftplan,circulantv,circulantv_fourier)
        !
        !
        !        !call fft_apply_plan(forward_fftplan,constant_factor,constant_factor_fourier)
        !
        !        if (minval(abs(circulantv_fourier))==0.0_f64) then
        !            print *, "WARNING: Calculating Eigenvalues of Circulant Matrix seed via fft: "
        !            print *, "Eigenvalues near machine precision, can't continue"
        !            print *, "Machine precision: ", epsilon(REALPART(circulantv_fourier(1)))
        !            print *, "Min Eigenvalue: ", minval(REALPART(circulantv_fourier)) ," + " , minval(IMAGPART(circulantv_fourier)) ,"i"
        !            !Install preconditioner
        !            !print *, REALPART(circulantv_fourier)
        !            !call sll_display(REALPART(circulantv_fourier), "(F14.14)")
        !            !call sll_display(IMAGPART(circulantv_fourier), "(F16.14)")
        !            !call sll_display(stiffn_matrix_first_line, "(F12.8)")
        !            stop
        !        endif
        !
        !
        !        data_complex=constant_factor_fourier/circulantv_fourier
        !
        !        SLL_DEALLOCATE_ARRAY(constant_factor_fourier,ierr)
        !        SLL_DEALLOCATE_ARRAY(circulantv_fourier,ierr)
        !
        !        call fft_apply_plan(backward_fftplan,data_complex,solution)
        !        SLL_DEALLOCATE_ARRAY(data_complex,ierr)
        !
        !        call fft_delete_plan(forward_fftplan)
        !        call fft_delete_plan(backward_fftplan)

    endfunction


    function bspline_fem_solver_1d_solve_matrix_equation( rightside  )  result(solution)
        sll_real64, dimension(:), intent(in)   :: rightside
        sll_real64, dimension(:) :: constant_factor(size(rightside))
        sll_real64 , dimension(:) :: solution(size(rightside))
        integer :: N

        !        type(sll_fft_plan), pointer :: forward_fftplan => null()
        !        type(sll_fft_plan), pointer :: backward_fftplan => null()

        !Since the input data is real, the data in fourier space satisfies
        !X_{N-k}=X_k^* (complex conjugate) which can be stored more efficiently,
        !so for the complex data there will be only allocated space for N/2 +1 values
        sll_comp64 , dimension(:),allocatable ::constant_factor_fourier
        sll_comp64 , dimension(:),allocatable :: data_complex
        !Determine dimension of problem
        N=size(rightside)
        SLL_ASSERT(N==n_cells)
        SLL_ASSERT(is_power_of_two(int(N,i64)))

        SLL_CLEAR_ALLOCATE(constant_factor_fourier(1:N/2+1), ierr)
        SLL_CLEAR_ALLOCATE(data_complex(1:N/2+1),ierr)
        constant_factor=0
        constant_factor=rightside *scale_matrix_equation

        !        forward_fftplan=>fft_new_plan(N,constant_factor,constant_factor_fourier,FFT_FORWARD+ FFT_NORMALIZE)
        !        backward_fftplan=>fft_new_plan(N,data_complex,solution, FFT_INVERSE  + FFT_NORMALIZE_INVERSE)

        call fft_apply_plan(forward_fftplan,constant_factor,constant_factor_fourier)

        !constant_factor_fourier(1)=0
        data_complex=constant_factor_fourier/(circulant_first_line_fourier)
        data_complex(1)=0
        SLL_DEALLOCATE_ARRAY(constant_factor_fourier,ierr)

        call fft_apply_plan(backward_fftplan,data_complex,solution)
        SLL_DEALLOCATE_ARRAY(data_complex,ierr)

        !        call fft_delete_plan(backward_fftplan)
        !        call fft_delete_plan(forward_fftplan)
        !!!
        !Somehow the normalization does not do what it should do:
        solution=solution/N
    endfunction




    function gen_mass_matrix_period (bspline_arbitrary_degree, knots ) &
            result( mass_matrix_period)
        type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline_arbitrary_degree
        sll_real64, dimension(:), intent(in)     :: knots
        sll_real64, dimension(:), allocatable    :: mass_matrix_period

        mass_matrix_period=gen_fem_matrix_period (bspline_arbitrary_degree, knots, b_splines_at_x )
    endfunction

    function gen_stiffn_matrix_period (bspline_arbitrary_degree, knots) &
            result( stiffn_matrix_period)
        type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline_arbitrary_degree
        sll_real64, dimension(:), intent(in)     :: knots
        sll_real64, dimension(:), allocatable    :: stiffn_matrix_period

        stiffn_matrix_period=gen_fem_matrix_period (bspline_arbitrary_degree, knots, b_spline_derivatives_at_x )
    endfunction


    !<Generates the FEM matrix for arbitrary spline evaluation functions
    function gen_fem_matrix_period (bspline_arbitrary_degree, knots, spline_evaluation, quadrature_points ) &
            result( fem_matrix_period)
        type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline_arbitrary_degree
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
                !print *, "Quadrature point values: ", i
                !call sll_display(bspline_qpoint_values(:,i), "(F8.4)")
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


        !call sll_display(knots(1:degree+2), "(F10.5)" )
        !call sll_display(quadrature_points_weights(1,:), "(F10.5)" )
        !call sll_display(quadrature_points_weights(2,:), "(F10.5)" )

        !call sll_display(bspline_qpoint_values, "(F10.5)" )
        !write(*,*)

        !!<SLL_DEALLOCATE_ARRAY( bspline_qpoint_values,ierr)
        SLL_DEALLOCATE_ARRAY( quadrature_points_weights,ierr)
    endfunction

    function ensure_boundaryc(x)  result(xout)
        sll_real64, dimension(:), intent(in) ::x
        sll_real64, dimension(size(x))  :: xout
        sll_real64 :: interval_a, interval_b,interval_length
        interval_a=knots_mesh(1)
        interval_b=knots_mesh(n_knots)
        interval_length=interval_b-interval_a
        SLL_ASSERT(interval_a < interval_b)

        selectcase(sll_bspline_fem_solver_boundary_type)
            case(SLL_PERIODIC)
                xout=x-interval_a

                do while (minval(xout)<0.0_f64)
                    xout=xout +interval_length
                enddo
                xout=mod(xout, interval_length)&
                    +interval_a

                SLL_ASSERT(minval(xout)>=interval_a)
                SLL_ASSERT(maxval(xout)<interval_b)
        endselect

    endfunction



end module
