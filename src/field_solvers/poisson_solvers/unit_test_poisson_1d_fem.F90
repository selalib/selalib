#ifndef SLL_CIRCULAR_MATRIX_SOLVER_TOL
#define SLL_CIRCULAR_MATRIX_SOLVER_TOL 1e-10
#endif
  module unit_test_poisson_1d_fem_module
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_constants
    sll_real64 :: testfunction_test_mode

    contains
    function sll_poisson_1d_fem_testfunction(x) result(y)
        sll_real64, dimension(:), intent(in) :: x
        !sll_real:: test_mode
        sll_real64, dimension(size(x)) :: y
        !test_mode=4
        y= ((testfunction_test_mode*sll_kx)**2 )*sin(testfunction_test_mode*sll_kx*x)
    endfunction

endmodule unit_test_poisson_1d_fem_module

program unit_test_poisson_1d_fem
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"
    use unit_test_poisson_1d_fem_module
    use sll_poisson_1d_fem
    use sll_logical_meshes
    use sll_constants
    implicit none
    integer :: ierr
    sll_int :: idx, jdx


    !call test_poisson_solver(10, 3)
    !stop
    do jdx=1,6
        do idx=jdx+1, 12
            print *, "--------------------------------------------------------------------------------------"

            print *, "Starting test at Spline degree ", jdx ," with ", 2**idx, " Cells"
            call test_poisson_solver(idx, jdx)
        enddo

    enddo
    stop
!    do idx=1, 2
!        call test_circulant_matrix_eq_solver(ierr, 18)
!    enddo

contains






!    subroutine test_poisson_solver_circulant_matrix_vector_product(test_power_two)
!        implicit none
!        integer :: ierr, idx
!        sll_int ,intent(in) :: test_power_two
!        sll_real64, dimension(2**test_power_two) ::  vector
!        sll_real64, dimension(2**test_power_two) :: circulant_matrix_first_line
!        sll_real64, dimension(2**test_power_two) :: solution
!
!
!        do idx=1,100
!            call  random_number(vector)
!            call  random_number(circulant_matrix_first_line)
!
!            vector=vector-0.5_f64
!            circulant_matrix_first_line=circulant_matrix_first_line
!
!            solution=bspline_fem_solver_1d_circulant_matrix_vector_product( circulant_matrix_first_line, vector )
!
!            if (dot_product(vector, solution)<0.0_f64 ) then
!                print *, "Error, Matrix not positiv definit", idx
!                stop
!            endif
!
!        enddo
!
!    endsubroutine


    subroutine test_poisson_solver(test_power_two, spline_degree)
        implicit none
        integer :: ierr
        sll_real64 :: interval_a=0, interval_b=1.0_f64
        type(sll_logical_mesh_1d), pointer :: mesh=>null()    !Finite Element mesh
        type(sll_logical_mesh_1d), pointer :: mesh_eval=>null()  !evaluation mesh
        class(poisson_1d_fem), pointer :: solver=>null()
        sll_int32, intent(in) :: test_power_two
        sll_int32, intent(in) :: spline_degree
        sll_real64, dimension(test_power_two) :: numerical_error
        sll_real64, dimension(test_power_two) :: numerical_error_deriv
        sll_real64, dimension(test_power_two) :: mode_error
        sll_real64 :: num_error_sum
        sll_int32 :: idx, jdx,node
        sll_int32 :: test_dimension

        sll_real64, dimension(:), allocatable :: equ_right_side
        sll_real64, dimension(:), allocatable :: solution

        !sll_real64 :: mesh_delta_t
        sll_real64, dimension(2**test_power_two+1 ) :: knots
        sll_real64, dimension(2)  :: knots_eval
        sll_real64, dimension(:), allocatable :: actual_solution
        sll_real64 :: H1seminorm, residual
        integer :: test_mode
        !Test Mass Matrix Assembler
        test_dimension=2**test_power_two


        mesh=>new_logical_mesh_1d( test_dimension, interval_a, interval_b )
        knots=sll_mesh_nodes(mesh)


        mesh_eval=>new_logical_mesh_1d( size(knots_eval)-1, interval_a, interval_b )
        knots_eval=sll_mesh_nodes(mesh_eval)

        solver=>new_poisson_1d_fem(mesh, spline_degree, SLL_PERIODIC, ierr)


        SLL_CLEAR_ALLOCATE(solution(1:size(knots_eval)),ierr)
        SLL_CLEAR_ALLOCATE(actual_solution(1:size(knots_eval)),ierr)

        numerical_error=0.0_f64
        do  test_mode=0,test_power_two-1,1
            print *, "--Test Mode: ", 2**test_mode,   " -----------------------------------------"

            testfunction_test_mode=1.0_f64*(2**(test_mode))
!            testfunction_test_mode=1.0_f64
!!              print *, sll_poisson_1d_fem_testfunction(knots)
!!              print *, sin(sll_kx*knots*testfunction_test_mode)*((testfunction_test_mode*sll_kx)**2 )
!!              stop
      !!  print *,sll_poisson_1d_fem_testfunction(1),  sll_poisson_1d_fem_testfunction(1)
            call solver%solve(sll_poisson_1d_fem_testfunction)

            call solver%eval_solution(knots_eval,solution)
            actual_solution=(((2.0_f64**test_mode)*sll_kx)**(-2))*sll_poisson_1d_fem_testfunction(knots_eval)
            numerical_error(test_mode+1)=maxval(abs( (solution - actual_solution)/maxval(abs(actual_solution)) ) )

            call solver%eval_solution_derivative(knots_eval,solution)
            actual_solution=(2.0_f64**test_mode)*sll_kx*cos((2.0_f64**test_mode)*sll_kx*knots_eval)
            numerical_error_deriv(test_mode+1)=maxval(abs( (solution - actual_solution)/maxval(abs(actual_solution)) ) )

            print *, "MAX",maxloc(abs( (solution - actual_solution)))

            H1seminorm=solver%H1seminorm_solution()
            residual=0
            !residual=sll_bspline_fem_solver_1d_calculate_residual()

            print *, "Relative Numerical Error: ", numerical_error(test_mode+1),"  ",numerical_error_deriv(test_mode+1)


            num_error_sum= sum(abs( (solution - actual_solution)))
            print *, "Summed Numerical Error: ", num_error_sum


            print *, "H1-Seminorm of solution: ",&
                H1seminorm, "Rel. Error: " , abs(H1seminorm - 0.5_f64*((2**test_mode)*sll_kx)**2)/(((2**test_mode)*sll_kx)**2)/2.0_f64

            print *, "Residual FEM Equation: ", residual


        enddo
!        do  test_mode=0,test_power_two-1,1
!             print*, numerical_error(test_mode+1) / (4**test_mode)
!        enddo
        do  test_mode=1,test_power_two-1,1
             print*, numerical_error(test_mode+1) / numerical_error(test_mode)
        enddo

        SLL_DEALLOCATE_ARRAY(solution,ierr)
        SLL_DEALLOCATE_ARRAY(actual_solution,ierr)


        call solver%delete(ierr)
        call delete(mesh)
        call delete(mesh_eval)
    endsubroutine

!
!    subroutine test_mass_matrix_assembler(test_power_two)
!        implicit none
!        integer :: ierr
!        sll_real64 :: interval_a=-1, interval_b=1
!        sll_int, intent(in) :: test_power_two
!        sll_real64 :: numerical_error
!        sll_int64 :: idx, jdx
!        sll_int64 :: test_dimension
!
!        sll_real64, dimension(:), allocatable :: circulant_matrix_first_line
!        sll_real64, dimension(:), allocatable :: equ_right_side
!        sll_real64, dimension(:), allocatable :: matrix_product
!        sll_real64, dimension(:), allocatable :: solution
!
!        sll_real64 :: mesh_delta_t
!        sll_real64, dimension(:), allocatable :: knots
!        type(arbitrary_degree_spline_1d), pointer :: arbitrarydegree1D
!        sll_real64, dimension(:) , allocatable:: stiffn_matrix_period
!        type(arb_deg_1d_interpolator) :: arbitrarydegree1D_interpolator
!        sll_real64, dimension(:), allocatable :: interpolant_original
!        sll_real64, dimension(:), allocatable :: interpolant
!        sll_real64,dimension(:), allocatable :: coeff_splines
!        !Initialize random seed
!        !call init_random_seed()
!        !Test Mass Matrix Assembler
!        test_dimension=2**test_power_two
!
!
!        mesh_delta_t = (interval_b - interval_a)/test_dimension
!        SLL_ALLOCATE(knots(test_dimension+1),ierr)
!        knots(1)=interval_a
!        do idx=2,size(knots)
!            knots(idx)=knots(idx-1) + mesh_delta_t
!        enddo
!
!        print *, "Cell Size: ", mesh_delta_t
!        arbitrarydegree1D=>new_arbitrary_degree_spline_1d( 3, knots, size(knots), PERIODIC_ARBITRARY_DEG_SPLINE)
!
!        stiffn_matrix_period=gen_mass_matrix_period (arbitrarydegree1D, knots)
!
!        !call sll_display(stiffn_matrix_period, "(F8.4)")
!
!        !Here the mass matrix is a circulant matrix so we use the solver we just tested
!        SLL_CLEAR_ALLOCATE(equ_right_side(test_dimension),ierr)
!        SLL_CLEAR_ALLOCATE(solution(test_dimension),ierr)
!        equ_right_side(test_dimension/2)=-1
!        equ_right_side(test_dimension/4)=-2
!
!
!
!        !    solution=sll_bspline_fem_solve_poisson_1d(arbitrarydegree1D, knots,equ_right_side)
!
!        !Test solver
!        call arbitrarydegree1D_interpolator%initialize( size(knots), &
!            interval_a, &
!            interval_b, &
!            SLL_PERIODIC, &
!            SLL_PERIODIC, &
!            4)
!
!
!
!        SLL_ALLOCATE(interpolant_original(test_dimension+1),ierr)
!        SLL_CLEAR_ALLOCATE(interpolant(test_dimension+1),ierr)
!        forall (idx=1:test_dimension+1)
!            interpolant_original(idx)=sin( (2*SLL_PI)*(knots(idx)-interval_a)/(interval_b - interval_a)*1_f64 )
!        endforall
!
!
!
!        call arbitrarydegree1D_interpolator%compute_interpolants( interpolant_original)
!
!
!
!
!        call arbitrarydegree1D_interpolator%compute_interpolants(  interpolant_original , knots,size(knots))
!
!
!        SLL_ALLOCATE(coeff_splines(size(knots)-1),ierr)
!
!        coeff_splines=arbitrarydegree1D_interpolator%coeff_splines(1:size(knots)-1)
!
!        !coeff_splines(1)=coeff_splines(1)+ arbitrarydegree1D_interpolator%coeff_splines(size(knots)+1)
!        call sll_display(coeff_splines, "(F10.4)")
!        call sll_display(interpolant_original, "(F10.4)")
!
!
!
!        !do idx=-3,3,1
!        !interpolant=bspline_basis_to_realvals (  arbitrarydegree1D, knots, &
!            ! cshift(coeff_splines,idx), knots )
!        !call sll_display(interpolant-interpolant_original, "(F10.4)")
!        !  print *, idx, ": " ,sum(abs(interpolant-interpolant_original))
!        !call sll_display(interpolant-interpolant_original, "(F10.4)")
!        !enddo
!
!
!
!    endsubroutine

!
!
!
!    !>Test Solver for Circulant Matrix equations
!    subroutine  test_circulant_matrix_eq_solver ( ierr,test_power_two)
!        implicit none
!        integer, intent(inout) :: ierr
!        sll_int, intent(in) :: test_power_two
!        sll_real64 :: interval_a=-1, interval_b=1
!        sll_real64 :: numerical_error
!        sll_int64 :: idx, jdx
!        sll_int64 :: test_dimension
!
!        sll_real64, dimension(:), allocatable :: circulant_matrix_first_line
!        sll_real64, dimension(:), allocatable :: equ_right_side
!        sll_real64, dimension(:), allocatable :: matrix_product
!        sll_real64, dimension(:), allocatable :: solution
!
!
!
!        !Init random seed
!        do idx=2, test_power_two
!            test_dimension=2**idx
!            SLL_CLEAR_ALLOCATE(circulant_matrix_first_line(1:test_dimension),ierr)
!            SLL_CLEAR_ALLOCATE(equ_right_side(1:test_dimension),ierr)
!            SLL_CLEAR_ALLOCATE(solution(1:test_dimension),ierr)
!
!
!            call random_number(circulant_matrix_first_line)
!            call random_number(equ_right_side)
!
!            !Test cases
!            !idx=3
!            !circulant_matrix_first_line=(/1,2,3,0,0,0,6,5/)
!            !equ_right_side=(/ 1,3,-1,0,5,3,0,0/)
!
!            solution= bspline_fem_solver_1d_solve_circulant_matrix_equation( circulant_matrix_first_line, &
!                equ_right_side  )
!            SLL_CLEAR_ALLOCATE(matrix_product(1:test_dimension),ierr)
!            do jdx=1, test_dimension
!                matrix_product(jdx)=dot_product( cshift(circulant_matrix_first_line, -(jdx-1)), solution  )
!            enddo
!
!            !call sll_display( matrix_product-equ_right_side , "(F10.5)")
!            numerical_error=sum(abs(matrix_product-equ_right_side))/sum(abs(equ_right_side))
!
!            if(numerical_error .lt. SLL_CIRCULAR_MATRIX_SOLVER_TOL ) then
!                print *, 'Dimension: ', test_dimension,  ' PASSED'
!            else
!                if (test_dimension<20) call sll_display(solution, "(F10.5)")
!                print *, 'FAILED'
!                print *, 'Error: ', numerical_error
!                print *, 'Dimension: ', test_dimension
!                stop(1)
!                stop
!            end if
!
!            SLL_DEALLOCATE_ARRAY(circulant_matrix_first_line,ierr)
!            SLL_DEALLOCATE_ARRAY(equ_right_side,ierr)
!            SLL_DEALLOCATE_ARRAY(solution,ierr)
!            SLL_DEALLOCATE_ARRAY(matrix_product,ierr)
!
!        enddo
!
!
!
!    endsubroutine
!
!
!
!    !>Subroutine taken from http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
!    subroutine init_random_seed()
!        implicit none
!        integer, allocatable :: seed(:)
!        integer :: i, n, un, istat, dt(8), pid, t(2), s
!        integer(8) :: count, tms
!        un=201
!        call random_seed(size = n)
!        allocate(seed(n))
!        ! First try if the OS provides a random number generator
!        open(unit=un, file="/dev/urandom", access="stream", &
!            form="unformatted", action="read", status="old", iostat=istat)
!        if (istat == 0) then
!            read(un) seed
!            close(un)
!        else
!            ! Fallback to XOR:ing the current time and pid. The PID is
!            ! useful in case one launches multiple instances of the same
!            ! program in parallel.
!            call system_clock(count)
!            if (count /= 0) then
!                t = transfer(count, t)
!            else
!                call date_and_time(values=dt)
!                tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
!                    + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
!                    + dt(3) * 24 * 60 * 60 * 60 * 1000 &
!                    + dt(5) * 60 * 60 * 1000 &
!                    + dt(6) * 60 * 1000 + dt(7) * 1000 &
!                    + dt(8)
!                t = transfer(tms, t)
!            end if
!            s = ieor(t(1), t(2))
!            pid = getpid() + 1099279 ! Add a prime
!            s = ieor(s, pid)
!            if (n >= 3) then
!                seed(1) = t(1) + 36269
!                seed(2) = t(2) + 72551
!                seed(3) = pid
!                if (n > 3) then
!                    seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
!                end if
!            else
!                seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
!            end if
!        end if
!        call random_seed(put=seed)
!    end subroutine init_random_seed
end program
