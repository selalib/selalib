#ifndef SLL_CIRCULAR_MATRIX_SOLVER_TOL
#define SLL_CIRCULAR_MATRIX_SOLVER_TOL 1e-10
#endif

program unit_test_poisson_1d_fd
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
    use sll_m_poisson_1d_fd
    use sll_m_cartesian_meshes
    use sll_m_constants

    implicit none

    !integer    :: ierr
    sll_int    :: idx, jdx
    sll_real64 :: testfunction_test_mode

    !call test_poisson_solver(10, 3)
    !stop
    do jdx=1,3
        do idx=jdx+2, 6
            print *, "--------------------------------------------------------------------------------------"

            print *, "Starting test at Spline degree ", jdx ," with ", 2**idx, " Cells"
            !call test_poisson_solver(idx, jdx,1000000)
            call test_momentum(idx, jdx,20000)
            !call test_monte_carlo(idx, jdx,500000)
        enddo

    enddo
    stop
!    do idx=1, 2
!        call test_circulant_matrix_eq_solver(ierr, 18)
!    enddo

contains

    function sll_poisson_1d_fd_testfunction(x) result(y)
        sll_real64, dimension(:), intent(in) :: x
        !sll_real:: test_mode
        sll_real64, dimension(size(x)) :: y
        !test_mode=4
        y= ((testfunction_test_mode*sll_kx)**2 )*sin(testfunction_test_mode*sll_kx*x)
    end function


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


    subroutine test_poisson_solver(test_power_two, spline_degree,npart)
        implicit none
        integer :: ierr
        sll_real64 :: interval_a=0.0_f64, interval_b=1.0_f64
        class(sll_cartesian_mesh_1d), pointer :: mesh=>null()    !Finite Element mesh
        class(sll_cartesian_mesh_1d), pointer :: mesh_eval=>null()  !evaluation mesh
        class(poisson_1d_fd), pointer :: solver=>null()
        sll_int32, intent(in) :: test_power_two
        sll_int32, intent(in) :: spline_degree
        sll_real64, dimension(test_power_two) :: numerical_error
        sll_real64, dimension(test_power_two) :: numerical_error_deriv
        !sll_real64, dimension(test_power_two) :: mode_error
        sll_real64 :: num_error_sum
        sll_int32, intent(in) ::  npart !<Number of particles
        !sll_int32 :: idx, jdx,node
        sll_int32 :: test_dimension

        sll_real64, dimension(:), allocatable :: rhs
        sll_real64, dimension(:), allocatable :: solution

        !sll_real64 :: mesh_delta_t
        sll_real64, dimension(2**test_power_two+1 ) :: knots
        sll_real64, dimension(2**test_power_two +1 )  :: knots_eval
        sll_real64, dimension(:), allocatable :: actual_solution
        sll_real64 :: H1seminorm, residual
        sll_real64, dimension(npart) :: xx
        sll_real64, dimension(npart) :: ww

        integer :: test_mode
        !Test Mass Matrix Assembler
        test_dimension=2**test_power_two

        mesh=>new_cartesian_mesh_1d( test_dimension, interval_a, interval_b )
        knots=mesh%eta1_nodes()

        mesh_eval=>new_cartesian_mesh_1d( test_dimension, interval_a, interval_b )
        knots_eval=mesh_eval%eta1_nodes()

        solver=>new_poisson_1d_fd(mesh, 2,spline_degree, SLL_PERIODIC, ierr)

        call random_number(xx)
        xx=xx*(interval_b- interval_a)+interval_a
        
        SLL_CLEAR_ALLOCATE(rhs(1:size(knots)),ierr)

        SLL_CLEAR_ALLOCATE(solution(1:size(knots_eval)),ierr)
        SLL_CLEAR_ALLOCATE(actual_solution(1:size(knots_eval)),ierr)

                
        numerical_error=0.0_f64
        do  test_mode=0,test_power_two-1,1
            print *, "--Test Mode: ", 2**test_mode,   " -----------------------------------------"

            testfunction_test_mode=1.0_f64*(2**(test_mode))

            ww=sll_poisson_1d_fd_testfunction(xx)/real(npart,8)/(interval_b- interval_a)
            
            
            rhs=solver%get_rhs_klimontovich(xx,ww)  !*(solver%bspline_degree+1)
              
            !rhs=solver%get_rhs_fun(sll_poisson_1d_fd_testfunction)
            
            call solver%set_solution(rhs)
            actual_solution=sll_poisson_1d_fd_testfunction(knots_eval)
            
            !call solver%solve(rhs)
            !actual_solution=sll_poisson_1d_fem_testfunction(knots_eval)/((testfunction_test_mode*sll_kx)**2)            
            call solver%eval_solution(knots_eval, solution)
            
            numerical_error(test_mode+1)=maxval(abs( (solution - actual_solution)/maxval(abs(actual_solution)) ) )
            
            !print *,"rhs=[", rhs,"];"
            !print *,"sol=[", solution,"];"

            
            !print *,"plot([", knots_eval(1: (size(knots_eval)-1):1),"], [", rhs,"]);hold on;"
            !print *,"plot([", knots_eval,"], [", actual_solution,"])"

            !print *,"plot([", rhs,"])"
            !print *,"plot([", solver%fd_solution,"])"
            !            testfunction_test_mode=1.0_f64
!!              print *, sll_poisson_1d_fem_testfunction(knots)
!!              print *, sin(sll_kx*knots*testfunction_test_mode)*((testfunction_test_mode*sll_kx)**2 )
!!              stop
      !!  print *,sll_poisson_1d_fem_testfunction(1),  sll_poisson_1d_fem_testfunction(1)
            !!call solver%solve(sll_poisson_1d_fem_testfunction)

            !call solver%eval_solution(knots_eval,solution)
            !actual_solution=(((2.0_f64**test_mode)*sll_kx)**(-2))*sll_poisson_1d_fem_testfunction(knots_eval)
            !numerical_error(test_mode+1)=maxval(abs( (solution - actual_solution)/maxval(abs(actual_solution)) ) )

            call solver%eval_solution_derivative(knots_eval,solution)
            actual_solution=(2.0_f64**test_mode)*sll_kx*cos((2.0_f64**test_mode)*sll_kx*knots_eval)
            numerical_error_deriv(test_mode+1)=maxval(abs( (solution - actual_solution)/maxval(abs(actual_solution)) ) )

            print *, "MAX",maxloc(abs( (solution - actual_solution)))

            H1seminorm=solver%H1seminorm_solution()
            residual=0.0_8
            !residual=sll_bspline_fem_solver_1d_calculate_residual()

            print *, "Relative Numerical Error: ", numerical_error(test_mode+1),"  ",numerical_error_deriv(test_mode+1)

            num_error_sum= sqrt(sum( (solution - actual_solution)**2))/sqrt(sum( actual_solution**2))
            print *, "Relative l2 Error: ", num_error_sum

            print *, "H1-Seminorm of solution: ",&
                H1seminorm, "Rel. Error: " , abs(H1seminorm - 0.5_f64*((2**test_mode)*sll_kx)**2)/(((2**test_mode)*sll_kx)**2)/2.0_f64

            print *, "Residual FEM Equation: ", residual

        enddo

        do  test_mode=0,test_power_two-1,1
             print*, numerical_error(test_mode+1)! / (4**test_mode)
        enddo
!         do  test_mode=1,test_power_two-1,1
!              print*, numerical_error(test_mode+1) / numerical_error(test_mode)
!         enddo

        SLL_DEALLOCATE_ARRAY(solution,ierr)
        SLL_DEALLOCATE_ARRAY(actual_solution,ierr)

        call solver%delete(ierr)
        !call delete(mesh)
        !call delete(mesh_eval)
    endsubroutine


    subroutine test_momentum(test_power_two, spline_degree,npart)
        implicit none
        integer :: ierr
        sll_real64 :: interval_a=0.0_f64, interval_b=1.0_f64
        class(sll_cartesian_mesh_1d), pointer :: mesh=>null()    !Finite Element mesh
        !class(sll_cartesian_mesh_1d), pointer :: mesh_eval=>null()  !evaluation mesh
        class(poisson_1d_fd), pointer :: solver=>null()
        sll_int32, intent(in) :: test_power_two
        sll_int32, intent(in) :: spline_degree
        sll_real64, dimension(test_power_two) :: numerical_error
        !sll_real64, dimension(test_power_two) :: numerical_error_deriv
        !sll_real64, dimension(test_power_two) :: mode_error
        !sll_real64 :: num_error_sum
        sll_int32, intent(in) ::  npart !<Number of particles
        !sll_int32 :: jdx
        sll_int32 :: test_dimension, test_mode
        sll_real64, dimension(npart) :: xx,ww, Phi, E
        sll_real64, dimension(2**test_power_two) :: rhs

        test_dimension=2**test_power_two

        mesh=>new_cartesian_mesh_1d( test_dimension, interval_a, interval_b )
        solver=>new_poisson_1d_fd(mesh, 2,spline_degree, SLL_PERIODIC, ierr)

        call random_number(xx)
        xx=xx*(interval_b- interval_a)+interval_a

        !print *, "Nodes:" , mesh%eta1_nodes()

        numerical_error=0.0_f64
        do  test_mode=0,test_power_two-2,1
            print *, "--Test Mode: ", 2**test_mode,   " -----------------------------------------"

            testfunction_test_mode=1.0_f64*(2**(test_mode))

            ww=sll_poisson_1d_fd_testfunction(xx)/(interval_b- interval_a)

            rhs=solver%get_rhs_klimontovich(xx(1:npart),ww(1:npart))/npart    

            !call solver%set_solution(rhs)
            call solver%solve(rhs)
            call solver%eval_solution_derivative(xx,E)
            call solver%eval_solution(xx,Phi)

            print *, "Int[E]= ", sum(E*ww)/real(npart,8),"	Int[Phi]= ", sum(Phi*ww)/real(npart,8)
        enddo

        call solver%delete(ierr)

    end subroutine


    subroutine test_monte_carlo(test_power_two, spline_degree,npartmax)
        implicit none
        integer :: ierr
        sll_real64 :: interval_a=0.0_f64, interval_b=1.0_f64
        class(sll_cartesian_mesh_1d), pointer :: mesh=>null()    !Finite Element mesh
        !class(sll_cartesian_mesh_1d), pointer :: mesh_eval=>null()  !evaluation mesh
        class(poisson_1d_fd), pointer :: solver=>null()
        sll_int32, intent(in) :: test_power_two
        sll_int32, intent(in) :: spline_degree
        sll_real64, dimension(test_power_two) :: numerical_error
        !sll_real64, dimension(test_power_two) :: numerical_error_deriv
        !sll_real64, dimension(test_power_two) :: mode_error
        !sll_real64 :: num_error_sum
        sll_int32, intent(in) ::  npartmax !<Number of particles
        sll_int32 :: npart
        !sll_int32 :: node
        sll_int32 :: test_dimension, test_mode
        sll_real64, dimension(npartmax) :: xx
        sll_real64, dimension(npartmax) :: ww
        sll_real64, dimension(2**test_power_two) :: rhs, rhsfun

        test_dimension=2**test_power_two

        mesh=>new_cartesian_mesh_1d( test_dimension, interval_a, interval_b )
        solver=>new_poisson_1d_fd(mesh, 2,spline_degree, SLL_PERIODIC, ierr)

        call random_number(xx)
        xx=xx*(interval_b- interval_a)+interval_a

        !print *, "Nodes:" , mesh%eta1_nodes()

        numerical_error=0.0_f64
        do  test_mode=0,test_power_two-2,1
            print *, "--Test Mode: ", 2**test_mode,   " -----------------------------------------"

            testfunction_test_mode=1.0_f64*(2**(test_mode))
            rhsfun=solver%get_rhs_fun(sll_poisson_1d_fd_testfunction)   
            
            ww=sll_poisson_1d_fd_testfunction(xx)/(interval_b- interval_a)
        
            npart=npartmax

            do while (npart > 100)
   
                rhs=solver%get_rhs_klimontovich(xx(1:npart),ww(1:npart))/npart    
   
                print *, "L2 Error (rhs vector): ", sqrt(sum( (rhs-rhsfun)**2))/sqrt(sum( rhsfun**2))

                npart= ceiling(npart/2.0_f64)

            enddo
        enddo
            
        call solver%delete(ierr)

    end subroutine

end program
