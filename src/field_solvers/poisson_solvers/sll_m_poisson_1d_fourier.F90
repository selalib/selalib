!**************************************************************
!  Author: Jakob Ameres, jakob.ameres@tum.de
!**************************************************************

!> Module to solve Poisson equation on one dimensional mesh using Finite Elements
module sll_m_poisson_1d_fourier
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_1d

  use sll_m_constants, only: &
    sll_p_i1, &
    sll_p_kx, &
    sll_p_pi

  implicit none

  public :: &
    sll_f_new_poisson_1d_fourier, &
    sll_t_poisson_1d_fourier

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
    type :: sll_t_poisson_1d_fourier
        class(sll_t_cartesian_mesh_1d), private, pointer  :: cartesian_mesh

        !> Fourier modes
        sll_int32 , private :: num_modes !<Number of fourier modes
        sll_real64, private :: Ilength
        sll_real64, private :: eta_min

        sll_comp64,  dimension(:), allocatable :: fourier_fmode


    contains
        procedure,  pass(self) :: initialize =>sll_initialize_poisson_1d_fourier
        !procedure,  pass(poisson) :: new =>initialize_poisson_1d_fourier
        !procedure,  pass(poisson) :: initialize =>initialize_poisson_1d_fourier
        procedure,  pass(self) :: delete =>sll_delete_poisson_1d_fourier
        procedure,  pass(self) :: solve =>solve_poisson_1d_fourier_rhs
        procedure,  pass(self) :: eval_solution=>poisson_1d_fourier_eval_solution
        procedure,  pass(self) :: eval_solution_derivative=>poisson_1d_fourier_eval_solution_derivative

        procedure, pass(self) :: H1seminorm_solution=>poisson_1d_fourier_H1seminorm_solution
        procedure, pass(self) :: L2norm_solution=>poisson_1d_fourier_L2norm_solution
        procedure, pass(self) ::  get_rhs_from_klimontovich_density_weighted=>&
                                            poisson_1d_fourier_get_rhs_from_klimontovich_density_weighted
    end type sll_t_poisson_1d_fourier

    !Interface for one dimensional function for right hand side
    abstract interface
        function poisson_1d_fourier_rhs_function(x) result(y)
            use sll_m_working_precision
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
    !         module procedure sll_f_new_poisson_1d_fourier
    !    endinterface

contains

    !>Destructor
    subroutine sll_delete_poisson_1d_fourier(self,ierr)
        class(sll_t_poisson_1d_fourier),intent(inout) :: self     !< Solver data structure
        sll_int32, intent(out)                :: ierr    !< error code
        SLL_DEALLOCATE_ARRAY(self%fourier_fmode,ierr)

    endsubroutine


    function  sll_f_new_poisson_1d_fourier(cartesian_mesh_1d,num_modes, bc_type,ierr) &
            result(solver)
        type(sll_t_poisson_1d_fourier), pointer :: solver     !< Solver data structure
        class(sll_t_cartesian_mesh_1d), intent(in),pointer  :: cartesian_mesh_1d !< Logical mesh
        sll_int32, intent(out)                :: ierr    !< error code
        sll_int32, intent(in)                :: num_modes !<Degree of the finite differences approximation
        sll_int32, intent(in)               :: bc_type !< type of boundary connditions

        SLL_ALLOCATE(solver,ierr)
        call sll_initialize_poisson_1d_fourier(solver,  cartesian_mesh_1d,num_modes,bc_type,ierr)
    endfunction


    subroutine sll_initialize_poisson_1d_fourier(self, cartesian_mesh_1d,num_modes,bc_type, ierr)
        class(sll_t_poisson_1d_fourier),intent(inout) :: self     !< Solver data structure
        class(sll_t_cartesian_mesh_1d), intent(in),pointer  :: cartesian_mesh_1d !< Logical mesh
        sll_int32, intent(out)                :: ierr    !< error code
        sll_int32, intent(in)                :: num_modes !<Degree of the bsplines
        sll_int32, intent(in)               :: bc_type !< type of boundary connditions
        ierr=0
        !self%boundarycondition=bc_type !sll_p_periodic

        self%cartesian_mesh=>cartesian_mesh_1d
        self%Ilength=self%cartesian_mesh%eta_max-self%cartesian_mesh%eta_min
        self%eta_min=self%cartesian_mesh%eta_min

        self%num_modes=num_modes

        !        selectcase(self%boundarycondition)
        !            case(sll_p_periodic)
        !
        !
        !            case(sll_p_dirichlet)
        !
        !        endselect
        SLL_ALLOCATE(self%fourier_fmode(self%num_modes),ierr)
        return

        print*,bc_type !sll_p_periodic

    endsubroutine

    !    subroutine poisson_1d_fourier_set_solution(self, solution_vector)
    !        class(sll_t_poisson_1d_fourier),intent(inout) :: self     !< Solver data structure
    !        sll_real64, dimension(:) :: solution_vector
    !        SLL_ASSERT(size(solution_vector)==self%num_cells)
    !
    !        self%fd_solution=solution_vector
    !    endsubroutine

    !    !<Calculates the inhomogenity rhs with given function f
    !    !<by Gauss Legendre integration
    !    function sll_poisson_1d_fourier_get_rhs_from_function(self, eval_function) &
        !            result( rhs )
    !        implicit none
    !        class(sll_t_poisson_1d_fourier),intent(in) :: self     !< Solver data structure
    !        procedure (poisson_1d_fourier_rhs_function) :: eval_function
    !        sll_real64, dimension(self%num_cells ) :: rhs !<Right hand side
    !        sll_int32 :: ierr=0
    !        sll_real64, dimension(self%num_cells+1) :: evalpoints
    !        evalpoints=eval_function(nodes_cartesian_mesh_1d(self%cartesian_mesh))
    !        !Map to right hand side according to boundary condition
    !
    !        rhs=evalpoints(1:self%num_cells)
    !        selectcase(self%boundarycondition)
    !            case(sll_p_periodic)
    !            case(sll_p_dirichlet)
    !                rhs(1)=0
    !        endselect
    !
    !    endfunction






    subroutine solve_poisson_1d_fourier_rhs_and_get(self, field, rhs)
        class(sll_t_poisson_1d_fourier),intent(inout) :: self
        sll_comp64, dimension(:), intent(in)      :: rhs
        sll_comp64, dimension(self%num_modes), intent(out)  :: field
        sll_real64 :: coeff
        coeff=2.0_f64*sll_p_pi/self%Ilength

        call  solve_poisson_1d_fourier_rhs(self, rhs)

        field=-self%fourier_fmode/coeff/ &
               cmplx(solve_poisson_1d_fourier_get_modes(self),0.0_f64,kind=f64) &
               /sll_p_i1

    end subroutine


    function solve_poisson_1d_fourier_get_modes(self)  result(modes)
        class(sll_t_poisson_1d_fourier),intent(inout) :: self
        sll_int32, dimension(self%num_modes) :: modes
        sll_int32 :: idx

        do idx=1,self%num_modes
            modes(idx)=idx
        enddo
    endfunction

    subroutine solve_poisson_1d_fourier_rhs(self, rhs)
        class(sll_t_poisson_1d_fourier),intent(inout) :: self
        sll_comp64, dimension(:), intent(in)      :: rhs
        SLL_ASSERT(size(rhs)==self%num_modes)
        SLL_ASSERT(size(rhs)==size(self%fourier_fmode))

        !Nothing to be done here
        !rhs=self%fourier_fmode
        self%fourier_fmode=rhs
    endsubroutine

    !
    !    !> @brief Solves the poisson equation for a given right hand side function
    !    !> @param self pointer to a sll_t_poisson_1d_fourier object.
    !    !> @param rhs right hand side function of type poisson_1d_fourier_rhs_function
    !    subroutine solve_poisson_1d_fourier_functionrhs(self, rhs_fun)
    !        implicit none
    !        class(sll_t_poisson_1d_fourier),intent(inout) :: self
    !        procedure (poisson_1d_fourier_rhs_function) :: rhs_fun
    !        sll_real64, dimension(self%num_cells) :: rhs
    !
    !        rhs=sll_poisson_1d_fourier_get_rhs_from_function(self, rhs_fun)
    !        call solve_poisson_1d_fourier_rhs(self,rhs  )
    !    endsubroutine




    !< Evaluates the first derivative of the solution at the given points knots_eval
    !< The Result is written into eval_solution
    subroutine poisson_1d_fourier_eval_solution(self, knots_eval, eval_solution)
        class(sll_t_poisson_1d_fourier),intent(inout) :: self
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        sll_int32, dimension(self%num_modes) :: fmode
        sll_int32 :: idx
        sll_real64 :: coeff
        coeff=2.0_f64*sll_p_pi/self%Ilength

        SLL_ASSERT(size(knots_eval)==size(eval_solution))
        do idx=1,self%num_modes
            fmode(idx)=idx
        enddo

        do idx=1, size(knots_eval)
                eval_solution(idx)= -2.0_f64* sum( &
                    real(  self%fourier_fmode/(coeff*fmode)**2 )*cos(knots_eval(idx)*fmode*coeff) &
                     -aimag( self%fourier_fmode/(coeff*fmode)**2 )*sin(knots_eval(idx)*fmode*coeff))

        enddo

    endsubroutine

    !< Evaluates the solution at the given points knots_eval
    !< The Result is written into eval_solution
    subroutine poisson_1d_fourier_eval_solution_derivative(self, knots_eval, eval_solution)
        class(sll_t_poisson_1d_fourier),intent(inout) :: self
        sll_real64, dimension(:), intent(in)     :: knots_eval
        sll_real64, dimension(:), intent(out)     :: eval_solution
        sll_int32, dimension(self%num_modes) :: fmode
        sll_int32 :: idx
        sll_real64 :: coeff
        coeff=2.0_f64*sll_p_pi/self%Ilength

        SLL_ASSERT(size(knots_eval)==size(eval_solution))

        do idx=1,self%num_modes
            fmode(idx)=idx
        enddo
        do idx=1, size(knots_eval)
                eval_solution(idx)= -2.0_f64* sum( &
                    real(   self%fourier_fmode/(coeff*fmode*sll_p_i1) )*cos(knots_eval(idx)*fmode*coeff) &
                    - aimag( self%fourier_fmode/(coeff*fmode*sll_p_i1) )*sin(knots_eval(idx)*fmode*coeff))

        enddo

    endsubroutine

    !    !<Gives the squared H1-seminorm of the solution $\Phi$: $|\nabla \Phi|^2$
    function  poisson_1d_fourier_H1seminorm_solution(self) result(seminorm)
        class(sll_t_poisson_1d_fourier),intent(inout) :: self
        sll_real64 :: seminorm
        sll_int32 :: fmode_a!, fmode_b
        sll_real64 :: coeff
        coeff=2.0_f64*sll_p_pi/self%Ilength

        seminorm=0.0_f64
        do fmode_a=1,self%num_modes
!!            do fmode_b=1,self%num_modes
!            seminorm=seminorm+ 2.0_f64*abs( (self%fourier_fmode(fmode_a)/coeff/fmode_a/sll_p_i1)&
!                                *(self%fourier_fmode(fmode_b)/coeff/fmode_b/sll_p_i1))
            seminorm=seminorm +  2.0_f64*abs( (self%fourier_fmode(fmode_a)/coeff/fmode_a/sll_p_i1))**2


            !seminorm=seminorm+ 2.0_f64*real((self%fourier_fmode(fmode_a)/coeff/fmode_a/sll_p_i1 )**2)/self%Ilength
!            seminorm=seminorm+  (real(self%fourier_fmode(fmode_a))**2 -aimag(self%fourier_fmode(fmode_a))**2) &
!                                  *(1.0_f64/fmode_a)**2/coeff**2/self%Ilength*2
!            seminorm=seminorm=
            !seminorm=seminorm+  (real(self%fourier_fmode(fmode_a))**2 +aimag(self%fourier_fmode(fmode_a))**2)
            !seminorm=seminorm+  2.0_f64*real(self%fourier_fmode(fmode_a)*conjg(self%fourier_fmode(fmode_a)))
        enddo
    endfunction
    !
    !<Gives the squared L2-norm of the solution $\Phi$: $|\Phi|^2$
    function  poisson_1d_fourier_L2norm_solution(self) result(l2norm)
        class(sll_t_poisson_1d_fourier),intent(inout) :: self
        sll_real64 :: l2norm
        sll_int32 :: fmode
        sll_real64 :: coeff
        coeff=2.0_f64*sll_p_pi/self%Ilength
        do fmode=1,self%num_modes
            l2norm=l2norm+ real((self%fourier_fmode(fmode)/(coeff*fmode*sll_p_i1)**2 )**2*2.0_f64)
        enddo
    endfunction


    function poisson_1d_fourier_get_rhs_from_klimontovich_density(self, &
            ppos)  result(rhs)
        class(sll_t_poisson_1d_fourier),intent(inout) :: self
        sll_real64, dimension(:), intent(in) ::ppos
        sll_comp64, dimension(self%num_modes) :: rhs
        sll_int32 ::fmode
        do fmode=1,self%num_modes
            rhs(fmode)=sum(exp(-sll_p_i1 * fmode*ppos*sll_p_kx/self%Ilength));
        enddo
    endfunction



    function poisson_1d_fourier_get_rhs_from_klimontovich_density_weighted( self,&
            ppos, pweight) result(rhs)
        class(sll_t_poisson_1d_fourier),intent(inout) :: self
        sll_real64, dimension(:), intent(in) ::ppos
        sll_real64, dimension(:), intent(in) ::pweight
        sll_comp64, dimension(self%num_modes) :: rhs
        sll_int32 :: fmode

        SLL_ASSERT(size(ppos)==size(pweight))
        do fmode=1,self%num_modes
            !Be careful here, the dot_product tends to complex conjugate stuff
            !which we don't want in self case
            !rhs(fmode)=dot_product(exp(-sll_p_i1*fmode*ppos*2.0_f64*sll_p_pi/self%Ilength), pweight )
            rhs(fmode)= &
          sum(exp(-fmode*sll_p_i1*ppos*sll_p_kx/self%Ilength)&
             * cmplx(pweight,0.0_f64,kind=f64))

        enddo
    endfunction



    !    function poisson_1d_fourier_calculate_residual(self) result(residuum)
    !                class(sll_t_poisson_1d_fourier),intent(inout) :: self
    !
    !        sll_real64 :: residuum
    !        residuum= sqrt(sum(( &
        !            poisson_1d_fourier_circulant_matrix_vector_product&
        !            ( self%fd_matrix_first_line, self%fd_solution)&
        !            - self%fem_inhomogenity &
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

end module sll_m_poisson_1d_fourier
