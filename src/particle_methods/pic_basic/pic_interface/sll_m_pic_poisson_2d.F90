!> @ingroup pic_interface
!> @author Katharina Kormann, IPP
!> @brief Factory method for Poisson solver for particle methods in 2d build from 2d Poisson solver and a kernel smoother.
!> @details 

module sll_m_pic_poisson_2d
#include "sll_assert.h"
#include "sll_working_precision.h"

  use sll_m_pic_poisson_base, only : &
       sll_c_pic_poisson
  use sll_m_kernel_smoother_base, only : &
       sll_c_kernel_smoother
  use sll_m_poisson_2d_base, only : &
       sll_c_poisson_2d_base, sll_f_function_of_position
  use sll_m_collective, only : &
       sll_v_world_collective, sll_o_collective_allreduce
  use sll_mpi, only: &
       MPI_SUM


  implicit none

  public :: sll_t_pic_poisson_2d, &
       sll_f_new_pic_poisson_2d, &
       sll_s_new_pic_poisson_2d

  private

  !> PIC Poisson solver 2d 
  type, public, extends(sll_c_pic_poisson) :: sll_t_pic_poisson_2d
     
     sll_int32 :: no_gridpts(2)
     sll_int32 :: no_dofs

     class(sll_c_kernel_smoother), pointer :: kernel
     class(sll_c_poisson_2d_base),   pointer :: solver
     sll_real64, allocatable               :: rho_dofs(:)
     sll_real64, allocatable               :: rho_dofs_local(:)
     sll_real64, allocatable               :: rho_analyt_dofs(:)
     sll_real64, allocatable               :: efield_dofs(:,:)
     sll_real64, allocatable               :: phi_dofs(:)
     sll_real64, allocatable               :: rho2d(:,:)
     sll_real64, allocatable               :: efield1(:,:)
     sll_real64, allocatable               :: efield2(:,:)
     sll_real64, allocatable               :: phi2d(:,:)

     logical                               :: rho_collected

     contains
       procedure :: add_charge_single => add_charge_single_2d !< Add contribution of one particle to the charge density.
       procedure :: reset => reset_2d !< Reset accumulated charge density to zero.
       procedure :: evaluate_rho_single => evaluate_rho_single_2d !< Evaluate charge density at given position.
       procedure :: evaluate_phi_single => evaluate_phi_single_2d !< Evaluate potential at given position.
       procedure :: evaluate_field_single => evaluate_field_single_2d !< Evaluate field components at given positions.
       procedure :: solve => solve_2d !< Solve Poisson's equation for potential and E-fields
       procedure :: solve_phi => solve_phi_2d !< Solve Poisson's equation for potential.
       procedure :: solve_fields => solve_fields_2d !< Solve for the electric field.
       procedure :: add_analytic_charge => add_analytic_charge_2d !< !< Set charge as linear combination of previously accumulated charge and previously set analytic charge.
       procedure :: set_analytic_charge => set_analytic_charge_2d  !< Set the value of the analytic charge contribution from a given function.
       procedure :: compute_field_energy => compute_field_energy_2d !< Compute the field energy.
       procedure :: init => init_pic_poisson_2d !< Initialize the type
       procedure :: free => free_pic_poisson_2d !< finalization


    end type sll_t_pic_poisson_2d

contains

  !> Add charge from one particle
  subroutine add_charge_single_2d(self, position, weight)
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object
    sll_real64,                intent( in ) :: position(self%dim) !< Position of the particle
    sll_real64,                intent( in ) :: weight !< Weight of the particle

    call self%kernel%add_charge( position, weight, self%rho_dofs_local )
       

  end subroutine add_charge_single_2d
  
  !> Evaluate charge density \a rho at one position
  subroutine evaluate_rho_single_2d(self, position, func_value)
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object
    sll_real64,                intent( in ) :: position(self%dim) !< Position of the particle
    sll_real64, intent(out) :: func_value !< Value of rho at given position

    if (self%rho_collected .EQV. .FALSE.) then
       self%rho_collected = .TRUE.
       self%rho_dofs = 0.0_f64
       call sll_o_collective_allreduce( sll_v_world_collective, self%rho_dofs_local, &
            self%no_dofs, MPI_SUM, self%rho_dofs)
    end if

    call self%kernel%evaluate( position, self%rho_dofs, func_value)

  end subroutine evaluate_rho_single_2d

  !> Evaluate potential \a phi at one position
  subroutine evaluate_phi_single_2d(self, position, func_value)
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object
    sll_real64,                intent( in ) :: position(self%dim) !< Position of the particle
    sll_real64, intent(out) :: func_value !< Value of phi at given position

    call self%kernel%evaluate( position, self%phi_dofs, func_value)

  end subroutine evaluate_phi_single_2d

  !> Evaluate components \a components of the electric field as one position
  subroutine evaluate_field_single_2d(self, position, components, func_value)
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object
    sll_int32, intent(in ) :: components(:)
    sll_real64,                intent( in ) :: position(self%dim) !< Position of the particle
    sll_real64, intent(out) :: func_value(:)

    call self%kernel%evaluate_multiple( position, components, self%efield_dofs, &
         func_value)

  end subroutine evaluate_field_single_2d


  !> Solve for phi and fields
  subroutine solve_2d(self)
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object

    call self%solve_phi()
    call self%solve_fields()
    
  end subroutine solve_2d

  !> Solve for potential
  subroutine solve_phi_2d(self)
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object

    if (self%rho_collected .EQV. .FALSE.) then
       self%rho_collected = .TRUE.
       self%rho_dofs = 0.0_f64
       call sll_o_collective_allreduce( sll_v_world_collective, self%rho_dofs_local, &
            self%no_dofs, MPI_SUM, self%rho_dofs)
    end if
    self%rho2d = reshape(self%rho_dofs, self%no_gridpts)
    call self%solver%compute_phi_from_rho(self%phi2d, self%rho2d)
    self%phi_dofs = reshape(self%phi2d, [self%no_dofs])

  end subroutine solve_phi_2d

  !> Solve efields from rho
  subroutine solve_fields_2d(self)
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object

    if (self%rho_collected .EQV. .FALSE.) then
       self%rho_collected = .TRUE.
       self%rho_dofs = 0.0_f64
       call sll_o_collective_allreduce( sll_v_world_collective, self%rho_dofs_local, &
            self%no_dofs, MPI_SUM, self%rho_dofs)
    end if
    self%rho2d = reshape(self%rho_dofs, self%no_gridpts)
    call self%solver%compute_E_from_rho(self%efield1, self%efield2, self%rho2d)
    self%efield_dofs(:,1) = reshape(self%efield1, [self%no_dofs])
    self%efield_dofs(:,2) = reshape(self%efield2, [self%no_dofs])

  end subroutine solve_fields_2d

  !> Reset charge to zero
  subroutine reset_2d(self)
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object

    self%rho_dofs_local = 0.0_f64
    self%rho_dofs = 0.0_f64
    self%rho_collected = .FALSE.

  end subroutine reset_2d

  !> Add analytic charge (set by \a set_analytic_charge ) to the accumulated charge
  subroutine add_analytic_charge_2d(self, factor_present, factor_analytic)   
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object
    sll_real64, intent( in ) :: factor_present !< Factor to multiply accumulated charge with
    sll_real64, intent( in ) :: factor_analytic !< Factor to multiply added analytic charge with

    self%rho_dofs = factor_present * self%rho_dofs + &
         factor_analytic * self%rho_analyt_dofs

  end subroutine add_analytic_charge_2d

  !> Set analytic charge defined by a function \a func obeying the interface \a sll_f_function_of_position
  subroutine set_analytic_charge_2d(self, func)
    class( sll_t_pic_poisson_2d ), intent( inout )    :: self !< PIC Poisson solver object.
    procedure(sll_f_function_of_position)                :: func !< Function to be projected.

    call self%solver%compute_rhs_from_function(func, self%rho_analyt_dofs)

  end subroutine set_analytic_charge_2d

  !> Compute the squared l2 norm of component \a component of the field
  function compute_field_energy_2d(self, component) result(energy)
    class (sll_t_pic_poisson_2d), intent( in ) :: self
    sll_int32, intent( in ) :: component !< Component of the electric field for which the energy should be computed
    sll_real64 :: energy !< L2 norm squarred of 


    if (component == 1) then
       energy = self%solver%l2norm_squared(self%efield1)
    elseif (component == 2) then
       energy = self%solver%l2norm_squared(self%efield2)
    end if

  end function compute_field_energy_2d

  !-------------------------------------------------------------------------------------------
  !< Constructor 
  subroutine init_pic_poisson_2d(self, no_gridpts, solver, kernel)
    class( sll_t_pic_poisson_2d), intent(out) :: self
    sll_int32, intent(in) :: no_gridpts(2)
    class( sll_c_poisson_2d_base), pointer, intent(in) :: solver
    class( sll_c_kernel_smoother), pointer, intent(in)   :: kernel !< kernel smoother object

    !local variables
    sll_int32 :: ierr

    self%dim = 1
    self%no_gridpts = no_gridpts
    self%no_dofs = product(no_gridpts)
    self%rho_collected = .FALSE.
    
    self%solver => solver
    self%kernel => kernel

    allocate(self%rho_dofs(self%no_dofs), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(self%rho_dofs_local(self%no_dofs), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(self%rho_analyt_dofs(self%no_dofs), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(self%efield_dofs(self%no_dofs, 2), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(self%phi_dofs(self%no_dofs), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(self%rho2d(self%no_gridpts(1), self%no_gridpts(2)), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(self%efield1(self%no_gridpts(1), self%no_gridpts(2)), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(self%efield2(self%no_gridpts(1), self%no_gridpts(2)), stat=ierr)  
    SLL_ASSERT( ierr == 0 )  
    allocate(self%phi2d(self%no_gridpts(1), self%no_gridpts(2)), stat=ierr)
    SLL_ASSERT( ierr == 0 )

  end subroutine init_pic_poisson_2d

  !> Destructor
  subroutine free_pic_poisson_2d(self)
    class( sll_t_pic_poisson_2d), intent(inout) :: self

    deallocate(self%rho_dofs)
    deallocate(self%rho_dofs_local)
    deallocate(self%rho_analyt_dofs)
    deallocate(self%efield_dofs)
    deallocate(self%phi_dofs)
    deallocate(self%rho2d)
    deallocate(self%efield1)
    deallocate(self%efield2)
    deallocate(self%phi2d)
    self%solver => null()
    self%kernel => null()

  end subroutine free_pic_poisson_2d


  !> Constructor (legacy version)
  function sll_f_new_pic_poisson_2d(no_gridpts, solver, kernel) result (poisson_solver)
    sll_int32, intent(in) :: no_gridpts(2)
    class( sll_c_poisson_2d_base), pointer, intent(in) :: solver
    class( sll_c_kernel_smoother), pointer, intent(in)   :: kernel !< kernel smoother object
    class( sll_t_pic_poisson_2d), pointer :: poisson_solver

    !local variables
    sll_int32 :: ierr
 
    allocate( poisson_solver, stat=ierr)

    call poisson_solver%init(no_gridpts, solver, kernel)

  end function sll_f_new_pic_poisson_2d


  !> Constructor for abstract type
  subroutine sll_s_new_pic_poisson_2d(poisson_solver, no_gridpts, solver, kernel)    
    class( sll_c_pic_poisson), pointer, intent(out) :: poisson_solver
    sll_int32, intent(in) :: no_gridpts(2)
    class( sll_c_poisson_2d_base), pointer, intent(in) :: solver
    class( sll_c_kernel_smoother), pointer, intent(in)   :: kernel !< kernel smoother object

    !local variables
    sll_int32 :: ierr
 
    allocate( sll_t_pic_poisson_2d :: poisson_solver, stat=ierr)
    select type( poisson_solver)
    type is ( sll_t_pic_poisson_2d )
       call poisson_solver%init( no_gridpts, solver, kernel)
    end select

  end subroutine sll_s_new_pic_poisson_2d

end module sll_m_pic_poisson_2d
