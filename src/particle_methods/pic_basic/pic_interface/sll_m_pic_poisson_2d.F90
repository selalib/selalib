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
       sll_f_new_pic_poisson_2d

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
       procedure :: initialize => initialize_pic_poisson_2d
       procedure :: delete => delete_pic_poisson_2d


    end type sll_t_pic_poisson_2d

contains

  subroutine add_charge_single_2d(this, position, weight)
    class(sll_t_pic_poisson_2d), intent( inout ) :: this !< Pic Poisson solver object
    sll_real64,                intent( in ) :: position(this%dim) !< Position of the particle
    sll_real64,                intent( in ) :: weight !< Weight of the particle

    call this%kernel%add_charge( position, weight, this%rho_dofs_local )
       

  end subroutine add_charge_single_2d
  

  subroutine evaluate_rho_single_2d(this, position, func_value)
    class(sll_t_pic_poisson_2d), intent( inout ) :: this !< Pic Poisson solver object
    sll_real64,                intent( in ) :: position(this%dim) !< Position of the particle
    sll_real64, intent(out) :: func_value !< Value of rho at given position

    if (this%rho_collected .EQV. .FALSE.) then
       this%rho_collected = .TRUE.
       this%rho_dofs = 0.0_f64
       call sll_o_collective_allreduce( sll_v_world_collective, this%rho_dofs_local, &
            this%no_dofs, MPI_SUM, this%rho_dofs)
    end if

    call this%kernel%evaluate( position, this%rho_dofs, func_value)

  end subroutine evaluate_rho_single_2d


  subroutine evaluate_phi_single_2d(this, position, func_value)
    class(sll_t_pic_poisson_2d), intent( inout ) :: this !< Pic Poisson solver object
    sll_real64,                intent( in ) :: position(this%dim) !< Position of the particle
    sll_real64, intent(out) :: func_value !< Value of phi at given position

    call this%kernel%evaluate( position, this%phi_dofs, func_value)

  end subroutine evaluate_phi_single_2d

  subroutine evaluate_field_single_2d(this, position, components, func_value)
    class(sll_t_pic_poisson_2d), intent( inout ) :: this !< Pic Poisson solver object
    sll_int32, intent(in ) :: components(:)
    sll_real64,                intent( in ) :: position(this%dim) !< Position of the particle
    sll_real64, intent(out) :: func_value(:)

    call this%kernel%evaluate_multiple( position, components, this%efield_dofs, &
         func_value)

  end subroutine evaluate_field_single_2d


  subroutine solve_2d(this)
    class(sll_t_pic_poisson_2d), intent( inout ) :: this !< Pic Poisson solver object

    call this%solve_phi()
    call this%solve_fields()
    
  end subroutine solve_2d

  subroutine solve_phi_2d(this)
    class(sll_t_pic_poisson_2d), intent( inout ) :: this !< Pic Poisson solver object

    if (this%rho_collected .EQV. .FALSE.) then
       this%rho_collected = .TRUE.
       this%rho_dofs = 0.0_f64
       call sll_o_collective_allreduce( sll_v_world_collective, this%rho_dofs_local, &
            this%no_dofs, MPI_SUM, this%rho_dofs)
    end if
    this%rho2d = reshape(this%rho_dofs, this%no_gridpts)
    call this%solver%compute_phi_from_rho(this%phi2d, this%rho2d)
    this%phi_dofs = reshape(this%phi2d, [this%no_dofs])

  end subroutine solve_phi_2d

  subroutine solve_fields_2d(this)
    class(sll_t_pic_poisson_2d), intent( inout ) :: this !< Pic Poisson solver object

    if (this%rho_collected .EQV. .FALSE.) then
       this%rho_collected = .TRUE.
       this%rho_dofs = 0.0_f64
       call sll_o_collective_allreduce( sll_v_world_collective, this%rho_dofs_local, &
            this%no_dofs, MPI_SUM, this%rho_dofs)
    end if
    this%rho2d = reshape(this%rho_dofs, this%no_gridpts)
    call this%solver%compute_E_from_rho(this%efield1, this%efield2, this%rho2d)
    this%efield_dofs(:,1) = reshape(this%efield1, [this%no_dofs])
    this%efield_dofs(:,2) = reshape(this%efield2, [this%no_dofs])

  end subroutine solve_fields_2d

  subroutine reset_2d(this)
    class(sll_t_pic_poisson_2d), intent( inout ) :: this !< Pic Poisson solver object

    this%rho_dofs_local = 0.0_f64
    this%rho_dofs = 0.0_f64
    this%rho_collected = .FALSE.

  end subroutine reset_2d

  subroutine add_analytic_charge_2d(this, factor_present, factor_analytic)   
    class(sll_t_pic_poisson_2d), intent( inout ) :: this !< Pic Poisson solver object
    sll_real64, intent( in ) :: factor_present !< Factor to multiply accumulated charge with
    sll_real64, intent( in ) :: factor_analytic !< Factor to multiply added analytic charge with

    this%rho_dofs = factor_present * this%rho_dofs + &
         factor_analytic * this%rho_analyt_dofs

  end subroutine add_analytic_charge_2d

  subroutine set_analytic_charge_2d(this, func)
    class( sll_t_pic_poisson_2d ), intent( inout )    :: this !< PIC Poisson solver object.
    procedure(sll_f_function_of_position)                :: func !< Function to be projected.

    call this%solver%compute_rhs_from_function(func, this%rho_analyt_dofs)

  end subroutine set_analytic_charge_2d

  function compute_field_energy_2d(this, component) result(energy)
    class (sll_t_pic_poisson_2d), intent( in ) :: this
    sll_int32, intent( in ) :: component !< Component of the electric field for which the energy should be computed
    sll_real64 :: energy !< L2 norm squarred of 


    if (component == 1) then
       energy = this%solver%l2norm_squared(this%efield1)
    elseif (component == 2) then
       energy = this%solver%l2norm_squared(this%efield2)
    end if

  end function compute_field_energy_2d

  !-------------------------------------------------------------------------------------------
  !< Constructor 
  subroutine initialize_pic_poisson_2d(this, no_gridpts, solver, kernel)
    class( sll_t_pic_poisson_2d), intent(out) :: this
    sll_int32, intent(in) :: no_gridpts(2)
    class( sll_c_poisson_2d_base), pointer, intent(in) :: solver
    class( sll_c_kernel_smoother), pointer, intent(in)   :: kernel !< kernel smoother object

    !local variables
    sll_int32 :: ierr

    this%dim = 1
    this%no_gridpts = no_gridpts
    this%no_dofs = product(no_gridpts)
    this%rho_collected = .FALSE.
    
    this%solver => solver
    this%kernel => kernel

    allocate(this%rho_dofs(this%no_dofs), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(this%rho_dofs_local(this%no_dofs), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(this%rho_analyt_dofs(this%no_dofs), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(this%efield_dofs(this%no_dofs, 2), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(this%phi_dofs(this%no_dofs), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(this%rho2d(this%no_gridpts(1), this%no_gridpts(2)), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(this%efield1(this%no_gridpts(1), this%no_gridpts(2)), stat=ierr)
    SLL_ASSERT( ierr == 0 )
    allocate(this%efield2(this%no_gridpts(1), this%no_gridpts(2)), stat=ierr)  
    SLL_ASSERT( ierr == 0 )  
    allocate(this%phi2d(this%no_gridpts(1), this%no_gridpts(2)), stat=ierr)
    SLL_ASSERT( ierr == 0 )

  end subroutine initialize_pic_poisson_2d

  
  subroutine delete_pic_poisson_2d(this)
    class( sll_t_pic_poisson_2d), intent(inout) :: this

    deallocate(this%rho_dofs)
    deallocate(this%rho_dofs_local)
    deallocate(this%rho_analyt_dofs)
    deallocate(this%efield_dofs)
    deallocate(this%phi_dofs)
    deallocate(this%rho2d)
    deallocate(this%efield1)
    deallocate(this%efield2)
    deallocate(this%phi2d)
    deallocate(this%solver)
    deallocate(this%kernel)

  end subroutine delete_pic_poisson_2d


  !< Constructor 
  function sll_f_new_pic_poisson_2d(no_gridpts, solver, kernel) result (poisson_solver)
    sll_int32, intent(in) :: no_gridpts(2)
    class( sll_c_poisson_2d_base), pointer, intent(in) :: solver
    class( sll_c_kernel_smoother), pointer, intent(in)   :: kernel !< kernel smoother object
    class( sll_t_pic_poisson_2d), pointer :: poisson_solver

    !local variables
    sll_int32 :: ierr
 
    allocate( poisson_solver, stat=ierr)

    call initialize_pic_poisson_2d(poisson_solver, no_gridpts, solver, kernel)

  end function sll_f_new_pic_poisson_2d

end module sll_m_pic_poisson_2d
