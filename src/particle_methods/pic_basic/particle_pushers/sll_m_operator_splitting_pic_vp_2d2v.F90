!> @ingroup particle_pushers
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on operator splitting for 2d2v Vlasov-Poisson.
!> @details MPI parallelization by domain cloning. Periodic boundaries.
module sll_m_operator_splitting_pic_vp_2d2v
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_errors.h"


 
  use sll_m_particle_group_base
  use sll_m_pic_poisson_base
  use sll_m_operator_splitting
  use sll_m_collective
  use sll_m_control_variate

  implicit none
  private

  public :: sll_new_hamiltonian_splitting_pic_vp_2d2v

  !> Operator splitting type for 2d2v Vlasov-Poisson
  type, public, extends(operator_splitting) :: sll_t_operator_splitting_pic_vp_2d2v
     class(sll_c_pic_poisson), pointer :: solver
!     class(poisson_2d_fft_solver), pointer    :: poisson_solver      !< Poisson solver (TODO: Use a base class here)
!     class(sll_kernel_smoother_base), pointer :: kernel_smoother  !< Kernel smoother
     class(sll_particle_group_base), pointer  :: particle_group    !< Particle group
!!$
!!$     sll_real64, allocatable :: rho_dofs(:)                      !< Degrees of freedom for kernel respresentation of the charge density. 
!!$     sll_real64, allocatable :: rho_dofs_local(:)                !< Processor local part of \a rho_dofs 
!!$     sll_real64, allocatable :: rho_2d(:,:)                      !< 2d representation of \a rho_dofs for use in Poisson solver.
!!$     sll_real64, allocatable :: efield1(:,:) !< 2d representation of electric field, first component for use in Poisson solver. 
!!$     sll_real64, allocatable :: efield2(:,:) !< 2d representation of electric field, second component for use in Poisson solver. 
!!$     sll_real64, pointer     :: efield_dofs(:,:)  !< Values of the electric field at grid points (1d representation).

     ! For version with control variate
     class(sll_t_control_variate), pointer :: control_variate => null()
     sll_int32 :: i_weight

   contains
     procedure :: operatorT => advection_x_pic_vp_2d2v  !< Operator for x advection
     procedure :: operatorV => advection_v_pic_vp_2d2v  !< Operator for v advection
     procedure :: strang_splitting => strang_splitting_pic_vp_2d2v !< Strang splitting

     procedure :: initialize => initialize_operator_splitting_pic_vp_2d2v
  end type sll_t_operator_splitting_pic_vp_2d2v

contains

  !---------------------------------------------------------------------------!
  !> Strang splitting
  subroutine strang_splitting_pic_vp_2d2v(this, dt)
    class(sll_t_operator_splitting_pic_vp_2d2v), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step


    call this%operatorT(0.5_f64*dt)
    call this%operatorV(dt)
    call this%operatorT(0.5_f64*dt)


  end subroutine strang_splitting_pic_vp_2d2v
  
  !---------------------------------------------------------------------------!


  !---------------------------------------------------------------------------!
  !> Push x 
  subroutine advection_x_pic_vp_2d2v(this, dt)
    class(sll_t_operator_splitting_pic_vp_2d2v), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step


    !local variables
    sll_int32 :: i_part
    sll_real64 :: x_new(3), vi(3)
    sll_real64 :: wi(this%particle_group%n_weights)

    ! X_new = X_old + dt * V
    do i_part=1,this%particle_group%n_particles
       x_new = this%particle_group%get_x(i_part) + dt * this%particle_group%get_v(i_part) 
       call this%particle_group%set_x(i_part, x_new)

       if (this%particle_group%n_weights == 3) then
          vi = this%particle_group%get_v(i_part)
          ! Update weights if control variate
          wi = this%particle_group%get_weights(i_part)
          wi(3) = this%control_variate%update_df_weight( x_new, vi, 0.0_f64, wi(1), wi(2))
          call this%particle_group%set_weights(i_part, wi)
       end if
    end do


  end subroutine advection_x_pic_vp_2d2v
  
  !---------------------------------------------------------------------------!
  ! Push v
  subroutine advection_v_pic_vp_2d2v(this, dt)
    class(sll_t_operator_splitting_pic_vp_2d2v), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: v_new(3)
    sll_real64 :: efield(2)
    sll_real64 :: qm
    sll_real64 :: xi(3)
    sll_real64 :: wi


    ! Assemble right-hand-side
    call this%solver%reset()
    do i_part = 1, this%particle_group%n_particles
       xi = this%particle_group%get_x(i_part)
       wi = this%particle_group%get_charge(i_part, this%i_weight)
       call this%solver%add_charge(xi(1:2), wi)
    end do

    call this%solver%solve_fields()

    ! V_new = V_old + dt *q/m* E
    qm = this%particle_group%species%q_over_m();
    do i_part=1,this%particle_group%n_particles
       ! Evaluate efields at particle position
       xi = this%particle_group%get_x(i_part)
       call this%solver%evaluate_field_single(xi(1:2), [1,2], efield)
       v_new = this%particle_group%get_v(i_part)
       v_new(1:2) = v_new(1:2) + dt * qm * efield
       call this%particle_group%set_v(i_part, v_new)

       ! Update particle weights
    end do
    

  end subroutine advection_v_pic_vp_2d2v
  
  !---------------------------------------------------------------------------!
  !> Initialization function
  subroutine initialize_operator_splitting_pic_vp_2d2v(this, solver, particle_group)
    class(sll_t_operator_splitting_pic_vp_2d2v), intent(out) :: this !< object 
    class(sll_c_pic_poisson),pointer, intent(in) :: solver !< poisson solver
    class(sll_particle_group_base),pointer, intent(in) :: particle_group !< particle group

    !local variables
    sll_int32 :: ierr

    this%solver => this%solver
    !this%poisson_solver => poisson_solver
    !this%kernel_smoother => kernel_smoother
    this%particle_group => particle_group

!!$    SLL_ALLOCATE(this%rho_dofs(this%kernel_smoother%n_dofs), ierr)
!!$    SLL_ALLOCATE(this%rho_2d(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)
!!$    SLL_ALLOCATE(this%efield1(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)   
!!$    SLL_ALLOCATE(this%efield2(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)
!!$    SLL_ALLOCATE(this%efield_dofs(this%kernel_smoother%n_dofs, 2), ierr)

  end subroutine initialize_operator_splitting_pic_vp_2d2v


  !---------------------------------------------------------------------------!
  !> Constructor.
  function sll_new_hamiltonian_splitting_pic_vp_2d2v &
       (solver, &
       particle_group, &
       control_variate, &
       i_weight) result(this)
    class(sll_t_operator_splitting_pic_vp_2d2v), pointer :: this !< time splitting object 
    class(sll_c_pic_poisson),pointer, intent(in) :: solver !< Poisson solver
    class(sll_particle_group_base),pointer, intent(in) :: particle_group !< Particle group
    class(sll_t_control_variate), optional, pointer, intent(in) :: control_variate
    sll_int32, optional, intent(in) :: i_weight
 
    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE(this, ierr)

    this%solver => solver
    !this%poisson_solver => poisson_solver
    !this%kernel_smoother => kernel_smoother
    this%particle_group => particle_group
!    this%efield_dofs => efield_dofs

!!$    SLL_ALLOCATE(this%rho_dofs(this%kernel_smoother%n_dofs), ierr)
!!$    SLL_ALLOCATE(this%rho_dofs_local(this%kernel_smoother%n_dofs), ierr)
!!$    SLL_ALLOCATE(this%rho_2d(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)
!!$    SLL_ALLOCATE(this%efield1(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)   
!!$    SLL_ALLOCATE(this%efield2(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)

    this%i_weight = 1
    if(present(i_weight)) this%i_weight = i_weight
    if(present(control_variate)) this%control_variate => control_variate

  end function sll_new_hamiltonian_splitting_pic_vp_2d2v




end module sll_m_operator_splitting_pic_vp_2d2v
