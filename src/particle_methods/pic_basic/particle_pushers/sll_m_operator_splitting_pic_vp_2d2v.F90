!> @ingroup particle_pushers
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on operator splitting for 2d2v Vlasov-Poisson.
!> @details MPI parallelization by domain cloning. Periodic boundaries.
module sll_m_operator_splitting_pic_vp_2d2v
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_collective_allreduce, &
    sll_world_collective

  use sll_m_kernel_smoother_base, only: &
    sll_kernel_smoother_base

  use sll_m_operator_splitting, only: &
    operator_splitting

  use sll_m_particle_group_base, only: &
    sll_particle_group_base

  use sll_m_poisson_2d_fft, only: &
    poisson_2d_fft_solver

  use sll_mpi, only: &
    mpi_sum

  implicit none

  public :: &
    sll_new_hamiltonian_splitting_pic_vp_2d2v, &
    sll_t_operator_splitting_pic_vp_2d2v

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Operator splitting type for 2d2v Vlasov-Poisson
  type, extends(operator_splitting) :: sll_t_operator_splitting_pic_vp_2d2v
     class(poisson_2d_fft_solver), pointer    :: poisson_solver      !< Poisson solver (TODO: Use a base class here)
     class(sll_kernel_smoother_base), pointer :: kernel_smoother  !< Kernel smoother
     class(sll_particle_group_base), pointer  :: particle_group    !< Particle group

     sll_real64, allocatable :: rho_dofs(:)                      !< Degrees of freedom for kernel respresentation of the charge density. 
     sll_real64, allocatable :: rho_dofs_local(:)                !< Processor local part of \a rho_dofs 
     sll_real64, allocatable :: rho_2d(:,:)                      !< 2d representation of \a rho_dofs for use in Poisson solver.
     sll_real64, allocatable :: efield1(:,:) !< 2d representation of electric field, first component for use in Poisson solver. 
     sll_real64, allocatable :: efield2(:,:) !< 2d representation of electric field, second component for use in Poisson solver. 
     sll_real64, pointer     :: efield_dofs(:,:)  !< Values of the electric field at grid points (1d representation).

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
    sll_real64 :: x_new(3)

    ! X_new = X_old + dt * V
    do i_part=1,this%particle_group%n_particles
       x_new = this%particle_group%get_x(i_part) + dt * this%particle_group%get_v(i_part) 
       call this%particle_group%set_x(i_part, x_new)
    end do
    
    ! Update shape factors since position is changed
    call this%kernel_smoother%compute_shape_factors(this%particle_group)

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


    ! Assemble right-hand-side
    !call this%kernel_smoother%compute_shape_factors(this%particle_group)
    this%rho_dofs_local = 0.0_f64
    call this%kernel_smoother%accumulate_rho_from_klimontovich(this%particle_group, &
         this%rho_dofs_local)

    ! MPI to sum up contributions from each processor
    this%rho_dofs = 0.0_f64
    call sll_collective_allreduce( sll_world_collective, this%rho_dofs_local, &
         this%kernel_smoother%n_dofs, MPI_SUM, this%rho_dofs)

    ! Bring to right shape for Poisson solver
    this%rho_2d(1:this%kernel_smoother%n_grid(1),1:this%kernel_smoother%n_grid(2)) = &
         -1.0_f64 +&
         reshape(this%rho_dofs, this%kernel_smoother%n_grid(1:2) )
    this%rho_2d(this%kernel_smoother%n_grid(1)+1,:) = this%rho_2d(1,:)
    this%rho_2d(:,this%kernel_smoother%n_grid(2)+1) = this%rho_2d(:,1)
    ! Solve Poisson problem
    call this%poisson_solver%compute_E_from_rho(this%efield1, this%efield2, this%rho_2d)

    ! Evaluate efield1 at particle positions
    this%efield_dofs(:,1) = reshape(this%efield1(1:this%kernel_smoother%n_grid(1),1:this%kernel_smoother%n_grid(2)), [this%kernel_smoother%n_dofs])
    !call this%kernel_smoother%evaluate_kernel_function_particles(this%particle_group, &
    !     this%rho_dofs, this%efield(:,1))
    
    ! Evaluate efield2 at particle positions
    this%efield_dofs(:,2) = reshape(this%efield2(1:this%kernel_smoother%n_grid(1),1:this%kernel_smoother%n_grid(2)), [this%kernel_smoother%n_dofs])
    !call this%kernel_smoother%evaluate_kernel_function_particles(this%particle_group, &
    !     this%rho_dofs, this%efield(:,2))


    ! V_new = V_old + dt *q/m* E
    qm = this%particle_group%species%q_over_m();
    do i_part=1,this%particle_group%n_particles
       ! Evaluate efields at particle position
       call this%kernel_smoother%evaluate_kernel_function_particle&
            (this%efield_dofs(:,1), i_part, efield(1))
       call this%kernel_smoother%evaluate_kernel_function_particle&
            (this%efield_dofs(:,2), i_part, efield(2))
       v_new = this%particle_group%get_v(i_part)
       v_new(1:2) = v_new(1:2) + dt * qm * efield
       call this%particle_group%set_v(i_part, v_new)
    end do
    

  end subroutine advection_v_pic_vp_2d2v
  
  !---------------------------------------------------------------------------!
  !> Initialization function
  subroutine initialize_operator_splitting_pic_vp_2d2v(this, poisson_solver, kernel_smoother, particle_group)
    class(sll_t_operator_splitting_pic_vp_2d2v), intent(out) :: this !< object 
    class(poisson_2d_fft_solver),pointer, intent(in) :: poisson_solver !< poisson solver
    class(sll_kernel_smoother_base),pointer, intent(in) :: kernel_smoother !< kernel smoother
    class(sll_particle_group_base),pointer, intent(in) :: particle_group !< particle group

    !local variables
    sll_int32 :: ierr

    this%poisson_solver => poisson_solver
    this%kernel_smoother => kernel_smoother
    this%particle_group => particle_group

    SLL_ALLOCATE(this%rho_dofs(this%kernel_smoother%n_dofs), ierr)
    SLL_ALLOCATE(this%rho_2d(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)
    SLL_ALLOCATE(this%efield1(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)   
    SLL_ALLOCATE(this%efield2(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)
    SLL_ALLOCATE(this%efield_dofs(this%kernel_smoother%n_dofs, 2), ierr)

  end subroutine initialize_operator_splitting_pic_vp_2d2v


  !---------------------------------------------------------------------------!
  !> Constructor.
  function sll_new_hamiltonian_splitting_pic_vp_2d2v(poisson_solver, kernel_smoother, particle_group, efield_dofs) result(this)
    class(sll_t_operator_splitting_pic_vp_2d2v), pointer :: this !< time splitting object 
    class(poisson_2d_fft_solver),pointer, intent(in) :: poisson_solver !< Poisson solver
    class(sll_kernel_smoother_base),pointer, intent(in) :: kernel_smoother !< Kernel smoother
    class(sll_particle_group_base),pointer, intent(in) :: particle_group !< Particle group
    sll_real64, pointer  :: efield_dofs(:,:)  !< Values of the electric field at grid points
 
    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE(this, ierr)

    this%poisson_solver => poisson_solver
    this%kernel_smoother => kernel_smoother
    this%particle_group => particle_group
    this%efield_dofs => efield_dofs

    SLL_ALLOCATE(this%rho_dofs(this%kernel_smoother%n_dofs), ierr)
    SLL_ALLOCATE(this%rho_dofs_local(this%kernel_smoother%n_dofs), ierr)
    SLL_ALLOCATE(this%rho_2d(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)
    SLL_ALLOCATE(this%efield1(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)   
    SLL_ALLOCATE(this%efield2(this%kernel_smoother%n_grid(1)+1, this%kernel_smoother%n_grid(2)+1), ierr)

  end function sll_new_hamiltonian_splitting_pic_vp_2d2v




end module sll_m_operator_splitting_pic_vp_2d2v
