!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on operator splitting for 2d2v Vlasov-Poisson.
!> @details MPI parallelization by domain cloning. Periodic boundaries.
module sll_m_operator_splitting_pic_vp_2d2v
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_o_collective_allreduce, &
    sll_v_world_collective

  use sll_m_control_variate, only: &
    sll_t_control_variate

  use sll_m_kernel_smoother_base, only: &
    sll_c_kernel_smoother

  use sll_m_operator_splitting, only: &
    sll_t_operator_splitting

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_pic_poisson_base, only: &
    sll_c_pic_poisson

  use sll_mpi, only: &
    mpi_sum

  implicit none

  public :: &
    sll_s_new_operator_splitting_pic_vp_2d2v, &
    sll_t_operator_splitting_pic_vp_2d2v

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Operator splitting type for 2d2v Vlasov-Poisson
  type, extends(sll_t_operator_splitting) :: sll_t_operator_splitting_pic_vp_2d2v
     class(sll_c_pic_poisson), pointer :: solver
     class(sll_c_particle_group_base), pointer  :: particle_group    !< Particle group


     ! For version with control variate
     class(sll_t_control_variate), pointer :: control_variate => null()
     sll_int32 :: i_weight

   contains
     procedure :: operatorT => advection_x_pic_vp_2d2v  !< Operator for x advection
     procedure :: operatorV => advection_v_pic_vp_2d2v  !< Operator for v advection
     procedure :: strang_splitting => strang_splitting_pic_vp_2d2v !< Strang splitting

     procedure :: init => initialize_operator_splitting_pic_vp_2d2v !>initializer
  end type sll_t_operator_splitting_pic_vp_2d2v

contains

  !---------------------------------------------------------------------------!
  !> Strang splitting
  subroutine strang_splitting_pic_vp_2d2v(this, dt)
    class(sll_t_operator_splitting_pic_vp_2d2v), intent(inout) :: this !< time splitting object 
    sll_real64,                                  intent(in)    :: dt   !< time step


    call this%operatorT(0.5_f64*dt)
    call this%operatorV(dt)
    call this%operatorT(0.5_f64*dt)


  end subroutine strang_splitting_pic_vp_2d2v
  
  !---------------------------------------------------------------------------!


  !---------------------------------------------------------------------------!
  !> Push x 
  subroutine advection_x_pic_vp_2d2v(this, dt)
    class(sll_t_operator_splitting_pic_vp_2d2v), intent(inout) :: this !< time splitting object 
    sll_real64,                                  intent(in)    :: dt   !< time step


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
          wi(3) = this%control_variate%update_df_weight( x_new(1:2), vi(1:2), 0.0_f64, wi(1), wi(2))
          call this%particle_group%set_weights(i_part, wi)
       end if
    end do


  end subroutine advection_x_pic_vp_2d2v
  
  !---------------------------------------------------------------------------!
  ! Push v
  subroutine advection_v_pic_vp_2d2v(this, dt)
    class(sll_t_operator_splitting_pic_vp_2d2v), intent(inout) :: this !< time splitting object 
    sll_real64,                                  intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: v_new(3)
    sll_real64 :: efield(2)
    sll_real64 :: qm
    sll_real64 :: xi(3)
    sll_real64 :: wi, wall(3)


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
       if (this%particle_group%n_weights == 3) then
          ! Update weights if control variate
          wall = this%particle_group%get_weights(i_part)          
          wall(3) = this%control_variate%update_df_weight( xi(1:2), v_new(1:2), 0.0_f64, wall(1), wall(2))
          call this%particle_group%set_weights(i_part, wall)
       end if
    end do
    

  end subroutine advection_v_pic_vp_2d2v
  
  !---------------------------------------------------------------------------!
  !> Initialization function
  subroutine initialize_operator_splitting_pic_vp_2d2v(&
       this, &
       solver, &
       particle_group, &
       control_variate, &
       i_weight)
    class(sll_t_operator_splitting_pic_vp_2d2v),    intent(out):: this !< time splitting object 
    class(sll_c_pic_poisson),pointer,               intent(in) :: solver !< Poisson solver
    class(sll_c_particle_group_base),pointer,       intent(in) :: particle_group !< Particle group
    class(sll_t_control_variate), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    sll_int32, optional,                            intent(in) :: i_weight !< Index of weight to be used by propagator

    this%solver => solver
    this%particle_group => particle_group

    this%i_weight = 1
    if(present(i_weight)) this%i_weight = i_weight
    if(present(control_variate)) this%control_variate => control_variate

  end subroutine initialize_operator_splitting_pic_vp_2d2v

  !---------------------------------------------------------------------------!
  !> Destructor
  subroutine delete_operator_splitting_pic_vp_2d2v( this )
    class(sll_t_operator_splitting_pic_vp_2d2v), intent(inout) :: this !< time splitting object 

    this%solver => null()
    this%particle_group => null()
    this%control_variate => null()

  end subroutine delete_operator_splitting_pic_vp_2d2v

  !---------------------------------------------------------------------------!
  !> Constructor for abstract type
  subroutine sll_s_new_operator_splitting_pic_vp_2d2v &
       (splitting, &
       solver, &
       particle_group, &
       control_variate, &
       i_weight) 
    class(sll_t_operator_splitting), pointer,        intent(out):: splitting !< time splitting object 
    class(sll_c_pic_poisson),pointer,                intent(in) :: solver !< Poisson solver
    class(sll_c_particle_group_base),pointer,        intent(in) :: particle_group !< Particle group
    class(sll_t_control_variate), optional, pointer, intent(in) :: control_variate !< Control variate (if delta f)
    sll_int32, optional,                             intent(in) :: i_weight !< Index of weight to be used by propagator
 
    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE( sll_t_operator_splitting_pic_vp_2d2v :: splitting , ierr)

    select type( splitting )
    type is ( sll_t_operator_splitting_pic_vp_2d2v )
       call splitting%init(solver, particle_group, control_variate, i_weight)
    end select

  end subroutine sll_s_new_operator_splitting_pic_vp_2d2v



end module sll_m_operator_splitting_pic_vp_2d2v
