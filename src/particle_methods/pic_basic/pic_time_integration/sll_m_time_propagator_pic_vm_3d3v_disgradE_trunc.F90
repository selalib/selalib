!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on antisymmetric splitting with disgradE for 3d3v Vlasov-Maxwell, modified by approximating particle mass by mass matrix
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_time_propagator_pic_vm_3d3v_disgradE_trunc
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
  
  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use sll_m_time_propagator_pic_vm_3d3v_disgradE, only : &
       sll_t_time_propagator_pic_vm_3d3v_disgradE

  use sll_m_time_propagator_pic_vm_3d3v_helper, only: &
       sll_t_time_propagator_pic_vm_3d3v_helper
  
  use sll_m_particle_mesh_coupling_spline_3d_feec, only: &
      sll_t_particle_mesh_coupling_spline_3d_feec
 
  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base, &
       sll_t_particle_array

  use sll_m_filter_base_3d, only: &
       sll_c_filter_base_3d
  
  
  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_3d3v_disgradE_trunc

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Time propagator for Vlasov-Maxwell 3d3v
  type, extends(sll_t_time_propagator_pic_vm_3d3v_disgradE) :: sll_t_time_propagator_pic_vm_3d3v_disgradE_trunc

   contains

     procedure :: lie_splitting => lie_splitting_pic_vm_3d3v_trunc !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_3d3v_trunc !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_3d3v_trunc !> Strang splitting propagator

     
  end type sll_t_time_propagator_pic_vm_3d3v_disgradE_trunc

contains

  
  !> Strang splitting
  subroutine strang_splitting_pic_vm_3d3v_trunc(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_disgradE_trunc), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    do i_step = 1, number_steps
       call self%helper%advect_x( dt*0.5_f64 )
       call self%helper%advect_vb( dt*0.5_f64 )
       call self%helper%advect_eb( dt*0.5_f64 )
       call self%helper%advect_e_trunc( dt )
       call self%helper%advect_eb( dt*0.5_f64 )
       call self%helper%advect_vb( dt*0.5_f64 )
       call self%helper%advect_x( dt*0.5_f64 )
       
    end do

  end subroutine strang_splitting_pic_vm_3d3v_trunc

  
  !> Lie splitting
  subroutine lie_splitting_pic_vm_3d3v_trunc(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_disgradE_trunc), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    do i_step = 1,number_steps
       call self%helper%advect_x( dt )
       call self%helper%advect_vb( dt )
       call self%helper%advect_eb( dt )
       call self%helper%advect_e_trunc( dt )
       
    end do


  end subroutine lie_splitting_pic_vm_3d3v_trunc

  
  !> Lie splitting (oposite ordering)
  subroutine lie_splitting_back_pic_vm_3d3v_trunc(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_disgradE_trunc), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    do i_step = 1,number_steps
       call self%helper%advect_e_trunc( dt )
       call self%helper%advect_eb( dt )
       call self%helper%advect_vb( dt )
       call self%helper%advect_x( dt )
    end do

  end subroutine lie_splitting_back_pic_vm_3d3v_trunc

  
end module sll_m_time_propagator_pic_vm_3d3v_disgradE_trunc
