!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on a variation of the energy and charge-conserving discrete gradient method with subcycling, implicit
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
!> Reference: Kormann, Sonnendrücker, Energy-conserving time propagation for a structure-preserving particle-in-cell Vlasov–Maxwell solver, Journal of Computational Physics 425, 109890, 2021
module sll_m_time_propagator_pic_vm_1d2v_disgradEC_sub
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_binomial_filter, only: &
       sll_t_binomial_filter

  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use  sll_m_time_propagator_pic_vm_1d2v_helper, only : &
       sll_t_time_propagator_pic_vm_1d2v_helper

  use sll_m_particle_mesh_coupling_base_1d, only: &
       sll_c_particle_mesh_coupling_1d

  use sll_m_maxwell_1d_base, only: &
       sll_c_maxwell_1d_base

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub

     type(sll_t_time_propagator_pic_vm_1d2v_helper) :: helper !< helper

   contains

     procedure :: lie_splitting => lie_splitting_pic_vm_1d2v_disgradEC !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_1d2v_disgradEC !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_1d2v_disgradEC !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_1d2v_disgradEC !> Initialize the type
     procedure :: free => delete_pic_vm_1d2v_disgradEC !> Finalization

     procedure :: reinit_fields => reinit_fields_disgradEC_sub

  end type sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub


contains


  !> Strang splitting
  subroutine strang_splitting_pic_vm_1d2v_disgradEC(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1, number_steps
       call self%helper%advect_vb(dt*0.5_f64)
       call self%helper%advect_eb(dt*0.5_f64)

       call self%helper%advect_e_sub(dt)
       call self%helper%advect_eb(dt*0.5_f64)
       call self%helper%advect_vb(dt*0.5_f64)
    end do

  end subroutine strang_splitting_pic_vm_1d2v_disgradEC


  !> Lie splitting
  subroutine lie_splitting_pic_vm_1d2v_disgradEC(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step
    sll_real64 :: transfer

    do i_step = 1, number_steps
       call self%helper%advect_eb(dt)
       call self%helper%advect_vb(dt)
       call self%helper%advect_e_sub(dt)
    end do

  end subroutine lie_splitting_pic_vm_1d2v_disgradEC


  !> Lie splitting
  subroutine lie_splitting_back_pic_vm_1d2v_disgradEC(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step
    sll_real64 :: transfer

    do i_step = 1, number_steps
       call self%helper%advect_e_sub(dt)
       call self%helper%advect_vb(dt)
       call self%helper%advect_eb(dt)
    end do

  end subroutine lie_splitting_back_pic_vm_1d2v_disgradEC


  !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_1d2v_disgradEC(&
       self, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       filter, &
       filename) 
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub), intent(out) :: self !< time splitting object 
    class(sll_c_maxwell_1d_base), target,          intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,      intent(in)  :: particle_group !< Particle group
    sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
    type(sll_t_binomial_filter), target,            intent(in) :: filter
    character(len=*), intent(in) :: filename !< name of nml-file


    call self%helper%init( maxwell_solver, &
         kernel_smoother_0, &
         kernel_smoother_1, &
         particle_group, &
         efield_dofs, &
         bfield_dofs, &
         x_min, &
         Lx, &
         filter, &
         filename) 

  end subroutine initialize_pic_vm_1d2v_disgradEC


  subroutine delete_pic_vm_1d2v_disgradEC(self)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub), intent( inout ) :: self !< time splitting object 

    call self%helper%free()

  end subroutine delete_pic_vm_1d2v_disgradEC


  subroutine reinit_fields_disgradEC_sub( self ) 
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub), intent(inout) :: self !< time splitting object 

    call self%helper%reinit_fields()

  end subroutine reinit_fields_disgradEC_sub


end module sll_m_time_propagator_pic_vm_1d2v_disgradEC_sub
