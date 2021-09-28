!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on energy and charge-conserving discrete gradient method, implicit
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
!> Reference: Kormann, Sonnendrücker, Energy-conserving time propagation for a structure-preserving particle-in-cell Vlasov–Maxwell solver, Journal of Computational Physics 425, 109890, 2021
module sll_m_time_propagator_pic_vm_1d2v_disgradEC
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_filter_base_1d, only: &
       sll_c_filter_base_1d

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
       sll_t_time_propagator_pic_vm_1d2v_disgradEC

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_1d2v_disgradEC

     type(sll_t_time_propagator_pic_vm_1d2v_helper) :: helper !< helper
     logical    :: electrostatic = .false. !< true for electrostatic simulation
   contains

     procedure :: lie_splitting => lie_splitting_pic_vm_1d2v_disgradEC !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_1d2v_disgradEC !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_1d2v_disgradEC !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_1d2v_disgradEC !> Initialize the type
     procedure :: init_from_file => initialize_file_pic_vm_1d2v_disgradEC !> Initialize the type
     procedure :: free => delete_pic_vm_1d2v_disgradEC !> Finalization

     procedure :: reinit_fields => reinit_fields_disgradEC

  end type sll_t_time_propagator_pic_vm_1d2v_disgradEC


contains

  !> Strang splitting
  subroutine strang_splitting_pic_vm_1d2v_disgradEC(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variable
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_ex( dt )
          call self%helper%advect_vb( dt*0.5_f64 )
       end do
    else
       do i_step = 1, number_steps
          call self%helper%advect_eb( dt*0.5_f64 )
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_ex( dt )
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_eb( dt*0.5_f64 )
       end do
    end if


  end subroutine strang_splitting_pic_vm_1d2v_disgradEC

  !> Lie splitting
  subroutine lie_splitting_pic_vm_1d2v_disgradEC(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variable
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1,number_steps
          call self%helper%advect_vb( dt )
          call self%helper%advect_ex( dt )
       end do
    else
       do i_step = 1,number_steps
          call self%helper%advect_eb( dt )
          call self%helper%advect_vb( dt )
          call self%helper%advect_ex( dt )
       end do
    end if


  end subroutine lie_splitting_pic_vm_1d2v_disgradEC


  !> Lie splitting
  subroutine lie_splitting_back_pic_vm_1d2v_disgradEC(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variable
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1,number_steps
          call self%helper%advect_ex( dt )
          call self%helper%advect_vb( dt )
       end do
    else
       do i_step = 1,number_steps
          call self%helper%advect_ex( dt )
          call self%helper%advect_vb( dt )
          call self%helper%advect_eb( dt )
       end do
    end if


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
       boundary_particles, &
       solver_tolerance, &
       iter_tolerance, max_iter, &
       force_sign, &
       betar, &
       electrostatic, &
       jmean)  
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC), intent(out) :: self !< time splitting object 
    class(sll_c_maxwell_1d_base), target,          intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,      intent(in)  :: particle_group !< Particle group
    sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
    class(sll_c_filter_base_1d), target,            intent(in) :: filter
    sll_int32, optional,                                  intent( in )  :: boundary_particles !< particle boundary conditions
    sll_real64, optional,                                 intent( in )  :: solver_tolerance !< solver tolerance
    sll_real64, optional,                          intent( in ) :: iter_tolerance !< iteration tolerance
    sll_int32,  optional,                          intent( in ) :: max_iter !< maximal number of iterations
    sll_real64, optional,                                 intent( in )  :: force_sign !< sign of particle force
    sll_real64, optional, intent(in) :: betar(2) !< reciprocal plasma beta
    logical, optional,                                    intent( in )  :: electrostatic !< true for electrostatic simulation
    logical, optional,                                    intent( in )  :: jmean !< logical for mean value of current
    !local variables
    sll_int32 :: max_iter_set
    sll_int32 :: boundary_particles_set
    sll_real64 :: solver_tolerance_set, iter_tolerance_set
    sll_real64 :: force_sign_set, betar_set(2)
    logical :: jmean_set

    if (present(boundary_particles)) then
       boundary_particles_set = boundary_particles
    else
       boundary_particles_set = 100
    end if

    if (present(solver_tolerance) )  then
       solver_tolerance_set = solver_tolerance
    else
       solver_tolerance_set = 1d-12
    end if

    if (present(iter_tolerance) )  then
       iter_tolerance_set = iter_tolerance
    else
       iter_tolerance_set = 1d-12
    end if

    if (present(max_iter) )  then
       max_iter_set = max_iter
    else
       max_iter_set = 10
    end if

    if( present(force_sign) )then
       force_sign_set = force_sign
    else
       force_sign_set = 1._f64
    end if

    if( present(betar) )then
       betar_set = betar
    else
       betar_set = 1._f64
    end if

    if( present(electrostatic) )then
       self%electrostatic = electrostatic
    end if

    if (present(jmean)) then
       jmean_set = jmean
    else
       jmean_set = .false.
    end if

    call self%helper%init( maxwell_solver, &
         kernel_smoother_0, &
         kernel_smoother_1, &
         particle_group, &
         efield_dofs, &
         bfield_dofs, &
         x_min, &
         Lx, &
         filter,&
         .false., &
         boundary_particles, &
         solver_tolerance_set, &
         betar=betar_set, &
         force_sign=force_sign_set,&
         jmean=jmean_set)

  end subroutine initialize_pic_vm_1d2v_disgradEC


  !> Constructor.
  subroutine initialize_file_pic_vm_1d2v_disgradEC(&
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
       filename, &
       boundary_particles, &
       force_sign, &
       betar, &
       electrostatic, &
       jmean)  
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_1d_base), target,                 intent( in )  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent( in )  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent( in )  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,                  intent( in )  :: particle_group !< Particle group
    sll_real64, target,                                   intent( in )  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                                   intent( in )  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                           intent( in )  :: x_min !< Lower bound of x domain
    sll_real64,                                           intent( in )  :: Lx !< Length of the domain in x direction.
    class(sll_c_filter_base_1d), target,            intent(in) :: filter
    character(len=*),                                     intent( in )  :: filename !< filename
    sll_int32, optional, intent( in )  :: boundary_particles !< particle boundary conditions
    sll_real64, optional, intent(in) :: force_sign !< sign of particle force
    sll_real64, optional,                          intent( in ) :: betar(2) !< reciprocal plasma beta
    logical, optional, intent(in)    :: electrostatic !< true for electrostatic simulation
    logical, optional, intent(in)    :: jmean !< logical for mean value of current
    !local variables
    sll_int32 :: boundary_particles_set
    sll_real64 :: force_sign_set, betar_set(2)
    logical :: jmean_set


    if (present(boundary_particles)) then
       boundary_particles_set = boundary_particles
    else
       boundary_particles_set = 100
    end if

    if( present(force_sign) )then
       force_sign_set = force_sign
    else
       force_sign_set = 1._f64
    end if

    if( present(betar) )then
       betar_set = betar
    else
       betar_set = 1._f64
    end if

    if( present(electrostatic) )then
       self%electrostatic = electrostatic
    end if


    if (present(jmean)) then
       jmean_set = jmean
    else
       jmean_set = .false.
    end if

    call self%helper%init_from_file( maxwell_solver, &
         kernel_smoother_0, &
         kernel_smoother_1, &
         particle_group, &
         efield_dofs, &
         bfield_dofs, &
         x_min, &
         Lx, &
         filter,&
         filename, &
         .false., &
         boundary_particles = boundary_particles_set, &
         betar=betar_set, &
         force_sign=force_sign_set, &
         jmean=jmean_set)


  end subroutine initialize_file_pic_vm_1d2v_disgradEC


  subroutine delete_pic_vm_1d2v_disgradEC(self)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC), intent( inout ) :: self !< time splitting object 

    call self%helper%free()

  end subroutine delete_pic_vm_1d2v_disgradEC


  subroutine reinit_fields_disgradEC( self ) 
    class(sll_t_time_propagator_pic_vm_1d2v_disgradEC), intent(inout) :: self !< time splitting object 

    call self%helper%reinit_fields()

  end subroutine reinit_fields_disgradEC

end module sll_m_time_propagator_pic_vm_1d2v_disgradEC
