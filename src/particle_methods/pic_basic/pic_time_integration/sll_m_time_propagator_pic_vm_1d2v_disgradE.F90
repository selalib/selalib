!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on energy (not charge-conserving) discrete gradient method, semi-implicit
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
!> Reference: Kormann, Sonnendrücker, Energy-conserving time propagation for a structure-preserving particle-in-cell Vlasov–Maxwell solver, Journal of Computational Physics 425, 109890, 2021
module sll_m_time_propagator_pic_vm_1d2v_disgradE
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_control_variate, only: &
       sll_t_control_variates

  use sll_m_filter_base_1d, only: &
       sll_c_filter_base_1d

  use sll_m_maxwell_1d_base, only: &
       sll_c_maxwell_1d_base

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use sll_m_particle_mesh_coupling_base_1d, only: &
       sll_c_particle_mesh_coupling_1d

  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use  sll_m_time_propagator_pic_vm_1d2v_helper, only : &
       sll_t_time_propagator_pic_vm_1d2v_helper


  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_1d2v_disgradE

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_1d2v_disgradE
     type(sll_t_time_propagator_pic_vm_1d2v_helper) :: helper !< helper
     logical    :: electrostatic = .false. !< true for electrostatic simulation

   contains
     procedure :: lie_splitting => lie_splitting_pic_vm_1d2v_disgradE !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_1d2v_disgradE !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_1d2v_disgradE !> Strang splitting propagator

     procedure :: reinit_fields

     procedure :: init => initialize_pic_vm_1d2v_disgradE !> Initialize the type
     procedure :: init_from_file => initialize_file_pic_vm_1d2v_disgradE !> Initialize the type
     procedure :: free => delete_pic_vm_1d2v_disgradE !> Finalization

  end type sll_t_time_propagator_pic_vm_1d2v_disgradE

contains
  subroutine reinit_fields( self ) 
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE), intent(inout) :: self !< time propagator object 

    call self%helper%reinit_fields()

  end subroutine reinit_fields


  !> Strang splitting
  subroutine strang_splitting_pic_vm_1d2v_disgradE(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%helper%advect_x( dt*0.5_f64 )
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_e( dt )
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_x( dt*0.5_f64 )
       end do
    else
       do i_step = 1, number_steps
          call self%helper%advect_eb( dt*0.5_f64 )
          call self%helper%advect_x( dt*0.5_f64 )
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_e( dt )
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_x( dt*0.5_f64 )
          call self%helper%advect_eb( dt*0.5_f64 )
       end do
    end if

  end subroutine strang_splitting_pic_vm_1d2v_disgradE


  !> Lie splitting
  subroutine lie_splitting_pic_vm_1d2v_disgradE(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1,number_steps
          call self%helper%advect_e( dt )
          call self%helper%advect_vb( dt )
          call self%helper%advect_x( dt )
       end do
    else
       do i_step = 1,number_steps
          call self%helper%advect_e( dt )
          call self%helper%advect_vb( dt )
          call self%helper%advect_x( dt )
          call self%helper%advect_eb( dt )
       end do
    end if

  end subroutine lie_splitting_pic_vm_1d2v_disgradE


  !> Lie splitting
  subroutine lie_splitting_back_pic_vm_1d2v_disgradE(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1,number_steps
          call self%helper%advect_e( dt )
          call self%helper%advect_vb( dt )
          call self%helper%advect_x( dt )
       end do
    else
       do i_step = 1,number_steps
          call self%helper%advect_e( dt )
          call self%helper%advect_vb( dt )
          call self%helper%advect_x( dt )
          call self%helper%advect_eb( dt )
       end do
    end if

  end subroutine lie_splitting_back_pic_vm_1d2v_disgradE


  !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_1d2v_disgradE(&
       self, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       phi_dofs, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       filter, &
       boundary_particles, &
       solver_tolerance, &
       force_sign, &
       control_variate, &
       i_weight, &
       betar, &
       electrostatic, &
       jmean) 
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE), intent(out) :: self !< time propagator object 
    class(sll_c_maxwell_1d_base), target,          intent(in)  :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,      intent(in)  :: particle_group !< Particle group
    sll_real64, target,                            intent(in)  :: phi_dofs(:) !< array for the coefficients of phi
    sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
    class( sll_c_filter_base_1d ), intent( in ), target :: filter !< filter
    sll_int32, optional,                                  intent( in )  :: boundary_particles !< particle boundary conditions
    sll_real64, intent(in), optional :: solver_tolerance !< solver tolerance
    sll_real64, optional, intent(in) :: force_sign !< sign of particle force
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    sll_int32, optional,                            intent(in) :: i_weight !< Index of weight to be used by propagator
    sll_real64, optional, intent(in) :: betar(2) !< reciprocal plasma beta
    logical, optional    :: electrostatic !< true for electrostatic simulation
    logical, optional, intent(in) :: jmean !< logical for mean value of current
    !local variables
    sll_int32 :: boundary_particles_set 
    sll_real64 :: solver_tolerance_set, force_sign_set, betar_set(2)
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

    if( present( control_variate) ) then
       call self%helper%init(maxwell_solver, &
            kernel_smoother_0, &
            kernel_smoother_1, &
            particle_group, &
            phi_dofs, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx, &
            filter, &
            boundary_particles=boundary_particles_set, &
            solver_tolerance=solver_tolerance_set, &
            force_sign=force_sign_set, &
            control_variate = control_variate,&
            i_weight=i_weight, &
            betar=betar_set, &
            jmean=jmean_set)
    else
       call self%helper%init(maxwell_solver, &
            kernel_smoother_0, &
            kernel_smoother_1, &
            particle_group, &
            phi_dofs, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx, &
            filter, &
            boundary_particles=boundary_particles_set, &
            solver_tolerance=solver_tolerance_set, &
            betar=betar_set, &
            force_sign=force_sign_set, &
            jmean=jmean_set)
    end if


  end subroutine initialize_pic_vm_1d2v_disgradE


  !> Constructor.
  subroutine initialize_file_pic_vm_1d2v_disgradE(&
       self, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       phi_dofs, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       filter, &
       filename, &
       boundary_particles, &
       force_sign, &
       control_variate, &
       i_weight, &
       betar, &
       electrostatic, &
       jmean) 
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE), intent(out) :: self !< time propagator object 
    class(sll_c_maxwell_1d_base), target,          intent(in)  :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,      intent(in)  :: particle_group !< Particle group
    sll_real64, target,                            intent(in)  :: phi_dofs(:) !< array for the coefficients of phi
    sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
    class( sll_c_filter_base_1d ), intent( in ), target :: filter !< filter
    character(len=*), intent(in) :: filename !< filename
    sll_int32, optional,  intent( in )  :: boundary_particles !< particle boundary conditions
    sll_real64, optional, intent(in) :: force_sign !< sign of particle force
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    sll_int32, optional,                            intent(in) :: i_weight !< Index of weight to be used by propagator
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
    else
       self%electrostatic = .false.
    end if

    if (present(jmean)) then
       jmean_set = jmean
    else
       jmean_set = .false.
    end if

    if( present( control_variate) ) then
       call self%helper%init_from_file(  maxwell_solver, &
            kernel_smoother_0, &
            kernel_smoother_1, &
            particle_group, &
            phi_dofs, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx, &
            filter, &
            filename, &
            boundary_particles = boundary_particles_set, &
            force_sign=force_sign_set, &
            control_variate = control_variate, &
            i_weight=i_weight, &
            betar=betar_set, &
            jmean=jmean_set)
    else
       call self%helper%init_from_file(  maxwell_solver, &
            kernel_smoother_0, &
            kernel_smoother_1, &
            particle_group, &
            phi_dofs, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx, &
            filter, &
            filename, &
            boundary_particles = boundary_particles_set, &
            betar=betar_set, &
            force_sign=force_sign_set, &
            jmean=jmean_set)
    end if


  end subroutine initialize_file_pic_vm_1d2v_disgradE


  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_1d2v_disgradE(self)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE), intent( inout ) :: self !< time propagator object
    call self%helper%free() 

  end subroutine delete_pic_vm_1d2v_disgradE


end module sll_m_time_propagator_pic_vm_1d2v_disgradE
