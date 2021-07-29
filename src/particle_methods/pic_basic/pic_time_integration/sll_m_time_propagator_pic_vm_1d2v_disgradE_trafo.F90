!> @ingroup pic_time_integration
!> @author Benedikt Perse, IPP
!> @brief Particle pusher based on antisymmetric splitting with AVF for 1d2v Vlasov-Poisson with coordinate transformation.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_time_propagator_pic_vm_1d2v_disgradE_trafo
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
       sll_f_get_collective_rank, &
       sll_o_collective_allreduce, &
       sll_v_world_collective

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_maxwell_1d_base, only: &
       sll_c_maxwell_1d_base

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use sll_m_particle_mesh_coupling_base_1d, only: &
       sll_c_particle_mesh_coupling_1d

  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use sll_m_time_propagator_pic_vm_1d2v_trafo_helper, only:&
       sll_t_time_propagator_pic_vm_1d2v_trafo_helper



  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo
     type(sll_t_time_propagator_pic_vm_1d2v_trafo_helper) :: helper !< helper
     logical    :: electrostatic = .false. !< true for electrostatic simulation


   contains
    
     procedure :: lie_splitting => lie_splitting_pic_vm_1d2v !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_1d2v !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_1d2v !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_1d2v !> Initialize the type
     procedure :: init_from_file => initialize_file_pic_vm_1d2v !> Initialize from nml file
     procedure :: free => delete_pic_vm_1d2v !> Finalization

  end type sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo

contains


  !> Strang splitting
  subroutine strang_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%helper%advect_x(dt*0.5_f64)
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_e(dt)
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_x(dt*0.5_f64)
       end do
    else
       do i_step = 1, number_steps
          call self%helper%advect_eb(dt*0.5_f64)
          call self%helper%advect_x(dt*0.5_f64)
          call self%helper%advect_vb(dt*0.5_f64)
          call self%helper%advect_e(dt)
          call self%helper%advect_vb(dt*0.5_f64)
          call self%helper%advect_x(dt*0.5_f64)
          call self%helper%advect_eb(dt*0.5_f64)
       end do
    end if

  end subroutine strang_splitting_pic_vm_1d2v


  !> Lie splitting
  subroutine lie_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%helper%advect_x(dt)
          call self%helper%advect_vb( dt )
          call self%helper%advect_e(dt)
       end do
    else
       do i_step = 1, number_steps
          call self%helper%advect_x(dt)
          call self%helper%advect_eb(dt)
          call self%helper%advect_vb(dt)
          call self%helper%advect_e(dt)
       end do
    end if

  end subroutine lie_splitting_pic_vm_1d2v


  !> Lie splitting
  subroutine lie_splitting_back_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%helper%advect_e(dt)
          call self%helper%advect_vb( dt )
          call self%helper%advect_x(dt)
       end do
    else
       do i_step = 1, number_steps
          call self%helper%advect_e(dt)
          call self%helper%advect_vb(dt)
          call self%helper%advect_eb(dt)
          call self%helper%advect_x(dt)
       end do
    end if

  end subroutine lie_splitting_back_pic_vm_1d2v

    
  !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_1d2v(&
       self, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       map, &
       boundary_particles, &
       solver_tolerance, &
       iter_tolerance, max_iter, &
       force_sign, &
       electrostatic, &
       jmean)  
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_1d_base), target,                 intent( in )  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent( in )  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent( in )  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,                  intent( in )  :: particle_group !< Particle group
    sll_real64, target,                                   intent( in )  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                                   intent( in )  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                           intent( in )  :: x_min !< Lower bound of x domain
    sll_real64,                                           intent( in )  :: Lx !< Length of the domain in x direction.
    type(sll_t_mapping_3d), target,                       intent( in )  :: map !< Coordinate transformation
    sll_int32, optional,                                  intent( in )  :: boundary_particles !< particle boundary conditions
    sll_real64, optional,                                 intent( in )  :: solver_tolerance !< solver tolerance
    sll_real64, optional,                          intent( in ) :: iter_tolerance !< iteration tolerance
    sll_int32,  optional,                          intent( in ) :: max_iter !< maximal number of iterations
    sll_real64, optional,                                 intent( in )  :: force_sign !< sign of particle force
    logical, optional,                                    intent( in )  :: electrostatic !< true for electrostatic simulation
    logical, optional,                                    intent( in )  :: jmean !< logical for mean value of current
    !local variables
    sll_int32 :: max_iter_set
    sll_int32 :: boundary_particles_set
    sll_real64 :: solver_tolerance_set, iter_tolerance_set
    sll_real64 :: force_sign_set
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
         map,&
         boundary_particles, &
         solver_tolerance_set, &
         iter_tolerance_set, &
         max_iter_set, &
         force_sign=force_sign_set,&
         jmean=jmean_set)

  end subroutine initialize_pic_vm_1d2v


  !> Constructor.
  subroutine initialize_file_pic_vm_1d2v(&
       self, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       map, &
       filename, &
       boundary_particles, &
       force_sign, &
       electrostatic, &
       jmean)  
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_1d_base), target,                 intent( in )  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent( in )  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent( in )  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,                  intent( in )  :: particle_group !< Particle group
    sll_real64, target,                                   intent( in )  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                                   intent( in )  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                           intent( in )  :: x_min !< Lower bound of x domain
    sll_real64,                                           intent( in )  :: Lx !< Length of the domain in x direction.
    type(sll_t_mapping_3d), target,                       intent( in )  :: map  !< Coordinate transformation  
    character(len=*),                                     intent( in )  :: filename !< filename
    sll_int32, optional, intent( in )  :: boundary_particles !< particle boundary conditions
    sll_real64, optional, intent(in) :: force_sign !< sign of particle force
    logical, optional, intent(in)    :: electrostatic !< true for electrostatic simulation
    logical, optional, intent(in)    :: jmean !< logical for mean value of current
    !local variables
    sll_int32 :: boundary_particles_set
    sll_real64 :: force_sign_set
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
         map,&
         filename, &
         boundary_particles = boundary_particles_set, &
         force_sign=force_sign_set, &
         jmean=jmean_set)


  end subroutine initialize_file_pic_vm_1d2v

  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_1d2v(self)
    class(sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo), intent( inout ) :: self !< time propagator object

    call self%helper%free()

  end subroutine delete_pic_vm_1d2v



end module sll_m_time_propagator_pic_vm_1d2v_disgradE_trafo
