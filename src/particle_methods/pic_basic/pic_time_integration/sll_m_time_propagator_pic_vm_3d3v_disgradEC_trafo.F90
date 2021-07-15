!> @ingroup pic_time_integration
!> @author Benedikt Perse, IPP
!> @brief Particle pusher based on antisymmetric splitting with discrete gradient method for 3d3v Vlasov-Maxwell with coordinate transformation.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_time_propagator_pic_vm_3d3v_disgradEC_trafo
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use sll_m_time_propagator_pic_vm_3d3v_trafo_helper, only:&
       sll_t_time_propagator_pic_vm_3d3v_trafo_helper

  use sll_m_particle_mesh_coupling_base_3d, only: &
       sll_c_particle_mesh_coupling_3d

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base, &
       sll_t_particle_array

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_collective, only: &
       sll_f_get_collective_rank, &
       sll_v_world_collective

  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Time propagator for Vlasov-Maxwell 3d3v with coordinate transformation
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo
     type(sll_t_time_propagator_pic_vm_3d3v_trafo_helper) :: helper
     logical    :: electrostatic = .false.

   contains
     procedure :: lie_splitting => lie_splitting_pic_vm_3d3v !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_3d3v !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_3d3v !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_3d3v !> Initialize the type
     procedure :: init_from_file => initialize_file_pic_vm_3d3v !> Initialize the type
     procedure :: free => delete_pic_vm_3d3v !> Finalization

  end type sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo

contains


  !> Strang splitting
  subroutine strang_splitting_pic_vm_3d3v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
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


  end subroutine strang_splitting_pic_vm_3d3v


  !> Lie splitting
  subroutine lie_splitting_pic_vm_3d3v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
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

  end subroutine lie_splitting_pic_vm_3d3v


  !> Lie splitting (oposite ordering)
  subroutine lie_splitting_back_pic_vm_3d3v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: dt   !< time step
    sll_int32,                                            intent( in    ) :: number_steps !< number of time steps
    !local variables
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

  end subroutine lie_splitting_back_pic_vm_3d3v

  !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_3d3v(&
       self, &
       maxwell_solver, &
       particle_mesh_coupling, &
       particle_group, &
       phi_dofs, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       map, &
       solver_tolerance, &
       iter_tolerance, max_iter, &
       betar, &
       electrostatic) 
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,          intent( in ) :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target, intent(in) :: particle_mesh_coupling !< Kernel smoother
    class(sll_t_particle_array), target,           intent( in ) :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:) !< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent( in ) :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent( in ) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                    intent( in ) :: x_min(3) !< Lower bound of x domain
    sll_real64,                                    intent( in ) :: Lx(3) !< Length of the domain in x direction.
    type(sll_t_mapping_3d), target,                intent( inout ) :: map !< Coordinate transformation
    sll_real64, optional,                          intent( in ) :: solver_tolerance !< Solver tolerance
    sll_real64, optional,                          intent( in ) :: iter_tolerance !< iteration tolerance
    sll_int32,  optional,                          intent( in ) :: max_iter !< maximal number of iterations
    sll_real64, optional,                          intent( in ) :: betar(2) !< reciprocal plasma beta
    logical, optional    :: electrostatic !< true for electrostatic simulation
    !local variables
    sll_int32 :: max_iter_set
    sll_real64 :: iter_tolerance_set
    sll_real64 :: betar_set(2)

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

    if (present(betar) )  then
       betar_set = betar
    else
       betar_set = 1d0
    end if

    if( present(electrostatic) )then
       self%electrostatic = electrostatic
    end if

    call self%helper%init(maxwell_solver, &
         particle_mesh_coupling, &
         particle_group, &
         phi_dofs, &
         efield_dofs, &
         bfield_dofs, &
         x_min, &
         Lx, &
         map,  &
         iter_tolerance=iter_tolerance_set, &
         max_iter=max_iter_set, &
         betar=betar_set)

  end subroutine initialize_pic_vm_3d3v


  !> Constructor.
  subroutine initialize_file_pic_vm_3d3v(&
       self, &
       maxwell_solver, &
       particle_mesh_coupling, &
       particle_group, &
       phi_dofs, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       map, &
       filename, &
       betar, &
       electrostatic) 
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,          intent( in ) :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target, intent(in) :: particle_mesh_coupling !< Kernel smoother
    class(sll_t_particle_array), target,           intent( in ) :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:) !< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent( in ) :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent( in ) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                    intent( in ) :: x_min(3) !< Lower bound of x domain
    sll_real64,                                    intent( in ) :: Lx(3) !< Length of the domain in x direction.
    type(sll_t_mapping_3d), target,                intent( inout ) :: map !< Coordinate transformation
    character(len=*),                              intent( in ) :: filename !< filename 
    sll_real64, optional,                          intent( in ) :: betar(2) !< reciprocal plasma beta
    logical, optional    :: electrostatic !< true for electrostatic simulation
    !local variables
    sll_int32 :: rank, file_id, io_stat
    sll_int32 :: input_file
    sll_real64 :: iter_tolerance
    sll_int32 :: max_iter
    sll_real64 :: betar_set(2)

    if( present(electrostatic) )then
       self%electrostatic = electrostatic
    end if

    if (present(betar) )  then
       betar_set = betar
    else
       betar_set = 1d0
    end if

    rank = sll_f_get_collective_rank(sll_v_world_collective)
    if (rank == 0 ) then
       open(newunit=file_id, file=trim(filename)//'_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
       write(file_id, *) 'electrostatic:', self%electrostatic
       close(file_id) 
    end if

    call self%helper%init_from_file(  maxwell_solver, &
         particle_mesh_coupling, &
         particle_group, &
         phi_dofs, &
         efield_dofs, &
         bfield_dofs, &
         x_min, &
         Lx, &
         map, &
         filename, &
         betar=betar_set)

  end subroutine initialize_file_pic_vm_3d3v


  !> Destructor.
  subroutine delete_pic_vm_3d3v(self)
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo), intent( inout ) :: self !< time propagator object

    call self%helper%free()

  end subroutine delete_pic_vm_3d3v


end module sll_m_time_propagator_pic_vm_3d3v_disgradEC_trafo
