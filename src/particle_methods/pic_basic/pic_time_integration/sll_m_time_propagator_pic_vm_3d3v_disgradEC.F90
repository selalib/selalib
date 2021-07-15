!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on antisymmetric splitting with discrete gradient method for 3d3v Vlasov-Maxwell.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_time_propagator_pic_vm_3d3v_disgradEC
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use sll_m_time_propagator_pic_vm_3d3v_helper, only: &
       sll_t_time_propagator_pic_vm_3d3v_helper

  use sll_m_particle_mesh_coupling_base_3d, only: &
       sll_c_particle_mesh_coupling_3d

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base, &
       sll_t_particle_array

  use sll_m_filter_base_3d, only: &
       sll_c_filter_base_3d

  use sll_m_control_variate, only: &
       sll_t_control_variates

  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_3d3v_disgradEC

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Time propagator for Vlasov-Maxwell 3d3v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_3d3v_disgradEC

     type(sll_t_time_propagator_pic_vm_3d3v_helper) :: helper
     logical    :: electrostatic = .false.

   contains
     procedure :: lie_splitting => lie_splitting_pic_vm_3d3v !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_3d3v !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_3d3v !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_3d3v !> Initialize the type
     procedure :: init_from_file => initialize_file_pic_vm_3d3v !> Initialize the type
     procedure :: free => delete_pic_vm_3d3v !> Finalization

  end type sll_t_time_propagator_pic_vm_3d3v_disgradEC

contains


  !> Strang splitting
  subroutine strang_splitting_pic_vm_3d3v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step
    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_ex ( dt )  
          call self%helper%advect_vb( dt*0.5_f64 )
       end do
    else
       do i_step = 1, number_steps
          call self%helper%advect_eb( dt*0.5_f64 )
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_ex ( dt )  
          call self%helper%advect_vb( dt*0.5_f64 )
          call self%helper%advect_eb( dt*0.5_f64 )
       end do
    end if

  end subroutine strang_splitting_pic_vm_3d3v


  !> Lie splitting
  subroutine lie_splitting_pic_vm_3d3v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

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
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

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
       filter, &
       solver_tolerance, &
       betar, &
       force_sign, &
       electrostatic, &
       control_variate, &
       jmean) 
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,          intent( in ) :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target,   intent(in)  :: particle_mesh_coupling !< Kernel smoother
    class(sll_t_particle_array), target,           intent( in ) :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:) !< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent( in ) :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent( in ) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                    intent( in ) :: x_min(3) !< Lower bound of x domain
    sll_real64,                                    intent( in ) :: Lx(3) !< Length of the domain in x direction.
    class(sll_c_filter_base_3d), target :: filter !< filter
    sll_real64, optional,                          intent( in ) :: solver_tolerance !< solver tolerance
    sll_real64, optional, intent(in) :: betar(2) !< reciprocal plasma beta
    sll_real64, optional, intent(in) :: force_sign !< sign of particle force
    logical, optional, intent(in)    :: electrostatic !< true for electrostatic simulation
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    logical, optional, intent(in)    :: jmean !< logical for mean value of current
    !local variables
    sll_real64 :: solver_tolerance_set
    sll_real64 :: betar_set(2), force_sign_set
    logical :: jmean_set

    if (present(solver_tolerance) )  then
       solver_tolerance_set = solver_tolerance
    else
       solver_tolerance_set = 1d-12
    end if

    if (present(betar) )  then
       betar_set = betar
    else
       betar_set = 1d0
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

    if( present( control_variate) ) then
       call self%helper%init(maxwell_solver, &
            particle_mesh_coupling, &
            particle_group, &
            phi_dofs, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx, &
            filter, &
            solver_tolerance=solver_tolerance_set, &
            betar=betar_set, &
            force_sign=force_sign_set, &
            control_variate = control_variate,&
            jmean=jmean_set)
    else
       call self%helper%init(maxwell_solver, &
            particle_mesh_coupling, &
            particle_group, &
            phi_dofs, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx, &
            filter, &
            solver_tolerance=solver_tolerance_set, &
            betar=betar_set, &
            force_sign=force_sign_set, &
            jmean=jmean_set)
    end if

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
       filter, &
       filename, &
       betar, &
       force_sign, &
       electrostatic, &
       control_variate, &
       jmean) 
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,          intent( in ) :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target, intent(in) :: particle_mesh_coupling !< Kernel smoother
    class(sll_t_particle_array), target,           intent( in ) :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:) !< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent( in ) :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent( in ) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                    intent( in ) :: x_min(3) !< Lower bound of x domain
    sll_real64,                                    intent( in ) :: Lx(3) !< Length of the domain in x direction.
    class(sll_c_filter_base_3d), target :: filter !< filter
    character(len=*),                              intent( in ) :: filename !< filename
    sll_real64, optional,                          intent( in ) :: betar(2) !< reciprocal plasma beta
    sll_real64, optional, intent(in) :: force_sign !< sign of particle force
    logical, optional    :: electrostatic !< true for electrostatic simulation
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    logical, optional, intent(in)    :: jmean !< logical for mean value of current
    !local variables
    sll_real64 :: betar_set(2), force_sign_set
    logical :: jmean_set

    if( present(force_sign) )then
       force_sign_set = force_sign
    else
       force_sign_set = 1._f64
    end if

    if( present(electrostatic) )then
       self%electrostatic = electrostatic
    end if

    if( present(betar) )then
       betar_set = betar
    else
       betar_set = 1d0
    end if

    if (present(jmean)) then
       jmean_set = jmean
    else
       jmean_set = .false.
    end if

    if( present( control_variate) ) then
       call self%helper%init_from_file(  maxwell_solver, &
            particle_mesh_coupling, &
            particle_group, &
            phi_dofs, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx, &
            filename, &
            filter, &
            betar=betar_set, &
            force_sign=force_sign_set, &
            control_variate = control_variate, &
            jmean=jmean_set)
    else
       call self%helper%init_from_file(  maxwell_solver, &
            particle_mesh_coupling, &
            particle_group, &
            phi_dofs, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx, &
            filename, &
            filter, &
            betar=betar_set, &
            force_sign=force_sign_set, &
            jmean=jmean_set)
    end if

  end subroutine initialize_file_pic_vm_3d3v


  !> Destructor.
  subroutine delete_pic_vm_3d3v(self)
    class(sll_t_time_propagator_pic_vm_3d3v_disgradEC), intent( inout ) :: self !< time propagator object 
    call self%helper%free()

  end subroutine delete_pic_vm_3d3v


end module sll_m_time_propagator_pic_vm_3d3v_disgradEC
