!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on the subcycling algorithm for the 1d2v Vlasov-Maxwell equation with splitting of the three H_p parts
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
!> Reference: Hirvijoki, Kormann, Zonta, Subcycling of particle orbits in variational, geometric electromagnetic particle-in-cell methods.
!> Control variate: Note the we do not account for the analytic j at the moment (TODO: control_variate for current)
module sll_m_time_propagator_pic_vm_1d2v_zigsub
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_o_collective_allreduce, &
    sll_v_world_collective

  use sll_m_binomial_filter, only: &
       sll_t_binomial_filter

  use sll_m_control_variate, only: &
    sll_t_control_variates

  use sll_m_time_propagator_base, only: &
    sll_c_time_propagator_base

  use sll_m_particle_mesh_coupling_base_1d, only: &
    sll_c_particle_mesh_coupling_1d

  use sll_m_maxwell_1d_base, only: &
    sll_c_maxwell_1d_base

  use sll_m_particle_group_base, only: &
    sll_t_particle_array

  use sll_mpi, only: &
       mpi_sum

  implicit none

  public :: &
    sll_s_new_time_propagator_pic_vm_1d2v_zigsub, &
    sll_s_new_time_propagator_pic_vm_1d2v_zigsub_ptr, &
    sll_t_time_propagator_pic_vm_1d2v_zigsub

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_1d2v_zigsub
     class(sll_c_maxwell_1d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_0  !< Kernel smoother (order p+1)
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_1  !< Kernel smoother (order p)
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

     sll_int32 :: spline_degree !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx !< Size of the domain
     sll_real64 :: x_min !< Lower bound for x domain
     sll_real64 :: delta_x !< Grid spacing

     sll_real64 :: cell_integrals_0(4) !< Integral over the spline function on each interval (order p+1)
     sll_real64 :: cell_integrals_1(3) !< Integral over the spline function on each interval (order p)


     sll_real64, pointer     :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
     sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the magnetic field
     sll_real64, allocatable :: j_dofs(:,:)      !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:,:)!< MPI-processor local part of one component of \a j_dofs
     sll_int32 :: n_species !< number of species

     logical :: electrostatic = .false. !< true for electrostatic simulation

     logical :: jmean = .false. !< logical for mean value of current
     sll_real64, allocatable     :: bfield_old(:)   !< DoFs describing the magnetic field

     type(sll_t_binomial_filter), pointer :: filter !< filter

     ! For version with control variate
     class(sll_t_control_variates), pointer :: control_variate !< control variate
     sll_int32 :: i_weight !< number of weights

     sll_int32 :: n_sub_iter !< number of subiterations
     
   contains
     procedure :: lie_splitting => lie_splitting_pic_vm_1d2v !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_1d2v !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_1d2v !> Strang splitting propagator
     procedure :: reinit_fields

     procedure :: init => initialize_pic_vm_1d2v !> Initialize the type
     procedure :: free => delete_pic_vm_1d2v !> Finalization

  end type sll_t_time_propagator_pic_vm_1d2v_zigsub

contains

  subroutine reinit_fields( self ) 
    class(sll_t_time_propagator_pic_vm_1d2v_zigsub), intent(inout) :: self !< time propagator object 
    
  end subroutine reinit_fields
  
  !> Strang splitting
  subroutine strang_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_zigsub), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call operator_all(self, dt)
       end do
    else
       do i_step = 1, number_steps
          call operator_all(self, dt)
       end do
    end if
    
  end subroutine strang_splitting_pic_vm_1d2v

  !> Lie splitting
  subroutine lie_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_zigsub), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step
    if(self%electrostatic) then
       do i_step = 1, number_steps
          call operator_all(self, dt)
       end do
    else
       do i_step = 1,number_steps
          call operator_all(self, dt)
       end do
    end if


  end subroutine lie_splitting_pic_vm_1d2v

  !> Lie splitting (oposite ordering)
  subroutine lie_splitting_back_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_zigsub), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call operator_all(self, dt)
       end do
    else
       do i_step = 1,number_steps
          call operator_all(self, dt)
       end do
    end if
 
  end subroutine lie_splitting_back_pic_vm_1d2v

  subroutine operator_all( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_zigsub), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_sp, i_part
    sll_real64 :: bfield, efield(2)
    sll_real64 :: dtau
    sll_real64 :: qoverm
    sll_real64 :: wi(1), vi(3), xi(3), wp(3), x_new(3)
    sll_int32  :: n_cells, iter
    
    dtau = dt/real(self%n_sub_iter, f64)
    
    n_cells = self%kernel_smoother_0%n_dofs

    ! Update bfield
    self%bfield_old = self%bfield_dofs
    call self%maxwell_solver%compute_B_from_E( &
         dt, self%efield_dofs(:,2), self%bfield_dofs)

    ! Particle loop with subcycling
    self%j_dofs_local = 0.0_f64
    ! Then update particle position:  X_new = X_old + dt * V
    do i_sp = 1,self%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          
          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

          ! First push v_1
          ! Evaluate bfield at particle position (splines of order p)
          call self%kernel_smoother_1%evaluate &
               (xi(1), self%bfield_old, bfield)
          vi(1) = vi(1) + dtau*qoverm*vi(2)*bfield
          ! Evaluate efields at particle position
          call self%kernel_smoother_1%evaluate &
               (xi(1), self%efield_dofs(:,1), efield(1))
          vi(1) = vi(1) + dt* qoverm * efield(1)
          
          call self%kernel_smoother_0%add_charge(xi(1:1), wi(1)*vi(2)*dtau, self%j_dofs_local(:,2)) 
                    
          if (self%particle_group%group(i_sp)%n_weights == 3) then
             ! Update weights if control variate
             wp = self%particle_group%group(i_sp)%get_weights(i_part)          
             wp(3) = self%control_variate%cv(i_sp)%update_df_weight( xi(1:1), vi(1:2), 0.0_f64, wp(1), wp(2))
             call self%particle_group%group(i_sp)%set_weights(i_part, wp)
          
             ! Get charge for accumulation of j
             wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)
          end if
          
          ! Then update particle position:  X_new = X_old + dt * V
          x_new = xi  + dtau * vi
          
          call self%kernel_smoother_1%add_current_update_v( xi, x_new, wi(1), qoverm, &
               self%bfield_dofs, vi, self%j_dofs_local(:,1))
          ! Accumulate rho for Poisson diagnostics
          !call self%kernel_smoother_0%add_charge( x_new, wi(1), &
          !     self%j_dofs_local(:,2))
          ! Evaluate efields at particle position
          call self%kernel_smoother_0%evaluate &
               (xi(1), self%efield_dofs(:,2), efield(2))
          vi(2) = vi(2) + dt* qoverm * efield(2)
          
          x_new(1) = modulo(x_new(1), self%Lx)
     
          if (self%particle_group%group(i_sp)%n_weights == 3) then
             ! Update weights if control variate
             wp = self%particle_group%group(i_sp)%get_weights(i_part)          
             wp(3) = self%control_variate%cv(i_sp)%update_df_weight( x_new(1:1), vi(1:2), 0.0_f64, wp(1), wp(2))
             call self%particle_group%group(i_sp)%set_weights(i_part, wp)
          end if

          do iter = 2, self%n_sub_iter
             xi = x_new
             ! First push v_1
             ! Evaluate bfield at particle position (splines of order p)
             call self%kernel_smoother_1%evaluate &
                  (xi(1), self%bfield_dofs, bfield)
             vi(1) = vi(1) + dtau*qoverm*vi(2)*bfield
          
             call self%kernel_smoother_0%add_charge(xi(1:1), wi(1)*vi(2)*dtau, self%j_dofs_local(:,2)) 
                    
             if (self%particle_group%group(i_sp)%n_weights == 3) then
                ! Update weights if control variate
                wp = self%particle_group%group(i_sp)%get_weights(i_part)          
                wp(3) = self%control_variate%cv(i_sp)%update_df_weight( xi(1:1), vi(1:2), 0.0_f64, wp(1), wp(2))
                call self%particle_group%group(i_sp)%set_weights(i_part, wp)
                
                ! Get charge for accumulation of j
                wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)
             end if
          
             ! Then update particle position:  X_new = X_old + dt * V
             x_new = xi  + dtau * vi
          
             call self%kernel_smoother_1%add_current_update_v( xi, x_new, wi(1), qoverm, &
                  self%bfield_dofs, vi, self%j_dofs_local(:,1))
             ! Accumulate rho for Poisson diagnostics
             !call self%kernel_smoother_0%add_charge( x_new, wi(1), &
             !     self%j_dofs_local(:,2))        
             x_new(1) = modulo(x_new(1), self%Lx)
             
             if (self%particle_group%group(i_sp)%n_weights == 3) then
                ! Update weights if control variate
                wp = self%particle_group%group(i_sp)%get_weights(i_part)          
                wp(3) = self%control_variate%cv(i_sp)%update_df_weight( x_new(1:1), vi(1:2), 0.0_f64, wp(1), wp(2))
                call self%particle_group%group(i_sp)%set_weights(i_part, wp)
             end if
          end do

          call self%particle_group%group(i_sp)%set_x(i_part, x_new)
          call self%particle_group%group(i_sp)%set_v(i_part, vi)
       end do
    end do

    ! Update the electric field
    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
         n_cells, MPI_SUM, self%j_dofs(:,1))
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
         n_cells, MPI_SUM, self%j_dofs(:,2))
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,1), 1, self%efield_dofs(:,1))
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,2), 2, self%efield_dofs(:,2))
    call self%maxwell_solver%compute_E_from_B(&
         dt, self%bfield_dofs, self%efield_dofs(:,2))


  end subroutine operator_all
 


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
       n_sub_iter, &
       filter, &
       jmean, &
       control_variate, &
       i_weight, &
       electrostatic) 
    class(sll_t_time_propagator_pic_vm_1d2v_zigsub), intent(out) :: self !< time propagator object 
    class(sll_c_maxwell_1d_base), target,          intent(in)  :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,           intent(in)  :: particle_group !< particle group
    sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
    sll_int32,                                     intent(in)  :: n_sub_iter !< number of subiterations
    type( sll_t_binomial_filter ), intent( in ), target :: filter !< filter
    logical, optional, intent(in) :: jmean !< logical for mean value of current
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    sll_int32, optional,                            intent(in) :: i_weight !< Index of weight to be used by propagator
    logical, optional    :: electrostatic !< true for electrostatic simulation
    !local variables
    sll_int32 :: ierr, j

    if( present(electrostatic) )then
       self%electrostatic = electrostatic
    end if

    self%maxwell_solver => maxwell_solver
    self%kernel_smoother_0 => kernel_smoother_0
    self%kernel_smoother_1 => kernel_smoother_1

    self%n_species = particle_group%n_species
    
    self%particle_group => particle_group
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs

    ! Check that n_dofs is the same for both kernel smoothers.
    SLL_ASSERT( self%kernel_smoother_0%n_dofs == self%kernel_smoother_1%n_dofs )

    SLL_ALLOCATE(self%j_dofs(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%j_dofs_local(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%bfield_old(self%kernel_smoother_0%n_dofs), ierr)

    self%spline_degree = self%kernel_smoother_0%spline_degree
    self%x_min = x_min
    self%Lx = Lx
    self%delta_x = self%Lx/self%kernel_smoother_1%n_dofs
    
    self%cell_integrals_1 = [0.5_f64, 2.0_f64, 0.5_f64]
    self%cell_integrals_1 = self%cell_integrals_1 / 3.0_f64

    self%cell_integrals_0 = [1.0_f64,11.0_f64,11.0_f64,1.0_f64]
    self%cell_integrals_0 = self%cell_integrals_0 / 24.0_f64

    
    if (present(jmean)) then
       self%jmean = jmean
    end if

    self%i_weight = 1
    if (present(i_weight)) self%i_weight = i_weight
    if(present(control_variate)) then
       allocate(self%control_variate )
       allocate(self%control_variate%cv(self%n_species) )
       self%control_variate => control_variate
    end if

    self%n_sub_iter = n_sub_iter
    
  end subroutine initialize_pic_vm_1d2v

  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_1d2v(self)
    class(sll_t_time_propagator_pic_vm_1d2v_zigsub), intent( inout ) :: self !< time propagator object 

    deallocate(self%j_dofs)
    deallocate(self%j_dofs_local)
    self%maxwell_solver => null()
    self%kernel_smoother_0 => null()
    self%kernel_smoother_1 => null()
    self%particle_group => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()

  end subroutine delete_pic_vm_1d2v


  !---------------------------------------------------------------------------!
  !> Constructor for allocatable abstract type.
  subroutine sll_s_new_time_propagator_pic_vm_1d2v_zigsub(&
       splitting, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       n_sub_iter, &
       filter, &
       jmean, &
       control_variate, &
       i_weight, &
       electrostatic) 
    class(sll_c_time_propagator_base), allocatable, intent(out) :: splitting !< time propagator object 
    class(sll_c_maxwell_1d_base), target,                intent(in)  :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,  intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,  intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,           intent(in)  :: particle_group !< Particle group
    sll_real64, target,                                  intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                                  intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                           intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                           intent(in)  :: Lx !< Length of the domain in x direction.
    sll_int32,                                     intent(in)  :: n_sub_iter !< number of subiterations
    type( sll_t_binomial_filter ), intent( in ), target :: filter !< filter
    logical, optional, intent(in) :: jmean !< Should jmean be substracted in Ampere's law?
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    sll_int32, optional,                            intent(in) :: i_weight !< Index of weight to be used by propagator
    logical, optional    :: electrostatic !< true for electrostatic simulation
    !local variables
    sll_int32 :: ierr
    logical :: jmean_val 

    SLL_ALLOCATE(sll_t_time_propagator_pic_vm_1d2v_zigsub :: splitting, ierr)

    if (present(jmean) ) then
       jmean_val = jmean
    else
       jmean_val = .false.
    end if
    
    select type (splitting)
    type is ( sll_t_time_propagator_pic_vm_1d2v_zigsub )
       if (present(control_variate) ) then
          call splitting%init(&
               maxwell_solver, &
               kernel_smoother_0, &
               kernel_smoother_1, &
               particle_group, &
               efield_dofs, &
               bfield_dofs, &
               x_min, &
               Lx, &
               n_sub_iter, &
               filter, &
               jmean_val, &
               control_variate, &
               i_weight, &
               electrostatic=electrostatic)
       else
          call splitting%init(&
               maxwell_solver, &
               kernel_smoother_0, &
               kernel_smoother_1, &
               particle_group, &
               efield_dofs, &
               bfield_dofs, &
               x_min, &
               Lx, &
               n_sub_iter, &
               filter, &
               jmean_val, &
               electrostatic=electrostatic)
       end if
    end select

  end subroutine sll_s_new_time_propagator_pic_vm_1d2v_zigsub

  !---------------------------------------------------------------------------!
  !> Constructor for pointer abstract type.
  subroutine sll_s_new_time_propagator_pic_vm_1d2v_zigsub_ptr(&
       splitting, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       n_sub_iter, &
       filter, &
       jmean, &
       electrostatic) 
    class(sll_c_time_propagator_base), pointer, intent(out) :: splitting !< time propagator object 
    class(sll_c_maxwell_1d_base), target,            intent(in)  :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,            intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,            intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,           intent(in)  :: particle_group !< Particle group
    sll_real64, target,                              intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                              intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                       intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                       intent(in)  :: Lx !< Length of the domain in x direction.
    sll_int32,                                     intent(in)  :: n_sub_iter !< number of subiterations
    type( sll_t_binomial_filter ), intent( in ), target :: filter !< filter
    logical, optional, intent(in) :: jmean !< logical for mean value of current
    logical, optional    :: electrostatic  !< true for electrostatic simulation
    !local variables
    sll_int32 :: ierr
    logical :: jmean_val

    SLL_ALLOCATE(sll_t_time_propagator_pic_vm_1d2v_zigsub :: splitting, ierr)


    if (present(jmean) ) then
       jmean_val = jmean
    else
       jmean_val = .false.
    end if
    
    select type (splitting)
    type is ( sll_t_time_propagator_pic_vm_1d2v_zigsub )
       call splitting%init(&
            maxwell_solver, &
            kernel_smoother_0, &
            kernel_smoother_1, &
            particle_group, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx, &
            n_sub_iter, &
            filter, &
            jmean_val, &
            electrostatic=electrostatic)
    end select

  end subroutine sll_s_new_time_propagator_pic_vm_1d2v_zigsub_ptr


end module sll_m_time_propagator_pic_vm_1d2v_zigsub
