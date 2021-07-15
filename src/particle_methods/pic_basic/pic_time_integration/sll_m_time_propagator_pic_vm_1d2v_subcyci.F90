!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on the subcycling algorithm for the 1d2v Vlasov-Maxwell equation
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
!> Reference: Hirvijoki, Kormann, Zonta, Subcycling of particle orbits in variational, geometric electromagnetic particle-in-cell methods.

!> Control variate: Note the we do not account for the analytic j at the moment (TODO: control_variate for current)
!> Explanation to the various versions:
!> 1.) first_step_modified2: This has to be ran in the version with Faraday's equation solved backwards in time
!> 2.) first_step_modified: This has to be ran in the version described in the notes (section 2.5 but with x_-1 instead of x_1.
!> 3.) first_step: This is the original version.
!> The flags step_one(_mod) are set to fit. Between the various version can be switched by
!> calling the corresponding first step function in the simulation.


module sll_m_time_propagator_pic_vm_1d2v_subcyci
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

  use sll_m_particle_mesh_coupling_spline_1d, only : &
       sll_t_particle_mesh_coupling_spline_1d

  use sll_m_maxwell_1d_base, only: &
       sll_c_maxwell_1d_base
  
  use sll_m_maxwell_1d_fem, only: &
       sll_t_maxwell_1d_fem

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use sll_mpi, only: &
       mpi_sum

  use sll_m_spline_fem_utilities, only : &
       sll_s_spline_fem_mass_line
  
  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr
 
  use sll_m_linear_solver_mgmres, only : &
       sll_t_linear_solver_mgmres
  
  implicit none

  public :: &
       sll_s_new_time_propagator_pic_vm_1d2v_subcyci, &
       sll_s_new_time_propagator_pic_vm_1d2v_subcyci_ptr, &
       sll_t_time_propagator_pic_vm_1d2v_subcyci

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_1d2v_subcyci
     class(sll_c_maxwell_1d_base), pointer :: maxwell_solver      !< Maxwell solver
     !class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_0  !< Kernel smoother (order p+1)
     !class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_1  !< Kernel smoother (order p)
     class(sll_t_particle_mesh_coupling_spline_1d), pointer :: kernel_smoother_0  !< Kernel smoother (order p+1)
     class(sll_t_particle_mesh_coupling_spline_1d), pointer :: kernel_smoother_1  !< Kernel smoother (order p)
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

     sll_int32 :: spline_degree !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx !< Size of the domain
     sll_real64 :: x_min !< Lower bound for x domain
     sll_real64 :: delta_x !< Grid spacing

     sll_real64 :: cell_integrals_0(4) !< Integral over the spline function on each interval (order p+1)
     sll_real64 :: cell_integrals_1(3) !< Integral over the spline function on each interval (order p)


     sll_real64, pointer     :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
     sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the magnetic field
     sll_real64, allocatable :: j1_dofs(:,:)      !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j1_dofs_local(:,:)!< MPI-processor local part of one component of \a j_dofs
     sll_real64, allocatable :: j2_dofs(:,:)      !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j2_dofs_local(:,:)!< MPI-processor local part of one component of \a j_dofs
     sll_real64, allocatable :: rho_dofs_local(:) 
     sll_real64, allocatable :: rho_dofs(:) 
     sll_int32 :: n_species

     logical :: jmean = .false.
     sll_real64, allocatable :: efield_filter(:,:) !< DoFs describing the two components of the electric field
     sll_real64, allocatable :: bfield_filter(:)   !< DoFs describing the magnetic field
     sll_real64, allocatable :: bfield_dofs_old(:) 
     sll_real64, allocatable :: bfield_linear_combination(:)
     sll_real64, allocatable :: efield_tmp(:)
     sll_real64, allocatable     :: bfield_dofs_even_older(:)
     sll_real64, allocatable     :: efield_dofs_old(:,:)
     sll_real64, allocatable     :: j1_dofs_old(:)
     sll_real64, allocatable     :: j2_dofs_old(:)

     sll_real64, allocatable :: xvec(:,:)
     sll_real64, allocatable :: vvec(:,:,:)

     type(sll_t_binomial_filter), pointer :: filter

     ! For version with control variate
     class(sll_t_control_variates), pointer :: control_variate
     sll_int32 :: i_weight

     sll_int32 :: n_sub_iter = 4
     sll_int32 :: n_max_iter = 50
     sll_int32 :: n_max_iter_field = 30
     sll_real64 :: tolerance = 1D-10
     sll_real64 :: tolerance_field = 1D-13
     sll_int32 :: file_id
     sll_int32 :: file_id_rho
     sll_int32 :: file_id_iter_field

     logical   :: step_one = .false.
     logical   :: step_one_mod = .false.

     ! For the implicit field solve
     type(sll_t_matrix_csr) :: lhs_eb
     type(sll_t_matrix_csr) :: rhs_eb
     type(sll_t_linear_solver_mgmres) :: linear_solver_eb

   contains
     procedure :: operatorall => operatorHp_pic_vm_1d2v_newton  !> Operator for H_p1 part
     procedure :: lie_splitting => lie_splitting_pic_vm_1d2v !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_1d2v !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_1d2v !> Strang splitting propagator
     procedure :: reinit_fields !> Apply the filter to set the dofs of the filtered fields

     procedure :: init => initialize_pic_vm_1d2v !> Initialize the type
     procedure :: free => delete_pic_vm_1d2v !> Finalization
     procedure :: first_step
     procedure :: first_step_modified
     procedure :: first_step_modified2

  end type sll_t_time_propagator_pic_vm_1d2v_subcyci

contains

  !> Apply the filter to set the dofs of the filtered fields
  subroutine reinit_fields( self ) 
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(inout) :: self !< time splitting object 

    call self%filter%apply( self%efield_dofs(:,1), self%efield_filter(:,1) ) 
    call self%filter%apply( self%efield_dofs(:,2), self%efield_filter(:,2) ) 
    call self%filter%apply( self%bfield_dofs, self%bfield_filter )

  end subroutine reinit_fields

  !> Strang splitting
  subroutine strang_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step


    do i_step = 1, number_steps
       call self%operatorall( dt )
    end do

  end subroutine strang_splitting_pic_vm_1d2v

  !> Lie splitting
  subroutine lie_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1,number_steps
       call self%operatorall( dt )
    end do


  end subroutine lie_splitting_pic_vm_1d2v

  !> Lie splitting (oposite ordering)
  subroutine lie_splitting_back_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1,number_steps
       call self%operatorall( dt )
    end do

  end subroutine lie_splitting_back_pic_vm_1d2v



  ! Implementation of the full time step with Newton iteration for the nonlinear part
  !---------------------------------------------------------------------------!
  subroutine operatorHp_pic_vm_1d2v_newton(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
  !  sll_int32 :: i_part
    sll_real64 :: x_new(3), v_new(3)
    sll_int32  :: n_cells, i_sp, i_part, n_iter_inner
    sll_real64 :: residual_field, residual_field_inner
    sll_int32  :: n_iter_field


    n_cells = self%kernel_smoother_0%n_dofs

    self%efield_dofs_old = self%efield_dofs
    ! Update the magnetic field
    if ( self%step_one .eqv. .false. ) then 
       self%bfield_dofs_even_older = self%bfield_dofs_old
       self%j1_dofs_old = self%j1_dofs(:,1)
       self%j2_dofs_old = self%j2_dofs(:,1)
       
    else
       self%bfield_dofs_even_older = self%bfield_dofs
       !   self%j1_dofs_old = self%j1_dofs(:,1)
       !   self%j2_dofs_old = self%j2_dofs(:,1)
    end if
    self%bfield_dofs_old = self%bfield_dofs
    call self%maxwell_solver%compute_B_from_E( &
         dt, self%efield_dofs(:,2), self%bfield_dofs )

    residual_field = 1.0_f64
    n_iter_field = 0
    do while( residual_field > self%tolerance_field .and. n_iter_field < self%n_max_iter_field )
       n_iter_field = n_iter_field + 1

       call particle_loop_subcycle( self, dt )

       ! Update the electric field
       self%efield_filter = self%efield_dofs
       self%efield_dofs = self%efield_dofs_old
       if ( self%step_one .eqv. .false. ) then
          self%j1_dofs(:,2) = self%j1_dofs_old + self%j1_dofs(:,2)
          self%j2_dofs(:,2) = self%j2_dofs_old + self%j2_dofs(:,2)
       else
          self%j1_dofs(:,2) = self%j1_dofs(:,1) + self%j1_dofs(:,2)
          self%j2_dofs(:,2) = self%j2_dofs(:,1) + self%j2_dofs(:,2)
       end if
       call self%maxwell_solver%compute_E_from_j(self%j1_dofs(:,2), 1, self%efield_dofs(:,1))
       call self%maxwell_solver%compute_E_from_j(self%j2_dofs(:,2), 2, self%efield_dofs(:,2))

!!$       ! Explicit field solve
!!$       call self%maxwell_solver%compute_E_from_B(&
!!$            dt, self%bfield_dofs_old, self%efield_dofs(:,2))
!!$       
!!$       self%bfield_dofs = self%bfield_dofs_old
!!$       call self%maxwell_solver%compute_B_from_E( &
!!$            dt, self%efield_dofs(:,2), self%bfield_dofs)
!!$       
!!$       ! Implicit field version (with explicit implementation)
!!$       self%bfield_linear_combination = 2.0_f64/3.0_f64 * self%bfield_dofs_old + &
!!$            ( self%bfield_dofs + self%bfield_dofs_even_older ) / 6.0_f64
!!$       call self%maxwell_solver%compute_E_from_B(&
!!$            dt, self%bfield_linear_combination, self%efield_dofs(:,2))      
!!$       self%bfield_dofs = self%bfield_dofs_old
!!$       call self%maxwell_solver%compute_B_from_E( &
!!$            dt, self%efield_dofs(:,2), self%bfield_dofs)
!!$
!!$       ! Implicit field version with Picard iteration
!!$       n_iter_inner = 0
!!$       residual_field_inner = 1.0_f64
!!$       do while (residual_field_inner > self%tolerance_field*0.1_f64 .and. n_iter_inner < self%n_max_iter_field )
!!$          self%efield_tmp = self%efield_dofs(:,2) 
!!$          self%bfield_linear_combination = 2.0_f64/3.0_f64 * self%bfield_dofs_old + &
!!$               ( self%bfield_dofs + self%bfield_dofs_even_older ) / 6.0_f64
!!$          call self%maxwell_solver%compute_E_from_B(&
!!$               dt, self%bfield_linear_combination, self%efield_tmp )
!!$          self%bfield_linear_combination = self%bfield_dofs ! Save the old value of the field for error control
!!$          self%bfield_dofs = self%bfield_dofs_old
!!$          call self%maxwell_solver%compute_B_from_E( &
!!$               dt, self%efield_tmp, self%bfield_dofs)
!!$          residual_field_inner = sum( (self%bfield_dofs - self%bfield_linear_combination)**2) * self%delta_x
!!$          n_iter_inner = n_iter_inner + 1
!!$       end do
!!$       self%efield_dofs(:,2) = self%efield_tmp
!!$       print*, 'Inner Picard', n_iter_inner, residual_field_inner

       ! Implicit version with Schur complement field solver
       self%bfield_linear_combination = 2.0_f64/3.0_f64 * self%bfield_dofs_old + &
            ( self%bfield_dofs_old + self%bfield_dofs_even_older ) / 6.0_f64
       select type( ms=> self%maxwell_solver )
       type is( sll_t_maxwell_1d_fem) 
          call ms%solve_e_b( dt, self%bfield_linear_combination, self%efield_dofs(:,2) )
       end select
       self%bfield_dofs = self%bfield_dofs_old
       call self%maxwell_solver%compute_b_from_e( dt, self%efield_dofs(:,2), &
            self%bfield_dofs )
       
       residual_field = sum( (self%efield_dofs - self%efield_filter)**2 )* self%delta_x

    end do
    write(self%file_id_iter_field,*) n_iter_field, residual_field
  
    ! Temporary check of Gauss law for the case with substepping in the electric field
    self%rho_dofs = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, self%rho_dofs_local, n_cells, &
         MPI_SUM, self%rho_dofs )
    call self%maxwell_solver%compute_rho_from_e( self%efield_dofs(:,1), self%rho_dofs_local )
    self%rho_dofs = self%rho_dofs/real(self%n_sub_iter,f64)
    self%rho_dofs = self%rho_dofs - sum(self%rho_dofs)/real(n_cells, f64)

    
    write(self%file_id_rho,*) maxval( abs(self%rho_dofs - self%rho_dofs_local ) )

    ! Set the new particle positions
    do i_sp = 1,self%n_species       
       do i_part=1,self%particle_group%group(i_sp)%n_particles
          x_new(1) = self%xvec(i_part, i_sp)
          v_new(1:2) = self%vvec(:,i_part, i_sp)
          call self%particle_group%group(i_sp)%set_x(i_part, x_new )
          call self%particle_group%group(i_sp)%set_v(i_part, v_new )
       end do
    end do

    self%step_one = .false.
    self%step_one_mod = .false.

  end subroutine operatorHp_pic_vm_1d2v_newton

  ! Implementation of the first step with the strategy where S_0 is modified to the explicit algorithm
  subroutine first_step_modified( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(inout) :: self !< time splitting object 
    sll_real64,                                             intent(in)    :: dt   !< time step

    sll_real64 :: dtau
    sll_int32  :: n_cells
    sll_real64 :: x_new(3), x_old(3), vi(3), wi(1)
    sll_int32  :: i_sp, i_part

    self%step_one = .false.
    self%step_one_mod = .true.
       
    dtau =  dt/real(self%n_sub_iter,f64)

    n_cells = self%kernel_smoother_0%n_dofs

    self%bfield_dofs_old = self%bfield_dofs
   ! call self%maxwell_solver%compute_B_from_E( &
   !     dtau, self%efield_dofs(:,2), self%bfield_dofs )

    ! Particle loop for x_1
    self%j1_dofs_local = 0.0_f64
    self%j2_dofs_local = 0.0_f64
    self%rho_dofs_local = 0.0_f64
    do i_sp = 1,self%n_species
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          x_new = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

          ! Then update particle position (first part):  X_new = X_old + dt * 0.5 * V
          x_old = x_new - dtau * vi

          call self%kernel_smoother_1%add_current( x_old(1), x_new(1), wi(1), self%j1_dofs_local(:,1) )
        !  call self%kernel_smoother_1%add_current_split( x_old(1), x_new(1), &
        !       0, 1, wi(1)*vi(1)*dtau, &
        !       self%j1_dofs_local )
          call self%kernel_smoother_0%add_charge_int( x_old(1), x_new(1), wi(1)*vi(2)*dtau, &
               self%j2_dofs_local(:,1) )

          call self%kernel_smoother_0%add_charge( x_new(1), wi(1), self%rho_dofs_local )

        !  call self%particle_group%group(i_sp)%set_x(i_part, x_new)
       end do
    end do

    self%j1_dofs = 0.0_f64
    self%j2_dofs = 0.0_f64
   ! self%j1_dofs_local(:,1) = self%j1_dofs_local(:,1) + self%j1_dofs_local(:,2)
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j1_dofs_local(:,1), &
         n_cells, MPI_SUM, self%j1_dofs_old)
    call sll_o_collective_allreduce( sll_v_world_collective, self%j2_dofs_local(:,1), &
         n_cells, MPI_SUM, self%j2_dofs_old)

    self%rho_dofs = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, self%rho_dofs_local, &
         n_cells, MPI_SUM, self%rho_dofs )

    call self%maxwell_solver%compute_E_from_rho( self%rho_dofs, self%efield_dofs(:,1) )
!    call self%maxwell_solver%compute_B_from_E( &
!        dtau, self%efield_dofs(:,2), self%bfield_dofs )
    
    self%efield_dofs_old = self%efield_dofs
!!$    ! Initial guess for e_1 and b_2
!!$    self%bfield_dofs_even_older = self%bfield_dofs_old
!!$    self%bfield_dofs_old = self%bfield_dofs
!!$    call self%maxwell_solver%compute_B_from_E( &
!!$        dt, self%efield_dofs(:,2), self%bfield_dofs )

 

  end subroutine first_step_modified

    ! Implementation of the first step with the implicit strategy for a first step that is as small as the substep
  subroutine first_step_modified2( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(inout) :: self !< time splitting object 
    sll_real64,                                             intent(in)    :: dt   !< time step

    sll_real64 :: dtau
    sll_int32  :: n_cells
    sll_real64 :: x_new(3), x_old(3), vi(3), wi(1)
    sll_int32  :: i_sp, i_part

    ! In this case, we prepare everything for a normal first step
    self%step_one_mod = .false.
    self%step_one = .false.
    
    dtau =  dt/real(self%n_sub_iter,f64)

    n_cells = self%kernel_smoother_0%n_dofs

    self%bfield_dofs_old = self%bfield_dofs
    call self%maxwell_solver%compute_B_from_E( &
        -dtau, self%efield_dofs(:,2), self%bfield_dofs_old )

    ! Particle loop for x_1
    self%j1_dofs_local = 0.0_f64
    self%j2_dofs_local = 0.0_f64
    self%rho_dofs_local = 0.0_f64
    do i_sp = 1,self%n_species
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          x_old = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

          ! Then update particle position (first part):  X_new = X_old + dt * 0.5 * V
          x_new = x_old + dtau * vi

     !     call self%kernel_smoother_1%add_current( x_old(1), x_new(1), wi(1), self%j1_dofs_local(:,1) )
          call self%kernel_smoother_1%add_current_split( x_old(1), x_new(1), &
               0, 1, wi(1)*vi(1)*dtau, &
               self%j1_dofs_local )
          call self%kernel_smoother_0%add_current_split( x_old(1), x_new(1), &
               0, 1, wi(1)*vi(2)*dtau, self%j2_dofs_local )
      !    call self%kernel_smoother_0%add_charge_int( x_old(1), x_new(1), wi(1)*vi(2)*dtau, &
      !         self%j2_dofs_local(:,1) )

      !    call self%kernel_smoother_0%add_charge( x_new(1), wi(1), self%rho_dofs_local )
          call self%kernel_smoother_0%add_charge_int ( x_old(1), x_new(1), wi(1), &
               self%rho_dofs_local )
          
          call self%particle_group%group(i_sp)%set_x(i_part, x_new)
       end do
    end do

    self%j1_dofs = 0.0_f64
    self%j2_dofs = 0.0_f64
   ! self%j1_dofs_local(:,1) = self%j1_dofs_local(:,1) + self%j1_dofs_local(:,2)
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j1_dofs_local(:,1), &
         n_cells, MPI_SUM, self%j1_dofs(:,1))
    call sll_o_collective_allreduce( sll_v_world_collective, self%j2_dofs_local(:,1), &
         n_cells, MPI_SUM, self%j2_dofs(:,1))

    self%rho_dofs = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, self%rho_dofs_local, &
         n_cells, MPI_SUM, self%rho_dofs )

    call self%maxwell_solver%compute_E_from_rho( self%rho_dofs, self%efield_dofs(:,1) )
    call self%maxwell_solver%compute_B_from_E( &
        dtau, self%efield_dofs(:,2), self%bfield_dofs )
    
    self%efield_dofs_old = self%efield_dofs
!!$    ! Initial guess for e_1 and b_2
!!$    self%bfield_dofs_even_older = self%bfield_dofs_old
!!$    self%bfield_dofs_old = self%bfield_dofs
!!$    call self%maxwell_solver%compute_B_from_E( &
!!$        dt, self%efield_dofs(:,2), self%bfield_dofs )

 

  end subroutine first_step_modified2
  
 ! Implementation of the first time step with dtau length and without push (can be changed to original version with full push)
  !---------------------------------------------------------------------------!
  subroutine first_step(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(inout) :: self !< time splitting object 
    sll_real64,                                             intent(in)    :: dt   !< time step

    !local variables
  !  sll_int32 :: i_part
    sll_real64 :: x_new(3), v_new(3)
    sll_int32  :: n_cells, i_sp, i_part
    sll_real64 :: residual_field
    sll_int32  :: n_iter_field

    self%step_one_mod = .false.
    self%step_one = .true.

    n_cells = self%kernel_smoother_0%n_dofs

    self%bfield_dofs_even_older = self%bfield_dofs

    self%efield_dofs_old = self%efield_dofs
    self%bfield_dofs_old = self%bfield_dofs
    call self%maxwell_solver%compute_B_from_E( &
         dt, self%efield_dofs(:,2), self%bfield_dofs )

    print*, 'Norm e', sum(self%efield_dofs(:,1)**2), self%maxwell_solver%L2norm_squared( self%efield_dofs(:,1), self%kernel_smoother_1%spline_degree )
    residual_field = 1.0_f64
    n_iter_field = 0
    do while( residual_field > self%tolerance_field .and. n_iter_field < self%n_max_iter_field )
       n_iter_field = n_iter_field + 1

       call particle_loop_subcycle_first2( self, dt )

       ! Update the electric field
       self%efield_filter = self%efield_dofs
       call self%maxwell_solver%compute_E_from_rho( self%rho_dofs, self%efield_dofs(:,1) )
       
       self%bfield_dofs = self%bfield_dofs_old
       call self%maxwell_solver%compute_B_from_E( &
            dt/real(self%n_sub_iter,f64), self%efield_dofs(:,2), self%bfield_dofs)
       
       residual_field = sum( (self%efield_dofs - self%efield_filter)**2)* self%delta_x
       write(32,*) self%efield_dofs(:,1)
       !     write(33,*) self%efield_filter(:,1)
       print*, 'Norm e', sum(self%efield_dofs(:,1)**2), self%maxwell_solver%L2norm_squared( self%efield_dofs(:,1), self%kernel_smoother_1%spline_degree )

    end do
    write(self%file_id_iter_field,*) n_iter_field, residual_field
 !   stop
    
    ! Temporary check of Gauss law for the case with substepping in the electric field
   ! self%rho_dofs = 0.0_f64
   ! call sll_o_collective_allreduce( sll_v_world_collective, self%rho_dofs_local, n_cells, &
   !      MPI_SUM, self%rho_dofs )
    call self%maxwell_solver%compute_rho_from_e( self%efield_dofs(:,1), self%rho_dofs_local )
    ! self%rho_dofs = self%rho_dofs/real(self%n_sub_iter,f64)
    write(72,*) self%rho_dofs
    self%rho_dofs = self%rho_dofs - sum(self%rho_dofs)/real(n_cells, f64)
    write(self%file_id_rho,*) maxval( abs(self%rho_dofs - self%rho_dofs_local ) )

    ! Set the new particle positions
    do i_sp = 1,self%n_species       
       do i_part=1,self%particle_group%group(i_sp)%n_particles
          x_new(1) = self%xvec(i_part, i_sp)
          v_new(1:2) = self%vvec(:,i_part, i_sp)
          call self%particle_group%group(i_sp)%set_x(i_part, x_new )
          call self%particle_group%group(i_sp)%set_v(i_part, v_new )
       end do
    end do    

  
  end subroutine first_step


  subroutine particle_loop_subcycle(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: x_new(3), vi(3), wi(1), x_old(3), wp(3), x_future(3), v_new(3), v_nnew(3)
    sll_int32  :: n_cells, i_sp
    sll_real64 :: qoverm
    sll_real64 :: efield(2),  efield_new(2), bfield_old(2), bfield_new(2)
    sll_real64 :: residual
    sll_int32  :: n_iter, n_iter_mean
    sll_int32  :: iter
    sll_real64 :: dtau
    sll_real64 :: gamma, residuum(2), gamma_e1, gamma_e2
    sll_real64 :: df11, df12, df21, df22, det
    sll_real64 :: kinetic_energy
    
    dtau = dt/real(self%n_sub_iter,f64)

    n_iter_mean = 0

    n_cells = self%kernel_smoother_0%n_dofs
    
    self%j1_dofs_local = 0.0_f64
    self%j2_dofs_local = 0.0_f64
    self%rho_dofs_local = 0.0_f64
  !  self%rho_dofs_local2 = 0.0_f64
       ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
       ! Then update particle position:  X_new = X_old + dt * V
       do i_sp = 1,self%n_species
          qoverm = self%particle_group%group(i_sp)%species%q_over_m();
          do i_part=1,self%particle_group%group(i_sp)%n_particles  
             ! Read out particle position and velocity
             x_new = self%particle_group%group(i_sp)%get_x(i_part)
             vi = self%particle_group%group(i_sp)%get_v(i_part)

             ! Get charge for accumulation of j
             wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

             ! Then update particle position (first part):  X_new = X_old + dt * 0.5 * V
             x_old = x_new - dtau * vi


             ! Now the first half of the velocity update
             ! First we evaluate the electric field
             if ( self%step_one_mod .eqv. .false. ) then
                call self%kernel_smoother_1%evaluate_int_quad(x_old(1), x_new(1), self%efield_dofs_old(:,1), efield_new(1))
                call self%kernel_smoother_0%evaluate_int_quad(x_old(1), x_new(1), self%efield_dofs_old(:,2), efield_new(2))
             else
                efield_new = 0.0_f64
             end if
             !  call self%kernel_smoother_1%evaluate &
             !       (x_new(1), self%efield_filter(:,1), efield(1))
             !  call self%kernel_smoother_0%evaluate &
             !       (x_new(1), self%efield_filter(:,2), efield(2))
             !  call self%kernel_smoother_1%evaluate_int_quad( x_old(1), x_new(1), self%bfield_dofs_old, bfield_old )
             if ( self%step_one_mod .eqv. .false. ) then
                call self%kernel_smoother_1%evaluate_int_linear_quad( x_old(1), x_new(1), self%n_sub_iter-1, self%n_sub_iter, self%bfield_dofs_even_older, self%bfield_dofs_old, bfield_old )
             else
                call self%kernel_smoother_1%evaluate_int_quad( x_old(1), x_new(1), self%bfield_dofs_old, bfield_old(1) )
              !  call self%kernel_smoother_1%evaluate_int_linear_quad( x_old(1), x_new(1), 0, 1, self%bfield_dofs_even_older, self%bfield_dofs_old, bfield_old )                
             end if
             ! Now we add up the fixed update into efield
             ! Note: Use dtau also for efield for the variant with orbit-averaged electric field
             efield(1) = qoverm * dtau * (  efield_new(1) +  vi(2) * bfield_old(1)  ) 
             efield(2) = qoverm * dtau * (  efield_new(2) -  vi(1) * bfield_old(1)  )

             v_new(1) = vi(1) + efield(1) + qoverm * dtau *  ( vi(2) * bfield_old(1) + efield_new(1) )
             v_new(2) = vi(2) + efield(2) + qoverm * dtau *  ( -vi(1) * bfield_old(1) + efield_new(2) )
             x_future = x_new + dtau * v_new

             call self%kernel_smoother_1%evaluate_int_linear_quad( x_new(1), x_future(1),  &
                  0, self%n_sub_iter, self%bfield_dofs_old, self%bfield_dofs, bfield_new )
             call self%kernel_smoother_1%evaluate_int_quad(x_future(1), x_new(1), &
                  self%efield_dofs(:,1), efield_new(1))
             call self%kernel_smoother_0%evaluate_int_quad(x_future(1), x_new(1), &
                  self%efield_dofs(:,2), efield_new(2))


             residuum(1) = v_new(1) - ( vi(1) + efield(1)  + &
                  qoverm * dtau * (v_new(2) * bfield_new(2) +  efield_new(1) ) )
             residuum(2) = v_new(2) - ( vi(2) + efield(2)  + &
                  qoverm * dtau * (- v_new(1) * bfield_new(2) + efield_new(2) ) )
             residual = 1.0_f64!maxval( abs(residuum) )

             n_iter = 0
             do while( residual > self%tolerance .and. n_iter < self%n_max_iter )
                n_iter = n_iter + 1 
                call self%kernel_smoother_1%evaluate_int_subnewton( x_new(1), x_future(1), &
                     self%efield_dofs(:,1), gamma_e1 )
                call self%kernel_smoother_0%evaluate_int_subnewton( x_new(1), x_future(1), &
                     self%efield_dofs(:,2), gamma_e2 )
                call self%kernel_smoother_1%evaluate_int_linear_quad_subnewton( x_new(1), x_future(1), &
                     0, self%n_sub_iter, self%bfield_dofs_old, self%bfield_dofs,gamma )
                gamma = gamma  *  qoverm
                
                df11 = 1.0_f64 - gamma * v_new(2) - gamma_e1 * qoverm
                df12 =  - qoverm * dtau * bfield_new(2)
                df21 =  ( qoverm * dtau * bfield_new(2) + gamma * v_new(1) - gamma_e2 * qoverm )
                df22 = 1.0_f64
                det = df11*df22- df12*df21

                v_nnew(1) = v_new(1) - ( df22 * residuum(1) - df12 * residuum(2) )/det
                v_nnew(2) = v_new(2) - ( - df21 * residuum(1) + df11 * residuum(2) )/det

                x_future = x_new + dtau * v_nnew
                v_new = v_nnew
                call self%kernel_smoother_1%evaluate_int_linear_quad( x_new(1), x_future(1),  &
                     0, self%n_sub_iter, self%bfield_dofs_old, self%bfield_dofs, bfield_new )
                call self%kernel_smoother_1%evaluate_int_quad(x_future(1), x_new(1), &
                     self%efield_dofs(:,1), efield_new(1))
                call self%kernel_smoother_0%evaluate_int_quad(x_future(1), x_new(1), &
                     self%efield_dofs(:,2), efield_new(2))

                residuum(1) = v_new(1) - ( vi(1) + efield(1) + qoverm * dtau * (v_new(2) * bfield_new(2) + efield_new(1) ) )
                residuum(2) = v_new(2) - ( vi(2) + efield(2) + qoverm * dtau * ( - v_new(1) * bfield_new(2) + efield_new(2) ) )
                residual = maxval( abs(residuum) )
             end do
             v_nnew(1) = vi(1) + efield(1) + qoverm * dtau * ( v_new(2) * bfield_new(2) + efield_new(1) ) 
             v_nnew(2) = vi(2) + efield(2) + qoverm * dtau * (-v_new(1) * bfield_new(2) + efield_new(2) )
             v_new = v_nnew
             x_future = x_new + dtau * v_nnew

             if ( n_iter .eq. self%n_max_iter ) then
                print*, 'First iteration did not converge. Residual:', residual
             end if
             n_iter_mean = n_iter_mean + n_iter
             ! End of the loop

             if (abs (v_new(1))> 1E-16_f64) then
                ! Now the charge deposition
                call self%kernel_smoother_1%add_current_split( x_new(1), x_future(1), &
                     0, self%n_sub_iter, wi(1)*v_new(1)*dtau, &
                     self%j1_dofs_local )

                call self%kernel_smoother_0%add_current_split( x_new(1), x_future(1), &
                     0, self%n_sub_iter, wi(1)*v_new(2)*dtau, self%j2_dofs_local )
                ! For rho check
                !call self%kernel_smoother_0%add_current( x_new, x_future, &
                !     wi(1)/v_new(1)/dtau, self%rho_dofs_local )
                !   else
              !  self%rho_dofs = 0.0_f64
              !  call self%kernel_smoother_0%add_charge( x_new, wi(1), self%rho_dofs )
             end if
             call self%kernel_smoother_0%add_charge_int( x_new(1), x_future(1), wi(1), &
                  self%rho_dofs_local )
             
      !       call self%kernel_smoother_0%add_current_split( x_new(1), x_future(1), &
      !            0,self%n_sub_iter, wi(1), self%rho_dofs_local2 )
         !    print*, self%rho_dofs_local
         !    print*, self%rho_dofs_local2(:,1)+ self%rho_dofs_local2(:,2)
         !    stop
            ! print*, 'charge', self%rho_dofs*dtau
            ! print*, 'current', self%rho_dofs_local
            ! stop

             do iter=2,self%n_sub_iter
                ! Then update particle position (second part):  X_new = X_old + dt * 0.5 * V
                x_old = x_new
                x_new = x_future
                vi = v_new

                ! Then the second half of the velocity update
                ! First the fixed old part again
                bfield_old(1) = bfield_new(1)
                call self%kernel_smoother_1%evaluate_int_quad(x_old(1), x_new(1), &
                     self%efield_dofs(:,1), efield_new(1))
                call self%kernel_smoother_0%evaluate_int_quad(x_old(1), x_new(1), &
                     self%efield_dofs(:,2), efield_new(2))

                ! Now we add up the fixed update into efield
                efield(1) = qoverm * dtau * (  vi(2) * bfield_old(1) + efield_new(1) )
                efield(2) = qoverm * dtau * (- vi(1) * bfield_old(1) + efield_new(2) )

                ! First guess for v_new = v_old (already set)
                v_new(1) = vi(1) + efield(1) + qoverm * dtau * ( vi(2) * bfield_old(1) + efield_new(1) ) 
                v_new(2) = vi(2) + efield(2) + qoverm * dtau * ( -vi(1) * bfield_old(1) + efield_new(2) )
                x_future = x_new + dtau * v_new

                ! Loop here for convergence
                ! call self%kernel_smoother_1%evaluate_int_quad( x_future(1), x_new(1), self%bfield_filter, bfield_new )
                call self%kernel_smoother_1%evaluate_int_linear_quad( x_new(1), x_future(1), &
                     iter-1, self%n_sub_iter, self%bfield_dofs_old, self%bfield_dofs, bfield_new )
                call self%kernel_smoother_1%evaluate_int_quad(x_future(1), x_new(1), &
                     self%efield_dofs(:,1), efield_new(1))
                call self%kernel_smoother_0%evaluate_int_quad(x_future(1), x_new(1), &
                     self%efield_dofs(:,2), efield_new(2))

                residuum(1) = v_new(1) - ( vi(1) + efield(1) + &
                     qoverm * dtau * (v_new(2) * bfield_new(2) + efield_new(1) ) )
                residuum(2) = v_new(2) - ( vi(2) + efield(2) +&
                     qoverm * dtau * (-v_new(1) * bfield_new(2)+ efield_new(2) ) )
                residual = 1.0_f64!maxval( abs(residuum) )
                n_iter = 0
                do while( residual > self%tolerance .and. n_iter < self%n_max_iter )
                   n_iter = n_iter + 1 
                   call self%kernel_smoother_1%evaluate_int_subnewton( x_new(1), x_future(1), &
                        self%efield_dofs(:,1), gamma_e1 )
                   call self%kernel_smoother_0%evaluate_int_subnewton( x_new(1), x_future(1), &
                        self%efield_dofs(:,2), gamma_e2 )
                   call self%kernel_smoother_1%evaluate_int_linear_quad_subnewton( x_new(1), x_future(1), &
                        iter-1, self%n_sub_iter, self%bfield_dofs_old, self%bfield_dofs, gamma )
                   gamma = gamma *  qoverm

                   df11 = 1.0_f64 - gamma * v_new(2) - gamma_e1 * qoverm
                   df12 =  - qoverm * dtau * bfield_new(2)
                   df21 =  ( qoverm * dtau * bfield_new(2) + gamma * v_new(1) - gamma_e2 * qoverm )
                   df22 = 1.0_f64
                   det = df11*df22- df12*df21

                   v_nnew(1) = v_new(1) - ( df22 * residuum(1) - df12 * residuum(2) )/det
                   v_nnew(2) = v_new(2) - ( - df21 * residuum(1) + df11 * residuum(2) )/det

                   x_future = x_new + dtau * v_nnew
                   v_new = v_nnew
                   !   call self%kernel_smoother_1%evaluate_int_quad( x_future(1), x_new(1), self%bfield_filter, bfield_new )
                   call self%kernel_smoother_1%evaluate_int_linear_quad( x_new(1), x_future(1), &
                        iter-1, self%n_sub_iter, self%bfield_dofs_old, self%bfield_dofs, bfield_new )
                   call self%kernel_smoother_1%evaluate_int_quad(x_future(1), x_new(1), &
                        self%efield_dofs(:,1), efield_new(1))
                   call self%kernel_smoother_0%evaluate_int_quad(x_future(1), x_new(1), &
                        self%efield_dofs(:,2), efield_new(2))

                   residuum(1) = v_new(1) - ( vi(1) + efield(1) + qoverm * dtau * (  v_new(2) * bfield_new(2) + efield_new(1) ) )
                   residuum(2) = v_new(2) - ( vi(2) + efield(2) + qoverm * dtau * ( -v_new(1) * bfield_new(2) + efield_new(2) ) )
                   residual = maxval( abs(residuum) )
                end do
                v_nnew(1) = vi(1) + efield(1) + qoverm * dtau * (  v_new(2) * bfield_new(2) + efield_new(1) ) 
                v_nnew(2) = vi(2) + efield(2) + qoverm * dtau * ( -v_new(1) * bfield_new(2) + efield_new(2) ) 
                v_new = v_nnew
                x_future = x_new + dtau * v_nnew
                if ( n_iter .eq. self%n_max_iter ) then
                   print*, 'Second iteration did not converge. Residual:', residual, self%tolerance
                end if
                n_iter_mean = n_iter_mean + n_iter
                ! End of the loop

                if (abs (v_new(1))> 1E-16_f64) then
                   ! Now the charge deposition
                   call self%kernel_smoother_1%add_current_split( x_new(1), x_future(1), &
                        iter-1,self%n_sub_iter, wi(1)*v_new(1)*dtau, self%j1_dofs_local )
                   call self%kernel_smoother_0%add_current_split( x_new(1), x_future(1), &
                        iter-1,self%n_sub_iter, wi(1)*v_new(2)*dtau, self%j2_dofs_local )
                   ! For rho check                  
                !   call self%kernel_smoother_0%add_current( x_new, x_future, &
                !        wi(1)/v_new(1)/dtau, self%rho_dofs_local )
               ! else
              !     call self%kernel_smoother_0%add_charge( x_new, wi(1), self%rho_dofs_local )
                end if
                call self%kernel_smoother_0%add_charge_int( x_new(1), x_future(1), wi(1), &
                     self%rho_dofs_local )
       !         call self%kernel_smoother_0%add_current_split( x_new(1), x_future(1), &
       !              iter-1,self%n_sub_iter, wi(1), self%rho_dofs_local2 )
             end do

             x_new(1) = modulo(x_future(1), self%Lx)
             self%vvec(:,i_part,i_sp) = v_new(1:2)
             self%xvec(i_part, i_sp) = x_new(1)
             ! call self%particle_group%group(i_sp)%set_x(i_part, x_new)
             ! call self%particle_group%group(i_sp)%set_v(i_part, v_new)

             if (self%particle_group%group(i_sp)%n_weights == 3) then
                ! Update weights if control variate
                wp = self%particle_group%group(i_sp)%get_weights(i_part)          
                wp(3) = self%control_variate%cv(i_sp)%update_df_weight( x_new(1:1), vi(1:2), 0.0_f64, wp(1), wp(2))
                call self%particle_group%group(i_sp)%set_weights(i_part, wp)
             end if
            ! x_old = self%particle_group%group(i_sp)%get_x(i_part)
            ! kinetic_energy = kinetic_energy + ( (x_future(1)- x_old(1))/dt )**2

          end do
       end do

       write(self%file_id,*) real(n_iter_mean, f64)/real(self%particle_group%group(1)%n_particles,f64)/real(self%n_sub_iter,f64)

       ! Temporary check of Gauss law for the case with substepping in the electric field
       ! self%j_dofs = 0.0_f64
       ! call sll_o_collective_allreduce( sll_v_world_collective, self%rho_dofs_local, n_cells, &
       !      MPI_SUM, self%j_dofs(:,1) )
       ! call self%maxwell_solver%compute_rho_from_e( self%efield_dofs(:,1), self%j_dofs(:,2) )

       ! self%j_dofs(:,1) = self%j_dofs(:,1)/real(self%n_sub_iter,f64)
       ! self%j_dofs(:,1) = self%j_dofs(:,1) - sum(self%j_dofs(:,1))/real(n_cells, f64)
       ! write(self%file_id_rho,*) maxval( abs(self%j_dofs(:,1) - self%j_dofs(:,2) ) )
       !    write(16,*) self%j1_dofs_local(:,1)
       !    write(17,*) self%j2_dofs_local(:,1)
       !    stop

       self%j1_dofs = 0.0_f64
       self%j2_dofs = 0.0_f64
       ! MPI to sum up contributions from each processor
       call sll_o_collective_allreduce( sll_v_world_collective, self%j1_dofs_local(:,1), &
            n_cells, MPI_SUM, self%j1_dofs(:,1))
       call sll_o_collective_allreduce( sll_v_world_collective, self%j1_dofs_local(:,2), &
            n_cells, MPI_SUM, self%j1_dofs(:,2))
       call sll_o_collective_allreduce( sll_v_world_collective, self%j2_dofs_local(:,1), &
            n_cells, MPI_SUM, self%j2_dofs(:,1))
       call sll_o_collective_allreduce( sll_v_world_collective, self%j2_dofs_local(:,2), &
            n_cells, MPI_SUM, self%j2_dofs(:,2))

  end subroutine particle_loop_subcycle


  subroutine particle_loop_subcycle_first2(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    
    !local variables
    sll_int32 :: i_part
    sll_real64 :: x_new(3), vi(3), wi(1), x_old(3), wp(3), x_future(3), v_new(3), v_nnew(3)
    sll_int32  :: n_cells, i_sp
    sll_real64 :: qoverm
    sll_real64 :: efield(2),  efield_new(2), bfield_old(2), bfield_new(2)
    sll_real64 :: residual
    sll_int32  :: n_iter, n_iter_mean
    sll_int32  :: iter
    sll_real64 :: dtau
    sll_real64 :: gamma, residuum(2), gamma_e1, gamma_e2
    sll_real64 :: df11, df12, df21, df22, det

    
    dtau = dt/real(self%n_sub_iter,f64)

    n_iter_mean = 0

    n_cells = self%kernel_smoother_0%n_dofs
    
    self%j1_dofs_local = 0.0_f64
    self%j2_dofs_local = 0.0_f64
    self%rho_dofs_local = 0.0_f64
    !   self%rho_dofs_local2 = 0.0_f64

    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
       ! Then update particle position:  X_new = X_old + dt * V
       do i_sp = 1,self%n_species
          qoverm = self%particle_group%group(i_sp)%species%q_over_m();
          do i_part=1,self%particle_group%group(i_sp)%n_particles  
             ! Read out particle position and velocity
             x_new = self%particle_group%group(i_sp)%get_x(i_part)
             v_new = self%particle_group%group(i_sp)%get_v(i_part)

             ! Get charge for accumulation of j
             wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

             ! Then update particle position (first part):  X_new = X_old + dt * 0.5 * V
             x_future = x_new + dtau * v_new

             if (abs (v_new(1))> 1E-16_f64) then
                ! Now the charge deposition
                call self%kernel_smoother_1%add_current_split( x_new(1), x_future(1), &
                     0, 1, wi(1)*v_new(1)*dtau, &
                     self%j1_dofs_local )

                call self%kernel_smoother_0%add_current_split( x_new(1), x_future(1), &
                     0, 1, wi(1)*v_new(2)*dtau, self%j2_dofs_local )
                ! For rho check
                !call self%kernel_smoother_0%add_current( x_new, x_future, &
                !     wi(1)/v_new(1)/dtau, self%rho_dofs_local )
                !   else
              !  self%rho_dofs = 0.0_f64
              !  call self%kernel_smoother_0%add_charge( x_new, wi(1), self%rho_dofs )
             end if
             call self%kernel_smoother_0%add_charge_int( x_new(1), x_future(1), wi(1), &
                  self%rho_dofs_local )
        !     call self%kernel_smoother_0%add_current_split( x_new(1), x_future(1), &
        !          iter-1,self%n_sub_iter, wi(1), self%rho_dofs_local2 )

             

             x_new(1) = modulo(x_future(1), self%Lx)
             self%vvec(:,i_part,i_sp) = v_new(1:2)
             self%xvec(i_part, i_sp) = x_new(1)
             ! call self%particle_group%group(i_sp)%set_x(i_part, x_new)
             ! call self%particle_group%group(i_sp)%set_v(i_part, v_new)

             if (self%particle_group%group(i_sp)%n_weights == 3) then
                ! Update weights if control variate
                wp = self%particle_group%group(i_sp)%get_weights(i_part)          
                wp(3) = self%control_variate%cv(i_sp)%update_df_weight( x_new(1:1), vi(1:2), 0.0_f64, wp(1), wp(2))
                call self%particle_group%group(i_sp)%set_weights(i_part, wp)
             end if

          end do
       end do

       write(self%file_id,*) real(n_iter_mean, f64)/real(self%particle_group%group(1)%n_particles,f64)/real(self%n_sub_iter,f64)


       self%j1_dofs = 0.0_f64
       self%j2_dofs = 0.0_f64
       ! MPI to sum up contributions from each processor
       call sll_o_collective_allreduce( sll_v_world_collective, self%j1_dofs_local(:,1), &
            n_cells, MPI_SUM, self%j1_dofs(:,1))
       call sll_o_collective_allreduce( sll_v_world_collective, self%j1_dofs_local(:,2), &
            n_cells, MPI_SUM, self%j1_dofs(:,2))
       call sll_o_collective_allreduce( sll_v_world_collective, self%j2_dofs_local(:,1), &
            n_cells, MPI_SUM, self%j2_dofs(:,1))
       call sll_o_collective_allreduce( sll_v_world_collective, self%j2_dofs_local(:,2), &
            n_cells, MPI_SUM, self%j2_dofs(:,2))

       self%rho_dofs = 0.0_f64
       call sll_o_collective_allreduce( sll_v_world_collective, self%rho_dofs_local, n_cells, &
            MPI_SUM, self%rho_dofs )
       self%rho_dofs = self%rho_dofs!/real(self%n_sub_iter,f64)

     end subroutine particle_loop_subcycle_first2

  subroutine particle_loop_subcycle_first(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    
    !local variables
    sll_int32 :: i_part
    sll_real64 :: x_new(3), vi(3), wi(1), x_old(3), wp(3), x_future(3), v_new(3), v_nnew(3)
    sll_int32  :: n_cells, i_sp
    sll_real64 :: qoverm
    sll_real64 :: efield(2),  efield_new(2), bfield_old(2), bfield_new(2)
    sll_real64 :: residual
    sll_int32  :: n_iter, n_iter_mean
    sll_int32  :: iter
    sll_real64 :: dtau
    sll_real64 :: gamma, residuum(2), gamma_e1, gamma_e2
    sll_real64 :: df11, df12, df21, df22, det

    
    dtau = dt/real(self%n_sub_iter,f64)

    n_iter_mean = 0

    n_cells = self%kernel_smoother_0%n_dofs
    
    self%j1_dofs_local = 0.0_f64
    self%j2_dofs_local = 0.0_f64
    self%rho_dofs_local = 0.0_f64
    !   self%rho_dofs_local2 = 0.0_f64

    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
       ! Then update particle position:  X_new = X_old + dt * V
       do i_sp = 1,self%n_species
          qoverm = self%particle_group%group(i_sp)%species%q_over_m();
          do i_part=1,self%particle_group%group(i_sp)%n_particles  
             ! Read out particle position and velocity
             x_new = self%particle_group%group(i_sp)%get_x(i_part)
             v_new = self%particle_group%group(i_sp)%get_v(i_part)

             ! Get charge for accumulation of j
             wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

             ! Then update particle position (first part):  X_new = X_old + dt * 0.5 * V
             x_future = x_new + dtau * v_new

             if (abs (v_new(1))> 1E-16_f64) then
                ! Now the charge deposition
                call self%kernel_smoother_1%add_current_split( x_new(1), x_future(1), &
                     0, self%n_sub_iter, wi(1)*v_new(1)*dtau, &
                     self%j1_dofs_local )

                call self%kernel_smoother_0%add_current_split( x_new(1), x_future(1), &
                     0, self%n_sub_iter, wi(1)*v_new(2)*dtau, self%j2_dofs_local )
                ! For rho check
                !call self%kernel_smoother_0%add_current( x_new, x_future, &
                !     wi(1)/v_new(1)/dtau, self%rho_dofs_local )
                !   else
              !  self%rho_dofs = 0.0_f64
              !  call self%kernel_smoother_0%add_charge( x_new, wi(1), self%rho_dofs )
             end if
             call self%kernel_smoother_0%add_charge_int( x_new(1), x_future(1), wi(1), &
                  self%rho_dofs_local )
        !     call self%kernel_smoother_0%add_current_split( x_new(1), x_future(1), &
        !          iter-1,self%n_sub_iter, wi(1), self%rho_dofs_local2 )

             
             call self%kernel_smoother_1%evaluate_int_linear_quad( x_new(1), x_future(1),  &
                  0, self%n_sub_iter, self%bfield_dofs_old, self%bfield_dofs, bfield_new )
             
             do iter=2,self%n_sub_iter
                ! Then update particle position (second part):  X_new = X_old + dt * 0.5 * V
                x_old = x_new
                x_new = x_future
                vi = v_new

                ! Then the second half of the velocity update
                ! First the fixed old part again
                bfield_old(1) = bfield_new(1)
                call self%kernel_smoother_1%evaluate_int_quad(x_old(1), x_new(1), &
                     self%efield_dofs(:,1), efield_new(1))
                call self%kernel_smoother_0%evaluate_int_quad(x_old(1), x_new(1), &
                     self%efield_dofs(:,2), efield_new(2))

                ! Now we add up the fixed update into efield
                efield(1) = qoverm * dtau * (  vi(2) * bfield_old(1) + efield_new(1) )
                efield(2) = qoverm * dtau * (- vi(1) * bfield_old(1) + efield_new(2) )

                ! First guess for v_new = v_old (already set)
                v_new(1) = vi(1) + efield(1) + qoverm * dtau * (  vi(2) * bfield_old(1) + efield_new(1) ) 
                v_new(2) = vi(2) + efield(2) + qoverm * dtau * ( -vi(1) * bfield_old(1) + efield_new(2) )
                x_future = x_new + dtau * v_new

                ! Loop here for convergence
                ! call self%kernel_smoother_1%evaluate_int_quad( x_future(1), x_new(1), self%bfield_filter, bfield_new )
                call self%kernel_smoother_1%evaluate_int_linear_quad( x_new(1), x_future(1), &
                     iter-1, self%n_sub_iter, self%bfield_dofs_old, self%bfield_dofs, bfield_new )
                call self%kernel_smoother_1%evaluate_int_quad(x_future(1), x_new(1), &
                     self%efield_dofs(:,1), efield_new(1))
                call self%kernel_smoother_0%evaluate_int_quad(x_future(1), x_new(1), &
                     self%efield_dofs(:,2), efield_new(2))


                residuum(1) = v_new(1) - ( vi(1) + efield(1) + &
                     qoverm * dtau * (v_new(2) * bfield_new(2) + efield_new(1) ) )
                residuum(2) = v_new(2) - ( vi(2) + efield(2) +&
                     qoverm * dtau * (-v_new(1) * bfield_new(2)+ efield_new(2) ) )
                residual = 1.0_f64!maxval( abs(residuum) )
                n_iter = 0
                do while( residual > self%tolerance .and. n_iter < self%n_max_iter )
                   n_iter = n_iter + 1 
                   call self%kernel_smoother_1%evaluate_int_subnewton( x_new(1), x_future(1), &
                        self%efield_dofs(:,1), gamma_e1 )
                   call self%kernel_smoother_0%evaluate_int_subnewton( x_new(1), x_future(1), &
                        self%efield_dofs(:,2), gamma_e2 )
                   call self%kernel_smoother_1%evaluate_int_linear_quad_subnewton( x_new(1), x_future(1), &
                        iter-1, self%n_sub_iter, self%bfield_dofs_old, self%bfield_dofs, gamma )
                   gamma = gamma *  qoverm

                   df11 = 1.0_f64 - gamma * v_new(2) - gamma_e1 * qoverm
                   df12 =  - qoverm * dtau * bfield_new(2)
                   df21 =  ( qoverm * dtau * bfield_new(2) + gamma * v_new(1) - gamma_e2 * qoverm )
                   df22 = 1.0_f64
                   det = df11*df22- df12*df21

                   v_nnew(1) = v_new(1) - ( df22 * residuum(1) - df12 * residuum(2) )/det
                   v_nnew(2) = v_new(2) - ( - df21 * residuum(1) + df11 * residuum(2) )/det

                   x_future = x_new + dtau * v_nnew
                   v_new = v_nnew
                   !   call self%kernel_smoother_1%evaluate_int_quad( x_future(1), x_new(1), self%bfield_filter, bfield_new )
                   call self%kernel_smoother_1%evaluate_int_linear_quad( x_new(1), x_future(1), &
                        iter-1, self%n_sub_iter, self%bfield_dofs_old, self%bfield_dofs, bfield_new )
                   call self%kernel_smoother_1%evaluate_int_quad(x_future(1), x_new(1), &
                        self%efield_dofs(:,1), efield_new(1))
                   call self%kernel_smoother_0%evaluate_int_quad(x_future(1), x_new(1), &
                        self%efield_dofs(:,2), efield_new(2))

                   residuum(1) = v_new(1) - ( vi(1) + efield(1) + qoverm * dtau * (  v_new(2) * bfield_new(2) + efield_new(1) ) )
                   residuum(2) = v_new(2) - ( vi(2) + efield(2) + qoverm * dtau * ( -v_new(1) * bfield_new(2) + efield_new(2) ) )
                   residual = maxval( abs(residuum) )
                end do
                v_nnew(1) = vi(1) + efield(1) + qoverm * dtau * (  v_new(2) * bfield_new(2) + efield_new(1) ) 
                v_nnew(2) = vi(2) + efield(2) + qoverm * dtau * ( -v_new(1) * bfield_new(2) + efield_new(2) ) 
                v_new = v_nnew
                x_future = x_new + dtau * v_nnew
                if ( n_iter .eq. self%n_max_iter ) then
                   print*, 'Second iteration did not converge. Residual:', residual, self%tolerance
                end if
                n_iter_mean = n_iter_mean + n_iter
                ! End of the loop

                if (abs (v_new(1))> 1E-16_f64) then
                   ! Now the charge deposition
                   call self%kernel_smoother_1%add_current_split( x_new(1), x_future(1), &
                        iter-1,self%n_sub_iter, wi(1)*v_new(1)*dtau, self%j1_dofs_local )
                   call self%kernel_smoother_0%add_current_split( x_new(1), x_future(1), &
                        iter-1,self%n_sub_iter, wi(1)*v_new(2)*dtau, self%j2_dofs_local )
                   ! For rho check                  
                !   call self%kernel_smoother_0%add_current( x_new, x_future, &
                !        wi(1)/v_new(1)/dtau, self%rho_dofs_local )
               ! else
              !     call self%kernel_smoother_0%add_charge( x_new, wi(1), self%rho_dofs_local )
                end if
                call self%kernel_smoother_0%add_charge_int( x_new(1), x_future(1), wi(1), &
                     self%rho_dofs_local )
         !       call self%kernel_smoother_0%add_current_split( x_new(1), x_future(1), &
         !            iter-1,self%n_sub_iter, wi(1), self%rho_dofs_local2 )
             end do

             x_new(1) = modulo(x_future(1), self%Lx)
             self%vvec(:,i_part,i_sp) = v_new(1:2)
             self%xvec(i_part, i_sp) = x_new(1)
             ! call self%particle_group%group(i_sp)%set_x(i_part, x_new)
             ! call self%particle_group%group(i_sp)%set_v(i_part, v_new)

             if (self%particle_group%group(i_sp)%n_weights == 3) then
                ! Update weights if control variate
                wp = self%particle_group%group(i_sp)%get_weights(i_part)          
                wp(3) = self%control_variate%cv(i_sp)%update_df_weight( x_new(1:1), vi(1:2), 0.0_f64, wp(1), wp(2))
                call self%particle_group%group(i_sp)%set_weights(i_part, wp)
             end if

          end do
       end do

       write(self%file_id,*) real(n_iter_mean, f64)/real(self%particle_group%group(1)%n_particles,f64)/real(self%n_sub_iter,f64)

       ! Temporary check of Gauss law for the case with substepping in the electric field
       ! self%j_dofs = 0.0_f64
       ! call sll_o_collective_allreduce( sll_v_world_collective, self%rho_dofs_local, n_cells, &
       !      MPI_SUM, self%j_dofs(:,1) )
       ! call self%maxwell_solver%compute_rho_from_e( self%efield_dofs(:,1), self%j_dofs(:,2) )

       ! self%j_dofs(:,1) = self%j_dofs(:,1)/real(self%n_sub_iter,f64)
       ! self%j_dofs(:,1) = self%j_dofs(:,1) - sum(self%j_dofs(:,1))/real(n_cells, f64)
       ! write(self%file_id_rho,*) maxval( abs(self%j_dofs(:,1) - self%j_dofs(:,2) ) )
       !    write(16,*) self%j1_dofs_local(:,1)
       !    write(17,*) self%j2_dofs_local(:,1)
       !    stop

       self%j1_dofs = 0.0_f64
       self%j2_dofs = 0.0_f64
       ! MPI to sum up contributions from each processor
       call sll_o_collective_allreduce( sll_v_world_collective, self%j1_dofs_local(:,1), &
            n_cells, MPI_SUM, self%j1_dofs(:,1))
       call sll_o_collective_allreduce( sll_v_world_collective, self%j1_dofs_local(:,2), &
            n_cells, MPI_SUM, self%j1_dofs(:,2))
       call sll_o_collective_allreduce( sll_v_world_collective, self%j2_dofs_local(:,1), &
            n_cells, MPI_SUM, self%j2_dofs(:,1))
       call sll_o_collective_allreduce( sll_v_world_collective, self%j2_dofs_local(:,2), &
            n_cells, MPI_SUM, self%j2_dofs(:,2))

       self%rho_dofs = 0.0_f64
       call sll_o_collective_allreduce( sll_v_world_collective, self%rho_dofs_local, n_cells, &
            MPI_SUM, self%rho_dofs )
       self%rho_dofs = self%rho_dofs/real(self%n_sub_iter,f64)

     end subroutine particle_loop_subcycle_first

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
       file_prefix, &
       filter, &
       jmean, &
       control_variate, &
       i_weight) 
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent(out) :: self !< time splitting object 
    class(sll_c_maxwell_1d_base), target,          intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,           intent(in)  :: particle_group
    sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
    sll_int32,                                      intent(in) :: n_sub_iter !< Number of subcycling steps
    character(*),                                   intent(in) :: file_prefix !< File prefix
    type( sll_t_binomial_filter ), intent( in ), target :: filter
    logical, optional, intent(in) :: jmean
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    sll_int32, optional,                            intent(in) :: i_weight !< Index of weight to be used by propagator

    !local variables
    sll_int32 :: ierr, i_sp
    sll_int32 :: n_particles_max

    print*, 'Subcycling propagator with efield'

    self%maxwell_solver => maxwell_solver

    ! self%kernel_smoother_0 => kernel_smoother_0
    ! self%kernel_smoother_1 => kernel_smoother_1
    select type ( ks0=>kernel_smoother_0)
    type is ( sll_t_particle_mesh_coupling_spline_1d )
       self%kernel_smoother_0 => ks0
    end select
    select type (ks1=>kernel_smoother_1)
    type is ( sll_t_particle_mesh_coupling_spline_1d )
       self%kernel_smoother_1 => ks1
    end select

    self%n_species = particle_group%n_species
    !allocate( sll_t_time_propagator_pic_vm_1d2v_subcyci :: self%particle_group(self%n_species) )
    !do j=1,self%n_species
    !   self%particle_group(j) => particle_group(j)
    !end do

    self%particle_group => particle_group
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs

    self%n_sub_iter = n_sub_iter

    ! Check that n_dofs is the same for both kernel smoothers.
    SLL_ASSERT( self%kernel_smoother_0%n_dofs == self%kernel_smoother_1%n_dofs )

    SLL_ALLOCATE(self%j1_dofs(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%j1_dofs_local(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%j2_dofs(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%j2_dofs_local(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%rho_dofs_local(self%kernel_smoother_0%n_dofs), ierr)
    SLL_ALLOCATE(self%rho_dofs(self%kernel_smoother_0%n_dofs), ierr)
    SLL_ALLOCATE(self%efield_filter(self%kernel_smoother_1%n_dofs,2), ierr)
    SLL_ALLOCATE(self%bfield_filter(self%kernel_smoother_1%n_dofs), ierr)
    SLL_ALLOCATE(self%bfield_dofs_old(self%kernel_smoother_1%n_dofs), ierr)
    SLL_ALLOCATE(self%bfield_dofs_even_older(self%kernel_smoother_1%n_dofs), ierr)
    SLL_ALLOCATE(self%bfield_linear_combination(self%kernel_smoother_1%n_dofs), ierr)
    SLL_ALLOCATE(self%efield_dofs_old(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%j1_dofs_old(self%kernel_smoother_0%n_dofs), ierr)
    SLL_ALLOCATE(self%j2_dofs_old(self%kernel_smoother_0%n_dofs), ierr)
    SLL_ALLOCATE(self%efield_tmp(self%kernel_smoother_0%n_dofs), ierr)
    
 !   SLL_ALLOCATE(self%rho_dofs_local2(self%kernel_smoother_0%n_dofs,2), ierr)

    n_particles_max = self%particle_group%group(1)%n_particles
    do i_sp = 2,self%particle_group%n_species       
       n_particles_max = max(n_particles_max, self%particle_group%group(i_sp)%n_particles )
    end do

    SLL_ALLOCATE( self%xvec(n_particles_max, self%particle_group%n_species), ierr )
    SLL_ALLOCATE( self%vvec(2,n_particles_max,self%particle_group%n_species), ierr )

    self%spline_degree = self%kernel_smoother_0%spline_degree
    self%x_min = x_min
    self%Lx = Lx
    self%delta_x = self%Lx/real(self%kernel_smoother_1%n_dofs, f64)

    self%cell_integrals_1 = [0.5_f64, 2.0_f64, 0.5_f64]
    self%cell_integrals_1 = self%cell_integrals_1 / 3.0_f64

    self%cell_integrals_0 = [1.0_f64,11.0_f64,11.0_f64,1.0_f64]
    self%cell_integrals_0 = self%cell_integrals_0 / 24.0_f64

    self%filter => filter

    call self%filter%apply( self%efield_dofs(:,1), self%efield_filter(:,1) ) 
    call self%filter%apply( self%efield_dofs(:,2), self%efield_filter(:,2) ) 
    call self%filter%apply( self%bfield_dofs, self%bfield_filter ) 

    if (present(jmean)) then
       self%jmean = jmean
    end if

    self%i_weight = 1
    if (present(i_weight)) self%i_weight = i_weight
    if(present(control_variate)) then
       allocate(self%control_variate )
       allocate(self%control_variate%cv(self%n_species) )
       self%control_variate => control_variate
       !do j=1,self%n_species
       !   self%control_variate%cv(j) => control_variate%cv(j)
       !end do
    end if

    ! File to write out no. of iterations
    open(newunit=self%file_id, file=trim(file_prefix)//"_subcycling_iterations.dat",action='write')

    open(newunit=self%file_id_rho, file=trim(file_prefix)//"_rho.dat",action='write')
    
    open(newunit=self%file_id_iter_field, file=trim(file_prefix)//"_subcycling_outer_iterations.dat",action='write')


  end subroutine initialize_pic_vm_1d2v

  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_1d2v(self)
    class(sll_t_time_propagator_pic_vm_1d2v_subcyci), intent( inout ) :: self !< time splitting object 

    deallocate(self%j1_dofs)
    deallocate(self%j1_dofs_local)
    deallocate(self%j2_dofs)
    deallocate(self%j2_dofs_local)
    self%maxwell_solver => null()
    self%kernel_smoother_0 => null()
    self%kernel_smoother_1 => null()
    self%particle_group => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()

  end subroutine delete_pic_vm_1d2v


  !---------------------------------------------------------------------------!
  !> Constructor for allocatable abstract type.
  subroutine sll_s_new_time_propagator_pic_vm_1d2v_subcyci(&
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
       file_prefix, &
       filter, &
       jmean, &
       control_variate, &
       i_weight) 
    class(sll_c_time_propagator_base), allocatable, intent(out) :: splitting !< time splitting object 
    class(sll_c_maxwell_1d_base), target,                intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,                intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,                intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,           intent(in)  :: particle_group
    !class(sll_c_particle_group_base),target,             intent(in)  :: particle_group(:) !< Particle group
    sll_real64, target,                                  intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                                  intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                           intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                           intent(in)  :: Lx !< Length of the domain in x direction.
    sll_int32,                                      intent(in) :: n_sub_iter !< Number of subcycling steps
    character(*),                                   intent(in) :: file_prefix !< File prefix
    type( sll_t_binomial_filter ), intent( in ), target :: filter
    logical, optional, intent(in) :: jmean !< Should jmean be substracted in Ampere's law?
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    sll_int32, optional,                            intent(in) :: i_weight !< Index of weight to be used by propagator

    !local variables
    sll_int32 :: ierr
    logical :: jmean_val 

    SLL_ALLOCATE(sll_t_time_propagator_pic_vm_1d2v_subcyci :: splitting, ierr)

    if (present(jmean) ) then
       jmean_val = jmean
    else
       jmean_val = .false.
    end if

    select type (splitting)
    type is ( sll_t_time_propagator_pic_vm_1d2v_subcyci )
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
               file_prefix, &
               filter, &
               jmean_val, &
               control_variate, &
               i_weight)
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
               file_prefix, &
               filter, &
               jmean_val)
       end if
    end select

  end subroutine sll_s_new_time_propagator_pic_vm_1d2v_subcyci

  !---------------------------------------------------------------------------!
  !> Constructor for pointer abstract type.
  subroutine sll_s_new_time_propagator_pic_vm_1d2v_subcyci_ptr(&
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
       file_prefix, &
       filter, &
       jmean) 
    class(sll_c_time_propagator_base), pointer, intent(out) :: splitting !< time splitting object 
    class(sll_c_maxwell_1d_base), target,            intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,            intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,            intent(in)  :: kernel_smoother_1  !< Kernel smoother
    !class(sll_c_particle_group_base),target,         intent(in)  :: particle_group(:) !< Particle group
    class(sll_t_particle_array), target,           intent(in)  :: particle_group
    sll_real64, target,                              intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                              intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                       intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                       intent(in)  :: Lx !< Length of the domain in x direction.
    sll_int32,                                      intent(in) :: n_sub_iter !< Number of subcycling steps
    character(*),                                   intent(in) :: file_prefix !< File prefix
    type( sll_t_binomial_filter ), intent( in ), target :: filter
    logical, optional, intent(in) :: jmean !< Should jmean be substracted in Ampere's law?



    !local variables
    sll_int32 :: ierr
    logical :: jmean_val

    SLL_ALLOCATE(sll_t_time_propagator_pic_vm_1d2v_subcyci :: splitting, ierr)


    if (present(jmean) ) then
       jmean_val = jmean
    else
       jmean_val = .false.
    end if

    select type (splitting)
    type is ( sll_t_time_propagator_pic_vm_1d2v_subcyci )
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
            file_prefix, &
            filter, &
            jmean_val)
    end select

  end subroutine sll_s_new_time_propagator_pic_vm_1d2v_subcyci_ptr


end module sll_m_time_propagator_pic_vm_1d2v_subcyci
