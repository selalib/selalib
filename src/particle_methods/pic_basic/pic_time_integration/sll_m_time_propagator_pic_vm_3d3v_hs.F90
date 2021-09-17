!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on Hamiltonian splitting for 3d3v Vlasov-Maxwell.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
!> Reference: Kraus, Kormann, Sonnendrücker, Morrison: GEMPIC: Geometric ElectroMagnetic Particle-In-Cell Methods
module sll_m_time_propagator_pic_vm_3d3v_hs
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
       sll_o_collective_allreduce, &
       sll_v_world_collective

  use sll_m_control_variate, only: &
       sll_t_control_variates

  use sll_m_filter_base_3d, only: &
       sll_c_filter_base_3d

  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use sll_m_time_propagator_pic_vm_3d3v_cl_helper, only: &
       sll_p_boundary_particles_periodic, &
       sll_p_boundary_particles_singular, &
       sll_p_boundary_particles_reflection, &
       sll_p_boundary_particles_absorption

  use sll_m_initial_distribution, only: &
       sll_t_params_cos_gaussian_screwpinch

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_mpi, only: &
       mpi_sum

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use sll_m_particle_mesh_coupling_base_3d, only: &
       sll_c_particle_mesh_coupling_3d

  use sll_m_profile_functions, only: &
       sll_t_profile_functions

  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_3d3v_hs

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 3d3v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_3d3v_hs
     class(sll_c_maxwell_3d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_particle_mesh_coupling_3d), pointer :: particle_mesh_coupling
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

     sll_int32 :: spline_degree(3) !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx(3) !< Size of the domain
     sll_real64 :: x_min(3) !< Lower bound for x domain
     sll_real64 :: x_max(3) !< Upper bound for x domain
     sll_int32 :: n_total0  !< total number of Dofs for 0form
     sll_int32 :: n_total1  !< total number of Dofs for 1form

     sll_real64 :: betar(2) !< reciprocal of plasma beta

     sll_real64, pointer     :: phi_dofs(:) !< DoFs describing the scalar potential
     sll_real64, pointer     :: efield_dofs(:) !< DoFs describing the three components of the electric field
     sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the three components of the magnetic field
     sll_real64, allocatable :: j_dofs(:)      !< DoFs for representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:)!< MPI-processor local part of one component of \a j_dofs

     sll_int32 :: boundary_particles = 100 !< particle boundary conditions
     sll_int32 :: counter_left = 0 !< boundary counter
     sll_int32 :: counter_right = 0 !< boundary counter
     sll_real64, pointer     :: rhob(:) => null() !< charge at the boundary

     sll_real64 :: force_sign = 1._f64 !< sign of particle force
     logical :: electrostatic = .false. !< true for electrostatic simulation
     logical :: jmean = .false. !< logical for mean value of current
     logical :: adiabatic_electrons = .false. !< true for simulation with adiabatic electrons
     logical :: lindf = .false. !< true for simulation with linear delta f method
     logical :: boundary = .false. !< true for non periodic boundary conditions

     !filtering
     sll_real64, allocatable     :: efield_filter(:) !< DoFs describing the two components of the electric field
     sll_real64, allocatable     :: bfield_filter(:)   !< DoFs describing the magnetic field
     class(sll_c_filter_base_3d), pointer :: filter  !< filter

     sll_int32 :: n_species !< number of species

     ! For control variate
     class(sll_t_control_variates), pointer :: control_variate => null()!< control variate
     sll_int32 :: i_weight = 1 !< number of weights


   contains
     procedure :: operatorHp1 => operatorHp1_pic_vm_3d3v_hs  !> Operator for H_p1 part
     procedure :: operatorHp2 => operatorHp2_pic_vm_3d3v_hs  !> Operator for H_p2 part
     procedure :: operatorHp3 => operatorHp3_pic_vm_3d3v_hs  !> Operator for H_p3 part
     procedure :: operatorHE => operatorHE_pic_vm_3d3v_hs  !> Operator for H_E part
     procedure :: operatorHB => operatorHB_pic_vm_3d3v_hs  !> Operator for H_B part
     procedure :: lie_splitting => lie_splitting_pic_vm_3d3v_hs !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_3d3v_hs !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_3d3v_hs !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_3d3v_hs !> Initialize the type
     procedure :: free => delete_pic_vm_3d3v_hs !> Finalization

  end type sll_t_time_propagator_pic_vm_3d3v_hs

contains


  !> Strang splitting
  subroutine strang_splitting_pic_vm_3d3v_hs(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%operatorHE(0.5_f64*dt)
          call self%operatorHp3(0.5_f64*dt)
          call self%operatorHp2(0.5_f64*dt)
          call self%operatorHp1(dt)
          call self%operatorHp2(0.5_f64*dt)
          call self%operatorHp3(0.5_f64*dt)
          call self%operatorHE(0.5_f64*dt)
       end do
    else
       do i_step = 1, number_steps
          call self%operatorHB(0.5_f64*dt)
          call self%operatorHE(0.5_f64*dt)
          call self%operatorHp3(0.5_f64*dt)
          call self%operatorHp2(0.5_f64*dt)
          call self%operatorHp1(dt)
          call self%operatorHp2(0.5_f64*dt)
          call self%operatorHp3(0.5_f64*dt)
          call self%operatorHE(0.5_f64*dt)
          call self%operatorHB(0.5_f64*dt)
       end do
    end if

  end subroutine strang_splitting_pic_vm_3d3v_hs


  !> Lie splitting
  subroutine lie_splitting_pic_vm_3d3v_hs(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1,number_steps
          call self%operatorHE(dt)
          call self%operatorHp1(dt)
          call self%operatorHp2(dt)
          call self%operatorHp3(dt)
       end do
    else
       do i_step = 1,number_steps
          call self%operatorHE(dt)
          call self%operatorHB(dt)
          call self%operatorHp1(dt)
          call self%operatorHp2(dt)
          call self%operatorHp3(dt)
       end do
    end if

  end subroutine lie_splitting_pic_vm_3d3v_hs


  !> Lie splitting (oposite ordering)
  subroutine lie_splitting_back_pic_vm_3d3v_hs(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1,number_steps
          call self%operatorHp3(dt)
          call self%operatorHp2(dt)
          call self%operatorHp1(dt)
          call self%operatorHE(dt)
       end do
    else
       do i_step = 1,number_steps
          call self%operatorHp3(dt)
          call self%operatorHp2(dt)
          call self%operatorHp1(dt)
          call self%operatorHB(dt)
          call self%operatorHE(dt)
       end do
    end if

  end subroutine lie_splitting_back_pic_vm_3d3v_hs


  !---------------------------------------------------------------------------!
  !> Push Hp1: Equations to solve are
  !> \partial_t f + v_1 \partial_{x_1} f = 0    -> xnew = xold + dt V_1
  !> V_new,2 = V_old,2 + \int_0 h V_old,1 B_old
  !> \partial_t E_1 = - \int v_1 f(t,x_1, v) dv -> E_{1,new} = E_{1,old} - \int \int v_1 f(t,x_1+s v_1,v) dv ds
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B = 0 => B_new = B_old 
  subroutine operatorHp1_pic_vm_3d3v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_3d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: xnew(3), vi(3), wi(1), xold(3), wall(3)
    sll_real64 :: qoverm
    sll_int32 :: i_weight
    sll_real64 :: work(self%n_total1+2*self%n_total0)

    i_weight = self%i_weight

    ! Here we have to accumulate j and integrate over the time interval.
    ! At each k=1,...,n_grid, we have for s \in [0,dt]:
    ! j_k(s) =  \sum_{i=1,..,N_p} q_i N((x_k+sv_{1,k}-x_i)/h) v_k,
    ! where h is the grid spacing and N the normalized B-spline
    ! In order to accumulate the integrated j, we normalize the values of x to the grid spacing, calling them y, we have
    ! j_k(s) = \sum_{i=1,..,N_p} q_i N(y_k+s/h v_{1,k}-y_i) v_k.
    ! Now, we want the integral 
    ! \int_{0..dt} j_k(s) d s = \sum_{i=1,..,N_p} q_i v_k \int_{0..dt} N(y_k+s/h v_{1,k}-y_i) ds =  \sum_{i=1,..,N_p} q_i v_k  \int_{0..dt}  N(y_k + w v_{1,k}-y_i) dw


    self%j_dofs_local = 0.0_f64

    self%bfield_filter = self%bfield_dofs
    call self%filter%apply_inplace(self%bfield_filter(1:self%n_total0))
    call self%filter%apply_inplace(self%bfield_filter(self%n_total0+1:self%n_total0+self%n_total1))
    call self%filter%apply_inplace(self%bfield_filter(self%n_total0+self%n_total1+1:self%n_total0+2*self%n_total1))

    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
    ! Then update particle position:  xnew = xold + dt * V
    do i_sp = 1, self%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part = 1, self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          xold = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          ! Then update particle position:  xnew = xold + dt * V
          xnew(1) = xold(1) + dt * vi(1)
          xnew(2:3) = xold(2:3)

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, i_weight)

          call compute_particle_boundary_current( self, xold, xnew, vi, wi, qoverm )
          call self%particle_group%group(i_sp)%set_x(i_part, xnew)
          call self%particle_group%group(i_sp)%set_v(i_part, vi)

          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             if( self%lindf .eqv. .false.) then
                wall = self%particle_group%group(i_sp)%get_weights(i_part)
                wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xnew, vi, 0.0_f64, wall(1), wall(2) )
                call self%particle_group%group(i_sp)%set_weights( i_part, wall )
             end if
          end if
       end do
    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         self%n_total1, MPI_SUM, self%j_dofs(1:self%n_total1))

    if ( self%jmean .eqv. .true. ) then
       self%j_dofs(1:self%n_total1) = self%j_dofs(1:self%n_total1) - sum(self%j_dofs(1:self%n_total1))/real(self%n_total1, f64)
    end if

    call self%filter%apply_inplace(self%j_dofs)
    ! Update the electric field.
    if( self%adiabatic_electrons) then
       work = 0._f64
       work(1:self%n_total1) = self%j_dofs(1:self%n_total1)
       call self%maxwell_solver%compute_phi_from_j( work, self%phi_dofs, self%efield_dofs )
    else
       call self%maxwell_solver%compute_E_from_j( self%force_sign*self%betar(2)*self%j_dofs(1:self%n_total1), self%efield_dofs(1:self%n_total1), 1)
    end if

  end subroutine operatorHp1_pic_vm_3d3v_hs


  !> Helper function for \a operatorHp1
  subroutine compute_particle_boundary_current( self, xold, xnew, vi, wi, qoverm  )
    class(sll_t_time_propagator_pic_vm_3d3v_hs), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( inout ) :: xold(3)
    sll_real64,                                           intent( inout ) :: xnew(3)
    sll_real64,                                           intent( inout ) :: vi(3)
    sll_real64,                                           intent( in    ) :: wi(1)
    sll_real64,                                           intent( in    ) :: qoverm
    !local variables
    sll_real64 :: xmid(3), xbar, dx

    if(xnew(1) < -self%x_max(1) .or.  xnew(1) > 2._f64*self%x_max(1) ) then
       print*, xnew
       SLL_ERROR("particle boundary", "particle out of bound")
    else if(xnew(1) < self%x_min(1) .or. xnew(1) > self%x_max(1) )then
       if(xnew(1) < self%x_min(1)  )then
          xbar = self%x_min(1)
          self%counter_left = self%counter_left+1
       else if(xnew(1) > self%x_max(1))then
          xbar = self%x_max(1)
          self%counter_right = self%counter_right+1
       end if

       dx = (xbar- xold(1))/(xnew(1)-xold(1))
       xmid = xold + dx * (xnew-xold)
       xmid(1) = xbar
       call self%particle_mesh_coupling%add_current_update_v_component1( xold, xmid(1), wi(1), qoverm, &
            self%bfield_filter, vi, self%j_dofs_local(1:self%n_total1))
       select case(self%boundary_particles)
       case(sll_p_boundary_particles_reflection)
          xnew(1) = 2._f64*xbar-xnew(1)
          vi(1) = -vi(1)
       case(sll_p_boundary_particles_absorption)
          call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
       case( sll_p_boundary_particles_periodic)
          xnew(1) = self%x_min(1) + modulo(xnew(1)-self%x_min(1), self%Lx(1))
          xmid(1) = self%x_max(1)+self%x_min(1)-xbar
       case default
          xnew(1) = self%x_min(1) + modulo(xnew(1)-self%x_min(1), self%Lx(1))
          xmid(1) = self%x_max(1)+self%x_min(1)-xbar
       end select
       if (xnew(1) >= self%x_min(1) .and. xnew(1) <= self%x_max(1)) then
          call self%particle_mesh_coupling%add_current_update_v_component1( xmid, xnew(1), wi(1), qoverm, &
               self%bfield_filter, vi, self%j_dofs_local(1:self%n_total1))
       end if
    else
       call self%particle_mesh_coupling%add_current_update_v_component1( xold, xnew(1), wi(1), qoverm, &
            self%bfield_filter, vi, self%j_dofs_local(1:self%n_total1))
    end if

  end subroutine compute_particle_boundary_current


  !---------------------------------------------------------------------------!
  !> Push Hp2: Equations to solve are
  subroutine operatorHp2_pic_vm_3d3v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_3d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: xnew(3), vi(3), wi(1), xold(3), wall(3)
    sll_real64 :: qoverm
    sll_int32 :: i_weight
    sll_real64 :: work(self%n_total1+2*self%n_total0)

    i_weight = self%i_weight

    self%j_dofs_local = 0.0_f64

    self%bfield_filter = self%bfield_dofs
    call self%filter%apply_inplace(self%bfield_filter(1:self%n_total0))
    call self%filter%apply_inplace(self%bfield_filter(self%n_total0+1:self%n_total0+self%n_total1))
    call self%filter%apply_inplace(self%bfield_filter(self%n_total0+self%n_total1+1:self%n_total0+2*self%n_total1))

    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
    ! Then update particle position:  xnew = xold + dt * V
    do i_sp = 1, self%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part = 1, self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          xold = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          ! Then update particle position:  xnew = xold + dt * V
          xnew(2) = xold(2) + dt * vi(2)
          xnew(1) = xold(1)
          xnew(3) = xold(3)

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, i_weight)

          call self%particle_mesh_coupling%add_current_update_v_component2( xold, xnew(2), wi(1), qoverm, &
               self%bfield_filter, vi, self%j_dofs_local)

          xnew(2) =  self%x_min(2) + modulo(xnew(2)-self%x_min(2), self%Lx(2)) 
          call self%particle_group%group(i_sp)%set_x(i_part, xnew)
          call self%particle_group%group(i_sp)%set_v(i_part, vi)

          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             if( self%lindf .eqv. .false.) then
                wall = self%particle_group%group(i_sp)%get_weights(i_part)
                wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xnew, vi, 0.0_f64, wall(1), wall(2) )
                call self%particle_group%group(i_sp)%set_weights( i_part, wall )
             end if
          end if
       end do
    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         self%n_total0, MPI_SUM, self%j_dofs)
    if ( self%jmean .eqv. .true. ) then
       self%j_dofs = self%j_dofs - sum(self%j_dofs)/real(self%n_total0, f64)
    end if

    call self%filter%apply_inplace(self%j_dofs)

    ! Update the electric field.
    if( self%adiabatic_electrons) then
       work = 0._f64
       work(self%n_total1+1:self%n_total0+self%n_total1) = self%j_dofs
       call self%maxwell_solver%compute_phi_from_j( work, self%phi_dofs, self%efield_dofs )
    else
       call self%maxwell_solver%compute_E_from_j(self%force_sign*self%betar(2)*self%j_dofs, self%efield_dofs(self%n_total1+1:self%n_total1+self%n_total0), 2)
    end if

  end subroutine operatorHp2_pic_vm_3d3v_hs


  !---------------------------------------------------------------------------!
  !> Push Hp3: Equations to solve are
  subroutine operatorHp3_pic_vm_3d3v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_3d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: xnew(3), vi(3), wi(1), xold(3), wall(3)
    sll_real64 :: qoverm
    sll_int32 :: i_weight
    sll_real64 :: work(self%n_total1+2*self%n_total0)

    i_weight = self%i_weight


    self%j_dofs_local = 0.0_f64

    self%bfield_filter = self%bfield_dofs
    call self%filter%apply_inplace(self%bfield_filter(1:self%n_total0))
    call self%filter%apply_inplace(self%bfield_filter(self%n_total0+1:self%n_total0+self%n_total1))
    call self%filter%apply_inplace(self%bfield_filter(self%n_total0+self%n_total1+1:self%n_total0+2*self%n_total1))

    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
    ! Then update particle position:  xnew = xold + dt * V
    do i_sp = 1, self%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part = 1, self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          xold = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          ! Then update particle position:  xnew = xold + dt * V
          xnew(3) = xold(3) + dt * vi(3)
          xnew(1:2) = xold(1:2)

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, i_weight)

          call self%particle_mesh_coupling%add_current_update_v_component3( xold, xnew(3), wi(1), qoverm, &
               self%bfield_filter, vi, self%j_dofs_local)

          xnew(3) =  self%x_min(3) + modulo(xnew(3)-self%x_min(3), self%Lx(3)) 
          call self%particle_group%group(i_sp)%set_x(i_part, xnew)
          call self%particle_group%group(i_sp)%set_v(i_part, vi)

          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             if( self%lindf .eqv. .false.) then
                wall = self%particle_group%group(i_sp)%get_weights(i_part)
                wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xnew, vi, 0.0_f64, wall(1), wall(2) )
                call self%particle_group%group(i_sp)%set_weights( i_part, wall )
             end if
          end if
       end do
    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         self%n_total0, MPI_SUM, self%j_dofs)
    if ( self%jmean .eqv. .true. ) then
       self%j_dofs = self%j_dofs - sum(self%j_dofs)/real(self%n_total0, f64)
    end if

    call self%filter%apply_inplace(self%j_dofs)
    ! Update the electric field.
    if( self%adiabatic_electrons) then
       work = 0._f64
       work(self%n_total1+self%n_total0+1:self%n_total1+2*self%n_total0) = self%j_dofs
       call self%maxwell_solver%compute_phi_from_j( work, self%phi_dofs, self%efield_dofs )
    else
       call self%maxwell_solver%compute_E_from_j(self%force_sign*self%betar(2)*self%j_dofs, self%efield_dofs(self%n_total0+self%n_total1+1:self%n_total1+2*self%n_total0), 3)
    end if

  end subroutine operatorHp3_pic_vm_3d3v_hs


  !---------------------------------------------------------------------------!
  !> Push H_E: Equations to be solved
  !> $V^{n+1}=V^n+\Delta t\mathbb{W}_{\frac{q}{m}} \mathbb{Lambda}^1(X^n) e^n$
  !> $b^{n+1}=b^n-\Delta t C e^n$
  subroutine operatorHE_pic_vm_3d3v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_3d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: v_new(3), xi(3), Rx(3), wall(3)
    sll_real64 :: efield(3)
    sll_real64 :: qoverm, q, m
    sll_int32 :: i_weight

    i_weight = self%i_weight

    self%efield_filter = self%efield_dofs
    call self%filter%apply_inplace(self%efield_filter(1:self%n_total1))
    call self%filter%apply_inplace(self%efield_filter(1+self%n_total1:self%n_total1+self%n_total0))
    call self%filter%apply_inplace(self%efield_filter(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0))

    do i_sp = 1, self%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       q = self%particle_group%group(i_sp)%species%q
       m = self%particle_group%group(i_sp)%species%m
       ! V_new = V_old + dt * E
       do i_part=1,self%particle_group%group(i_sp)%n_particles
          v_new = self%particle_group%group(i_sp)%get_v(i_part)
          ! Evaluate efields at particle position
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
               self%efield_filter(1:self%n_total1), efield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
               self%efield_filter(1+self%n_total1:self%n_total1+self%n_total0), efield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
               self%efield_filter(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0), efield(3))
          if( self%lindf .eqv. .false.) then
             v_new = v_new + dt* qoverm * efield
             call self%particle_group%group(i_sp)%set_v(i_part, v_new)

             if (self%particle_group%group(i_sp)%n_weights == 3 ) then
                wall = self%particle_group%group(i_sp)%get_weights(i_part)
                wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xi, v_new, 0.0_f64, wall(1), wall(2) )
                call self%particle_group%group(i_sp)%set_weights( i_part, wall )
             end if
          else             
             ! Update particle weights
             if (self%particle_group%group(i_sp)%n_weights == 3 ) then
                wall = self%particle_group%group(i_sp)%get_weights(i_part)
                select type(p => self%control_variate%cv(i_sp)%control_variate_distribution_params)
                type is(sll_t_params_cos_gaussian_screwpinch)
                   Rx = xi
                   Rx(1) = Rx(1) + m*v_new(2)/q
                   Rx = (Rx-self%x_min)/self%Lx

                   wall(3) = wall(3) + dt* (q/p%profile%T_i(Rx(1)) * sum(efield * v_new)-&
                        efield(2)*(p%profile%drho_0(Rx(1))/p%profile%rho_0(Rx(1))+(0.5_f64*m*sum(v_new**2)/p%profile%T_i(Rx(1)) - 1.5_f64)* p%profile%dT_i(Rx(1))/p%profile%T_i(Rx(1))  ) ) *&
                        p%eval_v_density(v_new, Rx, m)/p%eval_v_density(v_new, xi, m)*product(self%Lx)
                class default
                   wall(3) = wall(3) + dt* qoverm* sum(efield * v_new) *product(self%Lx)
                end select
                call self%particle_group%group(i_sp)%set_weights( i_part, wall )
             else
                v_new = v_new + dt* qoverm * efield
                call self%particle_group%group(i_sp)%set_v(i_part, v_new) 
             end if
          end if
       end do
    end do

    if(self%electrostatic .eqv. .false.) then
       ! Update bfield
       call self%maxwell_solver%compute_B_from_E( &
            dt, self%efield_dofs, self%bfield_dofs)
    end if

  end subroutine operatorHE_pic_vm_3d3v_hs


  !---------------------------------------------------------------------------!
  !> Push H_B: Equations to be solved
  !> $M_1 e^{n+1}=M_1 e^n+\Delta t C^\top M_2 b^n$
  subroutine operatorHB_pic_vm_3d3v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_3d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    ! Update efield2
    call self%maxwell_solver%compute_E_from_B(&
         dt*self%betar(1), self%bfield_dofs, self%efield_dofs)

  end subroutine operatorHB_pic_vm_3d3v_hs


  !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_3d3v_hs(&
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
       boundary_particles, &
       betar, &
       force_sign, &
       electrostatic, &
       rhob, &
       control_variate, &
       jmean, &
       lindf)
    class(sll_t_time_propagator_pic_vm_3d3v_hs), intent(out) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,           intent(in)  :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target,   intent(in)  :: particle_mesh_coupling !< Particle mesh coupling
    class(sll_t_particle_array), target,      intent(in)  :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:) !< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent(in)  :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min(3) !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx(3) !< Length of the domain in x direction. 
    class(sll_c_filter_base_3d), target :: filter !< filter
    sll_int32, optional,                           intent( in ) :: boundary_particles !< particle boundary conditions

    sll_real64, optional, intent(in) :: betar(2) !< reciprocal plasma beta
    sll_real64, optional, intent(in) :: force_sign !< sign of particle force
    logical, optional, intent(in)    :: electrostatic !< true for electrostatic simulation
    sll_real64, optional, target, intent( in ) :: rhob(:)  !< charge at the boundary
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    logical, optional, intent(in) :: jmean !< logical for mean value of current
    logical, optional, intent(in) :: lindf !< true for linear delta f method
    !local variables
    sll_int32 :: ierr

    if( present(electrostatic) )then
       self%electrostatic = electrostatic
    end if

    if( present(force_sign) )then
       self%force_sign = force_sign
    end if

    if (self%force_sign == 1._f64) then
       if( particle_group%group(1)%species%q > 0._f64) self%adiabatic_electrons = .true.
    end if

    if (present(boundary_particles) ) then
       self%boundary_particles = boundary_particles
    end if

    if( present(rhob) )then
       self%rhob => rhob
    end if


    self%maxwell_solver => maxwell_solver
    self%particle_mesh_coupling => particle_mesh_coupling
    self%particle_group => particle_group
    self%phi_dofs => phi_dofs
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs

    self%n_total0 = self%maxwell_solver%n_total0
    self%n_total1 = self%maxwell_solver%n_total1

    SLL_ALLOCATE(self%j_dofs(self%n_total0), ierr)
    SLL_ALLOCATE(self%j_dofs_local(self%n_total0), ierr)

    self%spline_degree = self%particle_mesh_coupling%spline_degree
    self%x_min = x_min
    self%x_max = x_min + Lx
    self%Lx = Lx

    if( self%particle_mesh_coupling%n_cells(1)+self%spline_degree(1) == self%maxwell_solver%n_dofs(1)   ) then
       self%boundary = .true.
    end if

    self%n_species = particle_group%n_species

    self%filter => filter

    if (present(control_variate)) then
       allocate(self%control_variate )
       allocate(self%control_variate%cv(self%n_species) )
       self%control_variate => control_variate
       self%i_weight = 3
       if (present(lindf)) then
          self%lindf = lindf
       end if
    end if

    if (present(betar)) then
       self%betar = betar
    else
       self%betar = 1.0_f64
    end if

  end subroutine initialize_pic_vm_3d3v_hs


  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_3d3v_hs(self)
    class(sll_t_time_propagator_pic_vm_3d3v_hs), intent( inout ) :: self !< time propagator object 

    deallocate(self%j_dofs)
    deallocate(self%j_dofs_local)
    self%maxwell_solver => null()
    self%particle_mesh_coupling => null()
    self%particle_group => null()
    self%phi_dofs => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()

  end subroutine delete_pic_vm_3d3v_hs


end module sll_m_time_propagator_pic_vm_3d3v_hs
