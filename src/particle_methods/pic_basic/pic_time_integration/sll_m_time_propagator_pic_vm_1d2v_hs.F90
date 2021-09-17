!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on Hamiltonian splitting for 1d2v Vlasov-Poisson.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
!> Reference: Kraus, Kormann, Sonnendrücker, Morrison: GEMPIC: Geometric ElectroMagnetic Particle-In-Cell Methods

!> Control variate: Note the we do not account for the analytic j at the moment (TODO: control_variate for current)
module sll_m_time_propagator_pic_vm_1d2v_hs
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

  use sll_m_filter_base_1d, only: &
       sll_c_filter_base_1d

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

  use sll_m_time_propagator_pic_vm_3d3v_cl_helper, only: &
       sll_p_boundary_particles_periodic, &
       sll_p_boundary_particles_singular, &
       sll_p_boundary_particles_reflection, &
       sll_p_boundary_particles_absorption

  implicit none

  public :: &
       sll_s_new_time_propagator_pic_vm_1d2v_hs, &
       sll_t_time_propagator_pic_vm_1d2v_hs

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_1d2v_hs
     class(sll_c_maxwell_1d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_0  !< Kernel smoother (order p+1)
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_1  !< Kernel smoother (order p)
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

     sll_int32 :: spline_degree !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx !< Size of the domain
     sll_real64 :: x_min !< Lower bound for x domain
     sll_real64 :: x_max
     sll_real64 :: delta_x !< Grid spacing
     sll_int32  :: n_cells
     sll_int32  :: n_total0
     sll_int32  :: n_total1
     sll_int32  :: boundary_particles = 100
     sll_int32 :: counter_left = 0
     sll_int32 :: counter_right = 0
     logical :: boundary = .false. !< true for non periodic boundary conditions
     sll_real64, pointer     :: rhob(:) => null()

     sll_real64 :: betar(2) !< reciprocal of plasma beta

     sll_real64, pointer     :: phi_dofs(:) !< DoFs describing then scalar potential phi
     sll_real64, pointer     :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
     sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the magnetic field
     sll_real64, allocatable :: j_dofs(:,:)      !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:,:)!< MPI-processor local part of one component of \a j_dofs

     sll_real64 :: force_sign = 1._f64
     logical :: electrostatic = .false.
     logical :: adiabatic_electrons = .false.

     logical :: jmean = .false.

     sll_real64, allocatable     :: efield_filter(:,:) !< DoFs describing the two components of the electric field
     sll_real64, allocatable     :: bfield_filter(:)   !< DoFs describing the magnetic field


     sll_real64, allocatable :: efield_to_val(:,:)
     sll_real64, allocatable :: bfield_to_val(:)

     class(sll_c_filter_base_1d), pointer :: filter

     sll_int32 :: n_species

     ! For version with control variate
     class(sll_t_control_variates), pointer :: control_variate => null()
     sll_int32 :: i_weight

   contains
     procedure :: operatorHp1 => operatorHp1_pic_vm_1d2v_hs  !> Operator for H_p1 part
     procedure :: operatorHp2 => operatorHp2_pic_vm_1d2v_hs  !> Operator for H_p2 part
     procedure :: operatorHE => operatorHE_pic_vm_1d2v_hs  !> Operator for H_E part
     procedure :: operatorHB => operatorHB_pic_vm_1d2v_hs  !> Operator for H_B part
     procedure :: lie_splitting => lie_splitting_pic_vm_1d2v_hs !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_1d2v_hs !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_1d2v_hs !> Strang splitting propagator
     procedure :: reinit_fields

     procedure :: init => initialize_pic_vm_1d2v_hs !> Initialize the type
     procedure :: free => delete_pic_vm_1d2v_hs !> Finalization

  end type sll_t_time_propagator_pic_vm_1d2v_hs

contains

  subroutine reinit_fields( self ) 
    class(sll_t_time_propagator_pic_vm_1d2v_hs), intent(inout) :: self !< time splitting object 

    call self%filter%apply( self%efield_dofs(:,1), self%efield_filter(:,1) ) 
    call self%filter%apply( self%efield_dofs(:,2), self%efield_filter(:,2) ) 
    call self%filter%apply( self%bfield_dofs, self%bfield_filter )

  end subroutine reinit_fields

  !> Strang splitting
  subroutine strang_splitting_pic_vm_1d2v_hs(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_hs), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%operatorHE(0.5_f64*dt)
          call self%operatorHp1(dt)
          call self%operatorHE(0.5_f64*dt)
       end do
    else
       do i_step = 1, number_steps
          call self%operatorHB(0.5_f64*dt)
          call self%operatorHE(0.5_f64*dt)
          call self%operatorHp2(0.5_f64*dt)
          call self%operatorHp1(dt)
          call self%operatorHp2(0.5_f64*dt)
          call self%operatorHE(0.5_f64*dt)
          call self%operatorHB(0.5_f64*dt)
       end do
    end if

  end subroutine strang_splitting_pic_vm_1d2v_hs

  !> Lie splitting
  subroutine lie_splitting_pic_vm_1d2v_hs(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_hs), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step
    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%operatorHE(dt)
          call self%operatorHp1(dt)
       end do
    else
       do i_step = 1,number_steps
          call self%operatorHE(dt)
          call self%operatorHB(dt)
          call self%operatorHp1(dt)
          call self%operatorHp2(dt)
       end do
    end if


  end subroutine lie_splitting_pic_vm_1d2v_hs

  !> Lie splitting (oposite ordering)
  subroutine lie_splitting_back_pic_vm_1d2v_hs(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_hs), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%operatorHp1(dt)
          call self%operatorHE(dt)
       end do
    else
       do i_step = 1,number_steps
          call self%operatorHp2(dt)
          call self%operatorHp1(dt)
          call self%operatorHB(dt)
          call self%operatorHE(dt)
       end do
    end if

  end subroutine lie_splitting_back_pic_vm_1d2v_hs


  !---------------------------------------------------------------------------!
  !> Push Hp1: Equations to solve are
  !> \partial_t f + v_1 \partial_{x_1} f = 0    -> X_new = X_old + dt V_1
  !> V_new,2 = V_old,2 + \int_0 h V_old,1 B_old
  !> \partial_t E_1 = - \int v_1 f(t,x_1, v) dv -> E_{1,new} = E_{1,old} - \int \int v_1 f(t,x_1+s v_1,v) dv ds
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B = 0 => B_new = B_old 
  subroutine operatorHp1_pic_vm_1d2v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_hs), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_part
    sll_real64 :: xnew(3), vi(3), wi(1), xold(3), wp(3)
    sll_int32  :: i_sp
    sll_real64 :: qoverm

    call self%maxwell_solver%transform_dofs( self%bfield_filter, self%bfield_to_val, 0 )


    ! Here we have to accumulate j and integrate over the time interval.
    ! At each k=1,...,n_grid, we have for s \in [0,dt]:
    ! j_k(s) =  \sum_{i=1,..,N_p} q_i N((x_k+sv_{1,k}-x_i)/h) v_k,
    ! where h is the grid spacing and N the normalized B-spline
    ! In order to accumulate the integrated j, we normalize the values of x to the grid spacing, calling them y, we have
    ! j_k(s) = \sum_{i=1,..,N_p} q_i N(y_k+s/h v_{1,k}-y_i) v_k.
    ! Now, we want the integral 
    ! \int_{0..dt} j_k(s) d s = \sum_{i=1,..,N_p} q_i v_k \int_{0..dt} N(y_k+s/h v_{1,k}-y_i) ds =  \sum_{i=1,..,N_p} q_i v_k  \int_{0..dt}  N(y_k + w v_{1,k}-y_i) dw


    self%j_dofs_local = 0.0_f64

    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
    ! Then update particle position:  X_new = X_old + dt * V
    do i_sp = 1,self%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          xold = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          ! Then update particle position:  X_new = X_old + dt * V
          xnew = xold + dt * vi

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

          call compute_particle_boundary_current( self, xold, xnew, vi, wi(1), qoverm )
          call self%particle_group%group(i_sp)%set_x(i_part, xnew)
          call self%particle_group%group(i_sp)%set_v(i_part, vi)

          if (self%particle_group%group(i_sp)%n_weights == 3) then
             ! Update weights if control variate
             wp = self%particle_group%group(i_sp)%get_weights(i_part)          
             wp(3) = self%control_variate%cv(i_sp)%update_df_weight( xnew(1:1), vi(1:2), 0.0_f64, wp(1), wp(2))
             call self%particle_group%group(i_sp)%set_weights(i_part, wp)
          end if

       end do
    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
         self%n_total1, MPI_SUM, self%j_dofs(1:self%n_total1,1))

    !call filter( self%j_dofs(:,1), self%j_dofs_local(:,1), n_cells )
    !call filter( self%j_dofs_local(:,1), self%j_dofs(:,1), n_cells )
    !self%j_dofs(:,1) = self%j_dofs_local(:,1)
    call self%filter%apply_inplace( self%j_dofs(:,1) )

    if ( self%jmean .eqv. .true. ) then
       self%j_dofs(:,1) = self%j_dofs(:,1) - sum(self%j_dofs(1:self%n_total1,1))/real(self%n_total1, f64)
    end if
    ! Update the electric field.
    if(self%adiabatic_electrons) then
       call self%maxwell_solver%compute_phi_from_j(self%j_dofs(1:self%n_total1,1), self%phi_dofs, self%efield_dofs(1:self%n_total1,1))
    else
       call self%maxwell_solver%compute_E_from_j(self%force_sign*self%betar(2)*self%j_dofs(1:self%n_total1,1), 1, self%efield_dofs(1:self%n_total1,1))
    end if
    call self%filter%apply( self%efield_dofs(:,1), self%efield_filter(:,1) )

  end subroutine operatorHp1_pic_vm_1d2v_hs


  subroutine compute_particle_boundary_current( self, xold, xnew, vi, wi, qoverm  )
    class(sll_t_time_propagator_pic_vm_1d2v_hs), intent( inout ) :: self !< time splitting object 
    sll_real64,                                           intent( inout ) :: xold(3)
    sll_real64,                                           intent( inout ) :: xnew(3)
    sll_real64,                                           intent( inout ) :: vi(3)
    sll_real64,                                           intent( in    ) :: wi(1)
    sll_real64,                                           intent( in    ) :: qoverm
    !local variables
    sll_real64 :: xmid(3), xbar

    if(xnew(1) < self%x_min .or. xnew(1) > self%x_max )then
       if(xnew(1) < self%x_min  )then
          xbar = self%x_min
          self%counter_left = self%counter_left+1
       else if(xnew(1) > self%x_max)then
          xbar = self%x_max
          self%counter_right = self%counter_right+1
       end if

       xmid = xnew
       xmid(1) = xbar

       call self%kernel_smoother_1%add_current_update_v( xold, xmid, wi(1), qoverm, &
            self%bfield_to_val, vi, self%j_dofs_local(1:self%n_total1,1))
       select case(self%boundary_particles)
       case(sll_p_boundary_particles_reflection)
          xnew(1) = 2._f64*xbar-xnew(1)
          vi(1) = -vi(1)
       case(sll_p_boundary_particles_absorption)
          call self%kernel_smoother_0%add_charge(xmid, wi(1), self%rhob)
       case( sll_p_boundary_particles_periodic)
          xnew(1) = self%x_min + modulo(xnew(1)-self%x_min, self%Lx)
          xmid = self%x_max+self%x_min-xbar
       case default
          xnew(1) = self%x_min + modulo(xnew(1)-self%x_min, self%Lx)
          xmid = self%x_max+self%x_min-xbar
       end select
       if (xnew(1) >= self%x_min .and. xnew(1) <= self%x_max) then
          call self%kernel_smoother_1%add_current_update_v( xmid, xnew, wi(1), qoverm, &
               self%bfield_to_val, vi, self%j_dofs_local(1:self%n_total1,1))
       end if
    else if(xnew(1) < -self%x_max .or.  xnew(1) > 2._f64*self%x_max ) then
       print*, xnew
       SLL_ERROR("particle boundary", "particle out of bound")
    else
       call self%kernel_smoother_1%add_current_update_v( xold, xnew, wi(1), qoverm, &
            self%bfield_to_val, vi, self%j_dofs_local(1:self%n_total1,1))
    end if

  end subroutine compute_particle_boundary_current


  !---------------------------------------------------------------------------!
  !> Push Hp2: Equations to solve are
  !> X_new = X_old
  !> V_new,1 = V_old,1 + \int_0 h V_old,2 B_old
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = - \int v_2 f(t,x_1, v) dv -> E_{2,new} = E_{2,old} - \int \int v_2 f(t,x_1+s v_1,v) dv ds
  !> \partial_t B = 0 => B_new = B_old
  subroutine operatorHp2_pic_vm_1d2v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_hs), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32  :: i_part, i_sp
    sll_real64 :: vi(3), xi(3), wi(1), wp(3)
    sll_real64 :: bfield
    sll_real64 :: qm


    self%j_dofs_local = 0.0_f64

    call self%maxwell_solver%transform_dofs( self%bfield_filter, self%bfield_to_val, 0 )

    do i_sp = 1, self%n_species
       qm = self%particle_group%group(i_sp)%species%q_over_m();
       ! Update v_1
       do i_part=1,self%particle_group%group(i_sp)%n_particles
          ! Evaluate bfield at particle position (splines of order p)
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          call self%kernel_smoother_1%evaluate &
               (xi(1), self%bfield_to_val, bfield)
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          vi(1) = vi(1) + dt*qm*vi(2)*bfield
          call self%particle_group%group(i_sp)%set_v(i_part, vi)

          xi = self%particle_group%group(i_sp)%get_x(i_part)

          ! Scale vi by weight to combine both factors for accumulation of integral over j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)*vi(2)

          call self%kernel_smoother_0%add_charge(xi(1:1), wi(1), self%j_dofs_local(:,2)) 

          if (self%particle_group%group(i_sp)%n_weights == 3) then
             ! Update weights if control variate
             wp = self%particle_group%group(i_sp)%get_weights(i_part)          
             wp(3) = self%control_variate%cv(i_sp)%update_df_weight( xi(1:1), vi(1:2), 0.0_f64, wp(1), wp(2))
             call self%particle_group%group(i_sp)%set_weights(i_part, wp)
          end if
       end do
    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
         self%n_total0, MPI_SUM, self%j_dofs(:,2))
    ! Update the electric field. Also, we still need to scale with 1/Lx ! TODO: Which scaling?
    if ( self%jmean .eqv. .true. ) then
       self%j_dofs(:,2) = self%j_dofs(:,2) - sum(self%j_dofs(:,2))/real(self%n_total0, f64)
    end if
    self%j_dofs(:,2) = self%j_dofs(:,2)*dt!/self%Lx


    call self%filter%apply_inplace( self%j_dofs(:,2) )

    if(self%adiabatic_electrons) then
       call self%maxwell_solver%compute_phi_from_j(self%j_dofs(:,2), self%phi_dofs, self%efield_dofs(:,1))
    else   
       call self%maxwell_solver%compute_E_from_j(self%force_sign*self%betar(2)*self%j_dofs(:,2), 2, self%efield_dofs(:,2))
    end if
    call self%filter%apply( self%efield_dofs(:,2), self%efield_filter(:,2) )

  end subroutine operatorHp2_pic_vm_1d2v_hs

  !---------------------------------------------------------------------------!
  !> Push H_E: Equations to be solved
  !> \partial_t f + E_1 \partial_{v_1} f + E_2 \partial_{v_2} f = 0 -> V_new = V_old + dt * E
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B + \partial_{x_1} E_2 = 0 => B_new = B_old - dt \partial_{x_1} E_2
  subroutine operatorHE_pic_vm_1d2v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_hs), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: v_new(3), xi(3), wp(3)
    sll_real64 :: efield(2)
    sll_real64 :: qm


    call self%maxwell_solver%transform_dofs( self%efield_filter(:,1), self%efield_to_val(:,1), 0 )    
    call self%maxwell_solver%transform_dofs( self%efield_filter(:,2), self%efield_to_val(:,2), 1 )

    do i_sp = 1,self%n_species
       qm = self%particle_group%group(i_sp)%species%q_over_m();
       ! V_new = V_old + dt * E
       do i_part=1,self%particle_group%group(i_sp)%n_particles
          v_new = self%particle_group%group(i_sp)%get_v(i_part)
          ! Evaluate efields at particle position
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          call self%kernel_smoother_1%evaluate &
               (xi(1), self%efield_to_val(1:self%n_total1,1), efield(1))
          call self%kernel_smoother_0%evaluate &
               (xi(1), self%efield_to_val(:,2), efield(2))
          v_new = self%particle_group%group(i_sp)%get_v(i_part)
          v_new(1:2) = v_new(1:2) + dt* qm * efield
          call self%particle_group%group(i_sp)%set_v(i_part, v_new)

          if (self%particle_group%group(i_sp)%n_weights == 3) then
             ! Update weights if control variate
             wp = self%particle_group%group(i_sp)%get_weights(i_part)          
             !wp(3) = wp(3) + dt*qm*efield(1)*v_new(1)*self%Lx
             wp(3) = self%control_variate%cv(i_sp)%update_df_weight( xi(1:1), v_new(1:2), 0.0_f64, wp(1), wp(2))
             call self%particle_group%group(i_sp)%set_weights(i_part, wp)
          end if
       end do
    end do

    if(self%electrostatic .eqv. .false.) then
       ! Update bfield
       call self%maxwell_solver%compute_B_from_E( &
            dt, self%efield_dofs(:,2), self%bfield_dofs)
       call self%filter%apply( self%bfield_dofs, self%bfield_filter )
    end if

  end subroutine operatorHE_pic_vm_1d2v_hs


  !---------------------------------------------------------------------------!
  !> Push H_B: Equations to be solved
  !> V_new = V_old
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old}
  !> \partial_t E_2 = - \partial_{x_1} B -> E_{2,new} = E_{2,old}-dt*\partial_{x_1} B
  !> \partial_t B = 0 -> B_new = B_old
  subroutine operatorHB_pic_vm_1d2v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_hs), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    ! Update efield2
    call self%maxwell_solver%compute_E_from_B(&
         dt*self%betar(1), self%bfield_dofs, self%efield_dofs(:,2))

    call self%filter%apply( self%efield_dofs(:,2), self%efield_filter(:,2) )

  end subroutine operatorHB_pic_vm_1d2v_hs


  !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_1d2v_hs(&
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
       force_sign, &
       jmean, &
       control_variate, &
       i_weight, &
       betar, &
       electrostatic, &
       rhob) 
    class(sll_t_time_propagator_pic_vm_1d2v_hs), intent(out) :: self !< time splitting object 
    class(sll_c_maxwell_1d_base), target,          intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,           intent(in)  :: particle_group
    sll_real64, target,                            intent(in)  :: phi_dofs(:) !< array for the coefficients of phi
    sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
    class( sll_c_filter_base_1d ), intent( in ), target :: filter
    sll_int32, optional, intent( in )  :: boundary_particles
    sll_real64, optional, intent(in) :: force_sign
    logical, optional, intent(in) :: jmean
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    sll_int32, optional,                            intent(in) :: i_weight !< Index of weight to be used by propagator
    sll_real64, optional, intent(in) :: betar(2) !< reciprocal plasma beta
    logical, optional    :: electrostatic
    sll_real64, optional, target,                  intent( in ) :: rhob(:)
    !local variables
    sll_int32 :: ierr

    if( present(electrostatic) )then
       self%electrostatic = electrostatic
    end if

    if( present(rhob) )then
       self%rhob => rhob
    end if

    if (present(betar)) then
       self%betar = betar!32.89_f64
    else
       self%betar = 1.0_f64
    end if

    if( particle_group%group(1)%species%q > 0._f64) self%adiabatic_electrons = .true.

    self%maxwell_solver => maxwell_solver
    self%kernel_smoother_0 => kernel_smoother_0
    self%kernel_smoother_1 => kernel_smoother_1
    self%n_species = particle_group%n_species
    self%particle_group => particle_group
    self%phi_dofs => phi_dofs
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs
    self%spline_degree = self%kernel_smoother_0%spline_degree
    self%x_min = x_min
    self%x_max = x_min + Lx
    self%Lx = Lx
    self%delta_x = self%Lx/real(self%maxwell_solver%n_cells,f64)
    self%n_cells = self%maxwell_solver%n_cells
    self%n_total0 = self%maxwell_solver%n_dofs0
    self%n_total1 = self%maxwell_solver%n_dofs1

    ! Check that n_dofs is the same for both kernel smoothers.
    SLL_ASSERT( self%kernel_smoother_0%n_cells == self%kernel_smoother_1%n_cells )
    SLL_ASSERT( self%kernel_smoother_0%n_cells == self%maxwell_solver%n_cells )

    SLL_ALLOCATE(self%j_dofs(self%n_total0,2), ierr)
    SLL_ALLOCATE(self%j_dofs_local(self%n_total0,2), ierr)
    SLL_ALLOCATE(self%efield_filter(self%n_total0,2), ierr)
    SLL_ALLOCATE(self%bfield_filter(self%n_total1), ierr)
    SLL_ALLOCATE(self%efield_to_val(self%n_total0,2), ierr)
    SLL_ALLOCATE(self%bfield_to_val(self%n_total1), ierr)

    if( self%n_cells+self%spline_degree == self%maxwell_solver%n_dofs0   ) then
       self%boundary = .true.
       self%boundary_particles = boundary_particles
    end if

    if( present(force_sign) )then
       self%force_sign = force_sign
    end if
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
    end if
  end subroutine initialize_pic_vm_1d2v_hs

  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_1d2v_hs(self)
    class(sll_t_time_propagator_pic_vm_1d2v_hs), intent( inout ) :: self !< time splitting object

    if( self%boundary ) then
       print*, 'left boundary', self%counter_left
       print*, 'right boundary', self%counter_right
    end if

    deallocate(self%j_dofs)
    deallocate(self%j_dofs_local)
    self%maxwell_solver => null()
    self%kernel_smoother_0 => null()
    self%kernel_smoother_1 => null()
    self%particle_group => null()
    self%phi_dofs => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()

  end subroutine delete_pic_vm_1d2v_hs


  !---------------------------------------------------------------------------!
  !> Constructor for allocatable abstract type.
  subroutine sll_s_new_time_propagator_pic_vm_1d2v_hs(&
       splitting, &
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
       force_sign, &
       jmean, &
       control_variate, &
       i_weight, &
       betar, &
       electrostatic, &
       rhob) 
    class(sll_c_time_propagator_base), allocatable, intent(out) :: splitting !< time splitting object 
    class(sll_c_maxwell_1d_base), target,                intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,                intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,                intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,           intent(in)  :: particle_group
    !class(sll_c_particle_group_base),target,             intent(in)  :: particle_group(:) !< Particle group
    sll_real64, target,                            intent(in)  :: phi_dofs(:) !< array for the coefficients of phi
    sll_real64, target,                                  intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                                  intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                           intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                           intent(in)  :: Lx !< Length of the domain in x direction.
    class( sll_c_filter_base_1d ), intent( in ), target :: filter
    sll_int32, optional,                                  intent( in )  :: boundary_particles
    sll_real64, optional, intent(in) :: force_sign
    logical, optional, intent(in) :: jmean !< Should jmean be substracted in Ampere's law?
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    sll_int32, optional,                            intent(in) :: i_weight !< Index of weight to be used by propagator
    sll_real64, intent( in ), optional :: betar(2)
    logical, optional    :: electrostatic
    sll_real64, optional, target,                  intent( in ) :: rhob(:)
    !local variables
    sll_int32 :: ierr, boundary_particles_set 
    sll_real64 :: force_sign_set
    logical :: jmean_set
    sll_real64 :: betar_set(2)
    logical :: electrostatic_set

    SLL_ALLOCATE(sll_t_time_propagator_pic_vm_1d2v_hs :: splitting, ierr)

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

    if (present(jmean) ) then
       jmean_set = jmean
    else
       jmean_set = .false.
    end if

    if (present(betar) ) then
       betar_set = betar
    else
       betar_set = 1._f64
    end if

    if (present(electrostatic) ) then
       electrostatic_set = electrostatic
    else
       electrostatic = .false.
    end if

    select type (splitting)
    type is ( sll_t_time_propagator_pic_vm_1d2v_hs )
       if (present(control_variate) ) then
          call splitting%init(&
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
               boundary_particles_set, &
               force_sign_set, &
               jmean_set, &
               control_variate, &
               i_weight, &
               betar = betar_set, &
               electrostatic=electrostatic_set, &
               rhob = rhob)
       else
          call splitting%init(&
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
               boundary_particles=boundary_particles_set, &
               force_sign=force_sign_set, &
               jmean = jmean_set, &
               betar = betar_set, &
               electrostatic=electrostatic_set, &
               rhob = rhob)
       end if
    end select

  end subroutine sll_s_new_time_propagator_pic_vm_1d2v_hs



end module sll_m_time_propagator_pic_vm_1d2v_hs
