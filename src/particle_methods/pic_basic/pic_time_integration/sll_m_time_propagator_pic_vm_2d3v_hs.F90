!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on Hamiltonian splitting for 2d3v Vlasov-Maxwell.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
!> Reference: Kraus, Kormann, Sonnendrücker, Morrison: GEMPIC: Geometric ElectroMagnetic Particle-In-Cell Methods
module sll_m_time_propagator_pic_vm_2d3v_hs
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_o_collective_allreduce, &
    sll_v_world_collective
  
  use sll_m_control_variate, only: &
    sll_t_control_variates

  use sll_m_time_propagator_base, only: &
    sll_c_time_propagator_base

  use sll_m_particle_mesh_coupling_spline_2d_feec, only: &
    sll_t_particle_mesh_coupling_spline_2d_feec

  use sll_m_maxwell_2d_fem_fft, only: &
    sll_t_maxwell_2d_fem_fft

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base, &
       sll_t_particle_array

  use mpi, only: &
       mpi_sum

  implicit none

  public :: &
    sll_s_new_time_propagator_pic_vm_2d3v_hs, &
    sll_s_new_time_propagator_pic_vm_2d3v_hs_ptr, &
    sll_t_time_propagator_pic_vm_2d3v_hs

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_2d3v_hs
     type(sll_t_maxwell_2d_fem_fft), pointer :: maxwell_solver      !< Maxwell solver
     type(sll_t_particle_mesh_coupling_spline_2d_feec), pointer :: particle_mesh_coupling
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

     sll_int32 :: spline_degree !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx(2) !< Size of the domain
     sll_real64 :: x_min(2) !< Lower bound for x domain
     !sll_real64 :: delta_x !< Grid spacing
     sll_int32 :: ntotal !< Total number of points on the tensor product grid

     sll_real64 :: betar !< reciprocal of plasma beta
     
     sll_real64, pointer     :: efield_dofs(:) !< DoFs describing the two components of the electric field
     sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the magnetic field
     sll_real64, allocatable :: j_dofs(:)      !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:)!< MPI-processor local part of one component of \a j_dofs

     ! For control variate
     class(sll_t_control_variates), pointer :: control_variate => null()!< control variate
     sll_int32 :: i_weight = 1 !< number of weights

   contains
     procedure :: operatorHp1 => operatorHp1_pic_vm_2d3v_hs  !> Operator for H_p1 part
     procedure :: operatorHp2 => operatorHp2_pic_vm_2d3v_hs  !> Operator for H_p2 part
     procedure :: operatorHp3 => operatorHp3_pic_vm_2d3v_hs  !> Operator for H_p3 part
     procedure :: operatorHE => operatorHE_pic_vm_2d3v_hs  !> Operator for H_E part
     procedure :: operatorHB => operatorHB_pic_vm_2d3v_hs  !> Operator for H_B part
     procedure :: lie_splitting => lie_splitting_pic_vm_2d3v_hs !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_2d3v_hs !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_2d3v_hs !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_2d3v_hs !> Initialize the type
     procedure :: free => delete_pic_vm_2d3v_hs !> Finalization

  end type sll_t_time_propagator_pic_vm_2d3v_hs

contains
  
  
  !> Strang splitting
  subroutine strang_splitting_pic_vm_2d3v_hs(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_2d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

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

  end subroutine strang_splitting_pic_vm_2d3v_hs

  !> Lie splitting
  subroutine lie_splitting_pic_vm_2d3v_hs(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_2d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1,number_steps
       call self%operatorHE(dt)
       call self%operatorHB(dt)
       call self%operatorHp1(dt)
       call self%operatorHp2(dt)
       call self%operatorHp3(dt)
    end do


  end subroutine lie_splitting_pic_vm_2d3v_hs

  !> Lie splitting (oposite ordering)
  subroutine lie_splitting_back_pic_vm_2d3v_hs(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_2d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1,number_steps
       call self%operatorHp3(dt)
       call self%operatorHp2(dt)
       call self%operatorHp1(dt)
       call self%operatorHB(dt)
       call self%operatorHE(dt)
    end do

  end subroutine lie_splitting_back_pic_vm_2d3v_hs


  !---------------------------------------------------------------------------!
  !> Push Hp1: Equations to solve are
  !> TODO: UPDATE
  !> \partial_t f + v_1 \partial_{x_1} f = 0    -> X_new = X_old + dt V_1
  !> V_new,2 = V_old,2 + \int_0 h V_old,1 B_old
  !> \partial_t E_1 = - \int v_1 f(t,x_1, v) dv -> E_{1,new} = E_{1,old} - \int \int v_1 f(t,x_1+s v_1,v) dv ds
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B = 0 => B_new = B_old 
  subroutine operatorHp1_pic_vm_2d3v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_2d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: x_new(3), vi(3), wi(1), x_old(3), wall(3)
    sll_int32  :: n_cells
    sll_real64 :: qoverm
    sll_int32 :: i_weight

    i_weight = self%i_weight
 
    n_cells = self%maxwell_solver%n_total

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
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          x_old = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          
          ! Then update particle position:  X_new = X_old + dt * V
          x_new(1) = x_old(1) + dt * vi(1)
          x_new(2:3) = x_old(2:3)
          
          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, i_weight)
          
          call self%particle_mesh_coupling%add_current_update_v_component1( x_old, x_new(1), wi(1), qoverm, &
               self%bfield_dofs, vi, self%j_dofs_local)
          
          x_new(1) = modulo(x_new(1), self%Lx(1))
          call self%particle_group%group(i_sp)%set_x(i_part, x_new)
          call self%particle_group%group(i_sp)%set_v(i_part, vi)
          
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( x_new, vi, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do
    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         n_cells, MPI_SUM, self%j_dofs)

    ! Update the electric field.
    call self%maxwell_solver%compute_E_from_j(self%j_dofs, 1, self%efield_dofs(1:self%ntotal))

  end subroutine operatorHp1_pic_vm_2d3v_hs
  
  !---------------------------------------------------------------------------!
  !> Push Hp2: Equations to solve are
  subroutine operatorHp2_pic_vm_2d3v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_2d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: x_new(3), vi(3), wi(1), x_old(3), wall(3)
    sll_int32  :: n_cells
    sll_real64 :: qoverm
    sll_int32 :: i_weight

    i_weight = self%i_weight
 
    n_cells = self%maxwell_solver%n_total

    self%j_dofs_local = 0.0_f64

    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
    ! Then update particle position:  X_new = X_old + dt * V
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          x_old = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          
          ! Then update particle position:  X_new = X_old + dt * V
          x_new(2) = x_old(2) + dt * vi(2)
          x_new(1) = x_old(1)
          
          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, i_weight)

          call self%particle_mesh_coupling%add_current_update_v_component2( x_old, x_new(2), wi(1), qoverm, &
               self%bfield_dofs, vi, self%j_dofs_local)
          
          x_new(2) = modulo(x_new(2), self%Lx(2))
          call self%particle_group%group(i_sp)%set_x(i_part, x_new)
          call self%particle_group%group(i_sp)%set_v(i_part, vi)
          
          
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( x_new, vi, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do

    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         n_cells, MPI_SUM, self%j_dofs)

    ! Update the electric field.
    call self%maxwell_solver%compute_E_from_j(self%j_dofs, 2, self%efield_dofs(self%ntotal+1:2*self%ntotal))

  end subroutine operatorHp2_pic_vm_2d3v_hs

  !---------------------------------------------------------------------------!
  !> Push Hp3: Equations to solve are
  subroutine operatorHp3_pic_vm_2d3v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_2d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: vi(3), wi(1), x_old(3), wall(3)
    sll_int32  :: n_cells
    sll_real64 :: qoverm
    sll_real64 :: bfield(2)
    sll_int32 :: i_weight

    i_weight = self%i_weight
 
    n_cells = self%maxwell_solver%n_total

    self%j_dofs_local = 0.0_f64

    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
    ! Then update particle position:  X_new = X_old + dt * V
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          x_old = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, i_weight)

          call self%particle_mesh_coupling%evaluate &
               (x_old(1:2), [self%spline_degree, self%spline_degree-1], &
               self%bfield_dofs(1:self%ntotal), bfield(1))
          call self%particle_mesh_coupling%evaluate &
               (x_old(1:2), [self%spline_degree-1, self%spline_degree], &
               self%bfield_dofs(self%ntotal+1:self%ntotal*2), bfield(2))

          vi(1) = vi(1) - dt *qoverm*vi(3)*bfield(2)
          vi(2) = vi(2) + dt *qoverm*vi(3)*bfield(1)

          ! was pp
          call self%particle_mesh_coupling%add_charge( x_old(1:2), wi(1)*vi(3), &
               [self%spline_degree, self%spline_degree], self%j_dofs_local )
          
          call self%particle_group%group(i_sp)%set_v(i_part, vi)
          
          
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( x_old, vi, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do

    end do
    
    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         n_cells, MPI_SUM, self%j_dofs)

    ! Update the electric field.
    self%j_dofs = self%j_dofs*dt
    call self%maxwell_solver%compute_E_from_j(self%j_dofs, 3, self%efield_dofs(2*self%ntotal+1:3*self%ntotal))

  end subroutine operatorHp3_pic_vm_2d3v_hs

  
  !---------------------------------------------------------------------------!
  !> Push H_E: Equations to be solved
  !> \partial_t f + E_1 \partial_{v_1} f + E_2 \partial_{v_2} f = 0 -> V_new = V_old + dt * E
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B + \partial_{x_1} E_2 = 0 => B_new = B_old - dt \partial_{x_1} E_2
  subroutine operatorHE_pic_vm_2d3v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_2d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: v_new(3), xi(3), wall(3)
    sll_real64 :: efield(3)
    sll_real64 :: qm
    sll_int32 :: i_weight

    i_weight = self%i_weight

    ! TODO: Optimize using evaluate multiple

    !write(25,*) self%efield_dofs
    
    do i_sp = 1, self%particle_group%n_species
       qm = self%particle_group%group(i_sp)%species%q_over_m();
       ! V_new = V_old + dt * E
       do i_part=1,self%particle_group%group(i_sp)%n_particles
          v_new = self%particle_group%group(i_sp)%get_v(i_part)
          ! Evaluate efields at particle position
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          
          !print*, self%efield_dofs
          call self%particle_mesh_coupling%evaluate &
               (xi(1:2), [self%spline_degree-1, self%spline_degree], &
               self%efield_dofs(1:self%ntotal), efield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi(1:2), [self%spline_degree, self%spline_degree-1], &
               self%efield_dofs(self%ntotal+1:2*self%ntotal), efield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi(1:2), [self%spline_degree, self%spline_degree], &
               self%efield_dofs(2*self%ntotal+1:3*self%ntotal), efield(3))

          v_new = v_new + dt* qm * efield
          call self%particle_group%group(i_sp)%set_v(i_part, v_new)
          
       
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xi, v_new, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do
    end do
    
    ! Update bfield
    call self%maxwell_solver%compute_B_from_E( &
         dt, self%efield_dofs, self%bfield_dofs)
        
  end subroutine operatorHE_pic_vm_2d3v_hs
  

  !---------------------------------------------------------------------------!
  !> Push H_B: Equations to be solved
  !> V_new = V_old
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old}
  !> \partial_t E_2 = - \partial_{x_1} B -> E_{2,new} = E_{2,old}-dt*\partial_{x_1} B
  !> \partial_t B = 0 -> B_new = B_old
  subroutine operatorHB_pic_vm_2d3v_hs(self, dt)
    class(sll_t_time_propagator_pic_vm_2d3v_hs), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    ! Update efield2
    call self%maxwell_solver%compute_E_from_B(&
         dt*self%betar, self%bfield_dofs, self%efield_dofs)
      
  end subroutine operatorHB_pic_vm_2d3v_hs


 !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_2d3v_hs(&
       self, &
       maxwell_solver, &
       particle_mesh_coupling, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       control_variate, &
       betar)
    class(sll_t_time_propagator_pic_vm_2d3v_hs), intent(out) :: self !< time propagator object 
    type(sll_t_maxwell_2d_fem_fft), pointer,           intent(in)  :: maxwell_solver !< Maxwell solver
    type(sll_t_particle_mesh_coupling_spline_2d_feec), pointer, intent(in)  :: particle_mesh_coupling !< Particle mesh coupling
    class(sll_t_particle_array), pointer,      intent(in)  :: particle_group !< Particle group
    sll_real64, pointer,                            intent(in)  :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, pointer,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min(2) !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx(2) !< Length of the domain in x direction.
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    sll_real64, optional, intent(in) :: betar !< reciprocal plasma beta
        !local variables
    sll_int32 :: ierr

    self%maxwell_solver => maxwell_solver
    self%particle_mesh_coupling => particle_mesh_coupling
    self%particle_group => particle_group
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs

    self%ntotal =self%particle_mesh_coupling%n_dofs 
    SLL_ALLOCATE(self%j_dofs(self%ntotal), ierr)
    SLL_ALLOCATE(self%j_dofs_local(self%ntotal), ierr)

    self%spline_degree = self%particle_mesh_coupling%spline_degree
    self%x_min = x_min
    self%Lx = Lx

    if (present(control_variate)) then
       allocate(self%control_variate )
       allocate(self%control_variate%cv(particle_group%n_species) )
       self%control_variate => control_variate
      self%i_weight = 3
    end if

    if (present(betar)) then
       self%betar = betar!32.89_f64
    else
       self%betar = 1.0_f64
    end if

  end subroutine initialize_pic_vm_2d3v_hs

  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_2d3v_hs(self)
    class(sll_t_time_propagator_pic_vm_2d3v_hs), intent( inout ) :: self !< time propagator object 

    deallocate(self%j_dofs)
    deallocate(self%j_dofs_local)
    self%maxwell_solver => null()
    self%particle_mesh_coupling => null()
    self%particle_group => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()

  end subroutine delete_pic_vm_2d3v_hs


  !---------------------------------------------------------------------------!
  !> Constructor for allocatable abstract type.
  subroutine sll_s_new_time_propagator_pic_vm_2d3v_hs(&
       splitting, &
       maxwell_solver, &
       particle_mesh_coupling, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       control_variate, &
       betar)
    class(sll_c_time_propagator_base), allocatable, intent(out) :: splitting !< time propagator object 
    type(sll_t_maxwell_2d_fem_fft), pointer,          intent(in)  :: maxwell_solver !< Maxwell solver
    type(sll_t_particle_mesh_coupling_spline_2d_feec), pointer, intent(in)  :: particle_mesh_coupling  !< Particle mesh coupling
    class(sll_t_particle_array),pointer,             intent(in)  :: particle_group !< Particle group
    sll_real64, pointer,                                  intent(in)  :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, pointer,                                  intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                           intent(in)  :: x_min(2) !< Lower bound of x domain
    sll_real64,                                           intent(in)  :: Lx(2) !< Length of the domain in x direction.
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    sll_real64, intent( in ), optional :: betar !< reciprocal of plasma beta
    !local variables
    sll_int32 :: ierr
    
    
    SLL_ALLOCATE(sll_t_time_propagator_pic_vm_2d3v_hs :: splitting, ierr)

    select type (splitting)
    type is ( sll_t_time_propagator_pic_vm_2d3v_hs )
       if (present(control_variate) ) then
          if (present(betar) ) then
             call splitting%init(&
                  maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  control_variate=control_variate, &
                  betar=betar)
          else
             call splitting%init(&
                  maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  control_variate=control_variate)
          end if
       else
          if (present(betar) ) then
             call splitting%init(&
                  maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  betar=betar)
          else
             call splitting%init(&
                  maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx)
          end if
       end if
          
    end select

  end subroutine sll_s_new_time_propagator_pic_vm_2d3v_hs

  !---------------------------------------------------------------------------!
  !> Constructor for pointer abstract type.
  subroutine sll_s_new_time_propagator_pic_vm_2d3v_hs_ptr(&
       splitting, &
       maxwell_solver, &
       particle_mesh_coupling, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx) 
    class(sll_c_time_propagator_base), pointer, intent(out) :: splitting !< time propagator object 
    type(sll_t_maxwell_2d_fem_fft), pointer,            intent(in)  :: maxwell_solver !< Maxwell solver
    type(sll_t_particle_mesh_coupling_spline_2d_feec), pointer, intent(in)  :: particle_mesh_coupling  !< Particle mesh coupling
    class(sll_t_particle_array),pointer,         intent(in)  :: particle_group !< Particle group
    sll_real64, pointer,                              intent(in)  :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, pointer,                              intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                       intent(in)  :: x_min(2) !< Lower bound of x domain
    sll_real64,                                       intent(in)  :: Lx(2) !< Length of the domain in x direction.
    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE(sll_t_time_propagator_pic_vm_2d3v_hs :: splitting, ierr)

    select type (splitting)
    type is ( sll_t_time_propagator_pic_vm_2d3v_hs )
       call splitting%init(&
            maxwell_solver, &
            particle_mesh_coupling, &
            particle_group, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx)
    end select

  end subroutine sll_s_new_time_propagator_pic_vm_2d3v_hs_ptr


end module sll_m_time_propagator_pic_vm_2d3v_hs
