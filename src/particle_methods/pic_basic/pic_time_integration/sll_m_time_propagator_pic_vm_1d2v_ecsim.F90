!> @ingroup pic_time_integration
!> @author Benedikt Perse, IPP
!> @brief Particle pusher based on Lapentas splitting in Ecsim for 1d2v Vlasov-Poisson.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_time_propagator_pic_vm_1d2v_ecsim
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
       sll_o_collective_allreduce, &
       sll_v_world_collective

  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use sll_m_particle_mesh_coupling_base_1d, only: &
       sll_c_particle_mesh_coupling_1d

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_gauss_legendre_integration, only : &
       sll_f_gauss_legendre_points_and_weights

  use sll_mpi, only: &
       mpi_sum

  use sll_m_linear_operator_ecsim, only : &
       sll_t_linear_operator_ecsim

  use sll_m_linear_solver_mgmres, only : &
       sll_t_linear_solver_mgmres

  use sll_m_spline_fem_utilities, only : &
       sll_s_spline_fem_mass_line

  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_1d2v_ecsim

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_1d2v_ecsim
  class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_0  !< Kernel smoother (order p+1)
  class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_1  !< Kernel smoother (order p)
  class(sll_t_particle_array), pointer  :: particle_group    !< Particle group
  sll_int32 :: spline_degree !< Degree of the spline for j,B. Here 3.
  sll_real64 :: Lx !< Size of the domain
  sll_real64 :: x_min !< Lower bound for x domain
  sll_real64 :: domain(2) !< domain given by x_min and x_max
  sll_real64 :: delta_x !< Grid spacing
  sll_int32 :: n_dofs !< Number of degrees of freedom 

  sll_real64 :: cell_integrals_0(4) !< Integral over the spline function on each interval (order p+1)
  sll_real64 :: cell_integrals_1(3) !< Integral over the spline function on each interval (order p)


  sll_real64, pointer     :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
  sll_real64, allocatable :: efield_dofs_local(:) !< DoFs describing the first components of the electric field for the computation of the Gauss law
  sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the magnetic field
  sll_real64, allocatable :: j_dofs(:,:)      !< DoFs for kernel representation of current density. 
  sll_real64, allocatable :: j_dofs_local(:,:)!< MPI-processor local part of one component of \a j_dofs
  sll_real64, allocatable :: rho_dofs(:)      !< DoFs for kernel representation of charge density. 
  sll_real64, allocatable :: rho_dofs_local(:)!< MPI-processor local part of one component of \a rho_dofs


  sll_real64, allocatable :: mass_line_0(:) !< Entries for the massmatrix M0
  sll_real64, allocatable :: mass_line_1(:) !< Entries for the massmatrix M1

  sll_real64, allocatable :: m1(:,:),m2(:,:),m4(:,:) !< Blockmatrices which build the particle-massmatrix
  sll_real64, allocatable :: m1_local(:,:),m2_local(:,:),m4_local(:,:) !< !< MPI-processor local part of the blockmatrices
  type(sll_t_linear_operator_ecsim) :: linear_operator !< Linear operator which defines the matrix-vectorproduct used to calculate the advancing fields
  type(sll_t_linear_solver_mgmres) :: linear_solver !< Iterator to solve the system of fieldequations
 
contains
  procedure :: advect_x => advect_x_pic_vm_1d2v_ecsim !> Advection of positionvariable
  procedure :: advect_e => advect_e_pic_vm_1d2v_ecsim !> Advection of the velocity- and the fieldvariables

  procedure :: lie_splitting => lie_splitting_pic_vm_1d2v_ecsim !> Lie splitting propagator
  procedure :: lie_splitting_back => lie_splitting_back_pic_vm_1d2v_ecsim !> Lie splitting propagator
  procedure :: strang_splitting => strang_splitting_pic_vm_1d2v_ecsim !> Strang splitting propagator

  procedure :: init => initialize_pic_vm_1d2v_ecsim !> Initialize the type
  procedure :: init_from_file => initialize_file_pic_vm_1d2v_ecsim !> Initialize the type
  procedure :: free => delete_pic_vm_1d2v_ecsim !> Finalization

end type sll_t_time_propagator_pic_vm_1d2v_ecsim


contains

  !> Strang splitting
subroutine strang_splitting_pic_vm_1d2v_ecsim(self,dt, number_steps)
 class(sll_t_time_propagator_pic_vm_1d2v_ecsim), intent(inout) :: self !< time splitting object 
 sll_real64,                                     intent(in)    :: dt   !< time step
 sll_int32,                                      intent(in)    :: number_steps !< number of time steps
 ! local variable
 sll_int32 :: i_step

 do i_step = 1, number_steps
    call self%advect_x(dt*0.5_f64)
    call self%advect_e(dt)
    call self%advect_x(dt*0.5_f64)
 end do



end subroutine strang_splitting_pic_vm_1d2v_ecsim

!> Lie splitting
subroutine lie_splitting_pic_vm_1d2v_ecsim(self,dt, number_steps)
 class(sll_t_time_propagator_pic_vm_1d2v_ecsim), intent(inout) :: self !< time splitting object 
 sll_real64,                                     intent(in)    :: dt   !< time step
 sll_int32,                                      intent(in)    :: number_steps !< number of time steps
 ! local variable
 sll_int32 :: i_step

 do i_step = 1, number_steps
    call self%advect_x(dt)
    call self%advect_e(dt)

 end do

end subroutine lie_splitting_pic_vm_1d2v_ecsim

!> Lie splitting
subroutine lie_splitting_back_pic_vm_1d2v_ecsim(self,dt, number_steps)
 class(sll_t_time_propagator_pic_vm_1d2v_ecsim), intent(inout) :: self !< time splitting object 
 sll_real64,                                     intent(in)    :: dt   !< time step
 sll_int32,                                      intent(in)    :: number_steps !< number of time steps
 ! local variable
 sll_int32 :: i_step

 do i_step = 1, number_steps
    call self%advect_e(dt)
    call self%advect_x(dt)
 end do

end subroutine lie_splitting_back_pic_vm_1d2v_ecsim

!---------------------------------------------------------------------------!
subroutine advect_x_pic_vm_1d2v_ecsim ( self, dt )
 class(sll_t_time_propagator_pic_vm_1d2v_ecsim), intent(inout) :: self !< time splitting object 
 sll_real64,                                     intent(in)    :: dt   !< time step
 ! local variables
 sll_int32 :: i_part, i_sp
 sll_real64 :: xi(3), vi(3) !< particle position and velocity

 do i_sp = 1,self%particle_group%n_species
    do i_part=1,self%particle_group%group(i_sp)%n_particles  
       ! Read out particle position and velocity
       xi = self%particle_group%group(i_sp)%get_x(i_part)
       vi = self%particle_group%group(i_sp)%get_v(i_part)

       xi(1) = xi(1) + dt * vi(1)
       xi(1) = modulo(xi(1), self%Lx)
       call self%particle_group%group(i_sp)%set_x ( i_part, xi )
    end do
 end do

end subroutine advect_x_pic_vm_1d2v_ecsim

!---------------------------------------------------------------------------!
subroutine advect_e_pic_vm_1d2v_ecsim ( self, dt )
 class(sll_t_time_propagator_pic_vm_1d2v_ecsim), intent(inout) :: self !< time splitting object 
 sll_real64,                                     intent(in)    :: dt   !< time step
 ! local variables
 sll_int32 :: i_part, i_sp
 sll_real64 :: vi(3),v_new(3), xi(3), wi !< particle velocity, position and charge times particle weight
 sll_real64 :: qm, beta  !< charge divided by mass, delta t/2 times qoverm 
 sll_real64 :: mu !< data for particle mass matrix
 sll_real64 :: efield(2), bfield !< scratch data for evaluated fields at one particle position
 sll_real64 :: scratch(3*self%n_dofs) !< scratch data
 sll_real64 :: rhs(3*self%n_dofs) !< righthandside of the equationsystem


 self%linear_operator%dt=dt
 self%linear_operator%dx=self%delta_x
 ! Set to zero
 self%j_dofs_local = 0.0_f64
 self%j_dofs = 0.0_f64
 self%rho_dofs_local = 0.0_f64
 self%rho_dofs = 0.0_f64
 self%m1=0.0_f64
 self%m2=0.0_f64 
 self%m4=0.0_f64
 self%m1_local=0.0_f64
 self%m2_local=0.0_f64 
 self%m4_local=0.0_f64
 rhs=0.0_f64
 scratch=0.0_f64
 bfield=0.0_f64
 efield=0.0_f64

 ! First particle loop
 do i_sp = 1, self%particle_group%n_species
    qm = self%particle_group%group(i_sp)%species%q_over_m();
    beta=qm*dt* 0.5_f64;
    do i_part = 1,self%particle_group%group(i_sp)%n_particles
       ! Read out particle position and velocity
       vi = self%particle_group%group(i_sp)%get_v(i_part)
       xi = self%particle_group%group(i_sp)%get_x(i_part)

       ! Get charge for accumulation of j
       wi = self%particle_group%group(i_sp)%get_charge(i_part)
       call self%kernel_smoother_1%evaluate &
            (xi(1), self%bfield_dofs, bfield)
       mu=(1+(beta*bfield)**2.0_f64)

       ! Accumulate jx
       call self%kernel_smoother_1%add_charge( xi(1), wi*(vi(1)+vi(2)*beta* bfield)/mu, self%j_dofs_local(:,1) )

       ! Accumulate jy
       call self%kernel_smoother_0%add_charge( xi(1), wi*(vi(2)-vi(1)*beta* bfield)/mu, self%j_dofs_local(:,2) )

       !Compute particle-mass matrices m1,m2,m4
       call self%kernel_smoother_1%add_particle_mass( xi(1), (beta*wi)/mu,self%m1_local)

       call add_particle_mass_mixed_spline_1d(self, self%kernel_smoother_0, self%spline_degree, self%spline_degree-1,  xi(1), -(beta**2.0_f64*wi*bfield)/mu,self%m2_local)

       call self%kernel_smoother_0%add_particle_mass( xi(1), (beta*wi)/mu,self%m4_local)

    end do
 end do



 ! MPI to sum up contributions from each processor
 call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
      self%kernel_smoother_1%n_dofs, MPI_SUM, self%j_dofs(:,1) )
 call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
      self%kernel_smoother_1%n_dofs, MPI_SUM, self%j_dofs(:,2) )
 call sll_o_collective_allreduce( sll_v_world_collective, self%m1_local, &
      self%kernel_smoother_1%n_dofs*self%spline_degree, MPI_SUM, self%m1 )
 call sll_o_collective_allreduce( sll_v_world_collective, self%m2_local, &
      self%kernel_smoother_0%n_dofs*2*self%spline_degree, MPI_SUM, self%m2 )
 call sll_o_collective_allreduce( sll_v_world_collective, self%m4_local, &
      self%kernel_smoother_0%n_dofs* (self%spline_degree+1), MPI_SUM, self%m4 )

 ! maxwell equations solved for fields at next timestep
 call righthandside( self, dt, rhs)
 call self%linear_solver%set_guess([self%efield_dofs(:,1),self%efield_dofs(:,2),self%bfield_dofs])
 call self%linear_solver%solve(rhs, scratch)


 self%efield_dofs(:,1)=0.5_f64*(scratch(1:self%n_dofs)+self%efield_dofs(:,1))
 self%efield_dofs(:,2)=0.5_f64*(scratch(self%n_dofs+1:2*self%n_dofs)+self%efield_dofs(:,2))
 bfield=0.0_f64
 efield=0.0_f64
 ! Second particle loop for v update
 do i_sp = 1, self%particle_group%n_species
    qm = self%particle_group%group(i_sp)%species%q_over_m();
    beta=qm*dt* 0.5_f64;
    do i_part = 1,self%particle_group%group(i_sp)%n_particles
       ! Read out particle position and velocity
       vi = self%particle_group%group(i_sp)%get_v(i_part)
       xi = self%particle_group%group(i_sp)%get_x(i_part)
       ! evaluate fields at particle position
       call self%kernel_smoother_1%evaluate &
            (xi(1), self%bfield_dofs, bfield)
       mu=1+(beta*bfield)**2.0_f64
       call self%kernel_smoother_1%evaluate &
            (xi(1), self%efield_dofs(:,1), efield(1))
       call self%kernel_smoother_0%evaluate &
            (xi(1), self%efield_dofs(:,2), efield(2))

       !calculate v**(n+1/2)
       v_new(1)=(vi(1)+beta*efield(1)+beta*bfield*(vi(2)+beta*efield(2)))/mu
       v_new(2)=(vi(2)+beta*efield(2)-beta*bfield*(vi(1)+beta*efield(1)))/mu


       v_new(1)=2.0_f64*v_new(1)-vi(1)
       v_new(2)=2.0_f64*v_new(2)-vi(2)

       call self%particle_group%group(i_sp)%set_v(i_part, v_new)
    end do
 end do

 ! update fields
 self%efield_dofs(:,1)= scratch(1:self%n_dofs)  
 self%efield_dofs(:,2)= scratch(self%n_dofs+1:2*self%n_dofs)
 self%bfield_dofs= scratch(2*self%n_dofs+1:3*self%n_dofs)


end subroutine advect_e_pic_vm_1d2v_ecsim


!> Add charge of one particle for splines of different degrees
subroutine add_particle_mass_mixed_spline_1d(self, kernel_smoother, degree1, degree2, position, marker_charge, particle_mass)
 class(sll_t_time_propagator_pic_vm_1d2v_ecsim), intent(inout) :: self !< Time splitting object 
 class(sll_c_particle_mesh_coupling_1d),      intent( in )    :: kernel_smoother !< Kernel smoother
 sll_int32,                                intent( in )    :: degree1 !< Degrees of first spline
 sll_int32,                                intent( in )    :: degree2 !< Degree of second spline
 sll_real64,                               intent( in )    :: position(self%kernel_smoother_1%dim) !< Position of the particle
 sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
 sll_real64,                               intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution
 !local variables
 sll_int32 :: i1, column
 sll_int32 :: index1d, index
 sll_real64 :: xi(1)
 sll_real64 :: spline_val1(degree1+1) !< scratch data for spline evaluation
 sll_real64 :: spline_val2(degree2+1) !< scratch data for spline evaluation

 SLL_ASSERT( size(particle_mass,1) == degree1+degree2+1 )
 SLL_ASSERT( size(particle_mass,2) == kernel_smoother%n_dofs )

 xi(1) = (position(1) - self%x_min)/self%delta_x
 index = ceiling(xi(1))
 xi(1) = xi(1) - real(index-1, f64)
 index = index - degree1


 call sll_s_uniform_bsplines_eval_basis(degree1, xi(1), spline_val1)
 call sll_s_uniform_bsplines_eval_basis(degree2, xi(1), spline_val2)

 do i1 = 1, degree1+1
    index1d = modulo(index+i1-2,kernel_smoother%n_cells)+1 !index-degree1:index
    do column = 1, degree2+1
       particle_mass(column+degree1+1-i1, index1d) = particle_mass(column+degree1+1-i1, index1d)+ &
            marker_charge * spline_val1(i1) * spline_val2(column)
    end do
 end do

end subroutine add_particle_mass_mixed_spline_1d
!> calculation of the righthandside of the system of maxwellequations for halftimestep

!> calculation of the righthandside of the system of maxwellequations 
subroutine righthandside (self,dt,rhs)
 class(sll_t_time_propagator_pic_vm_1d2v_ecsim), intent(inout) :: self !< Time splitting object
 sll_real64, intent(in) :: dt !< Timestep
 sll_real64, intent(out):: rhs(3*self%n_dofs) !< Righthandside of the equationsystem
 ! local variables
 sll_int32 :: row, column
 sll_real64 :: scratch(self%n_dofs+1) ! scratch data
 ! set to zero 
 scratch=0.0_f64
 do row=1, self%n_dofs
    ! multiply field-dofs with diagonal entries of mass and particle-mass matrix and substract current
    rhs(row) = (self%mass_line_1(1)-0.5_f64*dt*self%m1(1,row))*self%efield_dofs(row,1)-dt*self%j_dofs(row,1)
    rhs(self%n_dofs+row)= (self%mass_line_0(1)-0.5_f64*dt*self%m4(1,row))*self%efield_dofs(row,2)-dt*self%j_dofs(row,2)
    scratch(row)=self%mass_line_1(1)*self%bfield_dofs(row)
    do column = 2, self%spline_degree
       rhs(row) = rhs(row) +&
            self%mass_line_1(column) * &
            (self%efield_dofs(modulo(row+column-2,self%n_dofs)+1,1) +&
            self%efield_dofs(modulo(row-column,self%n_dofs)+1,1))-&
            0.5_f64*dt*self%m1(column,row) * &
            self%efield_dofs(modulo(row+column-2,self%n_dofs)+1,1)-& 
            0.5_f64*dt*self%m1(column,modulo(row-column,self%n_dofs)+1) * &
            self%efield_dofs(modulo(row-column,self%n_dofs)+1,1) 

       scratch(row)=scratch(row)+&
            self%mass_line_1(column) * &
            (self%bfield_dofs(modulo(row+column-2,self%n_dofs)+1) +&
            self%bfield_dofs(modulo(row-column,self%n_dofs)+1))
    end do
    do column = 2, self%spline_degree+1
       rhs(self%n_dofs+row) = rhs(self%n_dofs+row) +&
            self%mass_line_0(column) * &
            (self%efield_dofs(modulo(row+column-2,self%n_dofs)+1,2) +&
            self%efield_dofs(modulo(row-column,self%n_dofs)+1,2))-&
            0.5_f64*dt*self%m4(column,row)*&
            self%efield_dofs(modulo(row+column-2,self%n_dofs)+1,2)-&
            0.5_f64*dt*self%m4(column,modulo(row-column,self%n_dofs)+1)*&
            self%efield_dofs(modulo(row-column,self%n_dofs)+1,2)
    end do
    ! multiply efields with entries of m2 and m3
    do column = 1, 2*self%spline_degree
       rhs(row) = rhs(row) + &
            0.5_f64*dt*self%m2(7-column,modulo(column-4+row-1,self%n_dofs)+1) *&
            self%efield_dofs(modulo(row+column-self%spline_degree-2,self%n_dofs)+1,2) 
       rhs(row+self%n_dofs) = rhs(row+self%n_dofs) - &
            0.5_f64*dt*self%m2(column,row) *&
            self%efield_dofs(modulo(row+column-self%spline_degree-1,self%n_dofs)+1,1)
    end do
 end do
 scratch(self%n_dofs+1)=scratch(1)
 do row=1, self%n_dofs
    ! update efield with curl of bfield 
    rhs(self%n_dofs+row)=  rhs(self%n_dofs+row)+&
         0.5_f64*dt/self%delta_x*(scratch(row)-scratch(row+1))
    ! update bfield with curl of efield 
    rhs(2*self%n_dofs+row)=self%bfield_dofs(row)-&
         0.5_f64*dt/self%delta_x*(self%efield_dofs(row,2)-self%efield_dofs(modulo(row-2,self%n_dofs)+1,2))

 end do
end subroutine righthandside


!> Constructor.
subroutine initialize_file_pic_vm_1d2v_ecsim(&
    self, &
    kernel_smoother_0, &
    kernel_smoother_1, &
    particle_group, &
    efield_dofs, &
    bfield_dofs, &
    x_min, &
    Lx, &
    filename ) 
  class(sll_t_time_propagator_pic_vm_1d2v_ecsim), intent(out) :: self !< time splitting object 
  class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
  class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
  class(sll_t_particle_array), target,      intent(in)  :: particle_group !< Particle group
  sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
  sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
  sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
  sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
  character(len=*), intent(in) :: filename

  sll_int32 :: input_file
  sll_int32 :: io_stat
  sll_real64 :: maxwell_tolerance
     
  namelist /time_solver/ maxwell_tolerance 

  ! Read in solver tolerance
  open(newunit = input_file, file=filename, status='old',IOStat=io_stat)
  if (io_stat /= 0) then
     print*, 'sll_m_time_propagator_pic_vm_1d2v_ecsim: Input file does not exist. Set default tolerance.'
     call self%init( kernel_smoother_0, &
          kernel_smoother_1, &
          particle_group, &
          efield_dofs, &
          bfield_dofs, &
          x_min, &
          Lx )
  else       
     read(input_file, time_solver)
     call self%init( kernel_smoother_0, &
          kernel_smoother_1, &
          particle_group, &
          efield_dofs, &
          bfield_dofs, &
          x_min, &
          Lx, &
          solver_tolerance=maxwell_tolerance)
     close(input_file)
  end if
 

end subroutine initialize_file_pic_vm_1d2v_ecsim


!---------------------------------------------------------------------------!
!> Constructor.
subroutine initialize_pic_vm_1d2v_ecsim(&
    self, &
    kernel_smoother_0, &
    kernel_smoother_1, &
    particle_group, &
    efield_dofs, &
    bfield_dofs, &
    x_min, &
    Lx, &
    solver_tolerance ) 
 class(sll_t_time_propagator_pic_vm_1d2v_ecsim), intent(out) :: self !< time splitting object 
 class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
 class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
 class(sll_t_particle_array), target,      intent(in)  :: particle_group !< Particle group
 sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
 sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
 sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
 sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
 sll_real64, intent(in), optional                            :: solver_tolerance !< solver tolerance
 !local variable
 sll_int32 :: ierr


 self%kernel_smoother_0 => kernel_smoother_0
 self%kernel_smoother_1 => kernel_smoother_1
 self%particle_group => particle_group
 self%efield_dofs => efield_dofs
 self%bfield_dofs => bfield_dofs

 ! Check that n_dofs is the same for both kernel smoothers.
 SLL_ASSERT( self%kernel_smoother_0%n_dofs == self%kernel_smoother_1%n_dofs )
 SLL_ALLOCATE(self%j_dofs(self%kernel_smoother_0%n_dofs,2), ierr)
 SLL_ALLOCATE(self%j_dofs_local(self%kernel_smoother_0%n_dofs,2), ierr)

 self%spline_degree = self%kernel_smoother_0%spline_degree
 self%x_min = x_min
 self%Lx = Lx
 self%domain=[x_min,x_min+Lx]
 self%delta_x = self%Lx/real(self%kernel_smoother_1%n_dofs, f64)
 self%n_dofs=self%kernel_smoother_1%n_dofs

 self%cell_integrals_1 = [0.5_f64, 2.0_f64, 0.5_f64]
 self%cell_integrals_1 = self%cell_integrals_1 / 3.0_f64

 self%cell_integrals_0 = [1.0_f64,11.0_f64,11.0_f64,1.0_f64]
 self%cell_integrals_0 = self%cell_integrals_0 / 24.0_f64


 allocate( self%mass_line_0( self%spline_degree +1) )
 allocate( self%mass_line_1( self%spline_degree ) )
 call sll_s_spline_fem_mass_line ( self%spline_degree, self%mass_line_0 )
 call sll_s_spline_fem_mass_line ( self%spline_degree-1, self%mass_line_1 )

 ! Scale with dx
 self%mass_line_0 = self%mass_line_0 * self%delta_x
 self%mass_line_1 = self%mass_line_1 * self%delta_x


 SLL_ALLOCATE(self%rho_dofs(self%kernel_smoother_0%n_dofs), ierr)
 SLL_ALLOCATE(self%rho_dofs_local(self%kernel_smoother_0%n_dofs), ierr)
 SLL_ALLOCATE(self%efield_dofs_local(self%kernel_smoother_1%n_dofs), ierr)

 allocate(self%m1(self%spline_degree,self%kernel_smoother_1%n_dofs))
 allocate(self%m2(2*self%spline_degree,self%kernel_smoother_0%n_dofs))
 allocate(self%m4(self%spline_degree+1,self%kernel_smoother_0%n_dofs))
 allocate(self%m1_local(self%spline_degree,self%kernel_smoother_1%n_dofs))
 allocate(self%m2_local(2*self%spline_degree,self%kernel_smoother_0%n_dofs))
 allocate(self%m4_local(self%spline_degree+1,self%kernel_smoother_0%n_dofs))
 call self%linear_operator%create(self%n_dofs, self%spline_degree, self%mass_line_0,self%mass_line_1, self%m1, self%m2, self%m4)


 call self%linear_solver%create(self%linear_operator)
 if (present(solver_tolerance)) then
    self%linear_solver%rtol= solver_tolerance
    self%linear_solver%atol= solver_tolerance
 else
    self%linear_solver%rtol= 1d-12
    self%linear_solver%atol= 1d-12
 end if
 !self%linear_solver%verbose=.true.
end subroutine initialize_pic_vm_1d2v_ecsim

!---------------------------------------------------------------------------!
!> Destructor.
subroutine delete_pic_vm_1d2v_ecsim(self)
 class(sll_t_time_propagator_pic_vm_1d2v_ecsim), intent( inout ) :: self !< time splitting object 

 deallocate(self%j_dofs)
 deallocate(self%j_dofs_local)
 deallocate(self%m1)
 deallocate(self%m2)
 deallocate(self%m4)
 deallocate(self%m1_local)
 deallocate(self%m2_local)
 deallocate(self%m4_local)
 self%kernel_smoother_0 => null()
 self%kernel_smoother_1 => null()
 self%particle_group => null()
 self%efield_dofs => null()
 self%bfield_dofs => null()


 call self%linear_operator%free()
 call self%linear_solver%free()

end subroutine delete_pic_vm_1d2v_ecsim





end module sll_m_time_propagator_pic_vm_1d2v_ecsim
