!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on Hamiltonian splitting proposed by Crouseilles, Einkemmer, Faou for 1d2v Vlasov-Maxwell.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
!> Reference: Kraus, Kormann, Sonnendrücker, Morrison: GEMPIC: Geometric ElectroMagnetic Particle-In-Cell Methods
module sll_m_time_propagator_pic_vm_1d2v_cef
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

  use sll_m_maxwell_1d_base, only: &
       sll_c_maxwell_1d_base

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use mpi, only: &
       mpi_sum

  use sll_m_time_propagator_pic_vm_3d3v_cl_helper, only: &
       sll_p_boundary_particles_periodic, &
       sll_p_boundary_particles_singular, &
       sll_p_boundary_particles_reflection, &
       sll_p_boundary_particles_absorption

  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_1d2v_cef

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_1d2v_cef
     class(sll_c_maxwell_1d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_0  !< Kernel smoother (order p+1)
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_1  !< Kernel smoother (order p)
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

     sll_int32 :: spline_degree !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx !< Size of the domain
     sll_real64 :: x_min !< Lower bound for x domain
     sll_real64 :: x_max
     sll_real64 :: delta_x !< Grid spacing
     sll_int32  :: n_cells !< number of grid cells
     sll_int32  :: n_total0  !< number of Dofs for 0form
     sll_int32  :: n_total1 !< number of Dofs for 1form

     sll_int32 :: boundary_particles = 0 !< particle boundary conditions
     sll_int32 :: counter_left = 0 !< boundary counter
     sll_int32 :: counter_right = 0 !< boundary counter
     logical :: boundary = .false. !< true for non periodic boundary conditions
     sll_real64, pointer     :: rhob(:) => null() !< charge at the boundary

     logical :: electrostatic = .false. !< true for electrostatic simulation

     sll_real64, pointer     :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
     sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the magnetic field
     sll_real64, allocatable :: j_dofs(:,:)      !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:,:)!< MPI-processor local part of one component of \a j_dofs

   contains
     procedure :: operatorHp => operatorHp_pic_vm_1d2v  !> Operator for H_p part
     procedure :: operatorHE => operatorHE_pic_vm_1d2v  !> Operator for H_E part
     procedure :: operatorHB => operatorHB_pic_vm_1d2v  !> Operator for H_B part
     procedure :: lie_splitting => lie_splitting_pic_vm_1d2v !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_1d2v !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_1d2v !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_1d2v !> Initialize the type
     procedure :: free => delete_pic_vm_1d2v !> Finalization

  end type sll_t_time_propagator_pic_vm_1d2v_cef

contains


  !> Strang splitting
  subroutine strang_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_cef), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1, number_steps
       call self%operatorHB(0.5_f64*dt)
       call self%operatorHE(0.5_f64*dt)
       call self%operatorHp(dt)
       call self%operatorHE(0.5_f64*dt)
       call self%operatorHB(0.5_f64*dt)

    end do

  end subroutine strang_splitting_pic_vm_1d2v


  !> Lie splitting
  subroutine lie_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_cef), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1,number_steps
       call self%operatorHE(dt)
       call self%operatorHB(dt)
       call self%operatorHp(dt)
    end do

  end subroutine lie_splitting_pic_vm_1d2v


  !> Lie splitting (oposite ordering)
  subroutine lie_splitting_back_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_cef), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1,number_steps
       call self%operatorHp(dt)
       call self%operatorHB(dt)
       call self%operatorHE(dt)
    end do

  end subroutine lie_splitting_back_pic_vm_1d2v


  !---------------------------------------------------------------------------!
  !> Push Hp1: Equations to solve are
  !> \partial_t f + v_1 \partial_{x_1} f = 0    -> X_new = X_old + dt V_1
  !> V_new,2 = V_old,2 + \int_0 h V_old,1 B_old
  !> \partial_t E_1 = - \int v_1 f(t,x_1, v) dv -> E_{1,new} = E_{1,old} - \int \int v_1 f(t,x_1+s v_1,v) dv ds
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B = 0 => B_new = B_old 
  subroutine operatorHp_pic_vm_1d2v(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_cef), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: xnew(3), vi(3), wi(1), xi(3)

    self%j_dofs_local = 0.0_f64
    do i_sp = 1, self%particle_group%n_species
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          ! Then update particle position:  xnew = xold + dt * V
          xnew = xi + dt * vi

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part)
          call compute_particle_boundary_current( self, xi(1), xnew(1), vi(1:2), wi, dt )

          call self%particle_group%group(i_sp)%set_x(i_part, xnew)
       end do

    end do
    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
         self%n_total1, MPI_SUM, self%j_dofs(1:self%n_total1,1))
    ! Update the electric field.
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(1:self%n_total1,1), 1, self%efield_dofs(1:self%n_total1,1))

    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
         self%n_total0, MPI_SUM, self%j_dofs(:,2))
    ! Update the electric field.
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,2), 2, self%efield_dofs(:,2))

  end subroutine operatorHp_pic_vm_1d2v


  !> Helper function for operatorHp
  subroutine compute_particle_boundary_current( self, xold, xnew, vi, wi, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_cef), intent( inout ) :: self !< time splitting object 
    sll_real64,                                           intent( in    ) :: xold(1)
    sll_real64,                                           intent( inout ) :: xnew(1)
    sll_real64,                                           intent( inout ) :: vi(2)
    sll_real64,                                           intent( in    ) :: wi(1)
    sll_real64,                                           intent( in    ) :: dt
    !local variable
    sll_real64 :: xmid(1), vh, xbar

    if(xnew(1) < self%x_min .or. xnew(1) > self%x_max )then
       if(xnew(1) < self%x_min  )then
          xbar = self%x_min
          self%counter_left = self%counter_left+1
       else if(xnew(1) > self%x_max)then
          xbar = self%x_max
          self%counter_right = self%counter_right+1
       end if

       xmid = xbar
       vh = vi(2)/vi(1) 

       call self%kernel_smoother_1%add_current( xold, xmid, wi(1), self%j_dofs_local(1:self%n_total1,1))
       call self%kernel_smoother_0%add_current( xold, xmid, wi(1)*vh, self%j_dofs_local(:,2))

       select case(self%boundary_particles)
       case(sll_p_boundary_particles_reflection) 
          xnew = 2._f64*xbar-xnew
          vi(1) = -vi(1)
       case(sll_p_boundary_particles_absorption)
          call self%kernel_smoother_0%add_charge(xmid, wi(1), self%rhob)
       case( sll_p_boundary_particles_periodic)
          call self%kernel_smoother_0%add_charge(xmid, wi(1), self%rhob)
          xnew = self%x_min + modulo(xnew-self%x_min, self%Lx)
          xmid = self%x_max+self%x_min-xbar
          call self%kernel_smoother_0%add_charge(xmid, -wi(1), self%rhob)
       case default
          call self%kernel_smoother_0%add_charge(xmid, wi(1), self%rhob)
          xnew = self%x_min + modulo(xnew-self%x_min, self%Lx)
          xmid = self%x_max+self%x_min-xbar
          call self%kernel_smoother_0%add_charge(xmid, -wi(1), self%rhob)
       end select
       if (xnew(1) >= self%x_min .and. xnew(1) <= self%x_max) then
          vh = vi(2)/vi(1) 
          call self%kernel_smoother_1%add_current( xmid, xnew, wi(1), self%j_dofs_local(1:self%n_total1,1))
          call self%kernel_smoother_0%add_current( xmid, xnew, wi(1)*vh, self%j_dofs_local(:,2))
       end if
    else
       if ( abs(vi(1))> 1d-16 ) then
          ! Scale vi by weight to combine both factors for accumulation of integral over j
          call self%kernel_smoother_1%add_current( xold, xnew, wi(1), self%j_dofs_local(1:self%n_total1,1))
          call self%kernel_smoother_0%add_current( xold, xnew, wi(1)*vi(2)/vi(1), self%j_dofs_local(:,2))
       end if
    end if

  end subroutine compute_particle_boundary_current



  !---------------------------------------------------------------------------!
  !> Push H_E: Equations to be solved
  !> \partial_t f + E_1 \partial_{v_1} f + E_2 \partial_{v_2} f = 0 -> V_new = V_old + dt * E
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B + \partial_{x_1} E_2 = 0 => B_new = B_old - dt \partial_{x_1} E_2
  subroutine operatorHE_pic_vm_1d2v(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_cef), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: vi(3), xi(3)
    sll_real64 :: efield(2)
    sll_real64 :: qm

    do i_sp = 1, self%particle_group%n_species
       qm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part=1,self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          ! Evaluate efields at particle position
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          call self%kernel_smoother_1%evaluate &
               (xi(1), self%efield_dofs(1:self%n_total1,1), efield(1))
          call self%kernel_smoother_0%evaluate &
               (xi(1), self%efield_dofs(:,2), efield(2))
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          ! V_new = V_old + dt * E
          vi(1:2) = vi(1:2) + dt* qm * efield

          call self%particle_group%group(i_sp)%set_v(i_part, vi)
       end do
    end do

    if(self%electrostatic .eqv. .false.) then
       ! Update bfield
       call self%maxwell_solver%compute_B_from_E( &
            dt, self%efield_dofs(:,2), self%bfield_dofs)
    end if

  end subroutine operatorHE_pic_vm_1d2v


  !---------------------------------------------------------------------------!
  !> Push H_B: Equations to be solved
  !> V_new = V_old
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old}
  !> \partial_t E_2 = - \partial_{x_1} B -> E_{2,new} = E_{2,old}-dt*\partial_{x_1} B
  !> \partial_t B = 0 -> B_new = B_old
  subroutine operatorHB_pic_vm_1d2v(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_cef), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32  :: i_part, i_sp
    sll_real64 :: qmdt
    sll_real64 :: vi(3), v_new(3), xi(3)
    sll_real64 :: bfield, cs, sn

    ! VxB part
    do i_sp = 1, self%particle_group%n_species
       qmdt = self%particle_group%group(i_sp)%species%q_over_m()*dt;
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          call self%kernel_smoother_1%evaluate &
               (xi(1), self%bfield_dofs, bfield)

          bfield = qmdt*bfield

          cs = cos(bfield)
          sn = sin(bfield)

          v_new(1) = cs * vi(1) + sn * vi(2)
          v_new(2) = -sn* vi(1) + cs * vi(2)

          call self%particle_group%group(i_sp)%set_v( i_part, v_new )
       end do
    end do

    if(self%electrostatic .eqv. .false.) then
       ! Update efield2
       call self%maxwell_solver%compute_E_from_B(&
            dt, self%bfield_dofs, self%efield_dofs(:,2))
    end if

  end subroutine operatorHB_pic_vm_1d2v


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
       boundary_particles, &
       electrostatic, &
       rhob)
    class(sll_t_time_propagator_pic_vm_1d2v_cef), intent(out) :: self !< time splitting object 
    class(sll_c_maxwell_1d_base), target,          intent(in)  :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,      intent(in)  :: particle_group !< Particle group
    sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
    sll_int32, optional,                             intent( in )  :: boundary_particles !< particle boundary conditions
    logical, optional,                               intent( in )  :: electrostatic !< true for electrostatic simulation
    sll_real64, optional, target,                  intent( in ) :: rhob(:) !< charge at the boundary
    !local variables
    sll_int32 :: ierr

    self%maxwell_solver => maxwell_solver
    self%kernel_smoother_0 => kernel_smoother_0
    self%kernel_smoother_1 => kernel_smoother_1
    self%particle_group => particle_group
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs

    self%spline_degree = self%kernel_smoother_0%spline_degree
    self%x_min = x_min
    self%x_max = x_min + Lx
    self%Lx = Lx
    self%delta_x = self%Lx/real(self%maxwell_solver%n_cells, f64)

    self%n_cells = self%maxwell_solver%n_cells
    self%n_total0 = self%maxwell_solver%n_dofs0
    self%n_total1 = self%maxwell_solver%n_dofs1

    ! Check that n_dofs is the same for both kernel smoothers.
    SLL_ASSERT( self%kernel_smoother_0%n_cells == self%kernel_smoother_1%n_cells )
    SLL_ASSERT( self%kernel_smoother_0%n_cells == self%maxwell_solver%n_cells )

    SLL_ALLOCATE(self%j_dofs(self%n_total0, 2), ierr)
    SLL_ALLOCATE(self%j_dofs_local(self%n_total0, 2), ierr)

    if( self%n_cells+self%spline_degree == self%maxwell_solver%n_dofs0   ) then
       self%boundary = .true.
       self%boundary_particles = boundary_particles
    end if

    if( present(electrostatic) )then
       self%electrostatic = electrostatic
    end if

    if( present(rhob) )then
       self%rhob => rhob
    end if

  end subroutine initialize_pic_vm_1d2v

  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_1d2v(self)
    class(sll_t_time_propagator_pic_vm_1d2v_cef), intent( inout ) :: self !< time splitting object

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
    self%efield_dofs => null()
    self%bfield_dofs => null()

  end subroutine delete_pic_vm_1d2v


end module sll_m_time_propagator_pic_vm_1d2v_cef
