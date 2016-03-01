!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on Hamiltonian splitting for 1d2v Vlasov-Poisson.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_hamiltonian_splitting_pic_vm_1d2v
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_o_collective_allreduce, &
    sll_v_world_collective

  use sll_m_hamiltonian_splitting_base, only: &
    sll_c_hamiltonian_splitting_base

  use sll_m_kernel_smoother_base, only: &
    sll_c_kernel_smoother

  use sll_m_maxwell_1d_base, only: &
    sll_c_maxwell_1d_base

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_mpi, only: &
    mpi_sum

  implicit none

  public :: &
    sll_s_new_hamiltonian_splitting_pic_vm_1d2v, &
    sll_s_new_hamiltonian_splitting_pic_vm_1d2v_ptr, &
    sll_t_hamiltonian_splitting_pic_vm_1d2v

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type, extends(sll_c_hamiltonian_splitting_base) :: sll_t_hamiltonian_splitting_pic_vm_1d2v
     class(sll_c_maxwell_1d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_kernel_smoother), pointer :: kernel_smoother_0  !< Kernel smoother (order p+1)
     class(sll_c_kernel_smoother), pointer :: kernel_smoother_1  !< Kernel smoother (order p)
     class(sll_c_particle_group_base), pointer  :: particle_group    !< Particle group

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

   contains
     procedure :: operatorHp1 => operatorHp1_pic_vm_1d2v  !> Operator for H_p1 part
     procedure :: operatorHp2 => operatorHp2_pic_vm_1d2v  !> Operator for H_p2 part
     procedure :: operatorHE => operatorHE_pic_vm_1d2v  !> Operator for H_E part
     procedure :: operatorHB => operatorHB_pic_vm_1d2v  !> Operator for H_B part
     procedure :: lie_splitting => lie_splitting_pic_vm_1d2v !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_1d2v !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_1d2v !> Initialize the type
     procedure :: free => delete_pic_vm_1d2v !> Finalization

  end type sll_t_hamiltonian_splitting_pic_vm_1d2v

contains

  !> Strang splitting
  subroutine strang_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1, number_steps
       call self%operatorHp1(0.5_f64*dt)
       call self%operatorHp2(0.5_f64*dt)
       call self%operatorHE(0.5_f64*dt)
       call self%operatorHB(dt)
       call self%operatorHE(0.5_f64*dt)
       call self%operatorHp2(0.5_f64*dt)
       call self%operatorHp1(0.5_f64*dt)
    end do

  end subroutine strang_splitting_pic_vm_1d2v

  !> Lie splitting
  subroutine lie_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1, number_steps
       call self%operatorHp1(dt)
       call self%operatorHp2(dt)
       call self%operatorHE(dt)
       call self%operatorHB(dt)
    end do

  end subroutine lie_splitting_pic_vm_1d2v

  !---------------------------------------------------------------------------!
  !> Push Hp1: Equations to solve are
  !> \partial_t f + v_1 \partial_{x_1} f = 0    -> X_new = X_old + dt V_1
  !> V_new,2 = V_old,2 + \int_0 h V_old,1 B_old
  !> \partial_t E_1 = - \int v_1 f(t,x_1, v) dv -> E_{1,new} = E_{1,old} - \int \int v_1 f(t,x_1+s v_1,v) dv ds
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B = 0 => B_new = B_old 
  subroutine operatorHp1_pic_vm_1d2v(self, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: x_new(3), vi(3), wi(1), x_old(3)
    sll_int32  :: n_cells
    sll_real64 :: qoverm


    n_cells = self%kernel_smoother_0%n_dofs

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
    do i_part=1,self%particle_group%n_particles  
       ! Read out particle position and velocity
       x_old = self%particle_group%get_x(i_part)
       vi = self%particle_group%get_v(i_part)

       ! Then update particle position:  X_new = X_old + dt * V
       x_new = x_old + dt * vi

       ! Get charge for accumulation of j
       wi = self%particle_group%get_charge(i_part)
       qoverm = self%particle_group%species%q_over_m();

       call self%kernel_smoother_1%add_current_update_v( x_old, x_new, wi(1), qoverm, &
            self%bfield_dofs, vi, self%j_dofs_local(:,1))
      
       x_new(1) = modulo(x_new(1), self%Lx)
       call self%particle_group%set_x(i_part, x_new)
       call self%particle_group%set_v(i_part, vi)


    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
         n_cells, MPI_SUM, self%j_dofs(:,1))

    ! Update the electric field.
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,1), 1, self%efield_dofs(:,1))



 end subroutine operatorHp1_pic_vm_1d2v




 !---------------------------------------------------------------------------!
  !> Push Hp2: Equations to solve are
  !> X_new = X_old
  !> V_new,1 = V_old,1 + \int_0 h V_old,2 B_old
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = - \int v_2 f(t,x_1, v) dv -> E_{2,new} = E_{2,old} - \int \int v_2 f(t,x_1+s v_1,v) dv ds
  !> \partial_t B = 0 => B_new = B_old
  subroutine operatorHp2_pic_vm_1d2v(self, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32  :: i_part, n_cells
    sll_real64 :: vi(3), xi(3), wi(1)
    sll_real64 :: bfield
    sll_real64 :: qm
    
    n_cells = self%kernel_smoother_0%n_dofs

    self%j_dofs_local = 0.0_f64

    qm = self%particle_group%species%q_over_m();
    ! Update v_1
    do i_part=1,self%particle_group%n_particles
       ! Evaluate bfield at particle position (splines of order p)
       xi = self%particle_group%get_x(i_part)
       call self%kernel_smoother_1%evaluate &
            (xi(1), self%bfield_dofs, bfield)
       vi = self%particle_group%get_v(i_part)
       vi(1) = vi(1) + dt*qm*vi(2)*bfield
       call self%particle_group%set_v(i_part, vi)

       xi = self%particle_group%get_x(i_part)

       ! Scale vi by weight to combine both factors for accumulation of integral over j
       wi = self%particle_group%get_charge(i_part)*vi(2)

       call self%kernel_smoother_0%add_charge(xi(1:1), wi(1), self%j_dofs_local(:,2)) 

    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
         n_cells, MPI_SUM, self%j_dofs(:,2))

    ! Update the electric field. Also, we still need to scale with 1/Lx ! TODO: Which scaling?
    self%j_dofs(:,2) = self%j_dofs(:,2)*dt!/self%Lx

    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,2), 2, self%efield_dofs(:,2))
    

  end subroutine operatorHp2_pic_vm_1d2v
  
  !---------------------------------------------------------------------------!
  !> Push H_E: Equations to be solved
  !> \partial_t f + E_1 \partial_{v_1} f + E_2 \partial_{v_2} f = 0 -> V_new = V_old + dt * E
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B + \partial_{x_1} E_2 = 0 => B_new = B_old - dt \partial_{x_1} E_2
  subroutine operatorHE_pic_vm_1d2v(self, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: v_new(3), xi(3)
    sll_real64 :: efield(2)
    sll_real64 :: qm


    qm = self%particle_group%species%q_over_m();
    ! V_new = V_old + dt * E
    do i_part=1,self%particle_group%n_particles
       v_new = self%particle_group%get_v(i_part)
       ! Evaluate efields at particle position
       xi = self%particle_group%get_x(i_part)
       call self%kernel_smoother_1%evaluate &
            (xi(1), self%efield_dofs(:,1), efield(1))
       call self%kernel_smoother_0%evaluate &
            (xi(1), self%efield_dofs(:,2), efield(2))
       v_new = self%particle_group%get_v(i_part)
       v_new(1:2) = v_new(1:2) + dt* qm * efield
       call self%particle_group%set_v(i_part, v_new)
    end do
    
    ! Update bfield
    call self%maxwell_solver%compute_B_from_E( &
         dt, self%efield_dofs(:,2), self%bfield_dofs)
    

  end subroutine operatorHE_pic_vm_1d2v
  

  !---------------------------------------------------------------------------!
  !> Push H_B: Equations to be solved
  !> V_new = V_old
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old}
  !> \partial_t E_2 = - \partial_{x_1} B -> E_{2,new} = E_{2,old}-dt*\partial_{x_1} B
  !> \partial_t B = 0 -> B_new = B_old
  subroutine operatorHB_pic_vm_1d2v(self, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    
    ! Update efield2
    call self%maxwell_solver%compute_E_from_B(&
         dt, self%bfield_dofs, self%efield_dofs(:,2))
    

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
       Lx) 
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(out) :: self !< time splitting object 
    class(sll_c_maxwell_1d_base), pointer,          intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_kernel_smoother), pointer,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_kernel_smoother), pointer,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_c_particle_group_base), pointer,      intent(in)  :: particle_group !< Particle group
    sll_real64, pointer,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, pointer,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.

    !local variables
    sll_int32 :: ierr

    self%maxwell_solver => maxwell_solver
    self%kernel_smoother_0 => kernel_smoother_0
    self%kernel_smoother_1 => kernel_smoother_1
    self%particle_group => particle_group
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs

    ! Check that n_dofs is the same for both kernel smoothers.
    SLL_ASSERT( self%kernel_smoother_0%n_dofs == self%kernel_smoother_1%n_dofs )

    SLL_ALLOCATE(self%j_dofs(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%j_dofs_local(self%kernel_smoother_0%n_dofs,2), ierr)

    self%spline_degree = 3
    self%x_min = x_min
    self%Lx = Lx
    self%delta_x = self%Lx/self%kernel_smoother_1%n_dofs
    
    self%cell_integrals_1 = [0.5_f64, 2.0_f64, 0.5_f64]
    self%cell_integrals_1 = self%cell_integrals_1 / 3.0_f64

    self%cell_integrals_0 = [1.0_f64,11.0_f64,11.0_f64,1.0_f64]
    self%cell_integrals_0 = self%cell_integrals_0 / 24.0_f64

  end subroutine initialize_pic_vm_1d2v

  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_1d2v(self)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent( inout ) :: self !< time splitting object 

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
  subroutine sll_s_new_hamiltonian_splitting_pic_vm_1d2v(&
       splitting, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx) 
    class(sll_c_hamiltonian_splitting_base), allocatable, intent(out) :: splitting !< time splitting object 
    class(sll_c_maxwell_1d_base), pointer,                intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_kernel_smoother), pointer,                intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_kernel_smoother), pointer,                intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_c_particle_group_base),pointer,             intent(in)  :: particle_group !< Particle group
    sll_real64, pointer,                                  intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, pointer,                                  intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                           intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                           intent(in)  :: Lx !< Length of the domain in x direction.

    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE(sll_t_hamiltonian_splitting_pic_vm_1d2v :: splitting, ierr)

    select type (splitting)
    type is ( sll_t_hamiltonian_splitting_pic_vm_1d2v )
       call splitting%init(&
            maxwell_solver, &
            kernel_smoother_0, &
            kernel_smoother_1, &
            particle_group, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx)
    end select

  end subroutine sll_s_new_hamiltonian_splitting_pic_vm_1d2v

  !---------------------------------------------------------------------------!
  !> Constructor for pointer abstract type.
  subroutine sll_s_new_hamiltonian_splitting_pic_vm_1d2v_ptr(&
       splitting, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx) 
    class(sll_c_hamiltonian_splitting_base), pointer, intent(out) :: splitting !< time splitting object 
    class(sll_c_maxwell_1d_base), pointer,            intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_kernel_smoother), pointer,            intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_kernel_smoother), pointer,            intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_c_particle_group_base),pointer,         intent(in)  :: particle_group !< Particle group
    sll_real64, pointer,                              intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, pointer,                              intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                       intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                       intent(in)  :: Lx !< Length of the domain in x direction.

    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE(sll_t_hamiltonian_splitting_pic_vm_1d2v :: splitting, ierr)

    select type (splitting)
    type is ( sll_t_hamiltonian_splitting_pic_vm_1d2v )
       call splitting%init(&
            maxwell_solver, &
            kernel_smoother_0, &
            kernel_smoother_1, &
            particle_group, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx)
    end select

  end subroutine sll_s_new_hamiltonian_splitting_pic_vm_1d2v_ptr

end module sll_m_hamiltonian_splitting_pic_vm_1d2v
