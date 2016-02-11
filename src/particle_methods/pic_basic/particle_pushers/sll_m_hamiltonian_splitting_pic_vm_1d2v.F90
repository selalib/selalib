!> @ingroup particle_pushers
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
     class(sll_c_maxwell_1d_base), pointer  :: maxwell_solver      !< Maxwell solver
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
     procedure :: operatorHp1_pic_vm_1d2v_prim !> Alternative implementation of operator Hp1 using primitive function (only for degree 3)

     procedure :: init => initialize_pic_vm_1d2v !> Initialize the type
     procedure :: free => delete_pic_vm_1d2v !> Finalization

  end type sll_t_hamiltonian_splitting_pic_vm_1d2v

contains

  !> Strang splitting
  subroutine strang_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v) :: self !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step
    sll_int32, intent(in)  :: number_steps !< number of time steps

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
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v) :: self !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step
    sll_int32, intent(in)  :: number_steps !< number of time steps

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
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: x_new(3), vi(3), wi(1), x_old(3)
    sll_int32  :: n_cells
    sll_real64 :: qoverm
    !sll_real64 :: xi
    !sll_real64 :: r_new, r_old, qoverm, bfield
    !sll_int32 :: index_new, index_old, ind
    !sll_real64 :: primitive_1(3)
    !sll_real64 :: jnorm

    !sll_real64 :: efield_test(self%kernel_smoother_0%n_dofs), rho0(self%kernel_smoother_0%n_dofs), jk(self%kernel_smoother_0%n_dofs)

    n_cells = self%kernel_smoother_0%n_dofs

    ! Here we have to accumulate j and integrate over the time interval.
    ! At each k=1,...,n_grid, we have for s \in [0,dt]:
    ! j_k(s) =  \sum_{i=1,..,N_p} q_i N((x_k+sv_{1,k}-x_i)/h) v_k,
    ! where h is the grid spacing and N the normalized B-spline
    ! In order to accumulate the integrated j, we normalize the values of x to the grid spacing, calling them y, we have
    ! j_k(s) = \sum_{i=1,..,N_p} q_i N(y_k+s/h v_{1,k}-y_i) v_k.
    ! Now, we want the integral 
    ! \int_{0..dt} j_k(s) d s = \sum_{i=1,..,N_p} q_i v_k \int_{0..dt} N(y_k+s/h v_{1,k}-y_i) ds =  \sum_{i=1,..,N_p} q_i v_k  \int_{0..dt}  N(y_k + w v_{1,k}-y_i) dw



!!$    ! Test if efield_dofs(:,1) is the same when computed from Poisson
!!$    self%j_dofs_local(:,1) = 0.0_f64
!!$    do i_part=1,self%particle_group%n_particles
!!$       x_new = self%particle_group%get_x(i_part)
!!$       wi = self%particle_group%get_charge(i_part)
!!$       call self%kernel_smoother_0%add_charge(x_new(1), wi(1), self%j_dofs_local(:,1))
!!$    end do
!!$    self%j_dofs(:,1) = 0.0_f64
!!$    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
!!$         n_cells, MPI_SUM, self%j_dofs(:,1))
!!$    call self%maxwell_solver%compute_E_from_rho(efield_test,&
!!$         self%j_dofs(:,1))
!!$    !print*, 'Error efield: ',  maxval(abs(efield_test - self%efield_dofs(:,1))), maxval(abs(self%efield_dofs(:,1)))
!!$    !print*, self%j_dofs(:,1)
!!$    !print*, 'rho0########'
!!$    !print*, efield_test
!!$    !print*, '--------'
!!$    !print*, self%efield_dofs(:,1)

    self%j_dofs_local = 0.0_f64

    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
    ! Then update particle position:  X_new = X_old + dt * V
    do i_part=1,self%particle_group%n_particles  
       ! Read out particle position and velocity
       x_old = self%particle_group%get_x(i_part)
       vi = self%particle_group%get_v(i_part)

       ! Then update particle position:  X_new = X_old + dt * V
       x_new = x_old + dt * vi!modulo(x_old + dt * vi, self%Lx)
       !call self%particle_group%set_x(i_part, x_new)

       ! Get charge for accumulation of j
       wi = self%particle_group%get_charge(i_part)
       qoverm = self%particle_group%species%q_over_m();
       !print*, dt, vi(1)

       call self%kernel_smoother_1%add_current_update_v( x_old, x_new, wi(1), qoverm, &
            self%bfield_dofs, vi, self%j_dofs_local(:,1))
       !call self%kernel_smoother_1%add_current_update_v( x_old, x_old, wi(1), qoverm, &
       !     self%bfield_dofs, vi, self%j_dofs_local(:,1))
       !call self%kernel_smoother_1%evaluate &
       !     (x_old(1), self%bfield_dofs, bfield)
       !vi(2) = vi(2) - dt*qoverm*vi(1)*bfield
       !wi = wi*vi(1)
       !call self%kernel_smoother_1%add_charge(x_old(1:1), wi(1), self%j_dofs_local(:,1))
       !call self%kernel_smoother_1%add_charge(x_new(1:1), wi(1), self%j_dofs_local(:,1))

       x_new(1) = modulo(x_new(1), self%Lx)
       call self%particle_group%set_x(i_part, x_new)
       call self%particle_group%set_v(i_part, vi)

       !print*, x_old(1), x_new(1), wi, qoverm
       !print*, self%j_dofs_local(:,1)
       !print*, '--------'

    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
         n_cells, MPI_SUM, self%j_dofs(:,1))

    !jnorm = sum(self%j_dofs(:,1))/n_cells
    !print*, 'sum j', jnorm
    !self%j_dofs(:,1) = self%j_dofs(:,1)-jnorm

    ! Update the electric field. Also, we still need to scale with 1/Lx
    !self%j_dofs(:,1) = self%j_dofs(:,1)*dt/2!self%delta_x!/self%Lx
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,1), 1, self%efield_dofs(:,1))
!!$    print*, 'ea', self%efield_dofs(:,1)
!!$    print*, 'j++++++'
!!$
!!$    ! Test if efield_dofs(:,1) is the same when computed from Poisson
!!$    self%j_dofs_local(:,1) = 0.0_f64
!!$    do i_part=1,self%particle_group%n_particles
!!$       x_new = self%particle_group%get_x(i_part)
!!$       wi = self%particle_group%get_charge(i_part)
!!$       call self%kernel_smoother_0%add_charge(x_new(1), wi(1), self%j_dofs_local(:,1))
!!$       !print*, 'a', wi, x_new
!!$       !print*, self%j_dofs_local(:,1)
!!$    end do
!!$    self%j_dofs(:,1) = 0.0_f64
!!$    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
!!$         n_cells, MPI_SUM, self%j_dofs(:,1))
!!$    call self%maxwell_solver%compute_E_from_rho(efield_test,&
!!$         self%j_dofs(:,1))
!!$    !print*, self%efield_dofs(:,1)
!!$    print*, 'Error efield: ',  maxval(abs(efield_test - self%efield_dofs(:,1))), maxval(abs(self%efield_dofs(:,1)))
!!$    print*, self%j_dofs(:,1)
!!$    print*, 'rho1++++++'
!!$    print*, 'ep', efield_test


 end subroutine operatorHp1_pic_vm_1d2v

!!$ ! TODO: Self is hard coded for quadratic, cubic splines. Make general.
!!$ subroutine update_jv(self, lower, upper, index, weight, sign, vi)
!!$   class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: self !< time splitting object 
!!$   sll_real64, intent(in) :: lower
!!$   sll_real64, intent(in) :: upper
!!$   sll_int32,  intent(in) :: index
!!$   sll_real64, intent(in) :: weight(1)
!!$   sll_real64, intent(in) :: sign
!!$   sll_real64, intent(inout) :: vi
!!$
!!$   !Local variables
!!$   sll_real64 :: m, c, y1, y2, fy(self%spline_degree+1)
!!$   sll_int32  :: ind, i_grid, i_mod, n_cells
!!$   sll_real64 :: qm
!!$
!!$
!!$   qm = self%particle_group%species%q_over_m();
!!$   n_cells = self%kernel_smoother_0%n_dofs
!!$
!!$
!!$   m = 0.5_f64*(upper-lower)
!!$   c = 0.5_f64*(upper+lower)
!!$   y1 = -m/sqrt(3.0_f64)+c
!!$   y2 = m/sqrt(3.0_f64) +c
!!$   fy = sign*m*self%delta_x*&
!!$        (sll_f_uniform_b_splines_at_x(self%spline_degree-1, y1) +&
!!$        sll_f_uniform_b_splines_at_x(self%spline_degree-1, y2))
!!$
!!$
!!$   ind = 1
!!$   do i_grid = index - self%spline_degree + 1, index
!!$      i_mod = modulo(i_grid, n_cells ) + 1
!!$      self%j_dofs_local(i_mod,1) = self%j_dofs_local(i_mod,1) + &
!!$           weight(1)*fy(ind)
!!$      vi = vi - qm* fy(ind)*self%bfield_dofs(i_mod)
!!$      ind = ind + 1
!!$   end do
!!$
!!$
!!$ end subroutine update_jv

 !> Alternative implementation of the Hp1 operator based on the primal function of the spline
 subroutine operatorHp1_pic_vm_1d2v_prim(self, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: self !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part, i_grid, i_mod, ind, j
    sll_real64 :: xi, x_new(3), vi(3), wi(1)
    sll_int32  :: n_cells
    sll_real64 :: r_new, r_old
    sll_int32 :: index_new, index_old
    sll_real64 :: primitive_1(3)
    sll_real64 :: qm

    self%j_dofs_local = 0.0_f64
    n_cells = self%kernel_smoother_0%n_dofs    
    qm = self%particle_group%species%q_over_m()

    ! Here we have to accumulate j and integrate over the time interval.
    ! At each k=1,...,n_grid, we hav for s \in [0,dt]:
    ! j_k(s) =  \sum_{i=1,..,N_p} q_i N((x_k+sv_{1,k}-x_i)/h) v_k,
    ! where h is the grid spacing and N the normalized B-spline
    ! In order to accumulate the integrated j, we normalize the values of x to the grid spacing, calling them y, we have
    ! j_k(s) = \sum_{i=1,..,N_p} q_i N(y_k+s/h v_{1,k}-y_i) v_k.
    ! Now, we want the integral 
    ! \int_{0..dt} j_k(s) d s = \sum_{i=1,..,N_p} q_i v_k \int_{0..dt} N(y_k+s/h v_{1,k}-y_i) ds =  \sum_{i=1,..,N_p} q_i v_k  \int_{0..dt}  N(y_k + w v_{1,k}-y_i) dw


    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
    ! Then update particle position:  X_new = X_old + dt * V
    do i_part=1,self%particle_group%n_particles  
       ! Read out particle position and velocity
       x_new = self%particle_group%get_x(i_part)
       vi = self%particle_group%get_v(i_part)
       ! Compute index_old, the index of the last DoF on the grid the particle contributes to, and r_old, its position (normalized to cell size one).
       xi = (x_new(1) - self%x_min) /&
            self%delta_x
       index_old = floor(xi)
       r_old = xi - real(index_old,f64)

       ! Then update particle position:  X_new = X_old + dt * V
       x_new = modulo(self%particle_group%get_x(i_part) + dt * vi, self%Lx)
       call self%particle_group%set_x(i_part, x_new)

       ! Compute the new box index index_new and normalized position r_old.
       xi = (x_new(1) - self%x_min) /&
            self%delta_x
       index_new = floor(xi)
       r_new = xi - real(index_new ,f64) 

       ! Scale vi by weight to combine both factors for accumulation of integral over j
       wi = self%particle_group%get_charge(i_part)

       ! Compute the primitives at r_old for each interval
       primitive_1 = -primitive_uniform_quadratic_b_spline_at_x( r_old )
       ! If we are integrating to the right, we need to also use the value of the primitive at the right interval bound. If larger, we would need the lower but it is zero by our normalization.   
       if (index_old < index_new) then
          primitive_1 = primitive_1 + self%cell_integrals_1
       end if
       ! Now, we loop through the DoFs and add the contributions.
       ind = 1
       do i_grid = index_old - self%spline_degree + 1, index_old
          i_mod = modulo(i_grid, n_cells ) + 1
          self%j_dofs_local(i_mod,1) = self%j_dofs_local(i_mod,1) + &
               primitive_1(ind)*wi(1)
          vi(2) = vi(2) - qm*primitive_1(ind)*self%bfield_dofs(i_mod)*self%delta_x
          ind = ind + 1
       end do      
       ! Now contribution from r_new in the same way but different sign
       primitive_1 = primitive_uniform_quadratic_b_spline_at_x( r_new )
       if (index_old > index_new) then
          primitive_1 = primitive_1 - self%cell_integrals_1
       end if 

       ind = 1
       do i_grid = index_new - self%spline_degree + 1, index_new
          i_mod = modulo(i_grid, n_cells ) + 1
          self%j_dofs_local(i_mod,1) = self%j_dofs_local(i_mod,1) + &
               primitive_1(ind)*wi(1)
          vi(2) = vi(2) - qm*primitive_1(ind)*self%bfield_dofs(i_mod)*self%delta_x
          ind = ind + 1
       end do
       ! If |index_new - index_old|<2, we are done. Otherwise, we have to account for the piecewise definition of the spline and add the contributions from the intervals inbetween. These are always integrals over the whole integrals.
       ! First if index_old<index_new-1: Add the whole integrals.
       do j = index_old+1, index_new-1
          ind = 1
          do i_grid = j - self%spline_degree + 1, j
             i_mod = modulo(i_grid, n_cells ) + 1
             self%j_dofs_local(i_mod,1) = self%j_dofs_local(i_mod,1) + &
                  self%cell_integrals_1(ind)*wi(1)
             vi(2) = vi(2) - &
                  qm*self%cell_integrals_1(ind)*self%bfield_dofs(i_mod)*self%delta_x
             ind = ind + 1
          end do
       end do
       ! Then if index_old>index_new-1: Subtract the whole integrals.
       do j = index_new+1, index_old-1
          ind = 1
          do i_grid = j - self%spline_degree + 1, j
             i_mod = modulo(i_grid, n_cells ) + 1
             self%j_dofs_local(i_mod,1) = self%j_dofs_local(i_mod,1) - &
                  self%cell_integrals_1(ind)*wi(1)
             vi(2) = vi(2) + qm*self%cell_integrals_1(ind)*self%bfield_dofs(i_mod)*self%delta_x
             ind = ind + 1
          end do
       end do 

       call self%particle_group%set_v(i_part, vi)

    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
         n_cells, MPI_SUM, self%j_dofs(:,1))

    ! Update the electric field. Also, we still need to scale with 1/Lx
    self%j_dofs(:,1) = self%j_dofs(:,1)*self%delta_x
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,1), 1, self%efield_dofs(:,1))


 end subroutine operatorHp1_pic_vm_1d2v_prim






 !---------------------------------------------------------------------------!
  !> Push Hp2: Equations to solve are
  !> X_new = X_old
  !> V_new,1 = V_old,1 + \int_0 h V_old,2 B_old
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = - \int v_2 f(t,x_1, v) dv -> E_{2,new} = E_{2,old} - \int \int v_2 f(t,x_1+s v_1,v) dv ds
  !> \partial_t B = 0 => B_new = B_old
  subroutine operatorHp2_pic_vm_1d2v(self, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: self !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32  :: i_part, n_cells
    sll_real64 :: vi(3), xi(3), wi(1)
    !sll_real64 :: x_box, r_box
    !sll_real64 :: values_0(4)
    sll_real64 :: bfield
    sll_real64 :: qm
    !sll_real64 :: jnorm
    
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

    !jnorm = sum(self%j_dofs(:,2))/n_cells
    !print*, 'sum j2', jnorm
    !self%j_dofs(:,2) = self%j_dofs(:,2)-jnorm

    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,2), 2, self%efield_dofs(:,2))
    !!self%efield_dofs(:,2) = self%efield_dofs(:,2) - self%j_dofs(:,2)/self%Lx*dt
    

  end subroutine operatorHp2_pic_vm_1d2v
  
  !---------------------------------------------------------------------------!
  !> Push H_E: Equations to be solved
  !> \partial_t f + E_1 \partial_{v_1} f + E_2 \partial_{v_2} f = 0 -> V_new = V_old + dt * E
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B + \partial_{x_1} E_2 = 0 => B_new = B_old - dt \partial_{x_1} E_2
  subroutine operatorHE_pic_vm_1d2v(self, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: self !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

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
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    !sll_int32 :: i_part
    !sll_real64 :: v_old(3), v_new(3)

    
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
    class(sll_c_maxwell_1d_base), pointer, intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_kernel_smoother), pointer, intent(in) :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_kernel_smoother), pointer, intent(in) :: kernel_smoother_1  !< Kernel smoother
    class(sll_c_particle_group_base), pointer, intent(in) :: particle_group !< Particle group
    sll_real64, pointer, intent(in) :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, pointer, intent(in) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64, intent(in) :: x_min !< Lower bound of x domain
    sll_real64, intent(in) :: Lx !< Length of the domain in x direction.

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
    deallocate(self%efield_dofs)
    deallocate(self%bfield_dofs)

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
    class(sll_c_maxwell_1d_base), pointer, intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_kernel_smoother), pointer, intent(in) :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_kernel_smoother), pointer, intent(in) :: kernel_smoother_1  !< Kernel smoother
    class(sll_c_particle_group_base),pointer, intent(in) :: particle_group !< Particle group
    sll_real64, pointer, intent(in) :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, pointer, intent(in) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64, intent(in) :: x_min !< Lower bound of x domain
    sll_real64, intent(in) :: Lx !< Length of the domain in x direction.

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
    class(sll_c_maxwell_1d_base), pointer, intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_kernel_smoother), pointer, intent(in) :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_kernel_smoother), pointer, intent(in) :: kernel_smoother_1  !< Kernel smoother
    class(sll_c_particle_group_base),pointer, intent(in) :: particle_group !< Particle group
    sll_real64, pointer, intent(in) :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, pointer, intent(in) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64, intent(in) :: x_min !< Lower bound of x domain
    sll_real64, intent(in) :: Lx !< Length of the domain in x direction.

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



  !> Compute the primitive of the cubic B-spline in each intervall at x. Primitive function normalized such that it is 0 at x=0. Analogon to uniform_b_spline_at_x in arbitrary degree splines for primitive, but specific for cubic.
   function primitive_uniform_cubic_b_spline_at_x( x) result(primitive)
     sll_real64, intent(in)  :: x !< position where to evaluate the primitive
     sll_real64 :: primitive(4) !< value of the primitive for each of the four intervals.
 
     sll_real64 :: xx(3)
 
     xx(1) = x**2
     xx(2) = x*xx(1)
     xx(3) = x*xx(2)
 
     primitive(4) = xx(3)/24.0_f64
     primitive(3) = (x + 1.5_f64*xx(1) + xx(2) - 0.75_f64* xx(3))/6.0_f64
     primitive(2) = (4.0_f64*x - 2.0_f64* xx(2) + 0.75_f64* xx(3))/6.0_f64
     primitive(1) = (1.0_f64 - (1.0_f64-x)**4)/24.0_f64
 
   end function primitive_uniform_cubic_b_spline_at_x
 
   !> Compute the primitive of the quadratic B-spline in each intervall at x. Primitive function normalized such that it is 0 at x=0. Analogon to uniform_b_spline_at_x in arbitrary degree splines for primitive, but specific for quadratic.
   function primitive_uniform_quadratic_b_spline_at_x( x) result(primitive)
     sll_real64, intent(in)  :: x !< position where to evaluate the primitive
     sll_real64 :: primitive(3) !< value of the primitive for each of the three intervals.
 
     sll_real64 :: xx(2)

     xx(1) = x**2
     xx(2) = x*xx(1)
     
     primitive(3) = xx(2)/6.0_f64
     primitive(2) = (x + xx(1))*0.5_f64-xx(2)/3.0_f64
     primitive(1) = (1.0_f64 - (1.0_f64-x)**3)/6.0_f64
 
   end function primitive_uniform_quadratic_b_spline_at_x
 

end module sll_m_hamiltonian_splitting_pic_vm_1d2v
