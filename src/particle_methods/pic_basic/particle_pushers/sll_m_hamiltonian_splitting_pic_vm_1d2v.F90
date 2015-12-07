!> @ingroup particle_pushers
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on Hamiltonian splitting for 1d2v Vlasov-Poisson.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_hamiltonian_splitting_pic_vm_1d2v
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_arbitrary_degree_splines, only: &
    uniform_b_splines_at_x

  use sll_m_collective, only: &
    sll_collective_allreduce, &
    sll_world_collective

  use sll_m_hamiltonian_splitting_base, only: &
    sll_t_hamiltonian_splitting_base

  use sll_m_kernel_smoother_base, only: &
    sll_kernel_smoother_base

  use sll_m_maxwell_1d_base, only: &
    sll_maxwell_1d_base

  use sll_m_particle_group_base, only: &
    sll_particle_group_base

  use sll_mpi, only: &
    mpi_sum

  implicit none

  public :: &
    sll_new_hamiltonian_splitting_pic_vm_1d2v, &
    sll_t_hamiltonian_splitting_pic_vm_1d2v

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Operator splitting type for Vlasov-Maxwell 1d2v
  type, extends(sll_t_hamiltonian_splitting_base) :: sll_t_hamiltonian_splitting_pic_vm_1d2v
     class(sll_maxwell_1d_base), pointer  :: maxwell_solver      !< Maxwell solver
     class(sll_kernel_smoother_base), pointer :: kernel_smoother_0  !< Kernel smoother (order p+1)
     class(sll_kernel_smoother_base), pointer :: kernel_smoother_1  !< Kernel smoother (order p)
     class(sll_particle_group_base), pointer  :: particle_group    !< Particle group

     sll_int32 :: spline_degree !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx !< Size of the domain
     sll_real64 :: x_min !< Lower bound for x domain
     sll_real64 :: delta_x !< Grid spacing

     sll_real64 :: cell_integrals_0(4) !< Integral over the spline function on each interval (order p+1)
     sll_real64 :: cell_integrals_1(3) !< Integral over the spline function on each interval (order p)


     sll_real64, pointer     :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
     sll_real64, pointer     :: bfield_dofs(:) !< DoFs describing the magnetic field
     sll_real64, allocatable :: j_dofs(:,:)    !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:,:)!< MPI-processor local part of one component of \a j_dofs

   contains
     procedure :: operatorHp1 => operatorHp1_pic_vm_1d2v  !> Operator for H_p1 part
     procedure :: operatorHp2 => operatorHp2_pic_vm_1d2v  !> Operator for H_p2 part
     procedure :: operatorHE => operatorHE_pic_vm_1d2v  !> Operator for H_E part
     procedure :: operatorHB => operatorHB_pic_vm_1d2v  !> Operator for H_B part
     procedure :: lie_splitting => lie_splitting_pic_vm_1d2v !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_1d2v !> Strang splitting propagator
     procedure :: update_jv !> Helper function to compute the integral of j using Gauss quadrature

  end type sll_t_hamiltonian_splitting_pic_vm_1d2v

contains

  subroutine strang_splitting_pic_vm_1d2v(this,dt, number_steps)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step
    sll_int32, intent(in)  :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1, number_steps
       call this%operatorHp1(0.5_f64*dt)
       call this%operatorHp2(0.5_f64*dt)
       call this%operatorHE(0.5_f64*dt)
       call this%operatorHB(dt)
       call this%operatorHE(0.5_f64*dt)
       call this%operatorHp2(0.5_f64*dt)
       call this%operatorHp1(0.5_f64*dt)
    end do

  end subroutine strang_splitting_pic_vm_1d2v

  subroutine lie_splitting_pic_vm_1d2v(this,dt, number_steps)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step
    sll_int32, intent(in)  :: number_steps !< number of time steps

    sll_int32 :: i_step

    do i_step = 1, number_steps
       call this%operatorHp1(dt)
       call this%operatorHp2(dt)
       call this%operatorHE(dt)
       call this%operatorHB(dt)
    end do

  end subroutine lie_splitting_pic_vm_1d2v

  !---------------------------------------------------------------------------!
  !> Push Hp1: Equations to solve are
  !> \partial_t f + v_1 \partial_{x_1} f = 0    -> X_new = X_old + dt V_1
  !> V_new,2 = V_old,2 + \int_0 h V_old,1 B_old
  !> \partial_t E_1 = - \int v_1 f(t,x_1, v) dv -> E_{1,new} = E_{1,old} - \int \int v_1 f(t,x_1+s v_1,v) dv ds
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B = 0 => B_new = B_old 
  subroutine operatorHp1_pic_vm_1d2v(this, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part, ind
    sll_real64 :: xi, x_new(3), vi(3), wi(1)
    sll_int32  :: n_cells
    sll_real64 :: r_new, r_old
    sll_int32 :: index_new, index_old
    !sll_real64 :: primitive_1(3)

    this%j_dofs_local = 0.0_f64
    n_cells = this%kernel_smoother_0%n_dofs

    ! Here we have to accumulate j and integrate over the time interval.
    ! At each k=1,...,n_grid, we have for s \in [0,dt]:
    ! j_k(s) =  \sum_{i=1,..,N_p} q_i N((x_k+sv_{1,k}-x_i)/h) v_k,
    ! where h is the grid spacing and N the normalized B-spline
    ! In order to accumulate the integrated j, we normalize the values of x to the grid spacing, calling them y, we have
    ! j_k(s) = \sum_{i=1,..,N_p} q_i N(y_k+s/h v_{1,k}-y_i) v_k.
    ! Now, we want the integral 
    ! \int_{0..dt} j_k(s) d s = \sum_{i=1,..,N_p} q_i v_k \int_{0..dt} N(y_k+s/h v_{1,k}-y_i) ds =  \sum_{i=1,..,N_p} q_i v_k  \int_{0..dt}  N(y_k + w v_{1,k}-y_i) dw


    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
    ! Then update particle position:  X_new = X_old + dt * V
    do i_part=1,this%particle_group%n_particles  
       ! Read out particle position and velocity
       x_new = this%particle_group%get_x(i_part)
       vi = this%particle_group%get_v(i_part)
       ! Compute index_old, the index of the last DoF on the grid the particle contributes to, and r_old, its position (normalized to cell size one).
       xi = (x_new(1) - this%x_min) /&
            this%delta_x
       index_old = floor(xi)
       r_old = xi - real(index_old,f64)

       ! Then update particle position:  X_new = X_old + dt * V
       x_new = modulo(this%particle_group%get_x(i_part) + dt * vi, this%Lx)
       call this%particle_group%set_x(i_part, x_new)

       ! Compute the new box index index_new and normalized position r_old.
       xi = (x_new(1) - this%x_min) /&
            this%delta_x
       index_new = floor(xi)
       r_new = xi - real(index_new ,f64) 

       ! Scale vi by weight to combine both factors for accumulation of integral over j
       wi = this%particle_group%get_charge(i_part)

       if (index_old == index_new) then
          if (r_old < r_new) then
             call this%update_jv(r_old, r_new, index_old, wi, 1.0_f64, vi(2))
          else
             call this%update_jv(r_new, r_old, index_old, wi, -1.0_f64, vi(2))
          end if
       elseif (index_old < index_new) then
          call this%update_jv (r_old, 1.0_f64, index_old, wi, 1.0_f64, vi(2))
          call this%update_jv (0.0_f64, r_new, index_new, wi, 1.0_f64, vi(2))
          do ind = index_old+1, index_new-1
             call this%update_jv (0.0_f64, 1.0_f64, ind, wi, 1.0_f64, vi(2))
          end do
       else
          call this%update_jv (r_new, 1.0_f64, index_new, wi, -1.0_f64, vi(2))
          call this%update_jv (0.0_f64, r_old, index_old, wi, -1.0_f64, vi(2))
          do ind = index_new+1, index_old-1
             call this%update_jv (0.0_f64, 1.0_f64, ind, wi, -1.0_f64, vi(2))
          end do
       end if
             

       call this%particle_group%set_v(i_part, vi)

    end do

    this%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_collective_allreduce( sll_world_collective, this%j_dofs_local(:,1), &
         n_cells, MPI_SUM, this%j_dofs(:,1))

    ! Update the electric field. Also, we still need to scale with 1/Lx
    !this%j_dofs(:,1) = this%j_dofs(:,1)*this%delta_x!/this%Lx
    call this%maxwell_solver%compute_E_from_j(this%j_dofs(:,1), 1, this%efield_dofs(:,1))

    ! Finally, we recompute the shape factors for the new positions such that we can evaluate the electric and magnetic fields when calling the operators H_E and H_B
   call  this%kernel_smoother_0%compute_shape_factors(this%particle_group)
   call  this%kernel_smoother_1%compute_shape_factors(this%particle_group)

 end subroutine operatorHp1_pic_vm_1d2v

 ! TODO: This is hard coded for quadratic, cubic splines. Make general.
 subroutine update_jv(this, lower, upper, index, weight, sign, vi)
   class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: this !< time splitting object 
   sll_real64, intent(in) :: lower
   sll_real64, intent(in) :: upper
   sll_int32,  intent(in) :: index
   sll_real64, intent(in) :: weight(1)
   sll_real64, intent(in) :: sign
   sll_real64, intent(inout) :: vi

   !Local variables
   sll_real64 :: m, c, y1, y2, fy(this%spline_degree+1)
   sll_int32  :: ind, i_grid, i_mod, n_cells
   sll_real64 :: qm


   qm = this%particle_group%species%q_over_m();
   n_cells = this%kernel_smoother_0%n_dofs


   m = 0.5_f64*(upper-lower)
   c = 0.5_f64*(upper+lower)
   y1 = -m/sqrt(3.0_f64)+c
   y2 = m/sqrt(3.0_f64) +c
   fy = sign*m*this%delta_x*&
        (uniform_b_splines_at_x(this%spline_degree-1, y1) +&
        uniform_b_splines_at_x(this%spline_degree-1, y2))

   ind = 1
   do i_grid = index - this%spline_degree + 1, index
      i_mod = modulo(i_grid, n_cells ) + 1
      this%j_dofs_local(i_mod,1) = this%j_dofs_local(i_mod,1) + &
           weight(1)*fy(ind)
      vi = vi - qm* fy(ind)*this%bfield_dofs(i_mod)
      ind = ind + 1
   end do


 end subroutine update_jv

!PN DEFINED BUT NOT USED
!PN !> Alternative implementation of the Hp1 operator based on the primal function of the spline
!PN subroutine operatorHp1_pic_vm_1d2v_prim(this, dt)
!PN    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: this !< time splitting object 
!PN    sll_real64, intent(in) :: dt   !< time step
!PN
!PN    !local variables
!PN    sll_int32 :: i_part, i_grid, i_mod, ind, j
!PN    sll_real64 :: xi, x_new(3), vi(3), wi(1)
!PN    sll_int32  :: n_cells
!PN    sll_real64 :: r_new, r_old
!PN    sll_int32 :: index_new, index_old
!PN    sll_real64 :: primitive_1(3)
!PN    sll_real64 :: qm
!PN
!PN    this%j_dofs_local = 0.0_f64
!PN    n_cells = this%kernel_smoother_0%n_dofs    
!PN    qm = this%particle_group%species%q_over_m()
!PN
!PN    ! Here we have to accumulate j and integrate over the time interval.
!PN    ! At each k=1,...,n_grid, we have for s \in [0,dt]:
!PN    ! j_k(s) =  \sum_{i=1,..,N_p} q_i N((x_k+sv_{1,k}-x_i)/h) v_k,
!PN    ! where h is the grid spacing and N the normalized B-spline
!PN    ! In order to accumulate the integrated j, we normalize the values of x to the grid spacing, calling them y, we have
!PN    ! j_k(s) = \sum_{i=1,..,N_p} q_i N(y_k+s/h v_{1,k}-y_i) v_k.
!PN    ! Now, we want the integral 
!PN    ! \int_{0..dt} j_k(s) d s = \sum_{i=1,..,N_p} q_i v_k \int_{0..dt} N(y_k+s/h v_{1,k}-y_i) ds =  \sum_{i=1,..,N_p} q_i v_k  \int_{0..dt}  N(y_k + w v_{1,k}-y_i) dw
!PN
!PN
!PN    ! For each particle compute the index of the first DoF on the grid it contributes to and its position (normalized to cell size one). Note: j_dofs(_local) does not hold the values for j itself but for the integrated j.
!PN    ! Then update particle position:  X_new = X_old + dt * V
!PN    do i_part=1,this%particle_group%n_particles  
!PN       ! Read out particle position and velocity
!PN       x_new = this%particle_group%get_x(i_part)
!PN       vi = this%particle_group%get_v(i_part)
!PN       ! Compute index_old, the index of the last DoF on the grid the particle contributes to, and r_old, its position (normalized to cell size one).
!PN       xi = (x_new(1) - this%x_min) /&
!PN            this%delta_x
!PN       index_old = floor(xi)
!PN       r_old = xi - real(index_old,f64)
!PN
!PN       ! Then update particle position:  X_new = X_old + dt * V
!PN       x_new = modulo(this%particle_group%get_x(i_part) + dt * vi, this%Lx)
!PN       call this%particle_group%set_x(i_part, x_new)
!PN
!PN       ! Compute the new box index index_new and normalized position r_old.
!PN       xi = (x_new(1) - this%x_min) /&
!PN            this%delta_x
!PN       index_new = floor(xi)
!PN       r_new = xi - real(index_new ,f64) 
!PN
!PN       ! Scale vi by weight to combine both factors for accumulation of integral over j
!PN       wi = this%particle_group%get_charge(i_part)
!PN
!PN       ! Compute the primitives at r_old for each interval
!PN       primitive_1 = -primitive_uniform_quadratic_b_spline_at_x( r_old )
!PN       ! If we are integrating to the right, we need to also use the value of the primitive at the right interval bound. If larger, we would need the lower but it is zero by our normalization.   
!PN       if (index_old < index_new) then
!PN          primitive_1 = primitive_1 + this%cell_integrals_1
!PN       end if
!PN       ! Now, we loop through the DoFs and add the contributions.
!PN       ind = 1
!PN       do i_grid = index_old - this%spline_degree + 1, index_old
!PN          i_mod = modulo(i_grid, n_cells ) + 1
!PN          this%j_dofs_local(i_mod,1) = this%j_dofs_local(i_mod,1) + &
!PN               primitive_1(ind)*wi(1)
!PN          vi(2) = vi(2) - qm*primitive_1(ind)*this%bfield_dofs(i_mod)*this%delta_x
!PN          ind = ind + 1
!PN       end do      
!PN       ! Now contribution from r_new in the same way but different sign
!PN       primitive_1 = primitive_uniform_quadratic_b_spline_at_x( r_new )
!PN       if (index_old > index_new) then
!PN          primitive_1 = primitive_1 - this%cell_integrals_1
!PN       end if 
!PN
!PN       ind = 1
!PN       do i_grid = index_new - this%spline_degree + 1, index_new
!PN          i_mod = modulo(i_grid, n_cells ) + 1
!PN          this%j_dofs_local(i_mod,1) = this%j_dofs_local(i_mod,1) + &
!PN               primitive_1(ind)*wi(1)
!PN          vi(2) = vi(2) - qm*primitive_1(ind)*this%bfield_dofs(i_mod)*this%delta_x
!PN          ind = ind + 1
!PN       end do
!PN       ! If |index_new - index_old|<2, we are done. Otherwise, we have to account for the piecewise definition of the spline and add the contributions from the intervals inbetween. These are always integrals over the whole integrals.
!PN       ! First if index_old<index_new-1: Add the whole integrals.
!PN       do j = index_old+1, index_new-1
!PN          ind = 1
!PN          do i_grid = j - this%spline_degree + 1, j
!PN             i_mod = modulo(i_grid, n_cells ) + 1
!PN             this%j_dofs_local(i_mod,1) = this%j_dofs_local(i_mod,1) + &
!PN                  this%cell_integrals_1(ind)*wi(1)
!PN             vi(2) = vi(2) - &
!PN                  qm*this%cell_integrals_1(ind)*this%bfield_dofs(i_mod)*this%delta_x
!PN             ind = ind + 1
!PN          end do
!PN       end do
!PN       ! Then if index_old>index_new-1: Subtract the whole integrals.
!PN       do j = index_new+1, index_old-1
!PN          ind = 1
!PN          do i_grid = j - this%spline_degree + 1, j
!PN             i_mod = modulo(i_grid, n_cells ) + 1
!PN             this%j_dofs_local(i_mod,1) = this%j_dofs_local(i_mod,1) - &
!PN                  this%cell_integrals_1(ind)*wi(1)
!PN             vi(2) = vi(2) + qm*this%cell_integrals_1(ind)*this%bfield_dofs(i_mod)*this%delta_x
!PN             ind = ind + 1
!PN          end do
!PN       end do 
!PN
!PN       call this%particle_group%set_v(i_part, vi)
!PN
!PN    end do
!PN
!PN    this%j_dofs = 0.0_f64
!PN    ! MPI to sum up contributions from each processor
!PN    call sll_collective_allreduce( sll_world_collective, this%j_dofs_local(:,1), &
!PN         n_cells, MPI_SUM, this%j_dofs(:,1))
!PN
!PN    ! Update the electric field. Also, we still need to scale with 1/Lx
!PN    this%j_dofs(:,1) = this%j_dofs(:,1)*this%delta_x
!PN    call this%maxwell_solver%compute_E_from_j(this%j_dofs(:,1), 1, this%efield_dofs(:,1))
!PN
!PN    ! Finally, we recompute the shape factors for the new positions such that we can evaluate the electric and magnetic fields when calling the operators H_E and H_B
!PN   call  this%kernel_smoother_0%compute_shape_factors(this%particle_group)
!PN   call  this%kernel_smoother_1%compute_shape_factors(this%particle_group)
!PN
!PN end subroutine operatorHp1_pic_vm_1d2v_prim





 !---------------------------------------------------------------------------!
  !> Push Hp2: Equations to solve are
  !> X_new = X_old
  !> V_new,1 = V_old,1 + \int_0 h V_old,2 B_old
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = - \int v_2 f(t,x_1, v) dv -> E_{2,new} = E_{2,old} - \int \int v_2 f(t,x_1+s v_1,v) dv ds
  !> \partial_t B = 0 => B_new = B_old
  subroutine operatorHp2_pic_vm_1d2v(this, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32  :: i_part, n_cells
    sll_real64 :: vi(3), xi(3), wi(1)
    sll_real64 :: x_box, r_box
    sll_int32  :: index_box, ind, i_mod, i_grid
    sll_real64 :: values_0(4)
    sll_real64 :: bfield
    sll_real64 :: qm

    n_cells = this%kernel_smoother_0%n_dofs

    this%j_dofs_local = 0.0_f64

    qm = this%particle_group%species%q_over_m();
    ! Update v_1
    do i_part=1,this%particle_group%n_particles
       ! Evaluate bfield at particle position (splines of order p)
       call this%kernel_smoother_1%evaluate_kernel_function_particle&
            (this%bfield_dofs, i_part, bfield)
       vi = this%particle_group%get_v(i_part)
       vi(1) = vi(1) + dt*qm*vi(2)*bfield
       call this%particle_group%set_v(i_part, vi)

       xi = this%particle_group%get_x(i_part)
       ! Compute the new box index index_box and normalized position r_box.
       x_box = (xi(1) - this%x_min) /&
            this%delta_x
       index_box = floor(x_box)
       r_box = x_box - real(index_box ,f64) 

       ! Scale vi by weight to combine both factors for accumulation of integral over j
       wi = this%particle_group%get_charge(i_part)
       vi = vi*wi(1)

       ! Now, we loop through the DoFs and add the contributions.
       values_0 = uniform_b_splines_at_x(this%spline_degree, r_box)
       ind = 1
       do i_grid = index_box - this%spline_degree, index_box
          i_mod = modulo(i_grid, n_cells ) + 1
          this%j_dofs_local(i_mod,2) = this%j_dofs_local(i_mod,2) + &
               values_0(ind)*vi (2)
          ind = ind + 1
       end do     

    end do

    this%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_collective_allreduce( sll_world_collective, this%j_dofs_local(:,2), &
         n_cells, MPI_SUM, this%j_dofs(:,2))

    ! Update the electric field. Also, we still need to scale with 1/Lx ! TODO: Which scaling?
    this%j_dofs(:,2) = this%j_dofs(:,2)*dt!/this%Lx
    call this%maxwell_solver%compute_E_from_j(this%j_dofs(:,2), 2, this%efield_dofs(:,2))
    !!this%efield_dofs(:,2) = this%efield_dofs(:,2) - this%j_dofs(:,2)/this%Lx*dt
    

  end subroutine operatorHp2_pic_vm_1d2v
  
  !---------------------------------------------------------------------------!
  !> Push H_E: Equations to be solved
  !> \partial_t f + E_1 \partial_{v_1} f + E_2 \partial_{v_2} f = 0 -> V_new = V_old + dt * E
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B + \partial_{x_1} E_2 = 0 => B_new = B_old - dt \partial_{x_1} E_2
  subroutine operatorHE_pic_vm_1d2v(this, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: v_new(3)
    sll_real64 :: efield(2)
    sll_real64 :: qm


    qm = this%particle_group%species%q_over_m();
    ! V_new = V_old + dt * E
    do i_part=1,this%particle_group%n_particles
       v_new = this%particle_group%get_v(i_part)
       ! Evaluate efields at particle position
       call this%kernel_smoother_1%evaluate_kernel_function_particle&
            (this%efield_dofs(:,1), i_part, efield(1))
       call this%kernel_smoother_0%evaluate_kernel_function_particle&
            (this%efield_dofs(:,2), i_part, efield(2))
       v_new = this%particle_group%get_v(i_part)
       v_new(1:2) = v_new(1:2) + dt* qm * efield
       call this%particle_group%set_v(i_part, v_new)
    end do
    
    ! Update bfield
    call this%maxwell_solver%compute_B_from_E( &
         dt, this%efield_dofs(:,2), this%bfield_dofs)
    

  end subroutine operatorHE_pic_vm_1d2v
  

  !---------------------------------------------------------------------------!
  !> Push H_B: Equations to be solved
  !> V_new = V_old
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old}
  !> \partial_t E_2 = - \partial_{x_1} B -> E_{2,new} = E_{2,old}-dt*\partial_{x_1} B
  !> \partial_t B = 0 -> B_new = B_old
  subroutine operatorHB_pic_vm_1d2v(this, dt)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    !sll_int32 :: i_part
    !sll_real64 :: v_old(3), v_new(3)

    
    ! Update efield2
    call this%maxwell_solver%compute_E_from_B(&
         dt, this%bfield_dofs, this%efield_dofs(:,2))
    

  end subroutine operatorHB_pic_vm_1d2v



  !---------------------------------------------------------------------------!
  !> Constructor.
  function sll_new_hamiltonian_splitting_pic_vm_1d2v(&
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx) result(this)
    class(sll_t_hamiltonian_splitting_pic_vm_1d2v), pointer :: this !< time splitting object 
    class(sll_maxwell_1d_base), pointer, intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_kernel_smoother_base), pointer, intent(in) :: kernel_smoother_0  !< Kernel smoother
    class(sll_kernel_smoother_base), pointer, intent(in) :: kernel_smoother_1  !< Kernel smoother
    class(sll_particle_group_base),pointer, intent(in) :: particle_group !< Particle group
    sll_real64, pointer, intent(in) :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, pointer, intent(in) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64, intent(in) :: x_min !< Lower bound of x domain
    sll_real64, intent(in) :: Lx !< Length of the domain in x direction.

    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE(this, ierr)

    this%maxwell_solver => maxwell_solver
    this%kernel_smoother_0 => kernel_smoother_0
    this%kernel_smoother_1 => kernel_smoother_1
    this%particle_group => particle_group
    this%efield_dofs => efield_dofs
    this%bfield_dofs => bfield_dofs

    ! TODO: Check that n_dofs is the same for both kernel smoothers.

    SLL_ALLOCATE(this%j_dofs(this%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(this%j_dofs_local(this%kernel_smoother_0%n_dofs,2), ierr)

    this%spline_degree = 3
    this%x_min = x_min
    this%Lx = Lx
    this%delta_x = this%Lx/this%kernel_smoother_1%n_dofs
    
    this%cell_integrals_1 = [0.5_f64, 2.0_f64, 0.5_f64]
    this%cell_integrals_1 = this%cell_integrals_1 / 3.0_f64

    this%cell_integrals_0 = [1.0_f64,11.0_f64,11.0_f64,1.0_f64]
    this%cell_integrals_0 = this%cell_integrals_0 / 24.0_f64

  end function sll_new_hamiltonian_splitting_pic_vm_1d2v

!PN DEFINED BUT NOT USED
!PN   !> Compute the primitive of the cubic B-spline in each intervall at x. Primitive function normalized such that it is 0 at x=0. Analogon to uniform_b_spline_at_x in arbitrary degree splines for primitive, but specific for cubic.
!PN   function primitive_uniform_cubic_b_spline_at_x( x) result(primitive)
!PN     sll_real64, intent(in)  :: x !< position where to evaluate the primitive
!PN     sll_real64 :: primitive(4) !< value of the primitive for each of the four intervals.
!PN 
!PN     sll_real64 :: xx(3)
!PN 
!PN     xx(1) = x**2
!PN     xx(2) = x*xx(1)
!PN     xx(3) = x*xx(2)
!PN 
!PN     primitive(4) = xx(3)/24.0_f64
!PN     primitive(3) = (x + 1.5_f64*xx(1) + xx(2) - 0.75_f64* xx(3))/6.0_f64
!PN     primitive(2) = (4.0_f64*x - 2.0_f64* xx(2) + 0.75_f64* xx(3))/6.0_f64
!PN     primitive(1) = (1.0_f64 - (1.0_f64-x)**4)/24.0_f64
!PN 
!PN   end function primitive_uniform_cubic_b_spline_at_x
!PN 
!PN   !> Compute the primitive of the quadratic B-spline in each intervall at x. Primitive function normalized such that it is 0 at x=0. Analogon to uniform_b_spline_at_x in arbitrary degree splines for primitive, but specific for quadratic.
!PN   function primitive_uniform_quadratic_b_spline_at_x( x) result(primitive)
!PN     sll_real64, intent(in)  :: x !< position where to evaluate the primitive
!PN     sll_real64 :: primitive(3) !< value of the primitive for each of the three intervals.
!PN 
!PN     sll_real64 :: xx(2)
!PN 
!PN     xx(1) = x**2
!PN     xx(2) = x*xx(1)
!PN 
!PN     primitive(3) = xx(2)/6.0_f64
!PN     primitive(2) = (x + xx(1))*0.5_f64-xx(2)/3.0_f64
!PN     primitive(1) = (1.0_f64 - (1.0_f64-x)**3)/6.0_f64
!PN 
!PN   end function primitive_uniform_quadratic_b_spline_at_x
!PN 
!PN 
!PN 

end module sll_m_hamiltonian_splitting_pic_vm_1d2v
