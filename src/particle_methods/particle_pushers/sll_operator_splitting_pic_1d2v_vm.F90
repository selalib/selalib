!> @ingroup particle_methods
!> @author Katharina Kormann
!> @brief Particle pusher based on operator splitting for 1d2v Vlasov-Poisson.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_operator_splitting_pic_1d2v_vm
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_errors.h"


 
  use sll_module_pic_base
  use sll_m_kernel_smoother_base
  use sll_collective
  use sll_arbitrary_degree_splines
  
  use sll_m_maxwell_1d_base

  implicit none

  !> Operator splitting type for Vlasov-Maxwell 1d2v
  type :: sll_operator_splitting_pic_1d2v_vm
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


     sll_real64, allocatable :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
     sll_real64, allocatable :: bfield_dofs(:) !< DoFs describing the magnetic field
     sll_real64, allocatable :: efield(:,:)    !< Efield at particle positions
     sll_real64, allocatable :: bfield(:)      !< Bfield at particle positions
     sll_real64, allocatable :: j_dofs(:,:)    !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:,:)!< MPI-processor local part of one component of \a j_dofs

   contains
     procedure :: operatorHp1 => operatorHp1_pic_1d2v_vm_prim  !< Operator for H_p1 part
     procedure :: operatorHp2 => operatorHp2_pic_1d2v_vm  !< Operator for H_p2 part
     procedure :: operatorHE => operatorHE_pic_1d2v_vm  !< Operator for H_E part
     procedure :: operatorHB => operatorHB_pic_1d2v_vm  !< Operator for H_B part
     procedure :: strang_splitting => strang_splitting_pic_1d2v_vm
     procedure :: update_jv

  end type sll_operator_splitting_pic_1d2v_vm

contains

  subroutine strang_splitting_pic_1d2v_vm(this,dt)
    class(sll_operator_splitting_pic_1d2v_vm), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    call this%operatorHp1(0.5_f64*dt)
    call this%operatorHp2(0.5_f64*dt)
    call this%operatorHE(0.5_f64*dt)
    call this%operatorHB(dt)
    call this%operatorHE(0.5_f64*dt)
    call this%operatorHp2(0.5_f64*dt)
    call this%operatorHp1(0.5_f64*dt)

  end subroutine strang_splitting_pic_1d2v_vm

  !---------------------------------------------------------------------------!
  !> Push Hf1: Equations to solve are
  !> \partial_t f + v_1 \partial_{x_1} f = 0    -> X_new = X_old + dt V_1
  !> \partial_t E_1 = - \int v_1 f(t,x_1, v) dv -> E_{1,new} = E_{1,old} - \int \int v_1 f(t,x_1+s v_1,v) dv ds
  !> \partial_t E_2 = - \int v_2 f(t,x_1, v) dv -> E_{2,new} = E_{2,old} - \int \int v_2 f(t,x_1+s v_1,v) dv ds
  !> \partial_t B = 0 => B_new = B_old 
  subroutine operatorHp1_pic_1d2v_vm(this, dt)
    class(sll_operator_splitting_pic_1d2v_vm), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part, i_grid, i_mod, ind, j
    sll_real64 :: xi, x_new(3), vi(3), wi
    sll_int32  :: n_cells
    sll_real64 :: r_new, r_old
    sll_int32 :: index_new, index_old
    sll_real64 :: primitive_1(3)

    this%j_dofs_local = 0.0_f64
    n_cells = this%kernel_smoother_0%n_dofs

    ! Here we have to accumulate j and integrate over the time interval.
    ! At each k=1,...,n_grid, we have for s \in [0,dt]:
    ! j_k(s) = 1/L_x \sum_{i=1,..,N_p} q_i N((x_k+sv_{1,k}-x_i)/h)/h v_k,
    ! where h is the grid spacing and N the normalized B-spline
    ! In order to accumulate the integrated j, we normalize the values of x to the grid spacing, calling them y, we have
    ! j_k(s) = 1/L_x  \sum_{i=1,..,N_p} q_i N(y_k+s/h v_{1,k}-y_i)/h v_k.
    ! Now, we want the integral 
    ! \int_{0..dt} j_k(s) d s = 1/Lx \sum_{i=1,..,N_p} q_i v_k \int_{0..dt} N(y_k+s/h v_{1,k}-y_i)/h ds =  1/Lx \sum_{i=1,..,N_p} q_i v_k  \int_{0..dt/h}  N(y_k + w v_{1,k}-y_i) dw


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

 end subroutine operatorHp1_pic_1d2v_vm

 ! TODO: This is hard coded for quadratic, cubic splines. Make general.
! TODO: VZ unklar.
 subroutine update_jv(this, lower, upper, index, weight, sign, vi)
   class(sll_operator_splitting_pic_1d2v_vm), intent(inout) :: this !< time splitting object 
   sll_real64, intent(in) :: lower
   sll_real64, intent(in) :: upper
   sll_int32,  intent(in) :: index
   sll_real64, intent(in) :: weight
   sll_real64, intent(in) :: sign
   sll_real64, intent(inout) :: vi

   !Local variables
   sll_real64 :: m, c, y1, y2, fy(this%spline_degree+1)
   sll_int32  :: ind, i_grid, i_mod, n_cells


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
           weight*fy(ind)
      vi = vi - fy(ind)*this%bfield_dofs(i_mod)
      ind = ind + 1
   end do


 end subroutine update_jv

 subroutine operatorHp1_pic_1d2v_vm_prim(this, dt)
    class(sll_operator_splitting_pic_1d2v_vm), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part, i_grid, i_mod, ind, j
    sll_real64 :: xi, x_new(3), vi(3), wi
    sll_int32  :: n_cells
    sll_real64 :: r_new, r_old
    sll_int32 :: index_new, index_old
    sll_real64 :: primitive_1(3)

    this%j_dofs_local = 0.0_f64
    n_cells = this%kernel_smoother_0%n_dofs

    ! Here we have to accumulate j and integrate over the time interval.
    ! At each k=1,...,n_grid, we have for s \in [0,dt]:
    ! j_k(s) = 1/L_x \sum_{i=1,..,N_p} q_i N((x_k+sv_{1,k}-x_i)/h)/h v_k,
    ! where h is the grid spacing and N the normalized B-spline
    ! In order to accumulate the integrated j, we normalize the values of x to the grid spacing, calling them y, we have
    ! j_k(s) = 1/L_x  \sum_{i=1,..,N_p} q_i N(y_k+s/h v_{1,k}-y_i)/h v_k.
    ! Now, we want the integral 
    ! \int_{0..dt} j_k(s) d s = 1/Lx \sum_{i=1,..,N_p} q_i v_k \int_{0..dt} N(y_k+s/h v_{1,k}-y_i)/h ds =  1/Lx \sum_{i=1,..,N_p} q_i v_k  \int_{0..dt/h}  N(y_k + w v_{1,k}-y_i) dw


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

       ! Compute the primitives at r_old for each interval
       primitive_1 = -primitive_uniform_quadratic_b_spline_at_x( r_old )
       ! If we are integrating to the right, we need to also use the value of the primitive at the right interval bound. If larger, we would need the lower but it is zero by our normalization.   
       if (index_old < index_new) then
          primitive_1 = primitive_1 + this%cell_integrals_1
       end if
       ! Now, we loop through the DoFs and add the contributions.
       ind = 1
       do i_grid = index_old - this%spline_degree + 1, index_old
          i_mod = modulo(i_grid, n_cells ) + 1
          this%j_dofs_local(i_mod,1) = this%j_dofs_local(i_mod,1) + &
               primitive_1(ind)*wi
          vi(2) = vi(2) - primitive_1(ind)*this%bfield_dofs(i_mod)*this%delta_x
          ind = ind + 1
       end do      
       ! Now contribution from r_new in the same way but different sign
       primitive_1 = primitive_uniform_quadratic_b_spline_at_x( r_new )
       if (index_old > index_new) then
          primitive_1 = primitive_1 - this%cell_integrals_1
       end if 

       ind = 1
       do i_grid = index_new - this%spline_degree + 1, index_new
          i_mod = modulo(i_grid, n_cells ) + 1
          this%j_dofs_local(i_mod,1) = this%j_dofs_local(i_mod,1) + &
               primitive_1(ind)*wi
          vi(2) = vi(2) - primitive_1(ind)*this%bfield_dofs(i_mod)*this%delta_x
          ind = ind + 1
       end do
       ! If |index_new - index_old|<2, we are done. Otherwise, we have to account for the piecewise definition of the spline and add the contributions from the intervals inbetween. These are always integrals over the whole integrals.
       ! First if index_old<index_new-1: Add the whole integrals.
       do j = index_old+1, index_new-1
          ind = 1
          do i_grid = j - this%spline_degree + 1, j
             i_mod = modulo(i_grid, n_cells ) + 1
             this%j_dofs_local(i_mod,1) = this%j_dofs_local(i_mod,1) + &
                  this%cell_integrals_1(ind)*wi
             vi(2) = vi(2) - &
                  this%cell_integrals_1(ind)*this%bfield_dofs(i_mod)*this%delta_x
             ind = ind + 1
          end do
       end do
       ! Then if index_old>index_new-1: Subtract the whole integrals.
       do j = index_new+1, index_old-1
          ind = 1
          do i_grid = j - this%spline_degree + 1, j
             i_mod = modulo(i_grid, n_cells ) + 1
             this%j_dofs_local(i_mod,1) = this%j_dofs_local(i_mod,1) - &
                  this%cell_integrals_1(ind)*wi
             vi(2) = vi(2) + this%cell_integrals_1(ind)*this%bfield_dofs(i_mod)*this%delta_x
             ind = ind + 1
          end do
       end do 

       call this%particle_group%set_v(i_part, vi)

    end do

    this%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_collective_allreduce( sll_world_collective, this%j_dofs_local(:,1), &
         n_cells, MPI_SUM, this%j_dofs(:,1))

    ! Update the electric field. Also, we still need to scale with 1/Lx
    this%j_dofs(:,1) = this%j_dofs(:,1)*this%delta_x!/this%Lx
    call this%maxwell_solver%compute_E_from_j(this%j_dofs(:,1), 1, this%efield_dofs(:,1))

    ! Finally, we recompute the shape factors for the new positions such that we can evaluate the electric and magnetic fields when calling the operators H_E and H_B
   call  this%kernel_smoother_0%compute_shape_factors(this%particle_group)
   call  this%kernel_smoother_1%compute_shape_factors(this%particle_group)

 end subroutine operatorHp1_pic_1d2v_vm_prim





 !---------------------------------------------------------------------------!
  !> Push Hf2: Equations to solve are
  !> \partial_t f + v_1 \partial_{x_1} f = 0    -> X_new = X_old + dt V_1
  !> \partial_t E_1 = - \int v_1 f(t,x_1, v) dv -> E_{1,new} = E_{1,old} - \int \int v_1 f(t,x_1+s v_1,v) dv ds
  !> \partial_t E_2 = - \int v_2 f(t,x_1, v) dv -> E_{2,new} = E_{2,old} - \int \int v_2 f(t,x_1+s v_1,v) dv ds
  !> \partial_t B = 0 => B_new = B_old
  subroutine operatorHp2_pic_1d2v_vm(this, dt)
    class(sll_operator_splitting_pic_1d2v_vm), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32  :: i_part, n_cells
    sll_real64 :: vi(3), xi(3)
    sll_real64 :: x_box, r_box
    sll_int32  :: index_box, ind, i_mod, i_grid
    sll_real64 :: values_0(4)

    ! Evaluate bfield at particle positions (splines of order p)
    call this%kernel_smoother_1%evaluate_kernel_function(this%particle_group, &
         this%bfield_dofs, this%bfield)

    n_cells = this%kernel_smoother_0%n_dofs

    this%j_dofs_local = 0.0_f64

    ! Update v_1
    do i_part=1,this%particle_group%n_particles
       vi = this%particle_group%get_v(i_part)
       vi(1) = vi(1) + dt*vi(2)*this%bfield(i_part)
       call this%particle_group%set_v(i_part, vi)

       xi = this%particle_group%get_x(i_part)
       ! Compute the new box index index_box and normalized position r_box.
       x_box = (xi(1) - this%x_min) /&
            this%delta_x
       index_box = floor(x_box)
       r_box = x_box - real(index_box ,f64) 

       ! Scale vi by weight to combine both factors for accumulation of integral over j
       vi = vi*this%particle_group%get_charge(i_part)

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
    

 end subroutine operatorHp2_pic_1d2v_vm
  
  !---------------------------------------------------------------------------!
  !> Push H_E: Equations to be solved
  !> \partial_t f + E_1 \partial_{v_1} f + E_2 \partial_{v_2} f = 0 -> V_new = V_old + dt * E
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B + \partial_{x_1} E_2 = 0 => B_new = B_old - dt \partial_{x_1} E_2
  subroutine operatorHE_pic_1d2v_vm(this, dt)
    class(sll_operator_splitting_pic_1d2v_vm), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: v_new(3)

    ! Evaluate efield at particle positions (quadratic for E_x and cubic for E_y)
    call this%kernel_smoother_1%evaluate_kernel_function(this%particle_group, &
         this%efield_dofs(:,1), this%efield(:,1))
    call this%kernel_smoother_0%evaluate_kernel_function(this%particle_group, &
         this%efield_dofs(:,2), this%efield(:,2))

    ! V_new = V_old + dt * E
    do i_part=1,this%particle_group%n_particles
       v_new = this%particle_group%get_v(i_part)
       v_new(1:2) = v_new(1:2) + dt * this%efield(i_part,:) 
       call this%particle_group%set_v(i_part, v_new)
    end do
    
    ! Update bfield
    call this%maxwell_solver%compute_B_from_E( &
         dt, this%efield_dofs(:,2), this%bfield_dofs)
    

  end subroutine operatorHE_pic_1d2v_vm
  

  !---------------------------------------------------------------------------!
  !> Push H_B: Equations to be solved
  !> V_new = V_old
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old}
  !> \partial_t E_2 = - \partial_{x_1} B -> E_{2,new} = E_{2,old}-dt*\partial_{x_1} B
  !> \partial_t B = 0 -> B_new = B_old
  subroutine operatorHB_pic_1d2v_vm(this, dt)
    class(sll_operator_splitting_pic_1d2v_vm), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: v_old(3), v_new(3)

    
    ! Update efield2
    call this%maxwell_solver%compute_E_from_B(&
         dt, this%bfield_dofs, this%efield_dofs(:,2))
    

  end subroutine operatorHB_pic_1d2v_vm



  !---------------------------------------------------------------------------!
  !> Constructor.
  function sll_new_splitting_pic_1d2v_vm(&
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       x_min, &
       Lx) result(this)
    class(sll_operator_splitting_pic_1d2v_vm), pointer :: this !< time splitting object 
    class(sll_maxwell_1d_base), pointer, intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_kernel_smoother_base), pointer, intent(in) :: kernel_smoother_0  !< Kernel smoother
    class(sll_kernel_smoother_base), pointer, intent(in) :: kernel_smoother_1  !< Kernel smoother
    class(sll_particle_group_base),pointer, intent(in) :: particle_group !< Particle group
    sll_real64, intent(in) :: x_min !< Lower bound of x domain
    sll_real64, intent(in) :: Lx !< Length of the domain in x direction.

    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE(this, ierr)

    this%maxwell_solver => maxwell_solver
    this%kernel_smoother_0 => kernel_smoother_0
    this%kernel_smoother_1 => kernel_smoother_1
    this%particle_group => particle_group

    ! TODO: Check that n_dofs is the same for both kernel smoothers.

    SLL_ALLOCATE(this%efield_dofs(this%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(this%bfield_dofs(this%kernel_smoother_0%n_dofs), ierr)
    SLL_ALLOCATE(this%efield(this%particle_group%n_particles,2), ierr)
    SLL_ALLOCATE(this%bfield(this%particle_group%n_particles), ierr)
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

  end function sll_new_splitting_pic_1d2v_vm


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




end module sll_m_operator_splitting_pic_1d2v_vm
