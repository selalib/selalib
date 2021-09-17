!> @ingroup pic_time_integration
!> @author Benedikt Perse, IPP
!> @brief Particle pusher based on Hamiltonian splitting for 1d2v Vlasov-Maxwell with coordinate transformation.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
!> Reference: Kraus, Kormann, Sonnendrücker, Morrison: GEMPIC: Geometric ElectroMagnetic Particle-In-Cell Methods
module sll_m_time_propagator_pic_vm_1d2v_hs_trafo
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

  use sll_mpi, only: &
       mpi_sum

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_1d2v_hs_trafo

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_1d2v_hs_trafo
     class(sll_c_maxwell_1d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_0  !< Kernel smoother (order p+1)
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_1  !< Kernel smoother (order p)
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

     type( sll_t_mapping_3d ), pointer      :: map !< Coordinate transformation

     sll_int32 :: spline_degree !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx !< Size of the domain
     sll_real64 :: x_min !< Lower bound for x domain
     sll_real64 :: delta_x !< Grid spacing

     sll_real64, pointer     :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
     sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the magnetic field
     sll_real64, allocatable :: j_dofs(:,:)      !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:,:)!< MPI-processor local part of one component of \a j_dofs

     logical :: electrostatic = .false.  !< true for electrostatic simulation

     sll_int32 :: max_iter = 1 !< maximal amount of iterations, set 1 for ctest
     sll_real64 :: iter_tolerance = 1d-10 !< iteration tolerance

   contains
     procedure :: operatorHp => operatorHp_pic_vm_1d2v  !> Operator for H_p part
     procedure :: operatorHE => operatorHE_pic_vm_1d2v  !> Operator for H_E part
     procedure :: operatorHB => operatorHB_pic_vm_1d2v  !> Operator for H_B part
     procedure :: lie_splitting => lie_splitting_pic_vm_1d2v !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_1d2v !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_1d2v !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_1d2v !> Initialize the type
     procedure :: free => delete_pic_vm_1d2v !> Finalization

  end type sll_t_time_propagator_pic_vm_1d2v_hs_trafo

contains


  !> Strang splitting
  subroutine strang_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_hs_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%operatorHE(0.5_f64*dt)
          call self%operatorHp(dt)
          call self%operatorHE(0.5_f64*dt)
       end do
    else
       do i_step = 1, number_steps
          call self%operatorHB(0.5_f64*dt)
          call self%operatorHE(0.5_f64*dt)
          call self%operatorHp(dt)
          call self%operatorHE(0.5_f64*dt)
          call self%operatorHB(0.5_f64*dt)
       end do
    end if

  end subroutine strang_splitting_pic_vm_1d2v

  
  !> Lie splitting
  subroutine lie_splitting_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_hs_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%operatorHE(dt)
          call self%operatorHp(dt)
       end do
    else
       do i_step = 1,number_steps
          call self%operatorHE(dt)
          call self%operatorHB(dt)
          call self%operatorHp(dt)
       end do
    end if


  end subroutine lie_splitting_pic_vm_1d2v

  
  !> Lie splitting (oposite ordering)
  subroutine lie_splitting_back_pic_vm_1d2v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_1d2v_hs_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps

    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%operatorHp(dt)
          call self%operatorHE(dt)
       end do
    else
       do i_step = 1,number_steps
          call self%operatorHp(dt)
          call self%operatorHB(dt)
          call self%operatorHE(dt)
       end do
    end if

  end subroutine lie_splitting_back_pic_vm_1d2v


  !---------------------------------------------------------------------------!
  !> Push Hp1: Equations to solve are
  !> \partial_t f + v_1 \partial_{x_1} f = 0    -> X_new = X_old + dt V_1
  !> V_new,2 = V_old,2 + \int_0 h V_old,1 B_old
  !> \partial_t E_1 = - \int v_1 f(t,x_1, v) dv -> E_{1,new} = E_{1,old} - \int \int v_1 f(t,x_1+s v_1,v) dv ds
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B = 0 => B_new = B_old 
  subroutine operatorHp_pic_vm_1d2v(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_hs_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp, i
    sll_real64 :: x_new(3), v_new(3), vi(3), wi(1), xi(3), xbar(3), vbar(3)
    sll_int32  :: n_cells
    sll_real64 :: jmatrix(3,3)
    sll_real64 :: qmdt, err
    sll_real64 :: bfield, rhs(2)


    n_cells = self%kernel_smoother_0%n_dofs


    self%j_dofs_local = 0.0_f64

    do i_sp = 1, self%particle_group%n_species
       qmdt = self%particle_group%group(i_sp)%species%q_over_m()*dt;
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          !Predictor-Corrector with loop for corrector step
          jmatrix=self%map%jacobian_matrix_inverse([xi(1), 0._f64, 0._f64])
          !Transformation of the v coordinates, !VT = DF^{-1} * vi 
          !x^\star=\mod(x^n+dt*DF^{-1}(x^n)v^n,1)
          x_new(1) = xi(1) + dt * jmatrix(1,1)*vi(1)

          call self%kernel_smoother_1%evaluate &
               (xi(1), self%bfield_dofs, bfield)

          !rhs = vi + sign * DF^{-T} * bfield x DF^{-1} V
          v_new(1) = vi(1) + qmdt*jmatrix(1,1)*bfield*vi(2)
          v_new(2) = vi(2) - qmdt*bfield*jmatrix(1,1)*vi(1)


          err= max( abs(xi(1) - x_new(1)), maxval(abs(vi(1:2) - v_new(1:2))) ) 
          i = 0
          do while(i < self%max_iter .and. err > self%iter_tolerance)
             xbar(1) = modulo(0.5_f64*(x_new(1)+xi(1)), 1._f64)
             vbar = 0.5_f64*(v_new+vi)

             jmatrix = self%map%jacobian_matrix_inverse( [xbar(1),0._f64,0._f64] )

             call self%kernel_smoother_1%evaluate &
                  (xbar(1), self%bfield_dofs, bfield)

             xbar(1) = xi(1) + dt * jmatrix(1,1)*vbar(1)
             !rhs = vi + sign * DF^{-T} * bfield x DF^{-1} V
             rhs(1) = vi(1) + qmdt*jmatrix(1,1)*bfield*vbar(2)
             rhs(2) = vi(2) - qmdt*bfield*jmatrix(1,1)*vbar(1)

             err = max(abs(x_new(1) - xbar(1)), maxval(abs(v_new(1:2) - rhs(1:2))) )
             x_new = xbar
             v_new(1:2) = rhs
             i = i + 1
          end do
          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part)

          if ( abs(vbar(1))> 1d-16 ) then
             call self%kernel_smoother_1%add_current( xi(1), x_new(1), wi(1), self%j_dofs_local(:,1) )
             call self%kernel_smoother_0%add_current( xi(1), x_new(1), wi(1)*vbar(2)/(jmatrix(1,1)*vbar(1)), self%j_dofs_local(:,2) )
          else
             call self%kernel_smoother_0%add_charge( xi(1), wi(1)*dt*vbar(2), self%j_dofs_local(:,2) )
          end if

          x_new(1) = modulo(x_new(1), 1._f64)
          call self%particle_group%group(i_sp)%set_x(i_part, x_new)
          call self%particle_group%group(i_sp)%set_v( i_part, v_new )

       end do
    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
         n_cells, MPI_SUM, self%j_dofs(:,1))

    ! Update the electric field.
    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,1), 1, self%efield_dofs(:,1))
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
         n_cells, MPI_SUM, self%j_dofs(:,2))

    call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,2), 2, self%efield_dofs(:,2))

  end subroutine operatorHp_pic_vm_1d2v


  !---------------------------------------------------------------------------!
  !> Push H_E: Equations to be solved
  !> \partial_t f + E_1 \partial_{v_1} f + E_2 \partial_{v_2} f = 0 -> V_new = V_old + dt * E
  !> \partial_t E_1 = 0 -> E_{1,new} = E_{1,old} 
  !> \partial_t E_2 = 0 -> E_{2,new} = E_{2,old}
  !> \partial_t B + \partial_{x_1} E_2 = 0 => B_new = B_old - dt \partial_{x_1} E_2
  subroutine operatorHE_pic_vm_1d2v(self, dt)
    class(sll_t_time_propagator_pic_vm_1d2v_hs_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: vi(3), xi(3)
    sll_real64 :: efield(2)
    sll_real64 :: qoverm, jmat(3,3)

    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       ! V_new = V_old + dt * E
       do i_part=1,self%particle_group%group(i_sp)%n_particles
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          ! Evaluate efields at particle position
          call self%kernel_smoother_1%evaluate &
               (xi(1), self%efield_dofs(:,1), efield(1))
          call self%kernel_smoother_0%evaluate &
               (xi(1), self%efield_dofs(:,2), efield(2))

          jmat = self%map%jacobian_matrix_inverse_transposed( [xi(1), 0._f64, 0._f64] )

          vi(1) = vi(1) + dt* qoverm * jmat(1,1)* efield(1)
          vi(2) = vi(2) + dt* qoverm * jmat(2,2)*efield(2)
          call self%particle_group%group(i_sp)%set_v(i_part, vi)
       end do
    end do

    ! Update bfield
    if(self%electrostatic .eqv. .false.) then
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
    class(sll_t_time_propagator_pic_vm_1d2v_hs_trafo), intent(inout) :: self !< time propagator object 
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
       Lx, &
       map, &
       electrostatic) 
    class(sll_t_time_propagator_pic_vm_1d2v_hs_trafo), intent(out) :: self !< time propagator object 
    class(sll_c_maxwell_1d_base), target,          intent(in)  :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,      intent(in)  :: particle_group !< Particle group
    sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
    type(sll_t_mapping_3d), target,                intent( inout ) :: map !< Coordinate transformation
    logical, optional    :: electrostatic !< true for electrostatic simulation
    !local variables
    sll_int32 :: ierr

    if( present(electrostatic) )then
       self%electrostatic = electrostatic
    end if

    self%maxwell_solver => maxwell_solver
    self%kernel_smoother_0 => kernel_smoother_0
    self%kernel_smoother_1 => kernel_smoother_1
    self%particle_group => particle_group
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs
    self%map => map

    ! Check that n_dofs is the same for both kernel smoothers.
    SLL_ASSERT( self%kernel_smoother_0%n_dofs == self%kernel_smoother_1%n_dofs )

    SLL_ALLOCATE(self%j_dofs(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%j_dofs_local(self%kernel_smoother_0%n_dofs,2), ierr)

    self%spline_degree = self%kernel_smoother_0%spline_degree
    self%x_min = x_min
    self%Lx = Lx
    self%delta_x = 1._f64/real(self%kernel_smoother_1%n_dofs,f64)
    
  end subroutine initialize_pic_vm_1d2v


  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_1d2v(self)
    class(sll_t_time_propagator_pic_vm_1d2v_hs_trafo), intent( inout ) :: self !< time propagator object
    
    deallocate(self%j_dofs)
    deallocate(self%j_dofs_local)
    self%maxwell_solver => null()
    self%kernel_smoother_0 => null()
    self%kernel_smoother_1 => null()
    self%particle_group => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()
    self%map => null()

  end subroutine delete_pic_vm_1d2v


end module sll_m_time_propagator_pic_vm_1d2v_hs_trafo
