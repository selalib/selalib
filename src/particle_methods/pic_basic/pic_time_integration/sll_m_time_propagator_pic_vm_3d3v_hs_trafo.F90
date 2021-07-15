!> @ingroup pic_time_integration
!> @author Benedikt Perse, IPP
!> @brief Particle pusher based on Hamiltonian splitting for 3d3v Vlasov-Maxwell with coordinate transformation.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_time_propagator_pic_vm_3d3v_hs_trafo
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

  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use sll_m_time_propagator_pic_vm_3d3v_cl_helper, only: &
       sll_p_boundary_particles_periodic, &
       sll_p_boundary_particles_singular, &
       sll_p_boundary_particles_reflection, &
       sll_p_boundary_particles_absorption, &
       sll_s_compute_particle_boundary_simple, &
       sll_s_compute_particle_boundary_trafo, &
       sll_s_compute_particle_boundary_trafo_current

  use sll_m_initial_distribution, only: &
       sll_t_params_cos_gaussian_screwpinch

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_mpi, only: &
       mpi_sum

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base, &
       sll_t_particle_array

  use sll_m_particle_mesh_coupling_base_3d, only: &
       sll_c_particle_mesh_coupling_3d

  use sll_m_profile_functions, only: &
       sll_t_profile_functions



  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_3d3v_hs_trafo

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 3d3v with coordinate transformation
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_3d3v_hs_trafo
     class(sll_c_maxwell_3d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group
     class(sll_c_particle_mesh_coupling_3d), pointer :: particle_mesh_coupling !< Particle mesh coupling
     type( sll_t_mapping_3d ), pointer      :: map !< coordinate transformation

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

     logical :: electrostatic = .false. !< true for electrostatic simulation
     logical :: adiabatic_electrons = .false. !< true for simulation with adiabatic electrons
     logical :: lindf = .false. !< true for simulation with linear delta f method

     sll_real64 :: iter_tolerance !< iteration tolerance
     sll_int32 :: max_iter !< maximal number of iterations

     ! For control variate
     class(sll_t_control_variates), pointer :: control_variate => null() !< control variate
     sll_int32 :: i_weight = 1 !< number of weights


   contains
     procedure :: operatorHp => operatorHp_pic_vm_3d3v_hs_trafo !> Operator for H_p part
     procedure :: operatorHE => operatorHE_pic_vm_3d3v_hs_trafo !> Operator for H_E part
     procedure :: operatorHB => operatorHB_pic_vm_3d3v_hs_trafo !> Operator for H_B part
     procedure :: lie_splitting => lie_splitting_pic_vm_3d3v_hs_trafo !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_3d3v_hs_trafo !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_3d3v_hs_trafo !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_3d3v_hs_trafo !> Initialize the type
     procedure :: init_from_file => initialize_file_pic_vm_3d3v_hs_trafo !> Initialize the type
     procedure :: free => delete_pic_vm_3d3v_hs_trafo !> Finalization

  end type sll_t_time_propagator_pic_vm_3d3v_hs_trafo

contains


  !> Strang splitting
  subroutine strang_splitting_pic_vm_3d3v_hs_trafo(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_hs_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
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

  end subroutine strang_splitting_pic_vm_3d3v_hs_trafo

  !> Lie splitting
  subroutine lie_splitting_pic_vm_3d3v_hs_trafo(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_hs_trafo), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%operatorHE(dt)
          call self%operatorHp(dt)
       end do
    else
       do i_step = 1,number_steps
          call self%operatorHB(dt)
          call self%operatorHE(dt)
          call self%operatorHp(dt)
       end do
    end if

  end subroutine lie_splitting_pic_vm_3d3v_hs_trafo


  !> Lie splitting (oposite ordering)
  subroutine lie_splitting_back_pic_vm_3d3v_hs_trafo(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_hs_trafo), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: dt   !< time step
    sll_int32,                                            intent( in    ) :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    if(self%electrostatic) then
       do i_step = 1, number_steps
          call self%operatorHp(dt)
          call self%operatorHE(dt)
       end do
    else
       do i_step = 1,number_steps
          call self%operatorHp(dt)
          call self%operatorHE(dt)
          call self%operatorHB(dt)
       end do
    end if

  end subroutine lie_splitting_back_pic_vm_3d3v_hs_trafo


  !---------------------------------------------------------------------------!
  !> Push H_p: Equations to be solved
  !> $\Xi^{n+1}=\Xi^n+ \Delta t DF^{-1}(\overline{\Xi})\overline{V}$,
  !> $V^{n+1}=V^{n}+ \Delta t \mathbb {W}_{\frac{q}{m}} DF^{-\top}(\overline{\Xi}\right) \tBB(\overline{\Xi},b^n) DF^{-1}(\overline{\Xi})\overline{V}$,
  !>$M_1 e^{n+1}= M_1 e^n-\int^{t^{n+1}}_{t^n}   \mathbb{Lambda}^1(\Xi(\tau))^\top d\tau \mathbb{W}_q DF^{-1}(\overline{\Xi})\overline{V}$
  subroutine operatorHp_pic_vm_3d3v_hs_trafo(self, dt)
    class(sll_t_time_propagator_pic_vm_3d3v_hs_trafo), intent( inout ) :: self !< time propagator object 
    sll_real64,                                        intent( in    ) :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i, j, i_sp
    sll_real64 :: xi(3), xbar(3), xnew(3), xphys(3), vt(3), vi(3), vbar(3), vnew(3)  
    sll_real64 :: jmat(3,3)
    sll_real64 :: qoverm, wi(1), wall(3)
    sll_real64 :: err
    sll_real64 :: qmdt
    sll_real64 :: bfield(3), rhs(3)
    sll_real64 :: c(3)

    self%j_dofs_local = 0.0_f64
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       qmdt = qoverm*dt;
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          ! Read out particle position and velocity
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          ! Evaulate B at particle position  
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)-1], self%bfield_dofs(1:self%n_total0), bfield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)-1],self%bfield_dofs(self%n_total0+1:self%n_total0+self%n_total1), bfield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2)-1, self%spline_degree(3)],self%bfield_dofs(self%n_total0+self%n_total1+1:self%n_total0+2*self%n_total1), bfield(3))
          jmat=self%map%jacobian_matrix_inverse(xi)
          if( self%map%inverse)then
             xphys = self%map%get_x(xi)
             xbar = xphys + dt * vi
             xnew = self%map%get_xi(xbar)

             vt = (xnew - xi)/dt 
             !c = VT x bfield 
             c(1) = bfield(3)*vt(2) - bfield(2)*vt(3)
             c(2) = bfield(1)*vt(3) - bfield(3)*vt(1)
             c(3) = bfield(2)*vt(1) - bfield(1)*vt(2)

             !vnew = vi + sign * DF^{-T} * c
             do j=1,3
                vnew(j)= vi(j) + qmdt*(jmat(1,j)*c(1)+jmat(2,j)*c(2)+jmat(3,j)*c(3))
             end do
             err= maxval(abs(vi - vnew))

             i = 0
             do while(i < self%max_iter .and. err > self%iter_tolerance)
                xbar = 0.5_f64*(xnew+xi)
                vbar = 0.5_f64*(vnew+vi)
                call sll_s_compute_particle_boundary_simple( self%boundary_particles, self%counter_left, self%counter_right, xi, xbar  )
                if (xbar(1) >= 0._f64 .and. xbar(1) <= 1._f64) then
                   jmat = self%map%jacobian_matrix_inverse(xbar)
                   ! Evaulate B at particle position   
                   call self%particle_mesh_coupling%evaluate &
                        (xbar, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)-1], self%bfield_dofs(1:self%n_total0), bfield(1))
                   call self%particle_mesh_coupling%evaluate &
                        (xbar, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)-1],self%bfield_dofs(self%n_total0+1:self%n_total0+self%n_total1), bfield(2))
                   call self%particle_mesh_coupling%evaluate &
                        (xbar, [self%spline_degree(1)-1, self%spline_degree(2)-1, self%spline_degree(3)],self%bfield_dofs(self%n_total0+self%n_total1+1:self%n_total0+2*self%n_total1), bfield(3))

                   xbar = xphys + dt * vbar
                   xnew = self%map%get_xi(xbar)

                   vt = (xnew - xi)/dt 
                   !c = VT x bfield 
                   c(1) = bfield(3)*vt(2) - bfield(2)*vt(3)
                   c(2) = bfield(1)*vt(3) - bfield(3)*vt(1)
                   c(3) = bfield(2)*vt(1) - bfield(1)*vt(2)

                   !rhs = vi + sign * DF^{-T} * c
                   do j=1,3
                      rhs(j)= vi(j) + qmdt*(jmat(1,j)*c(1)+jmat(2,j)*c(2)+jmat(3,j)*c(3))
                   end do
                   err = maxval(abs(vnew - rhs))
                   vnew = rhs
                end if
                i = i + 1
             end do
             call sll_s_compute_particle_boundary_trafo( self%boundary_particles, self%counter_left, self%counter_right, self%map, xi, xnew, vnew  )
          else !Predictor-Corrector with loop for corrector step
             !Transformation of the v coordinates, !VT = DF^{-1} * vi 
             do j=1,3
                vt(j) = jmat(j,1)*vi(1) + jmat(j,2)*vi(2) + jmat(j,3)*vi(3)
             end do
             !c = VT x bfield  
             c(1)=bfield(3)*vt(2)-bfield(2)*vt(3)
             c(2)=bfield(1)*vt(3)-bfield(3)*vt(1)
             c(3)=bfield(2)*vt(1)-bfield(1)*vt(2)
             !rhs = vi + sign * DF^{-T} * c
             do j=1,3
                vnew(j)= vi(j) + qmdt*(jmat(1,j)*c(1)+jmat(2,j)*c(2)+jmat(3,j)*c(3))
             end do
             !x^\star=\mod(x^n+dt*DF^{-1}(x^n)v^n,1)
             xnew = xi + dt * vt

             err= max( maxval(abs(xi - xnew)), maxval(abs(vi - vnew)) ) 
             i = 0
             do while(i < self%max_iter .and. err > self%iter_tolerance)
                xbar = 0.5_f64*(xnew+xi)
                vbar = 0.5_f64*(vnew+vi)
                call sll_s_compute_particle_boundary_simple( self%boundary_particles, self%counter_left, self%counter_right, xi, xbar)

                call self%particle_mesh_coupling%evaluate &
                     (xbar, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)-1], self%bfield_dofs(1:self%n_total0), bfield(1))
                call self%particle_mesh_coupling%evaluate &
                     (xbar, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)-1],self%bfield_dofs(self%n_total0+1:self%n_total0+self%n_total1), bfield(2))
                call self%particle_mesh_coupling%evaluate &
                     (xbar, [self%spline_degree(1)-1, self%spline_degree(2)-1, self%spline_degree(3)],self%bfield_dofs(self%n_total0+self%n_total1+1:self%n_total0+2*self%n_total1), bfield(3))

                jmat = self%map%jacobian_matrix_inverse( xbar )
                !VT = DF^{-1} * vi
                do j = 1, 3
                   vt(j) = jmat(j,1)*vbar(1) + jmat(j,2)*vbar(2)+ jmat(j,3)*vbar(3)
                end do
                xbar = xi + dt * vt
                !c = VT x bfield
                c(1)=bfield(3)*vt(2)-bfield(2)*vt(3)
                c(2)=bfield(1)*vt(3)-bfield(3)*vt(1)
                c(3)=bfield(2)*vt(1)-bfield(1)*vt(2)
                !rhs = vi + sign * DF^{-T} * c
                do j=1,3
                   rhs(j)= vi(j) +  qmdt*(jmat(1,j)*c(1)+jmat(2,j)*c(2)+jmat(3,j)*c(3))
                end do

                err = max( maxval(abs(xnew - xbar)), maxval(abs(vnew - rhs)) )
                xnew = xbar 
                vnew = rhs
                i = i + 1
             end do
          end if
          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)
          call sll_s_compute_particle_boundary_trafo_current( self%boundary_particles, self%counter_left, self%counter_right, self%map, self%particle_mesh_coupling, self%j_dofs_local, self%spline_degree, self%rhob, xi, xnew, vnew, wi )


          call self%particle_group%group(i_sp)%set_x( i_part, xnew )
          call self%particle_group%group(i_sp)%set_v( i_part, vnew )
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             if( self%lindf .eqv. .false.) then
                wall = self%particle_group%group(i_sp)%get_weights(i_part)
                wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xnew, vnew, 0.0_f64, wall(1), wall(2) )
                call self%particle_group%group(i_sp)%set_weights( i_part, wall )
             end if
          end if
       end do
    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         self%n_total1+2*self%n_total0, MPI_SUM, self%j_dofs)

    if( self%adiabatic_electrons) then
       call self%maxwell_solver%compute_phi_from_j( self%j_dofs, self%phi_dofs, self%efield_dofs )
    else
       call self%maxwell_solver%compute_E_from_j( self%betar(2)*self%j_dofs, self%efield_dofs )
    end if

  end subroutine operatorHp_pic_vm_3d3v_hs_trafo


  !---------------------------------------------------------------------------!
  !> Push H_B: Equations to be solved
  !> $M_1 e^{n+1}=M_1 e^n+\Delta t C^\top M_2 b^n$
  subroutine operatorHB_pic_vm_3d3v_hs_trafo ( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_hs_trafo), intent( inout ) :: self !< time propagator object 
    sll_real64,                                               intent( in    ) :: dt   !< time step

    call self%maxwell_solver%compute_E_from_B( self%betar(1)*dt, self%bfield_dofs, self%efield_dofs)

  end subroutine operatorHB_pic_vm_3d3v_hs_trafo


  !---------------------------------------------------------------------------!
  !> Push H_E: Equations to be solved
  !> $V^{n+1}=V^n+\Delta t\mathbb{W}_{\frac{q}{m}} DF^{-\top}(\Xi^n) \mathbb{Lambda}^1(\Xi^n) e^n$
  !> $b^{n+1}=b^n-\Delta t C e^n$
  subroutine operatorHE_pic_vm_3d3v_hs_trafo ( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_hs_trafo), intent( inout ) :: self !< time propagator object 
    sll_real64,                                        intent( in    ) :: dt   !< time step
    !local variables
    sll_int32 :: i_part, j, i_sp
    sll_real64 :: xi(3), vi(3), jmat(3,3)
    sll_real64 :: efield(3), ephys(3), wall(3), Rx(3)
    sll_real64 :: qoverm, q, m

    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       q = self%particle_group%group(i_sp)%species%q
       m = self%particle_group%group(i_sp)%species%m
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          ! Read out particle position and velocity
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          ! Evaulate E at particle position and propagate v       
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], self%efield_dofs(1:self%n_total1), efield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)],self%efield_dofs(self%n_total1+1:self%n_total1+self%n_total0), efield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1],self%efield_dofs(self%n_total1+self%n_total0+1:self%n_total1+2*self%n_total0), efield(3))

          jmat=self%map%jacobian_matrix_inverse_transposed(xi)
          do j=1, 3
             ephys(j) = jmat(j,1)* efield(1)+jmat(j,2)* efield(2)+jmat(j,3)* efield(3)
          end do

          if( self%lindf .eqv. .false.) then
             ! velocity update 
             vi = vi + dt* qoverm * ephys
             call self%particle_group%group(i_sp)%set_v( i_part, vi )

             ! Update particle weights
             if (self%particle_group%group(i_sp)%n_weights == 3 ) then
                wall = self%particle_group%group(i_sp)%get_weights(i_part)
                wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xi, vi, 0.0_f64, wall(1), wall(2) )
                call self%particle_group%group(i_sp)%set_weights( i_part, wall )
             end if
          else
             ! Update particle weights
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             select type(p => self%control_variate%cv(i_sp)%control_variate_distribution_params)
             type is(sll_t_params_cos_gaussian_screwpinch)
                Rx = xi 
                Rx(1) = self%map%get_x1(xi) + m*vi(2)/q
                Rx(1) = (Rx(1)-self%x_min(1))/self%Lx(1)

                wall(3) = wall(3) + dt* (q/p%profile%T_i(Rx(1)) *sum(ephys*vi) -&
                     ephys(2)*(p%profile%drho_0(Rx(1))/p%profile%rho_0(Rx(1))+(0.5_f64*m*sum(vi**2)/p%profile%T_i(Rx(1)) - 1.5_f64)* p%profile%dT_i(Rx(1))/p%profile%T_i(Rx(1)) )  ) *&
                     p%eval_v_density(vi, Rx, m)/p%eval_v_density(vi, xi, m)* self%map%jacobian(xi)
             class default
                wall(3) = wall(3) + dt* qoverm* sum(ephys*vi) * self%map%jacobian(xi)
             end select
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do
    end do

    if(self%electrostatic .eqv. .false.) then
       call self%maxwell_solver%compute_B_from_E( dt, self%efield_dofs, self%bfield_dofs)
    end if

  end subroutine operatorHE_pic_vm_3d3v_hs_trafo


  !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_3d3v_hs_trafo(&
       self, &
       maxwell_solver, &
       particle_mesh_coupling, &
       particle_group, &
       phi_dofs, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       map, &
       boundary_particles, &
       iter_tolerance, max_iter, &
       betar, &
       electrostatic, &
       rhob, &
       control_variate, &
       lindf) 
    class(sll_t_time_propagator_pic_vm_3d3v_hs_trafo), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,          intent( in ) :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target, intent(in) :: particle_mesh_coupling !< Particle mesh coupling
    class(sll_t_particle_array), target,           intent( in ) :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:) !< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent( in ) :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent( in ) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                    intent( in ) :: x_min(3) !< Lower bound of x domain
    sll_real64,                                    intent( in ) :: Lx(3) !< Length of the domain in x direction.
    type(sll_t_mapping_3d), target,                intent( inout ) :: map !< Coordinate transformation
    sll_int32, optional,                           intent( in ) :: boundary_particles !< particle boundary conditions
    sll_real64, optional,                          intent( in ) :: iter_tolerance !< iteration tolerance
    sll_int32,  optional,                          intent( in ) :: max_iter !< maximal number of iterations
    sll_real64, optional,                          intent( in ) :: betar(2) !< reciprocal plasma beta
    logical, optional,                             intent( in ) :: electrostatic !< true for electrostatic simulation
    sll_real64, optional, target,                  intent( in ) :: rhob(:) !< charge at the boundary
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    logical, optional, intent(in) :: lindf !< true for linear delta f method
    !local variables
    sll_int32 :: ierr

    if (present(iter_tolerance) )  then
       self%iter_tolerance = iter_tolerance
       self%max_iter = max_iter
    else
       self%iter_tolerance = 1d-12
       self%max_iter = 10
    end if

    if( present(electrostatic) )then
       self%electrostatic = electrostatic
    end if

    if (present(boundary_particles) ) then
       self%boundary_particles = boundary_particles
    end if

    if( present(rhob) )then
       self%rhob => rhob
    end if

    if( particle_group%group(1)%species%q > 0._f64) self%adiabatic_electrons = .true.

    self%maxwell_solver => maxwell_solver
    self%particle_mesh_coupling => particle_mesh_coupling
    self%particle_group => particle_group
    self%phi_dofs => phi_dofs
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs
    self%n_total0 = self%maxwell_solver%n_total0
    self%n_total1 = self%maxwell_solver%n_total1
    self%spline_degree = self%particle_mesh_coupling%spline_degree
    self%x_min = x_min
    self%x_max = x_min + Lx
    self%Lx = Lx
    self%map => map

    SLL_ALLOCATE( self%j_dofs(1:self%n_total1+self%n_total0*2), ierr )
    SLL_ALLOCATE( self%j_dofs_local(1:self%n_total1+self%n_total0*2), ierr )
    self%j_dofs = 0.0_f64
    self%j_dofs_local = 0.0_f64

    if (present(betar)) then
       self%betar = betar
    else
       self%betar = 1.0_f64
    end if

    if (present(control_variate)) then
       allocate(self%control_variate )
       allocate(self%control_variate%cv(self%particle_group%n_species) )
       self%control_variate => control_variate
       self%i_weight = 3
       if (present(lindf)) then
          self%lindf = lindf
       end if
    end if

  end subroutine initialize_pic_vm_3d3v_hs_trafo


  !> Constructor.
  subroutine initialize_file_pic_vm_3d3v_hs_trafo(&
       self, &
       maxwell_solver, &
       particle_mesh_coupling, &
       particle_group, &
       phi_dofs, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       map, &
       filename, &
       boundary_particles, &
       betar, &
       electrostatic, &
       rhob, &
       control_variate, &
       lindf) 
    class(sll_t_time_propagator_pic_vm_3d3v_hs_trafo), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,          intent( in ) :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target, intent(in) :: particle_mesh_coupling !< Particle mesh coupling
    class(sll_t_particle_array), target,           intent( in ) :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:) !< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent( in ) :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent( in ) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                    intent( in ) :: x_min(3) !< Lower bound of x domain
    sll_real64,                                    intent( in ) :: Lx(3) !< Length of the domain in x direction.
    type(sll_t_mapping_3d), target,                intent( inout ) :: map !< Coordinate transformation
    character(len=*),                              intent( in ) :: filename !< file name
    sll_int32, optional,                           intent( in ) :: boundary_particles !< particle boundary conditions
    sll_real64, optional,                          intent( in ) :: betar(2) !< reciprocal plasma beta
    logical, optional    :: electrostatic !< true for electrostatic simulation
    sll_real64, optional, target,                  intent( in ) :: rhob(:) !< charge at the boundary
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    logical, optional, intent(in) :: lindf !< true for linear delta f method
    !local variables
    sll_int32 :: input_file
    sll_int32 :: io_stat
    sll_real64 :: iter_tolerance
    sll_int32 :: max_iter, boundary_particles_set
    logical :: electrostatic_set
    logical :: lindf_set
    sll_real64 :: betar_set(2)

    namelist /time_iterate/ iter_tolerance, max_iter

    if( present(boundary_particles) )then
       boundary_particles_set = boundary_particles
    else
       boundary_particles_set = 100
    end if

    if( present(electrostatic) )then
       electrostatic_set = electrostatic
    else
       electrostatic_set = .false.
    end if

    if( present(lindf) )then
       lindf_set = lindf
    else
       lindf_set = .false.
    end if

    if( present(betar) )then
       betar_set = betar
    else
       betar_set = 1._f64
    end if

    ! Read in solver tolerance
    open(newunit = input_file, file=filename, status='old',IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'sll_m_time_propagator_pic_vm_3d3v_hs_trafo: Input file does not exist. Set default tolerance.'
       if( present( control_variate ) )then
          call self%init( maxwell_solver, &
               particle_mesh_coupling, &
               particle_group, &
               phi_dofs,&
               efield_dofs, &
               bfield_dofs, &
               x_min, &
               Lx, &
               map, &
               boundary_particles_set, &
               betar=betar_set,&
               electrostatic = electrostatic_set, &
               rhob = rhob, &
               control_variate = control_variate, &
               lindf=lindf_set)
       else
          call self%init( maxwell_solver, &
               particle_mesh_coupling, &
               particle_group, &
               phi_dofs,&
               efield_dofs, &
               bfield_dofs, &
               x_min, &
               Lx, &
               map, &
               boundary_particles_set, &
               betar=betar_set,&
               electrostatic = electrostatic_set, &
               rhob = rhob)
       end if
    else       
       read(input_file, time_iterate,IOStat=io_stat)
       if (io_stat /= 0 ) then
          if( present( control_variate ) )then
             call self%init( maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs,&
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  map, &
                  boundary_particles_set, &
                  betar=betar_set,&
                  electrostatic = electrostatic_set, &
                  rhob = rhob, &
                  control_variate = control_variate, &
                  lindf=lindf_set)
          else
             call self%init( maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs,&
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  map, &
                  boundary_particles_set, &
                  betar=betar_set,&
                  electrostatic = electrostatic_set, &
                  rhob = rhob)
          end if
       else
          if( present( control_variate ) )then
             call self%init(   maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs,&
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  map, &
                  boundary_particles_set, &
                  iter_tolerance, max_iter, &
                  betar=betar_set,&
                  electrostatic = electrostatic_set, &
                  rhob = rhob, &
                  control_variate = control_variate, &
                  lindf=lindf_set)
          else
             call self%init(   maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs,&
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  map, &
                  boundary_particles_set, &
                  iter_tolerance, max_iter, &
                  betar=betar_set,&
                  electrostatic = electrostatic_set, &
                  rhob = rhob)
          end if
       end if
       close(input_file)
    end if


  end subroutine initialize_file_pic_vm_3d3v_hs_trafo


  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_3d3v_hs_trafo(self)
    class(sll_t_time_propagator_pic_vm_3d3v_hs_trafo), intent( inout ) :: self !< time propagator object 

    deallocate( self%j_dofs )
    deallocate( self%j_dofs_local )

    self%maxwell_solver => null()
    self%particle_mesh_coupling => null()
    self%particle_group => null()
    self%phi_dofs => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()
    self%map => null()



  end subroutine delete_pic_vm_3d3v_hs_trafo



end module sll_m_time_propagator_pic_vm_3d3v_hs_trafo
