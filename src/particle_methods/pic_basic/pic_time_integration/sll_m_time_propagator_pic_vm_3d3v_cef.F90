!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on Hamiltonian splitting using in Crouseilles, Einkemmer, Faou for 3d3v Vlasov-Maxwell.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
!> Reference: Kraus, Kormann, Sonnendrücker, Morrison: GEMPIC: Geometric ElectroMagnetic Particle-In-Cell Methods
module sll_m_time_propagator_pic_vm_3d3v_cef
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
       sll_p_boundary_particles_absorption

  use sll_m_initial_distribution, only: &
       sll_t_params_cos_gaussian_screwpinch, &
       sll_t_params_noise_gaussian

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_mpi, only: &
       mpi_sum

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use sll_m_particle_mesh_coupling_base_3d, only: &
       sll_c_particle_mesh_coupling_3d

  use sll_m_profile_functions, only: &
       sll_t_profile_functions


  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_3d3v_cef

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 3d3v
  type, extends(sll_c_time_propagator_base) :: sll_t_time_propagator_pic_vm_3d3v_cef
     class(sll_c_maxwell_3d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_particle_mesh_coupling_3d), pointer :: particle_mesh_coupling
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

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

     ! For control variate
     class(sll_t_control_variates), pointer :: control_variate => null()
     sll_int32 :: i_weight = 1 !< number of weights


   contains
     procedure :: operatorHp => operatorHp_pic_vm_3d3v  !> Operator for H_p part
     procedure :: operatorHE => operatorHE_pic_vm_3d3v  !> Operator for H_E part
     procedure :: operatorHB => operatorHB_pic_vm_3d3v  !> Operator for H_B part
     procedure :: lie_splitting => lie_splitting_pic_vm_3d3v !> Lie splitting propagator
     procedure :: lie_splitting_back => lie_splitting_back_pic_vm_3d3v !> Lie splitting propagator
     procedure :: strang_splitting => strang_splitting_pic_vm_3d3v !> Strang splitting propagator

     procedure :: init => initialize_pic_vm_3d3v !> Initialize the type
     procedure :: free => delete_pic_vm_3d3v !> Finalization

  end type sll_t_time_propagator_pic_vm_3d3v_cef

contains


  !> Strang splitting
  subroutine strang_splitting_pic_vm_3d3v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_cef), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    do i_step = 1, number_steps
       call self%operatorHB(0.5_f64*dt)
       call self%operatorHE(0.5_f64*dt)
       call self%operatorHp(dt)
       call self%operatorHE(0.5_f64*dt)
       call self%operatorHB(0.5_f64*dt)
    end do

  end subroutine strang_splitting_pic_vm_3d3v


  !> Lie splitting
  subroutine lie_splitting_pic_vm_3d3v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_cef), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    do i_step = 1,number_steps
       call self%operatorHB(dt)
       call self%operatorHE(dt)
       call self%operatorHp(dt)
    end do


  end subroutine lie_splitting_pic_vm_3d3v


  !> Lie splitting (oposite ordering)
  subroutine lie_splitting_back_pic_vm_3d3v(self,dt, number_steps)
    class(sll_t_time_propagator_pic_vm_3d3v_cef), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_int32,                                      intent(in)    :: number_steps !< number of time steps
    !local variables
    sll_int32 :: i_step

    do i_step = 1,number_steps
       call self%operatorHp(dt)
       call self%operatorHE(dt)
       call self%operatorHB(dt)
    end do

  end subroutine lie_splitting_back_pic_vm_3d3v


  !---------------------------------------------------------------------------!
  !> Push H_p: Equations to be solved
  !> $X^{n+1}=X^n+\Delta t V^n$
  !> $M_1 e^n= M_1 e^n- \int_{t^n}^{t^{n+1}} \mathbb{\Lambda}^1(X(\tau))^\top  d\tau \mathbb{W}_q  V^n$
  subroutine operatorHp_pic_vm_3d3v(self, dt)
    class(sll_t_time_propagator_pic_vm_3d3v_cef), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: xnew(3), vi(3), wi(1), xi(3), wall(3)

    self%j_dofs_local = 0.0_f64
    do i_sp = 1, self%particle_group%n_species
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          ! Read out particle position and velocity
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          ! Then update particle position:  X_new = X_old + dt * V
          xnew = xi + dt * vi

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

          call compute_particle_boundary( self, xi, xnew, vi, wi)

          call self%particle_group%group(i_sp)%set_x(i_part, xnew)
          call self%particle_group%group(i_sp)%set_v(i_part, vi)
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xnew, vi, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do
    end do
    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         self%n_total1+2*self%n_total0, MPI_SUM, self%j_dofs )


    if( self%adiabatic_electrons) then
       call self%maxwell_solver%compute_phi_from_j( self%j_dofs, self%phi_dofs, self%efield_dofs )
    else
       call self%maxwell_solver%compute_E_from_j( self%j_dofs, self%efield_dofs )
    end if
  end subroutine operatorHp_pic_vm_3d3v


  !> Helper function for \a operatorHp
  subroutine compute_particle_boundary( self, xold, xnew, vi, wi )
    class(sll_t_time_propagator_pic_vm_3d3v_cef), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: xold(3)
    sll_real64,                                           intent( inout ) :: xnew(3)
    sll_real64,                                           intent( inout ) :: vi(3)
    sll_real64,                                           intent( in    ) :: wi(1)
    !local variables
    sll_real64 :: xmid(3), xbar, dx, vh(3)

    if(xnew(1) < -self%x_max(1) .or.  xnew(1) > 2._f64*self%x_max(1) ) then
       print*, xnew
       SLL_ERROR("particle boundary", "particle out of bound")
    else if(xnew(1) < self%x_min(1) .or. xnew(1) > self%x_max(1) )then
       if(xnew(1) < self%x_min(1)  )then
          xbar = self%x_min(1)
          self%counter_left = self%counter_left+1
       else if(xnew(1) > self%x_max(1))then
          xbar = self%x_max(1)
          self%counter_right = self%counter_right+1
       end if

       dx = (xbar- xold(1))/(xnew(1)-xold(1))
       xmid = xold + dx * (xnew-xold)
       xmid(1) = xbar
       vh = (xmid - xold)*wi(1)
       call self%particle_mesh_coupling%add_current( xold, xmid, vh, self%j_dofs_local )
       select case(self%boundary_particles)
       case(sll_p_boundary_particles_reflection) 
          xnew(1) = 2._f64*xbar-xnew(1)
          vi(1) = -vi(1)
       case(sll_p_boundary_particles_absorption)
          call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
       case( sll_p_boundary_particles_periodic)
          call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
          xnew(1) = self%x_min(1) + modulo(xnew(1)-self%x_min(1), self%Lx(1))
          xmid(1) = self%x_max(1)+self%x_min(1)-xbar
          call self%particle_mesh_coupling%add_charge(xmid, -wi(1), self%spline_degree, self%rhob)
       case default
          xnew(1) = self%x_min(1) + modulo(xnew(1)-self%x_min(1), self%Lx(1))
          xmid(1) = self%x_max(1)+self%x_min(1)-xbar
       end select
    else
       xmid = xold
    end if
    vh = (xnew - xmid)*wi(1)
    call self%particle_mesh_coupling%add_current( xmid, xnew, vh, self%j_dofs_local )
    xnew(2:3) = self%x_min(2:3) + modulo(xnew(2:3)-self%x_min(2:3), self%Lx(2:3))

  end subroutine compute_particle_boundary


  !---------------------------------------------------------------------------!
  !> Push H_B: Equations to be solved
  !> $(\mathbb{I}-\Delta \frac{\Delta t q}{2 m}  \mathbb{B}(X^n,b^n) V^{n+1}=(\mathbb{I}+ \frac{\Delta t q}{2 m} \mathbb{B}(X^n,b^n) ) V^n$
  !> $M_1 e^{n+1}=M_1 e^n+\Delta t C^\top M_2 b^n$
  subroutine operatorHB_pic_vm_3d3v(self, dt)
    class(sll_t_time_propagator_pic_vm_3d3v_cef), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32  :: i_part, i_sp
    sll_real64 :: qmdt
    sll_real64 :: vi(3), vnew(3), xi(3)
    sll_real64 :: bfield(3),  bsq(3), denom, wall(3)

    do i_sp = 1, self%particle_group%n_species
       ! Advect dV/dt=VxB
       qmdt = self%particle_group%group(i_sp)%species%q_over_m()*dt*0.5_f64;
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)-1], self%bfield_dofs(1:self%n_total0), bfield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)-1],self%bfield_dofs(self%n_total0+1:self%n_total0+self%n_total1), bfield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2)-1, self%spline_degree(3)],self%bfield_dofs(self%n_total0+self%n_total1+1:self%n_total0+2*self%n_total1), bfield(3))

          bfield = qmdt*bfield
          bsq = bfield**2       
          denom = sum(bsq)+1.0_f64

          vnew(1) = ( (bsq(1)-bsq(2)-bsq(3)+1.0_f64)*vi(1) + &
               2.0_f64*(bfield(1)*bfield(2)+bfield(3))*vi(2) + &
               2.0_f64*(bfield(1)*bfield(3)-bfield(2))*vi(3) )/denom
          vnew(2) = ( 2.0_f64*(bfield(1)*bfield(2)-bfield(3))*vi(1) + &
               (-bsq(1)+bsq(2)-bsq(3)+1.0_f64)*vi(2) + &
               2.0_f64*(bfield(2)*bfield(3)+bfield(1))*vi(3) )/denom
          vnew(3) = ( 2.0_f64*(bfield(1)*bfield(3)+bfield(2))*vi(1) + &
               2.0_f64*(bfield(2)*bfield(3)-bfield(1))*vi(2) + &
               (-bsq(1)-bsq(2)+bsq(3)+1.0_f64)*vi(3)  )/denom

          call self%particle_group%group(i_sp)%set_v( i_part, vnew )
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xi, vnew, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do
    end do

    if(self%electrostatic .eqv. .false.) then
       ! Update efield2
       call self%maxwell_solver%compute_E_from_B(&
            self%betar(1)*dt, self%bfield_dofs, self%efield_dofs)
    end if
  end subroutine operatorHB_pic_vm_3d3v


  !---------------------------------------------------------------------------!
  !> Push H_E: Equations to be solved
  !> $V^{n+1}=V^n+\Delta t\mathbb{W}_{\frac{q}{m}} \mathbb{Lambda}^1(X^n) e^n$
  !> $b^{n+1}=b^n-\Delta t C e^n$
  subroutine operatorHE_pic_vm_3d3v(self, dt)
    class(sll_t_time_propagator_pic_vm_3d3v_cef), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: vnew(3), xi(3), wall(self%i_weight), Rx(3)
    sll_real64 :: efield(3)
    sll_real64 :: qoverm, q, m
    sll_int32 :: i_weight

    i_weight = self%i_weight

    ! TODO: Optimize using evaluate multiple
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       q = self%particle_group%group(i_sp)%species%q
       m = self%particle_group%group(i_sp)%species%m
       ! V_new = V_old + dt * E
       do i_part=1,self%particle_group%group(i_sp)%n_particles
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vnew = self%particle_group%group(i_sp)%get_v(i_part)
          ! Evaluate efields at particle position
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], self%efield_dofs(1:self%n_total1), efield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)],self%efield_dofs(self%n_total1+1:self%n_total1+self%n_total0), efield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1],self%efield_dofs(self%n_total1+self%n_total0+1:self%n_total1+2*self%n_total0), efield(3))
          if( self%lindf .eqv. .false.) then
             vnew = vnew + dt* qoverm * efield
             call self%particle_group%group(i_sp)%set_v(i_part, vnew)
             ! Update particle weights
             if (self%particle_group%group(i_sp)%n_weights == 3 ) then
                wall = self%particle_group%group(i_sp)%get_weights(i_part)
                wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xi, vnew, 0.0_f64, wall(1), wall(2) )
                call self%particle_group%group(i_sp)%set_weights( i_part, wall )
             end if
          else
             ! Update particle weights
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             select type(p => self%control_variate%cv(i_sp)%control_variate_distribution_params)
             type is(sll_t_params_cos_gaussian_screwpinch)
                Rx = xi
                Rx(1) = Rx(1) + m*vnew(2)/q
                Rx = (Rx-self%x_min)/self%Lx

                wall(1) = wall(1) + dt* (q/p%profile%T_i(Rx(1)) * sum(efield * vnew)-&
                     efield(2)*(p%profile%drho_0(Rx(1))/p%profile%rho_0(Rx(1))+(0.5_f64*m*sum(vnew**2)/p%profile%T_i(Rx(1)) - 1.5_f64)* p%profile%dT_i(Rx(1))/p%profile%T_i(Rx(1))  ) ) *&
                     self%control_variate%cv(i_sp)%control_variate(Rx, vnew, 0._f64)/p%eval_v_density(vnew, xi, m)*product(self%Lx)
             type is(sll_t_params_noise_gaussian)
                Rx = xi
                wall(1) = wall(1) + dt* (q/p%profile%T_i(Rx(1)) * sum(efield * vnew)-&
                     efield(2)*(p%profile%drho_0(Rx(1))/p%profile%rho_0(Rx(1))+(0.5_f64*m*sum(vnew**2)/p%profile%T_i(Rx(1)) - 1.5_f64)* p%profile%dT_i(Rx(1))/p%profile%T_i(Rx(1))  ) ) *product(self%Lx)
             class default
                wall(1) = wall(1) + dt* qoverm* sum(efield * vnew) *product(self%Lx)
             end select
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do
    end do

    if(self%electrostatic .eqv. .false.) then
       ! Update bfield
       call self%maxwell_solver%compute_B_from_E( &
            dt, self%efield_dofs, self%bfield_dofs)
    end if
  end subroutine operatorHE_pic_vm_3d3v


  !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_3d3v(&
       self, &
       maxwell_solver, &
       particle_mesh_coupling, &
       particle_group, &
       phi_dofs, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       boundary_particles, &
       betar, &
       electrostatic, &
       rhob, &
       control_variate)
    class(sll_t_time_propagator_pic_vm_3d3v_cef), intent(out) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,           intent(in)  :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target,   intent(in)  :: particle_mesh_coupling !< Particle mesh coupling
    class(sll_t_particle_array), target,      intent(in)  :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:) !< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent(in)  :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min(3) !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx(3) !< Length of the domain in x direction.
    sll_int32, optional,                           intent( in ) :: boundary_particles !< particle boundary conditions
    sll_real64, optional, intent(in) :: betar(2) !< reciprocal plasma beta
    logical, optional, intent(in)   :: electrostatic !< true for electrostatic simulation
    sll_real64, optional, target, intent( in ) :: rhob(:) !< charge at the boundary
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    !local variables
    sll_int32 :: ierr

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

    SLL_ALLOCATE( self%j_dofs(1:self%n_total1+self%n_total0*2), ierr )
    SLL_ALLOCATE( self%j_dofs_local(1:self%n_total1+self%n_total0*2), ierr )

    self%spline_degree = self%particle_mesh_coupling%spline_degree
    self%x_min = x_min
    self%x_max = x_min + Lx
    self%Lx = Lx


    if (present(control_variate)) then
       allocate(self%control_variate )
       allocate(self%control_variate%cv(self%particle_group%n_species) )
       self%control_variate => control_variate
       if(self%particle_group%group(1)%n_weights == 1 ) self%lindf = .true.
       self%i_weight = self%particle_group%group(1)%n_weights
    end if

    if (present(betar)) then
       self%betar = betar!32.89_f64
    else
       self%betar = 1.0_f64
    end if

  end subroutine initialize_pic_vm_3d3v


  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_3d3v(self)
    class(sll_t_time_propagator_pic_vm_3d3v_cef), intent( inout ) :: self !< time propagator object 

    deallocate(self%j_dofs)
    deallocate(self%j_dofs_local)
    self%maxwell_solver => null()
    self%particle_mesh_coupling => null()
    self%particle_group => null()
    self%phi_dofs => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()

  end subroutine delete_pic_vm_3d3v


end module sll_m_time_propagator_pic_vm_3d3v_cef
