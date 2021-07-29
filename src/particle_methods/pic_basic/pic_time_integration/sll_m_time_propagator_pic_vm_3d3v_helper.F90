!> @ingroup pic_time_integration
!> @author Katharina Kormann, IPP
!> @brief Particle pusher based on antisymmetric splitting with AVF for 3d3v Vlasov-Maxwell.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_time_propagator_pic_vm_3d3v_helper
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
       sll_f_get_collective_rank, &
       sll_o_collective_allreduce, &
       sll_v_world_collective

  use sll_m_control_variate, only: &
       sll_t_control_variates

  use sll_m_filter_base_3d, only: &
       sll_c_filter_base_3d

  use sll_m_time_propagator_pic_vm_3d3v_cl_helper, only: &
       sll_p_boundary_particles_periodic, &
       sll_p_boundary_particles_singular, &
       sll_p_boundary_particles_reflection, &
       sll_p_boundary_particles_absorption

  use sll_m_initial_distribution, only: &
       sll_t_params_cos_gaussian_screwpinch

  use sll_m_linear_operator_block, only : &
       sll_t_linear_operator_block

  use sll_m_particle_mass_3d_base, only: &
       sll_c_particle_mass_3d_base

  use sll_m_linear_operator_particle_mass_3d_diag, only : &
       sll_t_linear_operator_particle_mass_3d_diag

  use sll_m_linear_operator_particle_mass_cl_3d_diag, only : &
       sll_t_linear_operator_particle_mass_cl_3d_diag

  use sll_m_linear_operator_schur_ev_3d, only : &
       sll_t_linear_operator_schur_ev_3d

  use sll_m_linear_operator_schur_phiv_3d, only : &
       sll_t_linear_operator_schur_phiv_3d

  use sll_m_linear_solver_cg, only : &
       sll_t_linear_solver_cg

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_mpi, only: &
       mpi_sum, &
       mpi_max

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use sll_m_particle_mesh_coupling_base_3d, only: &
       sll_c_particle_mesh_coupling_3d

  use sll_m_preconditioner_fft, only : &
       sll_t_preconditioner_fft

  use sll_m_preconditioner_singular, only : &
       sll_t_preconditioner_singular

  use sll_m_profile_functions, only: &
       sll_t_profile_functions

  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_3d3v_helper

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Helper for implicit time propagators for Vlasov-Maxwell 3d3v
  type :: sll_t_time_propagator_pic_vm_3d3v_helper
     class(sll_c_maxwell_3d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_particle_mesh_coupling_3d), pointer :: particle_mesh_coupling
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

     type( sll_t_linear_operator_schur_phiv_3d ) :: linear_op_schur_phiv !< Schur complement operator for advect_e
     type( sll_t_linear_solver_cg )         :: linear_solver_schur_phiv !< Schur complement solver for advect_e
     type(sll_t_preconditioner_fft)            :: preconditioner_fft !< Preconditioner for Schur ev complement solver
     type(sll_t_preconditioner_singular) :: preconditioner1 !< preconditioner for mass matrices
     type( sll_t_linear_operator_schur_ev_3d ) :: linear_op_schur_ev !< Schur complement operator for advect_e
     type( sll_t_linear_solver_cg )            :: linear_solver_schur_ev  !< Schur complement solver for advect_e
     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_1 !< Particle mass
     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_2 !< Particle mass
     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_3 !< Particle mass
     type( sll_t_linear_operator_block )    :: particle_mass_op !< Particle mass operator


     sll_int32 :: spline_degree(3) !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx(3) !< Size of the domain
     sll_real64 :: x_min(3) !< Lower bound for x domain
     sll_real64 :: x_max(3) !< Upper bound for x domain
     sll_int32 :: n_total0 !< total number of Dofs for 0form
     sll_int32 :: n_total1 !< total number of Dofs for 1form
     sll_int32 :: nspan(3) !< Number of intervals where spline tensor product is non zero 
     sll_real64 :: betar(2) !< reciprocal of plasma beta

     sll_real64, pointer     :: phi_dofs(:) !< DoFs describing the scalar potential
     sll_real64, pointer     :: efield_dofs(:) !< DoFs describing the three components of the electric field
     sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the three components of the magnetic field
     sll_real64, allocatable :: j_dofs(:)      !< DoFs for representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:)!< MPI-processor local part of one component of \a j_dofs
     sll_real64, allocatable :: particle_mass_1_local(:,:) !< MPI-processor local part of one component of \a particle_mass
     sll_real64, allocatable :: particle_mass_2_local(:,:) !< MPI-processor local part of one component of \a particle_mass
     sll_real64, allocatable :: particle_mass_3_local(:,:) !< MPI-processor local part of one component of \a particle_mass

     sll_real64 :: solver_tolerance !< solver tolerance
     sll_real64 :: iter_tolerance  !< iteration tolerance 
     sll_int32  :: max_iter !< maximal amount of iterations

     sll_int32 :: boundary_particles = 100 !< particle boundary conditions
     sll_int32 :: counter_left = 0 !< boundary counter
     sll_int32 :: counter_right = 0 !< boundary counter
     sll_real64, pointer     :: rhob(:) => null() !< charge at the boundary
     sll_real64 :: force_sign = 1._f64 !< sign of particle force
     logical :: adiabatic_electrons = .false. !< true for simulation with adiabatic electrons
     logical :: jmean = .false. !< logical for mean value of current 
     logical :: lindf = .false. !< true for simulation with linear delta f method
     logical :: boundary = .false. !< true for non periodic boundary conditions

     sll_int32 :: n_particles_max ! Helper variable for size of temporary particle arrays
     sll_real64, allocatable :: xnew(:,:,:) !< extra data for particle position
     sll_real64, allocatable :: vnew(:,:,:) !< extra data for particle velocity
     sll_real64, allocatable :: efield_dofs_new(:) !< extra data for efield DoFs
     sll_real64, allocatable :: efield_dofs_work(:) !< extra data for efield DoFs
     sll_real64, allocatable :: phi_dofs_new(:) !< extra data for scalar potential DoFs
     sll_real64, allocatable :: phi_dofs_work(:) !< extra data for scalar potential DoFs

     sll_int32 :: n_failed = 0 !< number of failed iterations
     sll_int32 :: iter_counter = 0 !< number of iterations
     sll_int32 :: niter(100000) !< list of iteration numbers

     ! For control variate
     class(sll_t_control_variates), pointer :: control_variate => null() !< control variate
     sll_int32 :: i_weight = 1 !< number of weights

     !filtering
     sll_real64, allocatable     :: efield_filter(:) !< DoFs describing the two components of the electric field
     sll_real64, allocatable     :: bfield_filter(:)   !< DoFs describing the magnetic field
     class(sll_c_filter_base_3d), pointer :: filter !< filter

   contains
     procedure :: advect_x => advect_x_pic_vm_3d3v_avf !> Advect x-part
     procedure :: advect_eb => advect_eb_schur_pic_vm_3d3v_avf !> Solve Faraday and B-part of Ampere using Schur complement
     procedure :: advect_vb => advect_vb_pic_vm_3d3v_avf !> Push v, vxB-part only
     procedure :: advect_e => advect_e_pic_vm_3d3v_avf !> Advect ev-part
     procedure :: advect_e_trunc => advect_e_pic_vm_3d3v_avf_trunc !> Advect ev-part with truncated particle mass
     procedure :: advect_ex !> Advect ev-part together with x-part in nonlinear iteration

     procedure :: init => initialize_pic_vm_3d3v !> Initialize the type
     procedure :: init_from_file => initialize_file_pic_vm_3d3v !> Initialize the type
     procedure :: free => delete_pic_vm_3d3v !> Finalization

  end type sll_t_time_propagator_pic_vm_3d3v_helper

contains


  !---------------------------------------------------------------------------!
  !> advect_x: Equations to be solved
  !> $X^{n+1}=X^n + \Delta t V^n
  subroutine advect_x_pic_vm_3d3v_avf ( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_helper), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: xi(3), xnew(3), vi(3), wall(3)

    do i_sp = 1, self%particle_group%n_species
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          xnew = xi + dt * vi
          call compute_particle_boundary( self, xi, xnew, vi )


          call self%particle_group%group(i_sp)%set_x ( i_part, xnew )
          call self%particle_group%group(i_sp)%set_v ( i_part, vi )
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xnew, vi, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do
    end do

  end subroutine advect_x_pic_vm_3d3v_avf


  !> Helper function for \a advect_x
  subroutine compute_particle_boundary( self, xold, xnew, vi  )
    class(sll_t_time_propagator_pic_vm_3d3v_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( inout ) :: xold(3)
    sll_real64,                                           intent( inout ) :: xnew(3)
    sll_real64,                                           intent( inout ) :: vi(3)
    !local variables
    sll_real64 ::xbar

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

       select case(self%boundary_particles)
       case(sll_p_boundary_particles_reflection)
          xnew(1) = 2._f64*xbar-xnew(1)
          vi(1) = -vi(1)
       case(sll_p_boundary_particles_absorption)
       case( sll_p_boundary_particles_periodic)
          xnew(1) = self%x_min(1) + modulo(xnew(1)-self%x_min(1), self%Lx(1))
       case default
          xnew(1) = self%x_min(1) + modulo(xnew(1)-self%x_min(1), self%Lx(1))
       end select
    end if
    xnew(2:3) = self%x_min(2:3) + modulo(xnew(2:3)-self%x_min(2:3), self%Lx(2:3))

  end subroutine compute_particle_boundary


  !---------------------------------------------------------------------------!
  !> advect_vb: Equations to be solved
  !> $(\mathbb{I}-\Delta \frac{\Delta t q}{2 m}  \mathbb{B}(X^n,b^n) V^{n+1}=(\mathbb{I}+ \frac{\Delta t q}{2 m} \mathbb{B}(X^n,b^n) V^n$
  subroutine advect_vb_pic_vm_3d3v_avf ( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_helper), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    !local variables
    sll_int32  :: i_part, i_sp
    sll_real64 :: qmdt
    sll_real64 :: vi(3), vnew(3), xi(3), wall(3)
    sll_real64 :: bfield(3),  bsq(3), denom

    self%bfield_filter = self%bfield_dofs
    call self%filter%apply_inplace(self%bfield_filter(1:self%n_total0))
    call self%filter%apply_inplace(self%bfield_filter(self%n_total0+1:self%n_total0+self%n_total1))
    call self%filter%apply_inplace(self%bfield_filter(self%n_total0+self%n_total1+1:self%n_total0+2*self%n_total1))

    do i_sp = 1, self%particle_group%n_species
       qmdt = self%particle_group%group(i_sp)%species%q_over_m()*dt*0.5_f64;
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)-1], self%bfield_filter(1:self%n_total0), bfield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)-1],self%bfield_filter(self%n_total0+1:self%n_total0+self%n_total1), bfield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2)-1, self%spline_degree(3)],self%bfield_filter(self%n_total0+self%n_total1+1:self%n_total0+2*self%n_total1), bfield(3))


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

  end subroutine advect_vb_pic_vm_3d3v_avf


  !---------------------------------------------------------------------------!
  !> advect_eb: Equations to be solved
  !> Solution with Schur complement: $ S=M_1+\frac{\Delta t^2}{4} C^\top M_2 C $
  !> $ e^{n+1}=S^{-1}( (M_1-\frac{\Delta t^2}{4} C^\top M_2 C)e^n+\Delta t C^\top M_2 b^n) $
  !> $ b^{n+1}=b^n-\frac{\Delta t}{2} C(e^n+e^{n+1}) $
  subroutine advect_eb_schur_pic_vm_3d3v_avf ( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_helper), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    call self%maxwell_solver%compute_curl_part( dt, self%efield_dofs, self%bfield_dofs, self%betar(1) )

  end subroutine advect_eb_schur_pic_vm_3d3v_avf


  !---------------------------------------------------------------------------!
  !> advect_e: Equations to be solved
  !> Solution with Schur complement: $ S_{+}=M_1+\frac{\Delta t^2 q^2}{4 m} (\mathbb{\Lambda}^1)^T \mathbb{\Lambda}^1 $
  !> $e^{n+1}=S_{+}^{-1}\left(S_{-}e^n-\Delta t (\mathbb{\Lambda}^1)^\top \mathbb{W}_q V^n \right)$
  !> $V^{n+1}=V^n+\frac{\Delta t}{2} \mathbb{W}_{\frac{q}{m}} \mathbb{\Lambda}^1(e^{n+1}+e^n)$
  ! TODO: Here we need to merge at least the steps for the particle mass to improve the scaling
  subroutine advect_e_pic_vm_3d3v_avf( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_helper), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    ! local variables
    sll_int32 :: i_part, i_sp, j
    sll_real64 :: vi(3), xi(3), wi(1), wall(3)
    sll_real64 :: qoverm, factor
    sll_real64 :: efield(3)
    sll_real64 :: rhs(self%n_total1+2*self%n_total0)

    ! Set to zero
    self%j_dofs_local = 0.0_f64
    self%particle_mass_1_local = 0.0_f64
    self%particle_mass_2_local = 0.0_f64
    self%particle_mass_3_local = 0.0_f64

    self%efield_filter = self%efield_dofs
    call self%filter%apply_inplace(self%efield_filter(1:self%n_total1))
    call self%filter%apply_inplace(self%efield_filter(self%n_total1+1:self%n_total1+self%n_total0))
    call self%filter%apply_inplace(self%efield_filter(self%n_total1+self%n_total0+1:self%n_total1+2*self%n_total0))


    ! First particle loop
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       if( self%adiabatic_electrons ) then
          factor = dt**2 * 0.25_f64* qoverm
       else
          factor = dt**2 * 0.25_f64* qoverm*self%betar(2)
       end if
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

          ! Accumulate the particle mass matrix diagonals
          call self%particle_mesh_coupling%add_particle_mass( xi, wi(1) * factor, &
               [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
               self%particle_mass_1_local )
          call self%particle_mesh_coupling%add_particle_mass( xi, wi(1) * factor, &
               [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
               self%particle_mass_2_local )
          call self%particle_mesh_coupling%add_particle_mass( xi, wi(1) * factor, &
               [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
               self%particle_mass_3_local )

          ! Accumulate j
          call self%particle_mesh_coupling%add_charge( xi, wi(1)*vi(1), &
               [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
               self%j_dofs_local(1:self%n_total1) )
          call self%particle_mesh_coupling%add_charge( xi, wi(1)*vi(2), &
               [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
               self%j_dofs_local(1+self%n_total1:self%n_total1+self%n_total0) )
          call self%particle_mesh_coupling%add_charge( xi, wi(1)*vi(3), &
               [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
               self%j_dofs_local(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0) )

          ! Evaulate E at particle position and propagate v a half step       
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
               self%efield_filter(1:self%n_total1), efield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
               self%efield_filter(self%n_total1+1:self%n_total1+self%n_total0), efield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
               self%efield_filter(self%n_total1+self%n_total0+1:self%n_total1+2*self%n_total0), efield(3))

          ! velocity update
          vi = vi + dt* 0.5_f64* qoverm * efield
          call self%particle_group%group(i_sp)%set_v( i_part, vi )
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xi, vi, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if

       end do
    end do
    self%j_dofs = 0.0_f64
    self%particle_mass_1%particle_mass = 0.0_f64
    self%particle_mass_2%particle_mass = 0.0_f64
    self%particle_mass_3%particle_mass = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         self%n_total1+self%n_total0*2, MPI_SUM, self%j_dofs )

    if ( self%jmean ) then
       self%j_dofs(1:self%n_total1) = self%j_dofs(1:self%n_total1) - sum(self%j_dofs(1:self%n_total1))/real(self%n_total1, f64)
       self%j_dofs(1+self%n_total1:self%n_total1+self%n_total0) = self%j_dofs(1+self%n_total1:self%n_total1+self%n_total0) - sum(self%j_dofs(1+self%n_total1:self%n_total1+self%n_total0))/real(self%n_total0, f64)
       self%j_dofs(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0) = self%j_dofs(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0) - sum(self%j_dofs(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0))/real(self%n_total0, f64)
    end if

    call self%filter%apply_inplace(self%j_dofs(1:self%n_total1))
    call self%filter%apply_inplace(self%j_dofs(self%n_total1+1:self%n_total1+self%n_total0))
    call self%filter%apply_inplace(self%j_dofs(self%n_total1+self%n_total0+1:self%n_total1+2*self%n_total0))



    call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_1_local, &
         self%n_total1*self%nspan(1) , MPI_SUM, self%particle_mass_1%particle_mass)
    call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_2_local, &
         self%n_total0*self%nspan(2) , MPI_SUM, self%particle_mass_2%particle_mass)
    call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_3_local, &
         self%n_total0*self%nspan(3) , MPI_SUM, self%particle_mass_3%particle_mass)


    call self%particle_mass_op%set( 1, 1, self%particle_mass_1 )
    call self%particle_mass_op%set( 2, 2, self%particle_mass_2 )
    call self%particle_mass_op%set( 3, 3, self%particle_mass_3 )

    if( self%adiabatic_electrons ) then
       ! Compute rhs
       self%linear_op_schur_phiv%sign = -1._f64
       call self%linear_op_schur_phiv%dot(self%phi_dofs, rhs(1: self%n_total0 ))
       call self%maxwell_solver%multiply_gt( self%j_dofs, rhs(self%n_total0+1: 2*self%n_total0 ) )
       rhs(1: self%n_total0 ) = rhs(1: self%n_total0 ) + dt * rhs(self%n_total0+1: 2*self%n_total0 )

       ! Solve the Schur complement
       self%linear_op_schur_phiv%sign = 1._f64
       call self%linear_solver_schur_phiv%set_guess(self%phi_dofs)
       call self%linear_solver_schur_phiv%solve( rhs(1: self%n_total0 ), self%phi_dofs )
       if( self%boundary ) then !perfect conductor boundary
          do j = 1, self%maxwell_solver%n_dofs(2)*self%maxwell_solver%n_dofs(3) 
             self%phi_dofs(1+(j-1)*self%maxwell_solver%n_dofs(1)) = 0._f64
             self%phi_dofs(j*self%maxwell_solver%n_dofs(1)) = 0._f64
          end do
       end if
       call self%maxwell_solver%multiply_g( self%phi_dofs, self%efield_dofs )
       self%efield_dofs = -self%efield_dofs
    else
       ! Compute rhs
       self%linear_op_schur_ev%sign = -self%force_sign 
       call self%linear_op_schur_ev%dot(self%efield_dofs, rhs)
       rhs = rhs - self%force_sign *dt * self%betar(2) * self%j_dofs

       ! Solve the Schur complement
       self%linear_op_schur_ev%sign = self%force_sign 
       call self%linear_solver_schur_ev%set_guess(self%efield_dofs)
       call self%linear_solver_schur_ev%solve(rhs, self%efield_dofs)

    end if

    if( self%boundary ) then   !perfect conductor boundary
       do j = 1, self%maxwell_solver%n_dofs(2)*self%maxwell_solver%n_dofs(3)
          self%efield_dofs(self%n_total1+1+(j-1)*self%maxwell_solver%n_dofs(1)) = 0._f64
          self%efield_dofs(self%n_total1+j*self%maxwell_solver%n_dofs(1)) = 0._f64
          self%efield_dofs(self%n_total1+self%n_total0+1+(j-1)*self%maxwell_solver%n_dofs(1)) = 0._f64
          self%efield_dofs(self%n_total1+self%n_total0+j*self%maxwell_solver%n_dofs(1)) = 0._f64
       end do
    end if

    self%efield_filter = self%efield_dofs
    call self%filter%apply_inplace(self%efield_filter(1:self%n_total1))
    call self%filter%apply_inplace(self%efield_filter(self%n_total1+1:self%n_total1+self%n_total0))
    call self%filter%apply_inplace(self%efield_filter(self%n_total1+self%n_total0+1:self%n_total1+2*self%n_total0))



    ! Second particle loop (second half step of particle propagation)
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)

          ! Evaluate E at particle position and propagate v a half step
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
               self%efield_filter(1:self%n_total1), efield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
               self%efield_filter(self%n_total1+1:self%n_total1+self%n_total0), efield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
               self%efield_filter(self%n_total1+self%n_total0+1:self%n_total1+2*self%n_total0), efield(3))

          vi = vi + dt* 0.5_f64* qoverm * efield

          call self%particle_group%group(i_sp)%set_v( i_part, vi )
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             ! Update particle weights
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xi, vi, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )

          end if
       end do
    end do

  end subroutine advect_e_pic_vm_3d3v_avf


  !---------------------------------------------------------------------------!
  !> advect_ex: Equations to be solved
  !> $\frac{X^{n+1}-X^n}{\Delta t}=\frac{V^{n+1}+V^n}{2}$
  !> $\frac{V^{n+1}-V^n}{\Delta t}=\mathbb{W}_{\frac{q}{m}} \frac{1}{\Delta t}\int_{t^n}^{t^{n+1}} \mathbb{\Lambda}^1(X(\tau))  d\tau \frac{e^{n+1}+e^n}{2}$
  !> $\frac{M_1 e^{n+1}-M_1 e^n}{\Delta t} = - \frac{1}{\Delta t} \int_{t^n}^{t^{n+1}} \mathbb{\Lambda}^1(X(\tau))^\top  d\tau \mathbb{W}_q\frac{V^{n+1}+V^n}{2}$
  ! TODO: Here we need to merge at least the steps for the particle mass to improve the scaling
  subroutine advect_ex( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_helper), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    ! local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: vi(3), xi(3), wi(1), vnew(3), xnew(3), vbar(3),  xmid(3), xbar, dx
    sll_real64 :: qoverm, wall(3) 
    sll_real64 :: efield(3)
    sll_int32 :: niter
    sll_real64 :: residual(1), residual_local(1)
   
    self%efield_dofs_new = self%efield_dofs
    self%phi_dofs_new = self%phi_dofs
    do i_sp=1,self%particle_group%n_species
       do i_part = 1,self%particle_group%group(i_sp)%n_particles
          self%vnew(i_sp,:, i_part) = self%particle_group%group(i_sp)%get_v(i_part)
          self%xnew(i_sp,:, i_part) = self%particle_group%group(i_sp)%get_x(i_part)
       end do
    end do

    niter = 0
    residual = self%iter_tolerance + 1.0_f64
    do while ( (residual(1) > self%iter_tolerance) .and. niter < self%max_iter )
       niter = niter+1

       self%efield_dofs_work = 0.5_f64*( self%efield_dofs + self%efield_dofs_new )
       self%j_dofs_local = 0.0_f64

       ! Particle loop
       do i_sp=1,self%particle_group%n_species
          qoverm = self%particle_group%group(i_sp)%species%q_over_m();
          do i_part = 1,self%particle_group%group(i_sp)%n_particles
             vi = self%particle_group%group(i_sp)%get_v(i_part)
             xi = self%particle_group%group(i_sp)%get_x(i_part)

             ! Get charge for accumulation of j
             wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

             vnew = self%vnew(i_sp,:, i_part)
             vbar = 0.5_f64*(vi+vnew)

             xnew = xi + dt * vbar
             if(xnew(1) < -self%x_max(1) .or.  xnew(1) > 2._f64*self%x_max(1) ) then
                print*, xnew
                SLL_ERROR("particle boundary", "particle out of bound")
             else if(xnew(1) < self%x_min(1) .or. xnew(1) > self%x_max(1) )then
                if(xnew(1) < self%x_min(1)  )then
                   xbar = self%x_min(1)
                else if(xnew(1) > self%x_max(1))then
                   xbar = self%x_max(1)
                end if
                dx = (xbar- xi(1))/(xnew(1)-xi(1))
                xmid = xi + dx * (xnew-xi)
                xmid(1) = xbar

                vbar = (xmid - xi)*wi(1)
                call self%particle_mesh_coupling%add_current_evaluate( xi, xmid, vbar, self%efield_dofs_work, &
                     self%j_dofs_local, efield )
                vnew = vi + qoverm * dt* dx* efield
                select case(self%boundary_particles)
                case(sll_p_boundary_particles_reflection)
                   xnew(1) = 2._f64*xbar-xnew(1)
                   vnew(1) = - vnew(1)
                case(sll_p_boundary_particles_absorption)
                   call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
                   xnew(1) = xmid(1) + ((xbar-self%x_min(1))/self%Lx(1)-0.5_f64) * 1.9_f64*self%iter_tolerance 
                case( sll_p_boundary_particles_periodic)
                   xnew(1) = self%x_min(1) + modulo(xnew(1)-self%x_min(1), self%Lx(1))
                   xmid(1) = self%x_max(1)+self%x_min(1)-xbar
                case default
                   xnew(1) = self%x_min(1) + modulo(xnew(1)-self%x_min(1), self%Lx(1))
                   xmid(1) = self%x_max(1)+self%x_min(1)-xbar
                end select
                if( abs(xnew(1)-xmid(1)) > self%iter_tolerance ) then
                   vbar = (xnew - xmid)*wi(1)
                   call self%particle_mesh_coupling%add_current_evaluate( xmid, xnew, vbar, self%efield_dofs_work, &
                        self%j_dofs_local, efield )
                   vnew = vnew + qoverm * (1._f64-dx)*dt * efield
                   if(self%boundary_particles == sll_p_boundary_particles_reflection) then
                      vnew(1) = - vnew(1)
                   end if
                else
                   xnew(1) = xmid(1) 
                end if
             else
                vbar = (xnew - xi)*wi(1)

                call self%particle_mesh_coupling%add_current_evaluate( xi, xnew, vbar, self%efield_dofs_work, &
                     self%j_dofs_local, efield )
                vnew = vi + qoverm * dt * efield
             end if

             self%xnew(i_sp, :, i_part) = xnew
             self%vnew(i_sp, :, i_part) = vnew

          end do
       end do

       self%j_dofs = 0.0_f64
       ! MPI to sum up contributions from each processor
       call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
            self%n_total1+self%n_total0*2, MPI_SUM, self%j_dofs )

       if( self%adiabatic_electrons) then
          self%phi_dofs_work = self%phi_dofs
          call self%maxwell_solver%compute_phi_from_j( self%j_dofs, self%phi_dofs_work, self%efield_dofs_work )

          ! Compute residual based on phi
          residual_local = (sum((self%phi_dofs_work-self%phi_dofs_new)**2))*product(self%particle_mesh_coupling%delta_x)
          call sll_o_collective_allreduce( sll_v_world_collective, residual_local, 1, MPI_MAX, residual )
          residual = sqrt(residual)
       else
          self%efield_dofs_work = self%efield_dofs
          self%j_dofs=self%j_dofs*self%betar(2)
          call self%maxwell_solver%compute_E_from_j( self%j_dofs, self%efield_dofs_work )

          ! Compute residual based on e
          residual_local = (sum((self%efield_dofs_work-self%efield_dofs_new)**2))*product(self%particle_mesh_coupling%delta_x)
          call sll_o_collective_allreduce( sll_v_world_collective, residual_local, 1, MPI_MAX, residual )
          residual = sqrt(residual)
       end if
       self%phi_dofs_new = self%phi_dofs_work
       self%efield_dofs_new = self%efield_dofs_work
    end do

    if ( residual(1) > self%iter_tolerance ) then
       print*, 'Warning: Iteration no.', self%iter_counter+1 ,'did not converge.', residual, niter
       self%n_failed = self%n_failed+1
    end if
    
    
    self%efield_dofs_work = 0.5_f64*( self%efield_dofs + self%efield_dofs_new )
    self%j_dofs_local = 0.0_f64

    ! Particle loop
    do i_sp=1,self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part = 1,self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

          vnew = self%vnew(i_sp,:, i_part)
          vbar = 0.5_f64*(vi+vnew)

          xnew = xi + dt * vbar
          if(xnew(1) < -self%x_max(1) .or.  xnew(1) > 2._f64*self%x_max(1) ) then
             print*, xnew
             SLL_ERROR("particle boundary", "particle out of bound")
          else if(xnew(1) < self%x_min(1) .or. xnew(1) > self%x_max(1) )then
             if(xnew(1) < self%x_min(1)  )then
                xbar = self%x_min(1)
             else if(xnew(1) > self%x_max(1))then
                xbar = self%x_max(1)
             end if
             dx = (xbar- xi(1))/(xnew(1)-xi(1))
             xmid = xi + dx * (xnew-xi)
             xmid(1) = xbar

             vbar = (xmid - xi)*wi(1)
             call self%particle_mesh_coupling%add_current_evaluate( xi, xmid, vbar, self%efield_dofs_work, &
                  self%j_dofs_local, efield )
             vnew = vi + qoverm * dt* dx* efield
             select case(self%boundary_particles)
             case(sll_p_boundary_particles_reflection)
                xnew(1) = 2._f64*xbar-xnew(1)
                vnew(1) = - vnew(1)
             case(sll_p_boundary_particles_absorption)
                call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
             case( sll_p_boundary_particles_periodic)
                xnew(1) = self%x_min(1) + modulo(xnew(1)-self%x_min(1), self%Lx(1))
                xmid(1) = self%x_max(1)+self%x_min(1)-xbar
             case default
                xnew(1) = self%x_min(1) + modulo(xnew(1)-self%x_min(1), self%Lx(1))
                xmid(1) = self%x_max(1)+self%x_min(1)-xbar
             end select
             if( abs(xnew(1)-xmid(1)) > self%iter_tolerance ) then
                vbar = (xnew - xmid)*wi(1)
                call self%particle_mesh_coupling%add_current_evaluate( xmid, xnew, vbar, self%efield_dofs_work, &
                     self%j_dofs_local, efield )
                vnew = vnew + qoverm * (1._f64-dx)*dt * efield
             else
                xnew(1) = xmid(1) 
             end if
          else
             vbar = (xnew - xi)*wi(1)

             call self%particle_mesh_coupling%add_current_evaluate( xi, xnew, vbar, self%efield_dofs_work, &
                  self%j_dofs_local, efield )
             vnew = vi + qoverm * dt * efield
          end if

          xnew(2:3) = self%x_min(2:3) + modulo(xnew(2:3)-self%x_min(2:3), self%Lx(2:3))
          call self%particle_group%group(i_sp)%set_v( i_part, vnew )
          call self%particle_group%group(i_sp)%set_x( i_part, xnew )
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xnew, vnew, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do
    end do

    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         self%n_total1+self%n_total0*2, MPI_SUM, self%j_dofs )

    if( self%adiabatic_electrons) then
       call self%maxwell_solver%compute_phi_from_j( self%j_dofs, self%phi_dofs, self%efield_dofs )
    else
       self%j_dofs=self%j_dofs*self%betar(2)
       call self%maxwell_solver%compute_E_from_j( self%j_dofs, self%efield_dofs )

    end if
    
    self%iter_counter = self%iter_counter + 1
    self%niter(self%iter_counter) = niter

  end subroutine advect_ex


  !---------------------------------------------------------------------------!
  !> advect_e_trunc: This is a version of advect_e where the Particle Mass is
  !> approximated by the mass matrix
  subroutine advect_e_pic_vm_3d3v_avf_trunc( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_helper), intent(inout) :: self !< time propagator object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    ! local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: vi(3), xi(3), wi(1)
    sll_real64 :: qoverm
    sll_real64 :: efield(3)

    ! Set to zero
    self%j_dofs_local = 0.0_f64
    self%particle_mass_1_local = 0.0_f64
    self%particle_mass_2_local = 0.0_f64
    self%particle_mass_3_local = 0.0_f64


    ! First particle loop
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part = 1,self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

          ! Accumulate j
          call self%particle_mesh_coupling%add_charge( xi, wi(1)*vi(1), &
               [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
               self%j_dofs_local(1:self%n_total1) )
          call self%particle_mesh_coupling%add_charge( xi, wi(1)*vi(2), &
               [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
               self%j_dofs_local(1+self%n_total1:self%n_total1+self%n_total0) )
          call self%particle_mesh_coupling%add_charge( xi, wi(1)*vi(3), &
               [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
               self%j_dofs_local(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0) )


          ! Evaulate E at particle position and propagate v a half step       
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], self%efield_dofs(1:self%n_total1), efield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)],self%efield_dofs(self%n_total1+1:self%n_total1+self%n_total0), efield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1],self%efield_dofs(self%n_total1+self%n_total0+1:self%n_total1+2*self%n_total0), efield(3))

          vi = vi + dt* 0.5_f64* qoverm * efield
          call self%particle_group%group(i_sp)%set_v( i_part, vi )
       end do

    end do
    ! Update d_n
    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         self%n_total1+self%n_total0*2, MPI_SUM, self%j_dofs )


    ! Ex part
    ! Multiply E by mass matrix and add current
    call self%maxwell_solver%multiply_mass( [2,1,1], &
         self%efield_dofs(1:self%n_total1 ), self%j_dofs_local(1:self%n_total1) )
    self%j_dofs_local( 1:self%n_total1 ) = (1.0_f64-dt**2*0.25_f64*qoverm)*&
         self%j_dofs_local( 1:self%n_total1 ) - dt* self%j_dofs( 1:self%n_total1 )

    ! Ey part
    ! Multiply E by mass matrix and add current
    call self%maxwell_solver%multiply_mass( [1,2,1], &
         self%efield_dofs(1+self%n_total1:self%n_total1+self%n_total0 ), self%j_dofs_local(1+self%n_total1:self%n_total1+self%n_total0) )
    self%j_dofs_local( 1+self%n_total1:self%n_total1+self%n_total0 ) = (1.0_f64-dt**2*0.25_f64*qoverm)*&
         self%j_dofs_local( 1+self%n_total1:self%n_total1+self%n_total0 ) - dt* self%j_dofs( 1+self%n_total1:self%n_total1+self%n_total0 )

    ! Ez part
    ! Multiply E by mass matrix and add current
    call self%maxwell_solver%multiply_mass( [1,1,2], self%efield_dofs(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0 ), &
         self%j_dofs_local(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0) )
    self%j_dofs_local( 1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0 ) = (1.0_f64-dt**2*0.25_f64*qoverm)*&
         self%j_dofs_local( 1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0 ) - dt* self%j_dofs( 1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0 )

    ! Solve Schur complement
    ! Solve efield components together
    call self%maxwell_solver%multiply_mass_inverse( 1, self%j_dofs_local, self%efield_dofs )

    ! Second particle loop (second half step of particle propagation)
    do i_sp = 1, self%particle_group%n_species
       do i_part = 1, self%particle_group%group(i_sp)%n_particles

          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)

          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], self%efield_dofs(1:self%n_total1), efield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)],self%efield_dofs(self%n_total1+1:self%n_total1+self%n_total0), efield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1],self%efield_dofs(self%n_total1+self%n_total0+1:self%n_total1+2*self%n_total0), efield(3))

          vi = vi + dt* 0.5_f64* qoverm * efield
          call self%particle_group%group(i_sp)%set_v( i_part, vi )

       end do
    end do

  end subroutine advect_e_pic_vm_3d3v_avf_trunc


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
       filter, &
       boundary_particles, &
       solver_tolerance, &
       iter_tolerance, max_iter, &
       betar, &
       force_sign, &
       rhob, &
       control_variate, &
       jmean, &
       lindf)  
    class(sll_t_time_propagator_pic_vm_3d3v_helper), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,          intent( in ) :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target,   intent(in)  :: particle_mesh_coupling !< Kernel smoother
    class(sll_t_particle_array), target,           intent( in ) :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:) !< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent( in ) :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent( in ) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                    intent( in ) :: x_min(3) !< Lower bound of x domain
    sll_real64,                                    intent( in ) :: Lx(3) !< Length of the domain in x direction.
    class(sll_c_filter_base_3d), target :: filter !< filter
    sll_int32, optional,                           intent( in ) :: boundary_particles !< particle boundary conditions
    sll_real64, optional,                          intent( in ) :: solver_tolerance !< Solver tolerance
    sll_real64, optional,                          intent( in ) :: iter_tolerance !< iteration tolerance
    sll_int32,  optional,                          intent( in ) :: max_iter !< maximal number of iterations
    sll_real64, optional, intent(in) :: betar(2) !< reciprocal plasma beta
    sll_real64, optional, intent(in) :: force_sign !< sign of particle force
    sll_real64, optional, target, intent( in ) :: rhob(:) !< charge at the boundary
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    logical, optional, intent(in) :: jmean !< logical for mean value of current
    logical, optional, intent(in) :: lindf !< true for linear delta f method
    !local variables
    sll_int32 :: ierr, i_sp, j
    sll_real64, allocatable :: vec1(:)
    sll_real64, allocatable :: vec2(:)

    if (present(solver_tolerance) )  then
       self%solver_tolerance = solver_tolerance
    else
       self%solver_tolerance = 1d-12
    end if

    if (present(iter_tolerance) )  then
       self%iter_tolerance = iter_tolerance
       self%max_iter = max_iter
    else
       self%iter_tolerance = 1d-10
       self%max_iter = 20
    end if

    if( present(force_sign) )then
       self%force_sign = force_sign
    end if

    if (self%force_sign == 1._f64) then
       if( particle_group%group(1)%species%q > 0._f64) self%adiabatic_electrons = .true.
    end if

    if (present(jmean)) then
       self%jmean = jmean
    end if

    if (present(boundary_particles) ) then
       self%boundary_particles = boundary_particles
    end if

    if( present(rhob) )then
       self%rhob => rhob
    end if

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
    self%filter => filter

    self%nspan(1) =  (2*self%spline_degree(1)-1)*(2*self%spline_degree(2)+1)*(2*self%spline_degree(3)+1)
    self%nspan(2) =  (2*self%spline_degree(1)+1)*(2*self%spline_degree(2)-1)*(2*self%spline_degree(3)+1)
    self%nspan(3) =  (2*self%spline_degree(1)+1)*(2*self%spline_degree(2)+1)*(2*self%spline_degree(3)-1)

    SLL_ALLOCATE( vec1(1:self%n_total1+self%n_total0*2), ierr )
    SLL_ALLOCATE( vec2(1:self%n_total1+self%n_total0*2), ierr )
    vec1 = 0._f64
    vec2 = 0._f64

    SLL_ALLOCATE( self%j_dofs(1:self%n_total1+self%n_total0*2), ierr )
    SLL_ALLOCATE( self%j_dofs_local(1:self%n_total1+self%n_total0*2), ierr )

    SLL_ALLOCATE(self%particle_mass_1_local(self%nspan(1), self%n_total0), ierr )
    SLL_ALLOCATE(self%particle_mass_2_local(self%nspan(2), self%n_total0), ierr )
    SLL_ALLOCATE(self%particle_mass_3_local(self%nspan(3), self%n_total0), ierr )
    self%j_dofs = 0.0_f64
    self%j_dofs_local = 0.0_f64
    self%particle_mass_1_local = 0.0_f64
    self%particle_mass_2_local = 0.0_f64
    self%particle_mass_3_local = 0.0_f64

    if( self%particle_mesh_coupling%n_cells(1)+self%spline_degree(1) == self%maxwell_solver%n_dofs(1)   ) then
       self%boundary = .true.
       allocate( sll_t_linear_operator_particle_mass_cl_3d_diag  :: self%particle_mass_1 )
       allocate( sll_t_linear_operator_particle_mass_cl_3d_diag  :: self%particle_mass_2 )
       allocate( sll_t_linear_operator_particle_mass_cl_3d_diag  :: self%particle_mass_3 )
    else
       allocate( sll_t_linear_operator_particle_mass_3d_diag  :: self%particle_mass_1 )
       allocate( sll_t_linear_operator_particle_mass_3d_diag  :: self%particle_mass_2 )
       allocate( sll_t_linear_operator_particle_mass_3d_diag  :: self%particle_mass_3 )
    end if

    call self%particle_mass_1%create( self%particle_mesh_coupling%n_cells, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)] )
    call self%particle_mass_2%create( self%particle_mesh_coupling%n_cells, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)] )
    call self%particle_mass_3%create( self%particle_mesh_coupling%n_cells, [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1]  )

    call self%particle_mass_op%create( 3, 3 )
    if( self%adiabatic_electrons ) then
       call self%linear_op_schur_phiv%create( self%maxwell_solver, self%particle_mass_op )
       call self%linear_solver_schur_phiv%create( self%linear_op_schur_phiv )
       self%linear_solver_schur_phiv%atol = self%solver_tolerance
       !self%linear_solver_schur_phiv%verbose = .true.
    else
       call self%linear_op_schur_ev%create( self%maxwell_solver, self%particle_mass_op)
       call self%preconditioner_fft%init( self%maxwell_solver%Lx, self%maxwell_solver%n_cells, self%maxwell_solver%s_deg_0, self%boundary )
       if( self%boundary ) then
          vec1 = 1._f64
          call self%maxwell_solver%multiply_mass( [1], vec1, vec2 )
          do j = 1, self%maxwell_solver%n_dofs(2)*self%maxwell_solver%n_dofs(3)
             vec2(self%n_total1+1+(j-1)*self%maxwell_solver%n_dofs(1)) = 1._f64
             vec2(self%n_total1+j*self%maxwell_solver%n_dofs(1)) = 1._f64
             vec2(self%n_total1+self%n_total0+1+(j-1)*self%maxwell_solver%n_dofs(1)) = 1._f64
             vec2(self%n_total1+self%n_total0+j*self%maxwell_solver%n_dofs(1)) = 1._f64
          end do
          vec2 = 1._f64/sqrt(abs(vec2))

          call self%preconditioner1%create( self%preconditioner_fft%inverse_mass1_3d, vec2, 2*self%n_total0+self%n_total1 )

          call self%linear_solver_schur_ev%create( self%linear_op_schur_ev, self%preconditioner1)
       else
          call self%linear_solver_schur_ev%create( self%linear_op_schur_ev )
       end if
       self%linear_solver_schur_ev%atol = self%solver_tolerance
       !self%linear_solver_schur_ev%verbose = .true.
       self%linear_solver_schur_ev%n_maxiter = 2000

    end if

    if (present(betar)) then
       self%betar = betar
    else
       self%betar = 1.0_f64
    end if


    self%n_particles_max = self%particle_group%group(1)%n_particles
    do i_sp = 2,self%particle_group%n_species       
       self%n_particles_max = max(self%n_particles_max, self%particle_group%group(i_sp)%n_particles )
    end do

    SLL_ALLOCATE( self%xnew(self%particle_group%n_species,3,self%n_particles_max), ierr )
    SLL_ALLOCATE( self%vnew(self%particle_group%n_species,3,self%n_particles_max), ierr )
    SLL_ALLOCATE( self%efield_dofs_new(self%n_total1+2*self%n_total0), ierr )
    SLL_ALLOCATE( self%efield_dofs_work(self%n_total1+2*self%n_total0), ierr )
    SLL_ALLOCATE( self%phi_dofs_work(self%n_total0), ierr )
    SLL_ALLOCATE( self%phi_dofs_new(self%n_total0), ierr )
    SLL_ALLOCATE(self%efield_filter(self%n_total1+2*self%n_total0), ierr)
    SLL_ALLOCATE(self%bfield_filter(self%n_total0+2*self%n_total1), ierr)

    if (present(control_variate)) then
       allocate(self%control_variate )
       allocate(self%control_variate%cv(self%particle_group%n_species) )
       self%control_variate => control_variate
       self%i_weight = 3
       if (present(lindf)) then
          self%lindf = lindf
       end if
    end if

  end subroutine initialize_pic_vm_3d3v



  !> Constructor.
  subroutine initialize_file_pic_vm_3d3v(&
       self, &
       maxwell_solver, &
       particle_mesh_coupling, &
       particle_group, &
       phi_dofs, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       filename, &
       filter, &
       boundary_particles, &
       betar, &
       force_sign, &
       rhob, &
       control_variate, &
       jmean, &
       lindf)  
    class(sll_t_time_propagator_pic_vm_3d3v_helper), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,          intent( in ) :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target, intent(in) :: particle_mesh_coupling !< Kernel smoother
    class(sll_t_particle_array), target,           intent( in ) :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:)!< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent( in ) :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent( in ) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                    intent( in ) :: x_min(3) !< Lower bound of x domain
    sll_real64,                                    intent( in ) :: Lx(3) !< Length of the domain in x direction.
    character(len=*),                              intent( in ) :: filename !< filename
    class(sll_c_filter_base_3d), target, optional :: filter !< filter
    sll_int32, optional,                           intent( in ) :: boundary_particles !< particle boundary conditions
    sll_real64, optional,                          intent( in ) :: betar(2) !< reciprocal plasma beta
    sll_real64, optional, intent(in) :: force_sign !< sign of particle force
    sll_real64, optional, target,                  intent( in ) :: rhob(:) !< charge at the boundary
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    logical, optional, intent(in) :: jmean !< logical for mean value of current
    logical, optional, intent(in) :: lindf !< true for linear delta f method
    !local variables
    character(len=256) :: file_prefix
    sll_int32 :: input_file, rank
    sll_int32 :: io_stat, io_stat0, io_stat1, file_id
    sll_real64 :: maxwell_tolerance, iter_tolerance, betar_set(2), force_sign_set
    sll_int32 :: max_iter, boundary_particles_set
    logical :: jmean_set
    logical :: lindf_set

    namelist /output/ file_prefix
    namelist /time_solver/ maxwell_tolerance
    namelist /time_iterate/ iter_tolerance, max_iter

    rank = sll_f_get_collective_rank(sll_v_world_collective)

    if( present(boundary_particles) )then
       boundary_particles_set = boundary_particles
    else
       boundary_particles_set = 100
    end if

    if( present(betar) )then
       betar_set = betar
    else
       betar_set = 1._f64
    end if

    if( present(force_sign) )then
       force_sign_set = force_sign
    else
       force_sign_set = 1._f64
    end if

    if (present(jmean)) then
       jmean_set = jmean
    else
       jmean_set = .false.
    end if

    if (present(lindf)) then
       lindf_set = lindf
    else
       lindf_set = .false.
    end if

    if( present(control_variate) ) then
       ! Read in solver tolerance
       open(newunit = input_file, file=filename, status='old',IOStat=io_stat)
       if (io_stat /= 0) then
          if (rank == 0 ) then
             print*, 'sll_m_time_propagator_pic_vm_3d3v_helper: Input file does not exist. Set default tolerance.'
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'field solver tolerance:', 1d-12
             write(file_id, *) 'iter_tolerance:', 1d-10
             write(file_id, *) 'max_iter:', 20
             write(file_id, *) 'control_variate:', .true.
             close(file_id) 
          end if
          call self%init( maxwell_solver, &
               particle_mesh_coupling, &
               particle_group, &
               phi_dofs, &
               efield_dofs, &
               bfield_dofs, &
               x_min, &
               Lx, &
               filter, &
               boundary_particles = boundary_particles_set, &
               betar=betar_set, &
               force_sign=force_sign_set, &
               rhob = rhob, &
               control_variate = control_variate,&
               jmean=jmean_set, &
               lindf = lindf_set)
       else
          read(input_file, output,IOStat=io_stat0)
          read(input_file, time_solver,IOStat=io_stat)
          read(input_file, time_iterate,IOStat=io_stat1)
          if (io_stat /= 0.and. io_stat1 /= 0) then
             if (rank == 0 ) then
                print*, 'sll_m_time_propagator_pic_vm_3d3v_helper: Input parameter does not exist. Set default tolerance.'
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'field solver tolerance:', 1d-12
                write(file_id, *) 'iter_tolerance:', 1d-10
                write(file_id, *) 'max_iter:', 20
                write(file_id, *) 'control_variate:', .true.
                close(file_id) 
             end if
             call self%init( maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  filter, &
                  boundary_particles = boundary_particles_set, &
                  betar=betar_set, &
                  force_sign=force_sign_set, &
                  rhob = rhob, &
                  control_variate = control_variate,&
                  jmean=jmean_set)
          else if (io_stat == 0 .and. io_stat1 /= 0) then
             if (rank == 0 ) then
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'field solver tolerance:', maxwell_tolerance
                write(file_id, *) 'iter_tolerance:', 1d-10
                write(file_id, *) 'max_iter:', 20
                write(file_id, *) 'control_variate:', .true.
                close(file_id) 
             end if
             call self%init(   maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  filter, &
                  boundary_particles_set, &
                  maxwell_tolerance, &
                  betar=betar_set, &
                  force_sign=force_sign_set, &
                  rhob = rhob, &
                  control_variate = control_variate,&
                  jmean=jmean_set, &
                  lindf = lindf_set)
          else if (io_stat /= 0 .and. io_stat1 == 0) then
             if (rank == 0 ) then
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'field solver tolerance:', 1d-12
                write(file_id, *) 'iter_tolerance:', iter_tolerance
                write(file_id, *) 'max_iter:', max_iter
                write(file_id, *) 'control_variate:', .true.
                close(file_id) 
             end if
             call self%init(   maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  filter, &
                  boundary_particles = boundary_particles_set, &
                  iter_tolerance=iter_tolerance, max_iter=max_iter, &
                  betar=betar_set, &
                  force_sign=force_sign_set, &
                  rhob = rhob, &
                  control_variate = control_variate,&
                  jmean=jmean_set, &
                  lindf = lindf_set)
          else
             if (rank == 0 ) then
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'field solver tolerance:', maxwell_tolerance
                write(file_id, *) 'iter_tolerance:', iter_tolerance
                write(file_id, *) 'max_iter:', max_iter
                write(file_id, *) 'control_variate:', .true.
                close(file_id) 
             end if
             call self%init(   maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  filter, &
                  boundary_particles_set, &
                  maxwell_tolerance, &
                  iter_tolerance, max_iter, &
                  betar_set, &
                  force_sign=force_sign_set, &
                  rhob = rhob, &
                  control_variate = control_variate,&
                  jmean=jmean_set, &
                  lindf = lindf_set)
          end if
          close(input_file)
       end if

    else
       ! Read in solver tolerance
       open(newunit = input_file, file=filename, status='old',IOStat=io_stat)
       if (io_stat /= 0) then
          if (rank == 0 ) then
             print*, 'sll_m_time_propagator_pic_vm_3d3v_helper: Input file does not exist. Set default tolerance.'
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'field solver tolerance:', 1d-12
             write(file_id, *) 'iter_tolerance:', 1d-10
             write(file_id, *) 'max_iter:', 20
             close(file_id) 
          end if
          call self%init( maxwell_solver, &
               particle_mesh_coupling, &
               particle_group, &
               phi_dofs, &
               efield_dofs, &
               bfield_dofs, &
               x_min, &
               Lx, &
               filter, &
               boundary_particles = boundary_particles_set, &
               betar=betar_set, &
               force_sign=force_sign_set,&
               rhob = rhob, &
               jmean=jmean_set)
       else
          read(input_file, output,IOStat=io_stat0)
          read(input_file, time_solver,IOStat=io_stat)
          read(input_file, time_iterate,IOStat=io_stat1)
          if (io_stat /= 0.and. io_stat1 /= 0) then
             if (rank == 0 ) then
                print*, 'sll_m_time_propagator_pic_vm_3d3v_helper: Input parameter does not exist. Set default tolerance.'
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'field solver tolerance:', 1d-12
                write(file_id, *) 'iter_tolerance:', 1d-10
                write(file_id, *) 'max_iter:', 20
                close(file_id) 
             end if
             call self%init( maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  filter, &
                  boundary_particles = boundary_particles_set, &
                  betar=betar_set, &
                  force_sign=force_sign_set,&
                  rhob = rhob, &
                  jmean=jmean_set)
          else if (io_stat == 0 .and. io_stat1 /= 0) then
             if (rank == 0 ) then
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'field solver tolerance:', maxwell_tolerance
                write(file_id, *) 'iter_tolerance:', 1d-10
                write(file_id, *) 'max_iter:', 20
                close(file_id) 
             end if
             call self%init(   maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  filter, &
                  boundary_particles_set, &
                  maxwell_tolerance, &
                  betar=betar_set, &
                  force_sign=force_sign_set, &
                  rhob = rhob)
          else if (io_stat /= 0 .and. io_stat1 == 0) then
             if (rank == 0 ) then
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'field solver tolerance:', 1d-12
                write(file_id, *) 'iter_tolerance:', iter_tolerance
                write(file_id, *) 'max_iter:', max_iter
                close(file_id) 
             end if
             call self%init(   maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  filter, &
                  boundary_particles = boundary_particles_set, &
                  iter_tolerance=iter_tolerance, max_iter=max_iter, &
                  betar=betar_set, &
                  force_sign=force_sign_set,&
                  rhob = rhob, &
                  jmean=jmean_set)
          else
             if (rank == 0 ) then
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'solver tolerance:', maxwell_tolerance
                write(file_id, *) 'iter_tolerance:', iter_tolerance
                write(file_id, *) 'max_iter:', max_iter
                close(file_id) 
             end if
             call self%init(   maxwell_solver, &
                  particle_mesh_coupling, &
                  particle_group, &
                  phi_dofs, &
                  efield_dofs, &
                  bfield_dofs, &
                  x_min, &
                  Lx, &
                  filter, &
                  boundary_particles_set, &
                  maxwell_tolerance, &
                  iter_tolerance, max_iter, &
                  betar_set, &
                  force_sign=force_sign_set,&
                  rhob = rhob, &
                  jmean=jmean_set)
          end if
          close(input_file)
       end if
    end if

  end subroutine initialize_file_pic_vm_3d3v


  !> Destructor.
  subroutine delete_pic_vm_3d3v(self)
    class(sll_t_time_propagator_pic_vm_3d3v_helper), intent( inout ) :: self !< time propagator object
    !local variables
    sll_int32 :: file_id, j


    if(self%adiabatic_electrons) then
       call self%linear_solver_schur_phiv%free()
       call self%linear_op_schur_phiv%free()
    else
       call self%linear_solver_schur_ev%free()
       call self%linear_op_schur_ev%free()
    end if

    call self%particle_mass_op%free()
    call self%particle_mass_1%free()
    call self%particle_mass_2%free()
    call self%particle_mass_3%free()
    deallocate( self%particle_mass_1 )
    deallocate( self%particle_mass_2 )
    deallocate( self%particle_mass_3 )

    deallocate(self%j_dofs)
    deallocate(self%j_dofs_local)
    deallocate( self%particle_mass_1_local )
    deallocate( self%particle_mass_2_local )
    deallocate( self%particle_mass_3_local )
    self%maxwell_solver => null()
    self%particle_mesh_coupling => null()
    self%particle_group => null()
    self%phi_dofs => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()

    deallocate( self%vnew )
    deallocate( self%xnew )
    deallocate( self%efield_dofs_new )
    deallocate( self%efield_dofs_work )
    deallocate( self%phi_dofs_new )
    deallocate( self%phi_dofs_work )


    print*, 'No. of failed iterations:', self%n_failed

    open(newunit=file_id, file="outer-iters.dat", action='write')
    do j=1,self%iter_counter
       write(file_id,*) self%niter(j)
    end do

  end subroutine delete_pic_vm_3d3v


end module sll_m_time_propagator_pic_vm_3d3v_helper
