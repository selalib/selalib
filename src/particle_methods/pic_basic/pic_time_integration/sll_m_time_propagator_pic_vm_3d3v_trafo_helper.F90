!> @ingroup pic_time_integration
!> @author: Benedikt Perse, IPP
module sll_m_time_propagator_pic_vm_3d3v_trafo_helper
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

  use sll_m_time_propagator_pic_vm_3d3v_cl_helper, only: &
       sll_p_boundary_particles_periodic, &
       sll_p_boundary_particles_singular, &
       sll_p_boundary_particles_reflection, &
       sll_p_boundary_particles_absorption, &
       sll_s_compute_particle_boundary_simple, &
       sll_s_compute_matrix_inverse

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

  use sll_m_linear_operator_particle_mass_3d_od, only : &
       sll_t_linear_operator_particle_mass_3d_od

  use sll_m_linear_operator_particle_mass_cl_3d_od, only : &
       sll_t_linear_operator_particle_mass_cl_3d_od

  use sll_m_linear_operator_schur_ev_3d, only : &
       sll_t_linear_operator_schur_ev_3d

  use sll_m_linear_operator_schur_phiv_3d, only : &
       sll_t_linear_operator_schur_phiv_3d

  use sll_m_linear_solver_cg, only : &
       sll_t_linear_solver_cg

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

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
       sll_t_time_propagator_pic_vm_3d3v_trafo_helper

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Helper for implicit time propagator for 3d3v Vlasov-Maxwell with coordinate transformation 
  type :: sll_t_time_propagator_pic_vm_3d3v_trafo_helper
     class(sll_c_maxwell_3d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_particle_mesh_coupling_3d), pointer :: particle_mesh_coupling !< Particle mesh coupling
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_1 !< Particle mass
     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_2 !< Particle mass
     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_3 !< Particle mass
     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_12 !< Particle mass
     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_21 !< Particle mass
     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_13 !< Particle mass
     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_31 !< Particle mass
     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_23 !< Particle mass
     class(sll_c_particle_mass_3d_base), allocatable :: particle_mass_32 !< Particle mass

     type( sll_t_linear_operator_block )    :: particle_mass_op !< Particle mass operator
     type( sll_t_linear_operator_schur_ev_3d ) :: linear_op_schur_ev !< Schur complement operator for advect_e
     type( sll_t_linear_operator_schur_phiv_3d ) :: linear_op_schur_phiv !< Schur complement operator for advect_e
     type(sll_t_preconditioner_fft)         :: preconditioner_fft !< Preconditioner for Schur complement solver
     type(sll_t_preconditioner_singular) :: preconditioner1 !< preconditioner for mass matrices
     type( sll_t_linear_solver_cg )         :: linear_solver_schur_ev !< Schur complement solver for advect_e
     type( sll_t_linear_solver_cg )         :: linear_solver_schur_phiv !< Schur complement solver for advect_e

     type( sll_t_mapping_3d ), pointer      :: map !< coordinate transformation

     sll_int32 :: spline_degree(3) !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx(3) !< Size of the domain
     sll_real64 :: x_min(3) !< Lower bound for x domain
     sll_real64 :: x_max(3) !< Upper bound for x domain
     sll_int32 :: n_total0 !< total number of Dofs for 0form
     sll_int32 :: n_total1 !< total number of Dofs for 1form
     sll_int32 :: nspan_1(3), nspan_2(3) !< Number of intervals where spline tensor product is non zero 

     sll_real64 :: betar(2) !< reciprocal of plasma beta

     sll_real64, allocatable :: vec1(:) !< scratch data
     sll_real64, allocatable :: vec2(:) !< scratch data

     sll_real64, pointer     :: phi_dofs(:) !< DoFs describing the scalar potential
     sll_real64, pointer     :: efield_dofs(:) !< DoFs describing the three components of the electric field
     sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the three components of the magnetic field
     sll_real64, allocatable :: j_dofs(:)      !< DoFs for representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:) !< MPI-processor local part of one component of \a j_dofs
     sll_real64, allocatable :: particle_mass_1_local(:,:) !< MPI-processor local part of one component of \a particle_mass
     sll_real64, allocatable :: particle_mass_2_local(:,:) !< MPI-processor local part of one component of \a particle_mass
     sll_real64, allocatable :: particle_mass_3_local(:,:) !< MPI-processor local part of one component of \a particle_mass
     sll_real64, allocatable :: particle_mass_12_local(:,:) !< MPI-processor local part of one component of \a particle_mass
     sll_real64, allocatable :: particle_mass_21_local(:,:) !< MPI-processor local part of one component of \a particle_mass
     sll_real64, allocatable :: particle_mass_23_local(:,:) !< MPI-processor local part of one component of \a particle_mass
     sll_real64, allocatable :: particle_mass_32_local(:,:) !< MPI-processor local part of one component of \a particle_mass
     sll_real64, allocatable :: particle_mass_31_local(:,:) !< MPI-processor local part of one component of \a particle_mass
     sll_real64, allocatable :: particle_mass_13_local(:,:) !< MPI-processor local part of one component of \a particle_mass

     sll_real64 :: solver_tolerance, iter_tolerance !< solver and iteration tolerance
     sll_int32 :: max_iter !< maximal amount of iterations

     sll_int32 :: boundary_particles = 100 !< particle boundary conditions
     sll_int32 :: counter_left = 0 !< boundary counter
     sll_int32 :: counter_right = 0 !< boundary counter
     sll_real64, pointer     :: rhob(:) => null() !< charge at the boundary

     sll_real64 :: force_sign = 1._f64 !< sign of particle force
     logical :: adiabatic_electrons = .false. !< true for simulation with adiabatic electrons
     logical :: jmean = .false. !< logical for mean value of current 
     logical :: lindf = .false. !< true for simulation with linear delta f method
     logical :: boundary = .false. !< true for non periodic boundary conditions

     sll_int32 :: n_particles_max !< maximal number of particles
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

   contains
     procedure :: advect_x => advect_x_pic_vm_3d3v_trafo !> Advect x-part
     procedure :: advect_vb => advect_vb_pic_vm_3d3v_trafo !> Push v, vxB-part only
     procedure :: advect_eb => advect_eb_pic_vm_3d3v_trafo !> Solve Faraday and B-part of Ampere
     procedure :: advect_e => advect_e_pic_vm_3d3v_trafo !> Advect ev-part
     procedure :: advect_ex => advect_ex_pic_vm_3d3v_trafo !> Advect ev-part together with x-part in nonlinear iteration

     procedure :: init => initialize_pic_vm_3d3v_trafo !> Initialize the type
     procedure :: init_from_file => initialize_file_pic_vm_3d3v_trafo !> Initialize from nml file
     procedure :: free => delete_pic_vm_3d3v_trafo !> Finalization

  end type sll_t_time_propagator_pic_vm_3d3v_trafo_helper

contains


  !---------------------------------------------------------------------------!
  !> advect_x: Equations to be solved
  !> $\Xi^{n+1}=\Xi^n+ \frac{\Delta t}{2} (DF^{-1}(\Xi^{n+1})+DF^{-1}(\Xi^n)) V^n$
  subroutine advect_x_pic_vm_3d3v_trafo ( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_trafo_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i, j, i_sp, iter
    sll_real64 :: xi(3), xnew(3), xh(3), vi(3), vt(3), jmat(3,3), wall(3)
    sll_real64 :: err
    iter = 0
    do i_sp = 1, self%particle_group%n_species
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          ! Read out particle position and velocity
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          if( self%map%inverse)then
             xh = self%map%get_x(xi)
             xh = xh + dt * vi
             xnew = self%map%get_xi(xh)
             call sll_s_compute_particle_boundary( self, xi, xnew, vi )
          else
             !Predictor-Corrector with loop for corrector step
             jmat=self%map%jacobian_matrix_inverse(xi)
             !Transformation of the v coordinates 
             do j = 1, 3
                vt(j) = jmat(j,1)*vi(1) + jmat(j,2)*vi(2) + jmat(j,3)*vi(3)
             end do
             !x^\star=\mod(x^n+dt*DF^{-1}(x^n)v^n,1)
             xnew = xi + dt * vt
             i = 0
             err= maxval(abs(xi - xnew))
             do while(i < self%max_iter .and. err > self%iter_tolerance)
                call sll_s_compute_particle_boundary_simple( self%boundary_particles, self%counter_left, self%counter_right, xi, xnew )
                jmat=self%map%jacobian_matrix_inverse( xnew )
                do j = 1, 3
                   xh(j) = xi(j) + 0.5_f64 * dt * (vt(j) + jmat(j,1)*vi(1) + jmat(j,2)*vi(2)+ jmat(j,3)*vi(3))
                end do
                err = maxval(abs(xh - xnew))
                xnew = xh
                i = i + 1
             end do
             iter = iter + i - 1
             call sll_s_compute_particle_boundary( self, xi, xnew, vi )
          end if
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

  end subroutine advect_x_pic_vm_3d3v_trafo


  !> Helper function for \a advect_x
  subroutine sll_s_compute_particle_boundary( self, xold, xnew, vi  )
    class(sll_t_time_propagator_pic_vm_3d3v_trafo_helper), intent( inout ) :: self 
    sll_real64,                                           intent( inout ) :: xold(3)
    sll_real64,                                           intent( inout ) :: xnew(3)
    sll_real64,                                           intent( inout ) :: vi(3)
    !local variables
    sll_real64 :: xmid(3), xt(3), xbar, dx
    sll_real64 :: jmatrix(3,3)

    if(xnew(1) < 0._f64 .or. xnew(1) > 1._f64 )then
       if(xnew(1) < 0._f64  )then
          xbar = 0._f64
          self%counter_left = self%counter_left+1
       else if(xnew(1) > 1._f64)then
          xbar = 1._f64
          self%counter_right = self%counter_right+1
       end if
       dx = (xbar- xold(1))/(xnew(1)-xold(1))
       xmid = xold + dx * (xnew-xold)
       xmid(1) = xbar

       select case(self%boundary_particles)
       case(sll_p_boundary_particles_singular)
          if(xnew(1) < 0._f64 )then
             xnew(1) = -xnew(1)
             xnew(2) = xnew(2) + 0.5_f64
          else if(xnew(1) > 1._f64 )then
             jmatrix = self%map%jacobian_matrix_inverse_transposed(xmid)
             vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
             xnew(1) = 2._f64 - xnew(1)
             xnew(2) = 1._f64 - xnew(2)
          end if
       case(sll_p_boundary_particles_reflection)
          xt = xmid
          xt(2:3) = modulo(xt(2:3),1._f64)
          jmatrix = self%map%jacobian_matrix_inverse_transposed(xt)
          vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
          xnew(1) = 2._f64*xbar-xnew(1)
       case(sll_p_boundary_particles_absorption)
       case( sll_p_boundary_particles_periodic)
          xnew(1) = modulo(xnew(1), 1._f64)
       case default
          xnew(1) = modulo(xnew(1), 1._f64)
       end select
    else if(xnew(1) < -1._f64 .or.  xnew(1) > 2._f64 ) then
       print*, xnew
       SLL_ERROR("particle boundary", "particle out of bound")
    end if
    xnew(2:3) = modulo(xnew(2:3), 1._f64)

  end subroutine sll_s_compute_particle_boundary


  !---------------------------------------------------------------------------!
  !> advect_vb: Equations to be solved
  !> $(\mathbb{I}-\Delta \frac{\Delta t q}{2 m}  DF^{-\top} \mathbb{B}(\Xi^n,b^n) DF^{-1}) V^{n+1}=(\mathbb{I}+ \frac{\Delta t q}{2 m} DF^{-\top} \mathbb{B}(\Xi^n,b^n) DF^{-1}) V^n$
  subroutine advect_vb_pic_vm_3d3v_trafo ( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_trafo_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                               intent( in    ) :: dt   !< time step
    !local variables
    sll_int32  :: i_part, i_sp, j
    sll_real64 :: qmdt
    sll_real64 :: vi(3), xi(3), wall(3)
    sll_real64 :: bfield(3), jmatrix(3,3), rhs(3)
    sll_real64 :: vt(3), c(3)

    do i_sp = 1, self%particle_group%n_species
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

          jmatrix=self%map%jacobian_matrix_inverse(xi)

          !VT = DF^{-1} * vi
          do j=1,3
             vt(j)=jmatrix(j,1)*vi(1)+jmatrix(j,2)*vi(2)+jmatrix(j,3)*vi(3)
          end do
          !c = VT x bfield 
          c(1)=bfield(3)*vt(2)-bfield(2)*vt(3)
          c(2)=bfield(1)*vt(3)-bfield(3)*vt(1)
          c(3)=bfield(2)*vt(1)-bfield(1)*vt(2)
          !rhs = vi + sign * DF^{-T} * c
          do j=1,3
             rhs(j)= vi(j) + qmdt*(jmatrix(1,j)*c(1)+jmatrix(2,j)*c(2)+jmatrix(3,j)*c(3))
          end do

          call sll_s_compute_matrix_inverse(rhs, vi, bfield, jmatrix, qmdt)

          call self%particle_group%group(i_sp)%set_v( i_part, vi )
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xi, vi, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do
    end do

  end subroutine advect_vb_pic_vm_3d3v_trafo


  !---------------------------------------------------------------------------!
  !> advect_eb: Equations to be solved
  !> Solution with Schur complement: $ S=M_1+\frac{\Delta t^2}{4} C^\top M_2 C $
  !> $ e^{n+1}=S^{-1}( (M_1-\frac{\Delta t^2}{4} C^\top M_2 C)e^n+\Delta t C^\top M_2 b^n) $
  !> $ b^{n+1}=b^n-\frac{\Delta t}{2} C(e^n+e^{n+1}) $
  subroutine advect_eb_pic_vm_3d3v_trafo ( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_trafo_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                                  intent( in    ) :: dt   !< time step

    call self%maxwell_solver%compute_curl_part( dt, self%efield_dofs, self%bfield_dofs, self%betar(1) )

  end subroutine advect_eb_pic_vm_3d3v_trafo


  !---------------------------------------------------------------------------!
  !> advect_e: Equations to be solved
  !> Solution with Schur complement: $ S_{+}=M_1+\frac{\Delta t^2 q^2}{4 m} (\mathbb{\Lambda}^1)^T DF^{-1} DF^{-T} \mathbb{\Lambda}^1 $
  !> $e^{n+1}=S_{+}^{-1}\left(S_{-}e^n-\Delta t (\mathbb{\Lambda}^1)^\top DF^{-1}\mathbb{W}_q V^n \right)$
  !> $V^{n+1}=V^n+\frac{\Delta t}{2} \mathbb{W}_{\frac{q}{m}} DF^{-\top} \mathbb{\Lambda}^1(e^{n+1}+e^n)$
  subroutine advect_e_pic_vm_3d3v_trafo( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_trafo_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: dt   !< time step
    ! local variables
    sll_int32 :: i_part, j, i_sp
    sll_real64 :: vi(3), vt(3), xi(3), wi(1), wall(3)
    sll_real64 :: metric(3,3), jmat(3,3)
    sll_real64 :: factor, qoverm
    sll_real64 :: efield(3), ephys(3)
    sll_real64 :: rhs(self%n_total1+2*self%n_total0)

    ! Set to zero
    self%j_dofs_local = 0.0_f64
    self%j_dofs = 0.0_f64
    self%particle_mass_1_local = 0.0_f64
    self%particle_mass_2_local = 0.0_f64
    self%particle_mass_3_local = 0.0_f64
    if(self%map%flag2d)then
       self%particle_mass_12_local = 0.0_f64
       self%particle_mass_21_local = 0.0_f64
       if(self%map%flag3d)then
          self%particle_mass_23_local = 0.0_f64
          self%particle_mass_32_local = 0.0_f64
          self%particle_mass_31_local = 0.0_f64
          self%particle_mass_13_local = 0.0_f64
       end if
    end if

    self%particle_mass_1%particle_mass = 0.0_f64
    self%particle_mass_2%particle_mass = 0.0_f64
    self%particle_mass_3%particle_mass = 0.0_f64
    if(self%map%flag2d)then
       self%particle_mass_12%particle_mass = 0.0_f64
       self%particle_mass_21%particle_mass = 0.0_f64
       if(self%map%flag3d)then
          self%particle_mass_23%particle_mass = 0.0_f64
          self%particle_mass_32%particle_mass = 0.0_f64
          self%particle_mass_13%particle_mass = 0.0_f64
          self%particle_mass_31%particle_mass = 0.0_f64
       end if
    end if

    ! First particle loop
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       if( self%adiabatic_electrons ) then
          factor = dt**2 * 0.25_f64* qoverm
       else
          factor = dt**2 * self%betar(2) * 0.25_f64* qoverm
       end if
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part, self%i_weight)

          metric= self%map%metric_inverse( xi )
          jmat=self%map%jacobian_matrix_inverse( xi )
          do j=1,3
             vt(j)=vi(1)*jmat(j,1)+vi(2)*jmat(j,2)+vi(3)*jmat(j,3)
          end do

          ! Accumulate j
          call self%particle_mesh_coupling%add_charge( xi, wi(1)*vt(1), &
               [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
               self%j_dofs_local(1:self%n_total1) )
          call self%particle_mesh_coupling%add_charge( xi, wi(1)*vt(2), &
               [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
               self%j_dofs_local(1+self%n_total1:self%n_total1+self%n_total0) )
          call self%particle_mesh_coupling%add_charge( xi, wi(1)*vt(3), &
               [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
               self%j_dofs_local(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0) )


          ! Accumulate the particle mass matrix 
          call self%particle_mesh_coupling%add_particle_mass( xi, wi(1)*metric(1,1) * factor, &
               [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
               self%particle_mass_1_local )
          call self%particle_mesh_coupling%add_particle_mass( xi, wi(1)*metric(2,2) * factor, &
               [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
               self%particle_mass_2_local )
          call self%particle_mesh_coupling%add_particle_mass( xi, wi(1)*metric(3,3) * factor, &
               [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
               self%particle_mass_3_local )
          if(self%map%flag2d)then
             call self%particle_mesh_coupling%add_particle_mass_od( xi, wi(1)*metric(1,2) * factor, &
                  [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
                  [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)],&
                  self%particle_mass_12_local )
             call self%particle_mesh_coupling%add_particle_mass_od( xi, wi(1)*metric(2,1) * factor, &
                  [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)],&
                  [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
                  self%particle_mass_21_local )
             if(self%map%flag3d)then
                call self%particle_mesh_coupling%add_particle_mass_od( xi, wi(1)*metric(2,3) * factor, &
                     [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
                     [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
                     self%particle_mass_23_local )
                call self%particle_mesh_coupling%add_particle_mass_od( xi, wi(1)*metric(3,2) * factor, &
                     [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
                     [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
                     self%particle_mass_32_local )
                call self%particle_mesh_coupling%add_particle_mass_od( xi, wi(1)*metric(3,1) * factor, &
                     [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
                     [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
                     self%particle_mass_31_local )
                call self%particle_mesh_coupling%add_particle_mass_od( xi, wi(1)*metric(1,3) * factor, &
                     [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
                     [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
                     self%particle_mass_13_local )
             end if
          end if


          ! Evaulate E at particle position and propagate v a half step       
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
               self%efield_dofs(1:self%n_total1), efield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
               self%efield_dofs(1+self%n_total1:self%n_total1+self%n_total0), efield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
               self%efield_dofs(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0), efield(3))

          ! physical efield, jmatrix is not transposed
          do j=1, 3
             ephys(j) = jmat(1,j)* efield(1)+jmat(2,j)* efield(2)+jmat(3,j)* efield(3)
          end do

          ! velocity update 
          vi = vi + dt* 0.5_f64* qoverm * ephys

          call self%particle_group%group(i_sp)%set_v( i_part, vi )
          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xi, vi, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do
    end do

    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
         self%n_total1+self%n_total0*2, MPI_SUM, self%j_dofs )

    if ( self%jmean ) then
       self%j_dofs(1:self%n_total1) = self%j_dofs(1:self%n_total1) - sum(self%j_dofs(1:self%n_total1))/real(self%n_total1, f64)
       self%j_dofs(1+self%n_total1:self%n_total1+self%n_total0) = self%j_dofs(1+self%n_total1:self%n_total1+self%n_total0) - sum(self%j_dofs(1+self%n_total1:self%n_total1+self%n_total0))/real(self%n_total0, f64)
       self%j_dofs(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0) = self%j_dofs(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0) - sum(self%j_dofs(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0))/real(self%n_total0, f64)
    end if

    call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_1_local, &
         self%n_total1*self%nspan_1(1) , MPI_SUM, self%particle_mass_1%particle_mass)
    call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_2_local, &
         self%n_total0*self%nspan_1(2) , MPI_SUM, self%particle_mass_2%particle_mass)
    call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_3_local, &
         self%n_total0*self%nspan_1(3) , MPI_SUM, self%particle_mass_3%particle_mass)

    if(self%map%flag2d)then
       call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_12_local, &
            self%n_total1*self%nspan_2(1) , MPI_SUM, self%particle_mass_12%particle_mass)
       call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_21_local, &
            self%n_total0*self%nspan_2(1) , MPI_SUM, self%particle_mass_21%particle_mass)
       if(self%map%flag3d)then
          call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_23_local, &
               self%n_total0*self%nspan_2(2) , MPI_SUM, self%particle_mass_23%particle_mass)
          call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_31_local, &
               self%n_total0*self%nspan_2(3) , MPI_SUM, self%particle_mass_31%particle_mass)
          call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_32_local, &
               self%n_total0*self%nspan_2(2) , MPI_SUM, self%particle_mass_32%particle_mass)
          call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_13_local, &
               self%n_total1*self%nspan_2(3) , MPI_SUM, self%particle_mass_13%particle_mass)
       end if
    end if

    call self%particle_mass_op%set( 1, 1, self%particle_mass_1 )
    call self%particle_mass_op%set( 2, 2, self%particle_mass_2 )
    call self%particle_mass_op%set( 3, 3, self%particle_mass_3 )
    if(self%map%flag2d)then
       call self%particle_mass_op%set( 1, 2, self%particle_mass_12 )
       call self%particle_mass_op%set( 2, 1, self%particle_mass_21 )
       if(self%map%flag3d)then
          call self%particle_mass_op%set( 2, 3, self%particle_mass_23 )
          call self%particle_mass_op%set( 3, 2, self%particle_mass_32 )
          call self%particle_mass_op%set( 3, 1, self%particle_mass_31 )
          call self%particle_mass_op%set( 1, 3, self%particle_mass_13 )
       end if
    end if

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
       if( self%boundary ) then !perfect boundary conditions
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

    if( self%boundary ) then !perfect boundary conditions
       do j = 1, self%maxwell_solver%n_dofs(2)*self%maxwell_solver%n_dofs(3)
          self%efield_dofs(self%n_total1+1+(j-1)*self%maxwell_solver%n_dofs(1)) = 0._f64
          self%efield_dofs(self%n_total1+j*self%maxwell_solver%n_dofs(1)) = 0._f64
          self%efield_dofs(self%n_total1+self%n_total0+1+(j-1)*self%maxwell_solver%n_dofs(1)) = 0._f64
          self%efield_dofs(self%n_total1+self%n_total0+j*self%maxwell_solver%n_dofs(1)) = 0._f64
       end do
    end if
    ! Second particle loop (second half step of particle propagation)
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)

          ! Evaulate E at particle position and propagate v a half step
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
               self%efield_dofs(1:self%n_total1), efield(1))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
               self%efield_dofs(1+self%n_total1:self%n_total1+self%n_total0), efield(2))
          call self%particle_mesh_coupling%evaluate &
               (xi, [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
               self%efield_dofs(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0), efield(3))

          jmat = self%map%jacobian_matrix_inverse_transposed(xi)
          do j=1, 3
             ephys(j) = jmat(j,1)* efield(1)+jmat(j,2)* efield(2)+jmat(j,3)* efield(3)
          end do

          ! velocity update 
          vi = vi + dt* 0.5_f64* qoverm* ephys

          call self%particle_group%group(i_sp)%set_v( i_part, vi )

          ! Update particle weights
          if (self%particle_group%group(i_sp)%n_weights == 3 ) then
             wall = self%particle_group%group(i_sp)%get_weights(i_part)
             wall(3) = self%control_variate%cv(i_sp)%update_df_weight( xi, vi, 0.0_f64, wall(1), wall(2) )
             call self%particle_group%group(i_sp)%set_weights( i_part, wall )
          end if
       end do
    end do

  end subroutine advect_e_pic_vm_3d3v_trafo


  !---------------------------------------------------------------------------!
  !> advect_ex: Equations to be solved
  !> $\frac{\Xi^{n+1}-\Xi^n}{\Delta t}=\frac{DF^{-1}(\Xi^{n+1})+DF^{-1}(\Xi^n)}{2} \frac{V^{n+1}+V^n}{2}$
  !> $\frac{V^{n+1}-V^n}{\Delta t}=\mathbb{W}_{\frac{q}{m}} \frac{DF^{-\top}(\Xi^{n+1})+DF^{-\top}(\Xi^n)}{2} \frac{1}{\Delta t}\int_{t^n}^{t^{n+1}} \mathbb{\Lambda}^1(\Xi(\tau))  d\tau \frac{e^{n+1}+e^n}{2}$
  !> $\frac{M_1 e^{n+1}-M_1 e^n}{\Delta t} = - \frac{1}{\Delta t} \int_{t^n}^{t^{n+1}} \mathbb{\Lambda}^1(\Xi(\tau))^\top  d\tau  \frac{DF^{-1}(\Xi^{n+1})+DF^{-1}(\Xi^n)}{2} \mathbb{W}_q\frac{V^{n+1}+V^n}{2}$
  subroutine advect_ex_pic_vm_3d3v_trafo( self, dt )
    class(sll_t_time_propagator_pic_vm_3d3v_trafo_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: dt   !< time step
    ! local variables
    sll_int32 :: i_part, i_sp, j
    sll_real64 :: vi(3), vh(3), xi(3), xs(3), wi(1), xnew(3), vbar(3)
    sll_real64 :: qoverm, wall(3)
    sll_real64 :: jmat(3,3), jmatrix(3,3)
    sll_int32 :: niter
    sll_real64 :: residual(1), residual_local(1)
    sll_real64 :: xmid(3), xt(3), xbar, dx
    sll_real64 :: efield(3)

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
             vbar = 0.5_f64 * (self%vnew(i_sp,:, i_part)+vi)
             xnew = self%xnew(i_sp,:, i_part)

             if( self%map%inverse) then
                xs = self%map%get_x(xi)
                xs = xs + dt * vbar
                xnew = self%map%get_xi(xs)
             else
                jmat = self%map%jacobian_matrix_inverse_transposed( xi )
                jmatrix=self%map%jacobian_matrix_inverse_transposed( xnew )
                do j = 1, 3
                   vh(j) = 0.5_f64 * ((jmatrix(1,j)+jmat(1,j))*vbar(1) + (jmatrix(2,j)+jmat(2,j))*vbar(2) + (jmatrix(3,j)+ jmat(3,j))*vbar(3))
                end do
                xnew = xi + dt * vh
             end if
             
             if(xnew(1) < -1._f64 .or. xnew(1) > 2._f64)then
                print*, xnew
                SLL_ERROR("particle boundary", "particle out of bound")
             else if(xnew(1) < 0._f64 .or. xnew(1) > 1._f64 )then
                if(xnew(1) < 0._f64  )then
                   xbar = 0._f64
                   self%counter_left = self%counter_left+1
                else if(xnew(1) > 1._f64)then
                   xbar = 1._f64
                   self%counter_right = self%counter_right+1
                end if
                dx = (xbar- xi(1))/(xnew(1)-xi(1))
                xmid = xi + dx * (xnew-xi)
                xmid(1) = xbar

                vh = (xmid-xi) * wi(1) 
                call self%particle_mesh_coupling%add_current_evaluate( xi, xmid, vh, self%efield_dofs_work, &
                     self%j_dofs_local, efield )

                xt = xmid
                xt(2:3) = modulo(xt(2:3), 1._f64)
                jmatrix = self%map%jacobian_matrix_inverse_transposed(xt)
                do j = 1, 3
                   vi(j) = vi(j) + dx*dt*qoverm*0.5_f64*((jmatrix(j,1)+jmat(j,1))*efield(1) + (jmatrix(j,2)+jmat(j,2))*efield(2) + (jmatrix(j,3)+ jmat(j,3))*efield(3))
                end do
                select case(self%boundary_particles)
                case(sll_p_boundary_particles_singular)
                   call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
                   xmid(2) = xbar + (1._f64-2._f64*xbar)*xmid(2) + 0.5_f64-0.5_f64*xbar
                   xt(2:3) = modulo(xt(2:3), 1._f64)
                   jmatrix = self%map%jacobian_matrix_inverse_transposed(xt)
                   call self%particle_mesh_coupling%add_charge(xmid, -wi(1), self%spline_degree, self%rhob)
                   xnew(1) = 2._f64*xbar-xnew(1)
                   xnew(2) = xbar + (1._f64-2._f64*xbar)*xnew(2) + 0.5_f64-0.5_f64*xbar
                   if(xnew(1) > 1._f64 )then
                      vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
                   end if
                case(sll_p_boundary_particles_reflection)
                   xnew(1) = 2._f64*xbar-xnew(1)
                   vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
                case(sll_p_boundary_particles_absorption)
                   call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
                   xnew(1) = xmid(1) + (xbar-0.5_f64) * 1.9_f64* self%iter_tolerance 
                case( sll_p_boundary_particles_periodic)
                   xnew(1) = modulo(xnew(1), 1._f64)
                   xmid(1) = 1._f64-xbar
                case default
                   xnew(1) = modulo(xnew(1), 1._f64)
                   xmid(1) = 1._f64-xbar
                end select
                if( abs(xnew(1)-xmid(1)) > self%iter_tolerance ) then
                   vh = (xnew - xmid)*wi(1)
                   call self%particle_mesh_coupling%add_current_evaluate( xmid, xnew, vh, self%efield_dofs_work, &
                        self%j_dofs_local, efield )
                   xnew(2:3) = modulo(xnew(2:3), 1._f64)
                   jmat = self%map%jacobian_matrix_inverse_transposed(xnew)
                   do j = 1, 3
                      vi(j) = vi(j) + (1._f64-dx)*dt*qoverm*0.5_f64*((jmatrix(j,1)+jmat(j,1))*efield(1) + (jmatrix(j,2)+jmat(j,2))*efield(2) + (jmatrix(j,3)+ jmat(j,3))*efield(3))
                   end do
                   if(self%boundary_particles == sll_p_boundary_particles_reflection) then
                      do j = 1, 3
                         vh(j) = jmat(j,1) *vi(j)
                      end do
                      vh(1) = - vh(1)
                      jmat = self%map%jacobian_matrix(xnew)
                      do j = 1, 3
                         vi(j) = jmat(1,j) *vh(j)
                      end do
                   end if
                else
                   xnew(1) = xmid(1)   
                end if
             else   
                vh = (xnew-xi) * wi(1) 
                call self%particle_mesh_coupling%add_current_evaluate( xi, xnew, vh, self%efield_dofs_work, &
                     self%j_dofs_local, efield )
                xnew(2:3) =  modulo(xnew(2:3), 1._f64)
                jmatrix = self%map%jacobian_matrix_inverse_transposed(xnew)
                do j = 1, 3
                   vi(j) = vi(j) + dt*qoverm *0.5_f64*((jmatrix(j,1)+jmat(j,1))*efield(1) + (jmatrix(j,2)+jmat(j,2))*efield(2) + (jmatrix(j,3)+ jmat(j,3))*efield(3))
                end do
             end if
             
             self%xnew(i_sp, :, i_part) = xnew
             self%vnew(i_sp, :, i_part) = vi
          end do
       end do

       self%j_dofs = 0.0_f64
       ! MPI to sum up contributions from each processor
       call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local, &
            self%n_total1+2*self%n_total0, MPI_SUM, self%j_dofs )

       if( self%adiabatic_electrons) then
          self%phi_dofs_work = self%phi_dofs
          call self%maxwell_solver%compute_phi_from_j( self%j_dofs, self%phi_dofs_work, self%efield_dofs_work )

          ! Compute residual based on phi
          residual_local = (sum((self%phi_dofs_work-self%phi_dofs_new)**2))*product(self%particle_mesh_coupling%delta_x)
          call sll_o_collective_allreduce( sll_v_world_collective, residual_local, 1, MPI_MAX, residual )
          residual = sqrt(residual)
       else
          self%efield_dofs_work = self%efield_dofs
          call self%maxwell_solver%compute_E_from_j( self%betar(2)*self%j_dofs, self%efield_dofs_work )


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
          vbar = 0.5_f64 * (self%vnew(i_sp,:, i_part)+vi)
          xnew = self%xnew(i_sp,:, i_part)

          if( self%map%inverse) then
             xs = self%map%get_x(xi)
             xs = xs + dt * vbar
             xnew = self%map%get_xi(xs)
          else
             jmat = self%map%jacobian_matrix_inverse_transposed( xi )
             jmatrix=self%map%jacobian_matrix_inverse_transposed( xnew )
             do j = 1, 3
                vh(j) = 0.5_f64 * ((jmatrix(1,j)+jmat(1,j))*vbar(1) + (jmatrix(2,j)+jmat(2,j))*vbar(2) + (jmatrix(3,j)+ jmat(3,j))*vbar(3))
             end do
             xnew = xi + dt * vh
          end if

          call compute_particle_boundary_current_evaluate( self, xi, xnew, vi, wi, dt*qoverm )


          call self%particle_group%group(i_sp)%set_v( i_part, vi )
          call self%particle_group%group(i_sp)%set_x( i_part, xnew )
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
       call self%maxwell_solver%compute_E_from_j( self%betar(2)*self%j_dofs, self%efield_dofs )
    end if

    self%iter_counter = self%iter_counter + 1
    self%niter(self%iter_counter) = niter

  end subroutine advect_ex_pic_vm_3d3v_trafo


  !> Helper function for \a advect_ex
  subroutine compute_particle_boundary_current_evaluate( self, xi, xnew, vi, wi, sign )
    class(sll_t_time_propagator_pic_vm_3d3v_trafo_helper), intent( inout ) :: self !< time splitting object 
    sll_real64,                                           intent( in    ) :: xi(3)
    sll_real64,                                           intent( inout ) :: xnew(3)
    sll_real64,                                           intent( inout ) :: vi(3)
    sll_real64,                                           intent( in    ) :: wi(1)
    sll_real64,                                           intent( in    ) :: sign
    !local variables
    sll_real64 :: xmid(3), xt(3), vh(3), xbar, dx
    sll_real64 :: jmatrix(3,3), jmat(3,3), efield(3)
    sll_int32 :: j

    jmat = self%map%jacobian_matrix_inverse_transposed( xi )

    if(xnew(1) < -1._f64 .or. xnew(1) > 2._f64)then
       print*, xnew
       SLL_ERROR("particle boundary", "particle out of bound")
    else if(xnew(1) < 0._f64 .or. xnew(1) > 1._f64 )then
       if(xnew(1) < 0._f64  )then
          xbar = 0._f64
          self%counter_left = self%counter_left+1
       else if(xnew(1) > 1._f64)then
          xbar = 1._f64
          self%counter_right = self%counter_right+1
       end if
       dx = (xbar- xi(1))/(xnew(1)-xi(1))
       xmid = xi + dx * (xnew-xi)
       xmid(1) = xbar

       vh = (xmid-xi) * wi(1) 
       call self%particle_mesh_coupling%add_current_evaluate( xi, xmid, vh, self%efield_dofs_work, &
            self%j_dofs_local, efield )

       xt = xmid
       xt(2:3) = modulo(xt(2:3), 1._f64)
       jmatrix = self%map%jacobian_matrix_inverse_transposed(xt)
       do j = 1, 3
          vi(j) = vi(j) + dx*sign *0.5_f64*((jmatrix(j,1)+jmat(j,1))*efield(1) + (jmatrix(j,2)+jmat(j,2))*efield(2) + (jmatrix(j,3)+ jmat(j,3))*efield(3))
       end do
       select case(self%boundary_particles)
       case(sll_p_boundary_particles_singular)
          call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
          xmid(2) = xbar + (1._f64-2._f64*xbar)*xmid(2) + 0.5_f64-0.5_f64*xbar
          xt(2:3) = modulo(xt(2:3), 1._f64)
          jmatrix = self%map%jacobian_matrix_inverse_transposed(xt)
          call self%particle_mesh_coupling%add_charge(xmid, -wi(1), self%spline_degree, self%rhob)
          xnew(1) = 2._f64*xbar-xnew(1)
          xnew(2) = xbar + (1._f64-2._f64*xbar)*xnew(2) + 0.5_f64-0.5_f64*xbar
          if(xnew(1) > 1._f64 )then
             vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
          end if
       case(sll_p_boundary_particles_reflection)
          vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
          xnew(1) = 2._f64*xbar-xnew(1)
       case(sll_p_boundary_particles_absorption)
          call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
       case( sll_p_boundary_particles_periodic)
          xnew(1) = modulo(xnew(1), 1._f64)
          xmid(1) = 1._f64-xbar
       case default
          xnew(1) = modulo(xnew(1), 1._f64)
          xmid(1) = 1._f64-xbar
       end select
       if( abs(xnew(1)-xmid(1)) > self%iter_tolerance ) then
          vh = (xnew - xmid)*wi(1)
          call self%particle_mesh_coupling%add_current_evaluate( xmid, xnew, vh, self%efield_dofs_work, &
               self%j_dofs_local, efield )
          xnew(2:3) = modulo(xnew(2:3), 1._f64)
          jmat = self%map%jacobian_matrix_inverse_transposed(xnew)
          do j = 1, 3
             vi(j) = vi(j) + (1._f64-dx)*sign *0.5_f64*((jmatrix(j,1)+jmat(j,1))*efield(1) + (jmatrix(j,2)+jmat(j,2))*efield(2) + (jmatrix(j,3)+ jmat(j,3))*efield(3))
          end do
       else
          xnew(1) = xmid(1)   
       end if
    else   
       vh = (xnew-xi) * wi(1) 
       call self%particle_mesh_coupling%add_current_evaluate( xi, xnew, vh, self%efield_dofs_work, &
            self%j_dofs_local, efield )
       xnew(2:3) =  modulo(xnew(2:3), 1._f64)
       jmatrix = self%map%jacobian_matrix_inverse_transposed(xnew)
       do j = 1, 3
          vi(j) = vi(j) + sign *0.5_f64*((jmatrix(j,1)+jmat(j,1))*efield(1) + (jmatrix(j,2)+jmat(j,2))*efield(2) + (jmatrix(j,3)+ jmat(j,3))*efield(3))
       end do
    end if

  end subroutine compute_particle_boundary_current_evaluate


  !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_3d3v_trafo(&
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
       solver_tolerance, &
       iter_tolerance, max_iter, &
       betar,&
       force_sign, &
       rhob, &
       control_variate, &
       jmean, &
       lindf) 
    class(sll_t_time_propagator_pic_vm_3d3v_trafo_helper), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,          intent( in ) :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target, intent(in) :: particle_mesh_coupling !< Particle mesh coupling
    class(sll_t_particle_array), target,           intent( in ) :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:) !< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent( in ) :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent( in ) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                    intent( in ) :: x_min(3) !< Lower bound of x domain
    sll_real64,                                    intent( in ) :: Lx(3) !< Length of the domain in x direction.
    type(sll_t_mapping_3d), target,                intent( inout ) :: map !< Coordinate transformation
    sll_int32, optional,                           intent( in ) :: boundary_particles !< particle boundary conditions
    sll_real64, optional,                          intent( in ) :: solver_tolerance !< Solver tolerance
    sll_real64, optional,                          intent( in ) :: iter_tolerance !< iteration tolerance
    sll_int32,  optional,                          intent( in ) :: max_iter !< maximal number of iterations
    sll_real64, optional,                          intent( in ) :: betar(2) !< reciprocal plasma beta
    sll_real64, optional, intent(in) :: force_sign !< sign of particle force
    sll_real64, optional, target, intent( in ) :: rhob(:) !< charge at the boundary
    class(sll_t_control_variates), optional, target, intent(in) :: control_variate !< Control variate (if delta f)
    logical, optional, intent(in) :: jmean !< logical for mean value of current
    logical, optional, intent(in) :: lindf !< true for linear delta f method
    !local variables
    sll_int32 :: ierr, j

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
    self%map => map
    self%n_total0 = self%maxwell_solver%n_total0
    self%n_total1 = self%maxwell_solver%n_total1
    self%spline_degree = self%particle_mesh_coupling%spline_degree
    self%x_min = self%map%get_x([0._f64,0._f64,0._f64])
    self%Lx = self%map%Lx
    self%x_max = x_min + Lx


    self%nspan_1(1) =  (2*self%spline_degree(1)-1)*(2*self%spline_degree(2)+1)*(2*self%spline_degree(3)+1)
    self%nspan_1(2) =  (2*self%spline_degree(1)+1)*(2*self%spline_degree(2)-1)*(2*self%spline_degree(3)+1)
    self%nspan_1(3) =  (2*self%spline_degree(1)+1)*(2*self%spline_degree(2)+1)*(2*self%spline_degree(3)-1)

    self%nspan_2(1) =  (2*self%spline_degree(1))*(2*self%spline_degree(2))*(2*self%spline_degree(3)+1)
    self%nspan_2(2) =  (2*self%spline_degree(1)+1)*(2*self%spline_degree(2))*(2*self%spline_degree(3))
    self%nspan_2(3) =  (2*self%spline_degree(1))*(2*self%spline_degree(2)+1)*(2*self%spline_degree(3))

    SLL_ALLOCATE( self%vec1(1:self%n_total1+self%n_total0*2), ierr )
    SLL_ALLOCATE( self%vec2(1:self%n_total1+self%n_total0*2), ierr )
    self%vec1 = 0._f64
    self%vec2 = 0._f64
    SLL_ALLOCATE( self%j_dofs(1:self%n_total1+self%n_total0*2), ierr )
    SLL_ALLOCATE( self%j_dofs_local(1:self%n_total1+self%n_total0*2), ierr )
    SLL_ALLOCATE( self%particle_mass_1_local(1:self%nspan_1(1), 1:self%n_total1), ierr )
    SLL_ALLOCATE( self%particle_mass_2_local(1:self%nspan_1(2), 1:self%n_total0), ierr )
    SLL_ALLOCATE( self%particle_mass_3_local(1:self%nspan_1(3), 1:self%n_total0), ierr )
    SLL_ALLOCATE( self%particle_mass_12_local(1:self%nspan_2(1), 1:self%n_total1), ierr )
    SLL_ALLOCATE( self%particle_mass_23_local(1:self%nspan_2(2), 1:self%n_total0), ierr )
    SLL_ALLOCATE( self%particle_mass_31_local(1:self%nspan_2(3), 1:self%n_total0), ierr )
    SLL_ALLOCATE( self%particle_mass_21_local(1:self%nspan_2(1), 1:self%n_total0), ierr )
    SLL_ALLOCATE( self%particle_mass_32_local(1:self%nspan_2(2), 1:self%n_total0), ierr )
    SLL_ALLOCATE( self%particle_mass_13_local(1:self%nspan_2(3), 1:self%n_total1), ierr )

    if( self%particle_mesh_coupling%n_cells(1)+self%spline_degree(1) == self%maxwell_solver%n_dofs(1)   ) then
       self%boundary = .true.
       allocate( sll_t_linear_operator_particle_mass_cl_3d_diag :: self%particle_mass_1 )
       allocate( sll_t_linear_operator_particle_mass_cl_3d_diag :: self%particle_mass_2 )
       allocate( sll_t_linear_operator_particle_mass_cl_3d_diag :: self%particle_mass_3 )
       if(self%map%flag2d)then
          allocate( sll_t_linear_operator_particle_mass_cl_3d_od :: self%particle_mass_12 )
          allocate( sll_t_linear_operator_particle_mass_cl_3d_od :: self%particle_mass_21 )
          if(self%map%flag3d)then
             allocate( sll_t_linear_operator_particle_mass_cl_3d_od :: self%particle_mass_13 )
             allocate( sll_t_linear_operator_particle_mass_cl_3d_od :: self%particle_mass_31 )
             allocate( sll_t_linear_operator_particle_mass_cl_3d_od :: self%particle_mass_23 )
             allocate( sll_t_linear_operator_particle_mass_cl_3d_od :: self%particle_mass_32 )
          end if
       end if
    else
       allocate( sll_t_linear_operator_particle_mass_3d_diag :: self%particle_mass_1 )
       allocate( sll_t_linear_operator_particle_mass_3d_diag :: self%particle_mass_2 )
       allocate( sll_t_linear_operator_particle_mass_3d_diag :: self%particle_mass_3 )
       if(self%map%flag2d)then
          allocate( sll_t_linear_operator_particle_mass_3d_od :: self%particle_mass_12 )
          allocate( sll_t_linear_operator_particle_mass_3d_od :: self%particle_mass_21 )
          if(self%map%flag3d)then
             allocate( sll_t_linear_operator_particle_mass_3d_od :: self%particle_mass_13 )
             allocate( sll_t_linear_operator_particle_mass_3d_od :: self%particle_mass_31 )
             allocate( sll_t_linear_operator_particle_mass_3d_od :: self%particle_mass_23 )
             allocate( sll_t_linear_operator_particle_mass_3d_od :: self%particle_mass_32 )
          end if
       end if
    end if

    call self%particle_mass_1%create( self%particle_mesh_coupling%n_cells, [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)] )
    call self%particle_mass_2%create( self%particle_mesh_coupling%n_cells, [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)] )
    call self%particle_mass_3%create( self%particle_mesh_coupling%n_cells, [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1] )

    if(self%map%flag2d)then
       call self%particle_mass_12%create( self%particle_mesh_coupling%n_cells, &
            [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
            [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)] )
       call self%particle_mass_21%create( self%particle_mesh_coupling%n_cells, &
            [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
            [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)] )
       if(self%map%flag3d)then
          call self%particle_mass_23%create( self%particle_mesh_coupling%n_cells, &
               [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)], &
               [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1] )
          call self%particle_mass_32%create( self%particle_mesh_coupling%n_cells, &
               [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1], &
               [self%spline_degree(1), self%spline_degree(2)-1, self%spline_degree(3)] )
          call self%particle_mass_31%create( self%particle_mesh_coupling%n_cells, &
               [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1],  &
               [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)] )
          call self%particle_mass_13%create( self%particle_mesh_coupling%n_cells, &
               [self%spline_degree(1)-1, self%spline_degree(2), self%spline_degree(3)], &
               [self%spline_degree(1), self%spline_degree(2), self%spline_degree(3)-1] )
       end if
    end if
    call self%particle_mass_op%create( 3, 3 )

    if( self%adiabatic_electrons ) then
       call self%linear_op_schur_phiv%create( self%maxwell_solver, self%particle_mass_op )
       call self%linear_solver_schur_phiv%create( self%linear_op_schur_phiv )
       self%linear_solver_schur_phiv%atol = self%solver_tolerance
       !self%linear_solver_schur_phiv%verbose = .true.
    else
       call self%linear_op_schur_ev%create( self%maxwell_solver, self%particle_mass_op )
       call self%preconditioner_fft%init( self%maxwell_solver%Lx, self%maxwell_solver%n_cells, self%maxwell_solver%s_deg_0, self%boundary )

       if( self%boundary ) then
          self%vec1 = 1._f64
          call self%maxwell_solver%multiply_mass( [1], self%vec1, self%vec2 )
          do j = 1, self%maxwell_solver%n_dofs(2)*self%maxwell_solver%n_dofs(3)
             self%vec2(self%n_total1+1+(j-1)*self%maxwell_solver%n_dofs(1)) = 1._f64
             self%vec2(self%n_total1+j*self%maxwell_solver%n_dofs(1)) = 1._f64
             self%vec2(self%n_total1+self%n_total0+1+(j-1)*self%maxwell_solver%n_dofs(1)) = 1._f64
             self%vec2(self%n_total1+self%n_total0+j*self%maxwell_solver%n_dofs(1)) = 1._f64
          end do
          self%vec2 = 1._f64/sqrt(abs(self%vec2))
          call self%preconditioner1%create( self%preconditioner_fft%inverse_mass1_3d, self%vec2, 2*self%n_total0+self%n_total1 )
          call self%linear_solver_schur_ev%create( self%linear_op_schur_ev, self%preconditioner1)
       else
          call self%linear_solver_schur_ev%create( self%linear_op_schur_ev, self%preconditioner_fft%inverse_mass1_3d )
       end if
       self%linear_solver_schur_ev%atol = self%solver_tolerance/maxval(self%Lx)
       !self%linear_solver_schur_ev%verbose = .true.
    end if

    if (present(betar)) then
       self%betar = betar
    else
       self%betar = 1.0_f64
    end if

    self%n_particles_max = self%particle_group%group(1)%n_particles
    do j = 2,self%particle_group%n_species       
       self%n_particles_max = max(self%n_particles_max, self%particle_group%group(j)%n_particles )
    end do

    SLL_ALLOCATE( self%xnew(self%particle_group%n_species,3,self%n_particles_max), ierr )
    SLL_ALLOCATE( self%vnew(self%particle_group%n_species,3,self%n_particles_max), ierr )
    SLL_ALLOCATE( self%efield_dofs_new(self%n_total1+2*self%n_total0), ierr )
    SLL_ALLOCATE( self%efield_dofs_work(self%n_total1+2*self%n_total0), ierr )
    SLL_ALLOCATE( self%phi_dofs_work(self%n_total0), ierr )
    SLL_ALLOCATE( self%phi_dofs_new(self%n_total0), ierr )
    self%xnew = 0._f64
    self%vnew = 0._f64
    self%efield_dofs_new = 0._f64
    self%phi_dofs_new = 0._f64
    self%phi_dofs_work = 0._f64

    if (present(control_variate)) then
       allocate(self%control_variate )
       allocate(self%control_variate%cv(self%particle_group%n_species) )
       self%control_variate => control_variate
       self%i_weight = 3
       if (present(lindf)) then
          self%lindf = lindf
       end if
    end if

  end subroutine initialize_pic_vm_3d3v_trafo


  !> Constructor.
  subroutine initialize_file_pic_vm_3d3v_trafo(&
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
       force_sign, &
       rhob, &
       control_variate, &
       jmean, &
       lindf) 
    class(sll_t_time_propagator_pic_vm_3d3v_trafo_helper), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_3d_base), target,          intent( in ) :: maxwell_solver !< Maxwell solver
    class(sll_c_particle_mesh_coupling_3d), target, intent(in) :: particle_mesh_coupling !< Particle mesh coupling
    class(sll_t_particle_array), target,           intent( in ) :: particle_group !< Particle group
    sll_real64, target,                            intent( in ) :: phi_dofs(:) !< array for the coefficients of the scalar potential 
    sll_real64, target,                            intent( in ) :: efield_dofs(:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent( in ) :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                    intent( in ) :: x_min(3) !< Lower bound of x domain
    sll_real64,                                    intent( in ) :: Lx(3) !< Length of the domain in x direction.
    type(sll_t_mapping_3d), target,                intent( inout ) :: map !< Coordinate transformation
    character(len=*),                              intent( in ) :: filename !< filename
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
       betar_Set = 1._f64
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

    if( present( control_variate) ) then
       ! Read in solver tolerance
       open(newunit = input_file, file=filename, status='old',IOStat=io_stat)
       if (io_stat /= 0) then
          if (rank == 0 ) then
             print*, 'sll_m_time_propagator_pic_vm_3d3v_trafo_helper: Input file does not exist. Set default tolerance.'
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'solver tolerance:', 1d-12
             write(file_id, *) 'iter_tolerance:', 1d-12
             write(file_id, *) 'max_iter:', 10
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
               map, &
               boundary_particles_set, &
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
          if (io_stat /= 0 .and. io_stat1 /= 0) then
             if (rank == 0 ) then
                print*, 'sll_m_time_propagator_pic_vm_3d3v_trafo_helper: Input parameter does not exist. Set default tolerance.'
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'solver tolerance:', 1d-12
                write(file_id, *) 'iter_tolerance:', 1d-12
                write(file_id, *) 'max_iter:', 10
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
                  map, &
                  boundary_particles_set, &
                  betar=betar_set, &
                  force_sign=force_sign_set, &
                  rhob = rhob, &
                  control_variate = control_variate,&
                  jmean=jmean_set, &
                  lindf = lindf_set)
          else if (io_stat == 0 .and. io_stat1 /= 0) then
             if (rank == 0 ) then
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'solver tolerance:', maxwell_tolerance
                write(file_id, *) 'iter_tolerance:', 1d-12
                write(file_id, *) 'max_iter:', 10
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
                  map, &
                  boundary_particles_set, &
                  maxwell_tolerance, &
                  betar=betar_set, &
                  force_sign=force_sign_set, &
                  rhob = rhob, &
                  control_variate = control_variate, &
                  lindf = lindf_set)
          else if (io_stat /= 0 .and. io_stat1 == 0) then
             if (rank == 0 ) then
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'solver tolerance:', 1d-12
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
                  map, &
                  boundary_particles_set, &
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
                  map, &
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
             print*, 'sll_m_time_propagator_pic_vm_3d3v_trafo_helper: Input file does not exist. Set default tolerance.'
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'solver tolerance:', 1d-12
             write(file_id, *) 'iter_tolerance:', 1d-12
             write(file_id, *) 'max_iter:', 10
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
               map, &
               boundary_particles_set, &
               betar=betar_set, &
               force_sign=force_sign_set,&
               rhob = rhob, &
               jmean=jmean_set)
       else
          read(input_file, output,IOStat=io_stat0)
          read(input_file, time_solver,IOStat=io_stat)
          read(input_file, time_iterate,IOStat=io_stat1)
          if (io_stat /= 0 .and. io_stat1 /= 0) then
             if (rank == 0 ) then
                print*, 'sll_m_time_propagator_pic_vm_3d3v_trafo_helper: Input parameter does not exist. Set default tolerance.'
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'solver tolerance:', 1d-12
                write(file_id, *) 'iter_tolerance:', 1d-12
                write(file_id, *) 'max_iter:', 10
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
                  map, &
                  boundary_particles_set, &
                  betar=betar_set, &
                  force_sign=force_sign_set,&
                  rhob = rhob, &
                  jmean=jmean_set)
          else if (io_stat == 0 .and. io_stat1 /= 0) then
             if (rank == 0 ) then
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'solver tolerance:', maxwell_tolerance
                write(file_id, *) 'iter_tolerance:', 1d-12
                write(file_id, *) 'max_iter:', 10
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
                  map, &
                  boundary_particles_set, &
                  maxwell_tolerance, &
                  betar=betar_set, &
                  force_sign=force_sign_set,&
                  rhob = rhob, &
                  jmean=jmean_set)
          else if (io_stat /= 0 .and. io_stat1 == 0) then
             if (rank == 0 ) then
                open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
                write(file_id, *) 'solver tolerance:', 1d-12
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
                  map, &
                  boundary_particles_set, &
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
                  map, &
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


  end subroutine initialize_file_pic_vm_3d3v_trafo


  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_3d3v_trafo(self)
    class(sll_t_time_propagator_pic_vm_3d3v_trafo_helper), intent( inout ) :: self !< time propagator object
    !local variables
    sll_int32 :: file_id, j

    if(self%boundary)then
       print*, 'left boundary', self%counter_left
       print*, 'right boundary', self%counter_right
    end if
    if(self%adiabatic_electrons) then
       call self%linear_solver_schur_phiv%free()
       call self%linear_op_schur_phiv%free()
    else
       call self%preconditioner_fft%free()
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
    if(self%map%flag2d)then
       call self%particle_mass_12%free()
       call self%particle_mass_21%free()
       deallocate( self%particle_mass_12 )
       deallocate( self%particle_mass_21 )
       if(self%map%flag3d)then
          call self%particle_mass_23%free()
          call self%particle_mass_32%free()
          call self%particle_mass_31%free()
          call self%particle_mass_13%free()
          deallocate( self%particle_mass_13 )
          deallocate( self%particle_mass_31 )
          deallocate( self%particle_mass_23 )
          deallocate( self%particle_mass_32 )
       end if
    end if

    deallocate( self%j_dofs )
    deallocate( self%j_dofs_local )
    deallocate( self%particle_mass_1_local )
    deallocate( self%particle_mass_2_local )
    deallocate( self%particle_mass_3_local )
    deallocate( self%particle_mass_12_local )
    deallocate( self%particle_mass_23_local )
    deallocate( self%particle_mass_31_local )

    self%maxwell_solver => null()
    self%particle_mesh_coupling => null()
    self%particle_group => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()
    self%map => null()

    print*, 'No. of failed iterations:', self%n_failed
    open(newunit=file_id, file="outer-iters.dat", action='write')
    do j=1,self%iter_counter
       write(file_id,*) self%niter(j)
    end do

  end subroutine delete_pic_vm_3d3v_trafo


end module sll_m_time_propagator_pic_vm_3d3v_trafo_helper
