!> @ingroup pic_time_integration
!> @author Benedikt Perse, IPP
!> @brief Particle pusher based on antisymmetric splitting with AVF for 1d2v Vlasov-Poisson with coordinate transformation.
!> @details MPI parallelization by domain cloning. Periodic boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_time_propagator_pic_vm_1d2v_trafo_helper
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
       sll_f_get_collective_rank, &
       sll_o_collective_allreduce, &
       sll_v_world_collective

  use sll_m_gauss_legendre_integration, only : &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_maxwell_1d_base, only: &
       sll_c_maxwell_1d_base

  use sll_mpi, only: &
       mpi_sum, &
       mpi_max

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use sll_m_particle_mass_1d_base, only: &
       sll_c_particle_mass_1d_base

  use sll_m_particle_mesh_coupling_base_1d, only: &
       sll_c_particle_mesh_coupling_1d

  use sll_m_linear_operator_particle_mass_1d, only : &
       sll_t_linear_operator_particle_mass_1d

  use sll_m_linear_operator_particle_mass_cl_1d, only : &
       sll_t_linear_operator_particle_mass_cl_1d

  use sll_m_linear_operator_schur_ev_1d, only : &
       sll_t_linear_operator_schur_ev_1d

  use sll_m_linear_solver_cg, only : &
       sll_t_linear_solver_cg

  use sll_m_spline_fem_utilities, only : &
       sll_s_spline_fem_mass_line

  use sll_m_spline_fem_utilities_sparse, only : &
       sll_s_spline_fem_sparsity_mass

  use sll_m_time_propagator_pic_vm_3d3v_cl_helper, only: &
       sll_p_boundary_particles_periodic, &
       sll_p_boundary_particles_singular, &
       sll_p_boundary_particles_reflection, &
       sll_p_boundary_particles_absorption


  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_1d2v_trafo_helper

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type :: sll_t_time_propagator_pic_vm_1d2v_trafo_helper
     class(sll_c_maxwell_1d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_0  !< Kernel smoother (order p+1)
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_1  !< Kernel smoother (order p)
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

     class( sll_c_particle_mass_1d_base ), allocatable :: particle_mass_1 !< Particle mass
     class( sll_c_particle_mass_1d_base ), allocatable :: particle_mass_0 !< Particle mass
     type( sll_t_linear_operator_schur_ev_1d ) :: linear_operator_1 !< Schur operator for advect_ev 
     type(sll_t_linear_solver_cg)    :: linear_solver_1 !< Linear solver for Schur operator for advect_ev 
     type( sll_t_linear_operator_schur_ev_1d ) :: linear_operator_0 !< Schur operator for advect_ev 
     type(sll_t_linear_solver_cg)    :: linear_solver_0 !< Linear solver for Schur operator for advect_ev 

     type(sll_t_mapping_3d), pointer      :: map !< Coordinate transformation

     sll_int32 :: spline_degree !< Spline degree
     sll_real64 :: Lx !< Size of the domain
     sll_real64 :: x_min !< Lower bound for x domain
     sll_int32  :: n_cells !< number of grid cells
     sll_int32  :: n_dofs0  !< number of Dofs for 0form
     sll_int32  :: n_dofs1 !< number of Dofs for 1form

     sll_real64, pointer     :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
     sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the magnetic field
     sll_real64, allocatable :: j_dofs(:,:)      !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:,:)!< MPI-processor local part of one component of \a j_dofs
     sll_real64, allocatable :: particle_mass_0_local(:,:) !< Array to hold the 2*spline_degree+1 diagonals of the matrix A_0 M_p A_0^T
     sll_real64, allocatable :: particle_mass_1_local(:,:) !< Array to hold the 2*(spline_degree-1)+1 diagonals of the matrix A_1 M_p A_1^T
     
     sll_int32 :: n_particles_max !< maximal number of particles
     sll_real64, allocatable :: xnew(:,:) !< extra data for particle position
     sll_real64, allocatable :: vnew(:,:,:) !< extra data for particle velocity
     sll_real64, allocatable :: efield_dofs_new(:,:) !< extra data for efield DoFs
     sll_real64, allocatable :: efield_dofs_work(:,:) !< extra data for efield DoFs
     
     sll_int32  :: boundary_particles = 100 !< particle boundary conditions
     sll_int32 :: counter_left = 0 !< boundary counter
     sll_int32 :: counter_right = 0 !< boundary counter
     logical :: boundary = .false. !< true for non periodic boundary conditions
     sll_real64, pointer     :: rhob(:) => null() !< charge at the boundary
     sll_real64 :: force_sign = 1._f64 !< sign of particle force
     logical    :: jmean = .false. !< logical for mean value of current

     sll_real64 :: solver_tolerance !< solver tolerance


     sll_int32 :: max_iter = 1 !< maximal amount of iterations, set 1 for ctest
     sll_real64 :: iter_tolerance = 1d-10 !< iteration tolerance

     sll_int32 :: n_failed = 0 !< number of failed iterations
     sll_int32 :: iter_counter = 0 !< number of iterations
     sll_int32 :: niter(100000) !< list of iteration numbers

   contains
     procedure :: advect_x => advect_x_pic_vm_1d2v_trafo !> Advect x-part
     procedure :: advect_vb => advect_vb_pic_vm_1d2v_trafo !> Push v, vxB-part only
     procedure :: advect_eb => advect_eb_pic_vm_1d2v_trafo !> Solve Faraday and B-part of Ampere
     procedure :: advect_e => advect_e_pic_vm_1d2v_trafo !> Advect ev-part 
     procedure :: advect_ex => advect_ex_pic_vm_1d2v_trafo !> Advect ev-part together with x-part in nonlinear iteration
     
     procedure :: init => initialize_pic_vm_1d2v_trafo !> Initialize the type
     procedure :: init_from_file => initialize_file_pic_vm_1d2v_trafo !> Initialize from nml file
     procedure :: free => delete_pic_vm_1d2v_trafo !> Finalization

  end type sll_t_time_propagator_pic_vm_1d2v_trafo_helper

contains
  

  !---------------------------------------------------------------------------!
  !> advect_x: Equations to be solved
  !> $\Xi^{n+1}=\Xi^n+ \frac{\Delta t}{2} (DF^{-1}(\Xi^{n+1})+DF^{-1}(\Xi^n)) V^n$
  subroutine advect_x_pic_vm_1d2v_trafo ( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_trafo_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: dt   !< time step
    !local variables
    sll_int32 :: i_part, i_sp, i
    sll_real64 :: xi(3), xnew(3), vi(3), vt, matrix(3,3), xs
    sll_real64 :: err

    do i_sp=1,self%particle_group%n_species
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)

          matrix=self%map%jacobian_matrix_inverse( [xi(1), 0._f64, 0._f64] )
          !Transformation of the v_1 coordinate 
          vt = matrix(1,1) * vi(1)  
          !x^\star=\mod(x^n+dt*DF^{-1}(x^n)v^n,1)
          xs = xi(1) + dt * vt
          xnew = xi
          xnew(1) = xs
          i = 0
          err= abs(xi(1) - xs)
          do while(i < self%max_iter .and. err > self%iter_tolerance)
             call compute_particle_boundary( self, xi(1), xs, vi(2) )
             matrix=self%map%jacobian_matrix_inverse( [xs, 0._f64, 0._f64] )
             xs = xi(1) + 0.5_f64 * dt * (vt + matrix(1,1) * vi(1) )
             err = abs(xs - xnew(1))
             xnew(1) = xs
             i = i + 1
          end do
          call compute_particle_boundary( self, xi(1), xnew(1), vi(1) )

          call self%particle_group%group(i_sp)%set_x ( i_part, xnew )
          call self%particle_group%group(i_sp)%set_v ( i_part, vi )
       end do
    end do

  end subroutine advect_x_pic_vm_1d2v_trafo


  !> Helper function for \a advect_x
  subroutine compute_particle_boundary( self, xold, xnew, vi  )
    class(sll_t_time_propagator_pic_vm_1d2v_trafo_helper), intent( inout ) :: self  
    sll_real64,                                     intent( inout ) :: xold
    sll_real64,                                     intent( inout ) :: xnew
    sll_real64,                                     intent( inout ) :: vi
    !local variables
    sll_real64 :: xmid, xbar, dx

    if(xnew < 0._f64 .or. xnew > 1._f64 )then
       if(xnew < 0._f64  )then
          xbar = 0._f64
          self%counter_left = self%counter_left+1
       else if(xnew > 1._f64)then
          xbar = 1._f64
          self%counter_right = self%counter_right+1
       end if
       dx = (xbar- xold)/(xnew-xold)
       xmid = xold + dx * (xnew-xold)
       xmid = xbar

       select case(self%boundary_particles)
       case(sll_p_boundary_particles_reflection)
          vi = -vi
          xnew = 2._f64*xbar-xnew
       case(sll_p_boundary_particles_absorption)
       case( sll_p_boundary_particles_periodic)
          xnew = modulo(xnew, 1._f64)
       case default
          xnew = modulo(xnew, 1._f64)
       end select
    end if

  end subroutine compute_particle_boundary


  !---------------------------------------------------------------------------!
  !> advect_vb: Equations to be solved
  !> $V^{n+1}= (\cos(DF/J_F B)&\sin(DF/J_F B) \\ -\sin(DF/J_F B) &\cos(DF/J_F B) ) V^n
  subroutine advect_vb_pic_vm_1d2v_trafo ( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_trafo_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: dt   !< time step
    !local variables
    sll_int32  :: i_part, i_sp
    sll_real64 :: qmdt
    sll_real64 :: vi(3), vnew(3), xi(3), factor(3,3)
    sll_real64 :: bfield, cs, sn

    do i_sp=1,self%particle_group%n_species
       qmdt = self%particle_group%group(i_sp)%species%q_over_m()*dt;

       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)

          call self%kernel_smoother_1%evaluate(xi(1), self%bfield_dofs, bfield)

          factor=self%map%jacobian_matrix( [xi(1), 0._f64, 0._f64])/self%map%jacobian( [xi(1), 0._f64, 0._f64])*qmdt

          cs = cos(factor(3,3)*bfield)
          sn = sin(factor(3,3)*bfield)

          vnew(1) = cs * vi(1) + sn * vi(2)
          vnew(2) = -sn* vi(1) + cs * vi(2)

          call self%particle_group%group(i_sp)%set_v( i_part, vnew )
       end do
    end do

  end subroutine advect_vb_pic_vm_1d2v_trafo


  !---------------------------------------------------------------------------!
  !> advect_eb: Equations to be solved
  !> Solution with Schur complement: $ S=M_1+\frac{\Delta t^2}{4} D^\top M_2 D $
  !> $ e_2^{n+1}=S^{-1}( (M_1-\frac{\Delta t^2}{4} D^\top M_2 D)e_2^n+\Delta t D^\top M_2 b_3^n) $
  !> $ b_3^{n+1}=b_3^n-\frac{\Delta t}{2} C(e_2^n+e_2^{n+1}) $
  subroutine advect_eb_pic_vm_1d2v_trafo ( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_trafo_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: dt   !< time step

    call self%maxwell_solver%compute_curl_part( dt, self%efield_dofs(:,2), self%bfield_dofs )


  end subroutine advect_eb_pic_vm_1d2v_trafo


  !---------------------------------------------------------------------------!
  !> advect_e: Equations to be solved
  !> Solution with Schur complement: $ S_{+}=M_1+\frac{\Delta t^2 q^2}{4 m} (\mathbb{\Lambda}^1)^T DF^{-1} DF^{-T} \mathbb{\Lambda}^1 $
  !> $e^{n+1}=S_{+}^{-1}\left(S_{-}e^n-\Delta t (\mathbb{\Lambda}^1)^\top DF^{-1}\mathbb{W}_q V^n \right)$
  !> $V^{n+1}=V^n+\frac{\Delta t}{2} \mathbb{W}_{\frac{q}{m}} DF^{-\top} \mathbb{\Lambda}^1(e^{n+1}+e^n)$
  subroutine advect_e_pic_vm_1d2v_trafo ( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_trafo_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: dt   !< time step
    ! local variables
    sll_int32 :: i_part, i_sp, j
    sll_real64 :: vi(3), xi(3), wi(1)
    sll_real64 :: qoverm, jmat(3,3), metric, factor
    sll_real64 :: efield(2)
    sll_real64 :: rhs0(self%n_dofs0)
    sll_real64 :: rhs1(self%n_dofs1)

    ! Set to zero
    self%j_dofs_local = 0.0_f64
    self%particle_mass_1_local = 0.0_f64
    self%particle_mass_0_local = 0.0_f64

    ! First particle loop
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       factor = dt**2 * 0.25_f64* qoverm
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)

          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part)

          metric= self%map%metric_inverse_single( [xi(1), 0._f64, 0._f64], 1, 1 )
          jmat=self%map%jacobian_matrix_inverse( [xi(1), 0._f64, 0._f64] )

          ! Accumulate the particle mass matrix diagonals
          call self%kernel_smoother_1%add_particle_mass_full( xi(1), wi(1) * metric * factor, &
               self%particle_mass_1_local )
          call self%kernel_smoother_0%add_particle_mass_full( xi(1), wi(1) * factor , &
               self%particle_mass_0_local )

          ! Accumulate jx
          call self%kernel_smoother_1%add_charge( xi(1), wi(1)*jmat(1,1)*vi(1), &
               self%j_dofs_local(1:self%n_dofs1,1) )
          ! Accumulate jy
          call self%kernel_smoother_0%add_charge( xi(1), wi(1)*vi(2), &
               self%j_dofs_local(:,2) )

          ! Evaulate E_x and E_y at particle position and propagate v a half step
          call self%kernel_smoother_1%evaluate &
               ( xi(1), self%efield_dofs(1:self%n_dofs1,1), efield(1) )
          call self%kernel_smoother_0%evaluate &
               ( xi(1), self%efield_dofs(:,2), efield(2))

          ! velocity update with jacobian matrix, matrix is not transposed
          do j= 1, 2
             vi(j) = vi(j) + dt* 0.5_f64* qoverm *(jmat(1,j)* efield(1)+jmat(2,j)* efield(2))
          end do

          call self%particle_group%group(i_sp)%set_v( i_part, vi )
       end do
    end do

    self%j_dofs = 0.0_f64
    self%particle_mass_0%particle_mass = 0.0_f64
    self%particle_mass_1%particle_mass = 0.0_f64  

    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(1:self%n_dofs1,1), &
         self%n_dofs1, MPI_SUM, self%j_dofs(1:self%n_dofs1,1) )
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
         self%n_dofs0, MPI_SUM, self%j_dofs(:,2) )

    call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_1_local, &
         self%n_dofs1*(2*self%spline_degree-1), MPI_SUM, self%particle_mass_1%particle_mass)
    call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_0_local, &
         self%n_dofs0*(2*self%spline_degree+1), MPI_SUM, self%particle_mass_0%particle_mass)


    if ( self%jmean ) then
       self%j_dofs(:,1) = self%j_dofs(:,1) - sum(self%j_dofs(:,1))/real(self%n_dofs1, f64)
    end if

    ! Compute rhs
    self%linear_operator_1%sign = -self%force_sign 
    call self%linear_operator_1%dot(self%efield_dofs(1:self%n_dofs1,1) , rhs1 )
    rhs1 = rhs1 - dt * self%force_sign * self%j_dofs(1:self%n_dofs1,1)
    ! Solve the Schur complement
    self%linear_operator_1%sign = self%force_sign 
    call self%linear_solver_1%set_guess(self%efield_dofs(1:self%n_dofs1,1))
    call self%linear_solver_1%solve( rhs1 , self%efield_dofs(1:self%n_dofs1,1) )

    ! Compute rhs
    self%linear_operator_0%sign = -self%force_sign 
    call self%linear_operator_0%dot(self%efield_dofs(:,2) , rhs0 )
    rhs0 = rhs0 - self%force_sign * dt * self%j_dofs(:,2)

    ! Solve the Schur complement
    self%linear_operator_0%sign = self%force_sign 
    call self%linear_solver_0%set_guess(self%efield_dofs(:,2))
    call self%linear_solver_0%solve( rhs0 , self%efield_dofs(:,2) )

    if( self%boundary ) then!perfect conductor boundary
       self%efield_dofs(1,2) = 0._f64
       self%efield_dofs(self%n_dofs0,2) = 0._f64
    end if

    ! Second particle loop (second half step of particle propagation)
    do i_sp = 1, self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)

          ! Evaulate E_x and E_y at particle position and propagate v a half step
          call self%kernel_smoother_1%evaluate &
               ( xi(1), self%efield_dofs(1:self%n_dofs1,1), efield(1))
          call self%kernel_smoother_0%evaluate &
               ( xi(1), self%efield_dofs(:,2), efield(2))

          jmat = self%map%jacobian_matrix_inverse_transposed( [xi(1), 0._f64, 0._f64] )

          do j = 1, 2
             vi(j) = vi(j) + dt* 0.5_f64* qoverm *(jmat(j,1)*efield(1) + jmat(j,2)*efield(2))
          end do

          call self%particle_group%group(i_sp)%set_v( i_part, vi )
       end do
    end do

  end subroutine advect_e_pic_vm_1d2v_trafo


  !---------------------------------------------------------------------------!
  !> advect_ex: Equations to be solved
  !> $\frac{\Xi^{n+1}-\Xi^n}{\Delta t}=\frac{DF^{-1}(\Xi^{n+1})+DF^{-1}(\Xi^n)}{2} \frac{V^{n+1}+V^n}{2}$
  !> $\frac{V^{n+1}-V^n}{\Delta t}=\mathbb{W}_{\frac{q}{m}} \frac{DF^{-\top}(\Xi^{n+1})+DF^{-\top}(\Xi^n)}{2} \frac{1}{\Delta t}\int_{t^n}^{t^{n+1}} \mathbb{\Lambda}^1(\Xi(\tau))  d\tau \frac{e^{n+1}+e^n}{2}$
  !> $\frac{M_1 e^{n+1}-M_1 e^n}{\Delta t} = - \frac{1}{\Delta t} \int_{t^n}^{t^{n+1}} \mathbb{\Lambda}^1(\Xi(\tau))^\top  d\tau  \frac{DF^{-1}(\Xi^{n+1})+DF^{-1}(\Xi^n)}{2} \mathbb{W}_q\frac{V^{n+1}+V^n}{2}$
  subroutine advect_ex_pic_vm_1d2v_trafo ( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_trafo_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: dt   !< time step
    ! local variables
    sll_int32 :: i_part, i_sp, j
    sll_real64 :: vi(3), vh, xi(3), wi(1), xnew(3), vbar(2)
    sll_real64 :: qoverm
    sll_real64 :: efield(3), jmat(3,3), jmatrix(3,3)
    sll_int32 :: niter
    sll_real64 :: residual(1), residual_local(1)

    self%efield_dofs_new = self%efield_dofs
    do i_sp=1,self%particle_group%n_species
       do i_part = 1,self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xnew = self%particle_group%group(i_sp)%get_x(i_part)
          self%vnew(i_sp,:, i_part) = vi(1:2)
          self%xnew(i_sp, i_part) = xnew(1)
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
             wi = self%particle_group%group(i_sp)%get_charge(i_part)
             vbar = 0.5_f64 * (self%vnew(i_sp,:, i_part)+vi(1:2))
             xnew(1) = self%xnew(i_sp, i_part)

             jmat = self%map%jacobian_matrix_inverse_transposed( [xi(1), 0._f64, 0._f64] )
             jmatrix=self%map%jacobian_matrix_inverse_transposed( [xnew(1), 0._f64, 0._f64] )
             vh = 0.5_f64 * (jmatrix(1,1)+jmat(1,1))*vbar(1) 
             xnew(1) = xi(1) + dt * vh
   
             call compute_particle_boundary_current_evaluate( self, xi, xnew, vi, vbar, wi, qoverm, dt )
             
             self%xnew(i_sp,i_part) = xnew(1)
             self%vnew(i_sp,:,i_part) = vi(1:2)
          end do
       end do

       self%j_dofs = 0.0_f64
       ! MPI to sum up contributions from each processor
       call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
            self%n_dofs1, MPI_SUM, self%j_dofs(1:self%n_dofs1,1) )
       call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
            self%n_dofs0, MPI_SUM, self%j_dofs(:,2) )

       self%efield_dofs_work = self%efield_dofs
       call self%maxwell_solver%compute_E_from_j(self%j_dofs(1:self%n_dofs1,1), 1, self%efield_dofs_work(1:self%n_dofs1,1) )
       call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,2), 2, self%efield_dofs_work(:,2) )

       ! Compute residual based on e
       residual_local = (sum((self%efield_dofs_work-self%efield_dofs_new)**2))*self%maxwell_solver%delta_x
       call sll_o_collective_allreduce( sll_v_world_collective, residual_local, 1, MPI_MAX, residual )
       residual = sqrt(residual)
       self%efield_dofs_new = self%efield_dofs_work
    end do

    if ( residual(1) > self%iter_tolerance ) then
       print*, 'Warning: Iteration no.', self%iter_counter+1 ,'did not converge.', residual, niter
       self%n_failed = self%n_failed+1
    end if

      self%efield_dofs = self%efield_dofs_new
      do i_sp = 1, self%particle_group%n_species
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi(1:2) = self%vnew(i_sp,:,i_part)
          xnew(1) = self%xnew(i_sp,i_part)
          call self%particle_group%group(i_sp)%set_v( i_part, vi )
          call self%particle_group%group(i_sp)%set_x( i_part, xnew )
       end do
    end do

    self%iter_counter = self%iter_counter + 1
    self%niter(self%iter_counter) = niter
    
  end subroutine advect_ex_pic_vm_1d2v_trafo
  
  
  !> Helper function for \a advect_ex
  subroutine compute_particle_boundary_current_evaluate( self, xi, xnew, vi, vbar, wi, qoverm, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_trafo_helper), intent( inout ) :: self !< time splitting object 
    sll_real64,                                           intent( in    ) :: xi(1)
    sll_real64,                                           intent( inout ) :: xnew(1)
    sll_real64,                                           intent( inout ) :: vi(2)
    sll_real64,                                           intent( inout ) :: vbar(2)
    sll_real64,                                           intent( in    ) :: wi(1)
    sll_real64,                                           intent( in    ) :: qoverm
    sll_real64,                                           intent( in    ) :: dt
    !local variables
    sll_real64 :: xmid(1), xbar
    sll_real64 :: jmatrix(3,3), jmat(3,3), efield(2)
    sll_int32 :: j


    jmat = self%map%jacobian_matrix_inverse_transposed( [xi(1), 0._f64, 0._f64] )

    if(xnew(1) < 0._f64 .or. xnew(1) > 1._f64 )then
       if(xnew(1) < 0._f64  )then
          xbar = 0._f64
          self%counter_left = self%counter_left+1
       else if(xnew(1) > 1._f64)then
          xbar = 1._f64
          self%counter_right = self%counter_right+1
       end if

       xmid = xbar

       call self%kernel_smoother_1%add_current_evaluate &
            ( xi(1), xmid(1), wi(1), vbar(1), self%efield_dofs_work(1:self%n_dofs1,1), self%j_dofs_local(1:self%n_dofs1,1), &
            efield(1) )

       call self%kernel_smoother_0%add_current_evaluate &
            ( xi(1), xmid(1), wi(1)*vbar(2)/vbar(1), vbar(1), &
            self%efield_dofs_work(:,2), self%j_dofs_local(:,2), &
            efield(2) )

       jmatrix = self%map%jacobian_matrix_inverse_transposed( [xmid(1), 0._f64, 0._f64] )
       do j = 1, 2
          vi(j) = vi(j) + qoverm *0.5_f64*((jmatrix(j,1)+jmat(j,1))*efield(1) + (jmatrix(j,2)+jmat(j,2))*efield(2))
       end do

       select case(self%boundary_particles)
       case(sll_p_boundary_particles_reflection) 
          xnew = 2._f64*xbar-xnew
          vi(1) = -vi(1)
          vbar(1) = -vbar(1) !?
       case(sll_p_boundary_particles_absorption)
          call self%kernel_smoother_0%add_charge(xmid, wi(1), self%rhob)
       case( sll_p_boundary_particles_periodic)
          call self%kernel_smoother_0%add_charge(xmid, wi(1), self%rhob)
          xnew = modulo(xnew, 1._f64)
          xmid = 1._f64-xbar
          call self%kernel_smoother_0%add_charge(xmid, -wi(1), self%rhob)
       case default
          !call self%kernel_smoother_0%add_charge(xmid, wi(1), self%rhob)
          xnew = modulo(xnew, 1._f64)
          xmid = 1._f64-xbar
          call self%kernel_smoother_0%add_charge(xmid, -wi(1), self%rhob)
       end select

       call self%kernel_smoother_1%add_current_evaluate &
            ( xmid(1), xnew(1), wi(1), vbar(1), self%efield_dofs_work(1:self%n_dofs1,1), self%j_dofs_local(1:self%n_dofs1,1), &
            efield(1) )

       call self%kernel_smoother_0%add_current_evaluate &
            ( xmid(1), xnew(1), wi(1)*vbar(2)/vbar(1), vbar(1), &
            self%efield_dofs_work(:,2), self%j_dofs_local(:,2), &
            efield(2) )

    else
       if ( abs(vbar(1)) > 1.0D-16 ) then
          call self%kernel_smoother_1%add_current_evaluate &
               ( xi(1), xnew(1), wi(1), vbar(1), self%efield_dofs_work(1:self%n_dofs1,1), self%j_dofs_local(1:self%n_dofs1,1), &
               efield(1) )

          call self%kernel_smoother_0%add_current_evaluate &
               ( xi(1), xnew(1), wi(1)*vbar(2)/vbar(1), vbar(1), &
               self%efield_dofs_work(:,2), self%j_dofs_local(:,2), &
               efield(2) )
       else
          call self%kernel_smoother_1%evaluate &
               (xi(1), self%efield_dofs_work(1:self%n_dofs1,1), efield(1) )
          efield(1) = efield(1)*dt
          call self%kernel_smoother_0%add_charge( xi(1), &
               wi(1)* vbar(2)*dt, self%j_dofs_local(:,2) )
          call self%kernel_smoother_0%evaluate &
               (xi(1), self%efield_dofs_work(:,2), efield(2) )
          efield(2) = efield(2)*dt
       end if
    end if
    jmatrix = self%map%jacobian_matrix_inverse_transposed( [xnew(1), 0._f64, 0._f64] )
    do j = 1, 2
       vi(j) = vi(j) + qoverm *0.5_f64*((jmatrix(j,1)+jmat(j,1))*efield(1) + (jmatrix(j,2)+jmat(j,2))*efield(2))
    end do

  end subroutine compute_particle_boundary_current_evaluate
    
  !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_1d2v_trafo(&
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
       boundary_particles, &
       solver_tolerance, &
       iter_tolerance, max_iter, &
       force_sign, &
       rhob, &
       jmean)  
    class(sll_t_time_propagator_pic_vm_1d2v_trafo_helper), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_1d_base), target,                 intent( in )  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent( in )  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent( in )  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,                  intent( in )  :: particle_group !< Particle group
    sll_real64, target,                                   intent( in )  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                                   intent( in )  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                           intent( in )  :: x_min !< Lower bound of x domain
    sll_real64,                                           intent( in )  :: Lx !< Length of the domain in x direction.
    type(sll_t_mapping_3d), target,                       intent( in )  :: map !< Coordinate transformation
    sll_int32, optional,                                  intent( in )  :: boundary_particles !< particle boundary conditions
    sll_real64, optional,                                 intent( in )  :: solver_tolerance !< solver tolerance
    sll_real64, optional,                          intent( in ) :: iter_tolerance !< iteration tolerance
    sll_int32,  optional,                          intent( in ) :: max_iter !< maximal number of iterations
    sll_real64, optional,                                 intent( in )  :: force_sign !< sign of particle force
    sll_real64, optional, target,                  intent( in ) :: rhob(:) !< charge at the boundary
    logical, optional,                                    intent( in )  :: jmean !< logical for mean value of current
    !local variables
    sll_int32 :: ierr, j


    if (present(solver_tolerance) )  then
       self%solver_tolerance = solver_tolerance
    else
       self%solver_tolerance = 1d-12
    end if

    if( present(force_sign) )then
       self%force_sign = force_sign
    end if

    if (present(jmean)) then
       self%jmean = jmean
    end if

    if (present(iter_tolerance) )  then
       self%iter_tolerance = iter_tolerance
       self%max_iter = max_iter
    else
       self%iter_tolerance = 1d-10
       self%max_iter = 10
    end if

    self%maxwell_solver => maxwell_solver
    self%kernel_smoother_0 => kernel_smoother_0
    self%kernel_smoother_1 => kernel_smoother_1
    self%particle_group => particle_group
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs
    self%map => map

    ! Check that n_dofs is the same for both kernel smoothers.
    SLL_ASSERT( self%kernel_smoother_0%n_cells == self%kernel_smoother_1%n_cells )
    SLL_ASSERT( self%kernel_smoother_0%n_cells == self%maxwell_solver%n_cells )


    self%n_cells = self%maxwell_solver%n_cells
    self%n_dofs0 = self%maxwell_solver%n_dofs0
    self%n_dofs1 = self%maxwell_solver%n_dofs1
    self%spline_degree = self%kernel_smoother_0%spline_degree
 
    SLL_ALLOCATE( self%j_dofs(self%n_dofs0,2), ierr )
    SLL_ALLOCATE( self%j_dofs_local(self%n_dofs0,2), ierr )
    SLL_ALLOCATE( self%particle_mass_1_local(2*self%spline_degree-1, self%n_dofs1), ierr )
    SLL_ALLOCATE( self%particle_mass_0_local(2*self%spline_degree+1, self%n_dofs0), ierr )
    self%j_dofs = 0.0_f64
    self%j_dofs_local = 0.0_f64
    self%particle_mass_1_local = 0.0_f64
    self%particle_mass_0_local = 0.0_f64

    if( self%n_cells+self%spline_degree == self%maxwell_solver%n_dofs0   ) then
       self%boundary = .true.
       self%boundary_particles = boundary_particles
       allocate( sll_t_linear_operator_particle_mass_cl_1d  :: self%particle_mass_1 )
       select type ( q => self%particle_mass_1 )
       type is ( sll_t_linear_operator_particle_mass_cl_1d )
          call q%create( self%spline_degree-1, self%n_dofs1 )
       end select

       allocate( sll_t_linear_operator_particle_mass_cl_1d  :: self%particle_mass_0 )
       select type ( q => self%particle_mass_0 )
       type is ( sll_t_linear_operator_particle_mass_cl_1d )
          call q%create( self%spline_degree, self%n_dofs0 )
       end select
    else if ( self%n_cells == self%maxwell_solver%n_dofs0 ) then
       allocate( sll_t_linear_operator_particle_mass_1d  :: self%particle_mass_1 )
       select type ( q => self%particle_mass_1 )
       type is ( sll_t_linear_operator_particle_mass_1d )
          call q%create( self%spline_degree-1, self%n_dofs1 )
       end select

       allocate( sll_t_linear_operator_particle_mass_1d  :: self%particle_mass_0 )
       select type ( q => self%particle_mass_0 )
       type is ( sll_t_linear_operator_particle_mass_1d )
          call q%create( self%spline_degree, self%n_dofs0 )
       end select
    end if

    call self%linear_operator_1%create( self%maxwell_solver, self%particle_mass_1,  self%n_dofs1, self%maxwell_solver%s_deg_1 )
    call self%linear_solver_1%create( self%linear_operator_1 )
    self%linear_solver_1%atol = self%solver_tolerance/Lx
    !self%linear_solver_1%verbose = .true.

    call self%linear_operator_0%create( self%maxwell_solver, self%particle_mass_0,  self%n_dofs0, self%maxwell_solver%s_deg_0 )
    call self%linear_solver_0%create( self%linear_operator_0 )
    self%linear_solver_0%atol = self%solver_tolerance
    !self%linear_solver_0%verbose = .true.

    self%n_particles_max = self%particle_group%group(1)%n_particles
    do j = 2,self%particle_group%n_species       
       self%n_particles_max = max(self%n_particles_max, self%particle_group%group(j)%n_particles )
    end do
    
    SLL_ALLOCATE( self%xnew(self%particle_group%n_species,self%n_particles_max), ierr )
    SLL_ALLOCATE( self%vnew(self%particle_group%n_species,2,self%n_particles_max), ierr )
    SLL_ALLOCATE( self%efield_dofs_new(self%n_dofs0,2), ierr )
    SLL_ALLOCATE( self%efield_dofs_work(self%n_dofs0,2), ierr )

    if( present(rhob) )then
       self%rhob => rhob
    end if
    
  end subroutine initialize_pic_vm_1d2v_trafo


  !> Constructor.
  subroutine initialize_file_pic_vm_1d2v_trafo(&
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
       filename, &
       boundary_particles, &
       force_sign, &
       rhob, &
       jmean)  
    class(sll_t_time_propagator_pic_vm_1d2v_trafo_helper), intent( out ) :: self !< time propagator object 
    class(sll_c_maxwell_1d_base), target,                 intent( in )  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent( in )  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent( in )  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,                  intent( in )  :: particle_group !< Particle group
    sll_real64, target,                                   intent( in )  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                                   intent( in )  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                           intent( in )  :: x_min !< Lower bound of x domain
    sll_real64,                                           intent( in )  :: Lx !< Length of the domain in x direction.
    type(sll_t_mapping_3d), target,                       intent( in )  :: map  !< Coordinate transformation  
    character(len=*),                                     intent( in )  :: filename !< filename
    sll_int32, optional, intent( in )  :: boundary_particles !< particle boundary conditions
    sll_real64, optional, intent(in) :: force_sign !< sign of particle force
    sll_real64, optional, target,                  intent( in ) :: rhob(:) !< charge at the boundary
    logical, optional, intent(in)    :: jmean !< logical for mean value of current
    !local variables
    character(len=256) :: file_prefix
    sll_int32 :: input_file, rank
    sll_int32 :: io_stat, io_stat0, io_stat1, file_id, boundary_particles_set, max_iter
    sll_real64 :: maxwell_tolerance, force_sign_set, iter_tolerance
    logical :: jmean_set

    namelist /output/ file_prefix
    namelist /time_solver/ maxwell_tolerance
    namelist /time_iterate/ iter_tolerance, max_iter

    rank = sll_f_get_collective_rank(sll_v_world_collective)

    if (present(boundary_particles)) then
       boundary_particles_set = boundary_particles
    else
       boundary_particles_set = 100
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

    ! Read in solver tolerance
    open(newunit = input_file, file=filename, status='old',IOStat=io_stat)
    if (io_stat /= 0) then
       if (rank == 0 ) then
          print*, 'sll_m_time_propagator_pic_vm_1d2v_trafo_helper: Input file does not exist. Set default tolerance.'
          open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
          write(file_id, *) 'solver tolerance:', 1d-12
          write(file_id, *) 'iter_tolerance:', 1d-10
          write(file_id, *) 'max_iter:', 10
          write(file_id, *) 'force_sign:', force_sign_set
          close(file_id)
       end if
       call self%init( maxwell_solver, &
            kernel_smoother_0, &
            kernel_smoother_1, &
            particle_group, &
            efield_dofs, &
            bfield_dofs, &
            x_min, &
            Lx, &
            map,&
            boundary_particles = boundary_particles_set, &
            force_sign=force_sign_set, &
            rhob = rhob, &
            jmean=jmean_set)
    else
       read(input_file, output, IOStat=io_stat0)
       read(input_file, time_solver,IOStat=io_stat)
       read(input_file, time_iterate,IOStat=io_stat1)
       if (io_stat /= 0 .and. io_stat1 /= 0) then
          if (rank == 0 ) then
             print*, 'sll_m_time_propagator_pic_vm_1d2v_trafo_helper: Tolerance not given. Set default tolerance.'
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'solver tolerance:', 1d-12
             write(file_id, *) 'iter_tolerance:', 1d-10
             write(file_id, *) 'max_iter:', 10
             write(file_id, *) 'force_sign:', force_sign_set
             close(file_id)
          end if
          call self%init( maxwell_solver, &
               kernel_smoother_0, &
               kernel_smoother_1, &
               particle_group, &
               efield_dofs, &
               bfield_dofs, &
               x_min, &
               Lx, &
               map,&
               boundary_particles = boundary_particles_set, &
               force_sign=force_sign_set, &
               rhob = rhob, &
               jmean=jmean_set)
       else if (io_stat == 0 .and. io_stat1 /= 0) then
          if (rank == 0 ) then
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'solver tolerance:', maxwell_tolerance
             write(file_id, *) 'iter_tolerance:', 1d-10
             write(file_id, *) 'max_iter:', 10
             write(file_id, *) 'force_sign:', force_sign_set
             close(file_id) 
          end if
          call self%init( maxwell_solver, &
               kernel_smoother_0, &
               kernel_smoother_1, &
               particle_group, &
               efield_dofs, &
               bfield_dofs, &
               x_min, &
               Lx, &
               map, &
               boundary_particles = boundary_particles_set, &
               solver_tolerance=maxwell_tolerance,&
               force_sign=force_sign_set, &
               rhob = rhob, &
               jmean=jmean_set)
       else if (io_stat /= 0 .and. io_stat1 == 0) then
          if (rank == 0 ) then
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'solver tolerance:', 1d-12
             write(file_id, *) 'iter_tolerance:', iter_tolerance
             write(file_id, *) 'max_iter:', max_iter
             write(file_id, *) 'force_sign:', force_sign_set
             close(file_id) 
          end if
           call self%init( maxwell_solver, &
               kernel_smoother_0, &
               kernel_smoother_1, &
               particle_group, &
               efield_dofs, &
               bfield_dofs, &
               x_min, &
               Lx, &
               map, &
               boundary_particles = boundary_particles_set, &
               iter_tolerance=iter_tolerance, max_iter=max_iter, &
               force_sign=force_sign_set, &
               rhob = rhob, &
               jmean=jmean_set)
       else
          if (rank == 0 ) then
             open(newunit=file_id, file=trim(file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=io_stat)
             write(file_id, *) 'solver tolerance:', maxwell_tolerance
             write(file_id, *) 'force_sign:', force_sign_set
             write(file_id, *) 'iter_tolerance:', iter_tolerance
             write(file_id, *) 'max_iter:', max_iter
             close(file_id)
          end if
          call self%init( maxwell_solver, &
               kernel_smoother_0, &
               kernel_smoother_1, &
               particle_group, &
               efield_dofs, &
               bfield_dofs, &
               x_min, &
               Lx, &
               map, &
               boundary_particles = boundary_particles_set, &
               solver_tolerance=maxwell_tolerance,&
               iter_tolerance=iter_tolerance, max_iter=max_iter, &
               force_sign=force_sign_set, &
               rhob = rhob, &
               jmean=jmean_set)
       end if
       close(input_file)
    end if


  end subroutine initialize_file_pic_vm_1d2v_trafo

  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_1d2v_trafo(self)
    class(sll_t_time_propagator_pic_vm_1d2v_trafo_helper), intent( inout ) :: self !< time propagator object

    if( self%boundary ) then
       print*, 'left boundary', self%counter_left
       print*, 'right boundary', self%counter_right
    end if

    call self%linear_solver_1%free()
    call self%linear_operator_1%free()
    call self%particle_mass_1%free()

    call self%linear_solver_0%free()
    call self%linear_operator_0%free()
    call self%particle_mass_0%free()

    deallocate( self%j_dofs )
    deallocate( self%j_dofs_local )
    deallocate( self%particle_mass_1_local )
    deallocate( self%particle_mass_0_local )
    self%maxwell_solver => null()
    self%kernel_smoother_0 => null()
    self%kernel_smoother_1 => null()
    self%particle_group => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()
    self%map=> null()

  end subroutine delete_pic_vm_1d2v_trafo


end module sll_m_time_propagator_pic_vm_1d2v_trafo_helper
