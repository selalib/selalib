module sll_m_time_propagator_pic_vm_1d2v_helper
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
  
  use sll_m_collective, only: &
       sll_o_collective_allreduce, &
       sll_v_world_collective

  use sll_m_fft, only: &
       sll_t_fft, &
       sll_s_fft_init_r2r_1d, &
       sll_s_fft_exec_r2r_1d, &
       sll_s_fft_free, &
       sll_p_fft_forward, &
       sll_p_fft_backward

  use sll_m_filter_base_1d, only: &
       sll_c_filter_base_1d
   
  use sll_m_linear_solver_cg, only : &
       sll_t_linear_solver_cg

  use sll_m_linear_solver_mgmres, only : &
       sll_t_linear_solver_mgmres
  
  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_maxwell_1d_base, only: &
       sll_c_maxwell_1d_base

  use sll_m_maxwell_1d_fem, only: &
       sll_t_maxwell_1d_fem

  use sll_mpi, only: &
       mpi_sum, &
       mpi_max

  use sll_m_particle_group_base, only: &
       sll_t_particle_array

  use sll_m_particle_mesh_coupling_base_1d, only: &
    sll_c_particle_mesh_coupling_1d

  use sll_m_particle_mesh_coupling_spline_1d, only: &
    sll_t_particle_mesh_coupling_spline_1d

  use sll_m_particle_mesh_coupling_spline_smooth_1d, only: &
    sll_t_particle_mesh_coupling_spline_smooth_1d

  use sll_m_particle_mesh_coupling_spline_strong_1d, only: &
    sll_t_particle_mesh_coupling_spline_strong_1d

  use sll_m_spline_fem_utilities, only : &
       sll_s_spline_fem_mass_line
  
  use sll_m_spline_fem_utilities_sparse, only : &
       sll_s_spline_fem_sparsity_mass

  use sll_m_time_propagator_base, only: &
    sll_c_time_propagator_base

  
  implicit none

  public :: &
    sll_t_time_propagator_pic_vm_1d2v_helper

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Hamiltonian splitting type for Vlasov-Maxwell 1d2v
  type :: sll_t_time_propagator_pic_vm_1d2v_helper
     class(sll_c_maxwell_1d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_0  !< Kernel smoother (order p+1)
     class(sll_c_particle_mesh_coupling_1d), pointer :: kernel_smoother_1  !< Kernel smoother (order p)
     class(sll_t_particle_array), pointer  :: particle_group    !< Particle group

     sll_int32 :: spline_degree !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx !< Size of the domain
     sll_real64 :: x_min !< Lower bound for x domain
     sll_real64 :: delta_x !< Grid spacing

     sll_real64 :: cell_integrals_0(4) !< Integral over the spline function on each interval (order p+1)
     sll_real64 :: cell_integrals_1(3) !< Integral over the spline function on each interval (order p)


     sll_real64, pointer     :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
     sll_real64, allocatable :: efield_dofs_iter(:,:) !< DoFs describing the two components of the electric field
     sll_real64, pointer     :: bfield_dofs(:)   !< DoFs describing the magnetic field
     sll_real64, allocatable :: j_dofs(:,:)      !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:,:)!< MPI-processor local part of one component of \a j_dofs
     sll_real64, allocatable :: rho_dofs(:)      !< DoFs for kernel representation of charge density. 
     sll_real64, allocatable :: rho_dofs_local(:)!< MPI-processor local part of one component of \a rho_dofs

     sll_real64, allocatable :: particle_mass_0(:,:) !< Array to hold the 2*spline_degree+1 diagonals of the matrix A_0 M_p A_0^T
     sll_real64, allocatable :: particle_mass_1(:,:) !< Array to hold the 2*(spline_degree-1)+1 diagonals of the matrix A_1 M_p A_1^T
     sll_real64, allocatable :: particle_mass_0_local(:,:) !< Array to hold the 2*spline_degree+1 diagonals of the matrix A_0 M_p A_0^T
     sll_real64, allocatable :: particle_mass_1_local(:,:) !< Array to hold the 2*(spline_degree-1)+1 diagonals of the matrix A_1 M_p A_1^T
     sll_real64, allocatable :: residual(:)

     sll_real64, allocatable :: mass_line_0(:) !< mass line for degree+1
     sll_real64, allocatable :: mass_line_1(:) !< mass line for degree

     sll_real64, allocatable :: xnew(:,:)
     sll_real64, allocatable :: vnew(:,:,:)

     ! SPM object for solution of schur complement
     type(sll_t_matrix_csr) :: schur_0 !< Schur operator for advect_ev 
     type(sll_t_matrix_csr) :: schur_1 !< Schur operator for advect_ev 
     
     type(sll_t_linear_solver_cg) :: linear_solver_0 !< Linear solver for Schur operator for advect_ev 
     type(sll_t_linear_solver_cg) :: linear_solver_1 !< Linear solver for Schur operator for advect_ev 
     
     sll_real64, allocatable :: dof2_vec1(:)
     sll_real64, allocatable :: dof2_vec2(:)

     sll_int32 :: niter(100000)
     sll_int32 :: nevaliter(100000)
     sll_real64, allocatable :: nsubiter(:,:)
     sll_int32 :: iter_counter = 0
     sll_int32, allocatable :: sub_iter_counter(:)
     sll_int32 :: eval_counter = 0
     sll_int32 :: max_iter = 10
     sll_int32 :: max_iter_sub = 100
     sll_real64 :: tolerance = 1D-12
     sll_real64 :: tolerance_sub = 1D-10
     sll_int32 :: n_failed = 0


     sll_real64 :: solver_tolerance = 1D-12 !< solver tolerance

     sll_int32 :: n_sub_iter(2) = [8,2] !< number of subiterations

     character(len=256) :: output_prefix = "disgrad"

     sll_int32 :: n_particles_max ! Helper variable for size of temporary particle arrays

     !FFT filter
     type(sll_t_fft) :: filter_fw, filter_bw !< fft filter

     logical    :: smooth !< logical to store if smoothing is applied
     sll_int32  :: size_particle_mass_0
     sll_int32  :: size_particle_mass_1

     sll_real64, allocatable     :: efield_filter(:,:) !< DoFs describing the two components of the electric field
     sll_real64, allocatable     :: bfield_filter(:)   !< DoFs describing the magnetic field

     sll_real64, allocatable :: efield_to_val(:,:)
     sll_real64, allocatable :: bfield_to_val(:)

     class(sll_c_filter_base_1d), pointer :: filter !< filter

     
   contains
     procedure :: advect_x => advect_x_pic_vm_1d2v_helper
     procedure :: advect_vb => advect_vb_pic_vm_1d2v_helper
     procedure :: advect_eb => advect_eb_pic_vm_1d2v_helper
     
     procedure :: init => initialize_pic_vm_1d2v_helper !> Initialize the type
     procedure :: free => delete_pic_vm_1d2v_helper !> Finalization

     procedure :: advect_e_sub => advect_e_sub_pic_vm_1d2v_helper
     procedure :: advect_e_start_avf => advect_e_start_avf_pic_vm_1d2v_helper
     procedure :: avf_for_start

     procedure :: reinit_fields

  end type sll_t_time_propagator_pic_vm_1d2v_helper

contains

  !---------------------------------------------------------------------------!
  !> Advection of x part separately
  subroutine advect_x_pic_vm_1d2v_helper ( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_helper), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    sll_int32 :: i_part, i_sp
    sll_real64 :: xi(3), vi(3)

    do i_sp = 1, self%particle_group%n_species
       do i_part=1,self%particle_group%group(i_sp)%n_particles  
          ! Read out particle position and velocity
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          
          xi(1) = xi(1) + dt * vi(1)
          xi(1) = modulo(xi(1), self%Lx)
          call self%particle_group%group(i_sp)%set_x ( i_part, xi )
       end do
    end do
    
  end subroutine advect_x_pic_vm_1d2v_helper

  !---------------------------------------------------------------------------!
  !> vxB part (exact solution)
  subroutine advect_vb_pic_vm_1d2v_helper ( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_helper), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    !local variables
    sll_int32  :: i_part, i_sp
    sll_real64 :: qmdt
    sll_real64 :: vi(3), v_new(3), xi(3)
    sll_real64 :: bfield, cs, sn

    
    call self%maxwell_solver%transform_dofs( self%bfield_filter, self%bfield_to_val, 0 )
    
    do i_sp = 1, self%particle_group%n_species
       qmdt = self%particle_group%group(i_sp)%species%q_over_m()*dt;

       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          call self%kernel_smoother_1%evaluate &
               (xi(1), self%bfield_to_val, bfield)
          
          bfield = qmdt*bfield
          
          cs = cos(bfield)
          sn = sin(bfield)
          
          v_new(1) = cs * vi(1) + sn * vi(2)
          v_new(2) = -sn* vi(1) + cs * vi(2)
          
          call self%particle_group%group(i_sp)%set_v( i_part, v_new )
       end do
    end do


  end subroutine advect_vb_pic_vm_1d2v_helper

  !---------------------------------------------------------------------------!
  !> Solution of the curl part of Faraday and Ampere (solved as a system with inversion in Fourier space
  subroutine advect_eb_pic_vm_1d2v_helper ( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_helper), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step


    call self%maxwell_solver%compute_curl_part( dt, self%efield_dofs(:,2), self%bfield_dofs, 1.0_f64 )
    

    call self%filter%apply( self%efield_dofs(:,2), self%efield_filter(:,2) )
    call self%filter%apply( self%bfield_dofs, self%bfield_filter )
    
  end subroutine advect_eb_pic_vm_1d2v_helper

  
  
  !> Operator for e,x-part with subcycling
  subroutine advect_e_sub_pic_vm_1d2v_helper ( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_helper), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    ! local variables
    sll_real64 :: residuals(4)
    sll_real64 :: residuals_all(4)
    sll_real64 :: residual
    sll_real64 :: qoverm
    sll_real64 :: dt_sub(self%particle_group%n_species)
    sll_real64 :: xi(3), vi(3), wi(1)
    sll_real64 :: efield_dofs(self%kernel_smoother_1%n_dofs,2)
    sll_real64 :: efield_dofs_old(self%kernel_smoother_1%n_dofs,2)
    sll_int32 :: i_part, i_sub, niter, i_sp
    
    

    do i_sp = 1, self%particle_group%n_species
       dt_sub(i_sp) = dt/real( self%n_sub_iter(i_sp), f64 )
    end do
    ! Set to zero
    self%j_dofs_local = 0.0_f64

    niter = 0
    
    ! Set initial value for electric field
    ! Two options to initialize the efield: either by old time step value or by value from avf propagator
    !efield_dofs = self%efield_dofs
    call avf_for_start ( self, dt, efield_dofs )
    
    efield_dofs_old = efield_dofs
    residual = 1.0_f64
    do while ( (residual > self%tolerance) .and. niter < self%max_iter )
       niter = niter + 1
       residuals = 0.0_f64       
       efield_dofs = 0.5_f64*( self%efield_dofs + efield_dofs )
       ! Set to zero
       self%j_dofs_local = 0.0_f64
       do i_sp = 1, self%particle_group%n_species
          qoverm = self%particle_group%group(i_sp)%species%q_over_m();
          do i_part = 1,self%particle_group%group(i_sp)%n_particles
             vi = self%particle_group%group(i_sp)%get_v(i_part)
             xi = self%particle_group%group(i_sp)%get_x(i_part)
             
             ! Get charge for accumulation of j
             wi = self%particle_group%group(i_sp)%get_charge(i_part)
             
             do i_sub = 1, self%n_sub_iter(i_sp)
                call subcycle_xv( self, dt_sub(i_sp), qoverm, efield_dofs, wi, xi(1:1), &
                     vi(1:2), self%sub_iter_counter(i_sp) )
             end do

             residuals(1) = residuals(1) + (xi(1)-self%xnew(i_sp,i_part))**2*abs(wi(1))
             residuals(2) = residuals(2) + (vi(1)-self%vnew(i_sp,1,i_part))**2*abs(wi(1))
             residuals(3) = residuals(3) + (vi(2)-self%vnew(i_sp,2,i_part))**2*abs(wi(1))
             
             self%xnew(i_sp,i_part) = xi(1)
             self%vnew(i_sp,:,i_part) = vi(1:2)
          end do
       end do

       ! Update d_n
       self%j_dofs = 0.0_f64
       ! MPI to sum up contributions from each processor
       call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
            self%kernel_smoother_1%n_dofs, MPI_SUM, self%j_dofs(:,1) )
       call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
            self%kernel_smoother_1%n_dofs, MPI_SUM, self%j_dofs(:,2) )
       
       ! Update the electric field.
       efield_dofs = self%efield_dofs
       call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,1), 1, efield_dofs(:,1) )
       call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,2), 2, efield_dofs(:,2) )
       residuals(4) = (sum((efield_dofs_old-efield_dofs)**2))*self%delta_x!maxval( efield_dofs_old - efield_dofs ) ! Here the error norm can be changed
       
       call sll_o_collective_allreduce( sll_v_world_collective, residuals, 4, MPI_MAX, residuals_all )
       residuals_all = sqrt(residuals_all)
       residual = residuals_all(4)!maxval(residuals_all) ! Here the error criterion for the iteration can be changed
       efield_dofs_old = efield_dofs
    end do
    self%efield_dofs = efield_dofs
    do i_sp =1,self%particle_group%n_species
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vi(1:2) = self%vnew(i_sp,:,i_part)
          xi(1) = modulo(self%xnew(i_sp, i_part), self%Lx )
          call self%particle_group%group(i_sp)%set_v( i_part, vi )
          call self%particle_group%group(i_sp)%set_x( i_part, xi )
       end do
    end do

    self%iter_counter = self%iter_counter + 1
    self%niter(self%iter_counter) = niter
    do i_sp = 1,self%particle_group%n_species
       self%nsubiter(i_sp,self%iter_counter) = real(self%sub_iter_counter(i_sp), f64)/ &
            real( self%particle_group%group(i_sp)%n_particles  * niter * self%n_sub_iter(i_sp),f64)
    end do
    self%sub_iter_counter = 0
    
  end subroutine advect_e_sub_pic_vm_1d2v_helper

  

  !> Helper function for subcycle (using Picard iteration)
  subroutine subcycle_xv( self, dt, qoverm, efield_dofs, wi, position, velocity, sub_iter_counter )
    class(sll_t_time_propagator_pic_vm_1d2v_helper), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_real64,                                     intent(in)    :: qoverm
    sll_real64,                                     intent(in)    :: efield_dofs(:,:)
    sll_real64,                                     intent(in   ) :: wi(1)
    sll_real64,                                     intent(inout) :: position(1)
    sll_real64,                                     intent(inout) :: velocity(2)
    sll_int32,                                      intent(inout) :: sub_iter_counter
    
    ! local variables
    sll_real64 :: xnew(1), vnew(2), vbar, xold(1), vold(2)
    sll_real64 :: efield(2)
    sll_int32 :: niter
    sll_real64 :: residuals(3), residuals_all(3), residual
    
    !xnew = position + dt * velocity(1)
    vnew = velocity
    xnew = position
    
    niter = 0
    residual = 1.0_f64
    do while ( (residual > self%tolerance_sub) .and. niter < self%max_iter_sub )
       niter = niter + 1
       residuals = 0.0_f64  
       vbar = 0.5_f64*(velocity(1)+vnew(1))
       vold = vnew
       xold = xnew
       xnew = position + dt*vbar
       
       if ( abs(vbar) > 1.0D-16 ) then
             
          select type ( q=> self%kernel_smoother_1 )
          type is (sll_t_particle_mesh_coupling_spline_1d )
             call q%evaluate_int &
                  (  position, xnew,  vbar, efield_dofs(:,1),  &
                  efield(1) )
          type is (sll_t_particle_mesh_coupling_spline_smooth_1d )
             call q%evaluate_int &
                  (  position, xnew,  vbar, efield_dofs(:,1),  &
                  efield(1) )
          end select
          select type ( q0=> self%kernel_smoother_0 )
          type is (sll_t_particle_mesh_coupling_spline_1d )
             call q0%evaluate_int &
                  ( position, xnew, vbar, &
                  efield_dofs(:,2),  &
                  efield(2) )
          type is (sll_t_particle_mesh_coupling_spline_smooth_1d )
             call q0%evaluate_int &
                  ( position, xnew, vbar, &
                  efield_dofs(:,2),  &
                  efield(2) )
          end select
       else
          call self%kernel_smoother_1%evaluate &
               (position, efield_dofs(:,1), efield(1) )
          efield(1) = efield(1)*dt
          
          call self%kernel_smoother_0%evaluate &
               (position, efield_dofs(:,2), efield(2) )
          efield(2) = efield(2)*dt
       end if
       
       vnew = velocity + qoverm * efield(1:2)

       residuals(1) = abs(xnew(1)-xold(1))
       residuals(2) = abs(vnew(1)-vold(1))
       residuals(3) = abs(vnew(2)-vold(2))
       call sll_o_collective_allreduce( sll_v_world_collective, residuals, 3, MPI_MAX, residuals_all )
       residual = maxval(residuals_all)

    end do

    if ( residual > self%tolerance_sub ) then
       print*, 'Warning: Inner iteration of iteration no.', self%iter_counter+1 ,'did not converge.', residuals_all
       !self%n_failed = self%n_failed+1
    end if

    vbar = 0.5_f64*(velocity(1)+vnew(1))
    xnew = position + dt*vbar
    if ( abs(vbar) > 1.0D-16 ) then
             
       select type ( q=> self%kernel_smoother_1 )
       type is (sll_t_particle_mesh_coupling_spline_1d )
          call q%add_current_evaluate &
               ( position, xnew, wi(1), vbar, efield_dofs(:,1), self%j_dofs_local(:,1), &
               efield(1) )
       type is (sll_t_particle_mesh_coupling_spline_smooth_1d )
          call q%add_current_evaluate &
               ( position, xnew, wi(1), vbar, efield_dofs(:,1), self%j_dofs_local(:,1), &
               efield(1) )
       end select
       select type ( q0=> self%kernel_smoother_0 )
       type is (sll_t_particle_mesh_coupling_spline_1d )
          call q0%add_current_evaluate &
               ( position, xnew, wi(1)*0.5_f64*(velocity(2)+vnew(2))/vbar, vbar, &
               efield_dofs(:,2),  self%j_dofs_local(:,2), &
               efield(2) )
       type is (sll_t_particle_mesh_coupling_spline_smooth_1d )
          call q0%add_current_evaluate &
               ( position, xnew, wi(1)*0.5_f64*(velocity(2)+vnew(2))/vbar, vbar, &
               efield_dofs(:,2),  self%j_dofs_local(:,2), &
               efield(2) )
       end select
    else
       call self%kernel_smoother_1%evaluate &
            (position, efield_dofs(:,1), efield(1) )
       efield(1) = efield(1)*dt
       
       call self%kernel_smoother_0%add_charge( position, &
            wi(1)*0.5_f64*(velocity(2)+vnew(2))*dt, self%j_dofs_local(:,2) )
       call self%kernel_smoother_0%evaluate &
            (position, efield_dofs(:,2), efield(2) )
       efield(2) = efield(2)*dt
    end if

    vnew = velocity + qoverm*efield(1:2)
    
    velocity = vnew
    position = modulo(xnew, self%Lx)

    sub_iter_counter = sub_iter_counter + niter + 1
    
  end subroutine subcycle_xv


  
  subroutine filter( field_in, field_out, n_dofs )
    sll_real64,                            intent(in)  :: field_in(:,:) !< array for the coefficients of the efields 
    sll_real64,                            intent(out)  :: field_out(:,:) !< array for the coefficients of the efields
    sll_int32, intent(in) :: n_dofs

    sll_int32 :: j
    
    ! Filter
    field_out(1,:) = 0.25_f64*( field_in(1,:)*2.0_f64 + &
               field_in(n_dofs,:)+field_in(2,:))
    do j=2,n_dofs-1
       field_out(j,:) = 0.25_f64*( field_in(j,:)*2.0_f64 + &
            field_in(j-1,:)+field_in(j+1,:))
    end do
    field_out(n_dofs,:) = 0.25_f64*( field_in(n_dofs,:)*2.0_f64 + &
         field_in(n_dofs-1,:)+field_in(1,:))

  end subroutine filter



 !---------------------------------------------------------------------------!
  !> Constructor.
  subroutine initialize_pic_vm_1d2v_helper(&
       self, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       particle_group, &
       efield_dofs, &
       bfield_dofs, &
       x_min, &
       Lx, &
       filter, &
       filename, &
       build_particle_mass) 
    class(sll_t_time_propagator_pic_vm_1d2v_helper), intent(out) :: self !< time splitting object 
    class(sll_c_maxwell_1d_base), target,          intent(in)  :: maxwell_solver      !< Maxwell solver
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_0  !< Kernel smoother
    class(sll_c_particle_mesh_coupling_1d), target,          intent(in)  :: kernel_smoother_1  !< Kernel smoother
    class(sll_t_particle_array), target,      intent(in)  :: particle_group !< Particle group
    sll_real64, target,                            intent(in)  :: efield_dofs(:,:) !< array for the coefficients of the efields 
    sll_real64, target,                            intent(in)  :: bfield_dofs(:) !< array for the coefficients of the bfield
    sll_real64,                                     intent(in)  :: x_min !< Lower bound of x domain
    sll_real64,                                     intent(in)  :: Lx !< Length of the domain in x direction.
    class( sll_c_filter_base_1d ), intent( in ), target :: filter
    character(len=*), intent(in) :: filename
    logical, optional, intent(in) :: build_particle_mass

    
    !local variables
    sll_int32 :: ierr, i_sp, io_stat
    sll_real64 :: iter_tolerance, sub_iter_tolerance, maxwell_tolerance
    sll_int32 :: max_iter, max_iter_sub, n_sub_iter(2), input_file
    sll_int32 :: n_span_particle_mass
    character(len=256) :: file_prefix
    logical :: output_fields, output_particles
    logical :: build_particle_mass_loc = .true.
    
    namelist /output/ file_prefix, output_fields, output_particles
    namelist /time_solver/ maxwell_tolerance
    namelist /time_iterate/ iter_tolerance, sub_iter_tolerance, max_iter, max_iter_sub, n_sub_iter
!!$
    
    if ( present( build_particle_mass ))  build_particle_mass_loc = build_particle_mass
    
    self%maxwell_solver => maxwell_solver
    self%kernel_smoother_0 => kernel_smoother_0
    self%kernel_smoother_1 => kernel_smoother_1
    self%particle_group => particle_group
    self%efield_dofs => efield_dofs
    self%bfield_dofs => bfield_dofs

    ! Check that n_dofs is the same for both kernel smoothers.
    SLL_ASSERT( self%kernel_smoother_0%n_dofs == self%kernel_smoother_1%n_dofs )
!!$
    SLL_ALLOCATE(self%j_dofs(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%j_dofs_local(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%efield_dofs_iter(self%kernel_smoother_0%n_dofs,2), ierr)
    SLL_ALLOCATE(self%residual(self%kernel_smoother_0%n_dofs*2), ierr)

    self%spline_degree = self%kernel_smoother_0%spline_degree
    self%x_min = x_min
    self%Lx = Lx
    self%delta_x = self%Lx/real(self%kernel_smoother_1%n_dofs, f64)
    
    self%cell_integrals_1 = [0.5_f64, 2.0_f64, 0.5_f64]
    self%cell_integrals_1 = self%cell_integrals_1 / 3.0_f64

    self%cell_integrals_0 = [1.0_f64,11.0_f64,11.0_f64,1.0_f64]
    self%cell_integrals_0 = self%cell_integrals_0 / 24.0_f64


    if ( self%maxwell_solver%s_deg_0 > -1 .and. build_particle_mass_loc .eqv. .true. ) then
       allocate( self%mass_line_0( self%maxwell_solver%s_deg_0 +1) )
       allocate( self%mass_line_1( self%maxwell_solver%s_deg_0 ) )
       call sll_s_spline_fem_mass_line ( self%maxwell_solver%s_deg_0, self%mass_line_0 )
       call sll_s_spline_fem_mass_line ( self%maxwell_solver%s_deg_0-1, self%mass_line_1 )

       ! Scale with dx
       self%mass_line_0 = self%mass_line_0 * self%delta_x
       self%mass_line_1 = self%mass_line_1 * self%delta_x

       self%smooth = .false.
       n_span_particle_mass = self%spline_degree+1
       select type( pm=>kernel_smoother_0)
       type is ( sll_t_particle_mesh_coupling_spline_smooth_1d )
          self%smooth = .true.
          n_span_particle_mass = self%spline_degree+3
       end select
       self%size_particle_mass_0 = n_span_particle_mass* self%kernel_smoother_0%n_dofs
       self%size_particle_mass_1 = (n_span_particle_mass-1)* self%kernel_smoother_1%n_dofs
    
       allocate( self%particle_mass_0(n_span_particle_mass, self%kernel_smoother_0%n_dofs))
       allocate( self%particle_mass_0_local(n_span_particle_mass, self%kernel_smoother_0%n_dofs))
       allocate( self%particle_mass_1(n_span_particle_mass-1, self%kernel_smoother_1%n_dofs) )
       allocate( self%particle_mass_1_local(n_span_particle_mass-1, self%kernel_smoother_1%n_dofs) )
    
       call sll_s_spline_fem_sparsity_mass( n_span_particle_mass-1, self%kernel_smoother_0%n_dofs, &
            self%schur_0 ) 
       call sll_s_spline_fem_sparsity_mass( n_span_particle_mass-2, self%kernel_smoother_1%n_dofs, &
            self%schur_1 )
    end if
       
    SLL_ALLOCATE(self%rho_dofs(self%kernel_smoother_0%n_dofs), ierr)
    SLL_ALLOCATE(self%rho_dofs_local(self%kernel_smoother_0%n_dofs), ierr)
    
    self%n_particles_max = self%particle_group%group(1)%n_particles
    do i_sp = 2,self%particle_group%n_species       
       self%n_particles_max = max(self%n_particles_max, self%particle_group%group(i_sp)%n_particles )
    end do
    
    SLL_ALLOCATE( self%xnew(self%particle_group%n_species,self%n_particles_max), ierr )
    SLL_ALLOCATE( self%vnew(self%particle_group%n_species,2,self%n_particles_max), ierr )


    ! Read tolerance and max iter from nml file

    open(newunit = input_file, file=filename, status='old',IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'initialize_pic_vm_1d2v_helper() failed to open nml-file.'
    else
       read(input_file, output )
       self%output_prefix = trim(file_prefix)
       read(input_file, time_solver )
       self%solver_tolerance = maxwell_tolerance
       read(input_file, time_iterate )
       self%tolerance = iter_tolerance
       self%tolerance_sub = sub_iter_tolerance
       self%max_iter = max_iter
       self%max_iter_sub = max_iter_sub
       self%n_sub_iter = n_sub_iter
       close(input_file)
    end if

    ! FFT filter
     call sll_s_fft_init_r2r_1d( self%filter_fw, self%kernel_smoother_0%n_dofs, self%j_dofs(:,1), self%j_dofs_local(:,1), sll_p_fft_forward, normalized=.false. )
     call sll_s_fft_init_r2r_1d( self%filter_bw, self%kernel_smoother_0%n_dofs, self%j_dofs_local(:,1), self%j_dofs(:,1), sll_p_fft_backward, normalized=.true. )

     allocate(self%sub_iter_counter( self%particle_group%n_species ) )
     self%sub_iter_counter = 0
     allocate(self%nsubiter(self%particle_group%n_species,100000 ) )

     ! For binomial filter
     self%filter => filter
     SLL_ALLOCATE(self%efield_filter(self%kernel_smoother_0%n_dofs,2), ierr)
     SLL_ALLOCATE(self%bfield_filter(self%kernel_smoother_0%n_dofs), ierr)
     SLL_ALLOCATE(self%efield_to_val(self%kernel_smoother_0%n_dofs,2), ierr)
     SLL_ALLOCATE(self%bfield_to_val(self%kernel_smoother_0%n_dofs), ierr)
     
  end subroutine initialize_pic_vm_1d2v_helper

  !---------------------------------------------------------------------------!
  !> Destructor.
  subroutine delete_pic_vm_1d2v_helper(self)
    class(sll_t_time_propagator_pic_vm_1d2v_helper), intent( inout ) :: self !< time splitting object 

    sll_int32 :: file_id, j

    open(newunit=file_id, file=trim(self%output_prefix)//"-outer-iters.dat", action='write')
    do j=1,self%iter_counter
       write(file_id,*) self%niter(j)
    end do
    close(file_id)
    open(newunit=file_id, file=trim(self%output_prefix)//"-inner-iters.dat", action='write')
    do j=1, self%iter_counter
       write(file_id,*) self%nsubiter(:,j)
    end do
    close(file_id)
    print*, 'No. of failed iterations:', self%n_failed
    
    deallocate(self%j_dofs)
    deallocate(self%j_dofs_local)
    self%maxwell_solver => null()
    self%kernel_smoother_0 => null()
    self%kernel_smoother_1 => null()
    self%particle_group => null()
    self%efield_dofs => null()
    self%bfield_dofs => null()

    call self%schur_0%free()
    call self%schur_1%free()
  end subroutine delete_pic_vm_1d2v_helper



  subroutine assemble_schur ( n_dofs, degree, mass, particle_mass, schur )
    sll_int32,  intent( in    ) :: n_dofs
    sll_int32,  intent( in    ) :: degree
    sll_real64, intent( in    ) :: mass(degree+1)
    sll_real64, intent( in    ) :: particle_mass(degree+1,n_dofs)
    type(sll_t_matrix_csr), intent( inout ) :: schur

    !local variables
    sll_int32 :: row, column, ind

    ind = 1
    ! For the first degree rows we need to put the first part to the back due to periodic boundaries
    do row = 1, degree

       do column = -row+1, -1
          schur%arr_a(1,ind) = mass(1-column) + particle_mass(1-column ,row+column)
          ind = ind+1
       end do
       schur%arr_a(1,ind:ind+degree) = mass + particle_mass(:, row)
       ind = ind+degree+1

       do column = -degree, -row
          schur%arr_a(1,ind) = mass(1-column) + &
               particle_mass(1-column,  n_dofs+row+column)
          ind = ind+1
       end do
       
    end do
    
    do row = degree+1, n_dofs-degree

       do column = -degree, -1
          schur%arr_a(1,ind) = mass(1-column) + particle_mass(1-column ,row+column)
          ind = ind+1
       end do
       schur%arr_a(1,ind:ind+degree) = mass + particle_mass(:, row)
       ind = ind+degree+1
    end do
    
    ! For the last degree rows, we need to put the second part to the front due to periodic boundaries
    do row = n_dofs-degree+1, n_dofs

       do column = n_dofs-row+1, degree
          schur%arr_a(1,ind) = mass(column+1) + particle_mass(column+1, row)
          ind = ind+1
       end do
       
       do column = -degree, -1
          schur%arr_a(1,ind) = mass(1-column) + particle_mass(1-column ,row+column)
          ind = ind+1
       end do

       do column = 0, n_dofs-row
          schur%arr_a(1,ind) = mass(column+1) + particle_mass(column+1, row)
          ind = ind+1
       end do
    end do

    SLL_ASSERT( ind == schur%n_nnz+1)

    
  end subroutine assemble_schur


  

subroutine assemble_schur_smooth ( n_dofs, degree, mass, particle_mass, schur )
 sll_int32,  intent( in    ) :: n_dofs
 sll_int32,  intent( in    ) :: degree
 sll_real64, intent( in    ) :: mass(degree+1)
 sll_real64, intent( in    ) :: particle_mass(degree+3,n_dofs)
 type(sll_t_matrix_csr), intent( inout ) :: schur

 !local variables
 sll_int32 :: row, column, ind

 
 ind = 1
 ! For the first degree rows we need to put the first part to the back due to periodic boundaries
 do row = 1, degree+2
    if (row == degree+2) then
       schur%arr_a(1,ind) = particle_mass(degree+2,row-degree-1)!row-degree-1)
       ind = ind+1
    end if
    
    do column = max(-row+1,-degree), -1
       schur%arr_a(1,ind) = mass(1-column) + particle_mass(1-column ,row+column)
       ind = ind+1
    end do
    schur%arr_a(1,ind:ind+degree) = mass + particle_mass(1:degree+1, row)
    schur%arr_a(1,ind+degree+1:ind+degree+2) = particle_mass(degree+2:degree+3, row)
    ind = ind+degree+3
    
    schur%arr_a(1,ind) = particle_mass(degree+3,n_dofs+row-degree-2)
    ind = ind+1
    if (row<degree+2) then
       schur%arr_a(1,ind) = particle_mass(degree+2,n_dofs+row-degree-1)!row-degree-1)
       ind = ind+1
       do column = -degree, -row
          schur%arr_a(1,ind) = mass(1-column) + &
               particle_mass(1-column,  n_dofs+row+column)
          ind = ind+1
       end do
    end if

 end do

 do row = degree+3, n_dofs-degree-2
    schur%arr_a(1,ind) = particle_mass(degree+3,row-degree-2)
    ind = ind+1
    schur%arr_a(1,ind) = particle_mass(degree+2,row-degree-1)
    ind = ind+1
    do column = -degree, -1
       schur%arr_a(1,ind) = mass(1-column) + particle_mass(1-column ,row+column)
       ind = ind+1
    end do
    schur%arr_a(1,ind:ind+degree) = mass + particle_mass(1:degree+1, row)
    schur%arr_a(1,ind+degree+1:ind+degree+2) = particle_mass(degree+2:degree+3, row)
    ind = ind+degree+3
 end do

 ! For the last degree rows, we need to put the second part to the front due to periodic boundaries
 do row = n_dofs-degree-1, n_dofs

    if (row>n_dofs-degree-1) then
       do column = n_dofs-row+1, degree
          schur%arr_a(1,ind) = mass(column+1) + particle_mass(column+1, row)
          ind = ind+1
       end do
       schur%arr_a(1,ind) = particle_mass(degree+2, row)
       ind = ind+1
    end if
    schur%arr_a(1,ind) = particle_mass(degree+3, row)
    ind = ind+1
    
    schur%arr_a(1,ind) = particle_mass(degree+3,row-degree-2)
    ind = ind+1
    schur%arr_a(1,ind) = particle_mass(degree+2,row-degree-1)
    ind = ind+1
    do column = -degree, -1
       schur%arr_a(1,ind) = mass(1-column) + particle_mass(1-column ,row+column)
       ind = ind+1
    end do

    do column = 0, min(n_dofs-row,degree)
       schur%arr_a(1,ind) = mass(column+1) + particle_mass(column+1, row)
       ind = ind+1
    end do
    if (row == (n_dofs-degree-1) ) then
       schur%arr_a(1,ind) = particle_mass(degree+2, row)
       ind = ind+1
    end if
 end do

 SLL_ASSERT( ind == schur%n_nnz+1)


end subroutine assemble_schur_smooth

  
  !> Operator for first variant without subcycling (Picard iteration, started by AVF)
  subroutine advect_e_start_avf_pic_vm_1d2v_helper ( self, dt )
    class(sll_t_time_propagator_pic_vm_1d2v_helper), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step

    ! local variables
    sll_int32 :: i_part, i_sp
    sll_int32 :: niter
    sll_real64 :: vi(3), xi(3), wi(1), xnew(3), vnew(3), vbar
    sll_real64 :: qoverm
    sll_real64 :: efield(2)
    sll_real64 :: residual, residuals(4), residuals_all(4)
    sll_real64 :: efield_dofs(self%kernel_smoother_1%n_dofs,2)
    !sll_real64 :: efield_dofs2(self%kernel_smoother_0%n_dofs)
    sll_real64 :: efield_dofs_old(self%kernel_smoother_1%n_dofs,2)
    sll_int32  :: ndh

    ndh = self%kernel_smoother_0%n_dofs/2
    
    residuals = 0.0_f64

    ! Instead of a first particle loop as in the standard Picard iteration, we now
    ! perform the advect_x and advect_e part of the AVF scheme
    efield_dofs = self%efield_dofs
    !call avf_for_start ( self, dt, efield_dofs )
    do i_sp = 1, self%particle_group%n_species
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vnew = self%particle_group%group(i_sp)%get_v(i_part)
          xnew = self%particle_group%group(i_sp)%get_x(i_part)
          self%xnew(i_sp,i_part) = xnew(1)
          self%vnew(i_sp,:,i_part) = vnew(1:2)
       end do
    end do
    
    efield_dofs_old = efield_dofs
    residual = 1.0_f64
    niter = 1
    
    !print*, 'Iteration', self%iter_counter+1
    do while ( (residual > self%tolerance) .and. niter < self%max_iter )
       niter = niter + 1
       residuals = 0.0_f64
       call self%filter%apply( efield_dofs_old(:,1), efield_dofs(:,1) )
       call self%filter%apply( efield_dofs_old(:,2), efield_dofs(:,2) )
       efield_dofs = 0.5_f64*(self%efield_filter + efield_dofs )
       
       call self%maxwell_solver%transform_dofs( efield_dofs(:,1), self%efield_to_val(:,1), 0 )    
       call self%maxwell_solver%transform_dofs( efield_dofs(:,2), self%efield_to_val(:,2), 1 )

       ! Set to zero
       self%j_dofs_local = 0.0_f64
       ! First particle loop
       do i_sp = 1, self%particle_group%n_species
          qoverm = self%particle_group%group(i_sp)%species%q_over_m();
          do i_part = 1,self%particle_group%group(i_sp)%n_particles
             vi = self%particle_group%group(i_sp)%get_v(i_part)
             xi = self%particle_group%group(i_sp)%get_x(i_part)
             
             ! Get charge for accumulation of j
             wi = self%particle_group%group(i_sp)%get_charge(i_part)
             
             vnew(1:2) = self%vnew(i_sp,1:2, i_part)
             xnew(1) = xi(1) + dt * 0.5_f64 * ( vi(1) + vnew(1) )
             
          
             vbar = 0.5_f64*(vi(1)+vnew(1))
             if ( abs(vbar) > 1.0D-16 ) then
                
                select type ( q=> self%kernel_smoother_1 )
                type is (sll_t_particle_mesh_coupling_spline_1d )
                   call q%add_current_evaluate &
                        ( xi(1), xnew(1), wi(1), vbar, self%efield_to_val(:,1), self%j_dofs_local(:,1), &
                        efield(1) )
                type is (sll_t_particle_mesh_coupling_spline_smooth_1d )
                   call q%add_current_evaluate &
                        ( xi(1), xnew(1), wi(1), vbar, self%efield_to_val(:,1), self%j_dofs_local(:,1), &
                        efield(1) )
                type is (sll_t_particle_mesh_coupling_spline_strong_1d )
                   call q%add_current_evaluate &
                        ( xi(1), xnew(1), wi(1), vbar, self%efield_to_val(:,1), self%j_dofs_local(:,1), &
                        efield(1) )
                end select
                select type ( q0=> self%kernel_smoother_0 )
                type is (sll_t_particle_mesh_coupling_spline_1d )
                   call q0%add_current_evaluate &
                        ( xi(1), xnew(1), wi(1)*(vi(2)+vnew(2))/(vi(1)+vnew(1)), vbar, &
                        self%efield_to_val(:,2), self%j_dofs_local(:,2), &
                        efield(2) )
                type is (sll_t_particle_mesh_coupling_spline_smooth_1d )
                   call q0%add_current_evaluate &
                        ( xi(1), xnew(1), wi(1)*(vi(2)+vnew(2))/(vi(1)+vnew(1)), vbar, &
                        self%efield_to_val(:,2), self%j_dofs_local(:,2), &
                        efield(2) )
                type is (sll_t_particle_mesh_coupling_spline_strong_1d )
                   call q0%add_current_evaluate &
                        ( xi(1), xnew(1), wi(1)*(vi(2)+vnew(2))/(vi(1)+vnew(1)), vbar, &
                        self%efield_to_val(:,2), self%j_dofs_local(:,2), &
                        efield(2) )
                end select
             else
                call self%kernel_smoother_1%evaluate &
                     (xi(1), self%efield_to_val(:,1), efield(1) )
                efield(1) = efield(1)*dt
                call self%kernel_smoother_0%add_charge( xi(1), &
                     wi(1)* 0.5_f64*(vi(2)+vnew(2))*dt, self%j_dofs_local(:,2) )
                call self%kernel_smoother_0%evaluate &
                     (xi(1), self%efield_to_val(:,2), efield(2) )
                efield(2) = efield(2)*dt
             end if
             vnew(1:2) = vi(1:2) + qoverm * efield(1:2)

             ! Here yo can change the residual criterion
             residuals(1) = residuals(1) + (xnew(1)-self%xnew(i_sp,i_part))**2*abs(wi(1))!max( residuals(1), abs(xnew(1)-self%xnew(i_sp,i_part)) )
             residuals(2) = residuals(2) + (vnew(1)-self%vnew(i_sp,1,i_part))**2*abs(wi(1))!max( residuals(2), abs(self%vnew(i_sp,1,i_part)-vnew(1)) )
             residuals(3) = residuals(3) + (vnew(2)-self%vnew(i_sp,2,i_part))**2*abs(wi(1))!max( residuals(3), abs(self%vnew(i_sp,2,i_part)-vnew(2)) )
             self%xnew(i_sp,i_part) = xnew(1)
             self%vnew(i_sp,:,i_part) = vnew(1:2)
          end do
       end do

       ! Update d_n
       self%j_dofs = 0.0_f64
       ! MPI to sum up contributions from each processor
       call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
            self%kernel_smoother_1%n_dofs, MPI_SUM, self%j_dofs(:,1) )
       call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
              self%kernel_smoother_1%n_dofs, MPI_SUM, self%j_dofs(:,2) )

       call self%filter%apply_inplace( self%j_dofs(:,1) )
       call self%filter%apply_inplace( self%j_dofs(:,2) )
       
       ! Filter
       !call filter( self%j_dofs, self%j_dofs_local, self%kernel_smoother_1%n_dofs )
       !self%j_dofs = self%j_dofs_local

       ! FFT filter
       !write(15, *) self%j_dofs(:,1)
       !call sll_s_fft_exec_r2r_1d ( self%filter_fw, self%j_dofs(:,1), self%j_dofs_local(:,1) )
       !write(14,*) self%j_dofs_local(:,1)
       !self%j_dofs_local(ndh+1-5:ndh+1+5,1) = 0.0_f64
       !call sll_s_fft_exec_r2r_1d ( self%filter_bw, self%j_dofs_local(:,1), self%j_dofs(:,1) )
       !write(16,*) self%j_dofs(:,1)
       !stop
       
       ! Update the electric field.
       efield_dofs = self%efield_dofs
       call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,1), 1, efield_dofs(:,1) )
       call self%maxwell_solver%compute_E_from_j(self%j_dofs(:,2), 2, efield_dofs(:,2) )

       if ( self%maxwell_solver%strong_ampere .eqv. .true. ) then
          residuals(4) =  self%maxwell_solver%l2norm_squared( efield_dofs_old(:,1)-&
               efield_dofs(:,1), self%maxwell_solver%s_deg_0 ) +  &
               self%maxwell_solver%l2norm_squared( efield_dofs_old(:,2)-&
               efield_dofs(:,2), self%maxwell_solver%s_deg_0 -1 )
       else
          residuals(4) =  self%maxwell_solver%l2norm_squared( efield_dofs_old(:,1)-&
               efield_dofs(:,1), self%maxwell_solver%s_deg_0 - 1 ) +  &
               self%maxwell_solver%l2norm_squared( efield_dofs_old(:,2)-&
               efield_dofs(:,2), self%maxwell_solver%s_deg_0  )
       end if
       ! Here you can change the residual norm
       residuals(4) = (sum((efield_dofs_old-efield_dofs)**2))*self%delta_x!maxval( efield_dofs_old - efield_dofs )


       call sll_o_collective_allreduce( sll_v_world_collective, residuals, 4, MPI_MAX, residuals_all )
       
       residuals_all = sqrt(residuals_all)
      residual = residuals_all(4)!max(residuals_all(4)), maxval(residuals_all(1:3))*0.01_f64)!residuals_all(4)!maxval(residuals_all)!residuals_all(4)! maxval(residuals_all)
        !call filter( efield_dofs, self%j_dofs_local, self%kernel_smoother_1%n_dofs )
       !efield_dofs = self%j_dofs_local
       ! FFT filter
       !write(15, *) self%j_dofs(:,1)
       !call sll_s_fft_exec_r2r_1d ( self%filter_fw, efield_dofs(:,1), self%j_dofs_local(:,1) )
       !write(14,*) self%j_dofs_local(:,1)
       !self%j_dofs_local(ndh+1-5:ndh+1+5,1) = 0.0_f64
       !call sll_s_fft_exec_r2r_1d ( self%filter_bw, self%j_dofs_local(:,1), efield_dofs(:,1) )
       !write(16,*) self%j_dofs(:,1)
       !stop
       
       efield_dofs_old = efield_dofs
    end do

!!$    ! NEW VERSION WITH SWAPP TO OTHER ITERATION
!!$    if ( residual > self%tolerance ) then
!!$        print*, 'Warning: Iteration no.', self%iter_counter+1 ,'did not converge.', residuals_all, niter
!!$        self%n_failed = self%n_failed+1
!!$        call advect_e_newtonavf(self, dt )
    
!!$        self%iter_failed = .true.
!!$    else
!!$
!!$       !print*, maxval(abs(self%efield_dofs-efield_dofs))
!!$       self%efield_dofs = efield_dofs
!!$       do i_sp = 1, self%particle_group%n_species
!!$          do i_part = 1, self%particle_group%group(i_sp)%n_particles
!!$             vnew(1:2) = self%vnew(i_sp,:,i_part)
!!$             xnew(1) = modulo(self%xnew(i_sp,i_part), self%Lx)
!!$             call self%particle_group%group(i_sp)%set_v( i_part, vnew )
!!$             call self%particle_group%group(i_sp)%set_x( i_part, xnew )
!!$          end do
!!$       end do
!!$       
!!$       self%iter_counter = self%iter_counter + 1
!!$       self%niter(self%iter_counter) = niter
!!$    end if
!!$    ! OLD VERSION TAKING FAILED ITERATION
    if ( residual > self%tolerance ) then
        print*, 'Warning: Iteration no.', self%iter_counter+1 ,'did not converge.', residuals_all, niter
       self%n_failed = self%n_failed+1
    end if

    !print*, maxval(abs(self%efield_dofs-efield_dofs))
    self%efield_dofs = efield_dofs
    call self%filter%apply( self%efield_dofs(:,1), self%efield_filter(:,1) )
    call self%filter%apply( self%efield_dofs(:,2), self%efield_filter(:,2) )
    do i_sp = 1, self%particle_group%n_species
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          vnew(1:2) = self%vnew(i_sp,:,i_part)
          xnew(1) = modulo(self%xnew(i_sp,i_part), self%Lx)
          call self%particle_group%group(i_sp)%set_v( i_part, vnew )
          call self%particle_group%group(i_sp)%set_x( i_part, xnew )
       end do
    end do

    self%iter_counter = self%iter_counter + 1
    self%niter(self%iter_counter) = niter
    
  end subroutine advect_e_start_avf_pic_vm_1d2v_helper


  !> AVF as first guess
  subroutine avf_for_start ( self, dt, efield_dofs )
    class(sll_t_time_propagator_pic_vm_1d2v_helper), intent(inout) :: self !< time splitting object 
    sll_real64,                                     intent(in)    :: dt   !< time step
    sll_real64,                                     intent(out)   :: efield_dofs(:,:)


    ! local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: vi(3), xi(3), wi(1)
    sll_real64 :: qoverm
    sll_real64 :: efield(2)


    ! Set to zero
    self%j_dofs_local = 0.0_f64
    self%particle_mass_1_local = 0.0_f64
    self%particle_mass_0_local = 0.0_f64
    self%rho_dofs_local = 0.0_f64

    ! First particle loop
    do i_sp=1,self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part = 1,self%particle_group%group(i_sp)%n_particles
          vi = self%particle_group%group(i_sp)%get_v(i_part)
          xi = self%particle_group%group(i_sp)%get_x(i_part)
          xi(1) = xi(1) + 0.5_f64 * dt * vi(1)
          xi(1) = modulo(xi(1), self%Lx)
          
          
          ! Get charge for accumulation of j
          wi = self%particle_group%group(i_sp)%get_charge(i_part)
          
          ! Accumulate the particle mass matrix diagonals
          call self%kernel_smoother_1%add_particle_mass( xi(1), &
               wi(1) * (dt**2*0.25_f64) * qoverm, &
               self%particle_mass_1_local )
          call self%kernel_smoother_0%add_particle_mass( xi(1), &
               wi(1) * (dt**2*0.25_f64) * qoverm, &
               self%particle_mass_0_local )
          ! Accumulate jx
          call self%kernel_smoother_1%add_charge( xi(1), wi(1)*vi(1), &
               self%j_dofs_local(:,1) )
          ! Accumulate jy
          call self%kernel_smoother_0%add_charge( xi(1), wi(1)*vi(2), &
               self%j_dofs_local(:,2) )
          ! Evaulate E_x at particle position and propagate v a half step
          call self%kernel_smoother_1%evaluate &
               (xi(1), self%efield_dofs(:,1), efield(1) )
          ! Evaulate E_y at particle position and propagate v a half step
          call self%kernel_smoother_0%evaluate &
               (xi(1), self%efield_dofs(:,2), efield(2))   
          vi(1:2) = vi(1:2) + dt* 0.5_f64* qoverm * efield
          self%vnew(i_sp,:,i_part) = vi(1:2)
          self%xnew(i_sp,i_part) = xi(1)
       
       end do
    end do

    ! Update d_n
    self%j_dofs = 0.0_f64
    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,1), &
         self%kernel_smoother_1%n_dofs, MPI_SUM, self%j_dofs(:,1) )
    call sll_o_collective_allreduce( sll_v_world_collective, self%j_dofs_local(:,2), &
         self%kernel_smoother_1%n_dofs, MPI_SUM, self%j_dofs(:,2) )
    
    self%particle_mass_1 = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_1_local, &
         self%kernel_smoother_1%n_dofs*self%spline_degree, MPI_SUM, self%particle_mass_1)
    
    self%particle_mass_0 = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, self%particle_mass_0_local, &
         self%kernel_smoother_0%n_dofs*(self%spline_degree+1), MPI_SUM, &
         self%particle_mass_0)
    
   
    !print*, 'Solve 1'
    ! Prepare right hand side
    if ( self%smooth .eqv. .false. ) then
       call assemble_schur ( self%kernel_smoother_1%n_dofs, &
            self%spline_degree-1, self%mass_line_1, -self%particle_mass_1, self%schur_1 )
    else
       call assemble_schur_smooth ( self%kernel_smoother_1%n_dofs, &
            self%spline_degree-1, self%mass_line_1, -self%particle_mass_1, self%schur_1 )
    end if
       
    call self%schur_1%dot( self%efield_dofs(:,1) , self%j_dofs_local(:,1) )
    self%j_dofs_local(:,1) = self%j_dofs_local(:,1) - dt * self%j_dofs(:,1)

    
    ! Solve the Schur complement
    if ( self%smooth .eqv. .false. ) then
       call assemble_schur ( self%kernel_smoother_1%n_dofs, &
            self%spline_degree-1, self%mass_line_1, self%particle_mass_1, self%schur_1 )
    else
       call assemble_schur_smooth ( self%kernel_smoother_1%n_dofs, &
            self%spline_degree-1, self%mass_line_1, self%particle_mass_1, self%schur_1 )
    end if
    !call self%preconditioner_1%create( matrix=self%schur_1 )
    call self%linear_solver_1%create( self%schur_1)!, pc_left=self%preconditioner_1 )
    call self%linear_solver_1%set_guess( self%efield_dofs(:,1) )
    self%linear_solver_1%atol = self%solver_tolerance
    call self%linear_solver_1%solve ( self%j_dofs_local(:,1) , efield_dofs(:,1) )
    !call self%preconditioner_1%free()
    call self%linear_solver_1%free()
    
    !print*, 'Solve 2'
    ! Prepare right hand side
    if ( self%smooth .eqv. .false. ) then 
       call assemble_schur ( self%kernel_smoother_0%n_dofs, &
            self%spline_degree, self%mass_line_0, -self%particle_mass_0, self%schur_0 )
    else
       call assemble_schur_smooth ( self%kernel_smoother_0%n_dofs, &
            self%spline_degree, self%mass_line_0, -self%particle_mass_0, self%schur_0 )
    end if
    call self%schur_0%dot( self%efield_dofs(:,2) , self%j_dofs_local(:,2) )
    self%j_dofs_local(:,2) = self%j_dofs_local(:,2) - dt * self%j_dofs(:,2)
    
    ! Solve the Schur complement
    if ( self%smooth .eqv. .false. ) then 
       call assemble_schur ( self%kernel_smoother_0%n_dofs, &
            self%spline_degree, self%mass_line_0, self%particle_mass_0, self%schur_0 )
    else
       call assemble_schur_smooth ( self%kernel_smoother_0%n_dofs, &
            self%spline_degree, self%mass_line_0, self%particle_mass_0, self%schur_0 )
    end if
    !call self%preconditioner_0%create( matrix=self%schur_0 )
    call self%linear_solver_0%create( self%schur_0)!, pc_left=self%preconditioner_0 )
    call self%linear_solver_0%set_guess( self%efield_dofs(:,2) )
    self%linear_solver_0%atol = self%solver_tolerance
    call self%linear_solver_0%solve ( self%j_dofs_local(:,2) , efield_dofs(:,2) )

    !call self%preconditioner_0%free()
    call self%linear_solver_0%free()


    ! Second particle loop (second half step of particle propagation)
    do i_sp=1,self%particle_group%n_species
       qoverm = self%particle_group%group(i_sp)%species%q_over_m();
       do i_part = 1, self%particle_group%group(i_sp)%n_particles
          
          vi(1:2) = self%vnew(i_sp,:,i_part)
          xi(1) = self%xnew(i_sp,i_part)
          ! Evaulate E_x at particle position and propagate v a half step
          call self%kernel_smoother_1%evaluate &
               (xi(1), efield_dofs(:,1), efield(1))
          ! Evaulate E_y at particle position and propagate v a half step
          call self%kernel_smoother_0%evaluate &
               (xi(1), efield_dofs(:,2), efield(2))
          vi(1:2) = vi(1:2) + dt* 0.5_f64* qoverm * efield
          self%vnew(i_sp,:,i_part) = vi(1:2)
          self%xnew(i_sp,i_part) =  modulo(xi(1) + 0.5_f64*dt*vi(1), self%Lx)
       end do
    end do
    

  end subroutine avf_for_start

  

  !> Computes the filtered dofs
  subroutine reinit_fields( self ) 
    class(sll_t_time_propagator_pic_vm_1d2v_helper), intent(inout) :: self !< time splitting object 

    call self%filter%apply( self%efield_dofs(:,1), self%efield_filter(:,1) ) 
    call self%filter%apply( self%efield_dofs(:,2), self%efield_filter(:,2) ) 
    call self%filter%apply( self%bfield_dofs, self%bfield_filter )
    
  end subroutine reinit_fields


end module 
