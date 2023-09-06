module sll_m_sim_6d_utilities
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only : &
       sll_s_collective_allreduce_sum_3d_real64, &
       sll_o_collective_bcast, &
       sll_s_collective_bcast_3d_real64, &
       sll_s_collective_reduce_real64, &
       sll_f_get_collective_rank, &
       sll_v_world_collective, &
       sll_t_collective_t

  use sll_m_constants, only : &
       sll_p_twopi

  use sll_m_dg_interpolator_1d

  use sll_m_distribution_function_initializer_6d, only : &
       sll_t_array

  use sll_m_hdf5_io_parallel
  use sll_m_hdf5_io_serial
  use mpi, only : &
       mpi_bcast, &
       mpi_double_precision, &
       mpi_sum, &
       mpi_success, &
       mpi_sendrecv, &
       mpi_status_ignore

  use sll_m_utilities, only: &
       sll_s_int2string

  use sll_m_remapper, only : &
       sll_t_remap_plan_3d_real64, &
       sll_o_apply_remap_3d

  use sll_m_decomposition, only : &
       sll_t_decomposition_6d, &
       sll_t_decomposition_slim_3d, &
       sll_t_decomposition_slim_6d, &
       sll_t_cartesian_topology_3d, &
       sll_t_cartesian_topology_6d, &
       sll_f_apply_halo_exchange_slim_3d_real64,&
       sll_f_apply_halo_exchange,&
       sll_f_apply_halo_exchange_slim_6d_real64,&
       sll_f_new_cartesian_domain_decomposition_slim_6d

  use sll_m_decomposition_advanced

  
  use sll_m_advection_6d_lagrange_dd_slim , only : &
       sll_t_advection_6d_lagrange_dd_slim, &
       sll_s_advection_6d_lagrange_dd_slim_init, &
       sll_s_advection_6d_lagrange_dd_slim_free, &
       sll_s_advection_6d_clagrange_dd_slim_advect_eta1_d45, &
       sll_s_advection_6d_clagrange_dd_slim_advect_eta2_d45, &
       sll_s_advection_6d_lagrange_dd_slim_fadvect_eta3, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta1_dispeta45, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta2_dispeta45, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta1, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta2, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta3, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta4, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta5, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta6, &
       sll_s_advection_6d_lagrange_dd_slim_set_eta123, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta1_givenv, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta2_givenv, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta3_givenv
  
  use sll_m_timer

#ifdef _OPENMP
  use omp_lib
#define OMP_COLLAPSE collapse(2)
#define OMP_SCHEDULE schedule(static)
#endif

#ifdef USE_HALO_REAL32
#define HALO_DTYPE sll_real32
#else
#define HALO_DTYPE sll_real64
#endif

  implicit none

  public ::  sll_s_compute_charge_density_6d_dd_slim, &
     !  sll_s_plot_distribution_time_6d, &
       sll_s_time_history_diagnostics, &
     !  sll_s_time_history_diagnostics_dg, &
     !  sll_s_time_history_diagnostics_transp, &
     !  sll_s_additional_time_history_diagnostics, &
     !  sll_s_read_distribution_6d, &
       sll_s_check_diagnostics, &
     !  sll_s_compute_phi_qn, &
     !  sll_s_compute_momentum_energy_6d, &
     !  sll_s_plot_diagnostic_time, &
     !  sll_s_plot_diagnostic_time_dd, &
     !  sll_s_compute_free_energy_6d, &
     !  sll_s_plot_v_diagnostic_time, &
       sll_s_write_simulation_info, &
     !  sll_s_quadrature_6dto3d, &
     !  sll_s_uniform_to_gauss, &
     !  sll_s_uniform_to_gauss_staggered, &
     !  sll_s_field_uniform_to_gauss, &
       sll_f_check_triggered_shutdown, &
       sll_t_clocks, sll_t_stopwatch, &
       sll_s_init_clocks, sll_s_finalize_clocks, &
       sll_s_start_clock, sll_s_stop_clock!, &
     !  sll_s_add_source_term, &
     !  sll_s_collision_crook, & 
     !  sll_s_set_collision_mask_v, &
     !  sll_s_set_source_term_v, &
     !  sll_s_compute_momentum_energy_6d_dd ,&
     !  sll_s_additional_time_history_diagnostics_dd,&
     !  sll_s_plot_v_diagnostic_time_dd,&
     !  sll_s_add_noise, &
     !  advect_background,&
     !  sll_s_double_dim_distribution_6d,&
     !  sll_s_half_dim_distribution_6d, &
     !  sll_s_read_distribution_6d_interpolate,&
     !  sll_s_set_f0_v,&
     !  sll_s_plot_diagnostic_time_dd_add,&
     !  sll_s_compute_heatflux_stress_6d_dd

  private

  type :: sll_t_stopwatch
    type(sll_t_time_mark) :: t
    sll_real64 :: elapsed
  end type sll_t_stopwatch

  type :: sll_t_clocks
    type(sll_t_stopwatch) :: slot(ichar('/'):ichar('Z'),&
                                  ichar('/'):ichar('Z'))  ! cf. ASCII table
  end type sll_t_clocks


contains


  ! PROTOTYPE : Compute charge density from a domain-decomposed distribution function.
  subroutine sll_s_compute_charge_density_6d_dd_slim(f6d, decomp_6d, rho_3d, coll_3d_v, volume_v)
    type(sll_t_decomposition_slim_6d), target, intent(in) :: decomp_6d
    sll_real64, intent(in)  :: f6d(decomp_6d%local%mn(1):decomp_6d%local%mx(1), &
                                   decomp_6d%local%mn(2):decomp_6d%local%mx(2), &
                                   decomp_6d%local%mn(3):decomp_6d%local%mx(3), &
                                   decomp_6d%local%mn(4):decomp_6d%local%mx(4), &
                                   decomp_6d%local%mn(5):decomp_6d%local%mx(5), &
                                   decomp_6d%local%mn(6):decomp_6d%local%mx(6))
    type(sll_t_collective_t), pointer, intent(in) :: coll_3d_v
    sll_real64, intent(in)  :: volume_v
    sll_real64, intent(out) :: rho_3d(decomp_6d%local%mn(1):decomp_6d%local%mx(1), &
                                      decomp_6d%local%mn(2):decomp_6d%local%mx(2), &
                                      decomp_6d%local%mn(3):decomp_6d%local%mx(3))
!    sll_real64, allocatable :: sm(:)
!    sll_int32 :: i, j, k, l, m, n
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    loop_mn => decomp_6d%local%mn
    loop_mx => decomp_6d%local%mx

!  !$omp parallel default(shared) private(i,j,k,l,m,n,sm)
!  !   allocate(sm(loop_mn(1):loop_mx(1)))
!  !$omp do OMP_COLLAPSE OMP_SCHEDULE
!  !   do k=loop_mn(3),loop_mx(3)
!  !      do j=loop_mn(2),loop_mx(2)
!  !         ! --- cache blocking in i by using an sm array ---
!  !         sm(:) = 0.0_f64
!  !         do n=loop_mn(6),loop_mx(6)
!  !           do m=loop_mn(5),loop_mx(5)
!  !             do l=loop_mn(4),loop_mx(4)
!  !               do i=loop_mn(1),loop_mx(1)
!  !                 sm(i) = sm(i) + f6d(i,j,k,l,m,n)
!  !               end do
!  !             enddo
!  !           enddo
!  !         enddo
!  !         do i=loop_mn(1),loop_mx(1)
!  !           rho_3d(i,j,k) = sm(i) * volume_v
!  !         enddo
!  !      end do
!  !   end do
!  !$omp end do
!  !   deallocate(sm)
!  !$omp end parallel
    call sll_s_compute_charge_density_6d_core(f6d, loop_mn, loop_mx, &
                                              rho_3d, loop_mn, loop_mx, &
                                              loop_mn, loop_mx, volume_v)

    ! --- use an allreduce operation until we know where the result needs to go
    call sll_s_collective_allreduce_sum_3d_real64(coll_3d_v, rho_3d)
  end subroutine sll_s_compute_charge_density_6d_dd_slim



  ! Charge density computation core routine used by the implementation-specific callers above,
  ! the index arrays are necessary to map remap, dd, dd_slim layouts to the same code.
  ! Routine uses cache blocking which is greatly improving the performance.
  subroutine sll_s_compute_charge_density_6d_core(f6d, f6d_idx_lo, f6d_idx_hi, &
                                                  rho, rho_idx_lo, rho_idx_hi, &
                                                  loop_mn, loop_mx, volume_v)
    sll_int32, intent(in)   :: f6d_idx_lo(6), f6d_idx_hi(6), &
                               rho_idx_lo(:), rho_idx_hi(:), &
                               loop_mn(6), loop_mx(6)
    sll_real64, intent(in)  :: f6d(f6d_idx_lo(1):f6d_idx_hi(1), &
                                   f6d_idx_lo(2):f6d_idx_hi(2), &
                                   f6d_idx_lo(3):f6d_idx_hi(3), &
                                   f6d_idx_lo(4):f6d_idx_hi(4), &
                                   f6d_idx_lo(5):f6d_idx_hi(5), &
                                   f6d_idx_lo(6):f6d_idx_hi(6))
    sll_real64, intent(out) :: rho(rho_idx_lo(1):rho_idx_hi(1), &
                                   rho_idx_lo(2):rho_idx_hi(2), &
                                   rho_idx_lo(3):rho_idx_hi(3))
    sll_real64, intent(in)  :: volume_v
    sll_real64, allocatable :: sm(:)
    sll_int32 :: i, j, k, l, m, n
!$omp parallel default(shared) private(i,j,k,l,m,n,sm)
    allocate(sm(loop_mn(1):loop_mx(1)))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k=loop_mn(3),loop_mx(3)
       do j=loop_mn(2),loop_mx(2)
          ! --- cache blocking in i ---
          sm(:) = 0.0_f64
          do n=loop_mn(6),loop_mx(6)
            do m=loop_mn(5),loop_mx(5)
              do l=loop_mn(4),loop_mx(4)
                do i=loop_mn(1),loop_mx(1)
                  sm(i) = sm(i) - f6d(i,j,k,l,m,n)
                end do
              enddo
            enddo
          enddo
          do i=loop_mn(1),loop_mx(1)
            rho(i,j,k) = sm(i) * volume_v
          enddo
       end do
    end do
!$omp end do
    deallocate(sm)
!$omp end parallel
  end subroutine sll_s_compute_charge_density_6d_core



  subroutine sll_s_time_history_diagnostics(&
       index_min, &
       index_max, &
       data_index_min, &
       time, &
       dV, &
       dV_x, &
       tensor_grid,&
       f6d, &
       rho, &
       phi, &
       ex, &
       ey, &
       ez, &
       file_id, &
       topology_3d_velocity_coords, &
       nan_blowup)
    sll_int32,         intent(in) :: index_min(6)
    sll_int32,         intent(in) :: index_max(6)
    sll_int32,         intent(in) :: data_index_min(6)
    sll_real64,        intent(in) :: time
    sll_real64,        intent(in) :: dV
    sll_real64,        intent(in) :: dV_x
    type(sll_t_array), intent(in) :: tensor_grid(6)
    sll_real64,        intent(in) :: f6d(data_index_min(1):,data_index_min(2):,&
                                         data_index_min(3):,data_index_min(4):,&
                                         data_index_min(5):,data_index_min(6):)
    sll_real64,        intent(in) :: rho(:,:,:)
    sll_real64,        intent(in) :: phi(:,:,:)
    sll_real64,        intent(in) :: ex(:,:,:)
    sll_real64,        intent(in) :: ey(:,:,:)
    sll_real64,        intent(in) :: ez(:,:,:)
    sll_int32,         intent(in) :: file_id
    sll_int32,optional,intent(in) :: topology_3d_velocity_coords(3)
    logical, optional, intent(out) :: nan_blowup

    sll_real64 :: diagn_data(13), sm, diagn_data_all(13)
    sll_real64 :: smm, smsq, smp4, smk4, smp5, smk5, smp6, smk6
    sll_int32 :: i, j, k, l, m, n
    sll_int32 :: i_mx, j_mx, k_mx
    sll_int32 :: i_phi_mx, j_phi_mx, k_phi_mx
    logical :: do_electrostatics

    ! --- check if the present MPI rank is allowed to do the electrostatic calculations
    if (present(topology_3d_velocity_coords)) then
      if (all(topology_3d_velocity_coords == 0)) then
        do_electrostatics = .true.
      else
        do_electrostatics = .false.
      end if
    else
      do_electrostatics = .true.
    end if

    diagn_data = 0.0_f64

    ! --- calculate the maximum indices of the rho, ex, ey, ez arrays
    ! (plausibly assuming they are always the same)
    i_mx = size(rho, 1)
    j_mx = size(rho, 2)
    k_mx = size(rho, 3)
    ! --- phi has different extents
    i_phi_mx = size(phi, 1)
    j_phi_mx = size(phi, 2)
    k_phi_mx = size(phi, 3)

    diagn_data(:) = 0.0_f64
!$omp parallel default(shared) private(sm,smsq,smm,smp4,smp5,smp6,smk4,smk5,smk6,i,j,k,l,m,n) reduction(+ : diagn_data)

!    ! mass
!    diagn_data(1) = sum(f6d(index_min(1):index_max(1), index_min(2):index_max(2), &
!                            index_min(3):index_max(3), index_min(4):index_max(4), &
!                            index_min(5):index_max(5), index_min(6):index_max(6)))&
!                  * dV
    sm = 0.0_f64
    smsq = 0.0_f64
    smp4 = 0.0_f64
    smp5 = 0.0_f64
    smp6 = 0.0_f64
    smk4 = 0.0_f64
    smk5 = 0.0_f64
    smk6 = 0.0_f64
    smm = 0.0_f64
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=index_min(6),index_max(6)
      do m=index_min(5),index_max(5)
         do l=index_min(4),index_max(4)
            sm = 0.0_f64
            do k=index_min(3),index_max(3)
               do j=index_min(2),index_max(2)
                  do i=index_min(1),index_max(1)
                     sm = sm + f6d(i,j,k,l,m,n)
                     smsq = smsq + f6d(i,j,k,l,m,n)**2
                  enddo
               enddo
            enddo
            smp4 = smp4 + sm *  tensor_grid(4)%vals(l-index_min(4)+1) * dV
            smp5 = smp5 + sm * tensor_grid(5)%vals(m-index_min(5)+1) * dV
            smp6 = smp6 + sm * tensor_grid(6)%vals(n-index_min(6)+1) * dV
            smk4 = smk4 + sm * tensor_grid(4)%vals(l-index_min(4)+1)**2 * dV
            smk5 = smk5 + sm * tensor_grid(5)%vals(m-index_min(5)+1)**2 * dV
            smk6 = smk6 + sm * tensor_grid(6)%vals(n-index_min(6)+1)**2 * dV

            smm = smm + sm
        enddo
      enddo
    enddo
!$omp end do nowait
    diagn_data(1) = smm * dV
    diagn_data(2) = smsq * dV
    diagn_data(8) = smp4
    diagn_data(9) = smp5
    diagn_data(10) = smp6
    diagn_data(11) = smk4
    diagn_data(12) = smk5
    diagn_data(13) = smk6

!    ! l2 norm
!    diagn_data(2) = sum(f6d(index_min(1):index_max(1), index_min(2):index_max(2), &
!                            index_min(3):index_max(3), index_min(4):index_max(4), &
!                            index_min(5):index_max(5), index_min(6):index_max(6))**2)&
!                  * dV
!!$    sm = 0.0_f64
!!$!$omp do OMP_COLLAPSE OMP_SCHEDULE
!!$    do n=index_min(6),index_max(6)
!!$      do m=index_min(5),index_max(5)
!!$        do l=index_min(4),index_max(4)
!!$          do k=index_min(3),index_max(3)
!!$            do j=index_min(2),index_max(2)
!!$              do i=index_min(1),index_max(1)
!!$                sm = sm + f6d(i,j,k,l,m,n)**2
!!$              enddo
!!$            enddo
!!$          enddo
!!$        enddo
!!$      enddo
!!$    enddo
!!$!$omp end do nowait
!!$    diagn_data(2) = sm * dV


    if (do_electrostatics) then
!       ! l2 norm of rho
!       diagn_data(3) = sum(rho**2)*dV_x
      sm = 0.0_f64
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,k_mx
        do j=1,j_mx
          do i=1,i_mx
            sm = sm + rho(i,j,k)**2
          enddo
        enddo
      enddo
!$omp end do nowait
      diagn_data(3) = sm * dV_x

!       ! l2 norm of phi
!     diagn_data(4) = sum(phi**2)*dV_x
      sm = 0.0_f64
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,k_phi_mx
        do j=1,j_phi_mx
          do i=1,i_phi_mx
            sm = sm + phi(i,j,k)**2
          enddo
        enddo
      enddo
!$omp end do nowait
      diagn_data(4) = sm * dV_x

!       ! l2 norm of ex
!       diagn_data(5) = sum(ex**2)*dV_x
      sm = 0.0_f64
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,k_mx
        do j=1,j_mx
          do i=1,i_mx
            sm = sm + ex(i,j,k)**2
          enddo
        enddo
      enddo
!$omp end do nowait
      diagn_data(5) = sm * dV_x

!       ! l2 norm of ey
!       diagn_data(6) = sum(ey**2)*dV_x
      sm = 0.0_f64
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,k_mx
        do j=1,j_mx
          do i=1,i_mx
            sm = sm + ey(i,j,k)**2
          enddo
        enddo
      enddo
!$omp end do nowait
      diagn_data(6) = sm * dV_x

!       ! l2 norm of ez
!       diagn_data(7) = sum(ez**2)*dV_x
      sm = 0.0_f64
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,k_mx
        do j=1,j_mx
          do i=1,i_mx
            sm = sm + ez(i,j,k)**2
          enddo
        enddo
      enddo
!$omp end do nowait
      diagn_data(7) = sm * dV_x
    endif  ! do_electrostatics

!    sm = 0.0_f64
!    do j=index_min(4), index_max(4)
!       sm = sm + &
!            sum(f6d(index_min(1):index_max(1), index_min(2):index_max(2), &
!                    index_min(3):index_max(3), j, &
!                    index_min(5):index_max(5), index_min(6):index_max(6)) &
!                * tensor_grid(4)%vals(j-index_min(4)+1))
!    end do
!    diagn_data(8) = sm*dV
!!$!$omp do OMP_SCHEDULE
!!$    do l=index_min(4),index_max(4)
!!$      sm = 0.0_f64
!!$      do n=index_min(6),index_max(6)
!!$        do m=index_min(5),index_max(5)
!!$          do k=index_min(3),index_max(3)
!!$            do j=index_min(2),index_max(2)
!!$              do i=index_min(1),index_max(1)
!!$                sm = sm + f6d(i,j,k,l,m,n)
!!$              enddo
!!$            enddo
!!$          enddo
!!$        enddo
!!$      enddo
!!$      diagn_data(8) = diagn_data(8) + sm * tensor_grid(4)%vals(l-index_min(4)+1) * dV
!!$    enddo
!!$!$omp end do nowait

!    sm = 0.0_f64
!    do j=index_min(5), index_max(5)
!       sm = sm + &
!            sum(f6d(index_min(1):index_max(1), index_min(2):index_max(2), &
!                    index_min(3):index_max(3), index_min(4):index_max(4), &
!                    j, index_min(6):index_max(6)) &
!                * tensor_grid(5)%vals(j-index_min(5)+1))
!    end do
!    diagn_data(9) = sm*dV
!!$!$omp do OMP_SCHEDULE
!!$    do m=index_min(5),index_max(5)
!!$      sm = 0.0_f64
!!$      do n=index_min(6),index_max(6)
!!$        do l=index_min(4),index_max(4)
!!$          do k=index_min(3),index_max(3)
!!$            do j=index_min(2),index_max(2)
!!$              do i=index_min(1),index_max(1)
!!$                sm = sm + f6d(i,j,k,l,m,n)
!!$              enddo
!!$            enddo
!!$          enddo
!!$        enddo
!!$      enddo
!!$      diagn_data(9) = diagn_data(9) + sm * tensor_grid(5)%vals(m-index_min(5)+1) * dV
!!$    enddo
!!$!$omp end do nowait

!    sm = 0.0_f64
!    do j=index_min(6), index_max(6)
!       sm = sm + &
!            sum(f6d(index_min(1):index_max(1), index_min(2):index_max(2), &
!                    index_min(3):index_max(3), index_min(4):index_max(4), &
!                    index_min(5):index_max(5), j) &
!                * tensor_grid(6)%vals(j-index_min(6)+1))
!    end do
!    diagn_data(10) = sm*dV
!!$!$omp do OMP_SCHEDULE
!!$    do n=index_min(6),index_max(6)
!!$      sm = 0.0_f64
!!$      do m=index_min(5),index_max(5)
!!$        do l=index_min(4),index_max(4)
!!$          do k=index_min(3),index_max(3)
!!$            do j=index_min(2),index_max(2)
!!$              do i=index_min(1),index_max(1)
!!$                sm = sm + f6d(i,j,k,l,m,n)
!!$              enddo
!!$            enddo
!!$          enddo
!!$        enddo
!!$      enddo
!!$      diagn_data(10) = diagn_data(10) + sm * tensor_grid(6)%vals(n-index_min(6)+1) * dV
!!$    enddo
!!$!$omp end do nowait

!    sm = 0.0_f64
!    do j=index_min(4), index_max(4)
!       sm = sm + &
!            sum(f6d(index_min(1):index_max(1), index_min(2):index_max(2), &
!                    index_min(3):index_max(3), j, &
!                    index_min(5):index_max(5), index_min(6):index_max(6)) &
!                * tensor_grid(4)%vals(j-index_min(4)+1)**2)
!    end do
!    diagn_data(11) = sm*dV
!!$!$omp do OMP_SCHEDULE
!!$    do l=index_min(4),index_max(4)
!!$      sm = 0.0_f64
!!$      do n=index_min(6),index_max(6)
!!$        do m=index_min(5),index_max(5)
!!$          do k=index_min(3),index_max(3)
!!$            do j=index_min(2),index_max(2)
!!$              do i=index_min(1),index_max(1)
!!$                sm = sm + f6d(i,j,k,l,m,n)
!!$              enddo
!!$            enddo
!!$          enddo
!!$        enddo
!!$      enddo
!!$      diagn_data(11) = diagn_data(11) + sm * tensor_grid(4)%vals(l-index_min(4)+1)**2 * dV
!!$    enddo
!!$!$omp end do nowait

!    sm = 0.0_f64
!    do j=index_min(5), index_max(5)
!       sm = sm + &
!            sum(f6d(index_min(1):index_max(1), index_min(2):index_max(2), &
!                    index_min(3):index_max(3), index_min(4):index_max(4), &
!                    j, index_min(6):index_max(6)) &
!                * tensor_grid(5)%vals(j-index_min(5)+1)**2)
!    end do
!    diagn_data(12) = sm*dV
!!$!$omp do OMP_SCHEDULE
!!$    do m=index_min(5),index_max(5)
!!$      sm = 0.0_f64
!!$      do n=index_min(6),index_max(6)
!!$        do l=index_min(4),index_max(4)
!!$          do k=index_min(3),index_max(3)
!!$            do j=index_min(2),index_max(2)
!!$              do i=index_min(1),index_max(1)
!!$                sm = sm + f6d(i,j,k,l,m,n)
!!$              enddo
!!$            enddo
!!$          enddo
!!$        enddo
!!$      enddo
!!$      diagn_data(12) = diagn_data(12) + sm * tensor_grid(5)%vals(m-index_min(5)+1)**2 * dV
!!$    enddo
!!$!$omp end do nowait
!!$
!!$!    sm = 0.0_f64
!!$!    do j=index_min(6), index_max(6)
!!$!       sm = sm + &
!!$!            sum(f6d(index_min(1):index_max(1), index_min(2):index_max(2), &
!!$!                    index_min(3):index_max(3), index_min(4):index_max(4), &
!!$!                    index_min(5):index_max(5), j) &
!!$!                * tensor_grid(6)%vals(j-index_min(6)+1)**2)
!!$!    end do
!!$!    diagn_data(13) = sm*dV
!!$!$omp do OMP_SCHEDULE
!!$    do n=index_min(6),index_max(6)
!!$      sm = 0.0_f64
!!$      do m=index_min(5),index_max(5)
!!$        do l=index_min(4),index_max(4)
!!$          do k=index_min(3),index_max(3)
!!$            do j=index_min(2),index_max(2)
!!$              do i=index_min(1),index_max(1)
!!$                sm = sm + f6d(i,j,k,l,m,n)
!!$              enddo
!!$            enddo
!!$          enddo
!!$        enddo
!!$      enddo
!!$      diagn_data(13) = diagn_data(13) + sm * tensor_grid(6)%vals(n-index_min(6)+1)**2 * dV
!!$    enddo
!!$!$omp end do nowait

!$omp end parallel
    ! Implicit summation over diagn_data at the end of the parallel region.
    ! --> The array now contains the sum of the values over all threads.

    diagn_data_all = 0.0_f64
    call sll_s_collective_reduce_real64(sll_v_world_collective, diagn_data, 13, MPI_SUM, 0, diagn_data_all)

    if(sll_f_get_collective_rank(sll_v_world_collective)==0) then
      write(file_id,'(e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12)') &!'(f12.5,2g20.12)') &
        time, diagn_data_all
    endif

    ! check if NaNs exist in the diagnostics array
    if (present(nan_blowup)) then
      if (all((diagn_data_all == diagn_data_all) .eqv. .true.)) then
        nan_blowup = .false.
      else
        nan_blowup = .true.
      endif
    endif
  end subroutine sll_s_time_history_diagnostics



  subroutine sll_s_check_diagnostics(reffile, simfile, srcfile, tolerated_error)
    character(*), intent(in) :: reffile !< Name of reference file (stored in same folder as source file)
    character(*), intent(in) :: simfile !< Name of file with simulation results
    ! --- legacy argument, unused and to be removed ---
    character(*), intent(in), optional :: srcfile
    sll_real64, intent(in), optional :: tolerated_error

    sll_real64 :: error, tol_err
    sll_real64 :: data_sim(3,14)
    sll_real64 :: data_ref(3,14)
    sll_int32  :: file_id

    if (present(tolerated_error)) then
      tol_err = tolerated_error
    else
      tol_err = 5.0d-7
    endif

    ! Read simulation result
    open(newunit=file_id, file=simfile, status='old', action='read')
    read(unit=file_id,fmt=*) data_sim
    close(file_id)

    ! Read reference
    open(newunit=file_id, file=reffile, status='old', action='read')
    read(unit=file_id,fmt=*) data_ref
    close(file_id)

    ! Compare
    data_sim = data_sim - data_ref
    error = maxval(abs(data_sim))

    print*, 'Max error in time history diagnostics: ', error
    if (error < tol_err) then
       print*, 'PASSED.'
    else
       print*, 'FAILED.'
       stop
    end if
  end subroutine sll_s_check_diagnostics
  

  subroutine sll_s_write_simulation_info( file_prefix, nml_filename, &
       num_cells, time_step, &
       num_mpi_procs, wall_time )
    character(len=*), intent( in    ) :: file_prefix
    character(len=*), intent( in    ) :: nml_filename
    sll_int32,        intent( in    ) :: num_cells(:)
    sll_real64,       intent( in    ) :: time_step(:)
    sll_int32,        intent( in    ) :: num_mpi_procs(:)
    sll_real64,       intent( in    ) :: wall_time(:)

    type(sll_t_hdf5_ser_handle) :: file_handle
    sll_int32      :: ierr
    sll_int32      :: num_omp_procs(1)


    call sll_s_hdf5_ser_file_create( trim(file_prefix)//'.h5', &
         file_handle, ierr )

    call sll_o_hdf5_ser_write_array( file_handle, num_cells, '/num_cells',&
         ierr)
    call sll_o_hdf5_ser_write_array( file_handle, time_step, '/time_step',&
         ierr)
    call sll_o_hdf5_ser_write_array( file_handle, num_mpi_procs, '/num_mpi_procs', &
         ierr)
#ifdef _OPENMP
        num_omp_procs = omp_get_max_threads()
#else
        num_omp_procs = 1
#endif

    call sll_o_hdf5_ser_write_array( file_handle, num_omp_procs, '/num_openmp_procs', &
         ierr)
    call sll_o_hdf5_ser_write_array( file_handle, wall_time, '/wall_time',&
            ierr)
    call sll_s_hdf5_ser_write_file( file_handle, trim(nml_filename),  '/input_file', ierr)
    call sll_s_hdf5_ser_file_close( file_handle, ierr )
  end subroutine sll_s_write_simulation_info



  !> @brief  Check if a file named "stop" exists and communicate the result to all
  !>         MPI processes. Finally, remove the file.  Return if the file existed.
  function sll_f_check_triggered_shutdown()
    logical :: sll_f_check_triggered_shutdown
    logical :: file_exists(1)
    character(len=*), parameter :: file_name = "stop"

    if (sll_f_get_collective_rank(sll_v_world_collective)==0) then
      inquire(file="stop", exist=file_exists(1))
    endif
    call sll_o_collective_bcast(sll_v_world_collective, file_exists, 0)
    if (file_exists(1)) then
      sll_f_check_triggered_shutdown = .true.
    else
      sll_f_check_triggered_shutdown = .false.
    endif
    if ((sll_f_get_collective_rank(sll_v_world_collective)==0) .and. (sll_f_check_triggered_shutdown)) then
      write(*,*) "Shutdown trigger file 'stop' detected."
      call sll_s_rm(file_name)
    endif
  end function sll_f_check_triggered_shutdown


  !> @brief  Delete a file from the file system.
  subroutine sll_s_rm(file_name)
    character(len=*), intent(in) :: file_name
    integer, parameter :: fp = 4711
    integer :: stat
    open(unit=fp, iostat=stat, file=file_name, status='old')
    if (stat == 0) then
      close(fp, status='delete')
    endif
  end subroutine sll_s_rm


  subroutine sll_s_init_clocks(clocks)
    type(sll_t_clocks) :: clocks
    integer :: i, j
    do i = ichar('/'), ichar('Z')
      do j = ichar('/'), ichar('Z')
        clocks%slot(i,j)%elapsed = 0.0_f64
      enddo
    enddo
  end subroutine sll_s_init_clocks

  subroutine sll_s_finalize_clocks(clocks)
    type(sll_t_clocks) :: clocks
    integer :: i, j
    integer, parameter :: fd = 68
    character(len=*), parameter :: filename = "sll_clocks.txt"
    character(len=2) :: label
    open(unit=fd, file=filename, action='write')
    do i = ichar('/'), ichar('Z')
      do j = ichar('/'), ichar('Z')
        if (clocks%slot(i,j)%elapsed > 0.0_f64) then
          if (j == ichar('/')) then
            label = char(i)
          else
            label = char(i)//char(j)
          endif
          write(fd,*) trim(label), clocks%slot(i,j)%elapsed
        endif
      enddo
    enddo
    close(fd)
  end subroutine sll_s_finalize_clocks

  subroutine sll_s_start_clock(clocks, label)
    type(sll_t_clocks) :: clocks
    character(len=*) :: label
    integer :: key_1, key_2
    if (len_trim(label) == 1) then
      key_1 = ichar(label(1:1))
      key_2 = ichar('/')
    else
      key_1 = ichar(label(1:1))
      key_2 = ichar(label(2:2))
    endif
    call sll_s_set_time_mark(clocks%slot(key_1, key_2)%t)
  end subroutine sll_s_start_clock

  subroutine sll_s_stop_clock(clocks, label)
    type(sll_t_clocks) :: clocks
    character(len=*) :: label
    type(sll_t_time_mark) :: t
    integer :: key_1, key_2
    call sll_s_set_time_mark(t)
    if (len_trim(label) == 1) then
      key_1 = ichar(label(1:1))
      key_2 = ichar('/')
    else
      key_1 = ichar(label(1:1))
      key_2 = ichar(label(2:2))
    endif
    clocks%slot(key_1, key_2)%elapsed = clocks%slot(key_1, key_2)%elapsed + &
             sll_f_time_elapsed_between(clocks%slot(key_1, key_2)%t, t)
  end subroutine sll_s_stop_clock

end module sll_m_sim_6d_utilities
