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
  use sll_mpi, only : &
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

  public ::  sll_s_compute_charge_density_6d, &
       sll_s_compute_charge_density_6d_transp, &
       sll_s_compute_charge_density_6d_dd, &
       sll_s_compute_charge_density_6d_dd_slim, &
       sll_s_compute_charge_density_6d_dd_slim_dg, &
       sll_s_compute_charge_density_6d_dd_slim_dg_staggered, &
       sll_s_compute_charge_density_6d_dd_slim_overlap, &
       sll_s_plot_distribution_time_6d, &
       sll_s_time_history_diagnostics, &
       sll_s_time_history_diagnostics_dg, &
       sll_s_time_history_diagnostics_transp, &
       sll_s_additional_time_history_diagnostics, &
       sll_s_read_distribution_6d, &
       sll_s_check_diagnostics, &
       sll_s_compute_phi_qn, &
       sll_s_compute_momentum_energy_6d, &
       sll_s_plot_diagnostic_time, &
       sll_s_plot_diagnostic_time_dd, &
       sll_s_compute_free_energy_6d, &
       sll_s_plot_v_diagnostic_time, &
       sll_s_write_simulation_info, &
       sll_s_quadrature_6dto3d, &
       sll_s_uniform_to_gauss, &
       sll_s_uniform_to_gauss_staggered, &
       sll_s_field_uniform_to_gauss, &
       sll_f_check_triggered_shutdown, &
       sll_t_clocks, sll_t_stopwatch, &
       sll_s_init_clocks, sll_s_finalize_clocks, &
       sll_s_start_clock, sll_s_stop_clock, &
       sll_s_add_source_term, &
       sll_s_collision_crook, & 
       sll_s_set_collision_mask_v, &
       sll_s_set_source_term_v, &
       sll_s_compute_momentum_energy_6d_dd ,&
       sll_s_additional_time_history_diagnostics_dd,&
       sll_s_plot_v_diagnostic_time_dd,&
       sll_s_add_noise, &
       advect_background,&
       sll_s_double_dim_distribution_6d,&
       sll_s_half_dim_distribution_6d, &
       sll_s_read_distribution_6d_interpolate,&
       sll_s_set_f0_v,&
       sll_s_plot_diagnostic_time_dd_add,&
       sll_s_compute_heatflux_stress_6d_dd

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


  !> Compute charge density from a domain-decomposed distribution function.
  subroutine sll_s_compute_charge_density_6d_dd_slim_dg(f6d, decomp_6d, &
                rho_3d, coll_3d_v, topo_3d, decomp_3d, &
                qweights4, qweights5, qweights6, eval_dg_to_equi, degree)
    sll_real64, intent(in) :: f6d(:,:,:,:,:,:)
    type(sll_t_decomposition_slim_6d), target, intent(in) :: decomp_6d
    sll_real64, intent(out) :: rho_3d(:,:,:)
    type(sll_t_collective_t), pointer, intent(in) :: coll_3d_v
    type(sll_t_cartesian_topology_3d), target, intent(in) :: topo_3d
    type(sll_t_decomposition_slim_3d), target, intent(inout) :: decomp_3d
    sll_real64, intent( in ) :: qweights4(:)
    sll_real64, intent( in ) :: qweights5(:)
    sll_real64, intent( in ) :: qweights6(:)
    sll_real64, intent( in ) :: eval_dg_to_equi(:,:,:)
    sll_int32, intent( in ) :: degree(6)

    sll_int32 :: i, j, k, id, deg
    sll_real64 :: rho_old
    sll_real64, allocatable :: rho_cell(:)
    sll_real64, allocatable :: rho_old_1_tx(:,:), rho_old_2_tx(:,:), rho_old_3_tx(:,:)
    sll_real64, allocatable :: rho_old_1_rx(:,:), rho_old_2_rx(:,:), rho_old_3_rx(:,:)

    ! --- local index copies introduced during the transition to a domain decomposition
    sll_int32, pointer :: n_points(:), n_cells(:)

    n_points => decomp_6d%local%nw
    n_cells => decomp_6d%local%n_cells

    ! Quadrature over local parts of dimensions 4-6
    call sll_s_quadrature_6dto3d( f6d, n_points, qweights4, qweights5, qweights6, rho_3d )
    ! --- use an allreduce operation until we know where the result needs to go
    call sll_s_collective_allreduce_sum_3d_real64(coll_3d_v, rho_3d)

    ! Interpolate to uniform mesh

    i=n_points(1)
    j=n_points(2)
    k=n_points(3)
    allocate(rho_old_1_tx(j,k)); allocate(rho_old_1_rx(j,k))
    allocate(rho_old_2_tx(i,k)); allocate(rho_old_2_rx(i,k))
    allocate(rho_old_3_tx(i,j)); allocate(rho_old_3_rx(i,j))


    deg = degree(1)
    i = n_cells(1)
    id = deg+1
    do k=1,n_points(3)
       do j=1,n_points(2)
          rho_old_1_tx(j,k) = sum(eval_dg_to_equi(:,id,1) * rho_3d((i-1)*deg+1:i*deg,j,k))
       end do
    end do
    ! --- send rho_old_1(:,:) to the right neighbour
    call sll_s_exchange_charge_density_boundary_value_dd_slim_dg(topo_3d, decomp_3d, rho_old_1_tx, rho_old_1_rx, 1)

!$omp parallel default(shared) private(i, j, k, id, rho_old, rho_cell)
    allocate(rho_cell(deg+1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k=1,n_points(3)
       do j=1,n_points(2)
          rho_old = rho_old_1_rx(j,k)
          do i=1,n_cells(1)
             do id=1,deg+1
                rho_cell(id) = sum(eval_dg_to_equi(:,id,1) * rho_3d((i-1)*deg+1:i*deg,j,k))
             end do
             rho_3d((i-1)*deg+1,j,k) = 0.5_f64*(rho_cell(1)+rho_old)
             rho_3d((i-1)*deg+2:i*deg,j,k) = rho_cell(2:deg)
             rho_old = rho_cell(deg+1)
          end do
       end do
    end do
!$omp end do
    deallocate(rho_cell)
!$omp end parallel


    deg = degree(2)
    j = n_cells(2)
    id = deg+1
    do k=1,n_points(3)
       do i=1,n_points(1)
          rho_old_2_tx(i,k) = sum(eval_dg_to_equi(:,id,2) * rho_3d(i,(j-1)*deg+1:j*deg,k))
       end do
    end do
    ! --- send rho_old_2(:,:) to the right neighbour
    call sll_s_exchange_charge_density_boundary_value_dd_slim_dg(topo_3d, decomp_3d, rho_old_2_tx, rho_old_2_rx, 2)

!$omp parallel default(shared) private(i, j, k, id, rho_old, rho_cell)
    allocate(rho_cell(deg+1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k=1,n_points(3)
       do i=1,n_points(1)
          rho_old = rho_old_2_rx(i,k)
          do j=1,n_cells(2)
             do id=1,deg+1
                rho_cell(id) = sum(eval_dg_to_equi(:,id,2) * rho_3d(i,(j-1)*deg+1:j*deg,k))
             end do
             rho_3d(i,(j-1)*deg+1,k) = 0.5_f64*(rho_cell(1)+rho_old)
             rho_3d(i,(j-1)*deg+2:j*deg,k) = rho_cell(2:deg)
             rho_old = rho_cell(deg+1)
          end do
       end do
    end do
!$omp end do
    deallocate(rho_cell)
!$omp end parallel


    deg = degree(3)
    k = n_cells(3)
    id = deg+1
    do j=1,n_points(2)
       do i=1,n_points(1)
          rho_old_3_tx(i,j) = sum(eval_dg_to_equi(:,id,3) * rho_3d(i,j,(k-1)*deg+1:k*deg))
       end do
    end do
    ! --- send rho_old_3(:,:) to the right neighbour
    call sll_s_exchange_charge_density_boundary_value_dd_slim_dg(topo_3d, decomp_3d, rho_old_3_tx, rho_old_3_rx, 3)

!$omp parallel default(shared) private(i, j, k, id, rho_old, rho_cell)
    allocate(rho_cell(deg+1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do j=1,n_points(2)
       do i=1,n_points(1)
          rho_old = rho_old_3_rx(i,j)
          do k=1,n_cells(3)
             do id=1,deg+1
                rho_cell(id) = sum(eval_dg_to_equi(:,id,3) * rho_3d(i,j,(k-1)*deg+1:k*deg))
             end do
             rho_3d(i,j,(k-1)*deg+1) = 0.5_f64*(rho_cell(1)+rho_old)
             rho_3d(i,j,(k-1)*deg+2:k*deg) = rho_cell(2:deg)
             rho_old = rho_cell(deg+1)
          end do
       end do
    end do
!$omp end do
    deallocate(rho_cell)
!$omp end parallel

    deallocate(rho_old_1_tx); deallocate(rho_old_1_rx)
    deallocate(rho_old_2_tx); deallocate(rho_old_2_rx)
    deallocate(rho_old_3_tx); deallocate(rho_old_3_rx)
  end subroutine sll_s_compute_charge_density_6d_dd_slim_dg


 !> Compute charge density from a domain-decomposed distribution function.
subroutine sll_s_compute_charge_density_6d_dd_slim_dg_staggered(f6d, decomp_6d, rho_3d, coll_3d_v, &
                                        qweights4, qweights5, qweights6, eval_dg_to_equi, deg)
    type(sll_t_decomposition_slim_6d), target, intent(in) :: decomp_6d
    sll_real64, intent(in)  :: f6d(:,:,:,:,:,:)
    type(sll_t_collective_t), pointer, intent(in) :: coll_3d_v
    sll_real64, intent(out) :: rho_3d(:,:,:)
    sll_real64, intent( in ) :: qweights4(:)
    sll_real64, intent( in ) :: qweights5(:)
    sll_real64, intent( in ) :: qweights6(:)
    sll_real64, intent( in ) :: eval_dg_to_equi(:,:,:)
    sll_int32, intent( in ) :: deg
    sll_int32 :: i,j,k,id
    sll_real64, allocatable :: rho_cell(:)
    sll_int32, pointer :: n_points(:), n_cells(:)

    n_points => decomp_6d%local%nw
    n_cells => decomp_6d%local%n_cells

    ! Quadrature over local parts of dimensions 4-6
    call sll_s_quadrature_6dto3d( f6d, n_points, qweights4, qweights5, qweights6, rho_3d )
    ! --- use an allreduce operation until we know where the result needs to go
    call sll_s_collective_allreduce_sum_3d_real64(coll_3d_v, rho_3d)

    ! Interpolate to uniform mesh

!$omp parallel default(shared) private(i, j, k, id, rho_cell)
    allocate(rho_cell(1:deg))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k=1,n_points(3)
       do j=1,n_points(2)
          do i=1,n_cells(1)
             do id=1,deg
                rho_cell(id) = sum(eval_dg_to_equi(:,id,1) * rho_3d((i-1)*deg+1:i*deg,j,k))
             end do
             rho_3d((i-1)*deg+1:i*deg,j,k) = rho_cell
          end do
       end do
    end do
!$omp end do
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k=1,n_points(3)
       do i=1,n_points(1)
          do j=1,n_cells(2)
             do id=1,deg
                rho_cell(id) = sum(eval_dg_to_equi(:,id,2) * rho_3d(i,(j-1)*deg+1:j*deg,k))
             end do
             rho_3d(i,(j-1)*deg+1:j*deg,k) = rho_cell
          end do
       end do
    end do
!$omp end do
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do j=1,n_points(2)
       do i=1,n_points(1)
          do k=1,n_cells(3)
             do id=1,deg
                rho_cell(id) = sum(eval_dg_to_equi(:,id,3) * rho_3d(i,j,(k-1)*deg+1:k*deg))
             end do
             rho_3d(i,j,(k-1)*deg+1:k*deg) = rho_cell
          end do
       end do
    end do
!$omp end do
    deallocate(rho_cell)
!$omp end parallel
  end subroutine sll_s_compute_charge_density_6d_dd_slim_dg_staggered


  subroutine sll_s_exchange_charge_density_boundary_value_dd_slim_dg(topo, decomp, buf_tx, buf_rx, id)
    integer, parameter :: nd = 3
    type(sll_t_cartesian_topology_3d), intent(in) :: topo
    type(sll_t_decomposition_slim_3d), target, intent(inout) :: decomp
    sll_real64, intent(in) :: buf_tx(:,:)
    sll_real64, intent(out) :: buf_rx(:,:)
    sll_int32, intent(in) :: id
    integer :: ierr
    ! --- MPI communication-related variables and buffers
    sll_int32 :: nel
    integer, parameter :: mpi_precision = MPI_DOUBLE_PRECISION
    integer :: mpi_tag
    ! ------
    mpi_tag = 0
    nel = size(buf_tx)
    SLL_ASSERT_ALWAYS(nel == size(buf_rx))
    SLL_ASSERT_ALWAYS((id >= 0) .and. (id <= 3))
    ! --- copy buffer to the right neighbor
    if (nel > 0) then
      if (topo%procs(id) == 1) then
        buf_rx = buf_tx
      else
        call MPI_Sendrecv(buf_tx, nel, mpi_precision, topo%neighbors(2*id),   mpi_tag,&
                          buf_rx, nel, mpi_precision, topo%neighbors(2*id-1), mpi_tag,&
                          topo%comm, MPI_STATUS_IGNORE, ierr)
        SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
      endif
    endif
  end subroutine sll_s_exchange_charge_density_boundary_value_dd_slim_dg


  !> Compute charge density from distribution function (locally without data exchange)
  subroutine sll_s_compute_charge_density_6d( local_sizes, volume_v, f6d, rho)
    sll_int32,  intent(in)  :: local_sizes(6)
    sll_real64, intent(in)  :: volume_v
    sll_real64, intent(in)  :: f6d(1:,1:,1:,1:,1:,1:)
    sll_real64, intent(out) :: rho(1:, 1:, 1:)
!    sll_real64 :: sm
!    sll_int32 :: i,j,k,l,m,n,l_mx,m_mx,n_mx
    sll_int32 :: i, idx_mn(6), idx_mx(6)
! !$omp parallel default(shared) private(i,j,k,l,m,n)
! !$omp do OMP_COLLAPSE OMP_SCHEDULE
! !   do k=1,local_sizes(3)
! !      do j=1,local_sizes(2)
! !         do i=1,local_sizes(1)
! !            rho(i,j,k) = volume_v * sum(f6d(i,j,k,:,:,:))
! !         end do
! !      end do
! !   end do
! !$omp end do
! !$omp end parallel
    ! --- improved implementation/parallelization below ---
!    l_mx = size(f6d, 4)
!    m_mx = size(f6d, 5)
!    n_mx = size(f6d, 6)
!  !$omp parallel default(shared) private(i,j,k,l,m,n,sm)
!  !$omp do OMP_COLLAPSE OMP_SCHEDULE
!      do k=1,local_sizes(3)
!         do j=1,local_sizes(2)
!            do i=1,local_sizes(1)
!              sm = 0.0_f64
!              do n = 1,n_mx
!                do m = 1,m_mx
!                  do l = 1,l_mx
!                    sm = sm + f6d(i,j,k,l,m,n)
!                  enddo
!                enddo
!              enddo
!              rho(i,j,k) = sm * volume_v
!            end do
!         end do
!      end do
!  !$omp end do
!  !$omp end parallel

    idx_mn(:) = 1
    do i=1,3
      idx_mx(i) = local_sizes(i)
    enddo
    do i=4,6
      idx_mx(i) = size(f6d, i)
    enddo
    call sll_s_compute_charge_density_6d_core(f6d, idx_mn, idx_mx, &
                                              rho, idx_mn, idx_mx, &
                                              idx_mn, idx_mx, volume_v)
  end subroutine sll_s_compute_charge_density_6d


  ! Compute charge density from a domain-decomposed distribution function.
  subroutine sll_s_compute_charge_density_6d_dd(f6d, decomp_6d, rho_3d, coll_3d_v, volume_v)
    type(sll_t_decomposition_6d), target, intent(in) :: decomp_6d
    sll_real64, intent(in)  :: f6d(decomp_6d%local%lo(1):decomp_6d%local%hi(1), &
                                   decomp_6d%local%lo(2):decomp_6d%local%hi(2), &
                                   decomp_6d%local%lo(3):decomp_6d%local%hi(3), &
                                   decomp_6d%local%lo(4):decomp_6d%local%hi(4), &
                                   decomp_6d%local%lo(5):decomp_6d%local%hi(5), &
                                   decomp_6d%local%lo(6):decomp_6d%local%hi(6))
    type(sll_t_collective_t), pointer, intent(in) :: coll_3d_v
    sll_real64, intent(in)  :: volume_v
    sll_real64, intent(out) :: rho_3d(decomp_6d%local%mn(1):decomp_6d%local%mx(1), &
                                      decomp_6d%local%mn(2):decomp_6d%local%mx(2), &
                                      decomp_6d%local%mn(3):decomp_6d%local%mx(3))
!    sll_real64 :: sm
!    sll_int32 :: i, j, k, l, m, n
!    sll_int32, pointer :: loop_mn(:), loop_mx(:)
! !$omp parallel default(shared) private(i,j,k)
! !$omp do OMP_COLLAPSE OMP_SCHEDULE
!     do k=decomp_6d%local%mn(3), decomp_6d%local%mx(3)
!        do j=decomp_6d%local%mn(2), decomp_6d%local%mx(2)
!           do i=decomp_6d%local%mn(1), decomp_6d%local%mx(1)
!              rho_3d(i,j,k) = volume_v * sum(f6d(i,j,k,&
!                                                 decomp_6d%local%mn(4):decomp_6d%local%mx(4),&
!                                                 decomp_6d%local%mn(5):decomp_6d%local%mx(5),&
!                                                 decomp_6d%local%mn(6):decomp_6d%local%mx(6)))
!           end do
!        end do
!     end do
! !$omp end do
! !$omp end parallel

!      ! --- improved implementation/parallelization below ---
!      loop_mn => decomp_6d%local%mn
!      loop_mx => decomp_6d%local%mx
!  !$omp parallel default(shared) private(i,j,k,l,m,n,sm)
!  !$omp do OMP_COLLAPSE OMP_SCHEDULE
!      do k=loop_mn(3),loop_mx(3)
!         do j=loop_mn(2),loop_mx(2)
!            do i=loop_mn(1),loop_mx(1)
!              sm = 0.0_f64
!              do n=loop_mn(6),loop_mx(6)
!                do m=loop_mn(5),loop_mx(5)
!                  do l=loop_mn(4),loop_mx(4)
!                    sm = sm + f6d(i,j,k,l,m,n)
!                  enddo
!                enddo
!              enddo
!              rho_3d(i,j,k) = sm * volume_v
!            end do
!         end do
!      end do
!  !$omp end do
!  !$omp end parallel

    call sll_s_compute_charge_density_6d_core(f6d, decomp_6d%local%lo, decomp_6d%local%hi, &
                                              rho_3d, decomp_6d%local%mn, decomp_6d%local%mx, &
                                              decomp_6d%local%mn, decomp_6d%local%mx, volume_v)

    ! --- use an allreduce operation until we know where the result needs to go
    call sll_s_collective_allreduce_sum_3d_real64(coll_3d_v, rho_3d)
  end subroutine sll_s_compute_charge_density_6d_dd


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


  subroutine sll_s_compute_charge_density_6d_dd_slim_overlap(f6d, decomp_6d, rho_3d, coll_3d_v, volume_v)
    type(sll_t_decomposition), target, intent(in) :: decomp_6d
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
!    sll_real64 :: sm
!    sll_int32 :: i, j, k, l, m, n
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    loop_mn => decomp_6d%local%mn
    loop_mx => decomp_6d%local%mx

!  !$omp parallel default(shared) private(i,j,k,l,m,n,sm)
!  !$omp do OMP_COLLAPSE OMP_SCHEDULE
!      do k=loop_mn(3),loop_mx(3)
!         do j=loop_mn(2),loop_mx(2)
!            do i=loop_mn(1),loop_mx(1)
!              sm = 0.0_f64
!              do n=loop_mn(6),loop_mx(6)
!                do m=loop_mn(5),loop_mx(5)
!                  do l=loop_mn(4),loop_mx(4)
!                    sm = sm + f6d(i,j,k,l,m,n)
!                  enddo
!                enddo
!              enddo
!              rho_3d(i,j,k) = sm * volume_v
!            end do
!         end do
!      end do
!  !$omp end do
!  !$omp end parallel
    call sll_s_compute_charge_density_6d_core(f6d, loop_mn, loop_mx, &
                                              rho_3d, loop_mn, loop_mx, &
                                              loop_mn, loop_mx, volume_v)

    ! --- use an allreduce operation until we know where the result needs to go
    call sll_s_collective_allreduce_sum_3d_real64(coll_3d_v, rho_3d)
  end subroutine sll_s_compute_charge_density_6d_dd_slim_overlap


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


    
    subroutine sll_s_set_f0_v ( local_sizes, etas, f0 )
      sll_int32, intent( in ) :: local_sizes(6)
      type(sll_t_array), intent( in ) :: etas(6)
      sll_real64, intent( inout ) :: f0(1:,1:,1:)

      sll_int32 :: l,m,n
      sll_real64 :: v1, v2, v3, factor

      factor = 1.0_f64/(sll_p_twopi)**1.5_f64

!$omp parallel default(shared) private(l,m,n, v1, v2, v3)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do n=1, local_sizes(6)
         do m=1, local_sizes(5)
            v3 = etas(6)%vals(n)**2
            v2 = etas(5)%vals(m)**2
            do l=1, local_sizes(4)
               v1 = etas(4)%vals(l)**2
               f0(l,m,n) = &
                    exp( - 0.5_f64 * (v1 + v2 + v3 ) ) * factor 
          end do
       end do
    end do
!$omp end do
!$omp end parallel
      
      
    end subroutine sll_s_set_f0_v
  
  ! Charge density computation core routine used by the implementation-specific callers above,
  ! the index arrays are necessary to map remap, dd, dd_slim layouts to the same code.
  ! Routine uses cache blocking which is greatly improving the performance.
  subroutine sll_s_add_source_term ( f6d, sourcex, sourcev, factor, loop_mn, loop_mx )
    sll_int32, intent(in)   :: loop_mn(6), loop_mx(6)
    sll_real64, intent(inout)  :: f6d(loop_mn(1):loop_mx(1), &
         loop_mn(2):loop_mx(2), &
         loop_mn(3):loop_mx(3), &
         loop_mn(4):loop_mx(4), &
         loop_mn(5):loop_mx(5), &
         loop_mn(6):loop_mx(6))
    sll_real64, intent(in) :: sourcex(loop_mn(1):loop_mx(1), &
         loop_mn(2):loop_mx(2), &
         loop_mn(3):loop_mx(3))
    sll_real64, intent(in) :: sourcev(loop_mn(4):loop_mx(4), &
         loop_mn(5):loop_mx(5), &
         loop_mn(6):loop_mx(6))
    sll_real64, intent(in) :: factor
    
    sll_int32 :: i, j, k, l, m, n
    
!$omp parallel default(shared) private(i,j,k,l,m,n)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6),loop_mx(6)
       do m=loop_mn(5),loop_mx(5)
          do l=loop_mn(4),loop_mx(4)
             do k=loop_mn(3),loop_mx(3)
                do j=loop_mn(2),loop_mx(2)
                   do i=loop_mn(1),loop_mx(1)
                      f6d(i,j,k,l,m,n) = f6d(i,j,k,l,m,n) + &
                           factor * sourcex(i,j,k) * sourcev(l,m,n)
                   end do
                enddo
             enddo
          enddo
       end do
    end do
!$omp end do
!$omp end parallel
    
  end subroutine sll_s_add_source_term
  
  
  subroutine sll_s_collision_crook( f6d,etas, loop_mn, loop_mx, tau_c , f0, collision_mask)
    sll_int32, intent(in)   :: loop_mn(6), loop_mx(6)
    sll_real64, intent(inout)  :: f6d(loop_mn(1):loop_mx(1), &
         loop_mn(2):loop_mx(2), &
         loop_mn(3):loop_mx(3), &
         loop_mn(4):loop_mx(4), &
         loop_mn(5):loop_mx(5), &
         loop_mn(6):loop_mx(6))
         
    sll_real64, intent(in) :: f0(loop_mn(4):loop_mx(4), &
         loop_mn(5):loop_mx(5), &
         loop_mn(6):loop_mx(6))
      sll_real64, intent(in) :: collision_mask(loop_mn(4):loop_mx(4), &
         loop_mn(5):loop_mx(5), &
         loop_mn(6):loop_mx(6))
    type(sll_t_array), intent(in) :: etas(6)
    sll_real64, intent(in) :: tau_c

    sll_real64 :: r, factor, v1,v2, v3
    sll_int32 :: i, j, k, l, m, n
    
    factor = 1.0_f64/(sll_p_twopi)**1.5_f64

    
!$omp parallel default(shared) private(i,j,k,l,m,n,v1,v2,v3)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6),loop_mx(6)
       do m=loop_mn(5),loop_mx(5)
          do l=loop_mn(4),loop_mx(4)
             do k=loop_mn(3),loop_mx(3)
                do j=loop_mn(2),loop_mx(2)
                   do i=loop_mn(1),loop_mx(1)
                      f6d(i,j,k,l,m,n) = f6d(i,j,k,l,m,n) + tau_c* collision_mask(l,m,n) * (f6d(i,j,k,l,m,n) -f0(l,m,n))
                   end do
                enddo
             enddo
          enddo
       end do
    end do
!$omp end do
!$omp end parallel
    
  end subroutine sll_s_collision_crook
    
subroutine advect_background(local_sizes, etas, ev1, ev2, ez,  f_inout)
    sll_int32,                                intent( in ) :: local_sizes(6)
    type(sll_t_array),                        intent( in ) :: etas(6)
    sll_real64,                               intent(in)    :: ev1(:,:,:) !< displacement vector
    sll_real64,                               intent(in)    :: ev2(:,:,:) !< displacement vector
    sll_real64,                               intent(in)    :: ez(:,:,:) !< displacement vector
    sll_real64,                               intent(inout) :: f_inout(:,:,:,:,:,:) !< value of the function on input and advected function values on output
    sll_real64:: v1,v2,v3
    sll_int32 :: i,j,k,l,n, m
    sll_real64 :: factor
    

    factor = 1.0_f64/(sll_p_twopi)**1.5_f64
!$omp parallel default(shared) private(i,j,k,l,m,n,v1,v2,v3)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=1, local_sizes(6)
       do m= 1,local_sizes(5)
          do l= 1,local_sizes(4)
          v1 = etas(4)%vals(l)
          v2 = etas(5)%vals(m)
          v3 = etas(6)%vals(n)
!           print*, 'Hallo'
             do k= 1,local_sizes(3)
                do j =1, local_sizes(2)
                   do i= 1,local_sizes(1)                     
                        f_inout(i,j,k,l,m,n) = f_inout(i,j,k,l,m,n) - & 
                                                (v1*ev1(i,j,k) +&
                                                 v2*ev2(i,j,k) +&
                                                 v3*ez(i,j,k)) *&
                                                exp(-0.5_f64*(v1**2+v2**2+v3**2))*factor !-&
!                                                (-ev1(i,j,k)**2 &
!                                                 - ev2(i,j,k)**2 &
!                                                 - ez(i,j,k)**2 &
!                                                 + ev1(i,j,k)**2*v1**2 &
!                                                 + ev2(i,j,k)**2*v2**2 &
!                                                 + ez(i,j,k)**2*v3**2&
!                                                 + 2*ev1(i,j,k)*ev2(i,j,k)*v1*v2 &
!                                                 + 2*ev1(i,j,k)*ez(i,j,k)*v1*v3 &
!                                                 + 2*ev2(i,j,k)*ez(i,j,k)*v2*v3) &
!                                                 * exp(-0.5_f64*(v1**2+v2**2+v3**2))*facto
!                         f_inout(i,j,k,l,m,n)  = f_inout(i,j,k,l,m,n) + &
!                                                 factor *(exp(-0.5_f64*((v1+ev1(i,j,k))**2&
!                                                                     +(v2+ev2(i,j,k))**2&
!                                                                     +(v3+ez(i,j,k))**2)) - &
!                                                          exp(-0.5_f64*(v1**2+v2**2+v3**2)))
                   end do
                end do
             end do
          end do
       end do
    end do
!$omp end do
!$omp end parallel
  end subroutine advect_background
  
  subroutine sll_s_add_noise ( f6d,etas, loop_mn, loop_mx )
    sll_int32, intent(in)   :: loop_mn(6), loop_mx(6)
    sll_real64, intent(inout)  :: f6d(loop_mn(1):loop_mx(1), &
         loop_mn(2):loop_mx(2), &
         loop_mn(3):loop_mx(3), &
         loop_mn(4):loop_mx(4), &
         loop_mn(5):loop_mx(5), &
         loop_mn(6):loop_mx(6))
    type(sll_t_array), intent(in) :: etas(6)

    sll_real64 :: r, factor, v1,v2, v3
    sll_int32 :: i, j, k, l, m, n
    
    factor = 1.0_f64/(sll_p_twopi)**1.5_f64


!$omp parallel default(shared) private(i,j,k,l,m,n,v1,v2,v3)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6),loop_mx(6)
       do m=loop_mn(5),loop_mx(5)
          do l=loop_mn(4),loop_mx(4)
             do k=loop_mn(3),loop_mx(3)
                do j=loop_mn(2),loop_mx(2)
                   do i=loop_mn(1),loop_mx(1)
                        v3 = (etas(6)%vals(n-loop_mn(6)+1))**2
                        v2 = (etas(5)%vals(m-loop_mn(5)+1))**2
                        v1 = (etas(4)%vals(l-loop_mn(4)+1))**2
                      CALL RANDOM_NUMBER(r)
                      f6d(i,j,k,l,m,n) = f6d(i,j,k,l,m,n) + &
                           1e-7 *(r-0.5_f64) * factor * exp(-0.5*(v1+v2+v3))
                   end do
                enddo
             enddo
          enddo
       end do
    end do
!$omp end do
!$omp end parallel
    
  end subroutine sll_s_add_noise

  

  subroutine sll_s_compute_charge_density_6d_transp( local_sizes, volume_v, f6d, rho)
    sll_int32,  intent(in)  :: local_sizes(6)
    sll_real64, intent(in)  :: volume_v
    sll_real64, intent(in)  :: f6d(1:,1:,1:,1:,1:,1:)
    sll_real64, intent(out) :: rho(1:, 1:, 1:)

    sll_int32 :: i, j, k, l, m, n
    sll_int32 :: loop_mn(6), loop_mx(6)

    loop_mn = 1
    do i=4,6
      loop_mx(i) = local_sizes(i)
    enddo
    do i=1,3
      loop_mx(i) = size(f6d, i)
    enddo


!$omp parallel default(shared) private(i,j,k,l,m,n)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6),loop_mx(6)
       do m=loop_mn(5),loop_mx(5)
          do l=loop_mn(4),loop_mx(4)
             rho(l,m,n) = 0.0_f64
             do k=loop_mn(3),loop_mx(3)
                do j=loop_mn(2),loop_mx(2)
                   do i=loop_mn(1),loop_mx(1)
                      rho(l,m,n) = rho(l,m,n) + f6d(i,j,k,l,m,n)
                   end do
                end do
             enddo
             rho(l,m,n) = rho(l,m,n) * volume_v
          enddo
       enddo
    end do
!$omp end do
!$omp end parallel
  end subroutine sll_s_compute_charge_density_6d_transp

  

  !> Compute charge density from distribution function (locally without data exchange)
  subroutine sll_s_compute_momentum_energy_6d( local_sizes, volume_v, etas, vsq, f6d, f0v, momentum, energy)
    sll_int32,  intent(in)  :: local_sizes(6)
    sll_real64, intent(in)  :: volume_v
    type(sll_t_array), intent(in) :: etas(6)
    sll_real64, intent(in) :: vsq(1:,1:,1:)
    sll_real64, intent(in)  :: f6d(1:,1:,1:,1:,1:,1:)
    sll_real64, intent(in)  :: f0v(1:,1:,1:)
    sll_real64, intent(out) :: momentum(1:, 1:, 1:,1:)
    sll_real64, intent(out) :: energy(1:, 1:, 1:)

    sll_int32 :: i,j,k, l

    do k=1,local_sizes(3)
       do j=1,local_sizes(2)
          do i=1,local_sizes(1)
             do l=1,local_sizes(4)
                momentum(1,i,j,k) = momentum(1,i,j,k) + sum(etas(4)%vals(l)*(f6d(i,j,k,l,:,:)+f0v(l,:,:)))
             end do
             do l=1,local_sizes(5)
                momentum(2,i,j,k) = momentum(2,i,j,k) + sum(etas(5)%vals(l)*(f6d(i,j,k,:,l,:)+f0v(:,l,:)))
             end do
             do l=1,local_sizes(6)
                momentum(3,i,j,k) = momentum(3,i,j,k) + sum(etas(6)%vals(l)*(f6d(i,j,k,:,:,l)+f0v(:,:,l)))
             end do
             energy(i,j,k)     = sum(vsq*(f6d(i,j,k,:,:,:)+f0v))
          end do
       end do
    end do
    momentum = momentum * volume_v
    energy = energy * volume_v
  end subroutine sll_s_compute_momentum_energy_6d

  
    
    !> Compute charge density from distribution function (locally without data exchange)
  subroutine sll_s_compute_momentum_energy_6d_dd(f6d, decomp_6d, coll_3d_v, volume_v, etas, momentum_3d, energy_3d)
    type(sll_t_decomposition_slim_6d), target, intent(in) :: decomp_6d
    sll_real64, intent(in)  :: f6d(decomp_6d%local%mn(1):decomp_6d%local%mx(1), &
                                   decomp_6d%local%mn(2):decomp_6d%local%mx(2), &
                                   decomp_6d%local%mn(3):decomp_6d%local%mx(3), &
                                   decomp_6d%local%mn(4):decomp_6d%local%mx(4), &
                                   decomp_6d%local%mn(5):decomp_6d%local%mx(5), &
                                   decomp_6d%local%mn(6):decomp_6d%local%mx(6))
    type(sll_t_collective_t), pointer, intent(in) :: coll_3d_v
    sll_real64, intent(in)  :: volume_v
    type(sll_t_array), intent(in) :: etas(6)

    sll_real64, intent(out) :: momentum_3d(1:3,decomp_6d%local%mn(1):decomp_6d%local%mx(1), &
                                      decomp_6d%local%mn(2):decomp_6d%local%mx(2), &
                                      decomp_6d%local%mn(3):decomp_6d%local%mx(3))

    sll_real64, intent(out) :: energy_3d(decomp_6d%local%mn(1):decomp_6d%local%mx(1), &
                                      decomp_6d%local%mn(2):decomp_6d%local%mx(2), &
                                      decomp_6d%local%mn(3):decomp_6d%local%mx(3))
                                      
                                      
    sll_int32, pointer :: loop_mn(:), loop_mx(:)


    sll_real64, allocatable :: sme(:)
    sll_real64, allocatable :: smm(:,:)
    sll_int32 :: i, j, k, l, m, n
    
    
    
    loop_mn => decomp_6d%local%mn
    loop_mx => decomp_6d%local%mx
    !$omp parallel default(shared) private(i,j,k,l,m,n,sme, smm)
    allocate(sme(loop_mn(1):loop_mx(1))) 
    allocate(smm(3,loop_mn(1):loop_mx(1)))
    !$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k=loop_mn(3),loop_mx(3)
       do j=loop_mn(2),loop_mx(2)
          ! --- cache blocking in i ---
          sme(:) = 0.0_f64
          smm(:,:) = 0.0_f64
          do n=loop_mn(6),loop_mx(6)
            do m=loop_mn(5),loop_mx(5)
              do l=loop_mn(4),loop_mx(4)
                if ( etas(4)%vals(l-loop_mn(4)+1) .ne. -6.0 .and. etas(5)%vals(m-loop_mn(5)+1) .ne. -6.0 .and. etas(6)%vals(n-loop_mn(6)+1) .ne. -6.0 ) then

                    do i=loop_mn(1),loop_mx(1)
                    
                    sme(i) = sme(i) + (f6d(i,j,k,l,m,n)) * (etas(4)%vals(l-loop_mn(4)+1)**2 &
                                                        +etas(5)%vals(m-loop_mn(5)+1)**2 &
                                                        +etas(6)%vals(n-loop_mn(6)+1)**2)
                    smm(1,i) = smm(1,i) + (f6d(i,j,k,l,m,n)) * (etas(4)%vals(l-loop_mn(4)+1))
                    smm(2,i) = smm(2,i) + (f6d(i,j,k,l,m,n)) * (etas(5)%vals(m-loop_mn(5)+1))
                    smm(3,i) = smm(3,i) + (f6d(i,j,k,l,m,n)) * (etas(6)%vals(n-loop_mn(6)+1))
                    end do
                end if
              enddo
            enddo
          enddo
          do i=loop_mn(1),loop_mx(1)
            energy_3d(i,j,k) = sme(i) * volume_v
            momentum_3d(1,i,j,k) = smm(1,i) * volume_v
            momentum_3d(2,i,j,k) = smm(2,i) * volume_v
            momentum_3d(3,i,j,k) = smm(3,i) * volume_v
            
          enddo
       end do
    end do
    !$omp end do
    deallocate(sme)
    deallocate(smm)
!$omp end parallel

    call sll_s_collective_allreduce_sum_3d_real64(coll_3d_v, energy_3d)
    call sll_s_collective_allreduce_sum_3d_real64(coll_3d_v, momentum_3d(1,:,:,:))
    call sll_s_collective_allreduce_sum_3d_real64(coll_3d_v, momentum_3d(2,:,:,:))
    call sll_s_collective_allreduce_sum_3d_real64(coll_3d_v, momentum_3d(3,:,:,:))
        
  end subroutine sll_s_compute_momentum_energy_6d_dd
  
    !> Compute charge density from distribution function (locally without data exchange)
  subroutine sll_s_compute_heatflux_stress_6d_dd(f6d, decomp_6d, coll_3d_v, volume_v, etas, heatflux_3d, stress_3d)
    type(sll_t_decomposition_slim_6d), target, intent(in) :: decomp_6d
    sll_real64, intent(in)  :: f6d(decomp_6d%local%mn(1):decomp_6d%local%mx(1), &
                                   decomp_6d%local%mn(2):decomp_6d%local%mx(2), &
                                   decomp_6d%local%mn(3):decomp_6d%local%mx(3), &
                                   decomp_6d%local%mn(4):decomp_6d%local%mx(4), &
                                   decomp_6d%local%mn(5):decomp_6d%local%mx(5), &
                                   decomp_6d%local%mn(6):decomp_6d%local%mx(6))
    type(sll_t_collective_t), pointer, intent(in) :: coll_3d_v
    sll_real64, intent(in)  :: volume_v
    type(sll_t_array), intent(in) :: etas(6)

    sll_real64, intent(out) :: heatflux_3d(1:3,decomp_6d%local%mn(1):decomp_6d%local%mx(1), &
                                      decomp_6d%local%mn(2):decomp_6d%local%mx(2), &
                                      decomp_6d%local%mn(3):decomp_6d%local%mx(3))

    sll_real64, intent(out) :: stress_3d(1:3,1:3,&
                                      decomp_6d%local%mn(1):decomp_6d%local%mx(1), &
                                      decomp_6d%local%mn(2):decomp_6d%local%mx(2), &
                                      decomp_6d%local%mn(3):decomp_6d%local%mx(3))
                                      
                                      
    sll_int32, pointer :: loop_mn(:), loop_mx(:)

    sll_real64:: v2
    sll_real64, allocatable :: sme(:,:)
    sll_real64, allocatable :: smm(:,:,:)
    sll_int32 :: i, j, k, l, m, n, o, p
    
    
    
    loop_mn => decomp_6d%local%mn
    loop_mx => decomp_6d%local%mx
    !$omp parallel default(shared) private(i,j,k,l,m,n,o,p,sme, smm,v2)
    allocate(sme(1:3,loop_mn(1):loop_mx(1))) 
    allocate(smm(1:3,1:3,loop_mn(1):loop_mx(1)))
    !$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k=loop_mn(3),loop_mx(3)
       do j=loop_mn(2),loop_mx(2)
          ! --- cache blocking in i ---
          sme(:,:) = 0.0_f64
          smm(:,:,:) = 0.0_f64
          do n=loop_mn(6),loop_mx(6)
            do m=loop_mn(5),loop_mx(5)
              do l=loop_mn(4),loop_mx(4)
                if ( etas(4)%vals(l-loop_mn(4)+1) .ne. -4.0 .and. etas(5)%vals(m-loop_mn(5)+1) .ne. -4.0 .and. etas(6)%vals(n-loop_mn(6)+1) .ne. -4.0 ) then
        
                v2 = etas(4)%vals(l-loop_mn(4)+1)**2&
                    +etas(5)%vals(m-loop_mn(5)+1)**2&
                    +etas(6)%vals(n-loop_mn(6)+1)**2
                do i=loop_mn(1),loop_mx(1)
    
                    sme(1,i) = sme(1,i) + f6d(i,j,k,l,m,n) * etas(4)%vals(l-loop_mn(4)+1)*v2
                    sme(2,i) = sme(2,i) + f6d(i,j,k,l,m,n) * etas(5)%vals(m-loop_mn(5)+1)*v2
                    sme(3,i) = sme(3,i) + f6d(i,j,k,l,m,n) * etas(6)%vals(n-loop_mn(6)+1)*v2
                    
                    
                    smm(1,1,i) = smm(1,1,i) + f6d(i,j,k,l,m,n) * etas(4)%vals(l-loop_mn(4)+1)&
                                                            * etas(4)%vals(l-loop_mn(4)+1)
                    smm(1,2,i) = smm(1,2,i) + f6d(i,j,k,l,m,n) * etas(4)%vals(l-loop_mn(4)+1)&
                                                            * etas(5)%vals(m-loop_mn(5)+1)
                    smm(1,3,i) = smm(1,3,i) + f6d(i,j,k,l,m,n) * etas(4)%vals(l-loop_mn(4)+1)&
                                                            * etas(6)%vals(n-loop_mn(6)+1)

                    smm(2,1,i) = smm(2,1,i) + f6d(i,j,k,l,m,n) * etas(5)%vals(m-loop_mn(5)+1)&
                                                            * etas(4)%vals(l-loop_mn(4)+1)
                    smm(2,2,i) = smm(2,2,i) + f6d(i,j,k,l,m,n) * etas(5)%vals(m-loop_mn(5)+1)&
                                                            * etas(5)%vals(m-loop_mn(5)+1)
                    smm(2,3,i) = smm(2,3,i) + f6d(i,j,k,l,m,n) * etas(5)%vals(m-loop_mn(5)+1)&
                                                            * etas(6)%vals(n-loop_mn(6)+1)
                    
                    smm(3,1,i) = smm(3,1,i) + f6d(i,j,k,l,m,n) * etas(6)%vals(n-loop_mn(6)+1)&
                                                            * etas(4)%vals(l-loop_mn(4)+1)
                    smm(3,2,i) = smm(3,2,i) + f6d(i,j,k,l,m,n) * etas(6)%vals(n-loop_mn(6)+1)&
                                                            * etas(5)%vals(m-loop_mn(5)+1)
                    smm(3,3,i) = smm(3,3,i) + f6d(i,j,k,l,m,n) * etas(6)%vals(n-loop_mn(6)+1)&
                                                           * etas(6)%vals(n-loop_mn(6)+1)
                enddo
                end if
              enddo
            enddo
          enddo

          do i=loop_mn(1),loop_mx(1)
              do o = 1, 3
                  heatflux_3d(o,i,j,k) = sme(o,i) * volume_v
                  do p = 1,3
                      stress_3d(o,p,i,j,k) = smm(o,p,i) * volume_v
                  enddo
              enddo
          enddo
       end do
    end do
    !$omp end do
    deallocate(sme)
    deallocate(smm)
!$omp end parallel
            do o = 1, 3
              call sll_s_collective_allreduce_sum_3d_real64(coll_3d_v, heatflux_3d(o,:,:,:))
            do p = 1,3
                call sll_s_collective_allreduce_sum_3d_real64(coll_3d_v, stress_3d(o,p,:,:,:))
            enddo
          enddo

  end subroutine sll_s_compute_heatflux_stress_6d_dd
  

  subroutine sll_s_compute_free_energy_6d( local_sizes, indices_min, volume_x, vsq, f6d, fenergy, dfdivf0 )
    sll_int32,  intent(in)  :: local_sizes(6)
    sll_int32,  intent(in)  :: indices_min(6)
    sll_real64, intent(in)  :: volume_x
    sll_real64, intent(in)  :: vsq(1:,1:,1:)
    sll_real64, intent(in)  :: f6d(1:,1:,1:,1:,1:,1:)
    sll_real64, intent(out) :: fenergy(1:, 1:, 1: )
    sll_real64, intent(out) :: dfdivf0(1:, 1:, 1: )

    sll_int32 :: l,m,n
    sll_real64 :: factor
    sll_real64 :: f0

    factor = sqrt(sll_p_twopi)**3*volume_x

    do n=1,local_sizes(6)
       do m=1, local_sizes(5)
          do l=1, local_sizes(4)
             f0 = exp(0.5_f64*vsq(indices_min(4)+l-1,indices_min(5)+m-1,indices_min(6)+n-1))*factor
             dfdivf0(l,m,n) = sum(f6d(:,:,:,l,m,n))* f0
             fenergy(l,m,n) = sum(f6d(:,:,:,l,m,n)**2)* f0
          end do
       end do
    end do
  end subroutine sll_s_compute_free_energy_6d


  subroutine sll_s_compute_phi_qn (local_sizes_split, local_sizes_seq3, rmp_split2seq3, volume_v, C_tedn0, f6d, rho, phi)
    sll_int32,  intent(in)  :: local_sizes_split(6)
    sll_int32,  intent(in)  :: local_sizes_seq3(3)
    type(sll_t_remap_plan_3d_real64), pointer, intent(in) :: rmp_split2seq3
    sll_real64, intent(in)  :: volume_v
    sll_real64, intent(in)  :: c_tedn0
    sll_real64, intent(in)  :: f6d(1:,1:,1:,1:,1:,1:)
    sll_real64, intent(out) :: rho(1:, 1:, 1:)
    sll_real64, intent(out) :: phi(1:, 1:, 1:)

    sll_int32  :: i,j
    sll_real64 :: rhoz

    call sll_s_compute_charge_density_6d( local_sizes_split, volume_v, f6d, rho )
    call sll_o_apply_remap_3d ( rmp_split2seq3, rho, phi )
    phi = phi * c_tedn0
    do j=1,local_sizes_seq3(2)
       do i=1,local_sizes_seq3(1)
          rhoz = sum(phi(i,j,:))/real(local_sizes_seq3(3), f64)
          phi(i,j,:) = phi(i,j,:)-rhoz
       end do
    end do
  end subroutine sll_s_compute_phi_qn


  subroutine sll_s_time_history_diagnostics_dg(&
       local_sizes, &
       weights1, &
       weights2, &
       weights3, &
       weights4, &
       weights5, &
       weights6, &
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
    sll_int32,  intent( in ) :: local_sizes(6)
    sll_real64, intent( in ) :: weights1(:)
    sll_real64, intent( in ) :: weights2(:)
    sll_real64, intent( in ) :: weights3(:)
    sll_real64, intent( in ) :: weights4(:)
    sll_real64, intent( in ) :: weights5(:)
    sll_real64, intent( in ) :: weights6(:)
    sll_real64,        intent(in) :: time
    sll_real64,        intent(in) :: dV
    sll_real64,        intent(in) :: dV_x
    type(sll_t_array), intent(in) :: tensor_grid(6)
    sll_real64,        intent(in) :: f6d(:,:,:,:,:,:)
    sll_real64,        intent(in) :: rho(:,:,:)
    sll_real64,        intent(in) :: phi(:,:,:)
    sll_real64,        intent(in) :: ex(:,:,:)
    sll_real64,        intent(in) :: ey(:,:,:)
    sll_real64,        intent(in) :: ez(:,:,:)
    sll_int32,         intent(in) :: file_id
    sll_int32,optional,intent(in) :: topology_3d_velocity_coords(3)
    logical, optional, intent(out) :: nan_blowup

    sll_real64 :: diagn_data(13), sm, diagn_data_all(13)
    sll_int32 :: i, j, k
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

    ! calculate the maximum indices of the rho, ex, ey, ez arrays
    ! (plausibly assuming they are always the same)
    i_mx = size(rho, 1)
    j_mx = size(rho, 2)
    k_mx = size(rho, 3)
    ! --- phi has different extents
    i_phi_mx = size(phi, 1)
    j_phi_mx = size(phi, 2)
    k_phi_mx = size(phi, 3)

    diagn_data(:) = 0.0_f64
!$omp parallel default(shared) private(sm, i, j, k) reduction(+ : diagn_data)

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
!$omp end parallel


    ! mass
    call  sll_s_quadrature_6d( f6d, local_sizes, weights1, weights2, weights3, weights4, weights5, weights6, sm)
    diagn_data(1) = sm !* dV

    ! mass
    call  sll_s_l2_quadrature_6d( f6d, local_sizes, weights1, weights2, weights3, weights4, weights5, weights6, sm)
    diagn_data(2) = sm !* dV


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
  end subroutine sll_s_time_history_diagnostics_dg


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


  subroutine sll_s_time_history_diagnostics_transp(&
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
    sll_real64, allocatable :: sma(:)  ! array version of sm for cache blocking
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
!$omp parallel default(shared) private(sm,sma,smsq,smm,smp4,smp5,smp6,smk4,smk5,smk6,i,j,k,l,m,n) reduction(+ : diagn_data)

!    ! mass
!    diagn_data(1) = sum(f6d(index_min(1):index_max(1), index_min(2):index_max(2), &
!                            index_min(3):index_max(3), index_min(4):index_max(4), &
!                            index_min(5):index_max(5), index_min(6):index_max(6)))&
!                  * dV

    allocate(sma(index_min(1):index_max(1)))

    sm = 0.0_f64
    smsq = 0.0_f64
    smp4 = 0.0_f64
    smp5 = 0.0_f64
    smp6 = 0.0_f64
    smk4 = 0.0_f64
    smk5 = 0.0_f64
    smk6 = 0.0_f64
    smm = 0.0_f64


! --- non-cache-blocked code version, commented ---
! !$omp do OMP_COLLAPSE OMP_SCHEDULE
!     do k=index_min(3),index_max(3)
!        do j=index_min(2),index_max(2)
!           do i=index_min(1),index_max(1)
!             sm = 0.0_f64
!              do n=index_min(6),index_max(6)
!                 do m=index_min(5),index_max(5)
!                    do l=index_min(4),index_max(4)
!                      sm = sm + f6d(i,j,k,l,m,n)
!                      smsq = smsq + f6d(i,j,k,l,m,n)**2
!                   enddo
!                enddo
!             enddo
!             smp4 = smp4 + sm * tensor_grid(1)%vals(i-index_min(1)+1) * dV
!             smp5 = smp5 + sm * tensor_grid(2)%vals(j-index_min(2)+1) * dV
!             smp6 = smp6 + sm * tensor_grid(3)%vals(k-index_min(3)+1) * dV
!             smk4 = smk4 + sm * tensor_grid(1)%vals(i-index_min(1)+1)**2 * dV
!             smk5 = smk5 + sm * tensor_grid(2)%vals(j-index_min(2)+1)**2 * dV
!             smk6 = smk6 + sm * tensor_grid(3)%vals(k-index_min(3)+1)**2 * dV
!             smm = smm + sm
!         enddo
!       enddo
!     enddo
! !$omp end do nowait
!     diagn_data(1) = smm * dV
!     diagn_data(2) = smsq * dV
!     diagn_data(8) = smp4
!     diagn_data(9) = smp5
!     diagn_data(10) = smp6
!     diagn_data(11) = smk4
!     diagn_data(12) = smk5
!     diagn_data(13) = smk6


!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k=index_min(3),index_max(3)
      do j=index_min(2),index_max(2)
        sma(:) = 0.0_f64
        do n=index_min(6),index_max(6)
          do m=index_min(5),index_max(5)
            do l=index_min(4),index_max(4)
              ! cache blocking in i
              do i=index_min(1),index_max(1)
                sma(i) = sma(i) + f6d(i,j,k,l,m,n)
                smsq = smsq + f6d(i,j,k,l,m,n)**2
              enddo
            enddo
          enddo
        enddo

        do i=index_min(1),index_max(1)
          smp4 = smp4 + sma(i) * tensor_grid(1)%vals(i-index_min(1)+1)
          smp5 = smp5 + sma(i) * tensor_grid(2)%vals(j-index_min(2)+1)
          smp6 = smp6 + sma(i) * tensor_grid(3)%vals(k-index_min(3)+1)
          smk4 = smk4 + sma(i) * tensor_grid(1)%vals(i-index_min(1)+1)**2
          smk5 = smk5 + sma(i) * tensor_grid(2)%vals(j-index_min(2)+1)**2
          smk6 = smk6 + sma(i) * tensor_grid(3)%vals(k-index_min(3)+1)**2
          smm = smm + sma(i)
        enddo
      enddo
    enddo
!$omp end do nowait

    diagn_data(1) = smm * dV
    diagn_data(2) = smsq * dV
    diagn_data(8) = smp4 * dV
    diagn_data(9) = smp5 * dV
    diagn_data(10) = smp6 * dV
    diagn_data(11) = smk4 * dV
    diagn_data(12) = smk5 * dV
    diagn_data(13) = smk6 * dV

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

    deallocate(sma)
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
  end subroutine sll_s_time_history_diagnostics_transp



  subroutine sll_s_additional_time_history_diagnostics(&
       index_min, &
       index_max, &
       data_index_min, &
       time, &
       dV, &
       dV_x, &
       tensor_grid,&
       f6d, &
       f0v, &
       rho, &
       phi, &
       energy, &
       momentum, &
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
    sll_real64,        intent(in) :: f0v(data_index_min(4):,&
                                         data_index_min(5):,data_index_min(6):)
    sll_real64,        intent(in) :: rho(:,:,:)
    sll_real64,        intent(in) :: phi(:,:,:)
    sll_real64,        intent(in) :: energy(:,:,:)
    sll_real64,        intent(in) :: momentum(:,:,:,:)
    sll_real64,        intent(in) :: ex(:,:,:)
    sll_real64,        intent(in) :: ey(:,:,:)
    sll_real64,        intent(in) :: ez(:,:,:)
    sll_int32,         intent(in) :: file_id
    sll_int32,optional,intent(in) :: topology_3d_velocity_coords(3)
    logical, optional, intent(out) :: nan_blowup

    sll_real64 :: diagn_data(12), diagn_tmp1, diagn_data_all(12),diagn_tmp2,diagn_tmp3
    sll_int32 :: i, j, k, l, m, n
    sll_int32 :: i_mx, j_mx, k_mx
!    sll_int32 :: i_phi_mx, j_phi_mx, k_phi_mx
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

    if (do_electrostatics) then
!       !
!       diagn_data(1) = sum(rho*ex)*dV_x
!       diagn_data(2) = sum(rho*ey)*dV_x
!       diagn_data(3) = sum(rho*ez)*dV_x
       diagn_tmp1 = 0.0_f64
       diagn_tmp2 = 0.0_f64
       diagn_tmp3 = 0.0_f64
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,k_mx
        do j=1,j_mx
          do i=1,i_mx
            diagn_tmp1 = diagn_tmp1 + rho(i,j,k)*ex(i,j,k)
            diagn_tmp2 = diagn_tmp2 + rho(i,j,k)*ey(i,j,k)
            diagn_tmp3 = diagn_tmp3 + rho(i,j,k)*ez(i,j,k)
          enddo
        enddo
      enddo
!$omp end do nowait
      diagn_data(1) = diagn_tmp1 * dV_x
      diagn_data(2) = diagn_tmp2 * dV_x
      diagn_data(3) = diagn_tmp3 * dV_x

       diagn_tmp1 = 0.0_f64
       diagn_tmp2 = 0.0_f64
       diagn_tmp3 = 0.0_f64
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,k_mx
        do j=1,j_mx
          do i=1,i_mx
            diagn_tmp1 = diagn_tmp1 + energy(i,j,k)*ex(i,j,k)
            diagn_tmp2 = diagn_tmp2 + energy(i,j,k)*ey(i,j,k)
            diagn_tmp3 = diagn_tmp3 + energy(i,j,k)*ez(i,j,k)
          enddo
        enddo
      enddo
!$omp end do nowait
      diagn_data(4) = diagn_tmp1 * dV_x
      diagn_data(5) = diagn_tmp2 * dV_x
      diagn_data(6) = diagn_tmp3 * dV_x

   end if


!$omp do OMP_SCHEDULE
   do n=index_min(6),index_max(6)
      do m=index_min(5),index_max(5)
         do l=index_min(4),index_max(4)
            diagn_tmp1 = 0.0_f64
            do k=index_min(3),index_max(3)
               do j=index_min(2),index_max(2)
                  do i=index_min(1),index_max(1)
                     diagn_tmp1 = diagn_tmp1 + (f6d(i,j,k,l,m,n)+f0v(l,m,n))
                  enddo
               enddo
            enddo
            diagn_tmp1 = diagn_tmp1 * (tensor_grid(4)%vals(l-index_min(4)+1)**2+ tensor_grid(5)%vals(m-index_min(5)+1)**2+tensor_grid(6)%vals(n-index_min(6)+1)**2)
            diagn_data(7) = diagn_data(7) + diagn_tmp1 * &
                 tensor_grid(4)%vals(l-index_min(4)+1) * dV
            diagn_data(8) = diagn_data(8) + diagn_tmp1 * &
                 tensor_grid(5)%vals(m-index_min(5)+1) * dV
            diagn_data(9) = diagn_data(9) + diagn_tmp1 * &
                 tensor_grid(6)%vals(n-index_min(6)+1) * dV

         enddo
      enddo
   enddo
   !$omp end do nowait
   diagn_data(7:9) = diagn_data(7:9)!*dV


   diagn_tmp1 = 0.0_f64
   diagn_tmp2 = 0.0_f64
   diagn_tmp3 = 0.0_f64
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,k_mx
        do j=1,j_mx
          do i=1,i_mx
            diagn_tmp1 = diagn_tmp1 + momentum(1,i,j,k)*phi(i,j,k)
            diagn_tmp2 = diagn_tmp2 + momentum(2,i,j,k)*phi(i,j,k)
            diagn_tmp3 = diagn_tmp3 + momentum(3,i,j,k)*phi(i,j,k)
          enddo
        enddo
      enddo
!$omp end do nowait
      diagn_data(10) = diagn_tmp1 * dV_x
      diagn_data(11) = diagn_tmp2 * dV_x
      diagn_data(12) = diagn_tmp3 * dV_x
!!$!$omp do OMP_SCHEDULE
!!$   do n=index_min(6),index_max(6)
!!$      do m=index_min(5),index_max(5)
!!$         do l=index_min(4),index_max(4)
!!$            diagn_tmp1 = 0.0_f64
!!$            diagn_tmp2 = 0.0_f64
!!$            diagn_tmp3 = 0.0_f64
!!$            do k=index_min(3),index_max(3)
!!$               do j=index_min(2),index_max(2)
!!$                  do i=index_min(1),index_max(1)
!!$                     diagn_tmp1 = diagn_tmp1 + f6d(i,j,k,l,m,n) * ex(i,j,k)
!!$                     diagn_tmp2 = diagn_tmp2 + f6d(i,j,k,l,m,n) * ey(i,j,k)
!!$                     diagn_tmp3 = diagn_tmp3 + f6d(i,j,k,l,m,n) * ez(i,j,k)
!!$                  enddo
!!$               enddo
!!$            enddo
!!$            diagn_tmp1 = diagn_tmp1 * (tensor_grid(4)%vals(l-index_min(4)+1)) + &
!!$                  diagn_tmp2 * (tensor_grid(5)%vals(m-index_min(5)+1)) + &
!!$                  diagn_tmp3 * (tensor_grid(6)%vals(n-index_min(6)+1))
!!$
!!$            diagn_data(10) = diagn_data(10) + diagn_tmp1 * &
!!$                 tensor_grid(4)%vals(l-index_min(4)+1)
!!$            diagn_data(11) = diagn_data(11) + diagn_tmp1 * &
!!$                 tensor_grid(5)%vals(m-index_min(5)+1)
!!$            diagn_data(12) = diagn_data(12) + diagn_tmp1 * &
!!$                 tensor_grid(6)%vals(n-index_min(6)+1)
!!$
!!$         enddo
!!$      enddo
!!$   enddo
!!$   !$omp end do nowait
!!$   diagn_data(10:12) = diagn_data(10:12)*dV


   ! MPI collect the data
   diagn_data_all = 0.0_f64
    call sll_s_collective_reduce_real64(sll_v_world_collective, diagn_data, 12, MPI_SUM, 0, diagn_data_all)

    if(sll_f_get_collective_rank(sll_v_world_collective)==0) then
      write(file_id,'(e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12)') &!'(f12.5,2g20.12)') &
        time, diagn_data_all
    endif

 end subroutine sll_s_additional_time_history_diagnostics


 
  subroutine sll_s_additional_time_history_diagnostics_dd(&
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
       energy, &
       momentum, &
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
    sll_real64,        intent(in) :: energy(:,:,:)
    sll_real64,        intent(in) :: momentum(:,:,:,:)
    sll_real64,        intent(in) :: ex(:,:,:)
    sll_real64,        intent(in) :: ey(:,:,:)
    sll_real64,        intent(in) :: ez(:,:,:)
    sll_int32,         intent(in) :: file_id
    sll_int32,optional,intent(in) :: topology_3d_velocity_coords(3)
    logical, optional, intent(out) :: nan_blowup

    sll_real64 :: diagn_data(12), diagn_tmp1, diagn_data_all(12),diagn_tmp2,diagn_tmp3
    sll_int32 :: i, j, k, l, m, n
    sll_int32 :: i_mx, j_mx, k_mx
!    sll_int32 :: i_phi_mx, j_phi_mx, k_phi_mx
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

    if (do_electrostatics) then
!       !
!       diagn_data(1) = sum(rho*ex)*dV_x
!       diagn_data(2) = sum(rho*ey)*dV_x
!       diagn_data(3) = sum(rho*ez)*dV_x
       diagn_tmp1 = 0.0_f64
       diagn_tmp2 = 0.0_f64
       diagn_tmp3 = 0.0_f64
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,k_mx
        do j=1,j_mx
          do i=1,i_mx
            diagn_tmp1 = diagn_tmp1 + rho(i,j,k)*ex(i,j,k)
            diagn_tmp2 = diagn_tmp2 + rho(i,j,k)*ey(i,j,k)
            diagn_tmp3 = diagn_tmp3 + rho(i,j,k)*ez(i,j,k)
          enddo
        enddo
      enddo
!$omp end do nowait
      diagn_data(1) = diagn_tmp1 * dV_x
      diagn_data(2) = diagn_tmp2 * dV_x
      diagn_data(3) = diagn_tmp3 * dV_x

       diagn_tmp1 = 0.0_f64
       diagn_tmp2 = 0.0_f64
       diagn_tmp3 = 0.0_f64
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,k_mx
        do j=1,j_mx
          do i=1,i_mx
            diagn_tmp1 = diagn_tmp1 + energy(i,j,k)*ex(i,j,k)
            diagn_tmp2 = diagn_tmp2 + energy(i,j,k)*ey(i,j,k)
            diagn_tmp3 = diagn_tmp3 + energy(i,j,k)*ez(i,j,k)
          enddo
        enddo
      enddo
!$omp end do nowait
      diagn_data(4) = diagn_tmp1 * dV_x
      diagn_data(5) = diagn_tmp2 * dV_x
      diagn_data(6) = diagn_tmp3 * dV_x

   end if


!$omp do OMP_SCHEDULE
   do n=index_min(6),index_max(6)
      do m=index_min(5),index_max(5)
         do l=index_min(4),index_max(4)
            diagn_tmp1 = 0.0_f64
            do k=index_min(3),index_max(3)
               do j=index_min(2),index_max(2)
                  do i=index_min(1),index_max(1)
                     diagn_tmp1 = diagn_tmp1 + (f6d(i,j,k,l,m,n))
                  enddo
               enddo
            enddo
            if ( tensor_grid(4)%vals(l-index_min(4)+1) .ne. -6.0 .and. tensor_grid(5)%vals(m-index_min(5)+1) .ne. -6.0 .and. tensor_grid(6)%vals(n-index_min(6)+1) .ne. -6.0 ) then

                diagn_tmp1 = diagn_tmp1 * (tensor_grid(4)%vals(l-index_min(4)+1)**2+ tensor_grid(5)%vals(m-index_min(5)+1)**2+tensor_grid(6)%vals(n-index_min(6)+1)**2)
                diagn_data(7) = diagn_data(7) + diagn_tmp1 * &
                    tensor_grid(4)%vals(l-index_min(4)+1) * dV
                diagn_data(8) = diagn_data(8) + diagn_tmp1 * &
                    tensor_grid(5)%vals(m-index_min(5)+1) * dV
                diagn_data(9) = diagn_data(9) + diagn_tmp1 * &
                    tensor_grid(6)%vals(n-index_min(6)+1) * dV
                 
            endif 

         enddo
      enddo
   enddo
   !$omp end do nowait
   diagn_data(7:9) = diagn_data(7:9)!*dV


   diagn_tmp1 = 0.0_f64
   diagn_tmp2 = 0.0_f64
   diagn_tmp3 = 0.0_f64
   if (do_electrostatics ) then
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,k_mx
        do j=1,j_mx
          do i=1,i_mx
            diagn_tmp1 = diagn_tmp1 + momentum(1,i,j,k)*phi(i,j,k)
            diagn_tmp2 = diagn_tmp2 + momentum(2,i,j,k)*phi(i,j,k)
            diagn_tmp3 = diagn_tmp3 + momentum(3,i,j,k)*phi(i,j,k)
          enddo
        enddo
      enddo
      !$omp end do nowait
   end if
      diagn_data(10) = diagn_tmp1 * dV_x
      diagn_data(11) = diagn_tmp2 * dV_x
      diagn_data(12) = diagn_tmp3 * dV_x
!!$!$omp do OMP_SCHEDULE
!!$   do n=index_min(6),index_max(6)
!!$      do m=index_min(5),index_max(5)
!!$         do l=index_min(4),index_max(4)
!!$            diagn_tmp1 = 0.0_f64
!!$            diagn_tmp2 = 0.0_f64
!!$            diagn_tmp3 = 0.0_f64
!!$            do k=index_min(3),index_max(3)
!!$               do j=index_min(2),index_max(2)
!!$                  do i=index_min(1),index_max(1)
!!$                     diagn_tmp1 = diagn_tmp1 + f6d(i,j,k,l,m,n) * ex(i,j,k)
!!$                     diagn_tmp2 = diagn_tmp2 + f6d(i,j,k,l,m,n) * ey(i,j,k)
!!$                     diagn_tmp3 = diagn_tmp3 + f6d(i,j,k,l,m,n) * ez(i,j,k)
!!$                  enddo
!!$               enddo
!!$            enddo
!!$            diagn_tmp1 = diagn_tmp1 * (tensor_grid(4)%vals(l-index_min(4)+1)) + &
!!$                  diagn_tmp2 * (tensor_grid(5)%vals(m-index_min(5)+1)) + &
!!$                  diagn_tmp3 * (tensor_grid(6)%vals(n-index_min(6)+1))
!!$
!!$            diagn_data(10) = diagn_data(10) + diagn_tmp1 * &
!!$                 tensor_grid(4)%vals(l-index_min(4)+1)
!!$            diagn_data(11) = diagn_data(11) + diagn_tmp1 * &
!!$                 tensor_grid(5)%vals(m-index_min(5)+1)
!!$            diagn_data(12) = diagn_data(12) + diagn_tmp1 * &
!!$                 tensor_grid(6)%vals(n-index_min(6)+1)
!!$
!!$         enddo
!!$      enddo
!!$   enddo
!!$   !$omp end do nowait
!!$   diagn_data(10:12) = diagn_data(10:12)*dV


   ! MPI collect the data
   diagn_data_all = 0.0_f64
    call sll_s_collective_reduce_real64(sll_v_world_collective, diagn_data, 12, MPI_SUM, 0, diagn_data_all)

    if(sll_f_get_collective_rank(sll_v_world_collective)==0) then
      write(file_id,'(e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12,e20.12)') &!'(f12.5,2g20.12)') &
        time, diagn_data_all
    endif

 end subroutine sll_s_additional_time_history_diagnostics_dd


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


   subroutine sll_s_plot_diagnostic_time( data_rho, data_energy, data_moment, data_ex, data_ey, data_ez, offset, global_points, itime, filename)
    !use sll_m_collective
    !use hdf5
    !use sll_m_hdf5_io_parallel
    sll_real64, intent(in) :: data_rho(:,:,:)
    sll_real64, intent(in) :: data_ex(:,:,:)
    sll_real64, intent(in) :: data_ey(:,:,:)
    sll_real64, intent(in) :: data_ez(:,:,:)
    sll_real64, intent(in) :: data_energy(:,:,:)
    sll_real64, intent(in) :: data_moment(:,:,:,:)
    integer(i64),  intent(in) :: offset(3)
    integer(i64),  intent(in) :: global_points(3)
    sll_int32,  intent(in) :: itime
    character(*), intent(in) :: filename

    character(len=7)      :: ctime
    type(sll_t_hdf5_par_handle)     :: hdf_handle
    sll_int32  :: error

    call sll_s_int2string(itime, ctime)
    call sll_s_hdf5_par_file_create(filename//"-"//ctime//".h5", &
         sll_v_world_collective%comm, hdf_handle, error)


         
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, data_rho, &
         "rho", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, data_energy, &
         "energy", error )
    call sll_o_hdf5_par_write_array( hdf_handle, [3_8,global_points], [0_8,offset], data_moment, &
         "moments", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, data_ex, &
         "ex", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, data_ey, &
         "ey", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, data_ez, &
         "ez", error )
    call sll_s_hdf5_par_file_close(hdf_handle,error)

  end subroutine sll_s_plot_diagnostic_time

  
subroutine sll_s_plot_diagnostic_time_dd( data_rho, data_energy, data_moment, data_ex, data_ey, data_ez, offset, global_points, itime, filename,coll_3d_v)
    !use sll_m_collective
    !use hdf5
    !use sll_m_hdf5_io_parallel
    sll_real64, intent(in) :: data_rho(:,:,:)
    sll_real64, intent(in) :: data_ex(:,:,:)
    sll_real64, intent(in) :: data_ey(:,:,:)
    sll_real64, intent(in) :: data_ez(:,:,:)
    sll_real64, intent(in) :: data_energy(:,:,:)
    sll_real64, intent(in) :: data_moment(:,:,:,:)
    integer(i64),  intent(in) :: offset(3)
    integer(i64),  intent(in) :: global_points(3)
    sll_int32,  intent(in) :: itime
    character(*), intent(in) :: filename
    type(sll_t_collective_t), pointer, intent(in) :: coll_3d_v


    character(len=7)      :: ctime
    type(sll_t_hdf5_par_handle)     :: hdf_handle
    sll_int32  :: error

    call sll_s_int2string(itime, ctime)
    call sll_s_hdf5_par_file_create(filename//"-"//ctime//".h5", &
         coll_3d_v%comm, hdf_handle, error)



    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, data_rho, &
         "rho", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, data_energy, &
         "energy", error )
    call sll_o_hdf5_par_write_array( hdf_handle, [3_8,global_points], [0_8,offset], data_moment, &
         "moments", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, data_ex, &
         "ex", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, data_ey, &
         "ey", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, data_ez, &
         "ez", error )
    
    call sll_s_hdf5_par_file_close(hdf_handle,error)

  end subroutine sll_s_plot_diagnostic_time_dd
  
subroutine sll_s_plot_diagnostic_time_dd_add( data_heatflux, data_stress, offset, global_points, itime, filename,coll_3d_v)
    !use sll_m_collective
    !use hdf5
    !use sll_m_hdf5_io_parallel

    sll_real64, intent(in) :: data_heatflux(:,:,:,:)
    sll_real64, intent(in) :: data_stress(:,:,:,:,:)
    integer(i64),  intent(in) :: offset(3)
    integer(i64),  intent(in) :: global_points(3)
    sll_int32,  intent(in) :: itime
    character(*), intent(in) :: filename
    type(sll_t_collective_t), pointer, intent(in) :: coll_3d_v


    character(len=7)      :: ctime
    type(sll_t_hdf5_par_handle)     :: hdf_handle
    sll_int32  :: error

    call sll_s_int2string(itime, ctime)
    call sll_s_hdf5_par_file_create(filename//"_add"//"-"//ctime//".h5", &
         coll_3d_v%comm, hdf_handle, error)


    call sll_o_hdf5_par_write_array( hdf_handle, [3_8,3_8,global_points], [0_8,0_8,offset], data_stress, &
         "stress", error )
    call sll_o_hdf5_par_write_array( hdf_handle, [3_8,global_points], [0_8,offset], data_heatflux, &
         "heatflux", error )

    
    call sll_s_hdf5_par_file_close(hdf_handle,error)

  end subroutine sll_s_plot_diagnostic_time_dd_add


  subroutine sll_s_plot_v_diagnostic_time( fenergy, dfdivf0, f6d, offset, global_points, index_fmid, itime, filename)
    !use sll_m_collective
    !use hdf5
    !use sll_m_hdf5_io_parallel
    sll_real64, intent(in) :: fenergy(:,:,:)
    sll_real64, intent(in) :: dfdivf0(:,:,:)
    sll_real64, intent(in) :: f6d(:,:,:,:,:,:)
    integer(i64),  intent(in) :: offset(3)
    integer(i64),  intent(in) :: global_points(3)
    sll_int32, intent(in ) :: index_fmid(3)
    sll_int32,  intent(in) :: itime
    character(*), intent(in) :: filename

    character(len=5)      :: ctime
    type(sll_t_hdf5_par_handle)     :: hdf_handle
    sll_int32  :: error

    call sll_s_int2string(itime, ctime)
    call sll_s_hdf5_par_file_create(filename//"-"//ctime//".h5", &
         sll_v_world_collective%comm, hdf_handle, error)
         


    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, fenergy, &
         "fenergy", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, dfdivf0, &
         "dfdivf0", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, f6d(1,1,1,:,:,:), &
         "f111", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, &
         f6d(index_fmid(1), index_fmid(2), index_fmid(3),:,:,:), &
         "fmid", error )
         

    call sll_s_hdf5_par_file_close(hdf_handle,error)

  end subroutine sll_s_plot_v_diagnostic_time

  
  subroutine sll_s_plot_v_diagnostic_time_dd( f6d, offset, global_points, index_fmid, itime, filename)
    !use sll_m_collective
    !use hdf5
    !use sll_m_hdf5_io_parallel
    sll_real64, intent(in) :: f6d(:,:,:,:,:,:)
    integer(i64),  intent(in) :: offset(3)
    integer(i64),  intent(in) :: global_points(3)
    sll_int32, intent(in ) :: index_fmid(3)
    sll_int32,  intent(in) :: itime
    character(*), intent(in) :: filename

    character(len=7)      :: ctime
    type(sll_t_hdf5_par_handle)     :: hdf_handle
    sll_int32  :: error

    call sll_s_int2string(itime, ctime)
    call sll_s_hdf5_par_file_create(filename//"-"//ctime//".h5", &
         sll_v_world_collective%comm, hdf_handle, error)


    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, f6d(1,1,1,:,:,:), &
         "f111", error )
    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, &
         f6d(index_fmid(1), index_fmid(2), index_fmid(3),:,:,:), &
         "fmid", error )

    call sll_s_hdf5_par_file_close(hdf_handle,error)

  end subroutine sll_s_plot_v_diagnostic_time_dd
  
  subroutine sll_s_plot_distribution_time_6d( f6d, offset, global_points, itime, filename)
    !use sll_m_collective
    !use sll_m_hdf5_io_parallel
    sll_real64, intent(in) :: f6d(:,:,:,:,:,:)
    integer(i64),  intent(in) :: offset(6)
    integer(i64),  intent(in) :: global_points(6)
    sll_int32,  intent(in) :: itime
    character(*), intent(in) :: filename

    character(len=4)      :: ctime
    type(sll_t_hdf5_par_handle)     :: hdf_handle
    sll_int32  :: error

    call sll_s_int2string(itime, ctime)
    call sll_s_hdf5_par_file_create(filename//"-"//ctime//".h5", &
         sll_v_world_collective%comm, hdf_handle, error)

    call sll_o_hdf5_par_write_array( hdf_handle, global_points, offset, f6d, &
         "f", error )

    call sll_s_hdf5_par_file_close(hdf_handle,error)

  end subroutine sll_s_plot_distribution_time_6d


  subroutine sll_s_read_distribution_6d( f6d, offset, global_points, filename)
    sll_real64, intent(inout) :: f6d(:,:,:,:,:,:)
    integer(i64),  intent(in) :: offset(6)
    integer(i64),  intent(in) :: global_points(6)
    character(*), intent(in) :: filename

    type(sll_t_hdf5_par_handle)     :: hdf_handle
    sll_int32  :: error

    call sll_s_hdf5_par_file_open( trim(filename), sll_v_world_collective%comm,&
         hdf_handle, error )

    call sll_o_hdf5_par_read_array( hdf_handle, global_points, offset , f6d, &
         "f", error )
    SLL_ASSERT_ALWAYS( error == 0 )
    

    call sll_s_hdf5_par_file_close( hdf_handle, error )

  end subroutine sll_s_read_distribution_6d

  
  subroutine sll_s_read_distribution_6d_interpolate(topo_6d, f6d, dist_mod, offset, global_points, filename)
    sll_real64, intent(inout) :: f6d(:,:,:,:,:,:)
    sll_int32     :: hw(6)
    type(sll_t_cartesian_topology_6d), target, intent(inout), optional :: topo_6d
     type(sll_t_decomposition_slim_6d), pointer :: decomp_6d ! Decomposition object for 6d distribution function

    sll_real64,  intent(in) :: dist_mod(6)
    integer(i64),  intent(in) :: offset(6)
    integer(i64),  intent(in) :: global_points(6)
    character(*), intent(in) :: filename
    sll_real64, allocatable :: f6d1(:,:,:,:,:,:)
    sll_real64, allocatable :: f6d2(:,:,:,:,:,:)
    type(sll_t_hdf5_par_handle)     :: hdf_handle
    sll_int32  :: error, i, j 
    
    sll_int32 :: dim_in(6)
    sll_int32 :: dim_out(6)
    sll_int32 :: mod_log(6)
    sll_real64 :: factor(6)
    sll_real64 :: factor_base(6)
    

    factor = 1.0_f64
    
    dim_out = shape(f6d)
    
    dim_in = dim_out/dist_mod
    
    mod_log = int(log(dist_mod)/log(2.0_f64), i64)
    
    allocate(f6d1(dim_in(1),dim_in(2),dim_in(3),dim_in(4),dim_in(5),dim_in(6)))

    call sll_s_read_distribution_6d ( f6d1, &
                  int(offset/dist_mod, i64), &
                  int(global_points/dist_mod, i64), &
                  filename )
    
    do i = 1, 6 
        do j = 1, abs(mod_log(i))
            if (mod_log(i) .gt. 0) then 
                factor(i) = 2.0 *factor(i)
            else 
                factor(i) = 0.5 *factor(i)
            endif
            
            if ( allocated(f6d2)) deallocate(f6d2)
                        
            allocate(f6d2(  int(factor(1)*dim_in(1),i64),&
                            int(factor(2)*dim_in(2),i64),&
                            int(factor(3)*dim_in(3),i64),&
                            int(factor(4)*dim_in(4),i64),&
                            int(factor(5)*dim_in(5),i64),&
                            int(factor(6)*dim_in(6),i64)))
                            
            if( mod_log(i) .lt. 0) then
                call sll_s_half_dim_distribution_6d(f6d1, i, f6d2,shape(f6d2))
            else if (mod_log(i) .gt. 0) then
                decomp_6d => &
                    sll_f_new_cartesian_domain_decomposition_slim_6d(topo_6d, &
                    topo_6d%procs*shape(f6d1))
                call sll_s_double_dim_distribution_6d(topo_6d, decomp_6d, f6d1, i, f6d2,shape(f6d1))
            end if

            deallocate(f6d1)
            
            allocate(f6d1(  int(factor(1)*dim_in(1),i64),&
                            int(factor(2)*dim_in(2),i64),&
                            int(factor(3)*dim_in(3),i64),&
                            int(factor(4)*dim_in(4),i64),&
                            int(factor(5)*dim_in(5),i64),&
                            int(factor(6)*dim_in(6),i64)))
            f6d1 = f6d2    
        enddo       
    enddo
    
    f6d = f6d2
    if ( allocated(f6d1)) deallocate(f6d1)
    if ( allocated(f6d2)) deallocate(f6d2)
  end subroutine sll_s_read_distribution_6d_interpolate

  
  
  
  subroutine sll_s_double_dim_distribution_6d(topo_6d, decomp_6d, f6d, ndim,f6d_out,local_sizes)
    sll_real64, intent(inout) :: f6d(:,:,:,:,:,:)
    sll_int32, intent(in)     :: ndim 
    sll_int32     :: hw(6)
    type(sll_t_cartesian_topology_6d), target, intent(inout), optional :: topo_6d
    type(sll_t_decomposition_slim_6d), target, intent(inout), optional :: decomp_6d
    type(sll_t_advection_6d_lagrange_dd_slim) :: ladvector
    sll_int32,  intent( in ) :: local_sizes(6)
    sll_int32       ::      arraysize(6)
    
    sll_real64, allocatable, dimension(:)  :: delta
    sll_real64, allocatable, dimension(:,:,:)  :: delta2
    
    sll_real64, intent(out), dimension(:,:,:,:,:,:) :: f6d_out
    sll_int32   ::  factor(6) = (/1,1,1,1,1,1/)
    sll_int32   ::  shift(6) = (/0,0,0,0,0,0/)
    
     sll_int32 :: i,j,k,l,m,n
     
     factor = 1
     shift = 0

     
    call sll_s_advection_6d_lagrange_dd_slim_init( &
        ladvector, (/7,7/) )

    hw = (7-1)/2
     
    arraysize = shape(f6d)
    arraysize(ndim) = 2* arraysize(ndim)
    
    factor(ndim) = 2
    shift(ndim) = 1
    
    if (ndim < 4) then
        allocate(delta(arraysize(ndim+3)))
        delta =  0.5_f64
    else
        allocate(delta2(arraysize(1),arraysize(2),arraysize(3)))
        delta2 = 0.5_f64
    endif
    
    
!     allocate(f6d_out(arraysize(1),arraysize(2),arraysize(3),arraysize(4),arraysize(5),arraysize(6)))
!$omp parallel default(shared) private(i,j,k,l,m,n)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n = 1, local_sizes(6)
       do m = 1, local_sizes(5)
          do l = 1, local_sizes(4)
             do k = 1, local_sizes(3)
                do j = 1, local_sizes(2)
                   do i = 1, local_sizes(1)
                      f6d_out(factor(1)*i-shift(1),factor(2)*j-shift(2),factor(3)*k-shift(3),factor(4)*l-shift(4),factor(5)*m-shift(5),factor(6)*n-shift(6)) = f6d(i,j,k,l,m,n)
                   end do
                end do
             end do
          end do
       end do
    end do
!$omp end do
!$omp end parallel 
   
    if (ndim < 4) then
        call sll_advect_1d_x(topo_6d, decomp_6d, f6d, ndim, hw,ladvector, delta)
    else
        call sll_advect_1d_v(topo_6d, decomp_6d, f6d, ndim, hw,ladvector, delta2)
    endif    



!$omp parallel default(shared) private(i,j,k,l,m,n)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n = 1, local_sizes(6)
       do m = 1, local_sizes(5)
          do l = 1, local_sizes(4)
             do k = 1, local_sizes(3)
                do j = 1, local_sizes(2)
                   do i = 1, local_sizes(1)
                      f6d_out(factor(1)*i,factor(2)*j,factor(3)*k,factor(4)*l,factor(5)*m,factor(6)*n)= f6d(i,j,k,l,m,n)
                   end do
                end do
             end do
          end do
       end do
    end do
!$omp end do
!$omp end parallel 
    

          
  end subroutine sll_s_double_dim_distribution_6d
  
  subroutine sll_s_half_dim_distribution_6d(f6d, ndim, f6d_out,local_sizes)
    sll_real64, intent(inout) :: f6d(:,:,:,:,:,:)
    sll_int32, intent(in)     :: ndim 
    sll_int32,  intent( in ) :: local_sizes(6)
    sll_int32       ::      arraysize(6)
    
    sll_real64, intent(out), dimension(:,:,:,:,:,:) :: f6d_out
    sll_int32   ::  factor(6) = (/1,1,1,1,1,1/)
    sll_int32   ::  shift(6) = (/0,0,0,0,0,0/)
    
     sll_int32 :: i,j,k,l,m,n

    
    factor(ndim) = 2    
    shift(ndim) = 1
    

    

!     allocate(f6d_out(arraysize(1),arraysize(2),arraysize(3),arraysize(4),arraysize(5),arraysize(6)))
!$omp parallel default(shared) private(i,j,k,l,m,n)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n = 1, local_sizes(6)
       do m = 1, local_sizes(5)
          do l = 1, local_sizes(4)
             do k = 1, local_sizes(3)
                do j = 1, local_sizes(2)
                   do i = 1, local_sizes(1)
                      f6d_out(i,j,k,l,m,n) = f6d(factor(1)*i-shift(1),factor(2)*j-shift(2),factor(3)*k-shift(3),factor(4)*l-shift(4),factor(5)*m-shift(5),factor(6)*n-shift(6))
                   end do
                end do
             end do
          end do
       end do
    end do
!$omp end do
!$omp end parallel 
    

  end subroutine sll_s_half_dim_distribution_6d
  
  
  subroutine sll_advect_1d_x(topo_6d, decomp_6d, f6d, ndim, hw,ladvector, delta)
    sll_real64, intent(inout) :: f6d(:,:,:,:,:,:)
    sll_int32, intent(in)     :: ndim 
    sll_int32, intent(in)     :: hw(6)
    sll_real64, intent(in)     :: delta(:)
    type(sll_t_cartesian_topology_6d), target, intent(inout), optional :: topo_6d
    type(sll_t_decomposition_slim_6d), target, intent(inout), optional :: decomp_6d
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: ladvector
    call sll_f_apply_halo_exchange(topo_6d, &
                decomp_6d, &
                f6d, &
                ndim, &
                hw(ndim), &
                hw(ndim))
    select case (ndim)
      case (1) 
        call sll_s_advection_6d_lagrange_dd_slim_advect_eta1(ladvector, &
                    decomp_6d, &
                    delta, &
                    f6d)
      case (2)
        call sll_s_advection_6d_lagrange_dd_slim_advect_eta2(ladvector, &
                    decomp_6d, &
                    delta, &
                    f6d)
      case (3) 
        call sll_s_advection_6d_lagrange_dd_slim_advect_eta3(ladvector, &
                    decomp_6d, &
                    delta, &
                    f6d)
      case default
         print*, "Invalid Dimension" 
      
    end select

  end subroutine sll_advect_1d_x

  
  subroutine sll_advect_1d_v(topo_6d, decomp_6d, f6d, ndim, hw,ladvector, delta)
    sll_real64, intent(inout) :: f6d(:,:,:,:,:,:)
    sll_int32, intent(in)     :: ndim 
    sll_int32, intent(in)     :: hw(6)
    sll_real64, intent(in)     :: delta(:,:,:)
    type(sll_t_cartesian_topology_6d), target, intent(inout), optional :: topo_6d
    type(sll_t_decomposition_slim_6d), target, intent(inout), optional :: decomp_6d
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: ladvector
    
    call sll_f_apply_halo_exchange(topo_6d, &
                decomp_6d, &
                f6d, &
                ndim, &
                hw(ndim), &
                hw(ndim))
    select case (ndim)
      case (4) 
        call sll_s_advection_6d_lagrange_dd_slim_advect_eta4(ladvector, &
                    decomp_6d, &
                    delta, &
                    f6d)
      case (5)
        call sll_s_advection_6d_lagrange_dd_slim_advect_eta5(ladvector, &
                    decomp_6d, &
                    delta, &
                    f6d)
      case (6) 
        call sll_s_advection_6d_lagrange_dd_slim_advect_eta6(ladvector, &
                    decomp_6d, &
                    delta, &
                    f6d)
      case default
         print*, "Invalid Dimension" 
      
    end select

  end subroutine sll_advect_1d_v

  
  

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


  subroutine sll_s_quadrature_3d( f3d, local_sizes, weights1, weights2, weights3, val )
    sll_real64, intent( in ) :: f3d(:,:,:)
    sll_int32,  intent( in ) :: local_sizes(3)
    sll_real64, intent( in ) :: weights1(:)
    sll_real64, intent( in ) :: weights2(:)
    sll_real64, intent( in ) :: weights3(:)
    sll_real64, intent( out) :: val

    sll_int32 :: i,j,k
    sll_real64 :: sm

    sm = 0.0_f64
!$omp parallel default(shared) private(i,j,k) reduction(+ : sm)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k = 1, local_sizes(3)
       do j = 1, local_sizes(2)
          do i = 1, local_sizes(1)
             sm = sm + f3d(i,j,k) * weights1(i) * weights2(j) * weights3(k)
          end do
       end do
    end do
!$omp end do
!$omp end parallel

    val = sm
  end subroutine sll_s_quadrature_3d


  subroutine sll_s_l2_quadrature_3d( f3d, local_sizes, weights1, weights2, weights3, val )
    sll_real64, intent( in ) :: f3d(:,:,:)
    sll_int32,  intent( in ) :: local_sizes(3)
    sll_real64, intent( in ) :: weights1(:)
    sll_real64, intent( in ) :: weights2(:)
    sll_real64, intent( in ) :: weights3(:)
    sll_real64, intent( out) :: val

    sll_int32 :: i,j,k
    sll_real64 :: sm

    sm = 0.0_f64
!$omp parallel default(shared) private(i,j,k) reduction(+ : sm)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k = 1, local_sizes(3)
       do j = 1, local_sizes(2)
          do i = 1, local_sizes(1)
             sm = sm + f3d(i,j,k)**2 * weights1(i) * weights2(j) * weights3(k)
          end do
       end do
    end do
!$omp end do
!$omp end parallel

    val = sm
  end subroutine sll_s_l2_quadrature_3d



  subroutine sll_s_quadrature_6d( f6d, local_sizes, weights1, weights2, weights3, weights4, weights5, weights6, val )
    sll_real64, intent( in ) :: f6d(:,:,:,:,:,:)
    sll_int32,  intent( in ) :: local_sizes(6)
    sll_real64, intent( in ) :: weights1(:)
    sll_real64, intent( in ) :: weights2(:)
    sll_real64, intent( in ) :: weights3(:)
    sll_real64, intent( in ) :: weights4(:)
    sll_real64, intent( in ) :: weights5(:)
    sll_real64, intent( in ) :: weights6(:)
    sll_real64, intent( out) :: val

    sll_int32 :: i,j,k,l,m,n
    sll_real64 :: sm

    sm = 0.0_f64
!$omp parallel default(shared) private(i,j,k,l,m,n) reduction(+ : sm)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n = 1, local_sizes(6)
       do m = 1, local_sizes(5)
          do l = 1, local_sizes(4)
             do k = 1, local_sizes(3)
                do j = 1, local_sizes(2)
                   do i = 1, local_sizes(1)
                      sm = sm + f6d(i,j,k,l,m,n) * weights1(i) * weights2(j) * weights3(k)* &
                                                   weights4(l) * weights5(m) * weights6(n)
                   end do
                end do
             end do
          end do
       end do
    end do
!$omp end do
!$omp end parallel

    val = sm
  end subroutine sll_s_quadrature_6d



  subroutine sll_s_l2_quadrature_6d( f6d, local_sizes, weights1, weights2, weights3, weights4, weights5, weights6, val )
    sll_real64, intent( in ) :: f6d(:,:,:,:,:,:)
    sll_int32,  intent( in ) :: local_sizes(6)
    sll_real64, intent( in ) :: weights1(:)
    sll_real64, intent( in ) :: weights2(:)
    sll_real64, intent( in ) :: weights3(:)
    sll_real64, intent( in ) :: weights4(:)
    sll_real64, intent( in ) :: weights5(:)
    sll_real64, intent( in ) :: weights6(:)
    sll_real64, intent( out) :: val

    sll_int32 :: i,j,k,l,m,n
    sll_real64 :: sm

    sm = 0.0_f64
!$omp parallel default(shared) private(i,j,k,l,m,n) reduction(+ : sm)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n = 1, local_sizes(6)
       do m = 1, local_sizes(5)
          do l = 1, local_sizes(4)
             do k = 1, local_sizes(3)
                do j = 1, local_sizes(2)
                   do i = 1, local_sizes(1)
                      sm = sm + f6d(i,j,k,l,m,n)**2 * weights1(i) * weights2(j) * weights3(k)* &
                                                      weights4(l) * weights5(m) * weights6(n)
                   end do
                end do
             end do
          end do
       end do
    end do
!$omp end do
!$omp end parallel

    val = sm
  end subroutine sll_s_l2_quadrature_6d


  subroutine sll_s_quadrature_6dto3d( f6d, local_sizes, weights4, weights5, weights6, f3d )
    sll_real64, intent( in ) :: f6d(:,:,:,:,:,:)
    sll_int32,  intent( in ) :: local_sizes(6)
    sll_real64, intent( in ) :: weights4(:)
    sll_real64, intent( in ) :: weights5(:)
    sll_real64, intent( in ) :: weights6(:)
    sll_real64, intent( out) :: f3d(:,:,:)

    sll_int32 :: i,j,k,l,m,n
    sll_real64, allocatable :: sm(:)

!$omp parallel default(shared) private(i,j,k,l,m,n,sm)
    allocate(sm(1:local_sizes(1)))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k = 1, local_sizes(3)
       do j = 1, local_sizes(2)

          sm(:) = 0.0_f64
          do n = 1, local_sizes(6)
             do m = 1, local_sizes(5)
                do l = 1, local_sizes(4)
                   do i = 1, local_sizes(1)
                      sm(i) = sm(i) + f6d(i,j,k,l,m,n) *  weights4(l) * weights5(m) * weights6(n)
                   end do
                end do
             end do
          end do
          do i = 1, local_sizes(1)
             f3d(i,j,k) = sm(i)
          end do

       end do
    end do
!$omp end do
    deallocate(sm)
!$omp end parallel

  end subroutine sll_s_quadrature_6dto3d



  subroutine sll_s_gauss_quadrature_cells_6dto3d( f6d, local_size, local_cells, degree, weights, f3d )
    sll_real64, intent( in ) :: f6d(:,:,:,:,:,:)
    sll_int32,  intent( in ) :: local_size(6)
    sll_int32,  intent( in ) :: local_cells(6)
    sll_int32,  intent( in ) :: degree(6)
    sll_real64, intent( in ) :: weights(:,:,:)
    sll_real64, intent( out) :: f3d(:,:,:)

    sll_int32 :: i,j,k,l,m,n,ld,li,md,mi,nd,ni
    sll_real64, allocatable :: sm(:)

!$omp parallel default(shared) private(i,j,k,l,m,n,ld,li,md,mi,nd,ni,sm)
    allocate(sm(1:local_size(1)))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k = 1, local_size(3)
       do j = 1, local_size(2)

          sm(:) = 0.0_f64
          ni = 0
          do n = 1, local_cells(6)
             do nd = 1, degree(6)
                ni = ni+1
                mi = 0
                do m = 1, local_cells(5)
                   do md = 1, degree(5)
                      mi = mi+1
                      li = 0
                      do l = 1, local_cells(4)
                         do ld = 1, degree(4)
                            li = li+1
                            do i =1, local_size(1)
                               sm(i) = sm(i) + f6d(i,j,k,li, mi, ni ) * weights(ld, md, nd)
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
          do i= 1, local_size(1)
             f3d(i,j,k) = sm(i)
          end do

       end do
    end do
!$omp end do
    deallocate(sm)
!$omp end parallel

  end subroutine sll_s_gauss_quadrature_cells_6dto3d


  subroutine sll_s_gauss_quadrature_cells_6d( f6d, local_size, local_cells, degree, weights, val )
    sll_real64, intent( in ) :: f6d(:,:,:,:,:,:)
    sll_int32,  intent( in ) :: local_size(6)
    sll_int32,  intent( in ) :: local_cells(6)
    sll_int32,  intent( in ) :: degree(6)
    sll_real64, intent( in ) :: weights(:,:,:,:,:,:)
    sll_real64, intent( out) :: val

    sll_int32 :: i,j,k,l,m,n,ld,li,md,mi,nd,ni, ii,ji,ki,id,jd,kd
    sll_real64 :: sm

    sm = 0.0_f64
!$omp parallel default(shared) private(i,j,k,l,m,n,ii,ji,ki,li,mi,ni,id,jd,kd,ld,md,nd) reduction(+ : sm)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n = 1, local_cells(6)
       do m = 1, local_cells(5)

          do nd = 1, degree(6)
             ni = (n-1)*degree(6)+nd
             do md = 1, degree(5)
                mi = (m-1)*degree(5)+md
                do l = 1, local_cells(4)
                   do ld = 1, degree(4)
                      li = (l-1)*degree(4)+ld
                      do k =1, local_cells(3)
                         do kd = 1, degree(3)
                            ki = (k-1)*degree(3)+kd
                            do j = 1, local_cells(2)
                               do jd = 1, degree(2)
                                  ji = (j-1)*degree(2)+jd
                                  do i=1, local_cells(1)
                                     do id =1, degree(1)
                                        ii = (i-1)*degree(1)+id
                                        sm = sm + f6d(ii,ji,ki,li,mi,ni) * weights(id,jd,kd,ld,md,nd)
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do

       end do
    end do
!$omp end do
!$omp end parallel

    val = sm
  end subroutine sll_s_gauss_quadrature_cells_6d


  subroutine sll_s_gauss_quadrature_cells_3d( f3d, local_cells, degree, weights, val )
    sll_real64, intent( in ) :: f3d(:,:,:)
    sll_int32,  intent( in ) :: local_cells(3)
    sll_int32,  intent( in ) :: degree(3)
    sll_real64, intent( in ) :: weights(:,:,:)
    sll_real64, intent( out) :: val

    sll_int32 :: i,j,k,ii,ji,ki,id,jd,kd
    sll_real64 :: sm

    sm = 0.0_f64
!$omp parallel default(shared) private(i,j,k,ii,ji,ki,id,jd,kd) reduction(+ : sm)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do k =1, local_cells(3)
       do j = 1, local_cells(2)

          do kd = 1, degree(3)
             ki = (k-1)*degree(3)+kd
             do jd = 1, degree(2)
                ji = (j-1)*degree(2)+jd
                do i=1, local_cells(1)
                   do id =1, degree(1)
                      ii = (i-1)*degree(1)+id
                      sm = sm + f3d(ii,ji,ki) * weights(id,jd,kd)
                   end do
                end do
             end do
          end do

       end do
    end do
!$omp end do
!$omp end parallel

    val = sm
  end subroutine sll_s_gauss_quadrature_cells_3d


  subroutine sll_s_uniform_to_gauss( interpolator, ncell_uniform, cell, output_array )
    type(sll_t_dg_interpolator_1d), intent( in ) :: interpolator
    sll_int32, intent( in ) :: ncell_uniform
    sll_int32, intent( out ) :: cell(:)
    sll_real64, intent( out ) :: output_array(:,:)

    sll_real64, allocatable :: qpoints(:)
    sll_int32 :: j, index, degree, iinterp, ierr

    degree = interpolator%degree
    iinterp = ((ncell_uniform+1)/2)

    SLL_ALLOCATE(qpoints( -iinterp+1:ncell_uniform+iinterp-1), ierr)

    do j=-iinterp+1, ncell_uniform+iinterp-1
       qpoints(j) = real(j,f64)/real(ncell_uniform,f64)
    end do

    do j=1,degree
       cell(j) = floor(interpolator%quadrature_points(j)*real(ncell_uniform, f64))
       ! For Gauss-Lobatto, we have the point 1. Here, we shall pin it to the last interval of the cell instead of the first of the next cell)
       if ( cell(j) > ncell_uniform-1 ) then
          cell(j) = ncell_uniform-1
       end if
    end do

    do j=1,degree
       do index=1,2*iinterp
          output_array(index,j) = lagrange_poly( qpoints(cell(j)-iinterp+1:cell(j)+iinterp), index, interpolator%quadrature_points(j) )
       end do
    end do

    cell = cell+1

    SLL_DEALLOCATE_ARRAY(qpoints, ierr)
  end subroutine sll_s_uniform_to_gauss


  subroutine sll_s_uniform_to_gauss_staggered( interpolator, ncell_uniform, cell, output_array )
    type(sll_t_dg_interpolator_1d), intent( in ) :: interpolator
    sll_int32, intent( in ) :: ncell_uniform
    sll_int32, intent( out ) :: cell(:)
    sll_real64, intent( out ) :: output_array(:,:)

    sll_real64, allocatable :: qpoints(:)
    sll_real64 :: deltax
    sll_int32 :: j, index, degree, iinterp, ierr

    degree = interpolator%degree
    iinterp = ((ncell_uniform+1)/2)
    deltax = 1.0_f64/real(ncell_uniform,f64)

    SLL_ALLOCATE(qpoints( -iinterp:ncell_uniform+iinterp-1), ierr)

    do j=-iinterp, ncell_uniform+iinterp-1
       qpoints(j) = (real(j,f64)+0.5_f64)*deltax
    end do

    do j=1,degree
       cell(j) = floor((interpolator%quadrature_points(j)+deltax*0.5_f64)*real(ncell_uniform, f64))
       ! For Gauss-Lobatto, we have the point 1. Here, we shall pin it to the last interval of the cell instead of the first of the next cell)
       !if ( cell(j) > ncell_uniform-1 ) then
       !   cell(j) = ncell_uniform-1
       !end if
    end do

    do j=1,degree
       do index=1,2*iinterp
          output_array(index,j) = lagrange_poly( qpoints(cell(j)-iinterp:cell(j)+iinterp-1), index, interpolator%quadrature_points(j) )
       end do
    end do

    !cell = cell+1

    SLL_DEALLOCATE_ARRAY(qpoints, ierr)
  end subroutine sll_s_uniform_to_gauss_staggered


  !> Evaluate Lagrange polynomial
  function lagrange_poly( quadrature_points, index, x ) result(r)
    sll_real64, intent( in ) :: quadrature_points(:)
    sll_int32, intent( in ) :: index
    sll_real64, intent( in ) :: x
    sll_real64 :: r

    sll_int32 :: j, degree

    degree = size(quadrature_points)

    r = 1.0_f64
    do j=1,degree
       if (j .ne. index) then
          r = r * (x-quadrature_points(j))/(quadrature_points(index)-quadrature_points(j))
       end if
    end do

  end function lagrange_poly


  subroutine sll_s_field_uniform_to_gauss( field_in, field_out, &
                                           deg, ldeg, cell_uniform_to_gauss, eval_uniform_to_gauss, &
                                           topo_3d, decomp_3d, &
                                           topo_6d, &  ! optional argument for debug purposes only
                                           decomp_6d)  ! optional argument for debug purposes only
    sll_real64, intent( inout ) :: field_in(:,:,:)  ! input array, is used as a temporary buffer and overwritten!
    sll_real64, intent( inout ) :: field_out(:,:,:)
    sll_int32, intent( in ) :: deg(3)
    sll_int32, intent( in ) :: ldeg(3)
    sll_int32, intent( in ) :: cell_uniform_to_gauss(:,:)
    sll_real64, intent( in ) :: eval_uniform_to_gauss(:,:,:)
    type(sll_t_cartesian_topology_3d), target, intent(in) :: topo_3d
    type(sll_t_decomposition_slim_3d), target, intent(inout) :: decomp_3d
    type(sll_t_cartesian_topology_6d), target, intent(in), optional :: topo_6d
    type(sll_t_decomposition_slim_6d), target, intent(in), optional :: decomp_6d

    sll_int32 :: i, j, k, ii, jj, kk, npoints(3), ncells(3)
    sll_int32 :: cell, first, l, ind, hw_left, hw_right, id, dd
    sll_int32 :: stripe_mn, stripe_mx
    sll_real64, allocatable :: stripe(:)
    HALO_DTYPE, pointer :: l_buf(:,:,:), r_buf(:,:,:)

    logical, parameter :: use_mpi_implementation = .true.

    npoints(1) = size(field_out,1)
    npoints(2) = size(field_out,2)
    npoints(3) = size(field_out,3)
    ncells = npoints/deg

    ! --- runtime checks during MPI parallallelization work, may be disabled later ---
    do i=1,3
      SLL_ASSERT_ALWAYS(decomp_3d%local%nw(i) == npoints(i))
      SLL_ASSERT_ALWAYS(decomp_3d%local%n_cells(i) == ncells(i))
      if (present(topo_6d)) then
        SLL_ASSERT_ALWAYS(topo_3d%coords(i) == topo_6d%coords(i))
      endif
      if (present(decomp_6d)) then
        SLL_ASSERT_ALWAYS(decomp_3d%local%mn(i) == decomp_6d%local%mn(i))
        SLL_ASSERT_ALWAYS(decomp_3d%local%mx(i) == decomp_6d%local%mx(i))
        SLL_ASSERT_ALWAYS(decomp_3d%local%nw(i) == decomp_6d%local%nw(i))
      endif
    enddo

    if (use_mpi_implementation) then

      id = 1  ! dimension
      ! determine min stripe index (see loops below)
      i = 1
      cell = (i-1)*deg(id)
      dd = 1
      first = cell + cell_uniform_to_gauss(dd,id)-ldeg(id)/2
      stripe_mn = first + 1
      ! determine max stripe index (see loops below)
      i = ncells(id)
      cell = (i-1)*deg(id)
      dd = deg(id)
      first = cell + cell_uniform_to_gauss(dd,id)-ldeg(id)/2
      stripe_mx = first + ldeg(id)
      ! calculate halo widths based on the stripe extents
      hw_left = 1 - stripe_mn
      hw_right = stripe_mx - npoints(id)

      call sll_f_apply_halo_exchange_slim_3d_real64(topo_3d, decomp_3d, field_in, id, hw_left, hw_right)
      ! UGLY: use 1-based indexing instead of the global indexing we're using with domain decomposition
      l_buf(1:,1:,1:) => decomp_3d%local%halo_left%buf
      r_buf(1:,1:,1:) => decomp_3d%local%halo_right%buf

!$omp parallel default(shared) private(i, j, k, cell, dd, stripe, first)
      allocate(stripe(stripe_mn:stripe_mx))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,npoints(3)
         do j=1,npoints(2)
            stripe(    1-hw_left:0                   ) = l_buf(:,j,k)
            stripe(            1:npoints(id)         ) = field_in(:,j,k)
            stripe(npoints(id)+1:npoints(id)+hw_right) = r_buf(:,j,k)
            do i=1,ncells(id)
               cell = (i-1)*deg(id)
               do dd=1,deg(id)
                 first = cell + cell_uniform_to_gauss(dd,id)-ldeg(id)/2
                 field_out(cell+dd,j,k) = sum( eval_uniform_to_gauss(:,dd,id) * stripe(first+1:first+ldeg(id)) )
               end do
            end do
         end do
      end do
!$omp end do
      deallocate(stripe)
!$omp end parallel


      id = 2  ! dimension
      ! determine min stripe index (see loops below)
      i = 1
      cell = (i-1)*deg(id)
      dd = 1
      first = cell + cell_uniform_to_gauss(dd,id)-ldeg(id)/2
      stripe_mn = first + 1
      ! determine max stripe index (see loops below)
      i = ncells(id)
      cell = (i-1)*deg(id)
      dd = deg(id)
      first = cell + cell_uniform_to_gauss(dd,id)-ldeg(id)/2
      stripe_mx = first + ldeg(id)
      ! calculate halo widths based on the stripe extents
      hw_left = 1 - stripe_mn
      hw_right = stripe_mx - npoints(id)

      call sll_f_apply_halo_exchange_slim_3d_real64(topo_3d, decomp_3d, field_out, id, hw_left, hw_right)
      ! UGLY: use 1-based indexing instead of the global indexing we're using with domain decomposition
      l_buf(1:,1:,1:) => decomp_3d%local%halo_left%buf
      r_buf(1:,1:,1:) => decomp_3d%local%halo_right%buf

!$omp parallel default(shared) private(i, j, k, cell, dd, stripe, first)
      allocate(stripe(stripe_mn:stripe_mx))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do k=1,npoints(3)
         do i=1,npoints(1)
            stripe(    1-hw_left:0                   ) = l_buf(i,:,k)
            stripe(            1:npoints(id)         ) = field_out(i,:,k)
            stripe(npoints(id)+1:npoints(id)+hw_right) = r_buf(i,:,k)
            do j=1,ncells(id)
               cell = (j-1)*deg(id)
               do dd=1,deg(id)
                  first = cell + cell_uniform_to_gauss(dd,id)-ldeg(id)/2
                  field_in(i,cell+dd,k) = sum( eval_uniform_to_gauss(:,dd,id) * stripe(first+1:first+ldeg(id)) )
               end do
            end do
         end do
      end do
!$omp end do
      deallocate(stripe)
!$omp end parallel


      id = 3
      ! determine min stripe index (see loops below)
      i = 1
      cell = (i-1)*deg(id)
      dd = 1
      first = cell + cell_uniform_to_gauss(dd,id)-ldeg(id)/2
      stripe_mn = first + 1
      ! determine max stripe index (see loops below)
      i = ncells(id)
      cell = (i-1)*deg(id)
      dd = deg(id)
      first = cell + cell_uniform_to_gauss(dd,id)-ldeg(id)/2
      stripe_mx = first + ldeg(id)
      ! calculate halo widths based on the stripe extents
      hw_left = 1 - stripe_mn
      hw_right = stripe_mx - npoints(id)

      call sll_f_apply_halo_exchange_slim_3d_real64(topo_3d, decomp_3d, field_in, id, hw_left, hw_right)
      ! UGLY: use 1-based indexing instead of the global indexing we're using with domain decomposition
      l_buf(1:,1:,1:) => decomp_3d%local%halo_left%buf
      r_buf(1:,1:,1:) => decomp_3d%local%halo_right%buf

!$omp parallel default(shared) private(i, j, k, cell, dd, stripe, first)
      allocate(stripe(stripe_mn:stripe_mx))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do j=1,npoints(2)
         do i=1,npoints(1)
            stripe(    1-hw_left:0                   ) = l_buf(i,j,:)
            stripe(            1:npoints(id)         ) = field_in(i,j,:)
            stripe(npoints(id)+1:npoints(id)+hw_right) = r_buf(i,j,:)
            do k=1,ncells(id)
               cell = (k-1)*deg(id)
               do dd=1,deg(id)
                  first = cell + cell_uniform_to_gauss(dd,id)-ldeg(id)/2
                  field_out(i,j,cell+dd) = sum( eval_uniform_to_gauss(:,dd,id) * stripe(first+1:first+ldeg(id)) )
               end do
            end do
         end do
      end do
!$omp end do
      deallocate(stripe)
!$omp end parallel

      if (allocated(decomp_3d%local%halo_right%buf)) then
        deallocate(decomp_3d%local%halo_right%buf)
      endif
      if (allocated(decomp_3d%local%halo_left%buf)) then
        deallocate(decomp_3d%local%halo_left%buf)
      endif
      l_buf => null()
      r_buf => null()

    else  ! --- previous sequential implementation with modulo-based periodic BCs below ---
      do k=1,npoints(3)
         do j=1,npoints(2)
            do i=1,ncells(1)
               cell = (i-1)*deg(1)
               do ii=1,deg(1)
                  first = cell + cell_uniform_to_gauss(ii,1)-ldeg(1)/2
                  if ((first .le. 0) .OR.  (first+ldeg(1) > npoints(1))) then
                     field_out(cell+ii,j,k) = 0.0_f64
                     do l=1,ldeg(1)
                        ind = modulo(first+l-1,npoints(1))+1
                        field_out(cell+ii,j,k) = field_out(cell+ii,j,k) + &
                             eval_uniform_to_gauss(l,ii,1) * &
                             field_in(ind,j,k)
                     end do
                  else
                     field_out(cell+ii,j,k) = sum(eval_uniform_to_gauss(:,ii,1) * &
                          field_in(first+1:first+ldeg(1),j,k))
                  end if
               end do
            end do
         end do
      end do
      do k=1,npoints(3)
         do i=1,npoints(1)
            do j=1,ncells(2)
               cell = (j-1)*deg(2)
               do jj=1,deg(2)
                  first = cell + cell_uniform_to_gauss(jj,2)-ldeg(2)/2
                  if ( (first .le. 0) .OR. (first+ldeg(2) > npoints(2))) then
                     field_in(i,cell+jj,k) = 0.0_f64
                     do l=1,ldeg(2)
                        ind = modulo(first+l-1,npoints(2))+1
                        ! field_in(i,cell+jj,k) = field_out(i,cell+jj,k) +
                        field_in(i,cell+jj,k) = field_in(i,cell+jj,k) + &
                             eval_uniform_to_gauss(l,jj,2) * field_out(i,ind,k)
                     end do
                  else
                     field_in(i,cell+jj,k) = sum(eval_uniform_to_gauss(:,jj,2) * &
                          field_out(i,first+1:first+ldeg(2),k))
                  end if
               end do
            end do
         end do
      end do
      do j=1,npoints(2)
         do i=1,npoints(1)
            do k=1,ncells(3)
               cell = (k-1)*deg(3)
               do kk=1,deg(3)
                  first = cell + cell_uniform_to_gauss(kk,3)-ldeg(3)/2
                  if ( (first .le. 0) .OR. (first+ldeg(3) > npoints(3))) then
                     field_out(i,j,cell+kk) = 0.0_f64
                     do l=1,ldeg(3)
                        ind = modulo(first+l-1,npoints(3))+1
                        field_out(i,j,cell+kk) = field_out(i,j,cell+kk) + &
                             eval_uniform_to_gauss(l,kk,3) * field_in(i,j,ind)
                     end do
                  else
                     field_out(i,j,cell+kk) = sum(eval_uniform_to_gauss(:,kk,1) * &
                          field_in(i,j,first+1:first+ldeg(3)))
                  end if
               end do
            end do
         end do
      end do

    endif
  end subroutine sll_s_field_uniform_to_gauss


  
    subroutine sll_s_set_source_term_v ( local_sizes, etas, source_term_v )
      sll_int32, intent( in ) :: local_sizes(6)
      type(sll_t_array), intent( in ) :: etas(6)
      sll_real64, intent( inout ) :: source_term_v(1:,1:,1:)

      sll_int32 :: l,m,n
      sll_real64 :: v1, v2, v3, factor

      factor = 1.0_f64/(sll_p_twopi)**1.5_f64

!$omp parallel default(shared) private(l,m,n, v1, v2, v3)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do n=1, local_sizes(6)
         do m=1, local_sizes(5)
            v3 = etas(6)%vals(n)**2
            v2 = etas(5)%vals(m)**2
            do l=1, local_sizes(4)
               v1 = etas(4)%vals(l)**2
               source_term_v(l,m,n) = &
                    exp( - 0.5_f64 * (v1 + v2 + v3 ) ) * factor * &
                     ((v1+v2+v3)*0.5_f64-1.5_f64)
          end do
       end do
    end do
!$omp end do
!$omp end parallel
      
      
    end subroutine sll_s_set_source_term_v

    
    subroutine sll_s_set_collision_mask_v ( local_sizes, etas, collision_mask, vmax )
      sll_int32, intent( in ) :: local_sizes(6)
      type(sll_t_array), intent( in ) :: etas(6)
      sll_real64, intent( inout ) :: collision_mask(1:,1:,1:)
    sll_real64, intent(in) ::  vmax
      sll_int32 :: l,m,n
      sll_real64 :: v1, v2, v3, factor

      factor = 1.0_f64/(sll_p_twopi)**1.5_f64

!$omp parallel default(shared) private(l,m,n, v1, v2, v3)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
      do n=1, local_sizes(6)
         do m=1, local_sizes(5)
            v2 = etas(5)%vals(m)**2
            do l=1, local_sizes(4)
               v1 = etas(4)%vals(l)**2
               collision_mask(l,m,n) = &
                    0.5_f64* (erf(-10.0_f64* ((v1+v2)/(0.5*vmax)**2 - 1)) - 1.0_f64)
          end do
       end do
    end do
!$omp end do
!$omp end parallel
      
      
    end subroutine sll_s_set_collision_mask_v
  

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
