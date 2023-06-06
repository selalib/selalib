! ---
!  6D domain decomposition demo program.
! ---
!  2016
!    Klaus Reuter, khr@mpcdf.mpg.de
!    Max Planck Computing and Data Facility (MPCDF)
!    Garching, Germany
! ---

program test_decomposition_advanced
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_mpi, only: &
    mpi_land

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective, &
    sll_s_collective_reduce_logical

  use sll_m_decomposition, only: &
    sll_t_cartesian_topology_6d, &
    sll_f_new_cartesian_topology_6d

  use sll_m_decomposition_advanced

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! --- variable section

  sll_int32, parameter :: nd = 6
  sll_int32 :: procs_per_dimension(nd)
  logical :: periodic(nd), check(2*nd+1)
  type(sll_t_cartesian_topology_6d), pointer :: topology
  type(sll_t_decomposition), pointer :: decomposition
  sll_int32 :: ierr, i_buffer_left, i_buffer_right
  sll_int32 :: world_size, my_rank
  sll_int32 :: global_grid_points_per_dimension(nd)
  sll_int32 :: i, j, id, hw_left, hw_right, block_dim, n_blocks, n_buffers, n_buffer_pairs, width
  sll_real64, dimension(:,:,:,:,:,:), pointer :: f6d
  sll_real64 :: val_left, val_right
  logical :: chk_tmp, chk_all, chk_send_buf(1), chk_recv_buf(1)


  ! --- executable section


  call sll_s_boot_collective()
  world_size = sll_f_get_collective_size(sll_v_world_collective)
  my_rank = sll_f_get_collective_rank(sll_v_world_collective)


  !> Create a distribution of the MPI processes over the dimensions.
  !> There is a convenience function provided by MPI as well, however
  !> we do this by hand for the moment to have better control.
  select case(world_size)
  case(1)
    procs_per_dimension(:) = 1
  case(2)
    procs_per_dimension(1) = 2
    procs_per_dimension(2:6) = 1
  case(4)
    procs_per_dimension(1:2) = 2
    procs_per_dimension(3:6) = 1
  case(8)
    procs_per_dimension(1:3) = 2
    procs_per_dimension(4:6) = 1
  case(16)
    procs_per_dimension(1:4) = 2
    procs_per_dimension(5:6) = 1
  case(32)
    procs_per_dimension(1:5) = 2
    procs_per_dimension(6) = 1
  case(64)
    procs_per_dimension(:) = 2
  case(96)
    procs_per_dimension(1) = 3
    procs_per_dimension(2:6) = 2
  case(128)
    procs_per_dimension(1) = 4
    procs_per_dimension(2:6) = 2
  case default
    if (my_rank == 0) &
      write(*,*) "WARNING: No topology implemented for ", world_size, " processes."
    chk_all = .false.
    go to 100
  end select


  !> (1) Create a topology based on the distribution.  In addition, a nd boolean array
  !> is passed to indicate if there is a periodic BC associated to a certain dimension.
  periodic(:) = .true.
  topology => &
    sll_f_new_cartesian_topology_6d(sll_v_world_collective, procs_per_dimension, periodic)

  !> (2) Create a domain decomposition, create a global array distributed over the
  !> MPI processes.  The size of the local box of the array is passed via an nd array
  !> as well as the widh of the halo to be exchanged between neighboring processes.
  !> The decomposition object contains all the information necessary to allocate and
  !> access the local array (ie it has all the indices), see below.
  global_grid_points_per_dimension(:) = 8

  ! --- initial try-out phase ---
  if (.false.) then
    decomposition => &
      sll_f_new_cartesian_domain_decomposition(topology, global_grid_points_per_dimension, nd)
    ! --- define blocks logically
    id = 1
    n_blocks = 2
    block_dim = 6
    call sll_s_dd_define_blocks(decomposition, id, n_blocks, block_dim)
    ! ---
    write(*,*) "--- basic decomposition index values ---"
    write(*,*) "--- local ---"
    write(*,*) decomposition%local%mn
    write(*,*) decomposition%local%mx
    write(*,*) decomposition%local%nw
    write(*,*) "--- block 1 ---"
    write(*,*) decomposition%local%dimension(id)%block(1)%mn
    write(*,*) decomposition%local%dimension(id)%block(1)%mx
    write(*,*) decomposition%local%dimension(id)%block(1)%nw
    write(*,*) "--- block 2 ---"
    write(*,*) decomposition%local%dimension(id)%block(2)%mn
    write(*,*) decomposition%local%dimension(id)%block(2)%mx
    write(*,*) decomposition%local%dimension(id)%block(2)%nw
    ! ---

    n_buffers = 2
    width = 2
    do i=1, n_blocks
      do j=1, n_buffers
        call sll_s_dd_allocate_buffer(decomposition, id, i, j, width)
      end do
    end do

    write(*,*) "--- buffer index values ---"
    do i=1, n_blocks
      write(*,*) "--- block ", i
      do j=1, n_buffers
        write(*,*) "--- buffer ", j
        write(*,*) decomposition%local%dimension(id)%block(i)%buffer(j)%mn
        write(*,*) decomposition%local%dimension(id)%block(i)%buffer(j)%mx
        write(*,*) decomposition%local%dimension(id)%block(i)%buffer(j)%nw
      end do
    end do

    do i=1, n_blocks
      do j=1, n_buffers
        call sll_s_dd_deallocate_buffer(decomposition, id, i, j)
      end do
    end do

    call sll_s_deallocate_cartesian_domain_decomposition(decomposition)
  endif

  ! --- start over
  decomposition => &
    sll_f_new_cartesian_domain_decomposition(topology, global_grid_points_per_dimension, nd)

  !> Allocate the local array.
  allocate(f6d(decomposition%local%mn(1):decomposition%local%mx(1), &
               decomposition%local%mn(2):decomposition%local%mx(2), &
               decomposition%local%mn(3):decomposition%local%mx(3), &
               decomposition%local%mn(4):decomposition%local%mx(4), &
               decomposition%local%mn(5):decomposition%local%mx(5), &
               decomposition%local%mn(6):decomposition%local%mx(6)),&
               stat=ierr)
  SLL_ASSERT( ierr == 0 )

  !> Fill the array with the MPI rank.
  f6d = real(my_rank, kind=f64)

  !> (3) Exercise and check data exchange between next neighbors.
  n_blocks = 2
  block_dim = 6  ! dimension to be blocked over
  do id = 1, 3
    call sll_s_dd_define_blocks(decomposition, id, n_blocks, block_dim)
  enddo
  block_dim = 3  ! dimension to be blocked over
  do id = 4, 6
    call sll_s_dd_define_blocks(decomposition, id, n_blocks, block_dim)
  enddo
  chk_all = .true.
  width = 2
  n_buffer_pairs = 1
  do id = 1, 1  ! ###
    do i=1, n_blocks
      do j=1, n_buffer_pairs
        i_buffer_left = 2*j-1
        i_buffer_right = 2*j

        call sll_s_dd_allocate_buffer(decomposition, id, i, i_buffer_left, width)
        call sll_s_dd_allocate_buffer(decomposition, id, i, i_buffer_right, width)

        call sll_s_post_halo_exchange_real64(topology, decomposition, f6d, id, i, i_buffer_left)
        call sll_s_post_halo_exchange_real64(topology, decomposition, f6d, id, i, i_buffer_right)
        call sll_s_wait_halo_exchange(topology, decomposition, id, i, i_buffer_left)
        call sll_s_wait_halo_exchange(topology, decomposition, id, i, i_buffer_right)

        ! -- ADVECTION WOULD TAKE PLACE HERE --

        ! --- check halo buffer exchange
        val_left = real(topology%neighbors(2*id-1), kind=f64)
        val_right = real(topology%neighbors(2*id), kind=f64)
        chk_tmp = all(sll_f_get_mem_6d_from_buffer_obj(decomposition%local%dimension(id)%block(i)%buffer(i_buffer_left)) == val_left) .and. &
                  all(sll_f_get_mem_6d_from_buffer_obj(decomposition%local%dimension(id)%block(i)%buffer(i_buffer_right)) == val_right)
        if (.not. chk_tmp) then
          chk_all = .false.
        endif

        call sll_s_dd_deallocate_buffer(decomposition, id, i, i_buffer_left)
        call sll_s_dd_deallocate_buffer(decomposition, id, i, i_buffer_right)
      end do
    end do
  enddo

  chk_send_buf(1) = chk_all
  call sll_s_collective_reduce_logical(sll_v_world_collective, chk_send_buf, 1, MPI_LAND, 0, chk_recv_buf)
  chk_all = chk_recv_buf(1)

  SLL_DEALLOCATE(topology, ierr)
  call sll_s_deallocate_cartesian_domain_decomposition(decomposition)
  deallocate( f6d )

100 continue

  if (my_rank == 0) then
    write(*,"(a)",advance="no") "test_decomposition_advanced unit test : "
    if (chk_all .eqv. .true.) then
      write(*,*) "PASSED"
    else
      write(*,*) "FAILED"
    endif
  endif

101   call sll_s_halt_collective()
end program test_decomposition_advanced
