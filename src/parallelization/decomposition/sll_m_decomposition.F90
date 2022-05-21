!**************************************************************
!
!  Domain decomposition module for Selalib (1D, ..., 6d)
!
!  Some design aspects are borrowed from sll_remap.F90.
!
!  Note:
!  For memory reasons, this code is implemented separately
!  from the comm/port model provided by the existing module
!  located at "point_to_point_communications".
!  In 6d, we cannot afford allocating 2*12=24 buffers to
!  implement the non-blocking ghost cell exchange using
!  the comm/port model.
!
!**************************************************************
!> @ingroup decomposition
!> @brief
!> Module providing data structures and tools to implement domain decompositions.
!> @author
!> Klaus Reuter, Max Planck Computing and Data Facility (MPCDF)


! Preprocessor macro:  Use single precision for the halo exchange, ie.
! basically a poor man's lossy compression to speed up MPI communication.
! NOTE: The define should be passed by the build system, please see
! the build script <mpcdf.sh> for an example.
!!!#define USE_HALO_REAL32
#ifdef USE_HALO_REAL32
#define HALO_DTYPE sll_real32
#else
#define HALO_DTYPE sll_real64
#endif


module sll_m_decomposition
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_t_collective_t

  use sll_mpi, only: &
    mpi_cart_create, &
    mpi_cart_get, &
    mpi_cart_sub, &
    mpi_cart_shift, &
    mpi_double_precision, &
    mpi_real, &
    mpi_send, &
    mpi_recv, &
    mpi_sendrecv, &
    mpi_status_ignore, &
    mpi_success, &
    mpi_allreduce, &
    mpi_integer, &
    mpi_sum, &
    mpi_group_incl, &
    mpi_comm_create, &
    mpi_group_free, &
    mpi_barrier, &
    mpi_finalize, &
    mpi_byte, &
    mpi_max_processor_name, &
    mpi_get_processor_name, &
    mpi_scatter, &
    mpi_character, &
    mpi_comm_split, &
    mpi_comm_free

  use sll_m_utilities, only: sll_f_query_environment

#ifdef _OPENMP
  use omp_lib
#define OMP_COLLAPSE collapse(2)
#define OMP_SCHEDULE schedule(static)
!#define OMP_SCHEDULE schedule(dynamic)
#endif

#ifdef USE_FMEMPOOL
  use fmempool
#endif

  use sll_m_compression

  implicit none

  public :: &
    sll_o_apply_halo_exchange, &
    sll_t_cartesian_topology_3d, &
    sll_t_cartesian_topology_6d, &
    sll_s_copy_array_to_buffer_6d_real64, &
    sll_t_decomposition_3d, &
    sll_t_decomposition_6d, &
    sll_s_mpi_sendrecv_compressed_core, &
    sll_f_new_cartesian_domain_decomposition_3d, &
    sll_f_new_cartesian_domain_decomposition_6d, &
    sll_f_new_cartesian_topology_3d, &
    sll_f_new_cartesian_topology_6d, &
    sll_f_new_cartesian_topology_3d_from_6d, &
    sll_f_new_cartesian_topology_3d_orthogonal, &
    sll_f_select_dim, &
    sll_f_set_process_grid, &
    sll_t_decomposition_slim_6d, &
    sll_f_new_cartesian_domain_decomposition_slim_6d, &
    sll_f_apply_halo_exchange, &
    sll_s_apply_halo_exchange_slim_6d_real64, &
    sll_f_apply_halo_exchange_slim_6d_real64, &
    sll_s_allocate_bc_buffers_6d, &
    sll_s_allocate_bc_buffers_6d_part, &
    sll_s_deallocate_bc_buffers, &
    sll_s_apply_bc_exchange_slim_6d_real64

  public :: &
       sll_f_apply_halo_exchange_slim_3d_real64, &
       sll_f_new_cartesian_domain_decomposition_slim_3d, &
       sll_f_new_cartesian_cell_domain_decomposition_slim_6d, &
       sll_t_decomposition_slim_3d

  private

  !> @brief    Information on the 6D cartesian process topology.
  type :: sll_t_cartesian_topology_6d
    ! topology-associated MPI communicator
    sll_int32 :: comm
    ! MPI rank
    sll_int32 :: rank
    ! MPI size
    sll_int32 :: nprocs
    ! store optional information (for debugging purposes)
    sll_int32 :: info = -1
    ! array containing the number of MPI processes to use along the n-th dimension
    sll_int32 :: procs(6)
    ! array indicating if a periodic boundary condition exists at the n-th dimension
    logical   :: periodic(6)
    ! coordinates of the current MPI process within the cartesian topology
    sll_int32 :: coords(6)
    ! MPI ranks of the topological neighbors of the current MPI process
    sll_int32 :: neighbors(12)
  end type sll_t_cartesian_topology_6d

  !> @brief    Information on the 3D cartesian process topology.
  type :: sll_t_cartesian_topology_3d
    ! topology-associated MPI communicator
    sll_int32 :: comm
    ! MPI rank
    sll_int32 :: rank
    ! MPI size
    sll_int32 :: nprocs
    ! store optional information (for debugging purposes)
    sll_int32 :: info = -1
    ! array containing the number of MPI processes to use along the n-th dimension
    sll_int32 :: procs(3)
    ! array indicating if a periodic boundary condition exists at the n-th dimension
    logical   :: periodic(3)
    ! coordinates of the current MPI process within the cartesian topology
    sll_int32 :: coords(3)
    ! MPI ranks of the topological neighbors of the current MPI process
    sll_int32 :: neighbors(6)
  end type sll_t_cartesian_topology_3d

  !> @brief    6D decomposition, index limits local to an MPI process.
  type :: decomposition_local_6d
    ! --- primary: local array extents usable for computation, width of the halo
    sll_int32, dimension(:) :: mn(6)  ! min index, w/o halo
    sll_int32, dimension(:) :: mx(6)  ! max index, w/o halo
    sll_int32, dimension(:) :: hw(6)  ! halo width
    ! --- derived: total allocated extents of the local array, including halo cells
    sll_int32, dimension(:) :: lo(6)  ! min index, w/ halo
    sll_int32, dimension(:) :: hi(6)  ! max index, w/ halo
    ! --- derived: auxiliary variables, array sizes
    sll_int32, dimension(:) :: nw(6)  ! net width of the array, w/o halo
    sll_int32, dimension(:) :: gw(6)  ! gross width of the array, w/ halo
    ! --- derived: ranges to copy between send (tx) and receive (rx) buffers and arrays
    sll_int32, dimension(:) :: tx_lolo(6)  ! lower/left neighbor, lower bound
    sll_int32, dimension(:) :: tx_lohi(6)  ! lower/left neighbor, upper bound
    sll_int32, dimension(:) :: tx_hilo(6)  ! upper/right neighbor, lower bound
    sll_int32, dimension(:) :: tx_hihi(6)  ! upper/right neighbor, upper bound
    sll_int32, dimension(:) :: rx_lolo(6)
    sll_int32, dimension(:) :: rx_lohi(6)
    sll_int32, dimension(:) :: rx_hilo(6)
    sll_int32, dimension(:) :: rx_hihi(6)
  end type decomposition_local_6d


! --- "slim" halo redesign
  type :: halo_buffer_6d
    sll_int32, dimension(:) :: mn(6)  ! min index of buf
    sll_int32, dimension(:) :: mx(6)  ! max index of buf
    sll_int32, dimension(:) :: nw(6)  ! net width of buf
#ifdef USE_FMEMPOOL
    HALO_DTYPE, pointer :: buf(:,:,:,:,:,:) => null() ! halo buffer
#else
    HALO_DTYPE, dimension(:,:,:,:,:,:), allocatable :: buf
#endif
  end type halo_buffer_6d


  type :: halo_buffer_3d
    sll_int32, dimension(:) :: mn(3)  ! min index of buf
    sll_int32, dimension(:) :: mx(3)  ! max index of buf
    sll_int32, dimension(:) :: nw(3)  ! net width of buf
    HALO_DTYPE, dimension(:,:,:), allocatable :: buf ! halo buffer
  end type halo_buffer_3d


  !> @brief    6D decomposition, "slim" redesign with dynamic halo cells
  type :: decomposition_slim_local_6d
    ! --- primary: local array extents usable for computation, width of the halo
    sll_int32, dimension(:) :: mn(6)  ! min index, w/o halo
    sll_int32, dimension(:) :: mx(6)  ! max index, w/o halo
    sll_int32, dimension(:) :: nw(6)  ! net width of the array, w/o halo
    ! --- dynamic part
    sll_int32 :: id  ! dimension to which the halos currently apply (valid: 1, ..., 6; invalid: -1)
    type(halo_buffer_6d) :: halo_left
    type(halo_buffer_6d) :: halo_right
    ! --- boundary condition buffers for spline interpolation
    HALO_DTYPE, allocatable :: bc_left_send(:,:,:,:,:)
    HALO_DTYPE, allocatable :: bc_right_send(:,:,:,:,:)
    HALO_DTYPE, allocatable :: bc_left(:,:,:,:,:)
    HALO_DTYPE, allocatable :: bc_right(:,:,:,:,:)
    ! --- cell (DG) interpolation-specific entries
    sll_int32, dimension(:) :: mn_cell(6)  ! min cell
    sll_int32, dimension(:) :: mx_cell(6)  ! max cell
    sll_int32, dimension(:) :: n_cells(6)  ! local number of cells
  end type decomposition_slim_local_6d



  type :: decomposition_slim_local_3d
    ! --- primary: local array extents usable for computation, width of the halo
    sll_int32, dimension(:) :: mn(3)  ! min index, w/o halo
    sll_int32, dimension(:) :: mx(3)  ! max index, w/o halo
    sll_int32, dimension(:) :: nw(3)  ! net width of the array, w/o halo
    ! --- dynamic part
    sll_int32 :: id  ! dimension to which the halos currently apply (valid: 1, ..., 3; invalid: -1)
    type(halo_buffer_3d) :: halo_left
    type(halo_buffer_3d) :: halo_right
    ! --- cell (DG) interpolation-specific entries
    sll_int32, dimension(:) :: mn_cell(6)  ! min cell
    sll_int32, dimension(:) :: mx_cell(6)  ! max cell
    sll_int32, dimension(:) :: n_cells(6)  ! local number of cells
  end type decomposition_slim_local_3d



  !> @brief    3D decomposition, index limits local to an MPI process.
  type :: decomposition_local_3d
    ! --- primary: local array extents usable for computation, width of the halo
    sll_int32, dimension(:) :: mn(3)  ! min index, w/o halo
    sll_int32, dimension(:) :: mx(3)  ! max index, w/o halo
    sll_int32, dimension(:) :: hw(3)  ! halo width
    ! --- derived: total allocated extents of the local array, including halo cells
    sll_int32, dimension(:) :: lo(3)  ! min index, w/ halo
    sll_int32, dimension(:) :: hi(3)  ! max index, w/ halo
    ! --- derived: auxiliary variables, array sizes
    sll_int32, dimension(:) :: nw(3)  ! net width of the array, w/o halo
    sll_int32, dimension(:) :: gw(3)  ! gross width of the array, w/ halo
    ! --- derived: ranges to copy between send (tx) and receive (rx) buffers and arrays
    sll_int32, dimension(:) :: tx_lolo(3)  ! lower/left neighbor, lower bound
    sll_int32, dimension(:) :: tx_lohi(3)  ! lower/left neighbor, upper bound
    sll_int32, dimension(:) :: tx_hilo(3)  ! upper/right neighbor, lower bound
    sll_int32, dimension(:) :: tx_hihi(3)  ! upper/right neighbor, upper bound
    sll_int32, dimension(:) :: rx_lolo(3)
    sll_int32, dimension(:) :: rx_lohi(3)
    sll_int32, dimension(:) :: rx_hilo(3)
    sll_int32, dimension(:) :: rx_hihi(3)
  end type decomposition_local_3d

  !> @brief    6D decomposition, global array size information and local information.
  type :: sll_t_decomposition_6d
    sll_int32, dimension(:) :: global(6)
    type(decomposition_local_6d) :: local
  end type sll_t_decomposition_6d


  !> @brief    6D decomposition, slim redesign, global array size information and local information.
  type :: sll_t_decomposition_slim_6d
    sll_int32, dimension(:) :: global(6)
    sll_int32, dimension(:) :: n_cells(6)  ! local number of cells (DG)
    type(decomposition_slim_local_6d) :: local
  end type sll_t_decomposition_slim_6d


  !> @brief    3D decomposition, global array size information and local information.
  type :: sll_t_decomposition_3d
    sll_int32, dimension(:) :: global(3)
    type(decomposition_local_3d) :: local
  end type sll_t_decomposition_3d


  !> 3D slim decomposition
  type :: sll_t_decomposition_slim_3d
    sll_int32, dimension(:) :: global(3)
    type(decomposition_slim_local_3d) :: local
  end type sll_t_decomposition_slim_3d


  interface sll_o_new_cartesian_topology
    !module procedure sll_f_new_cartesian_topology_3d
    module procedure sll_f_new_cartesian_topology_6d
  end interface sll_o_new_cartesian_topology

  interface sll_o_new_cartesian_domain_decomposition
    module procedure sll_f_new_cartesian_domain_decomposition_6d
  end interface sll_o_new_cartesian_domain_decomposition

  ! --- DEPRECATED - should be sll_f_*, see below!
  interface sll_o_apply_halo_exchange
    module procedure sll_f_apply_halo_exchange_6d_real64
  end interface sll_o_apply_halo_exchange

  interface sll_f_apply_halo_exchange
    module procedure sll_f_apply_halo_exchange_6d_real64
    module procedure sll_f_apply_halo_exchange_slim_6d_real64
    module procedure sll_s_apply_halo_exchange_slim_6d_real64
  end interface sll_f_apply_halo_exchange


contains


  !> @brief  Returns a mpi rank table with the processes transposed.
  subroutine get_transposed_process_map(procs_per_dimension, rank_map)
    integer, parameter :: nd=6
    sll_int32, intent(in) :: procs_per_dimension(nd)
    sll_int32, intent(inout) :: rank_map(0:)
    integer :: i,j,k,l,m,n,rk
    integer, allocatable :: mpi_grid(:,:,:,:,:,:)

    allocate(mpi_grid(0:procs_per_dimension(1)-1,&
                      0:procs_per_dimension(2)-1,&
                      0:procs_per_dimension(3)-1,&
                      0:procs_per_dimension(4)-1,&
                      0:procs_per_dimension(5)-1,&
                      0:procs_per_dimension(6)-1))

    ! set up integer array in row major (C) order with MPI ranks ('rk'),
    ! mimicking a compact layout of MPI processes onto nodes
    rk = 0
    do i=0,procs_per_dimension(1)-1
      do j=0,procs_per_dimension(2)-1
        do k=0,procs_per_dimension(3)-1
          do l=0,procs_per_dimension(4)-1
            do m=0,procs_per_dimension(5)-1
              do n=0,procs_per_dimension(6)-1
                mpi_grid(i,j,k,l,m,n) = rk
                rk = rk + 1
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    ! fill rank map such that the original process grid is transposed
    rk = 0
    do n=0,procs_per_dimension(6)-1
      do m=0,procs_per_dimension(5)-1
        do l=0,procs_per_dimension(4)-1
          do k=0,procs_per_dimension(3)-1
            do j=0,procs_per_dimension(2)-1
              do i=0,procs_per_dimension(1)-1
                rank_map(rk) = mpi_grid(i,j,k,l,m,n)
                rk = rk + 1
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    deallocate(mpi_grid)
  end subroutine


  !> @brief  6D Cartesian topology constructor function
  function sll_f_new_cartesian_topology_6d(top_collective, procs_per_dimension, periodic)
    type(sll_t_cartesian_topology_6d), pointer :: sll_f_new_cartesian_topology_6d
    type(sll_t_collective_t), intent(in) :: top_collective
    integer, parameter :: nd=6
    sll_int32, intent(in) :: procs_per_dimension(nd)
    logical, intent(in) :: periodic(nd)
    ! reordering of MPI ranks with MPI_Cart_create()
    logical :: q_reorder
    sll_int32 :: i, ierr
    integer :: name_len, my_rank, num_ranks, fd, in1, in2, new_rank, comm_temp
    character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
    character(len=MPI_MAX_PROCESSOR_NAME), allocatable :: hostnames(:)
    integer, allocatable :: ranks_reordered(:)
    logical :: q_block, q_block_6d_py, q_transpose

    my_rank = top_collective%rank
    num_ranks = top_collective%size

    ! find out if and how we are supposed to modify the process grid,
    q_transpose = sll_f_query_environment("SLL_TRANSPOSE_PROCESS_GRID", .false.)
    q_block = sll_f_query_environment("SLL_USE_PROCESS_BLOCKING", .false.)
    q_reorder = sll_f_query_environment("SLL_MPI_CART_CREATE_REORDER", .false.)
    inquire(file="./block6d.py", exist=q_block_6d_py)

    ! catch non-useful logical cases
    if (((q_transpose).or.(q_block)).and.(q_reorder)) then
      q_reorder = .false.
    endif
    if ((q_transpose).and.(q_block)) then
      q_transpose = .false.
      q_block = .false.
    endif
    if ((q_block).and.(.not.(q_block_6d_py))) then
      if (my_rank == 0) then
        write(*,*) " MPI topology: disabled blocking due to missing <block6d.py>"
      endif
      q_block = .false.
    endif
    if ((q_block).and.(top_collective%size < 64)) then
      if (my_rank == 0) then
        write(*,*) " MPI topology: disabled blocking due to small number of processes"
      endif
      q_block = .false.
    endif

    SLL_ALLOCATE(sll_f_new_cartesian_topology_6d, ierr)

    sll_f_new_cartesian_topology_6d%procs = procs_per_dimension
    sll_f_new_cartesian_topology_6d%periodic = periodic

    if ((q_transpose).or.(q_block)) then
!       if (my_rank == 0) then
        allocate(ranks_reordered(0:num_ranks-1))
!       endif

      ! fill the 'ranks_reordered' array, either blocked or transposed
      if (q_block) then
        ! implementation transferred from './blocking_prototype/cart.f90'
        call MPI_Get_Processor_Name(hostname, name_len, ierr)
        if (my_rank == 0) then
          write(*,*) " MPI topology: blocked process topology"

          allocate(hostnames(0:num_ranks-1))

          ! collect hostname table on rank 0
          hostnames(0) = hostname
          do i=1,num_ranks-1
            call MPI_Recv(hostnames(i), MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                          i, 0, top_collective%comm, MPI_STATUS_IGNORE, ierr)
          enddo

          ! write dims and hostnames table to text file
          fd = 42
          open(unit=fd, file='rank_hostnames.dat', status='replace', action='write')
          write(fd,*) num_ranks
          write(fd,*) nd
          do i=1,nd
            write(fd,*) sll_f_new_cartesian_topology_6d%procs(i)
          enddo
          do i=0,num_ranks-1
            write(fd,'(I6,A,A)') i, " ", trim(hostnames(i))
          enddo
          close(fd)
          deallocate(hostnames)

          ! Reorder MPI ranks in a blocked way,
          ! this is currently done separately in a Python script.
          ! In case this turns out to be successful a Fortran implementation will be done.
          ! call system("python ./block6d.py")
          call execute_command_line("python ./block6d.py")

          ! read reordered rank numbers
          open(unit=fd, file='rank_reordered.dat', status='old', action='read')
          do i=0,num_ranks-1
            read(fd,*) in1, in2
            ranks_reordered(i) = in2
          enddo
          close(fd)
        else
          call MPI_Send(hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                        0, 0, top_collective%comm, ierr)
        endif
      else
        ! q_transpose branch
        if (my_rank == 0) then
          write(*,*) " MPI topology: transposed process topology"
          call get_transposed_process_map(procs_per_dimension, ranks_reordered)
        endif
      endif

      ! distribute the individual new rank to all the processes
      call MPI_Scatter(ranks_reordered, 1, MPI_INTEGER, &
                       new_rank,        1, MPI_INTEGER, &
                       0, top_collective%comm, ierr)
      if (my_rank == 0) then
        deallocate(ranks_reordered)
      endif

      ! renumber the MPI ranks
      call MPI_Comm_split(top_collective%comm, 0, new_rank, comm_temp, ierr)
      ! create a new cartesian topology based on the renumbered ranks
      call MPI_Cart_create(comm_temp, nd, &
                           sll_f_new_cartesian_topology_6d%procs,&
                           sll_f_new_cartesian_topology_6d%periodic,&
                           q_reorder,&
                           sll_f_new_cartesian_topology_6d%comm,&
                           ierr)
      call MPI_Comm_Free(comm_temp, ierr)
    else
      ! native Cartesian topology creation
      if ((my_rank == 0).and.(q_reorder)) then
        write(*,*) " MPI topology: MPI_Cart_create() _may_ reorder processes."
      endif
      call MPI_Cart_create(top_collective%comm, nd,&
                           sll_f_new_cartesian_topology_6d%procs,&
                           sll_f_new_cartesian_topology_6d%periodic,&
                           q_reorder,&
                           sll_f_new_cartesian_topology_6d%comm,&
                           ierr)
      SLL_ASSERT(ierr == MPI_SUCCESS)
    endif

    call MPI_Comm_rank(sll_f_new_cartesian_topology_6d%comm,&
                       sll_f_new_cartesian_topology_6d%rank,&
                       ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    call MPI_Comm_size(sll_f_new_cartesian_topology_6d%comm,&
                       sll_f_new_cartesian_topology_6d%nprocs,&
                       ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    ! print new rank mapping
    if ((q_transpose).or.(q_block)) then
      write(*,'(A,I6,A,I6)') "             : ", my_rank, " --> ", sll_f_new_cartesian_topology_6d%rank
    endif

    ! query the coordinates of the current process within the cartesian topology
    sll_f_new_cartesian_topology_6d%coords = -1
    call MPI_Cart_get(sll_f_new_cartesian_topology_6d%comm, nd,&
                      sll_f_new_cartesian_topology_6d%procs,&
                      sll_f_new_cartesian_topology_6d%periodic,&
                      sll_f_new_cartesian_topology_6d%coords,&
                      ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    ! determine the neighbors within the cartesian topology
    sll_f_new_cartesian_topology_6d%neighbors = -1
    do i=1,nd
       call MPI_Cart_shift(sll_f_new_cartesian_topology_6d%comm, i-1, 1, &
                           sll_f_new_cartesian_topology_6d%neighbors(2*i-1), &
                           sll_f_new_cartesian_topology_6d%neighbors(2*i), &
                           ierr)
       SLL_ASSERT(ierr == MPI_SUCCESS)
    enddo
  end function sll_f_new_cartesian_topology_6d

  !> @brief 6D Cartesian topology destructor
  subroutine sll_s_deallocate_cartesian_topology_6d()
  end subroutine sll_s_deallocate_cartesian_topology_6d

  !> @brief  3D Cartesian topology constructor function
  function sll_f_new_cartesian_topology_3d(top_collective, procs_per_dimension, periodic)
    type(sll_t_cartesian_topology_3d), pointer :: sll_f_new_cartesian_topology_3d
    type(sll_t_collective_t), intent(in) :: top_collective
    integer, parameter :: nd=3
    sll_int32, intent(in) :: procs_per_dimension(nd)
    logical, intent(in) :: periodic(nd)
    ! disallow reordering of MPI ranks with MPI_Cart_create()
    logical, parameter :: reorder = .false.
    sll_int32 :: i, ierr

    SLL_ALLOCATE(sll_f_new_cartesian_topology_3d, ierr)

    sll_f_new_cartesian_topology_3d%procs = procs_per_dimension
    sll_f_new_cartesian_topology_3d%periodic = periodic

    ! create a cartesian process topology, return a new communicator
    call MPI_Cart_create(top_collective%comm, nd,&
                         sll_f_new_cartesian_topology_3d%procs,&
                         sll_f_new_cartesian_topology_3d%periodic,&
                         reorder,&
                         sll_f_new_cartesian_topology_3d%comm,&
                         ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)


    call MPI_Comm_rank(sll_f_new_cartesian_topology_3d%comm,&
                       sll_f_new_cartesian_topology_3d%rank,&
                       ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    call MPI_Comm_size(sll_f_new_cartesian_topology_3d%comm,&
                       sll_f_new_cartesian_topology_3d%nprocs,&
                       ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)


    ! query the coordinates of the current process within the cartesian topology
    sll_f_new_cartesian_topology_3d%coords = -1
    call MPI_Cart_get(sll_f_new_cartesian_topology_3d%comm, nd,&
                      sll_f_new_cartesian_topology_3d%procs,&
                      sll_f_new_cartesian_topology_3d%periodic,&
                      sll_f_new_cartesian_topology_3d%coords,&
                      ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    ! determine the neighbors within the cartesian topology
    sll_f_new_cartesian_topology_3d%neighbors = -1
    do i=1,nd
       call MPI_Cart_shift(sll_f_new_cartesian_topology_3d%comm, i-1, 1, &
                           sll_f_new_cartesian_topology_3d%neighbors(2*i-1), &
                           sll_f_new_cartesian_topology_3d%neighbors(2*i), &
                           ierr)
       SLL_ASSERT(ierr == MPI_SUCCESS)
    enddo
  end function sll_f_new_cartesian_topology_3d

  !> @brief 3D Cartesian topology destructor
  subroutine sll_s_deallocate_cartesian_topology_3d()
  end subroutine sll_s_deallocate_cartesian_topology_3d


  !> @brief  6D-->3D topology mapper, creates a 3D sub-topology from a 6D topology.
  function sll_f_new_cartesian_topology_3d_from_6d(t6d, keep_dim)
    type(sll_t_cartesian_topology_3d), pointer :: sll_f_new_cartesian_topology_3d_from_6d, t3d
    type(sll_t_cartesian_topology_6d), pointer, intent(in) :: t6d
    logical, dimension(6), intent(in) :: keep_dim
    sll_int32, parameter :: nd = 3
    sll_int32 :: i, j, ierr

    ! assert that we actually pick three dimensions from the 6D topology
    j = 0
    do i=1,6
      if (keep_dim(i) .eqv. .true.) j = j + 1
    end do
    SLL_ASSERT(j == nd)

    SLL_ALLOCATE(sll_f_new_cartesian_topology_3d_from_6d, ierr)
    t3d => sll_f_new_cartesian_topology_3d_from_6d  ! create a convenient pointer alias

    ! Create a 3d topology by projecting the velocity coordinates (dims 4,5,6)
    ! down to the spatial coordinates (dims 1,2,3).  As a result,
    ! we obtain a set of new communicators with 3D Cartesian topology.
    call MPI_Cart_sub(t6d%comm, keep_dim, &
                      t3d%comm, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    call MPI_Comm_rank(t3d%comm,&
                       t3d%rank,&
                       ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    call MPI_Comm_size(t3d%comm,&
                       t3d%nprocs,&
                       ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    ! Debug: Construct an identifier that is unique among the processors of the new 3D group.
    i = t6d%rank
    call MPI_Allreduce(i, t3d%info, 1,&
                       mpi_integer, mpi_sum,&
                       t3d%comm, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    ! Query the coordinates of the current process within the cartesian topology
    t3d%coords = -1
    t3d%procs = -1
    t3d%periodic = .false.
    call MPI_Cart_get(t3d%comm, nd,&
                      t3d%procs,&
                      t3d%periodic,&
                      t3d%coords,&
                      ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    ! Determine the neighbors within the cartesian topology
    t3d%neighbors = -1
    do i=1,nd
       call MPI_Cart_shift(t3d%comm, i-1, 1, &
                           t3d%neighbors(2*i-1), &
                           t3d%neighbors(2*i), &
                           ierr)
       SLL_ASSERT(ierr == MPI_SUCCESS)
    enddo
  end function sll_f_new_cartesian_topology_3d_from_6d


  !> @brief
  function sll_f_new_cartesian_topology_3d_orthogonal(topo_6d, topo_3d)
    type(sll_t_cartesian_topology_3d), pointer :: sll_f_new_cartesian_topology_3d_orthogonal, &
                                                  topo_3d_o  ! just a convenient pointer alias
    type(sll_t_cartesian_topology_6d), pointer, intent(in) :: topo_6d
    type(sll_t_cartesian_topology_3d), pointer, intent(in) :: topo_3d
    sll_int32 :: i, j, n3d, ierr
    sll_int32, allocatable :: topo_3d_rank_table(:)
    sll_int32, allocatable :: topo_6d_rank_grouped_by_3d_rank(:)
    sll_int32 :: group_6d, group_3d

    ! build a table that contains the 3D-ranks of each process in the 6D communicator
    allocate(topo_3d_rank_table(0:topo_6d%nprocs-1))
    topo_3d_rank_table(:) = -1
    call MPI_Allgather(topo_3d%rank, 1, MPI_INTEGER,&
                       topo_3d_rank_table, 1, MPI_INTEGER,&
                       topo_6d%comm, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)
    !write(*,*) topo_3d_rank_table

    ! determine the number of co-existing 3d communicators
    n3d = 0
    do i=0, topo_6d%nprocs-1
      if (topo_3d%rank == topo_3d_rank_table(i)) then
        n3d = n3d + 1
      end if
    end do

    ! build a table of the 6d ranks that have identical 3d ranks
    allocate(topo_6d_rank_grouped_by_3d_rank(n3d))
    topo_6d_rank_grouped_by_3d_rank(:) = -1
    j = 1
    do i=0, topo_6d%nprocs-1
      if (topo_3d%rank == topo_3d_rank_table(i)) then
        topo_6d_rank_grouped_by_3d_rank(j) = i
        j = j + 1
      end if
    end do
    !write(*,"(8192(I))") topo_6d%rank, topo_6d_rank_grouped_by_3d_rank


    ! obtain the group belonging to the 6d communicator
    call MPI_Comm_group(topo_6d%comm, group_6d, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    ! TODO : handle (re)ordering of the ranks in order to be compatible
    !        with sll_remap (ie data <--> rank mapping must be handled)

    ! form new groups defined by the 6d rank tables with identical 3d ranks (see above)
    call MPI_Group_incl(group_6d, n3d, topo_6d_rank_grouped_by_3d_rank, group_3d, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    ! create the new topology
    SLL_ALLOCATE(topo_3d_o, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)
    call MPI_Comm_create(topo_6d%comm, group_3d, topo_3d_o%comm, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)
    !
    call MPI_Comm_rank(topo_3d_o%comm,&
                       topo_3d_o%rank,&
                       ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)
    !
    call MPI_Comm_size(topo_3d_o%comm,&
                       topo_3d_o%nprocs,&
                       ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    deallocate(topo_3d_rank_table)
    deallocate(topo_6d_rank_grouped_by_3d_rank)

    call MPI_Group_free(group_3d, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    call MPI_Group_free(group_6d, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)

    sll_f_new_cartesian_topology_3d_orthogonal => topo_3d_o
  end function sll_f_new_cartesian_topology_3d_orthogonal


  function sll_f_new_cartesian_domain_decomposition_6d(topology, grid_size, halo_width)
    type(sll_t_decomposition_6d), pointer :: sll_f_new_cartesian_domain_decomposition_6d
    type(sll_t_cartesian_topology_6d), pointer, intent(in) :: topology
    integer, parameter :: nd=6
    sll_int32, intent(in) :: grid_size(nd)
    sll_int32, intent(in) :: halo_width(nd)
    sll_int32 :: i, ierr
    sll_int32 :: lp, l0, l1

    ! --- initial checks
    do i=1,nd
       SLL_ASSERT( halo_width(i) >= 0 )
    end do
    do i=1,nd
       ! for the moment, the global grid must be divisible evenly among the MPI processes
       SLL_ASSERT( mod(grid_size(i), topology%procs(i)) == 0 )
    end do

    ! --- copy and calculate all the necessary index values ---
    SLL_ALLOCATE(sll_f_new_cartesian_domain_decomposition_6d, ierr)

    sll_f_new_cartesian_domain_decomposition_6d%global =  grid_size

    ! --- loop over dimensions and compute index values for each dimension
    do i=1,nd
       ! compute the local number of grid points
       lp = grid_size(i) / topology%procs(i)
       ! compute the lower local index bound (the coords array starts at zero)
       l0 = 1 + topology%coords(i) * lp
       ! compute the upper local index bound (the coords array starts at zero)
       l1 = (topology%coords(i) + 1) * lp

       SLL_ASSERT( lp/2 >= halo_width(i) )

       sll_f_new_cartesian_domain_decomposition_6d%local%mn(i) = l0
       sll_f_new_cartesian_domain_decomposition_6d%local%mx(i) = l1
       sll_f_new_cartesian_domain_decomposition_6d%local%hw(i) = halo_width(i)
       sll_f_new_cartesian_domain_decomposition_6d%local%lo(i) = l0 - halo_width(i)
       sll_f_new_cartesian_domain_decomposition_6d%local%hi(i) = l1 + halo_width(i)

       sll_f_new_cartesian_domain_decomposition_6d%local%tx_lolo(i) = l0
       sll_f_new_cartesian_domain_decomposition_6d%local%tx_lohi(i) = l0 + (halo_width(i) - 1)
       sll_f_new_cartesian_domain_decomposition_6d%local%tx_hilo(i) = l1 - (halo_width(i) - 1)
       sll_f_new_cartesian_domain_decomposition_6d%local%tx_hihi(i) = l1
       sll_f_new_cartesian_domain_decomposition_6d%local%rx_lolo(i) = l0 - halo_width(i)
       sll_f_new_cartesian_domain_decomposition_6d%local%rx_lohi(i) = l0 - 1
       sll_f_new_cartesian_domain_decomposition_6d%local%rx_hilo(i) = l1 + 1
       sll_f_new_cartesian_domain_decomposition_6d%local%rx_hihi(i) = l1 + halo_width(i)
    end do

    ! compute net array width
    sll_f_new_cartesian_domain_decomposition_6d%local%nw = 1 + &
       sll_f_new_cartesian_domain_decomposition_6d%local%mx - sll_f_new_cartesian_domain_decomposition_6d%local%mn

    ! compute gross array width
    sll_f_new_cartesian_domain_decomposition_6d%local%gw = 1 + &
       sll_f_new_cartesian_domain_decomposition_6d%local%hi - sll_f_new_cartesian_domain_decomposition_6d%local%lo

  end function sll_f_new_cartesian_domain_decomposition_6d


  !> @brief 6D Cartesian domain decomposition destructor
  subroutine sll_s_deallocate_cartesian_domain_decomposition_6d()
  end subroutine sll_s_deallocate_cartesian_domain_decomposition_6d


  ! --- allocator for redesigned ("slim") decomposition with dynamic halos
  function sll_f_new_cartesian_domain_decomposition_slim_6d(topology, grid_size)
    type(sll_t_decomposition_slim_6d), pointer :: sll_f_new_cartesian_domain_decomposition_slim_6d
    type(sll_t_cartesian_topology_6d), pointer, intent(in) :: topology
    integer, parameter :: nd=6
    sll_int32, intent(in) :: grid_size(nd)
    sll_int32 :: i, ierr
    sll_int32 :: lp, l0, l1

    ! --- initial checks
    do i=1,nd
       ! for the moment, the global grid must be divisible among the MPI processes
       SLL_ASSERT( mod(grid_size(i), topology%procs(i)) == 0 )
    end do

    ! --- copy and calculate all the necessary index values
    SLL_ALLOCATE(sll_f_new_cartesian_domain_decomposition_slim_6d, ierr)
    sll_f_new_cartesian_domain_decomposition_slim_6d%global = grid_size
    ! loop over dimensions and compute index values for each dimension
    do i=1,nd
       ! compute the local number of grid points
       lp = grid_size(i) / topology%procs(i)
       ! compute the lower local index bound (the coords array starts at zero)
       l0 = 1 + topology%coords(i) * lp
       ! compute the upper local index bound (the coords array starts at zero)
       l1 = (topology%coords(i) + 1) * lp

       sll_f_new_cartesian_domain_decomposition_slim_6d%local%mn(i) = l0
       sll_f_new_cartesian_domain_decomposition_slim_6d%local%mx(i) = l1
       sll_f_new_cartesian_domain_decomposition_slim_6d%local%nw(i) = lp
    end do

    sll_f_new_cartesian_domain_decomposition_slim_6d%local%id = -1

! #ifdef USE_FMEMPOOL
!     nullify(sll_f_new_cartesian_domain_decomposition_slim_6d%local%halo_left%buf)
!     nullify(sll_f_new_cartesian_domain_decomposition_slim_6d%local%halo_right%buf)
! #endif
  end function sll_f_new_cartesian_domain_decomposition_slim_6d


  !> @brief 6D Cartesian slim domain decomposition destructor
  subroutine sll_s_deallocate_cartesian_domain_decomposition_slim_6d()
  end subroutine sll_s_deallocate_cartesian_domain_decomposition_slim_6d


  function sll_f_new_cartesian_domain_decomposition_slim_3d(topology, grid_size)
    type(sll_t_decomposition_slim_3d), pointer :: sll_f_new_cartesian_domain_decomposition_slim_3d
    type(sll_t_cartesian_topology_3d), pointer, intent(in) :: topology
    integer, parameter :: nd=3
    sll_int32, intent(in) :: grid_size(nd)
    sll_int32 :: i, ierr
    sll_int32 :: lp, l0, l1

    ! --- initial checks
    do i=1,nd
       ! for the moment, the global grid must be divisible among the MPI processes
       SLL_ASSERT( mod(grid_size(i), topology%procs(i)) == 0 )
    end do

    ! --- copy and calculate all the necessary index values
    SLL_ALLOCATE(sll_f_new_cartesian_domain_decomposition_slim_3d, ierr)
    sll_f_new_cartesian_domain_decomposition_slim_3d%global = grid_size
    ! loop over dimensions and compute index values for each dimension
    do i=1,nd
       ! compute the local number of grid points
       lp = grid_size(i) / topology%procs(i)
       ! compute the lower local index bound (the coords array starts at zero)
       l0 = 1 + topology%coords(i) * lp
       ! compute the upper local index bound (the coords array starts at zero)
       l1 = (topology%coords(i) + 1) * lp

       sll_f_new_cartesian_domain_decomposition_slim_3d%local%mn(i) = l0
       sll_f_new_cartesian_domain_decomposition_slim_3d%local%mx(i) = l1
       sll_f_new_cartesian_domain_decomposition_slim_3d%local%nw(i) = lp
    end do

    sll_f_new_cartesian_domain_decomposition_slim_3d%local%id = -1
  end function sll_f_new_cartesian_domain_decomposition_slim_3d


  !> @brief 3D Cartesian slim domain decomposition destructor
  subroutine sll_s_deallocate_cartesian_domain_decomposition_slim_3d()
  end subroutine sll_s_deallocate_cartesian_domain_decomposition_slim_3d


  ! --- allocator for redesigned ("slim") decomposition with dynamic halos
  function sll_f_new_cartesian_cell_domain_decomposition_slim_6d(topology, n_cells, degree)
    type(sll_t_decomposition_slim_6d), pointer :: sll_f_new_cartesian_cell_domain_decomposition_slim_6d
    type(sll_t_cartesian_topology_6d), pointer, intent(in) :: topology
    integer, parameter :: nd=6
    sll_int32, intent(in) :: n_cells(nd)
    sll_int32, intent(in) :: degree(nd)

    sll_int32 :: i, ierr
    sll_int32 :: lp, l0, l1

    ! --- initial checks
    do i=1,nd
       ! for the moment, the global grid must be divisible among the MPI processes
       SLL_ASSERT( mod(n_cells(i), topology%procs(i)) == 0 )
    end do

    ! --- copy and calculate all the necessary index values
    SLL_ALLOCATE(sll_f_new_cartesian_cell_domain_decomposition_slim_6d, ierr)

    sll_f_new_cartesian_cell_domain_decomposition_slim_6d%global = n_cells * degree
    ! GRID indices: loop over dimensions and compute index values for each dimension
    do i=1,nd
       ! compute the local number of grid points
       lp = (n_cells(i) / topology%procs(i)) * degree(i)
       ! compute the lower local index bound (the coords array starts at zero)
       l0 = 1 + topology%coords(i) * lp
       ! compute the upper local index bound (the coords array starts at zero)
       l1 = (topology%coords(i) + 1) * lp

       sll_f_new_cartesian_cell_domain_decomposition_slim_6d%local%mn(i) = l0
       sll_f_new_cartesian_cell_domain_decomposition_slim_6d%local%mx(i) = l1
       sll_f_new_cartesian_cell_domain_decomposition_slim_6d%local%nw(i) = lp
    end do

    sll_f_new_cartesian_cell_domain_decomposition_slim_6d%n_cells = n_cells
    ! CELL indices: loop over dimensions and compute index values for each dimension
    do i=1,nd
       ! compute the local number of cells
       lp = n_cells(i) / topology%procs(i)
       ! compute the lower local index bound (the coords array starts at zero)
       l0 = 1 + topology%coords(i) * lp
       ! compute the upper local index bound (the coords array starts at zero)
       l1 = (topology%coords(i) + 1) * lp

       sll_f_new_cartesian_cell_domain_decomposition_slim_6d%local%mn_cell(i) = l0
       sll_f_new_cartesian_cell_domain_decomposition_slim_6d%local%mx_cell(i) = l1
       sll_f_new_cartesian_cell_domain_decomposition_slim_6d%local%n_cells(i) = lp
    end do

    sll_f_new_cartesian_cell_domain_decomposition_slim_6d%local%id = -1

! #ifdef USE_FMEMPOOL
!     nullify(sll_f_new_cartesian_cell_domain_decomposition_slim_6d%local%halo_left%buf)
!     nullify(sll_f_new_cartesian_cell_domain_decomposition_slim_6d%local%halo_right%buf)
! #endif
  end function sll_f_new_cartesian_cell_domain_decomposition_slim_6d

  !> @brief 6D Cartesian slim domain decomposition destructor
  subroutine sll_s_deallocate_cartesian_cell_domain_decomposition_slim_6d()
!#ifdef USE_FMEMPOOL
!    call mp_cleanup()
!#endif
  end subroutine sll_s_deallocate_cartesian_cell_domain_decomposition_slim_6d


  function sll_f_new_cartesian_domain_decomposition_3d(topology, grid_size, halo_width)
    type(sll_t_decomposition_3d), pointer :: sll_f_new_cartesian_domain_decomposition_3d
    type(sll_t_cartesian_topology_3d), pointer, intent(in) :: topology
    integer, parameter :: nd=3
    sll_int32, intent(in) :: grid_size(nd)
    sll_int32, intent(in) :: halo_width(nd)
    sll_int32 :: i, ierr
    sll_int32 :: lp, l0, l1

    ! --- initial checks
    do i=1,nd
       SLL_ASSERT( halo_width(i) >= 0 )
    end do
    do i=1,nd
       ! for the moment, the global grid must be divisible evenly among the MPI processes
       SLL_ASSERT( mod(grid_size(i), topology%procs(i)) == 0 )
    end do

    ! --- copy and calculate all the necessary index values ---
    SLL_ALLOCATE(sll_f_new_cartesian_domain_decomposition_3d, ierr)

    sll_f_new_cartesian_domain_decomposition_3d%global =  grid_size

    ! --- loop over dimensions and compute index values for each dimension
    do i=1,nd
       ! compute the local number of grid points
       lp = grid_size(i) / topology%procs(i)
       ! compute the lower local index bound (the coords array starts at zero)
       l0 = 1 + topology%coords(i) * lp
       ! compute the upper local index bound (the coords array starts at zero)
       l1 = (topology%coords(i) + 1) * lp

       SLL_ASSERT( lp/2 >= halo_width(i) )

       sll_f_new_cartesian_domain_decomposition_3d%local%mn(i) = l0
       sll_f_new_cartesian_domain_decomposition_3d%local%mx(i) = l1
       sll_f_new_cartesian_domain_decomposition_3d%local%hw(i) = halo_width(i)
       sll_f_new_cartesian_domain_decomposition_3d%local%lo(i) = l0 - halo_width(i)
       sll_f_new_cartesian_domain_decomposition_3d%local%hi(i) = l1 + halo_width(i)

       sll_f_new_cartesian_domain_decomposition_3d%local%tx_lolo(i) = l0
       sll_f_new_cartesian_domain_decomposition_3d%local%tx_lohi(i) = l0 + (halo_width(i) - 1)
       sll_f_new_cartesian_domain_decomposition_3d%local%tx_hilo(i) = l1 - (halo_width(i) - 1)
       sll_f_new_cartesian_domain_decomposition_3d%local%tx_hihi(i) = l1
       sll_f_new_cartesian_domain_decomposition_3d%local%rx_lolo(i) = l0 - halo_width(i)
       sll_f_new_cartesian_domain_decomposition_3d%local%rx_lohi(i) = l0 - 1
       sll_f_new_cartesian_domain_decomposition_3d%local%rx_hilo(i) = l1 + 1
       sll_f_new_cartesian_domain_decomposition_3d%local%rx_hihi(i) = l1 + halo_width(i)
    end do

    ! compute net array width
    sll_f_new_cartesian_domain_decomposition_3d%local%nw = 1 + &
       sll_f_new_cartesian_domain_decomposition_3d%local%mx - sll_f_new_cartesian_domain_decomposition_3d%local%mn

    ! compute gross array width
    sll_f_new_cartesian_domain_decomposition_3d%local%gw = 1 + &
       sll_f_new_cartesian_domain_decomposition_3d%local%hi - sll_f_new_cartesian_domain_decomposition_3d%local%lo
  end function sll_f_new_cartesian_domain_decomposition_3d


  !> @brief 3D Cartesian domain decomposition destructor
  subroutine sll_s_deallocate_cartesian_domain_decomposition_3d()
  end subroutine sll_s_deallocate_cartesian_domain_decomposition_3d


  subroutine sll_s_copy_array_to_buffer_6d_real64(arr, arr_lo, arr_hi, buf, ranges, n_threads)
    sll_int32, dimension(6), intent(in) :: arr_lo
    sll_int32, dimension(6), intent(in) :: arr_hi
    sll_real64, dimension(arr_lo(1):arr_hi(1), arr_lo(2):arr_hi(2), arr_lo(3):arr_hi(3),  &
                          arr_lo(4):arr_hi(4), arr_lo(5):arr_hi(5), arr_lo(6):arr_hi(6)), &
                              intent(in) :: arr
    HALO_DTYPE, dimension(:), intent(out) :: buf
    sll_int32, dimension(6,2), intent(in) :: ranges
    sll_int32, intent(in), optional :: n_threads
    sll_int32 :: i,j,k,l,m,n
    sll_int32 :: ii,ij,ik,il,im,in  ! original indices, computed from zero-based ones
    sll_int32 :: wi,wj,wk,wl,wm,wn  ! widths of the array, for zero-based indexing
    sll_int32 :: oj,ok,ol,om,on  ! offsets of the buffer, for zero-based indexing
    sll_int32 :: idx, n_omp_threads

#ifdef _OPENMP
    if (present(n_threads)) then
      n_omp_threads = n_threads
    else
      n_omp_threads = omp_get_max_threads()
    endif
#else
    n_omp_threads = 1
#endif

    wi = ranges(1,2) - ranges(1,1) + 1
    wj = ranges(2,2) - ranges(2,1) + 1
    wk = ranges(3,2) - ranges(3,1) + 1
    wl = ranges(4,2) - ranges(4,1) + 1
    wm = ranges(5,2) - ranges(5,1) + 1
    wn = ranges(6,2) - ranges(6,1) + 1
!$omp parallel num_threads(n_omp_threads) default(shared) private(i,j,k,l,m,n, ii,ij,ik,il,im,in, oj,ok,ol,om,on, idx)
!$omp do OMP_SCHEDULE OMP_COLLAPSE
    do n=0,wn-1
       do m=0,wm-1
          ! --- enable OMP collapse
          in = n + ranges(6,1)
          on = n * wi * wj * wk * wl * wm
          ! ---
          im = m + ranges(5,1)
          om = m * wi * wj * wk * wl
          do l=0,wl-1
             il = l + ranges(4,1)
             ol = l * wi * wj * wk
             do k=0,wk-1
                ik = k + ranges(3,1)
                ok = k * wi * wj
                do j=0,wj-1
                   ij = j + ranges(2,1)
                   oj = j * wi
                   do i=0,wi-1
                      ii = i + ranges(1,1)
                      ! linear index calculation
                      idx = 1 + i + oj + ok + ol + om + on
                      buf(idx) = arr(ii,ij,ik,il,im,in)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel
! --- previous version, a little slower, non-parallel, kept for documentation purposes
!      idx=1
!      do n=ranges(6,1),ranges(6,2)
!         do m=ranges(5,1),ranges(5,2)
!            do l=ranges(4,1),ranges(4,2)
!               do k=ranges(3,1),ranges(3,2)
!                  do j=ranges(2,1),ranges(2,2)
!                     do i=ranges(1,1),ranges(1,2)
!                        buf(idx) = arr(i,j,k,l,m,n)
!                        idx = idx + 1
!                     enddo
!                  enddo
!               enddo
!            enddo
!         enddo
!      enddo
  end subroutine sll_s_copy_array_to_buffer_6d_real64


  subroutine sll_s_copy_array_to_buffer_3d_real64(arr, arr_lo, arr_hi, buf, ranges, n_threads)
    sll_int32, dimension(3), intent(in) :: arr_lo
    sll_int32, dimension(3), intent(in) :: arr_hi
    sll_real64, dimension(arr_lo(1):arr_hi(1), arr_lo(2):arr_hi(2), arr_lo(3):arr_hi(3)), intent(in) :: arr
    HALO_DTYPE, dimension(:), intent(out) :: buf
    sll_int32, dimension(3,2), intent(in) :: ranges
    sll_int32, intent(in), optional :: n_threads

    sll_int32 :: i,j,k
    sll_int32 :: ii,ij,ik  ! original indices, computed from zero-based ones
    sll_int32 :: wi,wj,wk  ! widths of the array, for zero-based indexing
    sll_int32 :: oj,ok  ! offsets of the buffer, for zero-based indexing
    sll_int32 :: idx, n_omp_threads

#ifdef _OPENMP
    if (present(n_threads)) then
      n_omp_threads = n_threads
    else
      n_omp_threads = omp_get_max_threads()
    endif
#else
    n_omp_threads = 1
#endif

    wi = ranges(1,2) - ranges(1,1) + 1
    wj = ranges(2,2) - ranges(2,1) + 1
    wk = ranges(3,2) - ranges(3,1) + 1
!$omp parallel num_threads(n_omp_threads) default(shared) private(i,j,k,ii,ij,ik,oj,ok,idx)
!$omp do OMP_SCHEDULE OMP_COLLAPSE
    do k=0,wk-1
       do j=0,wj-1
          ik = k + ranges(3,1)
          ok = k * wi * wj
          ij = j + ranges(2,1)
          oj = j * wi
          do i=0,wi-1
             ii = i + ranges(1,1)
             idx = 1 + i + oj + ok
             buf(idx) = arr(ii,ij,ik)
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel
  end subroutine sll_s_copy_array_to_buffer_3d_real64



  subroutine copy_buffer_to_array_6d_real64(buf, arr, arr_lo, arr_hi, ranges, n_threads)
    HALO_DTYPE, dimension(:), intent(in) :: buf
    sll_int32, dimension(6), intent(in) :: arr_lo
    sll_int32, dimension(6), intent(in) :: arr_hi
    sll_real64, dimension(arr_lo(1):arr_hi(1), arr_lo(2):arr_hi(2), arr_lo(3):arr_hi(3), &
                          arr_lo(4):arr_hi(4), arr_lo(5):arr_hi(5), arr_lo(6):arr_hi(6)),&
                              intent(out) :: arr
    sll_int32, dimension(6,2), intent(in) :: ranges
    sll_int32, intent(in), optional :: n_threads
    sll_int32 :: i,j,k,l,m,n
    sll_int32 :: ii,ij,ik,il,im,in  ! original indices, computed from zero-based ones
    sll_int32 :: wi,wj,wk,wl,wm,wn  ! widths of the array, for zero-based indexing
    sll_int32 :: oj,ok,ol,om,on  ! offsets of the buffer, for zero-based indexing
    sll_int32 :: idx, n_omp_threads

#ifdef _OPENMP
    if (present(n_threads)) then
      n_omp_threads = n_threads
    else
      n_omp_threads = omp_get_max_threads()
    endif
#else
    n_omp_threads = 1
#endif

    wi = ranges(1,2) - ranges(1,1) + 1
    wj = ranges(2,2) - ranges(2,1) + 1
    wk = ranges(3,2) - ranges(3,1) + 1
    wl = ranges(4,2) - ranges(4,1) + 1
    wm = ranges(5,2) - ranges(5,1) + 1
    wn = ranges(6,2) - ranges(6,1) + 1
!$omp parallel num_threads(n_omp_threads) default(shared) private(i,j,k,l,m,n, ii,ij,ik,il,im,in, oj,ok,ol,om,on, idx)
!$omp do OMP_SCHEDULE OMP_COLLAPSE
    do n=0,wn-1
       do m=0,wm-1
          ! --- enable OMP collapse
          in = n + ranges(6,1)
          on = n * wi * wj * wk * wl * wm
          ! ---
          im = m + ranges(5,1)
          om = m * wi * wj * wk * wl
          do l=0,wl-1
             il = l + ranges(4,1)
             ol = l * wi * wj * wk
             do k=0,wk-1
                ik = k + ranges(3,1)
                ok = k * wi * wj
                do j=0,wj-1
                   ij = j + ranges(2,1)
                   oj = j * wi
                   do i=0,wi-1
                      ii = i + ranges(1,1)
                      ! linear index calculation
                      idx = 1 + i + oj + ok + ol + om + on
                      arr(ii,ij,ik,il,im,in) = buf(idx)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
!$omp end do
!$omp end parallel
! --- previous version, a little slower, non-parallel, kept for documentation purposes
!      idx=1
!      do n=ranges(6,1),ranges(6,2)
!         do m=ranges(5,1),ranges(5,2)
!            do l=ranges(4,1),ranges(4,2)
!               do k=ranges(3,1),ranges(3,2)
!                  do j=ranges(2,1),ranges(2,2)
!                     do i=ranges(1,1),ranges(1,2)
!                        arr(i,j,k,l,m,n) = buf(idx)
!                        idx = idx + 1
!                     enddo
!                  enddo
!               enddo
!            enddo
!         enddo
!      enddo
  end subroutine copy_buffer_to_array_6d_real64


  subroutine sll_f_apply_halo_exchange_6d_real64(topo, decomp, arr, dim_mask_in)
    integer, parameter :: nd = 6  ! we handle 6 dimensions
    ! --- arguments
    type(sll_t_cartesian_topology_6d), intent(in) :: topo
    type(sll_t_decomposition_6d), intent(in) :: decomp
    sll_real64, dimension(:,:,:,:,:,:), intent(inout) :: arr(decomp%local%lo(1):decomp%local%hi(1), &
                                                             decomp%local%lo(2):decomp%local%hi(2), &
                                                             decomp%local%lo(3):decomp%local%hi(3), &
                                                             decomp%local%lo(4):decomp%local%hi(4), &
                                                             decomp%local%lo(5):decomp%local%hi(5), &
                                                             decomp%local%lo(6):decomp%local%hi(6))
    logical, dimension(:), intent(in), optional :: dim_mask_in(nd)
    ! --- MPI communication-related variables and buffers
    sll_int64, save :: bufsize = 0
    sll_int64 :: nxc  ! total number of elements to be exchanged
    sll_int64 :: off  ! offset where to start from
    sll_int64 :: rem  ! remaining elements yet to be exchanged
    sll_int32 :: nel  ! number of elements to be exchanged at a single MPI call
    integer, dimension(:,:) :: r_rx(nd,2)  ! index ranges for the direct copy operations
    integer, dimension(:,:) :: r_tx(nd,2)
#ifdef USE_HALO_REAL32
    sll_int32, parameter :: nxc_max = 128000000
    integer, parameter :: mpi_precision = MPI_REAL
    integer, parameter :: word_size = 4
#else
    sll_int32, parameter :: nxc_max = 64000000
    integer, parameter :: mpi_precision = MPI_DOUBLE_PRECISION
    integer, parameter :: word_size = 8
#endif
    HALO_DTYPE, dimension(:), allocatable, target, save :: sendbuf, recvbuf
    integer :: mpi_tag
    ! ---
    integer :: id, jd
    integer :: ierr

!    integer :: nbytes_in
!    integer :: off_elements_thread, off_bytes_thread, nel_thread, omp_rank, omp_size
!    sll_int32 :: n_packets

    logical, save :: first_call = .true.
    logical, save :: sll_use_mpi_sendrecv = .true.
    integer :: i,j,k,l,m,n
    ! --- DEBUG
    integer, save :: dump_buffer_invocation_count = 0
    integer, save :: dump_f6d_invocation_count = 0
    integer, save :: dump_dd_info_invocation_count = 0
    ! --- dimension selection mask
    logical, dimension(:) :: dim_mask(nd)

    dim_mask = .true.  ! by default, we exchange all 6 dimensions
    if (present(dim_mask_in))  dim_mask = dim_mask_in

    mpi_tag = 0

    ! --- query environment variables
    if (first_call) then
#ifdef USE_HALO_REAL32
      if (topo%rank == 0) then
        write(*,*) "sll_m_decomposition::apply_halo_exchange() uses single precision messages"
      endif
#endif
      sll_use_mpi_sendrecv = sll_f_query_environment("SLL_USE_MPI_SENDRECV", .true.)
      if (topo%rank == 0) then
        if (sll_use_mpi_sendrecv) then
          write(*,*) "sll_m_decomposition::apply_halo_exchange() uses MPI_SENDRECV()"
        else
          write(*,*) "sll_m_decomposition::apply_halo_exchange() uses MPI_SEND() and MPI_RECV()"
        endif
      end if
      first_call = .false.
    endif

    ! --- loop over dimensions and exchange data between neighbors
    do id=1,nd
      if (.not. dim_mask(id)) cycle

      if (topo%procs(id) == 1) then
        ! --- we copy the ghost cells directly, assuming periodic BCs
        ! (1) copy to the left
        r_tx(:,1)  = decomp%local%mn(:)
        r_tx(:,2)  = decomp%local%mx(:)
        r_tx(id,1) = decomp%local%tx_lolo(id)
        r_tx(id,2) = decomp%local%tx_lohi(id)
        r_rx(:,1)  = decomp%local%mn(:)
        r_rx(:,2)  = decomp%local%mx(:)
        r_rx(id,1) = decomp%local%rx_hilo(id)
        r_rx(id,2) = decomp%local%rx_hihi(id)
!$omp parallel default(shared) private(i,j,k,l,m,n)
!$omp do OMP_SCHEDULE OMP_COLLAPSE
        do n = 0, r_rx(6,2)-r_rx(6,1)
          do m = 0, r_rx(5,2)-r_rx(5,1)
            do l = 0, r_rx(4,2)-r_rx(4,1)
              do k = 0, r_rx(3,2)-r_rx(3,1)
                do j = 0, r_rx(2,2)-r_rx(2,1)
                  do i = 0, r_rx(1,2)-r_rx(1,1)
                    arr(r_rx(1,1)+i, r_rx(2,1)+j, r_rx(3,1)+k, r_rx(4,1)+l, r_rx(5,1)+m, r_rx(6,1)+n) = &
                    arr(r_tx(1,1)+i, r_tx(2,1)+j, r_tx(3,1)+k, r_tx(4,1)+l, r_tx(5,1)+m, r_tx(6,1)+n)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
!$omp end do
!$omp end parallel
        ! (2) copy to the right
        r_tx(:,1)  = decomp%local%mn(:)
        r_tx(:,2)  = decomp%local%mx(:)
        r_tx(id,1) = decomp%local%tx_hilo(id)
        r_tx(id,2) = decomp%local%tx_hihi(id)
        r_rx(:,1)  = decomp%local%mn(:)
        r_rx(:,2)  = decomp%local%mx(:)
        r_rx(id,1) = decomp%local%rx_lolo(id)
        r_rx(id,2) = decomp%local%rx_lohi(id)
!$omp parallel default(shared) private(i,j,k,l,m,n)
!$omp do OMP_SCHEDULE OMP_COLLAPSE
        do n = 0, r_rx(6,2)-r_rx(6,1)
          do m = 0, r_rx(5,2)-r_rx(5,1)
            do l = 0, r_rx(4,2)-r_rx(4,1)
              do k = 0, r_rx(3,2)-r_rx(3,1)
                do j = 0, r_rx(2,2)-r_rx(2,1)
                  do i = 0, r_rx(1,2)-r_rx(1,1)
                    arr(r_rx(1,1)+i, r_rx(2,1)+j, r_rx(3,1)+k, r_rx(4,1)+l, r_rx(5,1)+m, r_rx(6,1)+n) = &
                    arr(r_tx(1,1)+i, r_tx(2,1)+j, r_tx(3,1)+k, r_tx(4,1)+l, r_tx(5,1)+m, r_tx(6,1)+n)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
!$omp end do
!$omp end parallel
      else
        ! --- calculate the number of items to be exchanged
        nxc = int(decomp%local%hw(id), i64)
        do jd=1,nd
          if (jd == id) then
            cycle
          else
            nxc = nxc * decomp%local%nw(jd)
          endif
        end do

        if (nxc > bufsize) then
          if (allocated(sendbuf)) deallocate(sendbuf)
          if (allocated(recvbuf)) deallocate(recvbuf)
          allocate(sendbuf(nxc))
          allocate(recvbuf(nxc))
          bufsize = nxc
        end if

        ! --- (1) copy halo cells to the left neighbor
        r_tx(:,1) = decomp%local%mn(:)
        r_tx(:,2) = decomp%local%mx(:)
        r_tx(id,1) = decomp%local%tx_lolo(id)
        r_tx(id,2) = decomp%local%tx_lohi(id)
        call sll_s_copy_array_to_buffer_6d_real64(arr, decomp%local%lo, decomp%local%hi, sendbuf, r_tx)

        ! --- split MPI communication into pieces smaller than nxc_max
        rem = nxc
        off = 1
        do while (rem > 0)
          if (rem > nxc_max) then
            nel = nxc_max
          else
            nel = int(rem, i32)
          endif
          if (sll_use_mpi_sendrecv) then
            call MPI_Sendrecv(sendbuf(off), nel, mpi_precision, topo%neighbors(2*id-1), mpi_tag,&
                              recvbuf(off), nel, mpi_precision, topo%neighbors(2*id),   mpi_tag,&
                              topo%comm, MPI_STATUS_IGNORE, ierr)
            SLL_ASSERT(ierr == MPI_SUCCESS)
          else  ! sll_use_mpi_sendrecv
            ! --- odd processor coordinate numbers send first to the left
            if (mod(topo%coords(id), 2) > 0) then
              call MPI_Send(sendbuf(off), nel, mpi_precision, topo%neighbors(2*id-1), mpi_tag, topo%comm, ierr)
              SLL_ASSERT(ierr == MPI_SUCCESS)
              call MPI_Recv(recvbuf(off), nel, mpi_precision, topo%neighbors(2*id), mpi_tag, topo%comm, MPI_STATUS_IGNORE, ierr)
              SLL_ASSERT(ierr == MPI_SUCCESS)
            else
              ! --- even processor coordinate numbers receive first from the right
              call MPI_Recv(recvbuf(off), nel, mpi_precision, topo%neighbors(2*id), mpi_tag, topo%comm, MPI_STATUS_IGNORE, ierr)
              SLL_ASSERT(ierr == MPI_SUCCESS)
              call MPI_Send(sendbuf(off), nel, mpi_precision, topo%neighbors(2*id-1), mpi_tag, topo%comm, ierr)
              SLL_ASSERT(ierr == MPI_SUCCESS)
            endif
          endif  ! sll_use_mpi_sendrecv
          off = off + nel
          rem = rem - nel
        enddo  ! while(rem>0)

        r_rx(:,1) = decomp%local%mn(:)
        r_rx(:,2) = decomp%local%mx(:)
        r_rx(id,1) = decomp%local%rx_hilo(id)
        r_rx(id,2) = decomp%local%rx_hihi(id)
        call copy_buffer_to_array_6d_real64(recvbuf, arr, decomp%local%lo, decomp%local%hi, r_rx)


        ! --- (2) copy halo cells to the right neighbor
        r_tx(:,1) = decomp%local%mn(:)
        r_tx(:,2) = decomp%local%mx(:)
        r_tx(id,1) = decomp%local%tx_hilo(id)
        r_tx(id,2) = decomp%local%tx_hihi(id)
        call sll_s_copy_array_to_buffer_6d_real64(arr, decomp%local%lo, decomp%local%hi, sendbuf, r_tx)
        rem = nxc
        off = 1
        do while (rem > 0)
          if (rem > nxc_max) then
            nel = nxc_max
          else
            nel = int(rem, i32)
          endif
          if (sll_use_mpi_sendrecv) then
            call MPI_Sendrecv(sendbuf(off), nel, mpi_precision, topo%neighbors(2*id), 1,&
                              recvbuf(off), nel, mpi_precision, topo%neighbors(2*id-1), 1,&
                              topo%comm, MPI_STATUS_IGNORE, ierr)
            SLL_ASSERT(ierr == MPI_SUCCESS)
          else  ! sll_use_mpi_sendrecv
            ! --- odd processor coordinate numbers send first to the right
            if (mod(topo%coords(id), 2) > 0) then
              call MPI_Send(sendbuf(off), nel, mpi_precision, topo%neighbors(2*id), &
                            mpi_tag, topo%comm, ierr)
              SLL_ASSERT(ierr == MPI_SUCCESS)
              call MPI_Recv(recvbuf(off), nel, mpi_precision, topo%neighbors(2*id-1), &
                            mpi_tag, topo%comm, MPI_STATUS_IGNORE, ierr)
              SLL_ASSERT(ierr == MPI_SUCCESS)
            else
            ! --- even processor coordinate numbers receive first from the left
              call MPI_Recv(recvbuf(off), nel, mpi_precision, topo%neighbors(2*id-1), &
                            mpi_tag, topo%comm, MPI_STATUS_IGNORE, ierr)
              SLL_ASSERT(ierr == MPI_SUCCESS)
              call MPI_Send(sendbuf(off), nel, mpi_precision, topo%neighbors(2*id), &
                            mpi_tag, topo%comm, ierr)
              SLL_ASSERT(ierr == MPI_SUCCESS)
            endif
          endif  ! sll_use_mpi_sendrecv
          off = off + nel
          rem = rem - nel
        enddo

        r_rx(:,1) = decomp%local%mn(:)
        r_rx(:,2) = decomp%local%mx(:)
        r_rx(id,1) = decomp%local%rx_lolo(id)
        r_rx(id,2) = decomp%local%rx_lohi(id)
        call copy_buffer_to_array_6d_real64(recvbuf, arr, decomp%local%lo, decomp%local%hi, r_rx)

      endif  ! if branch (topo%nprocs > 1)
    enddo  ! id (loop over dimensions)

  end subroutine sll_f_apply_halo_exchange_6d_real64


!    call MPI_Sendrecv(sendbuf,                     nel, mpi_precision, topo%neighbors(2*id-1), mpi_tag,&
!                      decomp%local%halo_right%buf, nel, mpi_precision, topo%neighbors(2*id),   mpi_tag,&
!                      topo%comm, MPI_STATUS_IGNORE, ierr)
  subroutine mpi_sendrecv_compressed_6d_real64(sendbuf, recvbuf, nel, rank_send,&
                                          rank_recv, mpi_comm, verbose, mpi_tag)
    use iso_c_binding
    HALO_DTYPE, pointer :: sendbuf(:)
    HALO_DTYPE, pointer :: recvbuf(:,:,:,:,:,:), recv_view_1d(:)
    sll_int32 :: nel, rank_send, rank_recv, mpi_comm
    logical, intent(in), optional :: verbose
    integer, intent(in), optional :: mpi_tag
    logical :: verb
    integer :: tag

    if (present(verbose)) then
      verb = verbose
    else
      verb = .false.
    endif

    if (present(mpi_tag)) then
      tag = mpi_tag
    else
      tag = 0
    endif

    call c_f_pointer(c_loc(recvbuf), recv_view_1d, [nel])
    call mpi_sendrecv_compressed(sendbuf, recv_view_1d, nel, rank_send, rank_recv, mpi_comm, verb, tag)
    nullify(recv_view_1d)
  end subroutine mpi_sendrecv_compressed_6d_real64

  subroutine mpi_sendrecv_compressed(sendbuf, recvbuf, nel, rank_send, rank_recv,&
                                     mpi_comm, verbose, mpi_tag)
    HALO_DTYPE, pointer :: sendbuf(:)
    HALO_DTYPE, pointer :: recvbuf(:)
    sll_int32 :: nel, rank_send, rank_recv, mpi_comm
    logical, intent(in), optional :: verbose
    integer, intent(in), optional :: mpi_tag
    sll_int32 :: ierr
#ifdef USE_HALO_REAL32
    integer, parameter :: mpi_precision = MPI_REAL
#else
    integer, parameter :: mpi_precision = MPI_DOUBLE_PRECISION
#endif
    logical :: do_compress
    type(sll_t_compressed_buffer) :: comp_send, comp_recv  ! data structures containing compressed data and offsets
    integer ::  omp_size
    logical :: verb
    integer :: tag

    if (present(verbose)) then
      verb = verbose
    else
      verb = .false.
    endif

    if (present(mpi_tag)) then
      tag = mpi_tag
    else
      tag = 0
    endif

#ifdef _OPENMP
    omp_size = omp_get_max_threads()
#else
    omp_size = 1
#endif

#ifdef USE_HALO_REAL32
      do_compress = .false.
#else
    if (modulo(nel, zfp_blocksize) /= 0) then
      write(*,*) "mpi_sendrecv_compressed() : disabled due to blocksize mismatch"
      do_compress = .false.
!    elseif (modulo(nel/zfp_blocksize, omp_size) /= 0)
!      write(*,*) "mpi_sendrecv_compressed() : disabled due to n_slices vs. omp_size mismatch"
!      do_compress = .false.
    else
      do_compress = .true.
    endif
#endif

    if (do_compress) then
#ifndef USE_HALO_REAL32
      ! compress halo data into the comp_send object
      call deflate_buffer_real64(sendbuf, comp_send, n_doubles=nel)

      ! ! Communicate the index part of the comp_send object, assuming that each MPI rank
      ! ! uses the same number of threads (i.e. comp%n_slices). We concatenate the integers
      ! ! and arrays into a contiguous buffer, for simplicity. Could be done better using a
      ! ! MPI data type, of course.
      ! n_idx = concatenate_index_arrays(comp_send, comp_idx_send)  ! allocates comp_idx_send internally!
      ! allocate(comp_idx_recv(0:n_idx-1))
      ! call MPI_Sendrecv(comp_idx_send, n_idx, MPI_INTEGER, rank_send, mpi_tag,&
      !                   comp_idx_recv, n_idx, MPI_INTEGER, rank_recv, mpi_tag,&
      !                   mpi_comm, MPI_STATUS_IGNORE, ierr)
      ! SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
      ! call decatenate_index_arrays(comp_recv, comp_idx_recv)
      ! deallocate(comp_idx_send); nullify(comp_idx_send)
      ! deallocate(comp_idx_recv); nullify(comp_idx_recv)
      !
      ! if (verb) then
      !   call print_compression_information(comp_send, .false.)
      ! endif
      !
      ! ! communicate the compressed buffer
      ! allocate(comp_recv%buffer(0:comp_recv%n_bytes_deflated_total-1))
      ! call MPI_Sendrecv(comp_send%buffer, comp_send%n_bytes_deflated_total, MPI_BYTE, rank_send, mpi_tag,&
      !                   comp_recv%buffer, comp_recv%n_bytes_deflated_total, MPI_BYTE, rank_recv, mpi_tag,&
      !                   mpi_comm, MPI_STATUS_IGNORE, ierr)
      ! SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
      call sll_s_mpi_sendrecv_compressed_core(comp_send, comp_recv, rank_send, rank_recv, mpi_comm, verb, tag)

      ! decompress halo data into revcbuf
      call inflate_buffer_real64(recvbuf, comp_recv)

      call deallocate_compressed_buffer_obj(comp_send)
      call deallocate_compressed_buffer_obj(comp_recv)
#endif
    else
      call MPI_Sendrecv(sendbuf, nel, mpi_precision, rank_send, tag,&
                        recvbuf, nel, mpi_precision, rank_recv, tag,&
                        mpi_comm, MPI_STATUS_IGNORE, ierr)
      SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
    endif
  end subroutine mpi_sendrecv_compressed


  !> MPI sendrecv functionality, wrapped for a compressed buffer.
  subroutine sll_s_mpi_sendrecv_compressed_core(comp_send, comp_recv, rank_send, rank_recv, mpi_comm, verbose, mpi_tag)
    type(sll_t_compressed_buffer) :: comp_send, comp_recv  ! data structures containing compressed data and offsets
    sll_int32 :: rank_send, rank_recv, mpi_comm
    logical, intent(in), optional :: verbose
    integer, intent(in), optional :: mpi_tag

    integer :: ierr, n_idx
    integer, pointer :: comp_idx_send(:), comp_idx_recv(:)
    logical :: verb
    integer :: tag

    if (present(verbose)) then
      verb = verbose
    else
      verb = .false.
    endif

    if (present(mpi_tag)) then
      tag = mpi_tag
    else
      tag = 0
    endif

    ! Communicate the index part of the comp_send object, assuming that each MPI rank
    ! uses the same number of threads (i.e. comp%n_slices). We concatenate the integers
    ! and arrays into a contiguous buffer, for simplicity. Could be done better using a
    ! MPI data type, of course.
    n_idx = concatenate_index_arrays(comp_send, comp_idx_send)  ! allocates comp_idx_send internally!
    allocate(comp_idx_recv(0:n_idx-1))
    call MPI_Sendrecv(comp_idx_send, n_idx, MPI_INTEGER, rank_send, mpi_tag,&
                      comp_idx_recv, n_idx, MPI_INTEGER, rank_recv, mpi_tag,&
                      mpi_comm, MPI_STATUS_IGNORE, ierr)
    SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
    call decatenate_index_arrays(comp_recv, comp_idx_recv)
    deallocate(comp_idx_send); nullify(comp_idx_send)
    deallocate(comp_idx_recv); nullify(comp_idx_recv)

    if (verb) then
      call print_compression_information(comp_send, .false.)
    endif

    ! communicate the compressed buffer
    allocate(comp_recv%buffer(0:comp_recv%n_bytes_deflated_total-1))
    call MPI_Sendrecv(comp_send%buffer, comp_send%n_bytes_deflated_total, MPI_BYTE, rank_send, tag,&
                      comp_recv%buffer, comp_recv%n_bytes_deflated_total, MPI_BYTE, rank_recv, tag,&
                      mpi_comm, MPI_STATUS_IGNORE, ierr)
    SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
  end subroutine sll_s_mpi_sendrecv_compressed_core


  ! --- halo exchange routine compatible with the "slim" redesing using dynamic buffers ---
  subroutine sll_f_apply_halo_exchange_slim_6d_real64(topo, decomp, arr, id, hw_left, hw_right)
    integer, parameter :: nd = 6  ! we handle 6 dimensions
    ! --- arguments
    type(sll_t_cartesian_topology_6d), intent(in) :: topo
    type(sll_t_decomposition_slim_6d), intent(inout) :: decomp
    sll_real64, dimension(:,:,:,:,:,:), intent(inout) :: arr(decomp%local%mn(1):decomp%local%mx(1), &
                                                             decomp%local%mn(2):decomp%local%mx(2), &
                                                             decomp%local%mn(3):decomp%local%mx(3), &
                                                             decomp%local%mn(4):decomp%local%mx(4), &
                                                             decomp%local%mn(5):decomp%local%mx(5), &
                                                             decomp%local%mn(6):decomp%local%mx(6))
    sll_int32, intent(in) :: id, hw_left, hw_right
    sll_int32 :: halo_block(6,2)
    halo_block(:,1) = decomp%local%mn
    halo_block(:,2) = decomp%local%mx
    call sll_s_apply_halo_exchange_slim_6d_real64( topo, decomp, arr, id, hw_left, hw_right, halo_block )
  end subroutine sll_f_apply_halo_exchange_slim_6d_real64


 ! --- halo exchange routine compatible with the "slim" redesing using dynamic buffers ---
  subroutine sll_s_apply_halo_exchange_slim_6d_real64(topo, decomp, arr, id, hw_left, hw_right, halo_block)
    integer, parameter :: nd = 6  ! we handle 6 dimensions
    ! --- arguments
    type(sll_t_cartesian_topology_6d), intent(in) :: topo
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomp
    sll_real64, intent(inout) :: arr(decomp%local%mn(1):decomp%local%mx(1), &
                                     decomp%local%mn(2):decomp%local%mx(2), &
                                     decomp%local%mn(3):decomp%local%mx(3), &
                                     decomp%local%mn(4):decomp%local%mx(4), &
                                     decomp%local%mn(5):decomp%local%mx(5), &
                                     decomp%local%mn(6):decomp%local%mx(6))
    sll_int32, intent(in) :: id, hw_left, hw_right
    sll_int32, intent(in) :: halo_block(6,2)
    integer :: jd, i, j, k, l, m, n
    integer :: ierr
    logical, save :: first_call = .true.

    ! --- MPI communication-related variables and buffers
    sll_int64, save :: bufsize = 0
    sll_int64 :: nxc  ! total number of elements to be exchanged
    sll_int32 :: nel  ! number of elements to be exchanged at a single MPI call
    integer, dimension(:,:) :: r_rx(nd,2)  ! index ranges for copy operations
    integer, dimension(:,:) :: r_tx(nd,2)
#ifdef USE_HALO_REAL32
    sll_int32, parameter :: nxc_max = 2147483647
    integer, parameter :: mpi_precision = MPI_REAL
    integer, parameter :: word_size = 4
#else
    sll_int32, parameter :: nxc_max = 2147483647
    integer, parameter :: mpi_precision = MPI_DOUBLE_PRECISION
    integer, parameter :: word_size = 8
#endif
    HALO_DTYPE, pointer, save :: sendbuf(:)
    HALO_DTYPE, pointer :: recvbuf(:,:,:,:,:,:)
    integer :: mpi_tag
    logical, save :: use_compression
    logical, save :: compression_verbose
    integer, save :: prec
    ! ------
    logical, parameter :: sendbuf_dump = .false.
    integer, save :: dump_counter = 0
    character(len=32) :: dump_filename

    mpi_tag = 0

    if (first_call) then
      first_call = .false.
      compression_verbose = .false.
#ifdef USE_HALO_REAL32
      if (topo%rank == 0) then
        write(*,*) "sll_m_decomposition::apply_halo_exchange() uses single precision messages"
      endif
      use_compression = .false.
#else
#ifdef USE_ZFP
      use_compression = sll_f_query_environment("SLL_USE_COMPRESSION", .false.)
      if (use_compression) then
        prec = sll_f_query_environment("SLL_ZFP_PRECISION", 32)
        call set_compression_precision(prec)
        if (topo%rank == 0) then
          write(*,*) "sll_m_decomposition::apply_halo_exchange() uses message compression"
          compression_verbose = sll_f_query_environment("SLL_COMPRESSION_VERBOSE", .false.)
        endif
      endif
#else
      use_compression = .false.
#endif
#endif
      nullify(sendbuf)
    endif


    ! --- copy halo cells to the left neighbor / into the right halo buffer
    !
    ! A one-dimensional example, :
    !
    ! |---|----------|**|  # left neighbor
    !
    !            |---|**--------|$$|  # current process
    !
    !                       |---|$$--------|%%|  # right neighbor
    !

    decomp%local%id = id
    decomp%local%halo_right%mn(:)  = halo_block(:,1)
    decomp%local%halo_right%mx(:)  = halo_block(:,2)
    decomp%local%halo_right%nw(:)  = decomp%local%halo_right%mx(:)-decomp%local%halo_right%mn(:)+1

    decomp%local%halo_right%mn(id) = decomp%local%mx(id) + 1
    decomp%local%halo_right%mx(id) = decomp%local%mx(id) + hw_right
    decomp%local%halo_right%nw(id) = hw_right

    if (hw_right > 0) then
      ! --- change the decomposition object state as required by the requested halo_width
#ifdef USE_FMEMPOOL
      if (associated(decomp%local%halo_right%buf)) then
        call mp_release(decomp%local%halo_right%buf)
      endif
      call mp_acquire(decomp%local%halo_right%buf, decomp%local%halo_right%mn, decomp%local%halo_right%mx)
#else
      if (allocated(decomp%local%halo_right%buf)) &
        deallocate(decomp%local%halo_right%buf)
      allocate(decomp%local%halo_right%buf( &
        decomp%local%halo_right%mn(1):decomp%local%halo_right%mx(1), &
        decomp%local%halo_right%mn(2):decomp%local%halo_right%mx(2), &
        decomp%local%halo_right%mn(3):decomp%local%halo_right%mx(3), &
        decomp%local%halo_right%mn(4):decomp%local%halo_right%mx(4), &
        decomp%local%halo_right%mn(5):decomp%local%halo_right%mx(5), &
        decomp%local%halo_right%mn(6):decomp%local%halo_right%mx(6) ))
#endif

!      call mp_statistics()
      recvbuf => decomp%local%halo_right%buf
!      write(*,*) "###", lbound(recvbuf), ubound(recvbuf), size(recvbuf)

      ! --- calculate rx and tx index ranges
      ! index ranges on the computation array to be sent
      r_tx(:,1)  = decomp%local%halo_right%mn(:)!decomp%local%mn(:)
      r_tx(:,2)  = decomp%local%halo_right%mx(:)!decomp%local%mx(:)
      r_tx(id,1) = decomp%local%mn(id)
      r_tx(id,2) = decomp%local%mn(id) + hw_right - 1
      ! index ranges on the buffer are simply its extents
      r_rx(:,1)  = decomp%local%halo_right%mn(:)
      r_rx(:,2)  = decomp%local%halo_right%mx(:)

      if (topo%procs(id) == 1) then
        ! --- we copy the ghost cells directly, assuming periodic BCs for the moment
!$omp parallel default(shared) private(i,j,k,l,m,n)
!$omp do OMP_SCHEDULE OMP_COLLAPSE
        do n = 0, r_rx(6,2)-r_rx(6,1)
          do m = 0, r_rx(5,2)-r_rx(5,1)
            do l = 0, r_rx(4,2)-r_rx(4,1)
              do k = 0, r_rx(3,2)-r_rx(3,1)
                do j = 0, r_rx(2,2)-r_rx(2,1)
                  do i = 0, r_rx(1,2)-r_rx(1,1)
                    recvbuf(r_rx(1,1)+i, r_rx(2,1)+j, r_rx(3,1)+k, &
                            r_rx(4,1)+l, r_rx(5,1)+m, r_rx(6,1)+n) = &
                    arr(r_tx(1,1)+i, r_tx(2,1)+j, r_tx(3,1)+k, &
                        r_tx(4,1)+l, r_tx(5,1)+m, r_tx(6,1)+n)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
!$omp end do
!$omp end parallel
      else
        ! calculate the total number of elements to be exchanged
        nxc = int(hw_right, i64)
        do jd=1,nd
          if (jd == id) then
            cycle
          else
            nxc = nxc * decomp%local%halo_right%nw(jd)
          endif
        end do
        if (nxc > bufsize) then
          if (associated(sendbuf)) deallocate(sendbuf)
          allocate(sendbuf(nxc))
          bufsize = nxc
        end if

        call sll_s_copy_array_to_buffer_6d_real64(arr, decomp%local%mn, decomp%local%mx, sendbuf, r_tx)

        SLL_ASSERT_ALWAYS(nxc <= nxc_max)  ! check if message size is within the allowed limit
        nel = nxc  ! 64 bit -> 32 bit

        if (use_compression) then
          if (sendbuf_dump) then
            write(dump_filename,'(a,i1.1,a,i4.4,a)') "L", id, "_", dump_counter, ".txt"
            write(*,*) trim(dump_filename)
            open(unit=88, file=trim(dump_filename), status='replace')
            write(88,'(E24.16)') (sendbuf(i), i=1,nel)
            close(88)
          endif
          call mpi_sendrecv_compressed_6d_real64(sendbuf, recvbuf, nel, &
                                                 topo%neighbors(2*id-1), topo%neighbors(2*id), topo%comm, &
                                                 compression_verbose)
        else
          call MPI_Sendrecv(sendbuf, nel, mpi_precision, topo%neighbors(2*id-1), mpi_tag,&
                            recvbuf, nel, mpi_precision, topo%neighbors(2*id),   mpi_tag,&
                            topo%comm, MPI_STATUS_IGNORE, ierr)
          SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
        endif
      endif  ! if branch (topo%nprocs > 1)
      nullify(recvbuf)
    else
      ! --- nothing to do
      continue
    endif  ! (hw_right > 0)


    ! --- copy halo cells to the right neighbor / into the left halo buffer
    !
    ! A one-dimensional example, :
    !
    ! |%%%|-------***|--|  # left neighbor
    !
    !            |***|-------$$$|--|  # current process
    !
    !                       |$$$|-------%%%|--|  # right neighbor
    !
    decomp%local%id = id
    decomp%local%halo_left%mn(:)  = halo_block(:,1)
    decomp%local%halo_left%mx(:)  = halo_block(:,2)
    decomp%local%halo_left%nw(:)  = decomp%local%halo_left%mx(:)-decomp%local%halo_left%mn(:)+1

    decomp%local%halo_left%mx(id) = decomp%local%mn(id) - 1
    decomp%local%halo_left%mn(id) = decomp%local%mn(id) - hw_left
    decomp%local%halo_left%nw(id) = hw_left

    if (hw_left > 0) then
      ! --- change the decomposition object's state as required by the requested halo width
#ifdef USE_FMEMPOOL
      if (associated(decomp%local%halo_left%buf)) then
        call mp_release(decomp%local%halo_left%buf)
      endif
      call mp_acquire(decomp%local%halo_left%buf, decomp%local%halo_left%mn, decomp%local%halo_left%mx)
#else
      if (allocated(decomp%local%halo_left%buf)) &
        deallocate(decomp%local%halo_left%buf)
      allocate(decomp%local%halo_left%buf( &
        decomp%local%halo_left%mn(1):decomp%local%halo_left%mx(1), &
        decomp%local%halo_left%mn(2):decomp%local%halo_left%mx(2), &
        decomp%local%halo_left%mn(3):decomp%local%halo_left%mx(3), &
        decomp%local%halo_left%mn(4):decomp%local%halo_left%mx(4), &
        decomp%local%halo_left%mn(5):decomp%local%halo_left%mx(5), &
        decomp%local%halo_left%mn(6):decomp%local%halo_left%mx(6) ))
#endif

      recvbuf => decomp%local%halo_left%buf

      ! --- calculate rx and tx index ranges
      ! index ranges on the computation array to be sent
      r_tx(:,1)  = decomp%local%halo_left%mn(:)!decomp%local%mn(:)
      r_tx(:,2)  = decomp%local%halo_left%mx(:)!decomp%local%mx(:)
      r_tx(id,1) = decomp%local%mx(id) - hw_left + 1
      r_tx(id,2) = decomp%local%mx(id)
      ! index ranges on the buffer are just its extents
      r_rx(:,1)  = decomp%local%halo_left%mn(:)
      r_rx(:,2)  = decomp%local%halo_left%mx(:)

      if (topo%procs(id) == 1) then
        ! --- we copy the ghost cells directly, assuming periodic BCs for the moment
!$omp parallel default(shared) private(i,j,k,l,m,n)
!$omp do OMP_SCHEDULE OMP_COLLAPSE
        do n = 0, r_rx(6,2)-r_rx(6,1)
          do m = 0, r_rx(5,2)-r_rx(5,1)
            do l = 0, r_rx(4,2)-r_rx(4,1)
              do k = 0, r_rx(3,2)-r_rx(3,1)
                do j = 0, r_rx(2,2)-r_rx(2,1)
                  do i = 0, r_rx(1,2)-r_rx(1,1)
                    recvbuf(r_rx(1,1)+i, r_rx(2,1)+j, r_rx(3,1)+k, &
                            r_rx(4,1)+l, r_rx(5,1)+m, r_rx(6,1)+n) = &
                    arr(r_tx(1,1)+i, r_tx(2,1)+j, r_tx(3,1)+k, &
                        r_tx(4,1)+l, r_tx(5,1)+m, r_tx(6,1)+n)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
!$omp end do
!$omp end parallel
      else
        ! calculate the total number of elements to be exchanged
        nxc = int(hw_left, i64)
        do jd=1,nd
          if (jd == id) then
            cycle
          else
            nxc = nxc * decomp%local%halo_left%nw(jd)
          endif
        end do
        if (nxc > bufsize) then
          if (associated(sendbuf)) deallocate(sendbuf)
          allocate(sendbuf(nxc))
          bufsize = nxc
        end if

        call sll_s_copy_array_to_buffer_6d_real64(arr, decomp%local%mn, decomp%local%mx, sendbuf, r_tx)

        SLL_ASSERT_ALWAYS(nxc <= nxc_max)
        nel = nxc

        if (use_compression) then
          if (sendbuf_dump) then
            write(dump_filename,'(a,i1.1,a,i4.4,a)') "R", id, "_", dump_counter, ".txt"
            write(*,*) trim(dump_filename)
            open(unit=88, file=dump_filename, status='replace')
            write(88,'(E24.16)') (sendbuf(i), i=1,nel)
            close(88)
          endif
          call mpi_sendrecv_compressed_6d_real64(sendbuf, recvbuf, nel, &
                                                 topo%neighbors(2*id), topo%neighbors(2*id-1), topo%comm, &
                                                 compression_verbose)
        else
          call MPI_Sendrecv(sendbuf, nel, mpi_precision, topo%neighbors(2*id),   mpi_tag,&
                            recvbuf, nel, mpi_precision, topo%neighbors(2*id-1), mpi_tag,&
                            topo%comm, MPI_STATUS_IGNORE, ierr)
          SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
        endif
      endif  ! if branch (topo%nprocs > 1)
      nullify(recvbuf)
    else
      ! --- nothing to do
      continue
    endif  ! (hw_left > 0)

    dump_counter = dump_counter + 1

#ifdef USE_FMEMPOOL
!    call mp_statistics()
#endif
  end subroutine sll_s_apply_halo_exchange_slim_6d_real64



! --- 3D HALO EXCHANGE BELOW ---



  ! --- halo exchange routine compatible with the "slim" redesing using dynamic buffers ---
  subroutine sll_f_apply_halo_exchange_slim_3d_real64(topo, decomp, arr, id, hw_left, hw_right)
    integer, parameter :: nd = 3  ! we handle 3 dimensions
    type(sll_t_cartesian_topology_3d), intent(in) :: topo
    type(sll_t_decomposition_slim_3d), intent(inout) :: decomp
    sll_real64, dimension(:,:,:), intent(inout) :: arr(decomp%local%mn(1):decomp%local%mx(1), &
                                                       decomp%local%mn(2):decomp%local%mx(2), &
                                                       decomp%local%mn(3):decomp%local%mx(3))
    sll_int32, intent(in) :: id, hw_left, hw_right
    sll_int32 :: halo_block(3,2)
    halo_block(:,1) = decomp%local%mn
    halo_block(:,2) = decomp%local%mx
    call sll_s_apply_halo_exchange_slim_3d_real64( topo, decomp, arr, id, hw_left, hw_right, halo_block )
  end subroutine sll_f_apply_halo_exchange_slim_3d_real64


 ! --- halo exchange routine compatible with the "slim" redesing using dynamic buffers ---
  subroutine sll_s_apply_halo_exchange_slim_3d_real64(topo, decomp, arr, id, hw_left, hw_right, halo_block)
    integer, parameter :: nd = 3  ! we handle 3 dimensions
    type(sll_t_cartesian_topology_3d), intent(in) :: topo
    type(sll_t_decomposition_slim_3d), target, intent(inout) :: decomp
    sll_real64, intent(inout) :: arr(decomp%local%mn(1):decomp%local%mx(1), &
                                     decomp%local%mn(2):decomp%local%mx(2), &
                                     decomp%local%mn(3):decomp%local%mx(3))
    sll_int32, intent(in) :: id, hw_left, hw_right
    sll_int32, intent(in) :: halo_block(3,2)
    integer :: jd, i, j, k
    integer :: ierr
    logical, save :: first_call = .true.

    ! --- MPI communication-related variables and buffers
    sll_int64, save :: bufsize = 0
    sll_int64 :: nxc  ! total number of elements to be exchanged
    sll_int32 :: nel  ! number of elements to be exchanged at a single MPI call
    integer, dimension(:,:) :: r_rx(nd,2)  ! index ranges for copy operations
    integer, dimension(:,:) :: r_tx(nd,2)
#ifdef USE_HALO_REAL32
    sll_int32, parameter :: nxc_max = 2147483647
    integer, parameter :: mpi_precision = MPI_REAL
    integer, parameter :: word_size = 4
#else
    sll_int32, parameter :: nxc_max = 2147483647
    integer, parameter :: mpi_precision = MPI_DOUBLE_PRECISION
    integer, parameter :: word_size = 8
#endif
!    HALO_DTYPE, dimension(:), allocatable, target, save :: sendbuf
    HALO_DTYPE, pointer, save :: sendbuf(:)
    HALO_DTYPE, pointer :: recvbuf(:,:,:)
    integer :: mpi_tag

    ! ------


    mpi_tag = 0

    if (first_call) then
#ifdef USE_HALO_REAL32
      if (topo%rank == 0) then
        write(*,*) "sll_m_decomposition::sll_s_apply_halo_exchange_slim_3d_real64() uses single precision messages"
      endif
#endif
      first_call = .false.
      nullify(sendbuf)
    endif


    decomp%local%id = id
    decomp%local%halo_right%mn(:)  = halo_block(:,1)
    decomp%local%halo_right%mx(:)  = halo_block(:,2)
    decomp%local%halo_right%nw(:)  = decomp%local%halo_right%mx(:)-decomp%local%halo_right%mn(:)+1

    decomp%local%halo_right%mn(id) = decomp%local%mx(id) + 1
    decomp%local%halo_right%mx(id) = decomp%local%mx(id) + hw_right
    decomp%local%halo_right%nw(id) = hw_right

    if (hw_right > 0) then
      ! --- change the decomposition object's state as required by the requested halo_width
      if (allocated(decomp%local%halo_right%buf)) &
        deallocate(decomp%local%halo_right%buf)

      allocate(decomp%local%halo_right%buf( &
        decomp%local%halo_right%mn(1):decomp%local%halo_right%mx(1), &
        decomp%local%halo_right%mn(2):decomp%local%halo_right%mx(2), &
        decomp%local%halo_right%mn(3):decomp%local%halo_right%mx(3)))

      recvbuf => decomp%local%halo_right%buf

      ! --- calculate rx and tx index ranges
      ! index ranges on the computation array to be sent
      r_tx(:,1)  = decomp%local%halo_right%mn(:)!decomp%local%mn(:)
      r_tx(:,2)  = decomp%local%halo_right%mx(:)!decomp%local%mx(:)
      r_tx(id,1) = decomp%local%mn(id)
      r_tx(id,2) = decomp%local%mn(id) + hw_right - 1
      ! index ranges on the buffer are simply its extents
      r_rx(:,1)  = decomp%local%halo_right%mn(:)
      r_rx(:,2)  = decomp%local%halo_right%mx(:)

      if (topo%procs(id) == 1) then
        ! --- we copy the ghost cells directly, assuming periodic BCs for the moment
!$omp parallel default(shared) private(i,j,k)
!$omp do OMP_SCHEDULE OMP_COLLAPSE
        do k = 0, r_rx(3,2)-r_rx(3,1)
          do j = 0, r_rx(2,2)-r_rx(2,1)
            do i = 0, r_rx(1,2)-r_rx(1,1)
              recvbuf(r_rx(1,1)+i, r_rx(2,1)+j, r_rx(3,1)+k) = &
              arr(r_tx(1,1)+i, r_tx(2,1)+j, r_tx(3,1)+k)
            enddo
          enddo
        enddo
!$omp end do
!$omp end parallel
      else
        ! calculate the total number of elements to be exchanged
        nxc = int(hw_right, i64)
        do jd=1,nd
          if (jd == id) then
            cycle
          else
            nxc = nxc * decomp%local%halo_right%nw(jd)
          endif
        end do
        if (nxc > bufsize) then
          if (associated(sendbuf)) deallocate(sendbuf)
          allocate(sendbuf(nxc))
          bufsize = nxc
        end if

        call sll_s_copy_array_to_buffer_3d_real64(arr, decomp%local%mn, decomp%local%mx, sendbuf, r_tx)

        SLL_ASSERT_ALWAYS(nxc <= nxc_max)  ! check if message size is within the allowed limit
        nel = nxc  ! 64 bit -> 32 bit
        call MPI_Sendrecv(sendbuf,                     nel, mpi_precision, topo%neighbors(2*id-1), mpi_tag,&
                          decomp%local%halo_right%buf, nel, mpi_precision, topo%neighbors(2*id),   mpi_tag,&
                          topo%comm, MPI_STATUS_IGNORE, ierr)
        SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
      endif  ! if branch (topo%nprocs > 1)
      nullify(recvbuf)
    else
      ! --- nothing to do
      continue
    endif  ! (hw_right > 0)


    ! --- copy halo cells to the right neighbor / into the left halo buffer

    decomp%local%id = id
    decomp%local%halo_left%mn(:)  = halo_block(:,1)
    decomp%local%halo_left%mx(:)  = halo_block(:,2)
    decomp%local%halo_left%nw(:)  = decomp%local%halo_left%mx(:)-decomp%local%halo_left%mn(:)+1

    decomp%local%halo_left%mx(id) = decomp%local%mn(id) - 1
    decomp%local%halo_left%mn(id) = decomp%local%mn(id) - hw_left
    decomp%local%halo_left%nw(id) = hw_left

    if (hw_left > 0) then
      ! --- change the decomposition object's state as required by the requested halo width
      if (allocated(decomp%local%halo_left%buf)) &
           deallocate(decomp%local%halo_left%buf)
      allocate(decomp%local%halo_left%buf( &
        decomp%local%halo_left%mn(1):decomp%local%halo_left%mx(1), &
        decomp%local%halo_left%mn(2):decomp%local%halo_left%mx(2), &
        decomp%local%halo_left%mn(3):decomp%local%halo_left%mx(3)))

      recvbuf => decomp%local%halo_left%buf

      ! --- calculate rx and tx index ranges
      ! index ranges on the computation array to be sent
      r_tx(:,1)  = decomp%local%halo_left%mn(:)!decomp%local%mn(:)
      r_tx(:,2)  = decomp%local%halo_left%mx(:)!decomp%local%mx(:)
      r_tx(id,1) = decomp%local%mx(id) - hw_left + 1
      r_tx(id,2) = decomp%local%mx(id)
      ! index ranges on the buffer are just its extents
      r_rx(:,1)  = decomp%local%halo_left%mn(:)
      r_rx(:,2)  = decomp%local%halo_left%mx(:)

      if (topo%procs(id) == 1) then
        ! --- we copy the ghost cells directly, assuming periodic BCs for the moment
!$omp parallel default(shared) private(i,j,k)
!$omp do OMP_SCHEDULE OMP_COLLAPSE
        do k = 0, r_rx(3,2)-r_rx(3,1)
          do j = 0, r_rx(2,2)-r_rx(2,1)
            do i = 0, r_rx(1,2)-r_rx(1,1)
              recvbuf(r_rx(1,1)+i, r_rx(2,1)+j, r_rx(3,1)+k) = &
              arr(r_tx(1,1)+i, r_tx(2,1)+j, r_tx(3,1)+k)
            enddo
          enddo
        enddo
!$omp end do
!$omp end parallel
      else
        ! calculate the total number of elements to be exchanged
        nxc = int(hw_left, i64)
        do jd=1,nd
          if (jd == id) then
            cycle
          else
            nxc = nxc * decomp%local%halo_left%nw(jd)
          endif
        end do
        if (nxc > bufsize) then
          if (associated(sendbuf)) deallocate(sendbuf)
          allocate(sendbuf(nxc))
          bufsize = nxc
        end if

        call sll_s_copy_array_to_buffer_3d_real64(arr, decomp%local%mn, decomp%local%mx, sendbuf, r_tx)

        SLL_ASSERT_ALWAYS(nxc <= nxc_max)
        nel = nxc
        call MPI_Sendrecv(sendbuf,                    nel, mpi_precision, topo%neighbors(2*id),   mpi_tag,&
                          decomp%local%halo_left%buf, nel, mpi_precision, topo%neighbors(2*id-1), mpi_tag,&
                          topo%comm, MPI_STATUS_IGNORE, ierr)
        SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
      endif  ! if branch (topo%nprocs > 1)
      nullify(recvbuf)
    else
      ! --- nothing to do
      continue
    endif  ! (hw_left > 0)
  end subroutine sll_s_apply_halo_exchange_slim_3d_real64


  subroutine sll_s_apply_bc_exchange_slim_6d_real64(topo, decomp, id)
    type(sll_t_cartesian_topology_6d), intent(in) :: topo
    type(sll_t_decomposition_slim_6d), intent(inout) :: decomp
    sll_int32, intent(in) :: id

    integer :: ierr
    sll_int32 :: nel  ! number of elements to be exchanged at a single MPI call
#ifdef USE_HALO_REAL32
    integer, parameter :: mpi_precision = MPI_REAL
#else
    integer, parameter :: mpi_precision = MPI_DOUBLE_PRECISION
#endif
    integer, parameter :: mpi_tag=1024   ! large tag number to enable overlap with halo exchange

    if (allocated(decomp%local%bc_left_send) .and. allocated(decomp%local%bc_right)) then
      nel = size(decomp%local%bc_left_send)
      call MPI_Sendrecv(decomp%local%bc_left_send, nel, mpi_precision, topo%neighbors(2*id),   mpi_tag,&
                        decomp%local%bc_right,     nel, mpi_precision, topo%neighbors(2*id-1), mpi_tag,&
                        topo%comm, MPI_STATUS_IGNORE, ierr)
      SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
    endif

    if (allocated(decomp%local%bc_right_send) .and. allocated(decomp%local%bc_left)) then
      nel = size(decomp%local%bc_right_send)
      call MPI_Sendrecv(decomp%local%bc_right_send, nel, mpi_precision, topo%neighbors(2*id-1), mpi_tag,&
                        decomp%local%bc_left,       nel, mpi_precision, topo%neighbors(2*id),   mpi_tag,&
                        topo%comm, MPI_STATUS_IGNORE, ierr)
      SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
    endif
  end subroutine sll_s_apply_bc_exchange_slim_6d_real64


  function sll_f_select_dim(id)
    logical, dimension(:) :: sll_f_select_dim(6)
    sll_int32, optional :: id

    if (present(id)) then
      sll_f_select_dim = .false.
      sll_f_select_dim(id) = .true.
    else
      sll_f_select_dim = .true.
    endif
  end function sll_f_select_dim


  subroutine sll_s_allocate_bc_buffers_6d(decomp, id)
    type(sll_t_decomposition_slim_6d), intent(inout) :: decomp
    sll_int32, intent(in) :: id  ! dimension to be skipped

    sll_int32 :: idx_mn(5), idx_mx(5)
    sll_int32 :: i, j

    j = 1
    do i=1, 6
      if (i == id) cycle
      idx_mn(j) = decomp%local%mn(i)
      idx_mx(j) = decomp%local%mx(i)
      j = j + 1
   enddo

   call sll_s_allocate_bc_buffers_6d_part(decomp, id, idx_mn, idx_mx)
  end subroutine sll_s_allocate_bc_buffers_6d

  ! Helper routine to allocate boundary condition buffers
  subroutine sll_s_allocate_bc_buffers_6d_part(decomp, id, idx_mn, idx_mx)
    type(sll_t_decomposition_slim_6d), intent(inout) :: decomp
    sll_int32, intent(in) :: id  ! dimension to be skipped

    sll_int32, intent(in) :: idx_mn(5), idx_mx(5)
!!$    sll_int32 :: i, j
!!$
!!$    j = 1
!!$    do i=1, 6
!!$      if (i == id) cycle
!!$      idx_mn(j) = decomp%local%mn(i)
!!$      idx_mx(j) = decomp%local%mx(i)
!!$      j = j + 1
!!$    enddo

    call sll_s_deallocate_bc_buffers(decomp)

    allocate(decomp%local%bc_right_send(idx_mn(1):idx_mx(1), &
                                        idx_mn(2):idx_mx(2), &
                                        idx_mn(3):idx_mx(3), &
                                        idx_mn(4):idx_mx(4), &
                                        idx_mn(5):idx_mx(5)))
    allocate(decomp%local%bc_left_send( idx_mn(1):idx_mx(1), &
                                        idx_mn(2):idx_mx(2), &
                                        idx_mn(3):idx_mx(3), &
                                        idx_mn(4):idx_mx(4), &
                                        idx_mn(5):idx_mx(5)))
    allocate(decomp%local%bc_right(     idx_mn(1):idx_mx(1), &
                                        idx_mn(2):idx_mx(2), &
                                        idx_mn(3):idx_mx(3), &
                                        idx_mn(4):idx_mx(4), &
                                        idx_mn(5):idx_mx(5)))
    allocate(decomp%local%bc_left(      idx_mn(1):idx_mx(1), &
                                        idx_mn(2):idx_mx(2), &
                                        idx_mn(3):idx_mx(3), &
                                        idx_mn(4):idx_mx(4), &
                                        idx_mn(5):idx_mx(5)))
  end subroutine sll_s_allocate_bc_buffers_6d_part


  ! safety/helper routine to deallocate boundary condition buffers
  subroutine sll_s_deallocate_bc_buffers(decomp)
    type(sll_t_decomposition_slim_6d), intent(inout) :: decomp

    if (allocated(decomp%local%bc_right_send)) deallocate(decomp%local%bc_right_send)
    if (allocated(decomp%local%bc_left_send))  deallocate(decomp%local%bc_left_send)
    if (allocated(decomp%local%bc_right))      deallocate(decomp%local%bc_right)
    if (allocated(decomp%local%bc_left))       deallocate(decomp%local%bc_left)
  end subroutine sll_s_deallocate_bc_buffers


  ! DEBUG output routine
  subroutine dump_ascii(file, arr)
    character(len=*), intent(in) :: file
#ifdef USE_HALO_REAL32
    sll_real32, dimension(:), intent(in) :: arr
#else
    sll_real64, dimension(:), intent(in) :: arr
#endif
    integer, parameter :: fd = 67
    integer :: i
    open(unit=fd, file=file)
    do i=1, size(arr)
      write(fd,*) arr(i)
    end do
    close(fd)
  end subroutine


  ! DEBUG output routine
  subroutine dump_binary(filename, array)
    character(len=*), intent(in) :: filename
    sll_real64, intent(in) :: array(:)
    integer, parameter :: iunit = 67
    open(iunit, file=TRIM(filename), status='replace', form='unformatted')
    write(iunit) array
    close(iunit)
  end subroutine dump_binary


  ! DEBUG output routine
  subroutine dump_ascii_6d(file, decomp, arr)
    character(len=*), intent(in) :: file
    type(sll_t_decomposition_6d), intent(in) :: decomp
    sll_real64, dimension(:,:,:,:,:,:), intent(inout) :: arr(decomp%local%lo(1):decomp%local%hi(1), &
                                                             decomp%local%lo(2):decomp%local%hi(2), &
                                                             decomp%local%lo(3):decomp%local%hi(3), &
                                                             decomp%local%lo(4):decomp%local%hi(4), &
                                                             decomp%local%lo(5):decomp%local%hi(5), &
                                                             decomp%local%lo(6):decomp%local%hi(6))
    integer, parameter :: fd = 67
    sll_int32 :: i,j,k,l,m,n
    open(unit=fd, file=file)
    do n=decomp%local%mn(6),decomp%local%mx(6)
       do m=decomp%local%mn(5),decomp%local%mx(5)
          do l=decomp%local%mn(4),decomp%local%mx(4)
             do k=decomp%local%mn(3),decomp%local%mx(3)
                do j=decomp%local%mn(2),decomp%local%mx(2)
                   do i=decomp%local%mn(1),decomp%local%mx(1)
                      write(fd,'(I3,I3,I3,I3,I3,I3,F20.16)') i, j, k, l, m, n, arr(i,j,k,l,m,n)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    close(fd)
  end subroutine


  subroutine dump_dd_information(file, topo, decomp)
    character(len=*), intent(in) :: file
    type(sll_t_cartesian_topology_6d), intent(in) :: topo
    type(sll_t_decomposition_6d), intent(in) :: decomp
    integer, parameter :: fd = 67
    open(unit=fd, file=file)
    write(fd,'(A,1I4)')  "topo%rank            : ", topo%rank
    write(fd,'(A,1I4)')  "topo%nprocs          : ", topo%nprocs
    write(fd,'(A,6I4)')  "topo%procs           : ", topo%procs
    write(fd,'(A,6L4)')  "topo%periodic        : ", topo%periodic
    write(fd,'(A,6I4)')  "topo%coords          : ", topo%coords
    write(fd,'(A,12I4)') "topo%neighbors       : ", topo%neighbors
    write(fd,'(A,6I4)')  "decomp%global        : ", decomp%global
    write(fd,'(A,6I4)')  "decomp%local%mn      : ", decomp%local%mn
    write(fd,'(A,6I4)')  "decomp%local%mx      : ", decomp%local%mx
    write(fd,'(A,6I4)')  "decomp%local%hw      : ", decomp%local%hw
    write(fd,'(A,6I4)')  "decomp%local%lo      : ", decomp%local%lo
    write(fd,'(A,6I4)')  "decomp%local%hi      : ", decomp%local%hi
    write(fd,'(A,6I4)')  "decomp%local%nw      : ", decomp%local%nw
    write(fd,'(A,6I4)')  "decomp%local%gw      : ", decomp%local%gw
    write(fd,'(A,6I4)')  "decomp%local%tx_lolo : ", decomp%local%tx_lolo
    write(fd,'(A,6I4)')  "decomp%local%tx_lohi : ", decomp%local%tx_lohi
    write(fd,'(A,6I4)')  "decomp%local%tx_hilo : ", decomp%local%tx_hilo
    write(fd,'(A,6I4)')  "decomp%local%tx_hihi : ", decomp%local%tx_hihi
    write(fd,'(A,6I4)')  "decomp%local%rx_lolo : ", decomp%local%rx_lolo
    write(fd,'(A,6I4)')  "decomp%local%rx_lohi : ", decomp%local%rx_lohi
    write(fd,'(A,6I4)')  "decomp%local%rx_hilo : ", decomp%local%rx_hilo
    write(fd,'(A,6I4)')  "decomp%local%rx_hihi : ", decomp%local%rx_hihi
    close(fd)
  end subroutine

  ! dummy routine to check correct linking, do not call
  subroutine dummy_mempool()
#ifdef USE_FMEMPOOL
    call mp_statistics()
#endif
  end subroutine dummy_mempool


  function sll_f_set_process_grid(mpi_world_size, process_grid_par) result(process_grid)
    integer :: process_grid(6), i, j
    integer, intent(in) :: mpi_world_size  ! number of MPI processes
    integer, intent(in), optional :: process_grid_par(6)  ! process grid set via parameter file, shall be initialized to zero otherwise
    logical :: swap_process_grid

    swap_process_grid = sll_f_query_environment("SLL_SWAP_PROCESS_GRID", .false.)

    if ((present(process_grid_par)) .and. (product(process_grid_par) == mpi_world_size)) then
      ! apply the process grid passed as a parameter (i.e. set via the parameter file)
      process_grid = process_grid_par
    else
      !> Create a distribution of the MPI processes over the dimensions.
      !> There is a convenience function provided by MPI as well, however
      !> we do this by hand for the moment to have better control.

      ! divide starting with the last dimension (may be swapped later, see below)
      select case(mpi_world_size)
      case(1)
         process_grid = [1,1,1,1,1,1]
      case(2)
         process_grid = [1,1,1,1,1,2]
      case(4)
         process_grid = [1,1,1,1,2,2]
      case(8)
         process_grid = [1,1,1,2,2,2]
      case(16)
         process_grid = [1,1,2,2,2,2]
      case(32)
         process_grid = [1,2,2,2,2,2]
      case(64)
         process_grid = [2,2,2,2,2,2]
      case(96)
         process_grid = [2,2,2,2,2,3]
      case(128)
         process_grid = [2,2,2,2,2,4]
      case(256)
         process_grid = [2,2,2,2,4,4]
      case(512)
         process_grid = [2,2,2,4,4,4]
      case(1024)
         process_grid = [2,2,4,4,4,4]
      case(2048)
         process_grid = [2,4,4,4,4,4]
      case(4096)
         process_grid = [4,4,4,4,4,4]
      case(8192)
         process_grid = [4,4,4,4,4,8]
      case(16384)
         process_grid = [4,4,4,4,8,8]
      case(32768)
         process_grid = [4,4,4,8,8,8]
      case(65536)
         process_grid = [4,4,8,8,8,8]
      case(131072)
         process_grid = [4,8,8,8,8,8]
      case(262144)
         process_grid = [8,8,8,8,8,8]
      case default
         write(*,*) "Error: No process topology implemented for ", mpi_world_size, " processes. STOP."
      end select

      if (swap_process_grid) then
        do i=1,3
          j = process_grid(7-i)
          process_grid(7-i) = process_grid(i)
          process_grid(i) = j
        enddo
      endif
    endif
  end function sll_f_set_process_grid

end module sll_m_decomposition
