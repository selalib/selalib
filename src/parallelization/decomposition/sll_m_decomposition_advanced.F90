!**************************************************************
!
!  Domain decomposition module for Selalib (1D, ..., 6d)
!
!  Extension of the base module to support overlapping communication.
!
!**************************************************************
!> @ingroup decomposition
!> @brief
!> Module providing data structures and tools to implement domain decompositions.
!> @author
!> Klaus Reuter, Max Planck Computing and Data Facility (MPCDF)

! The hierarchy of the classes in this module is as follows:
!
! sll_t_decomposition, contains a single element of:
!   sll_t_decomposition__local, contains an array of:
!     sll_t_decomposition__dimension, contains an array of:
!       sll_t_decomposition__block, contains an array of:
!         sll_t_decomposition__buffer, specializes to:
!           sll_t_decomposition__buffer_3d
!           sll_t_decomposition__buffer_6d
!
! This complex setup is necessary to implement multiple dimensions,
! blocks per dimension, buffers per block.

! Preprocessor macro:  Use single precision for the halo exchange, ie.
! basically a poor man's lossy compression to speed up MPI communication.
! NOTE: The define should be passed by the build system, please see
! the build script <compile_mpcdf.sh> for an example.
#ifdef USE_HALO_REAL32
#define HALO_DTYPE sll_real32
#else
#define HALO_DTYPE sll_real64
#endif


module sll_m_decomposition_advanced
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use mpi, only: &
    mpi_double_precision, &
    mpi_real, &
    mpi_sendrecv, &
    mpi_status_ignore, &
    mpi_success

 ! use sll_m_decomposition
  use sll_m_decomposition, only: &
     sll_t_cartesian_topology_6d, &
     sll_s_copy_array_to_buffer_6d_real64, &
     sll_s_mpi_sendrecv_compressed_core

#ifdef _OPENMP
  use omp_lib
#define OMP_COLLAPSE collapse(2)
#define OMP_SCHEDULE schedule(static)
!#define OMP_SCHEDULE schedule(dynamic)
#endif

! FMEMPOOL is currently not compatible with OpenMP nested parallelism!
#undef USE_FMEMPOOL

#ifdef USE_FMEMPOOL
  use fmempool
#endif

  use sll_m_compression

  use sll_m_utilities, only: &
       sll_f_query_environment

  implicit none

  public :: &
    sll_t_decomposition__buffer, &
    sll_t_decomposition__block, &
    sll_t_decomposition__dimension, &
    sll_t_decomposition__local, &
    sll_t_decomposition
    ! sll_t_decomposition__buffer_3d, &
    ! sll_t_decomposition__buffer_6d, &

  public :: &
    sll_f_new_cartesian_domain_decomposition, &
    sll_s_deallocate_cartesian_domain_decomposition, &
    sll_s_dd_define_blocks, &
    sll_s_dd_allocate_buffer, &
    sll_s_dd_deallocate_buffer, &
    sll_f_get_mem_6d_from_buffer_obj, &
    sll_s_post_halo_exchange_real64, &
    sll_s_wait_halo_exchange, &
    sll_s_blocking_halo_exchange, &
    sll_f_dd_get_n_blocks


  private


  ! Temporary note: Only the decomposition stuff is re-implemented here,
  ! the topology stuff from sll_m_decomposition can be reused for the moment.


  ! --- nested hierarchy of decomposition object types, supporting 3d and 6d decompositions ---
  !     names were shortened, nesting is indicated by the double underscore __
  type :: sll_t_decomposition__buffer
    character(len=64) :: label
    sll_int32, allocatable :: mn(:)  ! min index specific to the buffer
    sll_int32, allocatable :: mx(:)  ! max index specific to the buffer
    sll_int32, allocatable :: nw(:)  ! net width of the buffer
    logical :: valid
    ! --- buffers for MPI communication ---
    HALO_DTYPE, pointer :: sendbuf(:) => null()
    ! entries only used by blocking MPI_Sendrecv()
    sll_int32 :: nel
    sll_int32 :: mpi_source
    sll_int32 :: mpi_dest
    sll_int32 :: mpi_tag
    sll_int32 :: mpi_precision
    ! entry only used by non-blocking communication (isend,irecv)
    sll_int32 :: request(2)
    ! support for compressed mpi messages
    logical :: use_compression = .false.
    logical :: compression_verbose = .false.
    type(sll_t_compressed_buffer) :: comp_sendbuf, comp_recvbuf
    ! temporarily placed 6d array here because of 'select type' compiler bug with intel/17.0.4
    HALO_DTYPE, pointer :: mem(:,:,:,:,:,:) => null()
  end type sll_t_decomposition__buffer

  ! type, extends(sll_t_decomposition__buffer) :: sll_t_decomposition__buffer_3d
  !   HALO_DTYPE, pointer :: mem(:,:,:)
  ! end type sll_t_decomposition__buffer_3d
  !
  ! type, extends(sll_t_decomposition__buffer) :: sll_t_decomposition__buffer_6d
  !   HALO_DTYPE, pointer :: mem(:,:,:,:,:,:)
  ! end type sll_t_decomposition__buffer_6d

  type :: sll_t_decomposition__block
    sll_int32, allocatable :: mn(:)  ! min index specific to the block, w/o halo
    sll_int32, allocatable :: mx(:)  ! max index specific to the block, w/o halo
    sll_int32, allocatable :: nw(:)  ! net width of the block, w/o halo
    class(sll_t_decomposition__buffer), pointer :: buffer(:) => null()
  end type sll_t_decomposition__block

  type :: sll_t_decomposition__dimension
    sll_int32 :: block_dim  ! dimension in which the blocking takes place
    type(sll_t_decomposition__block), pointer :: block(:) => null()
  end type sll_t_decomposition__dimension

  type :: sll_t_decomposition__local
    sll_int32, allocatable :: mn(:)  ! min index, w/o halo
    sll_int32, allocatable :: mx(:)  ! max index, w/o halo
    sll_int32, allocatable :: nw(:)  ! net width of the array, w/o halo
    type(sll_t_decomposition__dimension), pointer :: dimension(:) => null()
  end type sll_t_decomposition__local

  type :: sll_t_decomposition
    sll_int32, allocatable :: global(:)
    type(sll_t_decomposition__local) :: local
  end type sll_t_decomposition

  ! --- maximum number of buffers per block, increase if necessary
  sll_int32, parameter :: max_buffers_per_block = 8

!#define VERBOSE(name) write(*,*) name // " called"
#define VERBOSE(name)

  sll_int32, parameter :: LEFT_BUFFER = 1
  sll_int32, parameter :: RIGHT_BUFFER = 0

contains

  subroutine sll_s_destruct_decomposition__buffer(self)
    class(sll_t_decomposition__buffer) :: self
    sll_int32 :: ierr
    VERBOSE("sll_s_destruct_decomposition__buffer")
    if (self%valid) then
!       select type(ptr => self)
!         class is (sll_t_decomposition__buffer_3d)
! #ifdef USE_FMEMPOOL
!           call mp_release(ptr%mem)
! #else
!           SLL_DEALLOCATE(ptr%mem, ierr)
! #endif
!         class is (sll_t_decomposition__buffer_6d)
#ifdef USE_FMEMPOOL
          call mp_release(self%mem)
#else
          SLL_DEALLOCATE(self%mem, ierr)
#endif
      !   class default
      !     SLL_ERROR("sll_s_destruct_decomposition__buffer", "unknown buffer data structure")
      ! end select
      SLL_DEALLOCATE_ARRAY(self%mn, ierr)
      SLL_DEALLOCATE_ARRAY(self%mx, ierr)
      SLL_DEALLOCATE_ARRAY(self%nw, ierr)
      self%valid = .false.
    endif
  end subroutine sll_s_destruct_decomposition__buffer

  subroutine sll_s_construct_decomposition__block(self, nd)
    type(sll_t_decomposition__block) :: self
    sll_int32, intent(in) :: nd
    sll_int32 :: i, ierr
    character(len=16) :: str
    VERBOSE("sll_s_construct_decomposition__block")
    ! select case (nd)
    !   case (3)
    !     allocate(sll_t_decomposition__buffer_3d::self%buffer(max_buffers_per_block))
    !   case (6)
    !     allocate(sll_t_decomposition__buffer_6d::self%buffer(max_buffers_per_block))
    !   case default
    !     write (str, *) nd
    !     SLL_ERROR("sll_s_construct_decomposition__block", "nd = " // trim(str) // " is not supported.")
    ! end select
    allocate(self%buffer(max_buffers_per_block))
    do i=1,max_buffers_per_block
      self%buffer(i)%valid = .false.
    enddo
    SLL_ALLOCATE(self%mn(nd), ierr)
    SLL_ALLOCATE(self%mx(nd), ierr)
    SLL_ALLOCATE(self%nw(nd), ierr)
  end subroutine sll_s_construct_decomposition__block

  subroutine sll_s_destruct_decomposition__block(self)
    type(sll_t_decomposition__block) :: self
    sll_int32 :: i, ierr
    VERBOSE("sll_s_destruct_decomposition__block")
    if (associated(self%buffer)) then
      do i=1,size(self%buffer)
        call sll_s_destruct_decomposition__buffer(self%buffer(i))
      enddo
      SLL_DEALLOCATE(self%buffer, ierr)
    endif
    if (allocated(self%mn)) then
      SLL_DEALLOCATE_ARRAY(self%mn, ierr)
    endif
    if (allocated(self%mx)) then
      SLL_DEALLOCATE_ARRAY(self%mx, ierr)
    endif
    if (allocated(self%nw)) then
      SLL_DEALLOCATE_ARRAY(self%nw, ierr)
    endif
  end subroutine sll_s_destruct_decomposition__block


  subroutine sll_s_construct_decomposition__dimension(self, nd)
    type(sll_t_decomposition__dimension) :: self
    sll_int32, intent(in) :: nd
    sll_int32 :: i, ierr
    sll_int32, parameter :: n_default_blocks = 0  ! disabled
    VERBOSE("sll_s_construct_decomposition__dimension")
    self%block_dim = 0
  end subroutine sll_s_construct_decomposition__dimension

  subroutine sll_s_destruct_decomposition__dimension(self)
    type(sll_t_decomposition__dimension) :: self
    sll_int32 :: i, ierr
    VERBOSE("sll_s_destruct_decomposition__dimension")
    if (associated(self%block)) then
      do i=1,size(self%block)
        call sll_s_destruct_decomposition__block(self%block(i))
      enddo
      SLL_DEALLOCATE(self%block, ierr)
    endif
    self%block_dim = 0
  end subroutine sll_s_destruct_decomposition__dimension


  subroutine sll_s_construct_decomposition__local(self, topology, grid_size, nd)
    type(sll_t_decomposition__local) :: self
    type(sll_t_cartesian_topology_6d), pointer, intent(in) :: topology  ! TODO: support generic nd topology
    sll_int32, intent(in) :: nd
    sll_int32, intent(in) :: grid_size(nd)
    sll_int32 :: i, ierr
    sll_int32 :: lp, l0, l1
    VERBOSE("sll_s_construct_decomposition__local")
    ! --- initialize index values
    SLL_ALLOCATE(self%mn(nd), ierr)
    SLL_ALLOCATE(self%mx(nd), ierr)
    SLL_ALLOCATE(self%nw(nd), ierr)
    ! loop over dimensions and compute index values for each dimension
    do i=1,nd
       ! compute the local number of grid points
       lp = grid_size(i) / topology%procs(i)
       ! compute the lower local index bound (the coords array starts at zero)
       l0 = 1 + topology%coords(i) * lp
       ! compute the upper local index bound (the coords array starts at zero)
       l1 = (topology%coords(i) + 1) * lp
       ! ---
       self%mn(i) = l0
       self%mx(i) = l1
       self%nw(i) = lp
    end do
    SLL_ALLOCATE(self%dimension(nd), ierr)
    do i=1,nd
      call sll_s_construct_decomposition__dimension(self%dimension(i), nd)
    end do
  end subroutine sll_s_construct_decomposition__local

  subroutine sll_s_destruct_decomposition__local(self)
    type(sll_t_decomposition__local) :: self
    sll_int32 :: i, ierr
    VERBOSE("sll_s_destruct_decomposition__local")
    if (associated(self%dimension)) then
      do i=1,size(self%dimension)
        call sll_s_destruct_decomposition__dimension(self%dimension(i))
      end do
      SLL_DEALLOCATE(self%dimension, ierr)
    endif
    if (allocated(self%mn)) then
      SLL_DEALLOCATE_ARRAY(self%mn, ierr)
    endif
    if (allocated(self%mx)) then
      SLL_DEALLOCATE_ARRAY(self%mx, ierr)
    endif
    if (allocated(self%nw)) then
      SLL_DEALLOCATE_ARRAY(self%nw, ierr)
    endif
  end subroutine sll_s_destruct_decomposition__local


  ! --- public functions below ---

  ! Recursively build up a domain decomposition object (including nested objects).
  function sll_f_new_cartesian_domain_decomposition(topology, grid_size, nd)
    type(sll_t_decomposition), pointer :: sll_f_new_cartesian_domain_decomposition
    sll_int32, intent(in) :: nd
    type(sll_t_cartesian_topology_6d), pointer, intent(in) :: topology  ! TODO: support generic nd topology
    sll_int32, intent(in) :: grid_size(nd)
    sll_int32 :: i, ierr
    type(sll_t_decomposition), pointer :: self  ! short convenience alias
    VERBOSE("sll_f_new_cartesian_domain_decomposition")
    ! --- initial consistency checks
    SLL_ASSERT_ALWAYS(size(topology%procs) == nd)
    do i=1,nd
       SLL_ASSERT(mod(grid_size(i), topology%procs(i)) == 0)
    end do
    ! --- initialize index values
    SLL_ALLOCATE(self, ierr)
    SLL_ALLOCATE(self%global(nd), ierr)
    self%global(:) = grid_size(:)
    ! --- recursively initialize decomposition object
    call sll_s_construct_decomposition__local(self%local, topology, grid_size, nd)
    sll_f_new_cartesian_domain_decomposition => self
    nullify(self)
  end function sll_f_new_cartesian_domain_decomposition


  ! Recursively tear down a domain decomposition object (including nested objects and memory).
  subroutine sll_s_deallocate_cartesian_domain_decomposition(decomposition)
    type(sll_t_decomposition), pointer :: decomposition
    sll_int32 :: ierr
    VERBOSE("sll_s_deallocate_cartesian_domain_decomposition")
    ! --- recursively deallocate/invalidata decomposition object
    call sll_s_destruct_decomposition__local(decomposition%local)
    SLL_DEALLOCATE_ARRAY(decomposition%global, ierr)
    SLL_DEALLOCATE(decomposition, ierr)
    nullify(decomposition)
  end subroutine sll_s_deallocate_cartesian_domain_decomposition


  !> Initialize blocking of halo cells for dimension \a id
  subroutine sll_s_dd_define_blocks(decomposition, id, n_blocks, block_dim)
    type(sll_t_decomposition), pointer, intent(in) :: decomposition !> Decomposition
    sll_int32, intent(in) :: id !> dimensions to be blocked
    sll_int32, intent(in) :: n_blocks !> number of blocks in this dimension
    sll_int32, intent(in) :: block_dim !> dimension over which the halos shall be blocked
    sll_int32 :: i, nd, ierr, nw, mn
    nd = size(decomposition%global)
    SLL_ASSERT_ALWAYS( id <= nd )
    SLL_ASSERT_ALWAYS( block_dim <= nd )
    SLL_ASSERT_ALWAYS( associated(decomposition%local%dimension) )
    ! clean up first, if blocks exist already
    if (associated(decomposition%local%dimension(id)%block)) then
      do i=1,size(decomposition%local%dimension(id)%block)
        call sll_s_destruct_decomposition__block(decomposition%local%dimension(id)%block(i))
      enddo
      SLL_DEALLOCATE(decomposition%local%dimension(id)%block, ierr)
    endif
    SLL_ASSERT_ALWAYS(decomposition%local%nw(block_dim) >= n_blocks)
    SLL_ASSERT_ALWAYS(mod(decomposition%local%nw(block_dim), n_blocks) == 0)
    SLL_ALLOCATE(decomposition%local%dimension(id)%block(n_blocks), ierr)
    do i=1,n_blocks
      ! allocate index arrays, the buffer array, but not the memory fields inside
      call sll_s_construct_decomposition__block(decomposition%local%dimension(id)%block(i), nd)
    enddo
    ! copy index arrays
    do i=1,n_blocks
      decomposition%local%dimension(id)%block(i)%mn = decomposition%local%mn
      decomposition%local%dimension(id)%block(i)%mx = decomposition%local%mx
      decomposition%local%dimension(id)%block(i)%nw = decomposition%local%nw
    enddo
    ! modify index arrays, sub-divide the blocks along block_dim
    nw = decomposition%local%nw(block_dim) / n_blocks
    mn = decomposition%local%mn(block_dim)
    do i=1,n_blocks
      decomposition%local%dimension(id)%block(i)%nw(block_dim) = nw
      decomposition%local%dimension(id)%block(i)%mn(block_dim) = mn + (i-1) * nw
      decomposition%local%dimension(id)%block(i)%mx(block_dim) = mn + i * nw - 1
    end do
    decomposition%local%dimension(id)%block_dim = block_dim
  end subroutine sll_s_dd_define_blocks

  !> Allocate the \a i_block th halo buffers for dimension \a id
  subroutine sll_s_dd_allocate_buffer(decomposition, id, i_block, i_buffer, width, label)
    type(sll_t_decomposition), pointer :: decomposition !> decomposition object
    sll_int32, intent(in) :: id !> dimension for which to allocate the halo cells
    sll_int32, intent(in) :: i_block !>  ????
    sll_int32, intent(in) :: i_buffer !>  ????
    sll_int32, intent(in) ::  width !> width of the halo cells
    character(len=*), optional, intent(in) :: label
    sll_int32 :: nd, ierr
    class(sll_t_decomposition__buffer), pointer :: BUFALIAS

    nd = size(decomposition%global)
    SLL_ASSERT_ALWAYS( width >= 0 )
    SLL_ASSERT_ALWAYS( id <= nd )
    SLL_ASSERT_ALWAYS( i_buffer <= max_buffers_per_block )
    SLL_ASSERT_ALWAYS( associated(decomposition%local%dimension(id)%block) )
    SLL_ASSERT_ALWAYS( associated(decomposition%local%dimension(id)%block(i_block)%buffer) )

    BUFALIAS => decomposition%local%dimension(id)%block(i_block)%buffer(i_buffer)

    ! select type(ptr => BUFALIAS)
    !   class is (sll_t_decomposition__buffer_3d)
    !     SLL_ASSERT_ALWAYS(nd == 3)
    !   class is (sll_t_decomposition__buffer_6d)
    !     SLL_ASSERT_ALWAYS(nd == 6)
    !   class is (sll_t_decomposition__buffer)
    !     SLL_ERROR("sll_s_dd_allocate_buffer", "buffer is of generic type")
    !   class default
    !     SLL_ERROR("sll_s_dd_allocate_buffer", "unknown buffer data structure")
    ! end select

    ! initialize index arrays
    SLL_ALLOCATE(BUFALIAS%mn(nd), ierr)
    SLL_ALLOCATE(BUFALIAS%mx(nd), ierr)
    SLL_ALLOCATE(BUFALIAS%nw(nd), ierr)
    BUFALIAS%mn = decomposition%local%dimension(id)%block(i_block)%mn
    BUFALIAS%mx = decomposition%local%dimension(id)%block(i_block)%mx
    BUFALIAS%nw = decomposition%local%dimension(id)%block(i_block)%nw
    if (get_buffer_side(i_buffer) == LEFT_BUFFER) then
      ! --- odd indices (1,3,5, ...) shall designate "left" buffers
      BUFALIAS%mx(id) = BUFALIAS%mn(id) - 1
      BUFALIAS%mn(id) = BUFALIAS%mx(id) - width + 1
      BUFALIAS%nw(id) = width
    else ! RIGHT_BUFFER
      ! --- even indices (2,4,6, ...) shall designate "right" buffers
      BUFALIAS%mn(id) = BUFALIAS%mx(id) + 1
      BUFALIAS%mx(id) = BUFALIAS%mn(id) + width - 1
      BUFALIAS%nw(id) = width
    endif

!     select type(ptr => BUFALIAS)
!       class is (sll_t_decomposition__buffer_3d)
! #ifdef USE_FMEMPOOL
!         call mp_acquire(ptr%mem, ptr%mn, ptr%mx)
! #else
!         allocate(ptr%mem(ptr%mn(1):ptr%mx(1), &
!                          ptr%mn(2):ptr%mx(2), &
!                          ptr%mn(3):ptr%mx(3)))
! #endif
!       class is (sll_t_decomposition__buffer_6d)
! #ifdef USE_FMEMPOOL
!         call mp_acquire(ptr%mem, ptr%mn, ptr%mx)
! #else
!         allocate(ptr%mem(ptr%mn(1):ptr%mx(1), &
!                          ptr%mn(2):ptr%mx(2), &
!                          ptr%mn(3):ptr%mx(3), &
!                          ptr%mn(4):ptr%mx(4), &
!                          ptr%mn(5):ptr%mx(5), &
!                          ptr%mn(6):ptr%mx(6)))
! #endif
!       class default
!         SLL_ERROR("sll_s_dd_allocate_buffer", "unknown buffer data structure")
!     end select

#ifdef USE_FMEMPOOL
        call mp_acquire(BUFALIAS%mem, BUFALIAS%mn, BUFALIAS%mx)
#else
        allocate(BUFALIAS%mem(BUFALIAS%mn(1):BUFALIAS%mx(1), &
                              BUFALIAS%mn(2):BUFALIAS%mx(2), &
                              BUFALIAS%mn(3):BUFALIAS%mx(3), &
                              BUFALIAS%mn(4):BUFALIAS%mx(4), &
                              BUFALIAS%mn(5):BUFALIAS%mx(5), &
                              BUFALIAS%mn(6):BUFALIAS%mx(6)))
#endif

    BUFALIAS%valid = .true.

    nullify(BUFALIAS)
  end subroutine sll_s_dd_allocate_buffer


  subroutine sll_s_dd_deallocate_buffer(decomposition, id, i_block, i_buffer)
    type(sll_t_decomposition), pointer :: decomposition
    sll_int32, intent(in) :: id, i_block, i_buffer
    sll_int32 :: nd
    class(sll_t_decomposition__buffer), pointer :: buffer  ! convenience shortcut

    nd = size(decomposition%global)
    SLL_ASSERT_ALWAYS( id <= nd )
    SLL_ASSERT_ALWAYS( i_buffer <= max_buffers_per_block )
    SLL_ASSERT_ALWAYS( associated(decomposition%local%dimension(id)%block) )
    buffer => decomposition%local%dimension(id)%block(i_block)%buffer(i_buffer)
    call sll_s_destruct_decomposition__buffer(buffer)
    nullify(buffer)
  end subroutine sll_s_dd_deallocate_buffer


  subroutine sll_s_post_halo_exchange_real64(topology, decomposition, f6d, id, i_block, i_buffer, n_threads, post_halo_exchange)
    integer, parameter :: nd = 6
    type(sll_t_cartesian_topology_6d), pointer, intent(in) :: topology
    type(sll_t_decomposition), pointer, intent(inout) :: decomposition
    sll_real64, dimension(:,:,:,:,:,:), intent(inout) :: &
                                        f6d(decomposition%local%mn(1):decomposition%local%mx(1), &
                                            decomposition%local%mn(2):decomposition%local%mx(2), &
                                            decomposition%local%mn(3):decomposition%local%mx(3), &
                                            decomposition%local%mn(4):decomposition%local%mx(4), &
                                            decomposition%local%mn(5):decomposition%local%mx(5), &
                                            decomposition%local%mn(6):decomposition%local%mx(6))
    sll_int32, intent(in) :: id, i_block, i_buffer
    sll_int32, intent(in), optional :: n_threads
    logical, intent(in), optional :: post_halo_exchange

    class(sll_t_decomposition__buffer), pointer :: buffer  ! convenience shortcut
    HALO_DTYPE, pointer :: mem(:,:,:,:,:,:)  ! convenience shortcut
    integer :: jd, i, j, k, l, m, n
    integer :: ierr, n_omp_threads
    logical, save :: first_call = .true.
    logical :: do_post

    ! --- MPI communication-related variables and buffers
    sll_int64 :: nxc  ! total number of elements to be exchanged
    sll_int32 :: nel  ! number of elements to be exchanged at a single MPI call
    integer, dimension(:,:) :: r_rx(nd,2)  ! index ranges for copy operations
    integer, dimension(:,:) :: r_tx(nd,2)
    sll_int32, parameter :: nxc_max = 2147483647
#ifdef USE_HALO_REAL32
    integer, parameter :: mpi_precision = MPI_REAL
    integer, parameter :: word_size = 4
#else
    integer, parameter :: mpi_precision = MPI_DOUBLE_PRECISION
    integer, parameter :: word_size = 8
#endif
    integer :: mpi_tag, mpi_dest, mpi_source

    logical, save :: use_compression
    integer, save :: prec
    logical, save :: compression_verbose

#ifdef _OPENMP
    if (present(n_threads)) then
      n_omp_threads = n_threads
    else
      n_omp_threads = omp_get_max_threads()
    endif
#else
    n_omp_threads = 1
#endif

    if (present(post_halo_exchange)) then
      do_post = post_halo_exchange
    else
      do_post = .true.
    endif

    if (first_call) then
      first_call = .false.
      compression_verbose = .false.
#ifdef USE_HALO_REAL32
      if (topology%rank == 0) then
        write(*,*) "sll_m_decomposition::post_halo_exchange() uses single precision messages"
      endif
      use_compression = .false.
#else
#ifdef USE_ZFP
      use_compression = sll_f_query_environment("SLL_USE_COMPRESSION", .false.)
      if (use_compression) then
        prec = sll_f_query_environment("SLL_ZFP_PRECISION", 32)
        call set_compression_precision(prec)
        if (topology%rank == 0) then
          write(*,*) "sll_m_decomposition::post_halo_exchange() uses message compression"
          compression_verbose = sll_f_query_environment("SLL_COMPRESSION_VERBOSE", .false.)
        endif
      endif
#else
      use_compression = .false.
#endif
#endif
    endif

    SLL_ASSERT_ALWAYS( size(decomposition%global) == nd )
    SLL_ASSERT_ALWAYS( id <= nd )
    SLL_ASSERT_ALWAYS( i_buffer <= max_buffers_per_block )
    SLL_ASSERT_ALWAYS( associated(decomposition%local%dimension(id)%block) )
    buffer => decomposition%local%dimension(id)%block(i_block)%buffer(i_buffer)

    if (buffer%nw(id) <= 0) then
      return
    endif

    if (get_buffer_side(i_buffer) == LEFT_BUFFER) then
      ! --- fill the left buffer, equivalent to copying from f6d to the right neighbor
      mpi_dest = topology%neighbors(2*id)
      mpi_source = topology%neighbors(2*id-1)
      ! --- calculate rx and tx index ranges
      ! index ranges on the computation array to be sent
      r_tx(:,1)  = decomposition%local%dimension(id)%block(i_block)%mn(:)
      r_tx(:,2)  = decomposition%local%dimension(id)%block(i_block)%mx(:)
      r_tx(id,1) = decomposition%local%dimension(id)%block(i_block)%mx(id) - buffer%nw(id) + 1
      r_tx(id,2) = decomposition%local%dimension(id)%block(i_block)%mx(id)
    else ! RIGHT_BUFFER
      ! --- fill the right buffer, equivalent to copying from f6d to the left neighbor
      mpi_dest = topology%neighbors(2*id-1)
      mpi_source = topology%neighbors(2*id)
      ! --- calculate rx and tx index ranges
      ! index ranges on the computation array to be sent
      r_tx(:,1)  = decomposition%local%dimension(id)%block(i_block)%mn(:)
      r_tx(:,2)  = decomposition%local%dimension(id)%block(i_block)%mx(:)
      r_tx(id,1) = decomposition%local%dimension(id)%block(i_block)%mn(id)
      r_tx(id,2) = decomposition%local%dimension(id)%block(i_block)%mn(id) + buffer%nw(id) - 1
    endif
    ! index ranges on the buffer are just its extents
    r_rx(:,1)  = decomposition%local%dimension(id)%block(i_block)%buffer(i_buffer)%mn(:)
    r_rx(:,2)  = decomposition%local%dimension(id)%block(i_block)%buffer(i_buffer)%mx(:)

    mem => sll_f_get_mem_6d_from_buffer_obj(buffer)

    if (topology%procs(id) == 1) then
      ! --- we copy the ghost cells directly, assuming periodic BCs for the moment
!        buffer%mem(r_rx(1,1):r_rx(1,2), r_rx(2,1):r_rx(2,2), &
!                   r_rx(3,1):r_rx(3,2), r_rx(4,1):r_rx(4,2), &
!                   r_rx(5,1):r_rx(5,2), r_rx(6,1):r_rx(6,2)) = &
!        f6d(r_tx(1,1):r_tx(1,2), r_tx(2,1):r_tx(2,2), r_tx(3,1):r_tx(3,2), &
!            r_tx(4,1):r_tx(4,2), r_tx(5,1):r_tx(5,2), r_tx(6,1):r_tx(6,2))
!$omp parallel num_threads(n_omp_threads) default(shared) private(i,j,k,l,m,n)
!$omp do OMP_SCHEDULE OMP_COLLAPSE
      do n = 0, r_rx(6,2)-r_rx(6,1)
        do m = 0, r_rx(5,2)-r_rx(5,1)
          do l = 0, r_rx(4,2)-r_rx(4,1)
            do k = 0, r_rx(3,2)-r_rx(3,1)
              do j = 0, r_rx(2,2)-r_rx(2,1)
                do i = 0, r_rx(1,2)-r_rx(1,1)
                  mem(r_rx(1,1)+i, r_rx(2,1)+j, r_rx(3,1)+k, &
                      r_rx(4,1)+l, r_rx(5,1)+m, r_rx(6,1)+n) = &
                  f6d(r_tx(1,1)+i, r_tx(2,1)+j, r_tx(3,1)+k, &
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
      nxc = int(buffer%nw(id), i64)
      do jd=1,nd
        if (jd == id) then
          cycle
        else
          nxc = nxc * buffer%nw(jd)
        endif
      end do
      SLL_ALLOCATE(buffer%sendbuf(nxc), ierr)

      call sll_s_copy_array_to_buffer_6d_real64(f6d, &
                                          decomposition%local%mn, &
                                          decomposition%local%mx, &
                                          buffer%sendbuf, r_tx, &
                                          n_omp_threads)

      if (use_compression) then
#ifndef USE_HALO_REAL32
        if (modulo(nxc, int(zfp_blocksize, i64)) /= 0) then
          if (topology%rank == 0) then
            write(*,*) "mpi_sendrecv_compressed() : disabled due to blocksize mismatch"
          endif
          buffer%use_compression = .false.
        else
          buffer%use_compression = .true.
          buffer%compression_verbose = compression_verbose
          i = int(nxc, i32)  ! i64 -> i32
          call deflate_buffer_real64(buffer%sendbuf, buffer%comp_sendbuf, n_doubles=i, &
                                     n_threads=n_omp_threads)
          SLL_DEALLOCATE(buffer%sendbuf, ierr)
        endif
#endif
      endif

      SLL_ASSERT_ALWAYS(nxc <= nxc_max)  ! check if message size is within the allowed limit
      nel = int(nxc, i32)  ! 64 bit -> 32 bit
      mpi_tag = 1000*id + 10*i_block + i_buffer

!      call MPI_Sendrecv(buffer%sendbuf, nel, mpi_precision, mpi_dest, mpi_tag, &
!                        mem, nel, mpi_precision, mpi_source, mpi_tag, &
!                        topology%comm, MPI_STATUS_IGNORE, ierr)
!      SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)

      if (do_post) then
        call MPI_Isend(buffer%sendbuf, nel, mpi_precision, mpi_dest,   mpi_tag, topology%comm, buffer%request(1), ierr)
        SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
        call MPI_Irecv(mem,            nel, mpi_precision, mpi_source, mpi_tag, topology%comm, buffer%request(2), ierr)
        SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
      else
        buffer%nel = nel
        buffer%mpi_source = mpi_source
        buffer%mpi_dest = mpi_dest
        buffer%mpi_tag = mpi_tag
        buffer%mpi_precision = mpi_precision
      endif

    endif
  end subroutine sll_s_post_halo_exchange_real64


  subroutine sll_s_wait_halo_exchange(topology, decomposition, id, i_block, i_buffer)
    integer, parameter :: nd = 6
    type(sll_t_cartesian_topology_6d), pointer, intent(in) :: topology
    type(sll_t_decomposition), pointer, intent(inout) :: decomposition
    sll_int32, intent(in) :: id, i_block, i_buffer
    class(sll_t_decomposition__buffer), pointer :: buffer  ! convenience shortcut
    integer :: ierr
    SLL_ASSERT_ALWAYS( size(decomposition%global) == nd )
    SLL_ASSERT_ALWAYS( id <= nd )
    SLL_ASSERT_ALWAYS( i_buffer <= max_buffers_per_block )
    SLL_ASSERT_ALWAYS( associated(decomposition%local%dimension(id)%block) )
    buffer => decomposition%local%dimension(id)%block(i_block)%buffer(i_buffer)
    if (buffer%nw(id) <= 0) then
      return
    endif
    if (topology%procs(id) == 1) then
      ! --- do nothing ---
    else
      call MPI_Waitall(2, buffer%request, MPI_STATUS_IGNORE, ierr)
      SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
      SLL_DEALLOCATE(buffer%sendbuf, ierr)
    endif
  end subroutine sll_s_wait_halo_exchange


  subroutine sll_s_blocking_halo_exchange(topology, decomposition, id, i_block, i_buffer)
    integer, parameter :: nd = 6
    type(sll_t_cartesian_topology_6d), pointer, intent(in) :: topology
    type(sll_t_decomposition), pointer, intent(inout) :: decomposition
    sll_int32, intent(in) :: id, i_block, i_buffer
    class(sll_t_decomposition__buffer), pointer :: buffer  ! convenience shortcut
    HALO_DTYPE, pointer :: mem(:,:,:,:,:,:)
    integer :: ierr
    SLL_ASSERT_ALWAYS( size(decomposition%global) == nd )
    SLL_ASSERT_ALWAYS( id <= nd )
    SLL_ASSERT_ALWAYS( i_buffer <= max_buffers_per_block )
    SLL_ASSERT_ALWAYS( associated(decomposition%local%dimension(id)%block) )
    buffer => decomposition%local%dimension(id)%block(i_block)%buffer(i_buffer)
    if (buffer%nw(id) <= 0) then
      return
    else
      if (topology%procs(id) == 1) then
        ! --- do nothing ---
      else
        if (buffer%use_compression) then
          call sll_s_mpi_sendrecv_compressed_core(buffer%comp_sendbuf, buffer%comp_recvbuf,&
                                            buffer%mpi_dest, buffer%mpi_source,&
                                            topology%comm, buffer%compression_verbose,&
                                            buffer%mpi_tag)
          call deallocate_compressed_buffer_obj(buffer%comp_sendbuf)
        else
          SLL_ASSERT_ALWAYS(associated(buffer%sendbuf))
          mem => sll_f_get_mem_6d_from_buffer_obj(buffer)
          call MPI_Sendrecv(buffer%sendbuf, buffer%nel, buffer%mpi_precision, buffer%mpi_dest,   buffer%mpi_tag, &
                            mem,            buffer%nel, buffer%mpi_precision, buffer%mpi_source, buffer%mpi_tag, &
                            topology%comm, MPI_STATUS_IGNORE, ierr)
          SLL_ASSERT_ALWAYS(ierr == MPI_SUCCESS)
          SLL_DEALLOCATE(buffer%sendbuf, ierr)
        endif
      endif
    endif
  end subroutine sll_s_blocking_halo_exchange


  function get_buffer_side(i_buffer)
    sll_int32 :: get_buffer_side, i_buffer
    get_buffer_side = mod(i_buffer, 2)
  end function get_buffer_side


  function sll_f_get_mem_6d_from_buffer_obj(self)
    HALO_DTYPE, pointer :: sll_f_get_mem_6d_from_buffer_obj(:,:,:,:,:,:)
    class(sll_t_decomposition__buffer) :: self
    VERBOSE("sll_f_get_mem_6d_from_buffer_obj")
    nullify(sll_f_get_mem_6d_from_buffer_obj)
    if (self%valid) then
      ! select type(ptr => self)
        ! class is (sll_t_decomposition__buffer_6d)
          sll_f_get_mem_6d_from_buffer_obj => self%mem
        ! class default
          ! SLL_ERROR("sll_s_get_ptr_6d_from_buffer", "buffer does not contain 6d array")
      ! end select
    else
      SLL_ERROR("sll_s_get_ptr_6d_from_buffer", "attempt to access invalid buffer")
    endif
  end function sll_f_get_mem_6d_from_buffer_obj


  function sll_f_dd_get_n_blocks(decomposition, id)
    sll_int32 :: sll_f_dd_get_n_blocks
    type(sll_t_decomposition), pointer :: decomposition
    sll_int32, intent(in) :: id
    sll_int32 :: nd, ierr
    nd = size(decomposition%global)
    SLL_ASSERT_ALWAYS( id <= nd )
    SLL_ASSERT_ALWAYS( associated(decomposition%local%dimension) )
    if (associated(decomposition%local%dimension(id)%block)) then
      sll_f_dd_get_n_blocks = size(decomposition%local%dimension(id)%block)
    else
      sll_f_dd_get_n_blocks = 0
    endif
  end function sll_f_dd_get_n_blocks


end module sll_m_decomposition_advanced
