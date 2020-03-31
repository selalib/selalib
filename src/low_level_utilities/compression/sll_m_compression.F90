!**************************************************************
!> @brief
!> Module providing an F90 interface to the ZFP compression library:
!> http://computation.llnl.gov/projects/floating-point-compression
!> In addition it provides simple threaded (de)compression routines.
!> Important: This module uses C-like 0-based indexing!
!> @author
!> Klaus Reuter, Max Planck Computing and Data Facility (MPCDF)
!**************************************************************


module sll_m_compression
#include "sll_working_precision.h"

  use iso_c_binding

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

#ifdef USE_ZFP
  interface
    integer (c_int) function sll_f_compressbound_zfp(buf_in, n_doubles_in, prec) &
      bind(c, name='zfp_maximum_size')
      use iso_c_binding
      implicit none
      type (c_ptr), value :: buf_in
      integer (c_int), value :: n_doubles_in
      integer (c_int), value :: prec
    end function sll_f_compressbound_zfp

    integer (c_int) function sll_f_deflate_zfp(buf_in, buf_out, n_doubles_in, n_bytes_out_max, prec) &
      bind(c, name='zfp_compress_safe')
      use iso_c_binding
      implicit none
      type (c_ptr), value :: buf_in
      type (c_ptr), value :: buf_out
      integer (c_int), value :: n_doubles_in
      integer (c_int), value :: n_bytes_out_max
      integer (c_int), value :: prec
    end function sll_f_deflate_zfp

    integer (c_int) function sll_f_inflate_zfp(buf_in, buf_out, n_bytes_in, n_doubles_out_max, prec) &
      bind(c, name='zfp_decompress_safe')
      use iso_c_binding
      implicit none
      type (c_ptr), value :: buf_in
      type (c_ptr), value :: buf_out
      integer (c_int), value :: n_bytes_in
      integer (c_int), value :: n_doubles_out_max
      integer (c_int), value :: prec
    end function sll_f_inflate_zfp
  end interface
#else
    ! add stub code here
#endif

  ! ZFP compresses in blocks of 64 double precision numbers
  integer, parameter :: zfp_blocksize = 64
  integer :: zfp_precision = 32

  !> data structure to support threaded ZFP compression and decompression
  type :: sll_t_compressed_buffer
    character, dimension(:), pointer :: buffer => null()          ! concatenated slices of compressed data
    ! NOTE -- integers below need to be taken care of during MPI communication (concatenate/decatenate)
    integer :: n_slices                                           ! number of data slices (== threads used during compression)
    integer :: n_bytes_deflated_total = 0                         ! total size of compressed data
    integer :: n_bytes_inflated_total = 0                         ! total size of decompressed data
    integer, dimension(:), pointer :: n_bytes_deflated => null()  ! length of the i-th slice of compressed data
    integer, dimension(:), pointer :: offset_deflated => null()   ! offset for the i-th slice of compressed data
    integer, dimension(:), pointer :: n_bytes_inflated => null()  ! length of the i-th slice of uncompressed data
    integer, dimension(:), pointer :: offset_inflated => null()   ! offset for the i-th slice of uncompressed data
  end type sll_t_compressed_buffer

!  public :: sll_f_compressbound_zfp, &
!            sll_f_deflate_zfp, &
!            sll_f_inflate_zfp, &
!            zfp_blocksize, &
!            zfp_accuracy


contains


  !> compress buffer
  subroutine deflate_buffer_real64(buf, comp, n_doubles, n_threads)
    sll_real64, intent(in), target :: buf(0:)  ! buffer to be compressed
    type(sll_t_compressed_buffer), intent(inout) :: comp  ! data structure containing compressed data and offsets
    integer, intent(in), optional :: n_doubles, n_threads  ! number of doubles to be compressed, number of threads to be used
#ifdef USE_ZFP
    integer :: n_omp_threads, omp_size, omp_rank
    integer :: n_total, n_per_thread  ! total and per-thread number of elements in 'buf'
    integer :: n_offset_thread  ! offset to indicate where a thread shall start compressing
    integer :: n_bytes_comp_thread_max  ! maximum possible size for compressed data (needed for allocation)
    integer :: i, off
    character, allocatable, target :: thread_buf(:)
    integer, parameter :: word_size = 8  ! double precision

    if (present(n_doubles)) then
      n_total = n_doubles
    else
      n_total = size(buf)
    endif

#ifdef _OPENMP
    if (present(n_threads)) then
      n_omp_threads = n_threads
    else
      n_omp_threads = omp_get_max_threads()
    endif
    ! n_omp_threads = 1
#else
    n_omp_threads = 1
#endif

    call deallocate_compressed_buffer_obj(comp)
    call allocate_compressed_buffer_index_arrays(comp, n_omp_threads)

!$omp parallel num_threads(n_omp_threads) default(shared) &
!$omp& private(n_per_thread, n_offset_thread, n_bytes_comp_thread_max, thread_buf) &
!$omp& private(omp_size, omp_rank, i, off)
#ifdef _OPENMP
    omp_size = omp_get_num_threads()
    omp_rank = omp_get_thread_num()
#else
    omp_size = 1
    omp_rank = 0
#endif

    ! distribute the input data among the threads, taking care not to violate the ZFP blocksize
    n_per_thread = n_total / zfp_blocksize / omp_size
    n_offset_thread = n_per_thread * omp_rank
    if (omp_rank == omp_size-1) then
      n_per_thread = (n_total / zfp_blocksize) - n_per_thread * omp_rank
    endif
    n_per_thread = n_per_thread * zfp_blocksize
    n_offset_thread = n_offset_thread * zfp_blocksize

    ! add size and offset information to 'comp' data structure
    comp%n_bytes_inflated(omp_rank) = word_size * n_per_thread
    comp%offset_inflated(omp_rank) = word_size * n_offset_thread

    n_bytes_comp_thread_max = sll_f_compressbound_zfp(c_loc(buf(n_offset_thread)), n_per_thread, zfp_precision)
    allocate(thread_buf(0:n_bytes_comp_thread_max-1))

    comp%n_bytes_deflated(omp_rank) = sll_f_deflate_zfp(c_loc(buf(n_offset_thread)), c_loc(thread_buf), &
                                                        n_per_thread, n_bytes_comp_thread_max, zfp_precision)

!$omp barrier
!$omp master
    comp%n_bytes_deflated_total = sum(comp%n_bytes_deflated)
    allocate(comp%buffer(0:comp%n_bytes_deflated_total-1))
    comp%offset_deflated(0) = 0
    do i=1, omp_size-1
      comp%offset_deflated(i) = comp%offset_deflated(i-1) + comp%n_bytes_deflated(i-1)
    enddo
!$omp end master
!$omp barrier

    ! cache thread's offset for the use in the subsequent loop
    off = comp%offset_deflated(omp_rank)
    ! copy the thread buffers in parallel into the buffer
    do i=0, comp%n_bytes_deflated(omp_rank)-1
      comp%buffer(off + i) = thread_buf(i)
    enddo

    deallocate(thread_buf)
!$omp end parallel

    comp%n_bytes_inflated_total = word_size * n_total
#else
    ! add stub code here
#endif
  end subroutine deflate_buffer_real64


  !> decompress buffer
  subroutine inflate_buffer_real64(buf, comp, n_threads)
    sll_real64, intent(inout), target :: buf(0:)  ! 1d view on buffer to be decompressed into
    type(sll_t_compressed_buffer), intent(inout) :: comp  ! data structure containing compressed data and offsets
    integer, intent(in), optional :: n_threads
#ifdef USE_ZFP
    integer :: n_omp_threads, omp_size, omp_rank
    integer :: i, off, n_el, ierr
    integer, parameter :: word_size = 8  ! double precision

#ifdef _OPENMP
    if (present(n_threads)) then
      n_omp_threads = n_threads
    else
      n_omp_threads = omp_get_max_threads()
    endif
    ! --- DEBUG ---
    ! use as many threads as were used during compression
    ! n_omp_threads = comp%n_slices
    ! n_omp_threads = 1
#else
    n_omp_threads = 1
#endif

!$omp parallel num_threads(n_omp_threads) default(shared) &
!$omp& private(omp_size, omp_rank, i, ierr, off, n_el)
#ifdef _OPENMP
    omp_size = omp_get_num_threads()
    omp_rank = omp_get_thread_num()
#else
    omp_size = 1
    omp_rank = 0
#endif

!$omp do schedule(static,1)
    do i=0,comp%n_slices-1
      ! offset and element count must be in units of double precision numbers (not raw bytes)
      off = comp%offset_inflated(i)/word_size
      n_el = comp%n_bytes_inflated(i)/word_size
      ierr = sll_f_inflate_zfp(c_loc(comp%buffer(comp%offset_deflated(i))), c_loc(buf(off)), &
                               comp%n_bytes_deflated(i), n_el, zfp_precision)
    enddo
!$omp end do

!$omp end parallel

#else
    ! add stub code here
#endif
  end subroutine inflate_buffer_real64


  subroutine deallocate_compressed_buffer_obj(comp)
    type(sll_t_compressed_buffer), intent(inout) :: comp  ! data structure containing compressed data and offsets
    if (associated(comp%buffer)) then
      deallocate(comp%buffer); nullify(comp%buffer)
    endif
    if (associated(comp%n_bytes_deflated)) then
      deallocate(comp%n_bytes_deflated); nullify(comp%n_bytes_deflated)
    endif
    if (associated(comp%offset_deflated)) then
      deallocate(comp%offset_deflated); nullify(comp%offset_deflated)
    endif
    if (associated(comp%n_bytes_inflated)) then
      deallocate(comp%n_bytes_inflated); nullify(comp%n_bytes_inflated)
    endif
    if (associated(comp%offset_inflated)) then
      deallocate(comp%offset_inflated); nullify(comp%offset_inflated)
    endif
    comp%n_slices = 0
    comp%n_bytes_inflated_total = 0
    comp%n_bytes_deflated_total = 0
  end subroutine deallocate_compressed_buffer_obj


  subroutine allocate_compressed_buffer_index_arrays(comp, n_slices)
    type(sll_t_compressed_buffer), intent(inout) :: comp
    integer, intent(in) :: n_slices
    allocate(comp%n_bytes_deflated(0:n_slices-1))
    allocate(comp%offset_deflated(0:n_slices-1))
    allocate(comp%n_bytes_inflated(0:n_slices-1))
    allocate(comp%offset_inflated(0:n_slices-1))
    comp%n_slices = n_slices
  end subroutine allocate_compressed_buffer_index_arrays


  subroutine print_compression_information(comp, verbose)
    type(sll_t_compressed_buffer), intent(in) :: comp
    logical, intent(in), optional :: verbose

    write(*,*) "--- sll_m_compression ---------------------"
    if ((present(verbose)).and.(verbose)) then
      write(*,*) "n_slices               =", comp%n_slices
      !write(*,*) "n_bytes_inflated_total =", comp%n_bytes_inflated_total
      write(*,*) "n_bytes_inflated       =", comp%n_bytes_inflated
      write(*,*) "offset_inflated        =", comp%offset_inflated
      write(*,*) "n_bytes_deflated_total =", comp%n_bytes_deflated_total
      write(*,*) "n_bytes_deflated       =", comp%n_bytes_deflated
      write(*,*) "offset_deflated        =", comp%offset_deflated
      write(*,*) "zfp_precision          =", zfp_precision
    endif
    if (comp%n_bytes_deflated_total > 0) then
      write(*,*) "n_bytes_inflated_total =", comp%n_bytes_inflated_total
      write(*,*) "ratio                  =", &
        real(comp%n_bytes_inflated_total)/real(comp%n_bytes_deflated_total)
    endif
    write(*,*) "-------------------------------------------"
  end subroutine print_compression_information


  subroutine set_compression_precision(prec)
    integer, intent(in) :: prec
    zfp_precision = prec
  end subroutine set_compression_precision


  !> allocate array, copy indices from comp into array, return
  function concatenate_index_arrays(comp, array) result(n_el)
    type(sll_t_compressed_buffer) :: comp
    integer, pointer :: array(:)
    integer :: i, n_el, idx

    n_el = 3 + 4 * comp%n_slices
    allocate(array(0:n_el-1))

    idx = 0
    array(idx) = comp%n_slices
    idx = idx + 1
    array(idx) = comp%n_bytes_deflated_total
    idx = idx + 1
    array(idx) = comp%n_bytes_inflated_total
    idx = idx + 1
    do i=0,comp%n_slices-1
      array(idx) = comp%n_bytes_deflated(i)
      idx = idx + 1
    enddo
    do i=0,comp%n_slices-1
      array(idx) = comp%offset_deflated(i)
      idx = idx + 1
    enddo
    do i=0,comp%n_slices-1
      array(idx) = comp%n_bytes_inflated(i)
      idx = idx + 1
    enddo
    do i=0,comp%n_slices-1
      array(idx) = comp%offset_inflated(i)
      idx = idx + 1
    enddo
  end function concatenate_index_arrays


  subroutine decatenate_index_arrays(comp, array)
    type(sll_t_compressed_buffer) :: comp
    integer, pointer :: array(:)
    integer :: i, n_el, idx
    idx = 0
    comp%n_slices = array(idx)
    idx = idx + 1
    comp%n_bytes_deflated_total = array(idx)
    idx = idx + 1
    comp%n_bytes_inflated_total = array(idx)
    idx = idx + 1
    call allocate_compressed_buffer_index_arrays(comp, comp%n_slices)
    do i=0,comp%n_slices-1
      comp%n_bytes_deflated(i) = array(idx)
      idx = idx + 1
    enddo
    do i=0,comp%n_slices-1
      comp%offset_deflated(i) = array(idx)
      idx = idx + 1
    enddo
    do i=0,comp%n_slices-1
      comp%n_bytes_inflated(i) = array(idx)
      idx = idx + 1
    enddo
    do i=0,comp%n_slices-1
      comp%offset_inflated(i) = array(idx)
      idx = idx + 1
    enddo
  end subroutine decatenate_index_arrays


  ! simple time function for the unit test code to be independent from the rest of selalib
  function get_time()
    double precision :: get_time
#ifdef _OPENMP
    get_time = omp_get_wtime()
#else
    get_time = 0.0
#endif
  end function

end module sll_m_compression
