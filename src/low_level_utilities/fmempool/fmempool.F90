
#include "fassert.inc"

!> @brief Plain Fortran implementation of a memory pool.
module fmempool
   use iso_c_binding
#ifdef _OPENMP
   use omp_lib
#endif
   use assert
   implicit none

   interface mp_acquire
      module procedure mp_acquire_int_1d
      module procedure mp_acquire_int_2d
      module procedure mp_acquire_int_3d
      module procedure mp_acquire_int_6d
      module procedure mp_acquire_real_1d
      module procedure mp_acquire_real_2d
      module procedure mp_acquire_real_3d
      module procedure mp_acquire_real_6d
      module procedure mp_acquire_double_1d
      module procedure mp_acquire_double_2d
      module procedure mp_acquire_double_3d
      module procedure mp_acquire_double_6d
   end interface mp_acquire

   interface mp_release
      module procedure mp_release_int_1d
      module procedure mp_release_int_2d
      module procedure mp_release_int_3d
      module procedure mp_release_int_6d
      module procedure mp_release_real_1d
      module procedure mp_release_real_2d
      module procedure mp_release_real_3d
      module procedure mp_release_real_6d
      module procedure mp_release_double_1d
      module procedure mp_release_double_2d
      module procedure mp_release_double_3d
      module procedure mp_release_double_6d
   end interface mp_release

   integer, parameter :: int_size = 4
   integer, parameter :: real_size = 4
   integer, parameter :: double_size = 2*real_size
   integer, parameter :: n_max_slices = 128

   type :: memslice
      logical :: acquired
      integer(kind=1), pointer :: mem(:)
   end type memslice

   type :: mempool
      type(memslice), pointer :: slice(:)
   end type mempool

   type(mempool), pointer, save :: pool(:) => null()
   logical, save :: verbose = .false.
   logical, save :: disabled = .false.

   public :: mp_init, mp_finalize, mp_statistics, &
             mp_compactify, mp_acquire, mp_release, &
             mp_disabled

   private

contains

   ! check if an environment variable is set to true or false
   function query_environment(env_var, default_val)
      implicit none
      logical :: query_environment
      character(len=*), intent(in) :: env_var
      logical, intent(in) :: default_val
      character(len=255) :: env_str
      query_environment = default_val
      call get_environment_variable(env_var, env_str)
      if (len_trim(env_str) > 0) then
         select case (trim(env_str))
         case ("1", "ON", "TRUE", "on", "true")
            query_environment = .true.
         case ("0", "OFF", "FALSE", "off", "false")
            query_environment = .false.
         end select
      end if
   end function query_environment

   ! get the number of elements from multidimensional array specifiers
   function get_n_elem(mn, mx) result(n_elem)
      integer, intent(in) :: mn(:)
      integer, intent(in) :: mx(:)
      integer :: i, nd
      integer(kind=8) :: n_elem
      nd = size(mn)
      i = 1
      n_elem = mx(i) - mn(i) + 1
      do i = 2, nd
         n_elem = n_elem*(mx(i) - mn(i) + 1)
      end do
   end function get_n_elem

   ! get the maximum available number of threads
   function get_omp_world_size()
      integer :: get_omp_world_size
#ifdef _OPENMP
      get_omp_world_size = omp_get_max_threads(); 
#else
      get_omp_world_size = 1; 
#endif
   end function get_omp_world_size

   ! get the rank of the current thread
   function get_omp_thread_idx()
      integer :: get_omp_thread_idx
#ifdef _OPENMP
      get_omp_thread_idx = omp_get_thread_num(); 
#else
      get_omp_thread_idx = 0; 
#endif
   end function get_omp_thread_idx

   subroutine mp_acquire_bytes(f_pointer, n_bytes)
      integer(kind=1), pointer :: f_pointer(:)
      integer(kind=8), intent(in) :: n_bytes
      integer :: it, j
      it = get_omp_thread_idx()
      do j = 1, n_max_slices
         if (.not. pool(it)%slice(j)%acquired) then
            exit
         end if
      end do
      ASSERT(j <= n_max_slices)  ! not enough slices available
      if (associated(pool(it)%slice(j)%mem)) then
         if (size(pool(it)%slice(j)%mem) < n_bytes) then
            deallocate (pool(it)%slice(j)%mem)
            allocate (pool(it)%slice(j)%mem(n_bytes))
         end if
      else
         allocate (pool(it)%slice(j)%mem(n_bytes))
      end if
      pool(it)%slice(j)%acquired = .true.
      f_pointer => pool(it)%slice(j)%mem
   end subroutine mp_acquire_bytes

   subroutine mp_release_bytes(c_pointer)
      type(c_ptr) :: c_pointer
      integer :: it, j
      it = get_omp_thread_idx()
      do j = 1, n_max_slices
         if (c_associated(c_loc(pool(it)%slice(j)%mem), c_pointer)) then
            exit
         end if
      end do
      ASSERT(j <= n_max_slices)  ! could not find pointer to be released
      pool(it)%slice(j)%acquired = .false.
      if (disabled) then
         deallocate (pool(it)%slice(j)%mem)
         nullify (pool(it)%slice(j)%mem)
      end if
   end subroutine mp_release_bytes

! --- public API functions below ---

   ! initialize the memory pool
   subroutine mp_init(min_threads, verbosity, disable)
      integer, intent(in), optional :: min_threads
      logical, intent(in), optional :: verbosity
      logical, intent(in), optional :: disable
      integer :: i, j, nt
      character(len=32) :: disabled_str
      if (present(min_threads)) then
         if (min_threads < get_omp_world_size()) then
            nt = get_omp_world_size()
         else
            nt = min_threads
         end if
      else
         nt = get_omp_world_size()
      end if
      allocate (pool(0:nt - 1))
      do i = 0, nt - 1
         allocate (pool(i)%slice(n_max_slices))
         do j = 1, n_max_slices
            pool(i)%slice(j)%acquired = .false.
            nullify (pool(i)%slice(j)%mem)
         end do
      end do
      if (present(verbosity)) then
         verbose = verbosity
      else
         verbose = query_environment("MP_VERBOSE", .false.)
      end if
      if (present(disable)) then
         disabled = disable
      else
         disabled = query_environment("MP_DISABLE", .false.)
      end if
      if (disabled) then
         disabled_str = ", disabled"
      else
         disabled_str = ""
      end if
      if (verbose) then
#ifdef _OPENMP
         write (*, *) "mempool: initialized (OpenMP)"//trim(adjustl(disabled_str))
#else
         write (*, *) "mempool: initialized (not threaded)"//trim(adjustl(disabled_str))
#endif
      end if
   end subroutine mp_init

   ! deallocate any memory used by the memory pool
   subroutine mp_finalize()
      integer :: i, j, nt
      nt = get_omp_world_size()
      do i = 0, nt - 1
         do j = 1, n_max_slices
            if (associated(pool(i)%slice(j)%mem)) then
               deallocate (pool(i)%slice(j)%mem)
               nullify (pool(i)%slice(j)%mem)
            end if
            pool(i)%slice(j)%acquired = .false.
         end do
         deallocate (pool(i)%slice)
         nullify (pool(i)%slice)
      end do
      deallocate (pool)
      nullify (pool)
      if (verbose) then
         write (*, *) "mempool: finalized"
      end if
   end subroutine mp_finalize

   ! deallocate any unused memory from the memory pool
   subroutine mp_compactify()
      integer :: i, j, it, nt
      it = get_omp_thread_idx()
      nt = get_omp_world_size()
      if (it == 0) then
         do i = 0, nt - 1
            do j = 1, n_max_slices
               if ((.not. pool(i)%slice(j)%acquired) .and. (associated(pool(i)%slice(j)%mem))) then
                  deallocate (pool(i)%slice(j)%mem)
                  nullify (pool(i)%slice(j)%mem)
               end if
            end do
         end do
      end if
   end subroutine mp_compactify

   ! print statistics
   subroutine mp_statistics()
      integer :: i, j, it, nt, n_slices, n_acquired, n_allocated
      integer(kind=8) :: n_bytes
      character(len=32) :: thread_str, acquired_str, slices_str, bytes_str, alloc_str
      it = get_omp_thread_idx()
      nt = get_omp_world_size()
      if ((verbose) .and. (it == 0)) then
         write (*, *) "mempool statistics"
         do i = 0, nt - 1
            n_bytes = 0
            n_acquired = 0
            n_allocated = 0
            do j = 1, n_max_slices
               if (pool(i)%slice(j)%acquired) then
                  n_acquired = n_acquired + 1
               end if
               if (associated(pool(i)%slice(j)%mem)) then
                  n_allocated = n_allocated + 1
                  n_bytes = n_bytes + size(pool(i)%slice(j)%mem)
               end if
            end do
            ! --- do some nice printing ---
            write (thread_str, *) i
            write (slices_str, *) n_max_slices
            write (acquired_str, *) n_acquired
            write (bytes_str, *) n_bytes
            write (alloc_str, *) n_allocated
            write (*, '(A)') "    "// &
               "pool["//trim(adjustl(thread_str))//"]: "// &
               "n_slices="//trim(adjustl(slices_str))//", "// &
               "n_acquired="//trim(adjustl(acquired_str))//", "// &
               "n_allocated="//trim(adjustl(alloc_str))//", "// &
               "n_bytes="//trim(adjustl(bytes_str))
         end do
      end if
   end subroutine mp_statistics

   ! return if the memory pool is disabled
   function mp_disabled()
      logical :: mp_disabled
      mp_disabled = disabled
   end function mp_disabled

   ! --- mp_acquire and mp_release interface implementations below ---

   ! double precision routines

   subroutine mp_acquire_double_1d(f_pointer, mn, mx)
      double precision, pointer :: f_pointer(:)
      integer, intent(in) :: mn(1)
      integer, intent(in) :: mx(1)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      double precision, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*double_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_double_1d

   subroutine mp_release_double_1d(f_pointer)
      double precision, pointer :: f_pointer(:)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_double_1d

   subroutine mp_acquire_double_2d(f_pointer, mn, mx)
      double precision, pointer :: f_pointer(:, :)
      integer, intent(in) :: mn(2)
      integer, intent(in) :: mx(2)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      double precision, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*double_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1), &
                mn(2):mx(2)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_double_2d

   subroutine mp_release_double_2d(f_pointer)
      double precision, pointer :: f_pointer(:, :)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_double_2d

   subroutine mp_acquire_double_3d(f_pointer, mn, mx)
      double precision, pointer :: f_pointer(:, :, :)
      integer, intent(in) :: mn(3)
      integer, intent(in) :: mx(3)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      double precision, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*double_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1), &
                mn(2):mx(2), &
                mn(3):mx(3)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_double_3d

   subroutine mp_release_double_3d(f_pointer)
      double precision, pointer :: f_pointer(:, :, :)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_double_3d

   subroutine mp_acquire_double_6d(f_pointer, mn, mx)
      double precision, pointer :: f_pointer(:, :, :, :, :, :)
      integer, intent(in) :: mn(6)
      integer, intent(in) :: mx(6)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      double precision, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*double_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1), &
                mn(2):mx(2), &
                mn(3):mx(3), &
                mn(4):mx(4), &
                mn(5):mx(5), &
                mn(6):mx(6)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_double_6d

   subroutine mp_release_double_6d(f_pointer)
      double precision, pointer :: f_pointer(:, :, :, :, :, :)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_double_6d

   ! single precision routines

   subroutine mp_acquire_real_1d(f_pointer, mn, mx)
      real, pointer :: f_pointer(:)
      integer, intent(in) :: mn(1)
      integer, intent(in) :: mx(1)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      real, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*real_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_real_1d

   subroutine mp_release_real_1d(f_pointer)
      real, pointer :: f_pointer(:)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_real_1d

   subroutine mp_acquire_real_2d(f_pointer, mn, mx)
      real, pointer :: f_pointer(:, :)
      integer, intent(in) :: mn(2)
      integer, intent(in) :: mx(2)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      real, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*real_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1), &
                mn(2):mx(2)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_real_2d

   subroutine mp_release_real_2d(f_pointer)
      real, pointer :: f_pointer(:, :)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_real_2d

   subroutine mp_acquire_real_3d(f_pointer, mn, mx)
      real, pointer :: f_pointer(:, :, :)
      integer, intent(in) :: mn(3)
      integer, intent(in) :: mx(3)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      real, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*real_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1), &
                mn(2):mx(2), &
                mn(3):mx(3)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_real_3d

   subroutine mp_release_real_3d(f_pointer)
      real, pointer :: f_pointer(:, :, :)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_real_3d

   subroutine mp_acquire_real_6d(f_pointer, mn, mx)
      real, pointer :: f_pointer(:, :, :, :, :, :)
      integer, intent(in) :: mn(6)
      integer, intent(in) :: mx(6)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      real, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*real_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1), &
                mn(2):mx(2), &
                mn(3):mx(3), &
                mn(4):mx(4), &
                mn(5):mx(5), &
                mn(6):mx(6)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_real_6d

   subroutine mp_release_real_6d(f_pointer)
      real, pointer :: f_pointer(:, :, :, :, :, :)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_real_6d

   ! 4-byte integer routines

   subroutine mp_acquire_int_1d(f_pointer, mn, mx)
      integer, pointer :: f_pointer(:)
      integer, intent(in) :: mn(1)
      integer, intent(in) :: mx(1)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      integer, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*int_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_int_1d

   subroutine mp_release_int_1d(f_pointer)
      integer, pointer :: f_pointer(:)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_int_1d

   subroutine mp_acquire_int_2d(f_pointer, mn, mx)
      integer, pointer :: f_pointer(:, :)
      integer, intent(in) :: mn(2)
      integer, intent(in) :: mx(2)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      integer, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*int_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1), &
                mn(2):mx(2)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_int_2d

   subroutine mp_release_int_2d(f_pointer)
      integer, pointer :: f_pointer(:, :)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_int_2d

   subroutine mp_acquire_int_3d(f_pointer, mn, mx)
      integer, pointer :: f_pointer(:, :, :)
      integer, intent(in) :: mn(3)
      integer, intent(in) :: mx(3)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      integer, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*int_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1), &
                mn(2):mx(2), &
                mn(3):mx(3)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_int_3d

   subroutine mp_release_int_3d(f_pointer)
      integer, pointer :: f_pointer(:, :, :)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_int_3d

   subroutine mp_acquire_int_6d(f_pointer, mn, mx)
      integer, pointer :: f_pointer(:, :, :, :, :, :)
      integer, intent(in) :: mn(6)
      integer, intent(in) :: mx(6)
      integer(kind=8) :: n_bytes, n_elem
      integer(kind=1), pointer :: ptr(:)
      integer, pointer :: flat_ptr(:)

      n_elem = get_n_elem(mn, mx)
      n_bytes = n_elem*int_size
      call mp_acquire_bytes(ptr, n_bytes)

      call c_f_pointer(c_loc(ptr), flat_ptr, [n_elem])
      f_pointer(mn(1):mx(1), &
                mn(2):mx(2), &
                mn(3):mx(3), &
                mn(4):mx(4), &
                mn(5):mx(5), &
                mn(6):mx(6)) => flat_ptr
      nullify (flat_ptr)
   end subroutine mp_acquire_int_6d

   subroutine mp_release_int_6d(f_pointer)
      integer, pointer :: f_pointer(:, :, :, :, :, :)

      call mp_release_bytes(c_loc(f_pointer))
      nullify (f_pointer)
   end subroutine mp_release_int_6d

end module fmempool
