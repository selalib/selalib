!----------------------------------------------------------------------------
! Simple tests and usage examples of the fmempool library.
!----------------------------------------------------------------------------

program main
  use fmempool
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none

  integer, pointer :: int_1d(:), int_2d(:,:)
  double precision, pointer :: dbl_1d(:), dbl_3d(:,:,:), dbl_private(:)
  integer :: min_1d(1), max_1d(1), min_2d(2), max_2d(2), min_3d(3), max_3d(3)
  integer :: i

  call mp_init(verbosity=.true.)

  min_1d(1) = 1
  max_1d(1) = 2*1024
  call mp_acquire(dbl_1d, min_1d, max_1d)

  call mp_statistics()

  call mp_release(dbl_1d)

  max_1d(1) = 1024

!$omp parallel default(shared) private(dbl_private)
  call mp_acquire(dbl_private, min_1d, max_1d)
!$omp barrier
  call mp_statistics()
!$omp barrier
  call mp_release(dbl_private)
!$omp barrier
  call mp_statistics()
!$omp barrier
!$omp end parallel

  call mp_statistics()
  call mp_compactify()
  call mp_statistics()

  min_1d(1) = 1
  max_1d(1) = 1024
  call mp_acquire(int_1d, min_1d, max_1d)
  int_1d = 0
  call mp_release(int_1d)

  min_2d = [1, 1]
  max_2d = [16, 16]
  call mp_acquire(int_2d, min_2d, max_2d)
  int_2d = 1
  call mp_release(int_2d)

  ! reuse the same memory, verify by summing the previously set elements
  min_1d(1) = 1
  max_1d(1) = 1024
  call mp_acquire(int_1d, min_1d, max_1d)
  if (.not. mp_disabled()) then
    if ( sum(int_1d) /= product(max_2d) ) then
      stop 1
    endif
  endif

  call mp_finalize()
end program main
