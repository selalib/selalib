program test_compression

   use sll_m_compression

   implicit none

   type(sll_t_compressed_buffer) :: comp  ! data structure containing compressed data and offsets
   double precision, pointer :: arr(:), brr(:), delta(:)
   double precision, parameter :: pi = 3.14159265359
   double precision :: t0, t1
   integer :: prec
   integer :: n_elem, i
   logical, parameter :: print_verbose = .true.

   prec = 32
   call set_compression_precision(prec)

   n_elem = 64*2**18
   allocate (arr(n_elem))

   do i = 1, n_elem
      arr(i) = 10.*sin(real(i)/real(n_elem)*2.*pi)
   end do

   t0 = get_time()
   call deflate_buffer_real64(arr, comp)
   t1 = get_time()
   write (*, '(A,F8.2)') " deflation rate [MB/s]   = ", comp%n_bytes_deflated_total/(t1 - t0)/real(2**20)
   call print_compression_information(comp, print_verbose)

   allocate (brr(n_elem))
   brr(:) = 0.0

   t0 = get_time()
   call inflate_buffer_real64(brr, comp)
   t1 = get_time()
   write (*, '(A,F8.2)') " inflation rate [MB/s]   = ", comp%n_bytes_deflated_total/(t1 - t0)/real(2**20)

   call deallocate_compressed_buffer_obj(comp)

   allocate (delta(n_elem))
   delta(:) = 0.0

   do i = 1, n_elem
      delta(i) = abs(arr(i) - brr(i))
   end do

!  if (maxval(delta) <= accuracy) then
!    write(*,*) "PASSED"
!  else
!    write(*,*) "FAILED"
!  endif

   deallocate (arr)
   deallocate (brr)
   deallocate (delta)
end program test_compression
