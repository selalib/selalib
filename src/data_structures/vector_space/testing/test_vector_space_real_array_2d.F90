program test_vector_space_real_array_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use sll_m_working_precision, only: f64

   use sll_m_vector_space_real_array_2d, only: sll_t_vector_space_real_array_2d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Working precision
   integer, parameter :: wp = f64

   integer :: i1, i2, n1, n2

   real(wp), parameter :: tol = 1.0e-14_wp

   ! Reference array to check results
   real(wp), allocatable :: r(:, :)

   type(sll_t_vector_space_real_array_2d) :: v, w, z, a(2)

   ! For CTest
   logical :: passed, success

   passed = .true.

   n1 = 3
   n2 = 4

   ! Allocate reference array
   allocate (r(n1, n2))

   ! Allocate and initialize vector space
   do i2 = 1, n2
      r(:, i2) = (/(real(i1, wp), i1=1, n1)/)
   end do
   allocate (v%array(n1, n2), source=r)

   ! Allocate and initialize auxiliary vector space for composite operations
   do i2 = 1, n2
      r(:, i2) = (/(0.0_wp, i1=1, n1)/)
   end do
   allocate (z%array(n1, n2), source=r)

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_2d % copy
   !-----------------------------------------------------------------------------

   allocate (w%array(size(v%array, 1), size(v%array, 2)))
   call w%copy(v)

   ! Check test
   if (all((v%array - w%array) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_2d % copy FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_2d % incr
   !-----------------------------------------------------------------------------

   call v%incr(w)

   ! Set reference array
   do i2 = 1, n2
      r(:, i2) = (/2.0_wp, 4.0_wp, 6.0_wp/)
   end do

   ! Check test
   if (all((v%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_2d % incr FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_2d % scal
   !-----------------------------------------------------------------------------

   call v%scal(0.5_wp)

   ! Set reference array
   do i2 = 1, n2
      r(:, i2) = (/1.0_wp, 2.0_wp, 3.0_wp/)
   end do

   ! Check test
   if (all((v%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_2d % scal FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_2d % add
   !-----------------------------------------------------------------------------

   call z%add(v, w)

   ! Set reference array
   do i2 = 1, n2
      r(:, i2) = (/2.0_wp, 4.0_wp, 6.0_wp/)
   end do

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_2d % add FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_2d % mult
   !-----------------------------------------------------------------------------

   call z%mult(0.5_wp, v)

   ! Set reference array
   do i2 = 1, n2
      r(:, i2) = (/0.5_wp, 1.0_wp, 1.5_wp/)
   end do

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_2d % mult FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_2d % mult_add
   !-----------------------------------------------------------------------------

   call z%mult_add(0.5_wp, v, w)

   ! Set reference array
   do i2 = 1, n2
      r(:, i2) = (/1.5_wp, 3.0_wp, 4.5_wp/)
   end do

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_2d % mult_add FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_2d % incr_mult
   !-----------------------------------------------------------------------------

   call z%incr_mult(-1.0_wp, v)

   ! Set reference array
   do i2 = 1, n2
      r(:, i2) = (/0.5_wp, 1.0_wp, 1.5_wp/)
   end do

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_2d % incr_mult FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_2d % lcmb
   !-----------------------------------------------------------------------------

   a(1) = v
   a(2) = w
   call z%lcmb((/0.5_wp, -1.0_wp/), a)

   ! Set reference array
   do i2 = 1, n2
      r(:, i2) = (/-0.5_wp, -1.0_wp, -1.5_wp/)
   end do

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_2d % lcmb FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_2d % incr_lcmb
   !-----------------------------------------------------------------------------

   a(1) = v
   a(2) = w
   call z%incr_lcmb((/0.5_wp, -1.0_wp/), a)

   ! Set reference array
   do i2 = 1, n2
      r(:, i2) = (/-1.0_wp, -2.0_wp, -3.0_wp/)
   end do

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_2d % incr_lcmb FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_2d % norm
   !-----------------------------------------------------------------------------

   ! Check test
   if (abs(v%norm() - sqrt(14.0_wp*n2)) < tol) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_2d % norm FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_2d % inner
   !-----------------------------------------------------------------------------

   ! Check test
   if (abs(v%inner(w) - 14.0_wp*n2) < tol) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_2d % inner FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Deallocate arrays
   !-----------------------------------------------------------------------------

   deallocate (v%array, w%array, z%array)

   !-----------------------------------------------------------------------------
   ! Check if all tests passed
   !-----------------------------------------------------------------------------
   if (passed) then
      write (*, '(/a/)') "PASSED"
   else
      write (*, *)
   end if

end program test_vector_space_real_array_2d
