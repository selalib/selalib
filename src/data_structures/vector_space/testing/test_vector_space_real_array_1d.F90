program test_vector_space_real_array_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use sll_m_working_precision, only: f64

   use sll_m_vector_space_real_array_1d, only: sll_t_vector_space_real_array_1d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Working precision
   integer, parameter :: wp = f64

   integer :: i, n

   real(wp), parameter :: tol = 1.0e-14_wp

   ! Reference array to check results
   real(wp), allocatable :: r(:)

   type(sll_t_vector_space_real_array_1d) :: v, w, z, a(2)

   ! For CTest
   logical :: passed, success

   passed = .true.

   n = 3

   ! Allocate reference array
   allocate (r(n))

   ! Allocate and initialize vector space
   allocate (v%array(n), source=(/(real(i, wp), i=1, n)/))

   ! Allocate and initialize auxiliary vector space for composite operations
   allocate (z%array(n), source=(/(0.0_wp, i=1, n)/))

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_1d % copy
   !-----------------------------------------------------------------------------

   allocate (w%array(size(v%array)))
   call w%copy(v)

   ! Check test
   if (all((v%array - w%array) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_1d % copy FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_1d % incr
   !-----------------------------------------------------------------------------

   call v%incr(w)

   ! Set reference array
   r(:) = (/2.0_wp, 4.0_wp, 6.0_wp/)

   ! Check test
   if (all((v%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_1d % incr FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_1d % scal
   !-----------------------------------------------------------------------------

   call v%scal(0.5_wp)

   ! Set reference array
   r(:) = (/1.0_wp, 2.0_wp, 3.0_wp/)

   ! Check test
   if (all((v%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_1d % scal FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_1d % add
   !-----------------------------------------------------------------------------

   call z%add(v, w)

   ! Set reference array
   r(:) = (/2.0_wp, 4.0_wp, 6.0_wp/)

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_1d % add FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_1d % mult
   !-----------------------------------------------------------------------------

   call z%mult(0.5_wp, v)

   ! Set reference array
   r(:) = (/0.5_wp, 1.0_wp, 1.5_wp/)

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_1d % mult FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_1d % mult_add
   !-----------------------------------------------------------------------------

   call z%mult_add(0.5_wp, v, w)

   ! Set reference array
   r(:) = (/1.5_wp, 3.0_wp, 4.5_wp/)

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_1d % mult_add FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_1d % incr_mult
   !-----------------------------------------------------------------------------

   call z%incr_mult(-1.0_wp, v)

   ! Set reference array
   r(:) = (/0.5_wp, 1.0_wp, 1.5_wp/)

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_1d % incr_mult FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_1d % lcmb
   !-----------------------------------------------------------------------------

   a(1) = v
   a(2) = w
   call z%lcmb((/0.5_wp, -1.0_wp/), a)

   ! Set reference array
   r(:) = (/-0.5_wp, -1.0_wp, -1.5_wp/)

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_1d % lcmb FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_1d % incr_lcmb
   !-----------------------------------------------------------------------------

   a(1) = v
   a(2) = w
   call z%incr_lcmb((/0.5_wp, -1.0_wp/), a)

   ! Set reference array
   r(:) = (/-1.0_wp, -2.0_wp, -3.0_wp/)

   ! Check test
   if (all((z%array - r) == 0.0_wp)) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_1d % incr_lcmb FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_1d % norm
   !-----------------------------------------------------------------------------

   ! Check test
   if (abs(v%norm() - sqrt(14.0_wp)) < tol) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_1d % norm FAILED"
   end if

   passed = passed .and. success

   !-----------------------------------------------------------------------------
   ! Test sll_t_vector_space_real_array_1d % inner
   !-----------------------------------------------------------------------------

   ! Check test
   if (abs(v%inner(w) - 14.0_wp) < tol) then
      success = .true.
   else
      success = .false.
      write (*, '(/a)') "Test of sll_t_vector_space_real_array_1d % inner FAILED"
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

end program test_vector_space_real_array_1d
