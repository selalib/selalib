!----------------------------------------------------------------------------
! C-like assert implementation for Fortran. To be used together with
! 'fassert.inc'.
!
! See README for details.
!----------------------------------------------------------------------------

module assert

   implicit none

contains

   ! Assertion routine designed to mimick the behavior of C's 'assert' from 'assert.h'.
   ! To be called from the ASSERT() macro from 'fassert.h', not directly.
   subroutine assertion(cond, file, line)
      character(len=*), intent(in) :: cond
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      character(len=6) :: line_str
      write (line_str, '(i6)') line
      write (*, '(a, a, a, a, a, a)') file, ":", trim(adjustl(line_str)), ": Assertion `", cond, "' failed."
      write (*, '(a)') "Aborted"
      stop 1
   end subroutine assertion

end module assert
