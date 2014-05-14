
real(8) function sol_bz( x1, x2, time)

   use sll_constants
   implicit none

   real(8), intent(in) :: x1
   real(8), intent(in) :: x2
   real(8), intent(in) :: time

   sol_bz =     cos(x1)    &
              * cos(x2)    &
              * cos(time)

   return

end function sol_bz

real(8) function sol_ex( x1, x2, time)

   use sll_constants
   implicit none

   real(8), intent(in) :: x1
   real(8), intent(in) :: x2
   real(8), intent(in) :: time

   sol_ex =   - cos(x1)    &
              * sin(x2)    &
              * sin(time)
   return

end function sol_ex

real(8) function sol_ey( x1, x2, time)

   use sll_constants
   implicit none

   real(8), intent(in) :: x1
   real(8), intent(in) :: x2
   real(8), intent(in) :: time

   sol_ey =   + sin(x1)    &
              * cos(x2)    &
              * sin(time)
   return

end function sol_ey

real(8) function gaussian( x1, x2, time)

   use sll_constants
   implicit none

   real(8), intent(in) :: x1
   real(8), intent(in) :: x2
   real(8), intent(in) :: time

   gaussian =   exp(-(x1*x1+x2*x2)) * cos(time)
   return

end function gaussian


real(8) function add( x1, x2, time)

   use sll_constants
   implicit none

   real(8), intent(in) :: x1
   real(8), intent(in) :: x2
   real(8), intent(in) :: time

   add =   x1+x2
   return

end function add


