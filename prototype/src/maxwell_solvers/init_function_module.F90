
real(8) function fcos( x1, x2, time)

   use sll_constants
   implicit none

   real(8), intent(in) :: x1
   real(8), intent(in) :: x2
   real(8), intent(in) :: time

   fcos =   - cos(sll_pi*x1)    &
            * cos(sll_pi*x2)    &
            * cos(time)
   return

end function fcos


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

   add =   x1!+x2
   return

end function add


