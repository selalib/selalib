
#define MODE_X 2
#define MODE_Y 2
#define OMEGA sqrt((MODE_X*sll_pi)**2+(MODE_Y*sll_pi)**2)

real(8) function sol_bz( x1, x2, time)

   use sll_constants
   implicit none

   real(8), intent(in) :: x1
   real(8), intent(in) :: x2
   real(8), intent(in) :: time



   sol_bz =   - cos(MODE_X*sll_pi*x1)    &
              * cos(MODE_Y*sll_pi*x2)    &
              * cos(OMEGA*time)

   return

end function sol_bz

real(8) function sol_ex( x1, x2, time)

   use sll_constants
   implicit none

   real(8), intent(in) :: x1
   real(8), intent(in) :: x2
   real(8), intent(in) :: time

   sol_ex =   + cos(MODE_X*sll_pi*x1)    &
              * sin(MODE_Y*sll_pi*x2)    &
              * sin(OMEGA*time) * MODE_Y*sll_pi/OMEGA
   return

end function sol_ex

real(8) function sol_ey( x1, x2, time)

   use sll_constants
   implicit none

   real(8), intent(in) :: x1
   real(8), intent(in) :: x2
   real(8), intent(in) :: time

   sol_ey =   - sin(MODE_X*sll_pi*x1)    &
              * cos(MODE_Y*sll_pi*x2)    &
              * sin(OMEGA*time) * MODE_X*sll_pi/OMEGA
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


real(8) function linear_x( x1, x2, time)

   use sll_constants
   implicit none

   real(8), intent(in) :: x1
   real(8), intent(in) :: x2
   real(8), intent(in) :: time

   linear_x =   x1
   return

end function linear_x


real(8) function linear_y( x1, x2, time)

   use sll_constants
   implicit none

   real(8), intent(in) :: x1
   real(8), intent(in) :: x2
   real(8), intent(in) :: time

   linear_y =   x2
   return

end function linear_y


