module m_maxwell_helper_functions
#include "sll_working_precision.h"

  use sll_m_constants, only : &
       sll_pi
     

#define MODE_X 2
#define MODE_Y 2
#define OMEGA sqrt((MODE_X*sll_pi)**2+(MODE_Y*sll_pi)**2)

implicit none

 contains

   sll_real64 function sol_bz( x1, x2, time)

     sll_real64, intent(in) :: x1
     sll_real64, intent(in) :: x2
     sll_real64, intent(in) :: time
     

     sol_bz =   - cos(MODE_X*sll_pi*x1)    &
          * cos(MODE_Y*sll_pi*x2)    &
          * cos(OMEGA*time)
     
     return

   end function sol_bz

   sll_real64 function sol_ex( x1, x2, time)

     sll_real64, intent(in) :: x1
     sll_real64, intent(in) :: x2
     sll_real64, intent(in) :: time

     sol_ex =   + cos(MODE_X*sll_pi*x1)    &
          * sin(MODE_Y*sll_pi*x2)    &
          * sin(OMEGA*time) * real(MODE_Y,f64)*sll_pi/OMEGA
     return
     
   end function sol_ex

   sll_real64 function sol_ey( x1, x2, time)

     sll_real64, intent(in) :: x1
     sll_real64, intent(in) :: x2
     sll_real64, intent(in) :: time
   
     sol_ey =   - sin(MODE_X*sll_pi*x1)    &
          * cos(MODE_Y*sll_pi*x2)    &
          * sin(OMEGA*time) * real(MODE_X,f64)*sll_pi/OMEGA
     return

   end function sol_ey

   sll_real64 function gaussian( x1, x2, time)
    
     sll_real64, intent(in) :: x1
     sll_real64, intent(in) :: x2
     sll_real64, intent(in) :: time
    
     gaussian =   exp(-(x1*x1+x2*x2)) * cos(time)
     return
    
   end function gaussian

end module m_maxwell_helper_functions
