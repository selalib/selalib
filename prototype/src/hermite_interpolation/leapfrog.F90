module compute_characteristic

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_constants

  implicit none

contains


  subroutine slb_compute_characteristic_leapfrog( &
       x1,x2,Uxn,Uyn,i,j,y1,y2,dt )

    sll_real64,dimension(:),intent(in):: Uxn, Uyn, Uxn_1, Uyn_1
    sll_real64, intent(in)  :: x1, x2, dt
    sll_real64, intent(out) :: y1, y2
    sll_int32, intent(in)   :: i, j
    sll_real64              :: Uxn1, Uyn1, d1x, d1y,d1j0, d1j1
    
    ! needs to be optimized

    d1x = dt * Uxn(i,j)
    d1y = dt * Uyn(i,j)

    dij0 = d1x - dt * ( d1x*dxUxn + d1y*dyUxn ) 
    dij1 = d1y - dt * ( d1y*dyUyn + d1x*dxUyn )

    y1 = x1 - 2._f64 * dij0
    y2 = x2 - 2._f64 * dij1
	  

  end subroutine slb_compute_characteristic_leapfrog

  subroutine slb_compute_characteristic_adams_moulton2( &
       x1,x2,Uxn,Uyn,Uxn_1,Uyn_1,i,j,y1,y2,dt )

    sll_real64,dimension(:),intent(in):: Uxn, Uyn, Uxn_1, Uyn_1
    sll_real64, intent(in)  :: x1, x2, dt
    sll_real64, intent(out) :: y1, y2
    sll_int32, intent(in)   :: i, j
    sll_real64              :: Uxn1, Uyn1, d1x, d1y,d1j0, d1j1
    
    ! needs to be optimized

    Uxn1 = 2._f64*Uxn(i,j) - Uxn_1(i,j)
    Uyn1 = 2._f64*Uyn(i,j) - Uyn_1(i,j)

    d1x = 0.5_f64*dt * (Uxn1 + Uxn_1(i,j))
    d1y = 0.5_f64*dt * (Uyn1 + Uyn_1(i,j))

    dij0 = d1x - 0.5_f64*dt *( d1x*dxUxn + d1y*dyUxn )
    dij1 = d1y - 0.5_f64*dt *( d1x*dxUyn + d1y*dyUyn )

    y1 = x1 - dij0
    y2 = x2 - dij1

  end subroutine slb_compute_characteristic_adams_moulton2

end module compute_characteristic
