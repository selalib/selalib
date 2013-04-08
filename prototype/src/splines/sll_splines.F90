!> \file sll_splines.F90
!> \namespace sll_splines
!> \brief  
!> The splines module provides capabilities for 1D data interpolation with 
!> cubic B-splines and different boundary conditions
!>
!> (at the time of this writing: periodic, hermite). The data to be 
!> interpolated is represented by a simple array.  The spline coefficients 
!> and other information are stored in a spline object, which is also used 
!> to interpolate the fitted data.
!> 
module sll_splines

use sll_cubic_splines
use sll_quintic_splines
use sll_odd_degree_splines
use sll_arbitrary_degree_splines

end module sll_splines
