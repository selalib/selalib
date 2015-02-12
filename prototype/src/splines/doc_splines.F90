!> @namespace sll_splines
!> @brief 
!> Library to use splines
!> @details
!>
!> <b> Modules available </b>
!>
!> - sll_cubic_splines
!> - sll_quintic_splines
!> - cubic_non_uniform_splines
!> - deboor_arbitrary_degree_splines
!> - sll_arbitrary_degree_splines
!> - sll_odd_degree_splines
!>
!> <b> How to use it </b>
!> - Header file : \code #include "sll_file_io.h" \endcode
!> - Link with   <code>-lsll_file_io</code> library
!>
!> <b> Examples </b>
!> \code
!>
!>  type(sll_cubic_spline_2d), pointer      :: spline
!>  sll_real64                              :: x1
!>  sll_real64                              :: x2
!>  sll_real64, dimension(npts1, npst2)     :: data_2d
!>
!>  spline => new_cubic_spline_2d( npts1, npts2,  &
!>                                 x1min, x1max, x2min, x2max, &
!>                                 SLL_PERIODIC, SLL_PERIODIC )
!>
!>  call compute_cubic_spline_2d( data_2d, spline )
!>
!>  x1 = call random_number()
!>  x1 = x1min + x1 * (x1max-x1min)
!>  x2 = call_random_number()
!>  x2 = x2min + x2 * (x2max-x2min)
!>
!>  print*, 'value at ', x1, x2, '=', interpolate_value_2D(x1,x2,spline)
!>
!> \endcode

