!> @brief  
!> Pre filters for quasi-interpolqtion for box splines
!> @Condat2006 :Three-directional box splines
!> 

module hex_pre_filters
#include "sll_working_precision.h"
#include "sll_splines.h"
use hex_mesh!, only:find_neighbour

implicit none

contains 
  ! Pre-filter to compute the box splines coefficients
  ! Reference : @Condat and Van De Ville (2007)
  !             "Quasi-interpolating spline models 
  !             for hexagonally-sampled data."
  function pre_filter_piir2(local_index, deg) result(weight)
      sll_int32, intent(in)     :: local_index
      sll_int32, intent(in)     :: deg
      sll_real64                :: weight
      sll_int32                 :: k1, k2
      k1 = from_global_index_k1(local_index)
      k2 = from_global_index_k2(local_index)


      if (deg .eq. 1) then 
          ! prefiltre PIIR2 for box-splines chi2
          ! with coefficients h0 = 11/12 and h1 = 1/24
          !             |h1 0  0 |   |0  0  0 |   |h1 0  0 |
          ! prefilter = |0  h0 0 | * |h1 h0 h1| * |h0 0  0 |
          !             |0  0  h1|   |0  0  0 |   |h1 0  0 |
          ! where '*' symbolizes the 2d convolution operator

          if (local_index .eq. 0) then
              weight = 1775._f64/2304._f64
          else if (local_index .lt. 7) then
              weight = 253._f64/6912._f64
          else if (local_index .lt. 19) then
              if (modulo(local_index, 2) .eq. 1) then
                  weight = 1._f64/13824._f64
              else
                  weight = 11._f64/6912._f64
              end if
          else
              weight = 0.
          end if

      else if (deg .eq. 2) then 
          ! prefiltre PIIR2 for box-splines chi4
          ! with coefficients h0 = 97/120, h1 = 1/10 and h2 = -1/240
          !             |h2 0  0  0  0 |   |0  0  0  0  0 |   |h2 0  0  0  0 |
          !             |0  h1 0  0  0 |   |0  0  0  0  0 |   |h1 0  0  0  0 |
          ! prefilter = |0  0  h0 0  0 | * |h2 h1 h0 h1 h2| * |h0 0  0  0  0 |
          !             |0  0  0  h1 0 |   |0  0  0  0  0 |   |h1 0  0  0  0 |
          !             |0  0  0  0  h2|   |0  0  0  0  0 |   |h2 0  0  0  0 |
          ! where '*' symbolizes the 2d convolution operator

          if (local_index .eq. 0) then
              ! origin
              weight = 244301._f64/460800._f64
          else if (local_index .lt. 7) then
              ! First hexagon
              weight = 42269._f64/576000._f64
          else if (local_index .lt. 19) then
              ! Second hexagon
              if (modulo(local_index, 2) .eq. 1) then
                  weight = -11809._f64/6912000._f64
              else
                  weight = 1067._f64/144000._f64
              end if
          else if (local_index .lt. 37) then
              !Third hexagon
              if ((k1.eq.0).or.(k2.eq.0).or.(k1.eq.k2)) then
                  weight = -23._f64/576000._f64
              else
                  weight = -109._f64/288000._f64
              end if
          else if (local_index .lt. 61) then
              ! Forth hexagon
              if ((k1.eq.0).or.(k2.eq.0).or.(k1.eq.k2)) then
                  weight = -1._f64/13824000._f64
              else if ((abs(k1).eq.2).or.(abs(k2).eq.2))then
                  weight = 97._f64/6912000._f64
              else 
                  weight = 1._f64/576000._f64
              end if
          else
              weight = 0.
          end if
      end if
   end function pre_filter_piir2

  ! Pre-filter to compute the box splines coefficients
  ! Reference : @Condat and Van De Ville (2007)
  !             "Quasi-interpolating spline models 
  !             for hexagonally-sampled data."
  function pre_filter_piir1(local_index, deg) result(weight)
      sll_int32, intent(in)     :: local_index
      sll_int32, intent(in)     :: deg
      sll_real64                :: weight
      sll_int32                 :: k1, k2
      k1 = from_global_index_k1(local_index)
      k2 = from_global_index_k2(local_index)


      if (deg .eq. 1) then 
          ! prefiltre PIIR2 for box-splines chi2
          ! with coefficients h0 = 3/4 and h1 = 1/24
          !             |h1 0  0 |   |0  0  0 |   |h1 0  0 |
          ! prefilter = |0  h0 0 | * |h1 h0 h1| * |h0 0  0 |
          !             |0  0  h1|   |0  0  0 |   |h1 0  0 |
          ! where '*' symbolizes the 2d convolution operator

          if (local_index .eq. 0) then
              weight = 2917._f64/6912._f64
          else if (local_index .lt. 7) then
              weight = 19._f64/768._f64
          else if (local_index .lt. 19) then
              if (modulo(local_index, 2) .eq. 1) then
                  weight = 1._f64/13824._f64
              else
                  weight = 1._f64/768._f64
              end if
          else
              weight = 0.
          end if

      else if (deg .eq. 2) then 
       print *, 'ERROR: pre_filter_piir1(...): ', &
            '     function not implemented for order 2 splines '
       print *, "Exiting..."
       STOP

      end if
   end function pre_filter_piir1


  ! Pre-filter to compute the box splines coefficients
  ! Reference : @Condat and Van De Ville (2007)
  !             "Quasi-interpolating spline models 
  !             for hexagonally-sampled data."
  function pre_filter_pfir(local_index, deg) result(weight)
      sll_int32, intent(in)     :: local_index
      sll_int32, intent(in)     :: deg
      sll_real64                :: weight
      sll_int32                 :: k1, k2
      k1 = from_global_index_k1(local_index)
      k2 = from_global_index_k2(local_index)


      if (deg .eq. 1) then 
          ! prefiltre PIIR2 for box-splines chi2
          ! with coefficients h0 = 5/4 and h1 = -1/24
          !             |h1 0  0 |   |0  0  0 |   |h1 0  0 |
          ! prefilter = |0  h0 0 | * |h1 h0 h1| * |h0 0  0 |
          !             |0  0  h1|   |0  0  0 |   |h1 0  0 |
          ! where '*' symbolizes the 2d convolution operator

          if (local_index .eq. 0) then
              weight = 2949._f64/1510._f64
          else if (local_index .lt. 7) then
              weight = -145._f64/2304._f64
          else if (local_index .lt. 19) then
              if (modulo(local_index, 2) .eq. 1) then
                  weight = -1._f64/13824._f64
              else
                  weight = 5._f64/2304._f64
              end if
          else
              weight = 0.
          end if

      else if (deg .eq. 2) then 
       print *, 'ERROR: pre_filter_pfir(...): ', &
            '     function not implemented for order 2 splines '
       print *, "Exiting..."
       STOP

      end if
   end function pre_filter_pfir

end module hex_pre_filters
