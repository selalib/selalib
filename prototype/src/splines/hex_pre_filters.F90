!**************************************************************
!  This module defines pre-filters for quasi-interpolation for 
!  box splines on a hexagonal mesh subdivided in equilateral triangles
!  Reference :
!     @Condat2006 "Three-directional box splines"
!  Author : 
!     Laura Mendoza (mela@ipp.mpg.de)
!************************************************************** 


module hex_pre_filters
#include "sll_working_precision.h"
use hex_mesh

implicit none

contains 
  ! Pre-filter to compute the box splines coefficients
  ! Reference : @Condat and Van De Ville (2007)
  !             "Quasi-interpolating spline models 
  !             for hexagonally-sampled data."
  function pre_filter_piir2(mesh, local_index, deg) result(weight)
      type(hex_mesh_2d) :: mesh 
      sll_int32, intent(in)      :: local_index
      sll_int32, intent(in)      :: deg
      sll_real64                 :: weight
      sll_int32                  :: k1, k2
      k1 = mesh%global_to_hex1(local_index)
      k2 = mesh%global_to_hex2(local_index)


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
              weight = 0._f64
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
              weight = 0._f64
          end if
      end if
   end function pre_filter_piir2

  ! Pre-filter to compute the box splines coefficients
  ! Reference : @Condat and Van De Ville (2007)
  !             "Quasi-interpolating spline models 
  !             for hexagonally-sampled data."
  function pre_filter_piir1(mesh, local_index, deg) result(weight)
    type(hex_mesh_2d) :: mesh 
      sll_int32, intent(in)     :: local_index
      sll_int32, intent(in)     :: deg
      sll_real64                :: weight
      sll_int32                 :: k1, k2
      k1 = mesh%global_to_hex1(local_index)
      k2 = mesh%global_to_hex2(local_index)


      if (deg .eq. 1) then 
          ! prefiltre PIIR1 for box-splines chi2
          ! with coefficients h0 = 3/4 and h1 = 1/24
          if (local_index .eq. 0) then
              weight = 3._f64/4._f64
          else if (local_index .lt. 7) then
              weight = 1._f64/24._f64
          else
              weight = 0._f64
          end if
       else if (deg .eq. 2) then 
          ! prefiltre PIIR1 for box-splines of deg =2 chi4
          ! with coefficients h0 = 29/60 h1 = 7/80 h2 = -1/720
          if (local_index .eq. 0) then
             weight = 29._f64/60._f64
          else if (local_index .lt. 7) then
             weight = 7._f64/80._f64
          else if ((local_index.lt.19).and.(modulo(local_index, 2).eq.0)) then
             weight = -1._f64/720._f64
          else
             weight = 0._f64
          end if
       ! print *, 'ERROR: pre_filter_piir1(...): ', &
       !      '     function not implemented for degree 2 splines '
       ! print *, "Exiting..."
       ! STOP

      end if
   end function pre_filter_piir1


  ! Pre-filter to compute the box splines coefficients
  ! Reference : @Condat and Van De Ville (2007)
  !             "Quasi-interpolating spline models 
  !             for hexagonally-sampled data."
   function pre_filter_pfir(mesh, local_index, deg) result(weight)
     type(hex_mesh_2d) :: mesh 
     sll_int32, intent(in)     :: local_index
     sll_int32, intent(in)     :: deg
     sll_real64                :: weight
     sll_int32                 :: k1, k2
     k1 = mesh%global_to_hex1(local_index)
     k2 = mesh%global_to_hex2(local_index)
     

      if (deg .eq. 1) then 
          ! prefiltre PFIR for box-splines chi2
          ! with coefficients h0 = 5/4 and h1 = -1/24

          if (local_index .eq. 0) then
              weight = 5._f64/4._f64
          else if (local_index .lt. 7) then
              weight = -1._f64/24._f64
          else
              weight = 0._f64
          end if
       else if (deg .eq. 2) then 
          ! prefiltre PFIR for box-splines of deg =2 chi4
          ! with coefficients h0 = 37/20 h1 =-41/240 h2 = 7/240
          if (local_index .eq. 0) then
             weight = 37._f64/20._f64
          else if (local_index .lt. 7) then
             weight = -41._f64/240._f64
          else if ((local_index.lt.19).and.(modulo(local_index, 2).eq.0)) then
             weight = 7._f64/240._f64
          else
             weight = 0._f64
          end if
      ! else if (deg .eq. 2) then 
      !  print *, 'ERROR: pre_filter_pfir(...): ', &
      !       '     function not implemented for order 2 splines '
      !  print *, "Exiting..."
      !  STOP
      end if
   end function pre_filter_pfir



  ! Pre-filter to compute the box splines coefficients
  ! Reference : @Condat and Van De Ville (2007)
  !             "Quasi-interpolating spline models 
  !             for hexagonally-sampled data."
  function pre_filter_int(mesh, local_index, deg) result(weight)
      type(hex_mesh_2d) :: mesh
      sll_int32, intent(in)     :: local_index
      sll_int32, intent(in)     :: deg
      sll_real64                :: weight
      sll_int32                 :: k1, k2
      k1 = mesh%global_to_hex1(local_index)
      k2 = mesh%global_to_hex2(local_index)


      if (deg .eq. 2) then 
          ! prefiltre int for box-splines chi2
          ! with coefficients h0 = 3/4 and h1 = 1/24
          if (local_index .eq. 0) then
              weight = 0.5_f64
          else if (local_index .lt. 7) then
              weight = 1._f64/12._f64
          else
              weight = 0._f64
          end if
      else if (deg .gt. 2) then 
       print *, 'ERROR: pre_filter_int(...): ', &
            '     function not implemented for degree > 2 splines '
       print *, "Exiting..."
       STOP

      end if
   end function pre_filter_int


end module hex_pre_filters
