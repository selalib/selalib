!> @ingroup splines
!> @author Laura Mendoza (IPP-Garching)
!> @brief Pre-filter for box-splines quasi interpolation
!> @details This module defines pre-filters for quasi-interpolation for
!> box splines on a hexagonal mesh subdivided in equilateral triangles
!> Reference : Condat2006 "Three-directional box splines"
module sll_m_hex_pre_filters
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_hexagonal_meshes, only: &
    sll_t_hex_mesh_2d

  implicit none

  public :: &
    sll_s_pre_filter_pfir, &
    sll_f_pre_filter_int, &
    sll_f_pre_filter_piir2, &
    sll_f_pre_filter_piir1


  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

  !---------------------------------------------------------------------------
  !> @brief Pre-filter PFIR to compute the box splines coefficients
  !> @details Pre-filter to compute the box splines coefficients.
  !> Reference : @Condat and Van De Ville (2007)
  !> "Quasi-interpolating spline models for hexagonally-sampled data."
  !> @param[IN] mesh sll_t_hex_mesh_2d: hex-mesh containing the mesh description
  !> @param[IN] deg integer: representing the degree of the spline
  !> @param[OUT] weight float: filter (aka weight) at the local index
  subroutine sll_s_pre_filter_pfir(mesh, deg, weight_tab)
    sll_real64, allocatable, intent(out) :: weight_tab(:)
    type(sll_t_hex_mesh_2d)     :: mesh
    sll_int32, intent(in)     :: deg
    sll_int32                 :: num_wei
    sll_int32                 :: ierr
    sll_int32                 :: index

    SLL_ASSERT(mesh%num_cells>0)
    num_wei = 3*deg*(deg+1) + 1
    SLL_ALLOCATE(weight_tab(num_wei), ierr)
    weight_tab(:) = 0.0_f64

    select case(deg)
    case(1)
       ! prefiltre PFIR for box-splines chi2
       ! with coefficients h0 = 5/4 and h1 = -1/24
       weight_tab(1) = 5._f64/4._f64
       weight_tab(2:7) = -1._f64/24._f64
    case(2)
       ! prefiltre PFIR for box-splines of deg =2 chi4
       ! with coefficients h0 = 37/20 h1 =-41/240 h2 = 7/240
       weight_tab(1)   = 37._f64/20._f64
       weight_tab(2:7) = -41._f64/240._f64
       do index=8,18,2
          weight_tab(index+1) = 7._f64/240._f64
       end do
    case(4)
       ! prefiltre PFIR for box-splines of deg =4 chi8
       ! with coefficients  h0 : -300538194444442.,
       ! h1 : 100179398148148.,   h2 : -100179398148149.0,
       ! h3 : 50089699074074.4, h4 : -0.0102843915343915
       weight_tab(1)   = -300538194444442._f64
       weight_tab(2:7) = 100179398148148.0_f64
       do index=8,18,2
          weight_tab(index+1) = -100179398148149.0_f64
          weight_tab(index)   =  50089699074074.4_f64
       end do
       do index=19,36
          if (modulo(index-1,3).ne.0) then
             print *, index
             weight_tab(index) = -0.0102843915343915_f64
          end if
       end do
       print *, ""
       print *, "sum =", SUM(weight_tab)
       print *, ""
    case(5)
       weight_tab(1)   = 631779569171640.0_f64
       weight_tab(2:7) = -192578390410100.0
       do index=8,18,2
          weight_tab(index+1) = 174563591096320.0
          weight_tab(index)   = -105296594861938.41
       end do
       do index=19,36
          if (modulo(index-1,3).ne.0) then
             print *, index
             weight_tab(index) = 18014799313778.578
          else
             weight_tab(index) = -18014799313778.602
          end if
       end do
       print *, ""
       print *, weight_tab
       print *, "sum =", SUM(weight_tab)
       print *, ""
    case default
       print *, 'ERROR: sll_s_pre_filter_pfir(...): ', &
            '     function not implemented for splines of degree > 4'
       print *, "Exiting..."
       STOP
    end select

  end subroutine sll_s_pre_filter_pfir


  !---------------------------------------------------------------------------
  !> @brief Pre-filter PIIR2 to compute the box splines coefficients
  !> @details Pre-filter PIIR2 to compute the box splines coefficients.
  !> Reference: Condat and Van De Ville (2007),
  !> "Quasi-interpolating spline models for hexagonally-sampled data."
  !> @param[IN] mesh sll_t_hex_mesh_2d: hex-mesh containing the mesh description
  !> @param[IN] local_index integer: local index of the point we want the filter
  !> @param[IN] deg integer: degree of the spline
  !> @param[OUT] weight float: filter (aka weight) at the local index
  function sll_f_pre_filter_piir2(mesh, local_index, deg) result(weight)
    type(sll_t_hex_mesh_2d)      :: mesh
    sll_int32, intent(in)      :: local_index
    sll_int32, intent(in)      :: deg
    sll_real64, allocatable    :: weights_tab(:)
    sll_real64                 :: weight
    sll_int32                  :: k1, k2
    sll_int32                  :: ierr
    sll_int32                  :: i

    if (deg .eq. 1) then
       ! prefiltre PIIR2 for box-splines chi2
       ! with coefficients h0 = 11/12 and h1 = 1/24
       !             |h1 0  0 |   |0  0  0 |   |h1 0  0 |
       ! prefilter = |0  h0 0 | * |h1 h0 h1| * |h0 0  0 |
       !             |0  0  h1|   |0  0  0 |   |h1 0  0 |
       ! where '*' symbolizes the 2d convolution operator
       SLL_ALLOCATE(weights_tab(19), ierr)
       weights_tab(1)   = 1775._f64/2304._f64
       weights_tab(2:7) = 253._f64/6912._f64
       do i=8,18,2
          weights_tab(i)   = 1.0_f64/13824._f64
          weights_tab(i+1) = 11._f64/6912._f64
       end do

       if (local_index .le. 19) then
          weight = weights_tab(local_index)
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
       if (local_index .eq. 1) then
          ! origin
          weight = 244301._f64/460800._f64
       else if (local_index .le. 7) then
          ! First hexagon
          weight = 42269._f64/576000._f64
       else if (local_index .le. 19) then
          ! Second hexagon
          if (modulo(local_index, 2) .eq. 0) then
             weight = -11809._f64/6912000._f64
          else
             weight = 1067._f64/144000._f64
          end if
       else if (local_index .le. 37) then
          !Third hexagon
          k1 = mesh%global_to_hex1(local_index)
          k2 = mesh%global_to_hex2(local_index)
          if ((k1.eq.0).or.(k2.eq.0).or.(k1.eq.k2)) then
             weight = -23._f64/576000._f64
          else
             weight = -109._f64/288000._f64
          end if
       else if (local_index .le. 61) then
          k1 = mesh%global_to_hex1(local_index)
          k2 = mesh%global_to_hex2(local_index)
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
    else if (deg .eq. 3) then 
       ! prefiltre PIIR2 for box-splines chi6
       ! with coefficients h0 = 173863/241920, h1 = 47309/322560,
       !                   h2 = -209/32256,    h3 = 457/967680
       !  piir2 = &
       !|h3 0  0  0  0  0  0 |   |0  0  0  0  0  0  0 |   |h3 0  0  0  0  0  0 |
       !|0  h2 0  0  0  0  0 |   |0  0  0  0  0  0  0 |   |h2 0  0  0  0  0  0 |
       !|0  0  h1 0  0  0  0 |   |0  0  0  0  0  0  0 |   |h1 0  0  0  0  0  0 |
       !|0  0  0  h0 0  0  0 | * |h3 h2 h1 h0 h1 h2 h3| * |h0 0  0  0  0  0  0 |
       !|0  0  0  0  h1 0  0 |   |0  0  0  0  0  0  0 |   |h1 0  0  0  0  0  0 |
       !|0  0  0  0  0  h2 0 |   |0  0  0  0  0  0  0 |   |h2 0  0  0  0  0  0 |
       !|0  0  0  0  0  0  h3|   |0  0  0  0  0  0  0 |   |h3 0  0  0  0  0  0 |
       ! where '*' symbolizes the 2d convolution operator
       if (local_index .eq. 1) then
          ! origin
          weight = 244301._f64/460800._f64
       else if (local_index .le. 7) then
          ! First hexagon
          weight = 42269._f64/576000._f64
       else if (local_index .le. 19) then
          ! Second hexagon
          if (modulo(local_index, 2) .eq. 0) then
             weight = -11809._f64/6912000._f64
          else
             weight = 1067._f64/144000._f64
          end if
       else if (local_index .le. 37) then
          k1 = mesh%global_to_hex1(local_index)
          k2 = mesh%global_to_hex2(local_index)
          !Third hexagon
          if ((k1.eq.0).or.(k2.eq.0).or.(k1.eq.k2)) then
             weight = -23._f64/576000._f64
          else
             weight = -109._f64/288000._f64
          end if
       else if (local_index .le. 61) then
          k1 = mesh%global_to_hex1(local_index)
          k2 = mesh%global_to_hex2(local_index)
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
  end function sll_f_pre_filter_piir2

  !---------------------------------------------------------------------------
  !> @brief Pre-filter PIIR1 to compute the box splines coefficients
  !> @details Pre-filter to compute the box splines coefficients.
  !> Reference : @Condat and Van De Ville (2007)
  !> "Quasi-interpolating spline models for hexagonally-sampled data."
  !> @param[IN] mesh sll_t_hex_mesh_2d: hex-mesh containing the mesh description
  !> @param[IN] local_index integer: local index of the point we want the filter
  !> @param[IN] deg integer: degree of the spline
  !> @param[OUT] weight float: filter (aka weight) at the local index
  function sll_f_pre_filter_piir1(mesh, local_index, deg) result(weight)
    type(sll_t_hex_mesh_2d)     :: mesh
    sll_int32, intent(in)     :: local_index
    sll_int32, intent(in)     :: deg
    sll_real64                :: weight
    sll_int32                 :: k1, k2
    k1 = mesh%global_to_hex1(local_index)
    k2 = mesh%global_to_hex2(local_index)

    if (deg .eq. 1) then
       ! prefiltre PIIR1 for box-splines chi2
       ! with coefficients h0 = 3/4 and h1 = 1/24
       if (local_index .eq. 1) then
          weight = 3._f64/4._f64
       else if (local_index .le. 7) then
          weight = 1._f64/24._f64
       else
          weight = 0._f64
       end if
    else if (deg .eq. 2) then
       ! prefiltre PIIR1 for box-splines of deg =2 chi4
       ! with coefficients h0 = 29/60 h1 = 7/80 h2 = -1/720
       if (local_index .eq. 1) then
          weight = 29._f64/60._f64
       else if (local_index .le. 7) then
          weight = 7._f64/80._f64
       else if ((local_index.le.19).and.(modulo(local_index, 2).eq.1)) then
          weight = -1._f64/720._f64
       else
          weight = 0._f64
       end if
       ! print *, 'ERROR: pre_filter_piir1(...): ', &
       !      '     function not implemented for degree 2 splines '
       ! print *, "Exiting..."
       ! STOP

    end if
  end function sll_f_pre_filter_piir1

  !---------------------------------------------------------------------------
  !> @brief Pre-filter PINT to compute the box splines coefficients
  !> @details Pre-filter to compute the box splines coefficients.
  !> Reference : @Condat and Van De Ville (2007)
  !> "Quasi-interpolating spline models for hexagonally-sampled data."
  !> @param[IN] mesh sll_t_hex_mesh_2d: hex-mesh containing the mesh description
  !> @param[IN] local_index integer: local index of the point we want the filter
  !> @param[IN] deg integer: degree of the spline
  !> @param[OUT] weight float: filter (aka weight) at the local index
  function sll_f_pre_filter_int(mesh, local_index, deg) result(weight)
    type(sll_t_hex_mesh_2d) :: mesh
    sll_int32, intent(in)     :: local_index
    sll_int32, intent(in)     :: deg
    sll_real64                :: weight

    SLL_ASSERT(mesh%num_cells>0)
    select case(deg)
    case(1)
       if (local_index.eq.1) then
          weight = 1._f64
       else
          weight = 0._f64
       end if
    case(2)
       ! prefiltre int for box-splines chi2
       if (local_index .eq. 1) then
          weight = 0.5_f64
       else if (local_index .le. 7) then
          weight = 1._f64/12._f64
       else
          weight = 0._f64
       end if
    case(3)
       ! prefiltre int for box-splines chi3
       ! with coefficients h0 = 3/4 and h1 = 1/24
       if (local_index .eq. 1) then
          weight = 5._f64/12._f64!1._f64/3._f64
       else if (local_index .le. 7) then
          weight = 1._f64/18._f64
       else if ((local_index.le.19).and.(modulo(local_index, 2).eq.1)) then
          weight = 1._f64/28._f64!1._f64/26._f64
       else
          weight = 0._f64
       end if
    case default
       print *, 'ERROR: sll_f_pre_filter_int(...): ', &
            '     function not implemented for degree > 2 splines '
       print *, "Exiting..."
       STOP
    end select
  end function sll_f_pre_filter_int


end module sll_m_hex_pre_filters
