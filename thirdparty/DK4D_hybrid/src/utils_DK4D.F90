!===========================================================================
!> Different useful types for 4D drift-kinetic hybrid simulation
!>
!> \date 2014-08-19
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module utils_DK4D_module
#include "sll_working_precision.h"

  !===========================================================================
  !> Type used to define boundary conditions in each of 
  !>  the fourth directions
  !---------------------------------------------------------------------------
  type, public :: boundary_conditions_4d_t

    sll_int32 :: left_eta1
    sll_int32 :: right_eta1
    sll_int32 :: left_eta2
    sll_int32 :: right_eta2
    sll_int32 :: left_eta3
    sll_int32 :: right_eta3
    sll_int32 :: left_vpar
    sll_int32 :: right_vpar
    
  end type boundary_conditions_4d_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Type used to define boundary conditions in the three
  !>  directions (eta1,eta2,eta3)
  !---------------------------------------------------------------------------
  type, public :: boundary_conditions_3d_t

    sll_int32 :: left_eta1
    sll_int32 :: right_eta1
    sll_int32 :: left_eta2
    sll_int32 :: right_eta2
    sll_int32 :: left_eta3
    sll_int32 :: right_eta3
    
  end type boundary_conditions_3d_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Type used to define boundary conditions in the two
  !>  directions (eta1,eta2)
  !---------------------------------------------------------------------------
  type, public :: boundary_conditions_2d_t

    sll_int32 :: left_eta1
    sll_int32 :: right_eta1
    sll_int32 :: left_eta2
    sll_int32 :: right_eta2
    
  end type boundary_conditions_2d_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Type used to define spline degrees in the four directions
  !---------------------------------------------------------------------------
  type, public :: spline_degree_4d_t

    sll_int32 :: eta1
    sll_int32 :: eta2
    sll_int32 :: eta3
    sll_int32 :: vpar
    
  end type spline_degree_4d_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Type used to define spline degrees in the three directions
  !>  eta1, eta2 and eta3
  !---------------------------------------------------------------------------
  type, public :: spline_degree_3d_t

    sll_int32 :: eta1
    sll_int32 :: eta2
    sll_int32 :: eta3
    
  end type spline_degree_3d_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Type used to define spline degrees in the two directions
  !>  eta1 and eta2 
  !---------------------------------------------------------------------------
  type, public :: spline_degree_2d_t

    sll_int32 :: eta1
    sll_int32 :: eta2
    
  end type spline_degree_2d_t
  !---------------------------------------------------------------------------

end module utils_DK4D_module
