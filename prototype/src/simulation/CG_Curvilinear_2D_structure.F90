module module_cg_curvi_structure

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"


  use sll_cubic_splines 
  use sll_poisson_2d_polar
!>type sll_plan_adv_polar
  !>type for advection with center-guide equations
  !>the field and other needed data/object are within
  type sll_plan_adv_curvilinear
     sll_real64 :: eta1_min,eta1_max,delta_eta1
     sll_real64 :: eta2_min,eta2_max,delta_eta2
     sll_int32 :: N_eta1,N_eta2,bc1_type,bc2_type
     sll_real64 :: dt
     type(sll_cubic_spline_2D), pointer :: spl_f
     sll_int32 :: carac_case
     sll_real64, dimension(:,:,:), allocatable :: field 
  end type sll_plan_adv_curvilinear

  !>type sll_SL_polar
  !>type for semi Lagrangian
  !>contains other types for the routines called in SL routines
  type sll_SL_curvilinear
     type(sll_plan_adv_curvilinear), pointer :: adv
     type(plan_curvilinear_op), pointer :: grad
     type(sll_plan_poisson_polar), pointer :: poisson
     sll_real64, dimension(:,:), allocatable :: phi
  end type sll_SL_curvilinear

  !>type plan_curvilinear_op
  type plan_curvilinear_op
     sll_int32 :: grad_case
     type(sll_cubic_spline_2D), pointer :: spl_phi
  end type plan_curvilinear_op


end module module_cg_curvi_structure



!type sll_geom_x 
   !  sll_real64 :: x1_min,x2_min,x1_max,x2_max
    ! sll_real64 :: delta_x1,delta_2
  !end type  sll_geom_x

  !type sll_geom_eta
     !sll_real64 :: eta1_min,eta1_max,delta_eta1
     !sll_real64 :: eta2_min,eta2_max,delta_eta2
     !sll_int32 :: N_eta1,N_eta2
  !end type  sll_geom_x
