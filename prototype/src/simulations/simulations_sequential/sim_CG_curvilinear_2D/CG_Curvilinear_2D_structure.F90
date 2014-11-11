module module_cg_curvi_structure

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "selalib.h"


! type sll_plan_adv_polar
  ! type for advection with center-guide equations
  ! the field and other needed data/object are within
  type sll_plan_adv_curvilinear
     sll_real64 :: eta1_min,eta1_max,delta_eta1
     sll_real64 :: eta2_min,eta2_max,delta_eta2
     sll_real64 :: dt
     sll_int32  :: N_eta1,N_eta2,bc1_type,bc2_type
     sll_int32  :: carac_case
     type(sll_cubic_spline_2D),    pointer     :: spl_f
     sll_real64, dimension(:,:,:), allocatable :: field 
  end type sll_plan_adv_curvilinear

  ! type sll_SL_polar
  ! type for semi Lagrangian
  ! contains other types for the routines called in SL routines
  type sll_SL_curvilinear
     type(sll_plan_adv_curvilinear), pointer     :: adv
     type(plan_curvilinear_op),      pointer     :: grad
     sll_real64, dimension(:,:),     allocatable :: phi
     type(sll_plan_poisson_polar),   pointer     :: poisson
  end type sll_SL_curvilinear

  ! type plan_curvilinear_op
  type plan_curvilinear_op
     sll_int32 :: grad_case
     type(sll_cubic_spline_2D), pointer :: spl_phi
  end type plan_curvilinear_op


end module module_cg_curvi_structure


