module sll_m_hermite_aligned_interpolation_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_f_new_hermite_aligned_interpolation_2d, &
    sll_t_hermite_aligned_interpolation_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type :: sll_t_hermite_aligned_interpolation_2d
  sll_real64 :: eta_min(2) !< eta1 min and eta2 min
  sll_real64 :: eta_max(2) !< eta1 max and eta2 max

end type sll_t_hermite_aligned_interpolation_2d 

contains

function sll_f_new_hermite_aligned_interpolation_2d() &
  result(interp)
    
  type(sll_t_hermite_aligned_interpolation_2d), pointer :: interp
  sll_int32 :: ierr
  
  SLL_ALLOCATE(interp,ierr)
  call initialize_hermite_aligned_interpolation_2d( interp )

end function sll_f_new_hermite_aligned_interpolation_2d

subroutine initialize_hermite_aligned_interpolation_2d( interp )
  type(sll_t_hermite_aligned_interpolation_2d) :: interp
! sll_int32 :: r1 
! sll_int32 :: s1 
! sll_int32 :: r2 
! sll_int32 :: s2
! sll_int32 :: num_cells_x1 
! sll_int32 :: num_cells_x2 
  return
  SLL_ASSERT(interp%eta_max(1)>=interp%eta_min(1))
  SLL_ASSERT(interp%eta_max(2)>=interp%eta_min(2))
   
end subroutine initialize_hermite_aligned_interpolation_2d 

end module sll_m_hermite_aligned_interpolation_2d
