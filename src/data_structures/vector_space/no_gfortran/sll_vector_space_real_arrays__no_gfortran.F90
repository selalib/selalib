module sll_m_vector_space_real_arrays

#include "sll_working_precision.h"
   use sll_m_vector_space_base, only: sll_vector_space_base

   implicit none

   public :: &
      sll_vector_space_real_1d, &
      sll_vector_space_real_2d, &
      sll_vector_space_real_3d

   private

   !============================================================================
   type, extends(sll_vector_space_base) :: sll_vector_space_real_1d(n)

   sll_int32, len :: n
   sll_real64     :: array(n)

contains
   procedure      :: copy => copy__r1d
   procedure      :: incr => incr__r1d
   procedure      :: scal => scal__r1d

   procedure      :: norm => norm__r1d
   procedure      :: inner => inner__r1d

   procedure, private :: source_single => source_single__r1d
   procedure, private :: source_array => source_array__r1d

   end type sll_vector_space_real_1d

end module sll_vector_space_real_arrays
