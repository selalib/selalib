!**************************************************************
!  Author: Jakob Ameres, jakob.ameres@tum.de
!**************************************************************

!Transform samples of random distributions such that certain moments are
!exactly matched
module sll_m_moment_matching
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use sll_m_constants


implicit none

 !> @brief !Match mean E[X] and E[X^2]  = first and second order moment
  interface match_moment_1D_linear
     module procedure match_moment_1D_linear_real64
  end interface

  
  !> @brief !Match mean E[W X] and E[W X^2]  = first and second order moment
  !>       modifies only X,  montecarlo estimate is: sum(W X)/N, N=len(X)
   interface match_moment_1D_weight_linear
      module procedure match_moment_1D_weight_linear_real64
   end interface

contains


!Match mean E[X] and E[X^2]  = first and second order moment
subroutine match_moment_1D_linear_real64( x, mean, mom2 )
sll_real64, dimension(:), intent(inout) :: x !random samples
sll_real64,  intent(in) :: mean, mom2 !desired first and second order moment
sll_real64 :: meanx, mom2x
sll_int32 :: num

num=size(x)
meanx=sum(x)/real(num,f64)
mom2x=sum(x**2)/real(num,f64)

x= (x-meanx)*sqrt(  (mom2 - mean**2)/(mom2x-meanx**2)) + mean
end subroutine match_moment_1D_linear_real64



subroutine match_moment_1D_weight_linear_real64( x, w, mean, mom2, num)
sll_real64, dimension(:), intent(inout) :: x !random samples
sll_real64, dimension(:), intent(in) :: w !weights
sll_int32, intent(in) :: num    !number of samples, doestn have to coincide with size of x
sll_real64,  intent(in) :: mean, mom2 !desired first and second order moment
sll_real64 :: meanx, mom2x, sumw

meanx=sum(w*x)/real(num,f64)
mom2x=sum(w*x**2)/real(num,f64)
sumw=sum(w)/real(num,f64)

x= (x- meanx/sumw)*sqrt((sumw*mom2 - mean**2)/(sumw*mom2x-meanx**2)) + mean/sumw
end subroutine match_moment_1D_weight_linear_real64
end module sll_m_moment_matching
