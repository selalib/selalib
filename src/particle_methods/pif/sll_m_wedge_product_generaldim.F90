!**************************************************************
!  Author: Jakob Ameres, jakob.ameres@tum.de
!**************************************************************
!This module should provide the wedge product for arbitrary dimensions
!such that it is possible to define Vlasov equation with a kind of
!v x B term for any dimension
module sll_m_wedge_product_generaldim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_f_cross_product_2d, &
    sll_f_cross_product_3d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 contains

!This is just an arbitrary definition
pure function sll_f_cross_product_2d(v,w) result(cross)
sll_real64, dimension(:,:), intent(in) :: v,w
 sll_real64, dimension(2,size(v,2)):: cross
  sll_real64, dimension(size(v,2)):: determinante
  
 !! SLL_ASSERT(size(v,1)==2)
  
   determinante(:)=v(1,:)*w(2,:) - v(2,:)*w(1,:)
!   cross=v*(v(1,:)*w(2,:) - v(2,:)*w(1,:))
  cross(1,:)=v(1,:)*determinante
  cross(2,:)=v(2,:)*determinante
end function


pure function sll_f_cross_product_3d( v, w) result(cross)
 sll_real64, dimension(:,:), intent(in) :: v,w
 sll_real64, dimension(3,size(v,2)):: cross
 !! SLL_ASSERT(size(v,1)==3)
  
  cross(1,:) = v(2,:) * w(3,:) - v(3,:) * w(2,:)
  cross(2,:) = v(3,:) * w(1,:) - v(1,:) * w(3,:)
  cross(3,:) = v(1,:) * w(2,:) - v(2,:) * w(1,:)
end function






end module