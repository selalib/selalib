!> sll_m_gaussian random generator
!>
!>\author 
!>\date created: 
module sll_m_gaussian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_f_gaussian_deviate, &
    sll_s_gaussian_deviate_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains 
  !>\author Sever Hirstoaga
  !>\date created: 2012-05-31

  !> Returns a value of a random variable in ONE dimension
  !> following the Gaussian probability density
  !> with zero mean and unit variance, using random generator
  !> for the uniform deviates.
  function sll_f_gaussian_deviate()
    sll_real64 :: sll_f_gaussian_deviate

    sll_real64 :: rsq, v1, v2
    sll_real64, save :: g
    logical, save :: g_stored=.false.
    
    if (g_stored) then
       sll_f_gaussian_deviate=g
       g_stored=.false.
    else
       do
          call random_number(v1)
          call random_number(v2)
          v1=2._f64*v1 - 1._f64
          v2=2._f64*v2 - 1._f64
          rsq = v1**2+v2**2
          if (rsq>0._f64.and.rsq<1._f64) exit
       enddo
       rsq = sqrt(-2._f64*log(rsq)/rsq)
       sll_f_gaussian_deviate=v1*rsq
       g=v2*rsq
       g_stored=.true.
    endif
  end function sll_f_gaussian_deviate

  
  !> Returns a value of a random variable in TWO dimensions
  !> following the Gaussian probability density
  !> with zero mean and unit variance, using random generator
  !> for the uniform deviates.
  subroutine sll_s_gaussian_deviate_2d(res)
    implicit none
    sll_real64,intent(out) :: res(1:2)
    sll_real64 :: rsq, v1, v2
    
    do
       call random_number(v1)
       call random_number(v2)
       v1=2._f64*v1 - 1._f64
       v2=2._f64*v2 - 1._f64
       rsq = v1**2+v2**2
       if (rsq>0._f64.and.rsq<1._f64) exit
    enddo
    rsq = sqrt(-2._f64*log(rsq)/rsq)
    res(1) = v1*rsq
    res(2) = v2*rsq
  end subroutine sll_s_gaussian_deviate_2d

end module sll_m_gaussian
