!>@internal [example]
program test_io_gnuplot
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_pi

  use sll_m_gnuplot, only: &
    sll_gnuplot_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_int32 :: i, j     
sll_int32 :: error 
sll_int32 :: n1, n2
sll_int32 :: iplot

sll_real64, allocatable, dimension(:)   :: eta1
sll_real64, allocatable, dimension(:)   :: eta2
sll_real64, allocatable, dimension(:,:) :: x1
sll_real64, allocatable, dimension(:,:) :: x2
sll_real64, allocatable, dimension(:,:) :: array

sll_real64 :: eta1_min = 1.0_f64, eta1_max = 2.0_f64
sll_real64 :: eta2_min = 0.0_f64, eta2_max = 2*sll_pi

n1 = 32
n2 = 64

SLL_ALLOCATE(eta1(n1),error)
SLL_ALLOCATE(eta2(n2),error)
SLL_ALLOCATE(x1(n1,n2),error)
SLL_ALLOCATE(x2(n1,n2),error)

!Construct the mesh
do j = 1, n2
   eta2(j) = (j-1) * (eta2_max-eta2_min)/(n2-1)
   do i = 1, n1
      eta1(i) =  eta1_min + (i-1) * (eta1_max-eta1_min) / real(n1-1,f64)
      x1(i,j) = eta1(i) * cos(eta2(j))
      x2(i,j) = eta1(i) * sin(eta2(j))
   end do
end do 

SLL_ALLOCATE(array(n1,n2),error)

array = cos(2.*x1)*exp(-x2*x2)

! Set plot index to 1, otherwise a new file id is not computed (dangerous!)
iplot = 1

call sll_gnuplot_2d(eta1_min, eta1_max, n1, &
                    eta2_min, eta2_max, n2, &
                    array, 'plot_1', iplot, error)

call sll_gnuplot_2d(n1, eta1, n2, eta2, array, 'plot_2', iplot, error)
 
call sll_gnuplot_2d(n1, n2, x1, x2, array, 'plot_3', iplot, error)
 
call sll_gnuplot_2d(n1, n2, x1, x2, 'mesh', error)
 
print*, 'PASSED'
end program test_io_gnuplot
!>@internal [example]
