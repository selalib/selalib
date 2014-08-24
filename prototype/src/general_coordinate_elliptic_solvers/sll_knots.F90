module sll_knots
#include "sll_working_precision.h"


  implicit none

  sll_int32, parameter :: KNOTS_PERIODIC = 0, KNOTS_DIRICHLET = 1

contains

  subroutine initialize_knots_source(&
       spline_degree, &
       num_cells, &
       eta_min, &
       eta_max, &
       knots)

    sll_int32,  intent(in)   :: spline_degree 
    sll_int32,  intent(in)   :: num_cells
    sll_real64, intent(in)   :: eta_min
    sll_real64, intent(in)   :: eta_max
    sll_real64, dimension(:) :: knots
    sll_int32                :: i
    sll_real64               :: delta

    delta = (eta_max - eta_min)/num_cells

    knots(1:spline_degree+1) = eta_min
    knots(num_cells+2:num_cells+1+spline_degree+1) = eta_max
    
    if (mod(spline_degree + 1, 2) == 0) then
       do i = spline_degree + 2, num_cells + 1
          knots(i) = eta_min +  &
               (i - (spline_degree +1)/2 - 1) * delta  
       end do
    else
       do i = spline_degree + 2, num_cells + 1
          knots(i) = 0.5_f64*(eta_min + (i - (spline_degree)/2 - 1) * delta + &
               eta_min +  ( i -1 - (spline_degree)/2 -1) * delta )
       end do
    end if
  end subroutine initialize_knots_source



  ! SUBROUTINE ASSEMBLING KNOTS ARRAY
  subroutine initialize_knots(&
       spline_degree,&
       num_cells, &
       eta_min, &
       eta_max, &
       bc_left, &
       bc_right, &
       knots)
    
    sll_int32,  intent(in)   :: spline_degree 
    sll_int32,  intent(in)   :: num_cells
    sll_real64, intent(in)   :: eta_min
    sll_real64, intent(in)   :: eta_max
    sll_int32,  intent(in)   :: bc_left
    sll_int32,  intent(in)   :: bc_right
    sll_real64, dimension(:) :: knots
    sll_int32                :: i
    sll_real64               :: eta
    sll_real64               :: delta

    ! This routine needs some error-checking...
    delta = (eta_max - eta_min)/num_cells

    if ( (bc_left  == KNOTS_PERIODIC) .and. &
         (bc_right == KNOTS_PERIODIC) ) then 
       ! it is intentional that eta_min is not used, one intends to consider
       ! only the [0,1] interval...
       do i = -spline_degree, spline_degree+1
          knots(i+spline_degree+1) = delta*i 
       end do
    else if ( (bc_left  == KNOTS_DIRICHLET) .and. &
              (bc_right == KNOTS_DIRICHLET) ) then 
       do i = 1, spline_degree + 1
          knots(i) = eta_min
       enddo
       eta = eta_min
       do i = spline_degree + 2, num_cells + 1 + spline_degree
          eta = eta + delta
          knots(i) = eta
       enddo
       do i = num_cells + spline_degree + 1, num_cells + 1 + 2*spline_degree
          knots(i) = eta_max
       enddo
    end if

  end subroutine initialize_knots
  
end module sll_knots
