module WENO_interp
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none

  type WENO_interp_1d
     sll_real64                          :: xmin
     sll_real64                          :: xmax
     sll_int32                           :: n_points ! size
     sll_real64                          :: delta    ! discretization step
     sll_real64                          :: rdelta   ! reciprocal of delta
     sll_real64, dimension(:), pointer   :: data     ! data for interpolation
  end type WENO_INTERP_1D
  
  interface delete
     module procedure delete_WENO_1D
  end interface

contains  ! ****************************************************************

  function new_WENO_1d(num_points, xmin, xmax)
    
    type(WENO_interp_1d)             :: new_WENO_1d
    sll_int32,  intent(in)        :: num_points
    sll_real64, intent(in)        :: xmin
    sll_real64, intent(in)        :: xmax
    sll_int32       :: ierr  ! error flag for memory allocation

    new_WENO_1D%xmin     = xmin
    new_WENO_1D%xmax     = xmax
    new_WENO_1D%n_points = num_points
    new_WENO_1D%delta    = (xmax - xmin)/real((num_points-1),f64)
    new_WENO_1D%rdelta   = 1.0_f64/new_WENO_1D%delta
    SLL_ALLOCATE(new_WENO_1D%data(num_points), ierr)
    if( num_points .le. 28 ) then
       print *, 'ERROR, new_WENO_1D: Because of the algorithm used, this function is meant to be used with arrays that are at least of size = 28'
       STOP 'new_WENO_1D()'
    end if
    if( xmin .gt. xmax ) then
       print *, 'ERROR, new_WENO_1D: xmin is greater than xmax, this would cause all sorts of errors.'
       STOP
    end if
  end function new_WENO_1d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  subroutine interpolate_WENO_1D(fWENO, num_points, data, coordinate_o, order, f_out)
    type(WENO_interp_1d)    :: fWENO   ! pointer with basic discretization parameter set up and data
    sll_int32, intent(in)             ::  num_points  ! size of data
    sll_real64, dimension(1:num_points), intent(in)   ::  data        ! data to be interpolated
    sll_real64, dimension(1:num_points), intent(in)   :: coordinate_o ! coordinate for data and interpolation locations
    sll_int32, intent(in)                             :: order ! order of interpolation
    sll_real64, dimension(1:num_points), intent(out)  :: f_out          ! interpolated f at coordinates_o

    sll_real64, dimension(-6: num_points+6)  :: f_temp
    sll_real64                               :: xi_temp
    sll_real64                               :: g1, g2, g3, smo1, smo2, smo3
    sll_real64                               :: w1, w2, w3, w
    sll_int32                                :: i_temp, ic_temp
    sll_real64                               :: eps                

    eps = 1.e-8
    if(.not.(fWENO%n_points.eq.num_points))then
       write(*,*) 'size of f does not agree; check n_f!'
    endif

!!! copy data to f_temp with periodic b.c.
    f_temp(1:num_points) = data(1:num_points)
    do i_temp = 1, order ! periodic b.c.
       f_temp(1-i_temp) = f_temp(num_points+1-i_temp)
       f_temp(num_points+i_temp) = f_temp(i_temp)
    enddo

!!! interpolation based on f_temp
    do ic_temp = 1, num_points
       i_temp  = ceiling((coordinate_o(ic_temp) -fWENO%xmin)*fWENO%rdelta)+1
       xi_temp = (coordinate_o(ic_temp)-(fWENO%xmin+(i_temp-1.0_f64)*fWENO%delta))*fWENO%rdelta
       if(xi_temp.ge.eps)then
          write(*,*) 'something wrong in locating i_temp'
       elseif(xi_temp.le.-1.-eps)then
          write(*,*) 'something wrong in locating i_temp'
       endif

       if(order.eq.2)then

          f_out(ic_temp) = -xi_temp*f_temp(i_temp-1) + (xi_temp+1.0_f64)*f_temp(i_temp)

       elseif(order.eq.4)then

          g1 = f_temp(i_temp-2)*(xi_temp+1.0_f64)*xi_temp/2.0_f64 &
               & - f_temp(i_temp-1)*(xi_temp+2.0_f64)*xi_temp &
               & + f_temp(i_temp)*(xi_temp+2.0_f64)*(xi_temp+1.0_f64)/2.0_f64
          g2 = f_temp(i_temp-1)*(xi_temp)*(xi_temp-1.0)/2.0_f64 &
               &  - f_temp(i_temp)*(xi_temp+1.0_f64)*(xi_temp-1.0_f64) &
               & + f_temp(i_temp+1)*(xi_temp+1.0_f64)*xi_temp/2.0_f64
          smo1 = f_temp(i_temp-2)-2.0_f64*f_temp(i_temp-1) + f_temp(i_temp)
          smo2 = f_temp(i_temp-1)-2.0_f64*f_temp(i_temp) + f_temp(i_temp+1)

          w1 = -(xi_temp-1.0_f64)/3.0_f64
          w2 = (xi_temp+2.0_f64)/3.0_f64
          w1 = w1/(eps+smo1)**2.0_f64
          w2 = w2/(eps+smo2)**2.0_f64
          w = w1 + w2
          w1 = w1/w
          w2 = w2/w

          f_out(ic_temp) = w1*g1 + w2*g2 

       elseif(order.eq.6)then
          g1 = f_temp(i_temp) &
               &  + (11.0_f64/6.0_f64*f_temp(i_temp)-3.0_f64*f_temp(i_temp-1)+1.5_f64*f_temp(i_temp-2)-f_temp(i_temp-3)/3.0_f64)*xi_temp &
               &     + (f_temp(i_temp)-2.5_f64*f_temp(i_temp-1)+2.0_f64*f_temp(i_temp-2)-.5_f64*f_temp(i_temp-3))*xi_temp**2.0_f64 &
               &     + (f_temp(i_temp)/6.0_f64-f_temp(i_temp-1)/2.0_f64+f_temp(i_temp-2)/2.0_f64-f_temp(i_temp-3)/6.0_f64)*xi_temp**3.0_f64
          g2 = f_temp(i_temp) &
               &  + (f_temp(i_temp+1)/3.0_f64+f_temp(i_temp)/2.0_f64-f_temp(i_temp-1)+f_temp(i_temp-2)/6.0_f64)*xi_temp &
               &     + (0.5_f64*f_temp(i_temp+1)-f_temp(i_temp)+.5_f64*f_temp(i_temp-1))*xi_temp**2.0_f64 &
               &     + (f_temp(i_temp+1)/6.0_f64-f_temp(i_temp)/2.0_f64+f_temp(i_temp-1)/2.0_f64-f_temp(i_temp-2)/6.0_f64)*xi_temp**3.0_f64
          g3 = f_temp(i_temp) &
               &  + (-f_temp(i_temp+2)/6.0_f64+f_temp(i_temp+1)-0.5_f64*f_temp(i_temp)-f_temp(i_temp-1)/3.0_f64)*xi_temp &
               &       +(0.5_f64*f_temp(i_temp+1)-f_temp(i_temp)+.5_f64*f_temp(i_temp-1))*xi_temp**2.0_f64 &
               &     + (f_temp(i_temp+2)/6.0_f64-f_temp(i_temp+1)/2.0_f64+f_temp(i_temp)/2.0_f64-f_temp(i_temp-1)/6.0_f64)*xi_temp**3.0_f64

          smo1 = -9.0_f64*f_temp(i_temp-3)*f_temp(i_temp-2) + 4.0_f64/3.0_f64*f_temp(i_temp-3)**2.0_f64 &
               &     - 11.0_f64/3.0_f64*f_temp(i_temp-3)*f_temp(i_temp) + 10.0_f64*f_temp(i_temp-3)*f_temp(i_temp-1) &
               &     +14.0_f64*f_temp(i_temp-2)*f_temp(i_temp) +22.0_f64*f_temp(i_temp-1)**2.0_f64 &
               &     -17.0_f64*f_temp(i_temp-1)*f_temp(i_temp) + 10.0_f64/3.0_f64*f_temp(i_temp)**2.0_f64 &
               &     +16.0_f64*f_temp(i_temp-2)**2.0_f64-37.0_f64*f_temp(i_temp-2)*f_temp(i_temp-1)

          smo2 =  -7.0_f64*f_temp(i_temp-2)*f_temp(i_temp-1) + 4.0_f64/3.0_f64*f_temp(i_temp-2)**2.0_f64 &
               &     - 5.0_f64/3.0_f64*f_temp(i_temp-2)*f_temp(i_temp+1) + 6.0_f64*f_temp(i_temp-2)*f_temp(i_temp) &
               &     +6.0_f64*f_temp(i_temp-1)*f_temp(i_temp+1) +10.0_f64*f_temp(i_temp)**2.0_f64 &
               &     -7.0_f64*f_temp(i_temp)*f_temp(i_temp+1) + 4.0_f64/3.0_f64*f_temp(i_temp+1)**2.0_f64 &
               &     +10.0_f64*f_temp(i_temp-1)**2.0_f64-19.0_f64*f_temp(i_temp-1)*f_temp(i_temp)

          smo3 = -17.0_f64*f_temp(i_temp-1)*f_temp(i_temp) + 10.0_f64/3.0_f64*f_temp(i_temp-1)**2.0_f64 &
               &     - 11.0_f64/3.0_f64*f_temp(i_temp-1)*f_temp(i_temp+2) + 14.0_f64*f_temp(i_temp-1)*f_temp(i_temp+1) &
               &     +10.0_f64*f_temp(i_temp)*f_temp(i_temp+2) +16.0_f64*f_temp(i_temp+1)**2.0_f64 &
               &     -9.0_f64*f_temp(i_temp+1)*f_temp(i_temp+2) + 4.0_f64/3.0_f64*f_temp(i_temp+2)**2.0_f64 &
               &     +22.0_f64*f_temp(i_temp)**2.0_f64-37.0_f64*f_temp(i_temp)*f_temp(i_temp+1)

          w1 = (xi_temp-1.0_f64)*(xi_temp-2.0_f64)/20.0_f64
          w2 = -(xi_temp-2.0_f64)*(xi_temp+3.0_f64)/10.0_f64
          w3 = (xi_temp+2.0_f64)*(xi_temp+3.0_f64)/20.0_f64

          w1 = w1/(eps+smo1)**2.0_f64
          w2 = w2/(eps+smo2)**2.0_f64
          w3 = w3/(eps+smo3)**2.0_f64
          w = w1 + w2+ w3
          w1 = w1/w
          w2 = w2/w
          w3 = w3/w 
          f_out(ic_temp) = w1*g1 + w2*g2 + w3*g3

       else
          print *, 'other order options not available.'
       endif
    enddo

  end subroutine interpolate_WENO_1D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine delete_WENO_1D( fWENO )
    type(WENO_interp_1d) :: fWENO
    sll_int32         :: ierr
   
    SLL_DEALLOCATE(fWENO%data, ierr )
  end subroutine delete_WENO_1D
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module WENO_interp
