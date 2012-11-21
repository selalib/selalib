!!! odd (first, third and fifth) order WENO interpolation with right biased stencils
!!! for example: u = u_j for x \in [x_{j-1}, x_j]

module WENO_interp
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none

  type WENO_interp_1D
     sll_real64                                    :: xmin
     sll_real64                                    :: xmax
     sll_int32                                     :: n_points ! size
     sll_real64                                    :: delta    ! discretization step
     sll_real64                                    :: rdelta   ! reciprocal of delta
     sll_real64, dimension(:), pointer  :: data     ! data for interpolation
     sll_int32                                     :: i_weno ! flag for weno or linear interpolation
     sll_int32                                     :: order   ! order of interpolation
  end type WENO_interp_1D

  interface delete
     module procedure delete_WENO_interp_1D
  end interface delete
  interface new
     module procedure new_WENO_interp_1D
  end interface new

contains  ! ****************************************************************

  function new_WENO_interp_1D(num_points, xmin, xmax, order, i_weno)

    type(WENO_interp_1D), pointer                    :: new_WENO_interp_1D
    sll_int32,  intent(in)                           :: num_points
    sll_real64, intent(in)                           :: xmin
    sll_real64, intent(in)                           :: xmax   
    sll_int32, intent(in)                            :: order
    sll_int32, intent(in)                            :: i_weno
    sll_int32                                        :: ierr  ! allows for an error code return in allocating memory
    
    SLL_ALLOCATE( new_WENO_interp_1D, ierr )

    new_WENO_interp_1D%xmin     = xmin
    new_WENO_interp_1D%xmax     = xmax
    new_WENO_interp_1D%n_points = num_points
    new_WENO_interp_1D%delta    = (xmax - xmin)/real(num_points,f64)
    new_WENO_interp_1D%rdelta   = 1.0_f64/new_WENO_interp_1D%delta
    new_WENO_interp_1D%order = order
    new_WENO_interp_1D%i_weno = i_weno
    SLL_ALLOCATE(new_WENO_interp_1D%data(num_points), ierr)
    if( num_points .le. 28 ) then
       print *, 'ERROR, new_WENO_interp_1D: Because of the algorithm used, this function is meant to be used with arrays that are at least of size = 28'
       STOP 'new_WENO_interp_1D()'
    end if
    if( xmin .gt. xmax ) then
       print *, 'ERROR, new_WENO_interp_1D: xmin is greater than xmax, this would cause all sorts of errors.'
       STOP
    end if
  end function new_WENO_interp_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine interpolate_WENO_1D(fWENO, num_points, data, coordinate_o, f_out)
    type(WENO_interp_1D), pointer         :: fWENO   ! pointer with basic discretization parameter set up and data
    sll_int32, intent(in)                             ::  num_points  ! size of data
    sll_real64, dimension(1:num_points), intent(in)   ::  data        ! data to be interpolated
    sll_real64, dimension(1:num_points), intent(in)   :: coordinate_o ! coordinate for data and interpolation locations
    sll_real64, dimension(1:num_points), intent(out)  :: f_out          ! interpolated f at coordinates_o

    sll_real64, dimension(-6: num_points+6)         :: f_temp
    sll_real64                                      :: xi_temp
    sll_real64                                      :: g1, g2, g3, smo1, smo2, smo3
    sll_real64                                      :: w1, w2, w3, w
    sll_int32                                       :: i_temp, ic_temp, i_weno, order
    sll_real64                                      :: eps                

    eps = 1.d-8

    if(.not.(fWENO%n_points.eq.num_points))then
       write(*,*) 'size of f does not agree; check n_f!'
    endif

    order = fWENO%order
    i_weno = fWENO%i_weno

!!! copy data to f_temp with periodic b.c.
    f_temp(1:num_points) = data(1:num_points)
    do i_temp = 1, order ! periodic b.c.
       f_temp(1-i_temp) = f_temp(num_points+1-i_temp)
       f_temp(num_points+i_temp) = f_temp(i_temp)
    enddo

!!! interpolation based on f_temp
    do ic_temp = 1, num_points
       i_temp  = ceiling((coordinate_o(ic_temp) -fWENO%xmin-fWENO%delta/2.0_f64)*fWENO%rdelta)+1
       xi_temp = (coordinate_o(ic_temp)-(fWENO%xmin+fWENO%delta/2.0_f64+(i_temp-1.0_f64)*fWENO%delta))*fWENO%rdelta
        if(xi_temp.ge.eps)then
          write(*,*) 'something wrong in locating i_temp'
       elseif(xi_temp.le.-1.-eps)then
          write(*,*) 'something wrong in locating i_temp'
       endif

       if(order.eq.1)then

          f_out(ic_temp) = f_temp(i_temp)

       elseif(order.eq.3)then

          g1 = f_temp(i_temp-1)*(-xi_temp) &
               & + f_temp(i_temp)*(xi_temp+1.0_f64) 
          g2 = f_temp(i_temp)*(-xi_temp+1.0_f64)  &
               & + f_temp(i_temp+1)*xi_temp

          smo1 = (f_temp(i_temp-1)-f_temp(i_temp))**2.0_f64
          smo2 = (f_temp(i_temp)-f_temp(i_temp+1))**2.0_f64

          w1 = (-xi_temp+1.0_f64)/2.0_f64
          w2 = (xi_temp+1.0_f64)/2.0_f64
          if(i_weno.eq.1)then
             w1 = w1/(eps+smo1)**2.0_f64
             w2 = w2/(eps+smo2)**2.0_f64
          endif
          w = w1 + w2
          w1 = w1/w
          w2 = w2/w

          f_out(ic_temp) = w1*g1 + w2*g2 

       elseif(order.eq.5)then
          g1 = f_temp(i_temp-2)*(xi_temp+1.0_f64)*xi_temp/2.0_f64 &
               &       - f_temp(i_temp-1)*(xi_temp+2.0_f64)*xi_temp &
               &       + f_temp(i_temp)*(xi_temp+2.0_f64)*(xi_temp+1.0_f64)/2.0_f64 

          g2 = f_temp(i_temp-1)*xi_temp*(xi_temp-1.0_f64)/2.0_f64 &
               &       - f_temp(i_temp)*(xi_temp+1.0_f64)*(xi_temp-1.0_f64) &
               &       + f_temp(i_temp+1)*(xi_temp+1.0_f64)*xi_temp/2.0_f64 

          g3 = f_temp(i_temp)*(xi_temp-1.0_f64)*(xi_temp-2.0_f64)/2.0_f64 &
               &       - f_temp(i_temp+1)*(xi_temp-2.0_f64)*xi_temp &
               &       + f_temp(i_temp+2)*(xi_temp-1.0_f64)*xi_temp/2.0_f64 


          smo1 = 13.0_f64/12.0_f64*(f_temp(i_temp-2)-2.0_f64*f_temp(i_temp-1)+f_temp(i_temp))**2.0_f64 &
               &         + (-f_temp(i_temp-1)+f_temp(i_temp))**2.0_f64             

          smo2 = 13.0_f64/12.0_f64*(f_temp(i_temp-1)-2.0_f64*f_temp(i_temp)+f_temp(i_temp+1))**2.0_f64 &
               &         + (f_temp(i_temp-1)-f_temp(i_temp))**2.0_f64    

          smo3 = 13.0_f64/12.0_f64*(f_temp(i_temp)-2.0_f64*f_temp(i_temp+1)+f_temp(i_temp+2))**2.0_f64 &
               &         + (-2.0_f64*f_temp(i_temp) + 3.0_f64*f_temp(i_temp+1) - f_temp(i_temp+2))**2.0_f64

          w1 = (xi_temp-1.0_f64)*(xi_temp-2.0_f64)/12.0_f64
          w3 = (xi_temp+2.0_f64)*(xi_temp+1.0_f64)/12.0_f64
          w2 = 1.0_f64 - w1 - w3

          if(i_weno.eq.1)then
             w1 = w1/(eps+smo1)**2.0_f64
             w2 = w2/(eps+smo2)**2.0_f64
             w3 = w3/(eps+smo3)**2.0_f64
          endif
          w = w1 + w2+ w3
          w1 = w1/w
          w2 = w2/w
          w3 = w3/w       
          f_out(ic_temp) = w1*g1 + w2*g2 + w3*g3   

       else
          print *, 'interpolate_WENO_1D: other order options not available.'
          stop
       endif
    enddo

  end subroutine interpolate_WENO_1D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine delete_WENO_interp_1D( fWENO )
    type(WENO_interp_1D), pointer :: fWENO
    sll_int32                    :: ierr
    ! Fixme: some error checking, whether the spline pointer is associated
    ! for instance
    SLL_ASSERT( associated(fWENO) )
    fWENO%data => null()
    SLL_DEALLOCATE(fWENO, ierr )
  end subroutine delete_WENO_interp_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module WENO_interp
