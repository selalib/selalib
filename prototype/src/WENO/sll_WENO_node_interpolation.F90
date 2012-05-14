module sll_WENO
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_interpolator_1d_base
  implicit none
  type, extends(sll_interpolator_1d_base)  :: sll_weno_1d
     sll_real64                        :: xmin
     sll_real64                        :: xmax
     sll_int32                         :: n_points ! size
     sll_real64                        :: delta    ! discretization step
     sll_real64                        :: rdelta   ! reciprocal of delta
     sll_real64, dimension(:), pointer :: data     ! data for interpolation
   contains
     procedure, pass :: interpolate_array => interpolate_WENO_1D  
     procedure, pass :: reconstruct_array => reconstruct_WENO_1D
  end type sll_weno_1d

  interface delete
     module procedure delete_weno_1d
  end interface delete

contains  ! ****************************************************************

  function new_WENO_1D(num_points, xmin, xmax)

    type(sll_WENO_1D)           :: new_WENO_1D
    sll_int32,  intent(in)      :: num_points
    sll_real64, intent(in)      :: xmin
    sll_real64, intent(in)      :: xmax
    sll_int32                   :: ierr  ! allows for an error code return in allocating memory
    sll_int32                   :: i_temp
    !SLL_ALLOCATE( new_WENO_1D, ierr )
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
  end function new_WENO_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function interpolate_WENO_1D(this, num_points, data, coordinates) result(res)
    class(sll_WENO_1D), intent(in)  :: this        ! basic discretization parameter set up and data
    sll_int32, intent(in)                   :: num_points  ! size of data
    sll_real64, dimension(:), intent(in)    :: data        ! data to be interpolated
    sll_real64, dimension(:), intent(in)    :: coordinates ! coordinate for data and interpolation locations
    sll_real64, dimension(num_points)       :: res         ! interpolated f at coordinates

    sll_real64, dimension(-6: num_points+6)         :: f_temp
    sll_real64                                      :: xi_temp
    sll_real64                                      :: g1, g2, g3, smo1, smo2, smo3
    sll_real64                                      :: w1, w2, w3, ww1, ww2, ww3, w
    sll_int32                                       :: i_temp, ic_temp
    sll_real64                                      :: eps                
    sll_int32                                       :: order

    eps = 1.e-8
    order = 6
    SLL_ASSERT(this%n_points.eq.num_points)

!!! copy data to f_temp with periodic b.c.
!!!ES modified Jingmei's version so that data are at grid points and not at cell centers so that points 1 and num_points are same point
    f_temp(1:num_points) = data(1:num_points)
    do i_temp = 1, order ! periodic b.c.
       f_temp(1-i_temp) = f_temp(num_points-i_temp)
       f_temp(num_points+i_temp) = f_temp(1+i_temp)
    enddo

!!! interpolation based on f_temp
    do ic_temp = 1, num_points
       i_temp  = ceiling((coordinates(ic_temp) -this%xmin)*this%rdelta)+1
       xi_temp = (coordinates(ic_temp)-(this%xmin+(i_temp-1.0_f64)*this%delta))*this%rdelta
       if(xi_temp.ge.eps)then
          write(*,*) 'something wrong in locating i_temp'
       elseif(xi_temp.le.-1.-eps)then
          write(*,*) 'something wrong in locating i_temp'
       endif
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
       w = w1 + w2 + w3
       w1 = w1/w
       w2 = w2/w
       w3 = w3/w
       res(ic_temp) = w1*g1 + w2*g2 + w3*g3
    enddo

  end function interpolate_WENO_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function reconstruct_WENO_1D(this, num_points, data) result(res)
    use sll_working_precision
    class(sll_WENO_1D), intent(in)     :: this
    sll_int32, intent(in)                 :: num_points    ! size of output array
    sll_real64, dimension(:), intent(in)  :: data          ! data to be interpolated 
    sll_real64, dimension(num_points)         :: res
    !! dummy function

  end function reconstruct_WENO_1D

  subroutine delete_WENO_1D( fWENO )
    type(sll_WENO_1D) :: fWENO
    sll_int32                    :: ierr
    ! Fixme: some error checking, whether the spline pointer is associated
    ! for instance
    !SLL_ASSERT( associated(fWENO) )
    fWENO%data => null()
    !SLL_DEALLOCATE(fWENO, ierr )
  end subroutine delete_WENO_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sll_WENO
