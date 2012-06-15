module WENO_recon
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none

  type WENO_recon_1D
     sll_real64                :: xmin
     sll_real64                :: xmax
     sll_int32                 :: n_points ! size
     sll_real64                :: delta    ! discretization step
     sll_real64                :: rdelta   ! reciprocal of delta
     sll_real64, dimension(:), pointer  :: data     ! data for interpolation
  end type WENO_recon_1D
  
  interface delete
     module procedure delete_WENO_recon_1D
  end interface

contains  ! ****************************************************************

  function new_WENO_recon_1D(num_points, xmin, xmax)
    
    type(WENO_recon_1D)            :: new_WENO_recon_1D
    sll_int32,  intent(in)         :: num_points
    sll_real64, intent(in)         :: xmin
    sll_real64, intent(in)         :: xmax
  
    new_WENO_recon_1D%xmin     = xmin
    new_WENO_recon_1D%xmax     = xmax
    new_WENO_recon_1D%n_points = num_points
    new_WENO_recon_1D%delta    = (xmax - xmin)/real((num_points),f64)
    new_WENO_recon_1D%rdelta   = 1.0_f64/new_WENO_recon_1D%delta
    
    if( num_points .le. 28 ) then
       print *, 'ERROR, new_WENO_recon_1D: Because of the algorithm used,', &
            'this function is meant to be used with arrays that are at least',&
            ' of size = 28'
       STOP 'new_WENO_recon_1D()'
    end if
    if( xmin .gt. xmax ) then
       print *, 'ERROR, new_WENO_recon_1D: xmin is greater than xmax,',&
            ' this would cause all sorts of errors.'
       STOP
    end if
  end function new_WENO_recon_1D
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  subroutine FD_WENO_recon_1D(fWENO, num_points, data, order, deri)
    type(WENO_recon_1D)      :: fWENO  
    sll_int32, intent(in)    ::  num_points  ! size of data
    sll_int32, intent(in)    ::  order
    ! point values of function
    sll_real64, dimension(1:num_points), intent(in)    ::  data 
    ! derivatives of f 
    sll_real64, dimension(1:num_points), intent(out)   ::  deri   
    ! local variables
    sll_real64, dimension(-6: num_points+6)   :: f_temp
    sll_real64, dimension(0:num_points)       :: flux, flux_l, flux_r
    sll_real64                                :: g1, g2, g3, smo1, smo2, smo3
    sll_real64                                :: w1, w2, w3, w
    sll_real64                                :: eps 
    sll_int32                                 :: ic_temp               

    eps = 1.e-8

    if(.not.(fWENO%n_points.eq.num_points))then
       write(*,*) 'size of f does not agree; check n_f!'
    endif

!!! copy data to f_temp with periodic b.c.
    f_temp(1:num_points) = data(1:num_points)

    ! periodic b.c.
    do ic_temp = 1, order 
       f_temp(1-ic_temp) = f_temp(num_points+1-ic_temp)
       f_temp(num_points+ic_temp) = f_temp(ic_temp)
    enddo

    ! WENO reconstruction
    ! 2nd order    
    if(order.eq.2)then
       do ic_temp = 1, num_points
          flux_l(ic_temp) = f_temp(ic_temp)
          flux_r(ic_temp) = f_temp(ic_temp+1)
       enddo

       ! fourth order
    elseif(order.eq.4)then

       do ic_temp = 1, num_points
          g1 = -0.5_f64*f_temp(ic_temp-1) + 1.5_f64*f_temp(ic_temp)
          g2 = 0.5_f64*f_temp(ic_temp) + 0.5_f64*f_temp(ic_temp+1)
          smo1 = (f_temp(ic_temp-1) - f_temp(ic_temp))**2.0_f64
          smo2 = (f_temp(ic_temp) - f_temp(ic_temp+1))**2.0_f64
          w1 = 1.0_f64/3.0_f64
          w2 = 2.0_f64/3.0_f64
          w  = w1+w2
          w1 = w1/w
          w2 = 1.0_f64-w1
          flux_l(ic_temp) = g1*w1+g2*w2

          g1 = 0.5_f64*f_temp(ic_temp) + 0.5_f64*f_temp(ic_temp+1)
          g2 = 1.5_f64*f_temp(ic_temp+1) - 0.5_f64*f_temp(ic_temp+2)
          smo1 =(f_temp(ic_temp) - f_temp(ic_temp+1))**2.0_f64
          smo2 =(f_temp(ic_temp+2) - f_temp(ic_temp+1))**2.0_f64 
          w1 = 2.0_f64/3.0_f64
          w2 = 1.0_f64/3.0_f64
          w  = w1+w2
          w1 = w1/w
          w2 = 1.0_f64-w1
          flux_r(ic_temp) = g1*w1+g2*w2    
       enddo

       ! sixth order   
    elseif(order.eq.6)then
       do ic_temp = 1, num_points
          w1=.1_f64
          w2=.6_f64
          w3=.3_f64
          g1=1.0_f64/3.0_f64*f_temp(ic_temp-2) &
               -7.0_f64/6.0_f64*f_temp(ic_temp-1) &
               +11.0_f64/6.0_f64*f_temp(ic_temp)
          smo1=13.0_f64/12.0_f64*(f_temp(ic_temp-2) &
               -2.0_f64*f_temp(ic_temp-1)+f_temp(ic_temp))**2.0_f64 &
               +1.0_f64/4.0_f64*(f_temp(ic_temp-2) &
               -4.0_f64*f_temp(ic_temp-1)+3.0_f64*f_temp(ic_temp))**2.0_f64

          g2=-1.0_f64/6.0_f64*f_temp(ic_temp-1) &
               +5.0_f64/6.0_f64*f_temp(ic_temp) &
               +1.0_f64/3.0_f64*f_temp(ic_temp+1)
          smo2=13.0_f64/12.0_f64*(f_temp(ic_temp-1)-2.0_f64*f_temp(ic_temp) &
               +f_temp(ic_temp+1))**2.0_f64 &
               +1.0_f64/4.0_f64*(f_temp(ic_temp+1) &
               -f_temp(ic_temp-1))**2.0_f64

          g3=1.0_f64/3.0_f64*f_temp(ic_temp) &
               +5.0_f64/6.0_f64*f_temp(ic_temp+1) &
               -1.0_f64/6.0_f64*f_temp(ic_temp+2)
          smo3=13.0_f64/12.0_f64*(f_temp(ic_temp) &
               -2.0_f64*f_temp(ic_temp+1)+f_temp(ic_temp+2))**2.0_f64 &
               +1.0_f64/4.0_f64*(3.0_f64*f_temp(ic_temp) &
               -4.0_f64*f_temp(ic_temp+1)+f_temp(ic_temp+2))**2.0_f64

          w1=w1/(eps+smo1)**2.0_f64
          w2=w2/(eps+smo2)**2.0_f64
          w3=w3/(eps+smo3)**2.0_f64
          w=(w1+w2+w3)
          w1=w1/w
          w2=w2/w
          w3=1-w1-w2
          flux_l(ic_temp)=w1*g1+w2*g2+w3*g3

          w1=.3
          w2=.6
          w3=.1
          g1=-1.0_f64/6.0_f64*f_temp(ic_temp-1) &
               +5.0_f64/6.0_f64*f_temp(ic_temp) &
               +1.0_f64/3.0_f64*f_temp(ic_temp+1)
          smo1=13.0_f64/12.0_f64*(f_temp(ic_temp-1) &
               -2.0_f64*f_temp(ic_temp)+f_temp(ic_temp+1))**2.0_f64 &
                  +1.0_f64/4.0_f64*(f_temp(ic_temp-1) &
                  -4.0_f64*f_temp(ic_temp)+3.0_f64*f_temp(ic_temp+1))**2.0_f64
          g2=1.0_f64/3*f_temp(ic_temp)+5.0_f64/6.0_f64*f_temp(ic_temp+1) &
               -1.0_f64/6.0_f64*f_temp(ic_temp+2) 
          smo2=13.0_f64/12.0_f64*(f_temp(ic_temp)-2*f_temp(ic_temp+1) &
               +f_temp(ic_temp+2))**2.0_f64 &
               +1.0_f64/4.0_f64*(f_temp(ic_temp)-f_temp(ic_temp+2))**2.0_f64
          g3=11.0_f64/6.0_f64*f_temp(ic_temp+1) &
               -7.0_f64/6.0_f64*f_temp(ic_temp+2) &
               +1.0_f64/3.0_f64*f_temp(ic_temp+3)
          smo3=13.0_f64/12.0_f64*(f_temp(ic_temp+1) &
               -2.0_f64*f_temp(ic_temp+2)+f_temp(ic_temp+3))**2.0_f64 &
               &    +1.0_f64/4.0_f64*(3.0_f64*f_temp(ic_temp+1)  &
               -4.0_f64*f_temp(ic_temp+2)+f_temp(ic_temp+3))**2.0_f64
          w1 = w1/(eps+smo1)**2.0_f64
          w2 = w2/(eps+smo2)**2.0_f64
          w3 = w3/(eps+smo3)**2.0_f64
          w = w1 + w2+ w3
          w1 = w1/w
          w2 = w2/w
          w3 = w3/w 
          flux_r(ic_temp) = w1*g1 + w2*g2 + w3*g3
       enddo
    else
       write(*,*) 'other order options not available.'
    endif
    do ic_temp=1, num_points
       flux(ic_temp) = (flux_l(ic_temp) + flux_r(ic_temp))*0.5_f64
    enddo
    flux(0) = flux(num_points)
    do ic_temp=1, num_points
       deri(ic_temp) = (flux(ic_temp)-flux(ic_temp-1))*fWENO%rdelta
    enddo

  end subroutine FD_WENO_recon_1D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine delete_WENO_recon_1D( fWENO )
    type(WENO_recon_1D)   :: fWENO
    sll_int32             :: ierr
    
    if (associated(fWENO%data)) then
       SLL_DEALLOCATE(fWENO%data, ierr )
    end if
  end subroutine delete_WENO_recon_1D
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module WENO_recon
