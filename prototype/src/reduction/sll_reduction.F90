module sll_reduction_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"


use sll_logical_meshes
implicit none

  abstract interface
     function sll_integration_discrete_1d( data, Npts, delta, params ) result(res)
       use sll_working_precision
       sll_real64                                  :: res
       sll_real64,dimension(:), intent(in)         :: data
       sll_int32, intent(in)                       :: Npts
       sll_real64,intent(in)                       :: delta
       sll_real64, dimension(:), intent(in), optional :: params
     end function sll_integration_discrete_1d
  end interface



contains

  !--------------------------------------------------
  ! Generic function which can be used for computing charge density
  ! by default we suppose a uniform mesh node based
  ! and trapezoid quadrature formula is used
  ! as an option we can use a way to integrate
  ! by providing a function whose signature is sll_integration_discrete_1d
  !---------------------------------------------------  
    

  subroutine compute_reduction_4d_to_3d_direction4(&
    data_4d, &
    data_3d, &
    Npts1, &
    Npts2, &
    Npts3, &
    Npts4, &
    delta4, &    
    integration_func, &
    integration_func_params)
    
    sll_real64, dimension(:,:,:,:), intent(in)    :: data_4d
    sll_real64, dimension(:,:,:)  , intent(out) :: data_3d
    sll_int32, intent(in)  :: Npts1
    sll_int32, intent(in)  :: Npts2
    sll_int32, intent(in)  :: Npts3
    sll_int32, intent(in)  :: Npts4
    sll_real64, intent(in) :: delta4
    procedure(sll_integration_discrete_1d), optional :: integration_func
    sll_real64, dimension(:), optional :: integration_func_params
    sll_int32  :: i1, i2, i3, i4
    sll_real64 :: tmp 
    
    
    if(Npts1>size(data_4d,1))then
      print *,'#Problem for size1 in compute_reduction_4d_to_3d_direction4'
      stop
    endif
    if(Npts2>size(data_4d,2))then
      print *,'#Problem for size2 in compute_reduction_4d_to_3d_direction4'
      stop
    endif
    if(Npts3>size(data_4d,3))then
      print *,'#Problem for size3 in compute_reduction_4d_to_3d_direction4'
      stop
    endif
    if(Npts4>size(data_4d,4))then
      print *,'#Problem for size3 in compute_reduction_4d_to_3d_direction4'
      stop
    endif
    if(.not.(present(integration_func)))then
      do i3 = 1,Npts3
        do i2 = 1,Npts2
          do i1 = 1,Npts1
            tmp = 0.5_f64*(data_4d(i1,i2,i3,1)&
              +data_4d(i1,i2,i3,Npts4))
            do i4 = 2,Npts4-1
              tmp = tmp + data_4d(i1,i2,i3,i4)
            end do
            data_3d(i1,i2,i3) = tmp*delta4
          end do
        end do
      end do
    else
      do i3 = 1,Npts3
        do i2 = 1,Npts2
          do i1 = 1,Npts1
            data_3d(i1,i2,i3) = integration_func( &
              data_4d(i1,i2,i3,:), &
              Npts4, &
              delta4, &
              integration_func_params )
          end do
        end do
      end do      
    endif  
  end subroutine compute_reduction_4d_to_3d_direction4
  
  function compute_integral_trapezoid_1d(data, Npts, delta, func_params) result(res)
    sll_real64, dimension(:), intent(in)    :: data
    sll_int32, intent(in) :: Npts
    sll_real64,intent(in) :: delta
    sll_real64, dimension(:), intent(in) ,optional :: func_params
    sll_real64 :: res
    sll_int32 :: i
    
    res = 0.5*(data(1)+data(Npts))
    do i=2,Npts-1
      res = res + data(i)
    enddo
    res = res*delta
  end function compute_integral_trapezoid_1d



end module sll_reduction_module