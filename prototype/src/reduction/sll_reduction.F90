module sll_reduction_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_remapper
  use sll_logical_meshes

contains

  !--------------------------------------------------
  ! Generic function for computing charge density
  ! we should also add a choice for the integration 
  ! also should be generalized for more complicated data
  ! as here array of values f_0,\dots,f_N 
  !---------------------------------------------------  
  
  ! function abstract interface
  ! input,Npts,param optionnels
   
  
  
  ! compute_reduction_4d_to_3d_direction4

  subroutine compute_reduction_4d_to_3d_direction4(&
    data_4d, &
    data_3d, &
    Npts1, &
    Npts2, &
    Npts3, &
    Npts4, &
    delta)
    sll_real64, dimension(:,:,:,:), intent(in)    :: data_4d
    sll_real64, dimension(:,:,:)  , intent(out) :: data_3d
    sll_int32, intent(in)  :: Npts1
    sll_int32, intent(in)  :: Npts2
    sll_int32, intent(in)  :: Npts3
    sll_int32, intent(in)  :: Npts4
    sll_real64, intent(in) :: delta
    sll_int32  :: i1, i2, i3, i4

    
    
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

    do i3 = 1,Npts3
      do i2 = 1,Npts2
        do i1 = 1,Npts1
          tmp = 0.5_f64*(data_4d(i1,i2,i3,1)&
            +data_4d(i1,i2,i3,Npts4))
          do i4 = 2,Npts4-1
            tmp = tmp + data_4d(i1,i2,i3,i4)
          end do
          data_3d(i1,i2,i3) = tmp*delta
        end do
      end do
    end do
  end subroutine compute_reduction_4d_to_3d_direction4
  
!  function compute_integral_trapezoid_1d(data, Npts, delta)
!    sll_real64, dimension(:), intent(in)    :: data
!    sll_int32, intent(in) :: Npts
!    sll_real64,intent(in) :: delta
!    sll_real64 :: tmp
!    sll_int32 :: i
!    
!    tmp = 0.5*(data(1)+data(Npts))
!    do i=2,Npts-1
!      tmp = t
!    enddo
!    
!  end function compute_integral_trapezoid



end module sll_reduction_module