module sll_ascii_io
#include "sll_working_precision.h"

implicit none

contains
  



   subroutine sll_ascii_write_array_1d(file_id,array,error)
   sll_int32 , intent(in)  :: file_id
   sll_int32 , intent(out) :: error
   sll_real64, intent(in)  :: array(:)
                                 
   write(file_id,*,IOSTAT=error) array
                                  
   end subroutine sll_ascii_write_array_1d


   subroutine sll_ascii_write_array_2d(file_id,array,error)
   sll_int32 , intent(in)  :: file_id
   sll_int32 , intent(out) :: error
   sll_real64, intent(in)  :: array(:,:)
                                 
   write(file_id,*,IOSTAT=error) array
                                  
   end subroutine sll_ascii_write_array_2d


   subroutine sll_ascii_write_array_3d(file_id,array,error)
   sll_int32 , intent(in)  :: file_id
   sll_int32 , intent(out) :: error
   sll_real64, intent(in)  :: array(:,:,:)
                                 
   write(file_id,*,IOSTAT=error) array
                                  
   end subroutine sll_ascii_write_array_3d


end module sll_ascii_io
