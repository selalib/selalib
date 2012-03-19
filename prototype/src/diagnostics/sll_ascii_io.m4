module sll_ascii_io
#include "sll_working_precision.h"

implicit none

contains
  
define(`NEW_FUNCTION',
`
   subroutine `$1'(file_id,array,error)
   sll_int32 , intent(in)  :: file_id
   sll_int32 , intent(out) :: error
   sll_real64, intent(in)  :: `$2'
                                 
   write(file_id,*,IOSTAT=error) array
                                  
   end subroutine `$1'
')

NEW_FUNCTION(sll_ascii_write_array_1d,array(:))
NEW_FUNCTION(sll_ascii_write_array_2d,array(:,:))
NEW_FUNCTION(sll_ascii_write_array_3d,array(:,:,:))

end module sll_ascii_io
