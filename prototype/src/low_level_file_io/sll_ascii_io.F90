module sll_ascii_io
#include "sll_working_precision.h"
implicit none


interface sll_ascii_write_array
   module procedure sll_ascii_write_array_1d
   module procedure sll_ascii_write_array_2d
   module procedure sll_ascii_write_array_3d
end interface sll_ascii_write_array

contains
  
#define NEW_FUNCTION(func_name, array_name_and_dims) \
subroutine func_name(file_id,array,error);           \
sll_int32 , intent(in)  :: file_id;                  \
sll_int32 , intent(out) :: error;                    \
sll_real64, intent(in)  :: array_name_and_dims;      \
                                                     \
write(file_id,*,IOSTAT=error) array;                 \
                                                     \
end subroutine func_name

NEW_FUNCTION(sll_ascii_write_array_1d,array(:))
NEW_FUNCTION(sll_ascii_write_array_2d,array(:,:))
NEW_FUNCTION(sll_ascii_write_array_3d,array(:,:,:))

end module sll_ascii_io
