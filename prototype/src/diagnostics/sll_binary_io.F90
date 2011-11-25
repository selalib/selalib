module sll_binary_io
#include "sll_working_precision.h"
#include "sll_assert.h"

implicit none

contains
  
subroutine sll_binary_file_create(filename,file_id,error)
character(len=*) , intent(in)  :: filename  
sll_int32        , intent(out) :: file_id   
sll_int32        , intent(out) :: error
logical                        :: lopen

error=0

do 100 file_id=20,99

   inquire(unit=file_id,opened=lopen)
   if(lopen) then
      cycle
   else
      open(file_id ,status='SCRATCH',err=100)
      close(file_id,status='DELETE' ,err=100)
      goto 200
   end if
 
100 continue
error=1
200 continue
error=0

! Create a new file using default properties
inquire(file=filename,opened=lopen)
SLL_ASSERT(.not. lopen)

open(file_id,FILE=filename,ACCESS="STREAM",FORM='UNFORMATTED',IOSTAT=error)

end subroutine sll_binary_file_create

subroutine sll_binary_file_close(file_id,error)
sll_int32, intent(in)  :: file_id
sll_int32, intent(out) :: error

close(file_id, IOSTAT=error)

end subroutine sll_binary_file_close

#define NEW_FUNCTION(func_name, dataspace_dimension, array_name_and_dims) \
subroutine func_name(file_id,array,error);                                \
sll_int32 , intent(in)       :: file_id;                                  \
sll_int32 , intent(out)      :: error;                                    \
sll_real64, intent(in)       :: array_name_and_dims;                      \
                                                                          \
write(file_id,IOSTAT=error) array;                                        \
                                                                          \
end subroutine func_name

NEW_FUNCTION(sll_binary_write_array_1d, 1, array(:))
NEW_FUNCTION(sll_binary_write_array_2d, 2, array(:,:))
NEW_FUNCTION(sll_binary_write_array_3d, 3, array(:,:,:))

end module sll_binary_io
