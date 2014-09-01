!**************************************************************
!  Copyright INRIA
!  Authors : 
!     Pierre Navaro 
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @brief
!> Implements the functions to write binary file to store heavy data
!> @details
!> This functions can be used only for sequential application.
!> one file = one dataset \n
!> If HDF5 is not installed you can use this module. \n
!> This is controlled by the variable <code>HDF5_ENABLE</code>.
!>
module sll_binary_io
#include "sll_working_precision.h"
#include "sll_assert.h"
  
  implicit none
  
  !> Write a nD array in a binary file
  !> \param[in] file_id file unit number
  !> \param[in] array array
  !> \param[out] error error code
  interface sll_binary_write_array
     module procedure sll_binary_write_array_1d
     module procedure sll_binary_write_array_2d
     module procedure sll_binary_write_array_3d
  end interface

contains
  
  !> Create binary file
  subroutine sll_binary_file_create(filename,file_id,error)
    character(len=*) , intent(in)  :: filename   !<file name
    sll_int32        , intent(out) :: file_id    !<file unit number
    sll_int32        , intent(out) :: error      !<error code
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
       
100    continue
       error=1
200    continue
       error=0

       ! Create a new file using default properties
       inquire(file=filename,opened=lopen)
       SLL_ASSERT(.not. lopen)

       open(file_id, &
            FILE=filename, &
            ACCESS="STREAM", &
            FORM='UNFORMATTED', &
            IOSTAT=error)

     end subroutine sll_binary_file_create

!> Open binary file
!subroutine sll_binary_file_open(file_id,error)
!sll_int32, intent(in)  :: file_id !<file unit number
!sll_int32, intent(out) :: error   !<error code

!open(file_id, IOSTAT=error)
     
!end subroutine sll_binary_file_open


!> Close binary file
subroutine sll_binary_file_close(file_id,error)
sll_int32, intent(in)  :: file_id !<file unit number
sll_int32, intent(out) :: error   !<error code

close(file_id, IOSTAT=error)
     
end subroutine sll_binary_file_close

!> Write a 0D array in the binary file file_id
subroutine sll_binary_write_array_0d(file_id,array,error)
sll_int32 , intent(in)       :: file_id  !< file unit number
sll_int32 , intent(out)      :: error    !< error code
sll_real64, intent(in)       :: array !< data array
write(file_id,IOSTAT=error) array
end subroutine


!> Write a 1D array in the binary file file_id
subroutine sll_binary_write_array_1d(file_id,array,error)
sll_int32 , intent(in)       :: file_id  !< file unit number
sll_int32 , intent(out)      :: error    !< error code
sll_real64, intent(in)       :: array(:) !< data array
write(file_id,IOSTAT=error) array
end subroutine

!> Write a 2D array in the binary file file_id
subroutine sll_binary_write_array_2d(file_id,array,error)
sll_int32 , intent(in)       :: file_id    !< file unit number
sll_int32 , intent(out)      :: error      !< error code
sll_real64, intent(in)       :: array(:,:) !< adata array
write(file_id,IOSTAT=error) array
end subroutine

!> Write a 3D array in the binary file file_id
subroutine sll_binary_write_array_3d(file_id,array,error)
sll_int32 , intent(in)       :: file_id      !< file unit number
sll_int32 , intent(out)      :: error        !< error code
sll_real64, intent(in)       :: array(:,:,:) !< data array
write(file_id,IOSTAT=error) array
end subroutine


!> Read a 0D array in the binary file file_id
subroutine sll_binary_read_array_0d(file_id,array,error)
sll_int32 , intent(in)       :: file_id  !< file unit number
sll_int32 , intent(out)      :: error    !< error code
sll_real64, intent(out)       :: array !< data array
read(file_id,IOSTAT=error) array
end subroutine

!> Read a 2D array in the binary file file_id
subroutine sll_binary_read_array_2d(file_id,array,error)
sll_int32 , intent(in)       :: file_id    !< file unit number
sll_int32 , intent(out)      :: error      !< error code
sll_real64, intent(out)       :: array(:,:) !< adata array
read(file_id,IOSTAT=error) array
end subroutine



end module sll_binary_io
