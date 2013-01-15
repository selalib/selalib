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

!> @author Pierre Navaro
!> @details use it for GNUplot
!> Module that contains routines to write data in file in ASCII format
module sll_ascii_io
#include "sll_working_precision.h"
#include "sll_assert.h"
  implicit none
  
  !> Write nD array in ascii file 
  !> \param[in] file_id file unit number
  !> \param[in] array array
  !> \param[out] error error code
  interface sll_ascii_write_array
     module procedure sll_ascii_write_array_1d 
     module procedure sll_ascii_write_array_2d
     module procedure sll_ascii_write_array_3d
  end interface
  
contains

!>Create ASCII file
subroutine sll_ascii_file_create(file_name, file_id, error)

   character(len=*), intent(in)         :: file_name !< file name
   sll_int32, intent(out)               :: error     !< error code
   sll_int32, intent(out)               :: file_id   !< file unit number
   logical                              :: lopen     
   
   error=0

   do 100 file_id=20,99

      inquire(unit=file_id,opened=lopen)
      if(lopen) then
         cycle
      else
         open(file_id,status='SCRATCH',err=100)
         close(file_id,status='DELETE',err=100)
         goto 200
      end if
 
   100 continue
   error=1
   200 continue
   error=0

   open(file_id,FILE=file_name,FORM='FORMATTED',IOSTAT=error)
   rewind(file_id)

end subroutine sll_ascii_file_create

#define NEW_FUNCTION(func_name, array_name_and_dims) \
subroutine func_name(file_id,array,error);           \
sll_int32 , intent(in)  :: file_id;                  \
sll_int32 , intent(out) :: error;                    \
sll_real64, intent(in)  :: array_name_and_dims;      \
write(file_id,*,IOSTAT=error) array;                 \
end subroutine func_name

NEW_FUNCTION(sll_ascii_write_array_1d,array(:))

NEW_FUNCTION(sll_ascii_write_array_2d,array(:,:))

NEW_FUNCTION(sll_ascii_write_array_3d,array(:,:,:))

end module sll_ascii_io
