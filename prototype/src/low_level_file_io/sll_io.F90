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
!> @brief
!> Implements the functions to display selalib objects
!> @details
module sll_io
  
#include "sll_working_precision.h"
#include "sll_assert.h"

integer, parameter :: SLL_IO_XDMF = 0     !< xdmf outout
integer, parameter :: SLL_IO_VTK  = 1     !< vtk output
integer, parameter :: SLL_IO_GNUPLOT = 2  !< gnuplot output
                      
contains

!> Get a file unit number free before creating a file
subroutine sll_new_file_id(file_id, error)

   sll_int32, intent(out) :: error     !< error code
   sll_int32, intent(out) :: file_id   !< file unit number
   logical                :: lopen    
   
   error=1

   do 100 file_id=20,99

      inquire(unit=file_id,opened=lopen)
      if(lopen) then
         cycle
      else
         open(file_id,status='SCRATCH',err=100)
         close(file_id,status='DELETE',err=100)
         error=0
         exit
      end if
 
   100 continue

   SLL_ASSERT(error == 0)
   

end subroutine sll_new_file_id

  
end module sll_io
