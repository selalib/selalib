!------------------------------------------------------------------------------
! Selalib
!------------------------------------------------------------------------------
! MODULE: sll_ascii_io
!
! DESCRIPTION:
!> @file sll_ascii_io.F90
!> @namespace sll_ascii_io
!> @author Pierre Navaro IRMA http://www-irma.u-strasbg.fr
!> @brief Module that contains routine to write data in file in ASCII 
!! format
!> @details 
!> Don't forget to open a file before using these routines
!>
!> \remark use it for GNUplot
!!         
!>
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------

module sll_ascii_io
#include "sll_working_precision.h"
  implicit none
  
  
  !> @brief interface sll_ascii_write_array
  interface sll_ascii_write_array
     module procedure sll_ascii_write_array_1d
     module procedure sll_ascii_write_array_2d
     module procedure sll_ascii_write_array_3d
  end interface
  
contains
  
#define NEW_FUNCTION(func_name, array_name_and_dims) \
  subroutine func_name(file_id,array,error);         \
    sll_int32 , intent(in)  :: file_id;              \
    sll_int32 , intent(out) :: error;                \
    sll_real64, intent(in)  :: array_name_and_dims;  \
    write(file_id,*,IOSTAT=error) array;             \
  end subroutine func_name

!>Write a 1d array in ASCII file
!! @param[in] file_d file unit
!! @param[in] array
  NEW_FUNCTION(sll_ascii_write_array_1d,array(:))

!>Write a 2d array in ASCII file
!! @param[in] file_d file unit
!! @param[in] array
  NEW_FUNCTION(sll_ascii_write_array_2d,array(:,:))

!>Write a 3d array in ASCII file
!! @param[in] file_d file unit
!! @param[in] array
  NEW_FUNCTION(sll_ascii_write_array_3d,array(:,:,:))

end module sll_ascii_io
