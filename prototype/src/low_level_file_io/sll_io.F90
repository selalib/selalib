!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_io
!
!> @author
!> Pierre Navaro
!> Edwin
!
! DESCRIPTION: 
!
!> @brief
!> Implements the functions to display selalib objects
!>
!>@details
!> 
!
! REVISION HISTORY:
! 30 03 2012 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_io
  
  use sll_xdmf
  use sll_gnuplot

#ifdef STDF95
    integer, parameter :: SLL_IO_XDMF = 0, &
                  SLL_IO_VTK  = 1, &
                  SLL_IO_GNUPLOT = 2
#else
  enum, bind(C)
    enumerator :: SLL_IO_XDMF = 0, &
                  SLL_IO_VTK  = 1, &
                  SLL_IO_GNUPLOT = 2
  end enum
#endif
  
end module sll_io
