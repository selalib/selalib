module pic_1d_distribution

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_file_io.h"

  use sll_visu_pic, only: compute_df_cic

  implicit none

!==============================================================================

  type, public :: pic1d_eulerian_distribution

    sll_float64, allocatable :: f(:,:)
!    sll_float64, allocatable :: x(:), v(:)
    sll_int32                :: nx, nv
    sll_float64              :: xmin, xmax
    sll_float64              :: vmin, vmax

  contains

    procedure :: initialize      => pic1d_ed__initialize
    procedure :: compute_f       => pic1d_ed__compute_f
    procedure :: compute_moments => pic1d_ed__compute_moments
    procedure :: print_f_to_file => pic1d_ed__print_f_to_file

  end type

!==============================================================================
contains
!==============================================================================

  subroutine pic1d_ed__initialize( self, xmin, xmax, nx, vmin, vmax, nv )
    class( pic1d_eulerian_distribution ), intent( inout ) :: self
    sll_float64                         , intent( in    ) :: xmin
    sll_float64                         , intent( in    ) :: xmax
    sll_int32                           , intent( in    ) :: nx
    sll_float64                         , intent( in    ) :: vmin
    sll_float64                         , intent( in    ) :: vmax
    sll_int32                           , intent( in    ) :: nv

    sll_int32 :: ierr

    self%xmin = xmin
    self%xmax = xmax
    self%nx   = nx
    self%vmin = vmin
    self%vmax = vmax
    self%nv   = nv

    SLL_ALLOCATE( self%f(nx,nv), ierr )

  end subroutine pic1d_ed__initialize

  !----------------------------------------------------------------------------
  subroutine pic1d_ed__compute_f( self, xp, vp, wp )
    class( pic1d_eulerian_distribution ), intent( inout ) :: self
    sll_float64                         , intent( in    ) :: xp(:)
    sll_float64                         , intent( in    ) :: vp(:)
    sll_float64                         , intent( in    ) :: wp(:)
    
    ! Use cloud-in-cell (CIC) algorithm to compute 1D-1V distribution function
    call compute_df_cic( xp, vp, wp, &
      self%xmin, self%xmax, self%nx, &
      self%vmin, self%vmax, self%nv, &
      self%f )

  end subroutine pic1d_ed__compute_f

  !----------------------------------------------------------------------------
  subroutine pic1d_ed__compute_moments( self )
    class( pic1d_eulerian_distribution ), intent( inout ) :: self

    character( len=256 ), parameter :: this_sub_name = &
      'pic1d_ed__compute_moments'

    SLL_WARNING( this_sub_name, 'Not implemented yet.' )

  end subroutine pic1d_ed__compute_moments
  
  !----------------------------------------------------------------------------
  subroutine pic1d_ed__print_f_to_file( self, plot_name, iplot )
    class( pic1d_eulerian_distribution ), intent( inout ) :: self
    character( len=* )                  , intent( in    ) :: plot_name
    sll_int32                           , intent( in    ) :: iplot

    sll_float64      :: dx, dv
    character(len=4) :: fin

    dx = (self%xmax-self%xmin)/(self%nx-1)
    dv = (self%vmax-self%vmin)/(self%nv-1)
    call int2string( iplot, fin )
    call sll_xdmf_corect2d_nodes( plot_name//'_'//fin, self%f, "f(x,v)", &
      self%xmin, dx, self%vmin, dv )

  end subroutine pic1d_ed__print_f_to_file


end module pic_1d_distribution

