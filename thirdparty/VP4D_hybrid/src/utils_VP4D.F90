!===========================================================================
!> Different useful types for 4D Vlasov-Poisson hybrid simulation
!>
!> \date 2015-02-26
!> \author V. Grandgirard, A. Back
!---------------------------------------------------------------------------
module utils_VP4D_module
#include "sll_working_precision.h"

  use sll_m_collective
  use sll_m_remapper

  !===========================================================================
  !> Type used to define boundary conditions in each of 
  !>  the fourth directions
  !---------------------------------------------------------------------------
  type, public :: boundary_conditions_4d_t

    sll_int32 :: left_eta1
    sll_int32 :: right_eta1
    sll_int32 :: left_eta2
    sll_int32 :: right_eta2
    sll_int32 :: left_vx
    sll_int32 :: right_vx
    sll_int32 :: left_vy
    sll_int32 :: right_vy
    
  end type boundary_conditions_4d_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Type used to define boundary conditions in the two
  !>  directions (eta1,eta2)
  !---------------------------------------------------------------------------
  type, public :: boundary_conditions_2d_t

    sll_int32 :: left_eta1
    sll_int32 :: right_eta1
    sll_int32 :: left_eta2
    sll_int32 :: right_eta2
    
  end type boundary_conditions_2d_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Type used to define spline degrees in the four directions
  !---------------------------------------------------------------------------
  type, public :: spline_degree_4d_t

    sll_int32 :: eta1
    sll_int32 :: eta2
    sll_int32 :: vx
    sll_int32 :: vy
    
  end type spline_degree_4d_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Type used to define spline degrees in the two directions
  !>  eta1 and eta2 
  !---------------------------------------------------------------------------
  type, public :: spline_degree_2d_t

    sll_int32 :: eta1
    sll_int32 :: eta2
    
  end type spline_degree_2d_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Return value = 1.0
  !---------------------------------------------------------------------------
  function func_one( eta1, eta2, params ) result(res)

    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res

    res = 1.0_f64

  end function func_one
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Return value = -1.0
  !---------------------------------------------------------------------------
  function func_minus_one( eta1, eta2, params ) result(res)

    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res

    res = -1.0_f64

  end function func_minus_one
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Return value = 0.0
  !---------------------------------------------------------------------------
  function func_zero( eta1, eta2, params ) result(res)

    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res

    res = 0.0_f64

  end function func_zero
  !---------------------------------------------------------------------------


  !===========================================================================
  !> \todo: this function should be a common function in SELALIB
  !---------------------------------------------------------------------------
  subroutine compute_displacements_array( layout, collective_size, disps )

    type(sll_t_layout_2d)               , pointer     :: layout
    sll_int32                     , intent(in)  :: collective_size
    sll_int32, &
        dimension(collective_size), intent(out) :: disps

    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j
    sll_int32 :: counter
    sll_int32 :: rank

    counter  = 0
    disps(1) = counter
    do rank = 1,collective_size-1
       imin    = sll_o_get_layout_i_min( layout, rank-1 )
       imax    = sll_o_get_layout_i_max( layout, rank-1 )
       jmin    = sll_o_get_layout_j_min( layout, rank-1 )
       jmax    = sll_o_get_layout_j_max( layout, rank-1 )
       size_i  = imax - imin + 1
       size_j  = jmax - jmin + 1
       counter = counter + size_i*size_j
       disps(rank+1) = counter
    end do

  end subroutine compute_displacements_array
  !---------------------------------------------------------------------------


  !===========================================================================
  !> \todo: this function should be a common function in SELALIB
  !---------------------------------------------------------------------------
  subroutine load_buffer( layout, data, buffer )

    type(sll_t_layout_2d)           , pointer     :: layout
    sll_real64, dimension(:,:), intent(in)  :: data
    sll_real64, dimension(:)  , intent(out) :: buffer

    sll_int32 :: myrank
    sll_int32 :: data_size
    type(sll_t_collective_t), pointer :: col
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j
    sll_int32 :: i,j
    sll_int32 :: counter

    col => sll_o_get_layout_collective( layout )
    myrank    = sll_f_get_collective_rank( col )
    data_size = size(data,1)*size(data,2)

    imin = sll_o_get_layout_i_min( layout, myrank )
    imax = sll_o_get_layout_i_max( layout, myrank )
    jmin = sll_o_get_layout_j_min( layout, myrank )
    jmax = sll_o_get_layout_j_max( layout, myrank )
    size_i = imax - imin + 1
    size_j = jmax - jmin + 1

    counter = 0
    do j = 1,size_j
      do i = 1,size_i
        counter = counter + 1
        buffer(counter) = data(i,j)
      end do
    end do

  end subroutine load_buffer
  !---------------------------------------------------------------------------


  !===========================================================================
  !> \todo: this function should be a common function in SELALIB
  !---------------------------------------------------------------------------
  function receive_counts_array( layout, n ) result(rc)

    type(sll_t_layout_2d), pointer      :: layout
    sll_int32      , intent(in)   :: n
    sll_int32      , dimension(n) :: rc

    sll_int32 :: i
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j

    do i = 0,n-1
       imin    = sll_o_get_layout_i_min( layout, i )
       imax    = sll_o_get_layout_i_max( layout, i )
       jmin    = sll_o_get_layout_j_min( layout, i )
       jmax    = sll_o_get_layout_j_max( layout, i )
       size_i  = imax - imin + 1
       size_j  = jmax - jmin + 1
       rc(i+1) = size_i*size_j
    end do

  end function receive_counts_array
  !---------------------------------------------------------------------------


  !===========================================================================
  !> \todo: this function should be a common function in SELALIB
  !---------------------------------------------------------------------------
  subroutine unload_buffer( layout, buffer, data )

    type(sll_t_layout_2d), pointer                :: layout
    sll_real64, dimension(:,:), intent(out) :: data
    sll_real64, dimension(:), intent(in)    :: buffer

    sll_int32 :: col_sz
    type(sll_t_collective_t), pointer :: col
    sll_int32 :: i, j
    sll_int32 :: box
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: pos            ! position in buffer

    col => sll_o_get_layout_collective( layout )
    col_sz = sll_f_get_collective_size( col )

    ! loop over all the boxes in the layout and fill the data array by chunks.
    pos = 1
    do box = 0,col_sz-1
       imin = sll_o_get_layout_i_min( layout, box )
       imax = sll_o_get_layout_i_max( layout, box )
       jmin = sll_o_get_layout_j_min( layout, box )
       jmax = sll_o_get_layout_j_max( layout, box )
       ! this will fill the data array in whatever order that the boxes
       ! are ordered.
       do j = jmin,jmax
          do i = imin,imax
             data(i,j) = buffer(pos)
             pos       = pos + 1
          end do
       end do
    end do

  end subroutine unload_buffer
  !---------------------------------------------------------------------------

end module utils_VP4D_module
!---------------------------------------------------------------------------
