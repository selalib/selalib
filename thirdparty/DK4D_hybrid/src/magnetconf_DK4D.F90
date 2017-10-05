!===========================================================================
!> Initialization of the magnetic configuration for 
!>  4D drift-kinetic hybrid simulation
!>
!> \date 2014-08-20
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module magnetconf_DK4D_module
#include "sll_working_precision.h"
#include "sll_memory.h"

  use mesh_DK4D_module

  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: magnetconf_DK4D_t

    !> Magnetic field
    sll_real64, dimension(:,:), pointer :: B_xy

  end type magnetconf_DK4D_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Magnetic configuration: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_magnetconf_DK4D( magnetconf, mesh4d )

    type(magnetconf_DK4D_t), intent(inout) :: magnetconf
    type(mesh_DK4D_t)      , intent(in)    :: mesh4d

    sll_int32 :: Nx, Ny
    sll_int32 :: ierr

    Nx = size(mesh4d%xgrid_2d,1)
    Ny = size(mesh4d%xgrid_2d,2)
    SLL_ALLOCATE( magnetconf%B_xy(Nx,Ny), ierr )

  end subroutine new_magnetconf_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Magnetic configuration: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_magnetconf_DK4D( magnetconf )

    type(magnetconf_DK4D_t), intent(inout) :: magnetconf

    sll_int32 :: ierr

    SLL_DEALLOCATE( magnetconf%B_xy, ierr )

  end subroutine delete_magnetconf_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Magnetic configuration: Initialization
  !---------------------------------------------------------------------------
  subroutine init_magnetconf_DK4D( magnetconf, mesh4d )

    type(magnetconf_DK4D_t), intent(inout) :: magnetconf
    type(mesh_DK4D_t)      , intent(in)    :: mesh4d

    sll_int32 :: ix, iy
    sll_int32 :: Nx, Ny
 
    !*** Initialization of B(x,y) ***
    Nx = size(mesh4d%xgrid_2d,1)
    Ny = size(mesh4d%xgrid_2d,2)
    do iy = 1,Ny
      do ix = 1,Nx
        magnetconf%B_xy(ix,iy) = 1._f64
      end do
    end do

  end subroutine init_magnetconf_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Magnetic configuration: Printing in HDF5 format
  !---------------------------------------------------------------------------
  subroutine print_magnetconf_DK4D( magnetconf, filename )

    use sll_m_collective
    use sll_m_hdf5_io_serial, only: sll_s_hdf5_ser_file_create, &
        sll_o_hdf5_ser_write_array, sll_s_hdf5_ser_file_close, sll_t_hdf5_ser_handle

    type(magnetconf_DK4D_t), intent(in) :: magnetconf
    character(len=*)       , intent(in) :: filename

    sll_int32 :: my_rank

    !--> For  HDF5 saving
    integer   :: file_err
    type(sll_t_hdf5_ser_handle) :: handle    !< file handle

    my_rank = sll_f_get_collective_rank(sll_v_world_collective)
    if (my_rank.eq.0) then
      call sll_s_hdf5_ser_file_create( trim(filename), handle, file_err )
      call sll_o_hdf5_ser_write_array( handle, magnetconf%B_xy, &
          'B_xy', file_err )
      call sll_s_hdf5_ser_file_close( handle, file_err )
    end if

  end subroutine print_magnetconf_DK4D
  !---------------------------------------------------------------------------


end module magnetconf_DK4D_module
!---------------------------------------------------------------------------
