!===========================================================================
!> Initialization of the equilibrium for 4D Vlasov-Poisson hybrid
!>  simulation
!>
!> \date 2015-02-26
!> \author V. Grandgirard, A. Back
!---------------------------------------------------------------------------
module equilibrium_VP4D_module
#include "sll_working_precision.h"
#include "sll_memory.h"

  use input_VP4D_module
  use mesh_VP4D_module
  use sll_m_boundary_condition_descriptors

  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: equilibrium_VP4D_t

    !> Equilibrium distribution function
    sll_real64, dimension(:,:), pointer :: feq_vxvy

  end type equilibrium_VP4D_t
  !---------------------------------------------------------------------------

  contains


  !===========================================================================
  !> Equilibrium: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_equilibrium_VP4D( equilibrium, mesh4d ) 

    type(equilibrium_VP4D_t), intent(inout) :: equilibrium
    type(mesh_VP4D_t)       , intent(in)    :: mesh4d

    sll_int32 :: Nvx, Nvy
    sll_int32 :: ierr

    !*** Array allocations ***
    !-> For equilibrium distribution function
    Nvx = size(mesh4d%vx_grid)
    Nvy = size(mesh4d%vy_grid)
    SLL_ALLOCATE( equilibrium%feq_vxvy(Nvx,Nvy), ierr )

  end subroutine new_equilibrium_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Equilibrium: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_equilibrium_VP4D( equilibrium )

    type(equilibrium_VP4D_t), intent(inout) :: equilibrium
    
    sll_int32 :: ierr

    SLL_DEALLOCATE( equilibrium%feq_vxvy, ierr )

  end subroutine delete_equilibrium_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Equilibrium: Initialization
  !---------------------------------------------------------------------------
  subroutine init_equilibrium_VP4D( equilibrium, mesh4d )

    type(equilibrium_VP4D_t), intent(inout) :: equilibrium
    type(mesh_VP4D_t)       , intent(in)    :: mesh4d

    !*** Initialization of the equilibrium distribution function ***
    call compute_fequilibrium_xy( &
        mesh4d%vx_grid, &
        mesh4d%vy_grid, &
        equilibrium%feq_vxvy )

  end subroutine init_equilibrium_VP4D


  !===========================================================================
  !> Equilibrium: Print in HDF5
  !---------------------------------------------------------------------------
  subroutine print_equilibrium_VP4D( equilibrium, filename )

    use sll_m_collective
    use sll_m_hdf5_io_serial, only: sll_s_hdf5_ser_file_create, &
         sll_o_hdf5_ser_write_array, &
         sll_s_hdf5_ser_file_close, &
         sll_t_hdf5_ser_handle

    type(equilibrium_VP4D_t), intent(in) :: equilibrium
    character(len=*)        , intent(in) :: filename

    sll_int32 :: my_rank

    !--> For  HDF5 saving
    integer   :: file_err
    type(sll_t_hdf5_ser_handle) :: handle    !< file handle

    my_rank = sll_f_get_collective_rank(sll_v_world_collective)
    if (my_rank.eq.0) then
      call sll_s_hdf5_ser_file_create( trim(filename), handle, file_err )
      call sll_o_hdf5_ser_write_array( handle, equilibrium%feq_vxvy, &
          'feq_vxvy', file_err )
      call sll_s_hdf5_ser_file_close( handle, file_err )
    end if

  end subroutine print_equilibrium_VP4D


  !===========================================================================
  ! Initialization of the 2D array for the equilibrium
  !  distribution function feq(vx,vy)
  !  feq(vx,vy) = 1/(2*pi)*exp(-0.5*(vx+vy)**2)
  !---------------------------------------------------------------------------
  subroutine compute_fequilibrium_xy( vx_grid, vy_grid, feq_vxvy )

    use sll_m_constants, only: sll_p_pi
    
    sll_real64, dimension(:)  , intent(in)  :: vx_grid
    sll_real64, dimension(:)  , intent(in)  :: vy_grid
    sll_real64, dimension(:,:), intent(out) :: feq_vxvy

    sll_int32  :: Nvx, Nvy
    sll_int32  :: ivx, ivy
    sll_real64 :: factor
    sll_real64 :: vx_tmp, vy_tmp 

    Nvx = size(vx_grid,1)
    Nvy = size(vy_grid,1)

    factor = 1._f64/(2._f64*sll_p_pi)
    do ivy = 1,Nvy
      vy_tmp = vy_grid(ivy)
      do ivx = 1,Nvx
        vx_tmp = vx_grid(ivx)
        feq_vxvy(ivx,ivy) =  factor * &
            exp(-0.5_f64*(vx_tmp**2+vy_tmp**2))
      end do
    end do

  end subroutine compute_fequilibrium_xy
  !---------------------------------------------------------------------------

end module equilibrium_VP4D_module
!---------------------------------------------------------------------------
