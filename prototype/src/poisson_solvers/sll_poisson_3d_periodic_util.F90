
!***************************************************************************
!
! Selalib 2012     
! Module: sll_poisson_3d_periodic.F90
!
!> @brief 
!> Selalib 3D poisson solver
!> Start date: Feb. 23, 2012
!> Last modification: March 20, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************

module sll_poisson_3d_periodic_util

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"
#include "sll_remap.h"

  use sll_fft
  use numeric_constants
  use sll_collective

  implicit none

  type poisson_3d_periodic_plan_seq
     sll_int32                   :: nx ! Number of points-1 in x-direction
     sll_int32                   :: ny ! Number of points-1 in y-direction
     sll_int32                   :: nz ! Number of points-1 in z-direction
     sll_real64                  :: Lx
     sll_real64                  :: Ly
     sll_real64                  :: Lz
     type(sll_fft_plan), pointer :: px
     type(sll_fft_plan), pointer :: py
     type(sll_fft_plan), pointer :: pz
     type(sll_fft_plan), pointer :: px_inv
     type(sll_fft_plan), pointer :: py_inv
     type(sll_fft_plan), pointer :: pz_inv
  end type poisson_3d_periodic_plan_seq

  type poisson_3d_periodic_plan_par
     sll_int32                   :: nx ! Number of points-1 in x-direction
     sll_int32                   :: ny ! Number of points-1 in y-direction
     sll_int32                   :: nz ! Number of points-1 in z-direction
     sll_real64                  :: Lx
     sll_real64                  :: Ly
     sll_real64                  :: Lz
     type(sll_fft_plan), pointer :: px
     type(sll_fft_plan), pointer :: py
     type(sll_fft_plan), pointer :: pz
     type(sll_fft_plan), pointer :: px_inv
     type(sll_fft_plan), pointer :: py_inv
     type(sll_fft_plan), pointer :: pz_inv
     type(layout_3D_t),  pointer :: layout_x
     type(layout_3D_t),  pointer :: layout_y
     type(layout_3D_t),  pointer :: layout_z
  end type poisson_3d_periodic_plan_par

   contains


     function new_poisson_3d_periodic_plan_seq(nx ,ny ,nz, Lx, Ly, Lz) result(plan)

       sll_comp64,                    dimension(nx) :: x
       sll_comp64,                    dimension(ny) :: y
       sll_comp64,                    dimension(nz) :: z
       sll_int32                                    :: nx, ny, nz
       sll_int32                                    :: ierr
       sll_real64                                   :: Lx, Ly, Lz
       type (poisson_3d_periodic_plan_seq), pointer :: plan

       SLL_ALLOCATE(plan, ierr)

       ! Geometry informations
       plan%nx = nx
       plan%ny = ny
       plan%nz = nz
       plan%Lx = Lx
       plan%Ly = Ly
       plan%Lz = Lz

       ! For FFTs (in each direction)
       plan%px => new_plan_c2c_1d( nx, x, x, FFT_FORWARD )
       plan%py => new_plan_c2c_1d( ny, y, y, FFT_FORWARD )
       plan%pz => new_plan_c2c_1d( nz, z, z, FFT_FORWARD )

       ! For inverse FFTs (in each direction)
       plan%px_inv => new_plan_c2c_1d( nx, x, x, FFT_INVERSE)
       plan%py_inv => new_plan_c2c_1d( ny, y, y, FFT_INVERSE )
       plan%pz_inv => new_plan_c2c_1d( nz, z, z, FFT_INVERSE )

     end function new_poisson_3d_periodic_plan_seq


     function new_poisson_3d_periodic_plan_par(start_layout, nx, ny, nz, Lx, Ly, Lz) result(plan)

       type(layout_3D_t),  pointer                  :: start_layout
       sll_int32                                    :: nx, ny, nz
       sll_comp64,                    dimension(nx) :: x
       sll_comp64,                    dimension(ny) :: y
       sll_comp64,                    dimension(nz) :: z
       sll_real64                                   :: Lx, Ly, Lz
       sll_int64                                    :: colsz ! collective size
       sll_int32                                    :: npx, npy, npz
       sll_int32                                    :: e
       sll_int32                                    :: ierr
       type (poisson_3d_periodic_plan_par), pointer :: plan

       SLL_ALLOCATE(plan, ierr)

       ! Geometry informations
       plan%nx = nx
       plan%ny = ny
       plan%nz = nz
       plan%Lx = Lx
       plan%Ly = Ly
       plan%Lz = Lz

       ! For FFTs (in each direction)
       plan%px => new_plan_c2c_1d( nx, x, x, FFT_FORWARD )
       plan%py => new_plan_c2c_1d( ny, y, y, FFT_FORWARD )
       plan%pz => new_plan_c2c_1d( nz, z, z, FFT_FORWARD )

       ! For inverse FFTs (in each direction)
       plan%px_inv => new_plan_c2c_1d( nx, x, x, FFT_INVERSE )
       plan%py_inv => new_plan_c2c_1d( ny, y, y, FFT_INVERSE )
       plan%pz_inv => new_plan_c2c_1d( nz, z, z, FFT_INVERSE )

       ! Layout and local sizes for FFTs in x-direction
       plan%layout_x => start_layout

       ! Layout and local sizes for FFTs in y-direction
       colsz  = sll_get_collective_size(sll_world_collective)
       e = int(log(real(colsz))/log(2.))
       plan%layout_y => new_layout_3D( sll_world_collective )
       npx = 2**(e/2)
       npy = 1
       npz = int(colsz)/npx
       call initialize_layout_with_distributed_3D_array( nx, ny, nz, npx, npy, npz, plan%layout_y )

       ! Layout and local sizes for FFTs in z-direction
       plan%layout_z => new_layout_3D( sll_world_collective )
       npy = npz
       npz = 1
       call initialize_layout_with_distributed_3D_array( nx, ny, nz, npx, npy, npz, plan%layout_z )

     end function new_poisson_3d_periodic_plan_par


     subroutine delete_poisson_3d_periodic_plan_seq(plan)

       type (poisson_3d_periodic_plan_seq), pointer :: plan
       sll_int32                                    :: ierr

       ! Fixme: some error checking, whether the poisson pointer is associated
       ! for instance
       SLL_ASSERT( associated(plan) )

       call delete_fft_plan1d(plan%px)
       call delete_fft_plan1d(plan%py)
       call delete_fft_plan1d(plan%pz)

       call delete_fft_plan1d(plan%px_inv)
       call delete_fft_plan1d(plan%py_inv)
       call delete_fft_plan1d(plan%pz_inv)

       SLL_DEALLOCATE(plan, ierr)

     end subroutine delete_poisson_3d_periodic_plan_seq


     subroutine delete_poisson_3d_periodic_plan_par(plan)

       type (poisson_3d_periodic_plan_par), pointer :: plan
       sll_int32                                    :: ierr

       ! Fixme: some error checking, whether the poisson pointer is associated
       ! for instance
       SLL_ASSERT( associated(plan) )

       call delete_fft_plan1d(plan%px)
       call delete_fft_plan1d(plan%py)
       call delete_fft_plan1d(plan%pz)

       call delete_fft_plan1d(plan%px_inv)
       call delete_fft_plan1d(plan%py_inv)
       call delete_fft_plan1d(plan%pz_inv)

       call delete_layout_3D( plan%layout_x )
       call delete_layout_3D( plan%layout_y )
       call delete_layout_3D( plan%layout_z )

       SLL_DEALLOCATE(plan, ierr)

     end subroutine delete_poisson_3d_periodic_plan_par


end module sll_poisson_3d_periodic_util

