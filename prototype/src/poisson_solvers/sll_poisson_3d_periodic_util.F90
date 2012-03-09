
!***************************************************************************
!
! Selalib 2012     
! Module: sll_poisson_3d_periodic.F90
!
!> @brief 
!> Selalib 3D poisson solver
!> Start date: Feb. 23, 2012
!> Last modification: March 09, 2012
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
     type(layout_3D_t),                  pointer :: layout_x
     type(layout_3D_t),                  pointer :: layout_y
     type(layout_3D_t),                  pointer :: layout_z
     type(layout_3D_t),                  pointer :: layout_kernel
     sll_int32,                   dimension(4,3) :: loc_sizes ! local sizes in the 4 layouts
     end type poisson_3d_periodic_plan_par

   contains


     function new_poisson_3d_periodic_plan_seq(array, Lx, Ly, Lz)

       sll_comp64, dimension(:,:,:)                 :: array
       sll_int32                                    :: nx, ny, nz
       sll_int32                                    :: ierr
       sll_real64                                   :: Lx, Ly, Lz
       type (poisson_3d_periodic_plan_seq), pointer :: new_poisson_3d_periodic_plan_seq

       nx = size(array,1)
       ny = size(array,2)
       nz = size(array,3)

       SLL_ALLOCATE(new_poisson_3d_periodic_plan_seq, ierr)

       ! Geometry informations
       new_poisson_3d_periodic_plan_seq%nx = nx
       new_poisson_3d_periodic_plan_seq%ny = ny
       new_poisson_3d_periodic_plan_seq%nz = nz
       new_poisson_3d_periodic_plan_seq%Lx = Lx
       new_poisson_3d_periodic_plan_seq%Ly = Ly
       new_poisson_3d_periodic_plan_seq%Lz = Lz

       ! For FFTs (in each direction)
       new_poisson_3d_periodic_plan_seq%px => new_plan_c2c_1d( nx, array(:,1,1), array(:,1,1), FFT_FORWARD )
       new_poisson_3d_periodic_plan_seq%py => new_plan_c2c_1d( ny, array(1,:,1), array(1,:,1), FFT_FORWARD )
       new_poisson_3d_periodic_plan_seq%pz => new_plan_c2c_1d( nz, array(1,1,:), array(1,1,:), FFT_FORWARD )

       ! For inverse FFTs (in each direction)
       new_poisson_3d_periodic_plan_seq%px_inv => new_plan_c2c_1d( nx, array(:,1,1), array(:,1,1), FFT_INVERSE)
       new_poisson_3d_periodic_plan_seq%py_inv => new_plan_c2c_1d( ny, array(1,:,1), array(1,:,1), FFT_INVERSE )
       new_poisson_3d_periodic_plan_seq%pz_inv => new_plan_c2c_1d( nz, array(1,1,:), array(1,1,:), FFT_INVERSE )

     end function new_poisson_3d_periodic_plan_seq


     function new_poisson_3d_periodic_plan_par(array, Lx, Ly, Lz)

       sll_comp64, dimension(:,:,:)                 :: array
       sll_real64                                   :: Lx, Ly, Lz
       type (poisson_3d_periodic_plan_par), pointer :: new_poisson_3d_periodic_plan_par
       sll_int64                                    :: colsz ! collective size
       sll_int32                                    :: nx, ny, nz
       sll_int32                                    :: npx, npy, npz
       sll_int32                                    :: e, ex, ey, ez
       sll_int32,                    dimension(4,3) :: loc_sizes ! local sizes in the 4 layouts
       sll_int32                                    :: ierr

       SLL_ALLOCATE(new_poisson_3d_periodic_plan_par, ierr)

       nx = size(array,1)
       ny = size(array,2)
       nz = size(array,3)

       ! Geometry informations
       new_poisson_3d_periodic_plan_par%nx = nx
       new_poisson_3d_periodic_plan_par%ny = ny
       new_poisson_3d_periodic_plan_par%nz = nz
       new_poisson_3d_periodic_plan_par%Lx = Lx
       new_poisson_3d_periodic_plan_par%Ly = Ly
       new_poisson_3d_periodic_plan_par%Lz = Lz

       ! For FFTs (in each direction)
       new_poisson_3d_periodic_plan_par%px => new_plan_c2c_1d( nx, array(:,1,1), array(:,1,1), FFT_FORWARD )
       new_poisson_3d_periodic_plan_par%py => new_plan_c2c_1d( ny, array(1,:,1), array(1,:,1), FFT_FORWARD )
       new_poisson_3d_periodic_plan_par%pz => new_plan_c2c_1d( nz, array(1,1,:), array(1,1,:), FFT_FORWARD )

       ! For inverse FFTs (in each direction)
       new_poisson_3d_periodic_plan_par%px_inv => new_plan_c2c_1d( nx, array(:,1,1), array(:,1,1), FFT_INVERSE)
       new_poisson_3d_periodic_plan_par%py_inv => new_plan_c2c_1d( ny, array(1,:,1), array(1,:,1), FFT_INVERSE )
       new_poisson_3d_periodic_plan_par%pz_inv => new_plan_c2c_1d( nz, array(1,1,:), array(1,1,:), FFT_INVERSE )

       colsz  = sll_get_collective_size(sll_world_collective)
       e = int(log(real(colsz))/log(2.))

       ! Layout and local sizes for FFTs in x-direction

       new_poisson_3d_periodic_plan_par%layout_x => new_layout_3D( sll_world_collective )
       npx = 1
       npy = 2**(e/2)
       npz = int(colsz)/npy
       call initialize_layout_with_distributed_3D_array( nx, ny, nz, npx, npy, npz, &
            new_poisson_3d_periodic_plan_par%layout_x )

       new_poisson_3d_periodic_plan_par%loc_sizes(1,1) = nx/npx
       new_poisson_3d_periodic_plan_par%loc_sizes(1,2) = ny/npy
       new_poisson_3d_periodic_plan_par%loc_sizes(1,3) = nz/npz

       ! Layout and local sizes for FFTs in y-direction

       new_poisson_3d_periodic_plan_par%layout_y => new_layout_3D( sll_world_collective )
       npx = npy
       npy = 1
       call initialize_layout_with_distributed_3D_array( nx, ny, nz, npx, npy, npz, &
            new_poisson_3d_periodic_plan_par%layout_y )

       new_poisson_3d_periodic_plan_par%loc_sizes(2,1) = nx/npx
       new_poisson_3d_periodic_plan_par%loc_sizes(2,2) = ny/npy
       new_poisson_3d_periodic_plan_par%loc_sizes(2,3) = nz/npz

       ! Layout and local sizes for FFTs in z-direction

       new_poisson_3d_periodic_plan_par%layout_z => new_layout_3D( sll_world_collective )
       npy = npz
       npz = 1
       call initialize_layout_with_distributed_3D_array( nx, ny, nz, npx, npy, npz, &
            new_poisson_3d_periodic_plan_par%layout_z )

       new_poisson_3d_periodic_plan_par%loc_sizes(3,1) = nx/npx
       new_poisson_3d_periodic_plan_par%loc_sizes(3,2) = ny/npy
       new_poisson_3d_periodic_plan_par%loc_sizes(3,3) = nz/npz

       ! Layout and local sizes for poisson solver kernel
       ex = e/3
       ey = (e-ex)/2
       ez = e - (ex+ey)
       npx = 2**ex
       npy = 2**ey
       npz = 2**ez
       new_poisson_3d_periodic_plan_par%layout_kernel  => new_layout_3D( sll_world_collective ) 
       call initialize_layout_with_distributed_3D_array( nx, ny, nz, npx, npy, npz, &
            new_poisson_3d_periodic_plan_par%layout_kernel )

       new_poisson_3d_periodic_plan_par%loc_sizes(4,1) = nx/npx
       new_poisson_3d_periodic_plan_par%loc_sizes(4,2) = ny/npy
       new_poisson_3d_periodic_plan_par%loc_sizes(4,3) = nz/npz

     end function new_poisson_3d_periodic_plan_par


     subroutine delete_poisson_3d_periodic_plan_seq(plan)

       type (poisson_3d_periodic_plan_seq), pointer :: plan
       sll_int32                                    :: ierr

       ! Fixme: some error checking, whether the poisson pointer is associated
       ! for instance
       SLL_ASSERT( associated(plan) )

       call delete(plan%px)
       call delete(plan%py)
       call delete(plan%pz)

       call delete(plan%px_inv)
       call delete(plan%py_inv)
       call delete(plan%pz_inv)

       SLL_DEALLOCATE(plan, ierr)

     end subroutine delete_poisson_3d_periodic_plan_seq


     subroutine delete_poisson_3d_periodic_plan_par(plan)

       type (poisson_3d_periodic_plan_par), pointer :: plan
       sll_int32                                    :: ierr

       ! Fixme: some error checking, whether the poisson pointer is associated
       ! for instance
       SLL_ASSERT( associated(plan) )

       call delete(plan%px)
       call delete(plan%py)
       call delete(plan%pz)

       call delete(plan%px_inv)
       call delete(plan%py_inv)
       call delete(plan%pz_inv)

       call delete_layout_3D( plan%layout_x )
       call delete_layout_3D( plan%layout_y )
       call delete_layout_3D( plan%layout_z )
       call delete_layout_3D( plan%layout_kernel )

       SLL_DEALLOCATE(plan, ierr)

     end subroutine delete_poisson_3d_periodic_plan_par

     end module sll_poisson_3d_periodic_util

