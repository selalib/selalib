
!***************************************************************************
!
! Selalib 2012     
! Module: sll_poisson_3d_periodic.F90
!
!> @brief 
!> Selalib 3D poisson solver
!> Start date: Feb. 23, 2012
!> Last modification: Feb. 23, 2012
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

  use sll_fft
  use numeric_constants

  implicit none

  type poisson_3d_periodic_plan
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
  end type poisson_3d_periodic_plan

contains

  function new_poisson_3d_periodic_plan(array, Lx, Ly, Lz)

    sll_comp64, dimension(:,:,:)             :: array
    sll_int32                                :: nx, ny, nz
    sll_int32                                :: ierr
    sll_real64                               :: Lx, Ly, Lz
    type (poisson_3d_periodic_plan), pointer :: new_poisson_3d_periodic_plan

    nx = size(array,1)
    ny = size(array,2)
    nz = size(array,3)

    SLL_ALLOCATE(new_poisson_3d_periodic_plan, ierr)

    new_poisson_3d_periodic_plan%nx = nx
    new_poisson_3d_periodic_plan%ny = ny
    new_poisson_3d_periodic_plan%nz = nz
    new_poisson_3d_periodic_plan%Lx = Lx
    new_poisson_3d_periodic_plan%Ly = Ly
    new_poisson_3d_periodic_plan%Lz = Lz

    new_poisson_3d_periodic_plan%px => new_plan_c2c_1d( nx, array(:,1,1), array(:,1,1), FFT_FORWARD )
    new_poisson_3d_periodic_plan%py => new_plan_c2c_1d( ny, array(1,:,1), array(1,:,1), FFT_FORWARD )
    new_poisson_3d_periodic_plan%pz => new_plan_c2c_1d( nz, array(1,1,:), array(1,1,:), FFT_FORWARD )
    new_poisson_3d_periodic_plan%px_inv => new_plan_c2c_1d( nx, array(:,1,1), array(:,1,1), FFT_INVERSE)
    new_poisson_3d_periodic_plan%py_inv => new_plan_c2c_1d( ny, array(1,:,1), array(1,:,1), FFT_INVERSE )
    new_poisson_3d_periodic_plan%pz_inv => new_plan_c2c_1d( nz, array(1,1,:), array(1,1,:), FFT_INVERSE )

  end function new_poisson_3d_periodic_plan

  subroutine delete_poisson_3d_periodic_plan(plan)

    type (poisson_3d_periodic_plan), pointer :: plan
    sll_int32                                :: ierr

    ! Fixme: some error checking, whether the poisson pointer is associated
    ! for instance
    SLL_ASSERT( associated(plan) )

    call fft_delete_plan(plan%px)
    call fft_delete_plan(plan%py)
    call fft_delete_plan(plan%pz)

    call fft_delete_plan(plan%px_inv)
    call fft_delete_plan(plan%py_inv)
    call fft_delete_plan(plan%pz_inv)

    SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_poisson_3d_periodic_plan

end module sll_poisson_3d_periodic_util

