!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
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

!> @file sll_logical_multipatch.F90
!> @namespace sll_logical_meshes_multipatch
!> @brief basic types to describe a collection of logical meshes
!> associated with the decomposition of a physical region with the 
!> multipatch approach.

module sll_logical_meshes_multipatch
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_logical_meshes
  implicit none

  !> @brief basic logical mesh multipatch object.
  type sll_logical_mesh_multipatch_2d
     sll_int32 :: number_patches
     type(sll_logical_mesh_2d_ptr), dimension(:), pointer :: meshes
   contains
     procedure, pass :: initialize_patch => initialize_patch_lmmp2d
     procedure, pass :: get_num_cells1   => get_num_cells1_lmmp2d
     procedure, pass :: get_num_cells2   => get_num_cells2_lmmp2d
     procedure, pass :: get_eta1_min     => get_eta1_min_lmmp2d
     procedure, pass :: get_eta1_max     => get_eta1_max_lmmp2d
     procedure, pass :: get_eta2_min     => get_eta2_min_lmmp2d
     procedure, pass :: get_eta2_max     => get_eta2_max_lmmp2d
     procedure, pass :: get_delta_eta1   => get_delta_eta1_lmmp2d
     procedure, pass :: get_delta_eta2   => get_delta_eta2_lmmp2d
     procedure, pass :: get_logical_mesh => get_logical_mesh_lmmp2d
  end type sll_logical_mesh_multipatch_2d

  interface sll_delete
     module procedure delete_logical_mesh_lmmp2d
  end interface sll_delete


contains

  !> @brief allocates the memory space for a new 2D multipatch logical mesh 
  !> in the heap and returns a pointer to the newly allocated object. 
  !> This object must be then initialized by the appropriate call to the 
  !> initialize_patch method.
  !> @param num_cells1 integer denoting the number of cells.
  !> @param eta1_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh.
  !> @param eta1_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh.
  !> return a pointer to the newly allocated object.
  function new_logical_mesh_multipatch_2d( num_patches ) result(m)
    type(sll_logical_mesh_multipatch_2d), pointer :: m
    sll_int32, intent(in) :: num_patches
    sll_int32 :: ierr
    SLL_ALLOCATE(m,ierr)
    SLL_ALLOCATE(m%meshes(num_patches),ierr)
    m%number_patches = num_patches
  end function new_logical_mesh_multipatch_2d


  !> @brief initializes a previously allocated multipatch 2d logical mesh.
  !> @param mp the multipatch object, passed implicitly.
  !> @param num_cells1 integer denoting the number of cells in first direction.
  !> @param num_cells2 integer denoting the number of cells in second direction.
  !> @param eta1_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh.
  !> @param eta1_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh.
  !> @param eta2_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh.
  !> @param eta2_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh.
  subroutine initialize_patch_lmmp2d( &
    mp, &
    patch, &
    num_cells1, &
    num_cells2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max )

    class(sll_logical_mesh_multipatch_2d), intent(inout) :: mp
    sll_int32, intent(in)  :: patch
    sll_int32, intent(in)  :: num_cells1
    sll_int32, intent(in)  :: num_cells2
    sll_real64, optional, intent(in) :: eta1_min
    sll_real64, optional, intent(in) :: eta1_max
    sll_real64, optional, intent(in) :: eta2_min
    sll_real64, optional, intent(in) :: eta2_max

    if( (patch < 0) .or. (patch > mp%number_patches - 1) ) then
       print *, 'ERROR, initialize_patch_2d(): the patch argument has to be ',&
            'between 0 and the number of patches - 1'
       stop
    end if
    print *, 'initializing patch: ', patch, size(mp%meshes)
    mp%meshes(patch+1)%lm => new_logical_mesh_2d( &
                                num_cells1, &
                                num_cells2, &
                                eta1_min, &
                                eta1_max, &
                                eta2_min, &
                                eta2_max )
  end subroutine initialize_patch_lmmp2d

#define MAKE_ACCESS_FUNC_MULTIPATCH( fname, obj_typ, slot, ret_typ ) \
  function fname( mp, patch ) result(res); \
    ret_typ :: res; \
    class(obj_typ), intent(in) :: mp; \
    sll_int32, intent(in) :: patch; \
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%number_patches) ); \
    res = mp%meshes(patch+1)%lm%slot; \
  end function fname

#define OBJECT sll_logical_mesh_multipatch_2d

MAKE_ACCESS_FUNC_MULTIPATCH(get_num_cells1_lmmp2d,OBJECT,num_cells1,sll_int32)
MAKE_ACCESS_FUNC_MULTIPATCH(get_num_cells2_lmmp2d,OBJECT,num_cells2,sll_int32)
MAKE_ACCESS_FUNC_MULTIPATCH(get_eta1_min_lmmp2d, OBJECT, eta1_min,sll_real64)
MAKE_ACCESS_FUNC_MULTIPATCH(get_eta1_max_lmmp2d, OBJECT, eta1_max,sll_real64)
MAKE_ACCESS_FUNC_MULTIPATCH(get_eta2_min_lmmp2d, OBJECT, eta2_min,sll_real64)
MAKE_ACCESS_FUNC_MULTIPATCH(get_eta2_max_lmmp2d, OBJECT, eta2_max,sll_real64)
MAKE_ACCESS_FUNC_MULTIPATCH(get_delta_eta1_lmmp2d,OBJECT,delta_eta1,sll_real64)
MAKE_ACCESS_FUNC_MULTIPATCH(get_delta_eta2_lmmp2d,OBJECT,delta_eta2,sll_real64)

#undef MAKE_ACCESS_FUNC_MULTIPATCH
#undef LMMP2D

  function get_logical_mesh_lmmp2d( mp, patch ) result(res)
    type(sll_logical_mesh_2d), pointer               :: res
    class(sll_logical_mesh_multipatch_2d), intent(in) :: mp
    sll_int32, intent(in)                            :: patch
    SLL_ASSERT( (patch >=0) .and. (patch < mp%number_patches) )
    res => mp%meshes(patch+1)%lm
  end function get_logical_mesh_lmmp2d

  subroutine delete_logical_mesh_lmmp2d( mp )
    type(sll_logical_mesh_multipatch_2d), pointer :: mp
    sll_int32 :: ierr
    SLL_DEALLOCATE(mp%meshes, ierr)
    SLL_DEALLOCATE(mp, ierr)
  end subroutine delete_logical_mesh_lmmp2d

end module sll_logical_meshes_multipatch
