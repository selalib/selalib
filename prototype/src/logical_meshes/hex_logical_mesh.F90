!**************************************************************
!  Copyright INRIA
!  Authors : 
!     Laura Mendoza
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


module hex_logical_meshes
#include "sll_working_precision.h"
#include "sll_memory.h"
  implicit none

  type hex_logical_mesh_2d
     sll_int32  :: num_cells  !< number of cells in three directions
     sll_real64 :: radius     !< maximum value of eta, direction 2
     sll_real64 :: center_x1  ! x coordinate of the origin
     sll_real64 :: center_x2  ! x coordinate of the origin
     sll_real64 :: delta_eta  !< cell spacing, direction 1
  end type hex_logical_mesh_2d

  type hex_logical_mesh_2d_ptr
     type(hex_logical_mesh_2d), pointer :: hm
  end type hex_logical_mesh_2d_ptr

  ! this should be sll_delete library-wide...
  interface delete
     module procedure delete_hex_logical_mesh_2d
  end interface delete

  interface sll_display
     module procedure display_logical_mesh_2d
  end interface sll_display

contains

#define TEST_PRESENCE_AND_ASSIGN_VAL( obj, arg, slot, default_val ) \
  if( present(arg) ) then ; \
    obj%slot = arg; \
  else; \
    obj%slot = default_val; \
end if

  function new_hex_logical_mesh_2d( &
    num_cells, &
    centerx1, &
    centerx2, &
    radius ) result(m)

    type(hex_logical_mesh_2d), pointer :: m
    sll_int32, intent(in)  :: num_cells
    sll_real64, intent(in) :: radius
    sll_real64, optional, intent(in) :: center_x1
    sll_real64, optional, intent(in) :: center_x2
    sll_int32 :: ierr

    SLL_ALLOCATE(m, ierr)
    call initialize_hex_logical_mesh_2d( &
         m, &
         num_cells, &
         radius, &
         center_x1, &
         center_x2)

  end function new_hex_logical_mesh_2d


  subroutine initialize_hex_logical_mesh_2d( &
    m, & 
    num_cells, &
    radius, &
    center_x1, &
    center_x2)

    type(hex_logical_mesh_2d), pointer :: m
    sll_int32, intent(in)  :: num_cells
    sll_real64, intent(in) :: radius
    sll_real64, optional, intent(in) :: center_x1
    sll_real64, optional, intent(in) :: center_x2

    m%num_cells = num_cells
    m%num_radius = radius

    TEST_PRESENCE_AND_ASSIGN_VAL( m, center_x1, center_x1, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, center_x2, center_x2, 0.0_f64 )

    m%delta_eta   = (2. * m%radius)/real(num_cells1,f64)

    if ( m%radius <= 0.) then
       print*,'ERROR, initialize_hex_logical_mesh_2d(): ', &
            'Problem to construct the mesh 2d '
       print*,'because radius <= 0.'
    end if
  end subroutine initialize_hex_logical_mesh_2d


  subroutine delete_hex_logical_mesh_2d( mesh )
    type(hex_logical_mesh_2d), pointer :: mesh
    sll_int32 :: ierr
    if(.not. associated(mesh))then
       print *, 'delete_hex_logical_mesh_2d, ERROR: passed argument is not ', &
            'associated. Crash imminent...'
    end if
    SLL_DEALLOCATE(mesh, ierr)
  end subroutine delete_hex_logical_mesh_2d

  
#undef TEST_PRESENCE_AND_ASSIGN_VAL

end module hex_logical_meshes
