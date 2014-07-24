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


module sll_particle_group_1d1v_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "particle_representation.h"

!  use sll_particle_representations
  use sll_logical_meshes
  implicit none

  type :: sll_particle_group_1d1v
     sll_int32  :: number_particles
     sll_int32  :: active_particles
     sll_int32  :: guard_list_size
     sll_real64 :: qoverm
     type(sll_logical_mesh_1d), pointer                 :: mesh
     class(sll_particle_1d1v), dimension(:), pointer       :: p_list
     class(sll_particle_1d1v_guard), dimension(:), pointer :: p_guard
    contains


    interface getX
!            sll_particle_1d1v_idxx_getX_pointer
!            sll_particle_1d1v_idxx_getX_pidx
!            sll_particle_1d1v_getX_mesh_pointer
    endinterface



  end type sll_particle_group_1d1v

  interface sll_delete
     module procedure delete_particle_2d_group
  end interface sll_delete



contains

  function new_particle_1d1v_group( &
       num_particles,       &
       particle_array_size, &
       guard_list_size,     &
       qoverm,              &
       mesh ) result(res)

    type(sll_particle_group_2d), pointer :: res
    sll_int32,  intent(in) :: num_particles
    sll_int32,  intent(in) :: particle_array_size
    sll_int32,  intent(in) :: guard_list_size
    sll_real64, intent(in) :: qoverm
    type(sll_logical_mesh_2d), pointer :: mesh
    sll_int32 :: ierr

    if( num_particles > particle_array_size ) then
       print *, 'new_particle_2d_group(): ERROR,  num_particles should not ', &
            'be greater than the memory size requested, particle_array_size.'
       STOP
    end if

    SLL_ALLOCATE( res, ierr )
    res%number_particles = num_particles
    res%active_particles = num_particles
    res%guard_list_size  = guard_list_size
    res%qoverm           = qoverm

    SLL_ALLOCATE( res%p_list(particle_array_size), ierr )
    SLL_ALLOCATE( res%p_guard(guard_list_size), ierr )

    if (.not.associated(mesh) ) then
       print*, 'error: passed mesh not associated'
    endif
    res%mesh => mesh

  end function new_particle_2d_group

  subroutine delete_particle_2d_group(p_group)
    type(sll_particle_group_2d), pointer :: p_group
    sll_int32 :: ierr

    if(.not. associated(p_group) ) then
       print *, 'delete_particle_group_2d(): ERROR, passed group was not ', &
            'associated.'
    end if
    SLL_DEALLOCATE(p_group%p_list, ierr)
    SLL_DEALLOCATE(p_group%p_guard, ierr)
    SLL_DEALLOCATE(p_group, ierr)
  end subroutine delete_particle_2d_group


!  subroutine set_group_particle_values( &
!                    p_group,   &
!                    p_index,   &
!                    x,  y,     &
!                    vx, vy, q  )
!
!    sll_real64, intent(in) :: x, y, vx, vy
!    sll_real32, intent(in) :: q
!    sll_int64,  intent(in) :: p_index
!    sll_real32 :: off_x, off_y
!    sll_real64 :: tmp1, tmp2, xmin, ymin, rdx, rdy
!    sll_int32  :: ic_x, ic_y, ncx
!    type(sll_particle_group_2d), pointer :: p_group
!
!    SLL_ASSERT( (p_index >= 1).and.(p_index <= p_group%number_particles) )
!
!    ncx  = group%mesh%num_cells1
!    xmin = group%mesh%eta1_min
!    ymin = group%mesh%eta2_min
!    rdx = 1._f64/group%mesh%delta_eta1
!    rdy = 1._f64/group%mesh%delta_eta2
!
!    SET_PARTICLE_VALUES(p_group%p_list(p_index),x,y,vx,vy,q,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
!!!$    call set_particle_values( p_group%p_list(p_index), &
!!!$                              x,  y,       &
!!!$                              vx, vy, q,   &
!!!$                              p_group%mesh )
!
!  end subroutine set_group_particle_values
!
end module sll_particle_group_1d1v_module
