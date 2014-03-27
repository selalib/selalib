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


module sll_particle_group_2d_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_particle_representations
  implicit none

  type :: sll_particle_group_2d
     sll_int64 :: number_particles! peut etre a faire en SLL_PRIV
     sll_int64 :: active_particles
     sll_int32 :: guard_list_size
     type(sll_particle_2d), dimension(:), pointer       :: p_list
     type(sll_particle_2d_guard), dimension(:), pointer :: p_guard
  end type sll_particle_group_2d

  interface sll_delete
     module procedure delete_particle_2d_group
  end interface sll_delete

contains

  function new_particle_2d_group( &
       num_particles, &
       particle_array_size, &
       guard_list_size ) result(res)

    type(sll_particle_group_2d), pointer :: res
    sll_int64, intent(in) :: num_particles
    sll_int64, intent(in) :: particle_array_size
    sll_int64, intent(in) :: guard_list_size
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
    SLL_ALLOCATE( res%p_list(particle_array_size), ierr )
    SLL_ALLOCATE( res%p_guard(guard_list_size), ierr )
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

end module sll_particle_group_2d_module
