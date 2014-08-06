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
#include "particle_representation.h"

!  use sll_particle_representations
  use sll_logical_meshes
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  type :: sll_particle_group_2d
     sll_int32  :: number_particles! peut etre a faire en SLL_PRIV
     sll_int32  :: active_particles! tout ça doit passer en 32
     sll_int32  :: guard_list_size! tout ça doit passer en 32
     ! an array indexed by the thread number, of the number of particles
     ! to post-process after the main loop
     sll_int32, dimension(:), pointer :: num_postprocess_particles
     sll_real64 :: qoverm 
     type(sll_logical_mesh_2d), pointer                 :: mesh
     type(sll_particle_2d), dimension(:), pointer       :: p_list
     type(sll_particle_2d_guard_ptr), dimension(:), pointer :: p_guard
  end type sll_particle_group_2d
  
  interface sll_delete
     module procedure delete_particle_2d_group
  end interface sll_delete

contains

  function new_particle_2d_group( &
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
    sll_int32 :: n_thread
    sll_int32 :: thread_id
    sll_int32 :: nn

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

    n_thread  = 1
    thread_id = 0

    !$omp parallel default(SHARED) PRIVATE(thread_id)
#ifdef _OPENMP
    thread_id = OMP_GET_THREAD_NUM()
    if (thread_id ==0) then
       n_thread  = OMP_GET_NUM_THREADS()
    endif
#endif
    !$omp end parallel

    nn = guard_list_size/n_thread
    SLL_ALLOCATE( res%p_guard(1:n_thread), ierr)
    SLL_ALLOCATE( res%num_postprocess_particles(1:n_thread), ierr)

    !$omp parallel default(SHARED) PRIVATE(thread_id)
#ifdef _OPENMP
    thread_id = OMP_GET_THREAD_NUM()
#endif
    SLL_ALLOCATE( res%p_guard(thread_id+1)%g_list(1:nn),ierr)
    !$omp end parallel
    
    if (.not.associated(mesh) ) then
       print*, 'error: passed mesh not associated'
    endif
    res%mesh => mesh

  end function new_particle_2d_group

  subroutine delete_particle_2d_group(p_group)
    type(sll_particle_group_2d), pointer :: p_group
    sll_int32 :: ierr
    sll_int32 :: thread_id

    if(.not. associated(p_group) ) then
       print *, 'delete_particle_group_2d(): ERROR, passed group was not ', &
            'associated.'
    end if

    thread_id = 0
#ifdef _OPENMP
    thread_id = OMP_GET_THREAD_NUM()
#endif
    SLL_DEALLOCATE(p_group%num_postprocess_particles, ierr  )
    SLL_DEALLOCATE(p_group%p_guard(thread_id+1)%g_list, ierr)
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
end module sll_particle_group_2d_module
