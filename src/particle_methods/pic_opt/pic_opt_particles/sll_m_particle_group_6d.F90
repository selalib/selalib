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


module sll_m_particle_group_6d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_3d

  use sll_m_particle_representations, only: &
    sll_t_particle_6d

!#ifdef _OPENMP
!  use omp_lib, only: &
!    omp_get_num_threads, &
!    omp_get_thread_num
!
!#endif
  implicit none

  public :: &
    sll_s_particle_6d_group_init, &
    sll_s_particle_6d_group_free, &
    sll_t_particle_group_6d

  private

  type :: sll_t_particle_group_6d
     sll_int32  :: number_particles
     sll_int32  :: active_particles
     sll_int32  :: guard_list_size
     sll_int32, allocatable :: num_postprocess_particles(:)
     sll_real64 :: qoverm 
     type(sll_t_cartesian_mesh_3d), pointer :: mesh
     type(sll_t_particle_6d),       pointer :: p_list(:)
  end type sll_t_particle_group_6d
  
contains

  subroutine sll_s_particle_6d_group_init( &
       res,                 &
       num_particles,       &
       particle_array_size, &
       guard_list_size,     &
       qoverm,              &
       mesh ) 

    type(sll_t_particle_group_6d)         :: res
    sll_int32,  intent(in)                :: num_particles
    sll_int32,  intent(in)                :: particle_array_size
    sll_int32,  intent(in)                :: guard_list_size
    sll_real64, intent(in)                :: qoverm
    type(sll_t_cartesian_mesh_3d), target :: mesh

    sll_int32 :: ierr
    sll_int32 :: n_thread
    sll_int32 :: thread_id
!    sll_int32 :: nn

    if( num_particles > particle_array_size ) then
       print *, 'sll_f_new_particle_6d_group(): ERROR,  num_particles should not ', &
            'be greater than the memory size requested, particle_array_size.'
       stop
    end if

    res%number_particles = num_particles
    res%active_particles = num_particles
    res%guard_list_size  = guard_list_size
    res%qoverm           = qoverm

    SLL_ALLOCATE( res%p_list(particle_array_size), ierr )

    n_thread  = 1
    thread_id = 0
!
!    !$omp parallel PRIVATE(thread_id)
!#ifdef _OPENMP
!    thread_id = OMP_GET_THREAD_NUM()
!    if (thread_id ==0) then
!       n_thread  = OMP_GET_NUM_THREADS()
!    endif
!#endif
!    !$omp end parallel
!
!    nn = guard_list_size/n_thread
!    SLL_ALLOCATE( res%p_guard(1:n_thread), ierr)
    SLL_ALLOCATE( res%num_postprocess_particles(1:n_thread), ierr)
!
!    !$omp parallel PRIVATE(thread_id)
!#ifdef _OPENMP
!    thread_id = OMP_GET_THREAD_NUM()
!#endif
!    SLL_ALLOCATE( res%p_guard(thread_id+1)%g_list(1:nn),ierr)
!    !$omp end parallel
    
    res%mesh => mesh

  end subroutine sll_s_particle_6d_group_init

  subroutine sll_s_particle_6d_group_free(p_group)
    type(sll_t_particle_group_6d), pointer :: p_group
!    sll_int32 :: ierr
!    sll_int32 :: thread_id

!    if(.not. associated(p_group) ) then
!       print *, 'delete_particle_group_2d(): ERROR, passed group was not ', &
!            'associated.'
!    end if

!    thread_id = 0
!#ifdef _OPENMP
!    thread_id = OMP_GET_THREAD_NUM()
!#endif
    deallocate(p_group%num_postprocess_particles)
!    SLL_DEALLOCATE(p_group%p_guard(thread_id+1)%g_list, ierr)
    deallocate(p_group%p_list)
!    SLL_DEALLOCATE(p_group%p_guard, ierr)
!    SLL_DEALLOCATE(p_group, ierr)
  end subroutine sll_s_particle_6d_group_free


end module sll_m_particle_group_6d
