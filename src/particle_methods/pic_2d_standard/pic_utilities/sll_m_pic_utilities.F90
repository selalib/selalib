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

module sll_m_pic_utilities
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_accumulators.h"
  use sll_m_accumulators
  use sll_m_particle_group_4d
  use sll_m_particle_group_2d

#ifdef _OPENMP
   use omp_lib
#endif
implicit none
#ifdef _OPENMP
   logical :: openmp_st
#endif
!!    !$ openmp_st = OMP_IN_PARALLEL()
!!    print*, 'USE of omp', openmp_st
!!    !$omp end parallel

contains

  subroutine sll_first_charge_accumulation_2d( p_group, q_accumulator )
    type(sll_particle_group_4d), pointer                       :: p_group
    type(sll_charge_accumulator_2d_ptr), dimension(:), pointer :: q_accumulator
    type(sll_particle_4d), dimension(:), pointer :: p
    type(sll_charge_accumulator_2d), pointer :: q_accum
    sll_int64  :: i
    sll_int64  :: num_particles
    sll_real64 :: tmp1
    sll_real64 :: tmp2
    sll_int32  :: thread_id
    
    SLL_ASSERT( associated(p_group) .and. associated(q_accumulator))
    num_particles =  int(p_group%number_particles,i64)
    p             => p_group%p_list

    !$omp parallel PRIVATE(thread_id, tmp1, tmp2, q_accum)
#ifdef _OPENMP
    thread_id = OMP_GET_THREAD_NUM()
#else
    thread_id = 0
#endif
    q_accum  => q_accumulator(thread_id+1)%q
    !$omp do
    do i=1,num_particles
       SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p(i),tmp1,tmp2)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine sll_first_charge_accumulation_2d

  subroutine sll_first_charge_accumulation_2d_CS( p_group, q_accumulator )
! ----  Remember : _CS is for use of Cubic Splines  ----
! ------------------------------------------------------
    type(sll_particle_group_4d), pointer         :: p_group
    type(sll_charge_accumulator_2d_CS_ptr), dimension(:), pointer  :: q_accumulator
    type(sll_particle_4d), dimension(:), pointer :: p
    type(sll_charge_accumulator_2d_CS), pointer :: q_accum
    sll_int64  :: i
    sll_int64  :: num_particles
    sll_real32 :: tmp(1:4,1:2), temp
    sll_int32  :: thread_id

    SLL_ASSERT( associated(p_group) .and. associated(q_accumulator))
    num_particles =  int(p_group%number_particles,i64)
    p             => p_group%p_list
    
    !$omp parallel default(SHARED) PRIVATE(thread_id, tmp, temp, q_accum)
#ifdef _OPENMP
    thread_id = OMP_GET_THREAD_NUM()
#else
    thread_id = 0
#endif
    q_accum  => q_accumulator(thread_id+1)%q
    !$omp do
    do i=1,num_particles
       SLL_ACCUMULATE_PARTICLE_CHARGE_CS(q_accum,p(i),tmp,temp)
    end do
    !$omp end do
    !$omp end parallel
    
  end subroutine sll_first_charge_accumulation_2d_CS


!!$ - - - -  for the GUIDING CENTER model   - - - -
  subroutine sll_first_gc_charge_accumulation_2d( p_group, q_accumulator )
    type(sll_particle_group_2d), pointer                       :: p_group
    type(sll_charge_accumulator_2d_ptr), dimension(:), pointer :: q_accumulator
    type(sll_particle_2d), dimension(:), pointer :: p
    type(sll_charge_accumulator_2d), pointer :: q_accum
    sll_int64  :: i
    sll_int64  :: num_particles
    sll_real64 :: tmp1
    sll_real64 :: tmp2
    sll_int32  :: thread_id
    
    SLL_ASSERT( associated(p_group) .and. associated(q_accumulator))
    num_particles =  int(p_group%number_particles,i64)
    p             => p_group%p_list

    !$omp parallel default(SHARED) PRIVATE(thread_id, tmp1, tmp2, q_accum)
#ifdef _OPENMP
    thread_id = OMP_GET_THREAD_NUM()
#else
    thread_id = 0
#endif
    q_accum  => q_accumulator(thread_id+1)%q
    !$omp do
    do i=1,num_particles
       SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p(i),tmp1,tmp2)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine sll_first_gc_charge_accumulation_2d

  subroutine sll_first_gc_charge_accumulation_2d_CS( p_group, q_accumulator )
    type(sll_particle_group_2d), pointer         :: p_group
    type(sll_charge_accumulator_2d_CS_ptr), dimension(:), pointer  :: q_accumulator
    type(sll_particle_2d), dimension(:), pointer :: p
    type(sll_charge_accumulator_2d_CS), pointer :: q_accum
    sll_int64  :: i
    sll_int64  :: num_particles
    sll_real32 :: tmp(1:4,1:2), temp
    sll_int32  :: thread_id

    SLL_ASSERT( associated(p_group) .and. associated(q_accumulator))
    num_particles =  int(p_group%number_particles,i64)
    p             => p_group%p_list
    
    !$omp parallel default(SHARED) PRIVATE(thread_id, tmp, temp, q_accum)
#ifdef _OPENMP
    thread_id = OMP_GET_THREAD_NUM()
#else
    thread_id = 0
#endif
    q_accum  => q_accumulator(thread_id+1)%q
    !$omp do
    do i=1,num_particles
       SLL_ACCUMULATE_PARTICLE_CHARGE_CS(q_accum,p(i),tmp,temp)
    end do
    !$omp end do
    !$omp end parallel
    
  end subroutine sll_first_gc_charge_accumulation_2d_CS



end module sll_m_pic_utilities
