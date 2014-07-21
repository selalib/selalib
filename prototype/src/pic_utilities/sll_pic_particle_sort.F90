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


module sll_particle_sort_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_logical_meshes
  use sll_particle_representations
  use sll_particle_group_2d_module
  implicit none

  type :: sll_particle_sorter_2d
     sll_int32, dimension(:), pointer   :: pa
     sll_int32, dimension(:), pointer   :: pa_save
     type(sll_logical_mesh_2d), pointer :: mesh
     sll_int32                          :: num_cells
  end type sll_particle_sorter_2d
  
  interface sll_delete
     module procedure delete_particle_sorter_2d
  end interface sll_delete

contains

  function sll_new_particle_sorter_2d( mesh ) result(res)
    type(sll_particle_sorter_2d), pointer :: res
    type(sll_logical_mesh_2d), pointer    :: mesh
    sll_int32 :: ierr
    sll_int32 :: ncx
    sll_int32 :: ncy

    if( .not. associated(mesh) ) then
       print *, 'ERROR, sll_new_particle_sorter_2d(): passed mesh is not ', &
            'associated.'
       stop
    end if

    ncx = mesh%num_cells1
    ncy = mesh%num_cells2

    SLL_ALLOCATE(res,ierr)
    SLL_ALLOCATE(res%pa(ncx*ncy+1),ierr)
    SLL_ALLOCATE(res%pa_save(ncx*ncy+1),ierr)

    res%mesh => mesh
    res%num_cells = ncx*ncy
  end function sll_new_particle_sorter_2d

  subroutine sll_sort_particles_2d( sorter, group )
    type(sll_particle_sorter_2d), pointer :: sorter
    type(sll_particle_group_2d), pointer  :: group
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: k
    sll_int32 :: N     ! number of particles
    sll_int32 :: num_cells
    sll_int32 :: current_cell
    sll_int32 :: index_in
    sll_int32 :: index_out
    sll_int32 :: index_stop
    type(sll_particle_2d), dimension(:), pointer :: p
    type(sll_particle_2d)                        :: p_tmp
    sll_int32, dimension(:), pointer             :: pa
    sll_int32, dimension(:), pointer             :: pa_save
    ! make sure that the meshes are the same
    if( .not. associated(sorter%mesh, target=group%mesh) ) then
       print *, 'ERROR, sll_sort_particles_2d(): mesh passed to sorter ', &
            'and particle group mesh are not the same. Code will not stop ', &
            'but bad things may happen...'
    end if

    N          = group%number_particles
    num_cells  = sorter%num_cells
    p          => group%p_list
    pa         => sorter%pa(:)
    pa_save    => sorter%pa_save(:)
    pa(:)      = 0
    pa_save(:) = 0
    ! Count how many particles are there in each cell, name this 'cell count'.
    do i=1,N
       pa(p(i)%ic) = pa(p(i)%ic) + 1
    end do

    ! Convert the 'cell count' into an allocation within pa and save a copy
    ! of this allocation in pa_save. The allocation is just the sum of 
    ! particles up to, and including the previous cell. Thus pa(i) stores
    ! the number of particles up to cell i-1.
    k = 0
    do i=1,num_cells+1
       j          = pa(i)
       pa_save(i) = k
       pa(i)      = k
       k          = k + j
    end do

    i = 1
    do while (i < num_cells + 1)
       if( pa(i) >= pa_save(i+1) ) then
          i = i + 1 ! current cell is done, process next cell.
       else
          ! The current cell still contain sunsorted particles, get the next
          ! unsorted particle in the current cell for the next sorting cycle.
          ! pa(i) contains the index in which the next particle in cell 'i'
          ! should be put.
          index_in   = pa(i) + 1 ! which particle to process
          index_stop = pa(i) + 1 ! when to stop the swaps
          do
             ! Figure out where to store the input particle. Update the
             ! allocation accordingly.
             current_cell     = p(index_in)%ic
             pa(current_cell) = pa(current_cell) + 1
             index_out        = pa(current_cell)
             if( index_out .ne. index_stop ) then
                p_tmp        = p(index_out)
                p(index_out) = p(index_in)
                p(index_in)  = p_tmp
             else
                exit
             end if
          end do
       end if
    end do
  end subroutine sll_sort_particles_2d


  subroutine delete_particle_sorter_2d( sorter )
    type(sll_particle_sorter_2d), pointer :: sorter
  end subroutine delete_particle_sorter_2d


end module sll_particle_sort_module
