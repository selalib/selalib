program remap_test
  use remapper
  use sll_collective
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
  implicit none

  ! Test of the 1D remapper takes a 1D array whose global size N,
  ! distributed among NP processors.
  integer, dimension(:), allocatable :: a
  integer, dimension(:,:,:), allocatable :: a3
  integer, dimension(:,:,:), allocatable :: b3

  ! Take a 3D array of dimensions 8X8X1
  integer, parameter                 :: total_sz_i = 8
  integer, parameter                 :: total_sz_j = 8
  integer, parameter                 :: total_sz_k = 1
  integer, parameter                 :: local_sz = 4 ! input, local size
  ! process mesh
  integer, parameter                 :: pi = 4
  integer, parameter                 :: pj = 4
  integer, parameter                 :: pk = 1

  ! Split it in  16 processes, each with a local chunk 2X2X1
  integer                            :: local_sz_i 
  integer                            :: local_sz_j 
  integer                            :: local_sz_k 
  integer                            :: ierr
  integer                            :: myrank
  integer                            :: colsz        ! collective size
  integer                            :: i,j,k
  integer                            :: i_min, i_max
  integer                            :: j_min, j_max
  integer                            :: k_min, k_max
  integer                            :: node
  integer, dimension(1:3)            :: gcoords
  ! Remap stuff
  type(layout_3D_t), pointer         :: conf3_init
  type(layout_3D_t), pointer         :: conf3_final

!  integer, dimension(:), allocatable :: sendbuf
  type(remap_plan_3D_t), pointer     :: rmp3

  print *, ' '
  print *, '--------------- REMAP test ---------------------'
  print *, ' '

  call flush()
  call sll_boot_collective()
  SLL_ALLOCATE( a(local_sz), ierr )

  local_sz_i = total_sz_i/pi
  local_sz_j = total_sz_j/pj
  local_sz_k = total_sz_k/pk
  SLL_ALLOCATE( a3(1:local_sz_i,1:local_sz_j,1:local_sz_k), ierr )
  SLL_ALLOCATE( b3(1:local_sz_i,1:local_sz_j,1:local_sz_k), ierr )
  myrank    = sll_get_collective_rank(sll_world_collective)
  colsz     = sll_get_collective_size(sll_world_collective)
!  total_size = colsz*local_sz!

!  conf_init  => new_layout_1D( sll_world_collective )
!  conf_final => new_layout_1D( sll_world_collective )

  conf3_init  => new_layout_3D( sll_world_collective )
  conf3_final => new_layout_3D( sll_world_collective )

  do k=0, pk-1
     do j=0, pj-1
        do i=0, pi-1
           node = i+pi*(j+pj*k) ! linear index of node
           i_min = i*local_sz_i + 1
           i_max = i*local_sz_i + local_sz_i
           j_min = j*local_sz_j + 1
           j_max = j*local_sz_j + local_sz_j
           k_min = k*local_sz_k + 1
           k_max = k*local_sz_k + local_sz_k
           call set_layout_i_min( conf3_init, node, i_min )
           call set_layout_i_max( conf3_init, node, i_max )
           call set_layout_j_min( conf3_init, node, j_min )
           call set_layout_j_max( conf3_init, node, j_max )
           call set_layout_k_min( conf3_init, node, k_min )
           call set_layout_k_max( conf3_init, node, k_max )
        end do
     end do
  end do
      
  call sll_view_lims_3D( conf3_init )

  ! Initialize the data. We use the information in the layout.
  do k=1, local_sz_k
     do j=1, local_sz_j
        do i=1, local_sz_i
           gcoords =  local_to_global_3D( conf3_init, (/i,j,k/) )
!           write (*,'(a,i4)') 'gcoords in rank: ', myrank
!           print *, gcoords(:)
!           call flush()
           a3(i,j,k) = gcoords(1) + &
                total_sz_i*((gcoords(2)-1) + total_sz_j*(gcoords(3)-1))
        end do
     end do
  end do

  write (*,'(a,i4)') 'From rank: ', myrank
  print *, a3(:,:,:)
  call flush()

  ! Initialize the final layout, in this case, just a transposition
  do k=0, pk-1
     do j=0, pj-1
        do i=0, pi-1
           node = i+pi*(j+pj*k) ! linear index of node
           i_min = i*local_sz_i + 1
           i_max = i*local_sz_i + local_sz_i
           j_min = j*local_sz_j + 1
           j_max = j*local_sz_j + local_sz_j
           k_min = k*local_sz_k + 1
           k_max = k*local_sz_k + local_sz_k
           call set_layout_i_min( conf3_final, node, j_min )
           call set_layout_i_max( conf3_final, node, j_max )
           call set_layout_j_min( conf3_final, node, i_min )
           call set_layout_j_max( conf3_final, node, i_max )
           call set_layout_k_min( conf3_final, node, k_min )
           call set_layout_k_max( conf3_final, node, k_max )
        end do
     end do
  end do
  
  call sll_view_lims_3D( conf3_final )      
  
  rmp3 => new_remap_plan_3D( conf3_init, conf3_final, INT32_SIZEOF(a3(1,1,1)) )
  call apply_remap_3D_int( rmp3, a3, b3 )
  print *, 'Remap operation completed.'
  write (*,'(a, i4)') 'the output data in rank: ', myrank
  print *, b3(:,:,:)
  call flush()
  
  call delete_layout_3D( conf3_init )
  call delete_layout_3D( conf3_final )

  call sll_collective_barrier(sll_world_collective)
  call sll_halt_collective()
  print *, 'TEST COMPLETE'




end program remap_test
