!**************************************************************
!  Copyright INRIA, CNRS
!  Authors : 
!     Pierre Navaro 
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


program test_layout_output
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"

  use sll_collective, only: sll_boot_collective, &
       sll_halt_collective
  
  use hdf5
  use sll_hdf5_io_parallel, only: sll_hdf5_file_create, &
       sll_hdf5_write_array, &
       sll_hdf5_file_close
  use sll_remapper


  implicit none

  ! ni, nj, nk: global sizes
  integer , parameter                       :: ni = 128
  integer , parameter                       :: nj = 64
  integer , parameter                       :: nk = 32
  ! Local sizes
  integer                                   :: loc_sz_i_init
  integer                                   :: loc_sz_j_init
  integer                                   :: loc_sz_k_init

  ! the process mesh
  integer                                   :: npi
  integer                                   :: npj
  integer                                   :: npk
  integer                                   :: error
  integer                                   :: myrank
  sll_int64                                 :: colsz        ! collective size

  type(layout_3D), pointer                  :: layout

  sll_real64                                :: tcpu1
  sll_real64                                :: tcpu2

  sll_int32                                 :: file_id
 
  character(len=9), parameter               :: filename = "layout.h5"

  integer, parameter :: rank = 3
  integer(HSIZE_T),  dimension(rank) :: dims = (/ni,nj,nk/)
  sll_int32, dimension(:,:,:), allocatable :: array

  integer(HSSIZE_T), dimension(rank) :: offset 

  ! Boot parallel environment
  call sll_boot_collective()
  colsz  = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)

  tcpu1 = MPI_WTIME()

  if( myrank .eq. 0) then
     print *, ' '
     print *, '--------------- layout output test ---------------------'
     print *, ' '
     print *, 'Running a test on ', colsz, 'processes'
     call flush(6)
  end if

  if (.not. is_power_of_two(colsz)) then     
     print *, 'This test needs to run in a number of processes which is ',&
              'a power of 2.'
     call sll_halt_collective()
     stop
  end if

  layout  => new_layout_3D( sll_world_collective )        
  call two_power_rand_factorization(colsz, npi, npj, npk)

  if( myrank .eq. 0 ) &
     print *, '3D layout configuration: ', npi, npj, npk

  call initialize_layout_with_distributed_3D_array( &
                    ni, nj, nk, npi, npj, npk, layout )
     
  call compute_local_sizes( layout, loc_sz_i_init, &
                                    loc_sz_j_init, &
                                    loc_sz_k_init )        

  ! initialize the local data    
  print *, myrank, 'Printing layout1: '
  call sll_view_lims_3D( layout )

  SLL_ALLOCATE(array(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init),error)
 
  array = myrank

  offset(1) = get_layout_3D_i_min( layout, myrank ) - 1
  offset(2) = get_layout_3D_j_min( layout, myrank ) - 1
  offset(3) = get_layout_3D_k_min( layout, myrank ) - 1

  call sll_hdf5_file_create('layout3d.h5',file_id,error)
  call sll_hdf5_write_array(file_id,dims,offset,dble(array),'array',error)
  call sll_hdf5_file_close(file_id,error)

  call delete_layout_3D( layout )
  SLL_DEALLOCATE_ARRAY(array, error)

  tcpu2 = MPI_WTIME()
  if (myrank == 0) &
   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*colsz

  call sll_halt_collective()
  
contains

  subroutine two_power_rand_factorization(n, n1, n2, n3)
    sll_int64, intent(in) :: n
    integer, intent(out)  :: n1, n2, n3
    integer               :: expo, expo1, expo2, expo3
    sll_real64            :: rand_real
    if (.not.is_power_of_two(colsz)) then   
       print*, 'The number of processors must be a power of 2'
       call sll_halt_collective()
       stop
    endif 
    expo = int(log(real(n))/log(2.))  
    call random_number(rand_real)
    expo1 = int(rand_real*expo)
    call random_number(rand_real)
    expo2 = int(rand_real*(expo-expo1))
    expo3 = expo - (expo1+expo2)
    n1 = 2**expo1
    n2 = 2**expo2
    n3 = 2**expo3
  end subroutine two_power_rand_factorization

end program test_layout_output
