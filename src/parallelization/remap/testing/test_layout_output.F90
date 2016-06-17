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
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use hdf5, only: &
    hid_t, &
    hsize_t, &
    hssize_t

  use iso_fortran_env, only: &
    output_unit

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_hdf5_io_parallel, only: &
    sll_o_hdf5_file_create, &
    sll_o_hdf5_write_array

  use sll_m_hdf5_io_serial, only: &
    sll_o_hdf5_file_close

  use sll_m_remapper, only: &
    sll_o_compute_local_sizes, &
    sll_o_get_layout_i_min, &
    sll_o_get_layout_j_min, &
    sll_o_get_layout_k_min, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_3d, &
    sll_f_new_layout_3d, &
    sll_o_delete, &
    sll_o_view_lims

  use sll_m_utilities, only: &
    sll_f_is_power_of_two

  use sll_mpi, only: &
    mpi_wtime

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! ni, nj, nk: global sizes
  sll_int32 , parameter                       :: ni = 128
  sll_int32 , parameter                       :: nj = 64
  sll_int32 , parameter                       :: nk = 32
  ! Local sizes
  sll_int32                                   :: loc_sz_i_init
  sll_int32                                   :: loc_sz_j_init
  sll_int32                                   :: loc_sz_k_init

  ! the process mesh
  sll_int32                                   :: npi
  sll_int32                                   :: npj
  sll_int32                                   :: npk
  sll_int32                                   :: error
  sll_int32                                   :: myrank
  sll_int32                                   :: comm
  sll_int64                                 :: colsz        ! collective size

  type(sll_t_layout_3d), pointer                  :: layout

  sll_real64                                :: tcpu1
  sll_real64                                :: tcpu2

  !sll_int32                                 :: file_id
  integer(hid_t)                            :: hfile_id
 
  character(len=9), parameter               :: filename = "layout.h5"

  sll_int32, parameter :: rank = 3
  integer(HSIZE_T),  dimension(rank) :: dims = (/int(ni,HSIZE_T),int(nj,HSIZE_T),int(nk,HSIZE_T)/)
  sll_int32, dimension(:,:,:), allocatable :: array

  integer(HSSIZE_T), dimension(rank) :: offset 

  ! Boot parallel environment
  call sll_s_boot_collective()
  colsz  = int(sll_f_get_collective_size(sll_v_world_collective),i64)
  myrank = sll_f_get_collective_rank(sll_v_world_collective)

  tcpu1 = MPI_WTIME()

  if( myrank .eq. 0) then
     print *, ' '
     print *, '--------------- layout output test ---------------------'
     print *, ' '
     print *, 'Running a test on ', colsz, 'processes'
     flush( output_unit )
  end if

  if (.not. sll_f_is_power_of_two(colsz)) then     
     print *, 'This test needs to run in a number of processes which is ',&
              'a power of 2.'
     call sll_s_halt_collective()
     stop
  end if

  layout  => sll_f_new_layout_3d( sll_v_world_collective )        
  call two_power_rand_factorization(colsz, npi, npj, npk)

  if( myrank .eq. 0 ) &
     print *, '3D layout configuration: ', npi, npj, npk

  call sll_o_initialize_layout_with_distributed_array( &
                    ni, nj, nk, npi, npj, npk, layout )
     
  call sll_o_compute_local_sizes( layout, loc_sz_i_init, &
                                    loc_sz_j_init, &
                                    loc_sz_k_init )        

  ! initialize the local data    
  print *, myrank, 'Printing layout1: '
  call sll_o_view_lims( layout )

  SLL_ALLOCATE(array(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init),error)
 
  array = myrank

  offset(1) = int(sll_o_get_layout_i_min( layout, myrank ) - 1,HSIZE_T)
  offset(2) = int(sll_o_get_layout_j_min( layout, myrank ) - 1,HSIZE_T)
  offset(3) = int(sll_o_get_layout_k_min( layout, myrank ) - 1,HSIZE_T)

  comm   = sll_v_world_collective%comm
  call sll_o_hdf5_file_create('layout3d.h5',comm,hfile_id,error)
  call sll_o_hdf5_write_array(hfile_id,dims,offset,dble(array),'array',error)
  call sll_o_hdf5_file_close(hfile_id,error)

  call sll_o_delete( layout )
  SLL_DEALLOCATE_ARRAY(array, error)

  tcpu2 = MPI_WTIME()
  if (myrank == 0) &
   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*colsz

  call sll_s_halt_collective()
  
contains

  subroutine two_power_rand_factorization(n, n1, n2, n3)
    sll_int64, intent(in) :: n
    sll_int32, intent(out)  :: n1, n2, n3
    sll_int32               :: expo, expo1, expo2, expo3
    sll_real64            :: rand_real
    if (.not.sll_f_is_power_of_two(colsz)) then   
       print*, 'The number of processors must be a power of 2'
       call sll_s_halt_collective()
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
