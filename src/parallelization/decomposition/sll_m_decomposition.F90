!**************************************************************
!
!  Domain decomposition module for Selalib (1D, ..., 6d)
!
!  Some design aspects are borrowed from sll_remap.F90.
!
!  Note:
!  For memory reasons, this code is implemented separately
!  from the comm/port model provided by the existing module
!  located at "point_to_point_communications".
!  In 6d, we cannot afford allocating 2*12=24 buffers to
!  implement the non-blocking ghost cell exchange using
!  the comm/port model.
!
!  Copyright 2015
!    Klaus Reuter, khr@mpcdf.mpg.de
!    Max Planck Computing and Data Facility (MPCDF)
!    Garching, Germany
!  All rights reserved.
!
!--------------------------------------------------------------
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
!
!**************************************************************


! Enable the SLL_ASSERT() calls locally during development.
!#define DEBUG


!> @ingroup decomposition
!> @brief
!> Module providing data structures and tools to implement domain decompositions.
!> @author
!> Klaus Reuter, Max Planck Computing and Data Facility (MPCDF)
module sll_m_decomposition
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
   use sll_m_collective
   use mpi

   implicit none


   !> @brief    Information on the cartesian process topology.
   type :: cartesian_topology_6d
      ! topology-associated MPI communicator
      sll_int32 :: comm
      ! array containing the number of MPI processes to use along the n-th dimension
      sll_int32 :: procs(6)
      ! array indicating if a periodic boundary condition exists at the n-th dimension
      logical   :: periodic(6)
      ! coordinates of the current MPI process within the cartesian topology
      sll_int32 :: coords(6)
      ! MPI ranks of the topological neighbors of the current MPI process
      sll_int32 :: neighbors(12)
   end type cartesian_topology_6d



   !> @brief    Index limits local to a MPI process.
   type :: decomposition_local_6d
      ! --- primary: local array extents usable for computation, width of the halo
      sll_int32, dimension(:) :: mn(6)  ! min index, w/o halo
      sll_int32, dimension(:) :: mx(6)  ! max index, w/o halo
      sll_int32, dimension(:) :: hw(6)  ! halo width
      ! --- derived: total allocated extents of the local array, including halo cells
      sll_int32, dimension(:) :: lo(6)  ! min index, w/ halo
      sll_int32, dimension(:) :: hi(6)  ! max index, w/ halo
      ! --- derived: auxiliary variables, array sizes
      sll_int32, dimension(:) :: nw(6)  ! net width of the array, w/o halo
      sll_int32, dimension(:) :: gw(6)  ! gross width of the array, w/ halo
      ! --- derived: ranges to copy between send (tx) and receive (rx) buffers and arrays
      sll_int32, dimension(:) :: tx_lolo(6)  ! lower/left neighbor, lower bound
      sll_int32, dimension(:) :: tx_lohi(6)  ! lower/left neighbor, upper bound
      sll_int32, dimension(:) :: tx_hilo(6)  ! upper/right neighbor, lower bound
      sll_int32, dimension(:) :: tx_hihi(6)  ! upper/right neighbor, upper bound
      sll_int32, dimension(:) :: rx_lolo(6)
      sll_int32, dimension(:) :: rx_lohi(6)
      sll_int32, dimension(:) :: rx_hilo(6)
      sll_int32, dimension(:) :: rx_hihi(6)
   end type decomposition_local_6d


   !> @brief    Global array size information and local information.
   type :: decomposition_6d
      sll_int32, dimension(:) :: global(6)
      type(decomposition_local_6d) :: local
   end type decomposition_6d


   interface new_cartesian_topology
      module procedure new_cartesian_topology_6d
   end interface new_cartesian_topology


   interface new_cartesian_domain_decomposition
      module procedure new_cartesian_domain_decomposition_6d
   end interface new_cartesian_domain_decomposition


   interface apply_halo_exchange
      module procedure apply_halo_exchange_6d_real64
   end interface apply_halo_exchange


   contains


   function new_cartesian_topology_6d(top_collective, procs_per_dimension, periodic)
      type(cartesian_topology_6d), pointer :: new_cartesian_topology_6d
      type(sll_collective_t), intent(in) :: top_collective
      integer, parameter :: nd=6
      sll_int32, intent(in) :: procs_per_dimension(nd)
      logical, intent(in) :: periodic(nd)
      ! disallow reordering of MPI ranks with MPI_Cart_create()
      logical, parameter :: reorder = .false.
      sll_int32 :: i, ierr

      SLL_ALLOCATE(new_cartesian_topology_6d, ierr)

      new_cartesian_topology_6d%procs = procs_per_dimension
      new_cartesian_topology_6d%periodic = periodic

      ! create a cartesian process topology, return a new communicator
      call MPI_Cart_create(top_collective%comm, nd,&
                           new_cartesian_topology_6d%procs,&
                           new_cartesian_topology_6d%periodic,&
                           reorder,&
                           new_cartesian_topology_6d%comm,&
                           ierr)
      SLL_ASSERT(ierr == MPI_SUCCESS)

      ! query the coordinates of the current process within the cartesian topology
      new_cartesian_topology_6d%coords = -1
      call MPI_Cart_get(new_cartesian_topology_6d%comm, nd,&
                        new_cartesian_topology_6d%procs,&
                        new_cartesian_topology_6d%periodic,&
                        new_cartesian_topology_6d%coords,&
                        ierr)
      SLL_ASSERT(ierr == MPI_SUCCESS)

      ! determine the neighbors within the cartesian topology
      new_cartesian_topology_6d%neighbors = -1
      do i=1,nd
         call MPI_Cart_shift(new_cartesian_topology_6d%comm, i-1, 1, &
                             new_cartesian_topology_6d%neighbors(2*i-1), &
                             new_cartesian_topology_6d%neighbors(2*i), &
                             ierr)
         SLL_ASSERT(ierr == MPI_SUCCESS)
      enddo
   end function


   function new_cartesian_domain_decomposition_6d(topology, grid_size, halo_width)
      type(decomposition_6d), pointer :: new_cartesian_domain_decomposition_6d
      type(cartesian_topology_6d), pointer, intent(in) :: topology
      integer, parameter :: nd=6
      sll_int32, intent(in) :: grid_size(nd)
      sll_int32, intent(in) :: halo_width(nd)
      sll_int32 :: i, ierr
      sll_int32 :: lp, l0, l1

      ! --- initial checks
      do i=1,nd
         SLL_ASSERT( halo_width(i) >= 0 )
      end do
      do i=1,nd
         ! for the moment, the global grid must be divisible evenly among the MPI processes
         SLL_ASSERT( mod(grid_size(i), topology%procs(i)) == 0 )
      end do

      ! --- copy and calculate all the necessary index values ---
      SLL_ALLOCATE(new_cartesian_domain_decomposition_6d, ierr)

      new_cartesian_domain_decomposition_6d%global =  grid_size

      ! --- loop over dimensions and compute index values for each dimension
      do i=1,nd
         ! compute the local number of grid points
         lp = grid_size(i) / topology%procs(i)
         ! compute the lower local index bound (the coords array starts at zero)
         l0 = 1 + topology%coords(i) * lp
         ! compute the upper local index bound (the coords array starts at zero)
         l1 = (topology%coords(i) + 1) * lp

         SLL_ASSERT( lp/2 >= halo_width(i) )

         new_cartesian_domain_decomposition_6d%local%mn(i) = l0
         new_cartesian_domain_decomposition_6d%local%mx(i) = l1
         new_cartesian_domain_decomposition_6d%local%hw(i) = halo_width(i)
         new_cartesian_domain_decomposition_6d%local%lo(i) = l0 - halo_width(i)
         new_cartesian_domain_decomposition_6d%local%hi(i) = l1 + halo_width(i)

         new_cartesian_domain_decomposition_6d%local%tx_lolo(i) = l0
         new_cartesian_domain_decomposition_6d%local%tx_lohi(i) = l0 + (halo_width(i) - 1)
         new_cartesian_domain_decomposition_6d%local%tx_hilo(i) = l1 - (halo_width(i) - 1)
         new_cartesian_domain_decomposition_6d%local%tx_hihi(i) = l1
         new_cartesian_domain_decomposition_6d%local%rx_lolo(i) = l0 - halo_width(i)
         new_cartesian_domain_decomposition_6d%local%rx_lohi(i) = l0 - 1
         new_cartesian_domain_decomposition_6d%local%rx_hilo(i) = l1 + 1
         new_cartesian_domain_decomposition_6d%local%rx_hihi(i) = l1 + halo_width(i)
      end do

      ! compute net array width
      new_cartesian_domain_decomposition_6d%local%nw = 1 + &
         new_cartesian_domain_decomposition_6d%local%mx - new_cartesian_domain_decomposition_6d%local%mn

      ! compute gross array width
      new_cartesian_domain_decomposition_6d%local%gw = 1 + &
         new_cartesian_domain_decomposition_6d%local%hi - new_cartesian_domain_decomposition_6d%local%lo
   end function


   subroutine copy_array_to_buffer_6d_real64(arr, arr_lo, arr_hi, buf, ranges)
      sll_int32, dimension(6), intent(in) :: arr_lo
      sll_int32, dimension(6), intent(in) :: arr_hi
      sll_real64, dimension(arr_lo(1):arr_hi(1), arr_lo(2):arr_hi(2), arr_lo(3):arr_hi(3),  &
                            arr_lo(4):arr_hi(4), arr_lo(5):arr_hi(5), arr_lo(6):arr_hi(6)), &
                                intent(in) :: arr
      sll_real64, dimension(:), intent(inout) :: buf
      sll_int32, dimension(6,2), intent(in) :: ranges
      sll_int32 :: idx,i,j,k,l,m,n

      idx=1
      do n=ranges(6,1),ranges(6,2)
         do m=ranges(5,1),ranges(5,2)
            do l=ranges(4,1),ranges(4,2)
               do k=ranges(3,1),ranges(3,2)
                  do j=ranges(2,1),ranges(2,2)
                     do i=ranges(1,1),ranges(1,2)
                        buf(idx) = arr(i,j,k,l,m,n)
                        idx = idx + 1
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

   end subroutine


   subroutine copy_buffer_to_array_6d_real64(buf, arr, arr_lo, arr_hi, ranges)
      sll_real64, dimension(:), intent(in) :: buf
      sll_int32, dimension(6), intent(in) :: arr_lo
      sll_int32, dimension(6), intent(in) :: arr_hi
      sll_real64, dimension(arr_lo(1):arr_hi(1), arr_lo(2):arr_hi(2), arr_lo(3):arr_hi(3), &
                            arr_lo(4):arr_hi(4), arr_lo(5):arr_hi(5), arr_lo(6):arr_hi(6)),&
                                intent(inout) :: arr
      sll_int32, dimension(6,2), intent(in) :: ranges
      sll_int32 :: idx,i,j,k,l,m,n

      idx=1
      do n=ranges(6,1),ranges(6,2)
         do m=ranges(5,1),ranges(5,2)
            do l=ranges(4,1),ranges(4,2)
               do k=ranges(3,1),ranges(3,2)
                  do j=ranges(2,1),ranges(2,2)
                     do i=ranges(1,1),ranges(1,2)
                        arr(i,j,k,l,m,n) = buf(idx)
                        idx = idx + 1
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine


   subroutine apply_halo_exchange_6d_real64(topo, decomp, arr)
      type(cartesian_topology_6d), intent(in) :: topo
      type(decomposition_6d), intent(in) :: decomp
      sll_real64, dimension(:,:,:,:,:,:), intent(inout) :: arr

      integer, save :: bufsize = 0
      sll_real64, dimension(:), allocatable, save :: sendbuf, recvbuf
      integer :: nxc

      integer, parameter :: nd = 6  ! we handle 6 dimensions
      integer, dimension(:,:) :: ranges(nd,2)  ! index ranges for the copy operations

      integer :: id, jd
      integer :: ierr

      ! --- loop over dimensions and exchange data between neighbors
      do id=1,nd
         ! calculate the number of items to be exchanged
         nxc = decomp%local%hw(id)
         do jd=1,nd
            if (jd == id) continue
            nxc = nxc * decomp%local%nw(jd)
         end do

         ! allocate memory if necessary
         if (nxc > bufsize) then
            if (allocated(sendbuf)) deallocate(sendbuf)
            if (allocated(recvbuf)) deallocate(recvbuf)
            allocate(sendbuf(nxc))
            allocate(recvbuf(nxc))
            bufsize = nxc
         end if

         ! --- (1) copy halo cells to the left neighbor
         ranges(:,1) = decomp%local%mn(:)
         ranges(:,2) = decomp%local%mx(:)
         ranges(id,1) = decomp%local%tx_lolo(id)
         ranges(id,2) = decomp%local%tx_lohi(id)
         call copy_array_to_buffer_6d_real64(arr, decomp%local%lo, decomp%local%hi, sendbuf, ranges)

         call MPI_Sendrecv(sendbuf, nxc, MPI_DOUBLE, topo%neighbors(2*id-1), 1, &
                           recvbuf, nxc, MPI_DOUBLE, topo%neighbors(2*id),   1, &
                           topo%comm, MPI_STATUS_IGNORE, ierr)

         ranges(:,1) = decomp%local%mn(:)
         ranges(:,2) = decomp%local%mx(:)
         ranges(id,1) = decomp%local%rx_hilo(id)
         ranges(id,2) = decomp%local%rx_hihi(id)
         call copy_buffer_to_array_6d_real64(recvbuf, arr, decomp%local%lo, decomp%local%hi, ranges)


         ! --- (2) copy halo cells to the right neighbor
         ranges(:,1) = decomp%local%mn(:)
         ranges(:,2) = decomp%local%mx(:)
         ranges(id,1) = decomp%local%tx_hilo(id)
         ranges(id,2) = decomp%local%tx_hihi(id)
         call copy_array_to_buffer_6d_real64(arr, decomp%local%lo, decomp%local%hi, sendbuf, ranges)

         call MPI_Sendrecv(sendbuf, nxc, MPI_DOUBLE, topo%neighbors(2*id),   1,&
                           recvbuf, nxc, MPI_DOUBLE, topo%neighbors(2*id-1), 1,&
                           topo%comm, MPI_STATUS_IGNORE, ierr)

         ranges(:,1) = decomp%local%mn(:)
         ranges(:,2) = decomp%local%mx(:)
         ranges(id,1) = decomp%local%rx_lolo(id)
         ranges(id,2) = decomp%local%rx_lohi(id)
         call copy_buffer_to_array_6d_real64(recvbuf, arr, decomp%local%lo, decomp%local%hi, ranges)

      end do

   end subroutine


end module sll_m_decomposition
