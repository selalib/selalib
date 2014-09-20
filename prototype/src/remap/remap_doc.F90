! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/prototype/documentation/build/html/doxygen/html/namespaces.html 
! The following lines will be read by doxygen to generate documentation:


!> @namespace sll_remap 
!> @brief 
!> provides capabilities for global data reconfigurations in a parallel 
!> machine. 
!> @details
!> Suppose that a given dataset is represented by a multi-dimensional
!> array which is distributed on multiple processors. The specifics of this
!> distribution (which portion of the global array is contained in which 
!> processor) are contained in an object called a 'layout'. To reconfigure
!> data, we need initial and final layouts, a remap 'plan' which uses the 
!> given layouts for its initialization and finally an application of the 
!> plan on the data described by the layouts. The data reconfiguration is an
!> out-of-place operation, so the module client is responsible for the
!> allocation of the appropriate arrays. Remap operates on multi-dimensional
!> arrays of several of the basic Fortran types. 
!>
!> @author 
!> Selalib team 
!>
!> <b> Headers file available </b>
!>  - sll_remap.h
!>
!> <b> Modules available </b>
!>  List fortran module available
!>  - sll_remap
!>
!> <b> How to use it </b>
!> - Header file : \code #include 'sll_remap.h' \endcode
!> - Link with   <code>-lsll_remap</code>
!> - Add <code> use sll_remap </code>
!>
!> <b> Examples </b>
!> \code
!!
!!  sll_real64, dimension(:,:),  pointer       :: f_eta1
!!  sll_real64, dimension(:,:),  pointer       :: f_eta2
!!  type(layout_2D), pointer                   :: layout_eta1
!!  type(layout_2D), pointer                   :: layout_eta2
!!  type(remap_plan_2D_real64), pointer        :: eta1_to_eta2 
!!  type(remap_plan_2D_real64), pointer        :: eta2_to_eta1
!!
!!  sll_int32  :: prank, comm
!!  sll_int32  :: loc_sz_i, loc_sz_j
!!  sll_int64  :: psize
!!  sll_int32  :: gi
!!  sll_int32  :: gj
!!  sll_int32  :: global_indices(2)
!!  sll_int32  :: error
!!  sll_real64 :: eta1
!!  sll_real64 :: eta2
!!
!!  call sll_boot_collective()
!!
!!  prank = sll_get_collective_rank(sll_world_collective)
!!  psize = sll_get_collective_size(sll_world_collective)
!!  comm  = sll_world_collective%comm
!!
!!  layout_eta1 => new_layout_2D( sll_world_collective )        
!!
!!  call initialize_layout_with_distributed_2D_array( &
!!             nc_eta1+1, nc_eta2+1, 1,int(psize,4),layout_eta1)
!!
!!  if ( prank == MPI_MASTER ) call sll_view_lims_2D( layout_eta1 )
!!
!!  call compute_local_sizes_2d(layout_eta1,loc_sz_i,loc_sz_j)        
!!  SLL_CLEAR_ALLOCATE(f_eta1(1:loc_sz_i,1:loc_sz_j),error)
!!
!!  layout_eta2 => new_layout_2D( sll_world_collective )
!!
!!  call initialize_layout_with_distributed_2D_array( &
!!              nc_eta1+1, nc_eta2+1, int(psize,4),1,layout_eta2)
!!
!! if ( prank == MPI_MASTER ) call sll_view_lims_2D( layout_eta2 )
!!  call flush(6)
!!
!!  call compute_local_sizes_2d(layout_eta2,loc_sz_i,loc_sz_j)        
!!  SLL_CLEAR_ALLOCATE(f_eta2(1:loc_sz_i,1:loc_sz_j),error)
!!
!!  eta1_to_eta2 => new_remap_plan( layout_eta1, layout_eta2, f_eta1)     
!!  eta2_to_eta1 => new_remap_plan( layout_eta2, layout_eta1, f_eta2)     
!!  
!!  do j=1,loc_sz_j
!!  do i=1,loc_sz_i
!!
!!     global_indices = local_to_global_2D(layout_eta2,(/i,j/)) 
!!     gi = global_indices(1)
!!     gj = global_indices(2)
!!
!!     eta1  = eta1_min+(gi-1)*delta_eta1
!!     eta2  = eta2_min+(gj-1)*delta_eta2
!!
!!     f_eta2(i,j)=exp(-.5*(eta1*eta1+eta2*eta2))
!!
!!  end do
!!  end do
!!
!!  call apply_remap_2D( eta2_to_eta1, f_eta2, f_eta1 )
!!
!!  call apply_remap_2D( eta1_to_eta2, f_eta1, f_eta2 )
!!
!!  call delete_layout_2D(layout_eta1)
!!  call delete_layout_2D(layout_eta2)
!!  SLL_DEALLOCATE_ARRAY(f_eta1, error)
!!  SLL_DEALLOCATE_ARRAY(f_eta2, error)
!!
!!  call sll_halt_collective()
!!
!> \endcode
!>
