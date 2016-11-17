integer(HID_T)   :: plist_id
integer(HID_T)   :: dset_id
integer(HID_T)   :: memspace
integer(HID_T)   :: filespace
integer(HSIZE_T) :: dimsfi(dspace_dims)
integer(HSIZE_T) :: block (dspace_dims)
integer(HSIZE_T) :: count (dspace_dims)
integer(HSIZE_T) :: stride(dspace_dims)
sll_int32        :: rank, i

rank   = dspace_dims
dimsfi = global_size
block  = [(int(size(array,i),HSIZE_T), i=1,rank)]
count  = 1
stride = 1

!
! Create dataspaces ('file space' and 'memory space').
!
call h5screate_simple_f( rank, global_size, filespace, error )
SLL_ASSERT(error==0)

call h5screate_simple_f( rank, block, memspace, error )
SLL_ASSERT(error==0)

!
! Open dataset.
!
call h5pcreate_f( H5P_DATASET_CREATE_F, plist_id, error )
SLL_ASSERT(error==0)

!    call h5pset_chunk_f( plist_id, rank, block, error )
!    SLL_ASSERT(error==0)

call h5dopen_f( file_id, dsetname, dset_id, error )
SLL_ASSERT(error==0)

call h5pclose_f( plist_id, error )
SLL_ASSERT(error==0)

call h5sclose_f( filespace, error )
SLL_ASSERT(error==0)

!
! Select hyperslab in the file.
!
call h5dget_space_f( dset_id, filespace, error )
SLL_ASSERT(error==0)

call h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, &
                            offset, count, error, stride, block )
SLL_ASSERT(error==0)

!
! Initialize data buffer with trivial data.
!
!
! Create property list for collective dataset write
!
call h5pcreate_f( H5P_DATASET_XFER_F, plist_id, error )
SLL_ASSERT(error==0)

call h5pset_dxpl_mpio_f( plist_id, H5FD_MPIO_COLLECTIVE_F, error )
SLL_ASSERT(error==0)

!
! Read the dataset collectively.
!
call h5dread_f( dset_id, DATATYPE, array, dimsfi, error, &
                file_space_id = filespace, mem_space_id = memspace, &
                xfer_prp = plist_id )
SLL_ASSERT(error==0)

!
! Close the property list.
!
call h5pclose_f( plist_id, error )
SLL_ASSERT(error==0)

!
! Close dataspaces.
!
call h5sclose_f( filespace, error )
SLL_ASSERT(error==0)

call h5sclose_f( memspace, error )
SLL_ASSERT(error==0)

!
! Close the dataset.
!
call h5dclose_f( dset_id, error )
SLL_ASSERT(error==0)
