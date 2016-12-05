integer(hid_t)   :: plist_id
integer(hid_t)   :: dset_id
integer(hid_t)   :: memspace
integer(hid_t)   :: filespace
integer(hsize_t) :: dimsfi(dspace_dims)
integer(hsize_t) :: block (dspace_dims)
integer(hsize_t) :: start (dspace_dims)
integer(hsize_t) :: count (dspace_dims)
integer(hsize_t) :: stride(dspace_dims)
sll_int32        :: rank

! Basic dataset parameters
rank   = dspace_dims
count  = 1
stride = 1

! Integer kind type conversions
dimsfi = int(  global_size, hsize_t )
block  = int( shape(array), hsize_t )
start  = int(       offset, hsize_t )

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
call h5pcreate_f( h5p_dataset_create_f, plist_id, error )
SLL_ASSERT(error==0)

!    call h5pset_chunk_f( plist_id, rank, block, error )
!    SLL_ASSERT(error==0)

call h5dopen_f( handle%file_id, dsetname, dset_id, error )
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

call h5sselect_hyperslab_f( filespace, h5s_select_set_f, &
                            start, count, error, stride, block )
SLL_ASSERT(error==0)

!
! Initialize data buffer with trivial data.
!
!
! Create property list for collective dataset write
!
call h5pcreate_f( h5p_dataset_xfer_f, plist_id, error )
SLL_ASSERT(error==0)

call h5pset_dxpl_mpio_f( plist_id, h5fd_mpio_collective_f, error )
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
