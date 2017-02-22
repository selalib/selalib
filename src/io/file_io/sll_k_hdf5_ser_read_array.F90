integer(hsize_t) :: array_dims(rank)
integer(hid_t)   :: dataset_id

array_dims = int( shape(array), hsize_t )

! Open existing dataset
call h5dopen_f( handle%file_id, dsetname, dataset_id, error )
SLL_ASSERT(error==0)

! Read dataset into array
call h5dread_f( dataset_id, DATATYPE, array, array_dims, error )
SLL_ASSERT(error==0)

! Close dataset
call h5dclose_f( dataset_id, error )
SLL_ASSERT(error==0)
