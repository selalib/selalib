integer(hsize_t) :: array_dims(rank)
integer(hid_t)   :: dataset_id
integer(hid_t)   :: dataspace_id

array_dims = int( shape(array), hsize_t )

! Create dataspace
call h5screate_simple_f( rank, array_dims, dataspace_id, error )
SLL_ASSERT(error==0)

! Create dataset
call h5dcreate_f( handle%file_id, dsetname, DATATYPE, dataspace_id, &
                  dataset_id, error )
SLL_ASSERT(error==0)

! Write dataset
call h5dwrite_f( dataset_id, DATATYPE, array, array_dims, error )
SLL_ASSERT(error==0)

! Close dataspace
call h5sclose_f( dataspace_id, error )
SLL_ASSERT(error==0)

! Close dataset
call h5dclose_f( dataset_id, error )
SLL_ASSERT(error==0)
