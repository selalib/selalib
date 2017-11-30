integer(hid_t)   ::   loc_id
integer(hid_t)   :: space_id
integer(hid_t)   ::  attr_id
integer(hsize_t) :: dims(1)

dims = int( [1], hsize_t )

! Open existing dataset
call h5dopen_f( handle%file_id, dsetname, loc_id, error )
SLL_ASSERT( error==0 )

! Create dataspace for new scalar attribute
call h5screate_f( H5S_SCALAR_F, space_id, error )
SLL_ASSERT( error==0 )

! Create attribute
call h5acreate_f( loc_id, attrname, DATATYPE, space_id, attr_id, error )
SLL_ASSERT( error==0 )

! Write attribute
call h5awrite_f( attr_id, DATATYPE, attrvalue, dims, error )
SLL_ASSERT( error==0 )

! Close dataspace
call h5sclose_f( space_id, error )
SLL_ASSERT( error==0 )

! Close dataset
call h5dclose_f( loc_id, error )
SLL_ASSERT( error==0 )
