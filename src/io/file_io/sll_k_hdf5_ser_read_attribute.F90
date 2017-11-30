integer(hid_t)   ::   loc_id
integer(hid_t)   ::  attr_id
integer(hsize_t) :: dims(1)

dims = int( [1], hsize_t )

! Open existing dataset
call h5dopen_f( handle%file_id, dsetname, loc_id, error )
SLL_ASSERT( error==0 )

! Open existing scalar attribute
call h5aopen_f( loc_id, attrname, attr_id, error )

! Read attribute
call h5aread_f( attr_id, DATATYPE, attrvalue, dims, error )
SLL_ASSERT( error==0 )

! Close dataset
call h5dclose_f( loc_id, error )
SLL_ASSERT( error==0 )
