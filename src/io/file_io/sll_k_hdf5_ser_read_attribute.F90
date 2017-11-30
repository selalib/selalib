integer(hid_t)   ::   loc_id
integer(hid_t)   ::  attr_id
integer(hsize_t) :: dims(1)

dims = int( [1], hsize_t )

! Set starting location to root group
loc_id = handle%file_id

! If path was provided, open existing dataset
if (len_trim(dsetname) > 0) then
  call h5dopen_f( handle%file_id, dsetname, loc_id, error )
  SLL_ASSERT( error==0 )
end if

! Open existing scalar attribute
call h5aopen_f( loc_id, attrname, attr_id, error )

! Read attribute
call h5aread_f( attr_id, DATATYPE, attrvalue, dims, error )
SLL_ASSERT( error==0 )

! Close dataset if needed
if (len_trim(dsetname) > 0) then
  call h5dclose_f( loc_id, error )
  SLL_ASSERT( error==0 )
end if
