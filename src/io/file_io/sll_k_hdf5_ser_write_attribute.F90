integer(hid_t)   ::   loc_id
integer(hid_t)   :: space_id
integer(hid_t)   ::  attr_id
integer(hsize_t) :: dims(1)
type(h5o_info_t) :: object_info

dims = int([1], hsize_t)

! Inquire about type of object in target path: group or dataset?
call h5oget_info_by_name_f(handle%file_id, objpath, object_info, error)
SLL_ASSERT_ALWAYS(error == 0)

! Open existing group or dataset
select case (object_info%type)
case (H5O_TYPE_GROUP)
   call h5gopen_f(handle%file_id, objpath, loc_id, error)
   SLL_ASSERT_ALWAYS(error == 0)
case (H5O_TYPE_DATASET)
   call h5dopen_f(handle%file_id, objpath, loc_id, error)
   SLL_ASSERT_ALWAYS(error == 0)
end select

! Create dataspace for new scalar attribute
call h5screate_f(H5S_SCALAR_F, space_id, error)
SLL_ASSERT_ALWAYS(error == 0)

! Create attribute
call h5acreate_f(loc_id, attrname, DATATYPE, space_id, attr_id, error)
SLL_ASSERT_ALWAYS(error == 0)

! Write attribute
call h5awrite_f(attr_id, DATATYPE, attrvalue, dims, error)
SLL_ASSERT_ALWAYS(error == 0)

! Close attribute
call h5aclose_f(attr_id, error)
SLL_ASSERT_ALWAYS(error == 0)

! Close dataspace
call h5sclose_f(space_id, error)
SLL_ASSERT_ALWAYS(error == 0)

! Close existing group or dataset
select case (object_info%type)
case (H5O_TYPE_GROUP)
   call h5gclose_f(loc_id, error)
   SLL_ASSERT_ALWAYS(error == 0)
case (H5O_TYPE_DATASET)
   call h5dclose_f(loc_id, error)
   SLL_ASSERT_ALWAYS(error == 0)
end select
