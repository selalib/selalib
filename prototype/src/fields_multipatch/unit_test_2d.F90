program unit_test_fields_multipatch
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_module_scalar_field_2d_multipatch
  use sll_coordinate_transformation_multipatch_module
  implicit none

  
  type(sll_coordinate_transformation_multipatch_2d), pointer :: T
  type(sll_scalar_field_multipatch_2d), pointer              :: F
  type(sll_logical_mesh_2d), pointer                         :: m
  sll_int32  :: ipatch
  sll_int32  :: i
  sll_int32  :: j
  sll_int32  :: num_patches
  sll_int32  :: num_pts1
  sll_int32  :: num_pts2
  sll_real64 :: val
  sll_real64 :: eta1
  sll_real64 :: eta2
  sll_real64 :: delta1
  sll_real64 :: delta2

  T => new_coordinate_transformation_multipatch_2d("identity_mp_info.nml")
  print *, 'initialized multipatch transformation'
  
  F => new_scalar_field_multipatch_2d("test_field_multipatch", T)
  print *, 'initialized scalar field multipatch'

  call F%allocate_memory()
  print *, 'allocated memory within the field'

  print *, 'initializing the multipatch field'
  ! loop over patches. This is necessary as in the general case we will not
  ! have a 'global' 2d array in which we could simply loop.
  num_patches = F%get_number_patches()
  do ipatch=0, num_patches-1
     m => F%get_logical_mesh(ipatch)
     num_pts1 = m%num_cells1+1
     num_pts2 = m%num_cells2+1
     delta1   = m%delta_eta1
     delta2   = m%delta_eta2
     do j=1,num_pts1
        eta2 = (j-1)*delta2
        do i=1,num_pts2
           ! here it is assumed that the eta_min's are = 0. This is supposed
           ! to be the case for NURBS transformations.
           eta1 = (i-1)*delta1
           val  = test_function(eta1,eta2)
           call F%set_value_at_indices( i, j, ipatch, val )     
        end do
     end do
  end do

  call F%update_interpolation_coefficients()
  print*, trim(F%field_name)
  call F%write_to_file(0)

  call delete_stmp2d_ptr(T)
  call sll_delete(F)
  print *, 'PASSED'
  

contains

  function test_function( x, y ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    res = sin(x)*cos(y)
  end function test_function

end program unit_test_fields_multipatch

