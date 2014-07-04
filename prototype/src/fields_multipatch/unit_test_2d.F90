program unit_test_fields_multipatch
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_module_scalar_field_2d_multipatch
  use sll_coordinate_transformation_multipatch_module
  implicit none

  
  type(sll_coordinate_transformation_multipatch_2d), pointer :: T
  class(sll_scalar_field_multipatch_2d), pointer              :: F
  type(sll_logical_mesh_2d), pointer                         :: m
  type(sll_coordinate_transformation_2d_nurbs), pointer      :: transf
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
  sll_real64 :: x1,x2

  T => new_coordinate_transformation_multipatch_2d("square_4p_n10")
  print *, 'initialized multipatch transformation'
  
  
  F => new_scalar_field_multipatch_2d("test_field_multipatch", T)
  print *, 'initialized scalar field multipatch'

  call F%allocate_memory()
  print *, 'allocated memory within the field'

  print *, 'initializing the multipatch field'
  ! loop over patches. This is necessary as in the general case we will not
  ! have a 'global' or 'abstract' 2d array we could use to reason. This is
  ! because the physical domain could be modelled by a collection of patches
  ! connected in some weird way.
  num_patches = F%get_number_patches()

  ! There is a lot of work here to initialize the multipatch. Could we do any
  ! better and abstract this? This could be made into a function that takes
  ! an analytical function initializer... but I guess that with the information
  ! coming from CAID, we could also assume that the initial data also comes
  ! from outside...
  do ipatch=0, num_patches-1
     ! Get rid of this 'fix' whenever it is decided that gfortran 4.6 is not
     ! supported by Selalib anymore
     !     m        => F%get_logical_mesh(ipatch)
     m => F%transf%transfs(ipatch+1)%t%mesh
     transf   => F%get_transformation(ipatch)
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
           x1 = transf%x1(eta1,eta2)
           x2 = transf%x2(eta1,eta2)
           val  = test_function(x1,x2)
           call F%set_value_at_indices( i, j, ipatch, val ) 
        end do
     end do
  
     print *, 'updating multipatch field coefficients in the boundary'
     call set_slope_mp(F,ipatch)

  end do

  print *, 'updating multipatch field interpolation coefficients...'
  call F%update_interpolation_coefficients()

  print *, 'writing to file...'
  call F%write_to_file(0)

  call sll_delete(T) 
  call sll_delete(F)
!  call delete_field_sfmp2d_ptr(F)
  print *, 'PASSED'
  

contains

  function test_function( x, y ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    res = sin(x)*cos(y)
  end function test_function

end program unit_test_fields_multipatch

