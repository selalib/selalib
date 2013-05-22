
module sll_bsl
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_splines.h"

implicit none

type :: mesh_1d
   sll_real64 :: eta_min, eta_max
   sll_real64 :: delta_eta
   sll_real64 :: nc_eta
end type

type :: mesh_2d
   sll_real64 :: eta1_min, eta1_max, eta2_min, eta2_max
   sll_real64 :: delta_eta1, delta_eta2
   sll_int32  :: nc_eta1, nc_eta2
end type

interface operator (*)
   function z(x, y)
      type(mesh_1d), intent(in)  :: x
      type(mesh_1d), intent(in)  :: y
      type(mesh_2d), intent(out) :: z

      z%eta1_min   = x%eta_min
      z%eta1_max   = x%eta_max
      z%eta2_min   = y%eta_min
      z%eta2_max   = y%eta_max
      z%delta_eta1 = x%delta_eta
      z%delta_eta2 = y%delta_eta
      z%nc_eta1 = x%nc_eta
      z%nc_eta2 = y%nc_eta
   end function
end interface


!type bsl_workspace_2d
!   type (sll_spline_1D), pointer :: spl_eta1
!   type (sll_spline_1D), pointer :: spl_eta2
!end type bsl_workspace_2d
!
!type bsl_workspace_4d
!   type (sll_spline_1D), pointer :: spl_eta1
!   type (sll_spline_1D), pointer :: spl_eta2
!   type (sll_spline_1D), pointer :: spl_eta3
!   type (sll_spline_1D), pointer :: spl_eta4
!end type bsl_workspace_4d

interface new_bsl_workspace
   module procedure new_bsl_workspace_1d
   module procedure new_bsl_workspace_2d
   module procedure new_bsl_workspace_4d
end interface

interface bsl_step
   module procedure bsl_step_1d
   module procedure bsl_step_2d
   module procedure bsl_step_4d
end interface

!---------------------------------------

contains

!---------------------------------------
  
function new_bsl_workspace_1d(advfield_1D)

type (bsl_workspace_1d), pointer :: new_bsl_workspace_1d
!type (field_1d_vec1), pointer  :: advfield_1D 
sll_real64, dimension(:) :: advfield_1D

sll_int32  :: ierr
sll_int32  :: nc_eta
sll_real64 :: eta_min
sll_real64 :: eta_max

! allocate pointer
SLL_ALLOCATE(new_bsl_workspace_1d,ierr)

! get dimensions
nc_eta    = advfield_1D%descriptor%nc_eta1
eta_min   = advfield_1D%descriptor%eta1_min
eta_max   = advfield_1D%descriptor%eta1_max

! initialize splines
new_bsl_workspace_1d%spl_eta => new_spline_1D( nc_eta+1,        &
                                               eta_min,         &
                                               eta_max,         &
                                               PERIODIC_SPLINE )

end function new_bsl_workspace_1d

function new_bsl_workspace_2d(dist_func_2D)

type (bsl_workspace_2d), pointer :: new_bsl_workspace_2d
type (sll_distribution_function_2D_t), pointer  :: dist_func_2D 

sll_int32  :: ierr
sll_int32  :: nc_eta1
sll_int32  :: nc_eta2
sll_real64 :: eta1_min
sll_real64 :: eta1_max
sll_real64 :: eta2_min
sll_real64 :: eta2_max
sll_int32  :: boundary1_type
sll_int32  :: boundary2_type

! allocate pointer
SLL_ALLOCATE(new_bsl_workspace_2d,ierr)

! get dimensions
nc_eta1    = get_df_nc_eta1( dist_func_2D ) 
eta1_min   = get_df_eta1_min( dist_func_2D )
eta1_max   = get_df_eta1_max( dist_func_2D )
nc_eta2    = get_df_nc_eta2( dist_func_2D ) 
eta2_min   = get_df_eta2_min( dist_func_2D )
eta2_max   = get_df_eta2_max( dist_func_2D )
boundary1_type = get_df_boundary1_type( dist_func_2D )
boundary2_type = get_df_boundary2_type( dist_func_2D )

! initialize splines
new_bsl_workspace_2d%spl_eta1 => new_spline_1D( nc_eta1+1,        &
                                             eta1_min,         &
                                             eta1_max,         &
                                             PERIODIC_SPLINE )

new_bsl_workspace_2d%spl_eta2 => new_spline_1D( nc_eta2+1,        &
                                             eta2_min,         &
                                             eta2_max,         &
                                             HERMITE_SPLINE )  
end function new_bsl_workspace_2d

function new_bsl_workspace_4d(dist_func_4D)

type (bsl_workspace_4d), pointer :: new_bsl_workspace_4d
type (sll_distribution_function_4D_t), pointer  :: dist_func_4D 

sll_int32  :: ierr
sll_int32  :: nc_eta1, nc_eta2, nc_eta3, nc_eta4
sll_real64 :: eta1_min, eta1_max
sll_real64 :: eta2_min, eta2_max
sll_real64 :: eta3_min, eta3_max
sll_real64 :: eta4_min, eta4_max

! allocate pointer
SLL_ALLOCATE(new_bsl_workspace_4d,ierr)

! get dimensions
nc_eta1  = dist_func_4D%field%descriptor_1%nc_eta1
eta1_min = dist_func_4D%field%descriptor_1%eta1_min
eta1_max = dist_func_4D%field%descriptor_1%eta1_max
nc_eta2  = dist_func_4D%field%descriptor_1%nc_eta2
eta2_min = dist_func_4D%field%descriptor_1%eta2_min
eta2_max = dist_func_4D%field%descriptor_1%eta2_max

! initialize splines
new_bsl_workspace_4d%spl_eta1 => new_spline_1D( nc_eta1+1,     &
                                             eta1_min,         &
                                             eta1_max,         &
                                             PERIODIC_SPLINE )

new_bsl_workspace_4d%spl_eta2 => new_spline_1D( nc_eta2+1,     &
                                             eta2_min,         &
                                             eta2_max,         &
                                             PERIODIC_SPLINE )  

! get dimensions
nc_eta3  = dist_func_4D%field%descriptor_2%nc_eta1
eta3_min = dist_func_4D%field%descriptor_2%eta1_min
eta3_max = dist_func_4D%field%descriptor_2%eta1_max
nc_eta4  = dist_func_4D%field%descriptor_2%nc_eta2
eta4_min = dist_func_4D%field%descriptor_2%eta2_min
eta4_max = dist_func_4D%field%descriptor_2%eta2_max

! initialize splines
new_bsl_workspace_4d%spl_eta3 => new_spline_1D( nc_eta3+1,        &
                                                eta3_min,         &
                                                eta3_max,         &
                                                HERMITE_SPLINE )

new_bsl_workspace_4d%spl_eta4 => new_spline_1D( nc_eta4+1,        &
                                                eta4_min,         &
                                                eta4_max,         &
                                                HERMITE_SPLINE )  
end function new_bsl_workspace_4d

subroutine delete_bsl_workspace(bsl_worksp)
type (bsl_workspace_2d), pointer :: bsl_worksp
sll_int32   :: ierr

if( .not. (associated(bsl_worksp))) then
   write (*,'(a)') 'ERROR: delete_bsl_workspace(), not associated argument.'
   STOP
end if
nullify(bsl_worksp%spl_eta1)
nullify(bsl_worksp%spl_eta2)
SLL_DEALLOCATE(bsl_worksp, ierr)
end subroutine delete_bsl_workspace

subroutine bsl_step_1d( bsl_work,       &
                        field_1D,   &
                        advfield,       &
                        delta_t )

type (bsl_workspace_1d),  pointer :: bsl_work
sll_real64, dimension(:)          :: field_1D
type (field_1d_vec1),     pointer :: advfield
sll_real64  ::  delta_t  
sll_int32   :: boundary_type

sll_real64, dimension(:), pointer  ::  eta_out 
sll_int32  :: i
sll_int32  :: ierr
sll_int32  :: nc_eta
sll_real64 :: delta_eta
sll_real64 :: eta_min
sll_real64 :: eta_max
sll_real64 :: eta

SLL_ASSERT(associated(bsl_work))
SLL_ASSERT(associated(advfield))

! get dimensions
nc_eta        = advfield%descriptor%nc_eta1
delta_eta     = advfield%descriptor%delta_eta1
eta_min       = advfield%descriptor%eta1_min
eta_max       = advfield%descriptor%eta1_max
boundary_type = advfield%descriptor%boundary_type

! allocation
SLL_ALLOCATE(eta_out(nc_eta+1),ierr)
    
call compute_spline_1D_periodic( field_1D, bsl_work%spl_eta )

eta = eta_min  
do i = 1, nc_eta+1
   eta_out(i) = eta + delta_t * advfield%data(i)
   if (eta_out(i) < eta_min) then
      eta_out(i) = eta_out(i) + eta_max - eta_min
   else if (eta_out(i) > eta_max) then
      eta_out(i) = eta_out(i) - eta_max + eta_min
   end if
   eta = eta + delta_eta
end do

call interpolate_array_values( eta_out, field_1D, &
                               nc_eta+1, bsl_work%spl_eta )
end subroutine bsl_step_1d

! subroutine bsl_step_2d
! Advances the distribution function on a time step deltat using a second 
! order time split (Strang splitting)
! conservative semi-Lagrangian scheme
subroutine bsl_step_2d( bsl_work,       &
                        dist_func_2D,   &
                        advfield,       &
                        delta_t )

type (bsl_workspace_2d), pointer :: bsl_work
type (sll_distribution_function_2D_t), pointer  :: dist_func_2D
type (field_2D_vec1), pointer  :: advfield
sll_real64  ::  delta_t  

sll_real64, dimension(:), pointer  ::  eta1_out 
sll_real64, dimension(:), pointer  ::  eta2_out
sll_int32  :: i1
sll_int32  :: i2
sll_int32  :: ierr
sll_int32  :: nc_eta1
sll_int32  :: nc_eta2
sll_real64 :: delta_eta1
sll_real64 :: delta_eta2
sll_real64 :: eta1_min
sll_real64 :: eta2_min
sll_real64 :: eta1_max
sll_real64 :: eta2_max
sll_int32  :: boundary1_type
sll_int32  :: boundary2_type
sll_real64 :: val
sll_real64 :: eta1
sll_real64 :: eta2

SLL_ASSERT(associated(bsl_work))
SLL_ASSERT(associated(dist_func_2D))
SLL_ASSERT(associated(advfield))

! get dimensions
nc_eta1        = get_df_nc_eta1(    dist_func_2D ) 
delta_eta1     = get_df_delta_eta1( dist_func_2D )
eta1_min       = get_df_eta1_min(   dist_func_2D )
eta1_max       = get_df_eta1_max(   dist_func_2D )
nc_eta2        = get_df_nc_eta2(    dist_func_2D ) 
delta_eta2     = get_df_delta_eta2( dist_func_2D )
eta2_min       = get_df_eta2_min(   dist_func_2D )
eta2_max       = get_df_eta2_max(   dist_func_2D )
boundary1_type = get_df_boundary1_type( dist_func_2D )
boundary2_type = get_df_boundary2_type( dist_func_2D )

! allocation
SLL_ALLOCATE(eta1_out(nc_eta1+1),ierr)
SLL_ALLOCATE(eta2_out(nc_eta2+1),ierr)
    
eta2 = eta2_min 
do i2 = 1, nc_eta2
get_filename_component(target ${_file} NAME_WE)
   call compute_spline_1D_periodic( dist_func_2D%field%data(:,i2), bsl_work%spl_eta1 )

   eta1 = eta1_min  
   do i1 = 1, nc_eta1+1
      eta1_out(i1) = eta1 + delta_t * advfield%data(i1,i2)
      if (eta1_out(i1) < eta1_min) then
         eta1_out(i1) = eta1_out(i1) + eta1_max - eta1_min
      else if (eta1_out(i1) > eta1_max) then
         eta1_out(i1) = eta1_out(i1) - eta1_max + eta1_min
      end if
      eta1 = eta1 + delta_eta1
   end do

   call interpolate_array_values( eta1_out, dist_func_2D%field%data(:,i2), &
                                  nc_eta1+1, bsl_work%spl_eta1 )
   eta2 = eta2 + delta_eta2
   
end do


eta1 = eta1_min 
do i1 = 1, nc_eta1

   call compute_spline_1D_periodic( dist_func_2D%field%data(i1,:), bsl_work%spl_eta2 )

   eta2 = eta2_min 
   do i2 = 1, nc_eta2+1
      eta2_out(i2) = eta2 +  delta_t * advfield%data(i1,i2)
      if (eta2_out(i2) < eta2_min) then
         eta2_out(i2) = eta2_out(i2) + eta2_max - eta2_min
      else if (eta2_out(i2) > eta2_max) then
         eta2_out(i2) = eta2_out(i2) - eta2_max + eta2_min
      end if
      eta2 = eta2 + delta_eta2
   end do

   call interpolate_array_values( eta2_out, dist_func_2D%field%data(i1,:), &
                                  nc_eta2+1, bsl_work%spl_eta2 )
   eta1 = eta1 + delta_eta1
   
end do

   
end subroutine bsl_step_2d

!> Advances the distribution function on a time step deltat using a 
!> conservative semi-Lagrangian scheme
subroutine bsl_step_4d( bsl_work,       &
                        dist_func_4D,   &
                        advfield,       &
                        deltat )

type (bsl_workspace_4d), pointer :: bsl_work
type (sll_distribution_function_4D_t), pointer  :: dist_func_4D
type (field_2D_vec2), pointer  :: advfield
sll_real64  ::  deltat 

sll_int32  :: order 
sll_real64, dimension(:), pointer  ::  eta1_out 
sll_real64, dimension(:), pointer  ::  eta2_out

sll_int32  :: i1, i2, i3, i4
sll_int32  :: ierr
sll_int32  :: ihalf
sll_int32  :: nc_eta1, nc_eta2
sll_real64 :: delta_eta1, delta_eta2
sll_real64 :: eta1_min, eta2_min, eta1_max, eta2_max
sll_int32  :: boundary1_type, boundary2_type
sll_real64 :: val
sll_real64 :: etget_filename_component(target ${_file} NAME_WE)a1, eta2

! order of scheme
order = 2
! parameter checking
!SLL_ASSERT(associated(bsl_work))
!SLL_ASSERT(associated(dist_func_4D))
!SLL_ASSERT(associated(advfield_old))
!SLL_ASSERT(associated(advfield_new))

! get physical space dimensions
nc_eta1    = dist_func_4D%field%descriptor_1%nc_eta1
delta_eta1 = dist_func_4D%field%descriptor_1%delta_eta1
eta1_min   = dist_func_4D%field%descriptor_1%eta1_min
eta1_max   = dist_func_4D%field%descriptor_1%eta1_max
nc_eta2    = dist_func_4D%field%descriptor_1%nc_eta2 
delta_eta2 = dist_func_4D%field%descriptor_1%delta_eta2
eta2_min   = dist_func_4D%field%descriptor_1%eta2_min
eta2_max   = dist_func_4D%field%descriptor_1%eta2_max

! get physical space dimensions

boundary1_type = dist_func_4D%field%descriptor_1%boundary1_type
boundary2_type = dist_func_4D%field%descriptor_1%boundary2_type

end subroutine bsl_step_4d

end module sll_bsl
