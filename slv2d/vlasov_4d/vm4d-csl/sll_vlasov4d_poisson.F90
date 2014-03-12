module sll_vlasov4d_poisson

#include "selalib-mpi.h"

 use geometry_module
 use sll_vlasov4d_base

 implicit none
 private
 public :: new, free
 public :: advection_x1, advection_x2, advection_x3, advection_x4

 type, public, extends(vlasov4d_base)       :: vlasov4d_poisson

   class(sll_interpolator_1d_base), pointer :: interp_x1
   class(sll_interpolator_1d_base), pointer :: interp_x2
   class(sll_interpolator_1d_base), pointer :: interp_x3
   class(sll_interpolator_1d_base), pointer :: interp_x4

 end type vlasov4d_poisson

 sll_int32, private :: i, j, k, l
 sll_int32, private :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l
 sll_int32, private :: global_indices(4), gi, gj, gk, gl
 sll_int32, private :: ierr

 interface new
   module procedure new_vlasov4d_poisson
 end interface new

 interface free
   module procedure free_vlasov4d_poisson
 end interface free

contains

 subroutine new_vlasov4d_poisson(this,geomx,geomv,interp_x1,interp_x2,interp_x3,interp_x4,error)

  use sll_hdf5_io
  class(vlasov4d_poisson),intent(inout)   :: this
  type(geometry),intent(in)               :: geomx
  type(geometry),intent(in)               :: geomv
  class(sll_interpolator_1d_base), target :: interp_x1
  class(sll_interpolator_1d_base), target :: interp_x2
  class(sll_interpolator_1d_base), target :: interp_x3
  class(sll_interpolator_1d_base), target :: interp_x4
  sll_int32                               :: error

  this%interp_x1 => interp_x1
  this%interp_x2 => interp_x2
  this%interp_x3 => interp_x3
  this%interp_x4 => interp_x4

  call new_vlasov4d_base(this,geomx,geomv,error)

  SLL_CLEAR_ALLOCATE(this%ex(1:geomx%nx,1:geomx%ny),error)
  SLL_CLEAR_ALLOCATE(this%ey(1:geomx%nx,1:geomx%ny),error)
  SLL_CLEAR_ALLOCATE(this%rho(1:geomx%nx,1:geomx%ny),error)

 end subroutine new_vlasov4d_poisson

 subroutine free_vlasov4d_poisson(this)

  class(vlasov4d_poisson),intent(inout) :: this

  call delete_layout_4D(this%layout_x)
  call delete_layout_4D(this%layout_v)
  SLL_DEALLOCATE_ARRAY(this%f, ierr)
  SLL_DEALLOCATE_ARRAY(this%ft, ierr)

 end subroutine free_vlasov4d_poisson

 subroutine advection_x1(this,dt)

  class(vlasov4d_poisson), intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64 :: alpha, x3_min, delta_x3

  ! verifier que la transposition est a jours
  SLL_ASSERT( .not. this%transposed) 

  x3_min   = this%geomv%x0
  delta_x3 = this%geomv%dx

  call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  do l=1,loc_sz_l
  do k=1,loc_sz_k
     global_indices = local_to_global_4D(this%layout_x,(/1,1,k,l/)) 
     gk = global_indices(3)
     alpha = (x3_min +(gk-1)*delta_x3)*dt
     do j=1,loc_sz_j
        this%f(:,j,k,l) = this%interp_x1%interpolate_array_disp(loc_sz_i,this%f(:,j,k,l),alpha)
     end do
  end do
  end do

 end subroutine advection_x1

 subroutine advection_x2(this,dt)

  class(vlasov4d_poisson),intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64 :: x4_min, delta_x4
  sll_real64 :: alpha

  SLL_ASSERT( .not. this%transposed)

  x4_min   = this%geomv%y0
  delta_x4 = this%geomv%dy
  call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  do l=1,loc_sz_l

    global_indices = local_to_global_4D(this%layout_x,(/1,1,1,l/)) 
    gl = global_indices(4)
    alpha = (x4_min +(gl-1)*delta_x4)*dt

    do k=1,loc_sz_k
    do i=1,loc_sz_i

       this%f(i,:,k,l) = this%interp_x2%interpolate_array_disp(loc_sz_j,this%f(i,:,k,l),alpha)

    end do
    end do

  end do

 end subroutine advection_x2

 subroutine advection_x3(this,dt)

  class(vlasov4d_poisson), intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64 :: alpha
  SLL_ASSERT(this%transposed) 
  call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  do l=1,loc_sz_l
  do j=1,loc_sz_j
  do i=1,loc_sz_i

     global_indices = local_to_global_4D(this%layout_v,(/i,j,1,l/)) 
     gi = global_indices(1)
     gj = global_indices(2)
     alpha = this%ex(gi,gj)*dt

     this%ft(i,j,:,l) = this%interp_x3%interpolate_array_disp(loc_sz_k,this%ft(i,j,:,l),alpha)

  end do
  end do
  end do

 end subroutine advection_x3

 subroutine advection_x4(this,dt)

  class(vlasov4d_poisson),intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  SLL_ASSERT(this%transposed) 
  call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  do k=1,loc_sz_k
  do j=1,loc_sz_j
  do i=1,loc_sz_i

     global_indices = local_to_global_4D(this%layout_v,(/i,j,1,1/)) 
     gi = global_indices(1)
     gj = global_indices(2)
     alpha = this%ey(gi,gj)*dt
     this%ft(i,j,k,:) = this%interp_x4%interpolate_array_disp(loc_sz_l,this%ft(i,j,k,:),alpha)

  end do
  end do
  end do

 end subroutine advection_x4

end module sll_vlasov4d_poisson
