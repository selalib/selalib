module sll_vlasov4d_poisson

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

 use sll_m_remapper
 use sll_m_interpolators_1d_base
 use sll_m_interpolators_2d_base

 use sll_vlasov4d_base

 implicit none
 private
 public :: initialize, free
 public :: advection_x1, advection_x2, advection_x3, advection_x4

 type, public, extends(vlasov4d_base)       :: vlasov4d_poisson

   class(sll_c_interpolator_1d), pointer :: interp_x1
   class(sll_c_interpolator_1d), pointer :: interp_x2
   class(sll_c_interpolator_1d), pointer :: interp_x3
   class(sll_c_interpolator_1d), pointer :: interp_x4

 end type vlasov4d_poisson

 sll_int32, private :: i, j, k, l
 sll_int32, private :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l
 sll_int32, private :: global_indices(4), gi, gj, gk, gl
 sll_int32, private :: ierr

 interface initialize
   module procedure initialize_vlasov4d_poisson
 end interface initialize

 interface free
   module procedure free_vlasov4d_poisson
 end interface free

contains

 subroutine initialize_vlasov4d_poisson(this,       &
                                        interp_x1,  &
                                        interp_x2,  &
                                        interp_x3,  &
                                        interp_x4,  &
                                        error)

  use sll_m_hdf5_io_serial

  class(vlasov4d_poisson),intent(inout)   :: this
  class(sll_c_interpolator_1d), target :: interp_x1
  class(sll_c_interpolator_1d), target :: interp_x2
  class(sll_c_interpolator_1d), target :: interp_x3
  class(sll_c_interpolator_1d), target :: interp_x4
  sll_int32                               :: error 

  this%interp_x1 => interp_x1
  this%interp_x2 => interp_x2
  this%interp_x3 => interp_x3
  this%interp_x4 => interp_x4

  call initialize_vlasov4d_base(this)

  SLL_CLEAR_ALLOCATE(this%ex(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%ey(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%rho(1:this%np_eta1,1:this%np_eta2),error)

 end subroutine initialize_vlasov4d_poisson

 subroutine free_vlasov4d_poisson(this)

  class(vlasov4d_poisson),intent(inout) :: this

  call sll_o_delete(this%layout_x)
  call sll_o_delete(this%layout_v)
  SLL_DEALLOCATE_ARRAY(this%f, ierr)
  SLL_DEALLOCATE_ARRAY(this%ft, ierr)

 end subroutine free_vlasov4d_poisson

 subroutine advection_x1(this,dt)

  class(vlasov4d_poisson), intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  SLL_ASSERT( .not. this%transposed) 

  call sll_o_compute_local_sizes(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)
  do l=1,loc_sz_l
  do k=1,loc_sz_k
     global_indices = sll_o_local_to_global(this%layout_x,(/1,1,k,l/)) 
     gk = global_indices(3)
     alpha = (this%eta3_min +(gk-1)*this%delta_eta3)*dt
     do j=1,loc_sz_j
           call this%interp_x1%interpolate_array_disp_inplace(loc_sz_i,this%f(:,j,k,l),alpha)
     end do
  end do
  end do

 end subroutine advection_x1

 subroutine advection_x2(this,dt)

  class(vlasov4d_poisson),intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  SLL_ASSERT( .not. this%transposed)

  call sll_o_compute_local_sizes(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)

  do l=1,loc_sz_l

    global_indices = sll_o_local_to_global(this%layout_x,(/1,1,1,l/)) 
    gl = global_indices(4)
    alpha = (this%eta4_min +(gl-1)*this%delta_eta4)*dt

    do k=1,loc_sz_k
    do i=1,loc_sz_i

          call this%interp_x2%interpolate_array_disp_inplace(loc_sz_j,this%f(i,:,k,l),alpha)

    end do
    end do

  end do

 end subroutine advection_x2

 subroutine advection_x3(this,dt)

  class(vlasov4d_poisson), intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64 :: alpha
  SLL_ASSERT(this%transposed) 
  call sll_o_compute_local_sizes(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)

  do l=1,loc_sz_l
  do j=1,loc_sz_j
  do i=1,loc_sz_i

     global_indices = sll_o_local_to_global(this%layout_v,(/i,j,1,l/)) 
     gi = global_indices(1)
     gj = global_indices(2)
     alpha = this%ex(gi,gj)*dt

        call this%interp_x3%interpolate_array_disp_inplace(loc_sz_k,this%ft(i,j,:,l),alpha)

  end do
  end do
  end do

 end subroutine advection_x3

 subroutine advection_x4(this,dt)

  class(vlasov4d_poisson),intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  SLL_ASSERT(this%transposed) 
  call sll_o_compute_local_sizes(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  do k=1,loc_sz_k
  do j=1,loc_sz_j
  do i=1,loc_sz_i

     global_indices = sll_o_local_to_global(this%layout_v,(/i,j,k,1/)) 
     gi = global_indices(1)
     gj = global_indices(2)
     alpha = this%ey(gi,gj)*dt
        call this%interp_x4%interpolate_array_disp_inplace(loc_sz_l,this%ft(i,j,k,:),alpha)

  end do
  end do
  end do

 end subroutine advection_x4

end module sll_vlasov4d_poisson
