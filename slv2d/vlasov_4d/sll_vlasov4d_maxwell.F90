module sll_vlasov4d_maxwell

#include "selalib-mpi.h"

 use sll_vlasov4d_base

 implicit none
 private

 public :: initialize, free
 public :: advection_x1, advection_x2
 public :: advection_x3x4

 type, public, extends(vlasov4d_base) :: vlasov4d_maxwell

   class(sll_interpolator_1d_base), pointer :: interp_x1
   class(sll_interpolator_1d_base), pointer :: interp_x2
   class(sll_interpolator_2d_base), pointer :: interp_x3x4

 end type vlasov4d_maxwell

 sll_int32, private :: i, j, k, l
 sll_int32, private :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l
 sll_int32, private :: global_indices(4), gi, gj, gk, gl
 sll_int32, private :: ierr
 
 interface initialize
   module procedure initialize_vlasov4d_maxwell
 end interface initialize

 interface free
   module procedure free_vlasov4d_maxwell
 end interface free

contains

 subroutine initialize_vlasov4d_maxwell(this,        &
                                        interp_x1,   &
                                        interp_x2,   &
                                        interp_x3x4, &
                                        error )

  class(vlasov4d_maxwell),intent(inout)   :: this
  class(sll_interpolator_1d_base), target :: interp_x1
  class(sll_interpolator_1d_base), target :: interp_x2
  class(sll_interpolator_2d_base), target :: interp_x3x4
  sll_int32                               :: error

  call initialize_vlasov4d_base(this)

  this%interp_x1   => interp_x1
  this%interp_x2   => interp_x2
  this%interp_x3x4 => interp_x3x4

  SLL_CLEAR_ALLOCATE(this%ex(1:this%nc_eta1+1,1:this%nc_eta2+1),error)
  SLL_CLEAR_ALLOCATE(this%ey(1:this%nc_eta1+1,1:this%nc_eta2+1),error)
  SLL_CLEAR_ALLOCATE(this%jx(1:this%nc_eta1+1,1:this%nc_eta2+1),error)
  SLL_CLEAR_ALLOCATE(this%jy(1:this%nc_eta1+1,1:this%nc_eta2+1),error)
  SLL_CLEAR_ALLOCATE(this%bz(1:this%nc_eta1+1,1:this%nc_eta2+1),error)

  SLL_CLEAR_ALLOCATE(this%rho(1:this%nc_eta1+1,1:this%nc_eta2+1),error)

 end subroutine initialize_vlasov4d_maxwell

 subroutine free_vlasov4d_maxwell(this)

  class(vlasov4d_maxwell),intent(inout) :: this

  call delete_layout_4D(this%layout_x)
  call delete_layout_4D(this%layout_v)
  SLL_DEALLOCATE_ARRAY(this%f, ierr)
  SLL_DEALLOCATE_ARRAY(this%ft, ierr)

 end subroutine free_vlasov4d_maxwell

 subroutine advection_x1(this,dt)
  class(vlasov4d_maxwell), intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64 :: alpha, x3_min, delta_x3

  ! verifier que la transposition est a jours
  SLL_ASSERT( .not. this%transposed) 

  x3_min   = this%eta3_min
  delta_x3 = this%delta_eta3

  call compute_local_sizes_4d(this%layout_x, &
                              loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  do l=1,loc_sz_l
  do k=1,loc_sz_k
     global_indices = local_to_global_4D(this%layout_x,(/1,1,k,l/)) 
     gk = global_indices(3)
     alpha = (x3_min +(gk-1)*delta_x3)*dt
     do j=1,loc_sz_j
        this%f(:,j,k,l) = this%interp_x1%interpolate_array_disp(loc_sz_i, &
                                                                this%f(:,j,k,l), &
                                                                alpha)
     end do
  end do
  end do

 end subroutine advection_x1

 subroutine advection_x2(this,dt)
  class(vlasov4d_maxwell),intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64 :: x4_min, delta_x4
  sll_real64 :: alpha

  SLL_ASSERT( .not. this%transposed)

  x4_min   = this%eta4_min
  delta_x4 = this%delta_eta4
  call compute_local_sizes_4d(this%layout_x, &
                              loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  do l=1,loc_sz_l

    global_indices = local_to_global_4D(this%layout_x,(/1,1,1,l/)) 
    gl = global_indices(4)
    alpha = (x4_min +(gl-1)*delta_x4)*dt

    do k=1,loc_sz_k
    do i=1,loc_sz_i

       this%f(i,:,k,l) = this%interp_x2%interpolate_array_disp(loc_sz_j, &
                                                               this%f(i,:,k,l), &
                                                               alpha)

    end do
    end do

  end do

 end subroutine advection_x2

 subroutine advection_x3x4(this,dt)

  class(vlasov4d_maxwell),intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64, dimension(this%np_eta3,this%np_eta4) :: alpha_x
  sll_real64, dimension(this%np_eta3,this%np_eta4) :: alpha_y
  sll_real64 :: px, py, ctheta, stheta, depvx, depvy
  sll_real64 :: x3_min, x3_max, x4_min, x4_max
  sll_real64 :: delta_x3, delta_x4

  x3_min   = this%eta3_min
  x3_max   = this%eta3_max
  delta_x3 = this%delta_eta3
  x4_min   = this%eta4_min 
  x4_max   = this%eta4_max
  delta_x4 = this%delta_eta4

  SLL_ASSERT(this%transposed) 
  call compute_local_sizes_4d(this%layout_v, &
                              loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  do i=1,loc_sz_i
  do j=1,loc_sz_j

     do k=1,loc_sz_k
     do l=1,loc_sz_l

        global_indices = local_to_global_4D(this%layout_v,(/i,j,k,l/)) 
        gi = global_indices(1)
        gj = global_indices(2)
        gk = global_indices(3)
        gl = global_indices(4)
        px = x3_min+(gk-1)*delta_x3
        py = x4_min+(gl-1)*delta_x4
        ctheta = cos(this%bz(gi,gj)*dt)
        stheta = sin(this%bz(gi,gj)*dt)
        depvx  = 0.5*dt*this%ex(gi,gj)
        depvy  = 0.5*dt*this%ey(gi,gj)
        alpha_x(k,l) = - (px - (depvx+(px+depvx)*ctheta-(py+depvy)*stheta))
        alpha_y(k,l) = - (py - (depvy+(px+depvx)*stheta+(py+depvy)*ctheta))

     end do
     end do

     this%ft(i,j,:,:) = this%interp_x3x4%interpolate_array_disp(loc_sz_k,loc_sz_l, &
                                                 this%ft(i,j,:,:),alpha_x,alpha_y)
  end do
  end do

 end subroutine advection_x3x4

end module sll_vlasov4d_maxwell
