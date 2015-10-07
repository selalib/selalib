module sll_vlasov2d_spectral

#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_memory.h"
use sll_m_interpolators_1d_base
use sll_collective
use sll_remapper
#include "sll_fftw.h"

use, intrinsic :: iso_c_binding
use sll_vlasov2d_base
use fftw3

implicit none
private
public :: initialize, free,             &
          spectral_advection_x,         &
          advection_x, advection_v

type, public, extends(vlasov2d_base) :: vlasov2d_spectral

  sll_real64, dimension(:),        pointer :: exn
  sll_real64, dimension(:),        pointer :: jx1,jx2,jx3
  sll_real64, dimension(:),        pointer :: jy1,jy2,jy3
  sll_real64, dimension(:),        pointer :: d_dx
  sll_real64, dimension(:),        pointer :: kx
  fftw_plan                                :: fwx
  fftw_plan                                :: bwx
  fftw_plan                                :: p_tmp_x
  fftw_comp,  dimension(:),        pointer :: tmp_x

  sll_real64, dimension(:,:),      pointer :: f_star
  sll_real64, dimension(:,:),      pointer :: ft_star

  class(sll_interpolator_1d_base), pointer :: interp_x
  class(sll_interpolator_1d_base), pointer :: interp_v

end type vlasov2d_spectral



interface initialize
   module procedure initialize_vlasov2d_spectral
end interface initialize
interface free
   module procedure free_vlasov2d_spectral
end interface free

contains

 subroutine initialize_vlasov2d_spectral(this,      &
                                         interp_x, &
                                         interp_v, &
                                         error)

  class(vlasov2d_spectral),intent(inout)   :: this
  sll_int32                                :: error

  sll_real64  :: kx0
  fftw_int    :: sz_tmp_x
  sll_int32   :: psize, prank, comm
  sll_int32   :: loc_sz_i,loc_sz_j
  sll_int32   :: i

  class(sll_interpolator_1d_base), target :: interp_x
  class(sll_interpolator_1d_base), target :: interp_v

  this%interp_x => interp_x
  this%interp_v => interp_v

  call initialize_vlasov2d_base(this)

  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  SLL_CLEAR_ALLOCATE(this%ex(1:this%np_eta1),error)
  SLL_CLEAR_ALLOCATE(this%rho(1:this%np_eta1),error)
  SLL_CLEAR_ALLOCATE(this%jx(1:this%np_eta1),error)

  SLL_CLEAR_ALLOCATE(this%exn(1:this%np_eta1),error)
  SLL_CLEAR_ALLOCATE(this%jx1(1:this%np_eta1),error)
  SLL_CLEAR_ALLOCATE(this%jx2(1:this%np_eta1),error)
  SLL_CLEAR_ALLOCATE(this%jx3(1:this%np_eta1),error)
  
  FFTW_ALLOCATE(this%tmp_x,this%nc_eta1/2+1,sz_tmp_x,this%p_tmp_x)
  SLL_CLEAR_ALLOCATE(this%d_dx(1:this%nc_eta1),error)

  NEW_FFTW_PLAN_R2C_1D(this%fwx, this%nc_eta1, this%d_dx,  this%tmp_x)
  NEW_FFTW_PLAN_C2R_1D(this%bwx, this%nc_eta1, this%tmp_x, this%d_dx)

  SLL_CLEAR_ALLOCATE(this%kx(1:this%nc_eta1/2+1), error)
   
  kx0 = 2._f64*sll_pi/(this%nc_eta1*this%delta_eta1)

  this%kx(1) = 1.0_f64
  do i=2,this%nc_eta1/2+1
     this%kx(i) = (i-1)*kx0
  end do

  call compute_local_sizes(this%layout_x,loc_sz_i,loc_sz_j)
  SLL_CLEAR_ALLOCATE(this%f_star(1:loc_sz_i,1:loc_sz_j),error)

  SLL_CLEAR_ALLOCATE(this%ex( 1:loc_sz_i),error)
  SLL_CLEAR_ALLOCATE(this%jx( 1:loc_sz_i),error)
  SLL_CLEAR_ALLOCATE(this%rho(1:loc_sz_i),error)

  call compute_local_sizes(this%layout_v,loc_sz_i,loc_sz_j)
  SLL_CLEAR_ALLOCATE(this%ft_star(1:loc_sz_i,1:loc_sz_j),error)


 end subroutine initialize_vlasov2d_spectral

 subroutine free_vlasov2d_spectral(this)

  class(vlasov2d_spectral) :: this
  sll_int32                :: error

  call sll_delete(this%layout_x)
  call sll_delete(this%layout_v)
  SLL_DEALLOCATE_ARRAY(this%f,error) 
  SLL_DEALLOCATE_ARRAY(this%ft, error)

#ifdef FFTW_F2003
  if (c_associated(this%p_tmp_x)) call fftw_free(this%p_tmp_x)
#endif

  call fftw_destroy_plan(this%fwx)
  call fftw_destroy_plan(this%bwx)

 end subroutine free_vlasov2d_spectral

 subroutine spectral_advection_x(this,dt)

  class(vlasov2d_spectral), intent(inout) :: this

  sll_real64, intent(in) :: dt
  sll_real64 :: v, v_min, delta_v
  sll_int32  :: loc_sz_i,loc_sz_j
  sll_int32  :: nc_x
  sll_int32  :: i, j
  sll_int32  :: gj
  sll_int32  :: global_indices(2)

  SLL_ASSERT( .not. this%transposed) 

  nc_x    = this%nc_eta1
  v_min   = this%eta2_min
  delta_v = this%delta_eta2

  call compute_local_sizes(this%layout_x,loc_sz_i,loc_sz_j)
  
  do j=1,loc_sz_j

     global_indices = local_to_global(this%layout_x,(/1,j/)) 
     gj = global_indices(2)
     v  = (v_min +(gj-1)*delta_v)*dt

     this%d_dx = this%f(1:nc_x,j)
     call fftw_execute_dft_r2c(this%fwx, this%d_dx,this%tmp_x)

     !exact : f = f^n exp(-i kx vx dt)

     do i = 2, nc_x/2+1
       this%tmp_x(i) = this%tmp_x(i)*exp(-cmplx(0.0_f64,this%kx(i),kind=f64)*v)
     end do

     call fftw_execute_dft_c2r(this%bwx, this%tmp_x, this%d_dx)
     this%f(1:nc_x,j)= this%d_dx / nc_x

  end do

  this%f(nc_x+1,:) = this%f(1,:)

 end subroutine spectral_advection_x

 subroutine advection_x(this,dt)

  class(vlasov2d_spectral),intent(inout) :: this
  sll_real64, intent(in)                 :: dt
  sll_real64                             :: alpha
  sll_int32                              :: loc_sz_i
  sll_int32                              :: loc_sz_j
  sll_int32                              :: j
  sll_int32                              :: gj
  sll_int32                              :: global_indices(2)

  SLL_ASSERT( .not. this%transposed)

  call compute_local_sizes(this%layout_x,loc_sz_i,loc_sz_j)

  do j=1,loc_sz_j

    global_indices = local_to_global(this%layout_x,(/1,j/)) 
    gj = global_indices(2)
    alpha = (this%eta2_min + (gj-1)*this%delta_eta2)*dt

    this%f(:,j) = this%interp_x%interpolate_array_disp(loc_sz_i,this%f(:,j),alpha)

  end do

 end subroutine advection_x

 subroutine advection_v(this,dt)

  class(vlasov2d_spectral),intent(inout) :: this
  sll_real64, intent(in)                 :: dt
  sll_real64                             :: alpha
  sll_int32                              :: loc_sz_i
  sll_int32                              :: loc_sz_j
  sll_int32                              :: i
  sll_int32                              :: gi
  sll_int32                              :: global_indices(2)

  SLL_ASSERT(this%transposed)

  call compute_local_sizes(this%layout_v,loc_sz_i,loc_sz_j)

  do i=1,loc_sz_i

    global_indices = local_to_global(this%layout_v,(/i,1/)) 
    gi = global_indices(1)
    alpha = this%ex(gi) *dt
    this%ft(i,:) = this%interp_v%interpolate_array_disp(loc_sz_j,this%ft(i,:),alpha)

  end do

 end subroutine advection_v

end module sll_vlasov2d_spectral
