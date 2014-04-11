#include "sll_fftw.h"

module sll_vlasov4d_spectral

#include "selalib-mpi.h"

 use sll_vlasov4d_base

 implicit none
 private
 public :: initialize, free
 public :: advection_x1, advection_x2, advection_x3x4

 type, public, extends(vlasov4d_base) :: vlasov4d_spectral

   sll_real64, dimension(:),   allocatable  :: d_dx
   sll_real64, dimension(:),   allocatable  :: d_dy
   sll_real64, dimension(:),   allocatable  :: kx
   sll_real64, dimension(:),   allocatable  :: ky
   fftw_plan                                :: fwx
   fftw_plan                                :: fwy
   fftw_plan                                :: bwx
   fftw_plan                                :: bwy
   fftw_plan                                :: p_tmp_x
   fftw_plan                                :: p_tmp_y
   fftw_comp, dimension(:),  pointer        :: tmp_x
   fftw_comp, dimension(:),  pointer        :: tmp_y
   class(sll_interpolator_2d_base), pointer :: interp_x3x4

 end type vlasov4d_spectral

 sll_int32, private :: i, j, k, l
 sll_int32, private :: global_indices(4), gi, gj, gk, gl
 sll_int32, private :: ierr

 interface initialize
    module procedure initialize_vlasov4d_spectral
 end interface initialize

 interface free
    module procedure free_vlasov4d_spectral
 end interface free

contains

 subroutine initialize_vlasov4d_spectral(this,interp_x3x4,error)

  use sll_hdf5_io

  class(vlasov4d_spectral),intent(inout)  :: this
  class(sll_interpolator_2d_base), target :: interp_x3x4
  sll_int32                               :: error

  sll_real64 :: kx0, ky0
  fftw_int   :: sz_tmp_x, sz_tmp_y
  sll_int32  :: psize, prank, comm

  this%interp_x3x4 => interp_x3x4

  call initialize_vlasov4d_base(this)

  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  SLL_CLEAR_ALLOCATE(this%ex(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%ey(1:this%np_eta1,1:this%np_eta2),error)

  SLL_CLEAR_ALLOCATE(this%bz(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%rho(1:this%np_eta1,1:this%np_eta2),error)

  SLL_CLEAR_ALLOCATE(this%jx(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%jy(1:this%np_eta1,1:this%np_eta2),error)

  FFTW_ALLOCATE(this%tmp_x,this%nc_eta1/2+1,sz_tmp_x,this%p_tmp_x)
  FFTW_ALLOCATE(this%tmp_y,this%nc_eta2/2+1,sz_tmp_y,this%p_tmp_y)
  SLL_CLEAR_ALLOCATE(this%d_dx(1:this%nc_eta1),error)
  SLL_CLEAR_ALLOCATE(this%d_dy(1:this%nc_eta2),error)

  NEW_FFTW_PLAN_R2C_1D(this%fwx,this%nc_eta1,this%d_dx,this%tmp_x) 
  NEW_FFTW_PLAN_C2R_1D(this%bwx,this%nc_eta1,this%tmp_x,this%d_dx) 
  NEW_FFTW_PLAN_R2C_1D(this%fwy,this%nc_eta2,this%d_dy,this%tmp_y) 
  NEW_FFTW_PLAN_C2R_1D(this%bwy,this%nc_eta2,this%tmp_y,this%d_dy) 

  SLL_CLEAR_ALLOCATE(this%kx(1:this%nc_eta1/2+1), error)
  SLL_CLEAR_ALLOCATE(this%ky(1:this%nc_eta2/2+1), error)
   
  kx0 = 2._f64*sll_pi/(this%nc_eta1*this%delta_eta1)
  ky0 = 2._f64*sll_pi/(this%nc_eta2*this%delta_eta2)

  do i=1,this%nc_eta1/2+1
     this%kx(i) = (i-1)*kx0
  end do
  this%kx(1) = 1.0_f64
  do j=1,this%nc_eta2/2+1
     this%ky(j) = (j-1)*ky0
  end do
  this%ky(1) = 1.0_f64

 end subroutine initialize_vlasov4d_spectral

 subroutine free_vlasov4d_spectral(this)

  class(vlasov4d_spectral) :: this

  call delete_layout_4D(this%layout_x)
  call delete_layout_4D(this%layout_v)
  SLL_DEALLOCATE_ARRAY(this%f, ierr)
  SLL_DEALLOCATE_ARRAY(this%ft, ierr)

#ifdef FFTW_F2003
  if (c_associated(this%p_tmp_x)) call fftw_free(this%p_tmp_x)
  if (c_associated(this%p_tmp_y)) call fftw_free(this%p_tmp_y)
#endif

  call fftw_destroy_plan(this%fwx)
  call fftw_destroy_plan(this%fwy)
  call fftw_destroy_plan(this%bwx)
  call fftw_destroy_plan(this%bwy)

 end subroutine free_vlasov4d_spectral

 subroutine advection_x1(this,dt)

  class(vlasov4d_spectral), intent(inout) :: this

  sll_real64, intent(in) :: dt
  sll_real64 :: vx, x3_min, delta_x3
  sll_int32  :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l
  sll_int32  :: nc_x1

  ! verifier que la transposition est a jours
  SLL_ASSERT( .not. this%transposed) 

  nc_x1    = this%nc_eta1
  x3_min   = this%eta3_min
  delta_x3 = this%delta_eta3

  call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  do l=1,loc_sz_l
  do k=1,loc_sz_k
     global_indices = local_to_global_4D(this%layout_x,(/1,1,k,l/)) 
     gk = global_indices(3)
     vx = (x3_min +(gk-1)*delta_x3)*dt
     do j=1,loc_sz_j
        print*, j
        call fftw_execute_dft_r2c(this%fwx, this%f(1:nc_x1,j,k,l),this%tmp_x)
        this%tmp_x = this%tmp_x*exp(-cmplx(0.0_f64,this%kx,kind=f64)*vx)
        call fftw_execute_dft_c2r(this%bwx, this%tmp_x, this%d_dx)
        this%f(1:nc_x1,j,k,l)= this%d_dx / nc_x1
     end do
  end do
  end do

  this%f(nc_x1+1,:,:,:) = this%f(1,:,:,:)

!Euler explicite
!        this%tmp_x = this%tmp_x*(1._f64-cmplx(0.0_f64,this%kx,kind=f64)*vx)
!Euler implicite
!        this%tmp_x = this%tmp_x/(1._f64+cmplx(0.0_f64,this%kx,kind=f64)*vx)
!crank-nicolson
!        this%tmp_x = this%tmp_x/(1._f64+cmplx(0.0_f64,this%kx,kind=f64)*vx*0.5_f64)
!        this%tmp_x = this%tmp_x*(1._f64-cmplx(0.0_f64,this%kx,kind=f64)*vx*0.5_f64)
!Euler cn modified
!        this%tmp_x=this%tmp_x*(1._f64-0.5*vx*cmplx(0.0_f64,this%kx,kind=f64))
 end subroutine advection_x1

 subroutine advection_x2(this,dt)

  class(vlasov4d_spectral),intent(inout) :: this

  sll_real64, intent(in) :: dt
  sll_real64 :: x4_min, delta_x4
  sll_real64 :: vy
  sll_int32  :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l
  sll_int32  :: nc_x2

  SLL_ASSERT( .not. this%transposed)

  nc_x2    = this%nc_eta2
  x4_min   = this%eta4_min
  delta_x4 = this%delta_eta4
  call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)

  do l=1,loc_sz_l

    global_indices = local_to_global_4D(this%layout_x,(/1,1,1,l/)) 
    gl = global_indices(4)
    vy = (x4_min +(gl-1)*delta_x4)*dt

    do k=1,loc_sz_k
    do i=1,loc_sz_i
       call fftw_execute_dft_r2c(this%fwy, this%f(i,1:nc_x2,k,l), this%tmp_y)
       this%tmp_y = this%tmp_y*exp(-cmplx(0.0_f64,this%ky,kind=f64)*vy)
       call fftw_execute_dft_c2r(this%bwy, this%tmp_y, this%d_dy)
       this%f(i,1:nc_x2,k,l) = this%d_dy / nc_x2
    end do
    end do

  end do

  this%f(:,nc_x2+1,:,:) = this%f(:,1,:,:)

!euler explicite
!       this%tmp_y = this%tmp_y*(1._f64-cmplx(0.0_f64,this%ky,kind=f64)*vy)
!euler implicite
!       this%tmp_y = this%tmp_y/(1._f64+cmplx(0.0_f64,this%ky,kind=f64)*vy)
!crank-nicolson
!       this%tmp_y = this%tmp_y/(1._f64+cmplx(0.0_f64,this%ky,kind=f64)*vy*0.5_f64)
!       this%tmp_y = this%tmp_y*(1._f64-cmplx(0.0_f64,this%ky,kind=f64)*vy*0.5_f64)
!Euler cn modified
!       this%tmp_y = this%tmp_y*(1._f64-cmplx(0.0_f64,this%ky,kind=f64)*vy-0.5_f64*(this%ky*vy)**2)

 end subroutine advection_x2

 subroutine advection_x3x4(this,dt)

  class(vlasov4d_spectral),intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64, dimension(this%np_eta3,this%np_eta4) :: alpha_x
  sll_real64, dimension(this%np_eta3,this%np_eta4) :: alpha_y
  sll_real64 :: px, py, ctheta, stheta, depvx, depvy
  sll_real64 :: x3_min, x3_max, x4_min, x4_max
  sll_real64 :: delta_x3, delta_x4
  sll_int32  :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l

  x3_min   = this%eta3_min
  x3_max   = this%eta3_max
  delta_x3 = this%delta_eta3
  x4_min   = this%eta4_min 
  x4_max   = this%eta4_max
  delta_x4 = this%delta_eta4

  SLL_ASSERT(this%transposed) 
  call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)

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


end module sll_vlasov4d_spectral
