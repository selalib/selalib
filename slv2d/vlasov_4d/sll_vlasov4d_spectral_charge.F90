#define FFTW_ALLOCATE(array,array_size,sz_array,p_array)  \
sz_array = int((array_size/2+1),C_SIZE_T);                \
p_array = fftw_alloc_complex(sz_array);                   \
call c_f_pointer(p_array, array, [array_size/2+1])        \

module sll_vlasov4d_spectral_charge

#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use sll_module_interpolators_1d_base
use sll_module_interpolators_2d_base
use sll_collective
use sll_remapper
use sll_constants


 use, intrinsic :: iso_c_binding
 use sll_vlasov4d_base

 implicit none
 private
 public :: initialize, free, densite_courantx, densite_couranty
 public :: advection_x1, advection_x2, advection_x3x4

 type, public, extends(vlasov4d_base) :: vlasov4d_spectral_charge

   sll_real64, dimension(:,:), pointer               :: exn
   sll_real64, dimension(:,:), pointer               :: eyn
   sll_real64, dimension(:,:), pointer               :: bzn
   sll_real64, dimension(:,:), pointer               :: jx1,jx2,jx3
   sll_real64, dimension(:,:), pointer               :: jy1,jy2,jy3
   sll_real64, dimension(:),   allocatable           :: d_dx
   sll_real64, dimension(:),   allocatable           :: d_dy
   sll_real64, dimension(:),   allocatable           :: kx
   sll_real64, dimension(:),   allocatable           :: ky
   type(C_PTR)                                       :: fwx, fwy
   type(C_PTR)                                       :: bwx, bwy
   type(C_PTR)                                       :: p_tmp_x, p_tmp_y
   complex(C_DOUBLE_COMPLEX), dimension(:),  pointer :: tmp_x, tmp_y
   class(sll_interpolator_2d_base), pointer          :: interp_x3x4
   type(remap_plan_2D_real64), pointer               :: x1_to_x2 
   type(remap_plan_2D_real64), pointer               :: x2_to_x1

   sll_real64, dimension(:,:,:,:),  pointer :: f_star
   sll_real64, dimension(:,:,:,:),  pointer :: ft_star


 end type vlasov4d_spectral_charge

 sll_int32, private :: i, j, k, l
 sll_int32, private :: global_indices(4), gi, gj, gk, gl
 sll_int32, private :: ierr

include 'fftw3.f03'

 interface initialize
    module procedure initialize_vlasov4d_spectral_charge
 end interface initialize
 interface free
    module procedure free_vlasov4d_spectral_charge
 end interface free

contains

 subroutine initialize_vlasov4d_spectral_charge(this,interp_x3x4,error)

  use sll_hdf5_io

  class(vlasov4d_spectral_charge),intent(inout)   :: this
  class(sll_interpolator_2d_base), target :: interp_x3x4
  sll_int32                               :: error

  sll_real64        :: kx0, ky0
  integer(C_SIZE_T) :: sz_tmp_x, sz_tmp_y
  sll_int32         :: psize, prank, comm
  sll_int32         :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l

  this%interp_x3x4 => interp_x3x4

  call initialize_vlasov4d_base(this)

  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  SLL_CLEAR_ALLOCATE(this%ex(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%ey(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%exn(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%eyn(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%bzn(1:this%np_eta1,1:this%np_eta2),error)

  SLL_CLEAR_ALLOCATE(this%bz(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%rho(1:this%np_eta1,1:this%np_eta2),error)

  SLL_CLEAR_ALLOCATE(this%jx(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%jy(1:this%np_eta1,1:this%np_eta2),error)

  SLL_CLEAR_ALLOCATE(this%jx1(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%jx2(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%jx3(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%jy1(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%jy2(1:this%np_eta1,1:this%np_eta2),error)
  SLL_CLEAR_ALLOCATE(this%jy3(1:this%np_eta1,1:this%np_eta2),error)
  
  FFTW_ALLOCATE(this%tmp_x,this%nc_eta1/2+1,sz_tmp_x,this%p_tmp_x)
  FFTW_ALLOCATE(this%tmp_y,this%nc_eta2/2+1,sz_tmp_y,this%p_tmp_y)
  SLL_CLEAR_ALLOCATE(this%d_dx(1:this%nc_eta1),error)
  SLL_CLEAR_ALLOCATE(this%d_dy(1:this%nc_eta2),error)

  this%fwx = fftw_plan_dft_r2c_1d(this%nc_eta1,this%d_dx, this%tmp_x,FFTW_ESTIMATE)
  this%bwx = fftw_plan_dft_c2r_1d(this%nc_eta1,this%tmp_x,this%d_dx, FFTW_ESTIMATE)
  this%fwy = fftw_plan_dft_r2c_1d(this%nc_eta2,this%d_dy, this%tmp_y,FFTW_ESTIMATE)
  this%bwy = fftw_plan_dft_c2r_1d(this%nc_eta2,this%tmp_y,this%d_dy, FFTW_ESTIMATE)

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

  call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)
  SLL_CLEAR_ALLOCATE(this%f_star(1:loc_sz_i,1:loc_sz_j,1:loc_sz_k,1:loc_sz_l),ierr)

  call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)
  SLL_CLEAR_ALLOCATE(this%ft_star(1:loc_sz_i,1:loc_sz_j,1:loc_sz_k,1:loc_sz_l),ierr)

 end subroutine initialize_vlasov4d_spectral_charge

 subroutine free_vlasov4d_spectral_charge(this)

  class(vlasov4d_spectral_charge) :: this

  call delete_layout_4D(this%layout_x)
  call delete_layout_4D(this%layout_v)
  SLL_DEALLOCATE_ARRAY(this%f, ierr)
  SLL_DEALLOCATE_ARRAY(this%ft, ierr)
  if (c_associated(this%p_tmp_x)) call fftw_free(this%p_tmp_x)
  if (c_associated(this%p_tmp_y)) call fftw_free(this%p_tmp_y)

  call fftw_destroy_plan(this%fwx)
  call fftw_destroy_plan(this%fwy)
  call fftw_destroy_plan(this%bwx)
  call fftw_destroy_plan(this%bwy)

 end subroutine free_vlasov4d_spectral_charge

 subroutine advection_x1(this,dt)

  class(vlasov4d_spectral_charge), intent(inout) :: this

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
           call fftw_execute_dft_r2c(this%fwx, this%f(1:nc_x1,j,k,l),this%tmp_x)
           !exact : f* = f^n exp(-i kx vx dt)
           !calcul du flux
           this%tmp_x = this%tmp_x &
                 * (1._f64-exp(-cmplx(0.0_f64,1,kind=f64)*vx*this%kx)) &
                 * cmplx(0.0_f64,-1._f64,kind=f64)/(dt*this%kx)
           call fftw_execute_dft_c2r(this%bwx, this%tmp_x, this%d_dx)
           this%f_star(1:nc_x1,j,k,l)= this%d_dx / nc_x1
        end do
     end do
  end do

  this%f_star(nc_x1+1,:,:,:) = this%f_star(1,:,:,:)

  call apply_remap_4D( this%x_to_v, this%f_star, this%ft_star) 

!calculer le courant avec la formule 
! f^* = f^n *exp(-ik vx dt) = f^n - vx * dt * ik f^n (1-exp(-ik vx dt))/(ik*dt*vx)
! jx^* = int ik f^n (1-exp(-ik vx dt))/(ik*dt) dvxdvy
  call densite_courantx(this, "*")

  do l=1,loc_sz_l
     do k=1,loc_sz_k
        global_indices = local_to_global_4D(this%layout_x,(/1,1,k,l/)) 
        gk = global_indices(3)
        vx = (x3_min +(gk-1)*delta_x3)*dt
        do j=1,loc_sz_j
           call fftw_execute_dft_r2c(this%fwx, this%f(1:nc_x1,j,k,l),this%tmp_x)
           !exact : f* = f^n exp(-i kx vx dt)
           this%tmp_x = this%tmp_x * exp(-cmplx(0.0_f64,this%kx,kind=f64)*vx)
           call fftw_execute_dft_c2r(this%bwx, this%tmp_x, this%d_dx)
           this%f(1:nc_x1,j,k,l)= this%d_dx / nc_x1
        end do
     end do
  end do

  this%f(nc_x1+1,:,:,:) = this%f(1,:,:,:)

 end subroutine advection_x1

 subroutine advection_x2(this,dt)

  class(vlasov4d_spectral_charge),intent(inout) :: this

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
           this%tmp_y = this%tmp_y &
                * (1._f64-exp(-cmplx(0.0_f64,1,kind=f64)*vy*this%ky)) &
                * cmplx(0.0_f64,-1._f64,kind=f64)/(dt*this%ky)
           call fftw_execute_dft_c2r(this%bwy, this%tmp_y, this%d_dy)
           this%f_star(i,1:nc_x2,k,l) = this%d_dy / nc_x2
        end do
     end do
  end do

  this%f_star(:,nc_x2+1,:,:) = this%f_star(:,1,:,:)

  call apply_remap_4D( this%x_to_v, this%f_star, this%ft_star) 
  call densite_couranty(this, "*")

  do l=1,loc_sz_l
     global_indices = local_to_global_4D(this%layout_x,(/1,1,1,l/)) 
     gl = global_indices(4)
     vy = (x4_min +(gl-1)*delta_x4)*dt
     do k=1,loc_sz_k
        do i=1,loc_sz_i
           call fftw_execute_dft_r2c(this%fwy, this%f(i,1:nc_x2,k,l), this%tmp_y)
           this%tmp_y = this%tmp_y * exp(-cmplx(0.0_f64,this%ky,kind=f64)*vy)
           call fftw_execute_dft_c2r(this%bwy, this%tmp_y, this%d_dy)
           this%f(i,1:nc_x2,k,l) = this%d_dy / nc_x2
        end do
     end do
  end do
  
  this%f(:,nc_x2+1,:,:) = this%f(:,1,:,:)

end subroutine advection_x2

subroutine advection_x3x4(this,dt)

  class(vlasov4d_spectral_charge),intent(inout) :: this
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

        this%ft(i,j,:,:) = &
             this%interp_x3x4%interpolate_array_disp(loc_sz_k,loc_sz_l, &
             this%ft(i,j,:,:),alpha_x,alpha_y)
     end do
  end do

 end subroutine advection_x3x4

 subroutine densite_courantx(this,star)

   class(vlasov4d_spectral_charge),intent(inout)  :: this
   character(len=1), intent(in), optional  :: star
   sll_real64, dimension(:,:,:,:), pointer :: df

   sll_int32  :: error
   sll_int32  :: c
   sll_int32  :: comm
   sll_real64 :: dxy
   sll_real64 :: vx 
   sll_real64, dimension(this%np_eta1,this%np_eta2) :: locjx
   sll_int32  :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l

   if( present(star)) then
      df => this%ft_star
   else
      df => this%ft
   end if
   
   dxy = this%delta_eta3*this%delta_eta4
   SLL_ASSERT(this%transposed)

   call compute_local_sizes_4d(this%layout_v, &
                               loc_sz_i,      &
                               loc_sz_j,      &
                               loc_sz_k,      &
                               loc_sz_l)        

   locjx = 0.0_f64
   do l=1,loc_sz_l
      do k=1,loc_sz_k
         do j=1,loc_sz_j
            do i=1,loc_sz_i
               global_indices = local_to_global_4D(this%layout_v,(/i,j,k,l/)) 
               gi = global_indices(1)
               gj = global_indices(2)
               gk = global_indices(3)
               gl = global_indices(4)
               vx = this%eta3_min+(gk-1)*this%delta_eta3
               locjx(gi,gj) = locjx(gi,gj) + dxy*df(i,j,k,l)*vx
            end do
         end do
      end do
   end do

   this%jx1 = 0._f64
   comm = sll_world_collective%comm
   c    = this%np_eta1*this%np_eta2
   
   call mpi_barrier(comm,error)
   call mpi_allreduce(locjx,this%jx1,c,MPI_REAL8,MPI_SUM,comm,error)
   
 end subroutine densite_courantx

 subroutine densite_couranty(this, star)

   class(vlasov4d_spectral_charge),intent(inout) :: this
   character(len=1), intent(in), optional  :: star

   sll_int32  :: error
   sll_real64 :: vy 
   sll_int32  :: c
   sll_int32  :: comm
   sll_real64 :: dxy
   sll_real64, dimension(this%np_eta1,this%np_eta2) :: locjy
   sll_real64, dimension(:,:,:,:), pointer :: df
   sll_int32  :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l

   if( present(star)) then
      df => this%ft_star
   else
      df => this%ft
   end if

   dxy = this%delta_eta3*this%delta_eta4
   SLL_ASSERT(this%transposed)
   call compute_local_sizes_4d(this%layout_v, &
                               loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

   locjy(:,:) = 0.
   do l=1,loc_sz_l
      do k=1,loc_sz_k
         do j=1,loc_sz_j
            do i=1,loc_sz_i
               global_indices = local_to_global_4D(this%layout_v,(/i,j,k,l/)) 
               gi = global_indices(1)
               gj = global_indices(2)
               gk = global_indices(3)
               gl = global_indices(4)
               vy = this%eta4_min+(gl-1)*this%delta_eta4
               locjy(gi,gj) = locjy(gi,gj) + dxy*df(i,j,k,l) *vy
            end do
         end do
      end do
   end do
   
   this%jy1(:,:) = 0._f64
   comm   = sll_world_collective%comm
   call mpi_barrier(comm,error)
   c=this%np_eta1*this%np_eta2
   call mpi_allreduce(locjy,this%jy1,c, MPI_REAL8,MPI_SUM,comm,error)

 end subroutine densite_couranty

end module sll_vlasov4d_spectral_charge
