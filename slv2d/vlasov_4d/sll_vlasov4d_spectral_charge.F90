#define FFTW_ALLOCATE(array,array_size,sz_array,p_array)  \
sz_array = int((array_size/2+1),C_SIZE_T);                \
p_array = fftw_alloc_complex(sz_array);                   \
call c_f_pointer(p_array, array, [array_size/2+1])        \

module sll_vlasov4d_spectral

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
 use used_precision
 use geometry_module
 use diagnostiques_module
 use sll_vlasov4d_base

 implicit none
 private
 public :: new, free, densite_courantx, densite_couranty
 public :: advection_x1, advection_x2, advection_x3x4

 type, public, extends(vlasov4d_base) :: vlasov4d_spectral

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
   type(layout_2D), pointer                          :: layout_x1
   type(layout_2D), pointer                          :: layout_x2
   type(remap_plan_2D_real64), pointer               :: x1_to_x2 
   type(remap_plan_2D_real64), pointer               :: x2_to_x1

   sll_real64, dimension(:,:,:,:),  pointer :: f_star
   sll_real64, dimension(:,:,:,:),  pointer :: ft_star


 end type vlasov4d_spectral

 sll_int32, private :: i, j, k, l
 sll_int32, private :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l
 sll_int32, private :: global_indices(4), gi, gj, gk, gl
 sll_int32, private :: ierr

include 'fftw3.f03'

 interface new
    module procedure new_vlasov4d_spectral
 end interface new
 interface free
    module procedure free_vlasov4d_spectral
 end interface free

contains

 subroutine new_vlasov4d_spectral(this,geomx,geomv,interp_x3x4,error)

  use sll_hdf5_io

  class(vlasov4d_spectral),intent(inout)   :: this
  type(geometry),intent(in)               :: geomx
  type(geometry),intent(in)               :: geomv
  class(sll_interpolator_2d_base), target :: interp_x3x4
  sll_int32                               :: error

  sll_int32         :: nc_x1, nc_x2, nc_x3, nc_x4
  sll_real64        :: dx, dy, kx0, ky0
  integer(C_SIZE_T) :: sz_tmp_x, sz_tmp_y
  sll_int32         :: psize, prank, comm

  this%interp_x3x4 => interp_x3x4

  call new_vlasov4d_base(this,geomx,geomv,error)

  nc_x1 = this%geomx%nx
  nc_x2 = this%geomx%ny
  nc_x3 = this%geomv%nx
  nc_x4 = this%geomv%ny

  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  this%layout_x1 => new_layout_2D( sll_world_collective )        
  call initialize_layout_with_distributed_2D_array( &
             geomx%nx,geomx%ny,1,int(psize,4),this%layout_x1)

!  call compute_local_sizes_2d(this%layout_x1,loc_sz_i,loc_sz_j)        
!  SLL_CLEAR_ALLOCATE(this%jx1(1:loc_sz_i,1:loc_sz_j),error)

  this%layout_x2 => new_layout_2D( sll_world_collective )
  call initialize_layout_with_distributed_2D_array( &
              geomx%nx,geomx%ny,int(psize,4),1,this%layout_x2)

!  call compute_local_sizes_2d(this%layout_x2,loc_sz_i,loc_sz_j)        
!  SLL_CLEAR_ALLOCATE(this%jx2(1:loc_sz_i,1:loc_sz_j),ierr)

!  this%x1_to_x2 => new_remap_plan( this%layout_x1, this%layout_x2, this%jx1)     
!  this%x2_to_x1 => new_remap_plan( this%layout_x2, this%layout_x1, this%jx2)     
  
  if(prank == MPI_MASTER) then

     print *,'Printing layout x1: '
     call sll_view_lims_2D(this%layout_x1)
     print *,'Printing layout x2: '
     call sll_view_lims_2D(this%layout_x2)

  end if

  SLL_CLEAR_ALLOCATE(this%ex(1:nc_x1,1:nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%ey(1:nc_x1,1:nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%exn(1:nc_x1,1:nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%eyn(1:nc_x1,1:nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%bzn(1:nc_x1,1:nc_x2),error)

  SLL_CLEAR_ALLOCATE(this%bz(1:nc_x1,1:nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%rho(1:nc_x1,1:nc_x2),error)

  SLL_CLEAR_ALLOCATE(this%jx(1:nc_x1,1:nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jy(1:nc_x1,1:nc_x2),error)

  SLL_CLEAR_ALLOCATE(this%jx1(1:nc_x1,1:nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jx2(1:nc_x1,1:nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jx3(1:nc_x1,1:nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jy1(1:nc_x1,1:nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jy2(1:nc_x1,1:nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jy3(1:nc_x1,1:nc_x2),error)
  
  FFTW_ALLOCATE(this%tmp_x,nc_x1/2+1,sz_tmp_x,this%p_tmp_x)
  FFTW_ALLOCATE(this%tmp_y,nc_x2/2+1,sz_tmp_y,this%p_tmp_y)
  SLL_CLEAR_ALLOCATE(this%d_dx(1:nc_x1),error)
  SLL_CLEAR_ALLOCATE(this%d_dy(1:nc_x2),error)

  this%fwx = fftw_plan_dft_r2c_1d(nc_x1,this%d_dx, this%tmp_x,FFTW_ESTIMATE)
  this%bwx = fftw_plan_dft_c2r_1d(nc_x1,this%tmp_x,this%d_dx, FFTW_ESTIMATE)
  this%fwy = fftw_plan_dft_r2c_1d(nc_x2,this%d_dy, this%tmp_y,FFTW_ESTIMATE)
  this%bwy = fftw_plan_dft_c2r_1d(nc_x2,this%tmp_y,this%d_dy, FFTW_ESTIMATE)

  SLL_CLEAR_ALLOCATE(this%kx(1:nc_x1/2+1), error)
  SLL_CLEAR_ALLOCATE(this%ky(1:nc_x2/2+1), error)
   
  dx = geomx%dx
  dy = geomx%dy

  kx0 = 2._f64*sll_pi/(nc_x1*dx)
  ky0 = 2._f64*sll_pi/(nc_x2*dy)

  do i=1,nc_x1/2+1
     this%kx(i) = (i-1)*kx0
  end do
  this%kx(1) = 1.0_f64
  do j=1,nc_x2/2+1
     this%ky(j) = (j-1)*ky0
  end do
  this%ky(1) = 1.0_f64

  call compute_local_sizes_4d(this%layout_x, &
                              loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(this%f_star(1:loc_sz_i,1:loc_sz_j,1:loc_sz_k,1:loc_sz_l),ierr)

  call compute_local_sizes_4d(this%layout_v, &
                              loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(this%ft_star(1:loc_sz_i,1:loc_sz_j,1:loc_sz_k,1:loc_sz_l),ierr)

 end subroutine new_vlasov4d_spectral

 subroutine free_vlasov4d_spectral(this)

  class(vlasov4d_spectral) :: this

  call delete_layout_4D(this%layout_x)
  call delete_layout_4D(this%layout_v)
  SLL_DEALLOCATE_ARRAY(this%f, ierr)
  SLL_DEALLOCATE_ARRAY(this%ft, ierr)
  if (c_associated(this%p_tmp_x)) call fftw_free(this%p_tmp_x)
  if (c_associated(this%p_tmp_y)) call fftw_free(this%p_tmp_y)
  call dfftw_destroy_plan(this%fwx)
  call dfftw_destroy_plan(this%fwy)
  call dfftw_destroy_plan(this%bwx)
  call dfftw_destroy_plan(this%bwy)

 end subroutine free_vlasov4d_spectral

 subroutine advection_x1(this,dt)

  class(vlasov4d_spectral), intent(inout) :: this

  sll_real64, intent(in) :: dt
  sll_real64 :: vx, x3_min, delta_x3

  ! verifier que la transposition est a jours
  SLL_ASSERT( .not. this%transposed) 

  x3_min   = this%geomv%x0
  delta_x3 = this%geomv%dx

  call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  
  do l=1,loc_sz_l
     do k=1,loc_sz_k
        global_indices = local_to_global_4D(this%layout_x,(/1,1,k,l/)) 
        gk = global_indices(3)
        vx = (x3_min +(gk-1)*delta_x3)*dt
        do j=1,loc_sz_j
           call fftw_execute_dft_r2c(this%fwx, this%f(:,j,k,l),this%tmp_x)
           !exact : f* = f^n exp(-i kx vx dt)
           !calcul du flux
           do i=2,this%geomx%nx
              this%tmp_x(i) = this%tmp_x(i)*(1._f64-exp(-cmplx(0.0_f64,1,kind=f64)*vx*this%kx(i)))*cmplx(0.0_f64,-1._f64,kind=f64)/(dt*this%kx(i))
           enddo
           this%tmp_x(1)=0._f64
           call fftw_execute_dft_c2r(this%bwx, this%tmp_x, this%d_dx)
           this%f_star(:,j,k,l)= this%d_dx / loc_sz_i
        end do
     end do
  end do

!calculer le courant avec la formule f^* = f^n *exp(-ik vx dt) = f^n - vx * dt * ik f^n (1-exp(-ik vx dt))/(ik*dt*vx)
!jx^* = int ik f^n (1-exp(-ik vx dt))/(ik*dt) dvxdvy
  call densite_courantx(this, "*")

  do l=1,loc_sz_l
     do k=1,loc_sz_k
        global_indices = local_to_global_4D(this%layout_x,(/1,1,k,l/)) 
        gk = global_indices(3)
        vx = (x3_min +(gk-1)*delta_x3)*dt
        do j=1,loc_sz_j
           call fftw_execute_dft_r2c(this%fwx, this%f(:,j,k,l),this%tmp_x)
           !exact : f* = f^n exp(-i kx vx dt)
           do i=2,this%geomx%nx
              this%tmp_x(i) = this%tmp_x(i)*exp(-cmplx(0.0_f64,this%kx(i),kind=f64)*vx)
           enddo
           call fftw_execute_dft_c2r(this%bwx, this%tmp_x, this%d_dx)
           this%f(:,j,k,l)= this%d_dx / loc_sz_i
        end do
     end do
  end do

 end subroutine advection_x1

 subroutine advection_x2(this,dt)

  class(vlasov4d_spectral),intent(inout) :: this

  sll_real64, intent(in) :: dt
  sll_real64 :: x4_min, delta_x4
  sll_real64 :: vy

  SLL_ASSERT( .not. this%transposed)

  x4_min   = this%geomv%y0
  delta_x4 = this%geomv%dy
  call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  do l=1,loc_sz_l
     global_indices = local_to_global_4D(this%layout_x,(/1,1,1,l/)) 
     gl = global_indices(4)
     vy = (x4_min +(gl-1)*delta_x4)*dt
     do k=1,loc_sz_k
        do i=1,loc_sz_i
           call fftw_execute_dft_r2c(this%fwy, this%f(i,:,k,l), this%tmp_y)
           do j=2,this%geomx%ny
              this%tmp_y(j) = this%tmp_y(j)*(1._f64-exp(-cmplx(0.0_f64,1,kind=f64)*vy*this%ky(j)))*cmplx(0.0_f64,-1._f64,kind=f64)/(dt*this%ky(j))
           enddo
           this%tmp_y(1)=0._f64
           call fftw_execute_dft_c2r(this%bwy, this%tmp_y, this%d_dy)
           this%f_star(i,:,k,l) = this%d_dy / loc_sz_j
        end do
     end do
  end do

  call densite_couranty(this, "*")

  
  do l=1,loc_sz_l
     global_indices = local_to_global_4D(this%layout_x,(/1,1,1,l/)) 
     gl = global_indices(4)
     vy = (x4_min +(gl-1)*delta_x4)*dt
     do k=1,loc_sz_k
        do i=1,loc_sz_i
           call fftw_execute_dft_r2c(this%fwy, this%f(i,:,k,l), this%tmp_y)
           do j=2,this%geomx%ny
              this%tmp_y(j) = this%tmp_y(j)*exp(-cmplx(0.0_f64,this%ky(j),kind=f64)*vy)
           enddo
           call fftw_execute_dft_c2r(this%bwy, this%tmp_y, this%d_dy)
           this%f(i,:,k,l) = this%d_dy / loc_sz_j
        end do
     end do
  end do
  
end subroutine advection_x2

subroutine advection_x3x4(this,dt)

  class(vlasov4d_spectral),intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64, dimension(this%geomv%nx,this%geomv%ny) :: alpha_x
  sll_real64, dimension(this%geomv%nx,this%geomv%ny) :: alpha_y
  sll_real64, dimension(this%geomv%nx,this%geomv%ny) :: alpha_x_tmp
  sll_real64, dimension(this%geomv%nx,this%geomv%ny) :: alpha_y_tmp
  sll_real64 :: px, py, ctheta, stheta, depvx, depvy
  sll_real64 :: x3_min, x3_max, x4_min, x4_max
  sll_real64 :: delta_x3, delta_x4

  x3_min   = this%geomv%x0
  x3_max   = this%geomv%x1
  delta_x3 = this%geomv%dx
  x4_min   = this%geomv%y0 
  x4_max   = this%geomv%y1
  delta_x4 = this%geomv%dy

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


 subroutine densite_courantx(this,star)

   class(vlasov4d_spectral),intent(inout)  :: this
   character(len=1), intent(in), optional  :: star
   sll_real64, dimension(:,:,:,:), pointer :: df

   sll_int32  :: error
   sll_int32  :: c
   sll_int32  :: comm
   sll_real64 :: dxy
   sll_real64 :: vx 
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjx

   if( present(star)) then
!      df => this%ft_star
      df => this%f_star
   else
!      df => this%ft
      df => this%f
   end if
   
   dxy = this%geomv%dx*this%geomv%dy
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
               vx = this%geomv%x0+(gk-1)*this%geomv%dx
               locjx(gi,gj) = locjx(gi,gj) + dxy*df(i,j,k,l) 
            end do
         end do
      end do
   end do

   this%jx1 = 0._f64
   comm = sll_world_collective%comm
   c    = this%geomx%nx*this%geomx%ny
   
   call mpi_barrier(comm,error)
   call mpi_allreduce(locjx,this%jx1,c,MPI_REAL8,MPI_SUM,comm,error)
   
 end subroutine densite_courantx

 subroutine densite_couranty(this, star)

   class(vlasov4d_spectral),intent(inout) :: this
   character(len=1), intent(in), optional  :: star

   sll_int32  :: error
   sll_real64 :: vy 
   sll_int32  :: c
   sll_int32  :: comm
   sll_real64 :: dxy
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjy
   sll_real64, dimension(:,:,:,:), pointer :: df

   if( present(star)) then
!      df => this%ft_star
      df => this%f_star
   else
!      df => this%ft
      df => this%f
   end if


   dxy = this%geomv%dx*this%geomv%dy
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
      vy = this%geomv%y0+(gl-1)*this%geomv%dy
      locjy(gi,gj) = locjy(gi,gj) + dxy*df(i,j,k,l) 
   end do
   end do
   end do
   end do

   this%jy1(:,:) = 0.
   comm   = sll_world_collective%comm
   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locjy,this%jy1,c, MPI_REAL8,MPI_SUM,comm,error)

 end subroutine densite_couranty

end module sll_vlasov4d_spectral
