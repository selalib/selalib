#define FFTW_ALLOCATE(array,array_size,sz_array,p_array)  \
sz_array = int((array_size/2+1),C_SIZE_T);                \
p_array = fftw_alloc_complex(sz_array);                   \
call c_f_pointer(p_array, array, [array_size/2+1])        \

module sll_vlasov4d_spectral

#include "selalib.h"

 use, intrinsic :: iso_c_binding
 use used_precision
 use geometry_module
 use diagnostiques_module
 use sll_module_interpolators_1d_base
 use sll_module_interpolators_2d_base
 use remapper
 use sll_vlasov4d, only: vlasov4d

 implicit none
 private
 public :: new, free, densite_courantx, densite_couranty
 public :: advection_x1, advection_x2

 type, public, extends(vlasov4d) :: vlasov4d_spectral

   sll_real64, dimension(:,:), pointer               :: exn
   sll_real64, dimension(:,:), pointer               :: eyn
   sll_real64, dimension(:,:), pointer               :: jx1,jx2
   sll_real64, dimension(:,:), pointer               :: jy1,jy2
   sll_real64, dimension(:), allocatable             :: d_dx
   sll_real64, dimension(:), allocatable             :: d_dy
   sll_real64, dimension(:), allocatable             :: kx
   sll_real64, dimension(:), allocatable             :: ky
   type(C_PTR)                                       :: fwx, fwy
   type(C_PTR)                                       :: bwx, bwy
   type(C_PTR)                                       :: p_tmp_x, p_tmp_y
   complex(C_DOUBLE_COMPLEX), dimension(:),  pointer :: tmp_x, tmp_y

 end type vlasov4d_spectral

 sll_int32, private :: i, j, k, l
 sll_int32, private :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l
 sll_int32, private :: global_indices(4), gi, gj, gk, gl
 sll_int32, private :: ierr

 interface new
   module procedure new_vlasov4d
 end interface

 interface free
   module procedure dealloc_vlasov4d
 end interface

include 'fftw3.f03'

contains

 subroutine new_vlasov4d(this,geomx,geomv,interp_x1,interp_x2,interp_x3,interp_x4,interp_x3x4,error)

  use sll_hdf5_io

  type(vlasov4d_spectral),intent(inout)            :: this
  type(geometry),intent(in)               :: geomx
  type(geometry),intent(in)               :: geomv
  class(sll_interpolator_1d_base), target :: interp_x1
  class(sll_interpolator_1d_base), target :: interp_x2
  class(sll_interpolator_1d_base), target :: interp_x3
  class(sll_interpolator_1d_base), target :: interp_x4
  class(sll_interpolator_2d_base), target :: interp_x3x4
  sll_int32                               :: error

  sll_int32  :: nc_x1, nc_x2, nc_x3, nc_x4
  sll_real64 :: x1_min, x2_min, x3_min, x4_min
  sll_real64 :: x1_max, x2_max, x3_max, x4_max
  sll_int32  :: psize, prank, comm, file_id
  sll_real64 :: dx, dy, kx0, ky0
  integer(C_SIZE_T) :: sz_tmp_x, sz_tmp_y

  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  if (.not. is_power_of_two(int(psize,i64))) then     
     print *, 'This test needs to run in a number of processes which is ',&
          'greater than 4 and a power of 2.'
     call sll_halt_collective()
     stop
  end if

  error = 0

  this%transposed = .false.

  nc_x1  = geomx%nx
  nc_x2  = geomx%ny
  nc_x3  = geomv%nx
  nc_x4  = geomv%ny

  x1_min = geomx%x0
  x1_max = geomx%x1
  x2_min = geomx%y0
  x2_max = geomx%y1
  x3_min = geomv%x0
  x3_max = geomv%x1
  x4_min = geomv%y0
  x4_max = geomv%y1

  this%geomx=geomx
  this%geomv=geomv

  this%interp_x1 => interp_x1
  this%interp_x2 => interp_x2
  this%interp_x3 => interp_x3
  this%interp_x4 => interp_x4
  this%interp_x3x4 => interp_x3x4

  this%layout_x => new_layout_4D( sll_world_collective )        
  call initialize_layout_with_distributed_4D_array( &
             geomx%nx,geomx%ny,geomv%nx, geomv%ny, &
             1,1,1,int(psize,4),this%layout_x)


  call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(this%f(loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l),ierr)

  this%layout_v => new_layout_4D( sll_world_collective )
  call initialize_layout_with_distributed_4D_array( &
              geomx%nx,geomx%ny,geomv%nx, geomv%ny, &
              1,int(psize,4),1,1,this%layout_v)

  call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(this%ft(loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l),ierr)

  this%x_to_v => new_remap_plan( this%layout_x, this%layout_v, this%f)     
  this%v_to_x => new_remap_plan( this%layout_v, this%layout_x, this%ft)     
  
  if(prank == MPI_MASTER) then

     print *,'Printing layout x: '
     call sll_view_lims_4D(this%layout_x)
     print *,'Printing layout v: '
     call sll_view_lims_4D(this%layout_v)

     call sll_hdf5_file_create("mesh4d.h5",file_id,error)
     call sll_hdf5_write_array(file_id,this%geomx%xgrid,"/x1",error)
     call sll_hdf5_write_array(file_id,this%geomx%ygrid,"/x2",error)
     call sll_hdf5_write_array(file_id,this%geomv%xgrid,"/x3",error)
     call sll_hdf5_write_array(file_id,this%geomv%ygrid,"/x4",error)
     call sll_hdf5_file_close(file_id, error)

  end if

  SLL_CLEAR_ALLOCATE(this%ex(nc_x1,nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%ey(nc_x1,nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%exn(nc_x1,nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%eyn(nc_x1,nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%bz(nc_x1,nc_x2),error); this%bz = 0.0_f64
  SLL_CLEAR_ALLOCATE(this%rho(nc_x1,nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jx(nc_x1,nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jx1(nc_x1,nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jx2(nc_x1,nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jy(nc_x1,nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jy1(nc_x1,nc_x2),error)
  SLL_CLEAR_ALLOCATE(this%jy2(nc_x1,nc_x2),error)
  
  FFTW_ALLOCATE(this%tmp_x,nc_x1/2+1,sz_tmp_x,this%p_tmp_x)
  FFTW_ALLOCATE(this%tmp_y,nc_x2/2+1,sz_tmp_y,this%p_tmp_y)
  SLL_CLEAR_ALLOCATE(this%d_dx(nc_x1), error)
  SLL_CLEAR_ALLOCATE(this%d_dy(nc_x2), error)

  this%fwx = fftw_plan_dft_r2c_1d(nc_x1, this%d_dx,  this%tmp_x, FFTW_ESTIMATE)
  this%bwx = fftw_plan_dft_c2r_1d(nc_x1, this%tmp_x, this%d_dx,  FFTW_ESTIMATE)
  this%fwy = fftw_plan_dft_r2c_1d(nc_x2, this%d_dy,  this%tmp_y, FFTW_ESTIMATE)
  this%bwy = fftw_plan_dft_c2r_1d(nc_x2, this%tmp_y, this%d_dy,  FFTW_ESTIMATE)

  SLL_CLEAR_ALLOCATE(this%kx(nc_x1/2+1), error)
  SLL_CLEAR_ALLOCATE(this%ky(nc_x2/2+1), error)
   
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

 end subroutine new_vlasov4d

 subroutine dealloc_vlasov4d(this)

  type(vlasov4d_spectral),intent(out) :: this

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

 end subroutine dealloc_vlasov4d

 subroutine advection_x1(this,dt)
  type(vlasov4d_spectral), intent(inout) :: this
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
!exact
!        this%tmp_x = this%tmp_x*exp(-cmplx(0.0_f64,this%kx,kind=f64)*vx)
!Euler explicite
!        this%tmp_x = this%tmp_x*(1._f64-cmplx(0.0_f64,this%kx,kind=f64)*vx)
!Euler implicite
!        this%tmp_x = this%tmp_x/(1._f64+cmplx(0.0_f64,this%kx,kind=f64)*vx)
!crank-nicolson
!        this%tmp_x = this%tmp_x/(1._f64+cmplx(0.0_f64,this%kx,kind=f64)*vx*0.5_f64)
!        this%tmp_x = this%tmp_x*(1._f64-cmplx(0.0_f64,this%kx,kind=f64)*vx*0.5_f64)
!Euler cn modified
        this%tmp_x = this%tmp_x*(1._f64-cmplx(0.0_f64,this%kx,kind=f64)*vx-0.5_f64*(this%kx*vx)**2)
        call fftw_execute_dft_c2r(this%bwx, this%tmp_x, this%d_dx)
        this%f(:,j,k,l)= this%d_dx / loc_sz_i
     end do
  end do
  end do

!comments NC: il faudrait calculer le courantx de ftmp donnee par 

 end subroutine advection_x1

 subroutine advection_x2(this,dt)
  type(vlasov4d_spectral),intent(inout) :: this
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
!exact
!       this%tmp_y = this%tmp_y*exp(-cmplx(0.0_f64,this%ky,kind=f64)*vy)
!euler explicite
!       this%tmp_y = this%tmp_y*(1._f64-cmplx(0.0_f64,this%ky,kind=f64)*vy)
!euler implicite
!       this%tmp_y = this%tmp_y/(1._f64+cmplx(0.0_f64,this%ky,kind=f64)*vy)
!crank-nicolson
!       this%tmp_y = this%tmp_y/(1._f64+cmplx(0.0_f64,this%ky,kind=f64)*vy*0.5_f64)
!       this%tmp_y = this%tmp_y*(1._f64-cmplx(0.0_f64,this%ky,kind=f64)*vy*0.5_f64)
!Euler cn modified
       this%tmp_y = this%tmp_y*(1._f64-cmplx(0.0_f64,this%ky,kind=f64)*vy-0.5_f64*(this%ky*vy)**2)
       call fftw_execute_dft_c2r(this%bwy, this%tmp_y, this%d_dy)
       this%f(i,:,k,l) = this%d_dy / loc_sz_j
    end do
    end do

  end do

 end subroutine advection_x2

 subroutine densite_courantx(this)

   type(vlasov4d_spectral),intent(inout) :: this
   sll_int32 :: error
   sll_real64 :: vx 
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjx
   sll_int32  :: c
   sll_int32  :: comm
   sll_real64 :: dxy

   dxy = this%geomv%dx*this%geomv%dy
   SLL_ASSERT(this%transposed)
   call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

   locjx(:,:) = 0.
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
      locjx(gi,gj) = locjx(gi,gj) + dxy*this%ft(i,j,k,l) * vx
   end do
   end do
   end do
   end do

   this%jx1(:,:) = 0.
   comm   = sll_world_collective%comm
   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locjx,this%jx1,c, MPI_REAL8,MPI_SUM,comm,error)

 end subroutine densite_courantx

 subroutine densite_couranty(this)

   type(vlasov4d_spectral),intent(inout) :: this
   sll_int32 :: error
   sll_real64 :: vy 
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjy
   sll_int32  :: c
   sll_int32  :: comm
   sll_real64 :: dxy

   dxy = this%geomv%dx*this%geomv%dy
   SLL_ASSERT(this%transposed)
   call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

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
      locjy(gi,gj) = locjy(gi,gj) + dxy*this%ft(i,j,k,l) * vy
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
