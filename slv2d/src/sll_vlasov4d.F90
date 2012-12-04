module sll_vlasov4d

#include "selalib.h"

 use used_precision
 use geometry_module
 use diagnostiques_module
 use sll_module_interpolators_1d_base
 use remapper

 implicit none
 private
 public :: new, free, densite_charge,densite_courant
 public :: transposexv, transposevx
 public :: advection_x1, advection_x2, advection_x3, advection_x4
 public :: thdiag

 type, public :: vlasov4d
   type (geometry)                          :: geomx
   type (geometry)                          :: geomv
   logical                                  :: transposed      
   sll_real64, dimension(:,:,:,:), pointer  :: f
   sll_real64, dimension(:,:,:,:), pointer  :: ft
   class(sll_interpolator_1d_base), pointer :: interp_x1
   class(sll_interpolator_1d_base), pointer :: interp_x2
   class(sll_interpolator_1d_base), pointer :: interp_x3
   class(sll_interpolator_1d_base), pointer :: interp_x4
   type(layout_4D), pointer                 :: layout_x
   type(layout_4D), pointer                 :: layout_v
   type(remap_plan_4D_real64), pointer      :: x_to_v 
   type(remap_plan_4D_real64), pointer      :: v_to_x
 end type vlasov4d

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

contains

 subroutine new_vlasov4d(this,geomx,geomv,interp_x1,interp_x2,interp_x3,interp_x4,error)

  type(vlasov4d),intent(out)              :: this
  type(geometry),intent(in)               :: geomx
  type(geometry),intent(in)               :: geomv
  sll_int32, intent(out)                  :: error
  class(sll_interpolator_1d_base), target :: interp_x1
  class(sll_interpolator_1d_base), target :: interp_x2
  class(sll_interpolator_1d_base), target :: interp_x3
  class(sll_interpolator_1d_base), target :: interp_x4

  sll_int32  :: nc_x1, nc_x2, nc_x3, nc_x4
  sll_real64 :: x1_min, x2_min, x3_min, x4_min
  sll_real64 :: x1_max, x2_max, x3_max, x4_max

  sll_int32 :: psize, prank, comm

  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  error = 0

  this%transposed=.false.

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

  this%layout_x => new_layout_4D( sll_world_collective )        
  call initialize_layout_with_distributed_4D_array( &
             geomx%nx,geomx%ny,geomv%nx, geomv%ny, &
             1,1,1,int(psize,4),this%layout_x)


  call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_ALLOCATE(this%f(loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l),ierr)

  this%layout_v => new_layout_4D( sll_world_collective )
  call initialize_layout_with_distributed_4D_array( &
              geomx%nx,geomx%ny,geomv%nx, geomv%ny, &
              1,int(psize,4),1,1,this%layout_v)

  call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_ALLOCATE(this%ft(loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l),ierr )

  this%x_to_v => new_remap_plan( this%layout_x, this%layout_v, this%f)     
  this%v_to_x => new_remap_plan( this%layout_v, this%layout_x, this%ft)     
  
  if(prank == MPI_MASTER) then
    print *,'Printing layout x: '
    call sll_view_lims_4D( this%layout_x )
    print *,'Printing layout v: '
    call sll_view_lims_4D( this%layout_v )
  end if

 end subroutine new_vlasov4d

 subroutine dealloc_vlasov4d(this)
  type(vlasov4d),intent(out) :: this

  call delete_layout_4D( this%layout_x )
  call delete_layout_4D( this%layout_v )
  SLL_DEALLOCATE_ARRAY(this%f, ierr)
  SLL_DEALLOCATE_ARRAY(this%ft, ierr)

 end subroutine dealloc_vlasov4d

 subroutine advection_x1(this,dt)
  type(vlasov4d), intent(inout) :: this
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
  type(vlasov4d),intent(inout) :: this
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

 subroutine advection_x3(this,dt,ex,ey,bz)

  type(vlasov4d), intent(inout)          :: this
  sll_real64, dimension(:,:), intent(in) :: ex
  sll_real64, dimension(:,:), optional   :: ey
  sll_real64, dimension(:,:), optional   :: bz
  sll_real64, intent(in) :: dt
  sll_real64 :: x3_min, delta_x3
  sll_real64 :: x4_min, delta_x4
  sll_real64 :: alpha
  sll_real64 :: depvx, depvy, px, py, ctheta, stheta

  SLL_ASSERT(this%transposed) 

  x3_min   = this%geomv%x0
  delta_x3 = this%geomv%dx
  x4_min   = this%geomv%y0
  delta_x4 = this%geomv%dy

  call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  if(present(bz) .and. present(ey)) then

     do l=1,loc_sz_l
     do k=1,loc_sz_k
     do j=1,loc_sz_j
     do i=1,loc_sz_i

        global_indices = local_to_global_4D(this%layout_v,(/i,j,k,l/)) 
        gi = global_indices(1)
        gj = global_indices(2)
        gk = global_indices(3)
        gl = global_indices(4)
        px = x3_min+(gk-1)*delta_x3
        py = x4_min+(gl-1)*delta_x4
        ctheta = cos(bz(gi,gj)*dt)
        stheta = sin(bz(gi,gj)*dt)
        depvx  = -0.5*dt*ex(gi,gj)
        depvy  = -0.5*dt*ey(gi,gj)
        alpha  = depvx+(px+depvx)*ctheta-(py+depvy)*stheta
        this%ft(i,j,k,l) = this%interp_x3%interpolate_array_value(this%ft(i,j,k,l),alpha)

     end do
     end do
     end do
     end do

  else

     do l=1,loc_sz_l
     do j=1,loc_sz_j
     do i=1,loc_sz_i

        global_indices = local_to_global_4D(this%layout_v,(/i,j,1,1/)) 
        gi = global_indices(1)
        gj = global_indices(2)
        alpha = ex(gi,gj)*dt
        this%ft(i,j,:,l) = this%interp_x3%interpolate_array_disp(loc_sz_k,this%ft(i,j,:,l),alpha)

     end do
     end do
     end do

 end if

 end subroutine advection_x3

 subroutine advection_x4(this,dt,ey,ex,bz)

  type(vlasov4d),intent(inout) :: this
  sll_real64, dimension(:,:), intent(in) :: ey
  sll_real64, dimension(:,:), optional   :: ex
  sll_real64, dimension(:,:), optional   :: bz
  sll_real64, intent(in) :: dt
  sll_real64 :: x3_min, delta_x3
  sll_real64 :: x4_min, delta_x4
  sll_real64 :: alpha
  sll_real64 :: depvx, depvy, px, py, ctheta, stheta

  SLL_ASSERT(this%transposed) 
  call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  if (present(bz) .and. present(ex)) then

     do l=1,loc_sz_l
     do k=1,loc_sz_k
     do j=1,loc_sz_j
     do i=1,loc_sz_i

        global_indices = local_to_global_4D(this%layout_v,(/i,j,k,l/)) 
        gi = global_indices(1)
        gj = global_indices(2)
        gk = global_indices(3)
        gl = global_indices(4)
        px = x3_min+(gk-1)*delta_x3
        py = x4_min+(gl-1)*delta_x4
        ctheta = cos(bz(gi,gj)*dt)
        stheta = sin(bz(gi,gj)*dt)
        depvx  = -0.5*dt*ex(gi,gj)
        depvy  = -0.5*dt*ey(gi,gj)
        alpha  = depvy+(px+depvx)*stheta+(py+depvy)*ctheta
        this%ft(i,j,k,l) = this%interp_x3%interpolate_array_value(this%ft(i,j,k,l),alpha)

     end do
     end do
     end do
     end do

  else

     do k=1,loc_sz_k
     do j=1,loc_sz_j
     do i=1,loc_sz_i

        global_indices = local_to_global_4D(this%layout_v,(/i,j,1,1/)) 
        gi = global_indices(1)
        gj = global_indices(2)
        alpha = ey(gi,gj)*dt
        this%ft(i,j,k,:) = this%interp_x4%interpolate_array_disp(loc_sz_l,this%ft(i,j,k,:),alpha)

    end do
    end do
    end do

  end if

 end subroutine advection_x4

 subroutine densite_charge(this, rho)

   type(vlasov4d),intent(inout) :: this
   sll_int32 :: error
   sll_real64, dimension(:,:), intent(out)  :: rho
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locrho
   sll_int32  :: c
   sll_int32  :: comm
   sll_real64 :: dxy

   dxy = this%geomv%dx*this%geomv%dy

   SLL_ASSERT(this%transposed)
   call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
   
   locrho(:,:) = 0.
   do j=1,loc_sz_j
   do i=1,loc_sz_i
      global_indices = local_to_global_4D(this%layout_v,(/i,j,1,1/)) 
      gi = global_indices(1)
      gj = global_indices(2)
      locrho(gi,gj) = sum(this%ft(i,j,:,:))*dxy 
   end do
   end do
   rho(:,:) = 0.
   comm  = sll_world_collective%comm
   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locrho,rho,c,MPI_REAL8,MPI_SUM,comm,error)

 end subroutine densite_charge

 subroutine densite_courant(this, jx, jy)

   type(vlasov4d),intent(inout) :: this
   sll_int32 :: error
   sll_real64, dimension(:,:), intent(out)  :: jx, jy
   sll_real64 :: vx, vy 
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjx
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjy
   sll_int32  :: c
   sll_int32  :: comm
   sll_real64 :: dxy

   dxy = this%geomv%dx*this%geomv%dy
   SLL_ASSERT(this%transposed)
   call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

   locjx(:,:) = 0.; locjy(:,:) = 0.
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
      vy = this%geomv%y0+(gl-1)*this%geomv%dy
      locjx(gi,gj) = locjx(gi,gj) + dxy*this%ft(i,j,k,l) * vx
      locjy(gi,gj) = locjy(gi,gj) + dxy*this%ft(i,j,k,l) * vy
   end do
   end do
   end do
   end do

   jx(:,:) = 0.; jy(:,:) = 0.
   comm   = sll_world_collective%comm
   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locjx,jx,c, MPI_REAL8,MPI_SUM,comm,error)
   call mpi_allreduce(locjy,jy,c, MPI_REAL8,MPI_SUM,comm,error)

 end subroutine densite_courant

 subroutine thdiag(this,nrj,t)

   type(vlasov4d),intent(inout) :: this
   sll_real64, intent(in) :: t,nrj  
   sll_real64,dimension(11) :: auxloc
   sll_int32 :: my_num, num_threads
   sll_real64,dimension(13) :: aux
   sll_real64,dimension(0:9) :: diag
   sll_int32 :: comm, error
   sll_real64 :: cell_volume

   comm   = sll_world_collective%comm
   my_num = sll_get_collective_rank(sll_world_collective)
   num_threads = sll_get_collective_size(sll_world_collective)

   auxloc = 0.0
   cell_volume = this%geomx%dx * this%geomx%dy * this%geomv%dx * this%geomv%dy
   auxloc(1) = cell_volume * sum(this%f) ! avg(f)
   auxloc(2) = cell_volume * sum(abs(this%f)) ! L1 norm
   auxloc(3) = cell_volume * sum(this%f*this%f) ! L2 norm
   
   call mpi_reduce(auxloc,aux(1:11),11,MPI_REAL8,MPI_SUM,MPI_MASTER,comm, error)

   if (my_num == MPI_MASTER) then
      diag=0.
      aux=0.
      aux(13)=t
      aux(12)=nrj
      write(*,"('time ', g8.3,' test nrj',f10.5)") t, nrj
      call time_history("thf","(13(1x,e15.6))",aux(1:13))
   end if

 end subroutine thdiag

 subroutine transposexv(this)

   type(vlasov4d),intent(inout) :: this

   call apply_remap_4D( this%x_to_v, this%f, this%ft )

 end subroutine transposexv

 subroutine transposevx(this)

   type(vlasov4d),intent(inout) :: this

   call apply_remap_4D( this%v_to_x, this%ft, this%f )

 end subroutine transposevx

end module sll_vlasov4d
