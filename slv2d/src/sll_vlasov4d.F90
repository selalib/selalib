module sll_vlasov4d

#include "selalib.h"

 use used_precision
 use geometry_module
 use diagnostiques_module
 use sll_module_interpolators_1d_base
 use sll_module_interpolators_2d_base
 use remapper

 implicit none
 private
 public :: new, free, densite_charge, densite_courant
 public :: transposexv, transposevx
 public :: advection_x1, advection_x2, advection_x3, advection_x4
 public :: advection_x3x4
 public :: thdiag, write_xmf_file

 type, public :: vlasov4d
   type(geometry)                           :: geomx
   type(geometry)                           :: geomv
   logical                                  :: transposed      
   sll_real64, dimension(:,:,:,:),  pointer :: f
   sll_real64, dimension(:,:,:,:),  pointer :: ft
   class(sll_interpolator_1d_base), pointer :: interp_x1
   class(sll_interpolator_1d_base), pointer :: interp_x2
   class(sll_interpolator_1d_base), pointer :: interp_x3
   class(sll_interpolator_1d_base), pointer :: interp_x4
   class(sll_interpolator_2d_base), pointer :: interp_x1x2
   class(sll_interpolator_2d_base), pointer :: interp_x3x4

   type(layout_4D), pointer                 :: layout_x
   type(layout_4D), pointer                 :: layout_v
   type(remap_plan_4D_real64), pointer      :: x_to_v 
   type(remap_plan_4D_real64), pointer      :: v_to_x

   sll_real64, dimension(:,:), pointer      :: ex
   sll_real64, dimension(:,:), pointer      :: ey
   sll_real64, dimension(:,:), pointer      :: jx
   sll_real64, dimension(:,:), pointer      :: jy
   sll_real64, dimension(:,:), pointer      :: bz
   sll_real64, dimension(:,:), pointer      :: rho
 end type vlasov4d

 sll_int32, private :: i, j, k, l
 sll_int32, private :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l
 sll_int32, private :: global_indices(4), gi, gj, gk, gl
 sll_int32, private :: ierr

 interface new
   module procedure new_vlasov4d_poisson
   module procedure new_vlasov4d_maxwell
 end interface

 interface free
   module procedure dealloc_vlasov4d
 end interface

contains

 subroutine new_vlasov4d_poisson(this,geomx,geomv,interp_x1,interp_x2,interp_x3,interp_x4,error)

  use sll_hdf5_io

  type(vlasov4d),intent(inout)            :: this
  type(geometry),intent(in)               :: geomx
  type(geometry),intent(in)               :: geomv
  class(sll_interpolator_1d_base), target :: interp_x1
  class(sll_interpolator_1d_base), target :: interp_x2
  class(sll_interpolator_1d_base), target :: interp_x3
  class(sll_interpolator_1d_base), target :: interp_x4
  sll_int32                               :: error

  sll_int32  :: nc_x1, nc_x2, nc_x3, nc_x4
  sll_real64 :: x1_min, x2_min, x3_min, x4_min
  sll_real64 :: x1_max, x2_max, x3_max, x4_max

  sll_int32  :: psize, prank, comm, file_id

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
  SLL_ALLOCATE(this%ft(loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l),ierr)

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

  SLL_ALLOCATE(this%ex(geomx%nx,geomx%ny),error)
  SLL_ALLOCATE(this%ey(geomx%nx,geomx%ny),error)
  SLL_ALLOCATE(this%rho(geomx%nx,geomx%ny),error)

  nullify(this%bz)
  nullify(this%jx)
  nullify(this%jy)

 end subroutine new_vlasov4d_poisson


 subroutine new_vlasov4d_maxwell(this,geomx,geomv,interp_x1,interp_x2,interp_x3,interp_x4,interp_x3x4,error)

  type(vlasov4d),intent(inout)            :: this
  type(geometry),intent(in)               :: geomx
  type(geometry),intent(in)               :: geomv
  class(sll_interpolator_1d_base), target :: interp_x1
  class(sll_interpolator_1d_base), target :: interp_x2
  class(sll_interpolator_1d_base), target :: interp_x3
  class(sll_interpolator_1d_base), target :: interp_x4
  class(sll_interpolator_2d_base), target :: interp_x3x4
  sll_int32                               :: error

  this%interp_x3x4 => interp_x3x4

  call new_vlasov4d_poisson(this,geomx,geomv,interp_x1,interp_x2,interp_x3,interp_x4,error)
  
  SLL_ALLOCATE(this%jx(geomx%nx,geomx%ny),error)
  SLL_ALLOCATE(this%jy(geomx%nx,geomx%ny),error)
  SLL_ALLOCATE(this%bz(geomx%nx,geomx%ny),error)

 end subroutine new_vlasov4d_maxwell

 subroutine dealloc_vlasov4d(this)

  type(vlasov4d),intent(out) :: this

  call delete_layout_4D(this%layout_x)
  call delete_layout_4D(this%layout_v)
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

 subroutine advection_x3(this,dt)

  type(vlasov4d), intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  SLL_ASSERT(this%transposed) 
  call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  print*,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l
  do l=1,loc_sz_l
  do j=1,loc_sz_j
  do i=1,loc_sz_i

     global_indices = local_to_global_4D(this%layout_v,(/i,j,1,1/)) 
     gi = global_indices(1)
     gj = global_indices(2)
     alpha = this%ex(gi,gj)*dt
     !print*,i,j,l
     this%ft(i,j,:,l) = this%interp_x3%interpolate_array_disp(loc_sz_k,this%ft(i,j,:,l),alpha)

  end do
  end do
  end do

 end subroutine advection_x3

 subroutine advection_x4(this,dt)

  type(vlasov4d),intent(inout) :: this
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

 subroutine densite_charge(this)

   type(vlasov4d),intent(inout) :: this
   sll_int32 :: error
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
   this%rho(:,:) = 0.
   comm  = sll_world_collective%comm
   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locrho,this%rho,c,MPI_REAL8,MPI_SUM,comm,error)

 end subroutine densite_charge

 subroutine densite_courant(this)

   type(vlasov4d),intent(inout) :: this
   sll_int32 :: error
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

   this%jx(:,:) = 0.; this%jy(:,:) = 0.
   comm   = sll_world_collective%comm
   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locjx,this%jx,c, MPI_REAL8,MPI_SUM,comm,error)
   call mpi_allreduce(locjy,this%jy,c, MPI_REAL8,MPI_SUM,comm,error)

 end subroutine densite_courant

 subroutine thdiag(this,nrj,t)

   type(vlasov4d),intent(inout) :: this
   sll_real64, intent(in) :: t,nrj  
   sll_real64,dimension(11) :: auxloc
   sll_int32 :: prank, psize
   sll_real64,dimension(13) :: aux
   sll_real64,dimension(0:9) :: diag
   sll_int32 :: comm, error
   sll_real64 :: cell_volume

   comm  = sll_world_collective%comm
   prank = sll_get_collective_rank(sll_world_collective)
   psize = sll_get_collective_size(sll_world_collective)

   auxloc = 0.0
   cell_volume = this%geomx%dx * this%geomx%dy * this%geomv%dx * this%geomv%dy
   auxloc(1) = cell_volume * sum(this%f) ! avg(f)
   auxloc(2) = cell_volume * sum(abs(this%f)) ! L1 norm
   auxloc(3) = cell_volume * sum(this%f*this%f) ! L2 norm
   
   call mpi_reduce(auxloc,aux(1:11),11,MPI_REAL8,MPI_SUM,MPI_MASTER,comm, error)

   if (prank == MPI_MASTER) then
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

   SLL_ASSERT(.not. this%transposed)
   call apply_remap_4D( this%x_to_v, this%f, this%ft )
   this%transposed = .true.

 end subroutine transposexv

 subroutine transposevx(this)

   type(vlasov4d),intent(inout) :: this

   SLL_ASSERT(this%transposed)
   call apply_remap_4D( this%v_to_x, this%ft, this%f )
   this%transposed = .false.

 end subroutine transposevx

 subroutine write_xmf_file(this, iplot)

  use hdf5
  use sll_hdf5_io

  type(vlasov4d),intent(inout) :: this
  sll_int32, intent(in)        :: iplot
  sll_int32                    :: error
  character(len=4)             :: cplot
  sll_int32                    :: prank
  sll_int32, parameter         :: one = 1
  sll_int32                    :: file_id
  sll_int32                    :: nx1, nx2, nx3, nx4

  call int2string(iplot,cplot)
  call write_fx1x2(this,cplot)
  call write_fx1x3(this,cplot)
  call write_fx2x4(this,cplot)
  call write_fx3x4(this,cplot)

  prank = sll_get_collective_rank(sll_world_collective)
  if (prank == MPI_MASTER) then

     
     if(associated(this%ex)) then
        call sll_hdf5_file_create('ex_'//cplot//".h5",file_id,error)
        call sll_hdf5_write_array(file_id,this%ex,"/values",error)
        call sll_hdf5_file_close(file_id, error)
     end if
     if(associated(this%ey)) then
        call sll_hdf5_file_create('ey_'//cplot//".h5",file_id,error)
        call sll_hdf5_write_array(file_id,this%ey,"/values",error)
        call sll_hdf5_file_close(file_id, error)
     end if
     if(associated(this%rho)) then
        call sll_hdf5_file_create('rho_'//cplot//".h5",file_id,error)
        call sll_hdf5_write_array(file_id,this%rho,"/values",error)
        call sll_hdf5_file_close(file_id, error)
     end if
     if(size(this%jx,1) > 1) then
        call sll_hdf5_file_create('jx_'//cplot//".h5",file_id,error)
        call sll_hdf5_write_array(file_id,this%jx,"/values",error)
        call sll_hdf5_file_close(file_id, error)
     end if
     if(associated(this%jy)) then
        call sll_hdf5_file_create('jy_'//cplot//".h5",file_id,error)
        call sll_hdf5_write_array(file_id,this%jy,"/values",error)
        call sll_hdf5_file_close(file_id, error)
     end if
     if(associated(this%bz)) then
        call sll_hdf5_file_create('bz_'//cplot//".h5",file_id,error)
        call sll_hdf5_write_array(file_id,this%bz,"/values",error)
        call sll_hdf5_file_close(file_id, error)
     end if

     nx1 = this%geomx%nx
     nx2 = this%geomv%nx
     nx3 = this%geomx%ny
     nx4 = this%geomv%ny

     call sll_xml_file_create("fvalues_"//cplot//".xmf",file_id,error)
     call write_grid(this,file_id,nx1,nx2,"x1","x2",cplot)
     call write_grid(this,file_id,nx1,nx3,"x1","x3",cplot)
     call write_grid(this,file_id,nx2,nx4,"x2","x4",cplot)
     call write_grid(this,file_id,nx3,nx4,"x3","x4",cplot)
     write(file_id,"(a)")"</Domain>"
     write(file_id,"(a)")"</Xdmf>"
     close(file_id)

  endif

 end subroutine write_xmf_file
 
 subroutine write_grid(this,file_id,nx,ny,xname,yname,cplot)

  type(vlasov4d),intent(in) :: this
  sll_int32 :: file_id, nx, ny
  character(len=*) :: cplot, xname, yname

  write(file_id,"(a)")"<Grid Name='"//xname//yname//"' GridType='Uniform'>"
  write(file_id, &
   "(a,2i5,a)")"<Topology TopologyType='2DRectMesh' NumberOfElements='",ny,nx,"'/>"
  write(file_id,"(a)")"<Geometry GeometryType='VXVY'>"
  write(file_id,"(a,i5,a)")"<DataItem Dimensions='",nx, &
                           "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")"mesh4d.h5:/"//xname
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a,i5,a)")"<DataItem Dimensions='",ny, &
                           "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")"mesh4d.h5:/"//yname
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a)")"</Geometry>"
  if (xname == "x1" .and. yname == "x2") then
     if(associated(this%ex)) &
     call write_attribute(file_id,nx,ny,"ex",xname,yname,cplot)
     if(associated(this%ey)) &
     call write_attribute(file_id,nx,ny,"ey",xname,yname,cplot)
     if(associated(this%rho)) &
     call write_attribute(file_id,nx,ny,"rho",xname,yname,cplot)
     if(associated(this%bz)) &
     call write_attribute(file_id,nx,ny,"bz",xname,yname,cplot)
     if(associated(this%jx)) &
     call write_attribute(file_id,nx,ny,"jx",xname,yname,cplot)
     if(associated(this%jy)) &
     call write_attribute(file_id,nx,ny,"jy",xname,yname,cplot)
  end if
  call write_attribute(file_id,nx,ny,"df",xname,yname,cplot)
  write(file_id,"(a)")"</Grid>"

 end subroutine write_grid

 subroutine write_attribute(file_id,nx,ny,fname,xname,yname,cplot)
 sll_int32                    :: file_id
 sll_int32                    :: nx
 sll_int32                    :: ny
 character(len=*), intent(in) :: fname
 character(len=*), intent(in) :: xname
 character(len=*), intent(in) :: yname
 character(len=*), intent(in) :: cplot

  write(file_id,"(a)")"<Attribute Name='"//fname//"' AttributeType='Scalar' Center='Node'>"
  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ny,nx, &
                            "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")"f"//xname//yname//"_"//cplot//".h5:/values"
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a)")"</Attribute>"

 end subroutine write_attribute

 subroutine advection_x3x4(this,dt)

  type(vlasov4d),intent(inout) :: this
  sll_real64, intent(in) :: dt
  sll_real64, dimension(this%geomv%nx,this%geomv%ny) :: alpha_x
  sll_real64, dimension(this%geomv%nx,this%geomv%ny) :: alpha_y
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
        alpha_x(k,l) = px - (depvx+(px+depvx)*ctheta-(py+depvy)*stheta)
        alpha_y(k,l) = py - (depvy+(px+depvx)*stheta+(py+depvy)*ctheta)

     end do
     end do

     this%ft(i,j,:,:) = this%interp_x3x4%interpolate_array_disp(loc_sz_k,loc_sz_l, &
                                                 this%ft(i,j,:,:),alpha_x,alpha_y)
  end do
  end do

 end subroutine advection_x3x4

 subroutine write_fx1x2(this,cplot)
 use sll_hdf5_io
 type(vlasov4d),intent(in)           :: this
 character(len=*)                    :: cplot
 sll_int32                           :: error
 sll_int32                           :: file_id
 sll_int32                           :: prank
 sll_int32                           :: comm
 sll_real64, dimension(:,:), pointer :: fij
 sll_real64                          :: sumloc

 prank = sll_get_collective_rank(sll_world_collective)
 comm  = sll_world_collective%comm
 call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
 SLL_ALLOCATE(fij(loc_sz_i,loc_sz_j),error)
 do j=1,loc_sz_j
    do i=1,loc_sz_i
       sumloc = sum(this%f(i,j,:,:))
       call mpi_reduce(sumloc,fij(i,j),1,MPI_REAL8,MPI_SUM,MPI_MASTER,comm,error)
    end do
 end do
 if (prank == MPI_MASTER) then
    call sll_hdf5_file_create('fx1x2_'//cplot//".h5",file_id,error)
    call sll_hdf5_write_array(file_id,fij,"/values",error)
    call sll_hdf5_file_close(file_id, error)
 end if
 end subroutine write_fx1x2

 subroutine write_fx1x3(this,cplot)
 use sll_hdf5_io
 type(vlasov4d),intent(in)           :: this
 character(len=*)                    :: cplot
 sll_int32                           :: error
 sll_int32                           :: file_id
 sll_real64, dimension(:,:), pointer :: fik
 sll_int32                           :: prank
 sll_int32                           :: comm
 sll_real64                          :: sumloc

 prank = sll_get_collective_rank(sll_world_collective)
 comm  = sll_world_collective%comm
 call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
 SLL_ALLOCATE(fik(loc_sz_i,loc_sz_k),error)
 do k=1,loc_sz_k
    do i=1,loc_sz_i
       sumloc= sum(this%f(i,:,k,:))
       call mpi_reduce(sumloc,fik(i,k),1,MPI_REAL8,MPI_SUM,MPI_MASTER,comm,error)
    end do
 end do
 if (prank == MPI_MASTER) then
    call sll_hdf5_file_create('fx1x3_'//cplot//".h5",file_id,error)
    call sll_hdf5_write_array(file_id,fik,"/values",error)
    call sll_hdf5_file_close(file_id, error)
 end if
 end subroutine write_fx1x3

 subroutine write_fx2x4(this,cplot)
 use hdf5
 use sll_hdf5_io_parallel
 type(vlasov4d),intent(in)           :: this
 character(len=*)                    :: cplot
 integer(HID_T)                      :: pfile_id
 integer(HSSIZE_T)                   :: offset(2)
 integer(HSIZE_T)                    :: global_dims(2)
 sll_int32                           :: error
 sll_int32                           :: prank
 sll_real64, dimension(:,:), pointer :: fjl

 prank = sll_get_collective_rank(sll_world_collective)
 call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
 SLL_ALLOCATE(fjl(loc_sz_j,loc_sz_l),error)
 do l=1,loc_sz_l
    do j=1,loc_sz_j
       fjl(j,l) = sum(this%f(:,j,:,l))
    end do
 end do
 global_dims = (/this%geomx%ny,this%geomv%ny/)
 offset(1) = get_layout_4D_j_min(this%layout_x,prank)-1
 offset(2) = get_layout_4D_l_min(this%layout_x,prank)-1
 call sll_hdf5_file_create('fx2x4_'//cplot//".h5",pfile_id,error)
 call sll_hdf5_write_array_2d(pfile_id,global_dims,offset,fjl,"/values",error)
 call sll_hdf5_file_close(pfile_id, error)

 end subroutine write_fx2x4

 subroutine write_fx3x4(this,cplot)
 use hdf5
 use sll_hdf5_io_parallel
 type(vlasov4d),intent(in)           :: this
 character(len=*)                    :: cplot
 integer(HID_T)                      :: pfile_id
 integer(HSSIZE_T)                   :: offset(2)
 integer(HSIZE_T)                    :: global_dims(2)
 sll_int32                           :: error
 sll_int32                           :: prank
 sll_real64, dimension(:,:), pointer :: fkl

 prank = sll_get_collective_rank(sll_world_collective)
 call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
 SLL_ALLOCATE(fkl(loc_sz_k,loc_sz_l),error)
 do l=1,loc_sz_l
    do k=1,loc_sz_k
       fkl(k,l) = sum(this%f(:,:,k,l))
    end do
 end do
 global_dims = (/this%geomv%nx,this%geomv%ny/)
 offset(1) = get_layout_4D_k_min(this%layout_x,prank)-1
 offset(2) = get_layout_4D_l_min(this%layout_x,prank)-1
 call sll_hdf5_file_create('fx3x4_'//cplot//".h5",pfile_id,error)
 call sll_hdf5_write_array_2d(pfile_id,global_dims,offset,fkl,"/values",error)
 call sll_hdf5_file_close(pfile_id, error)

 end subroutine write_fx3x4

end module sll_vlasov4d
