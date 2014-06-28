module vlasov4d_plot

#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"
use sll_collective
use used_precision  
use geometry_module
use sll_xml_io

implicit none
sll_real64, dimension(:,:),pointer :: fxy
sll_real64, dimension(:,:),pointer :: fxvx
sll_real64, dimension(:,:),pointer :: fyvy
sll_real64, dimension(:,:),pointer :: fvxvy
sll_int32, private :: i, j, k, l

enum, bind(C)
enumerator :: XY_VIEW = 0, XVX_VIEW  = 1, YVY_VIEW = 2, VXVY_VIEW = 3
end enum

contains

 subroutine plot_mesh4d(geomx, geomv, jstartx,jendx,jstartv,jendv)

  use sll_hdf5_io_serial
  type(geometry)        :: geomx, geomv
  sll_int32             :: file_id
  sll_int32, intent(in) :: jstartx, jendx, jstartv, jendv
  sll_int32             :: nx, nvx, ny, nvy
  sll_int32             :: error
  sll_int32             :: my_num

  nx  = geomx%nx; ny  = geomx%ny
  nvx = geomv%nx; nvy = geomv%ny
  SLL_ALLOCATE(fxy(nx,ny),error)
  SLL_ALLOCATE(fxvx(nx,nvx),error)
  SLL_ALLOCATE(fyvy(ny,jstartv:jendv),error)
  SLL_ALLOCATE(fvxvy(nvx,jstartv:jendv),error)

  my_num = sll_get_collective_rank(sll_world_collective)

  if (my_num == MPI_MASTER) then
     call sll_hdf5_file_create("grid4d.h5",file_id,error)
     call sll_hdf5_write_array(file_id,geomx%xgrid,"/x",error)
     call sll_hdf5_write_array(file_id,geomx%ygrid,"/y",error)
     call sll_hdf5_write_array(file_id,geomv%xgrid,"/vx",error)
     call sll_hdf5_write_array(file_id,geomv%ygrid,"/vy",error)
     call sll_hdf5_file_close(file_id, error)
  end if

 end subroutine plot_mesh4d

subroutine plot_df(f4d,iplot,geomx,geomv,jstartx,jendx,jstartv,jendv,choice)
  
  use sll_hdf5_io_serial
  type(geometry), intent(in) :: geomx, geomv
  sll_real64, intent(in) :: f4d(:,:,:,jstartv:)
  sll_int32, intent(in)  :: jstartx, jendx, jstartv, jendv
  
  sll_int32, intent(in) :: iplot
  character(len=4)      :: cplot
  sll_int32         :: file_id
  sll_real64        :: sumloc
  sll_int32         :: nx , nvx, ny, nvy
  sll_int32         :: my_num
  sll_int32         :: error
  sll_int32         :: comm
  sll_int32         :: choice
  
  my_num = sll_get_collective_rank(sll_world_collective)
  comm = sll_world_collective%comm
  
  nx  = geomx%nx; ny  = geomx%ny
  nvx = geomv%nx; nvy = geomv%ny
  
  call int2string(iplot,cplot)
  
  select case (choice)
  case(XY_VIEW)
     SLL_ASSERT(size(f4d,1) == size(fxy,1))
     SLL_ASSERT(size(f4d,2) == size(fxy,2))
     do j=1,ny
        do i=1,nx
           sumloc= sum(f4d(i,j,:,jstartv:jendv))
           call mpi_reduce(sumloc,fxy(i,j),1,MPI_REAL8,MPI_SUM,0,comm,error)
        end do
     end do
     if (my_num == 0) then
        call sll_hdf5_file_create('fxy'//cplot//".h5",file_id,error)
        call sll_hdf5_write_array(file_id,fxy,"/fxy",error)
        call write_xdmf("fxy"//cplot//".xmf", &
                        'fxy'//cplot//".h5","x","y","fxy",nx,ny)
        call sll_hdf5_file_close(file_id, error)
     end if
  case(XVX_VIEW)
     SLL_ASSERT(size(f4d,1) == size(fxvx,1))
     SLL_ASSERT(size(f4d,3) == size(fxvx,2))
     do k=1,nvx
        do i=1,nx
           sumloc= sum(f4d(i,:,k,jstartv:jendv))
           call mpi_reduce(sumloc,fxvx(i,k),1,MPI_REAL8,MPI_SUM,0,comm,error)
        end do
     end do
     if (my_num == 0) then
        call sll_hdf5_file_create('fxvx'//cplot//".h5",file_id,error)
        call sll_hdf5_write_array(file_id,fxvx,"/fxvx",error)
        call write_xdmf("fxvx"//cplot//".xmf", &
                        'fxvx'//cplot//".h5","x","vx","fxvx",nx,nvx)
        call sll_hdf5_file_close(file_id, error)
     end if
  case(YVY_VIEW)
     SLL_ASSERT(size(f4d,2) == size(fyvy,1))
     SLL_ASSERT(size(f4d,4) == size(fyvy,2))
     do l=jstartv,jendv
     do j=1,ny
        fyvy(j,l)= sum(f4d(:,j,:,l))
     end do
     end do
     call write_fyvy(ny,nvy,cplot,jstartv)
     if (my_num == 0) &
     call write_xdmf("fyvy"//cplot//".xmf", &
                     'fyvy'//cplot//".h5","y","vy","fyvy",ny,nvy)
  case(VXVY_VIEW)
     SLL_ASSERT(size(f4d,3) == size(fvxvy,1))
     SLL_ASSERT(size(f4d,4) == size(fvxvy,2))
     do l=jstartv,jendv
        do k=1,nvx
           fvxvy(k,l)=sum(f4d(:,:,k,l))
        end do
     end do
     call write_fvxvy(nvx,nvy,cplot,jstartv)
     if (my_num == 0) &
     call write_xdmf("fvxvy"//cplot//".xmf",'fvxvy'//cplot//".h5", &
                     "vx","vy","fvxvy",nvx,nvy)
  end select
  
end subroutine plot_df

 subroutine write_fyvy(ny,nvy,cplot,jstartv)

 use hdf5
 use sll_hdf5_io_parallel
 character(len=4)  :: cplot
 integer(HID_T)    :: pfile_id
 integer(HSSIZE_T) :: offset(2)
 integer(HSIZE_T)  :: global_dims(2)
 sll_int32, intent(in)  :: ny,nvy,jstartv
 sll_int32 :: error

 global_dims = (/ny,nvy/)
 offset      = (/0, jstartv-1/)
 call sll_hdf5_file_create('fyvy'//cplot//".h5",MPI_COMM_WORLD,pfile_id,error)
 call sll_hdf5_write_array_2d(pfile_id,global_dims,offset,fyvy,"/fyvy",error)
 call sll_hdf5_file_close(pfile_id, error)

 end subroutine write_fyvy

 subroutine write_fvxvy(nvx,nvy,cplot,jstartv)

 use hdf5
 use sll_hdf5_io_parallel

 character(len=4) :: cplot
 integer(HID_T)    :: pfile_id
 integer(HSSIZE_T) :: offset(2)
 integer(HSIZE_T)  :: global_dims(2)
 sll_int32, intent(in)  :: nvx, nvy, jstartv
 sll_int32 :: error

 global_dims = (/nvx,nvy/)
 offset      = (/0, jstartv-1/)
 call sll_hdf5_file_create('fvxvy'//cplot//".h5",MPI_COMM_WORLD,pfile_id,error)
 call sll_hdf5_write_array_2d(pfile_id,global_dims,offset,fvxvy,"/fvxvy",error)
 call sll_hdf5_file_close(pfile_id, error)

 end subroutine write_fvxvy

 subroutine write_xdmf(xdmffilename, datafilename, xname, yname, fname, nx, ny)

  character(len=*), intent(in) :: xdmffilename
  character(len=*), intent(in) :: datafilename
  character(len=*), intent(in) :: xname
  character(len=*), intent(in) :: yname
  character(len=*), intent(in) :: fname

  sll_int32, intent(in) :: nx
  sll_int32, intent(in) :: ny
  sll_int32 :: file_id
  sll_int32 :: error

  call sll_xml_file_create(xdmffilename,file_id,error)
  write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
  write(file_id, &
   "(a,2i5,a)")"<Topology TopologyType='2DRectMesh' NumberOfElements='",ny,nx,"'/>"
  write(file_id,"(a)")"<Geometry GeometryType='VXVY'>"
  write(file_id,"(a,i5,a)")"<DataItem Dimensions='",nx, &
                           "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")"grid4d.h5:/"//xname
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a,i5,a)")"<DataItem Dimensions='",ny, &
                           "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")"grid4d.h5:/"//yname
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a)")"</Geometry>"
  write(file_id,"(a)")"<Attribute Name='"//fname//"' AttributeType='Scalar' Center='Node'>"
  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ny,nx, &
                            "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")datafilename//":/"//fname
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a)")"</Attribute>"
  !write(file_id,"(a)")"</Grid>"
  call sll_xml_file_close(file_id,error)

 end subroutine write_xdmf

end module vlasov4d_plot
