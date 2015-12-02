!> @ingroup fields
module sll_parallel_array_output_module

#include "sll_working_precision.h"
#include "sll_memory.h"

  use sll_m_collective
  use sll_m_remapper
  use sll_m_xml_io
  use sll_m_cartesian_meshes
  use sll_m_utilities, only : &
       int2string

#ifndef NOHDF5
  use hdf5, only: hid_t, hssize_t, hsize_t
  use sll_m_hdf5_io_serial
  use sll_m_hdf5_io_parallel
#endif

#define MPI_MASTER 0

  implicit none

contains

subroutine write_mesh_4d(mesh)

  type(sll_cartesian_mesh_4d), intent(in) :: mesh
  sll_int32                               :: error
  sll_int32                               :: file_id
  sll_real64, dimension(:), allocatable   :: eta1
  sll_real64, dimension(:), allocatable   :: eta2
  sll_real64, dimension(:), allocatable   :: eta3
  sll_real64, dimension(:), allocatable   :: eta4

  sll_int32 :: i, j, k, l

  error = 0
 
  SLL_ALLOCATE(eta1(mesh%num_cells1+1),error)
  SLL_ALLOCATE(eta2(mesh%num_cells2+1),error)
  SLL_ALLOCATE(eta3(mesh%num_cells3+1),error)
  SLL_ALLOCATE(eta4(mesh%num_cells4+1),error)
 
  do i = 1, mesh%num_cells1+1
     eta1(i) = mesh%eta1_min + real(i-1,f64)*mesh%delta_eta1
  end do
  do j = 1, mesh%num_cells2+1
     eta2(j) = mesh%eta2_min + real(j-1,f64)*mesh%delta_eta2
  end do
  do k = 1, mesh%num_cells3+1
     eta3(k) = mesh%eta3_min + real(k-1,f64)*mesh%delta_eta3
  end do
  do l = 1, mesh%num_cells4+1
     eta4(l) = mesh%eta4_min + real(l-1,f64)*mesh%delta_eta4
  end do
 
#ifndef NOHDF5
  call sll_hdf5_file_create("mesh4d.h5",file_id,error)
  call sll_hdf5_write_array(file_id,eta1,"/x1",error)
  call sll_hdf5_write_array(file_id,eta2,"/x2",error)
  call sll_hdf5_write_array(file_id,eta3,"/x3",error)
  call sll_hdf5_write_array(file_id,eta4,"/x4",error)
  call sll_hdf5_file_close(file_id, error)
#endif

end subroutine write_mesh_4d

subroutine write_xmf_file(mesh, f, layout, iplot)

  class(sll_cartesian_mesh_4d),   intent(in) :: mesh
  sll_real64, dimension(:,:,:,:), intent(in) :: f
  type(layout_4d), pointer,       intent(in) :: layout
  sll_int32, intent(in)                      :: iplot
  sll_int32                                  :: error
  character(len=4)                           :: cplot
  sll_int32                                  :: prank
  sll_int32, parameter                       :: one = 1
  sll_int32                                  :: file_id
  sll_int32                                  :: nx1
  sll_int32                                  :: nx2
  sll_int32                                  :: nx3
  sll_int32                                  :: nx4

  call int2string(iplot,cplot)
  call write_fx1x2(f, layout, cplot)
  call write_fx1x3(f, layout, cplot)
  call write_fx2x4(f, layout, cplot)
  call write_fx3x4(f, layout, cplot)

  prank = sll_get_collective_rank(sll_world_collective)
  if ( prank == 0) then

     nx1 = mesh%num_cells1+1
     nx2 = mesh%num_cells2+1
     nx3 = mesh%num_cells3+1
     nx4 = mesh%num_cells4+1

     call sll_xml_file_create("fvalues_"//cplot//".xmf",file_id,error)
     call write_grid(file_id,nx1,nx2,"x1","x2",cplot)
     call write_grid(file_id,nx1,nx3,"x1","x3",cplot)
     call write_grid(file_id,nx2,nx4,"x2","x4",cplot)
     call write_grid(file_id,nx3,nx4,"x3","x4",cplot)
     write(file_id,"(a)")"</Domain>"
     write(file_id,"(a)")"</Xdmf>"
     close(file_id)

  endif

end subroutine write_xmf_file

subroutine write_grid(file_id,nx,ny,xname,yname,cplot)

  sll_int32                                :: file_id
  sll_int32                                :: nx
  sll_int32                                :: ny
  character(len=*)                         :: cplot
  character(len=*)                         :: xname
  character(len=*)                         :: yname

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
  call write_attribute(file_id,nx,ny,"f",cplot,xname,yname)
  write(file_id,"(a)")"</Grid>"

end subroutine write_grid

subroutine write_attribute(file_id,nx,ny,fname,cplot,xname,yname)

  sll_int32                    :: file_id
  sll_int32                    :: nx
  sll_int32                    :: ny
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: cplot
  character(len=*), optional   :: xname
  character(len=*), optional   :: yname

  write(file_id,"(a)") &
     "<Attribute Name='"//fname//"' AttributeType='Scalar' Center='Node'>"
  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ny,nx, &
                            "' NumberType='Float' Precision='8' Format='HDF'>"
  if( present(xname) .and. present(yname)) then
     write(file_id,"(a)")fname//xname//yname//"_"//cplot//".h5:/values"
  else
     write(file_id,"(a)")fname//"_"//cplot//".h5:/values"
  end if
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a)")"</Attribute>"

end subroutine write_attribute

subroutine write_fx1x2(f, layout, cplot)

  sll_real64, dimension(:,:,:,:)          :: f
  type(layout_4D), pointer                :: layout
  character(len=*)                        :: cplot
  sll_int32                               :: error
  sll_int32                               :: file_id
  sll_int32                               :: prank
  sll_int32                               :: comm
  sll_real64, dimension(:,:), pointer     :: fij
  sll_real64                              :: sumloc
 
  sll_int32 :: i, j
  sll_int32 :: loc_sz_i, loc_sz_j, loc_sz_k, loc_sz_l

  prank = sll_get_collective_rank(sll_world_collective)
  comm  = sll_world_collective%comm

  call compute_local_sizes(layout,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(fij(1:loc_sz_i,1:loc_sz_j),error)

  do j=1,loc_sz_j
     do i=1,loc_sz_i
        sumloc = sum(f(i,j,:,:))
        call mpi_reduce(sumloc,fij(i,j),1,MPI_REAL8,MPI_SUM,MPI_MASTER,comm,error)
     end do
  end do
#ifndef NOHDF5
  if (prank == MPI_MASTER) then
     call sll_hdf5_file_create('fx1x2_'//cplot//".h5",file_id,error)
     call sll_hdf5_write_array(file_id,fij,"/values",error)
     call sll_hdf5_file_close(file_id, error)
  end if
#endif

end subroutine write_fx1x2

subroutine write_fx1x3(f, layout, cplot)

  sll_real64, dimension(:,:,:,:)          :: f
  type(layout_4D), pointer                :: layout
  character(len=*)                        :: cplot
  sll_int32                               :: error
  sll_int32                               :: file_id
  sll_real64, dimension(:,:), pointer     :: fik
  sll_int32                               :: prank
  sll_int32                               :: comm
  sll_real64                              :: sumloc

  sll_int32 :: i, k
  sll_int32 :: loc_sz_i, loc_sz_j, loc_sz_k, loc_sz_l

  prank = sll_get_collective_rank(sll_world_collective)
  comm  = sll_world_collective%comm

  call compute_local_sizes(layout,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(fik(1:loc_sz_i,1:loc_sz_k),error)

  do k=1,loc_sz_k
     do i=1,loc_sz_i
        sumloc= sum(f(i,:,k,:))
        call mpi_reduce(sumloc,fik(i,k),1,MPI_REAL8,MPI_SUM,MPI_MASTER,comm,error)
     end do
  end do
#ifndef NOHDF5
  if (prank == MPI_MASTER) then
     call sll_hdf5_file_create('fx1x3_'//cplot//".h5",file_id,error)
     call sll_hdf5_write_array(file_id,fik,"/values",error)
     call sll_hdf5_file_close(file_id, error)
  end if
#endif

end subroutine write_fx1x3

subroutine write_fx2x4(f, layout, cplot)

  sll_real64, dimension(:,:,:,:)          :: f
  type(layout_4D), pointer                :: layout
  character(len=*)                        :: cplot
  integer(HID_T)                          :: pfile_id
  integer(HSSIZE_T)                       :: offset(2)
  integer(HSIZE_T)                        :: global_dims(2)
  sll_int32                               :: error
  sll_int32                               :: prank
  sll_int32                               :: comm
  sll_real64, dimension(:,:), pointer     :: fjl

  sll_int32 :: j, l
  sll_int32 :: loc_sz_i, loc_sz_j, loc_sz_k, loc_sz_l

  prank = sll_get_collective_rank(sll_world_collective)
  comm  = sll_world_collective%comm

  call compute_local_sizes(layout,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(fjl(1:loc_sz_j,1:loc_sz_l),error)
  do l=1,loc_sz_l
     do j=1,loc_sz_j
        fjl(j,l) = sum(f(:,j,:,l))
     end do
  end do

#ifdef NOHDF5
#define HSIZE_T  i32
#define HSSIZE_T i32
#endif

  global_dims(1) = int(get_layout_global_size_j(layout),HSIZE_T)
  global_dims(2) = int(get_layout_global_size_l(layout),HSIZE_T)
  offset(1) = int(get_layout_j_min(layout,prank)-1,HSSIZE_T)
  offset(2) = int(get_layout_l_min(layout,prank)-1,HSSIZE_T)

#ifndef NOHDF5
  call sll_hdf5_file_create('fx2x4_'//cplot//".h5",comm,pfile_id,error)
  call sll_hdf5_write_array(pfile_id,global_dims,offset,fjl,"/values",error)
  call sll_hdf5_file_close(pfile_id, error)
#endif

end subroutine write_fx2x4

subroutine write_fx3x4(f, layout, cplot)

  sll_real64, dimension(:,:,:,:)          :: f
  Type(layout_4D), pointer                :: layout
  character(len=*)                        :: cplot
  integer(HID_T)                          :: pfile_id
  integer(HSSIZE_T)                       :: offset(2)
  integer(HSIZE_T)                        :: global_dims(2)
  sll_int32                               :: error
  sll_int32                               :: prank
  sll_int32                               :: comm
  sll_real64, dimension(:,:), pointer     :: fkl

  sll_int32 :: k, l
  sll_int32 :: loc_sz_i, loc_sz_j, loc_sz_k, loc_sz_l

  prank = sll_get_collective_rank(sll_world_collective)
  comm  = sll_world_collective%comm

  call compute_local_sizes(layout,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(fkl(1:loc_sz_k,1:loc_sz_l),error)
  do l=1,loc_sz_l
     do k=1,loc_sz_k
        fkl(k,l) = sum(f(:,:,k,l))
     end do
  end do
  global_dims(1) = int(get_layout_global_size_k(layout), HSIZE_T)
  global_dims(2) = int(get_layout_global_size_l(layout), HSIZE_T)
  offset(1) = int(get_layout_k_min(layout,prank)-1,i64)
  offset(2) = int(get_layout_l_min(layout,prank)-1,i64)

#ifndef NOHDF5
  call sll_hdf5_file_create('fx3x4_'//cplot//".h5",comm,pfile_id,error)
  call sll_hdf5_write_array(pfile_id,global_dims,offset,fkl,"/values",error)
  call sll_hdf5_file_close(pfile_id, error)
#endif

end subroutine write_fx3x4

end module sll_parallel_array_output_module
