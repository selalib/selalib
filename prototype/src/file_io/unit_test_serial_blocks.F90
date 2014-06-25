program test_serial_blocks
#include "sll_working_precision.h"
#include "sll_utilities.h"
use mpi
use hdf5
use sll_hdf5_io_serial
use sll_xdmf_serial_blocks
implicit none

sll_int32 :: i, j
sll_int32 :: nx, ny, sx, ex, sy, ey
sll_int32 :: error
sll_real64 :: dimx, dimy
sll_real64 :: hx, hy
sll_real64, dimension(:,:), allocatable :: x 
sll_real64, dimension(:,:), allocatable :: y
sll_real64, dimension(:,:), allocatable :: z

sll_int32                      :: prank
sll_int32                      :: psize
sll_int32                      :: code
sll_int32                      :: comm2D
sll_int32, parameter           :: tag = 1111
sll_int32, parameter           :: ndims=2
sll_int32, dimension(ndims)    :: dims,coords
logical, dimension(ndims)      :: periods
sll_int32, dimension(8)        :: voisin
sll_int32, parameter           :: N =1, S =2, W =3, E =4
sll_real64                     :: tcpu1, tcpu2
sll_int32                      :: nxp, nyp
logical                        :: reorder
sll_int32, parameter           :: xmf = 77
character(len=72)              :: field_label, mesh_label
character(len=72)              :: field_name, mesh_name
character(len=4)               :: my_proc, cproc
character(len=2), dimension(2) :: coordNames

sll_int32         :: iproc
integer(hid_t)    :: file_id
integer(hsize_t)  :: data_dims(2)
sll_int32         :: iplot = 0
character(len=4)  :: cplot
character(len=22) :: file_name = "test_serial_blocks.xmf"
sll_int32         :: xmf_id


dimx = 2.0_f64
nx   = 128
dimy = 1.0_f64
ny   = 64

hx = dimx / (nx-1)
hy = dimy / (ny-1)

!Initialisation de MPI

call MPI_INIT(code)
call MPI_COMM_SIZE(MPI_COMM_WORLD,psize,code)
call MPI_COMM_RANK(MPI_COMM_WORLD,prank,code)
call MPI_BARRIER(MPI_COMM_WORLD,code)

tcpu1 = MPI_WTIME()

!Nombre de processus suivant x et y
dims(:) = 0
CALL MPI_DIMS_CREATE(psize,ndims,dims,code)
nxp = dims(1)
nyp = dims(2)
if (prank == 0) then

   print*," nxp, nyp = ", nxp, nyp

end if

!Creation de la grille 2D periodique en x et y
periods(1) = .false.
periods(2) = .false.
reorder    = .true.

CALL MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,periods,reorder,comm2d,code)

!Initialisation du tableau voisin 
voisin(:) = MPI_PROC_NULL

!Recherche de mes voisins Sud et Nord
CALL MPI_CART_SHIFT(comm2d,0,1,voisin(N),voisin(S),code)
!Recherche de mes voisins Ouest et Est
CALL MPI_CART_SHIFT(comm2d,1,1,voisin(W),voisin(E),code)
!Connaitre mes coordonnees dans la topologie
CALL MPI_COMM_RANK(comm2d,prank,code)
CALL MPI_CART_COORDS(comm2d,prank,ndims,coords,code)

if (prank == 0) then
   write(*,"(a,g12.3)")" largeur dimx          = ", dimx
   write(*,"(a,g12.3)")" longueur dimy         = ", dimy
   write(*,"(a,g12.3)")" nombre nx             = ", nx
   write(*,"(a,g12.3)")" nombre ny             = ", ny
end if

call mpe_decomp1d(nx,psize,prank,sx,ex)
call mpe_decomp1d(ny,psize,prank,sy,ey)

allocate(x(sx:ex,sy:ey))
allocate(y(sx:ex,sy:ey))
allocate(z(sx:ex,sy:ey))

do j=sy,ey
   do i=sx,ex
      x(i,j)=i*hx
      y(i,j)=j*hy
   end do
end do

z = prank

call int2string(prank,my_proc)
field_label = "xdmf2d"
mesh_label  = "mesh2d"

!Write the data file.
coordNames(1) = "/X"
coordNames(2) = "/Y"
    
!Write separate coordinate arrays for the x and y coordinates.
nx = ex-sx+1
ny = ey-sy+1
data_dims = (/nx,ny/)

iplot = iplot+1
call int2string(iplot,cplot)
mesh_name  = trim(mesh_label)//my_proc//".h5"
field_name = trim(field_label)//my_proc//"-"//cplot//".h5"

call sll_hdf5_file_create(trim(mesh_name),file_id,error)
call sll_hdf5_write_array(file_id,x(sx:ex,sy:ey),coordnames(1),error)
call sll_hdf5_write_array(file_id,y(sx:ex,sy:ey),coordnames(2),error)
call sll_hdf5_file_close(file_id, error)

call sll_hdf5_file_create(trim(field_name),file_id,error)
call sll_hdf5_write_array(file_id,z(sx:ex,sy:ey),"/Z",error)
call sll_hdf5_file_close(file_id, error)

if (prank == 0) then
   open(xmf,file="all_domains"//cplot//".xmf")
   write(xmf,'(a)')"<?xml version=""1.0"" ?>"
   write(xmf,'(a)')"<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
   write(xmf,'(a)')"<Xdmf Version=""2.0"">"
   write(xmf,'(a)')"<Domain Name=""mesh2d"">"
   write(xmf,'(a)')  &
   "<Grid Name=""Domain"" GridType=""Collection"" CollectionType=""Spatial"">"
   do iproc = 0, psize-1
      call int2string(iproc,cproc)
      write(xmf,'(a)')"<Grid Name=""SubDomain"" GridType=""Uniform"">"
      write(xmf,'(a,2i6,a)') &
      "<Topology TopologyType=""2DSMesh"" NumberOfElements='",ny,nx,"'/>"
      write(xmf,'(a)')"<Geometry GeometryType=""X_Y"">"
      write(xmf,'(a,2i6,a)')"<DataItem Dimensions='",ny,nx, &
      "' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
      write(xmf,'(a)')trim(mesh_label)//cproc//".h5:"//coordnames(1)
      write(xmf,'(a)')"</DataItem>"
      write(xmf,'(a,2i6,a)')"<DataItem Dimensions='",ny,nx, &
      "' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
      write(xmf,'(a)')trim(mesh_label)//cproc//".h5:"//coordnames(2)
      write(xmf,'(a)')"</DataItem>"
      write(xmf,'(a)')"</Geometry>"
      write(xmf,'(a)') &
      "<Attribute Name=""Z"" AttributeType=""Scalar"" Center=""Node"">"
      write(xmf,'(a,2i6,a)')"<DataItem Dimensions='",ny,nx, &
      "' NumberType=""Float"" Precision=""8"" Format=""HDF"">"
      write(xmf,'(a)')trim(field_label)//cproc//"-"//cplot//".h5:/Z"
      write(xmf,'(a)')"</DataItem>"
      write(xmf,'(a)')"</Attribute>"
      write(xmf,'(a)')"</Grid>"
   end do
   write(xmf,'(a)')"</Grid>"
   write(xmf,'(a)')"</Domain>"
   write(xmf,'(a)')"</Xdmf>"
   close(xmf)
end if

call sll_xdmf_open_serial_blocks(file_name, xmf_id, error)
tcpu2 = MPI_WTIME()

print*, tcpu2-tcpu1

call MPI_FINALIZE(code)

end program test_serial_blocks
