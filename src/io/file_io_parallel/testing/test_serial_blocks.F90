program test_serial_blocks
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use sll_m_hdf5_io_serial, only: &
      sll_t_hdf5_ser_handle, &
      sll_s_hdf5_ser_file_create, &
      sll_s_hdf5_ser_file_close, &
      sll_o_hdf5_ser_write_array

   use sll_m_utilities, only: &
      sll_s_int2string, &
      sll_s_mpe_decomp1d

   use sll_m_xdmf_serial_blocks, only: &
      sll_s_xdmf_array_2d_serial_blocks, &
      sll_s_xdmf_close_serial_blocks, &
      sll_s_xdmf_open_serial_blocks

   use mpi, only: &
      mpi_barrier, &
      mpi_cart_coords, &
      mpi_cart_create, &
      mpi_cart_shift, &
      mpi_comm_rank, &
      mpi_comm_size, &
      mpi_comm_world, &
      mpi_dims_create, &
      mpi_finalize, &
      mpi_init, &
      mpi_proc_null, &
      mpi_wtime

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_int32 :: i, j
   sll_int32 :: nx, ny, sx, ex, sy, ey
   sll_int32 :: error
   sll_real64 :: dimx, dimy
   sll_real64 :: hx, hy
   sll_real64, dimension(:, :), allocatable :: x
   sll_real64, dimension(:, :), allocatable :: y
   sll_real64, dimension(:, :), allocatable :: z

   sll_int32                      :: prank
   sll_int32                      :: psize
   sll_int32                      :: code
   sll_int32                      :: comm2D
   sll_int32, parameter           :: tag = 1111
   sll_int32, parameter           :: ndims = 2
   sll_int32, dimension(ndims)    :: dims, coords
   logical, dimension(ndims)      :: periods
   sll_int32, dimension(8)        :: voisin
   sll_int32, parameter           :: N = 1, S = 2, W = 3, E = 4
   sll_real64                     :: tcpu1, tcpu2
   sll_int32                      :: nxp, nyp
   logical                        :: reorder
   sll_int32, parameter           :: xmf = 77
   character(len=72)              :: field_label, mesh_label
   character(len=72)              :: field_name, mesh_name
   character(len=4)               :: my_proc, cproc
   character(len=2), dimension(2) :: coordNames

   sll_int32                      :: iproc
   type(sll_t_hdf5_ser_handle)    :: hfile_id
   sll_int32                      :: iplot = 0
   character(len=4)               :: cplot
   sll_int32                      :: xmf_id

   dimx = 2.0_f64
   nx = 127
   dimy = 1.0_f64
   ny = 64

   hx = dimx/(nx - 1)
   hy = dimy/(ny - 1)

!Initialisation de MPI

   call MPI_INIT(code)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, psize, code)
   call MPI_COMM_RANK(MPI_COMM_WORLD, prank, code)
   call MPI_BARRIER(MPI_COMM_WORLD, code)

   tcpu1 = MPI_WTIME()

!Nombre de processus suivant x et y
   dims(:) = 0
   CALL MPI_DIMS_CREATE(psize, ndims, dims, code)
   nxp = dims(1)
   nyp = dims(2)
   if (prank == 0) then
      print *, " nxp, nyp = ", nxp, nyp
   end if

!Creation de la grille 2D periodique en x et y
   periods(1) = .false.
   periods(2) = .false.
   reorder = .true.

   CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, comm2d, code)

!Initialisation du tableau voisin
   voisin(:) = MPI_PROC_NULL

!Recherche de mes voisins Sud et Nord
   CALL MPI_CART_SHIFT(comm2d, 0, 1, voisin(N), voisin(S), code)
!Recherche de mes voisins Ouest et Est
   CALL MPI_CART_SHIFT(comm2d, 1, 1, voisin(W), voisin(E), code)
!Connaitre mes coordonnees dans la topologie
   CALL MPI_COMM_RANK(comm2d, prank, code)
   CALL MPI_CART_COORDS(comm2d, prank, ndims, coords, code)

   if (prank == 0) then
      write (*, "(a,g12.3)") " largeur dimx          = ", dimx
      write (*, "(a,g12.3)") " longueur dimy         = ", dimy
      write (*, "(a,g12.3)") " nombre nx             = ", nx
      write (*, "(a,g12.3)") " nombre ny             = ", ny
   end if

   call sll_s_mpe_decomp1d(nx, dims(1), coords(1), sx, ex)
   call sll_s_mpe_decomp1d(ny, dims(2), coords(2), sy, ey)

   do iproc = 1, psize
      if (prank == iproc - 1) then
         print"('Rank ',i3,'[',2i5,'][',2i5,']')", prank, sx, ex, sy, ey
      end if
      call MPI_Barrier(MPI_COMM_WORLD, error)
   end do

   allocate (x(sx:ex, sy:ey))
   allocate (y(sx:ex, sy:ey))
   allocate (z(sx:ex, sy:ey))

   do j = sy, ey
      do i = sx, ex
         x(i, j) = i*hx
         y(i, j) = j*hy
      end do
   end do

   z = real(prank, f64)

   call sll_s_int2string(prank, my_proc)
   field_label = "xdmf2d"
   mesh_label = "mesh2d"

!Write the data file.
   coordNames(1) = "/X"
   coordNames(2) = "/Y"

!Write separate coordinate arrays for the x and y coordinates.
   nx = ex - sx + 1
   ny = ey - sy + 1

   iplot = iplot + 1
   call sll_s_int2string(iplot, cplot)
   mesh_name = trim(mesh_label)//my_proc//".h5"
   field_name = trim(field_label)//my_proc//"-"//cplot//".h5"

   call sll_s_hdf5_ser_file_create(trim(mesh_name), hfile_id, error)
   call sll_o_hdf5_ser_write_array(hfile_id, x(sx:ex, sy:ey), coordnames(1), error)
   call sll_o_hdf5_ser_write_array(hfile_id, y(sx:ex, sy:ey), coordnames(2), error)
   call sll_s_hdf5_ser_file_close(hfile_id, error)

   call sll_s_hdf5_ser_file_create(trim(field_name), hfile_id, error)
   call sll_o_hdf5_ser_write_array(hfile_id, z(sx:ex, sy:ey), "/Z", error)
   call sll_s_hdf5_ser_file_close(hfile_id, error)

   if (prank == 0) then
      open (xmf, file="all_domains"//cplot//".xmf")
      write (xmf, '(a)') "<?xml version=""1.0"" ?>"
      write (xmf, '(a)') "<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
      write (xmf, '(a)') "<Xdmf Version=""2.0"">"
      write (xmf, '(a)') "<Domain Name=""mesh2d"">"
      write (xmf, '(a)') &
         "<Grid Name=""Domain"" GridType=""Collection"" CollectionType=""Spatial"">"
      do iproc = 0, psize - 1
         call sll_s_int2string(iproc, cproc)
         write (xmf, '(a)') "<Grid Name=""SubDomain"" GridType=""Uniform"">"
         write (xmf, '(a,2i6,a)') &
            "<Topology TopologyType=""2DSMesh"" NumberOfElements='", ny, nx, "'/>"
         write (xmf, '(a)') "<Geometry GeometryType=""X_Y"">"
         write (xmf, '(a,2i6,a)') "<DataItem Dimensions='", ny, nx, &
            "' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
         write (xmf, '(a)') trim(mesh_label)//cproc//".h5:"//coordnames(1)
         write (xmf, '(a)') "</DataItem>"
         write (xmf, '(a,2i6,a)') "<DataItem Dimensions='", ny, nx, &
            "' NumberType=""Float"" Precision=""4"" Format=""HDF"">"
         write (xmf, '(a)') trim(mesh_label)//cproc//".h5:"//coordnames(2)
         write (xmf, '(a)') "</DataItem>"
         write (xmf, '(a)') "</Geometry>"
         write (xmf, '(a)') &
            "<Attribute Name=""Z"" AttributeType=""Scalar"" Center=""Node"">"
         write (xmf, '(a,2i6,a)') "<DataItem Dimensions='", ny, nx, &
            "' NumberType=""Float"" Precision=""8"" Format=""HDF"">"
         write (xmf, '(a)') trim(field_label)//cproc//"-"//cplot//".h5:/Z"
         write (xmf, '(a)') "</DataItem>"
         write (xmf, '(a)') "</Attribute>"
         write (xmf, '(a)') "</Grid>"
      end do
      write (xmf, '(a)') "</Grid>"
      write (xmf, '(a)') "</Domain>"
      write (xmf, '(a)') "</Xdmf>"
      close (xmf)
   end if

   call sll_s_xdmf_open_serial_blocks("data_serial_blocks", x, y, xmf_id, error)
   call sll_s_xdmf_array_2d_serial_blocks(xmf_id, "data_serial_blocks", z, "z_values", error)
   call sll_s_xdmf_close_serial_blocks(xmf_id, error)

   tcpu2 = MPI_WTIME()

   print *, tcpu2 - tcpu1

   call MPI_FINALIZE(code)

   print *, 'PASSED'
end program test_serial_blocks
