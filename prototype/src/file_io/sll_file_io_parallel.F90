!> @namespace sll_file_io_parallel
!> @brief 
!> Library to create parallel files readable by visualization softwares.
!> @details
!> - External links:
!>      - VisIt   https://wci.llnl.gov/codes/visit/
!>      - Xdmf    http://www.xdmf.org/index.php/Main_Page
!>      - HDF5    http://www.hdfgroup.org/HDF5/
!>
!> link with <code>sll_low_level_io_parallel</code> library
!>
!> - <b> Modules available </b>
!>      + sll_hdf5_io_parallel
!>      + sll_xdmf_parallel
!>      + sll_gnuplot_parallel
!>
!> - <b> How to use it : </b>
!>      + Import the module :
!>        @code #include "sll_file_io_parallel.h" @endcode
!>      + Link with the library
!>        @code sll_file_io_parallel @endcode
!> - <b> Examples </b>
!>      + Variables and array initialization from a layout_2d
!>
!>        \code
!>        sll_int32   :: nx, ny          !local sizes of fields
!>        sll_int32   :: file_id         !file unit number of xmf file
!>        integer(HSSIZE_T) :: offset(2) !array offset 
!>        integer(HSIZE_T)  :: global_dims(2) !global sizes
!>        SLL_ALLOCATE(xdata(nx,ny),error)
!>        SLL_ALLOCATE(ydata(nx,ny),error)
!>        SLL_ALLOCATE(zdata(nx,ny),error)
!>  
!>        do j = 1, my
!>           do i = 1, mx
!>             global_indices =  local_to_global_2D( layout, (/i, j/) )
!>             gi = global_indices(1)
!>             gj = global_indices(2)
!>             xdata(i,j) = float(gi-1)/(nx-1)
!>             ydata(i,j) = float(gj-1)/(ny-1)
!>             zdata(i,j) = (myrank+1) * xdata(i,j) * ydata(i,j)
!>           end do
!>        end do
!>  
!>        offset(1) =  get_layout_2D_i_min( layout, myrank ) - 1
!>        offset(2) =  get_layout_2D_j_min( layout, myrank ) - 1
!>        \endcode
!>
!>      + xdmf : Write some parallel distributed 2d 
!> node centered fields in a xmf file, h5file is the hdf5 data file name.
!> \code
!> call sll_xdmf_open("fields.xmf",h5file,nx,ny,file_id,error)
!> call sll_xdmf_write_array(h5file,global_dims,offset,x,'x1',error)
!> call sll_xdmf_write_array(h5file,global_dims,offset,y,'x2',error)
!> call sll_xdmf_write_array(h5file,global_dims,offset,rho,"rho",error,file_id,"Node")
!> call sll_xdmf_write_array(h5file,global_dims,offset,phi,"phi",error,file_id,"Node")
!> call sll_xdmf_write_array(h5file,global_dims,offset,ex ,"ex" ,error,file_id,"Node")
!> call sll_xdmf_write_array(h5file,global_dims,offset,ey ,"ey" ,error,file_id,"Node")
!> call sll_xdmf_close(file_id,error)
!> \endcode
!>
!>      + gnuplot : write data on a mesh
!> @code
!>  call sll_gnuplot_rect_2d_parallel(dble(offset(1)), dble(1), &
!>                                    dble(offset(2)), dble(1), &
!>                                    zdata, "rect_mesh", 1, error)  
!>
!>  call sll_gnuplot_curv_2d_parallel(xdata, ydata, zdata, "curv_mesh", 1, error)  
!> @endcode 
!>  
!>      + hdf5 : write data on a parallel hdf5 file
!> @code
!>  call sll_hdf5_file_create(xfile, file_id, error)
!>  call sll_hdf5_write_array(file_id, datadims,offset,xdata,xdset,error)
!>  call sll_hdf5_file_close(file_id,error)
!> @endcode 
