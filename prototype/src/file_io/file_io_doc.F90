!> @namespace sll_file_io
!> @brief 
!> Library to create files readable by visualization softwares.
!> @details
!> External links:
!> <table>
!> <tr><td> gnuplot </td><td>http://www.gnuplot.info                 </td></tr>
!> <tr><td> VisIt   </td><td>https://wci.llnl.gov/codes/visit/       </td></tr>
!> <tr><td> Xdmf    </td><td>http://www.xdmf.org/index.php/Main_Page </td></tr>
!> <tr><td> HDF5    </td><td>http://www.hdfgroup.org/HDF5/           </td></tr>
!> <tr><td> Plotmtv </td><td>http://www.phy.ornl.gov/csep/CSEP/CORNELL/TUTORIAL/PLOTMTV/OVERVIEW.html           </td></tr>
!> </table>
!>
!> <b> Modules available </b>
!>
!> - sll_xdmf
!> - sll_ascii_io
!> - sll_gnuplot
!> - sll_binary_io
!> - sll_xml_io
!> - sll_hdf5_io
!> - sll_plotmtv
!>
!> <b> How to use it </b>
!> - Header file : \code #include "sll_file_io.h" \endcode
!> - Link with   <code>-lsll_file_io</code> library
!>
!> <b> Examples </b>
!> - Plot a 2D field described on a cartesian mesh, axis are perpendicular and spacing is constant.
!> \code
!> call sll_xdmf_corect2d_nodes("file_name",df,"field_name",x_min,dx,y_min,dy,"HDF5")
!> \endcode
!> - Plot a 2D field described on a cartesian mesh, axis are perpendicular and spacing is define by x and y (1d arrays).
!> \code
!> call sll_xdmf_rect2d_nodes("file_name",df,"field_name",x,y,"HDF5") 
!> \endcode
!> - Plot a 2D field described on a structured curvilinear mesh, nodes coordinates are defined by x and y (2d arrays).
!> \code
!> call sll_xdmf_curv2d_nodes("file_name",df,"field_name",x,y,"HDF5") 
!> \endcode
!>
!> - If you want to plot several attributes then call sequence to 
!>   create an xdmf file is :
!> \code
!> call sll_xdmf_open("file_name.xmf","mesh_name",nnodes_x,nnodes_y, &
!!                    file_id,error)
!> call sll_xdmf_write_array("mesh_name",x,'x',error)
!> call sll_xdmf_write_array("mesh_name",y,'y',error)
!> call sll_xdmf_write_array("field_name",df,"NodeVal",error,file_id,"Node")
!> call sll_xdmf_write_array("field",df(1:ncells_x1,1:ncells_x2),"CellVal",error,file_id,"Cell")
!> call sll_xdmf_close(file_id,error)
!> \endcode
!> - To save some disk space, mesh coordinates can be written only at the initial step
!> of the simulation :
!> \code
!! do istep = 1, nstep
!!
!!    if (istep == 1) then
!!       call sll_hdf5_file_create("mesh_name-x1.h5",file_id,error)
!!       call sll_hdf5_write_array(file_id,x1,"/x1",error)
!!       call sll_hdf5_file_close(file_id, error)
!!       call sll_hdf5_file_create("mesh_name-x2.h5",file_id,error)
!!       call sll_hdf5_write_array(file_id,x2,"/x2",error)
!!       call sll_hdf5_file_close(file_id, error)
!!    end if
!! 
!!    call int2string(iplot,cplot)
!!    call sll_xdmf_open("field_name"//cplot//".xmf","mesh_name", &
!!                    nnodes_x1,nnodes_x2,file_id,error)
!!    call sll_xdmf_write_array("f"//cplot,f,"f_values", &
!!                    error,file_id,"Node")
!!    call sll_xdmf_close(file_id,error)
!!
!! end do
!> \endcode
!> - Write an array in a HDF5 file
!> \code
!> call sll_hdf5_file_create("file_name.h5",file_id,error)
!> call sll_hdf5_write_array(file_id,array,"dataset_name",error)
!> call sll_hdf5_file_close(file_id, error)
!> \endcode

