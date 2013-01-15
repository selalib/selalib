!> @namespace sll_low_level_file_io
!> @author Pierre Navaro
!> @brief 
!> Library to create files readable by visualization softwares.
!> @details
!> External links:
!> - gnuplot http://www.gnuplot.info
!> - VisIt   https://wci.llnl.gov/codes/visit/
!> - Xdmf    http://www.xdmf.org/index.php/Main_Page
!> - HDF5    http://www.hdfgroup.org/HDF5/
!>
!> link with sll_low_level_io_parallel library
!>
!> <b>Plot a 2D field on a curvilinear mesh</b> \n
!> import the module :
!> \code
!> use sll_xdmf
!> \endcode
!> call sequence to create an xdmf file :
!> \code
!> call sll_xdmf_open("file_name","mesh_name",nnodes_x,nnodes_y,file_id,error)
!> call sll_xdmf_write_array("mesh_name",x,'x',error)
!> call sll_xdmf_write_array("mesh_name",y,'y',error)
!> call sll_xdmf_write_array("field_name",df,"NodeVal",error,file_id,"Node")
!> call sll_xdmf_write_array("field",df(1:ncells_x1,1:ncells_x2),"CellVal",error,file_id,"Cell")
!> call sll_xdmf_close(file_id,error)
!> \endcode
!> <b>Plot a 2D field described on a cartesian mesh, 
!> axis are perpendicular and spacing is constant.</b>
!> \code
!> call sll_xdmf_corect2d_nodes( "file_name", df, "field_name", x_min, dx, y_min, dy, "HDF5") 
!> \endcode
!> <b> Plot a 2D field described on a cartesian mesh,
!> axis are perpendicular and spacing is define by x and y (1d arrays)
!> \code
!> call sll_xdmf_rect2d_nodes( "file_name", df, "field_name", x, y, "HDF5") 
!> \endcode
!> <b> Plot a 2D field described on a structured curvilinear mesh,
!> nodes coordinates are defined by x and y (2d arrays)
!> \code
!> call sll_xdmf_curv2d_nodes( "file_name", df, "field_name", x, y, "HDF5") 
!> \endcode
