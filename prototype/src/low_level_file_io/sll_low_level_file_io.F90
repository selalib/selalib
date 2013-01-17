!> @mainpage sll_low_level_file_io
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
!> link with <code>sll_low_level_io</code> library
!>
!> - import the module :
!> \code
!> use sll_xdmf
!> \endcode
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

