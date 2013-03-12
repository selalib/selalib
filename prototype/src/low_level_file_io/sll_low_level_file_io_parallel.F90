!> @namespace sll_low_level_file_io_parallel
!> @author Pierre Navaro
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
!> - Modules available
!>      + sll_hdf5_io_parallel
!>      + sll_xdmf_parallel
!> - 2d example, import the module :
!> \code
!> use sll_xdmf_parallel
!> \endcode
!> - Write some parallel distributed 2d fields on a xmf file, h5file is a string = data file name
!> \code
!> call sll_xdmf_open("fields.xmf","h5file",nx,ny,file_id,error)
!> call sll_xdmf_write_array("h5file",global_dims,offset,x,'x1',error)
!> call sll_xdmf_write_array("h5file",global_dims,offset,y,'x2',error)
!> call sll_xdmf_write_array("h5file",global_dims,offset,rho,"rho",error,file_id,"Node")
!> call sll_xdmf_write_array("h5file",global_dims,offset,phi,"phi",error,file_id,"Node")
!> call sll_xdmf_write_array("h5file",global_dims,offset,ex ,"ex" ,error,file_id,"Node")
!> call sll_xdmf_write_array("h5file",global_dims,offset,ey ,"ey" ,error,file_id,"Node")
!> call sll_xdmf_close(file_id,error)
!> \endcode
!> - nx and ny are local sizes of fields
!> - global_dims : array with global size of the distributed array
!> type declaration example in 2d case :
!> \code
!> sll_int32  :: file_id
!> integer(HSSIZE_T) :: offset(2)
!> integer(HSIZE_T)  :: global_dims(2)
!> \endcode

