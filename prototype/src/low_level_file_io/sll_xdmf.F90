!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_low_level_file_io
!
!> @author
!> Pierre Navaro
!>
!
! DESCRIPTION: 
!
!> @brief
!> Implements the functions to write xdmf file plotable by VisIt
!>
!>@details
!> In <b> XDMF </b> (eXtensible Data Model and Format) the description of the 
!> data is separate from the values themselves. Light data is stored using XML, 
!> Heavy data is stored using HDF5 or Binary files. 
!>
!> This is control by the variable <code>NOHDF5</code>.
!> HDF5 is set by default but il you prefer binary just add 
!>
!> <code> env.Append(CPPDEFINES=['NOHDF5']) </code>
!>
!> in your SCons script
!>
!> <h2>How to use this module: </h2>
!>
!> <p> Just use the module \a sll_low_level_file_io 
!> \code use sll_low_level_file_io \endcode
!>
!> External links:
!> - https://wci.llnl.gov/codes/visit/
!> - http://www.xdmf.org/index.php/Main_Page
!> - HDF5 file (http://www.hdfgroup.org/HDF5/)
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_xdmf
#include "sll_working_precision.h"
#include "sll_assert.h"
  
  use sll_hdf5_io
  use sll_binary_io
  use sll_ascii_io
  use sll_xml_io
  
  implicit none
  
  interface sll_xdmf_open
     module procedure sll_xdmf_open_2d
     module procedure sll_xdmf_open_3d
  end interface

  interface sll_xdmf_write_array
     module procedure sll_xdmf_array_2d
     module procedure sll_xdmf_array_3d
  end interface
  
  interface sll_xdmf_close
     module procedure sll_xml_file_close
  end interface
  
contains  
  
  subroutine sll_xdmf_open_2d(prefix,file_id,nnodes_x1,nnodes_x2,error)
    character(len=*), intent(in) :: prefix
    sll_int32, intent(out)       :: file_id
    sll_int32, intent(out)       :: error
    sll_int32                    :: nnodes_x1
    sll_int32                    :: nnodes_x2
    
    call sll_xml_file_create(trim(prefix)//".xmf",file_id,error)
    call sll_xml_grid_geometry(file_id, trim(prefix), nnodes_x1, nnodes_x2)
  end subroutine sll_xdmf_open_2d
  
  subroutine sll_xdmf_array_2d(prefix,array,array_name,error,xmffile_id,center)
    character(len=*), intent(in) :: prefix
    sll_real64, intent(in)       :: array(:,:)
    character(len=*), intent(in) :: array_name
    sll_int32, intent(out)       :: error
    sll_int32                    :: file_id
    sll_int32                    :: npoints_x1
    sll_int32                    :: npoints_x2
    sll_int32, intent(in), optional :: xmffile_id
    character(len=4), optional      :: center
    
    npoints_x1 = size(array,1)
    npoints_x2 = size(array,2)
#ifdef NOHDF5
    call sll_binary_file_create(trim(prefix)//"-"//trim(array_name)//".bin",file_id,error)
    call sll_binary_write_array(file_id,array,error)
    call sll_binary_file_close(file_id,error)
#else
    call sll_hdf5_file_create(trim(prefix)//"-"//trim(array_name)//".h5",file_id,error)
    call sll_hdf5_write_array(file_id,array,"/"//trim(array_name),error)
    call sll_hdf5_file_close(file_id, error)
#endif

    if ( present(xmffile_id) .and. present(center)) then
#ifdef NOHDF5
       call sll_xml_field(xmffile_id,trim(array_name),trim(prefix)//"-"//trim(array_name)//".bin", &
            npoints_x1,npoints_x2,'Binary',center)
#else
       call sll_xml_field( &
            xmffile_id, &
            trim(array_name), &
            trim(prefix)//"-"//trim(array_name)//".h5:/"//trim(array_name), &
            npoints_x1,&
            npoints_x2,&
            'HDF',&
            center)
#endif
  end if
  end subroutine sll_xdmf_array_2d

  subroutine sll_xdmf_open_3d( &
    prefix, &
    file_id, &
    nnodes_x1, &
    nnodes_x2, &
    nnodes_x3, &
    error)
    
    character(len=*), intent(in) :: prefix
    sll_int32, intent(out)       :: file_id
    sll_int32, intent(out)       :: error
    sll_int32                    :: nnodes_x1
    sll_int32                    :: nnodes_x2
    sll_int32                    :: nnodes_x3
    call sll_xml_file_create(trim(prefix)//".xmf",file_id,error)
    call sll_xml_grid_geometry(file_id, trim(prefix), nnodes_x1, nnodes_x2, nnodes_x3)
  end subroutine sll_xdmf_open_3d

  subroutine sll_xdmf_array_3d(prefix,array,array_name,error,xmffile_id,center)
    character(len=*), intent(in) :: prefix
    sll_real64, intent(in)       :: array(:,:,:)
    character(len=*), intent(in) :: array_name
    sll_int32, intent(out)       :: error
    sll_int32                    :: file_id
    sll_int32, intent(in), optional :: xmffile_id
    character(len=4), optional      :: center
    sll_int32                       :: npoints_x1
    sll_int32                       :: npoints_x2
    sll_int32                       :: npoints_x3
    
    npoints_x1 = size(array,1)
    npoints_x2 = size(array,2)
    npoints_x3 = size(array,3)

#ifdef NOHDF5
    call sll_binary_file_create(trim(prefix)//"-"//trim(array_name)//".bin",file_id,error)
    call sll_binary_write_array(file_id,array,error)
    call sll_binary_file_close(file_id,error)
#else
    call sll_hdf5_file_create(trim(prefix)//"-"//trim(array_name)//".h5",file_id,error)
    call sll_hdf5_write_array(file_id,array,"/"//trim(array_name),error)
    call sll_hdf5_file_close(file_id, error)
#endif


    if ( present(xmffile_id) .and. present(center)) then

#ifdef NOHDF5
    call sll_xml_field(xmffile_id,trim(array_name),trim(prefix)//"-"//trim(array_name)//".bin", &
         npoints_x1,npoints_x2,npoints_x3,'Binary',center)
#else
    call sll_xml_field( &
         xmffile_id, &
         trim(array_name), &
         trim(prefix)//"-"//trim(array_name)//".h5:/"//trim(array_name), &
         npoints_x1,&
         npoints_x2,&
         npoints_x3,&
         'HDF', &
         center)
#endif
   end if

  end subroutine sll_xdmf_array_3d

  !>Writes data for mesh plotting:
  !> \param[in] mesh_prefix is the filename prefix for the mesh
  !> \param[in] x1 is array with nodes coordinates  along direction 1
  !> \param[in] x2 is array with nodes coordinates  along direction 2
  !> \param[in] nnodes_x1 is the number of nodes number along direction 1
  !> \param[in] nnodes_x2 is the number of nodes number along direction 2
  !> \param[out] error is a parameter that should be 0
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Write XDMF file for 2D mesh plot
  !!
  !! Light data are in the file:
  !!  - <c>mesh_prefix.xmf</c>.
  !!
  !! Heavy data (mesh nodes coordinates) are in :
  !!  - <c>mesh_prefix.h5</c>
  !!  - <c>mesh_prefix-x1.bin</c> and <c>mesh_prefix-x2.bin</c>
  !!
  !! \param[in] x1 array with first coordinates
  !! \param[in] x2 array with second coordinates
  !! \param[in] nnodes_x1 nodes number along direction 1
  !! \param[in] nnodes_x2 nodes number along direction 2
  !! \param[in] mesh_prefix filename prefix for the mesh
  !!
  !! Example of a cartesian grid (xmin,xmax) x (vmin,vmax):
  !!\verbatim
  !!program grid
  !!#include "sll_memory.h"
  !!use sll_low_level_file_io
  !!sll_int64, dimension(:,:), allocatable :: x, v
  !!sll_int32 :: nx, nv
  !!
  !!nx = 32; xmin = 0.; xmax = 1.
  !!nv = 64; vmin = 0.; vmax = 1.
  !!SLL_ALLOCATE(x(nx,nv))
  !!SLL_ALLOCATE(v(nx,nv)) 
  !!do j = 1, nv                                             
  !!   do i = 1, nx                                         
  !!      x(i,j) = float(i-1) * (xmax-xmin) / (nx-1)       
  !!      v(i,j) = float(j-1) * (vmax-vmin) / (nv-1)      
  !!   end do                                            
  !!end do                                              
  !!call write_mesh(x,v,nx,nv,"mesh")
  !!end program
  !!\endverbatim
  !! A file <i>mesh.xmf</i> is created, readable by VisIt or Paraview. For non
  !!moving mesh, call it only one time, heavy data produced will be reused.
  !!
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Write XDMF file for vector 1d on 2d mesh with data in binary or HDF5 files
  !!
  !! Light data are in :
  !!  - <c>field_prefix.xmf</c>.
  !! Heavy data :
  !!  - mesh nodes coordinates are in :
  !!     - <c>mesh_prefix.h5</c>
  !!     - <c>mesh_prefix-x1.bin</c> and <c>mesh_prefix-x2.bin</c>
  !!  - scalar values are in :
  !!     - <c>field_prefix.h5</c>
  !!     - <c>field_prefix.bin</c>
  !!
  !! \param[in] vec_values_x1 Array with vector values
  !! \param[in] nnodes_x1 nodes number along direction 1
  !! \param[in] nnodes_x2 nodes number along direction 2
  !! \param[in] field_prefix filename prefix for the vector
  !! \param[in] mesh_prefix filename prefix for the mesh
  !! \param[in] icenter parameter to distinguish cells (0) or nodes (1) values
  !!
  !! Example of f(x,v) = xv on a cartesian grid (xmin,xmax) x (vmin,vmax):
  !!\verbatim
  !!program grid
  !!#include "sll_memory.h"
  !!#include "sll_working_precision.h"
  !!use sll_low_level_file_io
  !!sll_int64, dimension(:,:), allocatable :: x, v
  !!sll_int32 :: nx, nv, error
  !!
  !!nx = 32; xmin = 0.; xmax = 1.
  !!nv = 64; vmin = 0.; vmax = 1.
  !!SLL_ALLOCATE(x(nx,nv),error)
  !!SLL_ALLOCATE(v(nx,nv),error) 
  !!SLL_ALLOCATE(f(nx,nv),error) 
  !!do j = 1, nv                                             
  !!   do i = 1, nx                                         
  !!      x(i,j) = float(i-1) * (xmax-xmin) / (nx-1)       
  !!      v(i,j) = float(j-1) * (vmax-vmin) / (nv-1)      
  !!      f(i,j) = x(i,j) * v(i,j)
  !!   end do                                            
  !!end do                                              
  !!call write_mesh(x,v,nx,nv,"mesh")
  !!call write_vec1d(f,nx,nv,"f","mesh",0) !write cell centered values
  !!end program
  !!\endverbatim
  !! A file <i>f.xmf</i> is created, readable by VisIt or Paraview. 
  !!
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Write XDMF file for vector 2d on 2d mesh with data in binary or hdf5 files
  !!
  !! Light data are in :
  !!  - <c>field_prefix.xmf</c>.
  !! Heavy data :
  !!  - mesh nodes coordinates created by write_mesh call are in :
  !!     - <c>mesh_prefix.h5</c>
  !!     - <c>mesh_prefix-x1.bin</c> and <c>mesh_prefix-x2.bin</c>
  !!  - scalar values are in :
  !!     - <c>field_prefix.h5</c>
  !!     - <c>field_prefix-x1.bin</c> and <c>field_prefix-x2.bin</c>
  !!
  !! \param[in] vec_values_x1 Array with vector values
  !! \param[in] vec_values_x2 Array with vector values
  !! \param[in] nnodes_x1 nodes number along direction 1
  !! \param[in] nnodes_x2 nodes number along direction 2
  !! \param[in] field_prefix filename prefix for the vector
  !! \param[in] mesh_prefix filename prefix for the mesh
  !! \param[in] icenter parameter to distinguish cells values or nodes values
  !!
  !! Example of f1(x,v) = xv , 
  !! f2(x,v) = x+v on a cartesian grid (xmin,xmax) x (vmin,vmax):
  !!\verbatim
  !!program grid
  !!#include "sll_memory.h"
  !!#include "sll_working_precision.h"
  !!use sll_low_level_file_io
  !!sll_int64, dimension(:,:), allocatable :: x, v
  !!sll_int32 :: nx, nv, error
  !!
  !!nx = 32; xmin = 0.; xmax = 1.
  !!nv = 64; vmin = 0.; vmax = 1.
  !!SLL_ALLOCATE(x(nx,nv),error)
  !!SLL_ALLOCATE(v(nx,nv),error) 
  !!SLL_ALLOCATE(f1(nx,nv),error) 
  !!SLL_ALLOCATE(f2(nx,nv),error) 
  !!do j = 1, nv                                             
  !!   do i = 1, nx                                         
  !!      x(i,j) = float(i-1) * (xmax-xmin) / (nx-1)       
  !!      v(i,j) = float(j-1) * (vmax-vmin) / (nv-1)      
  !!      f1(i,j) = x(i,j) * v(i,j)
  !!      f2(i,j) = x(i,j) + v(i,j)
  !!   end do                                            
  !!end do                                              
  !!call write_mesh(x,v,nx,nv,"mesh")
  !!call write_vec2d(f1,f2,nx,nv,"f","mesh",0) !write cell centered values
  !!end program
  !!\endverbatim
  !! A file <i>f.xmf</i> is created, readable by VisIt or Paraview. 
  !!
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Write XDMF file for vector 3d on 3d mesh with data in binary or hdf5 files
  !!
  !! Light data are in :
  !!  - <c>field_prefix.xmf</c>.
  !! Heavy data :
  !!  - mesh nodes coordinates created by write_mesh call are in :
  !!     - <c>mesh_prefix.h5</c>
  !!     - <c>mesh_prefix-x123.bin</c>
  !!  - scalar values are in :
  !!     - <c>field_prefix.h5</c>
  !!     - <c>field_prefix-x123.bin</c>
  !!
  !! \param[in] vec_values_x1 Array with vector values
  !! \param[in] vec_values_x2 Array with vector values
  !! \param[in] vec_values_x3 Array with vector values
  !! \param[in] nnodes_x1 nodes number along direction 1
  !! \param[in] nnodes_x2 nodes number along direction 2
  !! \param[in] nnodes_x3 nodes number along direction 2
  !! \param[in] field_prefix filename prefix for the vector
  !! \param[in] mesh_prefix filename prefix for the mesh
  !! \param[in] icenter parameter to distinguish cells values or nodes values
  !!
  
  
  !----------------------------------------------------------------------------
  !> Outputs an error message:
  !>   - PRTFIL : unit number for print-out
  !>   - SEVRTY : 'W' - Warning 'F' - Fatal
  !>   - WHERE  : in which program or subroutine
  !>   - ErrMsg : error message
  subroutine errout( prtfil, sevrty, lwhere, ErrMsg )
    sll_int32, intent(in) ::  prtfil
    character(len=1),intent(in) :: sevrty 
    character(len=*),intent(in) :: lwhere , ErrMsg
    
    write( prtfil, * )
    select case ( sevrty )  !     *** Severity ***
    case ( 'W' )
       write(prtfil,"(/10x,a)") '*** WARNING ***'
    case ( 'F' )
       write(prtfil,"(/10x,a)") '*** FATAL ERROR ***'
    case default
       write(prtfil,"(/10x,a)") '*** FATAL ERROR ***'
       write(prtfil,"(/10x,a)") &
            'Error handler (ERROUT) called with unknown severity level: ', &
            SEVRTY
    end select
    write( prtfil,"(/10x,a)") &
         'Generated by program or subroutine: ', trim(lwhere)
    write( prtfil,"(/10x,a)") trim(ErrMsg)
    write( prtfil,"(/10x,a)")
    
    ! return or stop depending on severity
    if ( sevrty == 'W' ) then
       return
    else
       stop 'Fatal Error: See print file for details'
    end if
    
  end subroutine errout
  
end module sll_xdmf
