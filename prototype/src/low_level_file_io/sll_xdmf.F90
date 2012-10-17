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
  
#ifdef NOHDF5
  use sll_binary_io
#else
  use sll_hdf5_io
#endif

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
  
  !>Open a XDMF format file for a 2d plot
  subroutine sll_xdmf_open_2d(file_name,mesh_name,nnodes_x1,nnodes_x2,file_id,error)

    character(len=*), intent(in) :: file_name
    character(len=*), intent(in) :: mesh_name
    sll_int32, intent(out)       :: file_id
    sll_int32, intent(out)       :: error
    sll_int32                    :: nnodes_x1
    sll_int32                    :: nnodes_x2
    
    call sll_xml_file_create(trim(file_name),file_id,error)
    call sll_xml_grid_geometry(file_id, trim(mesh_name), nnodes_x1, nnodes_x2)

  end subroutine sll_xdmf_open_2d

  !>Open a XDMF format file for a 3d plot
  subroutine sll_xdmf_open_3d( &
    file_name, &
    mesh_name, &
    nnodes_x1, &
    nnodes_x2, &
    nnodes_x3, &
    file_id,   &
    error)
    
    character(len=*), intent(in) :: file_name
    character(len=*), intent(in) :: mesh_name
    sll_int32, intent(out)       :: file_id
    sll_int32, intent(out)       :: error
    sll_int32                    :: nnodes_x1
    sll_int32                    :: nnodes_x2
    sll_int32                    :: nnodes_x3

    call sll_xml_file_create(trim(file_name),file_id,error)
    call sll_xml_grid_geometry(file_id, trim(mesh_name), nnodes_x1, nnodes_x2, nnodes_x3)

  end subroutine sll_xdmf_open_3d
  
  !>Write 2d array in binary or hdf5 file and the matching line in XDMF file
  subroutine sll_xdmf_array_2d(mesh_name,array,array_name,error,xmffile_id,center)

    character(len=*), intent(in)    :: mesh_name
    sll_real64, intent(in)          :: array(:,:)
    character(len=*), intent(in)    :: array_name
    sll_int32, intent(out)          :: error
    sll_int32                       :: file_id
    sll_int32                       :: npoints_x1
    sll_int32                       :: npoints_x2
    sll_int32, intent(in), optional :: xmffile_id
    character(len=4), optional      :: center
    
    npoints_x1 = size(array,1)
    npoints_x2 = size(array,2)

#ifdef NOHDF5
    call sll_binary_file_create(trim(mesh_name)//"-"//trim(array_name)//".bin", &
                                file_id,error)
    call sll_binary_write_array(file_id,array,error)
    call sll_binary_file_close(file_id,error)
#else
    call sll_hdf5_file_create(trim(mesh_name)//"-"//trim(array_name)//".h5", &
                              file_id,error)
    call sll_hdf5_write_array(file_id,array,"/"//trim(array_name),error)
    call sll_hdf5_file_close(file_id, error)
#endif

    if ( present(xmffile_id) .and. present(center)) then
#ifdef NOHDF5
       call sll_xml_field(xmffile_id,trim(array_name), &
                          trim(mesh_name)//"-"//trim(array_name)//".bin", &
                          npoints_x1,npoints_x2,'Binary',center)
#else
       call sll_xml_field( &
            xmffile_id, &
            trim(array_name), &
            trim(mesh_name)//"-"//trim(array_name)//".h5:/"//trim(array_name), &
            npoints_x1,&
            npoints_x2,&
            'HDF',&
            center)
#endif
     end if
  end subroutine sll_xdmf_array_2d

  !>Write 3d array in binary or hdf5 file and the matching line in XDMF file
  subroutine sll_xdmf_array_3d(mesh_name,array,array_name,error,xmffile_id,center)

    character(len=*), intent(in)    :: mesh_name
    sll_real64, intent(in)          :: array(:,:,:)
    character(len=*), intent(in)    :: array_name
    sll_int32, intent(out)          :: error
    sll_int32                       :: file_id
    sll_int32, intent(in), optional :: xmffile_id
    character(len=4), optional      :: center
    sll_int32                       :: npoints_x1
    sll_int32                       :: npoints_x2
    sll_int32                       :: npoints_x3
    
    npoints_x1 = size(array,1)
    npoints_x2 = size(array,2)
    npoints_x3 = size(array,3)

#ifdef NOHDF5
    call sll_binary_file_create(trim(mesh_name)//"-"//trim(array_name)//".bin", &
                                file_id,error)
    call sll_binary_write_array(file_id,array,error)
    call sll_binary_file_close(file_id,error)
#else
    call sll_hdf5_file_create(trim(mesh_name)//"-"//trim(array_name)//".h5", &
                              file_id,error)
    call sll_hdf5_write_array(file_id,array,"/"//trim(array_name),error)
    call sll_hdf5_file_close(file_id, error)
#endif


    if ( present(xmffile_id) .and. present(center)) then

#ifdef NOHDF5
    call sll_xml_field(xmffile_id,trim(array_name), &
                       trim(mesh_name)//"-"//trim(array_name)//".bin", &
                       npoints_x1,npoints_x2,npoints_x3,'Binary',center)
#else
    call sll_xml_field( &
         xmffile_id, &
         trim(array_name), &
         trim(mesh_name)//"-"//trim(array_name)//".h5:/"//trim(array_name), &
         npoints_x1,&
         npoints_x2,&
         npoints_x3,&
         'HDF', &
         center)
#endif

     end if

  end subroutine sll_xdmf_array_3d

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
