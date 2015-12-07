program test_xml

  use sll_m_xml, only:  &
    sll_t_xml_document, &
    sll_t_xml_element

  use sll_m_io_utilities, only: &
    sll_f_check_equal_files,    &
    sll_f_check_empty_file,     &
    sll_s_remove_file

  implicit none

  !----------------------------------------------------------------------------
  ! VARIABLES DECLARATION
  !----------------------------------------------------------------------------

  type(sll_t_xml_document)         :: xml_doc
  type(sll_t_xml_element), pointer :: root, domain, grid1, grid2
  type(sll_t_xml_element), pointer :: time, topology, geometry, dataitem, field
  !character(len=256)               :: header

  character(len=256)               :: reference_filename
  logical                          :: file_exists, equal, empty

  !----------------------------------------------------------------------------
  ! PARSE INPUT
  !----------------------------------------------------------------------------

  ! Check that input argument was given
  !------------------------------------
  if (command_argument_count() /= 1 ) then
    write(*,*) "ERROR: exactly 1 input argument is required"
    stop
  end if

  ! Read name of reference file from input argument
  !------------------------------------------------
  call get_command_argument( 1, reference_filename )

  ! Check that file exists    
  !-----------------------
  inquire( file=trim( reference_filename ), exist=file_exists )
  if (.not. file_exists) then
    write(*,*) &
      "ERROR: reference file '"//trim( reference_filename )//"' does not exist"
    stop
  end if

  !----------------------------------------------------------------------------
  ! XML FILE CREATION
  !----------------------------------------------------------------------------

  ! Add header
  !-----------
  call xml_doc%add_header_line( "<?xml version='1.0' ?>" )
  call xml_doc%add_header_line( "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>" )

  ! Create root element (only one may exist!)
  !------------------------------------------
  root => xml_doc%new_element( 'Xdmf' )
  call root%add_attribute( 'Version', '2.0' )

  ! Add Domain element to root
  !---------------------------
  domain => root%new_element( 'Domain' )

  ! Add Grid 1 to Domain
  !---------------------
  grid1 => domain%new_element( 'Grid' )
  call grid1%add_attribute( 'Name', 'mesh_x2x3_cart' )
  call grid1%add_attribute( 'GridType', 'Uniform' )

  ! Add Grid 2 to Domain
  !---------------------
  grid2 => domain%new_element( 'Grid' )
  call grid2%add_attribute( 'Name', 'mesh_x1x2_polar' )
  call grid2%add_attribute( 'GridType', 'Uniform' )

  ! Add time to both grids
  !-----------------------
  time => grid1%new_element( 'Time' ); call time%add_attribute( 'Value','8.0' )
  time => grid2%new_element( 'Time' ); call time%add_attribute( 'Value','8.0' )

  ! Add topology to both grids
  !---------------------------
  topology => grid1%new_element( 'Topology' )
  call topology%add_attribute( 'TopologyType', '2DSMesh' )
  call topology%add_attribute( 'NumberOfElements', '33 33' )

  topology => grid2%new_element( 'Topology' )
  call topology%add_attribute( 'TopologyType', '2DSMesh' )
  call topology%add_attribute( 'NumberOfElements', '33 33' )

  ! Add Geometry to Grid 1
  !-----------------------
  geometry => grid1%new_element( 'Geometry' )
  call geometry%add_attribute( 'GeometryType', 'X_Y' )

  dataitem => geometry%new_element( 'DataItem' )
  call dataitem%add_attribute( 'Dimensions', '33 33' )
  call dataitem%add_attribute( 'NumberType', 'Float' )
  call dataitem%add_attribute( 'Precision' , '8'     )
  call dataitem%add_attribute( 'Format'    , 'HDF'   )
  call dataitem%add_chardata( 'mesh_x2x3_cart.h5:/x2' ) ! dataset x2-cart

  dataitem => geometry%new_element( 'DataItem' )
  call dataitem%add_attribute( 'Dimensions', '33 33' )
  call dataitem%add_attribute( 'NumberType', 'Float' )
  call dataitem%add_attribute( 'Precision' , '8'     )
  call dataitem%add_attribute( 'Format'    , 'HDF'   )
  call dataitem%add_chardata( 'mesh_x2x3_cart.h5:/x3' ) ! dataset x3-cart

  ! Add Geometry to Grid 2
  !-----------------------
  geometry => grid2%new_element( 'Geometry' )
  call geometry%add_attribute( 'GeometryType', 'X_Y' )

  dataitem => geometry%new_element( 'DataItem' )
  call dataitem%add_attribute( 'Dimensions', '33 33' )
  call dataitem%add_attribute( 'NumberType', 'Float' )
  call dataitem%add_attribute( 'Precision' , '8'     )
  call dataitem%add_attribute( 'Format'    , 'HDF'   )
  call dataitem%add_chardata( 'mesh_x1x2_polar.h5:/x1' ) ! dataset x1-polar

  dataitem => geometry%new_element( 'DataItem' )
  call dataitem%add_attribute( 'Dimensions', '33 33' )
  call dataitem%add_attribute( 'NumberType', 'Float' )
  call dataitem%add_attribute( 'Precision' , '8'     )
  call dataitem%add_attribute( 'Format'    , 'HDF'   )
  call dataitem%add_chardata( 'mesh_x1x2_polar.h5:/x2' ) ! dataset x2-polar

  ! Add 2D dataset to Grid 1
  !-------------------------
  field => grid1%new_element( 'Attribute' )
  call field%add_attribute( 'Name'         , 'f_x2x3' )
  call field%add_attribute( 'AttributeType', 'Scalar' )
  call field%add_attribute( 'Center'       , 'Node'   )

  dataitem => field%new_element( 'DataItem' )
  call dataitem%add_attribute( 'Dimensions', '33 33' )
  call dataitem%add_attribute( 'NumberType', 'Float' )
  call dataitem%add_attribute( 'Precision' , '8'     )
  call dataitem%add_attribute( 'Format'    , 'HDF'   )
  call dataitem%add_chardata( 'diag2d_0001.h5:/f_x2x3' ) ! dataset f_x2x3-cart

  ! Add 2D dataset to Grid 2
  !-------------------------
  field => grid2%new_element( 'Attribute' )
  call field%add_attribute( 'Name'         , 'f_x1x2' )
  call field%add_attribute( 'AttributeType', 'Scalar' )
  call field%add_attribute( 'Center'       , 'Node'   )

  dataitem => field%new_element( 'DataItem' )
  call dataitem%add_attribute( 'Dimensions', '33 33' )
  call dataitem%add_attribute( 'NumberType', 'Float' )
  call dataitem%add_attribute( 'Precision' , '8'     )
  call dataitem%add_attribute( 'Format'    , 'HDF'   )
  call dataitem%add_chardata( 'diag2d_0001.h5:/f_x1x2' ) ! dataset f_x1x2-cart

  ! Write XML file
  !---------------
  call xml_doc%write( 'out1.xml' )

  ! Deallocate/delete
  !------------------
  call xml_doc%delete()
  call xml_doc%write( 'out2.xml' )

  !----------------------------------------------------------------------------
  ! UNIT TESTING
  !----------------------------------------------------------------------------

  ! Compare to reference files
  !---------------------------
  equal = sll_f_check_equal_files( 'out1.xml', reference_filename )
  empty = sll_f_check_empty_file ( 'out2.xml' )

  ! Remove temporary files
  !-----------------------
  call sll_s_remove_file( 'out1.xml' )
  call sll_s_remove_file( 'out2.xml' )

  ! Check test
  !-----------
  if (equal .and. empty) then
    write(*,*) "PASSED"
  else
    if (.not. equal) write(*,*) "ERROR: output file 1 does not match reference"
    if (.not. empty) write(*,*) "ERROR: output file 2 is not empty"
  end if

end program test_xml
