program test_xdmf_parallel

#include "sll_working_precision.h"

  use sll_collective, only: &
    sll_world_collective,   &
    sll_boot_collective,    &
    sll_halt_collective,    &
    sll_get_collective_rank

  use sll_m_xdmf_parallel, only: &
    sll_t_xdmf_parallel_file

  implicit none

  !----------------------------------------------------------------------------
  ! VARIABLES DECLARATION
  !----------------------------------------------------------------------------

  type(sll_t_xdmf_parallel_file) :: xdmf_file
  sll_int32                      :: gid_cart, gid_polar
  character(len=256)             :: reference_filename
  logical                        :: file_exists

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
  ! START PARALLEL ENVIRONMENT
  !----------------------------------------------------------------------------

  call sll_boot_collective()

  !----------------------------------------------------------------------------
  ! XDMF FILE CREATION
  !----------------------------------------------------------------------------

  ! Initialize with time and MPI communicator
  !------------------------------------------
  call xdmf_file%init( time=8.0_f64, comm=sll_world_collective )

  ! Add Grid 1 to Domain
  !---------------------
  call xdmf_file%add_grid( &
    grid_name = 'mesh_x2x3_cart', &
    x1_path   = 'mesh_x2x3_cart.h5:/x2', &
    x2_path   = 'mesh_x2x3_cart.h5:/x3', &
    dims      = [33,33], &
    gid       = gid_cart )

  ! Add Grid 2 to Domain
  !---------------------
  call xdmf_file%add_grid( &
    grid_name = 'mesh_x1x2_polar', &
    x1_path   = 'mesh_x1x2_polar.h5:/x1', &
    x2_path   = 'mesh_x1x2_polar.h5:/x2', &
    dims      = [33,33], &
    gid       = gid_polar )

  ! Add 2D dataset to Grid 1
  !-------------------------
  call xdmf_file%add_field( &
    grid_id    = gid_cart, &
    field_name = 'f_x2x3', &
    field_path = 'diag2d_0001.h5:/f_x2x3' )

  ! Add 2D dataset to Grid 2
  !-------------------------
  call xdmf_file%add_field( &
    grid_id    = gid_polar, &
    field_name = 'f_x1x2', &
    field_path = 'diag2d_0001.h5:/f_x1x2' )

  ! Write XDMF file, then delete object
  !------------------------------------
  call xdmf_file%write( 'out.xdmf' )
  call xdmf_file%delete()

  !----------------------------------------------------------------------------
  ! UNIT TESTING
  !----------------------------------------------------------------------------

  if (sll_get_collective_rank( sll_world_collective ) == 0) then

    ! Compare to reference file
    !--------------------------
    if (equal_files( 'out.xdmf', trim( reference_filename ) )) then
      write(*,*) "PASSED"
    else
      write(*,*) "ERROR: output file does not match reference"
    end if

    ! Remove temporary file
    !----------------------
!    call remove_file( 'out.xdmf' )

  end if

  !----------------------------------------------------------------------------
  ! END PARALLEL ENVIRONMENT
  !----------------------------------------------------------------------------

  call sll_halt_collective()

!==============================================================================
contains
!==============================================================================

  subroutine remove_file( filename )
    character(len=*), intent(in) :: filename

    integer :: iunit

    open( newunit=iunit, file=filename, status='OLD' )
    close( iunit, status='DELETE' )

  end subroutine

  !----------------------------------------------------------------------------
  subroutine read_file( filename, str )
    character(len=*),              intent(in   ) :: filename
    character(len=:), allocatable, intent(  out) :: str

    integer :: iunit, istat, filesize
    character(len=1) :: c

    ! Open required file
    open( newunit=iunit, file=filename, status='OLD', &
      form='UNFORMATTED', access='STREAM' )

    ! How many characters in file
    inquire( file=filename, size=filesize )
    if (filesize < 0) then
      write(*,*) 'ERROR: negative file size'
      stop
    end if

    ! Read whole file into one string
    allocate( character(len=filesize) :: str )
    read( iunit, pos=1 ) str

    ! Make sure it was all read by trying to read more
    read( iunit, pos=filesize+1, iostat=istat ) c
    if (.not. IS_IOSTAT_END(istat)) then
      write(*,*) 'Error: file was not completely read'
      stop
    end if

    ! Close file
    close( iunit, iostat=istat )

  end subroutine read_file

  !----------------------------------------------------------------------------
  function equal_files( filename1, filename2 ) result( equal )
    character(len=*), intent(in   ) :: filename1
    character(len=*), intent(in   ) :: filename2
    logical :: equal

    character(len=:), allocatable :: str1
    character(len=:), allocatable :: str2

    call read_file( filename1, str1 )
    call read_file( filename2, str2 )

    equal = (str1 == str2)
    
  end function equal_files

  !----------------------------------------------------------------------------
  function empty_file( filename ) result( empty )
    character(len=*), intent(in) :: filename
    logical :: empty

    integer :: filesize

    inquire( file=filename, size=filesize )
    empty = (filesize == 0)

  end function empty_file

end program test_xdmf_parallel
