program test_xdmf_parallel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_io_utilities, only: &
    sll_f_check_equal_files

  use sll_m_xdmf_light_parallel, only: &
    sll_t_xdmf_parallel_file

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !----------------------------------------------------------------------------
  ! VARIABLES DECLARATION
  !----------------------------------------------------------------------------

  type(sll_t_xdmf_parallel_file) :: xdmf_file
  sll_int32                      :: send_rank, gid_cart, gid_polar
  character(len=256)             :: reference_filename
  logical                        :: file_exists, to_file

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
  inquire( file=reference_filename, exist=file_exists )
  if (.not. file_exists) then
    write(*,*) &
      "ERROR: reference file '"//trim( reference_filename )//"' does not exist"
    stop
  end if

  !----------------------------------------------------------------------------
  ! START PARALLEL ENVIRONMENT
  !----------------------------------------------------------------------------

  call sll_s_boot_collective()

  ! Choose which processor will send data to the XDMF file
  send_rank =  sll_f_get_collective_size( sll_v_world_collective ) -1
  to_file   = (sll_f_get_collective_rank( sll_v_world_collective ) == send_rank)

  !----------------------------------------------------------------------------
  ! XDMF FILE CREATION
  !----------------------------------------------------------------------------

  ! Initialize with time and MPI communicator
  !------------------------------------------
  call xdmf_file%init( time=8.0_f64, comm=sll_v_world_collective )

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
    field_path = 'diag2d_0001.h5:/f_x2x3', &
    to_file    = to_file )

  ! Add 2D dataset to Grid 2
  !-------------------------
  call xdmf_file%add_field( &
    grid_id    = gid_polar, &
    field_name = 'f_x1x2', &
    field_path = 'diag2d_0001.h5:/f_x1x2', &
    to_file    = to_file )

  ! Write XDMF file, then delete object
  !------------------------------------
  call xdmf_file%write( 'out.xdmf' )
  call xdmf_file%delete()

  !----------------------------------------------------------------------------
  ! UNIT TESTING
  !----------------------------------------------------------------------------

  if (sll_f_get_collective_rank( sll_v_world_collective ) == 0) then

    ! Compare to reference file
    !--------------------------
    if (sll_f_check_equal_files( 'out.xdmf', reference_filename )) then
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

  call sll_s_halt_collective()

end program test_xdmf_parallel
