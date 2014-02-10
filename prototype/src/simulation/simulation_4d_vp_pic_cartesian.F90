module sll_pic_simulation_4d_cartesian_module

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"

  use sll_constants
  use sll_simulation_base
  use sll_logical_meshes
  use sll_timer
  implicit none
  
  type, extends(sll_simulation_base_class) :: sll_pic_simulation_4d_cartesian
     ! Physics/numerical parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
   contains
     procedure, pass(sim) :: run => run_4d_pic_cart
     procedure, pass(sim) :: init_from_file => init_4d_pic_cartesian
  end type sll_pic_simulation_4d_cartesian

  interface sll_delete
     module procedure delete_4d_pic_cart
  end interface sll_delete

!!$  interface initialize
!!$     module procedure initialize_4d_qns_general
!!$  end interface initialize

contains
  
  ! Tentative function to initialize the simulation object 'manually'.
  subroutine initialize_4d_pic_cartesian( &
       sim )
    
    type(sll_pic_simulation_4d_cartesian), intent(inout)     :: sim
  end subroutine initialize_4d_pic_cartesian


  subroutine init_4d_pic_cartesian( sim, filename )
    intrinsic :: trim
    class(sll_pic_simulation_4d_cartesian), intent(inout) :: sim
    character(len=*), intent(in)                          :: filename
    sll_int32             :: IO_stat
    sll_real64            :: dt
    sll_int32             :: number_iterations
    sll_int32             :: num_cells_x1
    sll_int32             :: num_cells_x2
    sll_int32             :: num_cells_x3
    sll_int32             :: num_cells_x4
    sll_int32, parameter  :: input_file = 99

    namelist /sim_params/ dt, number_iterations
    namelist /grid_dims/ num_cells_x1, num_cells_x2, num_cells_x3, num_cells_x4
    ! Try to add here other parameters to initialize the mesh values like
    ! xmin, xmax and also for the distribution function initializer.
    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_vp4d_par_cart() failed to open file ', filename
       STOP
    end if
    read(input_file, sim_params)
    read(input_file,grid_dims)
    close(input_file)

    sim%dt = dt
    sim%num_iterations = number_iterations
  end subroutine init_4d_pic_cartesian

  ! Note that the following function has no local variables, which is silly...
  ! This just happened since the guts of the unit test were transplanted here
  ! directly, but this should be cleaned up.
  subroutine run_4d_pic_cart(sim)
    class(sll_pic_simulation_4d_cartesian), intent(inout) :: sim
  end subroutine run_4d_pic_cart


  subroutine delete_4d_pic_cart( sim )
    type(sll_pic_simulation_4d_cartesian) :: sim
    sll_int32 :: ierr
  end subroutine delete_4d_pic_cart


end module sll_pic_simulation_4d_cartesian_module
