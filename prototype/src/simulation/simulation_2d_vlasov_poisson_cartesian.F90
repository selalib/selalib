module sll_simulation_2d_vlasov_poisson_cartesian

! mimic of VP1D_deltaf_BSL in the simulation class
! intend to use MPI instead of openmp
! main application is for KEEN waves


#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
  use sll_collective
  use sll_remapper
  use sll_constants
  use sll_cubic_spline_interpolator_1d
  use sll_test_4d_initializer
  use sll_poisson_2d_periodic_cartesian_par
  use sll_cubic_spline_interpolator_1d
  use sll_simulation_base
  implicit none

  type, extends(sll_simulation_base_class) :: &
       sll_simulation_2d_vlasov_poisson_cart
   contains
     procedure, pass(sim) :: run => run_vp2d_cartesian
     procedure, pass(sim) :: init_from_file => init_vp2d_par_cart
  end type sll_simulation_2d_vlasov_poisson_cart

  interface delete
     module procedure delete_vp2d_par_cart
  end interface delete

contains

  subroutine init_vp2d_par_cart( sim, filename )
    intrinsic :: trim
    class(sll_simulation_2d_vlasov_poisson_cart), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
    sll_int32             :: IO_stat
    sll_real64            :: dt
    sll_int32             :: number_iterations
    sll_int32             :: num_cells_x1
    sll_int32             :: num_cells_x2
    sll_int32, parameter  :: input_file = 99

    namelist /sim_params/ dt, number_iterations
    namelist /grid_dims/ num_cells_x1, num_cells_x2
    ! Try to add here other parameters to initialize the mesh values like
    ! xmin, xmax and also for the distribution function initializer.
    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, '#init_vp2d_par_cart() failed to open file ', filename
       STOP
    end if
    read(input_file, sim_params)
    read(input_file,grid_dims)
    close(input_file)

    !sim%dt = dt
    !sim%num_iterations = number_iterations
    ! In this particular simulation, since the system is periodic, the number
    ! of points is the same as the number of cells in all directions.
    !sim%nc_x1 = num_cells_x1
    !sim%nc_x2 = num_cells_x2
  end subroutine init_vp2d_par_cart


  subroutine run_vp2d_cartesian(sim)
    class(sll_simulation_2d_vlasov_poisson_cart), intent(inout) :: sim
  end subroutine run_vp2d_cartesian

  subroutine delete_vp2d_par_cart( sim )
    class(sll_simulation_2d_vlasov_poisson_cart) :: sim
    sll_int32 :: ierr
  end subroutine delete_vp2d_par_cart


end module sll_simulation_2d_vlasov_poisson_cartesian
