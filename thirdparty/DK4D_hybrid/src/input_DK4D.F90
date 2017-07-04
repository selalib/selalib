!===========================================================================
!> Input data for 4D drift-kinetic hybrid simulation
!>
!> \date 2014-08-19
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module input_DK4D_module
#include "sll_working_precision.h"

  use utils_DK4D_module

  implicit none
  
  !> For BSL scheme type choice
  sll_int32, parameter, public :: BSL_POLAR  = 1
  sll_int32, parameter, public :: BSL_HYBRID = 2

  !> For time loop choice
  sll_int32, parameter, public :: TIME_LOOP_EULER  = 1
  sll_int32, parameter, public :: TIME_LOOP_PREDICTOR_CORRECTOR = 2


  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: input_DK4D_t

    !> Mesh parameters
    sll_int32  :: nc_x1
    sll_int32  :: nc_x2
    sll_int32  :: nc_x3
    sll_int32  :: nc_x4
    sll_real64 :: r_min
    sll_real64 :: r_max
    sll_real64 :: eta3_min
    sll_real64 :: eta3_max
    sll_real64 :: vpar_min
    sll_real64 :: vpar_max

    !> Equilibrium
    sll_real64 :: tau0      !-> tau0 = Ti(rpeak)/Te(rpeak)
    sll_real64 :: rho_peak    
    sll_real64 :: kappan   
    sll_real64 :: deltarn  
    sll_real64 :: kappaTi  
    sll_real64 :: deltarTi 
    sll_real64 :: kappaTe  
    sll_real64 :: deltarTe     

    !> Pertubation
    sll_int32  :: perturb_choice
    sll_int32  :: mmode
    sll_int32  :: nmode
    sll_real64 :: eps_perturb   

    !> Numerical parameters
    sll_int32                :: scheme_type
    sll_int32                :: time_loop_type
    sll_real64               :: dt
    sll_int32                :: nb_iter
    type(spline_degree_4d_t) :: spline_degree

    !> Diagnostics
    sll_real64 :: diag2D_step

    !> Tests
    logical :: advec_in_eta1eta2
    logical :: advec_in_eta3
    logical :: advec_in_vpar

    !> Boundary conditions (imposed and not read in the input data file)
    type(boundary_conditions_4d_t) :: bound_cond

  end type input_DK4D_t
  !---------------------------------------------------------------------------

  contains

  !===========================================================================
  !> Read the initial data file : sim4d_DK_hybrid_input.txt
  !---------------------------------------------------------------------------
  subroutine read_input_DK4D( input_data, filename )

    use sll_m_boundary_condition_descriptors, only : sll_p_periodic
    use sll_m_boundary_condition_descriptors, only : sll_p_dirichlet

    type(input_DK4D_t), intent(inout) :: input_data
    character(len=*)  , intent(in)    :: filename

    sll_int32 :: IO_stat
    sll_int32 :: ierr_tests_nml
    sll_int32, parameter :: input_file = 99

    !--> Mesh
    sll_int32  :: num_cells_x1
    sll_int32  :: num_cells_x2
    sll_int32  :: num_cells_x3
    sll_int32  :: num_cells_x4
    sll_real64 :: r_min
    sll_real64 :: r_max
    sll_real64 :: eta3_min
    sll_real64 :: eta3_max
    sll_real64 :: vpar_min
    sll_real64 :: vpar_max
    !--> Equilibrium
    sll_real64 :: tau0
    sll_real64 :: rho_peak    
    sll_real64 :: kappan   
    sll_real64 :: deltarn  
    sll_real64 :: kappaTi  
    sll_real64 :: deltarTi 
    sll_real64 :: kappaTe  
    sll_real64 :: deltarTe     
    !--> Pertubation
    sll_int32  :: perturb_choice
    sll_int32  :: mmode
    sll_int32  :: nmode
    sll_real64 :: eps_perturb   
    !--> Algorithm
    character(len=80) :: scheme_type
    character(len=80) :: time_loop_type
    sll_real64 :: dt
    sll_int32  :: number_iterations
    sll_int32  :: spline_degree_eta1
    sll_int32  :: spline_degree_eta2
    sll_int32  :: spline_degree_eta3
    sll_int32  :: spline_degree_vpar
    !--> Diagnostics
    sll_real64 :: diag2D_step
    !--> Tests
    logical :: advec_in_eta1eta2
    logical :: advec_in_eta3
    logical :: advec_in_vpar

    !*** Reading of the input data ****
    namelist /mesh/ num_cells_x1, num_cells_x2, &
      num_cells_x3, num_cells_x4, &
      r_min, r_max, eta3_min, eta3_max, &
      vpar_min, vpar_max
    namelist /equilibrium/ tau0, rho_peak, kappan, deltarn, &
      kappaTi, deltarTi, kappaTe, deltarTe
    namelist /perturbation/ perturb_choice, mmode, nmode, eps_perturb
    namelist /sim_params/ scheme_type, time_loop_type, dt, number_iterations, &
        spline_degree_eta1, spline_degree_eta2, &
        spline_degree_eta3, spline_degree_vpar
    namelist /diagnostics/ diag2D_step
    namelist /tests/ advec_in_eta1eta2, advec_in_eta3, advec_in_vpar

    !*** Initialization of the boundary conditions ***
    input_data%bound_cond%left_eta1  = sll_p_dirichlet
    input_data%bound_cond%right_eta1 = sll_p_dirichlet
    input_data%bound_cond%left_eta2  = sll_p_periodic
    input_data%bound_cond%right_eta2 = sll_p_periodic
    input_data%bound_cond%left_eta3  = sll_p_periodic
    input_data%bound_cond%right_eta3 = sll_p_periodic
    input_data%bound_cond%left_vpar  = sll_p_dirichlet
    input_data%bound_cond%right_vpar = sll_p_dirichlet


    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_4d_DK_hybrid() failed to open file ', filename
       STOP
    end if
    read(input_file,mesh)
    read(input_file,equilibrium)
    read(input_file,perturbation)
    read(input_file,sim_params)
    read(input_file,diagnostics)
    read(input_file,tests,iostat=ierr_tests_nml)
    close(input_file)

    !--> Mesh
    input_data%nc_x1    = num_cells_x1
    input_data%nc_x2    = num_cells_x2
    input_data%nc_x3    = num_cells_x3
    input_data%nc_x4    = num_cells_x4
    input_data%r_min    = r_min
    input_data%r_max    = r_max
    input_data%eta3_min = eta3_min
    input_data%eta3_max = eta3_max
    input_data%vpar_min = vpar_min
    input_data%vpar_max = vpar_max

    !--> Equilibrium
    input_data%tau0     = tau0
    input_data%rho_peak = rho_peak 
    input_data%kappan   = kappan
    input_data%deltarn  = deltarn
    input_data%kappaTi  = kappaTi
    input_data%deltarTi = deltarTi
    input_data%kappaTe  = kappaTe
    input_data%deltarTe = deltarTe

    !--> Pertubation
    input_data%perturb_choice = perturb_choice
    input_data%mmode          = mmode
    input_data%nmode          = nmode
    input_data%eps_perturb    = eps_perturb

    !--> Algorithm
    select case (scheme_type)
    case("BSL_POLAR")
      input_data%scheme_type = BSL_POLAR
    case("BSL_HYBRID")
      input_data%scheme_type = BSL_HYBRID
    end select
    select case (time_loop_type)
    case("TIME_LOOP_EULER")
      input_data%time_loop_type = TIME_LOOP_EULER
    case ("TIME_LOOP_PREDICTOR_CORRECTOR")
      input_data%time_loop_type = TIME_LOOP_PREDICTOR_CORRECTOR
    end select
    input_data%dt                 = dt
    input_data%nb_iter            = number_iterations
    input_data%spline_degree%eta1 = spline_degree_eta1
    input_data%spline_degree%eta2 = spline_degree_eta2
    input_data%spline_degree%eta3 = spline_degree_eta3
    input_data%spline_degree%vpar = spline_degree_vpar

    !--> Diagnostics
    input_data%diag2D_step   = diag2D_step

    !--> Tests
    if ( ierr_tests_nml==0 ) then
      input_data%advec_in_eta1eta2 = advec_in_eta1eta2
      input_data%advec_in_eta3     = advec_in_eta3
      input_data%advec_in_vpar     = advec_in_vpar
    else
      input_data%advec_in_eta1eta2 = .true.
      input_data%advec_in_eta3     = .true.
      input_data%advec_in_vpar     = .true.
    end if

  end subroutine read_input_DK4D
  !---------------------------------------------------------------------------

end module input_DK4D_module
