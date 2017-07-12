!===========================================================================
!> Input data for 4D Vlasov-Poisson hybrid simulation
!>
!> \date 2015-02-26
!> \author V. Grandgirard, A. Back
!---------------------------------------------------------------------------
module input_VP4D_module
#include "sll_working_precision.h"

  use utils_VP4D_module

  implicit none
  
  !> For finite element choice
  sll_int32, parameter, public :: BSPLINES     = 1
  sll_int32, parameter, public :: BOX_SPLINES  = 2

  !> For mapping type choice
  sll_int32, parameter, public :: ANALYTICAL = 1
  sll_int32, parameter, public :: CAID_FILE  = 2

  !> For choice of analytical mapping formula
  sll_int32, parameter, public :: VOID      = -1
  sll_int32, parameter, public :: IDENTITY  =  1
  sll_int32, parameter, public :: COLELLA   =  2

  !> For time loop choice
  sll_int32, parameter, public :: TIME_LOOP_EULER  = 1
  sll_int32, parameter, public :: TIME_LOOP_PREDICTOR_CORRECTOR = 2

  !> For advec2D in (eta1,eta2) scheme choice
  sll_int32, parameter :: SLL_EULER      = 0
  sll_int32, parameter :: SLL_RUNGEKUTTA = 1

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: input_VP4D_t

    !> Mesh parameters
    sll_int32          :: nc_x1
    sll_int32          :: nc_x2
    sll_int32          :: nc_x3
    sll_int32          :: nc_x4
    sll_real64         :: x1_min
    sll_real64         :: x1_max
    sll_real64         :: x2_min
    sll_real64         :: x2_max
    sll_real64         :: vx_min
    sll_real64         :: vx_max
    sll_real64         :: vy_min
    sll_real64         :: vy_max
    sll_int32          :: finite_element_type
    sll_int32          :: mapping_type
    sll_int32          :: analytical_formula
    character(len=80 ) :: mapping_filename
    sll_real64         :: colella_coeff

    !> Pertubation
    sll_real64 :: kx1
    sll_real64 :: eps_perturb   

    !> Numerical parameters
    sll_int32                :: time_loop_type
    sll_int32                :: advec2D_scheme
    sll_real64               :: dt
    sll_int32                :: nb_iter
    type(spline_degree_4d_t) :: spline_degree

    !> Diagnostics
    sll_real64 :: diag2D_step

    !> Tests
    logical :: advec_in_eta1eta2
    logical :: advec_in_vx
    logical :: advec_in_vy

    !> Boundary conditions (imposed and not read in the input data file)
    type(boundary_conditions_4d_t) :: bound_cond

  end type input_VP4D_t
  !---------------------------------------------------------------------------

  contains

  !===========================================================================
  !> Read the initial data file : sim4d_VP_hybrid_input.txt
  !---------------------------------------------------------------------------
  subroutine read_input_VP4D( input_data, filename )

    use sll_m_boundary_condition_descriptors, only : sll_p_periodic, sll_p_dirichlet

    type(input_VP4D_t), intent(inout) :: input_data
    character(len=*)  , intent(in)    :: filename

    sll_int32 :: IO_stat
    sll_int32 :: ierr_tests_nml
    sll_int32, parameter :: input_file = 99

    !--> Mesh
    sll_int32         :: num_cells_x1
    sll_int32         :: num_cells_x2
    sll_int32         :: num_cells_x3
    sll_int32         :: num_cells_x4
    sll_real64        :: x1_min
    sll_real64        :: x1_max
    sll_real64        :: x2_min
    sll_real64        :: x2_max
    sll_real64        :: vx_min
    sll_real64        :: vx_max
    sll_real64        :: vy_min
    sll_real64        :: vy_max
    character(len=80) :: finite_element_type
    character(len=80) :: mapping_type
    character(len=80) :: analytical_formula
    character(len=80) :: mapping_filename
    sll_real64        :: colella_coeff
    !--> Pertubation
    sll_real64 :: kx1
    sll_real64 :: eps_perturb   
    !--> Algorithm
    character(len=80) :: time_loop_type
    character(len=80) :: advec2D_scheme
    sll_real64 :: dt
    sll_int32  :: number_iterations
    sll_int32  :: spline_degree_eta1
    sll_int32  :: spline_degree_eta2
    sll_int32  :: spline_degree_vx
    sll_int32  :: spline_degree_vy
    !--> Diagnostics
    sll_real64 :: diag2D_step
    !--> Tests
    logical :: advec_in_eta1eta2
    logical :: advec_in_vx
    logical :: advec_in_vy

    !*** Reading of the input data ****
    namelist /mesh/ num_cells_x1, num_cells_x2, &
      num_cells_x3, num_cells_x4, &
      x1_min, x1_max, x2_min, x2_max, &
      vx_min, vx_max, vy_min, vy_max, &
      finite_element_type, mapping_type, &
      analytical_formula, mapping_filename, colella_coeff
    namelist /perturbation/ kx1, eps_perturb
    namelist /sim_params/ time_loop_type, advec2D_scheme, &
        dt, number_iterations, &
        spline_degree_eta1, spline_degree_eta2, &
        spline_degree_vx, spline_degree_vy
    namelist /diagnostics/ diag2D_step
    namelist /tests/ advec_in_eta1eta2, advec_in_vx, advec_in_vy

    !*** Initialization of the boundary conditions ***
    input_data%bound_cond%left_eta1  = sll_p_periodic
    input_data%bound_cond%right_eta1 = sll_p_periodic
    input_data%bound_cond%left_eta2  = sll_p_periodic
    input_data%bound_cond%right_eta2 = sll_p_periodic
    input_data%bound_cond%left_vx    = sll_p_dirichlet
    input_data%bound_cond%right_vx   = sll_p_dirichlet
    input_data%bound_cond%left_vy    = sll_p_dirichlet
    input_data%bound_cond%right_vy   = sll_p_dirichlet

    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_4d_VP_hybrid() failed to open file ', filename
       STOP
    end if
    read(input_file,mesh)
    read(input_file,perturbation)
    read(input_file,sim_params)
    read(input_file,diagnostics)
    read(input_file,tests,iostat=ierr_tests_nml)
    close(input_file)

    !--> Mesh
    input_data%nc_x1  = num_cells_x1
    input_data%nc_x2  = num_cells_x2
    input_data%nc_x3  = num_cells_x3
    input_data%nc_x4  = num_cells_x4
    input_data%x1_min = x1_min
    input_data%x1_max = x1_max
    input_data%x2_min = x2_min
    input_data%x2_max = x2_max
    input_data%vx_min = vx_min
    input_data%vx_max = vx_max
    input_data%vy_min = vy_min
    input_data%vy_max = vy_max

    select case (finite_element_type)
    case("BSPLINES")
      input_data%finite_element_type = BSPLINES
    end select

    select case (mapping_type)
    case("ANALYTICAL")
      input_data%mapping_type = ANALYTICAL
      select case (analytical_formula)
      case("IDENTITY")
        input_data%analytical_formula = IDENTITY
      case("COLELLA")
        input_data%analytical_formula = COLELLA
        input_data%colella_coeff      = colella_coeff
      end select
      input_data%mapping_filename = ""

    case("CAID_FILE")
      input_data%mapping_type       = CAID_FILE
      input_data%mapping_filename   = mapping_filename
      input_data%analytical_formula = VOID
    end select

    !--> Pertubation
    input_data%kx1         = kx1
    input_data%eps_perturb = eps_perturb

    !--> Algorithm
    select case (time_loop_type)
    case("TIME_LOOP_EULER")
      input_data%time_loop_type = TIME_LOOP_EULER
    case ("TIME_LOOP_PREDICTOR_CORRECTOR")
      input_data%time_loop_type = TIME_LOOP_PREDICTOR_CORRECTOR
    end select
    select case (advec2D_scheme)
    case("EULER")
      input_data%advec2D_scheme = SLL_EULER
    case("RUNGEKUTTA")
      input_data%advec2D_scheme = SLL_RUNGEKUTTA
    end select
    input_data%dt                 = dt
    input_data%nb_iter            = number_iterations
    input_data%spline_degree%eta1 = spline_degree_eta1
    input_data%spline_degree%eta2 = spline_degree_eta2
    input_data%spline_degree%vx   = spline_degree_vx
    input_data%spline_degree%vy   = spline_degree_vy

    !--> Diagnostics
    input_data%diag2D_step   = diag2D_step

    !--> Tests
    if ( ierr_tests_nml==0 ) then
      input_data%advec_in_eta1eta2 = advec_in_eta1eta2
      input_data%advec_in_vx       = advec_in_vx
      input_data%advec_in_vy       = advec_in_vy
    else
      input_data%advec_in_eta1eta2 = .true.
      input_data%advec_in_vx       = .true.
      input_data%advec_in_vy       = .true.
    end if

  end subroutine read_input_VP4D
  !---------------------------------------------------------------------------

end module input_VP4D_module
