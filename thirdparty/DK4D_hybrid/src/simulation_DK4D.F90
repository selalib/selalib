module simulation_DK4D_module

#include "sll_working_precision.h"
#include "sll_memory.h"

  use diagnostics_DK4D_module
  use magnetconf_DK4D_module
  use mesh_DK4D_module
  use equilibrium_DK4D_module
  use fdistribu_DK4D_module
  use field3d_DK4D_module
  use input_DK4D_module
  use elecpot_DK4D_hybrid_module
  use elecpot_DK4D_polar_module
  use QNsolver_DK4D_hybrid_module
  use QNsolver_DK4D_polar_module
  use sll_m_sim_base
  use utils_DK4D_module
  use vlasov_DK4D_hybrid_module
  use vlasov_DK4D_polar_module

  implicit none

#define PRINT_PLOTS 1

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, private ::  elecpot_DK4D_t

    type(elecpot_DK4D_hybrid_t) :: hybrid
    type(elecpot_DK4D_polar_t)  :: polar

  end type Elecpot_DK4D_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, private ::  QNsolver_DK4D_t

    type(Qnsolver_DK4D_hybrid_t) :: hybrid
    type(Qnsolver_DK4D_polar_t)  :: polar

  end type QNsolver_DK4D_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, private ::  vlasov_DK4D_t

    type(vlasov_DK4D_hybrid_t) :: hybrid
    type(vlasov_DK4D_polar_t)  :: polar

  end type vlasov_DK4D_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, extends(sll_c_simulation_base_class) :: sll_simulation_DK4D

    !-> For MPI parallelization
    sll_int32 :: world_size
    sll_int32 :: my_rank
    
    type(input_DK4D_t)       :: input_data
    type(mesh_DK4D_t)        :: mesh4d
    type(magnetconf_DK4D_t)  :: magnetconf
    type(equilibrium_DK4D_t) :: equilibrium
    type(fdistribu_DK4D_t)   :: fdistribu
    type(elecpot_DK4D_t)     :: elec_potential
    type(QNsolver_DK4D_t)    :: QNsolver
    type(vlasov_DK4D_t)      :: vlasov
    type(diagnostics_DK4D_t) :: diagnostics

    !*** Definition of the positions for the cross-section saving ***
    sll_int32 :: ix1_diag
    sll_int32 :: ix2_diag
    sll_int32 :: ieta3_diag
    sll_int32 :: ivpar_diag

  contains

    procedure, pass(sim) :: init_from_file => input_DK4D
    procedure, pass(sim) :: run => run_DK4D

  end type sll_simulation_DK4D
  !---------------------------------------------------------------------------

contains


  !===========================================================================
  !> Constructor for the drift-kinetic 4D simulation
  !---------------------------------------------------------------------------
  subroutine new_DK4D( sim )

    use sll_m_constants

    class(sll_simulation_DK4D), intent(inout) :: sim

    !-> Local variables
    sll_int32 :: Neta1, Neta2, Neta3, Nvpar
    sll_int32 :: Nr, Nx, Ny
    !--> For electrostatic potential boundary conditions and spline degree
    type(boundary_conditions_3d_t) :: elec_pot_bc
    type(spline_degree_3d_t)       :: elec_pot_spline_degree

    !*** Allocation of the 4D mesh ***
    Neta1 = sim%input_data%nc_x1+1 
    Neta2 = sim%input_data%nc_x2+1
    Neta3 = sim%input_data%nc_x3+1 
    Nvpar = sim%input_data%nc_x4+1 

    call new_mesh_DK4D( sim%mesh4d, &
        sim%input_data%nc_x1, &
        sim%input_data%nc_x2, &
        sim%input_data%nc_x3, &
        sim%input_data%nc_x4, & 
        eta1_min=sim%input_data%r_min, &
        eta1_max=sim%input_data%r_max, &
        eta2_min=0.0_f64, &
        eta2_max=2._f64*sll_p_pi, &
        eta3_min=sim%input_data%eta3_min, &
        eta3_max=sim%input_data%eta3_max, &
        vpar_min=sim%input_data%vpar_min, &
        vpar_max=sim%input_data%vpar_max, &
        scheme_type=sim%input_data%scheme_type )

    !*** Allocation of the magnetic configuration ***
    call new_magnetconf_DK4D( sim%magnetconf, sim%mesh4d )

    !*** Allocation of the equilibrium ***
    call new_equilibrium_DK4D( sim%equilibrium, sim%mesh4d )

    !*** Allocation of the distribution function ***
    call new_fdistribu_DK4D( &
        sim%fdistribu, sim%mesh4d, &
        sim%input_data%bound_cond,  &
        sim%input_data%spline_degree )

    !*** Allocation of the electrostatic potential and its derivatives ***
    !---> Initialize the boundary conditions for Phi
    elec_pot_bc%left_eta1  = sim%input_data%bound_cond%left_eta1
    elec_pot_bc%right_eta1 = sim%input_data%bound_cond%right_eta1
    elec_pot_bc%left_eta2  = sim%input_data%bound_cond%left_eta2
    elec_pot_bc%right_eta2 = sim%input_data%bound_cond%right_eta2
    elec_pot_bc%left_eta3  = sim%input_data%bound_cond%left_eta3
    elec_pot_bc%right_eta3 = sim%input_data%bound_cond%right_eta3
    !---> Initialize the spline degree for Phi
    elec_pot_spline_degree%eta1 = sim%input_data%spline_degree%eta1
    elec_pot_spline_degree%eta2 = sim%input_data%spline_degree%eta2
    elec_pot_spline_degree%eta3 = sim%input_data%spline_degree%eta3
    call new_elecpot_DK4D( &
        sim%elec_potential, &
        sim%mesh4d, &
        elec_pot_bc, &
        elec_pot_spline_degree, &
        sim%input_data%scheme_type )

    !*** Allocation of the quasi-neutrality solver ***
    call new_QNsolver_DK4D( &
        sim%QNsolver, &
        sim%mesh4d, &
        elec_pot_bc, &
        elec_pot_spline_degree, &
        sim%input_data%scheme_type )

    !*** Allocation of the valsov solver ***
    call new_vlasov_DK4D( &
        sim%vlasov, &
        sim%mesh4d, &
        sim%input_data%bound_cond, &
        sim%input_data%spline_degree, &
        sim%input_data%scheme_type )

    !*** Allocation for the diagnostics ***
    call new_diagnostics_DK4D( sim%diagnostics )

    !*** Definition of the positions for the cross-section saving ***
    sim%ix1_diag   = int(Neta1/2)
    sim%ix2_diag   = int(Neta2/3)
    sim%ieta3_diag = int(Neta3/4)
    sim%ivpar_diag = int(Nvpar/3)

  end subroutine new_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Delete all what has been used for the drift-kinetic 4D simulation
  !---------------------------------------------------------------------------
  subroutine delete_DK4D( sim )
    
    class(sll_simulation_DK4D), intent(inout) :: sim

    call delete_mesh_DK4D( sim%mesh4d )
    call delete_magnetconf_DK4D( sim%magnetconf )
    call delete_equilibrium_DK4D( sim%equilibrium )
    call delete_fdistribu_DK4D( sim%fdistribu )

    if ( sim%input_data%scheme_type==BSL_HYBRID ) then
      call delete_elecpot_DK4D_hybrid( sim%elec_potential%hybrid )
      call delete_QNsolver_DK4D_hybrid( sim%QNsolver%hybrid )
      call delete_vlasov_DK4D_hybrid( sim%vlasov%hybrid )
    else
      call delete_elecpot_DK4D_polar( sim%elec_potential%polar )
      call delete_QNsolver_DK4D_polar( sim%QNsolver%polar )
      call delete_vlasov_DK4D_polar( sim%vlasov%polar )
    end if


    call delete_diagnostics_DK4D( sim%diagnostics )

  end subroutine delete_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Run drift-kinetic 4D simulation
  !---------------------------------------------------------------------------
  subroutine run_DK4D( sim )
    
    class(sll_simulation_DK4D), intent(inout) :: sim

    sll_int32  :: iter, nb_iter, idiag_num
    sll_real64 :: iter_time, sim_dt

    !*** Allocation of the simulation ***
    call new_DK4D( sim )

    !*** Initialization of the simulation (mesh and equilibrium ) ***
    call initialize_DK4D( sim )

    !*** First iteration of the simulation ***
    call first_step_DK4D( sim )

    !*** Global scheme ***
    idiag_num = 0
    iter_time = 0.0_f64
    sim_dt    = sim%input_data%dt
    nb_iter   = sim%input_data%nb_iter
    do iter = 1,nb_iter
      iter_time = iter_time + sim_dt
      if (sim%my_rank.eq.0) &
          print*,' ===> ITERATION = ',iter, '/',nb_iter,'==> ', &
          iter_time

      !--> Vlasov solving
      call iteration_vlasov_DK4D( &
          sim%vlasov, &
          sim%mesh4d, &
          sim%magnetconf, &
          sim%equilibrium, &
          sim%fdistribu, &
          sim%QNsolver, &
          sim%elec_potential, &
          sim%input_data%dt, &
          sim%input_data%scheme_type, &
          sim%input_data%time_loop_type, &
          sim%input_data%advec_in_eta1eta2, &
          sim%input_data%advec_in_eta3, &
          sim%input_data%advec_in_vpar )

      !--> Quasi-neutrality solving
      call solve_QNsolver_DK4D( &
          sim%QNsolver, &
          sim%mesh4d, &
          sim%magnetconf, &
          sim%equilibrium, &
          sim%fdistribu, &
          l_QNfactorize=.false., &
          elec_potential=sim%elec_potential, &
          scheme_type=sim%input_data%scheme_type )

      !--> Computation of the derivatives of the electrostatic potential
      call derivatives_elecpot_DK4D( &
          sim%elec_potential, &
          sim%mesh4d, &
          sim%input_data%scheme_type )

      !--> Diagnostic computation and saving
      if ( mod(iter,int(sim%input_data%diag2D_step/sim_dt)) == 0) then
        idiag_num = idiag_num + 1        
        call print_DK4D( sim, &
            idiag_num, &
            iter_time )
      end if
    end do

    !*** Deleting of the simulation ***
    call delete_DK4D( sim )

  end subroutine run_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Read the initial data file : sim4d_DK_hybrid_input.txt
  !---------------------------------------------------------------------------
  subroutine input_DK4D( sim, filename )

    use sll_m_collective

    class(sll_simulation_DK4D), intent(inout) :: sim
    character(len=*)          , intent(in)    :: filename

    !*** Initialization of the parallel parameters ***
    sim%world_size = sll_f_get_collective_size( sll_v_world_collective )
    sim%my_rank    = sll_f_get_collective_rank( sll_v_world_collective )

    !*** Reading of the input file ***
    call read_input_DK4D( sim%input_data, filename )

  end subroutine input_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Initialization of the drift-kinetic 4D simulation
  !---------------------------------------------------------------------------
  subroutine initialize_DK4D( sim )

    use sll_m_constants

    class(sll_simulation_DK4D), intent(inout) :: sim

    sll_int32 :: Neta1, Neta2, Neta3, Nvpar

    !*** Initialization of the 4D mesh ***
    call init_mesh_DK4D( sim%mesh4d )
    call print_mesh_DK4D( sim%mesh4d, 'mesh.h5' )

    !*** Initialization of the magnetic configuration ***
    call init_magnetconf_DK4D( sim%magnetconf, sim%mesh4d )
    call print_magnetconf_DK4D( sim%magnetconf, 'magnetconf.h5' )

    !*** Initialization of the equilibrium ***
    call init_equilibrium_DK4D( sim%equilibrium, &
        sim%input_data, sim%mesh4d )
    call print_equilibrium_DK4D( sim%equilibrium, 'equilibrium.h5' )

    !*** Definition of the positions for the cross-section saving ***
    Neta1 = sim%input_data%nc_x1+1 
    Neta2 = sim%input_data%nc_x2+1
    Neta3 = sim%input_data%nc_x3+1 
    Nvpar = sim%input_data%nc_x4+1 

    sim%ix1_diag   = int(Neta1/2)
    sim%ix2_diag   = int(Neta2/3)
    sim%ieta3_diag = int(Neta3/4)
    sim%ivpar_diag = int(Nvpar/3)

  end subroutine initialize_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Run first iteration of the drift-kinetic 4D simulation
  !>   1) Initialization of the distribution function 
  !>   2) Solving of the quasi-neutrality equation
  !---------------------------------------------------------------------------
  subroutine first_step_DK4D( sim )
    
    use sll_m_timer

    class(sll_simulation_DK4D), intent(inout) :: sim

    !-> Local variables
    !--> For timers
    type(sll_t_time_mark) :: t0, t1 
    double precision    :: elaps_time
     
    call sll_s_set_time_mark(t0)

    !*** Initialization of the distribution function ***
    call init_fdistribu_DK4D( sim%fdistribu, &
      sim%input_data, sim%mesh4d, sim%equilibrium )
    
    !*** Solving of the quasi-neutrality equation ***
    call solve_QNsolver_DK4D( &
        sim%QNsolver, &
        sim%mesh4d, &
        sim%magnetconf, &
        sim%equilibrium, &
        sim%fdistribu, &
        l_QNfactorize=.true., &
        elec_potential=sim%elec_potential, &
        scheme_type=sim%input_data%scheme_type )

    !*** Computation of the electric field (required for the advections) ***
    call derivatives_elecpot_DK4D( &
        sim%elec_potential, &
        sim%mesh4d, &
        sim%input_data%scheme_type )

    !*** Diagnostic computation and saving at initial time ***
    call print_DK4D( sim, &
        idiag_num=0, &
        diag_time=0.0_f64 )

    call sll_s_set_time_mark(t1)
    elaps_time = sll_f_time_elapsed_between(t0,t1)
    if ( sim%my_rank.eq.0 ) &
      print*, ' First time = ', elaps_time

  end subroutine first_step_DK4D
  !---------------------------------------------------------------------------


  !************************************************************
  !  ELECTROSTATIC POTENTIAL AND ITS DERIVATIVES
  !************************************************************
  !===========================================================================
  !>  Constructor for the electrostatic potential and its derivatives 
  !>   depending on the BSL scheme type
  !---------------------------------------------------------------------------
  subroutine new_elecpot_DK4D( &
      elecpot, &
      mesh4d, &
      bound_cond, &
      spline_degree, &
      scheme_type )

    type(elecpot_DK4D_t)          , intent(inout) :: elecpot
    type(mesh_DK4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_3d_t), intent(in)    :: bound_cond
    type(spline_degree_3d_t)      , intent(in)    :: spline_degree
    sll_int32                     , intent(in)    :: scheme_type

    select case ( scheme_type )   
    case ( BSL_HYBRID )
      call new_elecpot_DK4D_hybrid( &
          elecpot%hybrid, &
          mesh4d, &
          bound_cond, &
          spline_degree )

    case ( BSL_POLAR )
      call new_elecpot_DK4D_polar( &
          elecpot%polar, &
          mesh4d, &
          bound_cond, &
          spline_degree )
    end select

  end subroutine new_elecpot_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> \todo : Move this computation in elec_potential -> derivatives
  !---------------------------------------------------------------------------
  subroutine derivatives_elecpot_DK4D( &
      elecpot, &
      mesh4d, &
      scheme_type )

    type(elecpot_DK4D_t), intent(inout) :: elecpot
    type(mesh_DK4D_t)   , intent(in)    :: mesh4d
    sll_int32           , intent(in)    :: scheme_type

    select case ( scheme_type )
    case ( BSL_HYBRID )
      call derivatives_elecpot_DK4D_hybrid( &
          elecpot%hybrid, &
          mesh4d )

    Case ( BSL_POLAR )
      call derivatives_elecpot_DK4D_polar( &
          elecpot%polar, &
          mesh4d )
    end select

  end subroutine derivatives_elecpot_DK4D
  !---------------------------------------------------------------------------


  !************************************************************
  !  QUASI-NEUTRALITY SOLVING
  !************************************************************

  !===========================================================================
  !>  Constructor for the quasi-neutrality solver depending on the BSL
  !>   scheme type
  !---------------------------------------------------------------------------
  subroutine new_QNsolver_DK4D( &
      QNsolver, &
      mesh4d, &
      bound_cond, &
      spline_degree, & 
      scheme_type )

    type(QNsolver_DK4D_t)         , intent(inout) :: QNsolver
    type(mesh_DK4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_3d_t), intent(in)    :: bound_cond
    type(spline_degree_3d_t)      , intent(in)    :: spline_degree
    sll_int32                     , intent(in)    :: scheme_type

    select case ( scheme_type )   
    case ( BSL_HYBRID )
      call new_QNsolver_DK4D_hybrid( &
          QNsolver%hybrid, &
          mesh4d, &
          bound_cond, &
          spline_degree )

    case ( BSL_POLAR )
      call new_QNsolver_DK4D_polar( &
          QNsolver%polar, &
          mesh4d, &
          bound_cond, &
          spline_degree )
    end select

  end subroutine new_QNsolver_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !>  Solving of the quasi-neutrality solver depending on the BSL
  !>   scheme type
  !---------------------------------------------------------------------------
  subroutine solve_QNsolver_DK4D( &
      QNsolver, &
      mesh4d, &
      magnetconf, &
      equilibrium, &
      fdistribu, &
      l_QNfactorize, &
      elec_potential, &
      scheme_type )

    type(QNsolver_DK4D_t)   , intent(inout) :: QNsolver
    type(mesh_DK4D_t)       , intent(in)    :: mesh4d
    type(magnetconf_DK4D_t) , intent(in)    :: magnetconf
    type(equilibrium_DK4D_t), intent(in)    :: equilibrium
    type(fdistribu_DK4D_t)  , intent(in)    :: fdistribu
    logical                 , intent(in)    :: l_QNfactorize
    type(elecpot_DK4D_t)    , intent(inout) :: elec_potential
    sll_int32               , intent(in)    :: scheme_type

    select case ( scheme_type )
    case ( BSL_HYBRID )
      !---> Initialization of the quasi-neutrality solver 
      call init_QNsolver_DK4D_hybrid( &
          QNsolver%hybrid, &
          mesh4d, &
          magnetconf, &
          equilibrium )
      !\todo: Understand the fact to call the factorization several
      !  times is wrong (probably terms of the matrix are not put
      !  equal to 0)
      if (l_QNfactorize) then
        call factorize_QNsolver_DK4D_hybrid( QNsolver%hybrid )
      end if

      !---> Solving
      call solve_QNsolver_DK4D_hybrid( &
          QNsolver%hybrid, &
          mesh4d, &
          equilibrium, &
          fdistribu, &
          elec_potential%hybrid%Phi )
    
    case ( BSL_POLAR )
      !---> Initialization of the quasi-neutrality solver 
      call init_QNsolver_DK4D_polar( &
          QNsolver%polar, &
          mesh4d, &
          equilibrium )

      !---> Solving
      call solve_QNsolver_DK4D_polar( &
          QNsolver%polar, &
          mesh4d, &
          equilibrium, &
          fdistribu, &
          elec_potential%polar%Phi )

    end select
      
  end subroutine solve_QNsolver_DK4D
  !---------------------------------------------------------------------------


  !************************************************************
  !  VLASOV SOLVING
  !************************************************************

  !===========================================================================
  !>  Constructor for the vlasov simulation depending on the BSL
  !>   scheme type
  !---------------------------------------------------------------------------
  subroutine new_vlasov_DK4D( &
      vlasov, &
      mesh4d, &
      bound_cond, &
      spline_degree, &
      scheme_type )
    
    type(vlasov_DK4D_t)           , intent(inout) :: vlasov
    type(mesh_DK4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_4d_t), intent(in)    :: bound_cond
    type(spline_degree_4d_t)      , intent(in)    :: spline_degree
    sll_int32                     , intent(in)    :: scheme_type

    select case ( scheme_type )
    case ( BSL_HYBRID )
      call new_vlasov_DK4D_hybrid( &
          vlasov%hybrid, &
          mesh4d, &
          bound_cond, &
          spline_degree )

    case ( BSL_POLAR )
      call new_vlasov_DK4D_polar( &
          vlasov%polar, &
          mesh4d, &
          bound_cond, &
          spline_degree )
    end select

  end subroutine new_vlasov_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !---------------------------------------------------------------------------
  subroutine iteration_vlasov_DK4D( &
      vlasov, &
      mesh4d, &
      magnetconf, &
      equilibrium, &
      fdistribu, &
      QNsolver, &
      elec_potential, &
      dt, &
      scheme_type, &
      time_loop_type, &
      advec_in_eta1eta2, &
      advec_in_eta3, &
      advec_in_vpar )

    type(vlasov_DK4D_t)     , intent(inout) :: vlasov
    type(mesh_DK4D_t)       , intent(in)    :: mesh4d
    type(magnetconf_DK4D_t) , intent(in)    :: magnetconf
    type(equilibrium_DK4D_t), intent(inout) :: equilibrium
    type(fdistribu_DK4D_t)  , intent(inout) :: fdistribu
    type(QNsolver_DK4D_t)   , intent(inout) :: QNsolver
    type(elecpot_DK4D_t)    , intent(inout) :: elec_potential
    sll_real64              , intent(in)    :: dt
    sll_int32               , intent(in)    :: scheme_type
    sll_int32               , intent(in)    :: time_loop_type
    logical                 , intent(in)    :: advec_in_eta1eta2
    logical                 , intent(in)    :: advec_in_eta3
    logical                 , intent(in)    :: advec_in_vpar

    !--> Local variable used for predictor-corrector
    sll_int32  :: ierr
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3, loc4d_sz_x4
    sll_int32  :: ix1, ix2, ix3, ivpar
    sll_real64, dimension(:,:,:,:), pointer :: f4d_store_seqx3x4
    !--> Local variable used to choose the advection which are performed
    sll_real64 :: advec_eta1eta2_OK
    sll_real64 :: advec_eta3_OK
    sll_real64 :: advec_vpar_OK

    !--> Initialization of the multiply coefficients used to decide if the 
    !-->  advections are taken into account or not
    advec_eta1eta2_OK = 1._f64
    advec_eta3_OK     = 1._f64
    advec_vpar_OK     = 1._f64
    if ( .not.advec_in_eta1eta2 ) advec_eta1eta2_OK = 0._f64 
    if ( .not.advec_in_eta3 ) advec_eta3_OK = 0._f64 
    if ( .not.advec_in_vpar ) advec_vpar_OK = 0._f64 

    !*** Choice between the different scheme for the global scheme iteration ***
    if ( time_loop_type == TIME_LOOP_PREDICTOR_CORRECTOR ) then
      call sll_o_compute_local_sizes( fdistribu%layout4d_seqx3x4, &
          loc4d_sz_x1, &
          loc4d_sz_x2, &
          loc4d_sz_x3, &
          loc4d_sz_x4 )    
      SLL_ALLOCATE( f4d_store_seqx3x4(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4), ierr )
    end if

    select case ( scheme_type )

    case ( BSL_HYBRID )
      if ( time_loop_type == TIME_LOOP_EULER ) then
        call iteration_vlasov_DK4D_hybrid( &
            vlasov%hybrid, &
            fdistribu, &
            mesh4d, &
            magnetconf, &
            elec_potential%hybrid, &
            dt, & 
            advec_vpar_step1_coef     = 0.5_f64*advec_vpar_OK, &
            advec_eta3_step1_coef     = 0.5_f64*advec_eta3_OK, &
            advec_eta1eta2_step1_coef = 1._f64*advec_eta1eta2_OK, &
            advec_eta3_step2_coef     = 0.5_f64*advec_eta3_OK, &
            advec_vpar_step2_coef     = 0.5_f64*advec_vpar_OK )
      else if ( time_loop_type == TIME_LOOP_PREDICTOR_CORRECTOR ) then
        !***  PREDICTOR ***
        do ivpar = 1,loc4d_sz_x4
          do ix3 = 1,loc4d_sz_x3
            do ix2 = 1,loc4d_sz_x2
              do ix1 = 1,loc4d_sz_x1
                f4d_store_seqx3x4(ix1,ix2,ix3,ivpar) = &
                    fdistribu%val4d_seqx3x4(ix1,ix2,ix3,ivpar)
              end do
            end do
          end do
        end do
        !-> vlasov solving
        call iteration_vlasov_DK4D_hybrid( &
            vlasov%hybrid, &
            fdistribu, &
            mesh4d, &
            magnetconf, &
            elec_potential%hybrid, &
            dt, & 
            advec_vpar_step1_coef     = 0.5_f64 * advec_vpar_OK, &
            advec_eta3_step1_coef     = 0.5_f64 * advec_eta3_OK, &
            advec_eta1eta2_step1_coef = 0.5_f64  * advec_eta1eta2_OK )
        !--> quasi-neutrality solving
        call solve_QNsolver_DK4D_hybrid( &
            QNsolver%hybrid, &
            mesh4d, &
            equilibrium, &
            fdistribu, &
            elec_potential%hybrid%Phi )
        !--> Computation of the derivative of the electrostatic potential
        call derivatives_elecpot_DK4D_hybrid( &
            elec_potential%hybrid, &
            mesh4d )
        
        !***  CORRECTOR ***
        do ivpar = 1,loc4d_sz_x4
          do ix3 = 1,loc4d_sz_x3
            do ix2 = 1,loc4d_sz_x2
              do ix1 = 1,loc4d_sz_x1
                fdistribu%val4d_seqx3x4(ix1,ix2,ix3,ivpar) = &
                    f4d_store_seqx3x4(ix1,ix2,ix3,ivpar)
              end do
            end do
          end do
        end do
        call iteration_vlasov_DK4D_hybrid( &
            vlasov%hybrid, &
            fdistribu, &
            mesh4d, &
            magnetconf, &
            elec_potential%hybrid, &
            dt, & 
            advec_vpar_step1_coef     = 0.5_f64*advec_vpar_OK, &
            advec_eta3_step1_coef     = 0.5_f64*advec_eta3_OK, &
            advec_eta1eta2_step1_coef = 1._f64*advec_eta1eta2_OK, &
            advec_eta3_step2_coef     = 0.5_f64*advec_eta3_OK, &
            advec_vpar_step2_coef     = 0.5_f64*advec_vpar_OK )
      end if

    case ( BSL_POLAR )
      if ( time_loop_type == TIME_LOOP_EULER ) then
        call iteration_vlasov_DK4D_polar( &
            vlasov%polar, &
            fdistribu, &
            mesh4d, &
            magnetconf, &
            elec_potential%polar, &
            dt, &
            advec_vpar_step1_coef     = 0.5_f64*advec_vpar_OK, &
            advec_eta3_step1_coef     = 0.5_f64*advec_eta3_OK, &
            advec_eta1eta2_step1_coef = 1._f64*advec_eta1eta2_OK, &
            advec_eta3_step2_coef     = 0.5_f64*advec_eta3_OK, &
            advec_vpar_step2_coef     = 0.5_f64*advec_vpar_OK )

      else if ( time_loop_type == TIME_LOOP_PREDICTOR_CORRECTOR ) then
        !***  PREDICTOR ***
        do ivpar = 1,loc4d_sz_x4
          do ix3 = 1,loc4d_sz_x3
            do ix2 = 1,loc4d_sz_x2
              do ix1 = 1,loc4d_sz_x1
                f4d_store_seqx3x4(ix1,ix2,ix3,ivpar) = &
                    fdistribu%val4d_seqx3x4(ix1,ix2,ix3,ivpar)
              end do
            end do
          end do
        end do
        !--> Vlasov solving
        call iteration_vlasov_DK4D_polar( &
            vlasov%polar, &
            fdistribu, &
            mesh4d, &
            magnetconf, &
            elec_potential%polar, &
            dt, & 
            advec_vpar_step1_coef     = 0.5_f64*advec_vpar_OK, &
            advec_eta3_step1_coef     = 0.5_f64*advec_eta3_OK, &
            advec_eta1eta2_step1_coef = 0.5_f64*advec_eta1eta2_OK )
        !--> Quasi-neutraliy solving 
        call solve_QNsolver_DK4D_polar( &
            QNsolver%polar, &
            mesh4d, &
            equilibrium, &
            fdistribu, &
            elec_potential%polar%Phi )
        !--> Computation of the derivative of the electrostatic potential
        call derivatives_elecpot_DK4D_polar( &
            elec_potential%polar, &
            mesh4d )

        !***  CORRECTOR ***
        do ivpar = 1,loc4d_sz_x4
          do ix3 = 1,loc4d_sz_x3
            do ix2 = 1,loc4d_sz_x2
              do ix1 = 1,loc4d_sz_x1
                fdistribu%val4d_seqx3x4(ix1,ix2,ix3,ivpar) = &
                    f4d_store_seqx3x4(ix1,ix2,ix3,ivpar)
              end do
            end do
          end do
        end do
        call iteration_vlasov_DK4D_polar( &
            vlasov%polar, &
            fdistribu, &
            mesh4d, &
            magnetconf, &
            elec_potential%polar, &
            dt, & 
            advec_vpar_step1_coef     = 0.5_f64*advec_vpar_OK, &
            advec_eta3_step1_coef     = 0.5_f64*advec_eta3_OK, &
            advec_eta1eta2_step1_coef = 1._f64*advec_eta1eta2_OK, &
            advec_eta3_step2_coef     = 0.5_f64*advec_eta3_OK, &
            advec_vpar_step2_coef     = 0.5_f64*advec_vpar_OK )
      else
        print*,' Problem with the choice of the time loop scheme'
      end if
    end select

    if ( time_loop_type == TIME_LOOP_PREDICTOR_CORRECTOR ) then
      SLL_DEALLOCATE( f4d_store_seqx3x4, ierr )
    end if

  end subroutine iteration_vlasov_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Diagnostics: Printing the distribution function, 
  !>  the electrostatic potential and all other diagnostics in 
  !>  HDF5 format
  !---------------------------------------------------------------------------
  subroutine print_DK4D( sim, &
      idiag_num, &
      diag_time )

    class(sll_simulation_DK4D), intent(inout) :: sim
    sll_int32 , intent(in) :: idiag_num
    sll_real64, intent(in) :: diag_time

!\todo : understand why the pointer does not work
!VG!    type(field3d_DK4D_t), pointer :: Phi
!VG!    type(field3d_DK4D_t), pointer :: elec_field_eta1_3d
!VG!    type(field3d_DK4D_t), pointer :: elec_field_eta2_3d
!VG!    type(field3d_DK4D_t), pointer :: elec_field_eta3_3d
!VG!
!VG!    select case ( sim%input_data%scheme_type )
!VG!    case ( BSL_HYBRID )
!VG!      Phi                => sim%elec_potential%hybrid%Phi
!VG!      elec_field_eta1_3d => sim%vlasov%hybrid%elec_field_eta1_3d
!VG!      elec_field_eta2_3d => sim%vlasov%hybrid%elec_field_eta2_3d
!VG!      elec_field_eta3_3d => sim%vlasov%hybrid%elec_field_eta3_3d
!VG!
!VG!    case ( BSL_POLAR )
!VG!      Phi                => sim%elec_potential%polar%Phi
!VG!      elec_field_eta1_3d => sim%vlasov%polar%elec_field_eta1_3d
!VG!      elec_field_eta2_3d => sim%vlasov%polar%elec_field_eta2_3d
!VG!      elec_field_eta3_3d => sim%vlasov%polar%elec_field_eta3_3d
!VG!    end select

    !--> Saving of cross-sections of the distribution function
    call print_fdistribu_DK4D( &
        sim%fdistribu, &
        sim%ix1_diag, &
        sim%ix2_diag, &
        sim%ieta3_diag, &
        sim%ivpar_diag, &
        idiag_num ) 

    select case ( sim%input_data%scheme_type )
    case ( BSL_HYBRID )
      !-> Saving of cross-sections of the electrostatic potential
      call print_field3d_DK4D( &
          sim%elec_potential%hybrid%Phi, &
          sim%ix1_diag, &
          sim%ix2_diag, &
          sim%ieta3_diag, &
          idiag_num )

      !-> Saving of cross-sections of the derivatives of the 
      !->  electrostatic potential
      call print_field3d_DK4D( &
          sim%elec_potential%hybrid%dPhi_deta1, &
          sim%ix1_diag, &
          sim%ix2_diag, &
          sim%ieta3_diag, &
          idiag_num )
      call print_field3d_DK4D( &
          sim%elec_potential%hybrid%dPhi_deta2, &
          sim%ix1_diag, &
          sim%ix2_diag, &
          sim%ieta3_diag, &
          idiag_num )
      call print_field3d_DK4D( &
          sim%elec_potential%hybrid%dPhi_deta3, &
          sim%ix1_diag, &
          sim%ix2_diag, &
          sim%ieta3_diag, &
          idiag_num )

      !-> Saving of 0D and 1D quantities
      call print_diagnostics_DK4D( sim%diagnostics, &
          idiag_num, &
          diag_time, &
          sim%mesh4d, &
          sim%equilibrium, &
          sim%fdistribu, &
          sim%elec_potential%hybrid%Phi )

    case ( BSL_POLAR )
      !-> Saving of cross-sections of the electrostatic potential
      call print_field3d_DK4D( &
          sim%elec_potential%polar%Phi, &
          sim%ix1_diag, &
          sim%ix2_diag, &
          sim%ieta3_diag, &
          idiag_num )

      !-> Saving of cross-sections of the derivatives of the 
      !->  electrostatic potential
      call print_field3d_DK4D( &
          sim%elec_potential%polar%dPhi_deta1, &
          sim%ix1_diag, &
          sim%ix2_diag, &
          sim%ieta3_diag, &
          idiag_num )
      call print_field3d_DK4D( &
          sim%elec_potential%polar%dPhi_deta2, &
          sim%ix1_diag, &
          sim%ix2_diag, &
          sim%ieta3_diag, &
          idiag_num )
      call print_field3d_DK4D( &
          sim%elec_potential%polar%dPhi_deta3, &
          sim%ix1_diag, &
          sim%ix2_diag, &
          sim%ieta3_diag, &
          idiag_num )

      !-> Saving of 0D and 1D quantities
      call print_diagnostics_DK4D( sim%diagnostics, &
          idiag_num, &
          diag_time, &
          sim%mesh4d, &
          sim%equilibrium, &
          sim%fdistribu, &
          sim%elec_potential%polar%Phi )
    end select

  end subroutine print_DK4D
!---------------------------------------------------------------------------

end module simulation_DK4D_module
!---------------------------------------------------------------------------
