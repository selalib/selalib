!===========================================================================
!> Simulation for 4D Vlasov-Poisson hybrid simulation
!>  f(x,y,vx,vy)
!>
!> \date 2015-02-26
!> \author V. Grandgirard, A. Back
!---------------------------------------------------------------------------
module simulation_VP4D_module

#include "sll_working_precision.h"
#include "sll_memory.h"

  use diagnostics_VP4D_module
  use equilibrium_VP4D_module
  use elecpot_Bsplines_VP4D_module
  use fdistribu_VP4D_module
  use field2d_VP4D_module
  use input_VP4D_module
  use mesh_VP4D_module
  use QNsolver_Bsplines_VP4D_module
  use sll_m_sim_base
  use utils_VP4D_module
  use vlasov_VP4D_module

  implicit none

#define PRINT_PLOTS 1

  !===========================================================================
  !> Rk: Structure prepared to add BoxSpline after
  !---------------------------------------------------------------------------
  type, private ::  elecpot_VP4D_t

    type(elecpot_Bsplines_VP4D_t) :: Bsplines
!VG!    type(elecpot_BOXsplines_VP4D_t) :: BOXsplines

  end type Elecpot_VP4D_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, private ::  QNsolver_VP4D_t

    type(Qnsolver_Bsplines_VP4D_t)   :: Bsplines
!VG!    type(Qnsolver_BOXsplines_VP4D_t) :: BOXsplines

  end type QNsolver_VP4D_t
  !---------------------------------------------------------------------------


  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, extends(sll_c_simulation_base_class) :: sll_simulation_VP4D

    !-> For MPI parallelization
    sll_int32 :: world_size
    sll_int32 :: my_rank
    
    type(input_VP4D_t)       :: input_data
    type(mesh_VP4D_t)        :: mesh4d
    type(equilibrium_VP4D_t) :: equilibrium
    type(fdistribu_VP4D_t)   :: fdistribu
    type(elecpot_VP4D_t)     :: elec_potential
    type(QNsolver_VP4D_t)    :: QNsolver
    type(diagnostics_VP4D_t) :: diagnostics

    !*** Definition of the positions for the cross-section saving ***
    sll_int32 :: ix1_diag
    sll_int32 :: ix2_diag
    sll_int32 :: ivx_diag
    sll_int32 :: ivy_diag

  contains

    procedure, pass(sim) :: init_from_file => input_VP4D
    procedure, pass(sim) :: run => run_VP4D

  end type sll_simulation_VP4D
  !---------------------------------------------------------------------------

contains


  !===========================================================================
  !> Constructor for the drift-kinetic 4D simulation
  !---------------------------------------------------------------------------
  subroutine new_VP4D( sim )

    use sll_m_constants

    class(sll_simulation_VP4D), intent(inout) :: sim

    !-> Local variables
    sll_int32 :: Neta1, Neta2, Nvx, Nvy

    !--> For electrostatic potential boundary conditions and spline degree
    type(boundary_conditions_2d_t) :: elec_pot_bc
    type(spline_degree_2d_t)       :: elec_pot_spline_degree

    !*** Allocation of the 4D mesh ***
    !BECAREFUL : Modify to put eta1_min(max) and eta2_min(max)
    call new_mesh_VP4D( sim%mesh4d, &
        sim%input_data%nc_x1, &
        sim%input_data%nc_x2, &
        sim%input_data%nc_x3, &
        sim%input_data%nc_x4, & 
        eta1_min=sim%input_data%x1_min, &
        eta1_max=sim%input_data%x1_max, &
        eta2_min=sim%input_data%x2_min, &
        eta2_max=sim%input_data%x2_max, &
        vx_min=sim%input_data%vx_min, &
        vx_max=sim%input_data%vx_max, &
        vy_min=sim%input_data%vy_min, &
        vy_max=sim%input_data%vy_max, &
        mapping_type=sim%input_data%mapping_type, &
        analytical_formula=sim%input_data%analytical_formula, &
        mapping_filename=sim%input_data%mapping_filename, &
        colella_coeff=sim%input_data%colella_coeff )

    !*** Allocation of the equilibrium ***
    call new_equilibrium_VP4D( sim%equilibrium, sim%mesh4d )

    !*** Allocation of the distribution function ***
    call new_fdistribu_VP4D( &
        sim%fdistribu, sim%mesh4d, &
        sim%input_data%bound_cond,  &
        sim%input_data%spline_degree )

    !*** Allocation of the electrostatic potential and its derivatives ***
    !---> Initialize the boundary conditions for Phi
    elec_pot_bc%left_eta1  = sim%input_data%bound_cond%left_eta1
    elec_pot_bc%right_eta1 = sim%input_data%bound_cond%right_eta1
    elec_pot_bc%left_eta2  = sim%input_data%bound_cond%left_eta2
    elec_pot_bc%right_eta2 = sim%input_data%bound_cond%right_eta2
    !---> Initialize the spline degree for Phi
    elec_pot_spline_degree%eta1 = sim%input_data%spline_degree%eta1
    elec_pot_spline_degree%eta2 = sim%input_data%spline_degree%eta2
    call new_elecpot_VP4D( &
        sim%elec_potential, &
        sim%mesh4d, &
        elec_pot_bc, &
        sim%input_data%finite_element_type, &
        elec_pot_spline_degree )

    !*** Allocation of the quasi-neutrality solver ***
    call new_QNsolver_VP4D( &
        sim%QNsolver, &
        sim%mesh4d, &
        elec_pot_bc, &
        sim%input_data%finite_element_type, &
        elec_pot_spline_degree )

    !*** Allocation for the diagnostics ***
    call new_diagnostics_VP4D( sim%diagnostics )

    !*** Definition of the positions for the cross-section saving ***
    Neta1 = sim%mesh4d%eta1_eta2_mesh2d%num_cells1+1 
    Neta2 = sim%mesh4d%eta1_eta2_mesh2d%num_cells2+1
    Nvx   = sim%mesh4d%vx_vy_mesh2d%num_cells1+1 
    Nvy   = sim%mesh4d%vx_vy_mesh2d%num_cells2+1 
    sim%ix1_diag = int((Neta1-1)/2)
    sim%ix2_diag = int((Neta2-1)/3)
    sim%ivx_diag = int((Nvx-1)/4)
    sim%ivy_diag = int((Nvy-1)/3)

  end subroutine new_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Delete all what has been used for the drift-kinetic 4D simulation
  !---------------------------------------------------------------------------
  subroutine delete_VP4D( sim )
    
    class(sll_simulation_VP4D), intent(inout) :: sim

    call delete_mesh_VP4D( sim%mesh4d )
    call delete_equilibrium_VP4D( sim%equilibrium )
    call delete_fdistribu_VP4D( sim%fdistribu )

    call delete_elecpot_Bsplines_VP4D( sim%elec_potential%Bsplines )
    call delete_QNsolver_Bsplines_VP4D( sim%QNsolver%Bsplines )

    call delete_diagnostics_VP4D( sim%diagnostics )

  end subroutine delete_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Run drift-kinetic 4D simulation
  !---------------------------------------------------------------------------
  subroutine run_VP4D( sim )
    
    class(sll_simulation_VP4D), intent(inout) :: sim

    sll_int32  :: iter, nb_iter, idiag_num
    sll_real64 :: iter_time, sim_dt

    !*** Allocation of the simulation ***
    call new_VP4D( sim )

    !*** Initialization of the simulation (mesh and equilibrium ) ***
    call initialize_VP4D( sim )

    !*** First iteration of the simulation ***
    call first_step_VP4D( sim )

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
      call iteration_vlasov_VP4D( &
          sim%mesh4d, &
          sim%equilibrium, &
          sim%fdistribu, &
          sim%QNsolver, &
          sim%elec_potential, &
          sim%input_data%dt, &
          sim%input_data%time_loop_type, &
          sim%input_data%advec2D_scheme, &
          sim%input_data%advec_in_eta1eta2, &
          sim%input_data%advec_in_vx, &
          sim%input_data%advec_in_vy )

      !--> Quasi-neutrality solving
      call solve_QNsolver_VP4D( &
          sim%QNsolver, &
          sim%mesh4d, &
          sim%equilibrium, &
          sim%fdistribu, &
          l_QNfactorize=.false., &
          elec_potential=sim%elec_potential, &
          finite_element_type=sim%input_data%finite_element_type )

      !--> Computation of the derivatives of the electrostatic potential
      call derivatives_elecpot_VP4D( &
          sim%elec_potential, &
          sim%mesh4d, &
          sim%input_data%finite_element_type )

      !--> Diagnostic computation and saving
      if ( mod(iter,int(sim%input_data%diag2D_step/sim_dt)) == 0) then
        idiag_num = idiag_num + 1        
        call print_VP4D( sim, &
            idiag_num, &
            iter_time )
      end if
    end do

    !*** Deleting of the simulation ***
    call delete_VP4D( sim )

  end subroutine run_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Read the initial data file 
  !---------------------------------------------------------------------------
  subroutine input_VP4D( sim, filename )

    use sll_m_collective

    class(sll_simulation_VP4D), intent(inout) :: sim
    character(len=*)          , intent(in)    :: filename

    !*** Initialization of the parallel parameters ***
    sim%world_size = sll_f_get_collective_size( sll_v_world_collective )
    sim%my_rank    = sll_f_get_collective_rank( sll_v_world_collective )

    !*** Reading of the input file ***
    call read_input_VP4D( sim%input_data, filename )

  end subroutine input_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Initialization of the drift-kinetic 4D simulation
  !---------------------------------------------------------------------------
  subroutine initialize_VP4D( sim )

    use sll_m_constants

    class(sll_simulation_VP4D), intent(inout) :: sim

    !*** Initialization of the 4D mesh ***
    call init_mesh_VP4D( sim%mesh4d )
    call print_mesh_VP4D( sim%mesh4d, 'mesh.h5' )

    !*** Initialization of the equilibrium ***
    call init_equilibrium_VP4D( sim%equilibrium, sim%mesh4d )
    call print_equilibrium_VP4D( sim%equilibrium, 'equilibrium.h5' )

  end subroutine initialize_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Run first iteration of the drift-kinetic 4D simulation
  !>   1) Initialization of the distribution function 
  !>   2) Solving of the quasi-neutrality equation
  !---------------------------------------------------------------------------
  subroutine first_step_VP4D( sim )
    
    use sll_m_timer

    class(sll_simulation_VP4D), intent(inout) :: sim

    !-> Local variables
    !--> For timers
    type(sll_t_time_mark) :: t0, t1 
    double precision    :: elaps_time
     
    call sll_s_set_time_mark(t0)

    !*** Initialization of the distribution function ***
    call init_fdistribu_VP4D( sim%fdistribu, &
      sim%input_data, sim%mesh4d, sim%equilibrium )
    
    !*** Solving of the quasi-neutrality equation ***
    call solve_QNsolver_VP4D( &
        sim%QNsolver, &
        sim%mesh4d, &
        sim%equilibrium, &
        sim%fdistribu, &
        l_QNfactorize=.true., &
        elec_potential=sim%elec_potential, &
        finite_element_type=sim%input_data%finite_element_type )

    !*** Computation of the electric field (required for the advections) ***
    call derivatives_elecpot_VP4D( &
        sim%elec_potential, &
        sim%mesh4d, &
        sim%input_data%finite_element_type )

    !*** Diagnostic computation and saving at initial time ***
    call print_VP4D( sim, &
        idiag_num=0, &
        diag_time=0.0_f64 )

    call sll_s_set_time_mark(t1)
    elaps_time = sll_f_time_elapsed_between(t0,t1)
    if ( sim%my_rank.eq.0 ) &
      print*, ' First time = ', elaps_time

  end subroutine first_step_VP4D
  !---------------------------------------------------------------------------


  !************************************************************
  !  ELECTROSTATIC POTENTIAL AND ITS DERIVATIVES
  !************************************************************
  !===========================================================================
  !>  Constructor for the electrostatic potential and its derivatives 
  !>   depending on the BSL scheme type
  !---------------------------------------------------------------------------
  subroutine new_elecpot_VP4D( &
      elecpot, &
      mesh4d, &
      bound_cond, &
      finite_element_type, &
      spline_degree )

    type(elecpot_VP4D_t)           , intent(inout) :: elecpot
    type(mesh_VP4D_t)              , intent(in)    :: mesh4d
    type(boundary_conditions_2d_t) , intent(in)    :: bound_cond
    sll_int32                      , intent(in)    :: finite_element_type
    type(spline_degree_2d_t), optional, intent(in) :: spline_degree

    select case ( finite_element_type )   
    case ( BSPLINES )
      call new_elecpot_Bsplines_VP4D( &
          elecpot%Bsplines, &
          mesh4d, &
          bound_cond, &
          spline_degree )
    end select

  end subroutine new_elecpot_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> \todo : Move this computation in elec_potential -> derivatives
  !---------------------------------------------------------------------------
  subroutine derivatives_elecpot_VP4D( &
      elecpot, &
      mesh4d, &
      finite_element_type )

    type(elecpot_VP4D_t), intent(inout) :: elecpot
    type(mesh_VP4D_t)   , intent(in)    :: mesh4d
    sll_int32           , intent(in)    :: finite_element_type

    select case ( finite_element_type )
    case ( BSPLINES )
      call derivatives_elecpot_Bsplines_VP4D( &
          elecpot%Bsplines, &
          mesh4d )
    end select

  end subroutine derivatives_elecpot_VP4D
  !---------------------------------------------------------------------------


  !************************************************************
  !  QUASI-NEUTRALITY SOLVING
  !************************************************************

  !===========================================================================
  !>  Constructor for the quasi-neutrality solver depending on the BSL
  !>   scheme type
  !---------------------------------------------------------------------------
  subroutine new_QNsolver_VP4D( &
      QNsolver, &
      mesh4d, &
      bound_cond, &
      finite_element_type, & 
      spline_degree )

    type(QNsolver_VP4D_t)         , intent(inout) :: QNsolver
    type(mesh_VP4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_2d_t), intent(in)    :: bound_cond
    sll_int32                     , intent(in)    :: finite_element_type
    type(spline_degree_2d_t)      , intent(in)    :: spline_degree

    select case ( finite_element_type )   
    case ( BSPLINES )
      call new_QNsolver_Bsplines_VP4D( &
          QNsolver%Bsplines, &
          mesh4d, &
          bound_cond, &
          spline_degree )
    end select

  end subroutine new_QNsolver_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !>  Solving of the quasi-neutrality solver depending on the BSL
  !>   scheme type
  !---------------------------------------------------------------------------
  subroutine solve_QNsolver_VP4D( &
      QNsolver, &
      mesh4d, &
      equilibrium, &
      fdistribu, &
      l_QNfactorize, &
      elec_potential, &
      finite_element_type )

    type(QNsolver_VP4D_t)   , intent(inout) :: QNsolver
    type(mesh_VP4D_t)       , intent(in)    :: mesh4d
    type(equilibrium_VP4D_t), intent(in)    :: equilibrium
    type(fdistribu_VP4D_t)  , intent(in)    :: fdistribu
    logical                 , intent(in)    :: l_QNfactorize
    type(elecpot_VP4D_t)    , intent(inout) :: elec_potential
    sll_int32               , intent(in)    :: finite_element_type

    select case ( finite_element_type )
    case ( BSPLINES )
      !\todo: Understand the fact to call the factorization several
      !  times is wrong (probably terms of the matrix are not put
      !  equal to 0)
      if (l_QNfactorize) then
        call factorize_QNsolver_Bsplines_VP4D( QNsolver%Bsplines )
      end if

      !---> Solving
      call solve_QNsolver_Bsplines_VP4D( &
          QNsolver%Bsplines, &
          mesh4d, &
          equilibrium, &
          fdistribu, &
          elec_potential%Bsplines%Phi )    
    end select
      
  end subroutine solve_QNsolver_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !---------------------------------------------------------------------------
  subroutine iteration_vlasov_VP4D( &
      mesh4d, &
      equilibrium, &
      fdistribu, &
      QNsolver, &
      elec_potential, &
      dt, &
      time_loop_type, &
      advec2D_scheme, &
      advec_in_eta1eta2, &
      advec_in_vx, &
      advec_in_vy )

    type(mesh_VP4D_t)       , intent(in)    :: mesh4d
    type(equilibrium_VP4D_t), intent(inout) :: equilibrium
    type(fdistribu_VP4D_t)  , intent(inout) :: fdistribu
    type(QNsolver_VP4D_t)   , intent(inout) :: QNsolver
    type(elecpot_VP4D_t)    , intent(inout) :: elec_potential
    sll_real64              , intent(in)    :: dt
    sll_int32               , intent(in)    :: time_loop_type
    sll_int32               , intent(in)    :: advec2D_scheme
    logical                 , intent(in)    :: advec_in_eta1eta2
    logical                 , intent(in)    :: advec_in_vx
    logical                 , intent(in)    :: advec_in_vy

    !--> Local variable used for predictor-corrector
    sll_int32  :: ierr
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3, loc4d_sz_x4
    sll_int32  :: ix1, ix2, ix3, ivy
    sll_real64, dimension(:,:,:,:), pointer :: f4d_store_seqx3x4
    !--> Local variable used to choose the advection which are performed
    sll_real64 :: advec_eta1eta2_OK
    sll_real64 :: advec_vx_OK
    sll_real64 :: advec_vy_OK

    !--> Initialization of the multiply coefficients used to decide if the 
    !-->  advections are taken into account or not
    advec_eta1eta2_OK = 1._f64
    advec_vx_OK     = 1._f64
    advec_vy_OK     = 1._f64
    if ( .not.advec_in_eta1eta2 ) advec_eta1eta2_OK = 0._f64 
    if ( .not.advec_in_vx ) advec_vx_OK = 0._f64 
    if ( .not.advec_in_vy ) advec_vy_OK = 0._f64 

    !*** Choice between the different scheme for the global scheme iteration ***
    select case ( time_loop_type )

    case ( TIME_LOOP_EULER )
      call iteration_vlasov_Bsplines_VP4D( &
          fdistribu, &
          mesh4d, &
          elec_potential%Bsplines, &
          dt, & 
          advec2D_scheme, &
          advec_vy_step1_coef     = 0.5_f64*advec_vy_OK, &
          advec_vx_step1_coef     = 0.5_f64*advec_vx_OK, &
          advec_eta1eta2_step1_coef = 1._f64*advec_eta1eta2_OK, &
          advec_vx_step2_coef     = 0.5_f64*advec_vx_OK, &
          advec_vy_step2_coef     = 0.5_f64*advec_vy_OK )

    case ( TIME_LOOP_PREDICTOR_CORRECTOR )
      call sll_o_compute_local_sizes( fdistribu%layout4d_seqx3x4, &
          loc4d_sz_x1, &
          loc4d_sz_x2, &
          loc4d_sz_x3, &
          loc4d_sz_x4 )    
      SLL_ALLOCATE( f4d_store_seqx3x4(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4), ierr )

      !***  PREDICTOR ***
      do ivy = 1,loc4d_sz_x4
        do ix3 = 1,loc4d_sz_x3
          do ix2 = 1,loc4d_sz_x2
            do ix1 = 1,loc4d_sz_x1
              f4d_store_seqx3x4(ix1,ix2,ix3,ivy) = &
                  fdistribu%val4d_seqx3x4(ix1,ix2,ix3,ivy)
            end do
          end do
        end do
      end do
      !-> vlasov solving
      call iteration_vlasov_Bsplines_VP4D( &
          fdistribu, &
          mesh4d, &
          elec_potential%Bsplines, &
          dt, & 
          advec2D_scheme, &
          advec_vy_step1_coef     = 0.5_f64 * advec_vy_OK, &
          advec_vx_step1_coef     = 0.5_f64 * advec_vx_OK, &
          advec_eta1eta2_step1_coef = 0.5_f64  * advec_eta1eta2_OK )
      !--> quasi-neutrality solving
      call solve_QNsolver_Bsplines_VP4D( &
          QNsolver%Bsplines, &
          mesh4d, &
          equilibrium, &
          fdistribu, &
          elec_potential%Bsplines%Phi )
      !--> Computation of the derivative of the electrostatic potential
      call derivatives_elecpot_Bsplines_VP4D( &
          elec_potential%Bsplines, &
          mesh4d )

      !***  CORRECTOR ***
      do ivy = 1,loc4d_sz_x4
        do ix3 = 1,loc4d_sz_x3
          do ix2 = 1,loc4d_sz_x2
            do ix1 = 1,loc4d_sz_x1
              fdistribu%val4d_seqx3x4(ix1,ix2,ix3,ivy) = &
                  f4d_store_seqx3x4(ix1,ix2,ix3,ivy)
            end do
          end do
        end do
      end do
      call iteration_vlasov_Bsplines_VP4D( &
          fdistribu, &
          mesh4d, &
          elec_potential%Bsplines, &
          dt, & 
          advec2D_scheme, &
          advec_vy_step1_coef     = 0.5_f64*advec_vy_OK, &
          advec_vx_step1_coef     = 0.5_f64*advec_vx_OK, &
          advec_eta1eta2_step1_coef = 1._f64*advec_eta1eta2_OK, &
          advec_vx_step2_coef     = 0.5_f64*advec_vx_OK, &
          advec_vy_step2_coef     = 0.5_f64*advec_vy_OK )

      SLL_DEALLOCATE( f4d_store_seqx3x4, ierr )

    case default
      print*,"Scheme in time = ", time_loop_type, " not treated."
      stop

    end select

  end subroutine iteration_vlasov_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Diagnostics: Printing the distribution function, 
  !>  the electrostatic potential and all other diagnostics in 
  !>  HDF5 format
  !---------------------------------------------------------------------------
  subroutine print_VP4D( sim, &
      idiag_num, &
      diag_time )

    class(sll_simulation_VP4D), intent(inout) :: sim
    sll_int32 , intent(in) :: idiag_num
    sll_real64, intent(in) :: diag_time

    !--> Saving of cross-sections of the distribution function
    call print_fdistribu_VP4D( &
        sim%fdistribu, &
        sim%mesh4d, &
        sim%ix1_diag, &
        sim%ix2_diag, &
        sim%ivx_diag, &
        sim%ivy_diag, &
        idiag_num ) 

    select case ( sim%input_data%finite_element_type )
    case ( BSPLINES )
      !-> Saving of cross-sections of the electrostatic potential
      call print_field2d_VP4D( &
          sim%elec_potential%Bsplines%Phi, &
          idiag_num )

      !-> Saving of cross-sections of the derivatives of the 
      !->  electrostatic potential
      call print_field2d_VP4D( &
          sim%elec_potential%Bsplines%dPhi_deta1, &
          idiag_num )
      call print_field2d_VP4D( &
          sim%elec_potential%Bsplines%dPhi_deta2, &
          idiag_num )

      !-> Saving of 0D and 1D quantities
      call print_diagnostics_VP4D( sim%diagnostics, &
          idiag_num, &
          diag_time, &
          sim%mesh4d, &
          sim%equilibrium, &
          sim%fdistribu, &
          sim%elec_potential%Bsplines%Phi, &
          sim%elec_potential%Bsplines%dPhi_deta1, &
          sim%elec_potential%Bsplines%dPhi_deta2 )
    end select

  end subroutine print_VP4D
!---------------------------------------------------------------------------

end module simulation_VP4D_module
!---------------------------------------------------------------------------
