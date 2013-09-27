! The idea of this simulation is to merge the functionalities of the qns-based
! simulation and the older, cartesian 4d simulation for the purposes of 
! debugging/understanding the behavior of the QNS one. Once this objective is
! fulfilled, this simulation can be deleted.

module sll_simulation_4d_DK_hybrid_module

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"

  use sll_collective
  use sll_remapper
!VG!  use sll_constants
!VG!  use sll_cubic_spline_interpolator_2d
!VG!  use sll_poisson_2d_periodic_cartesian_par
!VG!  use sll_cubic_spline_interpolator_1d
  use sll_simulation_base
  use sll_logical_meshes
  use sll_parallel_array_initializer_module
!VG!  use sll_coordinate_transformation_2d_base_module
!VG!  use sll_gnuplot_parallel
!VG!  use sll_general_coordinate_qn_solver_module
!VG!  use sll_module_scalar_field_2d_base
!VG!  use sll_module_scalar_field_2d_alternative

  implicit none

#define PRINT_PLOTS 1
  type, extends(sll_simulation_base_class) :: sll_simulation_4d_DK_hybrid
     ! Parallel environment parameters
     sll_int32  :: world_size
     sll_int32  :: my_rank
     sll_int32  :: power2 ! 2^power2 = number of processes available
     ! Processor mesh sizes
     sll_int32  :: nproc_x1
     sll_int32  :: nproc_x2
     sll_int32  :: nproc_x3
     sll_int32  :: nproc_x4 
     ! Mesh parameters
     sll_int32  :: nc_x1
     sll_int32  :: nc_x2
     sll_int32  :: nc_x3
     sll_int32  :: nc_x4
     ! Physics/numerical parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
     sll_int32  :: spline_degree_eta1, spline_degree_eta2
     sll_int32  :: spline_degree_eta3, spline_degree_eta4
     !--> Equilibrium
     sll_real64 :: tau0      !-> tau0 = Ti(rpeak)/Te(rpeak)
     sll_real64 :: rpeak    
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

     !--> 4D logical mesh (eta1,eta2,eta3,v)
     type(sll_logical_mesh_4d), pointer :: logical_mesh4d
     sll_real64, dimension(:), pointer  :: eta1_grid
     sll_real64, dimension(:), pointer  :: eta2_grid
     sll_real64, dimension(:), pointer  :: eta3_grid
     sll_real64, dimension(:), pointer  :: vpar_grid

     !--> For boundary conditions
     sll_int32 :: bc_left_eta1
     sll_int32 :: bc_right_eta1
     sll_int32 :: bc_left_eta2
     sll_int32 :: bc_right_eta2
     sll_int32 :: bc_left_eta3
     sll_int32 :: bc_right_eta3
     sll_int32 :: bc_left_vpar
     sll_int32 :: bc_right_vpar

     !--> Parallel decomposition of the mesh
     sll_int32 :: loc_sz_x1
     sll_int32 :: loc_sz_x2
     sll_int32 :: loc_sz_x3
     sll_int32 :: loc_sz_x4
     type(layout_4D), pointer :: sequential_x1x2
     type(layout_4D), pointer :: sequential_x3x4

     !--> Density and temperature profiles
     sll_real64, dimension(:), pointer :: n0
     sll_real64, dimension(:), pointer :: Ti
     sll_real64, dimension(:), pointer :: Te

     !--> Equilibrium distribution function
     sll_real64, dimension(:,:), pointer :: feq_2d

     !--> 4D distribution function 
     !----> sequential in (x1,x2) and parallel in (x3,x4)
     sll_real64, dimension(:,:,:,:), pointer :: f4d_x1x2 
     !----> parallel in (x1,x2) and sequential in (x3,x4) 
     sll_real64, dimension(:,:,:,:), pointer :: f4d_x3x4

   contains
     procedure, pass(sim) :: run => run_4d_DK_hybrid
     procedure, pass(sim) :: init_from_file => init_4d_DK_hybrid
  end type sll_simulation_4d_DK_hybrid

  interface delete
     module procedure delete_4d_DK_hybrid
  end interface delete

  interface initialize
     module procedure initialize_4d_DK_hybrid
  end interface initialize

contains

  !----------------------------------------------------------
  !  Read the initial data file : sim4d_DK_hybrid_input.txt
  !----------------------------------------------------------
  subroutine init_4d_DK_hybrid( sim, filename)
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    character(len=*)                  , intent(in)    :: filename

    sll_int32            :: IO_stat
    sll_int32, parameter :: input_file = 99
    !--> Mesh
    sll_int32  :: num_cells_x1
    sll_int32  :: num_cells_x2
    sll_int32  :: num_cells_x3
    sll_int32  :: num_cells_x4
    !--> Algorithm
    sll_real64 :: dt
    sll_int32  :: number_iterations
    sll_int32  :: spline_degree
    !--> Equilibrium
    sll_real64 :: tau0
    sll_real64 :: rpeak    
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

    namelist /grid_dims/ num_cells_x1, num_cells_x2, &
      num_cells_x3, num_cells_x4
    namelist /sim_params/ dt, number_iterations, spline_degree
    namelist /equilibrium/ tau0, rpeak, kappan, deltarn, &
      kappaTi, deltarTi, kappaTe, deltarTe
    namelist /perturbation/ perturb_choice, mmode, nmode, eps_perturb

    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_4d_DK_hybrid() failed to open file ', filename
       STOP
    end if
    read(input_file,grid_dims)
    read(input_file,sim_params)
    read(input_file,equilibrium)
    read(input_file,perturbation)
    close(input_file)

    !--> Mesh
    sim%nc_x1 = num_cells_x1
    sim%nc_x2 = num_cells_x2
    sim%nc_x3 = num_cells_x3
    sim%nc_x4 = num_cells_x4
    !--> Algorithm
    sim%dt                 = dt
    sim%num_iterations     = number_iterations
    sim%spline_degree_eta1 = spline_degree
    sim%spline_degree_eta2 = spline_degree
    sim%spline_degree_eta3 = spline_degree
    sim%spline_degree_eta4 = spline_degree
    !--> Equilibrium
    sim%rpeak    = rpeak 
    sim%kappan   = kappan
    sim%deltarn  = deltarn
    sim%kappaTi  = kappaTi
    sim%deltarTi = deltarTi
    sim%kappaTe  = kappaTe
    sim%deltarTe = deltarTe
    !--> Pertubation
    sim%perturb_choice = perturb_choice
    sim%mmode          = mmode
    sim%nmode          = nmode
    sim%eps_perturb    = eps_perturb
  end subroutine init_4d_DK_hybrid


  !---------------------------------------- 
  ! sech = cosh^-1 definition
  !---------------------------------------- 
  function sech(x)
    sll_real64, intent(in) :: x   
    sll_real64             :: sech
    
    sech = 1._f64/cosh(x)
  end function sech  


  !---------------------------------------- -----------------------
  ! Initialization of the radial density profiles 
  !---------------------------------------- -----------------------
  subroutine init_n0_eta1(sim,rpeak,kappan,deltarn,n0_eta1min,n0_1d)
    class(sll_simulation_4d_DK_hybrid), intent(in) :: sim
    sll_real64, intent(in) :: rpeak
    sll_real64, intent(in) :: kappan
    sll_real64, intent(in) :: deltarn
    sll_real64, intent(in) :: n0_eta1min
    sll_real64, dimension(:), intent(inout) :: n0_1d

    sll_int32  :: ierr, ieta1, Neta1
    sll_real64 :: delta_eta1, etai, etai_mid
    sll_real64 :: tmp, inv_Ln0
    sll_real64 :: n0norm_tmp

    Neta1      = sim%logical_mesh4d%num_cells1+1
    delta_eta1 = sim%logical_mesh4d%delta_eta1
    
    !*** compute ns0 solution of :                           ***
    !***  2/(n0(r)+n0(r-1))*(n0(r)-n0(r-1))/dr               ***
    !***                  = -(1/Ln0)*cosh^-2(r-rpeak/deltar) ***
    inv_Ln0  = kappan           !??? (1/Ln) = kappa_n0/R
    n0_1d(1) = n0_eta1min 
    do ieta1 = 2,Neta1
      etai     = sim%eta1_grid(ieta1-1)
      etai_mid = etai + delta_eta1*0.5_f64
      tmp      = -inv_Ln0 * &
        sech((etai_mid-rpeak)/deltarn)**2
      tmp          = 0.5_f64*delta_eta1*tmp
      n0_1d(ieta1) = (1._f64+tmp)/(1._f64-tmp)*n0_1d(ieta1-1)
    enddo

    !*** normalisation of the density at int(n0(r)rdr)/int(rdr) ***
    ! -> computation of int(n0(r)rdr)
    n0norm_tmp = 0._f64
    do ieta1 = 2,Neta1-1
      n0norm_tmp = n0norm_tmp + n0_1d(ieta1)*sim%eta1_grid(ieta1)
    enddo
    n0norm_tmp = n0norm_tmp + 0.5_f64 * &
      (n0_1d(1)*sim%eta1_grid(1) + n0_1d(Neta1)*sim%eta1_grid(Neta1))
    ! -> division by int(rdr)
    n0norm_tmp = n0norm_tmp*2._f64*delta_eta1 / & 
      (sim%eta1_grid(Neta1)**2-sim%eta1_grid(1)**2)

    n0_1d(1:Neta1) = n0_1d(1:Neta1)/n0norm_tmp
  end subroutine init_n0_eta1


  !---------------------------------------------------------------
  ! Initialization of the radial temperature profiles 
  !---------------------------------------------------------------
  subroutine init_T_eta1(sim,rpeak,kappaT,deltarT,T_eta1min,T_scal,T_1d)
    class(sll_simulation_4d_DK_hybrid), intent(in) :: sim
    sll_real64, intent(in) :: rpeak
    sll_real64, intent(in) :: kappaT
    sll_real64, intent(in) :: deltarT
    sll_real64, intent(in) :: T_eta1min
    sll_real64, intent(in) :: T_scal
    sll_real64, dimension(:), intent(inout) :: T_1d
    
    sll_int32  :: ierr, ieta1, Neta1
    sll_real64 :: delta_eta1, etai, etai_mid
    sll_real64 :: tmp, inv_LT
    sll_real64 :: w0, w1, Tnorm_tmp

    Neta1      = sim%logical_mesh4d%num_cells1+1
    delta_eta1 = sim%logical_mesh4d%delta_eta1
    
    !*** compute ns0 solution of :                           ***
    !***  2/(n0(r)+n0(r-1))*(n0(r)-n0(r-1))/dr               ***
    !***                  = -(1/Ln0)*cosh^-2(r-rpeak/deltar) ***
    inv_LT  = kappaT           !??? (1/Ln) = kappa_n0/R
    T_1d(1) = T_eta1min
    do ieta1 = 2,Neta1
      etai     = sim%eta1_grid(ieta1-1)
      etai_mid = etai + delta_eta1*0.5_f64
      tmp      = -inv_LT * &
        sech((etai_mid-rpeak)/deltarT)**2
      tmp           = 0.5_f64*delta_eta1*tmp
      T_1d(ieta1) = (1._f64+tmp)/(1._f64-tmp)*T_1d(ieta1-1)
    enddo

    !*** normalisation of the temperature to 1 at r=rpeak ***
    ieta1         = int((rpeak-sim%eta1_grid(1))/delta_eta1)
    w1            = (rpeak-sim%eta1_grid(ieta1))/delta_eta1
    w0            = 1._f64-w1
    Tnorm_tmp     = w0*T_1d(ieta1)+w1*T_1d(ieta1+1)
    T_1d(1:Neta1) = (T_1d(1:Neta1)/Tnorm_tmp)/T_scal
  end subroutine init_T_eta1


  !-------------------------------------------------------------------------------------
  ! Computation of the equilibrium distribution function for
  !  drift-kinetic 4D simulation
  !   feq(eta1,vpar) = n0(eta1)/(2*pi*Ti(eta1))**(1/2) * exp(-0.5*eta4**2/Ti(eta1)
  !  where n0 and Ti are respectively the initial density and temperature profiles
  !-------------------------------------------------------------------------------------
  function compute_feq_val( eta1, vpar, n0_eta1, Ti_eta1) &
    result(val)

    sll_real64             :: val      ! sll_DK_initializer_4d
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: vpar
    sll_real64, intent(in) :: n0_eta1
    sll_real64, intent(in) :: Ti_eta1

    val = n0_eta1*sqrt(2._f64*sll_pi*Ti_eta1)*exp(-0.5_f64*vpar**2/Ti_eta1)
  end function compute_feq_val


  !----------------------------------------------------
  ! Initialization of the 2D array for the equilibrium
  !  distribution function feq(eta1,vpar)
  !----------------------------------------------------
  subroutine init_fequilibrium(Neta1,Nvpar,eta1_grid,vpar_grid,n0_1d,Ti_1d,feq_2d)
    sll_int32, intent(in) :: Neta1
    sll_int32, intent(in) :: Nvpar
    sll_real64, dimension(:)  , intent(in)  :: eta1_grid
    sll_real64, dimension(:)  , intent(in)  :: vpar_grid
    sll_real64, dimension(:)  , intent(in)  :: n0_1d
    sll_real64, dimension(:)  , intent(in)  :: Ti_1d
    sll_real64, dimension(:,:), intent(out) :: feq_2d
    
    sll_int32  :: ieta1, ivpar
    sll_real64 :: eta1, vpar, n0_eta1, Ti_eta1

    do ivpar = 1,Nvpar
      vpar = vpar_grid(ivpar)
      do ieta1 = 1,Neta1
        eta1    = eta1_grid(ieta1)
        n0_eta1 = n0_1d(ieta1)
        Ti_eta1 = Ti_1d(ieta1)
        feq_2d(ieta1,ivpar) = compute_feq_val(eta1,vpar,n0_eta1,Ti_eta1)
      end do
    end do
  end subroutine init_fequilibrium


  !----------------------------------------------------
  ! Initialization of the drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine initialize_4d_DK_hybrid( sim, logical_mesh4d)
!VG!    use sll_hdf5_io_parallel, only: sll_hdf5_file_create, &
!VG!      sll_hdf5_write_array, sll_hdf5_file_close
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    type(sll_logical_mesh_4d)         , pointer       :: logical_mesh4d

    sll_int32  :: ierr
    sll_int32  :: ieta1, Neta1
    sll_real64 :: delta_eta1, eta1_min
    sll_int32  :: ivpar, Nvpar
    sll_real64 :: delta_vpar, vpar_min
    sll_real64 :: n0_eta1min, Ti_eta1min, Te_eta1min, Ti_scal, Te_scal

!VG!    !--> For initial profile HDF5 saving
!VG!    integer                     :: file_err
!VG!    sll_int32                   :: file_id
!VG!    character(len=9), parameter :: filename_prof = "init_prof.h5"

    !--> Logical mesh initialization
    sim%logical_mesh4d => logical_mesh4d
    !----> Init eta1 grid
    Neta1      = sim%logical_mesh4d%num_cells1+1
    SLL_ALLOCATE(sim%eta1_grid(Neta1),ierr)
    delta_eta1 = sim%logical_mesh4d%delta_eta1
    eta1_min   = sim%logical_mesh4d%eta1_min
    do ieta1 = 1,Neta1
      sim%eta1_grid(ieta1) = eta1_min + (ieta1-1)*delta_eta1
    end do
    !----> Init vpar grid
    Nvpar      = sim%logical_mesh4d%num_cells1+1
    SLL_ALLOCATE(sim%vpar_grid(Nvpar),ierr)
    delta_vpar = sim%logical_mesh4d%delta_eta4
    vpar_min   = sim%logical_mesh4d%eta4_min
    do ivpar = 1,Nvpar
      sim%vpar_grid(ivpar) = vpar_min + (ivpar-1)*delta_vpar
    end do

    !--> Initialization of the boundary conditions
    sim%bc_left_eta1  = SLL_DIRICHLET
    sim%bc_right_eta1 = SLL_DIRICHLET
    sim%bc_left_eta2  = SLL_PERIODIC
    sim%bc_right_eta2 = SLL_PERIODIC
    sim%bc_left_eta3  = SLL_PERIODIC
    sim%bc_right_eta3 = SLL_PERIODIC
    sim%bc_left_vpar  = SLL_DIRICHLET
    sim%bc_right_vpar = SLL_DIRICHLET

    !--> Radial profile initialisation
    SLL_ALLOCATE(sim%n0(sim%nc_x1+1),ierr)
    SLL_ALLOCATE(sim%Ti(sim%nc_x1+1),ierr)
    SLL_ALLOCATE(sim%Te(sim%nc_x1+1),ierr)
    sim%n0(:) = 0.0_f64
    sim%Ti(:) = 0.0_f64
    sim%Te(:) = 0.0_f64
    n0_eta1min = 10._f64**19
    Ti_eta1min = 1.3_f64
    Ti_scal    = 1._f64
    call init_n0_eta1(sim,sim%rpeak,sim%kappan,sim%deltarn,n0_eta1min,sim%n0)
    call init_T_eta1(sim,sim%rpeak,sim%kappaTi,sim%deltarTi,Ti_eta1min,Ti_scal,sim%Ti)
    Te_eta1min = sim%Ti(1)
    Te_scal    = sim%tau0
    call init_T_eta1(sim,sim%rpeak,sim%kappaTe,sim%deltarTe,Te_eta1min,Te_scal,sim%Te)
    
    !--> Initialization of the equilibrium distribution function
    Neta1 = sim%logical_mesh4d%num_cells1+1
    Nvpar = sim%logical_mesh4d%num_cells4+1
    SLL_ALLOCATE(sim%feq_2d(sim%nc_x1+1,sim%nc_x4+1),ierr)
    call init_fequilibrium(Neta1,Nvpar,sim%eta1_grid,sim%vpar_grid, &
      sim%n0,sim%Ti,sim%feq_2d)

!VG!    !*** Saving of the radial profiles in HDF5 file ***
!VG!    call sll_hdf5_file_create(trim(filename_prof),file_id,file_err)
!VG!    call sll_hdf5_write_array(file_id,sim%n0,'n0',file_err)
!VG!    call sll_hdf5_file_close(file_id,file_err)
  end subroutine initialize_4d_DK_hybrid


  !----------------------------------------------------
  ! Run drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine run_4d_DK_hybrid( sim )
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    !-->
    sll_real64, dimension(:), pointer :: send_buf
    sll_real64, dimension(:), pointer :: recv_buf
    sll_int32 , dimension(:), pointer :: recv_sz
    sll_int32 , dimension(:), pointer :: disps ! for allgatherv operation

    sll_int32 :: ierr
    sll_int32 :: itemp

    !--> Initialization of the 4D parallel layout
    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)
    SLL_ALLOCATE(recv_sz(sim%world_size),ierr)
    SLL_ALLOCATE(disps(sim%world_size),ierr)

    ! layout for sequential operations in x3 and x4. Make an even split for
    ! x1 and x2, or as close as even if the power of 2 is odd. This should 
    ! be packaged in some sort of routine and set up at initialization time.
    sim%power2 = int(log(real(sim%world_size))/log(2.0))
    !--> special case N = 1, so power2 = 0
    if(sim%power2 == 0) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 1
       sim%nproc_x4 = 1
    end if
    
    if(is_even(sim%power2)) then
       sim%nproc_x1 = 2**(sim%power2/2)
       sim%nproc_x2 = 2**(sim%power2/2)
       sim%nproc_x3 = 1
       sim%nproc_x4 = 1
    else 
       sim%nproc_x1 = 2**((sim%power2-1)/2)
       sim%nproc_x2 = 2**((sim%power2+1)/2)
       sim%nproc_x3 = 1
       sim%nproc_x4 = 1
    end if

    !--> Initialization of parallel layout of f4d in (x1,x2) directions
    !-->  (x1,x2) : parallelized layout
    !-->  (x3,x4) : sequential
    sim%sequential_x3x4  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%nc_x1+1, &
      sim%nc_x2+1, &
      sim%nc_x3+1, &
      sim%nc_x4+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%sequential_x3x4 )
        
    call compute_local_sizes_4d( sim%sequential_x3x4, &
      sim%loc_sz_x1, &
      sim%loc_sz_x2, &
      sim%loc_sz_x3, &
      sim%loc_sz_x4 )    
    SLL_ALLOCATE(sim%f4d_x3x4(sim%loc_sz_x1,sim%loc_sz_x2,sim%loc_sz_x3,sim%loc_sz_x4),ierr)

    !--> Initialization of parallel layout of f4d in (x3,x4) directions
    !-->  (x1,x2) : sequential
    !-->  (x3,x4) : parallelized layout
    sim%sequential_x1x2  => new_layout_4D( sll_world_collective )
    ! switch x1 and x3:
    itemp        = sim%nproc_x3
    sim%nproc_x3 = sim%nproc_x1
    sim%nproc_x1 = itemp
    ! switch x2 and x4
    itemp        = sim%nproc_x4
    sim%nproc_x4 = sim%nproc_x2 
    sim%nproc_x2 = itemp
    call initialize_layout_with_distributed_4D_array( &
      sim%nc_x1+1, & 
      sim%nc_x2+1, & 
      sim%nc_x3+1, &
      sim%nc_x4+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%sequential_x1x2 )
    
    ! Allocate the array needed to store the local chunk of the distribution
    ! function data. First compute the local sizes. Since the remap operations
    ! are out-of-place, we will allocate two different arrays, one for each
    ! layout.
    call compute_local_sizes_4d( sim%sequential_x1x2, &
      sim%loc_sz_x1, &
      sim%loc_sz_x2, &
      sim%loc_sz_x3, &
      sim%loc_sz_x4 )
    SLL_ALLOCATE(sim%f4d_x1x2(sim%loc_sz_x1,sim%loc_sz_x2,sim%loc_sz_x3,sim%loc_sz_x4),ierr)

    !--> Initialization of f4d_x1x2
    

    !--> Deallocation of temparory array
    SLL_DEALLOCATE(recv_sz,ierr)
    SLL_DEALLOCATE(disps,ierr)
  end subroutine run_4d_DK_hybrid


  !----------------------------------------------------
  ! Initialization of the drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine delete_4d_DK_hybrid( sim )
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32 :: ierr

    SLL_DEALLOCATE(sim%eta1_grid,ierr)
    SLL_DEALLOCATE(sim%vpar_grid,ierr) 
    SLL_DEALLOCATE(sim%n0,ierr)
    SLL_DEALLOCATE(sim%Ti,ierr)
    SLL_DEALLOCATE(sim%Te,ierr)
    SLL_DEALLOCATE(sim%feq_2d,ierr)
    SLL_DEALLOCATE(sim%f4d_x1x2,ierr)
    SLL_DEALLOCATE(sim%f4d_x3x4,ierr)
    call delete(sim%sequential_x1x2)
    call delete(sim%sequential_x3x4)
  end subroutine delete_4d_DK_hybrid


!VG!  !----------------------------------------------------
!VG!  ! Initialization of the generalized QN solver
!VG!  !----------------------------------------------------
!VG!  subroutine initialize_QN_general( &
!VG!    sim, &
!VG!    transformation_x, &
!VG!    a11_f, &
!VG!    a12_f, &
!VG!    a21_f, &
!VG!    a22_f, &
!VG!    c_f)
!VG!
!VG!   type(sll_simulation_4d_DK_hybrid), intent(inout)      :: sim
!VG!   class(sll_coordinate_transformation_2d_base), pointer :: transformation_x
!VG!   procedure(two_var_parametrizable_function) :: a11_f
!VG!   procedure(two_var_parametrizable_function) :: a12_f
!VG!   procedure(two_var_parametrizable_function) :: a21_f
!VG!   procedure(two_var_parametrizable_function) :: a22_f
!VG!   procedure(two_var_parametrizable_function) :: c_f
!VG!
!VG!   sim%mesh2d_x  => mesh2d_x
!VG!   sim%mesh2d_v  => mesh2d_v
!VG!   sim%transfx   => transformation_x
!VG!   sim%init_func => init_func
!VG!   sim%params    => params
!VG!   sim%a11_f     => a11_f
!VG!   sim%a12_f     => a12_f
!VG!   sim%a21_f     => a21_f
!VG!   sim%a22_f     => a22_f
!VG!   sim%c_f       => c_f
!VG!   sim%spline_degree_eta1 = spline_degre1
!VG!   sim%spline_degree_eta2 = spline_degre2
!VG!
!VG!   sim%bc_left   = bc_left
!VG!   sim%bc_right  = bc_right
!VG!   sim%bc_bottom = bc_bottom
!VG!   sim%bc_top    = bc_top
!VG!
!VG!   call sim%interp_phi%initialize( &
!VG!        sim%mesh2d_x%num_cells1 +1, &
!VG!        sim%mesh2d_x%num_cells2 +1, &
!VG!        sim%mesh2d_x%eta1_min, &
!VG!        sim%mesh2d_x%eta1_max, &
!VG!        sim%mesh2d_x%eta2_min, &
!VG!        sim%mesh2d_x%eta2_max, &
!VG!        sim%bc_left, &
!VG!        sim%bc_right, &
!VG!        sim%bc_bottom, &
!VG!        sim%bc_top, &
!VG!        sim%spline_degree_eta1, &
!VG!        sim%spline_degree_eta2)
!VG!
!VG!   call sim%interp_rho%initialize( &
!VG!        sim%mesh2d_x%num_cells1 +1, &
!VG!        sim%mesh2d_x%num_cells2 +1, &
!VG!        sim%mesh2d_x%eta1_min, &
!VG!        sim%mesh2d_x%eta1_max, &
!VG!        sim%mesh2d_x%eta2_min, &
!VG!        sim%mesh2d_x%eta2_max, &
!VG!        sim%bc_left, &
!VG!        sim%bc_right, &
!VG!        sim%bc_bottom, &
!VG!        sim%bc_top, &
!VG!        sim%spline_degree_eta1, &
!VG!        sim%spline_degree_eta2)
!VG!  end subroutine initialize_QN_general

end module sll_simulation_4d_DK_hybrid_module
