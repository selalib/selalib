!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
!
! Vlasov-Ampere simulation in 1Dx1D
!
! contact: Pierre Navaro (navaro@math.unistra.fr)
!
! current investigations:
!   BERK-BREIZMAN test case

module sll_m_sim_bsl_va_1d1v_cart_berk_breizman

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_advection_1d_ampere, only: &
    ampere_1d_advector_ptr, &
    new_ampere_1d_advector

  use sll_m_advection_1d_base, only: &
    sll_advection_1d_base_ptr

  use sll_m_advection_1d_periodic, only: &
    new_periodic_1d_advector

  use sll_m_advection_1d_spectral, only: &
    new_spectral_1d_advector

  use sll_m_ascii_io, only: &
    sll_ascii_file_close, &
    sll_ascii_file_create

  use sll_m_buffer_loader_utilities, only: &
    compute_displacements_array_2d, &
    receive_counts_array_2d

  use sll_m_cartesian_meshes, only: &
    get_node_positions, &
    new_cartesian_mesh_1d, &
    sll_cartesian_mesh_1d, &
    sll_cartesian_mesh_2d, &
    operator(*)

  use sll_m_collective, only: &
    sll_collective_allreduce, &
    sll_collective_gatherv_real64, &
    sll_get_collective_rank, &
    sll_get_collective_size, &
    sll_world_collective

  use sll_m_common_array_initializers, only: &
    sll_beam_initializer_2d, &
    sll_bump_on_tail_initializer_2d, &
    sll_landau_initializer_2d, &
    sll_scalar_initializer_2d, &
    sll_two_stream_instability_initializer_2d

  use sll_m_constants, only: &
    sll_pi

  use sll_m_fft, only: &
    fft_apply_plan_c2r_1d, &
    fft_apply_plan_r2c_1d, &
    fft_forward, &
    fft_new_plan_r2r_1d, &
    sll_fft_plan

  use sll_m_parallel_array_initializer, only: &
    sll_2d_parallel_array_initializer_cartesian

  use sll_m_periodic_interp, only: &
    lagrange, &
    spline, &
    trigo

  use sll_m_poisson_1d_base, only: &
    sll_poisson_1d_base

  use sll_m_poisson_1d_periodic_solver, only: &
    new_poisson_1d_periodic_solver

  use sll_m_remapper, only: &
    apply_remap_2d, &
    compute_local_sizes, &
    initialize_layout_with_distributed_array, &
    layout_2d, &
    local_to_global, &
    new_layout_2d, &
    new_remap_plan, &
    remap_plan_2d_real64, &
    sll_view_lims

  use sll_m_sim_base, only: &
    sll_simulation_base_class

  use sll_m_utilities, only: &
    compute_bloc, &
    compute_mesh_from_bloc, &
    int2string, &
    sll_new_file_id

  use sll_mpi, only: &
    mpi_sum

#ifdef _OPENMP
  use omp_lib, only: &
    omp_get_num_threads, &
    omp_get_thread_num

#endif
  implicit none

  public :: &
    new_va2d_par_cart, &
    sll_simulation_2d_vlasov_ampere_cart

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

integer, parameter :: SLL_ADVECTIVE    = 0
integer, parameter :: SLL_CONSERVATIVE = 1

type, extends(sll_simulation_base_class) :: sll_simulation_2d_vlasov_ampere_cart

 sll_int32                            :: num_threads
 type(sll_cartesian_mesh_2d), pointer :: mesh2d
 sll_int32                            :: num_dof_x2

 sll_real64, dimension(:),   pointer  :: x1_array
 sll_real64, dimension(:),   pointer  :: x2_array
 sll_real64, dimension(:,:), pointer  :: x2_array_omp
 sll_real64, dimension(:),   pointer  :: integration_weight
 sll_int32,  dimension(:),   pointer  :: every_x1
 sll_int32,  dimension(:),   pointer  :: every_x2
 sll_int32                            :: num_bloc_x1
 sll_int32                            :: num_bloc_x2
 sll_int32,  dimension(:),   pointer  :: bloc_index_x1
 sll_int32,  dimension(:),   pointer  :: bloc_index_x2
 
 sll_real64                           :: kx
 sll_real64                           :: eps

 procedure(sll_scalar_initializer_2d), nopass, pointer :: init_func
 sll_real64, dimension(:), pointer                     :: params

 sll_real64          :: nrj0
 sll_real64          :: dt
 sll_int32           :: num_iterations 
 sll_int32           :: freq_diag
 sll_int32           :: freq_diag_time
 sll_int32           :: nb_mode
 sll_real64          :: time_init
 character(len=256)  :: thdiag_filename
 logical             :: driven 
 sll_real64          :: t0
 sll_real64          :: twL
 sll_real64          :: twR
 sll_real64          :: tflat
 sll_real64          :: tL
 sll_real64          :: tR
 sll_real64          :: Edrmax
 sll_real64          :: omegadr
 logical             :: turn_drive_off

 type(sll_advection_1d_base_ptr), dimension(:), pointer :: advect_x1 
 type(sll_advection_1d_base_ptr), dimension(:), pointer :: advect_x2
 type(ampere_1d_advector_ptr),    dimension(:), pointer :: advect_ampere_x1

 sll_real64 :: factor_x1
 sll_real64 :: factor_x2_rho
 sll_real64 :: factor_x2_1

 class(sll_poisson_1d_base), pointer :: poisson 
 logical :: ampere = .false.

 sll_real64 :: L

 sll_real64, dimension(:), allocatable :: node_positions_x2
 sll_real64         :: gamma_d
 sll_real64         :: nu_a 
         
contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 procedure, pass(sim) :: run => run_va2d_cartesian
 procedure, pass(sim) :: init_from_file => init_va2d_fake !init_va2d_par_cart

end type sll_simulation_2d_vlasov_ampere_cart

interface delete
  module procedure delete_va2d_par_cart
end interface delete


sll_int32                               :: istep
sll_int32                               :: iplot
sll_int32                               :: adr_id
sll_int32                               :: edr_id
sll_int32                               :: t_id
sll_int32                               :: deltaf_id
sll_int32                               :: rhotot_id
sll_int32                               :: efield_id
sll_int32                               :: thdiag_id
sll_real64                              :: time

sll_int32,  dimension(:),   allocatable :: collective_displs
sll_int32,  dimension(:),   allocatable :: collective_recvcnts

sll_real64, dimension(:),   allocatable :: buf_fft
sll_real64, dimension(:),   allocatable :: f_x1_buf1d
sll_real64, dimension(:),   allocatable :: f_visu_buf1d
sll_real64, dimension(:,:), allocatable :: f_visu 
sll_real64, dimension(:,:), allocatable :: f1d_omp_in
sll_real64, dimension(:,:), allocatable :: f1d_omp_out

type(sll_fft_plan), pointer :: pfwd
logical                     :: MPI_MASTER

!*********************************************************************************
contains
!*********************************************************************************

function new_va2d_par_cart( filename, num_run ) result(sim)    

type(sll_simulation_2d_vlasov_ampere_cart), pointer :: sim    
character(len=*), intent(in), optional              :: filename
sll_int32, intent(in), optional                     :: num_run
sll_int32                                           :: ierr

SLL_ALLOCATE(sim, ierr)
call init_va2d_par_cart( sim, filename, num_run )
     
end function new_va2d_par_cart

subroutine init_va2d_par_cart( sim, filename, num_run )

class(sll_simulation_2d_vlasov_ampere_cart) :: sim

character(len=*),      optional :: filename
sll_int32, intent(in), optional :: num_run

character(len=*), parameter :: this_sub_name = 'init_va2d_par_cart'
character(len=291)          :: err_msg

character(len=256) :: mesh_case_x1
sll_int32          :: num_cells_x1
sll_real64         :: x1_min
sll_real64         :: x1_max
sll_int32          :: nbox_x1

character(len=256) :: mesh_case_x2
sll_int32          :: num_cells_x2
sll_real64         :: x2_min
sll_real64         :: x2_max
sll_real64         :: x2_fine_min
sll_real64         :: x2_fine_max
sll_real64         :: density_x2_min_to_x2_fine_min
sll_real64         :: density_x2_fine_min_to_x2_fine_max
sll_real64         :: density_x2_fine_max_to_x2_max
sll_int32          :: every_x2_min_to_x2_fine_min
sll_int32          :: every_x2_fine_min_to_x2_fine_max
sll_int32          :: every_x2_fine_max_to_x2_max
sll_int32          :: every_x1
sll_int32          :: every_x2

character(len=256) :: initial_function_case
sll_real64         :: kmode
sll_real64         :: eps
sll_real64         :: alpha_gaussian

sll_real64         :: dt
sll_int32          :: number_iterations
sll_int32          :: freq_diag
sll_int32          :: freq_diag_time
sll_int32          :: nb_mode
sll_real64         :: time_init
character(len=256) :: split_case

!advector
character(len=256) :: advector_x1
sll_int32          :: order_x1
character(len=256) :: advector_x2
sll_int32          :: order_x2
character(len=256) :: integration_case
sll_real64         :: factor_x1
sll_real64         :: factor_x2_rho
sll_real64         :: factor_x2_1

character(len=256) :: poisson_solver
character(len=256) :: ampere_solver
character(len=256) :: drive_type
sll_real64         :: keen_t0
sll_real64         :: keen_tL
sll_real64         :: keen_tR
sll_real64         :: keen_twL
sll_real64         :: keen_twR
sll_real64         :: keen_tflat
logical            :: keen_turn_drive_off
sll_real64         :: keen_Edrmax
sll_real64         :: keen_omegadr

sll_int32                            :: IO_stat
sll_int32                            :: input_file
type(sll_cartesian_mesh_1d), pointer :: mesh_x1
type(sll_cartesian_mesh_1d), pointer :: mesh_x2
sll_int32                            :: ierr
sll_int32, parameter                 :: param_out = 37
sll_int32, parameter                 :: param_out_drive = 40
sll_real64                           :: bloc_coord(2)
sll_int32                            :: bloc_index(3)
sll_int32                            :: i
sll_int32                            :: num_threads
sll_int32                            :: tid
character(len=256)                   :: str_num_run
character(len=256)                   :: filename_loc
sll_real64                           :: gamma_d = 0.0_f64
sll_real64                           :: nu_a = 0.0_f64

namelist /geometry/                   &
  mesh_case_x1,                       &
  num_cells_x1,                       &
  x1_min,                             &
  x1_max,                             &
  nbox_x1,                            &
  mesh_case_x2,                       &
  num_cells_x2,                       &
  x2_min,                             &
  x2_max,                             &
  x2_fine_min,                        &
  x2_fine_max,                        &
  density_x2_min_to_x2_fine_min,      &
  density_x2_fine_min_to_x2_fine_max, &
  density_x2_fine_max_to_x2_max,      &
  every_x2_min_to_x2_fine_min,        &
  every_x2_fine_min_to_x2_fine_max,   &
  every_x2_fine_max_to_x2_max,        &
  every_x1,                           &
  every_x2

namelist /initial_function/           &
  initial_function_case,              &
  kmode,                              &
  eps,                                &
  alpha_gaussian

namelist /time_iterations/            &
  dt,                                 &
  number_iterations,                  &
  freq_diag,                          &
  freq_diag_time,                     &
  nb_mode,                            &
  time_init,                          &
  split_case

namelist /advector/                   &
  advector_x1,                        &
  order_x1,                           &
  advector_x2,                        &
  order_x2,                           &
  factor_x1,                          &
  factor_x2_rho,                      &
  factor_x2_1,                        &
  integration_case

namelist /poisson/                    &
  ampere_solver,                      &
  poisson_solver

namelist /berk_breizman/              &
  gamma_d,                            &
  nu_a

namelist / drive /                    &
  drive_type,                         &
  keen_t0,                            &
  keen_twL,                           &
  keen_twR,                           &
  keen_tflat,                         &
  keen_tL,                            &
  keen_tR,                            &
  keen_turn_drive_off,                &
  keen_Edrmax,                        &
  keen_omegadr

num_threads = 1

if (sll_get_collective_rank(sll_world_collective)==0) then
  MPI_MASTER = .true.
else
  MPI_MASTER = .false.
end if
    
#ifdef _OPENMP
!$OMP PARALLEL SHARED(num_threads)
if (omp_get_thread_num()==0) then
  num_threads =  omp_get_num_threads()      
endif
!$OMP END PARALLEL
#endif

sim%num_threads = num_threads
print *,'#num_threads=',num_threads

! set default parameters

!geometry
mesh_case_x1 = "SLL_LANDAU_MESH"
num_cells_x1 = 32
x1_min       = 0.0_f64
nbox_x1      = 1
mesh_case_x2 = "SLL_CARTESIAN_MESH"
num_cells_x2 = 64
x2_min       = -6._f64
x2_max       =  6._f64

mesh_case_x2 = "SLL_TWO_GRID_MESH"
num_cells_x2 = 64
x2_min       = -6._f64
x2_max       = 6._f64
x2_fine_min  = 0.36_f64
x2_fine_max  = 2.28_f64

density_x2_min_to_x2_fine_min      = 1.0_f64
density_x2_fine_min_to_x2_fine_max = 1.0_f64
density_x2_fine_max_to_x2_max      = 1.0_f64
every_x2_min_to_x2_fine_min        = 1
every_x2_fine_min_to_x2_fine_max   = 1
every_x2_fine_max_to_x2_max        = 1
every_x1                           = 1
every_x2                           = 1

    
!initial_function
initial_function_case = "SLL_LANDAU"   !"SLL_BEAM"
kmode                 = 0.5_f64
eps                   = 0.001_f64
alpha_gaussian        = 0.2_f64

!time_iterations
dt                    = 0.1_f64
number_iterations     = 600
freq_diag             = 100
freq_diag_time        = 1
nb_mode               = 5
time_init             = 0._f64
split_case            = "SLL_STRANG_VTV" 

!advector
advector_x1           = "SLL_LAGRANGE"
order_x1              = 4
advector_x2           = "SLL_LAGRANGE"
order_x2              = 4
factor_x1             = 1._f64
factor_x2_rho         = 1._f64
factor_x2_1           = 1._f64

integration_case      = "SLL_TRAPEZOID" 
poisson_solver        = "SLL_FFT"
drive_type            = "SLL_NO_DRIVE"  

sim%thdiag_filename = "thdiag.dat"
call sll_ascii_file_create(sim%thdiag_filename, thdiag_id, ierr)

if (present(filename)) then

  filename_loc = filename
  filename_loc = adjustl(filename_loc)
  
  if(present(num_run)) then
    filename_loc = trim(filename)//"_"//trim(str_num_run)
  endif

  call sll_new_file_id(input_file, ierr)

  open(unit = input_file, file=trim(filename_loc)//'.nml',IOStat=IO_stat)
  if ( IO_stat /= 0 ) then
    err_msg = 'failed to open file '//trim(filename_loc)//'.nml'
    SLL_ERROR( this_sub_name, err_msg )
  end if

  if (MPI_MASTER) then
    print *,'#initialization with filename:'
    print *,'#',trim(filename_loc)//'.nml'
  endif

  read(input_file, geometry) 
  read(input_file, initial_function)
  read(input_file, time_iterations)
  read(input_file, advector)
  read(input_file, poisson)
  read(input_file, berk_breizman)
  read(input_file, drive)
  close(input_file)

else

  if (sll_get_collective_rank(sll_world_collective)==0) then
    print *,'#initialization with default parameters'
  endif      

endif

select case (mesh_case_x1) !geometry
  case ("SLL_LANDAU_MESH")
    x1_max = real(nbox_x1,f64) * 2._f64 * sll_pi / kmode
    mesh_x1 => new_cartesian_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)
    call get_node_positions( mesh_x1, sim%x1_array )
  case ("SLL_CARTESIAN_MESH")
    mesh_x1 => new_cartesian_mesh_1d(num_cells_x1,eta_min=x1_min, eta_max=x1_max)  
    call get_node_positions( mesh_x1, sim%x1_array )
  case default
    err_msg = '#mesh_case_x1 '//mesh_case_x1//' not implemented'
    SLL_ERROR( this_sub_name, err_msg )
end select

sim%num_bloc_x1 = 1
SLL_ALLOCATE(sim%every_x1(sim%num_bloc_x1),ierr)
SLL_ALLOCATE(sim%bloc_index_x1(sim%num_bloc_x1),ierr)

sim%every_x1(1) = every_x1
sim%bloc_index_x1(sim%num_bloc_x1) = num_cells_x1+1
    
select case (mesh_case_x2)

  case ("SLL_CARTESIAN_MESH")

    mesh_x2 => new_cartesian_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
    call get_node_positions( mesh_x2, sim%x2_array )
    SLL_ALLOCATE(sim%x2_array_omp(num_cells_x2+1,0:sim%num_threads-1),ierr)
    do i=0,sim%num_threads-1
      sim%x2_array_omp(:,i) = sim%x2_array(:)
    enddo
    sim%num_bloc_x2 = 1
    SLL_ALLOCATE(sim%every_x2(sim%num_bloc_x2),ierr)
    SLL_ALLOCATE(sim%bloc_index_x2(sim%num_bloc_x2),ierr)
    sim%every_x2(1) = every_x2
    sim%bloc_index_x2(sim%num_bloc_x2) = num_cells_x2+1

  case ("SLL_TWO_GRID_MESH")

    bloc_coord(1) = (x2_fine_min-x2_min)/(x2_max-x2_min)
    bloc_coord(2) = (x2_fine_max-x2_min)/(x2_max-x2_min)
    bloc_index(1) = floor(density_x2_min_to_x2_fine_min)
    bloc_index(2) = floor(density_x2_fine_min_to_x2_fine_max)
    bloc_index(3) = floor(density_x2_fine_max_to_x2_max)
            
    call compute_bloc(bloc_coord,bloc_index,num_cells_x2)
    SLL_ALLOCATE(sim%x2_array(num_cells_x2+1),ierr)
    call compute_mesh_from_bloc(bloc_coord,bloc_index,sim%x2_array)
    sim%x2_array = x2_min+sim%x2_array*(x2_max-x2_min)
    SLL_ALLOCATE(sim%x2_array_omp(num_cells_x2+1,0:sim%num_threads-1),ierr)
    do i=0,sim%num_threads-1
      sim%x2_array_omp(:,i) = sim%x2_array(:)
    enddo
    mesh_x2 => new_cartesian_mesh_1d(num_cells_x2,eta_min=x2_min, eta_max=x2_max)
    sim%num_bloc_x2 = 3
    SLL_ALLOCATE(sim%every_x2(sim%num_bloc_x2),ierr)
    SLL_ALLOCATE(sim%bloc_index_x2(sim%num_bloc_x2),ierr)

    sim%bloc_index_x2(1:3) = bloc_index(1:3)
    sim%every_x2(1) = every_x2_min_to_x2_fine_min
    sim%every_x2(2) = every_x2_fine_min_to_x2_fine_max
    sim%every_x2(3) = every_x2_fine_max_to_x2_max

  case default
    err_msg = '#mesh_case_x2 '//mesh_case_x2//' not implemented'
    SLL_ERROR( this_sub_name, err_msg )

end select

sim%mesh2d => mesh_x1 * mesh_x2 ! tensor product

!initial function
sim%nrj0 = 0._f64
sim%kx   = kmode
sim%eps  = eps

write(*,*) " kx  = ", sim%kx
write(*,*) " eps = ", sim%eps

select case (initial_function_case)

  case ("SLL_LANDAU")
    sim%init_func => sll_landau_initializer_2d
    SLL_ALLOCATE(sim%params(2),ierr)
    sim%params(1) = kmode
    sim%params(2) = eps
    sim%nrj0      = 0._f64  !compute the right value
    !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
      !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
    !for the moment
    sim%kx = kmode
    sim%eps = eps

  case ("SLL_BUMP_ON_TAIL")
    sim%init_func => sll_bump_on_tail_initializer_2d
    SLL_ALLOCATE(sim%params(2),ierr)
    sim%params(1) = kmode
    sim%params(2) = eps
    sim%nrj0 = 0._f64  !compute the right value
    !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
      !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
    !for the moment
    sim%kx = kmode
    sim%eps = eps

  case ("SLL_TWO_STREAM_INSTABILITY")
    sim%init_func => sll_two_stream_instability_initializer_2d
    SLL_ALLOCATE(sim%params(2),ierr)
    sim%params(1) = kmode
    sim%params(2) = eps
    sim%nrj0 = 0._f64  !compute the right value
    !(0.5_f64*eps*sll_pi)**2/(kmode_x1*kmode_x2) &
      !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
    !for the moment
    sim%kx = kmode
    sim%eps = eps

  case ("SLL_BEAM")  
    sim%init_func => sll_beam_initializer_2d
    SLL_ALLOCATE(sim%params(1),ierr)
    sim%params(1) = alpha_gaussian             

  case default
    err_msg = '#init_func_case not implemented'
    SLL_ERROR( this_sub_name, err_msg )

end select

!time iterations
sim%dt             = dt
sim%num_iterations = number_iterations
sim%freq_diag      = freq_diag
sim%freq_diag_time = freq_diag_time
sim%nb_mode        = nb_mode
sim%time_init      = time_init

if (sim%nb_mode<0) then
  write( err_msg,* ) '#bad value of nb_mode=', nb_mode, '; #should be >=0'
  SLL_ERROR( this_sub_name, err_msg )
endif

SLL_ALLOCATE(sim%advect_x1(num_threads),ierr)
SLL_ALLOCATE(sim%advect_x2(num_threads),ierr)

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(tid)

#ifdef _OPENMP
tid = omp_get_thread_num()+1
#else
tid = 1
#endif

print*,'Advection x1 :', advector_x1
print*,'Advection x2 :', advector_x2

select case (advector_x1)

  case ("SLL_SPLINES") ! arbitrary order periodic splines

    sim%advect_x1(tid)%ptr => new_periodic_1d_advector( &
      num_cells_x1,                                     &
      x1_min,                                           &
      x1_max,                                           &
      SPLINE,                                           &  
      order_x1) 

  case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation

    sim%advect_x1(tid)%ptr => new_periodic_1d_advector( &
      num_cells_x1,                                     &
      x1_min,                                           &
      x1_max,                                           &
      LAGRANGE,                                         & 
      order_x1)

  case("SLL_TRIGO") ! trigo periodic advection

    sim%advect_x1(tid)%ptr => new_periodic_1d_advector( &
      num_cells_x1,                                     &
      x1_min,                                           &
      x1_max,                                           &
      TRIGO,                                            &
      order_x1)

  case("SLL_SPECTRAL") ! spectral periodic advection

    sim%advect_x1(tid)%ptr => new_spectral_1d_advector( &
      num_cells_x1,                                     &
      x1_min,                                           &
      x1_max)

  case default

    err_msg = '#advector in x1 '//advector_x1//' not implemented'
    SLL_ERROR( this_sub_name, err_msg )

end select    

select case (advector_x2)

  case ("SLL_SPLINES") ! arbitrary order periodic splines

    sim%advect_x2(tid)%ptr => new_periodic_1d_advector( &
      num_cells_x2,                                     &
      x2_min,                                           &
      x2_max,                                           &
      SPLINE,                                           & 
      order_x2) 

  case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation

    sim%advect_x2(tid)%ptr => new_periodic_1d_advector( &
      num_cells_x2,                                     &
      x2_min,                                           &
      x2_max,                                           &
      LAGRANGE,                                         & 
      order_x2)

  case default

    err_msg = '#advector in x2 '//advector_x2//' not implemented'
    SLL_ERROR( this_sub_name, err_msg )

end select

!$OMP END PARALLEL

sim%num_dof_x2    = num_cells_x2+1

sim%factor_x1     = factor_x1
sim%factor_x2_rho = factor_x2_rho
sim%factor_x2_1   = factor_x2_1

SLL_ALLOCATE(sim%integration_weight(sim%num_dof_x2),ierr)

select case (integration_case)
  case ("SLL_RECTANGLE")
    do i=1,num_cells_x2
      sim%integration_weight(i) = sim%x2_array(i+1)-sim%x2_array(i)
    enddo
    sim%integration_weight(num_cells_x2+1) = 0._f64   
  case ("SLL_TRAPEZOID")
    sim%integration_weight(1)=0.5_f64*(sim%x2_array(2)-sim%x2_array(1))
    do i=2,num_cells_x2
      sim%integration_weight(i) = 0.5_f64*(sim%x2_array(i+1)-sim%x2_array(i-1))
    enddo  
    sim%integration_weight(num_cells_x2+1) = &
      0.5_f64*(sim%x2_array(num_cells_x2+1)-sim%x2_array(num_cells_x2))
  case ("SLL_CONSERVATIVE")
    do i=1,num_cells_x2
      sim%integration_weight(i)=sim%x2_array(i+1)-sim%x2_array(i)
    enddo  
  case default
    err_msg = '#integration_case not implemented'
    SLL_ERROR( this_sub_name, err_msg )
end select  

select case (poisson_solver)
  case ("SLL_FFT")
    sim%poisson => new_poisson_1d_periodic_solver( x1_min, x1_max, num_cells_x1)
  case default
    err_msg = '#poisson_solver '//poisson_solver//' not implemented'
    SLL_ERROR( this_sub_name, err_msg )
end select

select case (ampere_solver)
  case ("SLL_VLASOV_AMPERE")
    print*,'########################'
    print*,'# Vlasov-Ampere scheme #'
    print*,'########################'
    sim%ampere = .true.
    SLL_ALLOCATE(sim%advect_ampere_x1(num_threads),ierr)
    tid = 1
    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(tid)
    !$ tid = omp_get_thread_num()+1
    sim%advect_ampere_x1(tid)%ptr => new_ampere_1d_advector( &
      num_cells_x1, &
      x1_min,       &
      x1_max )
    !$OMP END PARALLEL
  case default
    continue
end select

sim%driven         = .false.
sim%t0             = 0.0_f64
sim%twL            = 0.0_f64
sim%twR            = 0.0_f64
sim%tflat          = 0.0_f64
sim%tL             = 0.0_f64
sim%tR             = 0.0_f64
sim%turn_drive_off = .false.
sim%Edrmax         = 0.0_f64
sim%omegadr        = 0.0_f64    

sim%gamma_d = gamma_d
sim%nu_a    = nu_a

write(*,*) " gamma_d = ", sim%gamma_d
write(*,*) " nu_a    = ", sim%nu_a

if (sll_get_collective_rank(sll_world_collective)==0) then
        
  open(unit = param_out, file = 'parameters.dat')

  write(param_out, *) real(number_iterations,f64)*dt !Tmax
  write(param_out, *) dt                             !dt
  write(param_out, *) number_iterations              !Nt
  write(param_out, *) kmode                          !k0
  write(param_out, *) sim%omegadr                    !omega0
  write(param_out, *) x1_max-x1_min                  !L
  write(param_out, *) nbox_x1                        !mbox
  write(param_out, *) sim%Edrmax                     !Edrmax
  write(param_out, *) freq_diag_time                 !nsave
  write(param_out, *) freq_diag                      !nsavef1
  write(param_out, *) freq_diag                      !nsavef2
  write(param_out, *) sim%tR+sim%tflat               !Tsetup
  write(param_out, *) x2_max                         !vxmax
  write(param_out, *) x2_min                         !vxmin
  write(param_out, *) num_cells_x1                   !N
  write(param_out, *) num_cells_x2                   !Nv
  ! KEEN drive
  write(param_out, *) sim%tL                         !tL
  write(param_out, *) sim%tR                         !tR
  write(param_out, *) sim%twL                        !twL
  write(param_out, *) sim%twR                        !twR
  write(param_out, *) sim%tflat                      !tflat
      
  close(param_out)
  
endif

sim%L = x1_max - x1_min

end subroutine init_va2d_par_cart

subroutine init_va2d_fake(sim, filename)

class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
character(len=*), intent(in)                               :: filename

character(len=*), parameter :: this_sub_name = 'init_va2d_fake'
character(len=288)          :: err_msg

print *,sim%dt
print *,filename
err_msg = '# Do not use the routine init_va2d_fake\n'// & 
          '# Use instead init_va2d_par_cart'
SLL_ERROR( this_sub_name, err_msg )

end subroutine init_va2d_fake

subroutine run_va2d_cartesian(sim)

class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
character(len=*), parameter :: this_sub_name = 'run_va2d_cartesian'
character(len=*), parameter :: this_prog_name = 'vlasov_ampere_2d'

sll_int32           :: ierr
sll_int32           :: i
sll_int32           :: j

procedure(sll_scalar_initializer_2d), pointer :: init_func

sll_real64, dimension(:,:), pointer :: f_x1
sll_real64, dimension(:,:), pointer :: f_x2

sll_real64, dimension(:),   pointer :: rho
sll_real64, dimension(:),   pointer :: efield
sll_real64, dimension(:),   pointer :: j0
sll_real64, dimension(:),   pointer :: j1
sll_real64, dimension(:),   pointer :: F0

type(layout_2D),            pointer :: layout_x1
type(layout_2D),            pointer :: layout_x2
type(remap_plan_2D_real64), pointer :: remap_plan_x1_x2
type(remap_plan_2D_real64), pointer :: remap_plan_x2_x1
sll_real64, dimension(:),   pointer :: f1d
sll_int32                           :: nc_x1
sll_int32                           :: nc_x2
sll_int32                           :: np_x1
sll_int32                           :: np_x2
sll_int32                           :: nproc_x1
sll_int32                           :: nproc_x2
sll_int32                           :: global_indices(2)
sll_int32                           :: local_size_x1
sll_int32                           :: local_size_x2


sll_int32  :: nb_mode 
sll_int32  :: num_dof_x2 
sll_real64 :: time_init
sll_real64 :: dx, dv
 
! for parallelization (output of distribution function in one single file)
sll_int32                              :: collective_size
character(len=4)                       :: cproc
character(len=4)                       :: cplot
sll_int32                              :: iproc

iplot = 1

nb_mode          = sim%nb_mode
time_init        = sim%time_init
nc_x1            = sim%mesh2d%num_cells1
nc_x2            = sim%mesh2d%num_cells2
np_x1            = sim%mesh2d%num_cells1+1
np_x2            = sim%mesh2d%num_cells2+1
num_dof_x2       = sim%num_dof_x2

collective_size = sll_get_collective_size(sll_world_collective)

if (MPI_MASTER) then

  print *,'#collective_size=', collective_size
  SLL_ALLOCATE(f_visu(np_x1,num_dof_x2),ierr)
  SLL_ALLOCATE(f_visu_buf1d(np_x1*num_dof_x2),ierr)

else

  SLL_ALLOCATE(f_visu(1:1,1:1),ierr)          
  SLL_ALLOCATE(f_visu_buf1d(1:1),ierr)

endif

SLL_ALLOCATE(collective_displs(collective_size),ierr)
SLL_ALLOCATE(collective_recvcnts(collective_size),ierr)

SLL_ALLOCATE(buf_fft(np_x1-1),ierr)
pfwd => fft_new_plan_r2r_1d(np_x1-1,buf_fft,buf_fft,FFT_FORWARD,normalized = .TRUE.)

layout_x1       => new_layout_2D( sll_world_collective )
layout_x2       => new_layout_2D( sll_world_collective )    
nproc_x1 = sll_get_collective_size( sll_world_collective )
nproc_x2 = 1
call initialize_layout_with_distributed_array( &
  np_x1, num_dof_x2, nproc_x1, nproc_x2, layout_x2 )
call initialize_layout_with_distributed_array( &
  np_x1, num_dof_x2, nproc_x2, nproc_x1, layout_x1 )

call sll_view_lims( layout_x1 )
call sll_view_lims( layout_x2 )

call compute_local_sizes( layout_x2, local_size_x1, local_size_x2 )
SLL_ALLOCATE(f_x2(local_size_x1,local_size_x2),ierr)

call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices(1:2) = local_to_global( layout_x1, (/1, 1/) )
SLL_ALLOCATE(f_x1(local_size_x1,local_size_x2),ierr)    
SLL_ALLOCATE(f_x1_buf1d(local_size_x1*local_size_x2),ierr)    

remap_plan_x1_x2 => NEW_REMAP_PLAN(layout_x1, layout_x2, f_x1)
remap_plan_x2_x1 => NEW_REMAP_PLAN(layout_x2, layout_x1, f_x2)

SLL_ALLOCATE(rho(np_x1),ierr)
SLL_ALLOCATE(j0(np_x1),ierr)
SLL_ALLOCATE(j1(np_x1),ierr)
SLL_ALLOCATE(efield(np_x1),ierr)
SLL_ALLOCATE(f1d(max(np_x1,np_x2)),ierr)
SLL_ALLOCATE(f1d_omp_in(max(np_x1,np_x2),sim%num_threads),ierr)
SLL_ALLOCATE(f1d_omp_out(max(np_x1,np_x2),sim%num_threads),ierr)
SLL_ALLOCATE(sim%node_positions_x2(num_dof_x2),ierr)
SLL_ALLOCATE(F0(np_x2),ierr)

sim%node_positions_x2(1:num_dof_x2) = sim%x2_array(1:num_dof_x2)
    
call sll_2d_parallel_array_initializer_cartesian( &
   layout_x1,                                     &
   sim%x1_array,                                  &
   sim%node_positions_x2,                         &
   f_x1,                                          &
   sim%init_func,                                 &
   sim%params)

iproc = sll_get_collective_rank(sll_world_collective)
call int2string(iproc, cproc)
call int2string(iplot, cplot)    

call compute_displacements_array_2d( layout_x1,       &
                                     collective_size, &
                                     collective_displs )

collective_recvcnts = receive_counts_array_2d( layout_x1, &
                                               collective_size )

call sll_collective_gatherv_real64( sll_world_collective,        &
                                    f_x1_buf1d,                  &
                                    local_size_x1*local_size_x2, &
                                    collective_recvcnts,         &
                                    collective_displs,           &
                                    0,                           &
                                    f_visu_buf1d )
    
f_visu = reshape(f_visu_buf1d, shape(f_visu))

call compute_rho(sim, layout_x1, f_x1, rho)
call sim%poisson%compute_E_from_rho( efield, rho )

efield(1) = 0.0_f64
dx = sim%mesh2d%delta_eta1
dv = sim%mesh2d%delta_eta2
do i = 1, nc_x1
  rho(i) = sum(f_x1(i,:))*dv
end do
rho(np_x1) = rho(1)
do i = 2, nc_x1
  efield(i) = efield(i-1)+0.5*(rho(i-1)+rho(i)-2.)*dx
end do
efield(1:nc_x1) = efield - sum(efield(1:nc_x1)) / real(nc_x1,f64)
efield(1) = efield(np_x1)
print*,'ee=', sum(efield(1:nc_x1)*efield(1:nc_x1)) / real(nc_x1,f64)
do i = 1, nc_x1
  write(11,*) sim%x1_array(i), efield(i), rho(i)
end do
    
! Ponderomotive force at initial time. We use a sine wave
! with parameters k_dr and omega_dr.
istep = 0

F0 = (0.9_f64*exp(-0.5_f64*sim%node_positions_x2**2) &
     +0.2_f64*exp(-0.5_f64*(sim%node_positions_x2-4.5_f64)**2/0.5**2))/sqrt(2.0_f64*sll_pi)

do j = 1, np_x2
  write(12,*) sim%node_positions_x2(j), F0(j)
end do

do j = 1, nc_x2
  do i = 1, nc_x1
    write(13,*) sim%x1_array(i), sim%x2_array(j), f_x1(i,j)
  end do
  write(13,*) 
end do

iplot = iplot+1  
call diagnostics(sim, layout_x1, f_x1, efield)

do istep = 1, sim%num_iterations

  !call compute_current(sim, layout_x1, f_x1, j0)
  !call solve_ampere(sim, efield, j0, 0.5_f64*sim%dt)
  !call advection_x( sim, layout_x1, f_x1, 0.5_f64*sim%dt)
  !call advection_poisson_x( sim, layout_x1, f_x1, efield, rho, 0.5_f64*sim%dt)
  call advection_ampere_x(sim, layout_x1, efield, f_x1, 0.5_f64*sim%dt)

  call collision( sim, layout_x1, F0, f_x1, 0.5_f64*sim%dt)


  if (mod(istep,sim%freq_diag_time)==0) then
    call diagnostics(sim, layout_x1, f_x1, efield)
    if (mod(istep,sim%freq_diag)==0) then          
      if (MPI_MASTER) then        
        print *,'#step=',istep,sim%time_init+real(istep,f64)*sim%dt,'iplot=',iplot
      endif
      iplot = iplot+1  
    endif
  end if

  call apply_remap_2D( remap_plan_x1_x2, f_x1, f_x2 )
  call advection_v(sim, layout_x2, f_x2, efield, sim%dt)
  call apply_remap_2D( remap_plan_x2_x1, f_x2, f_x1 )

  call collision( sim, layout_x1, F0, f_x1, 0.5_f64*sim%dt)

  !call compute_current(sim, layout_x1, f_x1, j0)
  !call solve_ampere(sim, efield, j0, 0.5_f64*sim%dt)
  !call advection_x( sim, layout_x1, f_x1, 0.5_f64*sim%dt)
  !call advection_poisson_x( sim, layout_x1, f_x1, efield, rho, 0.5_f64*sim%dt)
  call advection_ampere_x(sim, layout_x1, efield, f_x1, 0.5_f64*sim%dt)


enddo

if (MPI_MASTER) call sll_ascii_file_close(thdiag_id,ierr) 
  
end subroutine run_va2d_cartesian
    
subroutine collision(sim, layout_x1, F0, f_x1, delta_t)

class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
type(layout_2D), intent(in),    pointer        :: layout_x1
sll_real64,      intent(in),    dimension(:)   :: F0
sll_real64,      intent(inout), dimension(:,:) :: f_x1
sll_real64,      intent(in)                    :: delta_t

sll_int32  :: gj
sll_int32  :: j
sll_int32  :: global_indices(2)
sll_int32  :: local_size_x1
sll_int32  :: local_size_x2
sll_int32  :: np_x1, np_x2
sll_int32  :: tid
sll_real64 :: nu_a, coef

nu_a  = sim%nu_a
np_x1 = sim%mesh2d%num_cells1+1
np_x2 = sim%mesh2d%num_cells2+1
call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x1, (/1, 1/) )

coef = exp(-nu_a*delta_t)
tid  = 1

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(j,gj,tid) 
!$ tid = omp_get_thread_num()+1
!$OMP DO
do j = 1, local_size_x2

  gj                 = global_indices(2)-1+j
  f1d_omp_in(:,tid)  = f_x1(:,j)
  f1d_omp_out(:,tid) = coef*f1d_omp_in(:,tid)+(1.0_f64-coef) * F0(gj)
  f_x1(:,j)          = f1d_omp_out(:,tid)

end do
!$OMP END DO          
!$OMP END PARALLEL
    
end subroutine collision
    
subroutine advection_x(sim, layout_x1, f_x1, delta_t)

class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
type(layout_2D), intent(in),    pointer                    :: layout_x1
sll_real64,      intent(inout), dimension(:,:)             :: f_x1
sll_real64                                                 :: delta_t
sll_int32                                                  :: ig_omp
sll_int32                                                  :: i_omp
sll_int32                                                  :: global_indices(2)
sll_int32                                                  :: local_size_x1
sll_int32                                                  :: local_size_x2
sll_int32                                                  :: np_x1
sll_int32                                                  :: tid
sll_real64                                                 :: alpha_omp

np_x1 = sim%mesh2d%num_cells1+1
call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x1, (/1, 1/) )

tid=1          
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid) 
!advection in x
!$ tid = omp_get_thread_num()+1
!$OMP DO
do i_omp = 1, local_size_x2

  ig_omp    = i_omp+global_indices(2)-1
  alpha_omp = sim%factor_x1*sim%node_positions_x2(ig_omp)
  f1d_omp_in(1:np_x1,tid) = f_x1(1:np_x1,i_omp)
  
  call sim%advect_x1(tid)%ptr%advect_1d_constant(  &
    alpha_omp,                                     &
    delta_t,                                       &
    f1d_omp_in(1:np_x1,tid),                       &
    f1d_omp_out(1:np_x1,tid))

  f_x1(1:np_x1,i_omp)=f1d_omp_out(1:np_x1,tid)

end do
!$OMP END DO          
!$OMP END PARALLEL

end subroutine advection_x
    
subroutine advection_poisson_x(sim, layout_x1, f_x1, efield, rho, delta_t)

class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
type(layout_2D), intent(in),    pointer                    :: layout_x1
sll_real64,      intent(inout), dimension(:,:)             :: f_x1
sll_real64,      intent(out),   dimension(:)               :: efield
sll_real64,      intent(out),   dimension(:)               :: rho
sll_real64                                                 :: delta_t
sll_int32                                                  :: ig_omp
sll_int32                                                  :: i_omp
sll_int32                                                  :: global_indices(2)
sll_int32                                                  :: local_size_x1
sll_int32                                                  :: local_size_x2
sll_int32                                                  :: np_x1
sll_int32                                                  :: tid
sll_real64                                                 :: alpha_omp

np_x1 = sim%mesh2d%num_cells1+1
call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x1, (/1, 1/) )

tid=1          
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid) 
!advection in x
!$ tid = omp_get_thread_num()+1
!$OMP DO
do i_omp = 1, local_size_x2

  ig_omp    = i_omp+global_indices(2)-1
  alpha_omp = sim%factor_x1*sim%node_positions_x2(ig_omp)
  f1d_omp_in(1:np_x1,tid) = f_x1(1:np_x1,i_omp)
  
  call sim%advect_x1(tid)%ptr%advect_1d_constant(  &
    alpha_omp,                                     &
    delta_t,                                       &
    f1d_omp_in(1:np_x1,tid),                       &
    f1d_omp_out(1:np_x1,tid))

  f_x1(1:np_x1,i_omp)=f1d_omp_out(1:np_x1,tid)

end do
!$OMP END DO          
!$OMP END PARALLEL

call compute_rho(sim, layout_x1, f_x1, rho)
call sim%poisson%compute_E_from_rho( efield, rho )

!efield = efield*(1.0_f64 - sim%gamma_d*delta_t)
  
end subroutine advection_poisson_x
    
subroutine advection_ampere_x(sim, layout_x1, efield, f_x1, delta_t)

class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
type(layout_2d) , pointer :: layout_x1
sll_real64 :: efield(:)
sll_real64 :: f_x1(:,:)
sll_int32  :: local_size_x1
sll_int32  :: local_size_x2
sll_int32  :: global_indices(2)

sll_real64 :: delta_t
sll_int32  :: nc_x1
sll_int32  :: np_x1
sll_int32  :: tid, ig_omp, i, i_omp
sll_real64 :: alpha_omp
sll_real64 :: L
sll_real64 :: gamma_d

#ifdef _OPENMP
sll_comp64 :: s0, s1
#endif

call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x1, (/1, 1/) )

np_x1 = sim%mesh2d%num_cells1+1
nc_x1 = np_x1-1
tid   = 1          

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid) 
!advection in x
!$ tid = omp_get_thread_num()+1

sim%advect_ampere_x1(tid)%ptr%r0 = cmplx(0.0,0.0,kind=f64)
sim%advect_ampere_x1(tid)%ptr%r1 = cmplx(0.0,0.0,kind=f64)

!$OMP DO 
do i_omp = 1, local_size_x2

  ig_omp    = i_omp+global_indices(2)-1
  alpha_omp = sim%factor_x1*sim%node_positions_x2(ig_omp)*delta_t
  f1d_omp_in(1:np_x1,tid) = f_x1(1:np_x1,i_omp)
  
  sim%advect_ampere_x1(tid)%ptr%d_dx = f1d_omp_in(1:nc_x1,tid)

  call fft_apply_plan_r2c_1d(sim%advect_ampere_x1(tid)%ptr%fwx,  &
       sim%advect_ampere_x1(tid)%ptr%d_dx, &
       sim%advect_ampere_x1(tid)%ptr%fk)

  sim%advect_ampere_x1(tid)%ptr%r0(2:nc_x1/2+1) =    &
       sim%advect_ampere_x1(tid)%ptr%r0(2:nc_x1/2+1) &
     + sim%advect_ampere_x1(tid)%ptr%fk(2:nc_x1/2+1) * sim%integration_weight(ig_omp)

  do i = 2, nc_x1/2+1
    sim%advect_ampere_x1(tid)%ptr%fk(i) =  &
       sim%advect_ampere_x1(tid)%ptr%fk(i) & 
       * cmplx(cos(sim%advect_ampere_x1(tid)%ptr%kx(i)*alpha_omp), &
              -sin(sim%advect_ampere_x1(tid)%ptr%kx(i)*alpha_omp),kind=f64)
  end do

  sim%advect_ampere_x1(tid)%ptr%r1(2:nc_x1/2+1) =    &
       sim%advect_ampere_x1(tid)%ptr%r1(2:nc_x1/2+1) &
     + sim%advect_ampere_x1(tid)%ptr%fk(2:nc_x1/2+1) * sim%integration_weight(ig_omp)

  call fft_apply_plan_c2r_1d(sim%advect_ampere_x1(tid)%ptr%bwx, &
       sim%advect_ampere_x1(tid)%ptr%fk,  &
       sim%advect_ampere_x1(tid)%ptr%d_dx)

  f1d_omp_out(1:nc_x1, tid) = sim%advect_ampere_x1(tid)%ptr%d_dx/nc_x1
  f1d_omp_out(np_x1, tid)   = f1d_omp_out(1, tid) 

  f_x1(1:np_x1,i_omp)=f1d_omp_out(1:np_x1,tid)

end do
!$OMP END DO          

!$OMP END PARALLEL

sim%advect_ampere_x1(tid)%ptr%d_dx = efield(1:nc_x1)
call fft_apply_plan_r2c_1d(sim%advect_ampere_x1(1)%ptr%fwx,  &
     sim%advect_ampere_x1(1)%ptr%d_dx, &
     sim%advect_ampere_x1(1)%ptr%ek)

#ifdef _OPENMP
do i = 2, nc_x1/2+1
  s0 = cmplx(0.0,0.0,kind=f64)
  s1 = cmplx(0.0,0.0,kind=f64)
  do tid = 1, sim%num_threads
    s0 = s0 + sim%advect_ampere_x1(tid)%ptr%r0(i)
    s1 = s1 + sim%advect_ampere_x1(tid)%ptr%r1(i)
  end do
  sim%advect_ampere_x1(1)%ptr%r0(i) = s0
  sim%advect_ampere_x1(1)%ptr%r1(i) = s1
end do
#endif

L =  sim%L / (2.0_f64*sll_pi)
gamma_d = sim%gamma_d * delta_t

sim%advect_ampere_x1(1)%ptr%ek(1) = cmplx(0.,0.,f64)
do i = 2, nc_x1/2+1
  sim%advect_ampere_x1(1)%ptr%ek(i) = + cmplx(L,0.,f64) / cmplx(0.0_f64,real(i-1,f64),f64) * &
     (sim%advect_ampere_x1(1)%ptr%r1(i)-sim%advect_ampere_x1(1)%ptr%r0(i))     &
      +(1.0_f64-gamma_d) *  sim%advect_ampere_x1(1)%ptr%ek(i)
     
end do

call fft_apply_plan_c2r_1d(sim%advect_ampere_x1(1)%ptr%bwx, &
     sim%advect_ampere_x1(1)%ptr%ek,  &
     efield)

efield(np_x1) = efield(1)
efield        = efield/real(nc_x1,f64)

end subroutine advection_ampere_x
  
subroutine compute_rho(sim, layout_x1, f_x1, rho)

class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
type(layout_2d), pointer :: layout_x1
sll_real64      :: rho(:)
sll_real64      :: f_x1(:,:)
sll_int32       :: local_size_x1, local_size_x2
sll_int32       :: global_indices(2)
sll_int32       :: nc_x1
sll_int32       :: np_x1
sll_int32       :: i
sll_int32       :: ig
sll_int32       :: ierr

sll_real64, dimension(:), allocatable :: rho_loc

nc_x1 = sim%mesh2d%num_cells1
np_x1 = sim%mesh2d%num_cells1+1
SLL_ALLOCATE(rho_loc(nc_x1),ierr)

call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x1, (/1, 1/) )

rho_loc = 0._f64
ig = global_indices(2)-1
do i=1,nc_x1
  rho_loc(i)=rho_loc(i)+sum(f_x1(i,1:local_size_x2)                     &
                       *sim%integration_weight(1+ig:local_size_x2+ig))
end do
    
call sll_collective_allreduce( sll_world_collective, &
                               rho_loc,              &
                               nc_x1,                &
                               MPI_SUM,              &
                               rho )

rho(np_x1) = rho(1)
!rho = sim%factor_x2_1-sim%factor_x2_rho*rho
rho = rho - sum(rho)/real(np_x1,f64)


deallocate(rho_loc)
  
end subroutine compute_rho
  
subroutine advection_v(sim, layout_x2, f_x2, efield, delta_t)

class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
type(layout_2d), pointer :: layout_x2
sll_real64 :: f_x2(:,:)
sll_real64 :: efield(:)
sll_real64 :: delta_t, alpha_omp
sll_int32  :: i_omp, ig_omp
sll_int32  :: local_size_x1, local_size_x2
sll_int32  :: global_indices(2)
sll_int32  :: tid, np_x2

np_x2 = sim%mesh2d%num_cells2+1
call compute_local_sizes( layout_x2, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x2, (/1, 1/) )
tid = 1

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid) 
!advection in v
!$ tid = omp_get_thread_num()+1
!advection in v
!$OMP DO
do i_omp = 1,local_size_x1

  ig_omp=i_omp+global_indices(1)-1

  alpha_omp = efield(ig_omp)

  f1d_omp_in(1:sim%num_dof_x2,tid) = f_x2(i_omp,1:sim%num_dof_x2)

  call sim%advect_x2(tid)%ptr%advect_1d_constant(    &
    alpha_omp,                                       &
    delta_t,                                         &
    f1d_omp_in(1:sim%num_dof_x2,tid),                &
    f1d_omp_out(1:sim%num_dof_x2,tid))

  f_x2(i_omp,1:sim%num_dof_x2) = f1d_omp_out(1:sim%num_dof_x2,tid)
end do
!$OMP END DO          
!$OMP END PARALLEL
  
end subroutine advection_v

subroutine compute_current(sim, layout_x1, f_x1, current)
  
class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
type(layout_2d), pointer                                   :: layout_x1
sll_real64,                                  intent(in)    :: f_x1(:,:)
sll_real64,                                  intent(out)   :: current(:)

sll_int32               :: nc_x1
sll_int32               :: np_x1
sll_int32               :: i
sll_int32               :: ierr
sll_int32               :: global_indices(2)
sll_int32               :: local_size_x1
sll_int32               :: local_size_x2
sll_real64              :: v
sll_int32               :: j
sll_int32               :: gj
sll_real64, allocatable :: j_loc(:)

nc_x1 = sim%mesh2d%num_cells1
np_x1 = sim%mesh2d%num_cells1+1
SLL_ALLOCATE(j_loc(nc_x1),ierr)
call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x1, (/1, 1/) )

do i = 1,nc_x1
  j_loc(i) = 0._f64
  do j = 1,local_size_x2
    global_indices = local_to_global( layout_x1, (/i, j/) )
    gj = global_indices(2)
    v  = sim%node_positions_x2(gj)
    j_loc(i) = j_loc(i)+f_x1(i,j)*sim%integration_weight(gj)*v
  end do
end do

call sll_collective_allreduce( sll_world_collective, j_loc, nc_x1, MPI_SUM, current)

current(np_x1) = current(1)
current = current - sum(current) / real(np_x1,f64)
  
end subroutine compute_current

subroutine solve_ampere(sim, e, j, delta_t)
  
class(sll_simulation_2d_vlasov_ampere_cart), intent(inout) :: sim
sll_real64,                                  intent(inout) :: e(:)
sll_real64,                                  intent(in)    :: j(:)
sll_real64,                                  intent(in)    :: delta_t

sll_real64 :: a, g

g = sim%gamma_d
a = exp(-g*delta_t)

if (sim%gamma_d > 0.0_f64) then
  e = a * e + j * (1.0_f64-a)/g
else
  e = a * e + j * delta_t
end if

!a = 1.0_f64-0.5_f64*sim%gamma_d*delta_t
!b = 1.0_f64+0.5_f64*sim%gamma_d*delta_t
!
!efield = (a/b)*efield +delta_t/b*current

end subroutine solve_ampere
  
subroutine diagnostics(sim, layout_x1, f_x1, efield)
  
class(sll_simulation_2d_vlasov_ampere_cart), intent(in) :: sim

type(layout_2d), pointer :: layout_x1
sll_real64, intent(in)   :: f_x1(:,:)
sll_real64, intent(in)   :: efield(:)

sll_int32  :: local_size_x1, local_size_x2
sll_int32  :: global_indices(2)
sll_real64 :: tmp_loc(5), tmp(5)
sll_real64 :: mass, l1norm, l2norm, momentum
sll_real64 :: potential_energy, kinetic_energy
sll_int32  :: i, ig, np_x1
sll_real64, pointer :: vx(:)
sll_real64 :: dx

np_x1            = sim%mesh2d%num_cells1+1
call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )

global_indices   = local_to_global( layout_x1, (/1, 1/) )
time             = sim%time_init+real(istep,f64)*sim%dt
mass             = 0._f64
momentum         = 0._f64
l1norm           = 0._f64
l2norm           = 0._f64
kinetic_energy   = 0._f64
potential_energy = 0._f64
tmp_loc          = 0._f64
ig               = global_indices(2)-1               

dx = sim%mesh2d%delta_eta1
vx => sim%integration_weight(1+ig:local_size_x2+ig)

do i = 1, np_x1-1        
  tmp_loc(1) = tmp_loc(1)+sum(f_x1(i,1:local_size_x2)*vx)
  tmp_loc(2) = tmp_loc(2)+sum(abs(f_x1(i,1:local_size_x2))*vx)
  tmp_loc(3) = tmp_loc(3)+sum((f_x1(i,1:local_size_x2))**2*vx)
  tmp_loc(4) = tmp_loc(4)+sum(f_x1(i,1:local_size_x2)*sim%x2_array(ig+1:ig+local_size_x2)*vx)          
  tmp_loc(5) = tmp_loc(5)+sum(f_x1(i,1:local_size_x2)*sim%x2_array(ig+1:ig+local_size_x2)**2*vx)          
end do

call sll_collective_allreduce( sll_world_collective, tmp_loc, 5, MPI_SUM, tmp )

mass             = tmp(1)  * dx
l1norm           = tmp(2)  * dx
l2norm           = tmp(3)  * dx
momentum         = tmp(4)  * dx
kinetic_energy   = 0.5_f64 * tmp(5) * dx
potential_energy = sum(efield(1:np_x1-1)**2) * dx

if (MPI_MASTER) then                  

  write(thdiag_id,"(f12.5,7f20.12)") &
   time,                             &
   mass,                             &
   l1norm,                           &
   momentum,                         &
   l2norm,                           &
   kinetic_energy,                   &
   potential_energy,                 &
   kinetic_energy+potential_energy

endif
  
end subroutine diagnostics
  
subroutine delete_va2d_par_cart( sim )

class(sll_simulation_2d_vlasov_ampere_cart) :: sim
sll_int32 :: ierr

if(associated(sim%x1_array)) then
  SLL_DEALLOCATE(sim%x1_array,ierr)
  nullify(sim%x1_array)
endif
if(associated(sim%x2_array)) then
  SLL_DEALLOCATE(sim%x2_array,ierr)
  nullify(sim%x2_array)
endif
      
end subroutine delete_va2d_par_cart

end module sll_m_sim_bsl_va_1d1v_cart_berk_breizman
