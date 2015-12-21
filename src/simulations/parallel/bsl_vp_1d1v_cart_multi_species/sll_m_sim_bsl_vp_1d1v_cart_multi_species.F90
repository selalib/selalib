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
! Vlasov-Poisson simulation in 1Dx1D
! Landau and KEEN waves test cases
!
! contact: Michel Mehrenberger (mehrenbe@math.unistra.fr)
!
! current investigations:
!   High order splitting in time
!   KEEN waves with uniform and non uniform grid in velocity

module sll_m_sim_bsl_vp_1d1v_cart_multi_species

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_advection_1d_non_uniform_cubic_splines, only: &
    sll_f_new_non_uniform_cubic_splines_1d_advector

  use sll_m_advection_1d_periodic, only: &
    sll_f_new_periodic_1d_advector

  use sll_m_ascii_io, only: &
    sll_s_ascii_file_close, &
    sll_s_ascii_file_create

  use sll_m_cartesian_meshes, only: &
    sll_o_get_node_positions, &
    sll_f_new_cartesian_mesh_1d, &
    sll_t_cartesian_mesh_1d, &
    sll_f_tensor_product_1d_1d

  use sll_m_collective, only: &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_fft, only: &
    sll_s_fft_apply_plan_r2r_1d, &
    sll_p_fft_forward, &
    sll_f_fft_get_mode_r2c_1d, &
    sll_f_fft_new_plan_r2r_1d, &
    sll_t_fft_plan

  use sll_m_periodic_interp, only: &
    sll_p_lagrange, &
    sll_p_spline

  use sll_m_poisson_1d_base, only: &
    sll_c_poisson_1d_base

  use sll_m_poisson_1d_periodic_solver, only: &
    sll_f_new_poisson_1d_periodic_solver

  use sll_m_poisson_1d_polar_solver, only: &
    sll_f_new_poisson_1d_polar_solver

  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  use sll_m_species, only: &
    sll_s_advection_x1, &
    sll_s_advection_x2, &
    sll_s_compute_rho, &
    sll_s_diagnostics, &
    sll_s_initialize_species, &
    sll_s_read_restart_file, &
    sll_s_set_initial_function, &
    sll_p_advective, &
    sll_p_conservative, &
    sll_t_species, &
    sll_s_write_f, &
    sll_s_write_restart_file

  use sll_m_time_splitting_coeff, only: &
    sll_f_new_time_splitting_coeff, &
    sll_p_lie_tv, &
    sll_p_lie_vt, &
    sll_p_order6_vtv, &
    sll_p_order6vp_tvt, &
    sll_p_order6vp_vtv, &
    sll_p_order6vpnew1_vtv, &
    sll_p_order6vpnew2_vtv, &
    sll_p_order6vpnew_tvt, &
    sll_p_strang_tvt, &
    sll_p_strang_vtv, &
    sll_p_triple_jump_tvt, &
    sll_p_triple_jump_vtv, &
    sll_t_splitting_coeff

  use sll_m_utilities, only: &
    sll_s_int2string, &
    sll_s_pfenvelope, &
    sll_s_new_file_id

#ifdef _OPENMP
  use omp_lib, only: &
    omp_get_num_threads, &
    omp_get_thread_num

#endif
  implicit none

  public :: &
    sll_s_delete_vp2d_par_cart_multi_species, &
    sll_f_new_vp2d_par_cart_multi_species, &
    sll_t_simulation_2d_vlasov_poisson_cart_multi_species

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


type, extends(sll_c_simulation_base_class) :: &
  sll_t_simulation_2d_vlasov_poisson_cart_multi_species
   
  sll_int32 :: istep
  sll_int32 :: iplot
  sll_int32 :: num_threads
  sll_int32 :: num_bloc_x1

  sll_int32,  dimension(:), pointer :: bloc_index_x1
  sll_int32,  dimension(:), pointer :: every_x1
  sll_int32,  dimension(:), pointer :: every_x2
     
  character(len=256) :: restart_file
  logical :: time_init_from_restart_file
  
  type(sll_t_species) :: sp(2)
  
  !time_iterations
  sll_real64 :: dt
  sll_int32  :: num_iterations
  sll_int32  :: freq_diag
  sll_int32  :: freq_diag_time
  sll_int32  :: freq_diag_restart
  sll_int32  :: nb_mode
  sll_real64 :: time_init
  type(sll_t_splitting_coeff), pointer :: split
  character(len=256)      :: thdiag_filename
  
  !parameters for drive
  logical :: driven
  sll_real64 :: t0
  sll_real64 :: twL
  sll_real64 :: twR
  sll_real64 :: tflat
  sll_real64 :: tL
  sll_real64 :: tR
  sll_real64 :: Edrmax
  sll_real64 :: omegadr
  logical :: turn_drive_off


  !poisson solver
  class(sll_c_poisson_1d_base), pointer     :: poisson
  sll_real64                              :: mass_ratio
  sll_real64, dimension(:),   pointer     :: efield
  sll_real64, dimension(:),   pointer     :: e_app
  type(sll_t_fft_plan),         pointer     :: pfwd
  sll_real64, dimension(:),   allocatable :: buf_fft
  sll_comp64, dimension(:),   allocatable :: rho_mode
  sll_int32                               :: rhotote_id
  sll_int32                               :: rhototi_id
  sll_int32                               :: efield_id     
  sll_int32                               :: th_diag_id
   
contains

  procedure, pass(sim) :: run => run_vp2d_cartesian_multi_species
  procedure, pass(sim) :: init_from_file => init_vp2d_fake 

end type sll_t_simulation_2d_vlasov_poisson_cart_multi_species

interface sll_o_delete
  module procedure sll_s_delete_vp2d_par_cart_multi_species
end interface sll_o_delete

!###################################################################
contains

function sll_f_new_vp2d_par_cart_multi_species( filename, num_run) result(sim)    

  type(sll_t_simulation_2d_vlasov_poisson_cart_multi_species), pointer :: sim    
  character(len=*), intent(in), optional               :: filename
  sll_int32, intent(in), optional                      :: num_run
  sll_int32 :: ierr   

  SLL_ALLOCATE(sim,ierr)            
  call init_vp2d_par_cart_multi_species( sim, filename, num_run)
       
end function sll_f_new_vp2d_par_cart_multi_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_vp2d_par_cart_multi_species( sim, filename, num_run )

class(sll_t_simulation_2d_vlasov_poisson_cart_multi_species), &
  intent(inout) :: sim
character(len=*), intent(in), optional :: filename
sll_int32, intent(in), optional :: num_run

character(len=*), parameter :: this_sub_name &
  = 'init_vp2d_par_cart_multi_species'
character(len=128)          :: err_msg

character(len=256) :: mesh_case_x1
sll_int32          :: num_cells_x1
sll_real64         :: x1_min
sll_real64         :: x1_max
sll_int32          :: nbox_x1

character(len=256) :: mesh_case_x2_sp1
sll_int32          :: num_cells_x2_sp1
sll_real64         :: x2_min_sp1
sll_real64         :: x2_max_sp1

character(len=256) :: mesh_case_x2_sp2
sll_int32          :: num_cells_x2_sp2
sll_real64         :: x2_min_sp2
sll_real64         :: x2_max_sp2

!physical_params
sll_real64         :: mass_ratio
sll_real64         :: vd
sll_real64         :: alpha_sp1
sll_real64         :: alpha_sp2

!initial_function
character(len=256) :: initial_function_case_sp1
sll_real64         :: kmode_sp1
sll_real64         :: eps_sp1
sll_real64         :: sigma_sp1
sll_real64         :: v0_sp1
sll_real64         :: factor1_sp1
sll_real64         :: alpha_gaussian_sp1

character(len=256) :: initial_function_case_sp2
sll_real64         :: kmode_sp2
sll_real64         :: eps_sp2
sll_real64         :: sigma_sp2
sll_real64         :: v0_sp2
sll_real64         :: factor1_sp2
sll_real64         :: alpha_gaussian_sp2

!time_iterations
character(len=256) :: restart_file
logical            :: time_init_from_restart_file
sll_real64         :: dt
sll_int32          :: number_iterations
sll_int32          :: freq_diag
sll_int32          :: freq_diag_restart
sll_int32          :: freq_diag_time
sll_int32          :: nb_mode
sll_real64         :: time_init
character(len=256) :: split_case

!advector
character(len=256) :: advector_x1_sp1
sll_int32          :: order_x1_sp1
character(len=256) :: advector_x1_sp2
sll_int32          :: order_x1_sp2
character(len=256) :: advector_x2_sp1
sll_int32          :: order_x2_sp1
character(len=256) :: advector_x2_sp2
sll_int32          :: order_x2_sp2
character(len=256) :: advection_form_x2_sp1
character(len=256) :: advection_form_x2_sp2
character(len=256) :: integration_case
sll_real64         :: factor_x1

!poisson
character(len=256) :: poisson_solver

!drive
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

!local variables
sll_int32                            :: io_stat
sll_int32                            :: input_file
type(sll_t_cartesian_mesh_1d), pointer :: mesh_x1
type(sll_t_cartesian_mesh_1d), pointer :: mesh_x2_sp1
type(sll_t_cartesian_mesh_1d), pointer :: mesh_x2_sp2
sll_int32                            :: ierr
sll_int32, parameter                 :: param_out = 37
sll_int32, parameter                 :: param_out_drive = 40
sll_int32                            :: i
sll_int32                            :: num_threads
sll_int32                            :: tid
character(len=256)                   :: str_num_run
character(len=256)                   :: filename_loc
sll_int32                            :: psize
sll_int32                            :: prank
logical                              :: mpi_master
sll_int32                            :: np_x1, nc_x1
sll_int32                            :: np_x2
  
! namelists for data input
namelist /geometry/          &
  mesh_case_x1,              &
  num_cells_x1,              &
  x1_min,                    &
  x1_max,                    &
  nbox_x1,                   &
  mesh_case_x2_sp1,          &
  mesh_case_x2_sp2,          &
  num_cells_x2_sp1,          &
  num_cells_x2_sp2,          &
  x2_min_sp1,                &
  x2_max_sp1,                &
  x2_min_sp2,                &
  x2_max_sp2

namelist /physical_params/   &
     mass_ratio,             &
     vd,                     &
     alpha_sp1,              &
     alpha_sp2

namelist /initial_function/  &
  initial_function_case_sp1, &
  kmode_sp1,                 &
  eps_sp1,                   &
  sigma_sp1,                 &
  v0_sp1,                    &
  factor1_sp1,               &
  alpha_gaussian_sp1,        &
  initial_function_case_sp2, &
  kmode_sp2,                 &
  eps_sp2,                   &
  sigma_sp2,                 &
  v0_sp2,                    &
  factor1_sp2,               &
  alpha_gaussian_sp2,        &
  restart_file,              &
  time_init_from_restart_file

namelist /time_iterations/   &
  dt,                        &
  number_iterations,         &
  freq_diag,                 &
  freq_diag_time,            &
  freq_diag_restart,         &
  nb_mode,                   &
  time_init,                 &
  split_case

namelist /advector/          &
  advector_x1_sp1,           &
  order_x1_sp1,              &
  advector_x1_sp2,           &
  order_x1_sp2,              &
  advector_x2_sp1,           &
  order_x2_sp1,              &
  advector_x2_sp2,           &
  order_x2_sp2,              &
  advection_form_x2_sp1,     &
  advection_form_x2_sp2,     &
  factor_x1,                 &
  integration_case

namelist /poisson/           &
  poisson_solver

namelist /drive/             &
  drive_type,                &
  keen_t0,                   &
  keen_twL,                  &
  keen_twR,                  &
  keen_tflat,                &
  keen_tL,                   &
  keen_tR,                   &
  keen_turn_drive_off,       &
  keen_Edrmax,               &
  keen_omegadr

prank = sll_f_get_collective_rank( sll_v_world_collective )
psize = sll_f_get_collective_size( sll_v_world_collective )
mpi_master = merge(.true., .false., prank == 0)

!set default parameters
!geometry
mesh_case_x1                = "SLL_LANDAU_MESH"
num_cells_x1                = 32
x1_min                      = 0.0_f64
nbox_x1                     = 1
mesh_case_x2_sp1            = "SLL_CARTESIAN_MESH"
num_cells_x2_sp1            = 64
x2_min_sp1                  = -6._f64
x2_max_sp1                  = 6._f64
mesh_case_x2_sp2            = "SLL_CARTESIAN_MESH"
num_cells_x2_sp2            = 64
x2_min_sp2                  = -0.5_f64
x2_max_sp2                  = 0.5_f64

!physical_params
mass_ratio                  = 0.0005_f64
vd                          = 0._f64
alpha_sp1                   = 1._f64
alpha_sp2                   = 1._f64
    
restart_file                = "no_restart_file"
time_init_from_restart_file = .false.

!time_iterations
dt                          = 0.1_f64
number_iterations           = 600
freq_diag                   = 100
freq_diag_time              = 1
freq_diag_restart           = 5000
nb_mode                     = 5
time_init                   = 0._f64
split_case                  = "sll_p_strang_vtv" 

initial_function_case_sp1   = "SLL_LANDAU"
kmode_sp1                   = 0.5_f64
eps_sp1                     = 0.001_f64
sigma_sp1                   = 1._f64
v0_sp1                      = 0._f64
factor1_sp1                 = 1._f64/sqrt(2._f64*sll_p_pi)
alpha_gaussian_sp1          = 0.2_f64

initial_function_case_sp2   = "SLL_LANDAU"
kmode_sp2                   = 0.5_f64
eps_sp2                     = 0.001_f64
sigma_sp2                   = 1._f64
v0_sp2                      = 0._f64
factor1_sp2                 = 1._f64/sqrt(2._f64*sll_p_pi)

advector_x1_sp1             = "SLL_LAGRANGE"
order_x1_sp1                = 4
advector_x2_sp1             = "SLL_LAGRANGE"
order_x2_sp1                = 4
advection_form_x2_sp1       = "sll_p_advective"

advector_x1_sp2             = "SLL_LAGRANGE"
order_x1_sp2                = 4
advector_x2_sp2             = "SLL_LAGRANGE"
order_x2_sp2                = 4
advection_form_x2_sp2       = "sll_p_advective"

factor_x1                   = 1._f64

integration_case            = "SLL_TRAPEZOID"
poisson_solver              = "SLL_FFT"
drive_type                  = "SLL_NO_DRIVE"

if(present(num_run))then
  write(str_num_run, *) num_run
  str_num_run = adjustl(str_num_run) 
  sim%thdiag_filename = "thdiag_"//trim(str_num_run)//".dat"
else      
  sim%thdiag_filename = "thdiag.dat"
endif

if(present(filename))then

  filename_loc = filename
  filename_loc = adjustl(filename_loc)
  
  if(present(num_run)) then
    filename_loc = trim(filename)//"_"//trim(str_num_run)
  endif

  call sll_s_new_file_id(input_file, ierr)

  open(unit = input_file, file=trim(filename_loc)//'.nml',IOStat=IO_stat)
  if ( IO_stat /= 0 ) then
    err_msg = 'failed to open file '//trim(filename_loc)//'.nml'
    SLL_ERROR( this_sub_name, err_msg )
  end if

  if(mpi_master)then
    print *,'#initialization with filename:'
    print *,'#',trim(filename_loc)//'.nml'
  endif

  read(input_file, geometry)
  read(input_file, physical_params)
  read(input_file, initial_function)
  read(input_file, time_iterations)
  read(input_file, advector)
  read(input_file, poisson)
  read(input_file, drive)
  close(input_file)

else

  if(mpi_master)then
    print *,'#initialization with default parameters'
  endif      

endif

sim%dt                = dt
sim%num_iterations    = number_iterations
sim%freq_diag         = freq_diag
sim%freq_diag_restart = freq_diag_restart
sim%freq_diag_time    = freq_diag_time
sim%nb_mode           = nb_mode
sim%time_init         = time_init

SLL_ASSERT(sim%nb_mode >= 0)

select case (split_case)    
  case ("sll_p_lie_tv")
    sim%split => sll_f_new_time_splitting_coeff(sll_p_lie_tv)
  case ("sll_p_lie_vt") 
    sim%split => sll_f_new_time_splitting_coeff(sll_p_lie_vt)
  case ("sll_p_strang_tvt") 
    sim%split => sll_f_new_time_splitting_coeff(sll_p_strang_tvt)
  case ("sll_p_strang_vtv") 
    sim%split => sll_f_new_time_splitting_coeff(sll_p_strang_vtv)
  case ("sll_p_triple_jump_tvt") 
    sim%split => sll_f_new_time_splitting_coeff(sll_p_triple_jump_tvt)
  case ("sll_p_triple_jump_vtv") 
    sim%split => sll_f_new_time_splitting_coeff(sll_p_triple_jump_vtv)
  case ("sll_p_order6_vtv") 
    sim%split => sll_f_new_time_splitting_coeff(sll_p_order6_vtv)
  case ("sll_p_order6vp_tvt") 
    sim%split => sll_f_new_time_splitting_coeff(sll_p_order6vp_tvt,dt=dt)
  case ("sll_p_order6vp_vtv") 
    sim%split => sll_f_new_time_splitting_coeff(sll_p_order6vp_vtv,dt=dt)
  case ("sll_p_order6vpnew_tvt") 
    sim%split => sll_f_new_time_splitting_coeff(sll_p_order6vpnew_tvt,dt=dt)
  case ("sll_p_order6vpnew1_vtv") 
    sim%split => sll_f_new_time_splitting_coeff(sll_p_order6vpnew1_vtv,dt=dt)
  case ("sll_p_order6vpnew2_vtv") 
    sim%split => sll_f_new_time_splitting_coeff(sll_p_order6vpnew2_vtv,dt=dt)
  case default
    print *,'#split_case not defined'
    print *,'#in initialize_vlasov_par_poisson_seq_cart'
    stop       
end select


select case (mesh_case_x1)
  case ("SLL_LANDAU_MESH")
    x1_max = real(nbox_x1,f64) * 2._f64 * sll_p_pi / kmode_sp1
    mesh_x1 => sll_f_new_cartesian_mesh_1d(num_cells_x1, &
     eta_min=x1_min, eta_max=x1_max)
    call sll_o_get_node_positions( mesh_x1, sim%sp(1)%x1_array )
    call sll_o_get_node_positions( mesh_x1, sim%sp(2)%x1_array )
  case ("SLL_CARTESIAN_MESH")
    mesh_x1 => sll_f_new_cartesian_mesh_1d(num_cells_x1, &
     eta_min=x1_min, eta_max=x1_max)  
    call sll_o_get_node_positions( mesh_x1, sim%sp(1)%x1_array )
    call sll_o_get_node_positions( mesh_x1, sim%sp(2)%x1_array )
  case default
    print*,'#mesh_case_x1', mesh_case_x1, ' not implemented'
    print*,'#in init_vp2d_par_cart'
    stop 
end select

mesh_x2_sp1 => sll_f_new_cartesian_mesh_1d(num_cells_x2_sp1,  &
  eta_min=x2_min_sp1, eta_max=x2_max_sp1)
call sll_o_get_node_positions( mesh_x2_sp1, sim%sp(1)%x2_array )
sim%sp(1)%mesh2d => sll_f_tensor_product_1d_1d( mesh_x1, mesh_x2_sp1)
mesh_x2_sp2 => sll_f_new_cartesian_mesh_1d(num_cells_x2_sp2, &
  eta_min=x2_min_sp2, eta_max=x2_max_sp2)
call sll_o_get_node_positions( mesh_x2_sp2, sim%sp(2)%x2_array )
sim%sp(2)%mesh2d => sll_f_tensor_product_1d_1d( mesh_x1, mesh_x2_sp2)

nc_x1 = sim%sp(1)%mesh2d%num_cells1
np_x1 = sim%sp(1)%mesh2d%num_cells1+1
select case (poisson_solver)
case ("SLL_FFT")
  sim%poisson => sll_f_new_poisson_1d_periodic_solver(x1_min,x1_max,nc_x1)
case ("SLL_POLAR")
  sim%poisson => sll_f_new_poisson_1d_polar_solver(x1_min,x1_max,nc_x1)
case default
  print*,'#poisson_solver', poisson_solver, ' not implemented'
  print *,'#in init_vp2d_par_cart'
  stop 
end select

call sll_s_set_initial_function(sim%sp(1)  &
, initial_function_case_sp1        & 
, kmode_sp1                        &
, eps_sp1                          &
, sigma_sp1                        &
, v0_sp1                           &
, factor1_sp1                      &
, alpha_gaussian_sp1)       

call sll_s_set_initial_function(sim%sp(2)  &
, initial_function_case_sp2        & 
, kmode_sp2                        &
, eps_sp2                          &
, sigma_sp2                        &
, v0_sp2                           &
, factor1_sp2                      &
, alpha_gaussian_sp2)       

sim%time_init_from_restart_file = time_init_from_restart_file
sim%restart_file = restart_file

num_threads = 1
tid = 1
!$OMP PARALLEL
!$ num_threads = omp_get_num_threads()
!$OMP MASTER
print *,'#num_threads=',num_threads
!$OMP END MASTER
!$OMP END PARALLEL

sim%num_threads = num_threads

SLL_ALLOCATE(sim%sp(1)%advect_x1(num_threads),ierr)
SLL_ALLOCATE(sim%sp(1)%advect_x2(num_threads),ierr)
SLL_ALLOCATE(sim%sp(2)%advect_x1(num_threads),ierr)
SLL_ALLOCATE(sim%sp(2)%advect_x2(num_threads),ierr)
np_x1 = sim%sp(1)%mesh2d%num_cells1+1
np_x2 = sim%sp(1)%mesh2d%num_cells2+1
SLL_ALLOCATE(sim%sp(1)%f1d_omp_in (max(np_x1,np_x2),num_threads),ierr)
SLL_ALLOCATE(sim%sp(1)%f1d_omp_out(max(np_x1,np_x2),num_threads),ierr)
np_x1 = sim%sp(2)%mesh2d%num_cells1+1
np_x2 = sim%sp(2)%mesh2d%num_cells2+1
SLL_ALLOCATE(sim%sp(2)%f1d_omp_in (max(np_x1,np_x2),num_threads),ierr)
SLL_ALLOCATE(sim%sp(2)%f1d_omp_out(max(np_x1,np_x2),num_threads),ierr)


!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(tid)
!$ tid = omp_get_thread_num()+1
select case (advector_x1_sp1)
case ("SLL_SPLINES") ! arbitrary order periodic splines
  sim%sp(1)%advect_x1(tid)%ptr => sll_f_new_periodic_1d_advector( &
    num_cells_x1,                                         &
    x1_min,                                               &
    x1_max,                                               &
    sll_p_spline,                                               & 
    order_x1_sp1) 
case("SLL_LAGRANGE") ! arbitrary order sll_p_lagrange periodic interpolation
  sim%sp(1)%advect_x1(tid)%ptr => sll_f_new_periodic_1d_advector( &
    num_cells_x1,                                         &
    x1_min,                                               &
    x1_max,                                               &
    sll_p_lagrange,                                             & 
    order_x1_sp1)
case default
  print*,'#advector in x1 for sll_t_species 1 ', advector_x1_sp1, ' not implemented'
  stop 
end select
   
select case (advector_x1_sp2)
case ("SLL_SPLINES") ! arbitrary order periodic splines
  sim%sp(2)%advect_x1(tid)%ptr => sll_f_new_periodic_1d_advector( &
    num_cells_x1, x1_min, x1_max, sll_p_spline, order_x1_sp2) 
case("SLL_LAGRANGE") ! arbitrary order sll_p_lagrange periodic interpolation
  sim%sp(2)%advect_x1(tid)%ptr => sll_f_new_periodic_1d_advector( &
    num_cells_x1, x1_min, x1_max, sll_p_lagrange, order_x1_sp2)
case default
  print*,'#advector in x1 for sll_t_species 2 ', advector_x1_sp2, ' not implemented'
  stop 
end select    

select case (advector_x2_sp2)
case ("SLL_SPLINES") ! arbitrary order periodic splines
  sim%sp(2)%advect_x2(tid)%ptr => sll_f_new_periodic_1d_advector( &
     num_cells_x2_sp2,                                    &
     x2_min_sp2,                                          &
     x2_max_sp2,                                          &
     sll_p_spline,                                              &  
     order_x2_sp2) 
case("SLL_LAGRANGE") ! arbitrary order sll_p_lagrange periodic interpolation
  sim%sp(2)%advect_x2(tid)%ptr => sll_f_new_periodic_1d_advector( &
     num_cells_x2_sp2,                                    &
     x2_min_sp2,                                          &
     x2_max_sp2,                                          &
     sll_p_lagrange,                                            & 
     order_x2_sp2)
case("SLL_NON_UNIFORM_CUBIC_SPLINES") ! arbitrary order sll_p_lagrange periodic interpolation
  sim%sp(2)%advect_x2(tid)%ptr => sll_f_new_non_uniform_cubic_splines_1d_advector( &
     num_cells_x2_sp2,                                                     &
     x2_min_sp2,                                                           &
     x2_max_sp2,                                                           &
     order_x2_sp2,                                                         &
     sim%sp(2)%x2_array)           
case default
    print*,'#advector in x2 for sll_t_species 2', advector_x2_sp2, ' not implemented'
    stop 
end select
  
select case (advector_x2_sp1)
case ("SLL_SPLINES") ! arbitrary order periodic splines
  sim%sp(1)%advect_x2(tid)%ptr => sll_f_new_periodic_1d_advector( &
    num_cells_x2_sp1,                                     &
    x2_min_sp1,                                           &
    x2_max_sp1,                                           &
    sll_p_spline,                                               & 
    order_x2_sp1) 
case("SLL_LAGRANGE") ! arbitrary order sll_p_lagrange periodic interpolation
  sim%sp(1)%advect_x2(tid)%ptr => sll_f_new_periodic_1d_advector( &
    num_cells_x2_sp1,                                     &
    x2_min_sp1,                                           &
    x2_max_sp1,                                           &
    sll_p_lagrange,                                             & 
    order_x2_sp1)
case("SLL_NON_UNIFORM_CUBIC_SPLINES") ! arbitrary order sll_p_lagrange periodic interpolation
  sim%sp(1)%advect_x2(tid)%ptr => sll_f_new_non_uniform_cubic_splines_1d_advector( &
   num_cells_x2_sp1,                                                       &
   x2_min_sp1,                                                             &
   x2_max_sp1,                                                             &
   order_x2_sp1,                                                           &
   sim%sp(1)%x2_array)           
case default
  print*,'#advector in x2 for sll_t_species 1', advector_x2_sp1, ' not implemented'
  stop 
end select
!$OMP END PARALLEL


select case (advection_form_x2_sp1)
case ("sll_p_advective")
  sim%sp(1)%advection_form_x2 = sll_p_advective
  sim%sp(1)%num_dof_x2 = num_cells_x2_sp1+1
case ("sll_p_conservative")
  sim%sp(1)%advection_form_x2 = sll_p_conservative
  sim%sp(1)%num_dof_x2 = num_cells_x2_sp1
case default
  print*,'#advection_form_x2_sp1', advection_form_x2_sp1, ' not implemented'
  print *,'#in init_vp2d_par_cart'
  stop 
end select  

select case (advection_form_x2_sp2)
case ("sll_p_advective")
  sim%sp(2)%advection_form_x2 = sll_p_advective
  sim%sp(2)%num_dof_x2 = num_cells_x2_sp2+1
case ("sll_p_conservative")
  sim%sp(2)%advection_form_x2 = sll_p_conservative
  sim%sp(2)%num_dof_x2 = num_cells_x2_sp2
case default
  print*,'#advection_form_x2', advection_form_x2_sp2, ' not implemented'
  print *,'#in init_vp2d_par_cart'
  stop 
end select  

sim%sp(1)%factor_x1     = factor_x1
sim%sp(2)%factor_x1     = factor_x1

SLL_ALLOCATE(sim%sp(1)%integration_weight(sim%sp(1)%num_dof_x2),ierr)
SLL_ALLOCATE(sim%sp(2)%integration_weight(sim%sp(2)%num_dof_x2),ierr)

select case (integration_case)

case ("SLL_RECTANGLE")

  do i=1,num_cells_x2_sp1
    sim%sp(1)%integration_weight(i) = sim%sp(1)%x2_array(i+1)-sim%sp(1)%x2_array(i)
  enddo
  sim%sp(1)%integration_weight(num_cells_x2_sp1+1) = 0._f64   

  do i=1,num_cells_x2_sp2
    sim%sp(2)%integration_weight(i) = sim%sp(2)%x2_array(i+1)-sim%sp(2)%x2_array(i)
  enddo
  sim%sp(2)%integration_weight(num_cells_x2_sp2+1) = 0._f64

case ("SLL_TRAPEZOID")

  sim%sp(1)%integration_weight(1)=0.5_f64*(sim%sp(1)%x2_array(2)-sim%sp(1)%x2_array(1))
  do i=2,num_cells_x2_sp1
    sim%sp(1)%integration_weight(i) = 0.5_f64*(sim%sp(1)%x2_array(i+1)-sim%sp(1)%x2_array(i-1))
  enddo  
  sim%sp(1)%integration_weight(num_cells_x2_sp1+1) = &
    0.5_f64*(sim%sp(1)%x2_array(num_cells_x2_sp1+1)-sim%sp(1)%x2_array(num_cells_x2_sp1))

  sim%sp(2)%integration_weight(1)=0.5_f64*(sim%sp(2)%x2_array(2)-sim%sp(2)%x2_array(1))
  do i=2,num_cells_x2_sp2
    sim%sp(2)%integration_weight(i) = 0.5_f64*(sim%sp(2)%x2_array(i+1)-sim%sp(2)%x2_array(i-1))
  enddo  
  sim%sp(2)%integration_weight(num_cells_x2_sp2+1) = &
    0.5_f64*(sim%sp(2)%x2_array(num_cells_x2_sp2+1)-sim%sp(2)%x2_array(num_cells_x2_sp2))

case ("sll_p_conservative")

  do i=1,num_cells_x2_sp1
    sim%sp(1)%integration_weight(i)=sim%sp(1)%x2_array(i+1)-sim%sp(1)%x2_array(i)
  enddo  

  do i=1,num_cells_x2_sp2
    sim%sp(2)%integration_weight(i)=sim%sp(2)%x2_array(i+10)-sim%sp(2)%x2_array(i)
  enddo  

case default

  print *,'#integration_case not implemented'
  print *,'#in init_vp2d_par_cart'  
  stop      

end select  

sim%mass_ratio = mass_ratio
sim%sp(1)%alpha = alpha_sp1    
sim%sp(2)%alpha = alpha_sp2

 
select case (drive_type)
case ("SLL_NO_DRIVE")
  sim%driven = .false.
case("SLL_KEEN_DRIVE")
  sim%driven         = .true.
  sim%t0             = keen_t0
  sim%twL            = keen_twL
  sim%twR            = keen_twR
  sim%tflat          = keen_tflat
  sim%tL             = keen_tL
  sim%tR             = keen_tR
  sim%turn_drive_off = keen_turn_drive_off
  sim%Edrmax         = keen_Edrmax
  sim%omegadr        = keen_omegadr        
case default
  print*,'#drive_type', drive_type, ' not implemented'
  stop 
end select

call sll_s_initialize_species(sim%sp(1), "e", nb_mode)
call sll_s_initialize_species(sim%sp(2), "i", nb_mode)


if(mpi_master .and. sim%driven) then     

  open(unit = param_out, file = 'parameters.dat')
  write(param_out, *) real(number_iterations,f64)*dt !Tmax
  write(param_out, *) dt                             !dt
  write(param_out, *) number_iterations              !Nt
  write(param_out, *) kmode_sp1                      !k0
  write(param_out, *) sim%omegadr                    !omega0
  write(param_out, *) x1_max-x1_min                  !L
  write(param_out, *) nbox_x1                        !mbox
  write(param_out, *) sim%Edrmax                     !Edrmax
  write(param_out, *) freq_diag_time                 !nsave
  write(param_out, *) freq_diag                      !nsavef1
  write(param_out, *) freq_diag                      !nsavef2
  write(param_out, *) sim%tR+sim%tflat               !Tsetup
  write(param_out, *) x2_max_sp1                     !vxmax (sp1)
  write(param_out, *) x2_min_sp1                     !vxmin (sp1)
  write(param_out, *) x2_max_sp2                     !vxmax (sp2)
  write(param_out, *) x2_min_sp2                     !vxmin (sp2)
  write(param_out, *) num_cells_x1                   !N
  write(param_out, *) num_cells_x2_sp1               !Nv
  write(param_out, *) num_cells_x2_sp2               !Nv
  write(param_out, *) sim%tL                         !tL
  write(param_out, *) sim%tR                         !tR
  write(param_out, *) sim%twL                        !twL
  write(param_out, *) sim%twR                        !twR
  write(param_out, *) sim%tflat                      !tflat
  close(param_out)

endif


SLL_ALLOCATE(sim%buf_fft(nc_x1),ierr)
sim%pfwd => sll_f_fft_new_plan_r2r_1d(nc_x1,sim%buf_fft,sim%buf_fft,sll_p_fft_forward,normalized = .TRUE.)
SLL_ALLOCATE(sim%rho_mode(0:nb_mode),ierr)      
SLL_ALLOCATE(sim%efield(np_x1),ierr)
SLL_ALLOCATE(sim%e_app(np_x1),ierr)
sim%e_app = 0._f64

end subroutine init_vp2d_par_cart_multi_species

subroutine init_vp2d_fake(sim, filename)
class(sll_t_simulation_2d_vlasov_poisson_cart_multi_species), intent(inout) :: sim
character(len=*), intent(in)                                :: filename

print *,'# Do not use the routine init_vp2d_fake'
print *,'#use instead init_vp2d_par_cart'
print *,sim%dt
print *,filename
stop

end subroutine init_vp2d_fake

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_vp2d_cartesian_multi_species(sim)

class(sll_t_simulation_2d_vlasov_poisson_cart_multi_species), intent(inout) :: sim

sll_int32  :: ierr
sll_int32  :: istep
sll_int32  :: nb_mode 
sll_int32  :: split_istep

sll_real64 :: time
sll_real64 :: t_step
sll_real64 :: time_init
 
logical :: split_T

character(len=4)   :: cproc
sll_int32          :: psize
sll_int32          :: prank
logical            :: mpi_master

prank = sll_f_get_collective_rank( sll_v_world_collective )
call sll_s_int2string(prank,cproc)
psize = sll_f_get_collective_size( sll_v_world_collective )
mpi_master = merge(.true., .false., prank == 0)

istep          = sim%istep
nb_mode        = sim%nb_mode
time_init      = sim%time_init

if (istep == 1) then

  call sll_s_read_restart_file(sim%restart_file, sim%sp(1), time_init)
  call sll_s_read_restart_file(sim%restart_file, sim%sp(2), time_init)
  
  if(sim%time_init_from_restart_file) sim%time_init = time_init  
  
  sim%iplot = 1
  
  call sll_s_write_f(sim%sp(1), sim%iplot, "fe", time_init)
  call sll_s_write_f(sim%sp(2), sim%iplot, "fi", time_init)
      
  call sll_s_compute_rho(sim%sp(1))
  call sll_s_compute_rho(sim%sp(2))
  
  call sim%poisson%compute_E_from_rho( sim%efield, &
                                       sim%sp(2)%rho-sim%sp(1)%rho )

  if (sim%driven) call compute_e_app(sim,time_init)

  if (mpi_master) then

    call sll_s_ascii_file_create('thdiag.dat' , sim%th_diag_id, ierr)
    call sll_s_ascii_file_create('rhotote.dat', sim%rhotote_id, ierr)
    call sll_s_ascii_file_create('rhototi.dat', sim%rhototi_id, ierr)
    call sll_s_ascii_file_create('efield.dat' , sim%efield_id , ierr)

    write(sim%rhototi_id,*) sim%sp(1)%x1_array
    write(sim%rhotote_id,*) sim%sp(1)%x1_array
    write(sim%efield_id,*)  sim%sp(1)%x1_array

    print *,'#step=',0,time_init,'iplot=',sim%iplot

  endif

endif

if (mod(istep,sim%freq_diag)==0) then
  if (mpi_master) then        
    print *,'#step=',istep,time_init+real(istep,f64)*sim%dt,'iplot=',sim%iplot
  endif
endif  

time_init = sim%time_init
split_T   = sim%split%split_begin_T
t_step    = real(istep-1,f64)

do split_istep=1,sim%split%nb_split_step

  if(split_T) then

    call sll_s_advection_x1(sim%sp(1),                           &
                      sim%split%split_step(split_istep), & 
                      sim%dt)

    call sll_s_advection_x1(sim%sp(2),                           &
                      sim%split%split_step(split_istep), & 
                      sim%dt)

    t_step = t_step+sim%split%split_step(split_istep)

    call sll_s_compute_rho(sim%sp(1))
    call sll_s_compute_rho(sim%sp(2))
    call sim%poisson%compute_E_from_rho( sim%efield, sim%sp(2)%rho-sim%sp(1)%rho )
    
    t_step = t_step+sim%split%split_step(split_istep)
    
    if (sim%driven) call compute_e_app(sim,time_init+t_step*sim%dt)
        
  else

    call sll_s_advection_x2(sim%sp(1),                            &
                      -1.0_f64,                           &
                      sim%efield+sim%sp(1)%alpha*sim%e_app, &
                      sim%split%split_step(split_istep),  &
                      sim%dt)

    call sll_s_advection_x2(sim%sp(2),                            &
                      sim%mass_ratio,                     &
                      sim%efield-sim%sp(2)%alpha*sim%e_app, &
                      sim%split%split_step(split_istep),  &
                      sim%dt)
  endif

  split_T = .not.(split_T)

enddo
    
if (mod(istep,sim%freq_diag_time)==0) then

  time = time_init+real(istep,f64)*sim%dt

  call sll_s_diagnostics(sim%sp(1), sim%pfwd, sim%buf_fft, nb_mode)
  call sll_s_diagnostics(sim%sp(2), sim%pfwd, sim%buf_fft, nb_mode)

  if (mod(istep,sim%freq_diag_restart)==0) then          
    call sll_s_write_restart_file(sim%sp(1), istep, time)
    call sll_s_write_restart_file(sim%sp(2), istep, time)
  endif 

  if (mpi_master) then
    call write_time_history(sim, time)
  end if
    
  if (mod(istep,sim%freq_diag)==0) then          

    sim%iplot = sim%iplot+1  
    call sll_s_write_f( sim%sp(1), sim%iplot, "fe", time)
    call sll_s_write_f( sim%sp(2), sim%iplot, "fi", time)

  endif

end if

if (istep == sim%num_iterations) then
  if(mpi_master)then
    call sll_s_ascii_file_close(sim%th_diag_id,ierr) 
    call sll_s_ascii_file_close(sim%efield_id,ierr)
    call sll_s_ascii_file_close(sim%rhotote_id,ierr)
    call sll_s_ascii_file_close(sim%rhototi_id,ierr)
  endif
  print*, " 176.00010668708197, 820.34117552361215 "
  print"(2f20.14)", sum(sim%sp(1)%f_x1), sum(sim%sp(2)%f_x1)
end if

end subroutine run_vp2d_cartesian_multi_species

subroutine sll_s_delete_vp2d_par_cart_multi_species( sim )
class(sll_t_simulation_2d_vlasov_poisson_cart_multi_species) :: sim
  
print *,'#delete_vp2d_par_cart not implemented'
print *,sim%dt
    
end subroutine sll_s_delete_vp2d_par_cart_multi_species

elemental function f_equilibrium(v)
  sll_real64, intent(in) :: v
  sll_real64 :: f_equilibrium

  f_equilibrium = 1.0_f64/sqrt(2*sll_p_pi)*exp(-0.5_f64*v*v)
end function f_equilibrium


subroutine compute_e_app(sim,t)

class(sll_t_simulation_2d_vlasov_poisson_cart_multi_species), intent(in)  :: sim
sll_real64,                                                 intent(in)  :: t

sll_int32  :: i
sll_int32  :: np_x1
sll_real64 :: adr

sim%e_app = 0._f64
np_x1 = sim%sp(1)%mesh2d%num_cells1+1      
call sll_s_pfenvelope(adr, &
  t,                 &
  sim%tflat,         &
  sim%tL,            &
  sim%tR,            &
  sim%twL,           &
  sim%twR,           &
  sim%t0,            &
  sim%turn_drive_off)

do i = 1, np_x1
 sim%e_app(i) = sim%Edrmax*adr*sim%sp(1)%kx &
   * sin(sim%sp(1)%kx*real(i-1,f64)*sim%sp(1)%mesh2d%delta_eta1 - sim%omegadr*t)
enddo

end subroutine compute_e_app

subroutine write_time_history(sim, time)
class(sll_t_simulation_2d_vlasov_poisson_cart_multi_species) :: sim
sll_int32  :: nb_mode
sll_real64 :: time
sll_real64 :: potential_energy
sll_int32  :: i, k
sll_int32  :: nc_x1

nb_mode = sim%nb_mode
nc_x1 = sim%sp(1)%mesh2d%num_cells1

potential_energy = 0._f64
do i=1, nc_x1
  potential_energy = potential_energy+(sim%efield(i)+sim%e_app(i))**2
enddo
potential_energy = 0.5_f64*potential_energy* sim%sp(1)%mesh2d%delta_eta1

sim%buf_fft = sim%sp(1)%rho(1:nc_x1)-sim%sp(2)%rho(1:nc_x1)
call sll_s_fft_apply_plan_r2r_1d(sim%pfwd,sim%buf_fft,sim%buf_fft)
do k=0,nb_mode
  sim%rho_mode(k)=sll_f_fft_get_mode_r2c_1d(sim%pfwd,sim%buf_fft,k)
enddo  
write(sim%th_diag_id,'(f12.5,12g20.12)',advance='no') &
  time,                                               &
  sim%sp(1)%mass,                                       &
  sim%sp(1)%l1norm,                                     &
  sim%sp(1)%momentum,                                   &
  sim%sp(1)%l2norm,                                     &
  sim%sp(1)%kinetic_energy,                             &
  sim%sp(2)%mass,                                       &
  sim%sp(2)%l1norm,                                     &
  sim%sp(2)%momentum,                                   &
  sim%sp(2)%l2norm,                                     &
  sim%sp(2)%kinetic_energy,                             &
  potential_energy,                                   &
  sim%sp(1)%kinetic_energy+sim%sp(2)%kinetic_energy+potential_energy

do k=0,nb_mode
  write(sim%th_diag_id,'(g20.12)',advance='no') abs(sim%rho_mode(k))
enddo
do k=0,nb_mode-1
  write(sim%th_diag_id,'(2g20.12)',advance='no') &
       sim%sp(1)%f_hat_x2(k+1), &
       sim%sp(2)%f_hat_x2(k+1)
enddo
write(sim%th_diag_id,'(2g20.12)') &
     sim%sp(1)%f_hat_x2(nb_mode+1), &
     sim%sp(2)%f_hat_x2(nb_mode+1)

write(sim%efield_id,*)  sim%efield
write(sim%rhotote_id,*) sim%sp(1)%rho
write(sim%rhototi_id,*) sim%sp(2)%rho

end subroutine write_time_history

end module sll_m_sim_bsl_vp_1d1v_cart_multi_species
