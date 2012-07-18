program VP_1d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"

  use numeric_constants
  use sll_module_mapped_meshes_2d_cartesian
  use sll_cubic_spline_interpolator_1d
  use sll_landau_2d_initializer
  use sll_tsi_2d_initializer
  use distribution_function
  use sll_poisson_1d_periodic
  implicit none

  type(cubic_spline_1d_interpolator), target  :: interp_spline_x, interp_spline_v
  class(sll_interpolator_1d_base), pointer    :: interp_x, interp_v
  type(sll_mapped_mesh_2d_cartesian), target   :: mesh2d 
  class(sll_mapped_mesh_2d_base), pointer :: mesh2d_base
  type(init_landau_2d), target :: init_landau
  type(init_tsi_2d), target :: init_tsi
  class(scalar_field_2d_initializer_base), pointer    :: p_init_f
  type(sll_distribution_function_2d)   :: f
  type(poisson_1d_periodic)  :: poisson_1d
  sll_real64, dimension(:), allocatable :: rho
  sll_real64, dimension(:), allocatable :: efield
  sll_real64, dimension(:), pointer :: f1d
  sll_real64, dimension(:), allocatable :: v_array
  sll_int32  :: Ncx, Ncv   ! number of cells
  sll_int32, parameter  :: th_diag = 34, ex_diag = 35
  sll_real64 :: kmode
  sll_real64 :: xmin, xmax, vmin, vmax, delta_v
  sll_real64 :: alpha
  sll_real64 :: dt 
  sll_int32  :: nbiter
  sll_real64 :: mass, momentum, kinetic_energy, potential_energy
  character(len=32) :: fname, name
  sll_int32  :: istep
  character(len=4) :: cstep
  sll_real64 :: eps
  sll_int32 :: i, j
  sll_int32 :: ierr   ! error flag 

  ! define uniform cartesian mesh in x and v
  kmode = 0.2_f64    ! mode for Landau Damping
  Ncx = 128 
  Ncv = 128
  xmin = 0.0_f64
  xmax = 2.0_f64 * sll_pi / kmode
  vmin = -8.0_f64
  vmax = 8.0_f64
  delta_v = (vmax - vmin) / Ncv
  SLL_ALLOCATE(v_array(Ncv+1),ierr)
  do j = 1, Ncv + 1
     v_array(j) = vmin + (j-1)*delta_v
  end do
  
  call initialize_mesh_2d_cartesian( &
    mesh2d,           &
    "mesh2d_cart",       &
    xmin,         &
    xmax,         &
    Ncx+1,          &
    vmin,         &
    vmax,         &
    Ncv+1           &
   )
  mesh2d_base => mesh2d

  ! initialize interpolators
  call interp_spline_x%initialize( Ncx + 1, xmin, xmax, PERIODIC_SPLINE )
  call interp_spline_v%initialize( Ncv + 1, vmin, vmax, HERMITE_SPLINE )
  interp_x => interp_spline_x
  interp_v => interp_spline_v

  ! allocate rho and phi
  SLL_ALLOCATE(rho(Ncx+1),ierr)
  SLL_ALLOCATE(efield(Ncx+1),ierr)

  ! initialization of distribution_function
  eps = 0.05_f64  ! perturbation for Landau
  call init_landau%initialize(mesh2d_base, NODE_CENTERED_FIELD, eps)
  call init_tsi%initialize(mesh2d_base, NODE_CENTERED_FIELD, eps)
  !p_init_f => init_landau
  p_init_f => init_tsi

  fname = "landau"
  call initialize_distribution_function_2d( &
       f, &
       1.0_f64, &
       1.0_f64, &
       fname, &
       mesh2d_base, &
       NODE_CENTERED_FIELD, &
       p_init_f )
  ! write mesh and initial distribution function
  istep = 0
  call int2string(istep,cstep)
  name = trim(fname)//cstep
  call write_scalar_field_2d(f,output_file_name=name) 

  ! initialise Poisson
  call new(poisson_1d,xmin,xmax,Ncx,ierr)
  call solve(poisson_1d, efield, rho)

  ! open files for time history diagnostics
  open(unit = th_diag, file = 'thdiag.dat') 
  open(unit = ex_diag, file = 'exdiag.dat') 
  
  ! time loop
  !----------
  nbiter = 500
  dt = .1_f64
  ! half time step advection in v
  do istep = 1, nbiter
     do i = 1, Ncx+1
        alpha = -efield(i) * 0.5_f64 * dt
        f1d => FIELD_DATA(f) (i,:) 
        f1d = interp_v%interpolate_array_disp(Ncv+1, f1d, alpha)
     end do
     ! full time step advection in x
     do j = 1, Ncv+1
        alpha = (vmin + (j-1) * delta_v) * dt
        f1d => FIELD_DATA(f) (:,j) 
        f1d = interp_x%interpolate_array_disp(Ncx+1, f1d, alpha)
     end do
     ! compute rho and electric field
     rho = 1.0_f64 - sum(FIELD_DATA(f), DIM = 2)
     call solve(poisson_1d, efield, rho)
     ! half time step advection in v
     do i = 1, Ncx+1
        alpha = -efield(i) * 0.5_f64 * dt
        f1d => FIELD_DATA(f) (i,:) 
        f1d = interp_v%interpolate_array_disp(Ncv+1, f1d, alpha)
     end do
     ! diagnostics
     mass = sum(FIELD_DATA(f))
     momentum = sum(matmul(FIELD_DATA(f),v_array))
     kinetic_energy = 0.5_f64 * sum(matmul(FIELD_DATA(f),v_array**2))
     potential_energy = 0.5_f64 * sum(efield**2)
     write(th_diag,*) mass, momentum, kinetic_energy, potential_energy
     write(ex_diag,*) efield
     if (mod(istep,10)==0) then
        call int2string(istep,cstep)
        name = trim(fname)//cstep
        call write_scalar_field_2d(f,output_file_name=name) 
     end if
  end do

  close(th_diag)
  close(ex_diag)
  print*, 'VP1D_cart has exited normally'
end program 
