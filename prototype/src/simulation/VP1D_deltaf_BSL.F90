!> Vlasov-Poisson 1D on a uniform cartesian grid
!> using the Backward Semi-Lagrangian (BSL) method.
!> For increased accuracy we work here on delta_f = f-f_M
!> where f_M is the equilibrium Maxwellian function


program VP1d_deltaf
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
  sll_int32, parameter  :: input_file = 33, th_diag = 34, ex_diag = 35
  sll_real64 :: kmode
  sll_real64 :: xmin, xmax, vmin, vmax
  sll_real64 :: delta_x, delta_v
  sll_real64 :: alpha
  sll_real64 :: dt 
  sll_int32  :: nbiter
  sll_real64 :: time, mass, momentum, kinetic_energy, potential_energy
  character(len=32) :: fname, case
  sll_int32  :: istep
  sll_int32  :: nbox
  sll_real64 :: eps
  sll_real64 :: v, v0
  sll_int32  :: i, j
  sll_int32  :: ierr   ! error flag 

  ! namelists for data input
  namelist / geom / xmin, Ncx, nbox, vmin, vmax, Ncv
  namelist / time_iterations / dt, nbiter
  namelist / landau / kmode, eps 
  namelist / tsi / kmode, eps, v0 

  ! determine what case is being run
  call GET_COMMAND_ARGUMENT(1,case)
  ! open and read input file
  if (case == "landau") then
     open(unit = input_file, file = 'landau_input.txt')
     read(input_file, geom) 
     read(input_file, time_iterations)
     read(input_file,landau)
     close(input_file)
  else if (case == "tsi") then
     open(unit = input_file, file = 'tsi_input.txt')
     read(input_file, geom) 
     read(input_file, time_iterations)
     read(input_file,tsi)
     close(input_file)
  else
     print*, 'test case ', case, ' not defined'
     print*, '   usage: VP1D_cart test_case'
     print*, '     where test_case is on of:'
     print*, '        landau'
     print*, '        tsi'
     stop
  endif

  ! define uniform cartesian mesh in x and v
  xmax = nbox * 2 * sll_pi / kmode
  delta_x = (xmax - xmin) / Ncx
  delta_v = (vmax - vmin) / Ncv
  SLL_ALLOCATE(v_array(Ncv+1),ierr)
  do j = 1, Ncv + 1
     v_array(j) = vmin + (j-1)*delta_v
  end do

  ! print out run parameters
  if (case == "landau") then
     print*, '     ----------------------'
     print*, '     | Landau delta_f run |'
     print*, '     ----------------------'
     print*, '   k=', kmode
     print*, '   perturbation=', eps
  else if (case == "tsi") then
     print*, '     -----------'
     print*, '     | TSI run |'
     print*, '     -----------'
     print*, '   k=', kmode
     print*, '   perturbation=', eps
     print*, '   v0=', v0
  end if
  print*, 'geometry of computational domain:'
  print*, '   xmin=', xmin
  print*, '   xmax=', xmax
  print*, '   Ncx=', Ncx
  print*, '   vmin=', vmin
  print*, '   vmax=', vmax
  print*, '   Ncv=', Ncv
  print*, 'time iterations:'
  print*, '   dt=',dt
  print*, '   number of iterations=', nbiter
  print*, ' '
  
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
  call init_landau%initialize(mesh2d_base, NODE_CENTERED_FIELD, eps, kmode, .true.)
  call init_tsi%initialize(mesh2d_base, NODE_CENTERED_FIELD, eps, kmode)
  if (case == "landau") then
     p_init_f => init_landau
     fname = "landau"
  else if (case == "tsi") then
     p_init_f => init_tsi
     fname = "tsi"
  end if
  
  call initialize_distribution_function_2d( &
       f, &
       1.0_f64, &
       1.0_f64, &
       fname, &
       mesh2d_base, &
       NODE_CENTERED_FIELD, &
       p_init_f )
  ! write mesh and initial distribution function
  call write_scalar_field_2d(f) 

  ! initialise Poisson
  call new(poisson_1d,xmin,xmax,Ncx,ierr)
  call solve(poisson_1d, efield, rho)

  ! open files for time history diagnostics
  open(unit = th_diag, file = 'thdiag.dat') 
  open(unit = ex_diag, file = 'exdiag.dat') 
  
  ! time loop
  !----------
  ! half time step advection in v
  do istep = 1, nbiter
     do i = 1, Ncx+1
        alpha = -efield(i) * 0.5_f64 * dt
        f1d => FIELD_DATA(f) (i,:) 
        f1d = interp_v%interpolate_array_disp(Ncv+1, f1d, alpha)
        ! add equilibrium contribution
        do j=1, Ncv + 1
           v = vmin + (j-1) * delta_v
           f1d(j) = f1d(j) + f_equilibrium(v-alpha) - f_equilibrium(v)
        end do
     end do
     ! full time step advection in x
     do j = 1, Ncv+1
        alpha = (vmin + (j-1) * delta_v) * dt
        f1d => FIELD_DATA(f) (:,j) 
        f1d = interp_x%interpolate_array_disp(Ncx+1, f1d, alpha)
     end do
     ! compute rho and electric field
     rho = 1.0_f64 - delta_v * sum(FIELD_DATA(f), DIM = 2)
     call solve(poisson_1d, efield, rho)
     ! half time step advection in v
     do i = 1, Ncx+1
        alpha = -efield(i) * 0.5_f64 * dt
        f1d => FIELD_DATA(f) (i,:) 
        f1d = interp_v%interpolate_array_disp(Ncv+1, f1d, alpha)
        ! add equilibrium contribution
        do j=1, Ncv + 1
           v = vmin + (j-1) * delta_v
           f1d(j) = f1d(j) + f_equilibrium(v-alpha) - f_equilibrium(v)
        end do
     end do
     ! diagnostics
     time = istep*dt
     mass = delta_x * delta_v * sum(FIELD_DATA(f)(1:Ncx,:))
     momentum = delta_x * delta_v * sum(matmul(FIELD_DATA(f)(1:Ncx,:),v_array))
     kinetic_energy = delta_x * delta_v * 0.5_f64 * &
          sum(matmul(FIELD_DATA(f)(1:Ncx,:),v_array**2))
     potential_energy = delta_x * 0.5_f64 * sum(efield**2)
     write(th_diag,*) time, mass, momentum, kinetic_energy, potential_energy
     write(ex_diag,*) efield
     if (mod(istep,10)==0) then
        call write_scalar_field_2d(f) 
     end if
  end do

  close(th_diag)
  close(ex_diag)
  print*, 'VP1D_deltaf_cart has exited normally'
  contains
    function f_equilibrium(v)
      sll_real64 :: v
      sll_real64 :: f_equilibrium

      f_equilibrium = 1.0_f64/sqrt(2*sll_pi)*exp(-0.5_f64*v*v)
    end function f_equilibrium
end program 
