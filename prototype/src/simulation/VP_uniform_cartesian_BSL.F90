program VP_1d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"

  use sll_constants
!  use sll_module_mapped_meshes_2d_cartesian
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_cubic_spline_interpolator_1d
  use sll_landau_2d_initializer
  use sll_tsi_2d_initializer
  use distribution_function
  use sll_poisson_1d_periodic
  implicit none

  type(cubic_spline_1d_interpolator), target  :: interp_spline_x, interp_spline_v
  class(sll_interpolator_1d_base), pointer    :: interp_x, interp_v
  !type(sll_mapped_mesh_2d_cartesian), target   :: mesh2d 
  !class(sll_mapped_mesh_2d_base), pointer :: mesh2d_base
  type(sll_logical_mesh_2d), pointer :: mesh2d_cart
  class(sll_coordinate_transformation_2d_base), pointer   :: mesh2d_base
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
  sll_int32  :: freqdiag = 1
  sll_real64 :: time, mass, momentum, kinetic_energy, potential_energy
  sll_real64 :: l1norm, l2norm
  character(len=32) :: fname, case
  sll_int32  :: istep
  sll_int32  :: nbox
  sll_real64 :: eps
  sll_real64 :: v0
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
     print*, '     --------------'
     print*, '     | Landau run |'
     print*, '     --------------'
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
  
 mesh2d_cart => new_logical_mesh_2d( &
       Ncx,  &
       Ncv,  &
       xmin, &  
       xmax, &
       vmin, &
       vmax &
   )
  mesh2d_base => new_coordinate_transformation_2d_analytic( &
       "mesh2d_cart",  &
       mesh2d_cart,    &
       identity_x1,    &
       identity_x2,    &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22, &
       (/ 0.0_f64 /)) 

  ! initialize interpolators
  call interp_spline_x%initialize( Ncx + 1, xmin, xmax, SLL_PERIODIC )
  call interp_spline_v%initialize( Ncv + 1, vmin, vmax, SLL_HERMITE )
  interp_x => interp_spline_x
  interp_v => interp_spline_v

  ! allocate rho and phi
  SLL_ALLOCATE(rho(Ncx+1),ierr)
  SLL_ALLOCATE(efield(Ncx+1),ierr)

  ! initialization of distribution_function
  call init_landau%initialize(mesh2d_base, NODE_CENTERED_FIELD, eps, kmode)
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
       interp_x, &
       interp_v, &
       p_init_f )
  ! write mesh and initial distribution function
  call write_scalar_field_2d(f) 

  ! initialise Poisson
  call initialize(poisson_1d,xmin,xmax,Ncx,ierr)
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
     end do
     ! diagnostics
     time = istep*dt
     mass = delta_x * delta_v * sum(FIELD_DATA(f)(1:Ncx,:))
     l1norm = delta_x * delta_v * sum(abs(FIELD_DATA(f)(1:Ncx,:)))
     l2norm = delta_x * delta_v * sum(FIELD_DATA(f)(1:Ncx,:)*FIELD_DATA(f)(1:Ncx,:))
     momentum = delta_x * delta_v * sum(matmul(FIELD_DATA(f)(1:Ncx,:),v_array))
     kinetic_energy = delta_x * delta_v * 0.5_f64 * &
          sum(matmul(FIELD_DATA(f)(1:Ncx,:),v_array**2))
     potential_energy = delta_x * 0.5_f64 * sum(efield**2)
     write(th_diag,*) time, mass, l1norm, momentum, l2norm, &
          kinetic_energy, potential_energy, kinetic_energy + potential_energy
     write(ex_diag,*) efield
     if (mod(istep,freqdiag)==0) then
        print*, 'iteration: ', istep
        call write_scalar_field_2d(f) 
     end if
  end do

  close(th_diag)
  close(ex_diag)
  print*, 'VP1D_cart has exited normally'
end program 
