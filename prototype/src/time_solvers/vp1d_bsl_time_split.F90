!> Vlasov-Poisson 1D on a uniform cartesian grid
!> using the Backward Semi-Lagrangian (BSL) method.

program VP1d_BSL_time_split
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"

  use sll_constants
!  use sll_module_mapped_meshes_2d_cartesian
  use sll_common_coordinate_transformations
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_cubic_spline_interpolator_1d
  use sll_landau_2d_initializer
  use sll_tsi_2d_initializer
  use distribution_function
  use sll_poisson_1d_periodic
  use sll_vp_cartesian_2d
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
  type(vp_cartesian_2d), target :: ts_vp_cart_2d
  class(time_splitting), pointer :: time_split
  type(app_field_params) :: params
  sll_real64, dimension(:), allocatable :: rho
  sll_real64, dimension(:), allocatable :: efield
  sll_real64, dimension(:), allocatable :: e_app ! applied field
  sll_real64, dimension(:), pointer :: f1d
  sll_real64, dimension(:), allocatable :: f_maxwellian
  sll_real64, dimension(:), allocatable :: v_array
  sll_int32  :: Ncx, Ncv   ! number of cells
  sll_int32, parameter  :: input_file = 33, th_diag = 34, ex_diag = 35, rho_diag = 36
  sll_int32, parameter  :: param_out = 37, eapp_diag = 38, adr_diag = 39
  sll_real64 :: kmode, omegadr, omegadr0
  logical    :: driven
  sll_int32  :: is_delta_f
  sll_real64 :: xmin, xmax, vmin, vmax
  sll_real64 :: delta_x, delta_v
  sll_real64 :: alpha
  sll_real64 :: dt 
  sll_int32  :: nbiter
  sll_int32  :: freqdiag
  sll_real64 :: time, mass, momentum, kinetic_energy, potential_energy
  sll_real64 :: l1norm, l2norm
  character(len=32) :: fname, case
  sll_int32  :: istep
  sll_int32  :: nbox
  sll_real64 :: eps
  sll_real64 :: v, v0
  sll_int32  :: i, j
  sll_int32  :: ierr   ! error flag 
  sll_real64 :: t0, twL, twR, tstart, tflat, tL, tR
  sll_real64 :: Edrmax, adr
  logical    :: turn_drive_off

  ! namelists for data input
  namelist / geom / xmin, Ncx, nbox, vmin, vmax, Ncv
  namelist / time_iterations / dt, nbiter, freqdiag
  namelist / landau / kmode, eps, is_delta_f, driven 
  namelist / tsi / kmode, eps, v0 
  namelist / drive / t0, twL, twR, tstart, tflat, tL, tR, turn_drive_off, Edrmax, omegadr


  ! determine what case is being run
  call GET_COMMAND_ARGUMENT(1,case)
  ! open and read input file
  if (case == "landau") then
     open(unit = input_file, file = 'landau_input.txt')
     read(input_file, geom) 
     read(input_file, time_iterations)
     read(input_file, landau)
     if (driven) then
        read(input_file, drive)
        eps = 0.0  ! no initial perturbation for driven simulation
     end if
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

  ! save parameters of applied field in params derived type
  params =  app_field_params( Edrmax, tflat, tL, tR, twL, twR, t0, kmode, omegadr, &
       turn_drive_off, driven )
  ! define uniform cartesian mesh in x and v
  xmax = nbox * 2 * sll_pi / kmode
  delta_x = (xmax - xmin) / Ncx
  delta_v = (vmax - vmin) / Ncv
  SLL_ALLOCATE(v_array(Ncv+1),ierr)
  do j = 1, Ncv + 1
     v_array(j) = vmin + (j-1)*delta_v
  end do
  ! allocate f_maxwellian for diagnostics
  SLL_ALLOCATE(f_maxwellian(Ncv+1),ierr)
  if (is_delta_f == 0) then
     f_maxwellian = f_equilibrium(v_array)
  else 
     f_maxwellian = 0.0_f64
  end if
  ! print out run parameters
  if (case == "landau") then
     print*, '     ----------------------'
     print*, '     | Landau delta_f run |'
     print*, '     ----------------------'
     print*, '   k=', kmode
     print*, '   perturbation=', eps
     print*, '   driven=', driven
     if (driven) then
        print*, 'omegadr=', omegadr
     endif
  else if (case == "tsi") then
     print*, '     -----------'
     print*, '     | TSI run |'
     print*, '     -----------'
     print*, '   k=', kmode
     print*, '   perturbation=', eps
     print*, '   v0=', v0
  end if
  if (is_delta_f == 0) then
     print*, '   delta_f version'
  else
     print*, '   full_f version'
  end if
  print*, 'geometry of computational domain:'
  print*, '   xmin=', xmin
  print*, '   xmax=', xmax
  print*, '   Ncx=', Ncx
  print*, '   vmin=', vmin
  print*, '   vmax=', vmax
  print*, '   Ncv=', Ncv
  print*, 'time iterations:'
  print*, '   dt=', dt
  print*, '   number of iterations=', nbiter
  print*, ' '
  open(unit = param_out, file = 'param_out.dat') 
  write(param_out,*) trim(case), xmin, xmax, ncx, vmin, vmax, ncv, dt, nbiter, freqdiag, &
       is_delta_f
  close(param_out)

!!$  call initialize_mesh_2d_cartesian( &
!!$       mesh2d,           &
!!$       "mesh2d_cart",       &
!!$       xmin,         &
!!$       xmax,         &
!!$       Ncx+1,          &
!!$       vmin,         &
!!$       vmax,         &
!!$       Ncv+1           &
!!$       )
!!$  mesh2d_base => mesh2d

  mesh2d_cart => new_logical_mesh_2d( &
       Ncx, &
       Ncv,  &
       xmin,       &
       xmax,       &
       vmin,       &
       vmax       &
   )
  mesh2d_base => new_coordinate_transformation_2d_analytic( &
       "mesh2d_cart",      &
       mesh2d_cart,             &
       identity_x1,    &
       identity_x2,    &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22, &
       (/0.0_f64/) ) 

  ! initialize interpolators
  call interp_spline_x%initialize( Ncx + 1, xmin, xmax, SLL_PERIODIC )
  call interp_spline_v%initialize( Ncv + 1, vmin, vmax, SLL_HERMITE )
  interp_x => interp_spline_x
  interp_v => interp_spline_v

  ! allocate rho and phi
  SLL_ALLOCATE(rho(Ncx+1),ierr)
  SLL_ALLOCATE(efield(Ncx+1),ierr)
  SLL_ALLOCATE(e_app(Ncx+1),ierr)

  ! initialization of distribution_function
  call init_landau%initialize(mesh2d_base, NODE_CENTERED_FIELD, eps, kmode, is_delta_f)
  call init_tsi%initialize(mesh2d_base, NODE_CENTERED_FIELD, eps, kmode, v0, is_delta_f)
  if (case == "landau") then
     p_init_f => init_landau
  else if (case == "tsi") then
     p_init_f => init_tsi
  end if

  fname = 'dist_func'
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

  ! initialize Poisson
  call initialize(poisson_1d,xmin,xmax,Ncx,ierr)
  call solve(poisson_1d, efield, rho)
  ! Ponderomotive force at initial time. We use a sine wave
  ! with parameters k_dr and omega_dr.
  istep = 0
  if (driven) then
     call PFenvelope(adr, istep*dt, params)
     do i = 1, Ncx + 1
        e_app(i) = Edrmax * adr * kmode * sin(kmode * (i-1) * delta_x)
     enddo
  endif

  ! open files for time history diagnostics
  open(unit = th_diag, file = 'thdiag.dat') 
  open(unit = ex_diag, file = 'exdiag.dat')
  open(unit = rho_diag, file = 'rhodiag.dat') 
  open(unit = eapp_diag, file = 'eappdiag.dat')
  open(unit = adr_diag, file = 'adrdiag.dat') 

  ! write initial fields
  write(ex_diag,*) efield
  write(rho_diag,*) rho
  write(eapp_diag,*) e_app
  write(adr_diag,*) istep*dt, adr

  ! time loop
  !----------
  ! half time step advection in v
  do istep = 1, nbiter
     call time_split%lie_splitting(dt, 1)
     ! diagnostics
     time = istep*dt
     mass = 0.
     momentum = 0.
     l1norm = 0.
     l2norm = 0.
     kinetic_energy = 0.
     potential_energy = 0.
     do i = 1, Ncx 
        mass = mass + sum(FIELD_DATA(f)(i,:) + f_maxwellian)   
        l1norm = l1norm + sum(abs(FIELD_DATA(f)(i,:) + f_maxwellian))
        l2norm = l2norm + sum((FIELD_DATA(f)(i,:) + f_maxwellian)**2)
        momentum = momentum + sum(FIELD_DATA(f)(i,:)*v_array)
        kinetic_energy = kinetic_energy + 0.5_f64 * &
             sum((FIELD_DATA(f)(i,:) + f_maxwellian)*(v_array**2))
     end do
     mass = mass * delta_x * delta_v 
     l1norm = l1norm  * delta_x * delta_v
     l2norm = l2norm  * delta_x * delta_v
     momentum = momentum * delta_x * delta_v
     kinetic_energy = kinetic_energy * delta_x * delta_v
     potential_energy =   0.5_f64 * sum(efield**2) * delta_x
     write(th_diag,*) time, mass, l1norm, momentum, l2norm, &
          kinetic_energy, potential_energy, kinetic_energy + potential_energy
     write(ex_diag,*) efield
     write(rho_diag,*) rho
     write(eapp_diag,*) e_app
     write(adr_diag,*) istep*dt, adr
     if (mod(istep,freqdiag)==0) then
        print*, 'iteration: ', istep
        call write_scalar_field_2d(f) 
     end if
  end do

  close(th_diag)
  close(ex_diag)
  print*, 'VP1D_deltaf_cart has exited normally'

end program VP1d_BSL_time_split
