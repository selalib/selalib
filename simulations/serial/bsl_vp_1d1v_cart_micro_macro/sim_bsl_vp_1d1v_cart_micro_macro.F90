!> Vlasov-Poisson 1D on a uniform cartesian grid
!> using the Backward Semi-Lagrangian (BSL) method.
!> For increased accuracy we can work here on delta_f = f-f_M
!> where f_M is the equilibrium Maxwellian function
!> driven simulations (i.e. with an external force or non driven can
!> be performed

program sim_bsl_vp_1d1v_cart_micro_macro
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_hermite, &
      sll_p_periodic

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_cubic_spline_interpolator_1d, only: &
      sll_t_cubic_spline_interpolator_1d

   use sll_m_cubic_splines, only: &
      sll_s_cubic_spline_1d_compute_interpolant, &
      sll_s_cubic_spline_1d_eval_array, &
      sll_s_cubic_spline_1d_init, &
      sll_t_cubic_spline_1d

   use sll_m_hdf5_io_serial, only: &
      sll_t_hdf5_ser_handle, &
      sll_s_hdf5_ser_file_create, &
      sll_s_hdf5_ser_file_close, &
      sll_o_hdf5_ser_write_array

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   use sll_m_periodic_interp, only: &
      sll_p_lagrange, &
      sll_p_spline

   use sll_m_periodic_interpolator_1d, only: &
      sll_t_periodic_interpolator_1d

   use sll_m_poisson_1d_periodic, only: &
      sll_o_initialize, &
      sll_t_poisson_1d_periodic, &
      sll_o_solve

   use sll_m_utilities, only: &
      sll_s_int2string, &
      sll_s_pfenvelope

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!  type(sll_t_cubic_spline_interpolator_1d), target  ::  interp_spline_x
   type(sll_t_cubic_spline_1d)  :: interp_spline_v, interp_spline_vh, interp_spline_x
   type(sll_t_periodic_interpolator_1d), target      :: interp_per_x, interp_per_v
   type(sll_t_cubic_spline_interpolator_1d), target      :: interp_comp_v
   class(sll_c_interpolator_1d), pointer    :: interp_x, interp_v
   type(sll_t_poisson_1d_periodic)  :: poisson_1d
   sll_real64, dimension(:, :), allocatable, target :: f
   sll_real64, dimension(:, :), allocatable :: fg, ff, ff1, ff2
   sll_real64, dimension(:), allocatable :: rho
   sll_real64, dimension(:), allocatable :: efield
   sll_real64, dimension(:), allocatable :: e_app ! applied field
   sll_real64, dimension(:), allocatable :: fx
   !sll_real64, dimension(:), pointer     :: f1d
   sll_real64, dimension(:), allocatable :: f_maxwellian
   sll_real64, dimension(:), allocatable :: v_array, vg_array, vh_array, vhg_array
   sll_real64, dimension(:), allocatable :: x_array, xg_array
   sll_int32  :: Ncx, Ncv, Ncvh   ! number of cells
   sll_int32  :: iHmin, iHmax, iraf   ! indices of the coarse mesh to define fine mesh; iraf=dv/dvh
   sll_int32, parameter  :: input_file = 33, th_diag = 34, ex_diag = 35, rho_diag = 36
   sll_int32, parameter  :: param_out = 37, eapp_diag = 38, adr_diag = 39
   sll_int32, parameter  :: param_out_drive = 40
   sll_real64 :: kmode, omegadr!, omegadr0
   sll_int32  :: is_delta_f
   logical    :: driven
   sll_real64 :: xmin, xmax, vmin, vmax
   sll_real64 :: delta_x, delta_v, dvh
   sll_int32  :: interpol_x, interpol_v  ! type of interpolator
   sll_int32  :: order_x, order_v   ! order of interpolator
   sll_real64 :: alpha
   sll_real64 :: dt
   sll_int32  :: nbiter
   sll_int32  :: freqdiag = 1
   sll_real64 :: time, mass, momentum, kinetic_energy, potential_energy
   sll_real64 :: l1norm, l2norm
   !sll_real64 :: equilibrium_contrib
   character(len=32) :: case
   sll_int32  :: istep
   sll_int32  :: nbox
   sll_real64 :: eps
   sll_real64 :: v, v0
   sll_int32  :: i, j
   sll_int32  :: ierr   ! error flag
   sll_real64 :: t0, twL, twR, tstart, tflat, tL, tR
   sll_real64 :: adr, Edrmax
   logical    :: turn_drive_off
   !sll_int32  :: istartx, jstartv, jendv
   !sll_int32  :: num_threads, my_num
   !sll_int32  :: ipiece_size_x, ipiece_size_v
   !type(sll_time_mark) :: time0
   !sll_real64 :: time1
   sll_int32  :: error!, file_id
   type(sll_t_hdf5_ser_handle) :: hfile_id
   character(len=4) :: cstep
   !character(len=32) :: dsetname

   ! namelists for data input
   namelist /geom/ xmin, Ncx, nbox, vmin, vmax, Ncv, iHmin, iHmax, iraf
   namelist /interpolator/ interpol_x, order_x, interpol_v, order_v
   namelist /time_iterations/ dt, nbiter, freqdiag
   namelist /landau/ kmode, eps, is_delta_f, driven
   namelist /tsi/ kmode, eps, v0, is_delta_f
   namelist /drive/ t0, twL, twR, tstart, tflat, tL, tR, turn_drive_off, Edrmax, omegadr

   ! determine what case is being run
   call get_command_argument(1, case)
   ! open and read input file
   if (case == "landau") then
      open (unit=input_file, file='landau_input.nml')
      read (input_file, geom)
      read (input_file, interpolator)
      read (input_file, time_iterations)
      read (input_file, landau)
      if (driven) then
         read (input_file, drive)
         eps = 0.0_f64  ! no initial perturbation for driven simulation
      end if
      close (input_file)
   else if (case == "tsi") then
      open (unit=input_file, file='tsi_input.nml')
      read (input_file, geom)
      read (input_file, interpolator)
      read (input_file, time_iterations)
      read (input_file, tsi)
      close (input_file)
   elseif (case == "keen") then
      open (unit=input_file, file='ex-keen.nml')
      read (input_file, geom)
      read (input_file, interpolator)
      read (input_file, time_iterations)
      read (input_file, landau)
      if (driven) then
         read (input_file, drive)
         eps = 0.0_f64  ! no initial perturbation for driven simulation
      end if
      close (input_file)
   else
      print *, 'test case ', case, ' not defined'
      print *, '   usage: VP1D_cart test_case'
      print *, '     where test_case is on of:'
      print *, '        landau'
      print *, '        tsi'
      print *, '        keen'
      stop
   end if

   ! define uniform cartesian mesh in x and coarse mesh in v
   xmax = nbox*2*sll_p_pi/kmode
   delta_x = (xmax - xmin)/Ncx
   delta_v = (vmax - vmin)/Ncv

   ! define fine mesh in v
   dvh = delta_v/real(iraf, 8)
   Ncvh = (iHmax - iHmin + 1)*iraf

   SLL_ALLOCATE(x_array(Ncx + 1), ierr)
   SLL_ALLOCATE(xg_array(Ncx + 1), ierr)
   SLL_ALLOCATE(v_array(Ncv + 1), ierr)
   SLL_ALLOCATE(vg_array(Ncv + 1), ierr)
   SLL_ALLOCATE(vh_array(Ncvh + 1), ierr)
   SLL_ALLOCATE(vhg_array(Ncvh + 1), ierr)

   do i = 1, Ncx + 1
      x_array(i) = xmin + (i - 1)*delta_x
   end do
   do j = 1, Ncv + 1
      v_array(j) = vmin + (j - 1)*delta_v
   end do
   do j = 1, Ncvh + 1
      vh_array(j) = v_array(iHmin) + (j - 1)*dvh
   end do

   ! allocate f_maxwellian for diagnostics
   SLL_ALLOCATE(f_maxwellian(Ncv + 1), ierr)
   if (is_delta_f == 0) then
      f_maxwellian = f_equilibrium(v_array)
   else
      f_maxwellian = 0.0_f64
   end if
   ! print out run parameters
   if (case == "landau") then
      print *, '     ----------------------'
      print *, '     | Landau delta_f run |'
      print *, '     ----------------------'
      print *, '   k=', kmode
      print *, '   perturbation=', eps
      print *, '   driven=', driven
      if (driven) then
         print *, 'omegadr=', omegadr
      end if
   else if (case == "tsi") then
      print *, '     -----------'
      print *, '     | TSI run |'
      print *, '     -----------'
      print *, '   k=', kmode
      print *, '   perturbation=', eps
      print *, '   v0=', v0
   elseif (case == "keen") then
      print *, '     ----------------------'
      print *, '     | keen delta_f run |'
      print *, '     ----------------------'
      print *, '   k=', kmode
      print *, '   perturbation=', eps
      print *, '   driven=', driven
      if (driven) then
         print *, 'omegadr=', omegadr
      end if
   end if
   if (is_delta_f == 0) then
      print *, '   delta_f version'
   else
      print *, '   full_f version'
   end if
   print *, 'geometry of computational domain:'
   print *, '   xmin=', xmin
   print *, '   xmax=', xmax
   print *, '   Ncx=', Ncx
   print *, '   vmin=', vmin
   print *, '   vmax=', vmax
   print *, '   Ncv=', Ncv
   print *, 'time iterations:'
   print *, '   dt=', dt
   print *, '   number of iterations=', nbiter
   print *, ' '
   print *, 'interpolator in x', interpol_x, order_x
   print *, 'interpolator in v', interpol_v, order_v
   open (unit=param_out, file='param_out.dat')
   write (param_out, '(A6,2f10.3,I5,2f10.3,I5,f10.3,I8,I5,I2,f10.3)') &
      trim(case), xmin, xmax, ncx, vmin, vmax, ncv, &
      dt, nbiter, freqdiag, is_delta_f, kmode
   close (param_out)

   if (driven) then
      tstart = t0  ! is this parameter actually necessary ?
      open (unit=param_out_drive, file='param_out_drive.dat')
      write (param_out_drive, '(9g15.5)') t0, twL, twR, tstart, tflat, tL, tR, &
         Edrmax, omegadr
      close (param_out_drive)
   end if

   ! allocate rho and phi
   SLL_ALLOCATE(fx(Ncx + 1), ierr)
   SLL_ALLOCATE(rho(Ncx + 1), ierr)
   SLL_ALLOCATE(efield(Ncx + 1), ierr)
   SLL_ALLOCATE(e_app(Ncx + 1), ierr)

   ! open files for time history diagnostics
   open (unit=th_diag, file='thdiag.dat')
   open (unit=ex_diag, file='exdiag.dat')
   open (unit=rho_diag, file='rhodiag.dat')
   open (unit=eapp_diag, file='eappdiag.dat')
   open (unit=adr_diag, file='adrdiag.dat')

   ! initialization of sll_m_distribution_function
   SLL_ALLOCATE(f(Ncx + 1, Ncv + Ncvh - (iHmax - iHmin + 1) + 1), ierr)
   SLL_ALLOCATE(fg(Ncx + 1, Ncv + 1), ierr)
   SLL_ALLOCATE(ff(Ncx + 1, Ncvh + 1), ierr)
   SLL_ALLOCATE(ff1(Ncx + 1, Ncvh + 1), ierr)
   SLL_ALLOCATE(ff2(Ncx + 1, Ncvh + 1), ierr)

#ifdef __OPENMP
!$omp parallel default(shared) &
!$omp& private(i,alpha,v,j,f1d,my_num,istartx,iendx, jstartv, jendv,  &
!$omp& interp_x, interp_v, interp_spline_x, interp_spline_v, &
!$omp& interp_per_x, interp_per_v, time1, time0)

   my_num = omp_get_thread_num()
   num_threads = omp_get_num_threads()
   print *, 'running with openmp using ', num_threads, ' threads'
   ipiece_size_v = (Ncv + 1)/num_threads
   ipiece_size_x = (Ncx + 1)/num_threads

   istartx = my_num*ipiece_size_x + 1
   if (my_num < num_threads - 1) then
      iendx = istartx - 1 + ipiece_size_x
   else
      iendx = Ncx + 1
   end if
   jstartv = my_num*ipiece_size_v + 1
   if (my_num < num_threads - 1) then
      jendv = jstartv - 1 + ipiece_size_v
   else
      jendv = Ncv + 1
   end if
#endif

   ! sll_o_initialize interpolators
   select case (interpol_x)
   case (1) ! periodic cubic sll_p_spline
      call sll_s_cubic_spline_1d_init(interp_spline_x, Ncx + 1, xmin, xmax, sll_p_periodic)
!     call interp_spline_x%init( Ncx + 1, xmin, xmax, sll_p_periodic )
!     interp_x => interp_spline_x
   case (2) ! arbitrary order periodic splines
      call interp_per_x%init(Ncx + 1, xmin, xmax, sll_p_spline, order_x)
      interp_x => interp_per_x
   case (3) ! arbitrary order sll_p_lagrange periodic interpolation
      call interp_per_x%init(Ncx + 1, xmin, xmax, sll_p_lagrange, order_x)
      interp_x => interp_per_x
   case default
      print *, 'interpolation in x number ', interpol_x, ' not implemented'
   end select
   select case (interpol_v)
   case (1) ! hermite cubic sll_p_spline
      call sll_s_cubic_spline_1d_init(interp_spline_v, Ncv + 1, vmin, vmax, sll_p_hermite)
      call sll_s_cubic_spline_1d_init(interp_spline_vh, Ncvh + 1, vh_array(1), vh_array(Ncvh + 1), sll_p_hermite)
   case (2) ! arbitrary order periodic splines
      call interp_per_v%init(Ncv + 1, vmin, vmax, sll_p_spline, order_v)
      interp_v => interp_per_v
   case (3) ! arbitrary order sll_p_lagrange periodic interpolation
      call interp_per_v%init(Ncv + 1, vmin, vmax, sll_p_lagrange, order_v)
      interp_v => interp_per_v
   case (4) ! arbitrary order open sll_p_spline interpolation
      !PN remove odd_degree replaced by cubic_spline
      call interp_comp_v%init(Ncv + 1, vmin, vmax, sll_p_hermite)!order_v)
   case default
      print *, 'interpolation in x number ', interpol_v, ' not implemented'
   end select

#ifdef __OPENMP
!$omp barrier
!$omp single
   fname = 'dist_func'
#endif

   ! write mesh and initial distribution function
!  call write_scalar_field_2d(f) --> remplacer par une ecriture hdf5

!fg = d.f. on coarse mesh in v (v_array)
!ff = d.f. on fine mesh in v (vh_array)
!such that sum_{v_j in vh_array} [ ff(x_i, v_j) - fg(x_i, v_j) ] = 0 forall i

   ff1 = 0._f64; ff2 = 0._f64
   do j = 1, Ncv + 1
      v = v_array(j)
      do i = 1, Ncx + 1
         fg(i, j) = (1._f64/sqrt(2._f64*sll_p_pi))*exp(-0.5_f64*v*v)
      end do
   end do

   !compute deltaf such that int deltaf(v)dv = 0
   !--------------------------------------------
   do i = 1, Ncx + 1
      !compute splines coef associated to fg and evalute splines on the fine mesh vh_array --> ff1
      call sll_s_cubic_spline_1d_compute_interpolant(fg(i, :), interp_spline_v)
      call sll_s_cubic_spline_1d_eval_array(vh_array, ff1(i, :), Ncvh + 1, interp_spline_v)

      !compute ff:=deltaf on the fine mesh: ff(v_j)=f(v_j)-ff1(v_j), v_j\in vh_array
      mass = 0._f64
      do j = 1, Ncvh + 1
         v = vh_array(j)
         ff(i, j) = (1._f64/sqrt(2._f64*sll_p_pi))*exp(-0.5_f64*v*v) - ff1(i, j)
         mass = mass + ff(i, j)
      end do
      ff(i, :) = ff(i, :) - mass/real(Ncvh + 1, f64)

   end do
   !compute f on the fine mesh
   ff = ff1 + ff

   !plot of f on the fine mesh
   open (12, file="ffinit")
   do i = 1, Ncx + 1
      do j = 1, Ncvh + 1
         write (12, *) i, j, ff(i, j)
      end do
      write (12, *)
   end do
   close (12)

   ! sll_o_initialize Poisson
   call sll_o_initialize(poisson_1d, xmin, xmax, Ncx, ierr)
   if (is_delta_f == 0) then
      rho = -delta_v*sum(f, DIM=2)
   else
      !compute density on coarse grid
      rho = 1.0_f64 - delta_v*sum(fg, DIM=2)
   end if
   call sll_o_solve(poisson_1d, efield, rho)

   ! Ponderomotive force at initial time. We use a sine wave
   ! with parameters k_dr and omega_dr.
   istep = 0
   if (driven) then
      call sll_s_pfenvelope(adr, istep*dt, tflat, tL, tR, twL, twR, &
                            t0, turn_drive_off)
      do i = 1, Ncx + 1
         e_app(i) = Edrmax*adr*kmode*sin(kmode*(i - 1)*delta_x)
      end do
   end if

   ! write initial fields
   do i = 1, Ncx + 1
      write (ex_diag, "(g15.5)", advance="no") efield(i)
      write (rho_diag, "(g15.5)", advance="no") rho(i)
      write (eapp_diag, "(g15.5)", advance="no") e_app(i)
   end do
   write (ex_diag, *)
   write (rho_diag, *)
   write (eapp_diag, *)
   write (adr_diag, '(2g15.5)') istep*dt, adr
#ifdef __OPENMP
!$omp end single
#endif

   ! sll_o_initialize timer
   ! time loop
   !----------

   do istep = 1, nbiter
!     time0 => reset_time_mark(time0)

      ! half time step advection in v
      do i = 1, Ncx + 1 !istartx, iendx

         !compute splines coef associated to fg
         !and compute fg^{n+1}(v_j)=fg^n(v_j^*) (v_j on the coarse mesh) -> fg
         call sll_s_cubic_spline_1d_compute_interpolant(fg(i, :), interp_spline_v)
         alpha = -(efield(i) + e_app(i))*0.5_f64*dt
         do j = 1, Ncv + 1
            vg_array(j) = v_array(j) + alpha
            if (vg_array(j) < vmin) then
               vg_array(j) = vmin
            end if
            if (vg_array(j) > vmax) then
               vg_array(j) = vmax
            end if
         end do
         call sll_s_cubic_spline_1d_eval_array(vg_array, fg(i, :), Ncv + 1, interp_spline_v)

         !compute fg^{n+1}(v_j)=fg^n(v_j^*) (v_j on the fine mesh) -> ff2
         call sll_s_cubic_spline_1d_eval_array(vh_array + alpha, ff2(i, :), Ncvh + 1, interp_spline_v)

         !compute fg on the fine mesh -> ff1
         call sll_s_cubic_spline_1d_eval_array(vh_array, ff1(i, :), Ncvh + 1, interp_spline_v)

         !compute deltaf=ff1-ff on the fine mesh + zero average -> ff
         mass = 0._f64
         do j = 1, Ncvh + 1
            ff(i, j) = ff1(i, j) - ff(i, j)
            mass = mass + ff(i, j)
         end do
         ff(i, :) = ff(i, :) - mass/real(Ncvh + 1, f64)

         !compute splines coef associated to ff
         call sll_s_cubic_spline_1d_compute_interpolant(ff(i, :), interp_spline_vh)

         do j = 1, Ncvh + 1
            vhg_array(j) = vh_array(j) + alpha
            if (vhg_array(j) < vh_array(1)) then
               vhg_array(j) = vh_array(1)
            end if
            if (vhg_array(j) > vh_array(Ncvh + 1)) then
               vhg_array(j) = vh_array(Ncvh + 1)
            end if
         end do

         call sll_s_cubic_spline_1d_eval_array(vhg_array, ff(i, :), Ncvh + 1, interp_spline_vh)
         !update deltaf on the fine mesh: delta^{n+1}=ff2+ff-ff1
         !f^{n+1} = f^n(v*)= (ff2 + ff)(v*)
         ff(i, :) = ff2(i, :) + ff(i, :)

      end do

#ifdef __OPENMP
!$omp barrier
#endif

      ! advection in x
      do j = 1, Ncv + 1

         !compute splines coef associated to fg
         !and compute, for the coarse grid in v fg^{n+1}(x_i)=fg^n(x_i^*) -> fg
         call sll_s_cubic_spline_1d_compute_interpolant(fg(:, j), interp_spline_x)
         alpha = (vmin + (j - 1)*delta_v)*dt
         do i = 1, Ncx + 1
            xg_array(i) = x_array(i) - alpha
            if (xg_array(i) < xmin) then
               xg_array(i) = xg_array(i) + xmax
            end if
            if (xg_array(i) > xmax) then
               xg_array(i) = xg_array(i) - xmax
            end if
         end do
         call sll_s_cubic_spline_1d_eval_array(xg_array, fg(:, j), Ncx + 1, interp_spline_x)
      end do

      do j = 1, Ncvh + 1
         !compute splines coef associated to ff
         !and compute, for the fine grid in v ff^{n+1}(x_i)=ff^n(x_i^*) -> ff
         call sll_s_cubic_spline_1d_compute_interpolant(ff(:, j), interp_spline_x)
         alpha = vh_array(j)*dt
         do i = 1, Ncx + 1
            xg_array(i) = x_array(i) - alpha
            if (xg_array(i) < xmin) then
               xg_array(i) = xg_array(i) + xmax
            end if
            if (xg_array(i) > xmax) then
               xg_array(i) = xg_array(i) - xmax
            end if
         end do
         call sll_s_cubic_spline_1d_eval_array(xg_array, ff(:, j), Ncx + 1, interp_spline_x)

      end do

#ifdef __OPENMP
!$omp barrier

!$omp single
#endif

      ! compute rho and electric field
      if (is_delta_f == 0) then
         rho = -delta_v*sum(f, DIM=2)
      else
         !compute density on coarse grid
         rho = 1.0_f64 - delta_v*sum(fg, DIM=2)
      end if
      call sll_o_solve(poisson_1d, efield, rho)

      if (driven) then
         call sll_s_pfenvelope(adr, istep*dt, tflat, tL, tR, twL, twR, &
                               t0, turn_drive_off)
         do i = 1, Ncx + 1
            E_app(i) = Edrmax*adr*kmode*sin(kmode*(i - 1)*delta_x &
                                            - omegadr*istep*dt)
         end do
      end if

#ifdef __OPENMP
!$omp end single
#endif

      do i = 1, Ncx + 1 !istartx, iendx !

         !compute splines coef associated to fg
         !and compute fg^{n+1}(v_j)=fg^n(v_j^*) (v_j on the coarse mesh) -> fg
         call sll_s_cubic_spline_1d_compute_interpolant(fg(i, :), interp_spline_v)
         alpha = -(efield(i) + e_app(i))*0.5_f64*dt
         do j = 1, Ncv + 1
            vg_array(j) = v_array(j) + alpha
            if (vg_array(j) < vmin) then
               vg_array(j) = vmin
            end if
            if (vg_array(j) > vmax) then
               vg_array(j) = vmax
            end if
         end do
         call sll_s_cubic_spline_1d_eval_array(vg_array, fg(i, :), Ncv + 1, interp_spline_v)

         !compute fg^{n+1}(v_j)=fg^n(v_j^*) (v_j on the fine mesh) -> ff2
         call sll_s_cubic_spline_1d_eval_array(vh_array + alpha, ff2(i, :), Ncvh + 1, interp_spline_v)

         !compute fg on the fine mesh -> ff1
         call sll_s_cubic_spline_1d_eval_array(vh_array, ff1(i, :), Ncvh + 1, interp_spline_v)

         !compute deltaf=ff1-ff on the fine mesh + zero average -> ff
         mass = 0._f64
         do j = 1, Ncvh + 1
            ff(i, j) = ff1(i, j) - ff(i, j)
            mass = mass + ff(i, j)
         end do
         ff(i, :) = ff(i, :) - mass/real(Ncvh + 1, f64)

         !compute splines coef associated to ff
         call sll_s_cubic_spline_1d_compute_interpolant(ff(i, :), interp_spline_vh)

         do j = 1, Ncvh + 1
            vhg_array(j) = vh_array(j) + alpha
            if (vhg_array(j) < vh_array(1)) then
               vhg_array(j) = vh_array(1)
            end if
            if (vhg_array(j) > vh_array(Ncvh + 1)) then
               vhg_array(j) = vh_array(Ncvh + 1)
            end if
         end do

         call sll_s_cubic_spline_1d_eval_array(vhg_array, ff(i, :), Ncvh + 1, interp_spline_vh)

         !update deltaf on the fine mesh: delta^{n+1}=ff2+ff-ff1
         !f^{n+1} = f^n(v*)= (ff2 + ff)(v*)
         ff(i, :) = ff2(i, :) + ff(i, :)

      end do

#ifdef __OPENMP
!$omp barrier
!$omp single
#endif

      ! diagnostics
      if (mod(istep, freqdiag) == 0) then
         time = istep*dt
         print *, 'time', time
         mass = 0.0_f64
         momentum = 0.0_f64
         l1norm = 0.0_f64
         l2norm = 0.0_f64
         kinetic_energy = 0.0_f64
         potential_energy = 0.0_f64
         do i = 1, Ncx
            mass = mass + sum(fg(i, :))
            l1norm = l1norm + sum(abs(fg(i, :)))
            l2norm = l2norm + sum((fg(i, :))**2)
            momentum = momentum + sum(fg(i, :)*v_array)
            kinetic_energy = kinetic_energy + 0.5_f64* &
                             sum((fg(i, :))*(v_array**2))
         end do
         mass = mass*delta_x*delta_v
         l1norm = l1norm*delta_x*delta_v
         l2norm = l2norm*delta_x*delta_v
         momentum = momentum*delta_x*delta_v
         kinetic_energy = kinetic_energy*delta_x*delta_v
         potential_energy = 0.5_f64*sum(efield**2)*delta_x
         write (th_diag, '(f12.5,7g20.12)') time, mass, l1norm, momentum, l2norm, &
            kinetic_energy, potential_energy, kinetic_energy + potential_energy
         do i = 1, Ncx + 1
            write (ex_diag, "(g15.5)", advance="no") efield(i)
            write (rho_diag, "(g15.5)", advance="no") rho(i)
            write (eapp_diag, "(g15.5)", advance="no") e_app(i)
         end do
         write (ex_diag, *)
         write (rho_diag, *)
         write (eapp_diag, *)

         write (adr_diag, '(2g15.5)') istep*dt, adr
         print *, 'iteration: ', istep
         !call write_scalar_field_2d(f)
         call sll_s_int2string(istep, cstep)
         call sll_s_hdf5_ser_file_create("f"//cstep//".h5", hfile_id, error)
         call sll_o_hdf5_ser_write_array(hfile_id, f, "f", error)
         call sll_s_hdf5_ser_file_close(hfile_id, error)
      end if

#ifdef __OPENMP
!$omp end single
#endif

      print *, 'ITERATION', istep

      call sll_s_int2string(istep, cstep)
      call sll_s_hdf5_ser_file_create("ff"//cstep//".h5", hfile_id, error)
      call sll_o_hdf5_ser_write_array(hfile_id, ff, "ff", error)
      call sll_o_hdf5_ser_write_array(hfile_id, ff1, "ff1", error)
      call sll_o_hdf5_ser_write_array(hfile_id, fg, "fg", error)
      call sll_s_hdf5_ser_file_close(hfile_id, error)
   end do

   !compute fg on the fine mesh -> ff1 (for diagnostic)
   do i = 1, Ncx + 1
      call sll_s_cubic_spline_1d_compute_interpolant(fg(i, :), interp_spline_v)
      call sll_s_cubic_spline_1d_eval_array(vh_array, ff1(i, :), Ncvh + 1, interp_spline_v)
   end do

   open (12, file="ffinalh")
   do i = 1, Ncx + 1
      do j = 1, Ncvh + 1
         write (12, *) x_array(i), vh_array(j), ff(i, j), ff1(i, j)
      end do
      write (12, *)
   end do
   close (12)

   open (12, file="ffinal")
   do i = 1, Ncx + 1
      do j = 1, Ncv + 1
         write (12, *) x_array(i), v_array(j), fg(i, j)
      end do
      write (12, *)
   end do
   close (12)

#ifdef __OPENMP
!$omp end parallel
#endif

   close (th_diag)
   close (ex_diag)

   print *, 'VP1D_deltaf_cart has exited normally'
contains
   elemental function f_equilibrium(v)
      sll_real64, intent(in) :: v
      sll_real64 :: f_equilibrium

      f_equilibrium = 1.0_f64/sqrt(2*sll_p_pi)*exp(-0.5_f64*v*v)
   end function f_equilibrium

end program sim_bsl_vp_1d1v_cart_micro_macro
