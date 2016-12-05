program sim_bsl_ad_2d0v_hex_rot

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_ascii_io, only: &
    sll_s_ascii_file_create

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet

  use sll_m_box_splines, only: &
    sll_f_hex_interpolate_value, &
    sll_f_new_box_spline_2d, &
    sll_t_box_spline_2d

  use sll_m_constants, only: &
       sll_p_sqrt3, &
       sll_p_pi

  use sll_m_hexagonal_meshes, only: &
    sll_s_get_cell_vertices_index, &
    sll_f_new_hex_mesh_2d, &
    sll_t_hex_mesh_2d

  use sll_m_utilities, only: &
    sll_s_int2string

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: SLL_LINEAR_ADV = 1
  sll_int32, parameter :: SLL_CIRCULAR_ADV = 2
  sll_int32, parameter :: SLL_NOTHING = 0
  sll_int32, parameter :: input_file = 99

  type(sll_t_hex_mesh_2d),   pointer        :: mesh
  type(sll_t_box_spline_2d), pointer        :: spline

  sll_real64, dimension(:), allocatable :: dist_tn
  sll_real64, dimension(:), allocatable :: dist_tn1
  sll_real64, dimension(:), allocatable :: dist_exact

  logical :: inside
  logical :: use_num

  character(len = 256) :: num_method_case
  character(len = 256) :: input_filename
  character(len = 256) :: input_filename_loc
  character(len = 256) :: dist_name
  character(len = 256) :: dist_exact_name
  character(len = 256) :: dist_error_name
  character(len = 256) :: filename
  character(len = 256) :: str_num_run
  character(len = 256) :: thdiag_1d_filename
  character(len = 256) :: thdiag_0d_filename
  character(len = 4)   :: filenum

  sll_int32 :: num_cells
  sll_int32 :: spline_degree
  sll_int32 :: number_iterations
  sll_int32 :: freq_diag
  sll_int32 :: freq_diag_time
  sll_int32 :: n_points
  sll_int32 :: step
  sll_int32 :: i
  sll_int32 :: ierr
  sll_int32 :: num_method
  sll_int32 :: IO_stat
  sll_int32 :: count
  sll_int32 :: num_run
  sll_int32 :: thdiag_1d_id
  sll_int32 :: thdiag_0d_id

  sll_real64 :: center_mesh_x1
  sll_real64 :: center_mesh_x2
  sll_real64 :: radius
  sll_real64 :: delta
  sll_real64 :: t_init
  sll_real64 :: t_end
  sll_real64 :: dt
  sll_real64 :: gauss_x1
  sll_real64 :: gauss_x2
  sll_real64 :: gauss_sig
  sll_real64 :: gauss_amp
  sll_real64 :: x
  sll_real64 :: y
  sll_real64 :: xx
  sll_real64 :: yy
  sll_real64 :: r11, r12
  sll_real64 :: r21, r22
  sll_real64 :: h1, h2
  sll_real64 :: a0, a1
  sll_real64 :: r, th
  sll_real64 :: t
  sll_real64 :: linf_err
  sll_real64 :: linf_err_loc
  sll_real64 :: l1_err
  sll_real64 :: l2_err
  sll_real64 :: l1_err_loc
  sll_real64 :: l2_err_loc
  sll_real64 :: dist_min
  sll_real64 :: dist_min_loc
  sll_real64 :: time_compute_interpolant
  sll_real64 :: time_interpolate
  sll_real64 :: time_t0
  sll_real64 :: time_t1
  sll_real64 :: cosdt
  sll_real64 :: sindt


  namelist /geometry/ &
       center_mesh_x1, &
       center_mesh_x2, &
       radius, &
       num_cells

  namelist /initial_function/ &
       gauss_x1,  &
       gauss_x2,  &
       gauss_sig, &
       gauss_amp

  namelist /time_iterations/ &
       dt, &
       number_iterations, &
       freq_diag, &
       freq_diag_time


  namelist /interpolation/ &
       spline_degree, &
       num_method_case, &
       a0, &
       a1


  call cpu_time(t_init)
  time_interpolate = 0._f64
  time_compute_interpolant = 0._f64

  ! ----------------------------
  ! Default parameters
  ! ----------------------------
  linf_err = 0._f64
  linf_err_loc = 0._f64
  l1_err = 0._f64
  l1_err_loc = 0._f64
  l2_err = 0._f64
  l2_err_loc = 0._f64

  center_mesh_x1 = 0._f64
  center_mesh_x2 = 0._f64
  radius = 8._f64
  num_cells = 80
  number_iterations = 300
  dt = 0.1_f64
  freq_diag = 10
  freq_diag_time = 1

  gauss_x1  = 2._f64
  gauss_x2  = 2._f64
  gauss_sig = 1._f64/( 2._f64 * sqrt(2._f64))
  gauss_amp = 1.0_f64

  num_method_case = "SLL_CIRCULAR_ADV"
  num_method = 0
  spline_degree = 2
  a0 = 0._f64
  a1 = 0._f64

  ! ----------------------------
  ! Reading from file
  ! ----------------------------
  call get_command_argument(1, input_filename)

  if (len_trim(input_filename) == 0)then
     use_num = .false.
     print *,'#initialization with default parameters'
  else
     call get_command_argument(2, str_num_run)
     if(len_trim(str_num_run) == 0)then
        use_num = .false.
        input_filename_loc = trim(input_filename)
     else
        use_num = .true.
        read(str_num_run, *) num_run
        input_filename_loc = &
             trim(input_filename)//"_"//trim(str_num_run)//".nml"
     endif
     open(unit = input_file, file=trim(input_filename_loc),IOStat=IO_stat)
     if( IO_stat /= 0 ) then
        SLL_ERROR('advection_2d_hexagonal', 'cannot open file')
     end if
     print *,'#initialization with filename:'
     print *,'#',trim(input_filename_loc)
     read(input_file, geometry)
     read(input_file, initial_function)
     read(input_file, time_iterations)
     read(input_file, interpolation)
     close(input_file)
  endif

  print *,'#parameters are'
  print *,'#geometry:'
  print *,'  #center_mesh_x1=',center_mesh_x1
  print *,'  #center_mesh_x2=',center_mesh_x2
  print *,'  #radius=',radius
  print *,'  #num_cells=',num_cells

  print *,'#initial_function:'
  print *,'  #gauss_x1=',gauss_x1
  print *,'  #gauss_x2=',gauss_x2
  print *,'  #gauss_sig=',gauss_sig
  print *,'  #gauss_amp=',gauss_amp

  print *,'#time_iterations:'
  print *,'  #dt=',dt
  print *,'  #number_iterations=',number_iterations
  print *,'  #freq_diag=',freq_diag
  print *,'  #freq_diag_time=',freq_diag_time

  print *,'#interpolation:'
  print *,'  #num_method_case=',trim(num_method_case)
  print *,'  #spline_degree=',spline_degree
  print *,'  #a0, a1 = ', a0, a1


  if(use_num)then
     dist_name = "dist_"//trim(str_num_run)//"_"
     dist_error_name = "dist_error_"//trim(str_num_run)//"_"
     dist_exact_name = "dist_exact_"//trim(str_num_run)//"_"
     thdiag_1d_filename = "thdiag_1d_"//trim(str_num_run)//".dat"
     thdiag_0d_filename = "thdiag_0d_"//trim(str_num_run)//".dat"
  else
     dist_name = "dist_"
     dist_error_name = "dist_error_"
     dist_exact_name = "dist_exact_"
     thdiag_1d_filename = "thdiag_1d.dat"
     thdiag_0d_filename = "thdiag_0d.dat"
  endif

  ! ----------------------------
  ! Allocations and initialization of variables
  ! ----------------------------
  count = 0

  call sll_s_ascii_file_create(thdiag_1d_filename, thdiag_1d_id, ierr)

  select case (num_method_case)
  case ("SLL_LINEAR_ADV")
     num_method = SLL_LINEAR_ADV
  case ("SLL_CIRCULAR_ADV")
     num_method = SLL_CIRCULAR_ADV
  case ("SLL_NOTHING")
     num_method = SLL_NOTHING
  case default
     SLL_ERROR("advection_2d_hexagonal", "bad value of num_method_case")
  end select

  mesh => sll_f_new_hex_mesh_2d( &
       num_cells, &
       center_mesh_x1, &
       center_mesh_x2,&
       radius=radius)

  n_points = mesh%num_pts_tot
  delta = mesh%delta
  r11 = mesh%r1_x1
  r12 = mesh%r1_x2
  r21 = mesh%r2_x1
  r22 = mesh%r2_x2

  spline => sll_f_new_box_spline_2d(mesh, sll_p_dirichlet)

  SLL_ALLOCATE(dist_tn(n_points),ierr)
  SLL_ALLOCATE(dist_tn1(n_points),ierr)
  SLL_ALLOCATE(dist_exact(n_points),ierr)

  ! ----------------------------
  ! Initialization of solution
  ! ----------------------------

  call gaussian_at_t( &
       dist_tn, &
       0._f64, &
       mesh, &
       gauss_x1, &
       gauss_x2, &
       gauss_sig, &
       gauss_amp, a0, a1)

  dist_min_loc = minval(dist_tn)
  dist_min = dist_min_loc

  ! ----------------------------
  ! Time loop
  ! ----------------------------
  cosdt = cos(dt)
  sindt = sin(dt)

  do step=1,number_iterations

     ! ----------------------------
     ! Compute interpolants
     ! ----------------------------
     call cpu_time(time_t0)

     call spline%compute_coeff_box_spline_2d( dist_tn, spline_degree)

     call cpu_time(time_t1)

     time_compute_interpolant = time_compute_interpolant+time_t1-time_t0

     ! ----------------------------
     ! Interpolation at vertices
     ! ----------------------------
     call cpu_time(time_t0)

     do i=1, n_points
        x = mesh%cartesian_coord(1,i)
        y = mesh%cartesian_coord(2,i)

        select case (num_method)
        case(SLL_CIRCULAR_ADV)
           r = sqrt(x**2 + y**2)
           if (x .eq. 0._f64 .and. y .eq. 0._f64 ) then
              th = 0._f64
           else
              th = atan2(y, x)
           end if
           xx = r * cos(-2._f64*sll_p_pi*real(step, f64)*dt + th)
           yy = r * sin(-2._f64*sll_p_pi*real(step, f64)*dt + th)
           ! xx = x*cosdt - y*sindt
           ! yy = x*sindt + y*cosdt
        case(SLL_LINEAR_ADV)
           xx = x - a0 *  dt
           yy = y - a1 *  dt
        case default
           SLL_ERROR("advection_2d_hexagonal", "unknow num_method")
        end select

        inside = .true.
        h1 =  xx*r11 + yy*r12
        h2 =  xx*r21 + yy*r22
        if (abs(h1) > radius-delta .or. abs(h2) > radius-delta) inside = .false.
        if (abs(xx) > (radius-delta)*sll_p_sqrt3*0.5_f64) inside = .false.
        if ( inside ) then
           dist_tn1(i) = sll_f_hex_interpolate_value(mesh, &
                xx, yy, &
                spline, spline_degree)
        else
           dist_tn1(i) = 0._f64 ! dirichlet boundary condition
        endif
     enddo

     call cpu_time(time_t1)

     ! ----------------------------
     ! Update solution
     ! ----------------------------

     dist_tn = dist_tn1
     time_interpolate = time_interpolate+time_t1-time_t0


     ! ----------------------------
     ! Error computation
     ! ----------------------------

     t = real(step,f64)*dt

     if(modulo(step,freq_diag_time)==0)then

        dist_exact(1:n_points) = 0._f64
        call gaussian_at_t( &
             dist_exact, &
             t, &
             mesh, &
             gauss_x1, &
             gauss_x2, &
             gauss_sig, &
             gauss_amp, &
             a0, a1)

        linf_err_loc = maxval(abs(dist_exact-dist_tn))
        linf_err = max(linf_err,linf_err_loc)
        dist_min_loc = minval(dist_tn)
        dist_min = min(dist_min_loc,dist_min)
        l1_err_loc = sum(abs(dist_exact-dist_tn))*mesh%delta**2
        l1_err = max(l1_err_loc,l1_err)
        l2_err_loc = sqrt(sum((dist_exact-dist_tn)**2))*mesh%delta**2
        l2_err = max(l2_err_loc,l2_err)

        write(thdiag_1d_id,*) &
             t, &
             dist_min, &
             l1_err_loc, &
             l2_err_loc, &
             linf_err_loc

        count = count+1
        print *,'#time, diag_number, dist_min',t,count,dist_min
        print *,'#err (l1,l2,linf)', &
             l1_err_loc, &
             l2_err_loc, &
             linf_err_loc

        call sll_s_int2string(count,filenum)
        filename  = trim(dist_name)//trim(filenum)
        call mesh%sll_s_write_field_hex_mesh_xmf( dist_tn, trim(filename))
        filename  = trim(dist_error_name)//trim(filenum)
        call mesh%sll_s_write_field_hex_mesh_xmf(dist_tn - dist_exact, &
             trim(filename))
        filename  = trim(dist_exact_name)//trim(filenum)
        call mesh%sll_s_write_field_hex_mesh_xmf(dist_exact, &
             trim(filename))
     endif

  enddo ! End of time loop

  close(thdiag_1d_id)

  call cpu_time(t_end)

  print *, &
       '#cpu time,interpo-lant/late', &
       t_end-t_init, &
       time_compute_interpolant, &
       time_interpolate


  print *, &
       '#efficiency', &
       real(n_points,f64)*real(number_iterations,f64) &
       /(1.e6_f64*(t_end-t_init)), &
       real(n_points,f64)*real(number_iterations,f64) &
       /(1.e6_f64*(time_compute_interpolant+time_interpolate))

  print *, &
       '#max of err (l1,l2,linf)', &
       l1_err, &
       l2_err, &
       linf_err

  call sll_s_ascii_file_create(thdiag_0d_filename, thdiag_0d_id, ierr)

  write(thdiag_0d_id,*) &
       t_end-t_init, &
       dist_min, &
       l1_err, &
       l2_err, &
       linf_err, &
       time_compute_interpolant, &
       time_interpolate

  close(thdiag_0d_id)


  print *,"#PASSED"

contains

  subroutine gaussian_at_t( &
       f_t, &
       t, &
       mesh, &
       center_x1, &
       center_x2, &
       sigma, &
       amplitude, &
       a0, a1)

    implicit none
    type(sll_t_hex_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: f_t(:)
    sll_real64, intent(in) :: t
    sll_real64, intent(in) :: center_x1
    sll_real64, intent(in) :: center_x2
    sll_real64, intent(in) :: sigma
    sll_real64, intent(in) :: amplitude
    sll_real64, optional, intent(in) :: a0
    sll_real64, optional, intent(in) :: a1
    ! Local
    sll_real64 :: x
    sll_real64 :: y
    sll_int32  :: i
    sll_real64 :: xx
    sll_real64 :: yy
    sll_real64 :: r, th

    do i = 1,mesh%num_pts_tot
       x = mesh%cartesian_coord(1,i)
       y = mesh%cartesian_coord(2,i)

       select case (num_method)
        case(SLL_CIRCULAR_ADV)
           r = sqrt(x**2 + y**2)
           if (x .eq. 0._f64 .and. y .eq. 0._f64 ) then
              th = 0._f64
           else
              th = atan2(y, x)
           end if
           xx = r * cos(-2._f64*sll_p_pi*t + th)
           yy = r * sin(-2._f64*sll_p_pi*t + th)
        case(SLL_LINEAR_ADV)
           xx = x - a0 * t
           yy = y - a1 * t
        case default
           SLL_ERROR("advection_2d_hexagonal", "unknow num_method")
        end select

       ! xx = x*cos(t) - y*sin(t)
       ! yy = x*sin(t) + y*cos(t)
       f_t(i) = amplitude * exp(-0.5_f64* &
            ((xx-center_x1)**2 + (yy-center_x2)**2) / sigma**2 )
    enddo

  end subroutine gaussian_at_t


end program sim_bsl_ad_2d0v_hex_rot
