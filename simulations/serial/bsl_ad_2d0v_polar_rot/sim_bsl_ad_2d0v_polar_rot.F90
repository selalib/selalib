program sim_bsl_ad_2d0v_polar_rot

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_ascii_io, only: &
    sll_s_ascii_file_create

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_hermite, &
    sll_p_periodic

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d

  use sll_m_common_coordinate_transformations, only: &
    sll_f_polar_jac11, &
    sll_f_polar_jac12, &
    sll_f_polar_jac21, &
    sll_f_polar_jac22, &
    sll_f_polar_x1, &
    sll_f_polar_x2

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_coordinate_transformations_2d, only: &
    sll_f_new_coordinate_transformation_2d_analytic

  use sll_m_cubic_spline_interpolator_2d, only: &
    sll_f_new_cubic_spline_interpolator_2d

  use sll_m_hermite_interpolation_2d, only: &
    sll_p_hermite_c0, &
    sll_p_hermite_dirichlet, &
    sll_p_hermite_periodic

  use sll_m_hermite_interpolator_2d, only: &
    sll_f_new_hermite_interpolator_2d

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  use sll_m_xdmf, only: &
    sll_s_plot_f

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  sll_int32 :: num_cells1
  sll_int32 :: num_cells2
  sll_real64 :: t_init
  sll_real64 :: t_end
  !sll_int32 :: spline_degree
  sll_real64 :: dt
  sll_int32 :: number_iterations
  sll_int32 :: freq_diag
  sll_int32 :: freq_diag_time
  sll_real64 :: gauss_x1
  sll_real64 :: gauss_x2
  sll_real64 :: gauss_sig
  sll_real64 :: gauss_amp
  sll_int32 :: p
  
  
  sll_int32 :: ierr
  sll_real64, allocatable :: rho_tn(:,:)
  sll_real64, allocatable :: rho_tn1(:,:)
  sll_real64, allocatable :: rho_exact(:,:)
  sll_real64, allocatable :: charac_feet1(:,:)
  sll_real64, allocatable :: charac_feet2(:,:)
  sll_int32 :: step
  sll_int32 :: i
  !sll_real64 :: x
  !sll_real64 :: y
  !sll_real64 :: xx
  !sll_real64 :: yy
  !logical :: inside
  sll_real64 :: t
  sll_real64 :: linf_err
  character(len=256) :: num_method_case
  !sll_int32 :: num_method
  sll_int32 :: IO_stat
  sll_int32, parameter :: input_file = 99
  character(len = 256) :: input_filename
  character(len = 256) :: input_filename_loc
  sll_int32 :: count
  character(len=256) :: rho_name
  character(len=256) :: rho_error_name
  !character(len=256)  :: filename
  !character(len=4)   :: filenum
  logical :: use_num
  sll_int32 :: num_run
  character(len=256)  :: str_num_run
  sll_int32 :: thdiag_1d_id
  sll_int32 :: thdiag_0d_id
  character(len=256)  :: thdiag_1d_filename
  character(len=256)  :: thdiag_0d_filename
  character(len=256)  :: mesh_name
  sll_real64 :: linf_err_loc
  sll_real64 :: l1_err
  sll_real64 :: l2_err
  sll_real64 :: l1_err_loc
  sll_real64 :: l2_err_loc
  sll_real64 :: rho_min
  sll_real64 :: rho_min_loc
  sll_real64 :: time_compute_interpolant
  sll_real64 :: time_interpolate
  sll_real64 :: time_t0
  sll_real64 :: time_t1
  !sll_real64 :: cosdt
  !sll_real64 :: sindt
  sll_real64 :: delta1
  sll_real64 :: delta2
  sll_real64 :: rmin
  sll_real64 :: rmax
  class(sll_c_interpolator_2d), pointer :: interp2d
  class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
  type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
  sll_int32 :: j
  
  
  namelist /geometry/ &
   num_cells1, &
   num_cells2, &
   rmin, &
   rmax
   
   
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
    !spline_degree, &
    p, &
    num_method_case


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
  
  rmin = 1.e-5_f64
  rmax = 8._f64
  num_cells1 = 80
  num_cells2 = 80
  
  !spline_degree = 2
  number_iterations = 300
  dt = 0.1_f64
  freq_diag = 10
  freq_diag_time = 1
  p = 6
  num_method_case = "SLL_CUBIC_SPLINES"

  gauss_x1  = 2._f64
  gauss_x2  = 2._f64
  gauss_sig = 1._f64/( 2._f64 * sqrt(2._f64)) 
  gauss_amp = 1.0_f64

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
      input_filename_loc = trim(input_filename)//".nml"
    else
      use_num = .true. 
      read(str_num_run, *) num_run
      input_filename_loc = trim(input_filename)//"_"//trim(str_num_run)//".nml"        
    endif
    open(unit = input_file, file=trim(input_filename_loc),IOStat=IO_stat)
     if( IO_stat /= 0 ) then
       SLL_ERROR('rotation_2d_hexagonal_hermite,','can not open file')
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
  print *,'  #num_cells1/2=',num_cells1,num_cells2
  print *,'  #rmin/rmax=',rmin,rmax

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
  !print *,'  #spline_degree=',spline_degree
  print *,'  #p=',p
  
  
  if(use_num)then
    rho_name = "rho_"//trim(str_num_run)//"_"
    rho_error_name = "rho_error_"//trim(str_num_run)//"_"
    thdiag_1d_filename = "thdiag_1d_"//trim(str_num_run)//".dat"
    thdiag_0d_filename = "thdiag_0d_"//trim(str_num_run)//".dat"
    mesh_name = "polar_"//trim(str_num_run)
  else
    rho_name = "rho_"
    rho_error_name = "rho_error_"        
    thdiag_1d_filename = "thdiag_1d.dat"
    thdiag_0d_filename = "thdiag_0d.dat"
    mesh_name = "polar"
  endif

  ! ----------------------------
  ! Allocations and initialization of variables
  ! ----------------------------
  count = 0
  
  delta1 = (rmax-rmin)/real(num_cells1,f64)
  delta2 = 2._f64*sll_p_pi/real(num_cells2,f64)

  mesh_2d => sll_f_new_cartesian_mesh_2d( &
    num_cells1, &
    num_cells2, &
    eta1_min = rmin, &
    eta1_max = rmax, &
    eta2_min = 0._f64, &
    eta2_max = 2._f64*sll_p_pi)      


  transformation => sll_f_new_coordinate_transformation_2d_analytic( &
    "analytic_polar_transformation", &
    mesh_2d, &
    sll_f_polar_x1, &
    sll_f_polar_x2, &
    sll_f_polar_jac11, &
    sll_f_polar_jac12, &
    sll_f_polar_jac21, &
    sll_f_polar_jac22, &
    params=(/0._f64,0._f64,0._f64,0._f64/) )     

    call sll_plot_polar_init( &
      mesh_2d, &
      transformation, &
      mesh_name )    

    




  call sll_s_ascii_file_create(thdiag_1d_filename, thdiag_1d_id, ierr)


  select case (num_method_case)
    case ("SLL_CUBIC_SPLINES")
      interp2d => sll_f_new_cubic_spline_interpolator_2d( &
        num_cells1+1, &
        num_cells2+1, &
        rmin, &
        rmax, &
        0._f64, &
        2._f64*sll_p_pi, &
        sll_p_hermite, &
        sll_p_periodic, &
        const_eta1_min_slope = 0._f64, &
        const_eta1_max_slope = 0._f64)
    case ("sll_p_hermite")
      interp2d => sll_f_new_hermite_interpolator_2d( &
        num_cells1+1, &
        num_cells2+1, &
        rmin, &
        rmax, &
        0._f64, &
        2._f64*sll_p_pi, &
        p, &          
        p, &          
        sll_p_hermite_c0, &
        sll_p_hermite_c0, &
        sll_p_hermite_dirichlet, &
        sll_p_hermite_periodic)
    case default    
      SLL_ERROR("rotation_2d_polar", "bad value of num_method_case")  
  end select
  
  
  SLL_ALLOCATE(rho_tn(num_cells1+1,num_cells2+1),ierr)
  SLL_ALLOCATE(rho_tn1(num_cells1+1,num_cells2+1),ierr)
  SLL_ALLOCATE(rho_exact(num_cells1+1,num_cells2+1),ierr)
  SLL_ALLOCATE(charac_feet1(num_cells1+1,num_cells2+1),ierr)
  SLL_ALLOCATE(charac_feet2(num_cells1+1,num_cells2+1),ierr)
  
  

  ! ----------------------------
  ! Initialization of solution
  ! ----------------------------

  call gaussian_at_t_polar( &
    rho_tn, &
    0._f64, &
    num_cells1, &
    num_cells2, &
    rmin, &
    rmax, &
    gauss_x1, &
    gauss_x2, &
    gauss_sig, &
    gauss_amp )
  
  rho_min_loc = minval(rho_tn)
  rho_min = rho_min_loc

  ! ----------------------------
  ! compute characteristics (are the same)
  ! ----------------------------
    do j=1,num_cells2+1
      do i=1,num_cells1+1
        charac_feet1(i,j) = rmin+real(i-1,f64)*delta1
        charac_feet2(i,j) = (real(j-1,f64)*delta2+dt)/(2._f64*sll_p_pi)
        charac_feet2(i,j) = charac_feet2(i,j)-real(floor(charac_feet2(i,j)),f64)
        charac_feet2(i,j) = charac_feet2(i,j)*(2._f64*sll_p_pi)                
      enddo
    enddo

  
  ! ----------------------------
  ! Time loop
  ! ----------------------------
  
  
  do step=1,number_iterations

    ! ----------------------------
    ! Compute interpolants
    ! ----------------------------

    call cpu_time(time_t0)

    call cpu_time(time_t1)
    time_compute_interpolant = time_compute_interpolant+time_t1-time_t0 

    ! ----------------------------
    ! Interpolation at vertices
    ! ----------------------------
    call cpu_time(time_t0)
    
    !cosdt = cos(dt)
    !sindt = sin(dt)
    
!    do j=1,num_cells2+1
!      do i=1,num_cells1+1
!        x = r*cos
!        
!      enddo
!    enddo
!    do i=1, n_points
!      x = mesh%cartesian_coord(1,i)
!      y = mesh%cartesian_coord(2,i)
!      xx = x*cosdt - y*sindt
!      yy = x*sindt + y*cosdt
!    enddo 
    call interp2d%interpolate_array( &
      num_cells1+1, &
      num_cells2+1, &
      rho_tn, &
      charac_feet1, &
      charac_feet2, &
      rho_tn1)      


 

    call cpu_time(time_t1)

    ! ----------------------------
    ! Update solution
    ! ----------------------------

    rho_tn = rho_tn1
    time_interpolate = time_interpolate+time_t1-time_t0 


    ! ----------------------------
    ! Error computation
    ! ----------------------------
    
    t = real(step,f64)*dt
    
    if(modulo(step,freq_diag_time)==0)then

      call gaussian_at_t_polar( &
        rho_exact, &
        t, &
        num_cells1, &
        num_cells2, &
        rmin, &
        rmax, &
        gauss_x1, &
        gauss_x2, &
        gauss_sig, &
        gauss_amp)
      linf_err_loc = maxval(abs(rho_exact-rho_tn))
      linf_err = max(linf_err,linf_err_loc)
      rho_min_loc = minval(rho_tn)
      rho_min = min(rho_min_loc,rho_min)
      l1_err_loc = sum(abs(rho_exact-rho_tn))*delta1*delta2
      l1_err = max(l1_err_loc,l1_err)
      l2_err_loc = sqrt(sum((rho_exact-rho_tn)**2)*delta1*delta2)
      l2_err = max(l2_err_loc,l2_err)

      write(thdiag_1d_id,*) &
        t, &
        rho_min, &
        l1_err_loc, &
        l2_err_loc, &
        linf_err_loc
              
    endif

    if(modulo(step,freq_diag)==0)then      
      count = count+1
      call gaussian_at_t_polar( &
        rho_exact, &
        t, &
        num_cells1, &
        num_cells2, &
        rmin, &
        rmax, &
        gauss_x1, &
        gauss_x2, &
        gauss_sig, &
        gauss_amp)
      linf_err_loc = maxval(abs(rho_exact-rho_tn))
      linf_err = max(linf_err,linf_err_loc)
      rho_min_loc = minval(rho_tn)
      rho_min = min(rho_min_loc,rho_min)
      l1_err_loc = sum(abs(rho_exact-rho_tn))*delta1*delta2
      l1_err = max(l1_err_loc,l1_err)
      l2_err_loc = sqrt(sum((rho_exact-rho_tn)**2)*delta1*delta2)
      l2_err = max(l2_err_loc,l2_err)



      print *,'#time,count,rho_min',t,count,rho_min
      print *,'#err (l1,l2,linf)', &
        l1_err_loc, &
        l2_err_loc, &
        linf_err_loc
      call sll_s_plot_f( &
        count, &
        rho_tn, &  
        num_cells1+1, &
        num_cells2+1, &
        rho_name, &
        mesh_name, &
        t )    
      call sll_s_plot_f( &
        count, &
        rho_exact-rho_tn, &  
        num_cells1+1, &
        num_cells2+1, &
        rho_error_name, &
        mesh_name, &
        t )    
    endif

      
  enddo
  
  
  
  close(thdiag_1d_id)
  
  
  call cpu_time(t_end)
  
  print *, &
    '#cpu time,interpo-lant/late', &
    t_end-t_init, &
    time_compute_interpolant, &
    time_interpolate
  
  
  print *, &
    '#efficiency', &
    real(num_cells1*num_cells2,f64)*real(number_iterations,f64) &
    /(1.e6_f64*(t_end-t_init)), &
    real(num_cells1*num_cells2,f64)*real(number_iterations,f64) &
    /(1.e6_f64*(time_compute_interpolant+time_interpolate))

    
    
    !real(use_edge,f64),real(use_tri,f64)   
  print *, &
    '#max of err (l1,l2,linf)', &
    l1_err, &
    l2_err, &
    linf_err

  call sll_s_ascii_file_create(thdiag_0d_filename, thdiag_0d_id, ierr)

  write(thdiag_0d_id,*) &
    t_end-t_init, &
    rho_min, &
    l1_err, &
    l2_err, &
    linf_err, &
    time_compute_interpolant, &
    time_interpolate

  close(thdiag_0d_id)
  
  
  print *,"#PASSED"

contains

  subroutine gaussian_at_t_polar( &
    f_t, &
    t, &
    num_cells1, &
    num_cells2, &
    rmin, &
    rmax, &
    center_x1, &
    center_x2, &
    sigma, &
    amplitude)
    sll_real64, intent(inout) :: f_t(:,:)
    sll_real64, intent(in) :: t
    sll_real64, intent(in) :: rmin
    sll_real64, intent(in) :: rmax
    sll_int32, intent(in) :: num_cells1
    sll_int32, intent(in) :: num_cells2
    sll_real64, intent(in) :: center_x1
    sll_real64, intent(in) :: center_x2
    sll_real64, intent(in) :: sigma
    sll_real64, intent(in) :: amplitude
    
    sll_real64 :: x
    sll_real64 :: y
    sll_int32  :: i
    sll_int32  :: j
    sll_real64 :: xx
    sll_real64 :: yy
    sll_real64 :: r
    sll_real64 :: theta

    do j=1,num_cells2+1
      theta = 2._f64*sll_p_pi*real(j-1,f64)/real(num_cells2,f64) 
      do i=1,num_cells1+1
        r = rmin+real(i-1,f64)*(rmax-rmin)/real(num_cells1,f64)
        x = r*cos(theta)
        y = r*sin(theta)
        xx = x*cos(t) - y*sin(t)
        yy = x*sin(t) + y*cos(t)       
        f_t(i,j) = amplitude * exp(-0.5_f64* &
            ((xx-center_x1)**2 + (yy-center_x2)**2) / sigma**2 )
      enddo
    enddo
    
  end subroutine gaussian_at_t_polar

  subroutine sll_plot_polar_init( &
    mesh_2d, &
    transf, &
    mesh_name )
    
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    class(sll_c_coordinate_transformation_2d_base), pointer :: transf
    character(len=*), intent(in) :: mesh_name 
    
    sll_real64, allocatable :: x1(:,:)
    sll_real64, allocatable :: x2(:,:)
    sll_real64, allocatable :: f(:,:)
    sll_int32 :: ierr
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: num_pts1
    sll_int32 :: num_pts2
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: eta1
    sll_real64 :: eta2

    num_pts1 = mesh_2d%num_cells1+1
    num_pts2 = mesh_2d%num_cells2+1
    eta1_min = mesh_2d%eta1_min
    eta1_max = mesh_2d%eta1_max
    eta2_min = mesh_2d%eta2_min
    eta2_max = mesh_2d%eta2_max
    delta1 = mesh_2d%delta_eta1
    delta2 = mesh_2d%delta_eta2
    
    SLL_ALLOCATE(x1(num_pts1,num_pts2),ierr)
    SLL_ALLOCATE(x2(num_pts1,num_pts2),ierr)
    SLL_ALLOCATE(f(num_pts1,num_pts2),ierr)
    
    f = 0._f64
    
    do j=1,num_pts2
      eta2 = eta2_min+real(j-1,f64)*delta2
      do i=1,num_pts1
        eta1 = eta1_min+real(i-1,f64)*delta1
        x1(i,j) = transf%x1(eta1,eta2)
        x2(i,j) = transf%x2(eta1,eta2) 
      enddo
    enddo
    call sll_s_plot_f( &
      0, &
      f, &  
      num_pts1, &
      num_pts2,  &
      "f", & !dummy (for sll_plt_f, we should be able to
      !initialize only the mesh TODO)
      mesh_name, &
      0._f64, &
      x1, &
      x2)    
        
  end subroutine sll_plot_polar_init

  
  
end program sim_bsl_ad_2d0v_polar_rot
