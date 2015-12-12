program sim_bsl_ad_2d0v_hex_rot

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_ascii_io, only: &
    sll_ascii_file_create

  use sll_m_boundary_condition_descriptors, only: &
    sll_dirichlet

  use sll_m_box_splines, only: &
    hex_interpolate_value, &
    new_box_spline_2d, &
    sll_box_spline_2d

  use sll_m_constants, only: &
    sll_sqrt3

  use sll_m_hermite_interpolation_2d, only: &
    compute_w_hermite

  use sll_m_hexagonal_meshes, only: &
    get_cell_vertices_index, &
    new_hex_mesh_2d, &
    sll_hex_mesh_2d

  use sll_m_interpolation_hex_hermite, only: &
    der_finite_difference, &
    hermite_interpolation

  use sll_m_utilities, only: &
    int2string

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: SLL_HEX_SPLINES = 2
  sll_int32, parameter :: SLL_HEX_Z9 = 9
  sll_int32, parameter :: SLL_HEX_Z10 = 10
  sll_int32, parameter :: SLL_HEX_HCTR = 11
  sll_int32, parameter :: SLL_HEX_HCTC = 12
  sll_int32, parameter :: SLL_HEX_GANEV_DIMITROV = 15
  sll_int32, parameter :: SLL_HEX_NOTHING = 0
  sll_int32, parameter :: SLL_HEX_P1 = 1
  sll_int32, parameter :: SLL_HEX_P1_new = 3
  sll_int32, parameter :: SLL_HEX_MITCHELL = 4
  sll_int32, parameter :: SLL_HEX_MITCHELL_new = 40
  sll_int32, parameter :: SLL_HEX_Z9_new = 90

  type(sll_hex_mesh_2d),   pointer        :: mesh
  type(sll_box_spline_2d), pointer        :: spline
  sll_real64, allocatable :: deriv(:,:)
  !sll_real64, allocatable :: deriv_new(:,:)
  sll_real64, dimension(:), allocatable :: center_values_tn
  sll_real64, dimension(:), allocatable :: edge_values_tn
  sll_real64, dimension(:), allocatable :: center_values_tn1
  sll_real64, dimension(:), allocatable :: edge_values_tn1
  sll_real64, dimension(:), allocatable :: center_values_exact
  sll_real64, dimension(:), allocatable :: edge_values_exact
  sll_int32 :: num_cells
  sll_real64 :: center_mesh_x1
  sll_real64 :: center_mesh_x2
  sll_real64 :: radius
  sll_real64 :: t_init
  sll_real64 :: t_end
  sll_int32 :: spline_degree
  sll_real64 :: dt
  sll_int32 :: number_iterations
  sll_int32 :: freq_diag
  sll_int32 :: freq_diag_time
  sll_real64 :: gauss_x1
  sll_real64 :: gauss_x2
  sll_real64 :: gauss_sig
  sll_real64 :: gauss_amp
  sll_int32 :: p


  sll_int32 :: n_points
  sll_int32 :: ierr
  sll_real64, dimension(:), allocatable :: rho_tn
  sll_real64, dimension(:), allocatable :: rho_tn1
  sll_real64, dimension(:), allocatable :: rho_exact
  sll_int32 :: step
  sll_int32 :: i
  sll_real64 :: x
  sll_real64 :: y
  sll_real64 :: xx
  sll_real64 :: yy
  sll_real64 :: r11
  sll_real64 :: r12
  sll_real64 :: r21
  sll_real64 :: r22
  sll_real64 :: det
  logical :: inside
  sll_real64 :: h1
  sll_real64 :: h2
  sll_real64 :: t
  sll_real64 :: linf_err
  sll_real64 :: linf_err_center
  sll_real64 :: linf_err_edge
  character(len=256) :: num_method_case
  sll_int32 :: num_method
  sll_real64 :: aire
  sll_int32 :: EXTRA_TABLES
  logical :: use_edge
  logical :: use_center
  sll_int32 :: IO_stat
  sll_int32, parameter :: input_file = 99
  character(len = 256) :: input_filename
  character(len = 256) :: input_filename_loc
  sll_int32 :: count
  character(len=256) :: rho_name
  character(len=256) :: rho_error_name
  character(len=256)  :: filename
  character(len=4)   :: filenum
  logical :: use_num
  sll_int32 :: num_run
  character(len=256)  :: str_num_run
  sll_int32 :: thdiag_1d_id
  sll_int32 :: thdiag_0d_id
  character(len=256)  :: thdiag_1d_filename
  character(len=256)  :: thdiag_0d_filename
  sll_real64 :: linf_err_loc
  sll_real64 :: linf_err_center_loc
  sll_real64 :: linf_err_edge_loc
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
  sll_int32 :: total_num_pts
  sll_real64 :: cosdt
  sll_real64 :: sindt
  sll_int32, allocatable :: index1(:,:)
  sll_int32, allocatable :: bounds1(:,:)
  sll_int32, allocatable :: bounds2(:,:)
  sll_real64, allocatable :: positions(:,:)
  sll_real64 :: r_vec(2,2)
  sll_int32, allocatable :: hex_stencil(:,:)
  sll_int32 :: size_hex_stencil

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
       p, &
       num_method_case


  call cpu_time(t_init)
  time_interpolate = 0._f64
  time_compute_interpolant = 0._f64


  ! ----------------------------
  ! Default parameters
  ! ----------------------------

  linf_err = 0._f64
  linf_err_center = 0._f64
  linf_err_edge = 0._f64
  linf_err_loc = 0._f64
  linf_err_center_loc = 0._f64
  linf_err_edge_loc = 0._f64
  l1_err = 0._f64
  l1_err_loc = 0._f64
  l2_err = 0._f64
  l2_err_loc = 0._f64

  center_mesh_x1 = 0._F64
  center_mesh_x2 = 0._F64
  radius = 8._f64
  num_cells = 80
  spline_degree = 2
  number_iterations = 300
  dt = 0.1_f64
  freq_diag = 10
  freq_diag_time = 1
  p = 6
  num_method_case = "SLL_HEX_Z9"

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
        SLL_ERROR('rotation_2d_hexagonal_hermite, &
             &','can not open file')
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
  print *,'  #p=',p


  if(use_num)then
     rho_name = "rho_"//trim(str_num_run)//"_"
     rho_error_name = "rho_error_"//trim(str_num_run)//"_"
     thdiag_1d_filename = "thdiag_1d_"//trim(str_num_run)//".dat"
     thdiag_0d_filename = "thdiag_0d_"//trim(str_num_run)//".dat"
  else
     rho_name = "rho_"
     rho_error_name = "rho_error_"        
     thdiag_1d_filename = "thdiag_1d.dat"
     thdiag_0d_filename = "thdiag_0d.dat"
  endif

  ! ----------------------------
  ! Allocations and initialization of variables
  ! ----------------------------
  count = 0

  call sll_ascii_file_create(thdiag_1d_filename, thdiag_1d_id, ierr)


  select case (num_method_case)
  case ("SLL_HEX_SPLINES")
     num_method = SLL_HEX_SPLINES
  case ("SLL_HEX_Z9")
     num_method = SLL_HEX_Z9  
  case ("SLL_HEX_Z9_new")
     num_method = SLL_HEX_Z9_new  
  case ("SLL_HEX_Z10")
     num_method = SLL_HEX_Z10  
  case ("SLL_HEX_HCTR")
     num_method = SLL_HEX_HCTR  
  case ("SLL_HEX_HCTC")
     num_method = SLL_HEX_HCTC  
  case ("SLL_HEX_GANEV_DIMITROV")
     num_method = SLL_HEX_GANEV_DIMITROV  
  case ("SLL_HEX_NOTHING")
     num_method = SLL_HEX_NOTHING  
  case ("SLL_HEX_P1")
     num_method = SLL_HEX_P1 
  case ("SLL_HEX_P1_new")
     num_method = SLL_HEX_P1_new 
  case ("SLL_HEX_MITCHELL")
     num_method = SLL_HEX_MITCHELL 
  case ("SLL_HEX_MITCHELL_new")
     num_method = SLL_HEX_MITCHELL_new 
  case default    
     SLL_ERROR("rotation_2d_hexagonal_hermite", "bad value of num_method_case")  
  end select


  if(num_method==SLL_HEX_GANEV_DIMITROV)then
     EXTRA_TABLES = 1
  else
     EXTRA_TABLES = 0  
  endif

  use_center = .false.
  use_edge = .false.

  if(num_method==SLL_HEX_Z10)then
     use_center = .true.
  endif
  if(num_method==SLL_HEX_GANEV_DIMITROV)then
     use_edge = .true.
  endif


  mesh => new_hex_mesh_2d( &
       num_cells, &
       center_mesh_x1, &
       center_mesh_x2,&
       radius=radius, &
       EXTRA_TABLES = EXTRA_TABLES )
  !call display_hex_mesh_2d(mesh)
  r_vec(1,1) = mesh%r1_x1
  r_vec(1,2) = mesh%r1_x2
  r_vec(2,1) = mesh%r2_x1
  r_vec(2,2) = mesh%r2_x2
  n_points = mesh%num_pts_tot
  det = (mesh%r1_x1*mesh%r2_x2-mesh%r1_x2*mesh%r2_x1)/mesh%delta
  r11 = + mesh%r2_x2/det
  r12 = - mesh%r2_x1/det
  r21 = - mesh%r1_x2/det
  r22 = + mesh%r1_x1/det
  aire = mesh%delta**2*sqrt(3._f64)*0.25_f64


  if(num_method==SLL_HEX_SPLINES)then
     spline => new_box_spline_2d(mesh, SLL_DIRICHLET)
  else if(num_method==SLL_HEX_MITCHELL)then
     SLL_ALLOCATE(deriv(13,n_points),ierr)  
  else if(num_method==SLL_HEX_MITCHELL_new)then
     SLL_ALLOCATE(deriv(13,n_points),ierr)  
  else if(num_method==SLL_HEX_Z9_new)then
     SLL_ALLOCATE(deriv(7,n_points),ierr)  
  else
     SLL_ALLOCATE(deriv(6,n_points),ierr)  
  endif
  SLL_ALLOCATE(rho_tn(n_points),ierr)
  SLL_ALLOCATE(rho_tn1(n_points),ierr)
  SLL_ALLOCATE(rho_exact(n_points),ierr)
  SLL_ALLOCATE(positions(2,n_points),ierr)


  if(use_center)then
     SLL_ALLOCATE(center_values_tn(mesh%num_triangles),ierr)
     SLL_ALLOCATE(center_values_tn1(mesh%num_triangles),ierr)
     SLL_ALLOCATE(center_values_exact(mesh%num_triangles),ierr)
  endif
  if(use_edge)then
     SLL_ALLOCATE(edge_values_tn(mesh%num_edges),ierr)
     SLL_ALLOCATE(edge_values_tn1(mesh%num_edges),ierr)
     SLL_ALLOCATE(edge_values_exact(mesh%num_edges),ierr)
  endif


  ! ----------------------------
  ! Initialization of solution
  ! ----------------------------

  call gaussian_at_t( &
       rho_tn, &
       0._f64, &
       mesh, &
       gauss_x1, &
       gauss_x2, &
       gauss_sig, &
       gauss_amp, &
       use_center=use_center, &
       use_edge=use_edge, &
       center_values_t=center_values_tn, &
       edge_values_t=edge_values_tn)

  rho_min_loc = minval(rho_tn)
  rho_min = rho_min_loc


  ! ----------------------------
  ! Steps concerning the mesh
  ! ----------------------------

  SLL_ALLOCATE(index1(2*num_cells+1,2*num_cells+1),ierr) 
  call get_numerotation_new(mesh,index1)
  SLL_ALLOCATE(bounds1(2,2*num_cells+1),ierr)
  call compute_bounds1(index1,bounds1,num_cells) 
  SLL_ALLOCATE(bounds2(2,2*num_cells+1),ierr)
  call compute_bounds2(index1,bounds2,num_cells) 


  select case (num_method)
  case(SLL_HEX_MITCHELL)
     size_hex_stencil = compute_size_hex_stencil( &
          index1, &
          num_cells, &
          p)
     SLL_ALLOCATE(hex_stencil(size_hex_stencil,n_points),ierr)  
     call compute_hex_stencil( &
          index1, &
          num_cells, &
          p, &
          hex_stencil)  
  case(SLL_HEX_Z9_new)
     size_hex_stencil = compute_size_hex_stencil( &
          index1, &
          num_cells, &
          p)
     SLL_ALLOCATE(hex_stencil(size_hex_stencil,n_points),ierr)  
     call compute_hex_stencil( &
          index1, &
          num_cells, &
          p, &
          hex_stencil)  
  case(SLL_HEX_MITCHELL_new)
     size_hex_stencil = compute_size_hex_stencil( &
          index1, &
          num_cells, &
          p)
     SLL_ALLOCATE(hex_stencil(size_hex_stencil,n_points),ierr)  
     call compute_hex_stencil_new( &
          index1, &
          num_cells, &
          p, &
          hex_stencil)  
  case default
  end select

  ! ----------------------------
  ! Time loop
  ! ----------------------------
  cosdt = cos(dt)
  sindt = sin(dt)

  do i=1,n_points
     x = mesh%cartesian_coord(1,i)
     y = mesh%cartesian_coord(2,i)
     positions(1,i) = x*cosdt - y*sindt
     positions(2,i) = x*sindt + y*cosdt                
  enddo

  do step=1,number_iterations

     ! ----------------------------
     ! Compute interpolants
     ! ----------------------------

     !rho_tn = 0._f64
     !rho_tn(n_points/2) = 1._f64

     call cpu_time(time_t0)
     select case (num_method)
     case (SLL_HEX_SPLINES)
        call spline%compute_coeff_box_spline_2d( rho_tn, spline_degree)
     case (SLL_HEX_NOTHING)
     case (SLL_HEX_P1)
     case (SLL_HEX_P1_new)
     case (SLL_HEX_MITCHELL)
        !        call compute_derivative( &
        !          index1, &
        !          bounds1, &
        !          bounds2, &
        !          num_cells, &
        !          p, &
        !          rho_tn, &
        !          deriv)
        call compute_derivative_mitchell( &
             index1, &
             hex_stencil, &
             num_cells, &
             p, &
             rho_tn, &
             deriv)
     case (SLL_HEX_MITCHELL_new)
        !        call compute_derivative( &
        !          index1, &
        !          bounds1, &
        !          bounds2, &
        !          num_cells, &
        !          p, &
        !          rho_tn, &
        !          deriv)
        call compute_derivative_mitchell_new( &
             index1, &
             hex_stencil, &
             num_cells, &
             p, &
             rho_tn, &
             deriv)
     case (SLL_HEX_Z9_new)
        !        call compute_derivative( &
        !          index1, &
        !          bounds1, &
        !          bounds2, &
        !          num_cells, &
        !          p, &
        !          rho_tn, &
        !          deriv)
        call compute_derivative_new( &
             index1, &
             hex_stencil, &
             num_cells, &
             p, &
             rho_tn, &
             deriv)
     case default
        call  der_finite_difference( rho_tn, p, mesh%delta, mesh, deriv) 
     end select
     !    if(num_method==SLL_HEX_SPLINES)then
     !      call spline%compute_coeff_box_spline_2d( rho_tn, spline_degree)
     !    elseif(num_method == SLL_HEX_NOTHING)then
     !    elseif(num_method == SLL_HEX_P1)then
     !    else
     !      call  der_finite_difference( rho_tn, p, mesh%delta, mesh, deriv)
     !    endif
     call cpu_time(time_t1)
     time_compute_interpolant = time_compute_interpolant+time_t1-time_t0 

     ! ----------------------------
     ! Interpolation at vertices
     ! ----------------------------
     call cpu_time(time_t0)


     select case (num_method)
     case (SLL_HEX_SPLINES)
        do i=1, n_points
           x = mesh%cartesian_coord(1,i)
           y = mesh%cartesian_coord(2,i)
           xx = x*cosdt - y*sindt
           yy = x*sindt + y*cosdt      
           inside = .true.
           h1 =  xx*r11 + yy*r12
           h2 =  xx*r21 + yy*r22
           if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
           if ( abs(xx) > (radius-mesh%delta)*sll_sqrt3*0.5_f64) inside = .false.      
           if ( inside ) then
              rho_tn1(i) = hex_interpolate_value(mesh, xx, yy, spline, spline_degree)
           else
              rho_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif
        enddo

     case (SLL_HEX_NOTHING)
        do i=1, n_points
           x = mesh%cartesian_coord(1,i)
           y = mesh%cartesian_coord(2,i)
           xx = x*cosdt - y*sindt
           yy = x*sindt + y*cosdt      
           inside = .true.
           h1 =  xx*r11 + yy*r12
           h2 =  xx*r21 + yy*r22
           if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
           if ( abs(xx) > (radius-mesh%delta)*sll_sqrt3*0.5_f64) inside = .false.      
           if ( inside ) then
              rho_tn1(i) = rho_tn(i)
           else
              rho_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif
        enddo


     case (SLL_HEX_P1)
        do i=1, n_points
           x = mesh%cartesian_coord(1,i)
           y = mesh%cartesian_coord(2,i)
           xx = x*cosdt - y*sindt
           yy = x*sindt + y*cosdt      
           inside = .true.
           h1 =  xx*r11 + yy*r12
           h2 =  xx*r21 + yy*r22
           if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
           if ( abs(xx) > (radius-mesh%delta)*sll_sqrt3*0.5_f64) inside = .false.      
           if ( inside ) then
              call p1_interpolation( &
                   i, &
                   xx, &
                   yy, &
                   rho_tn, &
                   rho_tn1, &
                   mesh)
           else
              rho_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif
        enddo
     case (SLL_HEX_P1_new)
        call interpolate_p1_new( &
             radius, &
             num_cells, &
             index1, &
             rho_tn, &
             rho_tn1, &
             positions)
     case (SLL_HEX_MITCHELL)            
        call interpolate_mitchell( &
             radius, &
             num_cells, &
             index1, &
             deriv, &
             rho_tn1, &
             positions)
     case (SLL_HEX_MITCHELL_new)            
        call interpolate_mitchell_new( &
             radius, &
             num_cells, &
             index1, &
             deriv, &
             rho_tn1, &
             positions)
     case (SLL_HEX_Z9_new)            
        call interpolate_z9_new( &
             radius, &
             num_cells, &
             index1, &
             deriv, &
             rho_tn1, &
             positions)
     case default
        do i=1, n_points
           x = mesh%cartesian_coord(1,i)
           y = mesh%cartesian_coord(2,i)
           xx = x*cosdt - y*sindt
           yy = x*sindt + y*cosdt      
           inside = .true.
           h1 =  xx*r11 + yy*r12
           h2 =  xx*r21 + yy*r22
           if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
           if ( abs(xx) > (radius-mesh%delta)*sll_sqrt3*0.5_f64) inside = .false.      
           if ( inside ) then
              call hermite_interpolation( &
                   i, &
                   xx, &
                   yy, &
                   rho_tn, &
                   center_values_tn,&
                   edge_values_tn, &
                   rho_tn1, &
                   mesh, &
                   deriv, &
                   aire,& 
                   num_method)
           else
              rho_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif
        enddo

     end select
     !    do i=1, n_points
     !      x = mesh%cartesian_coord(1,i)
     !      y = mesh%cartesian_coord(2,i)
     !      xx = x*cosdt - y*sindt
     !      yy = x*sindt + y*cosdt      
     !      inside = .true.
     !      h1 =  xx*r11 + yy*r12
     !      h2 =  xx*r21 + yy*r22
     !      if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
     !      if ( abs(xx) > (radius-mesh%delta)*sll_sqrt3*0.5_f64) inside = .false.      
     !      if ( inside ) then
     !        if(num_method==SLL_HEX_SPLINES)then
     !          rho_tn1(i) = hex_interpolate_value(mesh, xx, yy, spline, spline_degree)
     !        elseif(num_method==SLL_HEX_NOTHING)then
     !          rho_tn1(i) = rho_tn(i)
     !        elseif(num_method==SLL_HEX_P1)then
     !          rho_tn1(i) = rho_tn(i)
     !          call p1_interpolation( &
     !          i, &
     !          xx, &
     !          yy, &
     !          rho_tn, &
     !          rho_tn1, &
     !          mesh)
     !        else
     !          call hermite_interpolation( &
     !          i, &
     !          xx, &
     !          yy, &
     !          rho_tn, &
     !          center_values_tn,&
     !          edge_values_tn, &
     !          rho_tn1, &
     !          mesh, &
     !          deriv, &
     !          aire,& 
     !          num_method)
     !          endif    
     !      else
     !        rho_tn1(i) = 0._f64 ! dirichlet boundary condition
     !      endif
     !    enddo





     ! ----------------------------
     ! Interpolation at center of triangles
     ! ----------------------------


     if(use_center)then
        do i=1, mesh%num_triangles
           x = mesh%center_cartesian_coord(1,i)
           y = mesh%center_cartesian_coord(2,i)
           xx = x*cos(dt) - y*sin(dt)
           yy = x*sin(dt) + y*cos(dt)

           inside = .true.
           h1 =  xx*r11 + yy*r12
           h2 =  xx*r21 + yy*r22

           if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
           if ( abs(xx) > (radius-mesh%delta)*sll_sqrt3*0.5_f64) inside = .false.

           if ( inside ) then
              !rho_tn1(i) = hex_interpolate_value(mesh, xx, yy, spline, spline_degree)

              call hermite_interpolation( &
                   i, &
                   xx, &
                   yy, &
                   rho_tn, &
                   center_values_tn,&
                   edge_values_tn, &
                   center_values_tn1, &
                   mesh, &
                   deriv, &
                   aire,& 
                   num_method)
           else
              center_values_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif
        enddo
     endif

     ! ----------------------------
     ! Interpolation at middle of edges
     ! ----------------------------


     if(use_edge)then
        do i=1, mesh%num_edges
           x = mesh%edge_center_cartesian_coord(1,i)
           y = mesh%edge_center_cartesian_coord(2,i)
           xx = x*cos(dt) - y*sin(dt)
           yy = x*sin(dt) + y*cos(dt)

           inside = .true.
           h1 =  xx*r11 + yy*r12
           h2 =  xx*r21 + yy*r22

           if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
           if ( abs(xx) > (radius-mesh%delta)*sll_sqrt3*0.5_f64) inside = .false.

           if ( inside ) then
              call hermite_interpolation( &
                   i, &
                   xx, &
                   yy, &
                   rho_tn, &
                   center_values_tn,&
                   edge_values_tn, &
                   edge_values_tn1, &
                   mesh, &
                   deriv, &
                   aire,& 
                   num_method)
           else
              edge_values_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif
        enddo
     endif
     call cpu_time(time_t1)

     ! ----------------------------
     ! Update solution
     ! ----------------------------

     rho_tn = rho_tn1
     if(use_center)then
        center_values_tn = center_values_tn1
     endif
     if(use_edge)then
        edge_values_tn = edge_values_tn1
     endif
     time_interpolate = time_interpolate+time_t1-time_t0 


     ! ----------------------------
     ! Error computation
     ! ----------------------------

     t = real(step,f64)*dt

     if(modulo(step,freq_diag_time)==0)then

        call gaussian_at_t( &
             rho_exact, &
             t, &
             mesh, &
             gauss_x1, &
             gauss_x2, &
             gauss_sig, &
             gauss_amp, &
             use_center=use_center, &
             use_edge=use_edge, &
             center_values_t=center_values_exact, &
             edge_values_t=edge_values_exact)
        linf_err_loc = maxval(abs(rho_exact-rho_tn))
        linf_err = max(linf_err,linf_err_loc)
        rho_min_loc = minval(rho_tn)
        rho_min = min(rho_min_loc,rho_min)
        l1_err_loc = sum(abs(rho_exact-rho_tn))*mesh%delta**2
        l1_err = max(l1_err_loc,l1_err)
        l2_err_loc = sqrt(sum((rho_exact-rho_tn)**2)*mesh%delta**2)
        l2_err = max(l2_err_loc,l2_err)


        if(use_center)then
           linf_err_center_loc = &
                maxval(abs(center_values_exact-center_values_tn))
           linf_err_center =  max( &
                linf_err_center, &
                linf_err_center_loc)
        endif
        if(use_edge)then
           linf_err_edge_loc = & 
                maxval(abs(edge_values_exact-edge_values_tn))
           linf_err_edge =  max( &
                linf_err_edge, &
                linf_err_edge_loc)
        endif

        write(thdiag_1d_id,*) &
             t, &
             rho_min, &
             l1_err_loc, &
             l2_err_loc, &
             linf_err_loc, &
             linf_err_center_loc, &
             linf_err_edge_loc

     endif

     if(modulo(step,freq_diag)==0)then
        count = count+1
        print *,'#time,count,rho_min',t,count,rho_min
        print *,'#err (l1,l2,linf)', &
             l1_err_loc, &
             l2_err_loc, &
             linf_err_loc
        call int2string(count,filenum)
        filename  = trim(rho_name)//trim(filenum)
        call mesh%write_field_hex_mesh_xmf( rho_tn, trim(filename))
        filename  = trim(rho_error_name)//trim(filenum)
        call mesh%write_field_hex_mesh_xmf( rho_exact-rho_tn, trim(filename))
     endif


  enddo



  close(thdiag_1d_id)


  call cpu_time(t_end)

  print *, &
       '#cpu time,interpo-lant/late', &
       t_end-t_init, &
       time_compute_interpolant, &
       time_interpolate

  total_num_pts = compute_num_tot_points( &
       num_cells, &
       use_edge=use_edge, &
       use_center=use_center)

  print *, &
       '#efficiency', &
       real(total_num_pts,f64)*real(number_iterations,f64) &
       /(1.e6_f64*(t_end-t_init)), &
       real(total_num_pts,f64)*real(number_iterations,f64) &
       /(1.e6_f64*(time_compute_interpolant+time_interpolate))



  !real(use_edge,f64),real(use_tri,f64)   
  print *, &
       '#max of err (l1,l2,linf)', &
       l1_err, &
       l2_err, &
       linf_err

  call sll_ascii_file_create(thdiag_0d_filename, thdiag_0d_id, ierr)

  write(thdiag_0d_id,*) &
       t_end-t_init, &
       rho_min, &
       l1_err, &
       l2_err, &
       linf_err, &
       linf_err_center, &
       linf_err_edge, &
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
       use_center, &
       use_edge, &
       center_values_t, &
       edge_values_t)
    implicit none
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: f_t(:)
    sll_real64, intent(in) :: t
    sll_real64, intent(in) :: center_x1
    sll_real64, intent(in) :: center_x2
    sll_real64, intent(in) :: sigma
    sll_real64, intent(in) :: amplitude
    logical, intent(in), optional :: use_center
    logical, intent(in), optional :: use_edge
    sll_real64, intent(inout), optional :: center_values_t(:)
    sll_real64, intent(inout), optional :: edge_values_t(:)

    sll_real64 :: x
    sll_real64 :: y
    sll_int32  :: i
    sll_real64 :: xx
    sll_real64 :: yy
    logical :: use_center_loc
    logical :: use_edge_loc

    use_center_loc = .false.
    use_edge_loc = .false.

    if(present(use_center))then
       use_center_loc = use_center
    endif
    if(present(use_edge))then
       use_edge_loc = use_edge
    endif

    do i = 1,mesh%num_pts_tot
       x = mesh%cartesian_coord(1,i)
       y = mesh%cartesian_coord(2,i)
       xx = x*cos(t) - y*sin(t)
       yy = x*sin(t) + y*cos(t)       
       f_t(i) = amplitude * exp(-0.5_f64* &
            ((xx-center_x1)**2 + (yy-center_x2)**2) / sigma**2 )
    enddo


    if(use_center_loc)then
       do i = 1,mesh%num_triangles
          x = mesh%center_cartesian_coord(1,i)
          y = mesh%center_cartesian_coord(2,i)
          xx = x*cos(t) - y*sin(t)
          yy = x*sin(t) + y*cos(t)       
          center_values_t(i) = amplitude * exp(-0.5_f64* &
               ((xx-center_x1)**2 + (yy-center_x2)**2) / sigma**2 )
       enddo
    endif

    if(use_edge_loc)then
       do i = 1,mesh%num_edges
          x = mesh%edge_center_cartesian_coord(1,i)
          y = mesh%edge_center_cartesian_coord(2,i)
          xx = x*cos(t) - y*sin(t)
          yy = x*sin(t) + y*cos(t)       
          edge_values_t(i) = amplitude * exp(-0.5_f64* &
               ((xx-center_x1)**2 + (yy-center_x2)**2) / sigma**2 )
       enddo
    endif

  end subroutine gaussian_at_t

  function compute_num_tot_points( &
       num_cells, &
       use_edge, &
       use_center) &
       result(res)
    implicit none

    sll_int32, intent(in) :: num_cells
    logical, intent(in), optional :: use_edge
    logical, intent(in), optional :: use_center
    sll_int32 :: res

    logical :: use_center_loc
    logical :: use_edge_loc
    sll_int32 :: use_center_int
    sll_int32 :: use_edge_int
    sll_int32 :: num_pts
    sll_int32 :: num_edges
    sll_int32 :: num_tri



    if(present(use_center))then
       use_center_loc = use_center
    else
       use_center_loc = .false.  
    endif
    if(present(use_edge))then
       use_edge_loc = use_edge
    else
       use_edge_loc = .false.  
    endif

    if(use_edge_loc)then
       use_edge_int = 1
    else
       use_edge_int = 0
    endif

    if(use_center_loc)then
       use_center_int = 1
    else
       use_center_int = 0
    endif


    num_pts=3*num_cells*(num_cells+1)+1
    num_edges=3*num_cells*(3*num_cells+1)
    num_tri=6*num_cells**2

    res = num_pts+use_center_int*num_tri
    res = res+use_edge_int*num_edges  
  end function compute_num_tot_points


  subroutine p1_interpolation( &
       num, &
       x, &
       y, &
       f_tn, &
       output_tn1, &
       mesh )
    implicit none
    sll_int32,intent(in) :: num
    sll_real64,intent(in) :: x
    sll_real64,intent(in) :: y
    sll_real64, intent(in) :: f_tn(:)
    sll_real64, intent(out) :: output_tn1(:)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_int32 :: i
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i3
    sll_real64 :: x1
    sll_real64 :: x2
    sll_real64 :: x3
    sll_real64 :: y1
    sll_real64 :: y2
    sll_real64 :: y3
    !sll_int32 :: k11
    !sll_int32 :: k12
    sll_real64 :: freedom(3)
    sll_real64 :: base(3)
    sll_real64 :: f
    sll_real64                 :: a2
    sll_real64                 :: x1x,x2x,x3x,y1y,y2y,y3y
    sll_real64 :: aire
    sll_real64                 :: l1, l2, l3

    aire = mesh%delta**2*sqrt(3._f64)*0.25_f64
    !return
    call get_cell_vertices_index( x, y, mesh, i1, i2, i3 )

    x1 = mesh%cartesian_coord(1,i1) 
    x2 = mesh%cartesian_coord(1,i2) 
    x3 = mesh%cartesian_coord(1,i3) 
    y1 = mesh%cartesian_coord(2,i1) 
    y2 = mesh%cartesian_coord(2,i2) 
    y3 = mesh%cartesian_coord(2,i3) 


    !k11 = mesh%hex_coord(1,i1) 
    !k12 = mesh%hex_coord(2,i1) 

    ! get the first 3 degrees of freedom

    ! values at the vertices of the triangle

    freedom(1) = f_tn(i1)
    freedom(2) = f_tn(i2)
    freedom(3) = f_tn(i3)

    a2  = 0.5_f64/aire
    y3y = y3 - y
    y2y = y2 - y
    y1y = y1 - y
    x3x = x3 - x
    x2x = x2 - x
    x1x = x1 - x

    l1   = a2 * abs( x2x*y3y - x3x*y2y )    ! barycentric coordinates
    l2   = a2 * abs( x1x*y3y - x3x*y1y ) 
    l3   = 1._f64 - l1 - l2


    base(1) = l1
    base(2) = l2
    base(3) = l3

    f = 0._f64

    do i = 1,3
       f = f + freedom(i)*base(i)
    enddo

    output_tn1(num) = f    



  end subroutine p1_interpolation

  subroutine get_numerotation(mesh,index1)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_int32 :: i
    sll_real64 :: err
    sll_real64 :: err_loc
    sll_int32, intent(out) :: index1(:,:)
    sll_int32, allocatable :: index2(:,:)
    sll_int32, allocatable :: index3(:,:)
    sll_int32, allocatable :: bounds1(:,:)
    !sll_int32, allocatable :: bounds2(:,:)
    !sll_int32, allocatable :: bounds3(:,:)
    sll_int32 :: ierr
    sll_int32 :: j
    sll_int32 :: mini
    !sll_int32 :: maxi
    sll_int32, allocatable :: check(:)
    sll_int32 :: num_cells

    num_cells = mesh%num_cells

    !check first the hexagonal coordinates
    print *,'#mesh'
    print *,'#center=',mesh%center_x1,mesh%center_x2
    print *,'#radius=',mesh%radius
    print *,'#num_cells=',mesh%num_cells
    print *,'#r1=',mesh%r1_x1/mesh%delta,mesh%r1_x2/mesh%delta
    print *,'#r2=',mesh%r2_x1/mesh%delta,mesh%r2_x2/mesh%delta
    print *,'#r3=',mesh%r3_x1/mesh%delta,mesh%r3_x2/mesh%delta
    print *,'#dot=',(mesh%r1_x1/mesh%delta)*(mesh%r2_x1/mesh%delta) &
         +(mesh%r1_x2/mesh%delta)*(mesh%r2_x2/mesh%delta)
    print *,'#det=',(mesh%r1_x1*mesh%r2_x2-mesh%r1_x2*mesh%r2_x1)/mesh%delta**2
    err = 0._f64
    err = max(err,abs((mesh%r1_x1*mesh%r2_x2-mesh%r1_x2*mesh%r2_x1)/mesh%delta**2 &
         -0.5_f64*sqrt(3._f64)))
    err = max(err,abs((mesh%r1_x1/mesh%delta)*(mesh%r2_x1/mesh%delta) &
         +(mesh%r1_x2/mesh%delta)*(mesh%r2_x2/mesh%delta)+0.5_f64))
    err = max(err,abs(sqrt((mesh%r1_x1/mesh%delta)**2 &
         +(mesh%r1_x2/mesh%delta)**2)-1._f64))
    err = max(err,abs(sqrt((mesh%r2_x1/mesh%delta)**2 &
         +(mesh%r2_x2/mesh%delta)**2)-1._f64))
    err = max(err,abs(mesh%delta-mesh%radius/mesh%num_cells))
    err = max(err,abs(mesh%r3_x1-mesh%r1_x1-mesh%r2_x1))
    err = max(err,abs(mesh%r3_x2-mesh%r1_x2-mesh%r2_x2))
    do i=1,mesh%num_pts_tot
       err = max(err, &
            abs(mesh%center_x1+mesh%r1_x1*mesh%hex_coord(1,i)+mesh%r2_x1*mesh%hex_coord(2,i) &
            -mesh%cartesian_coord(1,i)))
       err = max(err, &
            abs(mesh%center_x2+mesh%r1_x2*mesh%hex_coord(1,i)+mesh%r2_x2*mesh%hex_coord(2,i) &
            -mesh%cartesian_coord(2,i)))  
    enddo
    !print *,'#err=',err
    !we consider lines: 
    !  line1  spanned by r1
    !  line2  spanned by r2
    !  line3  spanned by -r3
    !  line1 inter line2 inter line3 = center = (0,0,0)
    !  line1 + k1*r2 = line1 - k1*(-r3) -> k1
    !  line2 + k2*(-r3) = line2 - k2*r1 -> k2
    !  line3 + k3*r1 = line3 - k3*r2    -> k3
    !  now (line1+k1*r2) inter (line2-k2*r1) = -k2*r1+k1*r2
    !  = -k2*r1+k1*(r3-r1) = (-k2-k1)*r1 - k1*(-r3)
    ! so k3 = -k2-k1 i.e. k1+k2+k3 = 0
    ! so -k2 = hex_k1, k1 = hex_k2, k3 = -k1-k2 = -hex_k2+hex_k1

    !    do i=1,mesh%num_pts_tot
    !      print *, &
    !        i, &
    !        mesh%hex_coord(2,i), &
    !        -mesh%hex_coord(1,i), &
    !        mesh%hex_coord(1,i)-mesh%hex_coord(2,i)
    !    enddo
    ! -2 -2 .. 0
    ! -1 -2 .. 1
    !  0 -2 .. 2
    !  1 -1 .. 2
    !  2  0 .. 2
    index1 = 0
    do i=1,mesh%num_pts_tot
       index1(num_cells+1+mesh%hex_coord(2,i),num_cells+1-mesh%hex_coord(1,i)) = i
    enddo


    !print *,'maxval=',maxval(index1(:,:))

    SLL_ALLOCATE(index2(-num_cells:num_cells,-num_cells:num_cells),ierr) 
    index2 = 0
    do i=1,mesh%num_pts_tot
       index2(-mesh%hex_coord(1,i),mesh%hex_coord(1,i)-mesh%hex_coord(2,i)) = i
    enddo


    SLL_ALLOCATE(index3(-num_cells:num_cells,-num_cells:num_cells),ierr) 
    index3 = 0
    do i=1,mesh%num_pts_tot
       index3(mesh%hex_coord(1,i)-mesh%hex_coord(2,i),mesh%hex_coord(2,i)) = i
    enddo


    SLL_ALLOCATE(bounds1(2,-num_cells:num_cells),ierr) 
    do j=-num_cells,num_cells
       i=-num_cells
       do while(index1(num_cells+1+i,num_cells+1+j)==0)
          i=i+1
       enddo
       bounds1(1,j) = i
       do while((i<=num_cells-1).and.(index1(num_cells+1+i,num_cells+1+j)/=0))
          i=i+1
       enddo
       if((i==num_cells).and.(index1(num_cells+1+i,num_cells+1+j)/=0))then
          i=i+1
       endif
       i=i-1
       bounds1(2,j) = i      
    enddo
    !print *,'#err=',err

    !check that we have all the points
    err_loc = 0._f64
    do i=-num_cells,num_cells
       err_loc = err_loc+real(bounds1(2,i)-bounds1(1,i)+1,f64)
    enddo
    err = max(err,abs(err_loc-mesh%num_pts_tot))
    !print *,'#err=',err

    mini = maxval(index1)
    do j=-num_cells,num_cells
       do i=bounds1(1,j),bounds1(2,j)
          mini = min(mini,index1(num_cells+1+i,num_cells+1+j))  
       enddo
    enddo
    !print *,'#mini/maxi=',mini,maxval(index1),mesh%num_pts_tot
    err = max(err,abs(real(mini-1,f64)))
    err = max(err,abs(real(maxval(index1)-mesh%num_pts_tot,f64)))
    !print *,'#err=',err
    SLL_ALLOCATE(check(mesh%num_pts_tot),ierr)

    check = 0
    do j=-num_cells,num_cells
       do i=bounds1(1,j),bounds1(2,j)
          check(index1(num_cells+1+i,num_cells+1+j)) = &
               check(index1(num_cells+1+i,num_cells+1+j))+1  
       enddo
    enddo

    check = check-1

    err = max(err,real(maxval(abs(check)),f64))

    print *,'#err on hexagonal coordinates=',err

  end subroutine get_numerotation

  subroutine interpolate_p1_new_old( &
       radius, &
       num_cells, &
       r, &
                                !mesh, &
       index1, &
       rho_tn, &
       rho_tn1, &
       positions )
    !type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64, intent(in) :: radius
    sll_int32, intent(in) :: num_cells
    sll_real64, intent(in) :: r(2,2)
    sll_int32, intent(in) :: index1(:,:)
    sll_real64, intent(in) :: rho_tn(:)
    sll_real64, intent(out) :: rho_tn1(:)
    sll_real64, intent(in) :: positions(:,:)
    !sll_real64, intent(in) :: dt

    sll_int32 :: i
    !sll_real64 :: x
    !sll_real64 :: y
    !logical :: true
    !sll_real64 :: cosdt
    !sll_real64 :: sindt
    sll_real64 :: xx
    sll_real64 :: yy
    sll_real64 :: det
    sll_real64 :: r11
    sll_real64 :: r12
    sll_real64 :: r21
    sll_real64 :: r22
    sll_real64 :: h1
    sll_real64 :: h2
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i3
    !sll_real64 :: x1
    !sll_real64 :: x2
    !sll_real64 :: x3
    !sll_real64 :: y1
    !sll_real64 :: y2
    !sll_real64 :: y3
    sll_int32 :: ii
    sll_int32 :: jj
    !sll_real64 :: xi
    sll_real64 :: tmp
    sll_real64 :: aire
    sll_real64 :: a2
    sll_int32 :: s
    sll_real64 :: freedom(3)
    sll_real64 :: base(3)
    sll_real64 :: f
    !sll_real64                 :: x1x,x2x,x3x,y1y,y2y,y3y
    sll_real64 :: val1
    sll_real64 :: val2
    sll_real64 :: r1_x1
    sll_real64 :: r1_x2
    sll_real64 :: r2_x1
    sll_real64 :: r2_x2
    sll_real64 :: delta
    sll_int32 :: num_pts_tot
    sll_real64 :: x2p
    sll_real64 :: x3p
    sll_real64 :: x2n
    sll_real64 :: x3n
    sll_real64 :: y2p
    sll_real64 :: y3p
    sll_real64 :: y2n
    sll_real64 :: y3n

    r1_x1 = r(1,1)
    r2_x2 = r(2,2)
    r1_x2 = r(1,2)
    r2_x1 = r(2,1)
    delta = radius/real(num_cells,f64)
    num_pts_tot = 3*num_cells*(num_cells+1)+1
    det = (r1_x1*r2_x2-r1_x2*r2_x1) !/mesh%delta
    r11 = + r2_x2/det
    r12 = - r2_x1/det
    r21 = - r1_x2/det
    r22 = + r1_x1/det
    !cosdt = cos(dt)
    !sindt = sin(dt)

    tmp = delta*sqrt(3._f64)*0.5_f64
    aire = delta**2*sqrt(3._f64)*0.25_f64
    a2  = 0.5_f64/aire
    val1 = radius/delta -1._f64
    val2 = (radius-delta)*sqrt(3._f64)*0.5_f64
    x2p = a2*r1_x1
    y2p = -a2*r1_x2
    x3p = -a2*(r1_x1+r2_x1)
    y3p = a2*(r1_x2+r2_x2)

    x3n = -a2*r2_x1
    y3n = a2*r2_x2
    x2n = a2*(r1_x1+r2_x1)
    y2n = -a2*(r1_x2+r2_x2)


    do i=1, num_pts_tot
       !x = mesh%cartesian_coord(1,i)
       !y = mesh%cartesian_coord(2,i)
       !xx = x*cosdt - y*sindt
       !yy = x*sindt + y*cosdt      
       xx = positions(1,i)
       yy = positions(2,i)
       inside = .true.
       h1 =  xx*r11 + yy*r12
       h2 =  xx*r21 + yy*r22
       if ( (abs(h1) >  val1) .or. (abs(h2) >  val1) .or. (abs(xx) > val2) ) then
          inside = .false.
       endif
       if ( inside ) then
          ii = floor(h1)
          jj = floor(h2)

          h1 = h1-real(ii,f64)
          h2 = h2-real(jj,f64)
          xx = r1_x1*h1+r2_x1*h2
          yy = r1_x2*h1+r2_x2*h2

          !ii = floor(h1/mesh%delta)
          !jj = floor(h2/mesh%delta)
          !ii = cart_to_hex1(mesh, xx, yy)
          !jj = cart_to_hex2(mesh, xx, yy)
          !jacob = r22 
          !jacob = mesh%r1_x1 * mesh%r2_x2 - mesh%r2_x1 * mesh%r1_x2
          !k1 = floor((mesh%r2_x2 * x1 - mesh%r2_x1 * x2)/jacob)
          !jacob = mesh%r1_x1 * mesh%r2_x2 - mesh%r2_x1 * mesh%r1_x2
          !k2 = floor((mesh%r1_x1 * x2 - mesh%r1_x2 * x1)/jacob)

          i1 = index1(num_cells+1+jj,num_cells+1-ii)
          !x1 = r1_x1*real(ii,f64)+r2_x1*real(jj,f64)
          !y1 = r1_x2*real(ii,f64)+r2_x2*real(jj,f64)
          !x1 = 0._f64
          !y1 = 0._f64

          !xi = ( real(ii,f64) - real(jj,f64) ) * tmp
          !xi = 0._f64
          !if ( xx >= xi ) then
          if ( xx >= 0._f64 ) then
             !i1 = hex_to_global(mesh,ii,jj) 
             i2 = index1(num_cells+1+jj,num_cells+1-ii-1)
             !x2 = r1_x1*real(ii+1,f64)+r2_x1*real(jj,f64)
             !y2 = r1_x2*real(ii+1,f64)+r2_x2*real(jj,f64)
             !x2 = a2*r1_x1
             !y2 = -a2*r1_x2
             i3 = index1(num_cells+1+jj+1,num_cells+1-ii-1)
             !x3 = r1_x1*real(ii+1,f64)+r2_x1*real(jj+1,f64)
             !y3 = r1_x2*real(ii+1,f64)+r2_x2*real(jj+1,f64)
             !x3 = -a2*(r1_x1+r2_x1)
             !y3 = a2*(r1_x2+r2_x2)
             base(2) = y3p*xx+x3p*yy
             base(3) = y2p*xx+x2p*yy
             !i2 = hex_to_global(mesh,ii+1,jj) 
             !i3 = hex_to_global(mesh,ii+1,jj+1) 
          else
             !i1 = hex_to_global(mesh,ii,jj) 
             i3 = index1(num_cells+1+jj+1,num_cells+1-ii)
             !x3 = r1_x1*real(ii,f64)+r2_x1*real(jj+1,f64)
             !y3 = r1_x2*real(ii,f64)+r2_x2*real(jj+1,f64)
             !x3 = -a2*r2_x1
             !y3 = a2*r2_x2
             i2 = index1(num_cells+1+jj+1,num_cells+1-ii-1)
             !x2 = r1_x1*real(ii+1,f64)+r2_x1*real(jj+1,f64)
             !y2 = r1_x2*real(ii+1,f64)+r2_x2*real(jj+1,f64)
             !x2 = a2*(r1_x1+r2_x1)
             !y2 = -a2*(r1_x2+r2_x2)
             !i2 = hex_to_global(mesh,ii,jj+1)
             !i3 = hex_to_global(mesh,ii+1,jj+1)
             base(2) = y3n*xx+x3n*yy
             base(3) = y2n*xx+x2n*yy
          endif


          !        if ( xx > xi ) then
          !          i1 = hex_to_global(mesh,ii,jj) 
          !          i2 = hex_to_global(mesh,ii+1,jj) 
          !          i3 = hex_to_global(mesh,ii+1,jj+1) 
          !        else if ( xx < xi ) then
          !          i1 = hex_to_global(mesh,ii,jj) 
          !          i2 = hex_to_global(mesh,ii,jj+1)
          !          i3 = hex_to_global(mesh,ii+1,jj+1)
          !        else if ( xx == xi ) then
          !          if (xx < 0) then
          !            i1 = hex_to_global(mesh,ii,jj) 
          !            i2 = hex_to_global(mesh,ii+1,jj) 
          !            i3 = hex_to_global(mesh,ii+1,jj+1) 
          !          elseif (xx >= 0) then
          !            i1 = hex_to_global(mesh,ii,jj)
          !            i2 = hex_to_global(mesh,ii,jj+1) 
          !            i3 = hex_to_global(mesh,ii+1,jj+1) 
          !          endif
          !        endif


          !call get_cell_vertices_index( xx, yy, mesh, i1, i2, i3 )

          !x1 = mesh%r1_x1*mesh%hex_coord(1,i1)+mesh%r2_x1*mesh%hex_coord(2,i1)
          !x2 = mesh%r1_x1*mesh%hex_coord(1,i2)+mesh%r2_x1*mesh%hex_coord(2,i2)
          !x3 = mesh%r1_x1*mesh%hex_coord(1,i3)+mesh%r2_x1*mesh%hex_coord(2,i3)
          !y1 = mesh%r1_x2*mesh%hex_coord(1,i1)+mesh%r2_x2*mesh%hex_coord(2,i1)
          !y2 = mesh%r1_x2*mesh%hex_coord(1,i2)+mesh%r2_x2*mesh%hex_coord(2,i2)
          !y3 = mesh%r1_x2*mesh%hex_coord(1,i3)+mesh%r2_x2*mesh%hex_coord(2,i3)

          !x1 = mesh%cartesian_coord(1,i1) 
          !x2 = mesh%cartesian_coord(1,i2) 
          !x3 = mesh%cartesian_coord(1,i3) 
          !y1 = mesh%cartesian_coord(2,i1) 
          !y2 = mesh%cartesian_coord(2,i2) 
          !y3 = mesh%cartesian_coord(2,i3) 


          freedom(1) = rho_tn(i1)
          freedom(2) = rho_tn(i2)
          freedom(3) = rho_tn(i3)
          !y3y = y3 - yy
          !y2y = y2 - yy
          !y1y = y1 - yy
          !x3x = x3 - xx
          !x2x = x2 - xx
          !x1x = x1 - xx

          ! barycentric coordinates

          !base(1)   = a2 * ( x2x*y3y - x3x*y2y )    
          !base(2)   = a2 * ( x3x*y1y - x1x*y3y ) 
          !base(3)   = 1._f64 - base(1) - base(2)

          !base(1) = a2*((x2-xx)*(y3-yy)-(x3-xx)*(y2-yy))
          !base(2) = y3*xx+x3*yy
          !base(3) = y2*xx+x2*yy
          base(1) = 1._f64-base(2)-base(3)
          !base(3) = 1._f64-base(1)-base(2)
          f = 0._f64
          do s = 1,3
             f = f + freedom(s)*base(s)
          enddo

          rho_tn1(i) = f    



          !call p1_interpolation( &
          !	  i, &
          !	  xx, &
          !	  yy, &
          !	  rho_tn, &
          !	  rho_tn1, &
          !	  mesh)
       else
          rho_tn1(i) = 0._f64 ! dirichlet boundary condition
       endif
    enddo

  end subroutine interpolate_p1_new_old

  subroutine get_numerotation_new(mesh,index1)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_int32, intent(out) :: index1(:,:)
    index1 = 0
    do i=1,mesh%num_pts_tot
       index1(num_cells+1+mesh%hex_coord(1,i),num_cells+1+mesh%hex_coord(2,i)) = i
    enddo

  end subroutine get_numerotation_new

  subroutine interpolate_p1_new( &
       radius, &
       num_cells, &
                                !mesh, &
       index1, &
       rho_tn, &
       rho_tn1, &
       positions, &
       r )
    sll_real64, intent(in) :: radius
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: index1(:,:)
    sll_real64, intent(in) :: rho_tn(:)
    sll_real64, intent(out) :: rho_tn1(:)
    sll_real64, intent(in) :: positions(:,:)
    sll_real64, intent(in), optional :: r(2,2)
    !sll_real64, intent(in) :: dt

    sll_int32 :: i
    sll_real64 :: xx
    sll_real64 :: yy
    sll_real64 :: det
    sll_real64 :: r11
    sll_real64 :: r12
    sll_real64 :: r21
    sll_real64 :: r22
    sll_real64 :: h1
    sll_real64 :: h2
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i3
    sll_int32 :: ii
    sll_int32 :: jj
    sll_real64 :: a2
    sll_int32 :: s
    sll_real64 :: freedom(3)
    sll_real64 :: base(3)
    sll_real64 :: f
    sll_real64 :: r1_x1
    sll_real64 :: r1_x2
    sll_real64 :: r2_x1
    sll_real64 :: r2_x2
    sll_real64 :: delta
    sll_int32 :: num_pts_tot
    sll_real64 :: x2p
    sll_real64 :: x3p
    sll_real64 :: x2n
    sll_real64 :: x3n
    sll_real64 :: y2p
    sll_real64 :: y3p
    sll_real64 :: y2n
    sll_real64 :: y3n



    delta = radius/real(num_cells,f64)
    r1_x1 = sqrt(3._f64) * 0.5_f64*delta
    r1_x2 = 0.5_f64*delta
    r2_x1 = -r1_x1
    r2_x2 = r1_x2

    if(present(r))then
       if(abs(r1_x1-r(1,1))>1.e-14)then
          print *,'#value for r1_x1 not supported'
          stop
       endif
       if(abs(r2_x2-r(2,2))>1.e-14)then
          print *,'#value for r2_x2 not supported'
          stop
       endif
       if(abs(r1_x2-r(1,2))>1.e-14)then
          print *,'#value for r1_x2 not supported'
          stop
       endif
       if(abs(r2_x1-r(2,1))>1.e-14)then
          print *,'#value for r2_x1 not supported'
          stop
       endif
    endif
    num_pts_tot = 3*num_cells*(num_cells+1)+1
    det = r1_x1*r2_x2-r1_x2*r2_x1
    r11 = r2_x2/det
    r12 = -r2_x1/det
    r21 = -r1_x2/det
    r22 = r1_x1/det

    a2  = 2._f64/(delta**2*sqrt(3._f64))
    x2p = a2*r1_x1
    y2p = -a2*r1_x2
    x3p = -a2*(r1_x1+r2_x1)
    y3p = a2*(r1_x2+r2_x2)

    x3n = -a2*r2_x1
    y3n = a2*r2_x2
    x2n = a2*(r1_x1+r2_x1)
    y2n = -a2*(r1_x2+r2_x2)


    do i=1,num_pts_tot
       xx = positions(1,i)
       yy = positions(2,i)
       h1 =  xx*r11 + yy*r12
       h2 =  xx*r21 + yy*r22
       ii = floor(h1)
       jj = floor(h2)        
       h1 = h1-real(ii,f64)
       h2 = h2-real(jj,f64)
       if(ii>=num_cells)then
          ii=num_cells-1
          h1 = 1._f64
       endif
       if(ii<=-num_cells)then
          ii=-num_cells
          h1 = 0._f64
       endif
       if(jj>=num_cells)then
          jj=num_cells-1
          h2 = 1._f64
       endif
       if(jj<=-num_cells)then
          jj=-num_cells
          h2 = 0._f64
       endif

       xx = r1_x1*h1+r2_x1*h2
       yy = r1_x2*h1+r2_x2*h2
       i1 = index1(num_cells+1+ii,num_cells+1+jj)
       if ( xx >= 0._f64 ) then
          i2 = index1(num_cells+1+ii+1,num_cells+1+jj)          
          i3 = index1(num_cells+1+ii+1,num_cells+1+jj+1)
          base(2) = y3p*xx+x3p*yy
          base(3) = y2p*xx+x2p*yy
       else
          i2 = index1(num_cells+1+ii+1,num_cells+1+jj+1)
          i3 = index1(num_cells+1+ii,num_cells+1+jj+1)
          base(2) = y3n*xx+x3n*yy
          base(3) = y2n*xx+x2n*yy
       endif
       if(i1/=0)then
          freedom(1) = rho_tn(i1)
       else
          freedom(1) = 0._f64
       endif
       if(i2/=0)then
          freedom(2) = rho_tn(i2)
       else
          freedom(2) = 0._f64
       endif
       if(i3/=0)then
          freedom(3) = rho_tn(i3)
       else
          freedom(3) = 0._f64
       endif
       base(1) = 1._f64-base(2)-base(3)
       f = 0._f64
       do s = 1,3
          f = f + freedom(s)*base(s)
       enddo
       rho_tn1(i) = f    

    enddo

  end subroutine interpolate_p1_new

  subroutine compute_bounds1( &
       index, &
       bounds1, &
       num_cells)
    sll_int32, intent(in) :: index(:,:)
    sll_int32, intent(out) :: bounds1(:,:)
    sll_int32, intent(in) :: num_cells
    sll_int32 :: i
    sll_int32 :: j

    do j=-num_cells,num_cells
       i=-num_cells
       do while(index1(num_cells+1+i,num_cells+1+j)==0)
          i=i+1
       enddo
       bounds1(1,j+num_cells+1) = i
       do while((i<=num_cells-1).and.(index1(num_cells+1+i,num_cells+1+j)/=0))
          i=i+1
       enddo
       if((i==num_cells).and.(index1(num_cells+1+i,num_cells+1+j)/=0))then
          i=i+1
       endif
       i=i-1
       bounds1(2,j+num_cells+1) = i      
    enddo



  end subroutine compute_bounds1


  subroutine compute_bounds2( &
       index, &
       bounds2, &
       num_cells)
    sll_int32, intent(in) :: index(:,:)
    sll_int32, intent(out) :: bounds2(:,:)
    sll_int32, intent(in) :: num_cells
    sll_int32 :: i
    sll_int32 :: j

    do i=-num_cells,num_cells
       j=-num_cells
       do while(index1(num_cells+1+i,num_cells+1+j)==0)
          j=j+1
       enddo
       bounds2(1,i+num_cells+1) = j
       do while((j<=num_cells-1).and.(index1(num_cells+1+i,num_cells+1+j)/=0))
          j=j+1
       enddo
       if((j==num_cells).and.(index1(num_cells+1+i,num_cells+1+j)/=0))then
          j=j+1
       endif
       j=j-1
       bounds2(2,i+num_cells+1) = j      
    enddo



  end subroutine compute_bounds2


  subroutine compute_hex_coord(index,num_cells,hex_coord)
    sll_int32, intent(in) :: index(:,:)
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(out) :: hex_coord(:,:)

    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: num_pts_tot
    sll_int32 :: ind

    num_pts_tot = 3*num_cells*(num_cells+1)+1


    do j=1,2*num_cells+1
       do i=1,2*num_cells+1
          ind = index(i,j)
          if(ind/=0)then
             hex_coord(1,ind) = i
             hex_coord(2,ind) = j
          endif
       enddo
    enddo


  end subroutine compute_hex_coord

  function compute_size_hex_stencil( &
       index, &
       num_cells, &
       p) &
       result(res)
    sll_int32, intent(in) :: index(:,:)
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: p
    sll_int32 :: res

    sll_int32 :: r_left
    sll_int32 :: s_left
    sll_int32 :: r_right
    sll_int32 :: s_right

    r_left=-p/2
    s_left=(p+1)/2
    r_right=(-p+1)/2
    s_right=p/2+1

    res = s_left-r_left+1        
    res = res+s_right-r_right+1
    res = res*3

  end function compute_size_hex_stencil

  subroutine compute_hex_stencil( &
       index, &
       num_cells, &
       p, &
       stencil)
    sll_int32, intent(in) :: index(:,:)
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: p
    sll_int32, intent(out) :: stencil(:,:)


    sll_int32 :: r_left
    sll_int32 :: s_left
    sll_int32 :: r_right
    sll_int32 :: s_right
    !sll_int32 :: r
    !sll_int32 :: s
    sll_int32 :: num_pts_tot
    sll_int32, allocatable :: bounds1(:,:)
    sll_int32, allocatable :: hex_coord(:,:)
    sll_int32 :: ierr
    sll_int32 :: hex1
    sll_int32 :: hex2
    sll_int32 :: ii 
    sll_int32 :: hex1_loc
    sll_int32 :: hex2_loc
    sll_int32 :: j
    sll_int32 :: num
    sll_int32 :: jj
    sll_int32 :: tmp
    sll_int32 :: k
    sll_int32 :: a1(6)
    sll_int32 :: a2(6)
    sll_int32 :: r_tab(6)
    sll_int32 :: s_tab(6)


    num_pts_tot = 3*num_cells*(num_cells+1)+1

    SLL_ALLOCATE(bounds1(2,2*num_cells+1),ierr) 
    call compute_bounds1(index,bounds1,num_cells) 
    SLL_ALLOCATE(hex_coord(2,num_pts_tot),ierr) 
    call compute_hex_coord(index,num_cells,hex_coord)

    r_left=-p/2
    s_left=(p+1)/2
    r_right=(-p+1)/2
    s_right=p/2+1
    !r = min(r_left,r_right)
    !s = max(s_left,s_right)

    a1(1) = 1
    a2(1) = 0
    r_tab(1) = r_left
    s_tab(1) = s_left
    a1(2) = 1
    a2(2) = 0
    r_tab(2) = r_right
    s_tab(2) = s_right
    a1(3) = 0
    a2(3) = 1
    r_tab(3) = r_left
    s_tab(3) = s_left
    a1(4) = 0
    a2(4) = 1
    r_tab(4) = r_right
    s_tab(4) = s_right
    a1(5) = 1
    a2(5) = 1
    r_tab(5) = r_left
    s_tab(5) = s_left
    a1(6) = 1
    a2(6) = 1
    r_tab(6) = r_right
    s_tab(6) = s_right


    do i=1,num_pts_tot
       hex1 = hex_coord(1,i)  
       hex2 = hex_coord(2,i)
       ii=0

       do k=1,6
          do j=r_tab(k),s_tab(k)
             hex1_loc = hex1+a1(k)*j
             if(hex1_loc<1)then
                hex1_loc = 1
             endif
             if(hex1_loc>2*num_cells+1)then
                hex1_loc = 2*num_cells+1
             endif
             hex2_loc = hex2+a2(k)*j
             if(hex2_loc<1)then
                hex2_loc = 1
             endif
             if(hex2_loc>2*num_cells+1)then
                hex2_loc = 2*num_cells+1
             endif
             ii=ii+1
             stencil(ii,i) = index(hex1_loc,hex2_loc)          
          enddo
          ii=ii-(s_tab(k)-r_tab(k)+1)
          num=0
          do while (stencil(ii,i)==0)
             ii=ii+1
             num=num+1
          enddo
          tmp = stencil(ii,i)
          !ii=ii-(s_tab(k)-r_tab(k)+1)
          do jj=1,num
             !stencil(ii+jj-1,i) = tmp
             ii=ii-1
             stencil(ii,i) = tmp
          enddo
          ii=ii+(s_tab(k)-r_tab(k)+1)
          num = 0
          do while (stencil(ii,i)==0)
             ii=ii-1
             num=num+1
          enddo
          tmp = stencil(ii,i)
          do jj=1,num
             ii=ii+1
             stencil(ii,i) = tmp
          enddo
       enddo
    enddo

    !    do j=-num_cells,num_cells
    !      ii=0
    !      do i=bounds1(1,j+num_cells+1),bounds1(2,j+num_cells+1)
    !        ii=ii+1
    !        bufin1(ii) = rho_tn(index1(num_cells+1+i,num_cells+1+j))  
    !      enddo
    !      do i=r,0
    !        bufin1(i) = bufin1(1)
    !      enddo
    !      do i=ii+1,ii+s
    !        bufin1(i) = bufin1(ii)
    !      enddo
    !      do i=1,ii
    !        tmp=0._f64
    !        do jj = r_left,s_left
    !          tmp = tmp+w_left(jj)*bufin1(i+jj)  
    !        enddo
    !        ind = index1(num_cells+1+i+bounds1(1,j+num_cells+1)-1,num_cells+1+j)
    !        deriv(1,ind) = bufin1(i)
    !        deriv(2,ind) = tmp
    !        tmp=0._f64
    !        do jj = r_right,s_right
    !          tmp = tmp+w_right(jj)*bufin1(i+jj)  
    !        enddo
    !        deriv(3,ind) = tmp        
    !      enddo      
    !    enddo



  end subroutine compute_hex_stencil


  subroutine compute_hex_stencil_new( &
       index, &
       num_cells, &
       p, &
       stencil)
    sll_int32, intent(in) :: index(:,:)
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: p
    sll_int32, intent(out) :: stencil(:,:)


    sll_int32 :: r_left
    sll_int32 :: s_left
    sll_int32 :: r_right
    sll_int32 :: s_right
    !sll_int32 :: r
    !sll_int32 :: s
    sll_int32 :: num_pts_tot
    sll_int32, allocatable :: bounds1(:,:)
    sll_int32, allocatable :: hex_coord(:,:)
    sll_int32 :: ierr
    sll_int32 :: hex1
    sll_int32 :: hex2
    sll_int32 :: ii 
    sll_int32 :: hex1_loc
    sll_int32 :: hex2_loc
    sll_int32 :: j
    !sll_int32 :: num
    !sll_int32 :: jj
    !sll_int32 :: tmp
    sll_int32 :: k
    sll_int32 :: a1(6)
    sll_int32 :: a2(6)
    logical :: inside



    num_pts_tot = 3*num_cells*(num_cells+1)+1

    SLL_ALLOCATE(bounds1(2,2*num_cells+1),ierr) 
    call compute_bounds1(index,bounds1,num_cells) 
    SLL_ALLOCATE(hex_coord(2,num_pts_tot),ierr) 
    call compute_hex_coord(index,num_cells,hex_coord)

    r_left=-p/2
    s_left=(p+1)/2
    r_right=(-p+1)/2
    s_right=p/2+1
    !r = min(r_left,r_right)
    !s = max(s_left,s_right)

    a1(1) = 1
    a2(1) = 0
    a1(2) = 1
    a2(2) = 1
    a1(3) = 0
    a2(3) = 1
    a1(4) = -1
    a2(4) = 0
    a1(5) = -1
    a2(5) = -1
    a1(6) = 0
    a2(6) = -1


    do i=1,num_pts_tot
       hex1 = hex_coord(1,i)  
       hex2 = hex_coord(2,i)
       ii=0      
       do k=1,6
          do j=r_left,s_left
             hex1_loc = hex1+a1(k)*j
             inside = .true.
             if(hex1_loc<1)then
                inside = .false.
             endif
             if(hex1_loc>2*num_cells+1)then
                inside = .false.
             endif
             hex2_loc = hex2+a2(k)*j
             if(hex2_loc<1)then
                inside = .false.
             endif
             if(hex2_loc>2*num_cells+1)then
                inside = .false.
             endif
             ii=ii+1
             if(inside)then
                stencil(ii,i) = index(hex1_loc,hex2_loc)
             else
                stencil(ii,i) = 0
             endif
          enddo
       enddo
    enddo
  end subroutine compute_hex_stencil_new



  subroutine compute_hex_stencil_old( &
       index, &
       num_cells, &
       p, &
       stencil)
    sll_int32, intent(in) :: index(:,:)
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: p
    sll_int32, intent(out) :: stencil(:,:)


    sll_int32 :: r_left
    sll_int32 :: s_left
    sll_int32 :: r_right
    sll_int32 :: s_right
    !sll_int32 :: r
    !sll_int32 :: s
    sll_int32 :: num_pts_tot
    sll_int32, allocatable :: bounds1(:,:)
    sll_int32, allocatable :: hex_coord(:,:)
    sll_int32 :: ierr
    sll_int32 :: hex1
    sll_int32 :: hex2
    sll_int32 :: ii 
    sll_int32 :: hex1_loc
    sll_int32 :: hex2_loc
    sll_int32 :: j
    sll_int32 :: num
    sll_int32 :: jj
    sll_int32 :: tmp
    sll_int32 :: k
    sll_int32 :: a1(6)
    sll_int32 :: a2(6)


    num_pts_tot = 3*num_cells*(num_cells+1)+1

    SLL_ALLOCATE(bounds1(2,2*num_cells+1),ierr) 
    call compute_bounds1(index,bounds1,num_cells) 
    SLL_ALLOCATE(hex_coord(2,num_pts_tot),ierr) 
    call compute_hex_coord(index,num_cells,hex_coord)

    r_left=-p/2
    s_left=(p+1)/2
    r_right=(-p+1)/2
    s_right=p/2+1
    !r = min(r_left,r_right)
    !s = max(s_left,s_right)

    a1(1) = 1
    a2(1) = 0
    a1(2) = 1
    a2(2) = 1
    a1(3) = 0
    a2(3) = 1
    a1(4) = -1
    a2(4) = 0
    a1(5) = -1
    a2(5) = -1
    a1(6) = 0
    a2(6) = -1


    do i=1,num_pts_tot
       hex1 = hex_coord(1,i)  
       hex2 = hex_coord(2,i)
       ii=0

       do k=1,6
          do j=r_left,s_left
             hex1_loc = hex1+a1(k)*j
             if(hex1_loc<1)then
                hex1_loc = 1
             endif
             if(hex1_loc>2*num_cells+1)then
                hex1_loc = 2*num_cells+1
             endif
             hex2_loc = hex2+a2(k)*j
             if(hex2_loc<1)then
                hex2_loc = 1
             endif
             if(hex2_loc>2*num_cells+1)then
                hex2_loc = 2*num_cells+1
             endif
             ii=ii+1
             stencil(ii,i) = index(hex1_loc,hex2_loc)          
          enddo
          ii=ii-(s_left-r_left+1)
          num=0
          do while (stencil(ii,i)==0)
             ii=ii+1
             num=num+1
          enddo
          tmp = stencil(ii,i)
          !ii=ii-(s_tab(k)-r_tab(k)+1)
          do jj=1,num
             !stencil(ii+jj-1,i) = tmp
             ii=ii-1
             stencil(ii,i) = tmp
          enddo
          ii=ii+(s_left-r_left+1)
          num = 0
          do while (stencil(ii,i)==0)
             ii=ii-1
             num=num+1
          enddo
          tmp = stencil(ii,i)
          do jj=1,num
             ii=ii+1
             stencil(ii,i) = tmp
          enddo
       enddo
    enddo
  end subroutine compute_hex_stencil_old






  subroutine compute_derivative( &
       index, &
       bounds1, &
       bounds2, &
       num_cells, &
       p, &
       rho_tn, &
       deriv)
    sll_int32, intent(in) :: index(:,:)
    sll_int32, intent(in) :: bounds1(:,:)
    sll_int32, intent(in) :: bounds2(:,:)
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: p
    sll_real64, intent(in) :: rho_tn(:)
    sll_real64, intent(out) :: deriv(:,:)

    sll_int32 :: i
    sll_int32 :: r_left
    sll_int32 :: s_left
    sll_int32 :: r_right
    sll_int32 :: s_right
    sll_int32 :: r
    sll_int32 :: s
    sll_int32 :: ierr
    !sll_int32 :: p2
    sll_real64, allocatable :: w_left(:)
    sll_real64, allocatable :: w_right(:)
    sll_int32 :: num_pts_tot
    !sll_int32, allocatable :: bounds1(:,:)
    !sll_int32, allocatable :: bounds2(:,:)
    !sll_int32, allocatable :: bounds3(:,:)
    sll_int32 :: maxi1
    sll_real64, allocatable :: bufin1(:)
    sll_real64, allocatable :: bufout1(:)
    sll_int32 :: j
    sll_int32 :: ii
    sll_int32 :: jj
    sll_real64 :: tmp
    sll_int32 :: ind



    r_left=-p/2
    s_left=(p+1)/2
    SLL_ALLOCATE( w_left(r_left:s_left),ierr )
    call compute_w_hermite(w_left,r_left,s_left)


    r_right=(-p+1)/2
    s_right=p/2+1
    SLL_ALLOCATE( w_right(r_right:s_right),ierr )
    if((2*(p/2)-p)==0)then
       w_right(r_right:s_right) = w_left(r_left:s_left)
    else
       w_right(r_right:s_right) = -w_left(s_left:r_left:-1)
    endif

    r = min(r_left,r_right)
    s = max(s_left,s_right)

    num_pts_tot = 3*num_cells*(num_cells+1)+1

    !    do i=1,num_pts_tot
    !      deriv(1,i) = rho_tn(i)
    !    enddo




    !on cartesian mesh, we have for p odd
    ! 9*num_cells storage
    !f(0,0)
    !fx(0,0)
    !fx(1,0)
    !fy(0,0)
    !fy(0,1)
    !fxy(0,0)
    !fxy(1,0)
    !fxy(0,1)
    !fxy(1,1)
    !dof to add: 7
    !f(0,1),f(1,0),f(1,1)
    !fx(0,1),fx(1,1),fy(1,0),fy(1,1)
    !on cartesian mesh, we have for p even
    ! 4*num_cells storage
    !f(0,0)
    !fx(0,0)
    !fy(0,0)
    !fxy(0,0)
    !dof to add: 12 
    !f(1,0),fx(1,0),fy(1,0),fxy(1,0)
    !f(0,1),fx(0,1),fy(0,1),fxy(0,1)
    !f(1,1),fx(1,1),fy(1,1),fxy(1,1)

    !on hexagonal grid, we have for p odd
    ! 13*num_cells storage
    ! 8 for p and 8 for n
    !f(0,0) p+n
    !fr1(0,0) p
    !fr1(1,0) p
    !fr3(0,0) p+n
    !fr3(1,1) p+n
    !frrp(0,0) p
    !frrp(1,0) p
    !frrp(1,1) p
    !fr2(0,0) n
    !fr2(0,1) n
    !frrn(0,0) n
    !frrn(0,1) n
    !frrn(1,1) n
    !dof to add:
    ! 4 for p and 4 for n
    !f(1,0) p
    !f(1,1) p+n
    !f(0,1) n
    !fr2(1,0) p
    !fr2(1,1) p
    !fr1(0,1) n
    !fr1(1,1) n
    !on hexagonal grid, we have for p even
    ! 7*num_cells storage
    ! 5 for p and 4 for n
    !f(0,0) p+n
    !fr1(0,0) p
    !fr3(0,0) p+n
    !frrp(0,0) p
    !frrp(1,0) p
    !fr2(0,0) n
    !frrn(0,0) n
    !dof to add: 7 for p and 8 for n
    !f(0,1) n
    !f(1,0) p
    !f(1,1) p+n
    !fr2(0,1) n
    !fr2(1,0) p
    !fr2(1,1) p
    !fr1(0,1) n
    !fr1(1,0) p
    !fr1(1,1) n
    !fr3(1,1) p+n
    !frrp(1,1) p
    !frrn(0,1) n
    !frrn(1,1) n
    !          x
    !      x       x
    !          x
    !      x       x 
    !          x
    ! frrp(0,0) : r13
    ! frrp(1,0) : -r12
    ! frrp(1,1) : r23
    ! frrn(0,0) : r23
    ! frrn(0,1) : -r12
    ! frrn(1,1) : r13

    !for the moment, we use the same storage for p odd
    ! and even


    !SLL_ALLOCATE(bounds1(2,2*num_cells+1),ierr) 

    !    maxi1 = 0
    !    do j=-num_cells,num_cells
    !      maxi1 = max(maxi1,bounds1(2,j)-bounds1(1,j)+1)
    !    enddo

    maxi1 = 2*num_cells+1

    !print *,'#maxi1=',maxi1

    !return

    SLL_ALLOCATE(bufin1(r:maxi1+s),ierr)
    SLL_ALLOCATE(bufout1(maxi1),ierr)

    !deriv = 0._f64


    do j=-num_cells,num_cells
       ii=0
       do i=bounds1(1,j+num_cells+1),bounds1(2,j+num_cells+1)
          ii=ii+1
          bufin1(ii) = rho_tn(index1(num_cells+1+i,num_cells+1+j))  
       enddo
       do i=r,0
          bufin1(i) = bufin1(1)
       enddo
       do i=ii+1,ii+s
          bufin1(i) = bufin1(ii)
       enddo
       do i=1,ii
          tmp=0._f64
          do jj = r_left,s_left
             tmp = tmp+w_left(jj)*bufin1(i+jj)  
          enddo
          ind = index1(num_cells+1+i+bounds1(1,j+num_cells+1)-1,num_cells+1+j)
          deriv(1,ind) = bufin1(i)
          deriv(2,ind) = tmp
          tmp=0._f64
          do jj = r_right,s_right
             tmp = tmp+w_right(jj)*bufin1(i+jj)  
          enddo
          deriv(3,ind) = tmp        
       enddo
    enddo


    do i=-num_cells,num_cells
       ii=0
       do j=bounds2(1,i+num_cells+1),bounds2(2,i+num_cells+1)
          ii=ii+1
          bufin1(ii) = rho_tn(index1(num_cells+1+i,num_cells+1+j))  
       enddo
       do j=r,0
          bufin1(j) = bufin1(1)
       enddo
       do j=ii+1,ii+s
          bufin1(j) = bufin1(ii)
       enddo
       do j=1,ii
          tmp=0._f64
          do jj = r_left,s_left
             tmp = tmp+w_left(jj)*bufin1(j+jj)  
          enddo
          ind = index1(num_cells+1+i,num_cells+1+j+bounds2(1,i+num_cells+1)-1)
          !deriv(1,ind) = bufin1(j)
          deriv(4,ind) = tmp
          tmp=0._f64
          do jj = r_right,s_right
             tmp = tmp+w_right(jj)*bufin1(j+jj)  
          enddo
          deriv(5,ind) = tmp        
       enddo
    enddo







  end subroutine compute_derivative


  subroutine compute_derivative_new( &
       index, &
       stencil, &
       num_cells, &
       p, &
       rho_tn, &
       deriv)
    sll_int32, intent(in) :: index(:,:)
    sll_int32, intent(in) :: stencil(:,:)
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: p
    sll_real64, intent(in) :: rho_tn(:)
    sll_real64, intent(out) :: deriv(:,:)

    sll_int32 :: i
    sll_int32 :: r_left
    sll_int32 :: s_left
    sll_int32 :: r_right
    sll_int32 :: s_right
    sll_int32 :: ierr
    sll_real64, allocatable :: w_left(:)
    sll_real64, allocatable :: w_right(:)
    sll_int32 :: num_pts_tot
    sll_int32 :: j
    sll_int32 :: ii
    !sll_int32 :: jj
    sll_real64 :: tmp
    !sll_int32 :: ind




    r_left=-p/2
    s_left=(p+1)/2
    SLL_ALLOCATE( w_left(r_left:s_left),ierr )
    call compute_w_hermite(w_left,r_left,s_left)


    r_right=(-p+1)/2
    s_right=p/2+1
    SLL_ALLOCATE( w_right(r_right:s_right),ierr )
    if((2*(p/2)-p)==0)then
       w_right(r_right:s_right) = w_left(r_left:s_left)
    else
       w_right(r_right:s_right) = -w_left(s_left:r_left:-1)
    endif

    !print *,'#left',r_left,s_left,w_left
    !print *,'#right',r_right,s_right,w_right

    num_pts_tot = 3*num_cells*(num_cells+1)+1





    !on cartesian mesh, we have for p odd
    ! 9*num_cells storage
    !f(0,0)
    !fx(0,0)
    !fx(1,0)
    !fy(0,0)
    !fy(0,1)
    !fxy(0,0)
    !fxy(1,0)
    !fxy(0,1)
    !fxy(1,1)
    !dof to add: 7
    !f(0,1),f(1,0),f(1,1)
    !fx(0,1),fx(1,1),fy(1,0),fy(1,1)
    !on cartesian mesh, we have for p even
    ! 4*num_cells storage
    !f(0,0)
    !fx(0,0)
    !fy(0,0)
    !fxy(0,0)
    !dof to add: 12 
    !f(1,0),fx(1,0),fy(1,0),fxy(1,0)
    !f(0,1),fx(0,1),fy(0,1),fxy(0,1)
    !f(1,1),fx(1,1),fy(1,1),fxy(1,1)

    !on hexagonal grid, we have for p odd
    ! 13*num_cells storage
    ! 8 for p and 8 for n
    !f(0,0) p+n
    !fr1(0,0) p
    !fr1(1,0) p
    !fr3(0,0) p+n
    !fr3(1,1) p+n
    !frrp(0,0) p
    !frrp(1,0) p
    !frrp(1,1) p
    !fr2(0,0) n
    !fr2(0,1) n
    !frrn(0,0) n
    !frrn(0,1) n
    !frrn(1,1) n
    !dof to add:
    ! 4 for p and 4 for n
    !f(1,0) p
    !f(1,1) p+n
    !f(0,1) n
    !fr2(1,0) p
    !fr2(1,1) p
    !fr1(0,1) n
    !fr1(1,1) n
    !on hexagonal grid, we have for p even
    ! 7*num_cells storage
    ! 5 for p and 4 for n
    !f(0,0) p+n
    !fr1(0,0) p
    !fr3(0,0) p+n
    !frrp(0,0) p
    !frrp(1,0) p
    !fr2(0,0) n
    !frrn(0,0) n
    !dof to add: 7 for p and 8 for n
    !f(0,1) n
    !f(1,0) p
    !f(1,1) p+n
    !fr2(0,1) n
    !fr2(1,0) p
    !fr2(1,1) p
    !fr1(0,1) n
    !fr1(1,0) p
    !fr1(1,1) n
    !fr3(1,1) p+n
    !frrp(1,1) p
    !frrn(0,1) n
    !frrn(1,1) n
    !          x
    !      x       x
    !          x
    !      x       x 
    !          x
    ! frrp(0,0) : r13
    ! frrp(1,0) : -r12
    ! frrp(1,1) : r23
    ! frrn(0,0) : r23
    ! frrn(0,1) : -r12
    ! frrn(1,1) : r13

    !for the moment, we use the same storage for p odd
    ! and even



    do i=1,num_pts_tot
       deriv(1,i) =  rho_tn(i)
       ii=0

       tmp=0._f64
       do j=r_left,s_left
          ii=ii+1
          tmp = tmp+w_left(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(2,i) = tmp
       tmp=0._f64
       do j=r_right,s_right
          ii=ii+1
          tmp = tmp+w_right(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(3,i) = tmp

       tmp=0._f64
       do j=r_left,s_left
          ii=ii+1
          tmp = tmp+w_left(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(4,i) = tmp
       tmp=0._f64
       do j=r_right,s_right
          ii=ii+1
          tmp = tmp+w_right(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(5,i) = tmp

       tmp=0._f64
       do j=r_left,s_left
          ii=ii+1
          tmp = tmp+w_left(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(6,i) = tmp
       tmp=0._f64
       do j=r_right,s_right
          ii=ii+1
          tmp = tmp+w_right(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(7,i) = tmp

    enddo




  end subroutine compute_derivative_new



  subroutine compute_derivative_mitchell( &
       index, &
       stencil, &
       num_cells, &
       p, &
       rho_tn, &
       deriv)
    sll_int32, intent(in) :: index(:,:)
    sll_int32, intent(in) :: stencil(:,:)
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: p
    sll_real64, intent(in) :: rho_tn(:)
    sll_real64, intent(out) :: deriv(:,:)

    sll_int32 :: i
    sll_int32 :: r_left
    sll_int32 :: s_left
    sll_int32 :: r_right
    sll_int32 :: s_right
    sll_int32 :: ierr
    sll_real64, allocatable :: w_left(:)
    sll_real64, allocatable :: w_right(:)
    sll_int32 :: num_pts_tot
    sll_int32 :: j
    sll_int32 :: ii
    !sll_int32 :: jj
    sll_real64 :: tmp
    !sll_int32 :: ind

    sll_int32, allocatable :: loc_stencil(:,:)



    r_left=-p/2
    s_left=(p+1)/2
    SLL_ALLOCATE( w_left(r_left:s_left),ierr )
    call compute_w_hermite(w_left,r_left,s_left)


    r_right=(-p+1)/2
    s_right=p/2+1
    SLL_ALLOCATE( w_right(r_right:s_right),ierr )
    if((2*(p/2)-p)==0)then
       w_right(r_right:s_right) = w_left(r_left:s_left)
    else
       w_right(r_right:s_right) = -w_left(s_left:r_left:-1)
    endif

    !print *,'#left',r_left,s_left,w_left
    !print *,'#right',r_right,s_right,w_right

    !stop

    num_pts_tot = 3*num_cells*(num_cells+1)+1





    !on cartesian mesh, we have for p odd
    ! 9*num_cells storage
    !f(0,0)
    !fx(0,0)
    !fx(1,0)
    !fy(0,0)
    !fy(0,1)
    !fxy(0,0)
    !fxy(1,0)
    !fxy(0,1)
    !fxy(1,1)
    !dof to add: 7
    !f(0,1),f(1,0),f(1,1)
    !fx(0,1),fx(1,1),fy(1,0),fy(1,1)
    !on cartesian mesh, we have for p even
    ! 4*num_cells storage
    !f(0,0)
    !fx(0,0)
    !fy(0,0)
    !fxy(0,0)
    !dof to add: 12 
    !f(1,0),fx(1,0),fy(1,0),fxy(1,0)
    !f(0,1),fx(0,1),fy(0,1),fxy(0,1)
    !f(1,1),fx(1,1),fy(1,1),fxy(1,1)

    !on hexagonal grid, we have for p odd
    ! 13*num_cells storage
    ! 8 for p and 8 for n
    !f(0,0) p+n
    !fr1(0,0) p
    !fr1(1,0) p
    !fr3(0,0) p+n
    !fr3(1,1) p+n
    !frrp(0,0) p
    !frrp(1,0) p
    !frrp(1,1) p
    !fr2(0,0) n
    !fr2(0,1) n
    !frrn(0,0) n
    !frrn(0,1) n
    !frrn(1,1) n
    !dof to add:
    ! 4 for p and 4 for n
    !f(1,0) p
    !f(1,1) p+n
    !f(0,1) n
    !fr2(1,0) p
    !fr2(1,1) p
    !fr1(0,1) n
    !fr1(1,1) n
    !on hexagonal grid, we have for p even
    ! 7*num_cells storage
    ! 5 for p and 4 for n
    !f(0,0) p+n
    !fr1(0,0) p
    !fr3(0,0) p+n
    !frrp(0,0) p
    !frrp(1,0) p
    !fr2(0,0) n
    !frrn(0,0) n
    !dof to add: 7 for p and 8 for n
    !f(0,1) n
    !f(1,0) p
    !f(1,1) p+n
    !fr2(0,1) n
    !fr2(1,0) p
    !fr2(1,1) p
    !fr1(0,1) n
    !fr1(1,0) p
    !fr1(1,1) n
    !fr3(1,1) p+n
    !frrp(1,1) p
    !frrn(0,1) n
    !frrn(1,1) n
    !          x
    !      x       x
    !          x
    !      x       x 
    !          x
    ! frrp(0,0) : r13
    ! frrp(1,0) : -r12
    ! frrp(1,1) : r23
    ! frrn(0,0) : r23
    ! frrn(0,1) : -r12
    ! frrn(1,1) : r13

    !for the moment, we use the same storage for p odd
    ! and even



    do i=1,num_pts_tot
       deriv(1,i) =  rho_tn(i)
       ii=0

       tmp=0._f64
       do j=r_left,s_left
          ii=ii+1
          tmp = tmp+w_left(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(2,i) = tmp
       tmp=0._f64
       do j=r_right,s_right
          ii=ii+1
          tmp = tmp+w_right(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(3,i) = tmp

       tmp=0._f64
       do j=r_left,s_left
          ii=ii+1
          tmp = tmp+w_left(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(4,i) = tmp
       tmp=0._f64
       do j=r_right,s_right
          ii=ii+1
          tmp = tmp+w_right(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(5,i) = tmp

       tmp=0._f64
       do j=r_left,s_left
          ii=ii+1
          tmp = tmp+w_left(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(6,i) = tmp
       tmp=0._f64
       do j=r_right,s_right
          ii=ii+1
          tmp = tmp+w_right(j)*rho_tn(stencil(ii,i))        
       enddo
       deriv(7,i) = tmp

    enddo


    SLL_ALLOCATE( loc_stencil(r_left:s_left,3),ierr )

    !we have to do a second loop to compute the mixed
    !derivatives
    do i=1,num_pts_tot      

       !first take the three needed stencil 
       ii=0
       do j=r_left,s_left
          ii=ii+1
          loc_stencil(j,1) = stencil(ii,i)
       enddo
       ii = ii+s_right-r_right+1
       do j=r_left,s_left
          ii=ii+1
          loc_stencil(j,2) = stencil(ii,i)
       enddo
       ii = ii+s_right-r_right+1
       do j=r_left,s_left
          ii=ii+1
          loc_stencil(j,3) = stencil(ii,i)
       enddo

       !now, we compute the 6 mixed derivatives
       ii=0

       tmp=0._f64
       do j=r_left,s_left
          tmp = tmp+w_left(j)*deriv(2,loc_stencil(j,3))        
       enddo
       deriv(8,i) = tmp

       tmp=0._f64
       do j=r_left,s_left
          tmp = tmp+w_left(j)*deriv(3,loc_stencil(j,2))        
       enddo
       deriv(9,i) = tmp

       tmp=0._f64
       do j=r_left,s_left
          tmp = tmp+w_left(j)*deriv(4,loc_stencil(j,3))        
       enddo
       deriv(10,i) = tmp

       tmp=0._f64
       do j=r_left,s_left
          tmp = tmp+w_left(j)*deriv(5,loc_stencil(j,1))        
       enddo
       deriv(11,i) = tmp

       tmp=0._f64
       do j=r_left,s_left
          tmp = tmp+w_left(j)*deriv(7,loc_stencil(j,1))        
       enddo
       deriv(12,i) = tmp

       tmp=0._f64
       do j=r_left,s_left
          tmp = tmp+w_left(j)*deriv(7,loc_stencil(j,2))        
       enddo
       deriv(13,i) = tmp


    enddo




  end subroutine compute_derivative_mitchell





  subroutine compute_derivative_mitchell_new( &
       index, &
       stencil, &
       num_cells, &
       p, &
       rho_tn, &
       deriv)
    sll_int32, intent(in) :: index(:,:)
    sll_int32, intent(in) :: stencil(:,:)
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: p
    sll_real64, intent(in) :: rho_tn(:)
    sll_real64, intent(out) :: deriv(:,:)

    sll_int32 :: i
    sll_int32 :: r_left
    sll_int32 :: s_left
    !sll_int32 :: r_right
    !sll_int32 :: s_right
    sll_int32 :: ierr
    sll_real64, allocatable :: w_left(:)
    !sll_real64, allocatable :: w_right(:)
    sll_int32 :: num_pts_tot
    sll_int32 :: j
    sll_int32 :: ii
    !sll_int32 :: jj
    sll_real64 :: tmp
    !sll_int32 :: ind

    sll_int32 :: k



    r_left=-p/2
    s_left=(p+1)/2
    SLL_ALLOCATE( w_left(r_left:s_left),ierr )
    call compute_w_hermite(w_left,r_left,s_left)



    !print *,'#left',r_left,s_left,w_left

    num_pts_tot = 3*num_cells*(num_cells+1)+1


    do i=1,num_pts_tot
       deriv(1,i) =  rho_tn(i)
       ii=0

       do k=1,6
          tmp=0._f64
          do j=r_left,s_left
             ii=ii+1
             if(stencil(ii,i)/=0)then
                tmp = tmp+w_left(j)*rho_tn(stencil(ii,i))
             endif
          enddo
          deriv(k+1,i) = tmp
       enddo
    enddo

    !we have to do a second loop to compute the mixed
    !derivatives
    do i=1,num_pts_tot      
       ii=0
       do k=1,5
          tmp=0._f64
          do j=r_left,s_left
             ii=ii+1
             if(stencil(ii,i)/=0)then
                tmp = tmp+w_left(j)*deriv(k+2,stencil(ii,i))
             endif
          enddo
          deriv(k+7,i) = tmp
       enddo

       tmp=0._f64
       do j=r_left,s_left
          ii=ii+1
          if(stencil(ii,i)/=0)then
             tmp = tmp+w_left(j)*deriv(2,stencil(ii,i))
          endif
       enddo
       deriv(13,i) = tmp

    enddo




  end subroutine compute_derivative_mitchell_new




  subroutine interpolate_mitchell( &
       radius, &
       num_cells, &
       index1, &
       deriv, &
       rho_tn1, &
       positions, &
       r )
    sll_real64, intent(in) :: radius
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: index1(:,:)
    sll_real64, intent(in) :: deriv(:,:)
    sll_real64, intent(out) :: rho_tn1(:)
    sll_real64, intent(in) :: positions(:,:)
    sll_real64, intent(in), optional :: r(2,2)
    !sll_real64, intent(in) :: dt

    sll_int32 :: i
    sll_real64 :: xx
    sll_real64 :: yy
    sll_real64 :: det
    sll_real64 :: r11
    sll_real64 :: r12
    sll_real64 :: r21
    sll_real64 :: r22
    sll_real64 :: h1
    sll_real64 :: h2
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i3
    sll_int32 :: ii
    sll_int32 :: jj
    sll_real64 :: a2
    sll_int32 :: s
    sll_real64 :: freedom(12)
    sll_real64 :: base(12)
    sll_real64 :: lam(3)
    sll_real64 :: f
    sll_real64 :: r1_x1
    sll_real64 :: r1_x2
    sll_real64 :: r2_x1
    sll_real64 :: r2_x2
    sll_real64 :: delta
    sll_int32 :: num_pts_tot
    sll_real64 :: x2p
    sll_real64 :: x3p
    sll_real64 :: x2n
    sll_real64 :: x3n
    sll_real64 :: y2p
    sll_real64 :: y3p
    sll_real64 :: y2n
    sll_real64 :: y3n

    !sll_real64 :: err_loc
    sll_real64 :: err

    err = 0._f64

    delta = radius/real(num_cells,f64)
    r1_x1 = sqrt(3._f64) * 0.5_f64*delta
    r1_x2 = 0.5_f64*delta
    r2_x1 = -r1_x1
    r2_x2 = r1_x2

    if(present(r))then
       if(abs(r1_x1-r(1,1))>1.e-14)then
          print *,'#value for r1_x1 not supported'
          stop
       endif
       if(abs(r2_x2-r(2,2))>1.e-14)then
          print *,'#value for r2_x2 not supported'
          stop
       endif
       if(abs(r1_x2-r(1,2))>1.e-14)then
          print *,'#value for r1_x2 not supported'
          stop
       endif
       if(abs(r2_x1-r(2,1))>1.e-14)then
          print *,'#value for r2_x1 not supported'
          stop
       endif
    endif
    num_pts_tot = 3*num_cells*(num_cells+1)+1
    det = r1_x1*r2_x2-r1_x2*r2_x1
    r11 = r2_x2/det
    r12 = -r2_x1/det
    r21 = -r1_x2/det
    r22 = r1_x1/det

    a2  = 2._f64/(delta**2*sqrt(3._f64))
    x2p = a2*r1_x1
    y2p = -a2*r1_x2
    x3p = -a2*(r1_x1+r2_x1)
    y3p = a2*(r1_x2+r2_x2)

    x3n = -a2*r2_x1
    y3n = a2*r2_x2
    x2n = a2*(r1_x1+r2_x1)
    y2n = -a2*(r1_x2+r2_x2)


    do i=1,num_pts_tot
       xx = positions(1,i)
       yy = positions(2,i)
       h1 =  xx*r11 + yy*r12
       h2 =  xx*r21 + yy*r22
       ii = floor(h1)
       jj = floor(h2)        
       h1 = h1-real(ii,f64)
       h2 = h2-real(jj,f64)
       if(ii>=num_cells)then
          ii=num_cells-1
          h1 = 1._f64
       endif
       if(ii<=-num_cells)then
          ii=-num_cells
          h1 = 0._f64
       endif
       if(jj>=num_cells)then
          jj=num_cells-1
          h2 = 1._f64
       endif
       if(jj<=-num_cells)then
          jj=-num_cells
          h2 = 0._f64
       endif

       xx = r1_x1*h1+r2_x1*h2
       yy = r1_x2*h1+r2_x2*h2
       i1 = index1(num_cells+1+ii,num_cells+1+jj)
       if ( xx >= 0._f64 ) then
          i2 = index1(num_cells+1+ii+1,num_cells+1+jj)          
          i3 = index1(num_cells+1+ii+1,num_cells+1+jj+1)
          lam(2) = y3p*xx+x3p*yy
          lam(3) = y2p*xx+x2p*yy
       else
          i2 = index1(num_cells+1+ii+1,num_cells+1+jj+1)
          i3 = index1(num_cells+1+ii,num_cells+1+jj+1)
          lam(2) = y3n*xx+x3n*yy
          lam(3) = y2n*xx+x2n*yy
       endif
       freedom = 0._f64

       if(xx>=0.)then

          if(i1/=0)then
             freedom(1) = deriv(1,i1)
             freedom(4) = deriv(2,i1)
             freedom(5) = -deriv(3,i1)
             freedom(8) = -deriv(7,i1)          
             freedom(9) = deriv(6,i1)          
             freedom(10) = deriv(8,i1)          
             freedom(11) = -deriv(9,i1)          
             freedom(12) = deriv(13,i1)          
          endif
          if(i2/=0)then
             freedom(2) = deriv(1,i2)
             freedom(6) = deriv(4,i2)
             freedom(7) = -deriv(5,i2)
          endif
          if(i3/=0)then
             freedom(3) = deriv(1,i3)
          endif



       else

          if(i1/=0)then
             freedom(1) = deriv(1,i1)
             freedom(4) = deriv(6,i1)
             freedom(5) = -deriv(7,i1)
             freedom(8) = -deriv(5,i1)
             freedom(9) = deriv(4,i1)
             freedom(10) = deriv(10,i1)
             freedom(11) = deriv(12,i1)
             freedom(12) = -deriv(11,i1)
          endif
          if(i2/=0)then
             freedom(2) = deriv(1,i2)
          endif
          if(i3/=0)then
             freedom(3) = deriv(1,i3)
             freedom(6) = -deriv(3,i3)
             freedom(7) = deriv(2,i3)
          endif


       endif

       !freedom(10:12) = 0._f64

       lam(1) = 1._f64-lam(2)-lam(3)

       base(1) =  lam(1)**2*(3._f64-2._f64*lam(1)+6._f64*lam(2)*lam(3))      
       base(2) =  lam(2)**2*(3._f64-2._f64*lam(2)+6._f64*lam(3)*lam(1))      
       base(3) =  lam(3)**2*(3._f64-2._f64*lam(3)+6._f64*lam(1)*lam(2))      

       base(4) = lam(1)**2*lam(2)*(1._f64+2._f64*lam(3))
       base(5) = lam(2)**2*lam(1)*(1._f64+2._f64*lam(3))

       base(6) = lam(2)**2*lam(3)*(1._f64+2._f64*lam(1))
       base(7) = lam(3)**2*lam(2)*(1._f64+2._f64*lam(1))

       base(8) = lam(3)**2*lam(1)*(1._f64+2._f64*lam(2))
       base(9) = lam(1)**2*lam(3)*(1._f64+2._f64*lam(2))

       base(10) = lam(1)**2*lam(2)*lam(3)
       base(11) = lam(2)**2*lam(3)*lam(1)
       base(12) = lam(3)**2*lam(1)*lam(2)

       !L[0]=p**2*(3-2*p+6*q*r) f(0,0) p+n -> deriv(1,i1)
       !L[1]=q**2*(3-2*q+6*r*p) f(1,0) p; f(1,1) n -> deriv(1,i2)
       !L[2]=r**2*(3-2*r+6*p*q) f(1,1) p; f(0,1) n -> deriv(i,i3)
       !L[3]=p**2*q*(1+2*r)     fr1(0,0) p; fr3(0,0) n 
       !L[4]=q**2*p*(1+2*r)     fr1(1,0) p; fr3(1,1) n
       !L[5]=p**2*q**2/24      0
       !L[6]=q**2*r*(1+2*p)  
       !L[7]=r**2*q*(1+2*p)
       !L[8]=q**2*r**2/24)      0
       !L[9]=r**2*p*(1+2*q)
       !L[10]=p**2*r*(1+2*q)
       !L[11]=r**2*p**2/24      0
       !L[12]=p**2*q*r          frrp(0,0) p; frrn(0,0) n
       !L[13]=q**2*r*p          frrp(1,0) p; frrn(1,1) n
       !L[14]=r**2*p*q          frrp(1,1) p; frrn(0,1) n



       f = 0._f64
       do s = 1,12
          f = f + freedom(s)*base(s)
       enddo
       rho_tn1(i) = f    
       !L[0]=p**2*(3-2*p+6*q*r) f(0,0) p+n -> deriv(1,i1)
       !L[1]=q**2*(3-2*q+6*r*p) f(1,0) p; f(1,1) n -> deriv(1,i2)
       !L[2]=r**2*(3-2*r+6*p*q) f(1,1) p; f(0,1) n -> deriv(i,i3)
       !L[3]-p**2*q*(1+2*r)     fr1(0,0) p; fr3(0,0) n 
       !L[4]=q**2*p*(1+2*r)     fr1(1,0) p; fr3(1,1) n
       !L[5]-=p**2*q**2/24      0
       !L[6]=q**2*r*(1+2*p)  
       !L[7]=r**2*q*(1+2*p)
       !L[8]=q**2*r**2/24)      0
       !L[9]=r**2*p*(1+2*q)
       !L[10]=p**2*r*(1+2*q)
       !L[11]=r**2*p**2/24      0
       !L[12]=p**2*q*r          frrp(0,0) p; frrn(0,0) n
       !L[13]=q**2*r*p          frrp(1,0) p; frrn(1,1) n
       !L[14]=r**2*p*q          frrp(1,1) p; frrn(0,1) n


       !check solution in P1 case
       !      err_loc = abs(rho_tn1(i) &
       !        -freedom(1)*lam(1)-freedom(2)*lam(2)-freedom(3)*lam(3))
       !      err = max(err,err_loc)
       !      
       !      !if((maxval(abs(freedom(1:3)))<1.e-10) .and. (err_loc>1e-4))then
       !      if(err_loc>1e-5)then
       !        print *,'#i,err_loc=',i,err_loc
       !        print *,'#lam=',lam
       !        print *,'#freedom=',freedom
       !        print *,'#i1,i2,i3=',i1,i2,i3
       !        print *,'#deriv(:,i1)=',deriv(:,i1)
       !        print *,'#deriv(:,i2)=',deriv(:,i2)
       !        print *,'#deriv(:,i3)=',deriv(:,i3)
       !        print *,'#xx=',xx
       !        stop
       !      endif





    enddo

    !print *,'#err with p1=',err


  end subroutine interpolate_mitchell


  subroutine interpolate_mitchell_new( &
       radius, &
       num_cells, &
       index1, &
       deriv, &
       rho_tn1, &
       positions, &
       r )
    sll_real64, intent(in) :: radius
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: index1(:,:)
    sll_real64, intent(in) :: deriv(:,:)
    sll_real64, intent(out) :: rho_tn1(:)
    sll_real64, intent(in) :: positions(:,:)
    sll_real64, intent(in), optional :: r(2,2)
    !sll_real64, intent(in) :: dt

    sll_int32 :: i
    sll_real64 :: xx
    sll_real64 :: yy
    sll_real64 :: det
    sll_real64 :: r11
    sll_real64 :: r12
    sll_real64 :: r21
    sll_real64 :: r22
    sll_real64 :: h1
    sll_real64 :: h2
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i3
    sll_int32 :: ii
    sll_int32 :: jj
    sll_real64 :: a2
    sll_int32 :: s
    sll_real64 :: freedom(12)
    sll_real64 :: base(12)
    sll_real64 :: lam(3)
    sll_real64 :: f
    sll_real64 :: r1_x1
    sll_real64 :: r1_x2
    sll_real64 :: r2_x1
    sll_real64 :: r2_x2
    sll_real64 :: delta
    sll_int32 :: num_pts_tot
    sll_real64 :: x2p
    sll_real64 :: x3p
    sll_real64 :: x2n
    sll_real64 :: x3n
    sll_real64 :: y2p
    sll_real64 :: y3p
    sll_real64 :: y2n
    sll_real64 :: y3n

    !sll_real64 :: err_loc
    sll_real64 :: err

    err = 0._f64

    delta = radius/real(num_cells,f64)
    r1_x1 = sqrt(3._f64) * 0.5_f64*delta
    r1_x2 = 0.5_f64*delta
    r2_x1 = -r1_x1
    r2_x2 = r1_x2

    if(present(r))then
       if(abs(r1_x1-r(1,1))>1.e-14)then
          print *,'#value for r1_x1 not supported'
          stop
       endif
       if(abs(r2_x2-r(2,2))>1.e-14)then
          print *,'#value for r2_x2 not supported'
          stop
       endif
       if(abs(r1_x2-r(1,2))>1.e-14)then
          print *,'#value for r1_x2 not supported'
          stop
       endif
       if(abs(r2_x1-r(2,1))>1.e-14)then
          print *,'#value for r2_x1 not supported'
          stop
       endif
    endif
    num_pts_tot = 3*num_cells*(num_cells+1)+1
    det = r1_x1*r2_x2-r1_x2*r2_x1
    r11 = r2_x2/det
    r12 = -r2_x1/det
    r21 = -r1_x2/det
    r22 = r1_x1/det

    a2  = 2._f64/(delta**2*sqrt(3._f64))
    x2p = a2*r1_x1
    y2p = -a2*r1_x2
    x3p = -a2*(r1_x1+r2_x1)
    y3p = a2*(r1_x2+r2_x2)

    x3n = -a2*r2_x1
    y3n = a2*r2_x2
    x2n = a2*(r1_x1+r2_x1)
    y2n = -a2*(r1_x2+r2_x2)


    do i=1,num_pts_tot
       xx = positions(1,i)
       yy = positions(2,i)
       h1 =  xx*r11 + yy*r12
       h2 =  xx*r21 + yy*r22
       ii = floor(h1)
       jj = floor(h2)        
       h1 = h1-real(ii,f64)
       h2 = h2-real(jj,f64)
       if(ii>=num_cells)then
          ii=num_cells-1
          h1 = 1._f64
       endif
       if(ii<=-num_cells)then
          ii=-num_cells
          h1 = 0._f64
       endif
       if(jj>=num_cells)then
          jj=num_cells-1
          h2 = 1._f64
       endif
       if(jj<=-num_cells)then
          jj=-num_cells
          h2 = 0._f64
       endif

       xx = r1_x1*h1+r2_x1*h2
       yy = r1_x2*h1+r2_x2*h2
       i1 = index1(num_cells+1+ii,num_cells+1+jj)
       if ( xx >= 0._f64 ) then
          i2 = index1(num_cells+1+ii+1,num_cells+1+jj)          
          i3 = index1(num_cells+1+ii+1,num_cells+1+jj+1)
          lam(2) = y3p*xx+x3p*yy
          lam(3) = y2p*xx+x2p*yy
       else
          i2 = index1(num_cells+1+ii+1,num_cells+1+jj+1)
          i3 = index1(num_cells+1+ii,num_cells+1+jj+1)
          lam(2) = y3n*xx+x3n*yy
          lam(3) = y2n*xx+x2n*yy
       endif
       freedom = 0._f64

       if(xx>=0.)then

          if(i1/=0)then
             freedom(1) = deriv(1,i1)
             freedom(4) = deriv(2,i1)
             freedom(9) = deriv(3,i1)
             freedom(10) = deriv(8,i1)
          endif
          if(i2/=0)then
             freedom(2) = deriv(1,i2)
             freedom(5) = deriv(5,i2)
             freedom(6) = deriv(4,i2)
             freedom(11) = deriv(10,i2)
          endif
          if(i3/=0)then
             freedom(3) = deriv(1,i3)
             freedom(7) = deriv(7,i3)
             freedom(8) = deriv(6,i3)
             freedom(12) = deriv(12,i3)
          endif



       else

          if(i1/=0)then
             freedom(1) = deriv(1,i1)
             freedom(4) = deriv(3,i1)
             freedom(9) = deriv(4,i1)
             freedom(10) = deriv(9,i1)
          endif
          if(i2/=0)then
             freedom(2) = deriv(1,i2)
             freedom(5) = deriv(6,i2)
             freedom(6) = deriv(5,i2)
             freedom(11) = deriv(11,i2)
          endif
          if(i3/=0)then
             freedom(3) = deriv(1,i3)
             freedom(7) = deriv(2,i3)
             freedom(8) = deriv(7,i3)
             freedom(12) = deriv(13,i3)
          endif


       endif

       !freedom(10:12) = 0._f64

       lam(1) = 1._f64-lam(2)-lam(3)

       base(1) =  lam(1)**2*(3._f64-2._f64*lam(1)+6._f64*lam(2)*lam(3))      
       base(2) =  lam(2)**2*(3._f64-2._f64*lam(2)+6._f64*lam(3)*lam(1))      
       base(3) =  lam(3)**2*(3._f64-2._f64*lam(3)+6._f64*lam(1)*lam(2))      

       base(4) = lam(1)**2*lam(2)*(1._f64+2._f64*lam(3))
       base(5) = lam(2)**2*lam(1)*(1._f64+2._f64*lam(3))

       base(6) = lam(2)**2*lam(3)*(1._f64+2._f64*lam(1))
       base(7) = lam(3)**2*lam(2)*(1._f64+2._f64*lam(1))

       base(8) = lam(3)**2*lam(1)*(1._f64+2._f64*lam(2))
       base(9) = lam(1)**2*lam(3)*(1._f64+2._f64*lam(2))

       base(10) = lam(1)**2*lam(2)*lam(3)
       base(11) = lam(2)**2*lam(3)*lam(1)
       base(12) = lam(3)**2*lam(1)*lam(2)

       !L[0]=p**2*(3-2*p+6*q*r) f(0,0) p+n -> deriv(1,i1)
       !L[1]=q**2*(3-2*q+6*r*p) f(1,0) p; f(1,1) n -> deriv(1,i2)
       !L[2]=r**2*(3-2*r+6*p*q) f(1,1) p; f(0,1) n -> deriv(i,i3)
       !L[3]=p**2*q*(1+2*r)     fr1(0,0) p; fr3(0,0) n 
       !L[4]=q**2*p*(1+2*r)     fr1(1,0) p; fr3(1,1) n
       !L[5]=p**2*q**2/24      0
       !L[6]=q**2*r*(1+2*p)  
       !L[7]=r**2*q*(1+2*p)
       !L[8]=q**2*r**2/24)      0
       !L[9]=r**2*p*(1+2*q)
       !L[10]=p**2*r*(1+2*q)
       !L[11]=r**2*p**2/24      0
       !L[12]=p**2*q*r          frrp(0,0) p; frrn(0,0) n
       !L[13]=q**2*r*p          frrp(1,0) p; frrn(1,1) n
       !L[14]=r**2*p*q          frrp(1,1) p; frrn(0,1) n



       f = 0._f64
       do s = 1,12
          f = f + freedom(s)*base(s)
       enddo
       rho_tn1(i) = f    
       !L[0]=p**2*(3-2*p+6*q*r) f(0,0) p+n -> deriv(1,i1)
       !L[1]=q**2*(3-2*q+6*r*p) f(1,0) p; f(1,1) n -> deriv(1,i2)
       !L[2]=r**2*(3-2*r+6*p*q) f(1,1) p; f(0,1) n -> deriv(i,i3)
       !L[3]-p**2*q*(1+2*r)     fr1(0,0) p; fr3(0,0) n 
       !L[4]=q**2*p*(1+2*r)     fr1(1,0) p; fr3(1,1) n
       !L[5]-=p**2*q**2/24      0
       !L[6]=q**2*r*(1+2*p)  
       !L[7]=r**2*q*(1+2*p)
       !L[8]=q**2*r**2/24)      0
       !L[9]=r**2*p*(1+2*q)
       !L[10]=p**2*r*(1+2*q)
       !L[11]=r**2*p**2/24      0
       !L[12]=p**2*q*r          frrp(0,0) p; frrn(0,0) n
       !L[13]=q**2*r*p          frrp(1,0) p; frrn(1,1) n
       !L[14]=r**2*p*q          frrp(1,1) p; frrn(0,1) n


       !check solution in P1 case
       !      err_loc = abs(rho_tn1(i) &
       !        -freedom(1)*lam(1)-freedom(2)*lam(2)-freedom(3)*lam(3))
       !      err = max(err,err_loc)
       !      
       !      !if((maxval(abs(freedom(1:3)))<1.e-10) .and. (err_loc>1e-4))then
       !      if(err_loc>1e-5)then
       !        print *,'#i,err_loc=',i,err_loc
       !        print *,'#lam=',lam
       !        print *,'#freedom=',freedom
       !        print *,'#i1,i2,i3=',i1,i2,i3
       !        print *,'#deriv(:,i1)=',deriv(:,i1)
       !        print *,'#deriv(:,i2)=',deriv(:,i2)
       !        print *,'#deriv(:,i3)=',deriv(:,i3)
       !        print *,'#xx=',xx
       !        stop
       !      endif
       !




    enddo

    !print *,'#err with p1=',err


  end subroutine interpolate_mitchell_new




  subroutine interpolate_z9_new( &
       radius, &
       num_cells, &
       index1, &
       deriv, &
       rho_tn1, &
       positions, &
       r )
    sll_real64, intent(in) :: radius
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: index1(:,:)
    sll_real64, intent(in) :: deriv(:,:)
    sll_real64, intent(out) :: rho_tn1(:)
    sll_real64, intent(in) :: positions(:,:)
    sll_real64, intent(in), optional :: r(2,2)
    !sll_real64, intent(in) :: dt

    sll_int32 :: i
    sll_real64 :: xx
    sll_real64 :: yy
    sll_real64 :: det
    sll_real64 :: r11
    sll_real64 :: r12
    sll_real64 :: r21
    sll_real64 :: r22
    sll_real64 :: h1
    sll_real64 :: h2
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i3
    sll_int32 :: ii
    sll_int32 :: jj
    sll_real64 :: a2
    sll_int32 :: s
    sll_real64 :: freedom(9)
    sll_real64 :: base(9)
    sll_real64 :: lam(3)
    sll_real64 :: f
    sll_real64 :: r1_x1
    sll_real64 :: r1_x2
    sll_real64 :: r2_x1
    sll_real64 :: r2_x2
    sll_real64 :: delta
    sll_int32 :: num_pts_tot
    sll_real64 :: x2p
    sll_real64 :: x3p
    sll_real64 :: x2n
    sll_real64 :: x3n
    sll_real64 :: y2p
    sll_real64 :: y3p
    sll_real64 :: y2n
    sll_real64 :: y3n

    sll_real64     :: phi, p05
    sll_real64     :: l12, l22, l32, ksi1, ksi2, ksi3 
    sll_real64     :: ksi12, ksi13, ksi21, ksi23, ksi31, ksi32 

    !sll_real64 :: err
    !sll_real64 :: err_loc



    delta = radius/real(num_cells,f64)
    r1_x1 = sqrt(3._f64) * 0.5_f64*delta
    r1_x2 = 0.5_f64*delta
    r2_x1 = -r1_x1
    r2_x2 = r1_x2

    if(present(r))then
       if(abs(r1_x1-r(1,1))>1.e-14)then
          print *,'#value for r1_x1 not supported'
          stop
       endif
       if(abs(r2_x2-r(2,2))>1.e-14)then
          print *,'#value for r2_x2 not supported'
          stop
       endif
       if(abs(r1_x2-r(1,2))>1.e-14)then
          print *,'#value for r1_x2 not supported'
          stop
       endif
       if(abs(r2_x1-r(2,1))>1.e-14)then
          print *,'#value for r2_x1 not supported'
          stop
       endif
    endif
    num_pts_tot = 3*num_cells*(num_cells+1)+1
    det = r1_x1*r2_x2-r1_x2*r2_x1
    r11 = r2_x2/det
    r12 = -r2_x1/det
    r21 = -r1_x2/det
    r22 = r1_x1/det

    a2  = 2._f64/(delta**2*sqrt(3._f64))
    x2p = a2*r1_x1
    y2p = -a2*r1_x2
    x3p = -a2*(r1_x1+r2_x1)
    y3p = a2*(r1_x2+r2_x2)

    x3n = -a2*r2_x1
    y3n = a2*r2_x2
    x2n = a2*(r1_x1+r2_x1)
    y2n = -a2*(r1_x2+r2_x2)

    !err = 0._f64


    do i=1,num_pts_tot
       xx = positions(1,i)
       yy = positions(2,i)
       h1 =  xx*r11 + yy*r12
       h2 =  xx*r21 + yy*r22
       ii = floor(h1)
       jj = floor(h2)        
       h1 = h1-real(ii,f64)
       h2 = h2-real(jj,f64)
       if(ii>=num_cells)then
          ii=num_cells-1
          h1 = 1._f64
       endif
       if(ii<=-num_cells)then
          ii=-num_cells
          h1 = 0._f64
       endif
       if(jj>=num_cells)then
          jj=num_cells-1
          h2 = 1._f64
       endif
       if(jj<=-num_cells)then
          jj=-num_cells
          h2 = 0._f64
       endif

       xx = r1_x1*h1+r2_x1*h2
       yy = r1_x2*h1+r2_x2*h2
       i1 = index1(num_cells+1+ii,num_cells+1+jj)
       if ( xx >= 0._f64 ) then
          i2 = index1(num_cells+1+ii+1,num_cells+1+jj)          
          i3 = index1(num_cells+1+ii+1,num_cells+1+jj+1)
          lam(2) = y3p*xx+x3p*yy
          lam(3) = y2p*xx+x2p*yy
       else
          i2 = index1(num_cells+1+ii+1,num_cells+1+jj+1)
          i3 = index1(num_cells+1+ii,num_cells+1+jj+1)
          lam(2) = y3n*xx+x3n*yy
          lam(3) = y2n*xx+x2n*yy
       endif
       freedom = 0._f64
       if(i1/=0)then
          freedom(1) = deriv(1,i1)
       else
          freedom(1) = 0._f64
       endif
       if(i2/=0)then
          freedom(2) = deriv(1,i2)
       else
          freedom(2) = 0._f64
       endif
       if(i3/=0)then
          freedom(3) = deriv(1,i3)
       else
          freedom(3) = 0._f64
       endif

       !      if(xx>0._f64)then
       !        if(i1/=0)then
       !          freedom(4) = deriv(2,i1)
       !          freedom(5) = deriv(3,i1)
       !          freedom(8) = deriv(6,i1) 
       !          freedom(9) = deriv(7,i1) 
       !        endif
       !        if(i2/=0)then  
       !          freedom(6) = deriv(4,i2) 
       !          freedom(7) = deriv(5,i2)
       !        endif  
       !      else
       !        if(i1/=0)then
       !          freedom(4) = deriv(6,i1)
       !          freedom(5) = deriv(7,i1)
       !          freedom(8) = deriv(4,i1) 
       !          freedom(9) = deriv(5,i1)       
       !        endif
       !        if(i3/=0)then
       !          freedom(6) = deriv(3,i3) 
       !          freedom(7) = deriv(2,i3)
       !        endif  
       !      endif      


       if(xx>=0._f64)then
          if(i1/=0)then
             freedom(4) = deriv(2,i1)
             freedom(5) = deriv(6,i1)
             freedom(7) = -deriv(3,i1)
             freedom(8) = -deriv(7,i1)            
          endif
          if(i2/=0)then  
             freedom(6) = deriv(4,i2)
             freedom(9) = -deriv(5,i2) 
          endif
       else
          !        if(i1/=0)then
          !          freedom(4) = deriv(4,i1)
          !          freedom(5) = deriv(6,i1)
          !          freedom(7) = -deriv(5,i1)
          !          freedom(8) = -deriv(7,i1) 
          !        endif
          !        if(i3/=0)then
          !          freedom(6) = -deriv(3,i3) 
          !          freedom(9) = deriv(4,i3)       
          !        endif  
          if(i1/=0)then
             freedom(4) = deriv(6,i1)
             freedom(5) = deriv(4,i1)
             freedom(7) = -deriv(7,i1)
             freedom(8) = -deriv(5,i1)
          endif
          if(i3/=0)then
             freedom(6) = -deriv(3,i3)
             freedom(9) = deriv(2,i3) 
          endif
       endif


       !      if(maxval(abs(freedom))>1e-1.and.i1/=0 .and. i2/=0 .and. i3/=0)then
       !      print *,xx,i1,i2,i3
       !      do s=1,9
       !        print *,s,freedom(s)
       !      enddo
       !      print *,deriv(:,i1)
       !      print *,deriv(:,i2)
       !      print *,deriv(:,i3)
       !       stop
       !      endif

       lam(1) = 1._f64-lam(2)-lam(3)

       phi = lam(1)*lam(2)*lam(3)  ! "bubble function"

       ksi1 = lam(1)**3 - phi
       ksi2 = lam(2)**3 - phi
       ksi3 = lam(3)**3 - phi

       ! optimization variables

       p05 = phi*0.5_f64
       l12 = lam(1)**2
       l22 = lam(2)**2
       l32 = lam(3)**2

       ksi12 = l12*lam(2) + p05
       ksi13 = l12*lam(3) + p05
       ksi21 = l22*lam(1) + p05
       ksi23 = l22*lam(3) + p05
       ksi31 = l32*lam(1) + p05
       ksi32 = l32*lam(2) + p05

       base(1) = 3._f64 * l12 - 2._f64*ksi1 
       base(2) = 3._f64 * l22 - 2._f64*ksi2 
       base(3) = 3._f64 * l32 - 2._f64*ksi3 
       base(4) = ksi12 
       base(5) = ksi13
       !base(6) = ksi21
       !base(7) = ksi23
       base(6) = ksi23
       base(7) = ksi21
       base(8) = ksi31
       base(9) = ksi32      

       !      base(4) = ksi12 
       !      base(5) = ksi21
       !      base(6) = ksi23
       !      base(7) = ksi32
       !      base(8) = ksi13
       !      base(9) = ksi31      


       f = 0._f64
       do s = 1,9
          f = f + freedom(s)*base(s)
       enddo
       rho_tn1(i) = f

       !check solution in P1 case
       !      err_loc = abs(rho_tn1(i) &
       !        -freedom(1)*lam(1)-freedom(2)*lam(2)-freedom(3)*lam(3))
       !      err = max(err,err_loc)

       !if((maxval(abs(freedom(1:3)))<1.e-10) .and. (err_loc>1e-4))then
       !      if(err_loc>1e-4)then
       !        print *,'#i,err_loc=',i,err_loc
       !        print *,'#lam=',lam
       !        print *,'#freedom=',freedom
       !        print *,'#i1,i2,i3=',i1,i2,i3
       !        print *,'#deriv(:,i1)=',deriv(:,i1)
       !        print *,'#deriv(:,i2)=',deriv(:,i2)
       !        print *,'#deriv(:,i3)=',deriv(:,i3)
       !        print *,'#xx=',xx
       !        !stop
       !      endif


    enddo

    !print *,'#err with p1=',err
  end subroutine interpolate_z9_new




end program sim_bsl_ad_2d0v_hex_rot
