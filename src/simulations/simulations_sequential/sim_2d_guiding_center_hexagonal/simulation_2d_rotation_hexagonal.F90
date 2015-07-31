program rotation_2d_hexagonal

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_errors.h"
use sll_ascii_io

  use sll_constants
  use sll_hex_meshes
  use sll_box_splines
  use sll_interpolation_hex_hermite
  implicit none
  sll_int32, parameter :: SLL_HEX_SPLINES = 1 
  sll_int32, parameter :: SLL_HEX_Z9 = 9 
  sll_int32, parameter :: SLL_HEX_Z10 = 10 
  sll_int32, parameter :: SLL_HEX_HCTR = 11 
  sll_int32, parameter :: SLL_HEX_HCTC = 12 
  sll_int32, parameter :: SLL_HEX_GANEV_DIMITROV = 15 

  type(sll_hex_mesh_2d),   pointer        :: mesh
  type(sll_box_spline_2d), pointer        :: spline
  sll_real64, dimension(:,:), allocatable :: deriv
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
       ','can not open file')
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
    case ("SLL_HEX_Z10")
      num_method = SLL_HEX_Z10  
    case ("SLL_HEX_HCTR")
      num_method = SLL_HEX_HCTR  
    case ("SLL_HEX_HCTC")
      num_method = SLL_HEX_HCTC  
    case ("SLL_HEX_GANEV_DIMITROV")
      num_method = SLL_HEX_GANEV_DIMITROV  
    case default    
      SLL_ERROR("rotation_2d_hexagonal_hermite&
      ", "bad value of num_method_case")  
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
  n_points = mesh%num_pts_tot
  det = (mesh%r1_x1*mesh%r2_x2-mesh%r1_x2*mesh%r2_x1)/mesh%delta
  r11 = + mesh%r2_x2/det
  r12 = - mesh%r2_x1/det
  r21 = - mesh%r1_x2/det
  r22 = + mesh%r1_x1/det
  aire = mesh%delta**2*sqrt(3._f64)*0.25_f64
  
  
  if(num_method==SLL_HEX_SPLINES)then
    spline => new_box_spline_2d(mesh, SLL_DIRICHLET)
  else
    SLL_ALLOCATE(deriv(6,n_points),ierr)  
  endif
  SLL_ALLOCATE(rho_tn(n_points),ierr)
  SLL_ALLOCATE(rho_tn1(n_points),ierr)
  SLL_ALLOCATE(rho_exact(n_points),ierr)
  
  
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
  ! Time loop
  ! ----------------------------
  
  
  do step=1,number_iterations

    ! ----------------------------
    ! Compute interpolants
    ! ----------------------------

    
    if(num_method==SLL_HEX_SPLINES)then
      call compute_coeff_box_spline_2d( rho_tn, spline_degree, spline )
    else
      call  der_finite_difference( rho_tn, p, mesh%delta, mesh, deriv)
    endif

    ! ----------------------------
    ! Interpolation at vertices
    ! ----------------------------

    do i=1, n_points
      x = mesh%cartesian_coord(1,i)
      y = mesh%cartesian_coord(2,i)
      xx = x*cos(dt) - y*sin(dt)
      yy = x*sin(dt) + y*cos(dt)

      inside = .true.
      h1 =  xx*r11 + yy*r12
      h2 =  xx*r21 + yy*r22

      if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
      if ( abs(xx) > (radius-mesh%delta)*sll_sqrt3*0.5_f64) inside = .false.

      if ( inside ) then
        if(num_method==SLL_HEX_SPLINES)then
          rho_tn1(i) = hex_interpolate_value(mesh, xx, yy, spline, spline_degree)
        else
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
		endif    
      else
        rho_tn1(i) = 0._f64 ! dirichlet boundary condition
      endif
    enddo


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
      call write_field_hex_mesh_xmf(mesh, rho_tn, trim(filename))                
      filename  = trim(rho_error_name)//trim(filenum)
      call write_field_hex_mesh_xmf(mesh, rho_exact-rho_tn, trim(filename))                
    endif

      
  enddo
  
  
  
  close(thdiag_1d_id)
  
  
  call cpu_time(t_end)
  
  print *, &
    '#cpu time,min of rho_min', &
    t_end-t_init, &
    rho_min    
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
    linf_err_edge

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


end program rotation_2d_hexagonal