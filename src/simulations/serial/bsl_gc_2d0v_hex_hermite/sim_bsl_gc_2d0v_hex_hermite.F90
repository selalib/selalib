program sim_bsl_gc_2d0v_hex_hermite
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_m_constants
  use sll_m_euler_2d_hex
  use sll_m_hexagonal_meshes, only : &
       sll_hex_mesh_2d
  use sll_m_hex_poisson
  use sll_m_pivotbande
  use sll_m_interpolation_hex_hermite
  use sll_m_utilities, only: int2string, sll_new_file_id

  implicit none

  type(sll_hex_mesh_2d),      pointer     :: mesh
  type(sll_hex_mesh_2d),      pointer     :: mesh2
  sll_real64, dimension(:,:), allocatable :: deriv
  sll_real64, dimension(:),   allocatable :: dxuxn, dyuxn
  sll_real64, dimension(:),   allocatable :: dxuyn, dyuyn
  sll_real64, dimension(:),   allocatable :: second_term
  sll_real64, dimension(:),   allocatable :: uxn, uxn_1
  sll_real64, dimension(:),   allocatable :: uyn, uyn_1
  sll_real64, dimension(:),   allocatable :: phi
  sll_real64, dimension(:),   allocatable :: phi_interm
  sll_real64, dimension(:,:), allocatable :: matrix_poisson, l, u
  sll_real64, dimension(:),   allocatable :: rho_tn_1 ! distribution at time n-1
  sll_real64, dimension(:),   allocatable :: rho_tn   ! distribution at time n
  sll_real64, dimension(:),   allocatable :: rho_tn1  ! distribution at time n+1
  sll_real64, dimension(:),   allocatable :: dxuxn_1,dyuxn_1,dxuyn_1,dyuyn_1
  sll_real64, dimension(:),   allocatable :: rho_tn_2,uxn_2,uyn_2
  sll_real64, dimension(:),   allocatable :: dxuxn_2,dyuxn_2,dxuyn_2,dyuyn_2
  ! centers' variables
  sll_real64, dimension(:),   allocatable :: rho_center_tn, rho_center_tn1
  sll_real64, dimension(:),   allocatable :: uxn_center,uyn_center, phi_center_interm
  sll_real64, dimension(:),   allocatable :: dxuxn_center,dyuxn_center,dxuyn_center,dyuyn_center
  sll_real64, dimension(:),   allocatable :: second_term_center, phi_center
  sll_real64, dimension(:,:), allocatable :: matrix_poisson_center, l_center, u_center
  sll_real64, dimension(:),   allocatable :: rho_center_tn_1,uxn_1_center,uyn_1_center
  sll_real64, dimension(:),   allocatable :: dxuxn_1_center,dyuxn_1_center,dxuyn_1_center,dyuyn_1_center
  sll_real64, dimension(:),   allocatable :: rho_center_tn_2,uxn_2_center,uyn_2_center
  sll_real64, dimension(:),   allocatable :: dxuxn_2_center,dyuxn_2_center,dxuyn_2_center,dyuyn_2_center
  ! middle of the edge's variables
  sll_real64, dimension(:),   allocatable :: rho_edge_tn, rho_edge_tn1
  sll_real64, dimension(:),   allocatable :: dxuxn_edge,dyuxn_edge,dxuyn_edge,dyuyn_edge
  sll_real64, dimension(:),   allocatable :: uxn_edge,uyn_edge
  sll_real64, dimension(:),   allocatable :: rho_edge_tn_1,uxn_1_edge,uyn_1_edge
  sll_real64, dimension(:),   allocatable :: dxuxn_1_edge,dyuxn_1_edge,dxuyn_1_edge,dyuyn_1_edge
  sll_real64, dimension(:),   allocatable :: rho_edge_tn_2,uxn_2_edge,uyn_2_edge
  sll_real64, dimension(:),   allocatable :: dxuxn_2_edge,dyuxn_2_edge,dxuyn_2_edge,dyuyn_2_edge
  ! 2 times more dense mesh's variables
  sll_real64, dimension(:),   allocatable :: second_term2, phi2,uxn2,uyn2,phi2_interm,rho2
  sll_real64, dimension(:),   allocatable :: dxuxn2,dyuxn2,dxuyn2,dyuyn2
  sll_real64, dimension(:,:), allocatable :: matrix_poisson2, l2, u2
  sll_int32                               :: width_band1_2,width_band2_2
  sll_int32                               :: n_points2

  sll_int32    :: spline_degree
  sll_int32    :: hermite_method
  sll_int32    :: i, k1, k2, index_tab, type
  sll_int32    :: width_band1,width_band2
  sll_int32    :: i1,i2,i3
  sll_int32    :: num_cells, n_points, n_triangle, n_edge
  sll_int32    :: cells_min, cells_max
  sll_int32    :: cells_stp
  sll_int32    :: nloops,count, ierr, EXTRA_TABLES
  sll_real64   :: center_mesh_x1, center_mesh_x2, radius
  sll_real64   :: epsilon
  sll_real64   :: gauss_x1, gauss_x2
  sll_real64   :: gauss_sig
  sll_real64   :: gauss_amp
  sll_real64   :: dt
  sll_real64   :: aire
  sll_real64   :: step
  sll_real64   :: tmax
  sll_real64   :: t
  sll_real64   :: t_init, t_end
  sll_real64   :: t3
  sll_real64   :: h1, h2, x ,y,xx, yy
  sll_real64   :: r11,r12,r21,r22,det
  logical      :: inside
  sll_int32    :: p = 6 ! degree of the approximation for the derivative
  sll_int32            :: IO_stat
  sll_int32, parameter :: input_file = 99
  character(len = 256) :: input_filename
  character(len = 256) :: model_name
  character(len = 50)  :: filename
  character(len = 4)   :: filenum
  character(len = 4)   :: degnum


  namelist /geometry/ &
       center_mesh_x1, &
       center_mesh_x2, &
       radius, &
       cells_min, &
       cells_max, &
       cells_stp

  namelist /test_case/ &
       model_name

  namelist /initial_function/ &
       epsilon,   &
       gauss_x1,  &
       gauss_x2,  &
       gauss_sig, &
       gauss_amp

  namelist /time_iterations/ &
       dt, &
       tmax

  namelist /interpolation/ &
       spline_degree, &
       hermite_method

  ! ----------------------------
  ! Setting default parameters
  ! ----------------------------
  ! Test case :
  model_name = "GC"  
  ! Mesh :
  center_mesh_x1 = 0._f64
  center_mesh_x2 = 0._f64
  radius = 14._f64
  cells_min = 20
  cells_max = 160
  cells_stp = 20
  ! Initial function :
  epsilon = 0.001_f64
  gauss_x1  = 2._f64
  gauss_x2  = 2._f64
  gauss_sig = 1._f64/( 2._f64 * sqrt(2._f64)) 
  gauss_amp = 1.0_f64
  ! Time iterations:
  tmax  = 20._f64
  dt    = 0.1_f64
  ! Interpolation
  spline_degree = 2
  hermite_method = 9
  
  ! ----------------------------
  ! Reading from file
  ! ----------------------------
  call get_command_argument(1, input_filename)
  if(len_trim(input_filename).gt.0)then
     open(unit = input_file, file=trim(input_filename),IOStat=IO_stat)
     if( IO_stat /= 0 ) then
        print *, '#simulation_2d_guiding_center_hexagonal_splines() failed to open file ', &
             trim(input_filename)
        STOP
     end if
     print *,'#initialization with filename:'
     print *,'#',trim(input_filename)
     read(input_file, geometry)
     read(input_file, test_case)
     read(input_file, initial_function)
     read(input_file, time_iterations)
     read(input_file, interpolation)
     close(input_file)
  else
     print *,'#initialization with default parameters'
  endif

  ! Initializing the extra-tables variable:
  if (hermite_method.eq.15) then
     EXTRA_TABLES = 1
  end if

  do num_cells = cells_min,cells_max,cells_stp

     call int2string(hermite_method,degnum)
     if (model_name.eq."GC") then
        print*, ""
        print*," ********************************* "
        print*,"     Guiding-Center Simulation"
        print*,"        on a Hexagonal mesh"
        print*,"   using hermite elements method # =", degnum
        call print_method(hermite_method)
        print*," ********************************* "
     elseif (model_name.eq."CIRCULAR") then
        print*, ""
        print*," ********************************* "
        print*,"     Circular Advection Simulation"
        print*,"        on a Hexagonal mesh"
        print*,"   using hermite elements method # =", degnum
        call print_method(hermite_method)
        print*," ********************************* "
     else
        print *, "ERROR : No recognized model name;"
        print *, "        Options are: 'GC' or 'CIRCULAR'."
        print *, "        Please change it in your namelist file."
     end if

     t = 0._f64
     nloops = 0
     count  = 0
     if (model_name.eq."CIRCULAR") then
        dt = 0.1_f64 * 20._f64/num_cells
     end if

     !*********************************************************
     !             allocation
     !*********************************************************

     ! Mesh creation -------------------------
     mesh => new_hex_mesh_2d( num_cells, center_mesh_x1, center_mesh_x2,&
         radius=radius, EXTRA_TABLES = EXTRA_TABLES )
     call display_hex_mesh_2d(mesh)

     if  ( hermite_method == 15 ) then
        mesh2 => new_hex_mesh_2d( num_cells*2, center_mesh_x1, center_mesh_x2,&
             radius=radius, EXTRA_TABLES = EXTRA_TABLES ) 
     endif

     ! Inintialization of some other variables:
     n_points   = mesh%num_pts_tot
     n_triangle = mesh%num_triangles
     n_edge     = mesh%num_edges

     step = radius / real(num_cells,f64)
     aire = step**2*sqrt(3._f64)*0.25_f64

     det = (mesh%r1_x1*mesh%r2_x2 - mesh%r1_x2*mesh%r2_x1)/mesh%delta

     r11 = + mesh%r2_x2/det
     r12 = - mesh%r2_x1/det
     r21 = - mesh%r1_x2/det
     r22 = + mesh%r1_x1/det
     ! ---------------------------------------

     width_band1 = 2*num_cells+1
     width_band2 = width_band1

     allocate( deriv(1:6,n_points) )
     
     SLL_ALLOCATE(rho_tn_1( n_points),ierr)
     SLL_ALLOCATE(rho_tn( n_points),ierr)
     SLL_ALLOCATE(rho_tn1( n_points ),ierr)
     
     SLL_ALLOCATE(uxn( n_points),ierr)
     SLL_ALLOCATE(uyn( n_points ),ierr)
     
     SLL_ALLOCATE(dxuxn( n_points),ierr)
     SLL_ALLOCATE(dxuyn( n_points ),ierr)
     SLL_ALLOCATE(dyuxn( n_points),ierr)
     SLL_ALLOCATE(dyuyn( n_points ),ierr)

     SLL_ALLOCATE(uxn_1( n_points),ierr)
     SLL_ALLOCATE(uyn_1( n_points ),ierr)
     SLL_ALLOCATE(dxuxn_1( n_points),ierr)
     SLL_ALLOCATE(dxuyn_1( n_points ),ierr)
     SLL_ALLOCATE(dyuxn_1( n_points),ierr)
     SLL_ALLOCATE(dyuyn_1( n_points ),ierr)

     SLL_ALLOCATE(rho_tn_2( n_points),ierr)
     SLL_ALLOCATE(uxn_2( n_points),ierr)
     SLL_ALLOCATE(uyn_2( n_points ),ierr)
     SLL_ALLOCATE(dxuxn_2( n_points),ierr)
     SLL_ALLOCATE(dxuyn_2( n_points ),ierr)
     SLL_ALLOCATE(dyuxn_2( n_points),ierr)
     SLL_ALLOCATE(dyuyn_2( n_points ),ierr)

     
     if (model_name.eq."GC") then
        ! variables only used in the guiding center model
        SLL_ALLOCATE(phi( n_points),ierr)
        SLL_ALLOCATE(phi_interm( n_points),ierr)
        
        SLL_ALLOCATE(uxn_1( n_points),ierr)
        SLL_ALLOCATE(uyn_1( n_points ),ierr)
        SLL_ALLOCATE(second_term( n_points),ierr)

        SLL_ALLOCATE(matrix_poisson( n_points,1 + 4*num_cells + 2 ) , ierr)
        SLL_ALLOCATE(l( n_points,1 + 4*num_cells + 2 ) , ierr)
        SLL_ALLOCATE(u( n_points,1 + 4*num_cells + 2), ierr)
     end if

     if  ( hermite_method == 10 )  then
        SLL_ALLOCATE(rho_center_tn ( n_triangle),ierr)
        SLL_ALLOCATE(rho_center_tn1( n_triangle),ierr)
        SLL_ALLOCATE(uxn_center( n_triangle),ierr)
        SLL_ALLOCATE(uyn_center( n_triangle),ierr)
        SLL_ALLOCATE(dxuxn_center( n_triangle),ierr)
        SLL_ALLOCATE(dxuyn_center( n_triangle),ierr)
        SLL_ALLOCATE(dyuxn_center( n_triangle),ierr)
        SLL_ALLOCATE(dyuyn_center( n_triangle),ierr)

        SLL_ALLOCATE(rho_center_tn_1 ( n_triangle),ierr)
        SLL_ALLOCATE(uxn_1_center( n_triangle),ierr)
        SLL_ALLOCATE(uyn_1_center( n_triangle),ierr)
        SLL_ALLOCATE(dxuxn_1_center( n_triangle),ierr)
        SLL_ALLOCATE(dxuyn_1_center( n_triangle),ierr)
        SLL_ALLOCATE(dyuxn_1_center( n_triangle),ierr)
        SLL_ALLOCATE(dyuyn_1_center( n_triangle),ierr)

        SLL_ALLOCATE(rho_center_tn_2 ( n_triangle),ierr)
        SLL_ALLOCATE(uxn_2_center( n_triangle),ierr)
        SLL_ALLOCATE(uyn_2_center( n_triangle),ierr)
        SLL_ALLOCATE(dxuxn_2_center( n_triangle),ierr)
        SLL_ALLOCATE(dxuyn_2_center( n_triangle),ierr)
        SLL_ALLOCATE(dyuxn_2_center( n_triangle),ierr)
        SLL_ALLOCATE(dyuyn_2_center( n_triangle),ierr)

        SLL_ALLOCATE(phi_center_interm( n_points),ierr)          
        SLL_ALLOCATE(phi_center( n_triangle),ierr)  
        SLL_ALLOCATE(second_term_center( n_points),ierr)
        SLL_ALLOCATE(matrix_poisson_center( n_triangle,1 + 4*num_cells + 2), ierr)
        SLL_ALLOCATE(l_center( n_triangle,1 + 4*num_cells + 2 ) , ierr)
        SLL_ALLOCATE(u_center( n_triangle,1 + 4*num_cells + 2), ierr)
     endif

     if  ( hermite_method == 15 )  then

        width_band1_2 = 4*num_cells+1
        width_band2_2 = width_band1_2
        n_points2   = 1 + 6 * num_cells * (2*num_cells + 1) 

        SLL_ALLOCATE(rho_edge_tn ( n_edge),ierr)
        SLL_ALLOCATE(rho_edge_tn1( n_edge),ierr)
        SLL_ALLOCATE(uxn_edge( n_edge),ierr)
        SLL_ALLOCATE(uyn_edge( n_edge),ierr) 
        SLL_ALLOCATE(dxuxn_edge( n_edge),ierr)
        SLL_ALLOCATE(dxuyn_edge( n_edge),ierr)
        SLL_ALLOCATE(dyuxn_edge( n_edge),ierr)
        SLL_ALLOCATE(dyuyn_edge( n_edge),ierr)

        SLL_ALLOCATE(rho_edge_tn_1( n_edge),ierr)
        SLL_ALLOCATE(uxn_1_edge( n_edge),ierr)
        SLL_ALLOCATE(uyn_1_edge( n_edge),ierr) 
        SLL_ALLOCATE(dxuxn_1_edge( n_edge),ierr)
        SLL_ALLOCATE(dxuyn_1_edge( n_edge),ierr)
        SLL_ALLOCATE(dyuxn_1_edge( n_edge),ierr)
        SLL_ALLOCATE(dyuyn_1_edge( n_edge),ierr)

        SLL_ALLOCATE(rho_edge_tn_2( n_edge),ierr)
        SLL_ALLOCATE(uxn_2_edge( n_edge),ierr)
        SLL_ALLOCATE(uyn_2_edge( n_edge),ierr) 
        SLL_ALLOCATE(dxuxn_2_edge( n_edge),ierr)
        SLL_ALLOCATE(dxuyn_2_edge( n_edge),ierr)
        SLL_ALLOCATE(dyuxn_2_edge( n_edge),ierr)
        SLL_ALLOCATE(dyuyn_2_edge( n_edge),ierr)

        SLL_ALLOCATE(second_term2( n_points2 ),ierr) 
        SLL_ALLOCATE(rho2( n_points2),ierr)     
        SLL_ALLOCATE(phi2( n_points2),ierr)           
        SLL_ALLOCATE(phi2_interm( n_points2),ierr) 
        SLL_ALLOCATE(uxn2( n_points2),ierr)
        SLL_ALLOCATE(uyn2( n_points2),ierr) 
        SLL_ALLOCATE(dxuxn2( n_points2),ierr)
        SLL_ALLOCATE(dxuyn2( n_points2),ierr)
        SLL_ALLOCATE(dyuxn2( n_points2),ierr)
        SLL_ALLOCATE(dyuyn2( n_points2),ierr)
        SLL_ALLOCATE(l2( n_points2,1 + 8*num_cells + 2 ) , ierr)
        SLL_ALLOCATE(u2( n_points2,1 + 8*num_cells + 2), ierr)
        SLL_ALLOCATE(matrix_poisson2( n_points2,1 + 8*num_cells + 2 ) , ierr)

     endif
     
     call cpu_time(t_init)

     !*********************************************************
     !  Distribution function & density initialization
     !*********************************************************


     ! Initial distribution ------------------
     if (model_name.eq."GC") then
        call init_distr_gc(rho_tn,rho_center_tn,rho_edge_tn,hermite_method, mesh,epsilon)
     elseif (model_name.eq."CIRCULAR") then
        call init_distr_circ(rho_tn, rho_center_tn,rho_edge_tn,hermite_method,mesh,gauss_x1,gauss_x2,gauss_sig,gauss_amp)
     end if
     ! ---------------------------------------

     ! Poisson solver ------------------------
     if (model_name.eq."GC") then
        if ( hermite_method /= 15 ) then
           call hex_matrix_poisson( matrix_poisson, mesh,type=1)
           call factolub_bande(matrix_poisson,l,u,n_points,width_band1,width_band2)
           call hex_second_terme_poisson( second_term, mesh, rho_tn )
           call solvlub_bande(l,u,phi_interm,second_term,n_points,width_band1,width_band2)        
           do i = 1, mesh%num_pts_tot
              k1 = mesh%hex_coord(1, i)
              k2 = mesh%hex_coord(2, i)
              call mesh%index_hex_to_global(k1, k2, index_tab)
              phi(i) = phi_interm(index_tab)
           enddo
           call compute_hex_fields(mesh,uxn,uyn,dxuxn,dyuxn,dxuyn,dyuyn,phi,type=1)
           if ( hermite_method == 10 ) then
              do i = 1,mesh%num_triangles
                 x = mesh%center_cartesian_coord(1,i)
                 y = mesh%center_cartesian_coord(2,i)
                 call get_cell_vertices_index( x, y, mesh, i1, i2, i3 )
                 uxn_center(i) = (uxn(i1)+uxn(i2)+uxn(i3))/3._f64
                 uyn_center(i) = (uyn(i1)+uyn(i2)+uyn(i3))/3._f64
              enddo
           end if
        end if
        if ( hermite_method == 15 ) then
           call hex_matrix_poisson( matrix_poisson2, mesh2,1)
           call factolub_bande(matrix_poisson2,l2,u2,n_points2,width_band1_2,width_band2_2)
           call assemble(rho_tn,rho_edge_tn,rho2,mesh,mesh2)
           call hex_second_terme_poisson( second_term2, mesh2, rho2)
           call solvlub_bande(l2,u2,phi2_interm,second_term2,n_points2,width_band1_2,width_band2_2)
           do i = 1, n_points2   ! need to re-index phi : 
              k1 = mesh2%hex_coord(1, i)
              k2 = mesh2%hex_coord(2, i)
              call mesh2%index_hex_to_global(k1, k2, index_tab)
              phi2(i) = phi2_interm(index_tab)
           enddo
           call compute_hex_fields(mesh2,uxn2,uyn2,dxuxn2,dyuxn2,dxuyn2,dyuyn2,phi2,type=1)
           call deassemble(uxn,uxn_edge,uxn2,mesh,mesh2)
           call deassemble(uyn,uyn_edge,uyn2,mesh,mesh2)
           call deassemble(dxuxn,dxuxn_edge,dxuxn2,mesh,mesh2)
           call deassemble(dxuyn,dxuyn_edge,dxuyn2,mesh,mesh2)
           call deassemble(dyuxn,dyuxn_edge,dyuxn2,mesh,mesh2)
           call deassemble(dyuyn,dyuyn_edge,dyuyn2,mesh,mesh2)
        endif
     end if
     ! ---------------------------------------

     if (model_name.eq."GC") then
        if ( hermite_method /= 15 ) then
           call hex_diagnostics_gc(rho_tn,t,mesh,uxn,uyn,nloops,hermite_method,tmax,cells_min,cells_max)
        else
           call hex_diagnostics_gc(rho2,t,mesh2,uxn2,uyn2,nloops,hermite_method,tmax,cells_min,cells_max)
        end if
     elseif (model_name.eq."CIRCULAR") then
        call hex_diagnostics_circ(rho_tn,t,mesh,nloops,hermite_method,tmax,cells_min,cells_max,gauss_x1, gauss_x2, gauss_sig, gauss_amp)
     end if


     !*********************************************************
     !                          Time loop
     !*********************************************************

     call cpu_time(t3)

     print*,"fin init",t3 - t_init

     do while (t .lt. tmax)

        t = t + dt
        nloops = nloops + 1
        count = count + 1


        !*********************************************************
        !                     interpolation
        !*********************************************************
        call  der_finite_difference( rho_tn, p, step, mesh, deriv )
        
        do i=1, n_points

           x = mesh%cartesian_coord(1,i)
           y = mesh%cartesian_coord(2,i)

           !*************************************************
           !       computation of the characteristics
           !*************************************************
           if (model_name.eq."CIRCULAR") then
              xx = x*cos(dt) - y*sin(dt)
              yy = x*sin(dt) + y*cos(dt)

           elseif (model_name.eq."GC") then
              ! We use Adams2 for solving ODE except for first step
              if ( t <= dt + 1e-6 ) then ! first step with euler
                 call compute_characteristic_euler_2d_hex( &
                      x,y,uxn,uyn,i,xx,yy,dt )
              else !the rest is done with Adams 2
                 call compute_characteristic_adams2_2d_hex( x,y,uxn,uyn,uxn_1,uyn_1,&
                      dxuxn,dyuxn,dxuyn,dyuyn,i,xx,yy,dt)
              endif
           end if
           
           inside = .true.
           h1 =  xx*r11 + yy*r12
           h2 =  xx*r21 + yy*r22

           if (abs(h1)>radius-mesh%delta .or. abs(h2)>radius-mesh%delta) inside = .false.
           if ( abs(xx) > (radius-mesh%delta)*sqrt(3._f64)*0.5_f64) inside = .false.

           if ( inside ) then
              call hermite_interpolation(i, xx, yy, rho_tn, rho_center_tn,&
                   rho_edge_tn, rho_tn1, mesh, deriv, aire,& 
                   hermite_method)
           else
              rho_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif

        end do ! end of the computation of the mesh points

        if ( hermite_method == 10 ) then 

           do i = 1,mesh%num_triangles

              x = mesh%center_cartesian_coord(1,i)
              y = mesh%center_cartesian_coord(2,i)
              if (t < 2*dt) then
                 call compute_characteristic_euler_2d_hex( &
                   x,y,uxn_center,uyn_center,i,xx,yy,dt )
              else
                 call compute_characteristic_adams2_2d_hex(  &
                      x,y,uxn_center,uyn_center,uxn_1_center,uyn_1_center, &
                      dxuxn_center,dyuxn_center,dxuyn_center,dyuyn_center, &
                      i,xx,yy,dt)
              endif

              inside = .true.
              h1 =  xx*r11 + yy*r12
              h2 =  xx*r21 + yy*r22 

              if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
              if ( abs(xx) > (radius-mesh%delta)*sqrt(3._f64)*0.5_f64) inside = .false.

              if ( inside ) then
                 call hermite_interpolation(i, x, y, rho_tn, rho_center_tn,&
                      rho_edge_tn, rho_center_tn1, mesh, deriv, aire,& 
                      hermite_method)
              else 
                 rho_center_tn1(i) = 0._f64 ! dirichlet boundary condition
              endif

           enddo! end of the computation of the center of the triangles
        endif

        if ( hermite_method == 15 ) then 

           ! besoin de la séparation ici pour avoir rho_tn et rho_edge_tn

           do i=1, n_edge

              x = mesh%edge_center_cartesian_coord(1,i)
              y = mesh%edge_center_cartesian_coord(2,i) 
              
              if (t < 2*dt) then
                 call compute_characteristic_euler_2d_hex( &
                      x,y,uxn_edge,uyn_edge,i,xx,yy,dt )
              else
                 call compute_characteristic_adams2_2d_hex( x,y,uxn_edge,uyn_edge,uxn_1_edge,uyn_1_edge,&
                      dxuxn_edge,dyuxn_edge,dxuyn_edge,dyuyn_edge,i,xx,yy,dt)
              endif

              inside = .true.
              h1 =  xx*r11 + yy*r12
              h2 =  xx*r21 + yy*r22 


              if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
              if ( abs(xx) > (radius-mesh%delta)*sqrt(3._f64)*0.5_f64) inside = .false.

              if ( inside ) then
                 call hermite_interpolation(i,xx,yy,rho_tn,rho_center_tn&
                      ,rho_edge_tn, rho_edge_tn1,mesh,deriv,aire,&
                      hermite_method)
              else 
                 rho_edge_tn1(i) = 0._f64 ! dirichlet boundary condition
              endif

           enddo ! end of the computation at the middles of the edges

        endif


        ! Updating the new field values ............
        rho_tn_2 = rho_tn_1
        rho_tn_1 = rho_tn
        rho_tn = rho_tn1

        if ( hermite_method == 10 ) then
           rho_center_tn_2 = rho_center_tn_1
           rho_center_tn_1 = rho_center_tn
           rho_center_tn   = rho_center_tn1
        endif
        if ( hermite_method == 15 ) then
           rho_edge_tn_2 = rho_edge_tn_1
           rho_edge_tn_1 = rho_edge_tn
           rho_edge_tn   = rho_edge_tn1
        endif
        ! ...........................................

        !*********************************************************
        !      computing the solution of the poisson equation 
        !*********************************************************
        if (model_name.eq."GC") then
           if ( hermite_method /= 15 ) then
              call hex_second_terme_poisson( second_term, mesh, rho_tn ) 
              call solvlub_bande(l,u,phi_interm,second_term,n_points,width_band1,width_band2)
           
              do i = 1, mesh%num_pts_tot    ! need to re-index phi :
                 k1 = mesh%hex_coord(1, i)
                 k2 = mesh%hex_coord(2, i)
                 call mesh%index_hex_to_global(k1, k2, index_tab)
                 phi(i) = phi_interm(index_tab)
              enddo

              call compute_hex_fields(mesh,uxn,uyn,dxuxn,dyuxn,dxuyn,dyuyn,phi,type=1)
              if ( hermite_method == 10 ) then
                 do i = 1,mesh%num_triangles
                    x = mesh%center_cartesian_coord(1,i)
                    y = mesh%center_cartesian_coord(2,i)
                    call get_cell_vertices_index( x, y, mesh, i1, i2, i3 )
                    uxn_center(i) = (uxn(i1)+uxn(i2)+uxn(i3))/3._f64
                    uyn_center(i) = (uyn(i1)+uyn(i2)+uyn(i3))/3._f64
                 enddo
              endif
           endif

           if ( hermite_method == 15 ) then
              call assemble(rho_tn,rho_edge_tn,rho2,mesh,mesh2)
              call hex_second_terme_poisson( second_term2, mesh2, rho2)
              call solvlub_bande(l2,u2,phi2_interm,second_term2,n_points2,width_band1_2,width_band2_2)
              
              do i = 1, n_points2  ! need to re-index phi : 
                 k1 = mesh2%hex_coord(1, i)
                 k2 = mesh2%hex_coord(2, i)
                 call mesh2%index_hex_to_global(k1, k2, index_tab)
                 phi2(i) = phi2_interm(index_tab)
              enddo
              call compute_hex_fields(mesh2,uxn2,uyn2,dxuxn2,dyuxn2,dxuyn2,dyuyn2,phi2,type=1)
              call deassemble(uxn,uxn_edge,uxn2,mesh,mesh2)
              call deassemble(uyn,uyn_edge,uyn2,mesh,mesh2)
              call deassemble(dxuxn,dxuxn_edge,dxuxn2,mesh,mesh2)
              call deassemble(dxuyn,dxuyn_edge,dxuyn2,mesh,mesh2)
              call deassemble(dyuxn,dyuxn_edge,dyuxn2,mesh,mesh2)
              call deassemble(dyuyn,dyuyn_edge,dyuyn2,mesh,mesh2)
           endif
    
           ! Updating the new field values ............
           uxn_2 = uxn_1
           uxn_1 = uxn
           uyn_2 = uyn_1
           uyn_1 = uyn
           
           dxuxn_2 = dxuxn_1
           dyuxn_2 = dyuxn_1
           dxuyn_2 = dxuyn_1
           dyuyn_2 = dyuyn_1
           
           dxuxn_1 = dxuxn
           dyuxn_1 = dyuxn
           dxuyn_1 = dxuyn
           dyuyn_1 = dyuyn

           if ( hermite_method == 10 ) then
              uxn_2_center = uxn_1_center
              uxn_1_center = uxn_center
              uyn_2_center = uyn_1_center
              uyn_1_center = uyn_center
              dxuxn_2_center = dxuxn_1_center
              dyuxn_2_center = dyuxn_1_center
              dxuyn_2_center = dxuyn_1_center
              dyuyn_2_center = dyuyn_1_center
              dxuxn_1_center = dxuxn_center
              dyuxn_1_center = dyuxn_center
              dxuyn_1_center = dxuyn_center
              dyuyn_1_center = dyuyn_center
           endif
           if ( hermite_method == 15 ) then
              uxn_2_edge = uxn_1_edge
              uxn_1_edge = uxn_edge
              uyn_2_edge = uyn_1_edge
              uyn_1_edge = uyn_edge
              dxuxn_2_edge = dxuxn_1_edge
              dyuxn_2_edge = dyuxn_1_edge
              dxuyn_2_edge = dxuyn_1_edge
              dyuyn_2_edge = dyuyn_1_edge
              dxuxn_1_edge = dxuxn_edge
              dyuxn_1_edge = dyuxn_edge
              dxuyn_1_edge = dxuyn_edge
              dyuyn_1_edge = dyuyn_edge
           endif
           ! ...........................................

        end if
        
        !*********************************************************
        !                  writing diagostics
        !*********************************************************
        if (model_name.eq."GC") then
           call hex_diagnostics_gc(rho_tn,t,mesh,uxn,uyn,nloops,hermite_method,tmax,cells_min,cells_max)
           if (count == 10.and.nloops<10000.and.num_cells == cells_max) then
              call int2string(nloops,filenum)
              filename  = "center_guide_rho"//trim(filenum)
              call mesh%write_field_hex_mesh_xmf(rho_tn1, trim(filename))
              if (model_name.eq."GC") then
                 filename  = "center_guide_phi"//trim(filenum)
                 call mesh%write_field_hex_mesh_xmf(phi, trim(filename))
              end if
              count = 0
           endif
        elseif (model_name.eq."CIRCULAR") then
           call hex_diagnostics_circ(rho_tn,t,mesh,nloops,hermite_method,tmax,cells_min,cells_max, gauss_x1, gauss_x2, gauss_sig, gauss_amp)
           if (count == 10.and.nloops<10000.and.num_cells == cells_max) then
              call int2string(nloops,filenum)
              filename  = "circular_advection_rho"//trim(filenum)
              call mesh%write_field_hex_mesh_xmf(rho_tn1, trim(filename))
              count = 0
           endif
        end if

     enddo ! end of time loop

     SLL_DEALLOCATE_ARRAY(rho_tn_1,ierr)
     SLL_DEALLOCATE_ARRAY(rho_tn,ierr)
     SLL_DEALLOCATE_ARRAY(rho_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(uxn,ierr)
     SLL_DEALLOCATE_ARRAY(uyn,ierr)
     SLL_DEALLOCATE_ARRAY(dxuxn,ierr)
     SLL_DEALLOCATE_ARRAY(dxuyn,ierr)
     SLL_DEALLOCATE_ARRAY(dyuxn,ierr)
     SLL_DEALLOCATE_ARRAY(dyuyn,ierr)

     SLL_DEALLOCATE_ARRAY(uxn_1,ierr)
     SLL_DEALLOCATE_ARRAY(uyn_1,ierr)
     SLL_DEALLOCATE_ARRAY(dxuxn_1,ierr)
     SLL_DEALLOCATE_ARRAY(dxuyn_1,ierr)
     SLL_DEALLOCATE_ARRAY(dyuxn_1,ierr)
     SLL_DEALLOCATE_ARRAY(dyuyn_1,ierr)

     SLL_DEALLOCATE_ARRAY(rho_tn_2,ierr)
     SLL_DEALLOCATE_ARRAY(uxn_2,ierr)
     SLL_DEALLOCATE_ARRAY(uyn_2,ierr)
     SLL_DEALLOCATE_ARRAY(dxuxn_2,ierr)
     SLL_DEALLOCATE_ARRAY(dxuyn_2,ierr)
     SLL_DEALLOCATE_ARRAY(dyuxn_2,ierr)
     SLL_DEALLOCATE_ARRAY(dyuyn_2,ierr)

     deallocate(deriv)
     
     if (model_name.eq."GC") then
        SLL_DEALLOCATE_ARRAY(uxn_1,ierr)
        SLL_DEALLOCATE_ARRAY(uyn_1,ierr)
        SLL_DEALLOCATE_ARRAY(second_term,ierr)
        SLL_DEALLOCATE_ARRAY(phi,ierr)
        SLL_DEALLOCATE_ARRAY(phi_interm,ierr)
        deallocate(matrix_poisson)
        deallocate(l)
        deallocate(u)
     end if

     if  ( hermite_method == 10 )  then
        SLL_DEALLOCATE_ARRAY(rho_center_tn,ierr)
        SLL_DEALLOCATE_ARRAY(rho_center_tn1,ierr)
        SLL_DEALLOCATE_ARRAY(uxn_center,ierr)
        SLL_DEALLOCATE_ARRAY(uyn_center,ierr)
        SLL_DEALLOCATE_ARRAY(dxuxn_center,ierr)
        SLL_DEALLOCATE_ARRAY(dxuyn_center,ierr)
        SLL_DEALLOCATE_ARRAY(dyuxn_center,ierr)
        SLL_DEALLOCATE_ARRAY(dyuyn_center,ierr)

        SLL_DEALLOCATE_ARRAY(rho_center_tn_1,ierr)
        SLL_DEALLOCATE_ARRAY(uxn_1_center,ierr)
        SLL_DEALLOCATE_ARRAY(uyn_1_center,ierr)
        SLL_DEALLOCATE_ARRAY(dxuxn_1_center,ierr)
        SLL_DEALLOCATE_ARRAY(dxuyn_1_center,ierr)
        SLL_DEALLOCATE_ARRAY(dyuxn_1_center,ierr)
        SLL_DEALLOCATE_ARRAY(dyuyn_1_center,ierr)

        SLL_DEALLOCATE_ARRAY(rho_center_tn_2,ierr)
        SLL_DEALLOCATE_ARRAY(uxn_2_center,ierr)
        SLL_DEALLOCATE_ARRAY(uyn_2_center,ierr)
        SLL_DEALLOCATE_ARRAY(dxuxn_2_center,ierr)
        SLL_DEALLOCATE_ARRAY(dxuyn_2_center,ierr)
        SLL_DEALLOCATE_ARRAY(dyuxn_2_center,ierr)
        SLL_DEALLOCATE_ARRAY(dyuyn_2_center,ierr)

        
        SLL_DEALLOCATE_ARRAY(phi_center_interm,ierr) 
        SLL_DEALLOCATE_ARRAY(second_term_center,ierr)    
        SLL_DEALLOCATE_ARRAY(phi_center,ierr) 
        deallocate(matrix_poisson_center,l_center,u_center)
     endif

     if  ( hermite_method == 15 )  then

        SLL_DEALLOCATE_ARRAY(rho_edge_tn,ierr)
        SLL_DEALLOCATE_ARRAY(rho_edge_tn1,ierr)
        SLL_DEALLOCATE_ARRAY(uxn_edge,ierr)
        SLL_DEALLOCATE_ARRAY(uyn_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dxuxn_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dxuyn_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dyuxn_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dyuyn_edge,ierr)

        SLL_DEALLOCATE_ARRAY(rho_edge_tn_1,ierr)
        SLL_DEALLOCATE_ARRAY(uxn_1_edge,ierr)
        SLL_DEALLOCATE_ARRAY(uyn_1_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dxuxn_1_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dxuyn_1_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dyuxn_1_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dyuyn_1_edge,ierr)

        SLL_DEALLOCATE_ARRAY(rho_edge_tn_2,ierr)
        SLL_DEALLOCATE_ARRAY(uxn_2_edge,ierr)
        SLL_DEALLOCATE_ARRAY(uyn_2_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dxuxn_2_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dxuyn_2_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dyuxn_2_edge,ierr)
        SLL_DEALLOCATE_ARRAY(dyuyn_2_edge,ierr)


        SLL_DEALLOCATE_ARRAY(second_term2,ierr)  
        SLL_DEALLOCATE_ARRAY(phi2,ierr)        
        SLL_DEALLOCATE_ARRAY(phi2_interm,ierr)  
        SLL_DEALLOCATE_ARRAY(rho2,ierr)  
        SLL_DEALLOCATE_ARRAY(uxn2,ierr)
        SLL_DEALLOCATE_ARRAY(uyn2,ierr) 
        SLL_DEALLOCATE_ARRAY(dxuxn2,ierr)
        SLL_DEALLOCATE_ARRAY(dxuyn2,ierr)
        SLL_DEALLOCATE_ARRAY(dyuxn2,ierr)
        SLL_DEALLOCATE_ARRAY(dyuyn2,ierr)
        deallocate(l2,u2,matrix_poisson2)

     endif


     call delete_hex_mesh_2d( mesh )

     call cpu_time(t_end)
     print*, "time used =", t_end - t_init

  end do ! end of loop on space step (end of all simulations)

  ! Some test should probably be put here:
  print*, 'PASSED'

contains

  !*********initialization**************

  subroutine init_distr_gc(f_tn,center_values_tn,edge_values_tn,hermite_method, mesh, epsilon)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:), intent(inout) :: f_tn
    sll_real64, dimension(:), intent(inout) :: center_values_tn
    sll_real64, dimension(:), intent(inout) :: edge_values_tn
    sll_real64, intent(in) :: epsilon
    sll_int32,  intent(in) :: hermite_method
    sll_real64 :: x, y
    sll_real64 :: r
    sll_int32  :: i

    do i = 1,mesh%num_pts_tot
       x = mesh%cartesian_coord(1,i)
       y = mesh%cartesian_coord(2,i)
       r = sqrt( x**2 + y**2 )

       if ( r <= 8._f64  .and. r >= 5._f64 ) then
          f_tn(i) = (1._f64 + epsilon * cos( 9._f64 * atan2(y,x)) )*&
               exp( -4._f64*(r-6.5_f64)**2)
       else
          f_tn(i) = 0._f64
       endif
    enddo
        if ( hermite_method == 10 ) then 

       do i = 1,mesh%num_triangles
          x = mesh%center_cartesian_coord(1,i)
          y = mesh%center_cartesian_coord(2,i)
          r = sqrt( x**2 + y**2 )

          if ( r <= 8._f64  .and. r >= 5._f64 ) then
             center_values_tn(i) = (1._f64 + epsilon * cos( 9._f64 * atan2(y,x)) )*&
                  exp( -4._f64*(r-6.5_f64)**2)
          else
             center_values_tn(i) = 0._f64
          endif

       enddo

    endif
    if ( hermite_method == 15 ) then 

       do i = 1,mesh%num_edges
          x = mesh%edge_center_cartesian_coord(1,i)
          y = mesh%edge_center_cartesian_coord(2,i)
          r = sqrt( x**2 + y**2 )

          if ( r <= 8._f64  .and. r >= 5._f64 ) then
             edge_values_tn(i) = (1._f64 + epsilon *cos( 9._f64 * atan2(y,x)))*&
                  exp( -4._f64*(r-6.5_f64)**2)
          else
             edge_values_tn(i) = 0._f64
          endif

       enddo

    endif

  end subroutine init_distr_gc

  subroutine init_distr_circ(f_tn,center_values_tn,edge_values_tn,hermite_method, mesh, center_x1, center_x2, sigma, amplitude)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:)       :: f_tn, center_values_tn, edge_values_tn
    sll_int32,  intent(in) :: hermite_method
    sll_real64, intent(in) :: center_x1
    sll_real64, intent(in) :: center_x2
    sll_real64, intent(in) :: sigma
    sll_real64, intent(in) :: amplitude
    sll_real64 :: x, y
    sll_real64 :: r
    sll_int32  :: i

    do i = 1,mesh%num_pts_tot
       x = mesh%cartesian_coord(1,i)
       y = mesh%cartesian_coord(2,i)

       f_tn(i) = amplitude * exp(-0.5_f64* &
            ((x-center_x1)**2 + (y-center_x2)**2) / sigma**2 )
    enddo

        if ( hermite_method == 10 ) then 

       do i = 1,mesh%num_triangles
          x = mesh%center_cartesian_coord(1,i)
          y = mesh%center_cartesian_coord(2,i)
          r = sqrt( x**2 + y**2 )

          center_values_tn(i) = amplitude * exp(-0.5_f64* &
            ((x-center_x1)**2 + (y-center_x2)**2) / sigma**2 )
       enddo
    elseif ( hermite_method == 15 ) then 

       do i = 1,mesh%num_edges
          x = mesh%edge_center_cartesian_coord(1,i)
          y = mesh%edge_center_cartesian_coord(2,i)
          r = sqrt( x**2 + y**2 )

          edge_values_tn(i) = amplitude * exp(-0.5_f64* &
            ((x-center_x1)**2 + (y-center_x2)**2) / sigma**2 )

       enddo

    endif

  end subroutine init_distr_circ

  !-------------------------------------------------------------------------
  !> @brief Writes diagnostics files
  !> @details Write two sort of documents: "diag_gc_hermite*_*.dat" and 
  !> "diag_gc_hermite*_nc.dat". Where important values (i.e. time of sim, errors, number
  !> of cells, etc) are written in order to compute diagnostics. The first file is for
  !> time evolutions, the second one is regarding the space discretization.
  !> @param rho real: contains the value of the density of the gc at time t
  !> @param t real: time of the simulation
  !> @param mesh sll_hex_mesh_2d: hexagonal mesh where the simulation is made
  !> @param uxn real: equals sum( y_i) where (xi,yi) are the mesh points
  !> @param uyn real: equals sum(-x_i) where (xi,yi) are the mesh points
  !> @param nloop int: number of loops done
  !> @param num_meth int: number of the method used for the interpolation method
  !> @param tmax: maximum time that the will simulation will run
  !> @param cells_min int: min number of cells the mesh will have during this simulation
  !> @param cells_max int: max number of cells the mesh will have during this simulation
  !> return 
  subroutine hex_diagnostics_gc(rho,t,mesh,uxn,uyn,nloop,num_meth,tmax,cells_min,cells_max)
    type(sll_hex_mesh_2d),  pointer  :: mesh
    sll_real64, dimension(:) :: rho
    sll_real64, dimension(:) :: uxn
    sll_real64, dimension(:) :: uyn
    sll_real64, intent(in)   :: t
    sll_real64, intent(in)   :: tmax
    sll_int32 , intent(in)   :: nloop
    sll_int32 , intent(in)   :: num_meth
    sll_int32 , intent(in)   :: cells_min, cells_max
    sll_real64 :: mass
    sll_real64 :: rho_min
    sll_real64 :: norm_l1
    sll_real64 :: norm_l2
    sll_real64 :: norm_linf
    sll_real64 :: energy
    sll_int32  :: i
    sll_int32  :: out_unit
    character(len = 50) :: filename
    character(len =  4) :: filenum
    character(len =  4) :: method

    energy    = 0._f64
    mass      = 0._f64
    rho_min   = rho(1)
    norm_l1   = 0._f64
    norm_l2   = 0._f64
    norm_linf = 0._f64

    ! --------------------------------------------------
    ! Writing file in respect to time...................

    if (mesh%num_cells == cells_max) then
       call int2string(mesh%num_cells,filenum)
       call int2string(num_meth, method)
       filename  = "diag_gc_hermite"//trim(method)//"_tmax"//trim(filenum)//".dat"
       
       call sll_new_file_id(out_unit, ierr)
       if (nloop == 0) then
          open(unit = out_unit, file=filename, action="write", status="replace")
       else
          open(unit = out_unit, file=filename, action="write", status="old",position = "append")
       endif

       do i = 1,mesh%num_pts_tot
          mass = mass + rho(i)
          norm_l1 = norm_l1 + abs(rho(i))
          norm_l2 = norm_l2 + rho(i)**2
          energy = energy + uxn(i)**2 + uyn(i)**2
          if ( abs(rho(i)) > norm_linf ) norm_linf = rho(i)
          if ( rho(i) < rho_min  ) rho_min  = rho(i)
       enddo

       print*,"diagnostic for t = ",t
       energy  = sqrt(energy * mesh%delta**2)
       mass    = mass * mesh%delta**2
       norm_l1 = norm_l1 * mesh%delta**2
       norm_l2 = sqrt(norm_l2 * mesh%delta**2)

       write(out_unit,"(7(g18.10,1x))") t, &
            mass, &
            rho_min, &
            norm_l1, &
            norm_l2, &
            norm_linf, &
            energy

       close(out_unit)

    end if
    ! --------------------------------------------------
    ! Writing file in respect to num_cells..............
    if (t.gt.tmax) then !We write on this file only if it is the last time step

       call int2string(num_meth,method)
       filename  = "diag_gc_hermite"//trim(method)//"_nc.dat"

       if ( mesh%num_cells == cells_min ) then
          call sll_new_file_id(out_unit, ierr)
          open(unit = out_unit, file=filename, action="write", status="replace")
       else
          call sll_new_file_id(out_unit, ierr)
          open(unit = out_unit, file=filename, action="write", status="old",position = "append")
       endif
       
       do i = 1,mesh%num_pts_tot
          mass = mass + rho(i)
          norm_l1 = norm_l1 + abs(rho(i))
          norm_l2 = norm_l2 + rho(i)**2
          energy = energy + uxn(i)**2 + uyn(i)**2
          if ( abs(rho(i)) > norm_linf ) norm_linf = rho(i)
          if ( rho(i) < rho_min  ) rho_min  = rho(i)
       enddo

       energy  = sqrt(energy * mesh%delta**2)
       mass    = mass * mesh%delta**2
       norm_l1 = norm_l1 * mesh%delta**2
       norm_l2 = sqrt(norm_l2 * mesh%delta**2)

       write(out_unit,"((i6,1x),7(g18.10,1x))") mesh%num_cells, &
            t, &
            mass, &
            rho_min, &
            norm_l1, &
            norm_l2, &
            norm_linf, &
            energy

       close(out_unit) 
    end if
  end subroutine hex_diagnostics_gc

  
  !-------------------------------------------------------------------------
  !> @brief Writes diagnostics files
  !> @details Write two sort of documents: "diag_circ_hermite*_*.dat" and 
  !> "diag_circ_hermite*_nc.dat". Where important values (i.e. time of sim, errors, number
  !> of cells, etc) are written in order to compute diagnostics. The first file is for
  !> time evolutions, the second one is regarding the space discretization.
  !> @param rho real: contains the value of the density of the circ at time t
  !> @param t real: time of the simulation
  !> @param mesh sll_hex_mesh_2d: hexagonal mesh where the simulation is made
  !> @param nloop int: number of loops done
  !> @param num_meth int: number of the method used for the interpolation method
  !> @param tmax: maximum time that the will simulation will run
  !> @param cells_min int: min number of cells the mesh will have during this simulation
  !> @param cells_max int: max number of cells the mesh will have during this simulation
  !> return 
  subroutine hex_diagnostics_circ(rho,t,mesh,nloop,num_meth,tmax,cells_min,cells_max, &
       gauss_x1, gauss_x2, gauss_sig, gauss_amp)
    type(sll_hex_mesh_2d),  pointer  :: mesh
    sll_real64, dimension(:) :: rho
    sll_real64, intent(in)   :: t
    sll_real64, intent(in)   :: gauss_x1
    sll_real64, intent(in)   :: gauss_x2
    sll_real64, intent(in)   :: gauss_sig
    sll_real64, intent(in)   :: gauss_amp
    sll_real64, intent(in)   :: tmax
    sll_int32 , intent(in)   :: nloop
    sll_int32 , intent(in)   :: num_meth
    sll_int32 , intent(in)   :: cells_min, cells_max
    sll_real64 :: mass
    sll_real64 :: rho_min
    sll_real64 :: norm_l1
    sll_real64 :: norm_l2
    sll_real64 :: norm_linf
    sll_real64 :: x, y, f_exact
    sll_int32  :: i
    sll_int32  :: out_unit
    character(len = 50) :: filename
    character(len =  4) :: filenum
    character(len =  4) :: method

    mass      = 0._f64
    rho_min   = rho(1)
    norm_l1   = 0._f64
    norm_l2   = 0._f64
    norm_linf = 0._f64

    ! --------------------------------------------------
    ! Writing file in respect to time...................

    if (mesh%num_cells == cells_max) then
       call int2string(mesh%num_cells,filenum)
       call int2string(num_meth, method)
       filename  = "diag_circ_hermite"//trim(method)//"_tmax"//trim(filenum)//".dat"
       
       call sll_new_file_id(out_unit, ierr)
       if (nloop == 0) then
          open(unit = out_unit, file=filename, action="write", status="replace")
          write(out_unit,*) &
               "time", ",",&
               "mass", ",", &
               "minimum value", ",", &
               "l_1 norm", ",", &
               "l_2 norm", ",",&
               "l_\inf norm"
       else
          open(unit = out_unit, file=filename, action="write", status="old",position = "append")
       endif

       do i = 1,mesh%num_pts_tot
          ! We compute the exact solution:
          x = mesh%cartesian_coord(1,i)*cos(t) - & 
               mesh%cartesian_coord(2,i)*sin(t)
          y = mesh%cartesian_coord(1,i)*sin(t) + & 
               mesh%cartesian_coord(2,i)*cos(t)
          f_exact = gauss_amp * exp(-0.5_f64* &
               ((x-gauss_x1)**2 + (y-gauss_x2)**2) / gauss_sig**2 )
          ! Then the errors:
          mass = mass + rho(i)
          norm_l1 = norm_l1 + abs(f_exact - rho(i))
          norm_l2 = norm_l2 + (f_exact - rho(i))**2
          if ( abs(rho(i)) > norm_linf ) norm_linf = rho(i)
          if ( rho(i) < rho_min  ) rho_min  = rho(i)
       enddo

       print*,"diagnostic for t = ",t
       mass    = mass * mesh%delta**2
       norm_l1 = norm_l1 * mesh%delta**2
       norm_l2 = sqrt(norm_l2 * mesh%delta**2)

       write(out_unit,"(5(g25.18,1a,1x),1(g25.18))") t, ",",&
            mass, ",", &
            rho_min, ",", &
            norm_l1, ",", &
            norm_l2, ",",&
            norm_linf
              
       close(out_unit)

    end if
    ! --------------------------------------------------
    ! Writing file in respect to num_cells..............
    if (t.gt.tmax) then !We write on this file only if it is the last time step

       call int2string(num_meth, method)
       filename  = "diag_circ_hermite"//trim(method)//"_nc.dat"

       if ( mesh%num_cells == cells_min ) then
          call sll_new_file_id(out_unit, ierr)
          open(unit = out_unit, file=filename, action="write", status="replace")
          write(out_unit,*) &
               "number of cells", ",",&
               "time", ",",&
               "mass", ",", &
               "minimum value", ",", &
               "l_1 norm", ",", &
               "l_2 norm", ",",&
               "l_\inf norm"
       else
          call sll_new_file_id(out_unit, ierr)
          open(unit = out_unit, file=filename, action="write", status="old",position = "append")
       endif
       
       do i = 1,mesh%num_pts_tot
          ! We compute the exact solution:
          x = mesh%cartesian_coord(1,i)*cos(t) - & 
               mesh%cartesian_coord(2,i)*sin(t)
          y = mesh%cartesian_coord(1,i)*sin(t) + & 
               mesh%cartesian_coord(2,i)*cos(t)
          f_exact = gauss_amp * exp(-0.5_f64* &
               ((x-gauss_x1)**2 + (y-gauss_x2)**2) / gauss_sig**2 )
          ! Then the errors:
          mass = mass + rho(i)
          norm_l1 = norm_l1 + abs(f_exact - rho(i))
          norm_l2 = norm_l2 + (f_exact - rho(i))**2
          if ( abs(rho(i)) > norm_linf ) norm_linf = rho(i)
          if ( rho(i) < rho_min  ) rho_min  = rho(i)
       enddo

       mass    = mass * mesh%delta**2
       norm_l1 = norm_l1 * mesh%delta**2
       norm_l2 = sqrt(norm_l2 * mesh%delta**2)

       write(out_unit,"((i6,a,1x),5(g25.18,a,1x),(g25.18))") &
            mesh%num_cells, ",", &
            t, ",", &
            mass, ",", &
            rho_min, ",", &
            norm_l1, ",", &
            norm_l2, ",", &
            norm_linf

       close(out_unit) 
    end if
  end subroutine hex_diagnostics_circ


  ! assembling the values at the middle of the edges + the value at 
  ! the vertices into one common array rho2
  subroutine assemble(rho,rho_edge,rho2,mesh,mesh2)
    type(sll_hex_mesh_2d), pointer :: mesh,mesh2
    sll_real64,dimension(:) :: rho,rho_edge,rho2
    sll_int32               :: i,h1,h2,m1,m2,k1,k2,ne,ns, index_tab
    sll_real64              :: eps = 1.e-6

    do i = 1,mesh2%num_pts_tot 
       h1 = mesh2%hex_coord(1,i)
       h2 = mesh2%hex_coord(2,i)
       m1 = mod(h1,2)
       m2 = mod(h2,2)
       k1 = int(h1/2)
       k2 = int(h2/2)

       if (m1<0) k1 = k1 - 1
       if (m2<0) k2 = k2 - 1

       call mesh%index_hex_to_global(k1, k2, index_tab)
       ns = mesh%global_indices(index_tab)

       if ( abs(m1) > eps .and. abs(m2) > eps ) then
          ne =  mesh%edge_center_index(2,ns)
          rho2(i) = rho_edge(ne) 
       elseif( abs(m1) < eps .and. abs(m2) > eps ) then
          ne =  mesh%edge_center_index(1,ns)
          rho2(i) = rho_edge(ne)
       elseif( abs(m1) > eps .and. abs(m2) < eps ) then 
          ne =  mesh%edge_center_index(3,ns)
          rho2(i) = rho_edge(ne)
       else
          rho2(i) = rho(ns)
       endif
    enddo

  end subroutine assemble

  ! deassembling the values at the middle of the edges + the value at 
  ! the vertices into 2 separate arrays
  subroutine deassemble(rho,rho_edge,rho2,mesh,mesh2)
    type(sll_hex_mesh_2d), pointer :: mesh,mesh2
    sll_real64,dimension(:) :: rho,rho_edge,rho2
    sll_int32               ::i,h1,h2,m1,m2,k1,k2,ne,ns
    sll_real64              :: eps = 1.e-6


    do i = 1,mesh2%num_pts_tot 
       h1 = mesh2%hex_coord(1,i)
       h2 = mesh2%hex_coord(2,i)
       m1 = mod(h1,2)
       m2 = mod(h2,2)
       k1 = int(h1/2)
       k2 = int(h2/2)

       if (m1<0) k1 = k1 - 1
       if (m2<0) k2 = k2 - 1

       call mesh%index_hex_to_global(k1, k2, index_tab)
       ns = mesh%global_indices(index_tab)

       if ( abs(m1) > eps .and. abs(m2) > eps ) then
          ne =  mesh%edge_center_index(2,ns)
          rho_edge(ne) = rho2(i)
       elseif( abs(m1) < eps .and. abs(m2) > eps ) then
          ne =  mesh%edge_center_index(1,ns)
          rho_edge(ne) = rho2(i)
       elseif( abs(m1) > eps .and. abs(m2) < eps ) then 
          ne =  mesh%edge_center_index(3,ns)
          rho_edge(ne) = rho2(i)
       else
          rho(ns) = rho2(i)
       endif
    enddo

  end subroutine deassemble
  

  subroutine write_center(mesh,rho,name)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64,dimension(:)        :: rho
    sll_int32                      :: i
    character(len = 50)            :: name
    sll_real64                     :: x,y

    open (unit = 11,file=name, action="write", status="replace")

    do i = 1,mesh%num_triangles
       x = mesh%center_cartesian_coord(1,i)
       y = mesh%center_cartesian_coord(2,i)
       write(11,*) x,y,rho(i)
    enddo

    close (11)

  end subroutine write_center

end program sim_bsl_gc_2d0v_hex_hermite
