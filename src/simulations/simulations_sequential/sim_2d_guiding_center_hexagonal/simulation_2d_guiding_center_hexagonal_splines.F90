program sim2d_gc_hex_splines

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_constants
  use euler_2d_hex
  use sll_hex_meshes
  use hex_poisson
  use pivotbande
  use sll_box_splines
  implicit none

  type(sll_hex_mesh_2d),   pointer        :: mesh
  type(sll_box_spline_2d), pointer        :: spline
  sll_real64, dimension(:),   allocatable :: dxuxn, dyuxn
  sll_real64, dimension(:),   allocatable :: dxuyn, dyuyn
  sll_real64, dimension(:),   allocatable :: second_term
  sll_real64, dimension(:),   allocatable :: phi
  sll_real64, dimension(:),   allocatable :: uxn, uxn_1
  sll_real64, dimension(:),   allocatable :: uyn, uyn_1
  sll_real64, dimension(:),   allocatable :: phi_interm
  sll_real64, dimension(:,:), allocatable :: matrix_poisson, l, u
  sll_real64, dimension(:),   allocatable :: rho_tn_1   ! distribution at time n-1
  sll_real64, dimension(:),   allocatable :: rho_tn   ! distribution at time n
  sll_real64, dimension(:),   allocatable :: rho_tn1  ! distribution at time n+1
  sll_real64, dimension(:),   allocatable :: x1_char
  sll_real64, dimension(:),   allocatable :: x2_char

  sll_int32    :: spline_degree
  sll_int32    :: i,j, k1, k2, index_tab, type
  sll_int32    :: l1,l2
  sll_int32    :: i1,i2,i3
  sll_int32    :: num_cells, n_points, n_triangle, n_points2
  sll_int32    :: cells_min, cells_max
  sll_int32    :: cells_stp
  sll_int32    :: nloops,count, ierr, EXTRA_TABLES = 0
  sll_real64   :: center_mesh_x1, center_mesh_x2, radius
  sll_real64   :: epsilon
  sll_real64   :: dt
  sll_real64   :: tmax
  sll_real64   :: t
  sll_real64   :: t_init, t_end
  sll_real64   :: t1,t2,t3,t4,t5,t6
  sll_real64   :: h1, h2, f_min, x ,y,xx, yy
  sll_real64   :: r11,r12,r21,r22,det
  logical      :: inside
  sll_int32    :: p = 6 ! degree of the approximation for the derivative
  sll_int32            :: IO_stat
  sll_int32, parameter :: input_file = 99
  character(len = 256) :: input_filename
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

  namelist /initial_function/ &
!      initial_function_case, &
!      kmode_x1, &
!      kmode_x2, &
      epsilon

    namelist /time_iterations/ &
      dt, &
      tmax
!      freq_diag, &
!      freq_diag_time, &
!      time_loop_case

    namelist /interpolation/ &
      spline_degree


    ! namelist /advector/ &
    !   advect2d_case, &   
    !   f_interp2d_case, &
    !   phi_interp2d_case, &
    !   charac2d_case, &
    !   A_interp_case, &
    !   charac1d_x1_case, &
    !   charac1d_x2_case, &
    !   advect1d_x1_case, &   
    !   advect1d_x2_case  

    ! namelist /poisson/ &
    !   poisson_case, &
    !   poisson_solver, &
    !   mudpack_method, &    
    !   spline_degree_eta1, &
    !   spline_degree_eta2    

  ! ----------------------------
  ! Setting default parameters
  ! ----------------------------
  ! Mesh :
  center_mesh_x1 = 0._f64
  center_mesh_x2 = 0._f64
  radius = 14._f64
  cells_min = 20
  cells_max = 160
  cells_stp = 20
  ! Initial function :
  epsilon = 0.001_f64
  ! Time iterations:
  tmax  = 20._f64
  dt    = 0.1_f64
  ! Interpolation
  spline_degree = 2

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
     read(input_file, initial_function)
     read(input_file, time_iterations)
     read(input_file, interpolation)
     ! read(input_file, advector)
     ! read(input_file, poisson)
     close(input_file)
  else
     print *,'#initialization with default parameters'
  endif


  do num_cells = cells_min,cells_max,cells_stp

     call int2string(spline_degree,degnum)
     print*, ""
     print*," ********************************* "
     print*,"     Guiding-Center Simulation"
     print*,"        on a Hexagonal mesh"
     print*,"   using boxsplines of deg =", degnum
     print*," ********************************* "

     t = 0._f64
     nloops = 0
     count  = 0

     !*********************************************************
     !             allocation
     !*********************************************************

     ! Mesh creation -------------------------
     mesh => new_hex_mesh_2d( num_cells, center_mesh_x1, center_mesh_x2,&
         radius=radius, EXTRA_TABLES = EXTRA_TABLES )
     call display_hex_mesh_2d(mesh)

     n_points   = mesh%num_pts_tot
     n_triangle = mesh%num_triangles

     det = (mesh%r1_x1*mesh%r2_x2 - mesh%r1_x2*mesh%r2_x1)/mesh%delta

     r11 = + mesh%r2_x2/det
     r12 = - mesh%r2_x1/det
     r21 = - mesh%r1_x2/det
     r22 = + mesh%r1_x1/det
     ! ---------------------------------------

     ! TODO : what are these  ?
     l1 = 2*num_cells+1
     l2 = l1

     SLL_ALLOCATE(rho_tn_1( n_points),ierr)
     SLL_ALLOCATE(rho_tn( n_points),ierr)
     SLL_ALLOCATE(rho_tn1( n_points ),ierr)
     SLL_ALLOCATE(x1_char( n_points ),ierr)
     SLL_ALLOCATE(x2_char( n_points ),ierr)

     SLL_ALLOCATE(uxn( n_points),ierr)
     SLL_ALLOCATE(uyn( n_points ),ierr)
     SLL_ALLOCATE(uxn_1( n_points),ierr)
     SLL_ALLOCATE(uyn_1( n_points ),ierr)

     SLL_ALLOCATE(dxuxn( n_points),ierr)
     SLL_ALLOCATE(dxuyn( n_points ),ierr)
     SLL_ALLOCATE(dyuxn( n_points),ierr)
     SLL_ALLOCATE(dyuyn( n_points ),ierr)

     SLL_ALLOCATE(second_term( n_points),ierr)
     SLL_ALLOCATE(phi( n_points),ierr)
     SLL_ALLOCATE(phi_interm( n_points),ierr)
     SLL_ALLOCATE(matrix_poisson( n_points,1 + 4*num_cells + 2 ) , ierr)
     SLL_ALLOCATE(l( n_points,1 + 4*num_cells + 2 ) , ierr)
     SLL_ALLOCATE(u( n_points,1 + 4*num_cells + 2), ierr)

     call cpu_time(t_init)

     !*********************************************************
     !  Distribution function & density initialization
     !*********************************************************

     ! Spline initialization -----------------
     spline => new_box_spline_2d(mesh, SLL_DIRICHLET)
     ! ---------------------------------------

     ! Initial distribution ------------------
     call init_distr(rho_tn,mesh,epsilon)
     ! ---------------------------------------

     ! Poisson solver ------------------------
     call hex_matrix_poisson( matrix_poisson, mesh,type=1)
     call factolub_bande(matrix_poisson,l,u,n_points,l1,l2)
     call hex_second_terme_poisson( second_term, mesh, rho_tn )
     call solvlub_bande(l,u,phi_interm,second_term,n_points,l1,l2)

     do i = 1, mesh%num_pts_tot
        k1 = mesh%hex_coord(1, i)
        k2 = mesh%hex_coord(2, i)
        call index_hex_to_global(mesh, k1, k2, index_tab)
        phi(i) = phi_interm(index_tab)
     enddo
     call compute_hex_fields(mesh,uxn,uyn,dxuxn,dyuxn,dxuyn,dyuyn,phi,type=1)
     ! ---------------------------------------

     call hex_diagnostics(rho_tn,t,mesh,uxn,uyn,nloops,spline_degree,tmax,cells_min,cells_max)


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

        call compute_coeff_box_spline_2d( rho_tn, spline_degree, spline )

        do i=1, n_points

           x = mesh%cartesian_coord(1,i)
           y = mesh%cartesian_coord(2,i)

           !*************************************************
           !       computation of the characteristics
           !*************************************************

           if ( t <= dt + 1e-6 ) then ! first step with euler

              call compute_characteristic_euler_2d_hex( &
                   x,y,uxn,uyn,i,xx,yy,dt )

           else !the rest is done with Adams 2


              call compute_characteristic_adams2_2d_hex( x,y,uxn,uyn,uxn_1,uyn_1,&
                   dxuxn,dyuxn,dxuyn,dyuyn,i,xx,yy,dt)

           endif

           
           inside = .true.
           h1 =  xx*r11 + yy*r12
           h2 =  xx*r21 + yy*r22

           if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
           if ( abs(xx) > (radius-mesh%delta)*sqrt(3._f64)*0.5_f64) inside = .false.

           if ( inside ) then
              rho_tn1(i) = hex_interpolate_value(mesh, xx, yy, spline, spline_degree)
           else
              rho_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif


        end do ! end of the computation of the mesh points


        ! Updating the new field values ............
        rho_tn_1 = rho_tn
        rho_tn = rho_tn1
        ! ...........................................

        !*********************************************************
        !      computing the solution of the poisson equation 
        !*********************************************************

        call hex_second_terme_poisson( second_term, mesh, rho_tn ) 

        call solvlub_bande(l,u,phi_interm,second_term,n_points,l1,l2)

        do i = 1, mesh%num_pts_tot    ! need to re-index phi :
           k1 = mesh%hex_coord(1, i)
           k2 = mesh%hex_coord(2, i)
           call index_hex_to_global(mesh, k1, k2, index_tab)
           phi(i) = phi_interm(index_tab)
        enddo

        call compute_hex_fields(mesh,uxn,uyn,dxuxn,dyuxn,dxuyn,dyuyn,phi,type=1)

        ! Updating the new field values ............
        uxn_1 = uxn
        uyn_1 = uyn
        ! ...........................................

        !*********************************************************
        !                  writing diagostics
        !*********************************************************
        call hex_diagnostics(rho_tn,t,mesh,uxn,uyn,nloops,spline_degree,tmax,cells_min,cells_max)
        if (count == 10.and.nloops<10000.and.num_cells == cells_max) then
           call int2string(nloops,filenum)
           filename  = "center_guide_rho"//trim(filenum)
           call write_field_hex_mesh_xmf(mesh, rho_tn1, trim(filename))
           filename  = "center_guide_phi"//trim(filenum)
           call write_field_hex_mesh_xmf(mesh, phi, trim(filename))
           count = 0
        endif


     enddo

     SLL_DEALLOCATE_ARRAY(rho_tn_1,ierr)
     SLL_DEALLOCATE_ARRAY(rho_tn,ierr)
     SLL_DEALLOCATE_ARRAY(rho_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(x1_char,ierr)
     SLL_DEALLOCATE_ARRAY(x2_char,ierr)
     SLL_DEALLOCATE_ARRAY(uxn,ierr)
     SLL_DEALLOCATE_ARRAY(uyn,ierr)
     SLL_DEALLOCATE_ARRAY(uxn_1,ierr)
     SLL_DEALLOCATE_ARRAY(uyn_1,ierr)
     SLL_DEALLOCATE_ARRAY(dxuxn,ierr)
     SLL_DEALLOCATE_ARRAY(dxuyn,ierr)
     SLL_DEALLOCATE_ARRAY(dyuxn,ierr)
     SLL_DEALLOCATE_ARRAY(dyuyn,ierr)
     SLL_DEALLOCATE_ARRAY(second_term,ierr)
     SLL_DEALLOCATE_ARRAY(phi,ierr)
     SLL_DEALLOCATE_ARRAY(phi_interm,ierr)
     deallocate(matrix_poisson)
     deallocate(l)
     deallocate(u)

     call delete_hex_mesh_2d( mesh )

     call cpu_time(t_end)
     print*, "time used =", t_end - t_init

  end do

  ! Some test should probably be put here:
  print*, 'PASSED'

contains

  !*********initialization**************

  subroutine init_distr(f_tn, mesh, epsilon)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:)       :: f_tn
    sll_real64, intent(in) :: epsilon
    sll_real64 :: x, y
    sll_real64 :: rho
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


  end subroutine init_distr


  !-------------------------------------------------------------------------
  !> @brief Writes diagnostics files
  !> @param Write two sort of documents: "diag_gc_spline*_*.dat" and 
  !> "diag_gc_spline*_nc.dat". Where important values (i.e. time of sim, errors, number
  !> of cells, etc) are written in order to compute diagnostics. The first file is for
  !> time evolutions, the second one is regarding the space discretization.
  !> @param rho real: contains the value of the density of the gc at time t
  !> @param t real: time of the simulation
  !> @param mesh sll_hex_mesh_2d: hexagonal mesh where the simulation is made
  !> @param uxn real:
  !> @param uyn real:
  !> @param nloop int: number of loops done
  !> @param deg int: degree of the splines used for the interpolation method
  !> @param tmax: maximum time that the will simulation will run
  !> @param cells_min int: min number of cells the mesh will have during this simulation
  !> @param cells_max int: max number of cells the mesh will have during this simulation
  !> return 
  subroutine hex_diagnostics(rho,t,mesh,uxn,uyn,nloop,deg,tmax,cells_min,cells_max)
    type(sll_hex_mesh_2d),  pointer  :: mesh
    sll_real64, dimension(:) :: rho
    sll_real64, dimension(:) :: uxn
    sll_real64, dimension(:) :: uyn
    sll_real64, intent(in)   :: t
    sll_real64, intent(in)   :: tmax
    sll_int32 , intent(in)   :: nloop
    sll_int32 , intent(in)   :: deg
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
    character(len =  4) :: splinedeg

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
       call int2string(deg,splinedeg)
       filename  = "diag_gc_spline"//trim(splinedeg)//"_tmax"//trim(filenum)//".dat"
       
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

       call int2string(deg,splinedeg)
       filename  = "diag_gc_spline"//trim(splinedeg)//"_nc.dat"

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


  end subroutine hex_diagnostics


  subroutine assemble(rho,rho_edge,rho2,mesh,mesh2)
    type(sll_hex_mesh_2d), pointer :: mesh,mesh2
    sll_real64,dimension(:) :: rho,rho_edge,rho2
    sll_int32               :: i,j,k,i1,i2

    rho2(1) = rho(1)

    rho2(2) = rho_edge(3)
    rho2(3) = rho_edge(2)
    rho2(4) = rho_edge(1)
    rho2(5) = rho_edge(15)
    rho2(6) = rho_edge(17)
    rho2(7) = rho_edge(19)

    ! à optimiser quand ce sera validé

    do i = 1,mesh%num_cells!mesh%num_pts_tot
       do k = 1,6*i
          j  = 6*i*(2*i-1)+2*k
          i1 = 1+(i-1)*6+k
          rho2(j)   = rho(i1)

       enddo

       do k = 1,i
          ! first edge
          j  = 6*i*(2*i-1)+2*k
          i1 = 1+3*i*(i-1)+k
          i2 = mesh%edge_center_index(1,i1)
          rho2(j+1) = rho_edge(i2)
          ! second edge
          j  = 6*i*(2*i-1)+2*(k+i)
          i1 = 1+3*i*(i-1)+k+  i+1
          i2 = mesh%edge_center_index(3,i1)
          rho2(j+1) = rho_edge(i2)
          ! third edge
          j  = 6*i*(2*i-1)+2*(k+2*i)
          i1 = 1+3*i*(i-1)+k+2*i+1
          i2 = mesh%edge_center_index(2,i1)
          rho2(j+1) = rho_edge(i2)
          ! fourth edge
          j  = 6*i*(2*i-1)+2*(k+3*i)
          i1 = 1+3*i*(i-1)+k+3*i+1
          i2 = mesh%edge_center_index(1,i1)
          rho2(j+1) = rho_edge(i2)
          ! fifth edge
          j  = 6*i*(2*i-1)+2*(k+4*i)
          i1 = 1+3*i*(i-1)+k+4*i
          i2 = mesh%edge_center_index(3,i1)
          rho2(j+1) = rho_edge(i2)
          ! sixth edge
          j  = 6*i*(2*i-1)+2*(k+5*i)
          i1 = 1+3*i*(i-1)+k+5*i
          i2 = mesh%edge_center_index(2,i1)
          rho2(j+1) = rho_edge(i2)

       enddo

       if (i<mesh%num_cells) then

          !corners

          j  = 6*i*(2*i-1)+2
          i2 = mesh%edge_center_index(3,i1)
          rho2(j) = rho_edge(i2)

          j  = 6*i*(2*i-1)+2+i
          i2 = mesh%edge_center_index(3,i1)
          rho2(j) = rho_edge(i2)


          do k = 2,i-1

             ! first edge
             j  = 6*i*(2*i-1)+2*k
             i1 = 1+3*i*(i-1)+k
             i2 = mesh%edge_center_index(3,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(2,i1)
             rho2(j+1) = rho_edge(i2)

             ! second edge
             j  = 6*i*(2*i-1)+(k+i)
             i1 = 1+3*i*(i-1)+k+  i+1
             i2 = mesh%edge_center_index(2,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(1,i1)
             rho2(j+1) = rho_edge(i2)
             ! third edge
             j  = 6*i*(2*i-1)+(k+2*i)
             i1 = 1+3*i*(i-1)+k+2*i+1
             i2 = mesh%edge_center_index(1,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(3,i1)
             rho2(j+1) = rho_edge(i2)
             ! fourth edge
             j  = 6*i*(2*i-1)+(k+3*i)
             i1 = 1+3*i*(i-1)+k+3*i+1
             i2 = mesh%edge_center_index(2,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(3,i1)
             rho2(j+1) = rho_edge(i2)
             ! fifth edge
             j  = 6*i*(2*i-1)+(k+4*i)
             i1 = 1+3*i*(i-1)+k+4*i
             i2 = mesh%edge_center_index(1,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(2,i1)
             rho2(j+1) = rho_edge(i2)
             ! sixth edge
             j  = 6*i*(2*i-1)+2*(k+5*i)
             i1 = 1+3*i*(i-1)+k+5*i
             i2 = mesh%edge_center_index(1,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(3,i1)
             rho2(j+1) = rho_edge(i2)

          enddo
       endif


    enddo


  end subroutine assemble

  subroutine deassemble(rho,rho_edge,rho2,mesh,mesh2)
    type(sll_hex_mesh_2d), pointer :: mesh,mesh2
    sll_real64,dimension(:) :: rho,rho_edge,rho2
    sll_int32               :: i,j,k,i1,i2,j1,j2

    rho(1) = rho2(1)

    do i = 1,mesh%num_cells!mesh2%num_pts_tot
       i1 = 2*i
       i2 = 2*i + 1
       ! for i even (hexagone of middle point)
       do k = 1,i2
          rho_edge(j) = rho(j1)
       enddo

       ! for i odd (hexagone mix)
       do k = 1,i
          rho(j1) = rho2(j)
          rho_edge(j2) = rho(j+1)

       enddo
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


end program sim2d_gc_hex_splines
