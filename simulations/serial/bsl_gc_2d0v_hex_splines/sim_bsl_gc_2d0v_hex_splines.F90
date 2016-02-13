program sim_bsl_gc_2d0v_hex_splines

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet

  use sll_m_box_splines, only: &
    sll_f_hex_interpolate_value, &
    sll_f_new_box_spline_2d, &
    sll_t_box_spline_2d

  use sll_m_constants, only: &
    sll_p_sqrt3

  use sll_m_euler_2d_hex, only: &
    sll_s_compute_characteristic_adams2_2d_hex, &
    sll_s_compute_characteristic_euler_2d_hex

  use sll_m_hex_poisson, only: &
    sll_s_compute_hex_fields, &
    sll_s_hex_matrix_poisson, &
    sll_s_hex_second_terme_poisson

  use sll_m_hexagonal_meshes, only: &
    sll_s_delete_hex_mesh_2d, &
    sll_s_display_hex_mesh_2d, &
    sll_f_new_hex_mesh_2d, &
    sll_t_hex_mesh_2d

  use sll_m_pivotbande, only: &
    sll_s_factolub_bande, &
    sll_s_solvlub_bande

  use sll_m_utilities, only: &
    sll_s_int2string, &
    sll_s_new_file_id

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_hex_mesh_2d),   pointer        :: mesh
  type(sll_t_box_spline_2d), pointer        :: spline
  sll_real64, dimension(:),   allocatable :: dxuxn, dyuxn
  sll_real64, dimension(:),   allocatable :: dxuyn, dyuyn
  sll_real64, dimension(:),   allocatable :: second_term
  sll_real64, dimension(:),   allocatable :: phi
  sll_real64, dimension(:),   allocatable :: uxn, uxn_1
  sll_real64, dimension(:),   allocatable :: uyn, uyn_1
  sll_real64, dimension(:),   allocatable :: phi_interm
  sll_real64, dimension(:,:), allocatable :: matrix_poisson, l, u
  sll_real64, dimension(:),   allocatable :: rho_tn   ! distribution at time n
  sll_real64, dimension(:),   allocatable :: rho_tn1  ! distribution at time n+1

  sll_int32    :: spline_degree
  sll_int32    :: hermite_method
  sll_int32    :: i, k1, k2, index_tab, type
  sll_int32    :: width_band1,width_band2
  sll_int32    :: num_cells, n_points
  sll_int32    :: cells_min, cells_max
  sll_int32    :: cells_stp
  sll_int32    :: nloops,count, ierr, EXTRA_TABLES = 0
  sll_real64   :: center_mesh_x1, center_mesh_x2, radius
  sll_real64   :: epsilon
  sll_real64   :: gauss_x1, gauss_x2
  sll_real64   :: gauss_sig
  sll_real64   :: gauss_amp
  sll_real64   :: dt
  sll_real64   :: tmax
  sll_real64   :: t
  sll_real64   :: t_init, t_end
  sll_real64   :: t3
  sll_real64   :: h1, h2, x ,y,xx, yy
  sll_real64   :: r11,r12,r21,r22,det
  logical      :: inside
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


  do num_cells = cells_min,cells_max,cells_stp

     call sll_s_int2string(spline_degree,degnum)
     if (model_name.eq."GC") then
        print*, ""
        print*," ********************************* "
        print*,"     Guiding-Center Simulation"
        print*,"        on a Hexagonal mesh"
        print*,"   using boxsplines of deg =", degnum
        print*," ********************************* "
     elseif (model_name.eq."CIRCULAR") then
        print*, ""
        print*," ********************************* "
        print*,"     Circular Advection Simulation"
        print*,"        on a Hexagonal mesh"
        print*,"   using boxsplines of deg =", degnum
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
        print *, " ----> dt =", dt 
     end if
     !*********************************************************
     !             allocation
     !*********************************************************

     ! Mesh creation -------------------------
     mesh => sll_f_new_hex_mesh_2d( num_cells, center_mesh_x1, center_mesh_x2,&
          radius=radius, EXTRA_TABLES = EXTRA_TABLES )
     call sll_s_display_hex_mesh_2d(mesh)

     n_points   = mesh%num_pts_tot
     det = (mesh%r1_x1*mesh%r2_x2 - mesh%r1_x2*mesh%r2_x1)/mesh%delta

     r11 = + mesh%r2_x2/det
     r12 = - mesh%r2_x1/det
     r21 = - mesh%r1_x2/det
     r22 = + mesh%r1_x1/det
     ! ---------------------------------------

     width_band1 = 2*num_cells+1
     width_band2 = width_band1

     SLL_ALLOCATE(rho_tn( n_points),ierr)
     SLL_ALLOCATE(rho_tn1( n_points ),ierr)

     SLL_ALLOCATE(uxn( n_points),ierr)
     SLL_ALLOCATE(uyn( n_points ),ierr)

     SLL_ALLOCATE(dxuxn( n_points),ierr)
     SLL_ALLOCATE(dxuyn( n_points ),ierr)
     SLL_ALLOCATE(dyuxn( n_points),ierr)
     SLL_ALLOCATE(dyuyn( n_points ),ierr)

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


     call cpu_time(t_init)

     !*********************************************************
     !  Distribution function & density initialization
     !*********************************************************

     ! Spline initialization -----------------
     spline => sll_f_new_box_spline_2d(mesh, sll_p_dirichlet)
     ! ---------------------------------------

     ! Initial distribution ------------------
     if (model_name.eq."GC") then
        call init_distr_gc(rho_tn,mesh,epsilon)
     elseif (model_name.eq."CIRCULAR") then
        call init_distr_circ(rho_tn,mesh,gauss_x1,gauss_x2,gauss_sig,gauss_amp)
     end if
     ! ---------------------------------------

     ! Poisson solver ------------------------
     if (model_name.eq."GC") then
        call sll_s_hex_matrix_poisson( matrix_poisson, mesh,type=1)
        call sll_s_factolub_bande(matrix_poisson,l,u,n_points,width_band1,width_band2)
        call sll_s_hex_second_terme_poisson( second_term, mesh, rho_tn )
        call sll_s_solvlub_bande(l,u,phi_interm,second_term,n_points,width_band1,width_band2)

        do i = 1, mesh%num_pts_tot
           k1 = mesh%hex_coord(1, i)
           k2 = mesh%hex_coord(2, i)
           call mesh%index_hex_to_global(k1, k2, index_tab)
           phi(i) = phi_interm(index_tab)
        enddo
        call sll_s_compute_hex_fields(mesh,uxn,uyn,dxuxn,dyuxn,dxuyn,dyuyn,phi,type=1)
     end if
     ! ---------------------------------------

     if (model_name.eq."GC") then
        call hex_diagnostics_gc(rho_tn,t,mesh,uxn,uyn,nloops,spline_degree,tmax,cells_min,cells_max)
        call sll_s_int2string(nloops,filenum)
        filename  = "guiding_center_rho"//trim(filenum)
        call mesh%write_field_hex_mesh_xmf( rho_tn, trim(filename))
        filename  = "guiding_center_phi"//trim(filenum)
        call mesh%write_field_hex_mesh_xmf( phi, trim(filename))
     elseif (model_name.eq."CIRCULAR") then
        call hex_diagnostics_circ(rho_tn, t, &
             mesh, nloops, spline_degree, &
             tmax, cells_min, cells_max, &
             gauss_x1, gauss_x2, gauss_sig, gauss_amp)
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

        call spline%compute_coeff_box_spline_2d( rho_tn, spline_degree)

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
                 call sll_s_compute_characteristic_euler_2d_hex( &
                      x,y,uxn,uyn,i,xx,yy,dt )
              else !the rest is done with Adams 2
                 call sll_s_compute_characteristic_adams2_2d_hex(x, y, &
                      uxn, uyn, uxn_1, uyn_1, &
                      dxuxn,dyuxn,dxuyn,dyuyn, &
                      i, xx, yy, dt)
              endif
           end if

           inside = .true.
           h1 =  xx*r11 + yy*r12
           h2 =  xx*r21 + yy*r22

           if ( abs(h1) > radius-mesh%delta .or. abs(h2) > radius-mesh%delta ) &
                inside = .false.
           if ( abs(xx) > (radius-mesh%delta)*sll_p_sqrt3*0.5_f64) &
                inside = .false.

           if ( inside ) then
              rho_tn1(i) = sll_f_hex_interpolate_value(mesh, xx, yy, &
                   spline, spline_degree)
           else
              rho_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif


        end do ! end of the computation of the mesh points

        ! Updating the new field values ............
        rho_tn = rho_tn1
        ! ...........................................

        !*********************************************************
        !      computing the solution of the poisson equation
        !*********************************************************
        if (model_name.eq."GC") then
           call sll_s_hex_second_terme_poisson( second_term, mesh, rho_tn )

           call sll_s_solvlub_bande(l,u,phi_interm,second_term,n_points,&
                width_band1,width_band2)

           do i = 1, mesh%num_pts_tot    ! need to re-index phi :
              k1 = mesh%hex_coord(1, i)
              k2 = mesh%hex_coord(2, i)
              call mesh%index_hex_to_global(k1, k2, index_tab)
              phi(i) = phi_interm(index_tab)
           enddo

           call sll_s_compute_hex_fields(mesh,uxn,uyn,&
                dxuxn,dyuxn,dxuyn,dyuyn,phi,type=1)

           ! Updating the new field values ............
           uxn_1 = uxn
           uyn_1 = uyn
           ! ...........................................
        end if

        !*********************************************************
        !                  writing diagostics
        !*********************************************************
        if (model_name.eq."GC") then
           call hex_diagnostics_gc(rho_tn,t,mesh,uxn,uyn,nloops,spline_degree,&
                tmax,cells_min,cells_max)
           if (count == 10.and.nloops<10000.and.num_cells == cells_max) then
              call sll_s_int2string(nloops,filenum)
              filename  = "guiding_center_rho"//trim(filenum)
              call mesh%write_field_hex_mesh_xmf(rho_tn1, trim(filename))
              if (model_name.eq."GC") then
                 filename  = "guiding_center_phi"//trim(filenum)
                 call mesh%write_field_hex_mesh_xmf(phi, trim(filename))
              end if
              count = 0
           endif
        elseif (model_name.eq."CIRCULAR") then
           call hex_diagnostics_circ(rho_tn,t,mesh,nloops,spline_degree,tmax,&
                cells_min,cells_max, gauss_x1, gauss_x2, gauss_sig, gauss_amp)
           if (count == 10.and.nloops<10000.and.num_cells == cells_max) then
              call sll_s_int2string(nloops,filenum)
              filename  = "circular_advection_rho"//trim(filenum)
              call mesh%write_field_hex_mesh_xmf(rho_tn1, trim(filename))
              count = 0
           endif
        end if

     enddo ! end of time loop

     SLL_DEALLOCATE_ARRAY(rho_tn,ierr)
     SLL_DEALLOCATE_ARRAY(rho_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(uxn,ierr)
     SLL_DEALLOCATE_ARRAY(uyn,ierr)
     SLL_DEALLOCATE_ARRAY(dxuxn,ierr)
     SLL_DEALLOCATE_ARRAY(dxuyn,ierr)
     SLL_DEALLOCATE_ARRAY(dyuxn,ierr)
     SLL_DEALLOCATE_ARRAY(dyuyn,ierr)
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

     call sll_s_delete_hex_mesh_2d( mesh )

     call cpu_time(t_end)
     print*, "time used =", t_end - t_init

  end do ! end of loop on space step (end of all simulations)

  ! Some test should probably be put here:
  print*, 'PASSED'

contains

  !*********initialization**************

  subroutine init_distr_gc(f_tn, mesh, epsilon)
    type(sll_t_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:)       :: f_tn
    sll_real64, intent(in) :: epsilon
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
  end subroutine init_distr_gc

  subroutine init_distr_circ(f_tn, mesh, center_x1, center_x2, sigma, amplitude)
    type(sll_t_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:)       :: f_tn
    sll_real64, intent(in) :: center_x1
    sll_real64, intent(in) :: center_x2
    sll_real64, intent(in) :: sigma
    sll_real64, intent(in) :: amplitude
    sll_real64 :: x, y
    sll_int32  :: i

    do i = 1,mesh%num_pts_tot
       x = mesh%cartesian_coord(1,i)
       y = mesh%cartesian_coord(2,i)

       f_tn(i) = amplitude * exp(-0.5_f64* &
            ((x-center_x1)**2 + (y-center_x2)**2) / sigma**2 )
    enddo
  end subroutine init_distr_circ

  !-------------------------------------------------------------------------
  !> @brief Writes diagnostics files
  !> @details Write two sort of documents: "diag_gc_spline*_*.dat" and 
  !> "diag_gc_spline*_nc.dat". Where important values (i.e. time of sim, errors, number
  !> of cells, etc) are written in order to compute diagnostics. The first file is for
  !> time evolutions, the second one is regarding the space discretization.
  !> @param rho real: contains the value of the density of the gc at time t
  !> @param t real: time of the simulation
  !> @param mesh sll_t_hex_mesh_2d: hexagonal mesh where the simulation is made
  !> @param uxn real: equals sum( y_i) where (xi,yi) are the mesh points
  !> @param uyn real: equals sum(-x_i) where (xi,yi) are the mesh points
  !> @param nloop int: number of loops done
  !> @param deg int: degree of the splines used for the interpolation method
  !> @param tmax: maximum time that the will simulation will run
  !> @param cells_min int: min number of cells the mesh will have during this simulation
  !> @param cells_max int: max number of cells the mesh will have during this simulation
  !> return 
  subroutine hex_diagnostics_gc(rho,t,mesh,uxn,uyn,nloop,deg,tmax,cells_min,cells_max)
    type(sll_t_hex_mesh_2d),  pointer  :: mesh
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
       call sll_s_int2string(mesh%num_cells,filenum)
       call sll_s_int2string(deg,splinedeg)
       filename  = "diag_gc_spline"//trim(splinedeg)//"_tmax"//trim(filenum)//".dat"

       call sll_s_new_file_id(out_unit, ierr)
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

       print *,"diagnostic for t = ",t
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

       call sll_s_int2string(deg,splinedeg)
       filename  = "diag_gc_spline"//trim(splinedeg)//"_nc.dat"

       if ( mesh%num_cells == cells_min ) then
          call sll_s_new_file_id(out_unit, ierr)
          open(unit = out_unit, file=filename, action="write", status="replace")
       else
          call sll_s_new_file_id(out_unit, ierr)
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
  !> @details Write two sort of documents: "diag_circ_spline*_*.dat" and 
  !> "diag_circ_spline*_nc.dat". Where important values (i.e. time of sim, errors, number
  !> of cells, etc) are written in order to compute diagnostics. The first file is for
  !> time evolutions, the second one is regarding the space discretization.
  !> @param rho real: contains the value of the density of the circ at time t
  !> @param t real: time of the simulation
  !> @param mesh sll_t_hex_mesh_2d: hexagonal mesh where the simulation is made
  !> @param nloop int: number of loops done
  !> @param deg int: degree of the splines used for the interpolation method
  !> @param tmax: maximum time that the will simulation will run
  !> @param cells_min int: min number of cells the mesh will have during this simulation
  !> @param cells_max int: max number of cells the mesh will have during this simulation
  !> @param gauss_x1 real: 1st coordinate of center of the bulb distribution function
  !> @param gauss_x2 real: 2nd coordinate of center of the bulb distribution function
  !> @param gauss_sig real: width of the bulb distribution function
  !> @param gauss_amp real: amplitude of the bulb distribution function
  !> return 
  subroutine hex_diagnostics_circ(rho,t,mesh,nloop,deg,tmax,cells_min,cells_max, &
       gauss_x1, gauss_x2, gauss_sig, gauss_amp)
    type(sll_t_hex_mesh_2d),  pointer  :: mesh
    sll_real64, dimension(:) :: rho
    sll_real64, intent(in)   :: t
    sll_real64, intent(in)   :: tmax
    sll_real64, intent(in)   :: gauss_x1
    sll_real64, intent(in)   :: gauss_x2
    sll_real64, intent(in)   :: gauss_sig
    sll_real64, intent(in)   :: gauss_amp
    sll_int32 , intent(in)   :: nloop
    sll_int32 , intent(in)   :: deg
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
    character(len =  4) :: splinedeg

    mass      = 0._f64
    rho_min   = rho(1)
    norm_l1   = 0._f64
    norm_l2   = 0._f64
    norm_linf = 0._f64

    ! --------------------------------------------------
    ! Writing file in respect to time...................
    if (mesh%num_cells == cells_max) then
       call sll_s_int2string(mesh%num_cells,filenum)
       call sll_s_int2string(deg,splinedeg)
       filename  = "diag_circ_spline"//trim(splinedeg)//"_tmax"//trim(filenum)//".txt"

       call sll_s_new_file_id(out_unit, ierr)
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

       call sll_s_int2string(deg,splinedeg)
       filename  = "diag_circ_spline"//trim(splinedeg)//"_nc.txt"

       if ( mesh%num_cells == cells_min ) then
          call sll_s_new_file_id(out_unit, ierr)
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
          call sll_s_new_file_id(out_unit, ierr)
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

  subroutine write_center(mesh,rho,name)
    type(sll_t_hex_mesh_2d), pointer :: mesh
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

end program sim_bsl_gc_2d0v_hex_splines
