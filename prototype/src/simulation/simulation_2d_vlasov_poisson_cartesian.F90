module sll_simulation_2d_vlasov_poisson_cartesian

! mimic of VP1D_deltaf_BSL in the simulation class
! intend to use MPI instead of openmp
! main application is for KEEN waves

!time mpirun -np 8 ./bin/test_2d_vp_cartesian keen

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"
  use sll_collective
  use sll_remapper
  use sll_buffer_loader_utilities_module
  use sll_constants
  !for mesh
  use sll_logical_meshes
  
  use sll_gnuplot_parallel
  !for initialization of distribution function
  !use sll_landau_2d_initializer
  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  use sll_parallel_array_initializer_module
  
  ! for interpolators
  !use sll_cubic_spline_interpolator_1d
  !use sll_periodic_interpolator_1d
  !use sll_odd_degree_spline_interpolator_1d
  use sll_module_advection_1d_periodic

  use sll_poisson_1d_periodic
  
  use sll_fft

  use sll_simulation_base
  implicit none

  type, extends(sll_simulation_base_class) :: &
       sll_simulation_2d_vlasov_poisson_cart


   !geometry
   type(sll_logical_mesh_2d), pointer :: mesh2d
      
   !interpolator
   !class(sll_interpolator_1d_base), pointer    :: interp_x
   !class(sll_interpolator_1d_base), pointer    :: interp_v

   !advector (should replace interpolator)
   class(sll_advection_1d_base), pointer    :: advect_x
   class(sll_advection_1d_base), pointer    :: advect_v


   
   !time_iterations
   sll_real64 :: dt
   sll_int32 :: first_step
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag,freq_diag_time
   sll_int32  :: split_case
  
   !parameters for initial function
   sll_real64  :: kx
   sll_real64  :: eps

   !parameters for drive
   sll_real64  :: t0, twL, twR, tstart, tflat, tL, tR, Edrmax, omegadr
   logical     :: driven, turn_drive_off
    
       
   contains
     procedure, pass(sim) :: run => run_vp2d_cartesian
     procedure, pass(sim) :: init_from_file => init_vp2d_par_cart
  end type sll_simulation_2d_vlasov_poisson_cart

  interface delete
     module procedure delete_vp2d_par_cart
  end interface delete

contains

  subroutine init_vp2d_par_cart( sim, filename )
    intrinsic :: trim
    class(sll_simulation_2d_vlasov_poisson_cart), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
    sll_int32             :: IO_stat
    sll_int32, parameter  :: input_file = 99
    sll_real64            :: xmin, vmin, vmax
    sll_int32             :: Ncx, nbox, Ncv
    sll_int32             :: interpol_x, order_x, interpol_v, order_v
    sll_real64            :: dt
    sll_int32             :: nbiter, freqdiag,freqdiag_time,first_step, split_case
    sll_real64            :: kmode, eps
    sll_int32             :: is_delta_f
    logical               :: driven
    sll_real64            :: t0, twL, twR, tstart, tflat, tL, tR, Edrmax, omegadr
    logical               :: turn_drive_off
    sll_int32, parameter  :: param_out = 37, param_out_drive = 40
    
    sll_real64            :: xmax
      
    ! namelists for data input
    namelist / geom / xmin, Ncx, nbox, vmin, vmax, Ncv
    namelist / interpolator / interpol_x, order_x, interpol_v, order_v
    namelist / time_iterations / dt, first_step, nbiter, freqdiag,freqdiag_time,split_case
    namelist / landau / kmode, eps, is_delta_f, driven 
    namelist / drive / t0, twL, twR, tstart, tflat, tL, tR, turn_drive_off, Edrmax, omegadr


    split_case = 1
    open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, '#init_vp2d_par_cart() failed to open file ', trim(filename)//'.nml'
       STOP
    end if
    read(input_file, geom) 
    read(input_file, interpolator)
    read(input_file, time_iterations)
    read(input_file, landau)
    if (driven) then
        read(input_file, drive)
        eps = 0.0_f64  ! no initial perturbation for driven simulation
    end if
    close(input_file)
    
    xmax = real(nbox,f64) * 2._f64 * sll_pi / kmode
    
    !we write the parameters in the simulation class
    sim%mesh2d => new_logical_mesh_2d(Ncx,Ncv,&
      eta1_min=xmin, eta1_max=xmax,&
      eta2_min=vmin, eta2_max=vmax)
        
    
!    select case (interpol_x)
!      case (1) ! periodic cubic spline
!        sim%interp_x => new_cubic_spline_1d_interpolator( Ncx + 1, xmin, xmax, SLL_PERIODIC)
!      case (2) ! arbitrary order periodic splines         
!         sim%interp_x => new_periodic_1d_interpolator( Ncx + 1, xmin, xmax, SPLINE, order_x)
!      case(3) ! arbitrary order Lagrange periodic interpolation
!         sim%interp_x => new_periodic_1d_interpolator( Ncx + 1, xmin, xmax, LAGRANGE, order_x)
!       case default
!         print*,'#interpolation in x number ', interpol_x, ' not implemented'
!         stop 
!    end select
!    select case (interpol_v)
!      case (1) ! hermite cubic spline
!       sim%interp_v => new_cubic_spline_1d_interpolator( Ncv + 1, vmin, vmax, SLL_HERMITE)
!      case (2) ! arbitrary order periodic splines
!        sim%interp_v => new_periodic_1d_interpolator( Ncv + 1, vmin, vmax, SPLINE, order_v)
!      case (3) ! arbitrary order Lagrange periodic interpolation
!        sim%interp_v => new_periodic_1d_interpolator( Ncv + 1, vmin, vmax, LAGRANGE, order_v)
!      !case(4) ! arbitrary order open spline interpolation   
!        !call sim%interp_comp_v%initialize( Ncv + 1, vmin, vmax, order_v)
!      case default
!        print*,'#interpolation in x number ', interpol_v, ' not implemented'
!        stop 
!    end select


    select case (interpol_x)
      case (2) ! arbitrary order periodic splines
        sim%advect_x => new_periodic_1d_advector( &
          Ncx, &
          xmin, &
          xmax, &
          SPLINE, & 
          order_x) 
      case(3) ! arbitrary order Lagrange periodic interpolation
        sim%advect_x => new_periodic_1d_advector( &
          Ncx, &
          xmin, &
          xmax, &
          LAGRANGE, & 
          order_x) 
       case default
         print*,'#interpolation in x number ', interpol_x, ' not implemented'
         stop 
    end select

    select case (interpol_v)
      case (2) ! arbitrary order periodic splines
        sim%advect_v => new_periodic_1d_advector( &
          Ncv, &
          vmin, &
          vmax, &
          SPLINE, & 
          order_v) 
      case(3) ! arbitrary order Lagrange periodic interpolation
        sim%advect_v => new_periodic_1d_advector( &
          Ncv, &
          vmin, &
          vmax, &
          LAGRANGE, & 
          order_v) 
       case default
         print*,'#interpolation in v number ', interpol_v, ' not implemented'
         stop 
    end select



    
     
    sim%dt=dt
    
    sim%first_step = first_step
    
    sim%num_iterations=nbiter
    sim%freq_diag=freqdiag
    sim%freq_diag_time=freqdiag_time
    sim%split_case = split_case
    sim%kx = kmode
    sim%eps = eps

    
    !sim%kmode=kmode
    !sim%eps=eps
    !sim%is_delta_f=is_delta_f
    
    sim%driven=driven
    if(driven) then
      sim%driven=driven
      sim%t0=t0
      sim%twL=twL
      sim%twR=twR
      sim%tstart=tstart
      sim%tflat=tflat
      sim%tL=tL
      sim%tR=tR
      sim%turn_drive_off=turn_drive_off
      sim%Edrmax=Edrmax
      sim%omegadr=omegadr
    endif
    
    if(sll_get_collective_rank(sll_world_collective)==0)then
      print *,'##geom'
      print *,'#xmin=',xmin
      print *,'#Ncx=',Ncx
      print *,'#Nbox=',Nbox
      print *,'#vmin=',vmin
      print *,'#vmax=',vmax
      print *,'#Ncv=',Ncv

      print *,'##interpolator'
      print *,'#interpol_x=',interpol_x
      print *,'#order_x=',order_x
      print *,'#interpol_v=',interpol_v
      print *,'#order_v=',order_v

      print *,'##time_iterations'
      print *,'#dt=',dt
      print *,'#first_step=',first_step
      print *,'#nbiter=',nbiter
      print *,'#freqdiag=',freqdiag
      print *,'#freqdiag_time=',freqdiag_time
      print *,'#split_case=',split_case

      print *,'##landau'
      print *,'#kmode=',kmode
      print *,'#eps=',eps
      print *,'#is_delta_f=',is_delta_f
      print *,'#driven=',driven

      if(driven)then
        print *,'##drive'
        print *,'#t0=',t0
        print *,'#twL=',twL
        print *,'#twR=',twR
        print *,'#tstart=',tstart
        print *,'#tflat=',tflat
        print *,'#tL=',tL
        print *,'#tR=',tR
        print *,'#turn_drive_off=',turn_drive_off
        print *,'#Edrmax=',Edrmax
        print *,'#omegadr=',omegadr
      endif 

!      open(unit = param_out, file = 'param_out.dat') 
!        write(param_out,'(A6,2f10.3,I5,2f10.3,I5,f10.3,I8,I5,I2,f10.3)') &
!         "landau", xmin, xmax, ncx, vmin, vmax, ncv, &
!         dt, nbiter, freqdiag, is_delta_f, kmode
!      close(param_out)
!
!      if (driven) then
!        open(unit = param_out_drive, file = 'param_out_drive.dat') 
!        write(param_out_drive,*) t0, twL, twR, tstart, tflat, tL, tR, &
!          Edrmax, omegadr
!        close(param_out_drive)
!      end if
      
      
      !to be compatible with VPpostprocessing_drive_KEEN.f
      if(driven) then     
      open(unit = param_out, file = 'parameters.dat')
      
      
        write(param_out, *) real(nbiter,f64)*dt !Tmax
        write(param_out, *) dt                  !dt
        write(param_out, *) nbiter              !Nt
        write(param_out, *) kmode               !k0
        write(param_out, *) omegadr             !omega0
        write(param_out, *) xmax-xmin           !L
        write(param_out, *) nbox                !mbox
        write(param_out, *) Edrmax              !Edrmax
        write(param_out, *) freqdiag_time       !nsave
        write(param_out, *) freqdiag            !nsavef1
        write(param_out, *) freqdiag            !nsavef2
        write(param_out, *) tR+tflat            !Tsetup
        write(param_out, *) vmax                !vxmax
        write(param_out, *) vmin                !vxmin
        write(param_out, *) Ncx                 !N
        write(param_out, *) Ncv                 !Nv
        ! added ES
        write(param_out, *) tL                  !tL
        write(param_out, *) tR                  !tR
        write(param_out, *) twL                 !twL
        write(param_out, *) twR                 !twR
        write(param_out, *) tflat               !tflat
      
      
      close(param_out)
      endif  
      

    endif
  end subroutine init_vp2d_par_cart


  subroutine run_vp2d_cartesian(sim)
    class(sll_simulation_2d_vlasov_poisson_cart), intent(inout) :: sim
    sll_real64,dimension(:,:),pointer :: f_x1,f_x2,f_x1_init!,buf_x1
    sll_real64,dimension(:),pointer :: rho,efield,e_app,rho_loc
    sll_int32, parameter  :: input_file = 33, th_diag = 34, ex_diag = 35, rho_diag = 36
    sll_int32, parameter  :: param_out = 37, eapp_diag = 38, adr_diag = 39
    sll_int32, parameter  :: param_out_drive = 40
    
    sll_int32 :: rhotot_id, efield_id, adr_id, Edr_id, deltaf_id,t_id,th_diag_id
    
    type(layout_2D), pointer :: layout_x1
    type(layout_2D), pointer :: layout_x2
    type(remap_plan_2D_real64), pointer :: remap_plan_x1_x2
    type(remap_plan_2D_real64), pointer :: remap_plan_x2_x1
    sll_real64, dimension(:), pointer     :: f1d
    sll_int32 :: np_x1,np_x2
    sll_int32 :: nproc_x1,nproc_x2
    sll_int32 :: global_indices(2)
    sll_int32 :: ierr
    sll_int32 :: local_size_x1,local_size_x2
    type(poisson_1d_periodic)  :: poisson_1d
    sll_real64 :: adr
    sll_real64::alpha
    sll_real64 ::tmp_loc(5),tmp(5)
    sll_int32  ::i,istep,ig,k
    
    sll_real64  ::   time, mass, momentum, l1norm, l2norm
    sll_real64  ::   kinetic_energy,potential_energy

    sll_real64, dimension(:), allocatable :: v_array
    sll_real64, dimension(:), allocatable :: v_array_middle
    sll_real64, dimension(:), allocatable :: x_array
    character(len=4)           :: fin   
    sll_int32                  :: file_id
    
    type(sll_fft_plan), pointer         :: pfwd
    sll_real64, dimension(:), allocatable :: buf_fft
    sll_comp64,dimension(:),allocatable :: rho_mode
    
    sll_int32 :: nb_mode = 5
    sll_real64, dimension(:), allocatable :: split_step
    sll_real64 :: t_step
    sll_int32 :: nb_split_step
    sll_int32 :: split_case 
    sll_int32 :: split_istep
    sll_int32 :: split_x
    sll_int32 :: split_x_init
    sll_int32, parameter :: SLL_STRANG_TVT           = 0 
    sll_int32, parameter :: SLL_STRANG_VTV           = 1 
    sll_int32, parameter :: SLL_TRIPLE_JUMP_TVT      = 2 
    sll_int32, parameter :: SLL_TRIPLE_JUMP_VTV      = 3 
    !sll_int32, parameter :: SLL_ORDER6_TVT           = 4 
    sll_int32, parameter :: SLL_ORDER6_VTV           = 5 
    sll_int32, parameter :: SLL_ORDER6VP_TVT         = 6 
    sll_int32, parameter :: SLL_ORDER6VP_VTV         = 7 
    sll_int32, parameter :: SLL_ORDER6VPnew_TVT      = 8 
    sll_int32, parameter :: SLL_ORDER6VPnew1_VTV     = 9 
    sll_int32, parameter :: SLL_ORDER6VPnew2_VTV     = 10 

    sll_int32 ::conservative_case
    
    
    ! for parallelization (output of distribution function in one single file)
    sll_int32, dimension(:), allocatable :: collective_displs
    sll_int32, dimension(:), allocatable :: collective_recvcnts
    sll_int32 :: collective_size
    sll_real64,dimension(:,:),pointer :: f_visu 
    sll_real64,dimension(:),pointer :: f_visu_buf1d
    sll_real64,dimension(:),pointer :: f_x1_buf1d

    np_x1 = sim%mesh2d%num_cells1+1
    np_x2 = sim%mesh2d%num_cells2+1

    if(sll_get_collective_rank(sll_world_collective)==0)then
      SLL_ALLOCATE(f_visu(np_x1,np_x2),ierr)
      SLL_ALLOCATE(f_visu_buf1d(np_x1*np_x2),ierr)
    else
      SLL_ALLOCATE(f_visu(1:0,1:0),ierr)          
      SLL_ALLOCATE(f_visu_buf1d(1:0),ierr)          
    endif

    collective_size = sll_get_collective_size(sll_world_collective)
    SLL_ALLOCATE(collective_displs(collective_size),ierr)
    SLL_ALLOCATE(collective_recvcnts(collective_size),ierr)

    

    !read from namelist file
    select case (sim%split_case)    
      case (0)
        split_case = SLL_STRANG_TVT
      case (1)
        split_case = SLL_STRANG_VTV
      case (2) 
        split_case = SLL_TRIPLE_JUMP_TVT
      case (3) 
        split_case = SLL_TRIPLE_JUMP_VTV
      case (4)
        split_case = SLL_ORDER6_VTV
      case (5)
        split_case = SLL_ORDER6VP_TVT
      case (6)
        split_case = SLL_ORDER6VP_VTV
      case (7)
        split_case = SLL_ORDER6VPnew_TVT
      case (8)
        split_case = SLL_ORDER6VPnew1_VTV
      case (9)
        split_case = SLL_ORDER6VPnew2_VTV      
    end select 

    !print *,'##split_case=',sim%split_case,split_case
    
    select case (split_case)    
      case (SLL_STRANG_TVT) ! Strang splitting TVT
        nb_split_step = 3
        SLL_ALLOCATE(split_step(nb_split_step),ierr)
        split_x_init = 1
        split_step(1) = 0.5_f64
        split_step(2) = 1._f64
        split_step(3) = split_step(1)
      case (SLL_STRANG_VTV) ! Strang splitting VTV
        nb_split_step = 3
        SLL_ALLOCATE(split_step(nb_split_step),ierr)
        split_x_init = 0
        split_step(1) = 0.5_f64
        split_step(2) = 1._f64
        split_step(3) = split_step(1)
      case (SLL_TRIPLE_JUMP_TVT) ! triple jump TVT
        nb_split_step = 7
        SLL_ALLOCATE(split_step(nb_split_step),ierr)
        split_x_init = 1
        split_step(1) = 0.675603595979829_f64
        split_step(2) = 1.351207191959658_f64
        split_step(3) = -0.17560359597982855_f64
        split_step(4) = -1.702414383919315_f64
        split_step(5) = split_step(3)
        split_step(6) = split_step(2)
        split_step(7) = split_step(1)
      case (SLL_TRIPLE_JUMP_VTV) ! triple jump VTV
        nb_split_step = 7
        SLL_ALLOCATE(split_step(nb_split_step),ierr)
        split_x_init = 0
        split_step(1) = 0.675603595979829_f64
        split_step(2) = 1.351207191959658_f64
        split_step(3) = -0.17560359597982855_f64
        split_step(4) = -1.702414383919315_f64
        split_step(5) = split_step(3)
        split_step(6) = split_step(2)
        split_step(7) = split_step(1)
      case (SLL_ORDER6_VTV) ! Order 6 VTV
        nb_split_step = 23
        SLL_ALLOCATE(split_step(nb_split_step),ierr)
        split_x_init = 0
        split_step(1) = 0.0414649985182624_f64
        split_step(2) = 0.123229775946271_f64
        split_step(3) = 0.198128671918067_f64
        split_step(4) = 0.290553797799558_f64
        split_step(5) = -0.0400061921041533_f64
        split_step(6) = -0.127049212625417_f64
        split_step(7) = 0.0752539843015807_f64          
        split_step(8) = -0.246331761062075_f64
        split_step(9) = -0.0115113874206879_f64
        split_step(10) = 0.357208872795928_f64
        split_step(11) = 0.23666992478693111_f64
        split_step(12) = 0.20477705429147008_f64
        split_step(13) = split_step(11)
        split_step(14) = split_step(10)
        split_step(15) = split_step(9)          
        split_step(16) = split_step(8)
        split_step(17) = split_step(7)
        split_step(18) = split_step(6)
        split_step(19) = split_step(5)
        split_step(20) = split_step(4)
        split_step(21) = split_step(3)
        split_step(22) = split_step(2)
        split_step(23) = split_step(1)          
      case (SLL_ORDER6VP_TVT) ! Order 6 for Vlasov-Poisson TVT 
        nb_split_step = 9
        SLL_ALLOCATE(split_step(nb_split_step),ierr)
        split_x_init = 1
        split_step(1) = 0.1095115577513980413559540_f64
        split_step(2) = 0.268722208204814693684441_f64&
          -2._f64*sim%dt**2*0.000805681667096178271312_f64&
          +4._f64*sim%dt**4*0.000017695766224036466792_f64
        split_step(3) = 0.4451715080955340951457244_f64
        split_step(4) = 0.2312777917951853063155588_f64&
          -2._f64*sim%dt**2*0.003955911930042478239763_f64&
          +4._f64*sim%dt**4*0.000052384078562246674986_f64
        split_step(5) = -0.1093661316938642730033570_f64
        split_step(6) = split_step(4)
        split_step(7) = split_step(3)
        split_step(8) = split_step(2)
        split_step(9) = split_step(1)
      case (SLL_ORDER6VP_VTV) ! Order 6 for Vlasov-Poisson VTV 
        nb_split_step = 9
        SLL_ALLOCATE(split_step(nb_split_step),ierr)
        split_x_init = 0
        split_step(1) = 0.359950808794143627485664_f64&
          -2._f64*sim%dt**2*(-0.01359558332625151635_f64)&
          +4._f64*sim%dt**4*(-8.562814848565929e-6_f64)
        split_step(2) = 1.079852426382430882456991_f64
        split_step(3) = -0.1437147273026540434771131_f64&
          -2._f64*sim%dt**2*(-0.00385637757897273261_f64)&
          +4._f64*sim%dt**4*(0.0004883788785819335822_f64)
        split_step(4) = -0.579852426382430882456991_f64
        split_step(5) = 0.567527837017020831982899_f64&
          -2._f64*sim%dt**2*(-0.03227361602037480885_f64)&
          +4._f64*sim%dt**4*0.002005141087312622342_f64
        split_step(6) = split_step(4)
        split_step(7) = split_step(3)
        split_step(8) = split_step(2)
        split_step(9) = split_step(1)
      case (SLL_ORDER6VPnew_TVT) ! Order 6 for Vlasov-Poisson TVT (new)
        nb_split_step = 9
        SLL_ALLOCATE(split_step(nb_split_step),ierr)
        split_x_init = 1
        split_step(1) = 0.1095115577513980413559540_f64
        split_step(2) = 0.268722208204814693684441_f64&
          -2._f64*sim%dt**2*0.000805681667096178271312_f64&
          +4._f64*sim%dt**4*(8.643923349886021963e-6_f64)&
          -8._f64*sim%dt**6*(1.4231479258353431522e-6_f64)
        split_step(3) = 0.4451715080955340951457244_f64
        split_step(4) = 0.2312777917951853063155588_f64&
          -2._f64*sim%dt**2*0.003955911930042478239763_f64&
          +4._f64*sim%dt**4*(0.000061435921436397119815_f64)
        split_step(5) = -0.1093661316938642730033570_f64
        split_step(6) = split_step(4)
        split_step(7) = split_step(3)
        split_step(8) = split_step(2)
        split_step(9) = split_step(1)
      case (SLL_ORDER6VPnew1_VTV) ! Order 6 for Vlasov-Poisson VTV (new1)
        nb_split_step = 11
        SLL_ALLOCATE(split_step(nb_split_step),ierr)
        split_x_init = 0
        split_step(1) = 0.0490864609761162454914412_f64&
          -2._f64*sim%dt**2*(0.0000697287150553050840999_f64)
        split_step(2) = 0.1687359505634374224481957_f64
        split_step(3) = 0.2641776098889767002001462_f64&
          -2._f64*sim%dt**2*(0.000625704827430047189169_f64)&
          +4._f64*sim%dt**4*(-2.91660045768984781644e-6_f64)
        split_step(4) = 0.377851589220928303880766_f64
        split_step(5) = 0.1867359291349070543084126_f64&
          -2._f64*sim%dt**2*(0.00221308512404532556163_f64)&
          +4._f64*sim%dt**4*(0.0000304848026170003878868_f64)&
          -8._f64*sim%dt**6*(4.98554938787506812159e-7_f64)
        split_step(6) = -0.0931750795687314526579244_f64
        split_step(7) = split_step(5)
        split_step(8) = split_step(4)
        split_step(9) = split_step(3)
        split_step(10) = split_step(2)
        split_step(11) = split_step(1)
      case (SLL_ORDER6VPnew2_VTV) ! Order 6 for Vlasov-Poisson VTV (new2)
        nb_split_step = 11
        SLL_ALLOCATE(split_step(nb_split_step),ierr)
        split_x_init = 0
        split_step(1) = 0.083335463273305120964507_f64&
          -2._f64*sim%dt**2*(-0.00015280483587048489661_f64)&
          +4._f64*sim%dt**4*( -0.0017675734111895638156_f64)&
          -8._f64*sim%dt**6*( 0.00021214072262165668039_f64)
        split_step(2) = 0.72431592569108212422250_f64
        split_step(3) = 0.827694857845135145869413_f64&
          -2._f64*sim%dt**2*(-0.010726848627286273332_f64)&
          +4._f64*sim%dt**4*(0.012324362982853212700_f64)
        split_step(4) = -0.4493507217041624582458844_f64
        split_step(5) = -0.4110303211184402668339201_f64&
          -2._f64*sim%dt**2*(0.014962337009932798678_f64)
        !  +4._f64*sim%dt**4*(-0.20213028317750837901_f64)
        split_step(6) = 0.4500695920261606680467717_f64
        split_step(7) = split_step(5)
        split_step(8) = split_step(4)
        split_step(9) = split_step(3)
        split_step(10) = split_step(2)
        split_step(11) = split_step(1)
      case default
        print *,'#split_case not defined'
        stop       
    end select
    
    
    
    if(sll_get_collective_rank(sll_world_collective)==0)then
      SLL_ALLOCATE(buf_fft(np_x1-1),ierr)
      pfwd => fft_new_plan(np_x1-1,buf_fft,buf_fft,FFT_FORWARD,FFT_NORMALIZE)
      SLL_ALLOCATE(rho_mode(0:nb_mode),ierr)      
    endif

    ! allocate and initialize the layouts...
    layout_x1       => new_layout_2D( sll_world_collective )
    layout_x2       => new_layout_2D( sll_world_collective )    
    nproc_x1 = sll_get_collective_size( sll_world_collective )
    nproc_x2 = 1
    call initialize_layout_with_distributed_2D_array( &
      np_x1, np_x2, nproc_x1, nproc_x2, layout_x2 )
    call initialize_layout_with_distributed_2D_array( &
      np_x1, np_x2, nproc_x2, nproc_x1, layout_x1 )
    !call sll_view_lims_2D( layout_x1 )

    
    !allocation of distribution functions f_x1 and f_x2
    call compute_local_sizes_2d( layout_x2, local_size_x1, local_size_x2 )
    SLL_ALLOCATE(f_x2(local_size_x1,local_size_x2),ierr)

    call compute_local_sizes_2d( layout_x1, local_size_x1, local_size_x2 )
    global_indices(1:2) = local_to_global_2D( layout_x1, (/1, 1/) )
    SLL_ALLOCATE(f_x1(local_size_x1,local_size_x2),ierr)    
    SLL_ALLOCATE(f_x1_init(local_size_x1,local_size_x2),ierr)    
    SLL_ALLOCATE(f_x1_buf1d(local_size_x1*local_size_x2),ierr)    


    !definition of remap
    remap_plan_x1_x2 => NEW_REMAP_PLAN(layout_x1, layout_x2, f_x1)
    remap_plan_x2_x1 => NEW_REMAP_PLAN(layout_x2, layout_x1, f_x2)

    
    !allocation of 1d arrays
    SLL_ALLOCATE(rho(np_x1),ierr)
    SLL_ALLOCATE(rho_loc(np_x1),ierr)
    SLL_ALLOCATE(efield(np_x1),ierr)
    SLL_ALLOCATE(e_app(np_x1),ierr)
    SLL_ALLOCATE(f1d(max(np_x1,np_x2)),ierr)
    SLL_ALLOCATE(v_array(np_x2),ierr)
    SLL_ALLOCATE(x_array(np_x1),ierr)
    SLL_ALLOCATE(v_array_middle(np_x2),ierr)


    do i = 1, np_x2
       v_array(i) = sim%mesh2d%eta2_min + real(i-1,f64)*sim%mesh2d%delta_eta2
    end do
    do i = 1, np_x1
       x_array(i) = sim%mesh2d%eta1_min + real(i-1,f64)*sim%mesh2d%delta_eta1
    end do

    do i = 1, np_x2-1
       v_array_middle(i) = 0.5_f64*(v_array(i)+v_array(i+1))
    end do
    v_array_middle(np_x2) = v_array_middle(1)+v_array(np_x2)-v_array(1)


    if(sim%driven)then
    if(sll_get_collective_rank(sll_world_collective)==0)then
  !open(unit=11, file='x.bdat', &
  !     form='unformatted')
  !write(11) x_array
  !close(11)
      call sll_binary_file_create("x.bdat", file_id, ierr)
      call sll_binary_write_array_1d(file_id,x_array(1:np_x1-1),ierr)
      call sll_binary_file_close(file_id,ierr)                    
      call sll_binary_file_create("v.bdat", file_id, ierr)
      call sll_binary_write_array_1d(file_id,v_array(1:np_x2-1),ierr)
      call sll_binary_file_close(file_id,ierr)                                             
    endif
    endif
    
    !initialize distribution function     

    call sll_2d_parallel_array_initializer_cartesian( &
       layout_x1, &
       sim%mesh2d, &
       f_x1_init, &
       sll_landau_initializer_2d, &
       (/ sim%kx,0._f64 /))
    if(sim%first_step/sim%freq_diag==0)then
      call sll_2d_parallel_array_initializer_cartesian( &
       layout_x1, &
       sim%mesh2d, &
       f_x1, &
       sll_landau_initializer_2d, &
       (/ sim%kx,sim%eps /))
       sim%first_step = 0
       istep = 0
    else
      print *,'#restart strategy not implemented yet'
      stop
    endif
    
    call compute_displacements_array_2d( &
      layout_x1, &
      collective_size, &
      collective_displs )
    collective_recvcnts = receive_counts_array_2d( &
      layout_x1, &
      collective_size )

    call load_buffer_2d( layout_x1, f_x1, f_x1_buf1d )
    call sll_collective_gatherv_real64( &
      sll_world_collective, &
      f_x1_buf1d, &
      sll_get_collective_rank(sll_world_collective), &
      collective_recvcnts, &
      collective_displs, &
      0, &
      f_visu_buf1d )
    f_visu = reshape(f_visu_buf1d, shape(f_visu))
    if(sll_get_collective_rank(sll_world_collective)==0)then
      !print *,'#begin f0.bdat'                    
      call sll_binary_file_create('f0.bdat', file_id, ierr)
      call sll_binary_write_array_2d(file_id,f_visu(1:np_x1-1,1:np_x2-1),ierr)
      call sll_binary_file_close(file_id,ierr)
      !print *,'#finish f0.bdat'
      !stop                    
    endif
    
    
    
    !call new(poisson_1d,sim%mesh2d%eta1_min,sim%mesh2d%eta1_max,np_x1-1,ierr)    
    call initialize(poisson_1d,sim%mesh2d%eta1_min,sim%mesh2d%eta1_max,np_x1-1,ierr)

    !computation of electric field
    rho_loc = 0._f64
    do i=1,np_x1
      rho_loc(i)=rho_loc(i)+sum(f_x1(i,1:local_size_x2))
    end do    
    call mpi_barrier(sll_world_collective%comm,ierr)
    call mpi_allreduce(rho_loc,rho,np_x1,MPI_REAL8,MPI_SUM,sll_world_collective%comm,ierr)
    rho = 1._f64-sim%mesh2d%delta_eta2*rho        
    call solve(poisson_1d, efield, rho)    
    ! Ponderomotive force at initial time. We use a sine wave
    ! with parameters k_dr and omega_dr.
    istep = sim%first_step
    e_app = 0._f64
    if (sim%driven) then
      call PFenvelope(adr, istep*sim%dt, sim%tflat, sim%tL, sim%tR, sim%twL, sim%twR, &
          sim%t0, sim%turn_drive_off)
      do i = 1, np_x1
        e_app(i) = sim%Edrmax*adr*sim%kx&
          *sin(sim%kx*real(i-1,f64)*sim%mesh2d%delta_eta1&
          -sim%omegadr*real(istep,f64)*sim%dt)
      enddo
    endif

    ! write initial fields
    if(sll_get_collective_rank(sll_world_collective)==0)then
      call sll_ascii_file_create('thdiag.dat', th_diag_id, ierr)
      call sll_binary_file_create('deltaf.bdat', deltaf_id, ierr)
      if(sim%driven)then
        call sll_binary_file_create('rhotot.bdat', rhotot_id, ierr)
        call sll_binary_file_create('efield.bdat', efield_id, ierr)
        call sll_binary_file_create('adr.bdat', adr_id, ierr)
        call sll_binary_file_create('Edr.bdat', Edr_id, ierr)
        call sll_binary_file_create('t.bdat', t_id, ierr)
        call sll_binary_write_array_1d(efield_id,efield(1:np_x1-1),ierr)
        call sll_binary_write_array_1d(rhotot_id,rho(1:np_x1-1),ierr)
        call sll_binary_write_array_1d(Edr_id,e_app(1:np_x1-1),ierr)
        call sll_binary_write_array_0d(adr_id,adr,ierr)
        call sll_binary_write_array_0d(t_id,real(istep,f64)*sim%dt,ierr)
      endif                    
    endif
    
    
    !write also initial deltaf function
    call load_buffer_2d( layout_x1, f_x1-f_x1_init, f_x1_buf1d )
    call sll_collective_gatherv_real64( &
      sll_world_collective, &
      f_x1_buf1d, &
      sll_get_collective_rank(sll_world_collective), &
      collective_recvcnts, &
      collective_displs, &
      0, &
      f_visu_buf1d )
    f_visu = reshape(f_visu_buf1d, shape(f_visu))


    if(sll_get_collective_rank(sll_world_collective)==0)then
      !call sll_binary_file_create('f0.bdat', file_id, ierr)
      call sll_binary_write_array_2d(deltaf_id,f_visu(1:np_x1-1,1:np_x2-1),ierr)
    endif




    

    do istep = sim%first_step+1, sim%num_iterations


      split_x = split_x_init
      t_step = real(istep-1,f64)
      do split_istep=1,nb_split_step
        
        if(split_x==1)then
          !! T ADVECTION 
          !advection in x
          do i = 1, local_size_x2
            ig=global_indices(2)
            alpha = (sim%mesh2d%eta2_min + real(i+ig-2,f64) * sim%mesh2d%delta_eta2) &
              * split_step(split_istep) ! *sim%dt
            !print *,'alpha=',alpha,i
            f1d(1:np_x1) = f_x1(1:np_x1,i)
            
            !f1d(1:np_x1) = sim%interp_x%interpolate_array_disp(np_x1, f1d(1:np_x1), alpha)
            
            call sim%advect_x%advect_1d_constant(&
              alpha, &
              sim%dt, &
              f1d(1:np_x1), &
              f1d(1:np_x1))
            
            f_x1(1:np_x1,i)=f1d(1:np_x1)
            !print *,'ok'
          end do
          t_step = t_step+split_step(split_istep)
          !computation of electric field
          rho_loc = 0._f64
          do i=1,np_x1
            rho_loc(i)=rho_loc(i)+sum(f_x1(i,1:local_size_x2))
            !call mpi_reduce(tmp(1),rho(i),1,MPI_REAL8,MPI_SUM,0,sll_world_collective%comm,ierr)
          end do    
          call mpi_barrier(sll_world_collective%comm,ierr)
          call mpi_allreduce(rho_loc,rho,np_x1,MPI_REAL8,MPI_SUM,sll_world_collective%comm,ierr)
          rho = 1._f64-sim%mesh2d%delta_eta2*rho
          call solve(poisson_1d, efield, rho)
          if (sim%driven) then
            call PFenvelope(adr, t_step*sim%dt, sim%tflat, sim%tL, sim%tR, sim%twL, sim%twR, &
              sim%t0, sim%turn_drive_off)
            do i = 1, np_x1
              e_app(i) = sim%Edrmax*adr*sim%kx&
              *sin(sim%kx*real(i-1,f64)*sim%mesh2d%delta_eta1&
              -sim%omegadr*t_step*sim%dt)
            enddo
          endif
        
        else
          !! V ADVECTION 
          !transposition
          call apply_remap_2D( remap_plan_x1_x2, f_x1, f_x2 )
          call compute_local_sizes_2d( layout_x2, local_size_x1, local_size_x2 )
          global_indices(1:2) = local_to_global_2D( layout_x2, (/1, 1/) )
          !advection in v
          do i = 1,local_size_x1
            ig=i+global_indices(1)-1
            alpha = -(efield(ig)+e_app(ig)) * split_step(split_istep) !* sim%dt
            !print *,'alpha=',alpha,i
            f1d(1:np_x2) = f_x2(i,1:np_x2)
                        
            !f1d(1:np_x2) = sim%interp_v%interpolate_array_disp(np_x2, f1d(1:np_x2), alpha)
            
            call sim%advect_v%advect_1d_constant(&
              alpha, &
              sim%dt, &
              f1d(1:np_x2), &
              f1d(1:np_x2))

            
            f_x2(i,1:np_x2) = f1d(1:np_x2)
          end do
          !transposition
          call apply_remap_2D( remap_plan_x2_x1, f_x2, f_x1 )
          call compute_local_sizes_2d( layout_x1, local_size_x1, local_size_x2 )
          global_indices(1:2) = local_to_global_2D( layout_x1, (/1, 1/) )

        endif
        split_x= 1-split_x
      
      enddo
      

     !!DIAGNOSTICS
     if (mod(istep,sim%freq_diag_time)==0) then
        time = real(istep,f64)*sim%dt
        mass = 0._f64
        momentum = 0._f64
        l1norm = 0._f64
        l2norm = 0._f64
        kinetic_energy = 0._f64
        potential_energy = 0._f64
        tmp_loc = 0._f64        
        do i = 1, np_x1-1        
          tmp_loc(1)= tmp_loc(1)+sum(f_x1(i,1:local_size_x2))
          tmp_loc(2)= tmp_loc(2)+sum(abs(f_x1(i,1:local_size_x2)))
          tmp_loc(3)= tmp_loc(3)+sum((f_x1(i,1:local_size_x2))**2)
          tmp_loc(4)= tmp_loc(4)+sum(f_x1(i,1:local_size_x2)*v_array(1:local_size_x2))
          tmp_loc(5)= tmp_loc(5)+sum(f_x1(i,1:local_size_x2)*v_array(1:local_size_x2)**2)          
        end do
        call mpi_barrier(sll_world_collective%comm,ierr)
        call mpi_allreduce(tmp_loc,tmp,5,MPI_REAL8,MPI_SUM,sll_world_collective%comm,ierr)
        mass = tmp(1) * sim%mesh2d%delta_eta1 * sim%mesh2d%delta_eta2
        l1norm = tmp(2)  * sim%mesh2d%delta_eta1 * sim%mesh2d%delta_eta2
        l2norm = tmp(3)  * sim%mesh2d%delta_eta1 * sim%mesh2d%delta_eta2
        momentum = tmp(4) * sim%mesh2d%delta_eta1 * sim%mesh2d%delta_eta2
        kinetic_energy = 0.5_f64 *tmp(5) * sim%mesh2d%delta_eta1 * sim%mesh2d%delta_eta2
        potential_energy = 0._f64
        do i=1, np_x1-1
          potential_energy = potential_energy+(efield(i)+e_app(i))**2
          !  0.5_f64 * sum((efield(1:np_x1-1)+e_app(1:np_x1-1))**2) * sim%mesh2d%delta_eta1
        enddo
        potential_energy = 0.5_f64*potential_energy* sim%mesh2d%delta_eta1
        if(sll_get_collective_rank(sll_world_collective)==0)then                  
          buf_fft = rho(1:np_x1-1)
          call fft_apply_plan(pfwd,buf_fft,buf_fft)
          do k=0,nb_mode
            rho_mode(k)=fft_get_mode(pfwd,buf_fft,k)
          enddo  
          write(th_diag_id,'(f12.5,13g20.12)') time, mass, l1norm, momentum, l2norm, &
             kinetic_energy, potential_energy, kinetic_energy + potential_energy, &
             abs(rho_mode(0)),abs(rho_mode(1)),abs(rho_mode(2)),abs(rho_mode(3)), &
             abs(rho_mode(4)),abs(rho_mode(5))
          if(sim%driven)then
            call sll_binary_write_array_1d(efield_id,efield(1:np_x1-1),ierr)
            call sll_binary_write_array_1d(rhotot_id,rho(1:np_x1-1),ierr)
            call sll_binary_write_array_1d(Edr_id,e_app(1:np_x1-1),ierr)
            call sll_binary_write_array_0d(adr_id,adr,ierr)
            call sll_binary_write_array_0d(t_id,real(istep,f64)*sim%dt,ierr)
          endif   
        endif
          
        if (mod(istep,sim%freq_diag)==0) then          
          !we substract f0
          !f_x1 = f_x1-f_x1_init                    
          !we gather in one file
          call load_buffer_2d( layout_x1, f_x1-f_x1_init, f_x1_buf1d )
          call sll_collective_gatherv_real64( &
            sll_world_collective, &
            f_x1_buf1d, &
            sll_get_collective_rank(sll_world_collective), &
            collective_recvcnts, &
            collective_displs, &
            0, &
            f_visu_buf1d )
          f_visu = reshape(f_visu_buf1d, shape(f_visu))
          if(sll_get_collective_rank(sll_world_collective)==0) then
            call sll_binary_write_array_2d(deltaf_id,f_visu(1:np_x1-1,1:np_x2-1),ierr)  
          endif
        endif
          

     end if


    
    
    enddo
    
    ! close files
    if(sll_get_collective_rank(sll_world_collective)==0)then
      call sll_ascii_file_close(th_diag_id,ierr) 
      call sll_binary_file_close(deltaf_id,ierr) 
      if(sim%driven)then
        call sll_binary_file_close(efield_id,ierr)
        call sll_binary_file_close(rhotot_id,ierr)
        call sll_binary_file_close(Edr_id,ierr)
        call sll_binary_file_close(adr_id,ierr)
        call sll_binary_file_close(t_id,ierr)
      endif   
    endif


    
    
  end subroutine run_vp2d_cartesian

  subroutine delete_vp2d_par_cart( sim )
    class(sll_simulation_2d_vlasov_poisson_cart) :: sim
    sll_int32 :: ierr
  end subroutine delete_vp2d_par_cart


  elemental function f_equilibrium(v)
    sll_real64, intent(in) :: v
    sll_real64 :: f_equilibrium

    f_equilibrium = 1.0_f64/sqrt(2*sll_pi)*exp(-0.5_f64*v*v)
  end function f_equilibrium

  subroutine PFenvelope(S, t, tflat, tL, tR, twL, twR, t0, &
       turn_drive_off)

    ! DESCRIPTION
    ! -----------
    ! S: the wave form at a given point in time. This wave form is 
    !    not scaled (its maximum value is 1).
    ! t: the time at which the envelope is being evaluated
    ! tflat, tL, tR, twL, twR, tstart, t0: the parameters defining the
    !    envelope, defined in the main portion of this program.
    ! turn_drive_off: 1 if the drive should be turned off after a time
    !    tflat, and 0 otherwise

    sll_real64, intent(in) :: t, tflat, tL, tR, twL, twR, t0
    sll_real64, intent(out) :: S
    logical, intent(in) :: turn_drive_off
    ! local variables
    sll_int32 :: i 
    sll_real64 :: epsilon

    ! The envelope function is defined such that it is zero at t0,
    ! rises to 1 smoothly, stay constant for tflat, and returns
    ! smoothly to zero.
    if(turn_drive_off) then
       epsilon = 0.5*(tanh((t0-tL)/twL) - tanh((t0-tR)/twR))
       S = 0.5*(tanh((t-tL)/twL) - tanh((t-tR)/twR)) - epsilon
       S = S / (1-epsilon)
    else
       epsilon = 0.5*(tanh((t0-tL)/twL) + 1)
       S = 0.5*(tanh((t-tL)/twL) + 1) - epsilon
       S = S / (1-epsilon)
    endif
    if(S<0) then
       S = 0.
    endif
    return
  end subroutine PFenvelope



end module sll_simulation_2d_vlasov_poisson_cartesian
