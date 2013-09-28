module sll_simulation_2d_vlasov_poisson_cartesian

! mimic of VP1D_deltaf_BSL in the simulation class
! intend to use MPI instead of openmp
! main application is for KEEN waves


#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"
  use sll_collective
  use sll_remapper
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
  use sll_cubic_spline_interpolator_1d
  use sll_periodic_interpolator_1d
  use sll_odd_degree_spline_interpolator_1d

  use sll_poisson_1d_periodic

  use sll_simulation_base
  implicit none

  type, extends(sll_simulation_base_class) :: &
       sll_simulation_2d_vlasov_poisson_cart

    ! Parallel environment parameters
    sll_int32  :: world_size
    sll_int32  :: my_rank
    sll_int32  :: power2 ! 2^power2 = number of processes available

   !geometry
   type(sll_logical_mesh_2d), pointer :: mesh2d
      
   !interpolator
   class(sll_interpolator_1d_base), pointer    :: interp_x
   class(sll_interpolator_1d_base), pointer    :: interp_v
   !type(cubic_spline_1d_interpolator):: interp_spline_x, interp_spline_v
    !type(cubic_spline_1d_interpolator)         :: interp_spline_x
    !type(cubic_spline_1d_interpolator)         :: interp_spline_v
    !type(per_1d_interpolator)                  :: interp_per_x
    !type(per_1d_interpolator)                  :: interp_per_v
    !type(odd_degree_spline_1d_interpolator)    :: interp_comp_v
   
   !time_iterations
   sll_real64 :: dt
   sll_int32  :: num_iterations
   sll_int32  :: freq_diag
   
   !parameters for initial function
   sll_real64  :: kx
   sll_real64  :: eps

   !parameters for drive
   sll_real64  :: t0, twL, twR, tstart, tflat, tL, tR, Edrmax, omegadr
   logical     :: driven, turn_drive_off

   
   
!    namelist / interpolator / interpol_x, order_x, interpol_v, order_v
!    namelist / time_iterations / dt, nbiter, freqdiag
!    namelist / landau / kmode, eps, is_delta_f, driven 
!    namelist / drive / t0, twL, twR, tstart, tflat, tL, tR, turn_drive_off, Edrmax, omegadr
    
       
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
    sll_int32             :: nbiter, freqdiag
    sll_real64            :: kmode, eps
    sll_int32             :: is_delta_f
    logical               :: driven
    sll_real64            :: t0, twL, twR, tstart, tflat, tL, tR, Edrmax, omegadr
    logical               :: turn_drive_off
    sll_int32, parameter  :: param_out = 37, param_out_drive = 40
    
    sll_real64            :: xmax
    !class(sll_interpolator_1d_base), pointer    :: interp_x, interp_v

    !type(init_landau_2d), target :: init_landau
    !class(scalar_field_2d_initializer_base), pointer    :: p_init_f
    class(sll_coordinate_transformation_2d_base), pointer :: tx



      
    ! namelists for data input
    namelist / geom / xmin, Ncx, nbox, vmin, vmax, Ncv
    namelist / interpolator / interpol_x, order_x, interpol_v, order_v
    namelist / time_iterations / dt, nbiter, freqdiag
    namelist / landau / kmode, eps, is_delta_f, driven 
    namelist / drive / t0, twL, twR, tstart, tflat, tL, tR, turn_drive_off, Edrmax, omegadr

    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, '#init_vp2d_par_cart() failed to open file ', filename
       STOP
    end if
    read(input_file, geom) 
    read(input_file, interpolator)
    read(input_file, time_iterations)
    read(input_file, landau)
    if (driven) then
        read(input_file, drive)
        eps = 0.0  ! no initial perturbation for driven simulation
    end if
    close(input_file)
    
    xmax = nbox * 2._f64 * sll_pi / kmode
    
    !we write the parameters in the simulation class
    sim%mesh2d => new_logical_mesh_2d(Ncx,Ncv,&
      eta1_min=xmin, eta1_max=xmax,&
      eta2_min=vmin, eta2_max=vmax)
    

    
    
    !call init_landau%initialize(mesh2d_base, NODE_CENTERED_FIELD, eps, kmode, &
    !  is_delta_f)
    !p_init_f =>init_landau
    
    
    select case (interpol_x)
      case (1) ! periodic cubic spline
        sim%interp_x => new_cubic_spline_1d_interpolator( Ncx + 1, xmin, xmax, SLL_PERIODIC)
        !call sim%interp_spline_x%initialize( Ncx + 1, xmin, xmax, SLL_PERIODIC)
        !sim%interp_x => sim%interp_spline_x
       case (2) ! arbitrary order periodic splines
         !call sim%interp_per_x%initialize( Ncx + 1, xmin, xmax, SPLINE, order_x)
         !sim%interp_x => sim%interp_per_x
       case(3) ! arbitrary order Lagrange periodic interpolation
         !call sim%interp_per_x%initialize( Ncx + 1, xmin, xmax, LAGRANGE, order_x)
         !sim%interp_x => sim%interp_per_x
       case default
         print*,'#interpolation in x number ', interpol_x, ' not implemented'
         stop 
    end select
    select case (interpol_v)
      case (1) ! hermite cubic spline
       sim%interp_v => new_cubic_spline_1d_interpolator( Ncv + 1, vmin, vmax, SLL_HERMITE)
       !call sim%interp_spline_v%initialize( Ncv + 1, vmin, vmax, SLL_HERMITE)
       !sim%interp_v => sim%interp_spline_v
      case (2) ! arbitrary order periodic splines
        !call sim%interp_per_v%initialize( Ncv + 1, vmin, vmax, SPLINE, order_v)
        !sim%interp_v => sim%interp_per_v
      case (3) ! arbitrary order Lagrange periodic interpolation
        !call sim%interp_per_v%initialize( Ncv + 1, vmin, vmax, LAGRANGE, order_v)
        !sim%interp_v => sim%interp_per_v
      case(4) ! arbitrary order open spline interpolation   
        !call sim%interp_comp_v%initialize( Ncv + 1, vmin, vmax, order_v)
      case default
        print*,'#interpolation in x number ', interpol_v, ' not implemented'
        stop 
    end select
    
     
    sim%dt=dt
    sim%num_iterations=nbiter
    sim%freq_diag=freqdiag

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
      print *,'#nbiter=',nbiter
      print *,'#freqdiag=',freqdiag

      print *,'##landau'
      print *,'#kmode=',kmode
      print *,'#eps=',eps
      print *,'#is_delta_f=',is_delta_f
      print *,'#driven=',driven


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
    

      open(unit = param_out, file = 'param_out.dat') 
        write(param_out,'(A6,2f10.3,I5,2f10.3,I5,f10.3,I8,I5,I2,f10.3)') &
         "landau", xmin, xmax, ncx, vmin, vmax, ncv, &
         dt, nbiter, freqdiag, is_delta_f, kmode
      close(param_out)

      if (driven) then
        open(unit = param_out_drive, file = 'param_out_drive.dat') 
        write(param_out_drive,*) t0, twL, twR, tstart, tflat, tL, tR, &
          Edrmax, omegadr
        close(param_out_drive)
      end if

    endif
  end subroutine init_vp2d_par_cart


  subroutine run_vp2d_cartesian(sim)
    class(sll_simulation_2d_vlasov_poisson_cart), intent(inout) :: sim
    sll_real64,dimension(:,:),pointer :: f_x1,f_x2 
    sll_real64,dimension(:),pointer :: rho,efield,e_app
    sll_int32, parameter  :: input_file = 33, th_diag = 34, ex_diag = 35, rho_diag = 36
    sll_int32, parameter  :: param_out = 37, eapp_diag = 38, adr_diag = 39
    sll_int32, parameter  :: param_out_drive = 40
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
    sll_real64 ::tmp
    sll_int32  ::i,istep,ig

    !interp_x=sim%interp_x
    !spline_x=sim%interp_spline_x
    
    !interp_x=>spline_x
    
    np_x1 = sim%mesh2d%num_cells1+1
    np_x2 = sim%mesh2d%num_cells2+1

    !call spline_x%initialize( np_x1, sim%mesh2d%eta1_min, &
    !  sim%mesh2d%eta1_max, SLL_PERIODIC)
    !call spline_v%initialize( np_x2, sim%mesh2d%eta2_min, &
    !  sim%mesh2d%eta2_max, SLL_PERIODIC)


    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)

    ! allocate the layouts...
    layout_x1       => new_layout_2D( sll_world_collective )
    layout_x2       => new_layout_2D( sll_world_collective )

    ! Make an even split for x1 and x2, or as close as even if the power of 2 is odd. 
    !This should be packaged in some sort of routine and set up at initialization time.

    nproc_x1 = sll_get_collective_size( sll_world_collective )
    nproc_x2 = 1
    
    
    
    call initialize_layout_with_distributed_2D_array( &
      np_x1, np_x2, nproc_x1, nproc_x2, layout_x2 )
    call initialize_layout_with_distributed_2D_array( &
      np_x1, np_x2, nproc_x2, nproc_x1, layout_x1 )

    !call sll_view_lims_2D( layout_x1 )

    call compute_local_sizes_2d( layout_x1, local_size_x1, local_size_x2 )
    SLL_ALLOCATE(f_x1(local_size_x1,local_size_x2),ierr)
    
    call compute_local_sizes_2d( layout_x2, local_size_x1, local_size_x2 )
    SLL_ALLOCATE(f_x2(local_size_x1,local_size_x2),ierr)

    SLL_ALLOCATE(rho(np_x1),ierr)
    SLL_ALLOCATE(efield(np_x1),ierr)
    SLL_ALLOCATE(e_app(np_x1),ierr)
    SLL_ALLOCATE(f1d(max(np_x1,np_x2)),ierr)
    
    print *,'#N=',max(np_x1,np_x2)
    
    
    !initialize distribution function
    
    call sll_2d_parallel_array_initializer_cartesian( &
       layout_x2, &
       sim%mesh2d, &
       f_x2, &
       sll_landau_initializer_2d, &
       (/ sim%kx,sim%eps /))
    
    global_indices(1:2) = local_to_global_2D( layout_x2, (/1, 1/) )
    
    call sll_gnuplot_rect_2d_parallel( &
      sim%mesh2d%eta1_min+(global_indices(1)-1)*sim%mesh2d%delta_eta1, &
      sim%mesh2d%delta_eta1, &
      sim%mesh2d%eta2_min+(global_indices(2)-1)*sim%mesh2d%delta_eta2, &
      sim%mesh2d%delta_eta2, &
      f_x2, &
      "f_x2", &
      0, &
      ierr )
     
    remap_plan_x1_x2 => NEW_REMAP_PLAN(layout_x1, layout_x2, f_x1)
    remap_plan_x2_x1 => NEW_REMAP_PLAN(layout_x2, layout_x1, f_x2)


    call apply_remap_2D( remap_plan_x2_x1, f_x2, f_x1 )
    call compute_local_sizes_2d( layout_x1, local_size_x1, local_size_x2 )
    
    
    global_indices(1:2) = local_to_global_2D( layout_x1, (/1, 1/) )
    call sll_gnuplot_rect_2d_parallel( &
      sim%mesh2d%eta1_min+(global_indices(1)-1)*sim%mesh2d%delta_eta1, &
      sim%mesh2d%delta_eta1, &
      sim%mesh2d%eta2_min+(global_indices(2)-1)*sim%mesh2d%delta_eta2, &
      sim%mesh2d%delta_eta2, &
      f_x1, &
      "f_x1", &
      0, &
      ierr )

    
    !computation of rho
    call new(poisson_1d,sim%mesh2d%eta1_min,sim%mesh2d%eta1_max,np_x1-1,ierr)
    do i=1,np_x1-1
      tmp= sum(f_x1(i,1:local_size_x2))
      call mpi_reduce(tmp,rho(i),1,MPI_REAL8,MPI_SUM,0,sll_world_collective%comm,ierr)
    end do
    rho(np_x1)=rho(1)
    rho = 1._f64-sim%mesh2d%delta_eta2*rho
    call solve(poisson_1d, efield, rho)


    ! Ponderomotive force at initial time. We use a sine wave
    ! with parameters k_dr and omega_dr.
    istep = 0
    if (sim%driven) then
      call PFenvelope(adr, istep*sim%dt, sim%tflat, sim%tL, sim%tR, sim%twL, sim%twR, &
          sim%t0, sim%turn_drive_off)
      do i = 1, np_x1
        e_app(i) = sim%Edrmax*adr*sim%kx&
          *sin(sim%kx*real(i-1,f64)*sim%mesh2d%delta_eta1)
      enddo
    endif

    ! write initial fields
    if(sll_get_collective_rank(sll_world_collective)==0)then
      write(rho_diag,*) rho
      write(eapp_diag,*) e_app
      write(adr_diag,*) istep*sim%dt, adr
    endif
    
    !transposition
    call apply_remap_2D( remap_plan_x1_x2, f_x1, f_x2 )
    call compute_local_sizes_2d( layout_x2, local_size_x1, local_size_x2 )
    global_indices(1:2) = local_to_global_2D( layout_x2, (/1, 1/) )

    do istep = 1, sim%num_iterations


      
      do i = 1,local_size_x1
        ig=i+global_indices(1)-1
        alpha = -(efield(ig)+e_app(ig)) * 0.5_f64 * sim%dt
        alpha = 0._f64
        f1d(1:np_x2) = f_x2(i,1:np_x2)
        !print *,'#np_x2=',np_x2
        f1d(1:np_x2) = f_x2(i,1:np_x2)
        f1d = sim%interp_v%interpolate_array_disp(np_x2, f1d, alpha)
        !f1d = sim%interp_v%interpolate_array_disp(np_x2, f1d, alpha)
        !f1d(1:np_x2) = sim%interp_v%interpolate_array_disp(np_x2, f1d(1:np_x2), alpha)
        f_x2(i,1:np_x2) = f1d(1:np_x2)
      end do
      
 return     
      
      !transposition
      call apply_remap_2D( remap_plan_x2_x1, f_x2, f_x1 )
      call compute_local_sizes_2d( layout_x1, local_size_x1, local_size_x2 )
      global_indices(1:2) = local_to_global_2D( layout_x1, (/1, 1/) )

 
      do i = 1, local_size_x2
        ig=global_indices(2)
        alpha = (sim%mesh2d%eta2_min + real(i+ig-2,f64) * sim%mesh2d%delta_eta2) * sim%dt
        f1d(1:np_x1) = f_x1(1:np_x1,i)
        !f1d = sim%interp_x%interpolate_array_disp(np_x1, f1d, alpha)
        f_x1(1:np_x1,i)=f1d(1:np_x1)
      end do

      !transposition
      call apply_remap_2D( remap_plan_x1_x2, f_x1, f_x2 )
      call compute_local_sizes_2d( layout_x2, local_size_x1, local_size_x2 )
      global_indices(1:2) = local_to_global_2D( layout_x2, (/1, 1/) )

    
    
    enddo
    


    
    
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
