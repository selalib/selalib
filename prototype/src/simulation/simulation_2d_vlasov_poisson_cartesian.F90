module sll_simulation_2d_vlasov_poisson_cartesian

! mimic of VP1D_deltaf_BSL in the simulation class
! intend to use MPI instead of openmp
! main application is for KEEN waves


#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
  use sll_collective
  use sll_remapper
  use sll_constants
  use sll_cubic_spline_interpolator_1d
  use sll_test_4d_initializer
  use sll_poisson_2d_periodic_cartesian_par
  use sll_cubic_spline_interpolator_1d
  use sll_simulation_base
  implicit none

  type, extends(sll_simulation_base_class) :: &
       sll_simulation_2d_vlasov_poisson_cart
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


      
    ! namelists for data input
    namelist / geom / xmin, Ncx, nbox, vmin, vmax, Ncv
    namelist / interpolator / interpol_x, order_x, interpol_v, order_v
    namelist / time_iterations / dt, nbiter, freqdiag
    namelist / landau / kmode, eps, is_delta_f, driven 
    namelist / drive / t0, twL, twR, tflat, tL, tR, turn_drive_off, Edrmax, omegadr

    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, '#init_vp2d_par_cart() failed to open file ', filename
       STOP
    end if
    read(input_file, geom) 
    read(input_file, time_iterations)
    read(input_file, landau)
    if (driven) then
        read(input_file, drive)
        eps = 0.0  ! no initial perturbation for driven simulation
    end if
    close(input_file)
    
    xmax = nbox * 2._f64 * sll_pi / kmode
    
    
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
        write(param_out_drive,*) t0, twL, twR, tflat, tL, tR, &
          Edrmax, omegadr
        !! warning we have skipped tstart  in
        !! t0, twL, twR, tstart, tflat, tL, tR            
        close(param_out_drive)
      end if

    endif


    
    !sim%dt = dt
    !sim%num_iterations = number_iterations
    ! In this particular simulation, since the system is periodic, the number
    ! of points is the same as the number of cells in all directions.
    !sim%nc_x1 = num_cells_x1
    !sim%nc_x2 = num_cells_x2
  end subroutine init_vp2d_par_cart


  subroutine run_vp2d_cartesian(sim)
    class(sll_simulation_2d_vlasov_poisson_cart), intent(inout) :: sim
  end subroutine run_vp2d_cartesian

  subroutine delete_vp2d_par_cart( sim )
    class(sll_simulation_2d_vlasov_poisson_cart) :: sim
    sll_int32 :: ierr
  end subroutine delete_vp2d_par_cart


end module sll_simulation_2d_vlasov_poisson_cartesian
