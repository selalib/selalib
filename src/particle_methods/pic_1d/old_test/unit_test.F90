program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

     use sll_m_collective , only :       sll_v_world_collective , sll_collective_barrier ,&
     sll_s_boot_collective
    
    
    use sll_m_pic_1d  , only : new_sll_pic_1d, sll_solver_fourier, sll_solver_fd, sll_solver_fem, sll_pic1d_testcase_quiet ,&
    sll_pic1d_testcase_landau, sll_pic1d_testcase_ionbeam, sll_pic1d_testcase_bumpontail, sll_pic1d_ppusher_verlet ,&     
    enable_deltaf, interval_a, interval_b, landau_alpha, landau_mode, pic1d_testcase, root_path, sll_p_pi ,&
    sll_pic1d_ppusher_euler, sll_pic1d_ppusher_shift, sll_pic1d_ppusher_rk2, sll_pic1d_ppusher_rk4, sll_pic1d_ppusher_none ,&
    sll_pic1d_ppusher_merson, sll_pic1d_ppusher_leapfrog_v, set_loading_parameters ,&
    sll_pic1d_ppusher_leapfrog, sll_pic1d_ppusher_heun, sll_solver_spectral, sll_pic_1d_run, destroy_sll_pic_1d
  
    
    implicit none




    !new_sll_pic_1d(mesh_cells_user, spline_degree_user, numberofparticles_user,&
        !            timesteps_user, timestepwidth_user )



    character(len=255) :: string  ! stores command line argument
    character(len=255) :: message           ! use for I/O error messages
    integer :: ios, IO_stat

    !---------------------------------------------------------------------------
    !------------------DEFAULT VALUES-------------------------------------------------------
    sll_int32 :: tsteps = 100     !Number of timesteps
    sll_real64 :: tstepw=0.1_f64   !stepwidth
    sll_int32 :: nmark=40000        !Number of marker particles
    sll_real64 :: final_time
    sll_int32 :: sdeg=3            !spline degree

    !    LOGICAL :: gnuplot_inline_output=.TRUE. !write gnuplot output during the simulation

    !Landau damping
    sll_real64 :: lalpha!=0.1_f64
    sll_real64 :: lmode=0.5_f64

    !
    sll_real64 :: boxlenpi=4.0_f64

    character(len=255) :: ppusher!='lfrog'  !Particle pusher  euler
    character(len=255) :: psolver='fem'
    character(len=255) :: scenario='landau'

    !Streams
    sll_int32 :: nstreams=1
    logical  :: gnuplot_inline_output_user=.FALSE.
    LOGICAL :: deltaf=.TRUE.

  !------------------ END DEFAULT VALUES --------------------------------------------------

    !    LOGICAL :: gnuplot_inline_output=.TRUE. !write gnuplot output during the simulation

    !Landau damping
   


      !Particle pusher  euler 
     !Poisson solver
      !Which testcase to use

    !Choices are rk4, euler, verlet, lfrog, v_lfrog
    sll_int32 :: ppusher_int=0, psolver_int=0

    character(len=255) :: path="./"
    character(len=256) :: filename
    integer :: gpinline=0
    logical :: pi_unit
    sll_int32, parameter  :: input_file =99
    !input   
     sll_int32 :: femp=7            !exponent for number of mesh cells
     
     
    namelist /landau_params/ lalpha,lmode,pi_unit,interval_a,interval_b
    
    
    
    namelist /params/ nmark,final_time  , ppusher, scenario, psolver,gnuplot_inline_output_user,deltaf,nstreams
    
    namelist /numerical_params/sdeg,tstepw,femp
    
    
    
    
NAMELIST /cmd/ tsteps, tstepw ,  scenario, nstreams, femp, sdeg, lalpha, lmode,&
        ppusher,gpinline,path,boxlenpi,deltaf, psolver


!**************************************************************
    ! Read input data from file
    call get_command_argument(1,filename)
    open(unit=input_file, file=trim(filename), IOStat=IO_stat)
    if( IO_stat /= 0 ) then
        print *, 'init_file() failed to open file ', filename
        STOP
    end if
    read(input_file,landau_params)
    read(input_file,numerical_params)
    read(input_file,params)
    close(input_file)
!**************************************************************

tsteps=final_time/tstepw

    !Add input variables to namelist
    

    call sll_s_boot_collective()

    call sll_collective_barrier(sll_v_world_collective)

    !---------------------------------------------------------------------------
    CALL get_command_argument(1,string)                ! get command line arguments
    string="&cmd "//trim(string)//" /"            ! add NAMELIST prefix and terminator
    READ(string,NML=cmd,iostat=ios,iomsg=message) ! internal read of NAMELIST
    WRITE(*,NML=cmd)

    ppusher=trim(ppusher)
    scenario=trim(scenario)
    root_path=trim(path)

    if(ios.ne.0)then
        write(*,*)'NML-Error Message: ',ios
        write(*,*)message

    endif

        enable_deltaf=deltaf
 

    if(pi_unit)then
    interval_a=interval_a*sll_p_pi
    interval_b=interval_b*sll_p_pi
    endif
    selectcase (scenario)
        case("landau")
            pic1d_testcase = SLL_PIC1D_TESTCASE_LANDAU
        case("ionbeam")
            pic1d_testcase = SLL_PIC1D_TESTCASE_IONBEAM
            interval_a=0
            interval_b=200  !0.0022_f64 !20mm
        case("quiet")
            pic1d_testcase = SLL_PIC1D_TESTCASE_QUIET
        case("bump")
            pic1d_testcase = SLL_PIC1D_TESTCASE_BUMPONTAIL
            landau_alpha=0.001_f64
            landau_mode=0.5_f64
    end select



    call set_loading_parameters(lalpha,lmode,nstreams)




    selectcase (ppusher)
        case("rk4")
            ppusher_int=SLL_PIC1D_PPUSHER_RK4
        case("verlet")
            ppusher_int=SLL_PIC1D_PPUSHER_VERLET
        case("euler")
            ppusher_int=SLL_PIC1D_PPUSHER_EULER
        case("lfrog_v")
            ppusher_int=SLL_PIC1D_PPUSHER_LEAPFROG_V
        case("lfrog")
            ppusher_int=SLL_PIC1D_PPUSHER_LEAPFROG
        case("rk2")
            ppusher_int=SLL_PIC1D_PPUSHER_RK2
        case("rk3")
            ppusher_int=SLL_PIC1D_PPUSHER_RK2
        case("merson")
            ppusher_int=SLL_PIC1D_PPUSHER_MERSON
        case("heun")
            ppusher_int=SLL_PIC1D_PPUSHER_HEUN
        case("none")
            ppusher_int=SLL_PIC1D_PPUSHER_NONE
        case("shift ")
            ppusher_int=SLL_PIC1D_PPUSHER_SHIFT
        case default
            ppusher_int=SLL_PIC1D_PPUSHER_NONE
    end select

    selectcase (psolver)
        case("fem")
            psolver_int = SLL_SOLVER_FEM
        case("fd")
            psolver_int = SLL_SOLVER_FD
        case("fourier")
            psolver_int = SLL_SOLVER_FOURIER
        case("spec")
            psolver_int = SLL_SOLVER_SPECTRAL
        case default
            psolver_int = SLL_SOLVER_FEM

    end select

    call  new_sll_pic_1d(2**femp, sdeg, nmark, tsteps , tstepw,ppusher_int , psolver_int)


    !call  new_sll_pic_1d(2**, 3, 10000, 1000 , 0.001_f64 )

    !Nonlinear Landau damping
    !call  new_sll_pic_1d(2**12, 3, 20000, 1000 , 0.001_f64 )

    call sll_pic_1d_run(gnuplot_inline_output_user)


    call destroy_sll_pic_1d
   write(*,*) 'nber of steps  =',tsteps
    
end program unit_test


