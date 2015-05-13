program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"

    use sll_collective
    use sll_pic_1d
    implicit none




    !new_sll_pic_1d(mesh_cells_user, spline_degree_user, numberofparticles_user,&
        !            timesteps_user, timestepwidth_user )



    CHARACTER(LEN=255) :: string  ! stores command line argument
    CHARACTER(LEN=255) :: message           ! use for I/O error messages
    integer :: ios, IO_stat

    !---------------------------------------------------------------------------
    !------------------DEFAULT VALUES-------------------------------------------------------
    sll_int32 :: tsteps = 100     !Number of timesteps
    sll_real64 :: tstepw=0.1_f64   !stepwidth
    sll_int32 :: nmark=40000        !Number of marker particles
    sll_int32 :: femp=7            !exponent for number of mesh cells
    sll_int32 :: sdeg=3            !spline degree

    !    LOGICAL :: gnuplot_inline_output=.TRUE. !write gnuplot output during the simulation

    !Landau damping
    sll_real64 :: lalpha=0.1_f64
    sll_real64 :: lmode=0.5_f64

    !
    sll_real64 :: boxlenpi=4.0_f64

    CHARACTER(LEN=255) :: ppusher!='lfrog'  !Particle pusher  euler
    CHARACTER(LEN=255) :: psolver='fem'
    CHARACTER(LEN=255) :: scenario='landau'

    !Streams
    sll_int32 :: nstreams=1
    logical  :: gnuplot_inline_output_user=.FALSE.

    namelist /paramss/ tsteps, tstepw, nmark, femp, sdeg, ppusher, scenario, psolver, lalpha, gnuplot_inline_output_user

  !------------------ END DEFAULT VALUES --------------------------------------------------

    !    LOGICAL :: gnuplot_inline_output=.TRUE. !write gnuplot output during the simulation

    !Landau damping
   


      !Particle pusher  euler 
     !Poisson solver
      !Which testcase to use

    !Choices are rk4, euler, verlet, lfrog, v_lfrog
    sll_int32 :: ppusher_int=0, psolver_int=0

    CHARACTER(LEN=255) :: path="./"
    character(len=256) :: filename
    integer :: gpinline=0,deltaf=0
    sll_int32, parameter  :: input_file =99
NAMELIST /cmd/ tsteps, tstepw , nmark, scenario, nstreams, femp, sdeg, lalpha, lmode,&
        ppusher,gpinline,path,boxlenpi,deltaf, psolver

!**************************************************************
    ! Read input data from file
    call getarg(1,filename)
    open(unit=input_file, file=trim(filename), IOStat=IO_stat)
    if( IO_stat /= 0 ) then
        print *, 'init_file() failed to open file ', filename
        STOP
    end if
    read(input_file,paramss)
    close(input_file)
!**************************************************************

    !Add input variables to namelist
    

    call sll_boot_collective()

    call sll_collective_barrier(sll_world_collective)

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
    if (deltaf/=0) then
        enable_deltaf=.TRUE.
    else
        enable_deltaf=.FALSE.

    endif

    if (gpinline/=0) then
        gnuplot_inline_output_user=.TRUE.

    endif

    interval_a=0
    interval_b=boxlenpi*sll_pi

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

end program unit_test

