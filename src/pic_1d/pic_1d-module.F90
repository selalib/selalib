module sll_pic_1d_Class
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"

    use sll_constants
    !use sll_logical_meshes
    !use sll_poisson_1d_periodic
    !use sll_poisson_solvers
    use sll_collective !Parallel operations

    !  use hdf5, only: HID_T,HSIZE_T,HSSIZE_T
    !  ! #ifdef HDF5_PARALLEL
    !          use sll_hdf5_io_parallel
    !        use sll_xdmf_parallel
    ! #endif

    use sll_arbitrary_degree_splines
    use gauss_legendre_integration
    use sll_bspline_fem_solver_1d
    use sll_visu_pic !Visualization with gnuplot
    use pic_1d_particle_loading !should be integrated here
    use sll_timer
    use pic_postprocessing


    implicit none

    integer :: ierr
    integer, private :: i, j


    type particle
        sll_real64, DIMENSION(:), allocatable:: particleposition
        sll_real64, DIMENSION(:), allocatable :: particlespeed
        sll_real64, DIMENSION(:), allocatable :: particleweight          ! cf. y_k in delta f and c_k in full f pic

        endtype


    !1d1v data for electrostatic PIC
    sll_real64, DIMENSION(:), allocatable:: particleposition
    sll_real64, DIMENSION(:), allocatable :: particlespeed
    sll_real64, DIMENSION(:), allocatable :: particleweight          ! cf. y_k in delta f and c_k in full f pic
    sll_real64, DIMENSION(:), allocatable :: particleweight_constant !cf. c_k
    sll_real64, DIMENSION(:), allocatable:: steadyparticleposition !Steady non-moving objects aka Ions
    !Constants
    sll_real64 :: particle_mass=1.0_f64, particle_charge=1.0_f64

    !CHARACTER(LEN=255) :: particle_pusher
    sll_int32 :: particle_pusher

    integer, private :: timestep

    !Data to be collected throughout a run
    sll_real64, dimension(:), allocatable :: fieldenergy, kineticenergy, impulse, thermal_velocity_estimate

    !sll_real64 :: initial_fieldenergy , initial_kineticenergy , initial_total_energy
    sll_real64 :: initial_fem_inhom , fem_inhom
    sll_int64    ::  fastest_particle_idx
    sll_real64 , DIMENSION(:), allocatable:: eval_solution


    sll_int32 :: timesteps = 100
    sll_real64 :: timestepwidth=0.1_f64
    sll_int32 :: nparticles = 1000  !GLOBAL, marker particles
    sll_int32 :: particles_system= 1 !Number of the particles in the system
    LOGICAL :: gnuplot_inline_output=.FALSE. !write gnuplot output during the simulation
    sll_real64,parameter :: plasma_frequency = 1.0_f64 != sqrt( sll_e_charge**2/sll_e_mass/sll_epsilon_0)
    sll_real64 ,parameter:: thermal_velocity = 1.0_f64

    !Mesh parameters
    sll_int32 :: mesh_cells = 2**8
    sll_int32                   :: spline_degree=3

    sll_real :: interval_a=0, interval_b=4.0_f64*sll_pi
    sll_real64 :: mesh_delta_t
    sll_real64, dimension(:), allocatable, private   :: knots
    sll_int32                  :: num_pts

    integer, private ::  phasespace_file_id=-1

    CHARACTER(LEN=255) :: root_path
    !          interface deferred
    !            function sll_pic_1d_electric_field(x) result(E)
    !        real, dimension(:), intent(in) :: x
    !        real, dimension(:) :: E
    !        endfunction
    !              endinterface
    procedure(sll_pic_1d_electric_field), pointer:: electric_field_external

    !type(sll_logical_mesh_1d), pointer :: mesh1d
    !    type(poisson_1d_periodic), pointer :: poisson1dsolver


    sll_real64, dimension(:), allocatable ::electricpotential_interp
    !For parallelization
    sll_int32, private :: coll_rank, coll_size

    ! mesh1d=> new_logical_mesh_1d(mesh_cells, interval_a,interval_b)
    !call sll_display(mesh1d)

    !Definitions for different particle pushers, enumerator
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_NONE       = 0
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_EULER   = 1    !Explicit EULER
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_VERLET   = 2
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_RK2       = 3 !Runge Kutta 2
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_RK4       = 4   !Runge Kutta 4
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_HEUN       = 5
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_SHIFT       = 6 !Shift by velocity
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_LEAPFROG       = 7
    sll_int32, parameter :: SLL_PIC1D_PPUSHER_LEAPFROG_V       = 8 !Variational Leapfrog

    sll_int32, parameter :: SLL_PIC1D_FULLF       = 1
    sll_int32, parameter :: SLL_PIC1D_DELTAF       = 2

contains

    !<
    subroutine new_sll_pic_1d(mesh_cells_user, spline_degree_user, numberofparticles_user,&
            timesteps_user, timestepwidth_user,particle_pusher_user )
        implicit none
        sll_int32, intent(in) :: mesh_cells_user
        sll_int32 , intent(in):: spline_degree_user
        sll_int32 , intent(in):: numberofparticles_user
        sll_int32, intent(in) :: timesteps_user
        sll_real64 , intent(in):: timestepwidth_user
        sll_int32 :: idx
        !character(len=255), intent(in):: particle_pusher_user
        sll_int32:: particle_pusher_user
        SLL_ASSERT( is_power_of_two( int( mesh_cells,i64)))
        mesh_cells=mesh_cells_user
        spline_degree=spline_degree_user
        nparticles=numberofparticles_user
        timesteps=timesteps_user
        timestepwidth=timestepwidth_user
        !particle_pusher=trim(particle_pusher_user)
        particle_pusher=particle_pusher_user
        !############## PARALLEL ##################
        call sll_boot_collective()

        !really global variables
        coll_rank = sll_get_collective_rank( sll_world_collective )
        coll_size = sll_get_collective_size( sll_world_collective )
        call sll_collective_barrier(sll_world_collective)


        !The Mesh is needed on every Node, it is global
        mesh_delta_t = (interval_b - interval_a)/mesh_cells
        SLL_CLEAR_ALLOCATE(knots(1:mesh_cells+1),ierr)
        SLL_CLEAR_ALLOCATE(electricpotential_interp(1:mesh_cells),ierr)



        knots(1)=interval_a
        do i=2,size(knots)
            knots(i)=  knots(1) + (i-1)*mesh_delta_t*1.0_f64
        enddo
        knots(size(knots))=interval_b
        SLL_ASSERT(knots(size(knots))==interval_b)


        if (coll_rank==0) then
            print *, "Size of MPI-Collective: ", coll_size
        endif

        call sll_collective_barrier(sll_world_collective)


        electric_field_external=>sll_pic_1d_electric_field
    endsubroutine new_sll_pic_1d

    !<Destructor
    subroutine destroy_sll_pic_1d()
        call sll_halt_collective()

    endsubroutine


    subroutine sll_pic1d_write_phasespace(timestep)
        sll_int32, intent(in) :: timestep
        sll_int32:: grid=100
        if (coll_rank==0) then
            !            if (phasespace_file_id==-1) then
            !                        call sll_new_file_id(phasespace_file_id, ierr)
            !                    open(phasespace_file_id, file = trim(root_path)//'phasespace.dat')
            !            endif
            !
            !            write(phasespace_file_id,*)  particleposition, particlespeed
           ! if (phasespace_file_id==-1) then
                            SLL_ASSERT(    minval(particleposition) >=interval_a)
                            SLL_ASSERT(    maxval(particleposition) <=interval_b)
                            grid=floor(sqrt(nparticles*1.0_f64)/2.0_f64)
                            call distribution_xdmf("phasespace", particleposition, particlespeed, particleweight, &
                                        interval_a-1.0_f64, interval_b+1.0_f64, grid, minval(particlespeed)-1.0_f64, maxval(particlespeed)+1.0_f64, grid, timestep)
                                  !          phasespace_file_id=0
            !endif
        endif
    endsubroutine

    function ionbeam(x) result(y)
        sll_real64, dimension(:),intent(in) :: x
        sll_real64, dimension(size(x)) :: y
        y=0
        where (x>=interval_a .and. x-interval_a<=(interval_b-interval_a)/3.0_f64) y=3.0_f64
    end function ionbeam


    !<Calculates Error estimates for different Testcases using preknown analytical
    !<solutions. Call after the first field solve.
    subroutine sll_pic_1d_initial_error_estimates()
        sll_real64 :: num_err_mean, num_err_noise_max, num_err_noise_l2, num_err_seminorm
        sll_real64, dimension(size(knots)-1) :: analytical_solution
        sll_real64 :: flag=1.0_f64
        selectcase(pic1d_testcase)
            case(SLL_PIC1D_TESTCASE_LANDAU)
                !The first field solve is a Monte Carlo estimator
                if (coll_rank==0) then
                    !analytical_solution=landau_alpha/((interval_b-interval_a)*landau_mode)*&
                    if (landau_alpha/=0) then
                        analytical_solution=landau_alpha/(landau_mode)*&
                            sin(landau_mode*(knots(1:mesh_cells)+(sll_pi/6)*(knots(2)-knots(1))))
                        num_err_noise_max= maxval(abs((eval_solution(1:mesh_cells)-analytical_solution)))/&
                            maxval(abs(analytical_solution))
                        num_err_noise_l2=sum(((eval_solution(1:mesh_cells)-analytical_solution))**2)/&
                            sum(analytical_solution**2)
                    else
                        num_err_noise_max= maxval(abs((eval_solution(1:mesh_cells))))
                    endif


                    print *, "Error in MC estimate of E-field: ", num_err_noise_max, num_err_noise_l2
                endif

                    !num_err_seminorm=landau_alpha**2/(2*(interval_b-interval_a)*landau_mode**2)
                    if (landau_alpha/=0) then
                        num_err_seminorm=landau_alpha**2/(2*landau_mode**2)*(interval_b-interval_a)
                        num_err_seminorm=(sll_pic1d_calc_fieldenergy()-num_err_seminorm)/num_err_seminorm
                    else
                        num_err_seminorm=sll_pic1d_calc_fieldenergy()
                    endif
                    if (coll_rank==0) print *, "Error in MC estimate of E-Energy: ", num_err_seminorm
        endselect

    endsubroutine



    subroutine sll_pic_1d_run(gnuplot_inline_output_user)
        logical, intent(in), optional ::gnuplot_inline_output_user
        sll_real64, dimension(size(knots)-1) :: analytical_solution
        sll_int32 :: idx,jdx
        !Timers
        type(sll_time_mark)  :: tstart, tstop

        gnuplot_inline_output=gnuplot_inline_output_user
        !For Parallelization we take the number of particles nparticles and distribute these
        !particles on the different cores
        !this also means we only plot the particles on one core
        !energy calculation and everything has to be paralellized
        nparticles=nparticles/coll_size
        print*, "#Core ", coll_rank, " handles particles", coll_rank*nparticles +1, "-", (coll_rank+1)*nparticles
        call sll_collective_barrier(sll_world_collective)
        if (coll_rank==0) print*, "#Total Number of particles: ", nparticles*coll_size


        !Scalar values to be recorded
        SLL_CLEAR_ALLOCATE(fieldenergy(1:timesteps+1), ierr)
        SLL_CLEAR_ALLOCATE(kineticenergy(1:timesteps+1), ierr)
        SLL_CLEAR_ALLOCATE(impulse(1:timesteps+1), ierr)
        SLL_CLEAR_ALLOCATE(thermal_velocity_estimate(1:timesteps+1), ierr)


        SLL_CLEAR_ALLOCATE( steadyparticleposition(1:nparticles),ierr)
        SLL_CLEAR_ALLOCATE(particleposition(1:nparticles),ierr)
        SLL_CLEAR_ALLOCATE(particlespeed(1:nparticles),ierr)
        SLL_CLEAR_ALLOCATE(particleweight(1:nparticles),ierr)
        SLL_CLEAR_ALLOCATE(particleweight_constant(1:nparticles),ierr)

        SLL_CLEAR_ALLOCATE(eval_solution(1:size(knots)),ierr)



        particle_mass=1.0_f64

        particle_charge=1.0_f64
        if (coll_rank==0) write(*,*) "#PIC1D: Loading particles..."
        call load_particles (nparticles, interval_a, interval_b,&
            steadyparticleposition, particleposition, &
            particlespeed, particleweight, particleweight_constant)

        call sll_collective_barrier(sll_world_collective)
        if (coll_rank==0) write(*,*) "#PIC1D: Particles loaded..."

        !Determine index of fastest particle
        fastest_particle_idx=MAXLOC(particlespeed, dim=1)



        !Do precalculations, especially for the FEM solver
        !call sll_bspline_fem_solver_1d_initialize(knots, spline_degree, steadyparticleposition, -(sll_e_charge/sll_epsilon_0) )
        if (coll_rank==0) print *, "#Number of FEM-Cells ", mesh_cells
        !call sll_bspline_fem_solver_1d_initialize(knots, spline_degree, knots, 1.0_f64)
        !call sll_bspline_fem_solver_1d_initialize(knots, spline_degree, steadyparticleposition, 1.0_f64, sll_world_collective)
        call sll_bspline_fem_solver_1d_initialize(knots, spline_degree, steadyparticleposition, (interval_b-interval_a), sll_world_collective)

        !call sll_bspline_fem_solver_1d_initialize(knots, spline_degree)
        call sll_bspline_fem_solver_1d_set_weights(particleweight)


        !load Ions constant
        !call sll_bspline_fem_solver_1d_set_inhomogenity_constant(1.0_f64/(1.0_f64*mesh_cells*(interval_b-interval_a)) )
        !call sll_bspline_fem_solver_1d_set_inhomogenity_constant(1.0_f64*nparticles*coll_size/(1.0_f64*real(mesh_cells)))

        !  if  (pic1d_testcase == SLL_PIC1D_TESTCASE_LANDAU)

        selectcase(pic1d_testcase)
            case(SLL_PIC1D_TESTCASE_QUIET)
                call sll_bspline_fem_solver_1d_set_inhomogenity_constant(0.0_f64 )
            case(SLL_PIC1D_TESTCASE_LANDAU)
                !load constant ion background
                call sll_bspline_fem_solver_1d_set_inhomogenity_constant(0.0_f64 )
        endselect


        if (enable_deltaf .eqv. .TRUE.) then

            call sll_bspline_fem_solver_1d_set_weights(particleweight_constant)

            call sll_bspline_fem_solver_1d_set_inhomogenity_constant(0.0_f64)
            call sll_bspline_fem_solver_1d_add_inhomogenity_function(sll_pic_1d_landaudamp_PDF,10  )

        endif


        !Start with PIC Method
        !Do the initial solve of the field
        particleposition=sll_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)
        call sll_bspline_fem_solver_1d_solve(particleposition)
        print *, "#Initial field solve"


        print *, "#Test initial field, and loading condition for landau damping"
        call sll_bspline_fem_solver_1d_eval_solution_derivative(knots(1:mesh_cells)+ (sll_pi/6)*(knots(2)-knots(1)), eval_solution(1:mesh_cells))
        call  sll_pic_1d_initial_error_estimates()




        !Information on initial field

        print *, "#Initial Field average: " , sum( sll_bspline_fem_solver_1d_get_solution())
        print *, "#Initial Field variance: " , sum((sum( sll_bspline_fem_solver_1d_get_solution())/mesh_cells &
            - sll_bspline_fem_solver_1d_get_solution())**2)/mesh_cells

        !eval_solution=sll_bspline_fem_solver_1d_get_solution()
        !initial_fem_inhom=sum(sll_bspline_fem_solver_1d_get_inhomogenity())



        !---------------------------Initial Plots------------------------------------------------------------
        if (coll_rank==0) then


            !            call distribution_xv_gnuplot( "phasespace", particleposition, particlespeed, &
                !                                                 minval(particleposition), maxval(particleposition),&
                !                                                  size(particleposition), &
                !                                                  minval(particlespeed), maxval(particlespeed),&
                !                                                  size(particlespeed), &
                !                                                1, 1.0_f64)

!            call sll_bspline_fem_solver_1d_eval_solution(knots(1:mesh_cells), eval_solution(1:mesh_cells))
!
!            open(unit=20, file=trim(root_path)//"initial_field.txt")
!            write (20,*) eval_solution(1:mesh_cells)
!            close(20)
!
!            !Particle density
!            open(unit=20, file=trim(root_path)//"initial_electrondensity.txt")
!
!            write (20,*) bspline_particle_density(particleposition)
!            close(20)
!
!            !Particle density
!            open(unit=20, file=trim(root_path)//"initial_iondensity.txt")
!            write (20,*) bspline_particle_density(steadyparticleposition)
!            close(20)

            !            open(unit=20, file=trim(root_path)//"initial_phasespace.dat")
            !            do idx=1, size(particleposition)
            !                write (20,*) particleposition(idx), particlespeed(idx), particleweight(idx)
            !            enddo
            !            close(20)

            !            open(unit=20, file=trim(root_path)//"initial_electronposition.txt")
            !            do idx=1, size(particleposition)
            !                write (20,*) particleposition(idx)
            !            enddo
            !            close(20)
            !
            !            open(unit=20, file=trim(root_path)//"initial_electronspeed.txt")
            !            do idx=1, size(particleposition)
            !                write (20,*) particlespeed(idx)
            !            enddo
            !            close(20)
            !        print *, sum(abs(particleposition - sll_pic1d_ensure_periodicity( particleposition + 1.0D-3, &
                !            interval_a, interval_b)))

            !                call xv_particles_center_gnuplot('/tmp/initial_phasedensity', particleposition, particlespeed, &
                !                            interval_a, interval_b,&
                !                            -1.1_f64*maxval(particlespeed), 1.1_f64*maxval(particlespeed),&
                !                            ierr, 1.0_f64*timestep )


            print *, "#Field Energy                  " , &
                "Kinetic Energy                  ", "Total Energy                      "
        endif

        !Write initial phasespace
        do jdx=0, coll_size
            call sll_collective_barrier(sll_world_collective)
            if( coll_rank==jdx) then
                if (jdx==0) then
                    open(unit=20, file=trim(root_path)//"initial_phasespace.dat")
                else
                    open(unit=20, file=trim(root_path)//"initial_phasespace.dat",position="append")
                endif
                do idx=1, size(particleposition)
                    write (20,*) particleposition(idx), particlespeed(idx), particleweight(idx)
                enddo
                close(20)
            endif
        enddo

        !----------------------------------------------------------------------------------------------------------

        !open file

        call sll_set_time_mark(tstart)
        !Open HDF5 file for parallel write


        !---------------------------### Main Loop ###----------------------------------------------------
        do timestep=1, timesteps
            call sll_pic1d_write_phasespace(timestep)

            !print *, sum(abs((sll_bspline_fem_solver_1d_get_inhomogenity())))
            !            call xv_particles_center_gnuplot('/tmp/phasedensity', particleposition, particlespeed, &
                !                interval_a, interval_b,&
                !                -1.1_f64*maxval(particlespeed), 1.1_f64*maxval(particlespeed),&
                !                ierr, 1.0_f64*timestep )


            if (nparticles*coll_size<=10000 .AND. (gnuplot_inline_output .eqv. .true.) .AND. (coll_rank==0) .AND. mod(timestep-1,timesteps/100)==0) then
                call particles_center_gnuplot_inline(particleposition, particlespeed, &
                    interval_a, interval_b,&
                    -1.1_f64*maxval(particlespeed), 1.1_f64*maxval(particlespeed),&
                    timestepwidth*(timestep-1) )

            end if



            kineticenergy(timestep)=sll_pic1d_calc_kineticenergy(particle_mass,  particlespeed , particleweight )
            fieldenergy(timestep)=sll_pic1d_calc_fieldenergy()
            impulse(timestep)=sll_pic1d_calc_impulse(particle_mass,  particlespeed , particleweight)
            thermal_velocity_estimate(timestep)=sll_pic1d_calc_thermal_velocity(particlespeed, particleweight)

            if ( (gnuplot_inline_output.eqv. .true.) .AND. coll_rank==0 .AND. mod(timestep-1,timesteps/100)==0  ) then
                call energies_electrostatic_gnuplot_inline(kineticenergy(1:timestep), fieldenergy(1:timestep), impulse(1:timestep),timestepwidth)
            else

            endif


            if ( (gnuplot_inline_output.eqv. .true.) .AND. coll_rank==0 .AND. mod(timestep-1,timesteps/100)==0 ) then
                call sll_bspline_fem_solver_1d_eval_solution(knots(1:mesh_cells), electricpotential_interp)
                call electricpotential_gnuplot_inline( electricpotential_interp,knots(1:mesh_cells) )

            endif

            if ( coll_rank==0) then
                !Stop if Error ist catastrophal
                if ( (kineticenergy(timestep)+fieldenergy(timestep) -(kineticenergy(1)+fieldenergy(1)))&
                        /(kineticenergy(1)+fieldenergy(1)) > 10000.0_f64) then
                    print *, "Its Over, Giving Up!"
                    stop
                endif
            endif


            if ((gnuplot_inline_output .eqv. .FALSE.) .AND. coll_rank==0) then
                print *, timestep, fieldenergy(timestep),kineticenergy(timestep)

            endif


            !Push all particles at once
            !
            !




            selectcase (particle_pusher)
                case(SLL_PIC1D_PPUSHER_RK4)
                    call  sll_pic_1d_rungekutta4_step_array( particleposition, &
                        particlespeed, timestepwidth)
                case(SLL_PIC1D_PPUSHER_VERLET)
                    call sll_pic_1d_Verlet_scheme( particleposition, &
                        particlespeed, timestepwidth)
                case(SLL_PIC1D_PPUSHER_EULER)
                    call sll_pic_1d_explicit_euler( particleposition, &
                        particlespeed, timestepwidth)
                case(SLL_PIC1D_PPUSHER_LEAPFROG_V)
                    call  sll_pic_1d_variational_leap_frog( particleposition, &
                        particlespeed, timestepwidth)
                case(SLL_PIC1D_PPUSHER_LEAPFROG)
                    call  sll_pic_1d_leap_frog( particleposition, &
                        particlespeed, timestepwidth)
                case(SLL_PIC1D_PPUSHER_RK2)
                    call sll_pic_1d_rungekutta2( particleposition, &
                        particlespeed, timestepwidth)
                case(SLL_PIC1D_PPUSHER_HEUN)
                    call sll_pic_1d_heun( particleposition, &
                        particlespeed, timestepwidth)
                case(SLL_PIC1D_PPUSHER_NONE)
                    print *, "No Particle pusher choosen"
                case(SLL_PIC1D_PPUSHER_SHIFT)
                    particleposition=particleposition+timestepwidth*(interval_b-interval_a)
                    particleposition=sll_pic1d_ensure_periodicity( particleposition, &
                        interval_a, interval_b)

                    call sll_bspline_fem_solver_1d_solve(particleposition)

            end select




            if (coll_rank==0 .AND. timestep==1) then
                call sll_set_time_mark(tstop)
                !Calculate remaining Time
                print *, "Remaining Time: " , (sll_time_elapsed_between(tstart,tstop)/timestep)*real(timesteps-timestep,i64)
            endif

        enddo

        !Last update of relevant data
        kineticenergy(timestep)=sll_pic1d_calc_kineticenergy(particle_mass,  particlespeed , particleweight )
        fieldenergy(timestep)=sll_pic1d_calc_fieldenergy()
        impulse(timestep)=sll_pic1d_calc_impulse(particle_mass,  particlespeed , particleweight)
        thermal_velocity_estimate(timestep)=sll_pic1d_calc_thermal_velocity(particlespeed, particleweight)


        close(25)

        !Save Results to file
        call sll_pic1d_write_result("pic1dresult", kineticenergy, fieldenergy, impulse )


        if (coll_rank==0) then
            call sll_set_time_mark(tstop)
            !Calculate remaining Time
            print *, "Overall Time: " , sll_time_elapsed_between(tstart,tstop)
        endif
        if (coll_rank==0) call det_landau_damping((/ ( timestep*timestepwidth, timestep = 0, timesteps) /),fieldenergy)

        call sll_bspline_fem_solver_1d_destroy

    endsubroutine

    !Dummy function for solving the potential Phi
    !Input x \in [a, b] is the position in the electric field
    !Result is $\frac{DPhi}{dx}(x)$
    !This evaluates with the solution of the weak Poisson equation
    !the solution vector at the
    !    function sll_pic_1d_DPhidx(x, bspline, knots_mesh, fem_Phi_solution ) &
        !            result(DPhidx)
    !        sll_real64, dimension(:), intent(in) :: x
    !        type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline
    !        sll_real64, dimension(:), intent(in)     :: knots_mesh
    !        sll_real64, dimension(:), intent(in) :: fem_Phi_solution
    !        sll_real :: DPhidx(size(x))
    !        sll_real, dimension(:) :: x_trunc(size(x))
    !        !Account for periodic geometry
    !        x_trunc=sll_pic1d_ensure_periodicity(x,knots_mesh(1), knots_mesh(size(knots_mesh)))
    !
    !        DPhidx=bspline_basis_to_realvals(bspline, knots_mesh, fem_Phi_solution, x_trunc, b_spline_derivatives_at_x)
    !    endfunction



    !    function sll_pic_1d_solve_DPhidx(bspline, knots_mesh, fem_Phi_solution, particleposition, fem_inhomogenity_steady ) &
        !        sll_real64, dimension(:), intent(in) :: eval_points
    !        type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline
    !        sll_real64, dimension(:), intent(in)     :: knots_mesh
    !        sll_real64, dimension(:), intent(out) :: fem_Phi_solution
    !        sll_real64, dimension(:), intent(inout) particleposition
    !        sll_real :: DPhidx(size(x))
    !        sll_real64, dimension(:), intent(in) :: fem_inhomogenity_steady
    !        sll_real64, dimension(:), intent(inout) :: fem_inhomogenity
    !
    !
    !        !Account for periodic geometry
    !        x_trunc=sll_pic1d_ensure_periodicity(particleposition,knots_mesh(1), knots_mesh(size(knots_mesh)))
    !
    !        fem_inhomogenity=(sll_e_charge/sll_epsilon_0)   &
        !                *( fem_inhomogenity_steady + &
        !                interpolate_particles_bsplines( arbitrarydegree1D, knots, particleposition, b_splines_at_x))
    !
    !
    !            fem_Phi_solution=sll_bspline_fem_solve_poisson_1d(arbitrarydegree1D, knots,fem_inhomogenity)
    !
    !    endfunction




    subroutine sll_pic_1d_rungekutta4_step_array(x_0, v_0, h)
        sll_real64, intent(in):: h
        !procedure  (sll_pic_1d_DPhidx):: f_ode
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        !sll_real64 :: x_h,v_h
        sll_real64, dimension(:) :: k_x1(size(x_0)), k_x2(size(x_0)),&
            k_x3(size(x_0)), k_x4(size(x_0)), k_v1(size(x_0)), &
            k_v2(size(x_0)), k_v3(size(x_0)), k_v4(size(x_0)), stage_DPhidx(size(x_0))

        x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)

        !--------------------Stage 1-------------------------------------------------
        !call sll_bspline_fem_solver_1d_solve(x_0)
        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+electric_field_external(x_0))*plasma_frequency*thermal_velocity
        !--------------------Stage 2-------------------------------------------------
        k_x2= h*(v_0  + 0.5_f64  *k_v1)
        !call sll_bspline_fem_solver_1d_solve(sll_pic1d_ensure_periodicity(x_0 + (k_x1 + 2.0_f64*k_x2)/6.0_f64,  interval_a, interval_b))
        call sll_bspline_fem_solver_1d_solve(x_0 + 0.5_f64 *k_x1)
        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0 + 0.5_f64 *k_x1, stage_DPhidx)
        k_v2= h*( stage_DPhidx+electric_field_external(sll_pic1d_ensure_periodicity(x_0 + 0.5_f64 *k_x1,  interval_a, interval_b)))&
                                            *plasma_frequency*thermal_velocity
        !--------------------Stage 3-------------------------------------------------
        k_x3= h*(v_0+ 0.5_f64 * k_v2)
        !call sll_bspline_fem_solver_1d_solve(sll_pic1d_ensure_periodicity(x_0 + (1*k_x1 + 2.0_f64*k_x3)/6.0_f64,  interval_a, interval_b))
        call sll_bspline_fem_solver_1d_solve(x_0 + 0.5_f64 *k_x2)
        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0 + 0.5_f64 *k_x2, stage_DPhidx)
        k_v3= h*(stage_DPhidx+electric_field_external(x_0 + 0.5_f64 *k_x2))*plasma_frequency*thermal_velocity
        !--------------------Stage 4-------------------------------------------------
        k_x4= h*(v_0+ k_v3)
        !call sll_bspline_fem_solver_1d_solve(sll_pic1d_ensure_periodicity(x_0 + (k_x1 + 2.0_f64*k_x2 + 2.0_f64*k_x3 +k_x4  )/6.0_f64,  interval_a, interval_b))
        call sll_bspline_fem_solver_1d_solve(x_0 +  k_x3)
        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0 +  k_x3, stage_DPhidx)
        k_v4= h*(stage_DPhidx+electric_field_external(sll_pic1d_ensure_periodicity(x_0 +  k_x3,  interval_a, interval_b)))&
                    *plasma_frequency*thermal_velocity
        !Perform step---------------------------------------------------------------
        x_0= x_0 + (  k_x1 +  2.0_f64 *k_x2 +  2.0_f64*k_x3 + k_x4 )/6.0_f64
        v_0= v_0 + (  k_v1 +  2.0_f64 *k_v2 +  2.0_f64*k_v3 + k_v4)/6.0_f64

        x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_bspline_fem_solver_1d_solve(x_0)

    endsubroutine

    subroutine sll_pic_1d_rungekutta2(x_0, v_0, h)
        sll_real64, intent(in):: h
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(:) :: k_x1(size(x_0)), k_x2(size(x_0)),&
            k_x3(size(x_0)), k_x4(size(x_0)), k_v1(size(x_0)), &
            k_v2(size(x_0)), stage_DPhidx(size(x_0))

        x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)

        !--------------------Stage 1-------------------------------------------------
        !call sll_bspline_fem_solver_1d_solve(x_0)
        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx)*plasma_frequency*thermal_velocity
        !--------------------Stage 2-------------------------------------------------
        k_x2= h*(v_0  + 0.5_f64  *k_v1)
        call sll_bspline_fem_solver_1d_solve(x_0 + 0.5_f64 *k_x1)
        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0 + 0.5_f64 *k_x1, stage_DPhidx)
        k_v2= h*( stage_DPhidx)*plasma_frequency*thermal_velocity
        !Perform step---------------------------------------------------------------
        x_0= x_0 + k_x2
        v_0= v_0 + k_v2

        x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_bspline_fem_solver_1d_solve(x_0)

    endsubroutine

    subroutine sll_pic_1d_heun(x_0, v_0, h)
        sll_real64, intent(in):: h
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(:) :: k_x1(size(x_0)), k_x2(size(x_0)), &
            k_v1(size(x_0)), k_v2(size(x_0)),stage_DPhidx(size(x_0))

        x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)

        !--------------------Stage 1-------------------------------------------------
        !call sll_bspline_fem_solver_1d_solve(x_0)
        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0, stage_DPhidx)
        k_x1= h*v_0
        k_v1= h*(stage_DPhidx+electric_field_external(x_0))*plasma_frequency*thermal_velocity
        !--------------------Stage 2-------------------------------------------------
        call sll_bspline_fem_solver_1d_solve(x_0 + k_x1)
        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0 + k_x1, stage_DPhidx)
        k_x2= h*(v_0  + k_v1)
        k_v2= h*( stage_DPhidx+electric_field_external(x_0 + k_x1))*plasma_frequency*thermal_velocity

        !Perform step---------------------------------------------------------------
        x_0= x_0 + 0.5_f64 *(k_x1+k_x2   )
        v_0= v_0 + 0.5_f64 *(k_v1+k_v2  )

        x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_bspline_fem_solver_1d_solve(x_0)

    endsubroutine



    subroutine sll_pic_1d_explicit_euler(x_0, v_0, h)
        sll_real64, intent(in):: h
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(:) :: stage_DPhidx(size(x_0))

        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)

        !--------------------Stage 1-------------------------------------------------
        !call sll_bspline_fem_solver_1d_solve(x_0)
        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0, stage_DPhidx)

        x_0=x_0 + h*v_0
        v_0=v_0 + h*stage_DPhidx

        x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_bspline_fem_solver_1d_solve(x_0)
    endsubroutine

    !<From Erics lecture notes on Vaslov Equations
    subroutine sll_pic_1d_Verlet_scheme(x_0, v_0, h)
        sll_real64, intent(in):: h
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(size(x_0)) ::  DPhidx, x_1, v_1, v_05

        !x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)

        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0, DPhidx)
        v_05=v_0+ 0.5_f64*h*(DPhidx+electric_field_external(x_0))*plasma_frequency*thermal_velocity
        x_1=x_0 + h*v_05
        !print *, maxval(abs(v_05))
        x_1=sll_pic1d_ensure_periodicity(x_1,  interval_a, interval_b)
        call sll_bspline_fem_solver_1d_solve(x_1)

        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_1, DPhidx)
        v_1=v_05 + 0.5_f64*h*(DPhidx+electric_field_external(x_1))*plasma_frequency*thermal_velocity
        !Commit Push
        x_0=x_1
        v_0=v_1

    endsubroutine

    subroutine sll_pic_1d_variational_leap_frog(x_0, v_0, h)
        sll_real64, intent(in):: h
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(size(x_0)) ::  DPhidx_0, DPhidx_1, x_1, v_1


        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0, DPhidx_0)
        x_1=x_0+ h*v_0 - ((h**2)/2.0_f64) *DPhidx_0
        x_1=sll_pic1d_ensure_periodicity(x_1,  interval_a, interval_b)

        call sll_bspline_fem_solver_1d_solve(x_1)
        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_1, DPhidx_1)

        v_1=v_0 + 0.5_f64*h* (DPhidx_0 + DPhidx_1)

        !Push
        x_0=x_1
        v_0=v_1
    endsubroutine

    subroutine sll_pic_1d_leap_frog(x_0, v_0, h)
        sll_real64, intent(in):: h
        sll_real64, dimension(:) ,intent(inout) :: x_0
        sll_real64, dimension(:) ,intent(inout) :: v_0
        sll_real64, dimension(size(x_0)) ::  DPhidx_0, DPhidx_1, x_1, v_1

        x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_bspline_fem_solver_1d_eval_solution_derivative(x_0, DPhidx_0)
        !DPhidx_0=0

        v_1= v_0 - 0.5_f64 *DPhidx_0*h
        v_1=v_1+  DPhidx_0*h
        !v_1=v_0+ 0.5_f64* DPhidx_0*h
        x_1=x_0+ v_1*h

        !Check the CFL timestep criterion
        !if  ((maxval(v_1) *h)/2>1) print *, "CFL!!"

        !        x_1=x_0 + (2.0/(128.0))
        !        v_1=v_0
        !Push
        x_0=x_1
        v_0=v_1

        x_0=sll_pic1d_ensure_periodicity(x_0,  interval_a, interval_b)
        call sll_bspline_fem_solver_1d_solve(x_0)
    endsubroutine

    !<Nonrelativistic kinetic energy 0.5*m*v^2
    function sll_pic1d_calc_kineticenergy(particlemass,  particlespeed , particleweight ) &
            result(energy)
        sll_real64, DIMENSION(:), intent(in) :: particlespeed
        sll_real64, DIMENSION(:), intent(in):: particleweight
        sll_real64 , intent(in) ::particlemass
        sll_real64 :: energy

        energy=0.5_f64*particlemass*dot_product(particleweight,particlespeed**2)
        call sll_collective_globalsum(sll_world_collective, energy, 0)
    endfunction

    !<overall impulse 0.5*m*v
    function sll_pic1d_calc_impulse(particlemass,  particlespeed , particleweight) &
            result(impulse)
        sll_real64, DIMENSION(:), intent(in) :: particlespeed
        sll_real64, DIMENSION(:), intent(in):: particleweight
        sll_real64 , intent(in) ::particlemass
        sll_real64:: impulse

        impulse=particlemass*dot_product(particleweight,particlespeed)
        call sll_collective_globalsum(sll_world_collective, impulse, 0)
    endfunction

    !<Thermal velocity, mostly defined as the variance of the particle velocity
    function sll_pic1d_calc_thermal_velocity(particlespeed, particleweight) &
            result(vth)
        sll_real64, DIMENSION(:), intent(in) :: particlespeed
        sll_real64, DIMENSION(:), intent(in):: particleweight
        sll_real64 :: vth
        sll_real64 :: mean
        mean=dot_product(particlespeed,particleweight)
        call sll_collective_globalsum(sll_world_collective, mean, 0)

        vth=dot_product((particlespeed-mean)**2,particleweight)
        call sll_collective_globalsum(sll_world_collective, vth, 0)
    endfunction


    !<Nonrelativistic field energy, namely $|\nabla \Phi$|^2$
    !<which is nothing else than the H1 seminorm of the self generated field
    !<additionally we add the potential energy by the external field E
    function sll_pic1d_calc_fieldenergy() &
            result(energy)
        !sll_real64, DIMENSION(:), intent(in) :: particleposition
        sll_real64 :: energy
        energy=0.0_f64
        !Only for constant case
        energy=dot_product(electric_field_external(particleposition),particleposition)
        call sll_collective_globalsum(sll_world_collective, energy, 0)


        if (coll_rank==0)  energy=energy+0.5_f64*bspline_fem_solver_1d_H1seminorm_solution()/(interval_b-interval_a) !/2.0_f64

        !if (coll_rank==0)  energy=sqrt(bspline_fem_solver_1d_L2norm_solution())!  +bspline_fem_solver_1d_H1seminorm_solution()
    endfunction


    !Write all results of the pic simulation to file
    subroutine sll_pic1d_write_result(filename, kinetic_energy, electrostatic_energy, impulse )
        character(len=*), intent(in) :: filename
        sll_real64, dimension(:), intent(in) :: kinetic_energy, electrostatic_energy, impulse
        integer :: k,file_id,file_id_err
        sll_real64,  dimension(size(kinetic_energy)) :: total_energy


        if (coll_rank==0) then

            SLL_ASSERT(size(kineticenergy)==timesteps+1)
            SLL_ASSERT(size(fieldenergy)==timesteps+1)
            SLL_ASSERT(size(impulse)==timesteps+1)

            total_energy=kinetic_energy +electrostatic_energy



            call sll_new_file_id(file_id, ierr)

            !Write Data File
            !            open(file_id, file = plot_name//"_"//fin//'.dat' )
            open(file_id, file = trim(root_path)//filename//'.dat')

            write (file_id, *)  "#Full 1d1v Electrostatic PIC"
            write (file_id,*)  "#Time steps:", timesteps
            write (file_id,*)  "#Time stepwidth:", timestepwidth
            write (file_id,*)  "#Marko particles:", coll_size*nparticles
            write (file_id,*)  "#Particle Pusher:", particle_pusher
            write (file_id,*)  "#Finite Elements: 2^(", log(real(mesh_cells,i64))/log(2.0_f64),")"
            write (file_id,*)  "#Size of MPI Collective: ", coll_size

            write (file_id,*)  "time  ", "kineticenergy  ", "electrostaticenergy  ", "impulse  ", "vthermal"
            do k=1,timesteps+1
                write (file_id,*)  timestepwidth*(k-1), kineticenergy(k), fieldenergy(k), impulse(k),thermal_velocity_estimate(k)
            enddo
            close(file_id)

            call sll_new_file_id(file_id_err, ierr)

            open(file_id_err, file = trim(root_path)//filename//'-errors.dat')

            if (impulse(1)/=0 .AND. total_energy(1)/=0) then
                !Write the relative errors
                do k=1,timesteps+1
                    write (file_id_err,*)  timestepwidth*(k-1),&
                        abs((total_energy(k)-total_energy(1))/total_energy(1)), &
                        abs((  impulse(k)-impulse(1))/impulse(1))
                enddo
            endif
            close(file_id_err)


            !Write Gnuplot file
            open(file_id, file = trim(root_path)//filename//'.gnu')

            !First Plot
            write (file_id, *) "set term x11 1"

            write(file_id,"(A17,G10.3,A4,G10.3,A1)")"set title 'alpha=",landau_alpha,", k=",landau_mode,"'"
            write(file_id, "(A14,G10.5,A1)") "set xrange [0:",timestepwidth*timesteps,"]"
            write (file_id,*)  "set logscale y"
            write (file_id,*)  "unset logscale x"
            write (file_id,*)  "l_alpha=",landau_alpha
            write (file_id,"(A21,G4.3,A18,G10.3,A1)")  "set xlabel 'Time t in", timesteps, &
                " steps, stepwidth=", timestepwidth, "'"
            write (file_id,*)  "set ylabel 'Energy'"


            !Plot expected damping curve if known
            if (landau_alpha==0.001_f64 .AND. landau_mode==0.4_f64) then

                write (file_id,*) "E(x)=abs(0.002*0.424666*exp(-0.0661*x)*cos(1.2850*x-0.3357725))**2"
                write (file_id,*)  "plot '"//filename//".dat' using 1:3 with lines,\"
                write (file_id,*) "E(x) with lines linestyle 2;"
            elseif (landau_alpha/=0) then
                write (file_id,*) "E(x)=0.1*exp(-0.1533*2*x)"
                write (file_id,*)  "plot '"//filename//".dat' using 1:3 with lines,\"
                write (file_id,*) "E(x) with lines linestyle 2;"

            else

                write (file_id,*)  "plot '"//filename//".dat' using 1:3 with lines"


            endif
            write (file_id,*) "set autoscale x; unset logscale;unset xlabel;unset ylabel;unset xrange;unset yrange;"




            !Plot Kernel density estimates for initial distribution
            !!!Enter v_thermal and two stream here
            write (file_id, *) "set term x11 2"
            write (file_id, *)  "set multiplot layout 2,1 rowsfirst title 'Distribution Kernel Density Estimates'"
            write (file_id, *) "set autoscale x; set autoscale y"
            write (file_id,*)  "vel_pdf(x)=1/sqrt(2*pi)*exp(-0.5*(x**2))"
            write (file_id,*)  "set title 'Velocity Distribution'"
            write (file_id,*)  "plot '"//trim(root_path)//"initial_phasespace.dat' using 2:3 smooth kdensity, vel_pdf(x)"
            write (file_id,*)  "set title 'Spacial Distribution'"
            write (file_id,*)  "plot '"//trim(root_path)//"initial_phasespace.dat' using 1:3 smooth kdensity"
            write (file_id,*) "unset multiplot"


            write (file_id, *) "set term x11 3"
            write (file_id, *) "set autoscale x; set autoscale y"
            write (file_id,*)  "set title 'Thermal velocity estimate'"
            write (file_id,*)  "set xlabel 'time'"
            write (file_id,*)  "set ylabel 'v_th'"
            write (file_id,*)  "plot '"//filename//".dat' using 1:5 with lines"

            write (file_id, *) "set term x11 4"
            write (file_id, *) "set autoscale x; set autoscale y; unset logscale y"
            write (file_id,*)  "set title 'Absolute Energies'"
            write (file_id,*)  "set xlabel 'time'"
            write (file_id,*)  "set ylabel 'energy'"
            write (file_id,*)  "plot '"//filename//".dat' using 1:($2+$3) with lines\"
            write (file_id,*)  "title 'total (kinetic+electrostatic)',\"
            write (file_id,*)  "'"//filename//".dat' using 1:2 with lines \"
            write (file_id,*)  "title 'kinetic'"


            write (file_id, *) "set term x11 5"
            write (file_id, *) "set autoscale x; set logscale y"
            write (file_id,*)  "set title 'Relative Energy Error'"
            write (file_id,*)  "set xlabel 'time'"
            write (file_id,*)  "set ylabel 'rel. error'"
            write (file_id,*)  "plot '"//filename//"-errors.dat' using 1:2 with lines"

            close(file_id)
        endif
    endsubroutine

    !    !This function should not have any side-effects
    !    function sll_pic1d_ensure_periodicity( particle_position, &
        !                    interval_a, interval_b) result( particle_position_out)
    !        sll_real, intent(in)::interval_a, interval_b
    !        sll_real, intent(in) :: particle_position
    !        sll_real :: interval_length
    !        sll_real :: particle_position_out
    !        SLL_ASSERT(interval_a < interval_b)
    !
    !        interval_length=interval_b - interval_a
    !        particle_position_out=modulo(particle_position- interval_a, interval_length)&
        !                            +interval_a
    !
    !
    !        SLL_ASSERT(particle_position>=interval_a)
    !        SLL_ASSERT(particle_position<=interval_b)
    !    endfunction

    !> \brief Check whether timestepwidth is accurate and scale down, when particles are to fast
    subroutine sll_pic1d_check_timescale( particle_position, &
            interval_a, interval_b, timestepwidth)
        sll_real64, intent(in)::interval_a, interval_b
        sll_real64, dimension(:) ,intent(in) :: particle_position
        sll_real64, intent(inout) :: timestepwidth
        sll_real64 :: max_offset

        !Calculate maximum relative offset
        max_offset = maxval(max( abs(particle_position- interval_a),  abs(particle_position- interval_b)))&
            /(interval_b-interval_a)

        if (1.0_f64/max_offset> 2_f64) then
            print *, "Maximal offset is:",max_offset
            timestepwidth=timestepwidth*0.5_f64
            print *, "Scaling down timestepwidth to:", timestepwidth
        endif
    endsubroutine


    subroutine sll_pic1d_adjustweights(xold, xnew, vold, vnew)
        sll_real64, dimension(:),intent(in) ::vold, xold
        sll_real64, dimension(:),intent(out) ::vnew, xnew

        !Adjust weights
        if (enable_deltaf .eqv. .TRUE.) then
            particleweight = sll_pic_1d_landaudamp_PDFxv(particleposition, vnew)/ sll_pic_1d_landaudamp_PDFxv(particleposition, vold)
            particleweight = particleweight_constant - (particleweight_constant)
        endif


    endsubroutine


    !<Dummy function for Electric field
    function sll_pic_1d_electric_field(x) result(E)
        sll_real64, dimension(:), intent(in) :: x
        !sll_real64, intent(in) :: t
        sll_real64, dimension(size(x)) :: E
        E=0
        if (pic1d_testcase == SLL_PIC1D_TESTCASE_IONBEAM)         E=1.0_f64
    endfunction



endmodule  sll_pic_1d_Class


