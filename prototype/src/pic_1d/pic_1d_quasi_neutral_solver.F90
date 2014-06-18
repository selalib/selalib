
#ifndef FEM_MASS_MATRIX_TOL
#define  FEM_MASS_MATRIX_TOL 0.0000000001_f64
#endif



module sll_pic_1d_quasi_neutral_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"
    !  use sll_arbitrary_degree_spline_interpolator_1d_module
    use sll_arbitrary_degree_splines
    !!use gauss_lobatto_integration
    use sll_fft
    use sll_constants
    use sll_collective
    use sll_module_poisson_1d_periodic_solver
    use sll_poisson_1d_fem !Finite Element Bspline
    use sll_logical_meshes
    use sll_poisson_1d_fd !Finite difference solver
    use sll_poisson_1d_periodic
    use sll_particle_1d_description
    use sll_poisson_1d_fourier
    implicit none
    !initialize
    !solve
    !interpolate_solution
    !destroy
    sll_int32, parameter :: SLL_SOLVER_FEM = 1
    sll_int32, parameter :: SLL_SOLVER_FD = 2
    sll_int32, parameter :: SLL_SOLVER_SPECTRAL = 3
    sll_int32, parameter :: SLL_SOLVER_FOURIER = 4


    type  :: pic_1d_quasi_neutral_solver
        class(poisson_1d_fem), private , pointer ::  femsolver
        type(poisson_1d_periodic), pointer :: spectralsolver
        class(poisson_1d_fd), pointer :: fdsolver
        class(poisson_1d_fourier), pointer :: fouriersolver


        !This should not be here. It belongs to the solver, we just
        !keep that here because the spectral solver can't hold in its own solution
        sll_real64, dimension(:), allocatable :: poisson_solution


        type(sll_logical_mesh_1d),pointer ::  mesh

        sll_real64, private :: scalar_inhomogenity
        sll_real64, dimension(:), allocatable :: inhomogenity
        sll_real64, dimension(:), allocatable :: inhomogenity_steady
        sll_real64, dimension(:), allocatable :: inhomogenity_scndmom

        sll_comp64, dimension(:), allocatable :: inhomogenity_comp
        sll_comp64, dimension(:), allocatable :: inhomogenity_comp_steady



        sll_int32, private :: poisson_solver
        sll_int32, private :: boundary_type

        sll_int32, private :: problemsize

        sll_int32, private :: num_fourier_modes

        logical :: variance_calculation

        !Collective
        sll_int32, private :: coll_rank, coll_size
        type(sll_collective_t), pointer, private :: collective
    contains
        procedure,  pass(this) :: initialize =>pic_1d_quasi_neutral_solver_initialize
        procedure, pass(this), public ::  delete=>pic_1d_quasi_neutral_solver_delete


        !Solving
        procedure, pass(this), public ::  BC=>pic_1d_quasi_neutral_solver_boundary_conditions

        !Setting up solve
        procedure, pass(this), public ::  set_ions_constant=>pic_1d_quasi_neutral_solver_set_ions_constant
        procedure, pass(this), public :: set_ions_constant_particles=>pic_1d_quasi_neutral_solver_set_ions_constant_particles
        procedure, pass(this), public :: set_ions_constant_function=>pic_1d_quasi_neutral_solver_set_ions_constant_function
        procedure, pass(this), public :: calc_variance_rhs=>pic_1d_quasi_neutral_solver_calc_variance_rhs

        procedure,  pass(this) :: solve =>pic_1d_quasi_neutral_solver_solve

        procedure, pass(this), public :: reset_particles=>pic_1d_quasi_neutral_solver_reset_particles
        procedure, pass(this), public :: add_species=>pic_1d_quasi_neutral_solver_add_species

        procedure, pass(this), public :: set_electrons_only_weighted=>pic_1d_quasi_neutral_solver_set_electrons_only_weighted

        procedure, pass(this), public :: set_charged_allparticles_weighted=>&
            pic_1d_quasi_neutral_solver_set_charged_allparticles_weighted


        procedure, pass(this), public :: add_particles_weighted=>&
            pic_1d_quasi_neutral_solver_add_particles_weighted



        procedure, pass(this), public :: fieldenergy=>pic_1d_quasi_neutral_solver_fieldenergy
        !procedure,  pass(this) :: initialize =>pic_1d_quasi_neutral_solver_initialize

        procedure, pass(this), public :: get_problemsize=>pic_1d_quasi_neutral_solver_get_problemsize

        procedure, pass(this), public :: set_bg_particles=>pic_1d_quasi_neutral_solver_set_background_particles

        procedure, pass(this), public :: evalE=>pic_1d_quasi_neutral_solver_eval_electricfield
        procedure, pass(this), public :: evalPhi=>pic_1d_quasi_neutral_solver_eval_potential

        procedure, pass(this), public :: getPotential=>pic_1d_quasi_neutral_solver_get_potential

        !Interpolatin schemes
        procedure, pass(this), private :: get_rhs_cic=>pic_1d_quasi_neutral_solver_get_rhs_cic
    end type pic_1d_quasi_neutral_solver


    integer, private :: ierr !!!FIX THIS , DECIDE THIS



    !    sll_int, private :: n_knots
    !    sll_int, private :: n_cells
    !    sll_int, protected :: n_steadyparticles
    !    type(arbitrary_degree_spline_1d), pointer, private :: bspline_arbitrary_degree
    !
    !    sll_int, protected :: n_particles !number of particles
    !    sll_real64, dimension(:), allocatable, protected :: fem_inhomogenity_steady
    !    sll_real64, dimension(:), allocatable, private :: fem_inhomogenity
    !    !sll_real64, dimension(:), allocatable, private :: particleweights
    !
    !    sll_real64, private :: scale_matrix_equation    !<Scale for the stiffnes matrix and the inhomogenity
    !
    !    !Boundary Description handeled by the solver
    !    sll_int32, private :: sll_bspline_fem_solver_boundary_type = SLL_PERIODIC
    !
    !
    !    !Poisson Finite Differences Solver
    !    type (poisson_1d_periodic),pointer,private           :: poissonsolverFD => null()
    !



    !    interface bspline_fem_solver_1d_cell_number
    !        module procedure         bspline_fem_solver_1d_cell_number_single
    !        module procedure     bspline_fem_solver_1d_cell_number_array
    !    endinterface


    !  integer, private :: startsbefore=0

contains

    function new_pic_1d_quasi_neutral_solver(eta_min, eta_max, &
            spline_degree, num_cells, poisson_solver_type, collective ,boundary_type ) &
            result(qn_solver)
        sll_int32, intent(in):: spline_degree
        sll_int32, intent(in)::  poisson_solver_type
        sll_int32, intent(in)::  num_cells
        sll_real64, intent(in) :: eta_min, eta_max
        sll_int32, intent(in) :: boundary_type
        type(sll_collective_t), pointer , intent(in):: collective

        class(pic_1d_quasi_neutral_solver), pointer :: qn_solver

        SLL_ALLOCATE(qn_solver,ierr)

        call   pic_1d_quasi_neutral_solver_initialize( qn_solver, eta_min, eta_max, &
            spline_degree, num_cells, poisson_solver_type, collective,boundary_type  )
    endfunction

    subroutine pic_1d_quasi_neutral_solver_initialize( this, eta_min, eta_max, &
            spline_degree, num_cells, poisson_solver_type, collective ,boundary_type )
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        sll_int32, intent(in):: spline_degree
        sll_int32, intent(in)::  poisson_solver_type
        sll_int32, intent(in)::  num_cells
        sll_real64, intent(in) :: eta_min, eta_max
        sll_int32, intent(in) :: boundary_type

        type(sll_logical_mesh_1d), pointer :: mesh
        type(sll_collective_t), pointer , intent(in):: collective

        !Pepare collective
        this%collective=>collective
        this%coll_size=sll_get_collective_size( collective )
        this%coll_rank=sll_get_collective_rank( collective )

        !Boundary Conditions
        this%boundary_type=boundary_type

        !Scale the right side correctly
        this%scalar_inhomogenity=(eta_max -eta_min)
        this%poisson_solver=poisson_solver_type
        this%variance_calculation=.false.




        !Generate logical mesh
        mesh=>new_logical_mesh_1d( num_cells, eta_min, eta_max )
        this%mesh=>mesh
        selectcase (poisson_solver_type)
            case(SLL_SOLVER_FEM)
                this%femsolver=>new_poisson_1d_fem(this%mesh, spline_degree,this%boundary_type, ierr)
                this%problemsize=num_cells

            case(SLL_SOLVER_FD)
                this%problemsize=num_cells
                this%fdsolver=>new(this%mesh,1,this%boundary_type,ierr)
            case(SLL_SOLVER_FOURIER)
                this%num_fourier_modes=spline_degree
                this%problemsize=this%num_fourier_modes
                this%fouriersolver=>new_poisson_1d_fourier(this%mesh,   this%num_fourier_modes,this%boundary_type, ierr)
                SLL_CLEAR_ALLOCATE(this%inhomogenity_comp(this%num_fourier_modes),ierr)
                SLL_CLEAR_ALLOCATE(this%inhomogenity_comp_steady(this%num_fourier_modes),ierr)
            case(SLL_SOLVER_SPECTRAL)
                !All Knots, for boundary set first and last knot to zero
                this%problemsize=num_cells+1 !FIXED
                SLL_CLEAR_ALLOCATE(this%poisson_solution(this%problemsize),ierr)
                this%spectralsolver=>new(eta_min,eta_max, num_cells,ierr)
        endselect

        SLL_CLEAR_ALLOCATE(this%inhomogenity(this%problemsize),ierr)
        SLL_CLEAR_ALLOCATE(this%inhomogenity_steady(this%problemsize),ierr)
        !if (this%variance_calculation .eqv. .true.) then
        SLL_CLEAR_ALLOCATE(this%inhomogenity_scndmom(this%problemsize),ierr)
        !endif


    endsubroutine

    subroutine pic_1d_quasi_neutral_solver_delete(this)
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        !        print *, this%inhomogenity
        !!SLL_DEALLOCATE_ARRAY(this%inhomogenity_steady, ierr  )
        !!SLL_DEALLOCATE_ARRAY(this%inhomogenity, ierr  )
        !!if (allocated(this%poisson_solution) SLL_DEALLOCATE_ARRAY(this%poisson_solution,ierr)

        if (associated(this%femsolver))  then
            !call delete(this%mesh)
            !call this%femsolver%delete(ierr)
        endif
        !SLL_DEALLOCATE(this,ierr)


    endsubroutine

    function pic_1d_quasi_neutral_solver_get_problemsize(this) &
            result(problemsize)
        class(pic_1d_quasi_neutral_solver), intent(in) :: this
        sll_int32 :: problemsize
        problemsize=this%problemsize
    endfunction

    !
    !    function pic_1d_quasi_neutral_solver_fourier_rhs(this,  ppos, pweight  )
    !        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
    !        sll_real64, dimension(:) ,intent(in) :: ppos
    !        sll_real64, dimension(:) ,intent(in) :: pweight
    !        SLL_ASSERT(size(ppos)==size(pweight))
    !
    !
    !
    !
    !    endfunction
    function pic_1d_quasi_neutral_solver_calc_variance_rhs(this)&
            result(variance)
        class(pic_1d_quasi_neutral_solver), intent(in) :: this
        sll_real64, dimension(this%problemsize) :: variance_v
        sll_real64 :: variance

        selectcase(this%poisson_solver)
            case(SLL_SOLVER_FEM)
                !Biased estimate
                variance_v = this%inhomogenity_scndmom !  - this%inhomogenity**2

                variance=sum(variance_v)
            case default
                variance=-1.0_f64
        endselect

    endfunction



    subroutine pic_1d_quasi_neutral_solver_solve(this)
        class(pic_1d_quasi_neutral_solver) , intent(inout) :: this

        !        if  (.NOT. allocated(this%inhomogenity))  then
        !            SLL_CLEAR_ALLOCATE(this%inhomogenity(this%problemsize),ierr)
        !        endif

        selectcase(this%poisson_solver)
            case(SLL_SOLVER_FEM)
                !Here one would also have the option to solve on each core and
                !omit the broadcast
                call sll_collective_globalsum(this%collective, this%inhomogenity, 0)
                this%inhomogenity=this%scalar_inhomogenity*(this%inhomogenity_steady+this%inhomogenity)
                !On machines with more than 6 cores, we solve on each core
                !We do this, because a HexaCore is the biggest PC used in development
                if (this%coll_rank==0 .OR. this%coll_size>6) then

                    call this%femsolver%solve(this%inhomogenity)

                endif
                ! call sll_collective_barrier(sll_world_collective)
                !Broadcast result so it lies on every core for interpolation
                if (this%coll_size<=6) then
                    call sll_collective_bcast_real64(this%collective, this%femsolver%fem_solution,&
                        this%problemsize, 0 )

                endif
            case(SLL_SOLVER_FD)
                call sll_collective_globalsum(this%collective, this%inhomogenity, 0)
                this%inhomogenity=this%scalar_inhomogenity*(this%inhomogenity_steady+this%inhomogenity)
                !On machines with more than 6 cores, we solve on each core
                !We do this, because a HexaCore is the biggest PC used in development
                if (this%coll_rank==0 .OR. this%coll_size>6) then
                    call this%fdsolver%solve(this%inhomogenity)
                endif
                ! call sll_collective_barrier(sll_world_collective)
                !Broadcast result so it lies on every core for interpolation
                if (this%coll_size<=6) then
                    call sll_collective_bcast_real64(this%collective, this%fdsolver%fd_solution,&
                        this%problemsize, 0 )

                endif
            case(SLL_SOLVER_SPECTRAL)
                call sll_collective_globalsum(this%collective, this%inhomogenity, 0)
                this%inhomogenity=this%scalar_inhomogenity*(this%inhomogenity_steady+this%inhomogenity)
                call solve(this%spectralsolver, this%poisson_solution, this%inhomogenity)
            case(SLL_SOLVER_FOURIER)
                !Distribute fourier modes on each Core
                call sll_collective_globalsum(this%collective, this%inhomogenity_comp)
                call this%fouriersolver%solve(this%inhomogenity_comp)
        endselect

    endsubroutine

    !<@brief Takes negative particles with charge one and custom weights
    !>@param this pointer to a pic_1d_quasi_neutral_solver object.
    !>@param ppos spatial position of the negative ion/electron
    !>@param pweight corresponding weight of the negative ion/electron
    subroutine pic_1d_quasi_neutral_solver_set_electrons_only_weighted(this,  ppos, pweight  )
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        sll_real64, dimension(:) ,intent(in) :: ppos
        sll_real64, dimension(:) ,intent(in) :: pweight
        SLL_ASSERT(size(ppos)==size(pweight))

        selectcase(this%poisson_solver)
            case(SLL_SOLVER_FEM)
                this%inhomogenity=-this%femsolver%get_rhs_from_klimontovich_density_weighted( &
                    this%BC(ppos),pweight)
                if (this%variance_calculation .eqv. .true.) then
                    !this%inhomogenity will be overwritten later
                    this%inhomogenity_scndmom=abs(this%femsolver%get_rhs_from_klimontovich_density_moment( &
                        this%BC(ppos),pweight,2)) - this%inhomogenity**2

                endif
            case(SLL_SOLVER_FD)
                this%inhomogenity=-this%get_rhs_cic( this%BC(ppos),pweight)
            case(SLL_SOLVER_SPECTRAL)
                this%inhomogenity=-this%get_rhs_cic( this%BC(ppos),pweight)
            case(SLL_SOLVER_FOURIER)
                this%inhomogenity_comp=this%fouriersolver%get_rhs_from_klimontovich_density_weighted(&
                    this%BC(ppos),-pweight)

        endselect

    endsubroutine


    !<@brief Takes negative and positive particles in one array, charge has to be put in weights by user
    !>@param this pointer to a pic_1d_quasi_neutral_solver object.
    !>@param ppos spatial position of the negative ion/electron
    !>@param pweight corresponding weight of particle
    subroutine pic_1d_quasi_neutral_solver_set_charged_allparticles_weighted(this,&
            ppos, pweight)
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        sll_real64, dimension(:),intent(in) :: ppos
        sll_real64, dimension(:),intent(in)  :: pweight
        SLL_ASSERT(size(ppos)==size(pweight))
        this%inhomogenity=0
        call pic_1d_quasi_neutral_solver_add_particles_weighted(this,&
            ppos, pweight)
    endsubroutine


    !<@brief Takes negative and positive particles in one array, charge has to be put in weights by user
    !>@param this pointer to a pic_1d_quasi_neutral_solver object.
    !>@param ppos spatial position of the negative ion/electron
    !>@param pweight corresponding weight of particle
    subroutine pic_1d_quasi_neutral_solver_add_particles_weighted(this,&
            ppos, pweight)
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        sll_real64, dimension(:),intent(in) :: ppos
        sll_real64, dimension(:),intent(in)  :: pweight
        SLL_ASSERT(size(ppos)==size(pweight))

        selectcase(this%poisson_solver)
            case(SLL_SOLVER_FEM)
                this%inhomogenity= this%inhomogenity+ this%femsolver%get_rhs_from_klimontovich_density_weighted(&
                    this%BC(ppos),pweight)
                if (this%variance_calculation .eqv. .true.) then
                    this%inhomogenity_scndmom=  this%inhomogenity_scndmom+ this%femsolver%&
                        get_rhs_from_klimontovich_density_weighted(this%BC(ppos),pweight)
                endif
            case(SLL_SOLVER_FD)
                this%inhomogenity=  this%inhomogenity +this%get_rhs_cic( this%BC(ppos),pweight)
            case(SLL_SOLVER_SPECTRAL)
                this%inhomogenity=    this%inhomogenity +this%get_rhs_cic( this%BC(ppos),pweight)
            case(SLL_SOLVER_FOURIER)
                this%inhomogenity_comp=this%inhomogenity_comp +this%fouriersolver%get_rhs_from_klimontovich_density_weighted(&
                    this%BC(ppos),pweight)

        endselect
    endsubroutine





    subroutine pic_1d_quasi_neutral_solver_reset_particles(this)
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this

        selectcase(this%poisson_solver)
            case(SLL_SOLVER_FOURIER)
                this%inhomogenity_comp=0.0_f64
                    !!!TODO variance calculation
            case default
                this%inhomogenity=0.0_f64
                if (this%variance_calculation .eqv. .true.) then
                    this%inhomogenity_scndmom=0
                endif
        endselect
    endsubroutine

    !<@brief Takes negative and positive particles in one array, charge has to be put in weights by user
    !>@param this pointer to a pic_1d_quasi_neutral_solver object.
    subroutine pic_1d_quasi_neutral_solver_add_species(this,&
            species)
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        type(sll_particle_1d_group), intent(in) :: species

        selectcase(this%poisson_solver)
            case(SLL_SOLVER_FEM)
                this%inhomogenity=this%inhomogenity + this%femsolver%get_rhs_from_klimontovich_density_weighted(&
                    this%BC(species%particle%dx), sign(species%particle%weight,species%qm) )

                if (this%variance_calculation .eqv. .true.) then
                    this%inhomogenity_scndmom=this%inhomogenity_scndmom + &
                        this%femsolver%get_rhs_from_klimontovich_density_weighted( &
                        this%BC(species%particle%dx),sign(species%particle%weight,species%qm) )
                endif
            case(SLL_SOLVER_FD)
                this%inhomogenity=this%inhomogenity + &
                    this%get_rhs_cic( this%BC(species%particle%dx),sign(species%particle%weight,species%qm)  )
            case(SLL_SOLVER_SPECTRAL)
                this%inhomogenity=this%inhomogenity + &
                    this%get_rhs_cic( this%BC(species%particle%dx),sign(species%particle%weight,species%qm)  )
            case(SLL_SOLVER_FOURIER)
                this%inhomogenity_comp=this%inhomogenity_comp+ &
                this%fouriersolver%get_rhs_from_klimontovich_density_weighted(&
                this%BC(species%particle%dx), sign(species%particle%weight,species%qm) )
                print  *,  this%inhomogenity_comp

        endselect
    endsubroutine

    !<@brief Takes negative and positive particles with charge one and custom weights
    !>@param this pointer to a pic_1d_quasi_neutral_solver object.
    !>@param ppos_neg spatial position of the negative ion/electron
    !>@param pweight_neg corresponding weight of the negative ion/electron
    !>@param ppos_pos spatial position of the positive ion
    !>@param pweight_pos corresponding weight of the positive ion
    subroutine pic_1d_quasi_neutral_solver_set_posneg_particles_weighted(this,&
            ppos_pos, pweight_pos,  ppos_neg, pweight_neg )
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        sll_real64, dimension(:),intent(in) :: ppos_pos
        sll_real64, dimension(:),intent(in)  :: pweight_pos
        sll_real64, dimension(:) ,intent(in) :: ppos_neg
        sll_real64, dimension(:),intent(in)  :: pweight_neg
        sll_real64, dimension(size(ppos_pos)*2):: ppos
        sll_real64, dimension(size(pweight_pos)*2):: pweight
        sll_int32 :: npart

        npart=size(ppos_pos)
        SLL_ASSERT(size(ppos_neg)==size(pweight_neg))
        SLL_ASSERT(size(ppos_pos)==size(pweight_pos))

        !!!This can be solved differently
        ppos(1:npart)=ppos_neg
        pweight(1:npart)=pweight_pos
        ppos(npart+1:2*npart)=ppos_neg
        pweight(npart+1:2*npart)=-pweight_neg
        call this%set_charged_allparticles_weighted(ppos,pweight)

    endsubroutine

    subroutine pic_1d_quasi_neutral_solver_set_ions_constant_particles(this,  ppos)
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        sll_real64, dimension(:),intent(in)  :: ppos
        sll_int32 :: np
        np=size(ppos)

        this%inhomogenity_steady=&
            this%femsolver%get_rhs_from_klimontovich_density(this%BC(ppos)) / (np*this%coll_size  )
        call sll_collective_globalsum_array_real64(this%collective,this%inhomogenity_steady )
    endsubroutine

    !<Set the constant background, that will be added to each poisson solve
    subroutine pic_1d_quasi_neutral_solver_set_background_particles(this, ppos, pweight)
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        sll_real64, dimension(:),intent(in)  :: ppos  !< Particle Position
        sll_real64, dimension(:),intent(in)  :: pweight !< Particleweight
        SLL_ASSERT(size(ppos)==size(pweight))


        selectcase(this%poisson_solver)
        case(SLL_SOLVER_FOURIER)
                this%inhomogenity_comp_steady=this%inhomogenity_comp_steady+ &
                this%fouriersolver%get_rhs_from_klimontovich_density_weighted&
                        (this%BC(ppos), pweight  )

             call sll_collective_globalsum_array_comp64(this%collective,this%inhomogenity_comp_steady)
        case(SLL_SOLVER_FEM)
            this%inhomogenity_steady= this%femsolver%get_rhs_from_klimontovich_density_weighted&
            (this%BC(ppos), pweight  )
             call sll_collective_globalsum_array_real64(this%collective,this%inhomogenity_steady )

        case default

        endselect



    endsubroutine

    subroutine pic_1d_quasi_neutral_solver_set_ions_constant(this, const)
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        sll_real64, intent(in) :: const

         selectcase(this%poisson_solver)
        case(SLL_SOLVER_FOURIER)
            this%inhomogenity_comp_steady=0.0_f64
        case default
        this%inhomogenity_steady=const

        endselect
    endsubroutine

    subroutine pic_1d_quasi_neutral_solver_set_ions_constant_function(this,  distribution_function)
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        procedure (poisson_1d_fem_rhs_function) :: distribution_function
        this%inhomogenity_steady=&
            this%femsolver%get_rhs_from_function( distribution_function,10)
    endsubroutine

    subroutine pic_1d_quasi_neutral_solver_eval_electricfield(this,eval_points,eval_solution)

        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        sll_real64, dimension(:), intent(in)     :: eval_points
        sll_real64, dimension(:), intent(out)     :: eval_solution

        SLL_ASSERT(size(eval_points)==size(eval_solution))

        selectcase (this%poisson_solver)
            case(SLL_SOLVER_FEM)
                call this%femsolver%eval_solution_derivative&
                    (this%BC(eval_points) ,eval_solution)
            case(SLL_SOLVER_FD)
                call this%fdsolver%eval_solution_derivative&
                    (this%BC(eval_points) ,eval_solution)
            case(SLL_SOLVER_SPECTRAL)
                !!!not implemented by spectral solver
            case(SLL_SOLVER_FOURIER)
                call this%fouriersolver%eval_solution_derivative&
                    (this%BC(eval_points) ,eval_solution)

        endselect
    endsubroutine

    subroutine pic_1d_quasi_neutral_solver_eval_potential(this,eval_points,eval_solution)
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        sll_real64, dimension(:), intent(in)     :: eval_points
        sll_real64, dimension(:), intent(out)     :: eval_solution

        SLL_ASSERT(size(eval_points)==size(eval_solution))

        selectcase (this%poisson_solver)
            case(SLL_SOLVER_FEM)
                call this%femsolver%eval_solution&
                    (this%BC(eval_points) ,eval_solution)
            case(SLL_SOLVER_FD)
                call this%fdsolver%eval_solution&
                    (this%BC(eval_points) ,eval_solution)
            case(SLL_SOLVER_SPECTRAL)
                !!!not implemented by spectral solver
            case(SLL_SOLVER_FOURIER)
                    call this%fouriersolver%eval_solution&
                    (this%BC(eval_points) ,eval_solution)
        endselect
    endsubroutine





    function pic_1d_quasi_neutral_solver_get_potential(this) &
            result(potential)
        class(pic_1d_quasi_neutral_solver), intent(inout) :: this
        sll_real64, dimension(this%problemsize) :: potential

        selectcase (this%poisson_solver)
            case(SLL_SOLVER_FEM)
                potential=this%femsolver%fem_solution
            case(SLL_SOLVER_FD)
                potential=this%fdsolver%fd_solution
        endselect
    endfunction

    function  pic_1d_quasi_neutral_solver_get_solution(this) result(solution)
        class(pic_1d_quasi_neutral_solver), intent(in) :: this

        sll_real64, dimension(this%problemsize) :: solution

        solution=this%femsolver%fem_solution
    endfunction

    function  pic_1d_quasi_neutral_solver_get_inhomogenity(this) result(inhomogenity)
        class(pic_1d_quasi_neutral_solver), intent(in) :: this
        sll_real64, dimension(this%problemsize) :: inhomogenity

        inhomogenity=this%inhomogenity
    endfunction

    function pic_1d_quasi_neutral_solver_boundary_conditions(this, x)  result(xout)
        class(pic_1d_quasi_neutral_solver) , intent(in) :: this
        sll_real64, dimension(:), intent(in) ::x
        sll_real64, dimension(size(x))  :: xout
        sll_real64 :: interval_a, interval_b,interval_length

        interval_a=this%mesh%eta_min
        interval_b=this%mesh%eta_max
        interval_length=interval_b-interval_a
        SLL_ASSERT(interval_a < interval_b)

        selectcase(this%boundary_type)
            case(SLL_PERIODIC)
                xout=x-interval_a

                do while (minval(xout)<0.0_f64)
                    xout=xout +interval_length
                enddo
                xout=mod(xout, interval_length)&
                    +interval_a
                SLL_ASSERT(minval(xout)>=interval_a)
                SLL_ASSERT(maxval(xout)<interval_b)
            case(SLL_DIRICHLET)
                xout=x
                where (x<interval_a) xout=interval_a
                where (x>interval_b) xout=interval_b



        endselect


        SLL_ASSERT(minval(xout)>=interval_a)
        SLL_ASSERT(maxval(xout)<=interval_b)
    endfunction


    function  pic_1d_quasi_neutral_solver_fieldenergy(this) &
            result(energy)
        class(pic_1d_quasi_neutral_solver) , intent(inout) :: this
        sll_real64 :: energy
        energy=0
        selectcase(this%poisson_solver)
            case(SLL_SOLVER_FEM)
                energy=this%femsolver%H1seminorm_solution()
                energy=0.5_f64*energy/(this%mesh%eta_max-this%mesh%eta_min)
            case(SLL_SOLVER_FD)
                energy=this%fdsolver%H1seminorm_solution()
                energy=0.5_f64*energy/(this%mesh%eta_max-this%mesh%eta_min)
             case(SLL_SOLVER_FOURIER)
                energy=this%fouriersolver%H1seminorm_solution()
                energy=0.5_f64*energy/(this%mesh%eta_max-this%mesh%eta_min)
        endselect

    endfunction



    !< Cloud in Cell scheme for logical mesh
    function  pic_1d_quasi_neutral_solver_get_rhs_cic( this, ppos, pweight) result(rhs)
        class(pic_1d_quasi_neutral_solver),intent(inout) :: this
        sll_real64, dimension(:), intent(in) ::ppos
        sll_real64, dimension(:), intent(in) ::pweight
        sll_real64, dimension(this%problemsize) :: rhs
        sll_real64, dimension(size(ppos)) ::   pw
        sll_int32, dimension(size(ppos)) ::   pidx
        sll_real64, dimension(this%mesh%num_cells +1) ::   knotsvals
        sll_real64 :: eta_min, eta_max
        sll_int32 :: num_cells
        sll_int32 :: idx

        eta_min=this%mesh%eta_min
        eta_max=this%mesh%eta_max
        num_cells=this%mesh%num_cells


        pw=(ppos-eta_min)*(num_cells)/(eta_max-eta_min)

        knotsvals=0

        !Left side
        pidx=floor(pw)
        pidx=pidx+1  !first knot has index 1

        SLL_ASSERT(maxval(pidx)<=num_cells+1)
        SLL_ASSERT(minval(pidx)>=1)


        !knotsvals(pidx)=knotsvals(pidx)+ (1.0_f64-mod(pw,1.0_f64))*pweight
        do idx=1,size(pidx)
            knotsvals(pidx(idx))=knotsvals(pidx(idx))+ (1.0_f64-mod(pw(idx),1.0_f64))*pweight(idx)
        enddo

        !Right side
        where(pidx==num_cells +1)
            pw=0.0_f64
            pidx=num_cells
        endwhere
        pidx=pidx+1
        SLL_ASSERT(maxval(pidx)<=num_cells+1)
        SLL_ASSERT(minval(pidx)>=1)

        do idx=1,size(pidx)
            knotsvals(pidx(idx))=knotsvals(pidx(idx))+ mod(pw(idx),1.0_f64)*pweight(idx)
        enddo

        selectcase(this%boundary_type )
            case(SLL_PERIODIC)
                !Last knot is identical with the first knot
                SLL_ASSERT(this%problemsize==num_cells)
                rhs=knotsvals(1:num_cells)
                rhs(1)=rhs(1)+knotsvals(num_cells+1)

            case(SLL_DIRICHLET)
                SLL_ASSERT(this%problemsize==num_cells)
                rhs=knotsvals(1:num_cells)
                rhs(1)=0.0_f64 !homogenous dirichlet
                !works only because we have first order scheme
        endselect
    endfunction



    !
    !    subroutine sll_bspline_fem_solver_1d_set_weights(weights_user )
    !        sll_real64, dimension(:), allocatable, intent(in):: weights_user
    !        if (.NOT. allocated(weights)) then
    !            SLL_CLEAR_ALLOCATE(weights(1:size(weights_user)),ierr)
    !        endif
    !        weights=weights_user
    !    endsubroutine

    !    subroutine sll_bspline_fem_solver_1d_get_weights(weights_user )
    !        sll_real64, dimension(:), allocatable, intent(out):: weights_user
    !        weights_user=weights
    !    endsubroutine




    !<Solves the given FEM-Matrix, with right side as sum of deltafunctions at particleposition
    !<with sign (-1), so particleposition are the positions of the negative charges
    !<If none given solves the
    !< solves  - \laplace u = f See test case or more help
    !<
    !    subroutine sll_bspline_fem_solver_1d_solve( particleposition)
    !        sll_real64, DIMENSION(:), intent(in),optional :: particleposition
    !        !sll_real64, DIMENSION(:), intent(in),optional   :: particleweight_user
    !        !SLL_ALLOCATE(fem_inhomogenity(n_cells),ierr)
    !
    !
    !
    !        selectcase(sll_pic1d_poisson_solver)
    !            case(SLL_SOLVER_FEM)
    !                if (associated(collective)) then
    !                    fem_inhomogenity= interpolate_particles_bsplines( bspline_arbitrary_degree, knots_mesh,&
        !                        particleposition, b_splines_at_x,weights )
    !
    !                    !Here one would also have the option to solve on each core and
    !                    !omit the broadcast
    !                    call sll_collective_globalsum(collective, fem_inhomogenity, 0)
    !
    !                    if (collective_rank==0) then
    !                        !fem_inhomogenity=scalar_inhomogenity *( fem_inhomogenity_steady-fem_inhomogenity)
    !                        !print *, "Inhom:",  sum(fem_inhomogenity)
    !
    !                        !fem_inhomogenity=scalar_inhomogenity*(1.0_f64/n_cells-(fem_inhomogenity/fem_inhomogenity_steady))
    !                        fem_inhomogenity=scalar_inhomogenity*(fem_inhomogenity_steady-fem_inhomogenity)
    !                        !   print *,  sum(fem_inhomogenity)
    !                        !fem_inhomogenity=scalar_inhomogenity*(1.0_f64/n_cells-fem_inhomogenity&
        !                            !           *(knots_mesh(n_knots)-knots_mesh(1)) )
    !                        !Delta f, where fem_inhomogenity_steady is delta f
    !                        !                fem_inhomogenity=scalar_inhomogenity*(1.0_f64/n_cells - &
        !                            !                                (fem_inhomogenity-fem_inhomogenity_steady)*n_cells)
    !                        !                fem_inhomogenity=scalar_inhomogenity*(fem_inhomogenity*n_cells -1.0_f64/n_cells)
    !
    !                        fem_solution=bspline_fem_solver_1d_solve_matrix_equation(fem_inhomogenity)
    !                    endif
    !                    ! call sll_collective_barrier(sll_world_collective)
    !                    !Broadcast result so it lies on every core for interpolation
    !                    call sll_collective_bcast_real64( collective, fem_solution, n_cells, 0 )
    !                else
    !
    !                    if (present(particleposition)) then
    !                        !Calculate Electric Potential (gather)
    !                        !> \frac{q}{\epsilon_0}(n_i - n_e)
    !                        if (size(particleposition)==0) then
    !                            fem_inhomogenity=0
    !                        else
    !                            fem_inhomogenity=scalar_inhomogenity   &
        !                                *( fem_inhomogenity_steady - &
        !                                interpolate_particles_bsplines( bspline_arbitrary_degree, knots_mesh,&
        !                                particleposition, b_splines_at_x,weights ))
    !                        endif
    !                    else
    !                        fem_inhomogenity=scalar_inhomogenity*fem_inhomogenity_steady
    !                    endif
    !                    !Solve
    !                    fem_solution=bspline_fem_solver_1d_solve_matrix_equation(fem_inhomogenity)
    !
    !                endif
    !            case(SLL_SOLVER_FD)
    !                if (present(particleposition)) then
    !                    fem_inhomogenity=sll_cloudincell_1d(particleposition,weights)
    !                    if (associated(collective)) call sll_collective_globalsum(collective, fem_inhomogenity, 0)
    !                endif
    !
    !                if ( ( .NOT. associated(collective) ) .OR. collective_rank==0) then
    !                    fem_inhomogenity=scalar_inhomogenity*(fem_inhomogenity_steady-fem_inhomogenity)
    !
    !                    call solve(poissonsolverFD,fem_solution, fem_inhomogenity)
    !                endif
    !
    !                if (associated(collective)) call sll_collective_bcast_real64( collective, fem_solution, n_knots, 0 )
    !        endselect
    !
    !
    !
    !        !fem_solution=fem_solution-sum(fem_solution)/n_cells
    !        !print *, sum(fem_solution), sum (fem_inhomogenity)
    !        !SLL_DEALLOCATE_ARRAY(fem_inhomogenity, ierr)
    !    endsubroutine



    !<Linear interpolator for particles, which corresponds to a cloud in cell scheme in 1d
    !    function sll_cloudincell_1d(ppos, pweight) result(rho)
    !        sll_real64, dimension(:), intent(in) :: ppos
    !        sll_real64, dimension(:), intent(in) :: pweight
    !        sll_real64, dimension(n_knots)  :: rho
    !        sll_int32, dimension(size(ppos)) :: cell
    !        !!SLL_ASSERT(size(pweight)==size(ppos))
    !        rho=0
    !        cell=bspline_fem_solver_1d_cell_number(knots_mesh,ppos)
    !        rho(cell)=(ppos-knots_mesh(cell+1)  )/(knots_mesh(cell+1)-knots_mesh(cell))*pweight
    !        rho(cell+1)=(knots_mesh(cell)-ppos)/(knots_mesh(cell+1)-knots_mesh(cell))*pweight
    !
    !
    !        !periodic
    !        !rho=rho(1:end-1)
    !    endfunction
    !





    !    !<Mostly used to set the heavy ion contribution to a constant level
    !    subroutine sll_bspline_fem_solver_1d_set_inhomogenity_constant( average )
    !        sll_real64, intent(in) :: average
    !        fem_inhomogenity_steady=average
    !    endsubroutine
    !
    !
    !
    !    !<First Cell is in analogy to fortran arrays Cell Number 1
    !    function  bspline_fem_solver_1d_cell_number_array(knots_mesh, particle_position) result(cell)
    !        sll_real64, dimension(:), intent(in)     :: knots_mesh
    !        sll_real64, dimension(:),  intent(in) ::particle_position
    !        sll_int32, dimension(size(particle_position)) :: cell
    !
    !        cell=floor((particle_position - knots_mesh(1))/(knots_mesh(2)-knots_mesh(1)))+1
    !    endfunction
    !
    !
    !    function  bspline_fem_solver_1d_cell_number_single(knots_mesh, particle_position) result(cell)
    !        sll_real64, dimension(:), intent(in)     :: knots_mesh
    !        sll_real64, intent(in) ::particle_position
    !        sll_int32 :: cell
    !        SLL_ASSERT(size(knots_mesh)>=3)
    !        SLL_ASSERT(particle_position>=knots_mesh(1))
    !        SLL_ASSERT(particle_position<=knots_mesh(size(knots_mesh)))
    !
    !        cell=floor((particle_position - knots_mesh(1))/(knots_mesh(2)-knots_mesh(1)))+1
    !
    !
    !
    !
    !        !Periodicity Exception: last knot belongs to last cell
    !        if (cell==size(knots_mesh))  then
    !            cell=cell-1
    !            !cell=1
    !            print *, "Particle sitting on Last knot of Mesh", particle_position
    !            !stop
    !        else
    !            do while (particle_position>=knots_mesh(cell+1))
    !                cell=cell +1
    !            enddo
    !            do while (particle_position<knots_mesh(cell))
    !                cell=cell -1
    !            enddo
    !            SLL_ASSERT( particle_position < knots_mesh(cell+1))
    !        endif
    !
    !
    !        SLL_ASSERT(cell>=1)
    !        SLL_ASSERT(cell<=n_cells)
    !        SLL_ASSERT( particle_position >= knots_mesh(cell))
    !
    !    endfunction











end module
