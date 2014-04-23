program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_utilities.h"
    use sll_constants
    use sll_logical_meshes
    !use sll_poisson_1d_periodic
    !use sll_poisson_solvers
    use sll_arbitrary_degree_spline_interpolator_1d_module
    use sll_arbitrary_degree_splines
    use gauss_legendre_integration
    use pic_1d_particle_loading
    use sll_bspline_fem_solver_1d
    implicit none


    integer :: ierr
    integer :: i
    integer:: nparticles =10
    sll_real :: interval_a=-1, interval_b=9
    sll_real, DIMENSION(:), allocatable:: particleposition
    sll_real, DIMENSION(:), allocatable :: particlespeed



    !Mesh parameters
    sll_int32 :: mesh_cells = 20
    sll_real :: mesh_delta_t

    type(sll_logical_mesh_1d), pointer :: mesh1d
    !    type(poisson_1d_periodic), pointer :: poisson1dsolver


    type(arb_deg_1d_interpolator) :: ad1d
    type(arbitrary_degree_spline_1d), pointer :: arbitrarydegree1D

    sll_int32                   :: degree=3
    sll_real64, dimension(:), allocatable    :: knots
    sll_int32                  :: num_pts
    sll_int32                  :: spline_idx
    sll_real64, dimension(:), allocatable    :: mass_matrix_first_line
    sll_real64, dimension(:), allocatable    :: mass_matrix_period
    sll_real64, dimension(:), allocatable    :: bspline_knot_values

    sll_real64, dimension(:,:), allocatable    :: quadrature_point_vals
    integer                     :: quadrature_npoints
    sll_real64, dimension(:,:), allocatable   :: quadrature_points_weights
    sll_real64, dimension(:,:), allocatable    :: bspline_qpoint_values
    integer :: cell







    !ad1d=>new(interval_a, interval_b, SLL_PERIODIC , SLL_PERIODIC,
    mesh1d=> new_logical_mesh_1d(mesh_cells, interval_a,interval_b)
    call sll_display(mesh1d)

    !SLL_ASSERT(mesh_cells>1)

    mesh_delta_t = (interval_b - interval_a)/mesh_cells
    SLL_ALLOCATE(knots(mesh_cells+1),ierr)

    knots(1)=interval_a
    do i=2,size(knots)
        knots(i)=knots(i-1) + mesh_delta_t
    enddo

    !print *, knots

    arbitrarydegree1D=>new_arbitrary_degree_spline_1d( degree, knots, size(knots), SLL_PERIODIC)


    mass_matrix_period=gen_mass_matrix_period(arbitrarydegree1D, knots )



stop

    !Get Gauss Legendre points and weights to be exact for the selected spline degree
    !Note a Bspline is a piecewise polynom
    !quadrature_npoints=ceiling(0.5_f64*real(degree+1))+1
    quadrature_npoints=2*degree

    SLL_ALLOCATE(quadrature_points_weights(2,1:quadrature_npoints),ierr)
    quadrature_points_weights=gauss_legendre_points_and_weights(quadrature_npoints, knots(1), knots(degree+2))
    !print *,quadrature_points_weights
    !Get spline values at
    SLL_ALLOCATE(bspline_qpoint_values(1:degree+1,1:quadrature_npoints), ierr)

    do i=1,quadrature_npoints
        cell=floor((quadrature_points_weights(1,i ) - knots(1))/ mesh_delta_t )+1
        cell=min(cell, degree+1)
        bspline_qpoint_values(:,i)=b_splines_at_x(arbitrarydegree1D, cell, quadrature_points_weights(1,i ))
        !reorientate to the right side
        bspline_qpoint_values(:,i) = bspline_qpoint_values((degree+1):1:-1,i)
        !shift in order to account for intersecting and limited support
        bspline_qpoint_values(:,i) = eoshift(bspline_qpoint_values(:,i), (cell-1))
    enddo


call sll_display(knots(1:degree+2), "(F10.5)" )
call sll_display(quadrature_points_weights(1,:), "(F10.5)" )
call sll_display(quadrature_points_weights(2,:), "(F10.5)" )

call sll_display(bspline_qpoint_values, "(F10.5)" )

print *,dot_product(bspline_qpoint_values(1, :), quadrature_points_weights(2,:)  ), "##"

    !Do Gauss Legendre integral for each spline combination
    !Loop over first half of splines with intersecting support
    !spline index corresponds to the lowest index of the support cell
 SLL_ALLOCATE(mass_matrix_period(degree+1),ierr)
    do i=1, degree+1
        !print *,bspline_qpoint_values(1, :)*bspline_qpoint_values(i, :),"##"
        mass_matrix_period(i)=dot_product(bspline_qpoint_values(1, :)*bspline_qpoint_values(i,:), quadrature_points_weights(2,:)  )
    enddo

write(*,*)
call sll_display(mass_matrix_period, "(F10.5)" )

SLL_CLEAR_ALLOCATE(mass_matrix_first_line(mesh_cells),ierr)
mass_matrix_first_line(1:size(mass_matrix_period))=mass_matrix_period
mass_matrix_first_line(mesh_cells)=1.0_f64

call sll_display(mass_matrix_period(size(mass_matrix_period):2:-1), "(F10.5)" )
!call sll_display(mass_matrix_first_line(mesh_cells- (size(mass_matrix_period) -1) +1: mesh_cells:1 ), "(F10.5)" )
mass_matrix_first_line(mesh_cells- (size(mass_matrix_period) -1) +1: mesh_cells)=mass_matrix_period(size(mass_matrix_period):2:-1)

call sll_display(mass_matrix_first_line, "(F10.5)" )






!firstline= mass_matrix_period +

!print *, knots(1), knots(degree+2), "##"
    !SLL_DEALLOCATE_ARRAY(bspline_knot_values, ierr)
!print *, mass_matrix_period, "##"

    !bspline_knot_values=bspline_knot_values(degree+1:1:-1)



stop





    !print *, knots(1),";"
    !print *, eoshift(b_splines_at_x(arbitrarydegree1D, 1, knots(1)  ),1 ), ","

    SLL_ALLOCATE(bspline_knot_values(degree+1),ierr)
    bspline_knot_values=b_splines_at_x(arbitrarydegree1D, 1, knots(1))
    bspline_knot_values=bspline_knot_values(degree+1:1:-1)

    SLL_ALLOCATE(quadrature_point_vals( degree+1 ,degree+2),ierr )
    !Loop over first half of splines with intersecting support
    do i=0, degree+1
        print *, eoshift(bspline_knot_values, -i),";"
        quadrature_point_vals(:,i+1)= dot_product(bspline_knot_values , eoshift(bspline_knot_values, -i))

    enddo


    !eoshif(bspline_knot_values, i , bspline_knot_values)

    print *, bspline_knot_values
    quadrature_point_vals(:,1)= b_splines_at_x(arbitrarydegree1D, 1, knots(1) )


    do i=2,degree+2
        !Integrate with all splines
        !get quadratur points
        !quadrature_point_vals
        quadrature_point_vals(:,i)=b_splines_at_x(arbitrarydegree1D, i, knots(i) )
        !print *, knots(i),";"
        !print *, b_splines_at_x(arbitrarydegree1D, i, knots(i) ), ","
    enddo


    print *,quadrature_point_vals


    !print *, find_index( x, knots, num_pts )
    !cell =






    !arbitrarydegree1D=>new_arbitrary_degree_spline_1d( degree, knots, num_pts, SLL_PERIODIC)

    !sll_arbitrary_degree_splines



    SLL_ALLOCATE(particleposition(nparticles),ierr)
    SLL_ALLOCATE(particlespeed(nparticles),ierr)


    call load_particles (nparticles, interval_a, interval_b, particleposition, particlespeed )


    !Setup mesh for quasineutrality equation


    !poisson1dsolver=>new(interval_a,interval_b,100,ierr)
    !poisson1dsolver=>solve(poisson1dsolver, )




    !solve quasineutrality equation
    !sll_poisson_1d_periodic


    !Splines





    !call random_number(  particleposition )

    !CALL init_random_seed()
    !CALL RANDOM_NUMBER(particleposition)
    !CALL RANDOM_NUMBER(particlespeed)


    !    print *, "Hello Wodfrld!"

contains


    subroutine  load_particles (nparticles, interval_a, interval_b, particleposition, particlespeed)
#include "sll_working_precision.h"
#include "sll_memory.h"
        use sll_constants
        implicit none
        integer, intent(in) :: nparticles
        sll_real, intent(in) ::interval_a
        sll_real, intent(in) :: interval_b
        sll_real, DIMENSION(:), allocatable, intent(inout):: particleposition
        sll_real, DIMENSION(:), allocatable, intent(inout) :: particlespeed
        sll_real :: mu, sigma
        integer :: i

        !uniform random
        call random_number(particleposition)
        particleposition=interval_a + (interval_b -interval_a)*particleposition

        !Gaussian speed distribution
        mu=interval_a + (interval_b - interval_a )/2
        sigma=1.0
        do i=1,size(particlespeed)
            particlespeed(i)=gaussianrnd(mu, sigma  )
        end do
        !call random_number(particlespeed)
        !particlespeed=-1 + particlespeed*2
        return
    end subroutine load_particles


    sll_real function maxwellboltzmann (m ,T)
#include "sll_working_precision.h"
    use sll_constants
    IMPLICIT NONE

    sll_real, intent(in) :: T !< temperature in K
    sll_real, intent(in) :: m !< particle mass in kg

    !> @param Boltzmann constant (def) J/K
    sll_real64, parameter :: sll_kb = 1.3806488D-23
    sll_real ans

    call random_number(ans)

endfunction




end program unit_test
