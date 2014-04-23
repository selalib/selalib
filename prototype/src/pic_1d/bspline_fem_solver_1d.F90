module sll_bspline_fem_solver_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"
    use sll_arbitrary_degree_spline_interpolator_1d_module
    use sll_arbitrary_degree_splines
    use gauss_legendre_integration
    use sll_fft
    implicit none


contains



function sll_bspline_fem_solve_poisson_1d(bspline_arbitrary_degree, knots,b)  result(x)
 sll_real64, dimension(:), intent(in)   :: b
 sll_real64, dimension(:) :: stiffn_matrix_period
 sll_real64 ::stiffn_matrix_first_line(size(knots)-1)


stiffn_matrix_first_line=0.0_f64
stiffn_matrix_period=gen_stiffn_matrix_period (bspline_arbitrary_degree, knots)
stiffn_matrix_first_line(mesh_cells- (size(stiffn_matrix_period) -1) +1: size(knots)-1)=stiffn_matrix_period(size(stiffn_matrix_period):2:-1)





stiffn_matrix_first_line(1:size(stiffn_matrix_period))=stiffn_matrix_period

call sll_display(mass_matrix_period(size(mass_matrix_period):2:-1), "(F10.5)" )



endfunction



!<Solve Cx=b
function bspline_fem_solver_1d_solve_circulant_matrix_equation( c ,b  )  result(x)

    sll_real64, dimension(:), intent(in)   :: c
    sll_real64, dimension(:), intent(in)   :: b
    integer :: N
    sll_real64, dimension(:),allocatable   :: fft_computational_data
    sll_real64, dimension(:), allocatable   :: x
    integer :: ierr
    N =size(c)
SLL_ASSERT(size(c)==size(x))
SLL_ASSERT(size(c)==size(b))


!Initialize discrete fourier transformation with periodicity N
SLL_ALLOCATE( fft_computational_data(2*N+15), ierr)
call dffti(N,fft_computational_data)

call dfftf(N, b, fft_computational_data)
call dfftf(N, c, fft_computational_data)

!Backward transformation
SLL_ALLOCATE( x(N), ierr)
x=b/c
call dfftb(N, x, fft_computational_data)


SLL_DEALLOCATE_ARRAY(fft_computational_data,ierr)
endfunction




  function gen_mass_matrix_period (bspline_arbitrary_degree, knots ) &
        result( mass_matrix_period)
    type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline_arbitrary_degree
    sll_real64, dimension(:), intent(in)     :: knots
    sll_real64, dimension(:), allocatable    :: mass_matrix_period

       mass_matrix_period=gen_fem_matrix_period (bspline_arbitrary_degree, knots, b_splines_at_x )
  endfunction

  function gen_stiffn_matrix_period (bspline_arbitrary_degree, knots) &
        result( stiffn_matrix_period)
    type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline_arbitrary_degree
    sll_real64, dimension(:), intent(in)     :: knots
    sll_real64, dimension(:), allocatable    :: stiffn_matrix_period

      stiffn_matrix_period=gen_fem_matrix_period (bspline_arbitrary_degree, knots, b_spline_derivatives_at_x )
  endfunction



    function gen_fem_matrix_period (bspline_arbitrary_degree, knots, spline_evaluation ) &
        result( fem_matrix_period)
    type(arbitrary_degree_spline_1d), pointer, intent(in) :: bspline_arbitrary_degree
    sll_real64, dimension(:), intent(in)     :: knots
    sll_real64, dimension(:), allocatable    :: fem_matrix_period
    procedure (b_splines_at_x) :: spline_evaluation

    sll_int32                   :: degree


    sll_real64, dimension(:,:), allocatable   :: quadrature_point_vals
    sll_real64, dimension(:,:), allocatable   :: quadrature_points_weights
    sll_real64, dimension(:,:), allocatable   :: bspline_qpoint_values

    integer                                   :: quadrature_npoints
    integer :: cell
    integer :: i,j
    integer :: ierr
    degree=bspline_arbitrary_degree%degree

    !Get Gauss Legendre points and weights to be exact for the selected spline degree
    !Note a Bspline is a piecewise polynom
    quadrature_npoints=ceiling(0.5_f64*real(degree+1))

    SLL_ALLOCATE(quadrature_points_weights(2,1:quadrature_npoints),ierr)

    !print *,quadrature_points_weights
    !Get spline values at
    SLL_ALLOCATE(bspline_qpoint_values(1:degree+1,1:quadrature_npoints), ierr)

    SLL_ASSERT(size(knots) >degree+1)


    SLL_CLEAR_ALLOCATE(fem_matrix_period(degree+1),ierr)

    !Quadrature for each support cell
    do cell=1,degree+1
        quadrature_points_weights=gauss_legendre_points_and_weights(quadrature_npoints, knots(cell), knots(cell+1))
        do i=1,quadrature_npoints
            bspline_qpoint_values(:,i)=spline_evaluation(bspline_arbitrary_degree, cell, quadrature_points_weights(1,i ))
            !reorientate to the right side
            bspline_qpoint_values(:,i) = bspline_qpoint_values((degree+1):1:-1,i)
            !shift in order to account for intersecting and limited support
            bspline_qpoint_values(:,i) = eoshift(bspline_qpoint_values(:,i), (cell-1))
        enddo
        do j=1, degree+1
            fem_matrix_period(j)= fem_matrix_period(j) + &
                dot_product(bspline_qpoint_values(1, :)*bspline_qpoint_values(j,:), quadrature_points_weights(2,:)  )
        enddo
    enddo



    call sll_display(knots(1:degree+2), "(F10.5)" )
    call sll_display(quadrature_points_weights(1,:), "(F10.5)" )
    call sll_display(quadrature_points_weights(2,:), "(F10.5)" )

    call sll_display(bspline_qpoint_values, "(F10.5)" )
    write(*,*)
    call sll_display(fem_matrix_period, "(F10.5)" )

!!<SLL_DEALLOCATE( bspline_qpoint_values,ierr)
!!<SLL_DEALLOCATE( quadrature_points_weights,ierr)
endfunction



end module
