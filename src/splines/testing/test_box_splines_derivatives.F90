!*****************************************************************************
!> @brief
!> sll_m_box_splines unit test
!> @author
!> Laura S. Mendoza
!*****************************************************************************

program test_box_splines_derivatives

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet

  use sll_m_box_splines, only: &
    sll_f_boxspline_x1_derivative, &
    sll_f_boxspline_x2_derivative, &
    sll_f_compute_box_spline, &
    sll_f_new_box_spline_2d, &
    sll_t_box_spline_2d, &
    sll_f_hex_interpolate_value

  use sll_m_hex_pre_filters, only: &
    sll_s_pre_filter_pfir

  use sll_m_hexagonal_meshes, only: &
    sll_f_new_hex_mesh_2d, &
    sll_s_write_field_hex_mesh_xmf, &
    sll_t_hex_mesh_2d

  use sll_m_constants, only : &
       sll_p_pi

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Variable to check is successful. By default it is the case.
  logical :: passed_test = .true.
  ! Variables for the test:
  type(sll_t_hex_mesh_2d),   pointer   :: mesh
  type(sll_t_box_spline_2d), pointer   :: spline
  sll_int32,  parameter :: max_deg = 5
  sll_real64, parameter :: criterion = 1.0e-14_f64
  sll_int32    :: ierr
  sll_int32    :: num_cells
  sll_int32    :: degree
  sll_int32    :: rule
  sll_int32    :: i
  sll_real64   :: x1, x1p
  sll_real64   :: x2, x2p
  sll_real64   :: somme
  sll_real64   :: val
  sll_real64   :: error_current
  sll_real64   :: error_linf
  sll_real64   :: error_l2
  sll_real64   :: appx_f
  sll_real64   :: randnum
  sll_real64, dimension(:), allocatable :: filters
  sll_real64, dimension(:), allocatable :: dist
  sll_real64, dimension(:), allocatable :: f
  sll_real64, dimension(:), allocatable :: dxf
  sll_real64, dimension(:), allocatable :: dyf
  sll_real64, dimension(:), allocatable :: f2
  sll_real64, dimension(:), allocatable :: dxf2
  sll_real64, dimension(:), allocatable :: dyf2
  !sll_real64, dimension(:), allocatable :: splines_on_support

  print *, " ************************************************** "
  print *, "       Testing arbitrary degree box-splines"
  print *, "            degree = 1,..., ", max_deg
  print *, " ************************************************** "

  print *, "  Test will be made on the following mesh:"
  ! Mesh initialization
  num_cells = 20
  mesh => sll_f_new_hex_mesh_2d(num_cells, 0._f64, 0._f64, radius = 2._f64)
  call mesh%display()

  ! Spline initialization
  spline => sll_f_new_box_spline_2d(mesh, sll_p_dirichlet)

  ! Allocations for boxsplines and derivatives :
  SLL_ALLOCATE(f(mesh%num_pts_tot),ierr)
  SLL_ALLOCATE(dxf(mesh%num_pts_tot),ierr)
  SLL_ALLOCATE(dyf(mesh%num_pts_tot),ierr)
  ! Allocation for degree 2 splines:
  SLL_ALLOCATE(f2(mesh%num_pts_tot),ierr)
  SLL_ALLOCATE(dxf2(mesh%num_pts_tot),ierr)
  SLL_ALLOCATE(dyf2(mesh%num_pts_tot),ierr)
  ! Allocation for a distribution function
  SLL_ALLOCATE(dist(mesh%num_pts_tot),ierr)

  ! Initialization of distribution function
  do i=1, mesh%num_pts_tot
     x1 = mesh%global_to_x1(i)
     x2 = mesh%global_to_x2(i)
     dist(i) = cos(x1*2.0*sll_p_pi)
  end do

  do degree=1,max_deg
     ! Computing box-splines coefficients:
     call spline%compute_coeff_box_spline_2d(dist, degree)
     ! Reseting errors:
     error_linf = 0._f64
     error_l2   = 0._f64
     do i=1, mesh%num_pts_tot
        x1 = mesh%global_to_x1(i)
        x2 = mesh%global_to_x2(i)
        if (i .lt. 3*(num_cells-4)*(num_cells-5)+1) then
           call random_number(randnum)
           x1p = x1 - randnum/num_cells
           x2p = x2 - randnum/num_cells
           ! Interpolation on mesh points computation:
           appx_f = sll_f_hex_interpolate_value(mesh, x1p, x2p, spline, degree)
           ! Computing errors:
           error_current = ABS(cos(x1p*2.0*sll_p_pi) - appx_f)
           error_l2   = error_l2 + error_current**2
           error_linf = max(error_linf, error_current)
        end if
        ! Computing values of box-splines and derivatives
        ! Note: This is not needed for interpolation
        f(i) = sll_f_compute_box_spline(spline, x1, x2, degree)
        dxf(i) = sll_f_boxspline_x1_derivative(x1, x2, degree)
        dyf(i) = sll_f_boxspline_x2_derivative(x1, x2, degree)
     end do
     ! Storing degree two results:
     if (degree .eq. 2) then
        f2(:) = f(:)
        dxf2(:) = dxf(:)
        dyf2(:) = dyf(:)
     end if
     ! Printing errors of interpolation and values of box-splines :
     print *, " * Degree =", degree
     print *, "    --> Error in L_inf     = ", error_linf
     print *, "    --> Error in L_2       = ", SQRT(error_l2)
     ! print *, "    --> sum boxspline      = ", sum(f)
     ! print *, "    --> sum 1st derivative = ", sum(dxf)
     ! print *, "    --> sum 2nd derivative = ", sum(dyf)
     if (abs(sum(f)-1.0_f64).gt.criterion) then
        passed_test = .false.
        print *, " Test FAILED for box-spline of degree = ", degree
        STOP
     end if
  end do


  !Wrtting on docs:
  call sll_s_write_field_hex_mesh_xmf(mesh,   f2, "boxspline2")
  call sll_s_write_field_hex_mesh_xmf(mesh, dxf2, "der1_boxspline2")
  call sll_s_write_field_hex_mesh_xmf(mesh, dyf2, "der2_boxspline2")

  SLL_DEALLOCATE_ARRAY(dist, ierr)

  SLL_DEALLOCATE_ARRAY(f, ierr)
  SLL_DEALLOCATE_ARRAY(dxf, ierr)
  SLL_DEALLOCATE_ARRAY(dyf, ierr)

  SLL_DEALLOCATE_ARRAY(f2, ierr)
  SLL_DEALLOCATE_ARRAY(dxf2, ierr)
  SLL_DEALLOCATE_ARRAY(dyf2, ierr)


  ! ! Computing non null splines on one cell:
  ! SLL_ALLOCATE(splines_on_support(degree*degree*3),ierr)
  ! splines_on_support = non_zeros_splines(mesh, 2, degree)

  ! print *, "Non null splines of degree 1 at cell#2 (expected 1,2,3) =", &
  !      splines_on_support(1), &
  !      splines_on_support(2), &
  !      splines_on_support(3)

  ! SLL_DEALLOCATE_ARRAY(splines_on_support, ierr)
  ! call sll_o_delete(spline) !also deletes the mesh


  ! Testing degree 3 boxsplines:
  ! print *, ""
  ! print *, "------------ testing degree 3 -----------"
  ! num_cells = 3
  ! mesh => sll_f_new_hex_mesh_2d(num_cells, 0._f64, 0._f64, radius = 1._f64)
  ! call mesh%display()
  ! degree = 3
  ! spline => sll_f_new_box_spline_2d(mesh, sll_p_dirichlet)
  ! SLL_ALLOCATE(f(mesh%num_pts_tot),ierr)
  ! do i=1, mesh%num_pts_tot
  !    x1 = mesh%global_to_x1(i)
  !    x2 = mesh%global_to_x2(i)
  !    f(i) = sll_f_compute_box_spline(spline, x1, x2, degree)
  ! end do
  ! print *, sum(f)
  ! call sll_s_write_field_hex_mesh_xmf(mesh, f, "chi3")


  ! Testing pre-filter:
  degree = 1
  print *, ""
  print *, " ********** Testing the prefilters ********"
  print *, " ***************** PFIR *******************"
  print *, " Degree = ", degree
  somme = 0._f64
  call sll_s_pre_filter_pfir(mesh, degree, filters)
  do i=1, 3*(degree+1)*(degree)+1
     val = filters(i)
     somme = somme + val
     print *, "     ", i, " ---->", val
  end do
  print *, "     sum =", somme

end program test_box_splines_derivatives
