#define MPI_MASTER 0

! Sample computation with the following characteristics:
! - Simple advection 
! - Polar mesh with analytic coordinate transformation
! - Initial field is a gaussain function
!

program polar_advection
#include "selalib.h"

implicit none

type(cubic_spline_2d_interpolator), target  :: spline_xy

type(sll_logical_mesh_2d), pointer                    :: logical_mesh
class(sll_coordinate_transformation_2d_base), pointer :: transfx

sll_int32  :: itime
sll_int32  :: nbiter = 1000
sll_int32  :: i
sll_int32  :: j

sll_int32,  parameter :: nc_eta1  = 64 
sll_int32,  parameter :: nc_eta2  = 64
sll_real64, parameter :: eta1_min = 2.0_f64 
sll_real64, parameter :: eta1_max = 8.0_f64 
sll_real64, parameter :: eta2_min = -3.141592654_f64
sll_real64, parameter :: eta2_max =  3.141592654_f64
sll_real64, parameter :: deltat   = 0.01_f64

sll_int32  :: error
sll_real64 :: eta1, eta2, eta3, eta4
sll_real64 :: alpha1, alpha2
sll_real64 :: inv_j(2,2)
sll_real64, dimension(:,:), allocatable :: f
sll_real64, dimension(:,:), allocatable :: x1
sll_real64, dimension(:,:), allocatable :: x2
sll_int32 :: nc_x1, nc_x2

logical_mesh => new_logical_mesh_2d(nc_eta1,  nc_eta2,  &
                                    eta1_min, eta1_max, & 
                                    eta2_min, eta2_max)

call spline_xy%initialize( nc_eta1+1, nc_eta2+1, &
                           eta1_min, eta1_max,   &
                           eta2_min, eta2_max,   &
                           SLL_PERIODIC, SLL_PERIODIC )

transfx => new_coordinate_transformation_2d_analytic( &
       "analytic_polar_transformation", &
       logical_mesh, &
       polar_x1, &
       polar_x2, &
       polar_jac11, &
       polar_jac12, &
       polar_jac21, &
       polar_jac22, (/0.0_f64/) )

nc_x1 = nc_eta1
nc_x2 = nc_eta2
SLL_CLEAR_ALLOCATE(x1(1:nc_x1+1,1:nc_x2+1),error)
SLL_CLEAR_ALLOCATE(x2(1:nc_x1+1,1:nc_x2+1),error)
SLL_CLEAR_ALLOCATE(f(1:nc_eta1+1,1:nc_eta2+1),error)

do j=1,nc_eta2+1
   do i=1,nc_eta1+1
      x1(i,j) = transfx%x1_at_node(i,j)
      x2(i,j) = transfx%x2_at_node(i,j)
      f(i,j)  = exp(-0.5_f64*((x1(i,j)-5.)**2+(x2(i,j))**2))
   end do
end do


do itime = 1, nbiter

  call spline_xy%compute_interpolants(f)

   do j=1,nc_eta2+1
      do i=1,nc_eta1+1

         eta1 = eta1_min + (i-1)*logical_mesh%delta_eta1
         eta2 = eta2_min + (j-1)*logical_mesh%delta_eta2
         eta3 =  eta1*sin(eta2)
         eta4 = -eta1*cos(eta2)

         inv_j  = transfx%inverse_jacobian_matrix(eta1,eta2)
         alpha1 = -deltat*(inv_j(1,1)*eta3 + inv_j(1,2)*eta4)
         alpha2 = -deltat*(inv_j(2,1)*eta3 + inv_j(2,2)*eta4)

         eta1 = eta1_min + modulo(eta1 - eta1_min - alpha1, eta1_max-eta1_min)
         eta2 = eta2_min + modulo(eta2 - eta2_min - alpha2, eta2_max-eta2_min)

         f(i,j) = spline_xy%interpolate_value(eta1,eta2)

      end do
   end do

   call sll_gnuplot_2d(nc_x1+1,nc_x2+1,x1,x2,f,"f_polar_advection",itime,error)

end do


end program polar_advection


