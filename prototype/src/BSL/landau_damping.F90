program landau_damping
#include "sll_assert.h"
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"
#include "sll_poisson_solvers.h"

use numeric_constants
use distribution_function
use sll_splines
use sll_diagnostics
use sll_poisson_2d_periodic

implicit none
  

!Geometry
sll_real64 :: eta1_min, eta1_max, eta2_min, eta2_max
sll_real64 :: eta3_min, eta3_max, eta4_min, eta4_max
type(geometry_2D), pointer :: geom_x, geom_v

!Grids
sll_int32  :: nc_eta1, nc_eta2, nc_eta3, nc_eta4
sll_real64 :: delta_eta1, delta_eta2, delta_eta3, delta_eta4
type(mesh_descriptor_2D), pointer :: mesh_x, mesh_v

!Time domain
sll_int32  :: i_step, n_step
sll_real64 :: delta_t
sll_real32 :: time

!Distribution function 4D
type(sll_distribution_function_4D_t), pointer :: dist_func

!Poisson solver
type(poisson_2d_periodic), pointer :: poisson
type(field_2D_vec2),       pointer :: exy
type(field_2D_vec1),       pointer :: rho

!Semi lagrangian scheme
type (sll_spline_1D), pointer :: spl_eta1
type (sll_spline_1D), pointer :: spl_eta2
type (sll_spline_1D), pointer :: spl_eta3
type (sll_spline_1D), pointer :: spl_eta4

sll_real64, dimension(:), pointer  ::  eta1_out , eta1
sll_real64, dimension(:), pointer  ::  eta2_out , eta2
sll_real64, dimension(:), pointer  ::  eta3_out , eta3
sll_real64, dimension(:), pointer  ::  eta4_out , eta4

!Diagnostics and errors
sll_int32                          :: error
sll_real32 :: nrj
character(len=4) :: counter

!Local indices
sll_int32  :: i1, i2, i3, i4

!x domain
eta1_min =  0.0_f64; eta1_max =  2.0_f64 * sll_pi
eta2_min =  0.0_f64; eta2_max =  2.0_f64 * sll_pi

geom_x => new_geometry_2D('cartesian')

nc_eta1 = 128; nc_eta2 = 32

mesh_x => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
          PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom_x)

call write_mesh_2D(mesh_x,"mesh_x")

!v domain
eta3_min = -6.0_f64; eta3_max =  6.0_f64 
eta4_min = -6.0_f64; eta4_max =  6.0_f64 

geom_v => new_geometry_2D('cartesian')

nc_eta3 = 64; nc_eta4 = 64

mesh_v => new_mesh_descriptor_2D(eta3_min, eta3_max, nc_eta3, &
          PERIODIC, eta4_min, eta4_max, nc_eta4, PERIODIC, geom_v)


rho     => new_field_2D_vec1(mesh_x)
exy     => new_field_2D_vec2(mesh_x)

poisson   => new_poisson_2d_periodic(exy)

dist_func => sll_new_distribution_function_4D(mesh_x, mesh_v, NODE_CENTERED_DF, 'f')

call sll_init_distribution_function_4D( dist_func, LANDAU)

delta_eta1 = GET_MESH_DELTA_ETA1(mesh_x)
delta_eta2 = GET_MESH_DELTA_ETA2(mesh_x)
delta_eta3 = GET_MESH_DELTA_ETA1(mesh_v)
delta_eta4 = GET_MESH_DELTA_ETA2(mesh_v)

!call write_distribution_function(dist_func)

! initialize splines
spl_eta1 => new_spline_1D( nc_eta1+1,eta1_min,eta1_max,PERIODIC_SPLINE )
spl_eta2 => new_spline_1D( nc_eta2+1,eta2_min,eta2_max,PERIODIC_SPLINE )  
spl_eta3 => new_spline_1D( nc_eta3+1,eta3_min,eta3_max,HERMITE_SPLINE )
spl_eta4 => new_spline_1D( nc_eta4+1,eta4_min,eta4_max,HERMITE_SPLINE )  

SLL_ALLOCATE(eta1(nc_eta1+1),error)
SLL_ALLOCATE(eta2(nc_eta2+1),error)
SLL_ALLOCATE(eta3(nc_eta3+1),error)
SLL_ALLOCATE(eta4(nc_eta4+1),error)

SLL_ALLOCATE(eta1_out(nc_eta1+1),error)
SLL_ALLOCATE(eta2_out(nc_eta2+1),error)
SLL_ALLOCATE(eta3_out(nc_eta3+1),error)
SLL_ALLOCATE(eta4_out(nc_eta4+1),error)

do i1 = 1, nc_eta1+1
   eta1(i1) = eta1_min + (i1-1) * delta_eta1
end do
do i2 = 1, nc_eta2+1
   eta2(i2) = eta2_min + (i2-1) * delta_eta2
end do
do i3 = 1, nc_eta3+1
   eta3(i3) = eta3_min + (i3-1) * delta_eta3
end do
do i4 = 1, nc_eta4+1
   eta4(i4) = eta4_min + (i4-1) * delta_eta4
end do

n_step = 1000
delta_t = .01_f64
time = 0.0_f32

do i_step = 1, n_step !Loop over time

   if (i_step == 1 ) then
      call advection_x1(0.5_f64*delta_t)
   else
      call advection_x1(1.0_f64*delta_t)
   end if

   call advection_x2(0.5_f64*delta_t)

   call compute_rho(dist_func, rho)

   call solve_poisson_2d_periodic(poisson,exy,rho,error)

   call advection_v1(0.5_f64*delta_t)

   call advection_v2(1.0_f64*delta_t)

   call advection_v1(0.5_f64*delta_t)

   call advection_x2(0.5_f64*delta_t)

   time  = time + delta_t

   nrj = sum(exy%data%v1*exy%data%v1+exy%data%v2*exy%data%v2) 
   nrj = nrj*delta_eta1*delta_eta2
   nrj = 0.5_f64*log(nrj)
   
   write(11,*) time, nrj
   write(*,*) time, nrj

   call int2string(dist_func%plot_counter,counter)
   !call write_distribution_function(dist_func)
   !call write_vec1d(rho%data,mesh_x%nc_eta1+1,mesh_x%nc_eta2+1, &
   !                 "rho"//counter,"mesh_x",0)
   !call write_vec2d(exy%data%v1,exy%data%v2,mesh_x%nc_eta1+1, &
   !                 mesh_x%nc_eta2+1,"exy"//counter,"mesh_x",0)

end do !next time step

call delete_poisson_2d_periodic(poisson)
call delete_field_2D_vec1( rho )
call delete_field_2D_vec2( exy )

contains

subroutine advection_x1(dt)
sll_real64, intent(in) :: dt

   do i4 = 1, nc_eta4+1
   do i3 = 1, nc_eta3+1
   do i2 = 1, nc_eta2+1

      call sl_step( dist_func%field%data(:,i2,i3,i4),    &
                    eta1, nc_eta1, eta3(i3), eta1_min,   &
                    eta1_max, spl_eta1, dt )

   end do
   end do
   end do

end subroutine advection_x1

subroutine advection_x2(dt)
sll_real64, intent(in) :: dt

   do i4 = 1, nc_eta4+1
   do i3 = 1, nc_eta3+1
   do i1 = 1, nc_eta1+1

      call sl_step( dist_func%field%data(i1,:,i3,i4),    &
                    eta2, nc_eta2, eta4(i4), eta2_min,   &  
                    eta2_max, spl_eta2, dt )

   end do
   end do
   end do


end subroutine advection_x2

subroutine advection_v1(dt)
sll_real64, intent(in) :: dt

   do i4 = 1, nc_eta4+1
   do i2 = 1, nc_eta2+1
   do i1 = 1, nc_eta1+1

      call sl_step( dist_func%field%data(i1,i2,:,i4),     &
                    eta3, nc_eta3, exy%data(i1,i2)%v1,&
                    eta3_min, eta3_max, spl_eta3, dt )

   end do
   end do
   end do

end subroutine advection_v1

subroutine advection_v2(dt)
sll_real64, intent(in) :: dt

   do i3 = 1, nc_eta3+1
   do i2 = 1, nc_eta2+1
   do i1 = 1, nc_eta1+1

      call sl_step( dist_func%field%data(i1,i2,i3,:),      &
                    eta4, nc_eta4, exy%data(i1,i2)%v2, &
                    eta4_min, eta4_max, spl_eta4, dt )

   end do
   end do
   end do

end subroutine advection_v2

subroutine sl_step( array_1d, eta_in, nc_eta, &
                    flux, eta_min, eta_max, spl_eta, dt )

sll_real64, dimension(:)        :: array_1d
sll_real64, dimension(nc_eta+1) :: eta_in
sll_real64, dimension(nc_eta+1) :: eta_out
sll_real64                      :: flux
sll_int32                       :: i
sll_int32                       :: nc_eta
sll_real64                      :: eta_min
sll_real64                      :: eta_max
sll_real64                      :: dt
type (sll_spline_1D), pointer   :: spl_eta

!Periodic boundary conditions (both spaces)
do i = 1, nc_eta+1
   eta_out(i) = eta_in(i) - flux*dt
   if (eta_out(i) < eta_min) then
      eta_out(i) = eta_out(i) + eta_max - eta_min
   else if (eta_out(i) > eta_max) then
      eta_out(i) = eta_out(i) - eta_max + eta_min
   end if
end do

call compute_spline_1D_periodic( array_1d, spl_eta )
call interpolate_array_values( eta_out, array_1d, nc_eta+1, spl_eta )

end subroutine sl_step
  
end program landau_damping
