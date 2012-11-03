program landau_4d
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
use sll_bsl

implicit none
  

!Geometry
sll_real64 :: eta1_min, eta1_max, eta2_min, eta2_max
sll_real64 :: eta3_min, eta3_max, eta4_min, eta4_max

!Grids
sll_int32  :: nc_eta1, nc_eta2, nc_eta3, nc_eta4
sll_real64 :: delta_eta1, delta_eta2, delta_eta3, delta_eta4

!Time domain
sll_int32  :: i_step, n_step, j_step
sll_real64 :: delta_t
sll_real64 :: time

!Distribution function 4D
sll_real64, dimension(:,:,:,:), allocatable :: dist_func
sll_real64, dimension(:,:), allocatable :: ex
sll_real64, dimension(:,:), allocatable :: ey
sll_real64, dimension(:,:), allocatable :: rho

!Poisson solver
type(poisson_2d_periodic), pointer :: poisson

sll_real64, dimension(:), pointer  ::  eta1_out 
sll_real64, dimension(:), pointer  ::  eta2_out
sll_real64, dimension(:), pointer  ::  eta3_out
sll_real64, dimension(:), pointer  ::  eta4_out

!Diagnostics and errors
sll_int32                             :: error
sll_real64, dimension(:), allocatable :: nrj
character(len=4) :: cplot

!Local indices
sll_int32  :: i1, i2, i3, i4

!x domain
eta1_min =  0.0_f64; eta1_max =  4.0_f64 * sll_pi
eta2_min =  0.0_f64; eta2_max =  4.0_f64 * sll_pi

nc_eta1 = 31; nc_eta2 = 31

delta_eta1 = (eta1_max-eta1_min)/nc_eta1
delta_eta2 = (eta2_max-eta2_min)/nc_eta2

!v domain
eta3_min = -6.0_f64; eta3_max = 6.0_f64 
eta4_min = -6.0_f64; eta4_max = 6.0_f64 

nc_eta3 = 31; nc_eta4 = 31

delta_eta3 = (eta3_max-eta3_min)/nc_eta3
delta_eta4 = (eta4_max-eta4_min)/nc_eta4

SLL_ALLOCATE(rho(nc_eta1+1,nc_eta2+1)
SLL_ALLOCATE(ex(nc_eta1+1,nc_eta2+1)
SLL_ALLOCATE(ey(nc_eta1+1,nc_eta2+1)

poisson   => new_poisson_2d_periodic(exy)

dist_func => sll_new_distribution_function_4D(mesh_x, mesh_v, NODE_CENTERED_DF, 'f')
bsl       => new_bsl_workspace(dist_func)

call sll_init_distribution_function_4D( dist_func, LANDAU)

call compute_rho(dist_func, rho)

call solve_poisson_2d_periodic(poisson,exy,rho,error)

!call write_distribution_function(dist_func)

SLL_ALLOCATE(eta1_out(nc_eta1+1),error)
SLL_ALLOCATE(eta2_out(nc_eta2+1),error)
SLL_ALLOCATE(eta3_out(nc_eta3+1),error)
SLL_ALLOCATE(eta4_out(nc_eta4+1),error)

n_step = 1000
SLL_ALLOCATE(nrj(n_step), error)
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

   nrj(i_step) = sum(exy%data%v1*exy%data%v1+exy%data%v2*exy%data%v2) 
   nrj(i_step) = nrj(i_step)*delta_eta1*delta_eta2
   nrj(i_step) = 0.5_f64*log(nrj(i_step))
   
   write(11,*) time, nrj(i_step)

   call int2string(dist_func%plot_counter,counter)
   !call write_distribution_function(dist_func)
   !call write_vec1d(rho%data,mesh_x%nc_eta1+1,mesh_x%nc_eta2+1, &
   !                 "rho"//counter,"mesh_x",0)
   !call write_vec2d(exy%data%v1,exy%data%v2,mesh_x%nc_eta1+1, &
   !                 mesh_x%nc_eta2+1,"exy"//counter,"mesh_x",0)

   write(*,100) .0,10.,-15.,1.
   do j_step = 1, i_step
      print*, (j_step-1)*delta_t, nrj(j_step)
   end do
   print*, 'e'


end do !next time step

call delete_poisson_2d_periodic(poisson)
call delete_field_2D_vec1( rho )
call delete_field_2D_vec2( exy )

100 format('p [',f5.1,':',f5.1,'][',f6.1,':',f6.1,'] ''-'' w l')

contains

subroutine advection_x1(dt)

   sll_real64, intent(in) :: dt
   sll_real64             :: eta3

   do i4 = 1, nc_eta4+1
      eta3 = eta3_min
      do i3 = 1, nc_eta3+1
         do i2 = 1, nc_eta2+1
            call sl_step( dist_func%field%data(:,i2,i3,i4), eta3, bsl%spl_eta1, dt )
         end do
         eta3 = eta3 + delta_eta3
      end do
   end do

end subroutine advection_x1

subroutine advection_x2(dt)

   sll_real64, intent(in) :: dt
   sll_real64             :: eta4

   eta4 = eta4_min
   do i4 = 1, nc_eta4+1
      do i3 = 1, nc_eta3+1
         do i1 = 1, nc_eta1+1
            call sl_step( dist_func%field%data(i1,:,i3,i4), eta4, bsl%spl_eta2, dt )
         end do
      end do
      eta4 = eta4 + delta_eta4
   end do


end subroutine advection_x2

subroutine advection_v1(dt)
   sll_real64, intent(in) :: dt

   do i4 = 1, nc_eta4+1
   do i2 = 1, nc_eta2+1
   do i1 = 1, nc_eta1+1

      call sl_step( dist_func%field%data(i1,i2,:,i4), exy%data(i1,i2)%v1,&
                    bsl%spl_eta3, dt )

   end do
   end do
   end do

end subroutine advection_v1

subroutine advection_v2(dt)

   sll_real64, intent(in) :: dt

   do i3 = 1, nc_eta3+1
   do i2 = 1, nc_eta2+1
   do i1 = 1, nc_eta1+1

      call sl_step( dist_func%field%data(i1,i2,i3,:), exy%data(i1,i2)%v2, &
                    bsl%spl_eta4, dt )

   end do
   end do
   end do

end subroutine advection_v2

subroutine sl_step( array_1d, flux, spl_eta, dt )

   implicit none
   sll_real64, dimension(:)          :: array_1d
   sll_real64                        :: flux
   sll_int32                         :: i
   sll_real64                        :: dt
   type (sll_spline_1D), pointer     :: spl_eta
   sll_real64, dimension(:), pointer :: eta_out

   allocate(eta_out(spl_eta%n_points)) 
   !Periodic boundary conditions (both spaces)
   do i = 1, spl_eta%n_points
      eta_out(i) = spl_eta%xmin +(i-1)*spl_eta%delta - flux*dt
      if (eta_out(i) < spl_eta%xmin) then
         eta_out(i) = eta_out(i) + spl_eta%xmax - spl_eta%xmin
      else if (eta_out(i) > spl_eta%xmax) then
         eta_out(i) = eta_out(i) - spl_eta%xmax + spl_eta%xmin
      end if
   end do
   
   call compute_spline_1D_periodic( array_1d, spl_eta )
   call interpolate_array_values( eta_out, array_1d, size(array_1d), spl_eta )

end subroutine sl_step
  
end program landau_4d
