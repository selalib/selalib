program landau_2d
#include "sll_assert.h"
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"
#include "sll_poisson_solvers.h"

use numeric_constants
use distribution_function
use sll_splines
use sll_diagnostics
use sll_poisson_1d_periodic

implicit none
  
sll_int32  :: i_step, n_step, i1, i2, j_step
sll_int32  :: nc_eta1, nc_eta2

sll_real64 :: eta1_min, eta1_max
sll_real64 :: eta2_min, eta2_max
sll_real64 :: eta1, eta2
sll_real64 :: delta_eta1, delta_eta2
sll_real64 :: delta_t
sll_real64 :: time

sll_real64, dimension(:), allocatable :: nrj
sll_real64, dimension(:), allocatable :: eta1_out
sll_real64, dimension(:), allocatable :: eta2_out

type(mesh_descriptor_1D),             pointer :: mesh_x
type(mesh_descriptor_1D),             pointer :: mesh_v
type(mesh_descriptor_2D),             pointer :: mesh_xv
type(sll_distribution_function_2D_t), pointer :: df_2d

type (sll_spline_1D), pointer      :: spl_eta1
type (sll_spline_1D), pointer      :: spl_eta2

type (field_1D_vec1), pointer      :: ex
type (field_1D_vec1), pointer      :: rho
type (poisson_1d_periodic)         :: poisson

sll_real64 :: eps, kx
sll_int32  :: error

eta1_min = 0.0_f64
eta1_max = 4.0_f64*sll_pi
nc_eta1  = 256

mesh_x => new_mesh_descriptor_1D(eta1_min, eta1_max, nc_eta1, PERIODIC)

eta2_min = -6.0_f64
eta2_max = +6.0_f64
nc_eta2  = 256

mesh_v => new_mesh_descriptor_1D(eta2_min, eta2_max, nc_eta2, PERIODIC)

rho    => new_field_1D_vec1( mesh_x )
ex     => new_field_1D_vec1( mesh_x )

mesh_xv => mesh_x * mesh_v

df_2d => sll_new_distribution_function_2D(mesh_xv, NODE_CENTERED_DF, 'df')

delta_eta1 = mesh_xv%delta_eta1
delta_eta2 = mesh_xv%delta_eta2

call new(poisson, nc_eta1, error) 

call write_mesh_2D(mesh_xv)

!Initialize distribution function
eps = 0.05_f64
kx  = 0.5_f64 
eta2 = eta2_min
do i2=1, nc_eta2+1
   eta1 = eta1_min
   do i1=1, nc_eta1+1
      df_2d%field%data(i1,i2)=(1.0_f64+eps*cos(kx*eta1))/(2.0_f64*sll_pi) &
                              * exp(-0.5_f64*eta2*eta2)
      eta1 = eta1 + delta_eta1
   end do
   eta2 = eta2 + delta_eta2
end do

call write_vec1d(df_2d%field%data,nc_eta1+1,nc_eta2+1,"df","mesh",0)

! run BSL method using 10 time steps and second order splitting
time = 0.0_f64
n_step = 1000
delta_t = 0.05_f64
SLL_ALLOCATE(nrj(n_step), error)
print*, "set title'", delta_eta1, delta_eta2,"'"
print*, 'set term x11'

spl_eta1 => new_spline_1D(nc_eta1+1,eta1_min,eta1_max,PERIODIC_SPLINE)
spl_eta2 => new_spline_1D(nc_eta2+1,eta2_min,eta2_max,HERMITE_SPLINE)  

! allocation
SLL_ALLOCATE(eta1_out(nc_eta1+1),error)
SLL_ALLOCATE(eta2_out(nc_eta2+1),error)

do i_step = 1, n_step

   eta2 = eta2_min 
   do i2 = 1, nc_eta2+1
      call sl_step( df_2d%field%data(:,i2), eta2, spl_eta1, 0.5_f64*delta_t )
      eta2 = eta2 + delta_eta2
   end do

   do i1 = 1, nc_eta1+1
      rho%data(i1) = sum(df_2d%field%data(i1,:))*delta_eta2
   end do
   
   call solve(poisson, ex , rho)

   do i1 = 1, nc_eta1+1
      call sl_step( df_2d%field%data(i1,:), ex%data(i1), spl_eta2, delta_t )
   end do

   eta2 = eta2_min 
   do i2 = 1, nc_eta2+1
      call sl_step( df_2d%field%data(:,i2), eta2, spl_eta1, 0.5_f64*delta_t )
      eta2 = eta2 + delta_eta2
   end do

   nrj(i_step) = 0.5_f64*log(sum(ex%data*ex%data)*delta_eta1)
   
   time = time + delta_t

   write(*,100) .0,n_step*delta_t,-19.5,0.5
   do j_step = 1, i_step
      print*, (j_step-1)*delta_t, nrj(j_step)
   end do
   print*, 'e'

end do

100 format('p [',f5.1,':',f5.1,'][',f6.1,':',f6.1,'] ''-'' w l')

contains

subroutine sl_step( array_1d, flux, spl_eta, dt )

   implicit none
   sll_real64, dimension(:)          :: array_1d
   sll_real64                        :: flux
   sll_int32                         :: i
   sll_real64                        :: dt
   type (sll_spline_1D), pointer     :: spl_eta
   sll_real64, dimension(:), pointer :: eta_out

   allocate(eta_out(spl_eta%n_points)) 

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
  
end program landau_2d
