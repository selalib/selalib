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
use sll_bsl

implicit none
  
sll_real64 :: eta1_min, eta1_max, eta2_min, eta2_max
sll_int32  :: nc_eta1, nc_eta2
sll_real64 :: eta1, eta2, delta_eta1, delta_eta2
sll_int32  :: i_step, n_step, i1, i2, j_step
sll_real64 :: delta_t
sll_real64 :: time
sll_real64, dimension(:), allocatable :: nrj

type(mesh_descriptor_1D),             pointer :: mesh_x
type(mesh_descriptor_1D),             pointer :: mesh_v
type(mesh_descriptor_2D),             pointer :: mesh_xv
type(sll_distribution_function_2D_t), pointer :: df_2d
type(field_1D_vec1),                  pointer :: adv_field_x
type(bsl_workspace_1d),               pointer :: bsl_work_x
type(bsl_workspace_1d),               pointer :: bsl_work_v

type (field_1D_vec1), pointer      :: ex
type (field_1D_vec1), pointer      :: rho
type (poisson_1d_periodic)         :: poisson
sll_int32   :: error
sll_real64 :: eps, kx

eta1_min =  0.0_f64; eta1_max = 4*sll_pi;
eta2_min = -6.0_f64; eta2_max = 6.0_f64
nc_eta1 = 128
nc_eta2 = 128

mesh_x => new_mesh_descriptor_1D(eta1_min, eta1_max, nc_eta1, PERIODIC)
mesh_v => new_mesh_descriptor_1D(eta2_min, eta2_max, nc_eta2, PERIODIC)

rho    => new_field_1D_vec1( mesh_x )
ex     => new_field_1D_vec1( mesh_x )

call new(poisson, nc_eta1, error) 

mesh_xv => mesh_x * mesh_v

df_2d => sll_new_distribution_function_2D(mesh_xv, NODE_CENTERED_DF, 'one_d')

delta_eta1 = mesh_x%delta_eta1
delta_eta2 = mesh_v%delta_eta1

call write_mesh_2D(mesh_xv)

eps = 0.05_f64
kx  = 2.0_f64*sll_pi/(nc_eta1*delta_eta1)
eta2 = eta2_min
do i2=1, nc_eta2+1
   eta1 = eta1_min
   do i1=1, nc_eta1
      df_2d%field%data(i1,i2) = (1+eps*cos(kx*eta1))/(2.*sll_pi)*exp(-.5_f64*eta2*eta2)
      eta1=eta1+delta_eta1
   end do
   eta2 = eta2+delta_eta2
end do

call write_distribution_function ( df_2d )

adv_field_x => new_field_1D_vec1(mesh_v)

bsl_work_x => new_bsl_workspace(adv_field_x)
bsl_work_v => new_bsl_workspace(ex)

! run BSL method using 10 time steps and second order splitting
time = 0.0_f64
n_step = 1000
delta_t = 0.01_f64
SLL_ALLOCATE(nrj(n_step), error)
print*, 'set term x11'

do i_step = 1, n_step

   eta2 = eta2_min
   do i2 = 1, nc_eta2
      adv_field_x%data(:) = eta2
      call bsl_step_1d( bsl_work_x, df_2d%field%data(:,i2), &
                        adv_field_x, delta_t )
      eta2 = eta2 + delta_eta2
   end do

   do i1 = 1, nc_eta1
      rho%data(i1) = sum(df_2d%field%data(i1,:))
   end do

   call solve(poisson, ex, rho)

   do i1 = 1, nc_eta1
      call bsl_step_1d( bsl_work_v, df_2d%field%data(i1,:), &
                        ex, delta_t )
   end do

   time = time + delta_t

   nrj(i_step) = 0.5_f64*log(sum(ex%data*ex%data)*delta_eta1)

   write(*,100) .0,10.,-15.5,1.5
   do j_step = 1, i_step
      print*, (j_step-1)*delta_t, nrj(j_step)
   end do
   print*, 'e'

end do

100 format('p [',f5.1,':',f5.1,'][',f6.1,':',f6.1,'] ''-'' w l')
    
end program landau_2d
