#define MPI_MASTER 0

! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D: x, y, vx, vy (or x1, x2, x3, x4) with arbitrary coordinate 
!   transformation
!   in the x,y variables.
! - parallel

program vlasov_poisson_4d_polar
#include "selalib-mpi.h"
use sll_vlasov4d_polar
use sll_poisson_polar_parallel

implicit none

type(vlasov4d_polar)                :: sim
type(sll_logical_mesh_2d), pointer  :: mx
type(sll_logical_mesh_2d), pointer  :: mv
type(sll_poisson_polar)             :: poisson 

class(sll_coordinate_transformation_2d_base), pointer :: transformation

type(cubic_spline_2d_interpolator), target :: spl_x1x2
type(cubic_spline_1d_interpolator), target :: spl_x3
type(cubic_spline_1d_interpolator), target :: spl_x4

sll_int32 :: k, l
sll_int32 :: loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4
sll_int32 :: ntime

sll_real64 :: dt

call sll_boot_collective() ! Wrap this up somewhere else

#define nc_eta1 64
#define nc_eta2 64
#define nc_eta3 32
#define nc_eta4 32

eta1_min= 2.0_f64 ; eta1_max= 8.0_f64
eta2_min=.0_f64   ; eta2_max=2.0_f64*sll_pi

eta3_min=-6.0_f64 ; eta3_max=6.0_f64
eta4_min=-6.0_f64 ; eta4_max=6.0_f64

! logical mesh for velocity coordinates
mv => new_logical_mesh_2d( nc_eta3, nc_eta4,    &
      eta1_min=-6.0_f64, eta1_max=6.0_f64,      &
      eta2_min=-6.0_f64, eta2_max=6.0_f64)

! logical mesh for space coordinates
mx => new_logical_mesh_2d( nc_eta1, nc_eta2,    & 
      eta1_min= 2.0_f64, eta1_max= 8.0_f64,     &
      eta2_min=.0_f64,   eta2_max=2.0_f64*sll_pi)

! logical mesh for velocity coordinates
mv => new_logical_mesh_2d( nc_eta3, nc_eta4,     &
      eta1_min=-6.0_f64, eta1_max=6.0_f64,       &
      eta2_min=-6.0_f64, eta2_max=6.0_f64)


call spl_x1x2%initialize( nc_eta1+1, nc_eta2+1,  &
                          eta1_min, eta1_max,    &
                          eta2_min, eta2_max,    &
                          SLL_PERIODIC, SLL_PERIODIC )

call spl_x3%initialize(nc_eta3+1, eta3_min, eta3_max, SLL_PERIODIC)

call spl_x4%initialize(nc_eta4+1, eta4_min, eta4_max, SLL_PERIODIC)

! coordinate transformation associated with space coordinates
transformation => new_coordinate_transformation_2d_analytic( &
       "analytic_polar_transformation", &
       mx, &
       polar_x1, &
       polar_x2, &
       polar_jac11, &
       polar_jac12, &
       polar_jac21, &
       polar_jac22 )

call initialize_vp4d_polar( sim, mx, mv, transformation, &
                            spl_x1x2, spl_x3,spl_x4)

call initialize( poisson,       &
                 sim%layout_x1, &
                 sim%layout_x2, &
                 eta1_min,      &
                 eta2_max,      &
                 sim%nc_x1,     &
                 sim%nc_x2,     &
                 SLL_DIRICHLET, &
                 SLL_DIRICHLET)


call compute_local_sizes_4d(sim%layout_v, &
                            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4)

do l=1,loc_sz_x4 
do k=1,loc_sz_x3

   sim%ft(:,:,k,l)= exp(-0.5_f64*((sim%x1-5.)**2+(sim%x2)**2))

end do
end do

ntime = 100
dt    = 0.01
do itime = 1, ntime

   if(prank == MPI_MASTER) &
      print *, 'Starting iteration ', itime, ' of ', ntime

   call plot_f(sim)

   call apply_remap_4D( sim%v_to_x, sim%ft, sim%f )

   call advection_x1x2(sim, dt)

!   call compute_charge_density( sim )
!
!   call plot_rho(sim)
!
!   call solve(sim%poisson, sim%rho, sim%phi_x2)
!
!   call plot_phi(sim)
!
!   call plot_ft(sim)

   call apply_remap_4D( sim%x_to_v, sim%f, sim%ft )

   !call apply_remap_2D( sim%rmp_x2x1, sim%phi_x2, sim%phi_x1 )
   !call compute_electric_fields_eta1( sim )
   !call apply_remap_2D( sim%rmp_x1x2, sim%phi_x1, sim%phi_x2 )
   !call compute_electric_fields_eta2( sim )
   !call advection_x3(sim,sim%dt)
   !call advection_x4(sim,sim%dt)


end do ! next time step 

call delete( spl_x1x2 )
call delete( spl_x3 )
call delete( spl_x4 )

call sll_halt_collective()

end program vlasov_poisson_4d_polar
