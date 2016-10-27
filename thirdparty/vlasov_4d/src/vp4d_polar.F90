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

type(vlasov4d_polar)    :: sim
type(sll_poisson_polar) :: poisson 

type(sll_t_cubic_spline_interpolator_2d), target :: spl_x1x2
type(sll_cubic_spline_interpolator_1d), target :: spl_x3
type(sll_cubic_spline_interpolator_1d), target :: spl_x4

sll_int32 :: i, j, k, l
sll_int32 :: loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4
sll_int32 :: error

call sll_s_boot_collective() 

call read_input_file(sim)

call spl_x1x2%initialize( sim%nc_eta1+1, sim%nc_eta2+1, &
                          sim%eta1_min, sim%eta1_max,   &
                          sim%eta2_min, sim%eta2_max,   &
                          sll_p_periodic, sll_p_periodic )

call spl_x3%initialize(sim%nc_eta3+1,sim%eta3_min,sim%eta3_max,sll_p_periodic)

call spl_x4%initialize(sim%nc_eta4+1,sim%eta4_min,sim%eta4_max,sll_p_periodic)

call initialize_vp4d_polar( sim, spl_x1x2, spl_x3, spl_x4)

!call initialize( poisson,       &
!                 sim%layout_x1, &
!                 sim%layout_x2, &
!                 sim%eta1_min,  &
!                 sim%eta1_max,  &
!                 sim%nc_eta1,   &
!                 sim%nc_eta2,   &
!                 SLL_DIRICHLET, &
!                 SLL_DIRICHLET)

call sll_o_compute_local_sizes(sim%layout_v, &
                         loc_sz_x1,    &
                         loc_sz_x2,    &
                         loc_sz_x3,    &
                         loc_sz_x4)

do l=1,loc_sz_x4 
do k=1,loc_sz_x3
do j=1,loc_sz_x2 
do i=1,loc_sz_x1 

   sim%ft(i,j,k,l)= exp(-0.5_f64*((sim%x1(i,j)-5.)**2+(sim%x2(i,j))**2))

end do
end do
end do
end do



do itime = 1, sim%nbiter

   if(prank == MPI_MASTER) &
      print *, itime, ' of ', sim%nbiter, sim%dt

   call plot_f(sim)

   call sll_o_apply_remap_4d( sim%v_to_x, sim%ft, sim%f )

   call advection_x1x2(sim, sim%dt)

   !call compute_charge_density(sim)
   !call plot_rho(sim)
   !call solve(sim%poisson,sim%rho,sim%phi_x2)
   !call plot_phi(sim)
   !call plot_ft(sim)

   call sll_o_apply_remap_4d(sim%x_to_v,sim%f,sim%ft)

   !call apply_remap_2D(sim%rmp_x2x1,sim%phi_x2,sim%phi_x1)
   !call compute_electric_fields_eta1(sim)
   !call apply_remap_2D(sim%rmp_x1x2,sim%phi_x1,sim%phi_x2)
   !call compute_electric_fields_eta2(sim)
   !call advection_x3(sim,sim%dt)
   !call advection_x4(sim,sim%dt)

end do

call sll_o_delete( spl_x1x2 )
call sll_o_delete( spl_x3 )
call sll_o_delete( spl_x4 )

call sll_s_halt_collective()

end program vlasov_poisson_4d_polar
