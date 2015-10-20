program vm2d_spectral

#define MPI_MASTER 0
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"

use sll_m_boundary_condition_descriptors
use sll_m_cubic_spline_interpolator_1d
use sll_m_constants, only : &
     sll_pi

#include "sll_fftw.h"

use sll_m_vlasov2d_base
use sll_m_vlasov2d_spectral
use sll_m_poisson_1d_periodic  
use sll_m_poisson_1d_periodic_solver
use sll_m_ampere_vlasov_1d
use sll_m_collective
use sll_m_remapper
use sll_m_init_functions, only: landau_1d, tsi, TSI_CASE, LANDAU_X_CASE

implicit none

type(vlasov2d_spectral)                         :: vlasov 
type(sll_cubic_spline_interpolator_1d), target  :: spl_x1
type(sll_cubic_spline_interpolator_1d), target  :: spl_x2
class(sll_poisson_1d_base),             pointer :: poisson
type(sll_ampere_1d),                    pointer :: ampere 


sll_int32  :: iter 
sll_real64 :: tcpu1, tcpu2

sll_int32  :: prank, comm
sll_int64  :: psize

sll_int32  :: loc_sz_i, loc_sz_j

call sll_boot_collective()
prank = sll_get_collective_rank(sll_world_collective)
psize = sll_get_collective_size(sll_world_collective)
comm  = sll_world_collective%comm

tcpu1 = MPI_WTIME()
if (prank == MPI_MASTER) then
   print*,'MPI Version of slv2d running on ',psize, ' processors'
end if

call initlocal()

call transposexv(vlasov)
call compute_charge(vlasov)
call poisson%compute_E_from_rho( vlasov%ex, vlasov%rho )
call transposevx(vlasov)
call spectral_advection_x(vlasov, 0.5*vlasov%dt)

do iter=1,vlasov%nbiter

   if (iter ==1 .or. mod(iter,vlasov%fdiag) == 0) then 
#ifndef __INTEL_COMPILER
      call write_xmf_file(vlasov,iter/vlasov%fdiag)
#endif
   end if

   call transposexv(vlasov)

   call compute_current(vlasov)

   call solve_ampere(vlasov%dt)

   call advection_v(vlasov, vlasov%dt)

   call transposevx(vlasov)

   call spectral_advection_x(vlasov, vlasov%dt)

   if (mod(iter,vlasov%fthdiag) == 0) then 
      call write_energy(vlasov, iter*vlasov%dt)
   endif

end do

tcpu2 = MPI_WTIME()
if (prank == MPI_MASTER) then
     write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize
end if

!PN call delete(poisson)
call sll_halt_collective()

print*,'PASSED'

!####################################################################################

contains

!####################################################################################

subroutine initlocal()

  use sll_m_init_functions
  
  sll_real64 :: x
  sll_real64 :: vx
  sll_real64 :: v2
  sll_int32  :: i,j
  sll_real64 :: kx
  sll_int32  :: gi, gj
  sll_int32  :: global_indices(2)
  sll_int32  :: psize
  sll_int32  :: error

  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  call read_input_file(vlasov)

  call spl_x1%initialize(vlasov%np_eta1,  &
                         vlasov%eta1_min, &
                         vlasov%eta1_max, &
                         SLL_PERIODIC)

  call spl_x2%initialize(vlasov%np_eta2,  &
                         vlasov%eta2_min, &
                         vlasov%eta2_max, &
                         SLL_PERIODIC)


  call initialize(vlasov,spl_x1,spl_x2,error)

  call compute_local_sizes(vlasov%layout_x,loc_sz_i,loc_sz_j)        

  kx  = 2_f64*sll_pi/(vlasov%nc_eta1*vlasov%delta_eta1)

  do j=1,loc_sz_j
     do i=1,loc_sz_i
        
        global_indices = local_to_global(vlasov%layout_x,(/i,j/)) 
        gi = global_indices(1)
        gj = global_indices(2)
        
        x  = vlasov%eta1_min+(gi-1)*vlasov%delta_eta1
        vx = vlasov%eta2_min+(gj-1)*vlasov%delta_eta2
        
        v2 = vx*vx

        select case(vlasov%num_case)
        case(LANDAU_X_CASE)
            vlasov%f(i,j)= landau_1d(vlasov%eps,kx,x,v2)
        case(TSI_CASE)
            vlasov%f(i,j)= tsi(vlasov%eps,kx,x,vx,v2)
        end select
        
     end do
  end do
  
  ampere => new_ampere_1d( vlasov%eta1_min,  &
                           vlasov%eta1_max,  &
                           vlasov%nc_eta1)
  
  poisson => new_poisson_1d_periodic_solver( vlasov%eta1_min, &
                                             vlasov%eta1_max, &
                                             vlasov%nc_eta1)
  
  
end subroutine initlocal

subroutine solve_ampere(dt)

  sll_real64, intent(in)    :: dt
  
  call sll_solve(ampere, dt, vlasov%jx, vlasov%ex)

end subroutine solve_ampere

end program vm2d_spectral
