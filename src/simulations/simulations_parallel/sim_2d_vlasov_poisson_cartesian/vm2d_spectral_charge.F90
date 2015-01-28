program vm2d_spectral_charge

#define MPI_MASTER 0
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_constants.h"
#include "sll_interpolators.h"
#include "sll_fftw.h"

use init_functions
use sll_vlasov2d_base
use sll_vlasov2d_spectral_charge
use sll_poisson_1d_periodic  
use sll_module_poisson_1d_periodic_solver
use sll_module_ampere_1d_pstd
use sll_collective
use sll_remapper

implicit none

type(vlasov2d_spectral_charge)                  :: vlasov 
type(sll_cubic_spline_interpolator_1d), target  :: spl_x1
type(sll_cubic_spline_interpolator_1d), target  :: spl_x2
class(sll_poisson_1d_base),             pointer :: poisson
type(sll_ampere_1d_pstd),               pointer :: ampere 


sll_int32  :: iter 
sll_real64 :: tcpu1, tcpu2, mass0

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

!f --> ft
call transposexv(vlasov)
call compute_charge(vlasov)
call poisson%compute_E_from_rho( vlasov%ex, vlasov%rho )

!vlasov%exn=vlasov%ex

!ft --> f
call transposevx(vlasov)
call spectral_advection_x(vlasov, 0.5*vlasov%dt)

mass0=sum(vlasov%rho)*vlasov%delta_eta1

print *,'mass init',mass0

do iter=1,vlasov%nbiter

   if (iter ==1 .or. mod(iter,vlasov%fdiag) == 0) then 
      call write_xmf_file(vlasov,iter/vlasov%fdiag)
   end if

   call transposexv(vlasov)

   call compute_charge(vlasov)

   call poisson%compute_E_from_rho( vlasov%ex, vlasov%rho )

   call advection_v(vlasov, vlasov%dt)

   call transposevx(vlasov)

   call spectral_advection_charge_x(vlasov, vlasov%dt)


!   if ( vlasov%va == VA_VALIS .or. vlasov%va == VA_CLASSIC) then 
!
!      !f --> ft, current (this%jx,this%jy), ft-->f
!      call transposexv(vlasov)
!
!      !compute this%jx, this%jy (zero average) at time tn
!      call compute_current(vlasov)
!
!      call transposevx(vlasov)
!
!      !compute vlasov%bz=B^{n+1/2} from Ex^n, Ey^n, B^{n-1/2}  
!      !!!!Attention initialisation B^{-1/2}
!
!      !compute vlasov%bzn=B^n=0.5(B^{n+1/2}+B^{n-1/2})          
!      vlasov%exn=vlasov%ex
!      
!      !compute (vlasov%ex,vlasov%ey)=E^{n+1/2} from vlasov%bzn=B^n
!      call solve_ampere(0.5_f64*vlasov%dt) 
!
!      if (vlasov%va == VA_CLASSIC) then 
!
!         vlasov%jx3=vlasov%jx
!
!      endif
!
!   endif
!
!   !advec x + compute this%jx1
!   call advection_x1(vlasov,0.5_f64*vlasov%dt)
!
!   !advec y + compute this%jy1
!   call advection_x2(vlasov,0.5_f64*vlasov%dt)
!
!
!   if (vlasov%va == VA_OLD_FUNCTION) then 
!
!      !compute rho^{n+1}
!      call compute_charge(vlasov)
!      !compute E^{n+1} via Poisson
!      call poisson%compute_E_from_rho( vlasov%ex, vlasov%rho )
!
!   endif
!
!   call transposevx(vlasov)
!
!   !copy jy^{**}
!   call advection_x2(vlasov,0.5_f64*vlasov%dt)
!   
!   !copy jx^*
!   vlasov%jx=vlasov%jx1
!   !advec x + compute this%jx1
!   call advection_x1(vlasov,0.5_f64*vlasov%dt)
!
!   if (vlasov%va == VA_VALIS) then 
!
!      !compute the good jx current
!      vlasov%jx=0.5_f64*(vlasov%jx+vlasov%jx1)
!      print *,'sum jx ', &
!         sum(vlasov%jx)*vlasov%delta_eta1*vlasov%delta_eta2, &
!         maxval(vlasov%jx)
!
!      !compute E^{n+1} from B^{n+1/2}, vlasov%jx, vlasov%jy, E^n
!      vlasov%ex  = vlasov%exn
!      
!      call solve_ampere(vlasov%dt) 
!
!      !copy ex and ey at t^n for the next loop
!      vlasov%exn = vlasov%ex
!
!   else if (vlasov%va == VA_CLASSIC) then 
!
!      !f --> ft, current (this%jx,this%jy), ft-->f
!      call transposexv(vlasov)
!      !compute this%jx, this%jy (zero average) at time tn
!      call compute_current(vlasov)
!
!      call transposevx(vlasov)
!
!      !compute J^{n+1/2}=0.5*(J^n+J^{n+1})
!      vlasov%jx=0.5_f64*(vlasov%jx+vlasov%jx3)
!
!      !compute E^{n+1} from B^{n+1/2}, vlasov%jx, vlasov%jy, E^n
!      vlasov%ex  = vlasov%exn
!      
!      call solve_ampere(vlasov%dt) 
!
!      !copy ex and ey at t^n for the next loop
!      vlasov%exn = vlasov%ex
!
!   else if (vlasov%va == VA_VLASOV_POISSON) then 
!
!      call transposexv(vlasov)
!      !compute rho^{n+1}
!      call compute_charge(vlasov)
!      call transposevx(vlasov)
!
!      !compute E^{n+1} via Poisson
!      call poisson%compute_E_from_rho( vlasov%ex, vlasov%rho )
!
!      !print *,'verif charge conservation', &
!      !             maxval(vlasov%exn-vlasov%ex), &
!      !             maxval(vlasov%eyn-vlasov%ey)
!   
!   else if (vlasov%va==VA_OLD_FUNCTION) then 
!
!      !recompute the electric field at time (n+1) for diagnostics
!      call transposexv(vlasov)
!      !compute rho^{n+1}
!      call compute_charge(vlasov)
!      call transposevx(vlasov)
!      !compute E^{n+1} via Poisson
!      call poisson%compute_E_from_rho( vlasov%ex, vlasov%rho )
!   endif

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

  use init_functions
  
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
  
  ampere => new_ampere_1d_pstd( vlasov%eta1_min,  &
                                vlasov%eta1_max,  &
                                vlasov%nc_eta1)
  
  poisson => new_poisson_1d_periodic_solver( vlasov%eta1_min, &
                                             vlasov%eta1_max, &
                                             vlasov%nc_eta1)
  
  
end subroutine initlocal

subroutine solve_ampere(dt)

  sll_real64, intent(in)    :: dt
  
  print*, 'solve ampere'
  call sll_solve(ampere, dt, vlasov%jx, vlasov%ex)

end subroutine solve_ampere

end program vm2d_spectral_charge
