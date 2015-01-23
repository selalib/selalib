program vm2d_spectral_charge

#define MPI_MASTER 0
#include "selalib-mpi.h"

use init_functions
use sll_vlasov2d_base
use sll_vlasov2d_spectral_charge
use sll_poisson_1d_periodic  
use sll_module_poisson_1d_periodic_solver
use sll_module_ampere_1d_pstd

implicit none

type(vlasov2d_spectral_charge)                 :: vlasov2d 
type(sll_cubic_spline_interpolator_1d), target :: spl_x2
class(sll_poisson_1d_base), pointer            :: poisson
type(sll_ampere_1d_pstd),   pointer            :: ampere 


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
call transposexv(vlasov2d)

call compute_charge(vlasov2d)
call poisson%compute_E_from_rho( vlasov2d%ex, vlasov2d%rho )

vlasov2d%exn=vlasov2d%ex

!ft --> f
call transposevx(vlasov2d)

mass0=sum(vlasov2d%rho)*vlasov2d%delta_eta1

print *,'mass init',mass0

do iter=1,vlasov2d%nbiter

   if (iter ==1 .or. mod(iter,vlasov2d%fdiag) == 0) then 
      call write_xmf_file(vlasov2d,iter/vlasov2d%fdiag)
   end if

   if ( vlasov2d%va == VA_VALIS .or. vlasov2d%va == VA_CLASSIC) then 

      !f --> ft, current (this%jx,this%jy), ft-->f
      call transposexv(vlasov2d)

      !compute this%jx, this%jy (zero average) at time tn
      call compute_current(vlasov2d)

      call transposevx(vlasov2d)

      !compute vlasov2d%bz=B^{n+1/2} from Ex^n, Ey^n, B^{n-1/2}  
      !!!!Attention initialisation B^{-1/2}

      !compute vlasov2d%bzn=B^n=0.5(B^{n+1/2}+B^{n-1/2})          
      vlasov2d%exn=vlasov2d%ex
      
      !compute (vlasov2d%ex,vlasov2d%ey)=E^{n+1/2} from vlasov2d%bzn=B^n
      call solve_ampere(0.5_f64*vlasov2d%dt) 

      if (vlasov2d%va == VA_CLASSIC) then 

         vlasov2d%jx3=vlasov2d%jx

      endif

   endif

   !advec x + compute this%jx1
   call advection_x1(vlasov2d,0.5_f64*vlasov2d%dt)

   !advec y + compute this%jy1
   call advection_x2(vlasov2d,0.5_f64*vlasov2d%dt)

   call transposexv(vlasov2d)

   if (vlasov2d%va == VA_OLD_FUNCTION) then 

      !compute rho^{n+1}
      call compute_charge(vlasov2d)
      !compute E^{n+1} via Poisson
      call poisson%compute_E_from_rho( vlasov2d%ex, vlasov2d%rho )

   endif

   call transposevx(vlasov2d)

   !copy jy^{**}
   call advection_x2(vlasov2d,0.5_f64*vlasov2d%dt)
   
   !copy jx^*
   vlasov2d%jx=vlasov2d%jx1
   !advec x + compute this%jx1
   call advection_x1(vlasov2d,0.5_f64*vlasov2d%dt)

   if (vlasov2d%va == VA_VALIS) then 

      !compute the good jx current
      vlasov2d%jx=0.5_f64*(vlasov2d%jx+vlasov2d%jx1)
      print *,'sum jx ', &
         sum(vlasov2d%jx)*vlasov2d%delta_eta1*vlasov2d%delta_eta2, &
         maxval(vlasov2d%jx)

      !compute E^{n+1} from B^{n+1/2}, vlasov2d%jx, vlasov2d%jy, E^n
      vlasov2d%ex  = vlasov2d%exn
      
      call solve_ampere(vlasov2d%dt) 

      !copy ex and ey at t^n for the next loop
      vlasov2d%exn = vlasov2d%ex

   else if (vlasov2d%va == VA_CLASSIC) then 

      !f --> ft, current (this%jx,this%jy), ft-->f
      call transposexv(vlasov2d)
      !compute this%jx, this%jy (zero average) at time tn
      call compute_current(vlasov2d)

      call transposevx(vlasov2d)

      !compute J^{n+1/2}=0.5*(J^n+J^{n+1})
      vlasov2d%jx=0.5_f64*(vlasov2d%jx+vlasov2d%jx3)

      !compute E^{n+1} from B^{n+1/2}, vlasov2d%jx, vlasov2d%jy, E^n
      vlasov2d%ex  = vlasov2d%exn
      
      call solve_ampere(vlasov2d%dt) 

      !copy ex and ey at t^n for the next loop
      vlasov2d%exn = vlasov2d%ex

   else if (vlasov2d%va == VA_VLASOV_POISSON) then 

      call transposexv(vlasov2d)
      !compute rho^{n+1}
      call compute_charge(vlasov2d)
      call transposevx(vlasov2d)

      !compute E^{n+1} via Poisson
      call poisson%compute_E_from_rho( vlasov2d%ex, vlasov2d%rho )

      !print *,'verif charge conservation', &
      !             maxval(vlasov2d%exn-vlasov2d%ex), &
      !             maxval(vlasov2d%eyn-vlasov2d%ey)
   
   else if (vlasov2d%va==VA_OLD_FUNCTION) then 

      !recompute the electric field at time (n+1) for diagnostics
      call transposexv(vlasov2d)
      !compute rho^{n+1}
      call compute_charge(vlasov2d)
      call transposevx(vlasov2d)
      !compute E^{n+1} via Poisson
      call poisson%compute_E_from_rho( vlasov2d%ex, vlasov2d%rho )
   endif

   if (mod(iter,vlasov2d%fthdiag) == 0) then 
      call write_energy(vlasov2d, iter*vlasov2d%dt)
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

  call read_input_file(vlasov2d)

  call spl_x2%initialize(vlasov2d%np_eta2,  &
                         vlasov2d%eta2_min, &
                         vlasov2d%eta2_max, &
                         SLL_PERIODIC)


  call initialize(vlasov2d,spl_x2,error)

  call compute_local_sizes(vlasov2d%layout_x,loc_sz_i,loc_sz_j)        

  kx  = 2_f64*sll_pi/(vlasov2d%nc_eta1*vlasov2d%delta_eta1)

  do j=1,loc_sz_j
     do i=1,loc_sz_i
        
        global_indices = local_to_global(vlasov2d%layout_x,(/i,j/)) 
        gi = global_indices(1)
        gj = global_indices(2)
        
        x  = vlasov2d%eta1_min+(gi-1)*vlasov2d%delta_eta1
        vx = vlasov2d%eta2_min+(gj-1)*vlasov2d%delta_eta2
        
        v2 = vx*vx

        select case(vlasov2d%num_case)
        case(LANDAU_X_CASE)
            vlasov2d%f(i,j)= landau_1d(vlasov2d%eps,kx,x,v2)
        case(TSI_CASE)
            vlasov2d%f(i,j)= tsi(vlasov2d%eps,kx,x,vx,v2)
        end select
        
     end do
  end do
  
   ampere => new_ampere_1d_pstd( vlasov2d%eta1_min,  &
                                 vlasov2d%eta1_max,  &
                                 vlasov2d%nc_eta1)
  
   poisson => new_poisson_1d_periodic_solver( vlasov2d%eta1_min, &
                                              vlasov2d%eta1_max, &
                                              vlasov2d%nc_eta1)
  
  
end subroutine initlocal

subroutine solve_ampere(dt)

  sll_real64, intent(in)    :: dt
  
  print*, 'solve ampere'
  call sll_solve(ampere, vlasov2d%ex, dt, vlasov2d%jx)

end subroutine solve_ampere

end program vm2d_spectral_charge
