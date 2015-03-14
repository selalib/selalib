program pif_fieldsolver_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"


use sll_pif_fieldsolver
use sll_timer
use sll_sobol
use sll_prob
use sll_collective
use sll_visu_pic

implicit none


sll_int32, parameter :: SOBOL_OFFSET   = 10    !Default sobol offset, skip zeros


!------type--------
sll_real64, dimension(:,:), allocatable :: particle
sll_real64, dimension(:),allocatable :: rk_d, rk_c
type(pif_fieldsolver) :: SOLVER
sll_real64 :: Efield
!stencils
sll_int32, dimension(:),allocatable :: maskx,maskv,maskxw
sll_int32 :: maskw
sll_int32 :: tstep

!----Problem parameter----
sll_int32 :: dimx
sll_real64 :: boxlen
sll_real64 :: qm

!----Solver parameters----
sll_real64 :: dt
sll_int32 :: tsteps
sll_int32 :: npart, npart_loc


!results
sll_real64, dimension(:),allocatable :: kineticenergy,fieldenergy, energy,energy_error,weight_sum,weight_var,moment_error
sll_real64, dimension(:,:),allocatable ::moment,moment_var
sll_comp64, dimension(:), allocatable :: rhs, solution
!--------------

!------MPI---------
sll_int32 :: coll_rank, coll_size


!------MPI----------



sll_int32 ::ierr 


!-----------------
npart=5e6
dt=0.1
tsteps=100
qm=-1 !electrons

dimx=3
boxlen=4*sll_pi


call sll_boot_collective()
call sll_collective_barrier(sll_world_collective)
 !really global variables
 coll_rank = sll_get_collective_rank( sll_world_collective )
 coll_size = sll_get_collective_size( sll_world_collective )
 call sll_collective_barrier(sll_world_collective)



!--------------------------------------
call init_particle_masks()


SOLVER%dimx=dimx
npart_loc=npart
!Load some particles


call SOLVER%set_box_len(boxlen)
call SOLVER%init(2)

call set_symplectic_rungekutta_coeffs(3)
call init_diagnostics()

call init_particle()
call init_particle_prior_maxwellian()

!Set weights for strong landau damping
particle(maskw,:)=(1-0.4*sum(cos(0.5*particle(maskx,:) +sll_pi/3.0),1))*boxlen**dimx
  
  
!Initial field solv


SLL_ALLOCATE(rhs(SOLVER%problemsize()),ierr)  
SLL_ALLOCATE(solution(SOLVER%problemsize()),ierr)
  
print *, "# TIME                          IMPULS ERR.(abs.)       ENERGY ERROR(rel.)"
do tstep=1,tsteps
  call symplectic_rungekutta()
        print *,"#", dt*(tstep-1), moment_error(tstep), energy_error(tstep)

      !if ( (gnuplot_inline_output.eqv. .true.) .AND. coll_rank==0 .AND. mod(timestep-1,timesteps/100)==0  ) then
!                    call energies_electrostatic_gnuplot_inline(kineticenergy(1:tstep), fieldenergy(1:tstep),&
!   			  moment_error(1:tstep),dt)
       !     else

        !    endif
end do


call sll_halt_collective()
!--------------------------------------------------------------
contains

!Fills particle vector with random numbers, of users choice
subroutine init_particle
sll_int32 :: idx
sll_int64 :: seed
SLL_ALLOCATE(particle(2*dimx+1,npart_loc),ierr)

!Generate random numbers
seed=SOBOL_OFFSET
do idx=1,npart_loc
!call i8_sobol_generate ( int(dimx,8) , npart, SOBOL_OFFSET , particle(1:2*dimx,:))            
call i8_sobol( int(2*dimx,8), seed, particle(1:2*dimx,idx))
end do
end subroutine

!allocates space
subroutine init_diagnostics
SLL_CLEAR_ALLOCATE(kineticenergy(tsteps),ierr)
SLL_CLEAR_ALLOCATE(fieldenergy(tsteps),ierr)
SLL_CLEAR_ALLOCATE(energy(tsteps),ierr)
SLL_CLEAR_ALLOCATE(energy_error(tsteps),ierr)
SLL_CLEAR_ALLOCATE(moment(dimx,tsteps),ierr)
SLL_CLEAR_ALLOCATE(moment_error(tsteps),ierr)
SLL_CLEAR_ALLOCATE(moment_var(dimx,tsteps),ierr)
SLL_CLEAR_ALLOCATE(weight_sum(tsteps),ierr)
SLL_CLEAR_ALLOCATE(weight_var(tsteps),ierr)
end subroutine


!reads tstep
subroutine calculate_diagnostics
sll_int32 :: idx

kineticenergy(tstep)=sum(sum(particle(maskv,:)**2,1)*particle(maskw,:))/npart
fieldenergy(tstep)=abs(dot_product(solution,rhs))
! fieldenergy(tstep)=abs(Efield)
energy(tstep)=kineticenergy(tstep)+fieldenergy(tstep)
energy_error(tstep)=abs(energy(2)-energy(tstep))/abs(energy(1))

do idx=1,size(particle,2)
moment(:,tstep)=moment(:,tstep)+(particle(maskv,idx)*particle(maskw,idx))
end do
moment_error(tstep)=sqrt(sum((moment(:,1)-moment(:,tstep))**2))

end subroutine

subroutine init_particle_prior_maxwellian()
sll_int32 :: idx,jdx

!scale to boxlength
particle(maskx,:)=particle(maskx,:)*boxlen
!load gaussian profile
do idx=1,size(particle,2)
 do jdx=1, size(maskv)
  call normal_cdf_inv( particle(maskv(jdx),idx), 0.0_f64 , 1.0_f64, particle(maskv(jdx),idx))
 end do
end do
end subroutine




!Only for separable lagrangian
subroutine symplectic_rungekutta()
  sll_int32 :: rkidx
  sll_real64 :: t !time
  SLL_ASSERT(size(rk_d)==size(rk_c))
  
    
  !Loop over all stages
    do rkidx=1, size(rk_d);
     
     !Charge assignement, get right hand side
     rhs=SOLVER%get_rhs_particle(particle(maskxw,:))/npart
     !mpi allreduce
     
     solution=SOLVER%solve_poisson(rhs)
     
     if (rkidx==1) then
      call calculate_diagnostics()
     endif

     particle(maskv,:)=particle(maskv,:) -  &
          rk_c(rkidx)*dt*qm*(-SOLVER%eval_gradient(particle(maskx,:),solution));
     particle(maskx,:)=particle(maskx,:) +  rk_d(rkidx)*dt*particle(maskv,:);
        
     !   xx=mod(xx,repmat(L,1,npart));
        
      !  CV_new=CV(xx,vx);
       ! weight=forward_weight_CV(weight, CV_old, CV_new);CV_old=CV_new;
        
         t=(tstep-1)*dt+sum(rk_d(1:rkidx))*dt;
    end do
     t=(tstep)*dt;

end subroutine


subroutine set_symplectic_rungekutta_coeffs(rk_order)
sll_int32, intent(in) :: rk_order
sll_real64 :: rk4sx=((2**(1/3) +2**(-1/3)-1)/6)

SLL_ALLOCATE(rk_c(rk_order),ierr)
SLL_ALLOCATE(rk_d(rk_order),ierr)

SELECT CASE (rk_order)
   CASE (1)
    !euler not symplectic
    rk_d=1
    rk_c=1
   CASE (2)
      rk_d=(/0.5, 0.5 /)
      rk_c=(/0.0, 1.0/) 
   CASE (3)
     rk_d(:)=(/ 2.0/3.0, -2.0/3.0, 1.0 /) 
     rk_c(:)=(/ 7.0/24.0,3/4.0,-1.0/24.0 /)  
   CASE (4)
      rk_d=(/ 2.0*rk4sx+1.0 , -4.0*rk4sx-1.0, 2.0*rk4sx+1.0, 0.0_f64/) 
      rk_c=(/ rk4sx + 0.5 , -rk4sx, -rk4sx, rk4sx +0.5 /) 
END SELECT

print *, "Symplectic Runge Kutta of order:", rk_order
print *, "D_COEFFS:", rk_d
print *, "C_COEFFS:", rk_c

end subroutine set_symplectic_rungekutta_coeffs


subroutine init_particle_masks()
 sll_int32 :: idx

 !allocate stencils
SLL_ALLOCATE(maskx(dimx),ierr)
SLL_ALLOCATE(maskxw(dimx+1),ierr)
SLL_ALLOCATE(maskv(dimx),ierr)

maskx=(/( idx,idx=1,dimx)/)
maskv=(/( idx,idx=dimx+1,2*dimx)/)
maskw=2*dimx+1
maskxw(1:dimx)=maskx
maskxw(dimx+1)=maskw

end subroutine init_particle_masks


end program