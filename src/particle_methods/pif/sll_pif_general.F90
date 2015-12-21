!**************************************************************
!  Author: Jakob Ameres, jakob.ameres@tum.de
!**************************************************************

program sll_pif_general
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"


use sll_m_pif_fieldsolver
use sll_m_timer
use sll_m_sobol
use sll_m_prob
use sll_m_collective
use sll_m_pic_visu
use sll_m_moment_matching
use sll_m_pic_utilities
use sll_m_particle_method_descriptors

implicit none


sll_int32, parameter :: SOBOL_OFFSET   = 10    !Default sobol offset, skip zeros


!------type--------
sll_real64, dimension(:,:), allocatable :: particle
sll_real64, dimension(:), allocatable :: weight_const, prior_weight !initial weight if needed
sll_real64, dimension(:),allocatable :: rk_d, rk_c
type(pif_fieldsolver) :: SOLVER
sll_real64 :: Efield
!stencils
sll_int32, dimension(:),allocatable :: maskx,maskv,maskxw,maskxv
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
sll_int32 :: collisions=SLL_COLLISIONS_NONE
sll_int32 :: controlvariate=SLL_CONTROLVARIATE_NONE      !SLL_CONTROLVARIATE_MAXWELLIAN




!results
sll_real64, dimension(:),allocatable :: kineticenergy,fieldenergy, energy,energy_error,weight_sum,weight_var,moment_error
sll_real64, dimension(:,:),allocatable ::moment,moment_var
sll_comp64, dimension(:), allocatable :: rhs, solution
!--------------

!------MPI---------
sll_int32 :: coll_rank, coll_size


!------MPI----------



sll_int32 ::ierr , idx, jdx


!-----------------
npart=1e6
dt=0.1
tsteps=100
qm=-1 !electrons

dimx=2
boxlen=4*sll_p_pi


call sll_s_boot_collective()
call sll_collective_barrier(sll_v_world_collective)
 !really global variables
 coll_rank = sll_f_get_collective_rank( sll_v_world_collective )
 coll_size = sll_f_get_collective_size( sll_v_world_collective )
 call sll_collective_barrier(sll_v_world_collective)


if (coll_rank==0) then
        print *, "Size of MPI-Collective: ", coll_size
endif


!Determine MPI particle distribution scheme


npart_loc=npart/coll_size

!print*, "#Core ", coll_rank, " handles particles", coll_rank*nparticles +1, "-", (coll_rank+1)*nparticles
!call sll_collective_barrier(sll_v_world_collective)
if (coll_rank==0) print*, "#Total Number of particles: ", npart_loc*coll_size


 
 
!--------------------------------------
call init_particle_masks()


SOLVER%dimx=dimx
!Load some particles


call SOLVER%set_box_len(boxlen)
call SOLVER%init(10)

call set_symplectic_rungekutta_coeffs(3)
call init_diagnostics()

call init_particle()
call init_particle_prior_maxwellian()



!Set weights for strong landau damping
particle(maskw,:)=(1-0.4*sum(cos(0.5*particle(maskx,:) +sll_p_pi/3.0),1))
particle(maskw,:)=particle(maskw,:)/prior_weight

SLL_ALLOCATE(weight_const(npart_loc),ierr)
weight_const=particle(maskw,:)

!Match some moments
do idx=1, size(maskv)
 call match_moment_1D_weight_linear_real64(particle(maskv(idx),:), particle(maskw,:),0.0_f64, 1.0_f64,npart)
end do


  
!Initial field solve
SLL_ALLOCATE(rhs(SOLVER%problemsize()),ierr)  
SLL_ALLOCATE(solution(SOLVER%problemsize()),ierr)
  
if (coll_rank==0)print *, "# TIME                          IMPULSE ERR.(abs.)       ENERGY ERROR(rel.)"
do tstep=1,tsteps
  call symplectic_rungekutta()
        if (coll_rank==0) print *,"#", dt*(tstep-1),abs(moment(1,tstep))/npart, energy_error(tstep)

      !if ( (gnuplot_inline_output.eqv. .true.) .AND. coll_rank==0 .AND. mod(timestep-1,timesteps/100)==0  ) then
!                    call energies_electrostatic_gnuplot_inline(kineticenergy(1:tstep), fieldenergy(1:tstep),&
!   			  moment_error(1:tstep),dt)
       !     else

        !    endif
end do

call write_result()


call sll_s_halt_collective()
!--------------------------------------------------------------
contains

!Fills particle vector with random numbers, of users choice
subroutine init_particle
sll_int32 :: idx
sll_int64 :: seed

SLL_ALLOCATE(particle(2*dimx+1,npart_loc),ierr)

!Generate random numbers
seed=SOBOL_OFFSET + coll_rank*npart_loc
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
call sll_collective_globalsum(sll_v_world_collective, kineticenergy)

fieldenergy(tstep)=abs(dot_product(solution,rhs))
! fieldenergy(tstep)=abs(Efield)
energy(tstep)=kineticenergy(tstep)+fieldenergy(tstep)
energy_error(tstep)=abs(energy(2)-energy(tstep))/abs(energy(1))

do idx=1,size(particle,2)
moment(:,tstep)=moment(:,tstep)+(particle(maskv,idx)*particle(maskw,idx))
end do
call sll_collective_globalsum(sll_v_world_collective, moment(:,tstep))

moment_error(tstep)=sqrt(sum((moment(:,1)-moment(:,tstep))**2))

weight_sum(tstep)=sum(particle(maskw,:))/npart
call sll_collective_globalsum(sll_v_world_collective, weight_sum,0)


end subroutine

subroutine init_particle_prior_maxwellian()
sll_int32 :: idx,jdx

!scale to boxlength
particle(maskx,:)=particle(maskx,:)*boxlen
!load sll_m_gaussian profile
do idx=1,size(particle,2)
 do jdx=1, size(maskv)
  call normal_cdf_inv( particle(maskv(jdx),idx), 0.0_f64 , 1.0_f64, particle(maskv(jdx),idx))
 end do
end do

!set temperature and impulse
 do idx=1, size(maskv)
! print *, sum(particle(maskv(idx),:))/npart
 call match_moment_1D_linear_real64(particle(maskv(idx),:), 0.0_f64, 1.0_f64)
! print *, sum(particle(maskv(idx),:))/npart
  end do

  
SLL_ALLOCATE(prior_weight(1:npart_loc),ierr)
prior_weight=1.0/boxlen**dimx


end subroutine




!Only for separable lagrangian
subroutine symplectic_rungekutta()
  sll_int32 :: rkidx
  sll_real64 :: t !time
  SLL_ASSERT(size(rk_d)==size(rk_c))
  
    
  !Loop over all stages
    do rkidx=1, size(rk_d);
     
     !Set control variate if used
     call update_weight()
     
     !Charge assignement, get right hand side
     rhs=SOLVER%get_rhs_particle(particle(maskxw,:))/npart
     !mpi allreduce
     call sll_collective_globalsum(sll_v_world_collective, rhs)

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
SLL_ALLOCATE(maskxv(dimx*2),ierr)

maskx=(/( idx,idx=1,dimx)/)
maskv=(/( idx,idx=dimx+1,2*dimx)/)
maskw=2*dimx+1
maskxw(1:dimx)=maskx
maskxw(dimx+1)=maskw
maskxv(1:dimx)=maskx
maskxv(dimx+1:2*dimx)=maskv

end subroutine init_particle_masks


!Updates the weights for use with control variate
subroutine update_weight()


if (controlvariate /= 0) then

  if (collisions/=0) then
!particle(maskw,:)=initial_prior -  (particle(maskw,:)-initial   control_variate(particle(maskxv,:))
print *, "not implemented"
else
!no collisions, characteristics are conserved
 particle(maskw,:)=weight_const(:) - control_variate(particle(maskxv,:))/prior_weight(:)
 endif

 endif
end subroutine





function control_variate(particle) result(cv)
sll_real64, dimension(:,:), intent(in) :: particle !(x,v)
sll_real64, dimension(size(particle,2)) :: cv

!Standard maxwellian control variate
SELECT CASE (controlvariate)
   CASE (SLL_CONTROLVARIATE_NONE)
     !This should not happen
      cv=1
    CASE (SLL_CONTROLVARIATE_STANDARD)
      cv=sqrt(2*sll_p_pi)**(-size(maskv))*exp(-0.5*sum(particle(maskv,:),1)**2)
    CASE (SLL_CONTROLVARIATE_MAXWELLIAN)
      cv=sqrt(2*sll_p_pi)**(-size(maskv))*exp(-0.5*sum(particle(maskv,:),1)**2)
END SELECT

end function 



subroutine write_result()
        !character(len=*), intent(in) :: filename
        character(len=*),parameter :: filename="pif"
        
!         sll_real64, dimension(:), intent(in) :: kinetic_energy, electrostatic_energy, impulse,&
!             particleweight_mean,particleweight_var, perror_mean  , perror_var
        integer :: idx,file_id,file_id_err
!         sll_real64,  dimension(size(kinetic_energy)) :: total_energy


        if (coll_rank==0) then

             call sll_new_file_id(file_id, ierr)

            !Write Data File
            !            open(file_id, file = plot_name//"_"//fin//'.dat' )
            open(file_id, file = './'//filename//'.csv')

!             write (file_id, *)  "#Full 1d1v Electrostatic PIC"
!             write (file_id,*)  "#Time steps:", timesteps
!             write (file_id,*)  "#Time stepwidth:", timestepwidth
!             write (file_id,*)  "#Marko particles:", coll_size*nparticles
!             write (file_id,*)  "#Particle Pusher:", particle_pusher
!             write (file_id,*)  "#Finite Elements: 2^(", log(real(mesh_cells,i64))/log(2.0_f64),")"
!             write (file_id,*)  "#Size of MPI Collective: ", coll_size

            write (file_id,*)  ' \"time\", \"kineticenergy\", \"fieldenergy\", \"energy_error\", \"moment_error\"'
				
!             "impulse  ", &
!                 "vthermal  ", "weightmean   ", "weightvar  ", "perrormean  ", "perrorvar   "
            do idx=1,tsteps
                write (file_id,*) (idx-1)*dt,",",kineticenergy(idx),",",fieldenergy(idx),",", &
                        fieldenergy(idx),",", energy_error(idx),",", moment_error(idx)
            enddo
            close(file_id)
            endif
end subroutine

end program sll_pif_general
