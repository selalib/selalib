!**************************************************************
!  Author: Jakob Ameres, jakob.ameres@tum.de
!**************************************************************

module sll_general_vlasov_poisson_pif
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
use sll_moment_matching
use sll_pic_utilities
use sll_particle_method_descriptors
use sll_simulation_base
use sll_wedge_product_generaldim
use sll_descriptors
implicit none

! abstract interface
!         function electric_field_general(x,t) result(E)
!             use sll_working_precision
!             sll_real64, dimension(:,:),intent(in) :: x 
!             sll_real64, dimension(:),intent(in) :: t
!             sll_real64, dimension(size(x,1),size(x,2)) :: E
!         endfunction
! endinterface
! 
!  
! !Dummy zero field
! function zero_electric_field_general(x,t) result(E)
! 	    sll_real64, dimension(:,:),intent(in) :: x 
!             sll_real64, dimension(:),intent(in) :: t
!             sll_real64, dimension(size(x,1),size(x,2)) :: E
!             E=0
! end function
!  
 
 
 sll_int32, parameter :: PIF_INTEGRATOR_RK1S=1
 sll_int32, parameter :: PIF_INTEGRATOR_RK2S=2
 sll_int32, parameter :: PIF_INTEGRATOR_RK3S=3
 sll_int32, parameter :: PIF_INTEGRATOR_RK4S=4
 sll_int32, parameter :: PIF_INTEGRATOR_BORIS1=4
 
 sll_int32, parameter :: PIF_POISSON=1
 sll_int32, parameter :: PIF_QUASINEUTRAL_RHO=2
 
 
 
type, extends(sll_simulation_base_class) :: &
                sll_simulation_general_vlasov_poisson_pif
     
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
sll_real64, dimension(:), allocatable :: boxlen !L
sll_real64 :: qm
sll_real64, dimension(:), allocatable :: kmode !k
sll_real64 :: eps !size of disturbance
sll_real64, dimension(:), allocatable :: B0 !constant magentic field vector

!----Solver parameters----
sll_real64 :: dt
sll_int32 :: tsteps
sll_int32 :: npart, npart_loc
sll_int32 :: time_integrator_order
sll_int32 :: collisions=SLL_COLLISIONS_NONE
sll_int32 :: controlvariate=SLL_CONTROLVARIATE_NONE      !SLL_CONTROLVARIATE_MAXWELLIAN
sll_int32 :: momentmatch=SLL_MOMENT_MATCH_NONE
type(sll_vlasovpoisson_sim) :: testcase=SLL_LANDAU_SUM
sll_int32 :: RND_OFFSET   = 10    !Default sobol offset, skip zeros
sll_int32 :: num_modes = 1
!results
sll_real64, dimension(:),allocatable :: kineticenergy,fieldenergy, energy,energy_error,&
               weight_sum,weight_var,moment_error, l2potential
sll_real64, dimension(:,:),allocatable ::moment,moment_var
sll_comp64, dimension(:), allocatable :: rhs, solution
!--------------

!------MPI---------
sll_int32 :: coll_rank, coll_size
!------MPI---------

   contains
     
     procedure, pass(sim) :: init_particle => init_particle_generalvp_pif
     procedure, pass(sim) :: init_from_file => init_file_generalvp_pif
     procedure, pass(sim) :: write_result=>write_result_generalvp_pif
     procedure, pass(sim) :: run=>run_generalvp_pif
     procedure, pass(sim) :: update_weight=>update_weight_generalvp_pif
     procedure, pass(sim) :: init_diagnostics=>init_diagnostics_generalvp_pif
     procedure, pass(sim) :: init_particle_masks=>init_particle_masks_generalvp_pif
     procedure, pass(sim) :: init_particle_prior_maxwellian=>init_particle_prior_maxwellian_generalvp_pif
     procedure, pass(sim) :: set_symplectic_rungekutta_coeffs=>set_symplectic_rungekutta_coeffs_generalvp_pif
     procedure, pass(sim) :: symplectic_rungekutta=>symplectic_rungekutta_generalvp_pif
     procedure, pass(sim) :: YSun_g2h=>YSun_g2h_generalvp_pif
     procedure, pass(sim) :: calculate_diagnostics=>calculate_diagnostics_generalvp_pif
     procedure, pass(sim) :: control_variate=>control_variate_generalvp_pif
       procedure, pass(sim) :: vxB=>v_cross_B_generalvp_pif
     !Magnetic field, can be changed
     procedure, pass(sim) :: B=>Bstandard_generalvp_pif
     !>Unit vector of B
!      procedure, pass(sim) :: Bunit=
      procedure, pass(sim) :: E=>Ezero_generalvp_pif !Electric field
     
     !Loading
     procedure, pass(sim) :: load_landau_sum=>load_landau_sum_generalvp_pif
     procedure, pass(sim) :: load_landau_prod=>load_landau_prod_generalvp_pif
     procedure, pass(sim) :: load_landau_diag=>load_landau_diag_generalvp_pif
end type

!--------------------------------------------------------------

contains

!Gives back the v cross B
function v_cross_B_generalvp_pif(sim, x,v, time) result(vxB)
 class(sll_simulation_general_vlasov_poisson_pif), intent(in) :: sim
 sll_real64, dimension(:,:), intent(in) :: x !position for B
 sll_real64, dimension(:,:), intent(in) :: v !
 sll_real64, intent(in) :: time
 sll_real64, dimension(size(v,1),size(v,2)) :: vxB
 
 SELECT CASE (sim%dimx)
   CASE (1)
    vxB=0
   CASE (2)
      !det( [v;B])*v
      vxB=cross_product_2D(v, sim%B(x,time))      
   CASE (3)   
      vxB=cross_product_3D(v, sim%B(x,time))     
   CASE default
   vxB=0
END SELECT
 
end function


!Define standard zero fields
function Bzero_generalvp_pif(sim, x, time) result(B)
 class(sll_simulation_general_vlasov_poisson_pif), intent(in) :: sim
 sll_real64, dimension(:,:), intent(in) :: x !position for B
 sll_real64, intent(in) :: time
 sll_real64, dimension(size(x,1),size(x,2)) :: B
 B=0
end function

function Ezero_generalvp_pif(sim, x, time) result(E)
 class(sll_simulation_general_vlasov_poisson_pif), intent(in) :: sim
 sll_real64, dimension(:,:), intent(in) :: x !position for B
 sll_real64, intent(in) :: time
 sll_real64, dimension(size(x,1),size(x,2)) :: E
 E=0
end function


function E_KEENWAVE_2D(sim, x, time) result(E)
 class(sll_simulation_general_vlasov_poisson_pif), intent(in) :: sim
 sll_real64, dimension(:,:), intent(in) :: x !position for B
 sll_real64, intent(in) :: time
 sll_real64, dimension(size(x,1),size(x,2)) :: E

 sll_real64 :: t0=0
 sll_real64 :: tL=69
 sll_real64 :: tR=307
 sll_real64 :: tomegaL=20
 sll_real64 :: tomegaR=20
 sll_real64 :: omega=0.37
 sll_real64 :: Emax=0.2
 
 sll_real64 :: k=0.26
 sll_real64 :: L
 sll_real64 :: KEENeps
 tomegaR=tomegaL;
 L=2*sll_pi/k;
 KEENeps=0.5*(tanh( (t0-tL)/tomegaL) - tanh( (t0-tR)/tomegaR))
 E=Emax*k*(0.5*(tanh( (time-tL)/tomegaL) -tanh( (time-tR)/tomegaR)) -KEENeps)*sin(k*x - omega*time)/(1-KEENeps);
end function




function Bstandard_generalvp_pif(sim, x, time) result(B)
 class(sll_simulation_general_vlasov_poisson_pif), intent(in) :: sim
 sll_real64, dimension(:,:), intent(in) :: x !position for B
 sll_real64, intent(in) :: time
 sll_real64, dimension(size(x,1),size(x,2)) :: B
sll_int32 :: idx
 if (allocated(sim%B0)) then
   do idx=1,size(x,2)
     B(:,idx)=sim%B0(:)
   end do
 else
   B=0
 endif
 
!  !leads in three dimensions to (0,0,1) 
!  B(sim%dimx,:)=1
end function

subroutine run_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
  sll_int32 ::ierr , idx, tstep

call sll_collective_barrier(sll_world_collective)
 !really global variables
 sim%coll_rank = sll_get_collective_rank( sll_world_collective )
 sim%coll_size = sll_get_collective_size( sll_world_collective )
 call sll_collective_barrier(sll_world_collective)


if (sim%coll_rank==0) then
        print *, "Size of MPI-Collective: ", sim%coll_size
endif


!Determine MPI particle distribution scheme
sim%npart_loc=sim%npart/sim%coll_size

!print*, "#Core ", coll_rank, " handles particles", coll_rank*nparticles +1, "-", (coll_rank+1)*nparticles
!call sll_collective_barrier(sll_world_collective)
if (sim%coll_rank==0) print*, "#Total Number of particles: ", sim%npart_loc*sim%coll_size



!--------------------------------------
call sim%init_particle_masks()


sim%SOLVER%dimx=sim%dimx
!Load some particles

call sim%SOLVER%set_box_lens(sim%boxlen)
call sim%SOLVER%init(sim%num_modes)

if  (sim%coll_rank==0) call sim%SOLVER%visu_info()

call sim%set_symplectic_rungekutta_coeffs(sim%time_integrator_order)
call sim%init_diagnostics()

call sim%init_particle()
call sim%init_particle_prior_maxwellian()


!Load particles
SELECT CASE ( sim%testcase%id)
 CASE(SLL_LANDAU_SUM%id)
   call sim%load_landau_sum()
 CASE(SLL_LANDAU_PROD%id)
   call sim%load_landau_prod()
 CASE(SLL_LANDAU_DIAG%id)
   call sim%load_landau_diag()
 CASE DEFAULT
END SELECT

if (sim%coll_rank==0) print *, "TESTCASE: ", trim(sim%testcase%name())

SLL_ALLOCATE(sim%weight_const(sim%npart_loc),ierr)
sim%weight_const=sim%particle(sim%maskw,:)

!Match some moments
if (sim%momentmatch==SLL_MOMENT_MATCH_INITIAL) then
do idx=1, size(sim%maskv)
 call match_moment_1D_weight_linear_real64(sim%particle(sim%maskv(idx),:), &
                sim%particle(sim%maskw,:),0.0_f64, 1.0_f64,sim%npart)
end do
endif
  
!Initial field solve
SLL_ALLOCATE(sim%rhs(sim%SOLVER%problemsize()),ierr)  
SLL_ALLOCATE(sim%solution(sim%SOLVER%problemsize()),ierr)

call sll_collective_barrier(sll_world_collective)
if (sim%coll_rank==0) print *, "# TIME        |IMPULSE ERR.(abs.) | ENERGY ERROR(rel.) | FIELDENERGY | MOMENTUM"


do tstep=1,sim%tsteps
 sim%tstep=tstep
  
  SELECT CASE (sim%dimx)
   CASE(3)
!  call sim%symplectic_rungekutta()
   call sim%YSun_g2h()
   CASE DEFAULT
  call sim%symplectic_rungekutta()
  end SELECT
  
        if (sim%coll_rank==0) write(*,'(A2, G10.4,  G20.8,  G20.8, G20.8, G20.8)') '#', sim%dt*(sim%tstep-1), &
                           abs(sim%moment(1,sim%tstep))/sim%npart, sim%energy_error(sim%tstep), &
                           sim%fieldenergy(sim%tstep), sqrt(sum(sim%moment(:,tstep)**2))
      !if ( (gnuplot_inline_output.eqv. .true.) .AND. coll_rank==0 .AND. mod(timestep-1,timesteps/100)==0  ) then
!                    call energies_electrostatic_gnuplot_inline(kineticenergy(1:tstep), fieldenergy(1:tstep),&
!   			  moment_error(1:tstep),dt)
       !     else

        !    endif
end do

call sim%write_result()


end subroutine run_generalvp_pif

subroutine load_landau_sum_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim 
  sim%particle(sim%maskw,:)=(1.0+sim%eps*sum(cos(diag_dot_matrix_real64(sim%kmode,sim%particle(sim%maskx,:))),1))
  sim%particle(sim%maskw,:)=sim%particle(sim%maskw,:)/sim%prior_weight(:)  
  print *, "Load landau sum"
end subroutine

subroutine load_landau_prod_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
  sll_int32 :: idx
  do idx=1, sim%npart_loc
    sim%particle(sim%maskw,idx)=1+sim%eps*product(cos(sim%kmode(:)*sim%particle(sim%maskx,idx)),1)  
  end do
    sim%particle(sim%maskw,:)=sim%particle(sim%maskw,:)/sim%prior_weight  
end subroutine

subroutine load_landau_diag_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
!     sim%particle(sim%maskw,:)=(1+sim%eps*product(cos(0.5*sim%particle(sim%maskx,:))))
!     sim%particle(sim%maskw,:)=sim%particle(sim%maskw,:)/sim%prior_weight  
end subroutine


 subroutine init_file_generalvp_pif( sim, filename )
  intrinsic :: trim
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
  character(len=*), intent(in)                                :: filename
  sll_real64 :: dt
 sll_int32 :: NUM_TIMESTEPS, NUM_MODES, NUM_PARTICLES, DIMENSION, TIME_INTEGRATOR_ORDER
 sll_real64 :: QoverM, EPSILON
 sll_real64, dimension(10) :: K=0,L=0, B0=0
 sll_int32 :: CONTROLVARIATE, RND_OFFSET
 character(len=32) :: TESTCASE
     
     sll_int32, parameter  :: input_file = 99
     sll_int32             :: IO_stat,ierr,idx
 
     namelist /sim_params/ dt, NUM_PARTICLES, DIMENSION, NUM_TIMESTEPS, &
                          NUM_MODES, QoverM, EPSILON, CONTROLVARIATE, &
                          TIME_INTEGRATOR_ORDER,RND_OFFSET,K,L,B0,&
                          TESTCASE
     
     open(unit = input_file, file=trim(filename),IOStat=IO_stat)
     if( IO_stat /= 0 ) then
        print *, 'init_file_generalvp_pif() failed to open file ', filename
        STOP
     end if
     read(input_file, sim_params)
     close(input_file)

     !Fill simulation with parameters
     sim%dt = dt     
     sim%npart=NUM_PARTICLES
     sim%tsteps=NUM_TIMESTEPS
     sim%dimx=DIMENSION
     sim%qm=QoverM !electrons
     sim%controlvariate=CONTROLVARIATE
     sim%eps=EPSILON
     sim%time_integrator_order=TIME_INTEGRATOR_ORDER
     
     sim%num_modes=NUM_MODES
     sim%RND_OFFSET=RND_OFFSET

     call sim%testcase%parse(TESTCASE)
     
     SLL_ALLOCATE(sim%kmode(sim%dimx),ierr)
     SLL_ALLOCATE(sim%boxlen(sim%dimx),ierr)
 
     sim%kmode=K(1:sim%dimx)
     sim%boxlen=L(1:sim%dimx)
     do idx=1,sim%dimx
       if (sim%boxlen(idx)==0 .and. K(idx)/=0) then
        sim%boxlen=2*sll_pi/K(1:sim%dimx)
       endif
     end do
     
     if (sum(abs(B0(1:sim%dimx)))/=0) then
          SLL_ALLOCATE(sim%B0(sim%dimx),ierr)
         sim%B0=B0(1:sim%dimx)
         print *, "Constant B-Field set to", sim%B0
     endif
     
     
     
!      elseif (size(K)==0) then
!         print *, "k-vector not given"   
!         sim%kmode=0.5_f64
!        else
!        print *, "k-vector has wrong size"   
!         sim%kmode=1
!         sim%kmode(1:size(K))=K   
!      endif
!      
!      print *, K
     !      print *, sim%boxlen
     
  end subroutine init_file_generalvp_pif

  
!Fills particle vector with random numbers, of users choice
subroutine init_particle_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim

sll_int32 :: idx, ierr
sll_int64 :: seed

SLL_ALLOCATE(sim%particle(2*sim%dimx+1,sim%npart_loc),ierr)

!Generate random numbers
seed=10+ sim%RND_OFFSET + sim%coll_rank*sim%npart_loc
do idx=1,sim%npart_loc
!call i8_sobol_generate ( int(dimx,8) , npart, RND_OFFSET , particle(1:2*dimx,:))            
call i8_sobol( int(2*sim%dimx,8), seed, sim%particle(1:2*sim%dimx,idx))
end do
end subroutine init_particle_generalvp_pif

!allocates space
subroutine init_diagnostics_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
sll_int32 :: ierr
  SLL_CLEAR_ALLOCATE(sim%kineticenergy(1:sim%tsteps),ierr)
SLL_CLEAR_ALLOCATE(sim%fieldenergy(1:sim%tsteps),ierr)
SLL_CLEAR_ALLOCATE(sim%energy(1:sim%tsteps),ierr)
SLL_CLEAR_ALLOCATE(sim%energy_error(1:sim%tsteps),ierr)
SLL_CLEAR_ALLOCATE(sim%moment(1:sim%dimx,1:sim%tsteps),ierr)
SLL_CLEAR_ALLOCATE(sim%moment_error(1:sim%tsteps),ierr)
SLL_CLEAR_ALLOCATE(sim%moment_var(1:sim%dimx,1:sim%tsteps),ierr)
SLL_CLEAR_ALLOCATE(sim%weight_sum(1:sim%tsteps),ierr)
SLL_CLEAR_ALLOCATE(sim%weight_var(1:sim%tsteps),ierr)
SLL_CLEAR_ALLOCATE(sim%l2potential(1:sim%tsteps),ierr)

end subroutine init_diagnostics_generalvp_pif
 
!reads tstep
subroutine calculate_diagnostics_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim

sll_int32 :: idx

sim%kineticenergy(sim%tstep)=sum(sum(sim%particle(sim%maskv,:)**2,1)*sim%particle(sim%maskw,:))/sim%npart
call sll_collective_globalsum(sll_world_collective, sim%kineticenergy(sim%tstep))

sim%fieldenergy(sim%tstep)=abs(dot_product(sim%solution,sim%rhs))
sim%l2potential(sim%tstep)=sim%SOLVER%l2norm(sim%solution)

! fieldenergy(tstep)=abs(Efield)
sim%energy(sim%tstep)=sim%kineticenergy(sim%tstep)+sim%fieldenergy(sim%tstep)
sim%energy_error(sim%tstep)=abs(sim%energy(2)-sim%energy(sim%tstep))/abs(sim%energy(1))

do idx=1,size(sim%particle,2)
 sim%moment(:,sim%tstep)=sim%moment(:,sim%tstep)+(sim%particle(sim%maskv,idx)*sim%particle(sim%maskw,idx))
end do
! normalize
sim%moment(:,sim%tstep)=sim%moment(:,sim%tstep)/sim%npart
call sll_collective_globalsum(sll_world_collective, sim%moment(:,sim%tstep))

sim%moment_error(sim%tstep)=sqrt(sum((sim%moment(:,1)-sim%moment(:,sim%tstep))**2))

sim%weight_sum(sim%tstep)=sum(sim%particle(sim%maskw,:))/sim%npart
call sll_collective_globalsum(sll_world_collective, sim%weight_sum(sim%tstep))

sim%weight_var(sim%tstep)=sum( (sim%particle(sim%maskw,:)-sim%weight_sum(sim%tstep))**2)/sim%npart
call sll_collective_globalsum(sll_world_collective, sim%weight_var(sim%tstep))

end subroutine calculate_diagnostics_generalvp_pif



subroutine init_particle_prior_maxwellian_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim

sll_int32 :: idx,jdx,ierr

!scale to boxlength
do idx=1,sim%npart_loc
sim%particle(sim%maskx,idx)=sim%boxlen*sim%particle(sim%maskx,idx)
end do

!load gaussian profile
do idx=1,size(sim%particle,2)
 do jdx=1, size(sim%maskv)
  call normal_cdf_inv( sim%particle(sim%maskv(jdx),idx), 0.0_f64 , 1.0_f64, sim%particle(sim%maskv(jdx),idx))
 end do
end do

!set temperature and impulse
 do idx=1, size(sim%maskv)
! print *, sum(sim%particle(sim%maskv(idx),:))/sim%npart
 !call match_moment_1D_linear_real64(sim%particle(sim%maskv(idx),:), 0.0_f64, 1.0_f64)
! print *, sum(sim%particle(sim%maskv(idx),:))/sim%npart
  end do

 !print *,minval(sim%particle(sim%maskx,:)), maxval(sim%particle(sim%maskx,:))
 !  print *,minval(sim%particle(sim%maskv,:)), maxval(sim%particle(sim%maskv,:))
  
SLL_ALLOCATE(sim%prior_weight(1:sim%npart_loc),ierr)
sim%prior_weight=1.0/product(sim%boxlen)

end subroutine init_particle_prior_maxwellian_generalvp_pif


!>Only for separable lagrangian, no magnetic field
subroutine symplectic_rungekutta_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim

  sll_int32 :: rkidx
  sll_real64 :: t !time
  SLL_ASSERT(size(sim%rk_d)==size(sim%rk_c))
  
    
  !Loop over all stages
   t=(sim%tstep-1)*sim%dt
    do rkidx=1, size(sim%rk_d);
     
     !Set control variate if used
     call sim%update_weight()
     
     !Charge assignement, get right hand side
     sim%rhs=sim%SOLVER%get_rhs_particle(sim%particle(sim%maskxw,:))/sim%npart
     !mpi allreduce
     call sll_collective_globalsum(sll_world_collective, sim%rhs)

     !sim%solution=sim%SOLVER%solve_poisson(sim%rhs)
     sim%solution=sim%SOLVER%solve_quasineutral(sim%rhs)
     
     if (rkidx==1) then
      call sim%calculate_diagnostics()
     endif

     sim%particle(sim%maskv,:)=sim%particle(sim%maskv,:) -  &
          sim%rk_c(rkidx)*sim%dt*sim%qm*&
         (sim%E(sim%particle(sim%maskx,:),t)     + &                          ! Electric field external
         (-sim%SOLVER%eval_gradient(sim%particle(sim%maskx,:),sim%solution)) ); !Electric field selfconsistent
     sim%particle(sim%maskx,:)=sim%particle(sim%maskx,:) +  sim%rk_d(rkidx)*sim%dt*sim%particle(sim%maskv,:);
        
     ! xx=mod(xx,repmat(L,1,npart));
         t=(sim%tstep-1)*sim%dt+sum(sim%rk_d(1:rkidx))*sim%dt;
    end do
     t=(sim%tstep)*sim%dt;

end subroutine symplectic_rungekutta_generalvp_pif


subroutine YSun_g2h_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
  sll_real64, dimension(sim%dimx, sim%npart_loc) :: E,B
  sll_int32 :: rkidx
  sll_real64 :: t, h !time
  sll_real64, parameter :: gamma1=1.0_f64/(2.0_f64-2.0_f64**(1.0_f64/3.0_f64))
  sll_real64, parameter :: gamma0=1.0_f64-2.0_f64/(2.0_f64-2.0_f64**(1.0_f64/3.0_f64))
  sll_real64, dimension(3) :: gamma
     
    gamma(1)=gamma1
    gamma(2)=gamma0
    gamma(3)=gamma1
         
  !Loop over all stages
    t=(sim%tstep-1)*sim%dt
     do rkidx=1, size(gamma);
     !Set control variate if used
     call sim%update_weight()
     
     !Charge assignement, get right hand side
     sim%rhs=sim%SOLVER%get_rhs_particle(sim%particle(sim%maskxw,:))/sim%npart
     !mpi allreduce
     call sll_collective_globalsum(sll_world_collective, sim%rhs)
     !sim%solution=sim%SOLVER%solve_poisson(sim%rhs)
     sim%solution=sim%SOLVER%solve_quasineutral(sim%rhs)
    
       if (rkidx==1) then
         call sim%calculate_diagnostics()
       endif

      h=gamma(rkidx)*sim%dt
    
      E= sim%E(sim%particle(sim%maskx,:),t) + &  ! Electric field external
           (sim%SOLVER%eval_gradient(sim%particle(sim%maskx,:),sim%solution)) !Electric field selfconsistent
          
      B=sim%B(sim%particle(sim%maskx,:), t) !External magnetic field
      
      
         sim%particle(sim%maskv,:)=exp_skew_product2( B, &
                         sim%particle(sim%maskv,:) + h*sim%qm/2*E, h*(-sim%qm)*l2norm(B) ) + h*sim%qm/2*E

!        sim%particle(sim%maskv,:)=sim%particle(sim%maskv,:)- h*sim%qm*E
       
!        sim%particle(sim%maskv,:)=exp_skew_product(h*sim%qm*sim%B(sim%particle(sim%maskx,:), t) , &
!                        sim%particle(sim%maskv,:) + h*sim%qm/2*E) + h*sim%qm/2*E
     !sim%particle(sim%maskv,:)=sim%particle(sim%maskv,:) +  h*(sim%qm)*(sim%vxB(sim%particle(sim%maskx,:),sim%particle(sim%maskv,:),t)+  E);
    
    !x_{k+1}= x_k + dt*v_{k+1}
    sim%particle(sim%maskx,:)=sim%particle(sim%maskx,:) +  h*sim%particle(sim%maskv,:);
      
          t=(sim%tstep-1)*sim%dt+sum(gamma(1:rkidx))*sim%dt;
     end do
     t=(sim%tstep)*sim%dt;
     
end subroutine YSun_g2h_generalvp_pif
 
 function l2norm(x) result(norm)
   sll_real64, dimension(:,:), intent(in) :: x
   sll_real64, dimension(size(x,2)) :: norm
 norm=sqrt(sum(x**2,1))
 end function
 
! expm(omega*(v \times) ) w
  function exp_skew_product2( v, w , omega) result(c)
  sll_real64, dimension(:,:), intent(in) :: v, w
  sll_real64, dimension(3,size(v,2)) :: c
  sll_real64, dimension(:), intent(in) :: omega
     
   c(1,:) = w(3,:)*(v(3,:)*sin(omega) + (v(1,:)*v(3,:)*sin(omega/2)**2)/8) - w(2,:)*(v(3,:)*sin(omega) - (v(1,:)*v(3,:)*sin(omega/2)**2)/8) - w(1,:)*((sin(omega/2)**2*(v(3,:)**2 + v(3,:)**2))/8 - 1)
   c(2,:)=w(1,:)*(v(3,:)*sin(omega) + (v(1,:)*v(3,:)*sin(omega/2)**2)/8) - w(2,:)*((sin(omega/2)**2*(v(1,:)**2 + v(3,:)**2))/8 - 1) - w(3,:)*(v(1,:)*sin(omega) - (v(3,:)*v(3,:)*sin(omega/2)**2)/8)
   c(3,:)=w(2,:)*(v(1,:)*sin(omega) + (v(3,:)*v(3,:)*sin(omega/2)**2)/8) - w(1,:)*(v(3,:)*sin(omega) - (v(1,:)*v(3,:)*sin(omega/2)**2)/8) - w(3,:)*((sin(omega/2)**2*(v(1,:)**2 + v(3,:)**2))/8 - 1) 
     
!         c(1,:) = -w(1,:)*(sin(omega*(1.0D0/2.0D0))**2*(1.0D0/8.0D0)-1.0D0)-w(2,:)*sin(omega)
!         c(2,:) = -w(2,:)*(sin(omega*(1.0D0/2.0D0))**2*(1.0D0/8.0D0)-1.0D0)+w(1,:)*sin(omega)
!         c(3,:) = w(3,:)
!       
!      
!      c(1,:) = w(1,:)*cos(omega)-w(2,:)*sin(omega)
!       c(2,:) = w(2,:)*cos(omega)+w(1,:)*sin(omega)
!       c(3,:) = w(3,:)
!   
  end function


  function exp_skew_product( v, w) result(c)
  sll_real64, dimension(:,:), intent(in) :: v, w
  sll_real64, dimension(3,size(v,2)) :: c
  sll_comp64, dimension(size(v,2)) :: s
  sll_real64, dimension(size(v,2)) :: vv
  sll_comp64, dimension(size(v,2)) :: sqrtvv
  
  
  vv=-(v(1,:)**2+v(2,:)**2+v(3,:)**2)
  sqrtvv=sll_i1*sqrt(-vv)
  
     
     c(1,:) = real((w(1,:)*exp(-sqrtvv)*(v(1,:)**2*exp(sqrtvv)*2.0D0+v(3,:)**2*exp(sqrtvv*2.0D0)+v(3,:)**2*exp(sqrtvv*2.0D0)+v(3,:)**2+v(3,:)**2)*&
     (-1.0D0/2.0D0))/vv+1.0D0/sqrtvv**(3.0D0)*w(2,:)*exp(-sqrtvv)*(exp(sqrtvv)-1.0D0)*&
     (v(3,:)**3*exp(sqrtvv)+v(1,:)**2*v(3,:)+v(3,:)**2*v(3,:)+v(3,:)**3+v(1,:)**2*v(3,:)*&
     exp(sqrtvv)+v(3,:)**2*v(3,:)*exp(sqrtvv)-sqrtvv*v(1,:)*v(3,:)+sqrtvv*v(1,:)*v(3,:)*exp(sqrtvv))*(1.0D0/2.0D0)+(w(3,:)*&
     exp(-sqrtvv)*(exp(sqrtvv)-1.0D0)*(sqrtvv*v(3,:)-v(1,:)*v(3,:)+sqrtvv*v(3,:)*exp(sqrtvv)+v(1,:)*v(3,:)*exp(sqrtvv))*(1.0D0/2.0D0))/vv)
  
     c(2,:) =real( (w(2,:)*exp(-sqrtvv)*(v(3,:)**2*exp(sqrtvv)*2.0D0+v(1,:)**2*exp(sqrtvv*2.0D0)+v(3,:)**2*exp(sqrtvv*2.0D0)+v(1,:)**2+v(3,:)**2)*&
     (-1.0D0/2.0D0))/vv-1.0D0/sqrtvv**(3.0D0)*w(1,:)*exp(-sqrtvv)*(exp(sqrtvv)-1.0D0)*(v(3,:)**3*&
     exp(sqrtvv)+v(1,:)**2*v(3,:)+v(3,:)**2*v(3,:)+v(3,:)**3+v(1,:)**2*v(3,:)*exp(sqrtvv)+v(3,:)**2*v(3,:)*&
     exp(sqrtvv)+sqrtvv*v(1,:)*v(3,:)-sqrtvv*v(1,:)*v(3,:)*exp(sqrtvv))*(1.0D0/2.0D0)-&
     (w(3,:)*exp(-sqrtvv)*(exp(sqrtvv)-1.0D0)*(sqrtvv*v(1,:)+v(3,:)*v(3,:)+sqrtvv*v(1,:)*&
     exp(sqrtvv)-v(3,:)*v(3,:)*exp(sqrtvv))*(1.0D0/2.0D0))/vv)
     
     
     c(3,:) = real((w(3,:)*exp(-sqrtvv)*(v(3,:)**2*exp(sqrtvv)*2.0D0+v(1,:)**2*exp(sqrtvv*2.0D0)+v(3,:)**2*exp(sqrtvv*2.0D0)+v(1,:)**2+v(3,:)**2)*(-1.0D0/2.0D0))/&
     vv-1.0D0/sqrtvv**(3.0D0)*w(2,:)*exp(-sqrtvv)*(exp(sqrtvv)-1.0D0)*(v(1,:)**3*exp(sqrtvv)+v(1,:)*v(3,:)**2+v(1,:)*v(3,:)**2+v(1,:)**3+v(1,:)&
     *v(3,:)**2*exp(sqrtvv)+v(1,:)*v(3,:)**2*exp(sqrtvv)+sqrtvv*v(3,:)*v(3,:)-sqrtvv*v(3,:)*v(3,:)*exp(sqrtvv))*(1.0D0/2.0D0)+&
     1.0D0/sqrtvv**(3.0D0)*w(1,:)*exp(-sqrtvv)*(exp(sqrtvv)-1.0D0)*(v(3,:)**3*exp(sqrtvv)+v(1,:)**2*v(3,:)+v(3,:)*v(3,:)**2+&
     v(3,:)**3+v(1,:)**2*v(3,:)*exp(sqrtvv)+v(3,:)*v(3,:)**2*exp(sqrtvv)-sqrtvv*v(1,:)*v(3,:)+sqrtvv*v(1,:)*v(3,:)*exp(sqrtvv))*(1.0D0/2.0D0))
  
  
!    c(1,:) = w(1,:)*5.403023058681397D-1-w(2,:)*8.414709848078966D-1
!    c(2,:) = w(1,:)*8.414709848078965D-1+w(2,:)*5.403023058681398D-1
!    c(3,:) = w(3,:)
!     c=0
!      s(:)=sll_i1*sqrt(v(1,:)**2+v(2,:)**2+v(3,:)**2)  
!    c(1,:) =real( (w(1,:)*exp(-s)*(v(1,:)**2*exp(s)*2.0D0+v(2,:)**2*exp(s*2.0D0)+v(3,:)**2*exp(s*2.0D0)+v(2,:)**2+v(3,:)**2)*&
!            (1.0D0/2.0D0))/(v(1,:)**2+v(2,:)**2+v(3,:)**2)+w(2,:)*exp(s)*(exp(s)-1.0D0)*1.0D0/ &
!            sll_i1**3/sqrt(v(1,:)**2+v(2,:)**2+ v(3,:)**2)**(3.0D0)*(v(3,:)**3*exp(s)+v(1,:)**2*v(3,:)+v(2,:)**2*v(3,:)+v(3,:)**3+v(1,:)**2*v(3,:)*exp(s)+&
!            v(2,:)**2*v(3,:)*exp(s)-s*v(1,:)*v(2,:)+s*v(1,:)*v(2,:)*exp(s))*(1.0D0/2.0D0)-(w(3,:)*exp(s)*(exp(s)-1.0D0)*(s*v(2,:)&
!            -v(1,:)*v(3,:)+s*v(2,:)*exp(s)+v(1,:)*v(3,:)*exp(s))*(1.0D0/2.0D0))/(v(1,:)**2+v(2,:)**2+v(3,:)**2))

!   c(2,:) =real( (w(2,:)*exp(-sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*(v(1,:)**2*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2)*2.0D0)+v(3,:)**2*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2)*2.0D0)+v(1,:)**2+v(3,:)**2+v(2,:)**2*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*2.0D0)*(1.0D0/2.0D0))&
!   /(v(1,:)**2+v(2,:)**2+v(3,:)**2)+(w(3,:)*exp(-sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*&
!   (exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))-1.0D0)*(v(1,:)*sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2)+v(2,:)*v(3,:)-v(2,:)*&
!      v(3,:)*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))+v(1,:)*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*&
!      sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*(1.0D0/2.0D0))/(v(1,:)**2+v(2,:)**2+v(3,:)**2)-w(1,:)*&
!      exp(-sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*(exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))-1.0D0)&
!      *1.0D0/(-v(1,:)**2-v(2,:)**2-v(3,:)**2)**(3.0D0/2.0D0)*(v(1,:)**2*v(3,:)+v(2,:)**2*v(3,:)+v(3,:)**3&
!      +v(3,:)**3*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))+v(1,:)**2*v(3,:)*exp(sqrt(-v(1,:)**2-v(2,:)**2&
!      -v(3,:)**2))+v(2,:)**2*v(3,:)*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))+v(1,:)*v(2,:)*sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2)-&
!      v(1,:)*v(2,:)*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*(1.0D0/2.0D0))
!           
!  c(3,:) = real((w(3,:)*exp(-sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*(v(1,:)**2*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2)*2.0D0)+&
!  v(2,:)**2*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2)*2.0D0)+v(1,:)**2+&
!  v(2,:)**2+v(3,:)**2*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*2.0D0)*(1.0D0/2.0D0))&
!  /(v(1,:)**2+v(2,:)**2+v(3,:)**2)-w(2,:)*exp(-sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*&
!  (exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))-1.0D0)*1.0D0/(-v(1,:)**2-v(2,:)**2-v(3,:)**2)**(3.0D0/2.0D0)*(v(1,:)*v(2,:)**2+&
!  v(1,:)*v(3,:)**2+v(1,:)**3+v(1,:)**3*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))+&
!  v(1,:)*v(2,:)**2*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))+v(1,:)*v(3,:)**2*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))+v(2,:)*v(3,:)*&
!  sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2)-v(2,:)*v(3,:)*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*&
!  sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*(1.0D0/2.0D0)+w(1,:)*exp(-sqrt(-v(1,:)**2-v(2,:)**2-&
!  v(3,:)**2))*(exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))-1.0D0)*1.0D0/(-v(1,:)**2-v(2,:)**2-v(3,:)**2)**(3.0D0/2.0D0)*(v(1,:)**2*v(2,:)+v(2,:)*v(3,:)**2+v(2,:)**3+v(2,:)**3*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))+v(1,:)**2*v(2,:)*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))+v(2,:)*v(3,:)**2*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))-v(1,:)*v(3,:)*sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2)+v(1,:)*v(3,:)*exp(sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*sqrt(-v(1,:)**2-v(2,:)**2-v(3,:)**2))*(1.0D0/2.0D0))
 
!  print *, sum(c(1,:))
 
 
 
  end function
 
 
subroutine set_symplectic_rungekutta_coeffs_generalvp_pif(sim,rk_order)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
sll_int32, intent(in) :: rk_order
sll_int32 :: ierr
sll_real64,parameter :: rk4sx=((2**(1/3) +2**(-1/3)-1)/6)

SLL_ALLOCATE(sim%rk_c(rk_order),ierr)
SLL_ALLOCATE(sim%rk_d(rk_order),ierr)

SELECT CASE (rk_order)
   CASE (1)
    !euler not symplectic
    sim%rk_d=1
    sim%rk_c=1
   CASE (2)
      sim%rk_d=(/0.5, 0.5 /)
      sim%rk_c=(/0.0, 1.0/) 
   CASE (3)
     sim%rk_d(:)=(/ 2.0/3.0, -2.0/3.0, 1.0 /) 
     sim%rk_c(:)=(/ 7.0/24.0,3/4.0,-1.0/24.0 /)  
   CASE (4)
      sim%rk_d=(/ 2.0*rk4sx+1.0 , -4.0*rk4sx-1.0, 2.0*rk4sx+1.0, 0.0_f64/) 
      sim%rk_c=(/ rk4sx + 0.5 , -rk4sx, -rk4sx, rk4sx +0.5 /) 
END SELECT

print *, "Symplectic Runge Kutta of order:", rk_order
print *, "D_COEFFS:", sim%rk_d
print *, "C_COEFFS:", sim%rk_c

end subroutine set_symplectic_rungekutta_coeffs_generalvp_pif


subroutine init_particle_masks_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
 sll_int32 :: idx,ierr

 !allocate stencils
SLL_ALLOCATE(sim%maskx(sim%dimx),ierr)
SLL_ALLOCATE(sim%maskxw(sim%dimx+1),ierr)
SLL_ALLOCATE(sim%maskv(sim%dimx),ierr)
SLL_ALLOCATE(sim%maskxv(sim%dimx*2),ierr)

sim%maskx=(/( idx,idx=1,sim%dimx)/)
sim%maskv=(/( idx,idx=sim%dimx+1,2*sim%dimx)/)
sim%maskw=2*sim%dimx+1
sim%maskxw(1:sim%dimx)=sim%maskx
sim%maskxw(sim%dimx+1)=sim%maskw
sim%maskxv(1:sim%dimx)=sim%maskx
sim%maskxv(sim%dimx+1:2*sim%dimx)=sim%maskv

end subroutine init_particle_masks_generalvp_pif


!Updates the weights for use with control variate
subroutine update_weight_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim

if (sim%controlvariate /= 0) then

  if (sim%collisions/=0) then
!particle(maskw,:)=initial_prior -  (particle(maskw,:)-initial   control_variate(particle(maskxv,:))
print *, "not implemented"
else
!no collisions, characteristics are conserved
 sim%particle(sim%maskw,:)=sim%weight_const(:) - &
            sim%control_variate(sim%particle(sim%maskxv,:))/sim%prior_weight(:)
 endif

 endif
end subroutine update_weight_generalvp_pif



function control_variate_generalvp_pif(sim, particle) result(cv)
  class(sll_simulation_general_vlasov_poisson_pif), intent(in) :: sim
sll_real64, dimension(:,:), intent(in) :: particle !(x,v)
sll_real64, dimension(size(particle,2)) :: cv

!Standard maxwellian control variate
SELECT CASE (sim%controlvariate)
   CASE (SLL_CONTROLVARIATE_NONE)
     !This should not happen
      cv=1
    CASE (SLL_CONTROLVARIATE_STANDARD)
#pragma vector always
      cv=sqrt(2*sll_pi)**(-size(sim%maskv))*exp(-0.5*sum(sim%particle(sim%maskv,:),1)**2)
    CASE (SLL_CONTROLVARIATE_MAXWELLIAN)
      cv=sqrt(2*sll_pi)**(-size(sim%maskv))*exp(-0.5*sum(sim%particle(sim%maskv,:),1)**2)
END SELECT

end function control_variate_generalvp_pif



subroutine write_result_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(in) :: sim       
        integer :: idx,file_id,file_id_err,ierr


        if (sim%coll_rank==0) then

             call sll_new_file_id(file_id, ierr)

            !Write Data File
            !            open(file_id, file = plot_name//"_"//fin//'.dat' )
             open(file_id, file = 'pif_result.csv')
            !             open(file_id, file = './'//filename//'.csv')
!             write (file_id, *)  "#Full 1d1v Electrostatic PIC"
!             write (file_id,*)  "#Time steps:", timesteps
!             write (file_id,*)  "#Time stepwidth:", timestepwidth
!             write (file_id,*)  "#Marko particles:", coll_size*nparticles
!             write (file_id,*)  "#Particle Pusher:", particle_pusher
!             write (file_id,*)  "#Finite Elements: 2^(", log(real(mesh_cells,i64))/log(2.0_f64),")"
!             write (file_id,*)  "#Size of MPI Collective: ", coll_size
            write (file_id,*)  ' \"time\", \"kineticenergy\", \"fieldenergy\", \"energy_error\", \"moment_error\", \"l2potential\"'
            do idx=1,sim%tsteps
                write (file_id,*) (idx-1)*sim%dt,",",sim%kineticenergy(idx),",",sim%fieldenergy(idx),",", &
                        sim%fieldenergy(idx),",", sim%energy_error(idx),",", sim%moment_error(idx),",",&
                        sim%l2potential(idx)
            enddo
            close(file_id)
        
          open(file_id, file = 'pif_result.gnu')
       
          write(file_id,*) "set term x11 1"
          write(file_id,*) "set logscale y"
          write(file_id,*) "plot '-' using 1:2 title 'Absolute Impulse Error' with lines"
          do idx = 1, sim%tsteps
            write(file_id,*)  (idx-1)*sim%dt,  sim%moment_error(idx)
          end do
          write(file_id,*) "e"
          write(file_id,*) "   "
          
          
          write(file_id,*) "set term x11 2"
          write(file_id,*) "set logscale y"
          write(file_id,*) "plot '-' using 1:2 title 'Relative Energy Error' with lines"
          do idx = 1, sim%tsteps
            write(file_id,*)  (idx-1)*sim%dt,  sim%energy_error(idx)
          end do
          write(file_id,*) "e"
          write(file_id,*) "   "
       
          write(file_id,*) "set term x11 3"
          write(file_id,*) "set logscale y"
          write(file_id,*) "plot '-' using 1:2 title 'Electrostatic Energy' with lines"
          do idx = 1, sim%tsteps
            write(file_id,*)  (idx-1)*sim%dt,  sim%fieldenergy(idx)
          end do
          write(file_id,*) "e"
          write(file_id,*) "   "
       
          write(file_id,*) "set term x11 4"
          write(file_id,*) "set logscale y"
          write(file_id,*) "plot '-' using 1:2 title 'Kinetic Energy' with lines"
          do idx = 1, sim%tsteps
            write(file_id,*)  (idx-1)*sim%dt,  sim%kineticenergy(idx)
          end do
          write(file_id,*) "e"
          write(file_id,*) "   "
          
          write(file_id,*) "set term x11 5"
          write(file_id,*) "set logscale y"
          write(file_id,*) "plot '-' using 1:2 title 'Sum of weights' with lines"
          do idx = 1, sim%tsteps
            write(file_id,*)  (idx-1)*sim%dt,  sim%weight_sum(idx)
          end do
          write(file_id,*) "e"
          write(file_id,*) "   "
          close(file_id)
        
          write(file_id,*) "set term x11 6"
          write(file_id,*) "set logscale y"
          write(file_id,*) "plot '-' using 1:2 title 'Electrostatic Potential' with lines"
          do idx = 1, sim%tsteps
            write(file_id,*)  (idx-1)*sim%dt, sim%l2potential(idx)
          end do
          write(file_id,*) "e"
          write(file_id,*) "   "
        
        
        endif
end subroutine write_result_generalvp_pif





end module sll_general_vlasov_poisson_pif