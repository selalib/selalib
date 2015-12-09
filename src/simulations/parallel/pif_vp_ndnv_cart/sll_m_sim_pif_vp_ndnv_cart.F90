!*************************************************************
!  Author: Jakob Ameres, jakob.ameres@tum.de
!**************************************************************

module sll_m_sim_pif_vp_ndnv_cart
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use hdf5, only: &
    hid_t

  use sll_m_collective, only: &
    sll_collective_barrier, &
    sll_collective_globalsum, &
    sll_get_collective_rank, &
    sll_get_collective_size, &
    sll_world_collective

  use sll_m_constants, only: &
    sll_i1, &
    sll_pi

  use sll_m_descriptors, only: &
    sll_landau_diag, &
    sll_landau_prod, &
    sll_landau_sum, &
    sll_vlasovpoisson_sim

  use sll_m_hdf5_io_serial, only: &
    sll_hdf5_file_close, &
    sll_hdf5_file_create, &
    sll_hdf5_write_array

  use sll_m_moment_matching, only: &
    match_moment_1d_weight_linear_real64

  use sll_m_particle_method_descriptors, only: &
    sll_collisions_none, &
    sll_controlvariate_maxwellian, &
    sll_controlvariate_none, &
    sll_controlvariate_standard, &
    sll_moment_match_initial, &
    sll_moment_match_none

  use sll_m_pic_visu, only: &
    plot_format_points3d

  use sll_m_pic_visu_parallel, only: &
    distribution_xdmf_coll

  use sll_m_pif_fieldsolver, only: &
    diag_dot_matrix_real64, &
    pif_fieldsolver

  use sll_m_prob, only: &
    normal_cdf_inv

  use sll_m_sim_base, only: &
    sll_simulation_base_class

  use sll_m_sobol, only: &
    i8_sobol

  use sll_m_time_composition, only: &
    comp_coeff_sym_sym

  use sll_m_utilities, only: &
    int2string, &
    sll_new_file_id

  use sll_m_wedge_product_generaldim, only: &
    cross_product_2d, &
    cross_product_3d

  implicit none

  public :: &
    sll_simulation_general_vlasov_poisson_pif

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! abstract interface
!         function electric_field_general(x,t) result(E)
!             use sll_m_working_precision
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
 
 sll_int32, parameter :: PIF_NO_E_FIELD=0
 sll_int32, parameter :: PIF_POISSON=1
 sll_int32, parameter :: PIF_QUASINEUTRAL_RHO=2
 sll_int32, parameter :: PIF_QUASINEUTRAL_RHO_WO_ZONALFLOW=3
 sll_int32, parameter :: PIF_QUASINEUTRAL=4
 sll_int32, parameter :: PIF_QUASINEUTRAL_WO_ZONALFLOW=5

 
 
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
sll_int32 :: FIELDSOLVER=PIF_QUASINEUTRAL_RHO_WO_ZONALFLOW !PIF_POISSON
type(sll_vlasovpoisson_sim) :: testcase=SLL_LANDAU_SUM
sll_int32 :: RND_OFFSET   = 10    !Default sobol offset, skip zeros
sll_int32 :: num_modes = 1

!Plots
sll_int32, dimension(2) :: plot2d_idx ! Contains 2 indicies of the dimensions to be plotted
sll_int32, dimension(2) :: plot2d_bin ! Contains number of bins for each dimension
sll_int32 :: num_visu_particles = 0
character(len=64) :: prefix='pif'
sll_int32 :: WRITE_PHI = 0

!results
sll_real64, dimension(:),allocatable :: kineticenergy,fieldenergy, energy,energy_error,&
               weight_sum,weight_var,moment_error, l2potential
sll_real64, dimension(:,:),allocatable ::moment,moment_var
sll_comp64, dimension(:), allocatable :: rhs, solution
sll_real64 :: ExB=0
!--------------

!------MPI---------
sll_int32 :: coll_rank, coll_size
!------MPI---------

   contains
     
     procedure, pass(sim) :: init_particle => init_particle_generalvp_pif
     procedure, pass(sim) :: init_from_file => init_file_generalvp_pif
     procedure, pass(sim) :: delete=>delete_generalvp_pif
     procedure, pass(sim) :: write_result=>write_result_generalvp_pif
     procedure, pass(sim) :: run=>run_generalvp_pif
     procedure, pass(sim) :: update_weight=>update_weight_generalvp_pif
     procedure, pass(sim) :: init_diagnostics=>init_diagnostics_generalvp_pif
     procedure, pass(sim) :: init_particle_masks=>init_particle_masks_generalvp_pif
     procedure, pass(sim) :: init_particle_prior_maxwellian=>init_particle_prior_maxwellian_generalvp_pif
     procedure, pass(sim) :: set_symplectic_rungekutta_coeffs=>set_symplectic_rungekutta_coeffs_generalvp_pif
     procedure, pass(sim) :: symplectic_rungekutta=>symplectic_rungekutta_generalvp_pif
     procedure, pass(sim) :: YSun_g2h=>YSun_g2h_generalvp_pif
     procedure, pass(sim) :: heun=>heun_generalvp_pif
     procedure, pass(sim) :: calculate_diagnostics=>calculate_diagnostics_generalvp_pif
     procedure, pass(sim) :: control_variate=>control_variate_generalvp_pif
     procedure, pass(sim) :: vxB=>v_cross_B_generalvp_pif
     !Magnetic field, can be changed
     procedure, pass(sim) :: B=>Bstandard_generalvp_pif
     !>Unit vector of B
!      procedure, pass(sim) :: Bunit=
      procedure, pass(sim) :: E=>Ezero_generalvp_pif !Electric field
     
     !wrapper for Field solve
     procedure, pass(sim) :: solve_field=>solve_field_sum_generalvp_pif
     
     !Loading
     procedure, pass(sim) :: load_landau_sum=>load_landau_sum_generalvp_pif
     procedure, pass(sim) :: load_landau_prod=>load_landau_prod_generalvp_pif
     procedure, pass(sim) :: load_landau_diag=>load_landau_diag_generalvp_pif
     
     
     !Visualization
     procedure, pass(sim) :: visu_bound_low=>visu_bound_low
     procedure, pass(sim) :: visu_bound_up=>visu_bound_up
     procedure, pass(sim) :: visu_phasespace=>visu_phasespace
end type

!--------------------------------------------------------------

contains

!>Destructor, disallocates everything
subroutine delete_generalvp_pif(sim)
 class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
 sll_int32 :: ierr
if (allocated(sim%particle)) then
  SLL_DEALLOCATE_ARRAY(sim%particle,ierr)
endif
 
if (allocated(sim%weight_const)) then
  SLL_DEALLOCATE_ARRAY(sim%weight_const,ierr)
endif

if (allocated(sim%prior_weight)) then
  SLL_DEALLOCATE_ARRAY(sim%prior_weight,ierr)
endif

if (allocated(sim%rk_d)) then
  SLL_DEALLOCATE_ARRAY(sim%rk_d,ierr)
endif
if (allocated(sim%rk_c)) then
  SLL_DEALLOCATE_ARRAY(sim%rk_c,ierr)
endif

!type(pif_fieldsolver) :: SOLVER


  SLL_DEALLOCATE_ARRAY(sim%maskx,ierr)
  SLL_DEALLOCATE_ARRAY(sim%maskv,ierr)
  SLL_DEALLOCATE_ARRAY(sim%maskxw,ierr)
  SLL_DEALLOCATE_ARRAY(sim%maskxv,ierr)
  
  SLL_DEALLOCATE_ARRAY(sim%boxlen,ierr)
  SLL_DEALLOCATE_ARRAY(sim%kmode,ierr)

if (allocated(sim%B0)) then
  SLL_DEALLOCATE_ARRAY(sim%B0,ierr)
endif


SLL_DEALLOCATE_ARRAY(sim%kineticenergy,ierr)
SLL_DEALLOCATE_ARRAY(sim%fieldenergy,ierr)
SLL_DEALLOCATE_ARRAY(sim%energy_error,ierr)
SLL_DEALLOCATE_ARRAY(sim%energy,ierr)
SLL_DEALLOCATE_ARRAY(sim%weight_sum,ierr)
SLL_DEALLOCATE_ARRAY(sim%weight_var,ierr)
SLL_DEALLOCATE_ARRAY(sim%moment_error,ierr)
SLL_DEALLOCATE_ARRAY(sim%l2potential,ierr)
SLL_DEALLOCATE_ARRAY(sim%moment,ierr)
SLL_DEALLOCATE_ARRAY(sim%moment_var,ierr)

SLL_DEALLOCATE_ARRAY(sim%rhs,ierr)
SLL_DEALLOCATE_ARRAY(sim%solution,ierr)
end subroutine


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




function solve_field_sum_generalvp_pif(sim, rhs) result(solution)
 class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
 sll_comp64, dimension(:), intent(in) :: rhs
 sll_comp64, dimension(size(rhs)) :: solution
 
 SELECT CASE (sim%FIELDSOLVER)
    CASE(PIF_NO_E_FIELD)
     solution=0
    CASE(PIF_POISSON)
      solution=sim%SOLVER%solve_poisson(rhs)
    CASE(PIF_QUASINEUTRAL_RHO)
      solution=sim%SOLVER%solve_quasineutral(rhs)
    CASE(PIF_QUASINEUTRAL_RHO_WO_ZONALFLOW)
      solution=sim%SOLVER%solve_qn_rho_wo_zonalflow(rhs)
 END SELECT
 
end function solve_field_sum_generalvp_pif


subroutine run_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
  sll_int32 ::ierr , idx, tstep, phi_file_id
  character(len=4) :: timestr
  integer(hid_t) :: hdata_id    !< HDF5 data file for each timestep

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
if (allocated(sim%B0)) then
 if (sim%coll_rank==0) print *, "Constant B-Field set to", sim%B0
endif



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


!open file to write the Electric potential
if (sim%WRITE_PHI==1 .and. sim%coll_rank==0 ) then 
   call sll_new_file_id(phi_file_id, ierr)
   open(phi_file_id, file = trim(sim%prefix)//'_phi_unitmodes.dat')
     do idx=1,size(sim%SOLVER%allmodes,2)
       write (phi_file_id,*) sim%SOLVER%allmodes(:,idx)
     end do
   close(phi_file_id)
endif





do tstep=1,sim%tsteps
 sim%tstep=tstep
  
  SELECT CASE (sim%dimx)
   CASE(3)

   if (allocated(sim%B0)) then
   !  call sim%symplectic_rungekutta()
   !call sim%heun()
   call sim%YSun_g2h()
   else
   call sim%symplectic_rungekutta()
   endif
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
  !periodic Boundary conditions
  do idx=1, sim%npart_loc
   sim%particle(sim%maskx,idx)=modulo(sim%particle(sim%maskx,idx), sim%boxlen)
  end do
  
  
!   if (sim%WRITE_PHI==1 .and. sim%coll_rank==0 ) then 
!     write (phi_file_id, *) abs(sim%solution)
!   endif

  if (sim%coll_rank==0) then
   if (sim%WRITE_PHI==1) then
   call int2string(sim%tstep, timestr)
   call sll_hdf5_file_create(trim(sim%prefix)//'_data_'//timestr//'.h5',hdata_id,ierr)
   call sll_hdf5_write_array(hdata_id,abs(sim%solution),"/phi",ierr)
   call sll_hdf5_write_array(hdata_id,abs(sim%rhs),"/rho",ierr)
   call sll_hdf5_file_close(hdata_id, ierr)
   endif
  endif
  
  call sim%visu_phasespace()
end do


call sim%write_result()




end subroutine run_generalvp_pif

subroutine visu_phasespace(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(in) :: sim 
  
  if  (all(sim%plot2d_idx/=0)) then
    call distribution_xdmf_coll(trim(sim%prefix)//'_plot2d',  sim%particle(sim%plot2d_idx(1),:), sim%particle(sim%plot2d_idx(2),:), &
              sim%particle(sim%maskw,:)/sim%npart, &
              sim%visu_bound_low(sim%plot2d_idx(1)),sim%visu_bound_up(sim%plot2d_idx(1)), sim%plot2d_bin(1), &
              sim%visu_bound_low(sim%plot2d_idx(2)),sim%visu_bound_up(sim%plot2d_idx(2)), sim%plot2d_bin(2), sim%tstep,&
              sll_world_collective, 0)
  end if
  
  
  !characteristics
  if (sim%coll_rank==0 .and. sim%num_visu_particles/=0 ) then
  SELECT CASE(sim%dimx)
  CASE(1)
  call plot_format_points3d(trim(sim%prefix)//'_xw', &
               sim%particle(sim%maskx(1),1:sim%num_visu_particles), &
               sim%particle(sim%maskx(2),1:sim%num_visu_particles), &
               sim%particle(sim%maskw,1:sim%num_visu_particles), &
               sim%tstep)
  CASE(2)
  call plot_format_points3d(trim(sim%prefix)//'pif_xyw', &
               sim%particle(sim%maskx(1),1:sim%num_visu_particles), &
               sim%particle(sim%maskx(2),1:sim%num_visu_particles), &
               sim%particle(sim%maskw,1:sim%num_visu_particles), &
               sim%tstep)
  CASE(3)
  call plot_format_points3d(trim(sim%prefix)//'_xyzw', &
               sim%particle(sim%maskx(1),1:sim%num_visu_particles), &
               sim%particle(sim%maskx(2),1:sim%num_visu_particles), &
               sim%particle(sim%maskx(3),1:sim%num_visu_particles), &
               sim%particle(sim%maskw,1:sim%num_visu_particles), &
               sim%tstep)
  END SELECT
  endif
  
end subroutine 


function visu_bound_low(sim, dim)  result(bound)
  class(sll_simulation_general_vlasov_poisson_pif), intent(in) :: sim 
  sll_int32, intent(in) :: dim
  sll_real64 :: bound
  
  if (dim<=sim%dimx) then
    bound=0
  elseif (dim>sim%dimx .and. dim <=2*sim%dimx) then
    bound=-4.0_f64
  endif
end function 

function visu_bound_up(sim, dim) result(bound)
  class(sll_simulation_general_vlasov_poisson_pif), intent(in) :: sim 
  sll_int32, intent(in) :: dim
  sll_real64 :: bound
  
  if (dim<=sim%dimx) then
    bound=sim%boxlen(dim)
  elseif (dim>sim%dimx .and. dim <=2*sim%dimx) then
    bound=4.0_f64
  endif
end function 

subroutine load_landau_sum_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim 
  sim%particle(sim%maskw,:)=(1.0+sim%eps*sum(cos(diag_dot_matrix_real64(sim%kmode,sim%particle(sim%maskx,:))),1))
  sim%particle(sim%maskw,:)=sim%particle(sim%maskw,:)/sim%prior_weight(:)  
  if (sim%coll_rank==0) print *, "Load landau sum"
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
  sll_int32 :: idx
  do idx=1, sim%npart_loc
    sim%particle(sim%maskw,idx)=1+sim%eps*cos(sum(sim%kmode(:)*sim%particle(sim%maskx,idx),1))  
  end do
  sim%particle(sim%maskw,:)=sim%particle(sim%maskw,:)/sim%prior_weight 
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
 
 sll_int32, dimension(2) :: PLOT2D_IDX=0, PLOT2D_BIN=50
 sll_int32 :: NUM_VISU_PARTICLES=0 , WRITE_PHI=0
 character(len=32) :: TESTCASE
 character(len=64) :: PREFIX
     
     sll_int32, parameter  :: input_file = 99
     sll_int32             :: IO_stat,ierr,idx
 
     namelist /sim_params/ dt, NUM_PARTICLES, DIMENSION, NUM_TIMESTEPS, &
                          NUM_MODES, QoverM, EPSILON, CONTROLVARIATE, &
                          TIME_INTEGRATOR_ORDER,RND_OFFSET,K,L,B0,&
                          TESTCASE,&
                          PLOT2D_IDX,PLOT2D_BIN, NUM_VISU_PARTICLES, PREFIX,WRITE_PHI
     
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

     if (any(B0(1:sim%dimx)/=0)) then
         SLL_ALLOCATE(sim%B0(sim%dimx),ierr)
         sim%B0=B0(1:sim%dimx)
     endif

     !Output
     if (len(trim(PREFIX))/=0) then
       sim%PREFIX=trim(PREFIX)
     endif
     
     sim%WRITE_PHI=WRITE_PHI
     
     if (all(PLOT2D_IDX/=0)) then
         if (maxval(PLOT2D_IDX)>sim%dimx*2+1) then
           print *, "PLOT2D_IDX is out of bound:",  PLOT2D_IDX
           else
           sim%plot2d_idx=PLOT2D_IDX
           sim%plot2d_bin=PLOT2D_BIN
         endif
     endif
     
     sim%num_visu_particles=NUM_VISU_PARTICLES
     
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

  SLL_ALLOCATE( sim%particle(2*sim%dimx+1,sim%npart_loc), ierr )

  !Generate random numbers
  seed = 10 + sim%RND_OFFSET + sim%coll_rank*sim%npart_loc
  do idx=1,sim%npart_loc
    !call i8_sobol_generate ( int(dimx,8) , npart, RND_OFFSET , particle(1:2*dimx,:))
    call i8_sobol( int(2*sim%dimx,8), seed, sim%particle(1:2*sim%dimx,idx) )
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

  sll_int32  :: idx,jdx,ierr
  sll_real64 :: v_cdf

  !scale to boxlength
  do idx=1,sim%npart_loc
    sim%particle(sim%maskx,idx) = sim%boxlen * sim%particle(sim%maskx,idx)
  end do

  !load sll_m_gaussian profile
  do idx=1,size(sim%particle,2)
   do jdx=1, size(sim%maskv)
    v_cdf = sim%particle(sim%maskv(jdx),idx)
    call normal_cdf_inv( v_cdf, 0.0_f64 , 1.0_f64, sim%particle(sim%maskv(jdx),idx) )
   end do
  end do

  !set temperature and impulse
  ! do idx=1, size(sim%maskv)
  ! print *, sum(sim%particle(sim%maskv(idx),:))/sim%npart
   !call match_moment_1D_linear_real64(sim%particle(sim%maskv(idx),:), 0.0_f64, 1.0_f64)
  ! print *, sum(sim%particle(sim%maskv(idx),:))/sim%npart
  ! end do
  
   !print *,minval(sim%particle(sim%maskx,:)), maxval(sim%particle(sim%maskx,:))
   !  print *,minval(sim%particle(sim%maskv,:)), maxval(sim%particle(sim%maskv,:))

  SLL_ALLOCATE( sim%prior_weight(1:sim%npart_loc), ierr )
  sim%prior_weight = 1.0 / product(sim%boxlen)

end subroutine init_particle_prior_maxwellian_generalvp_pif


!>Only for separable lagrangian, no magnetic field
subroutine symplectic_rungekutta_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim

  sll_int32 :: rkidx
  sll_real64 :: t !time
  SLL_ASSERT(size(sim%rk_d)==size(sim%rk_c))
  
    
  !Loop over all stages
   t=(sim%tstep-1)*sim%dt
    do rkidx=1, size(sim%rk_d)
     
     !Set control variate if used
     call sim%update_weight()
     
     !Charge assignement, get right hand side
     sim%rhs=sim%SOLVER%get_rhs_particle(sim%particle(sim%maskxw,:))/sim%npart
     !mpi allreduce
     call sll_collective_globalsum(sll_world_collective, sim%rhs)

     sim%solution=sim%solve_field(sim%rhs)
     
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
  sll_int32 :: stage
  sll_real64 :: t, h !time
  type(comp_coeff_sym_sym)  :: compc
  
  SELECT CASE( sim%time_integrator_order  )
    CASE(1)
     call compc%init(1,1)
    CASE(2:4)
     call compc%init(2,3)
    CASE(5:8)
     call compc%init(6,9) !(6,7)
    CASE(9:10)
     call compc%init(8,17) !(8,15)
    CASE(11:12)
     call compc%init(10,35)
   CASE DEFAULT
    call compc%init(2,3)
  END SELECT
         
  !Loop over all stages
    t=(sim%tstep-1)*sim%dt
     do stage=1, compc%stages;
     !Set control variate if used
      call sim%update_weight()
     
     !Charge assignement, get right hand side
      sim%rhs=sim%SOLVER%get_rhs_particle(sim%particle(sim%maskxw,:))/sim%npart
     !mpi allreduce
     call sll_collective_globalsum(sll_world_collective, sim%rhs)
     sim%solution=sim%solve_field(sim%rhs)
    
       if (stage==1) then
         call sim%calculate_diagnostics()
       endif

      h=compc%gamma(stage)*sim%dt
    
     E= sim%E(sim%particle(sim%maskx,:),t) + &  ! Electric field external
            (sim%SOLVER%eval_gradient(sim%particle(sim%maskx,:),sim%solution)) !Electric field selfconsistent
      B=sim%B(sim%particle(sim%maskx,:), t) !External magnetic field
      
     sim%ExB=0 
      
     sim%particle(sim%maskv,:)=exp_skew_product2( normalize(B), &
                      sim%particle(sim%maskv,:) + h*sim%qm/2.0*E, h*(-sim%qm)*l2norm(B) ) + h*sim%qm/2.0*E    
    !x_{k+1}= x_k + dt*v_{k+1}
    sim%particle(sim%maskx,:)=sim%particle(sim%maskx,:) +  h*sim%particle(sim%maskv,:);
      
      t=(sim%tstep-1)*sim%dt+sum(compc%gamma(1:stage))*sim%dt;
     end do
     t=(sim%tstep)*sim%dt;
     
end subroutine YSun_g2h_generalvp_pif


function normalize( v ) result(vn)
 sll_real64, dimension(:,:), intent(in) :: v
 sll_real64, dimension(size(v,1),size(v,2)) :: vn
 sll_real64, dimension(size(v,2)) :: normv
 sll_int32 :: idx
 normv=l2norm(v)
 
 do idx=1,size(v,2)
   if (normv(idx)/=0) then
     vn(:,idx)=v(:,idx)/normv(idx)
   endif
 end do
end function normalize


subroutine heun_generalvp_pif(sim)
    class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
    sll_real64 :: t
    sll_real64, dimension(2*sim%dimx,sim%npart_loc) :: xv0
    sll_real64, dimension(sim%dimx, sim%npart_loc) :: E
     t=(sim%tstep-1)*sim%dt;

     xv0(1:sim%dimx,:)=sim%particle(sim%maskx,:)
     xv0(sim%dimx+1:2*sim%dimx,:)=sim%particle(sim%maskv,:)
     
    
     call sim%update_weight()
     sim%rhs=sim%SOLVER%get_rhs_particle(sim%particle(sim%maskxw,:))/sim%npart
     call sll_collective_globalsum(sll_world_collective, sim%rhs)
     sim%solution=sim%solve_field(sim%rhs)
     E= sim%E(sim%particle(sim%maskx,:),t) + (sim%SOLVER%eval_gradient(sim%particle(sim%maskx,:),sim%solution)) 
     call sim%calculate_diagnostics()

     !Stage 1
     sim%particle(sim%maskv,:)=sim%particle(sim%maskv,:)+ sim%dt*sim%qm*(sim%vxB(sim%particle(sim%maskx,:),sim%particle(sim%maskv,:),t) + E) 
     sim%particle(sim%maskx,:)=sim%particle(sim%maskx,:)+ sim%dt*sim%particle(sim%maskv,:)
     !stage 2
     call sim%update_weight()
     sim%rhs=sim%SOLVER%get_rhs_particle(sim%particle(sim%maskxw,:))/sim%npart
     call sll_collective_globalsum(sll_world_collective, sim%rhs)
     sim%solution=sim%solve_field(sim%rhs)
     E= sim%E(sim%particle(sim%maskx,:),t) + (sim%SOLVER%eval_gradient(sim%particle(sim%maskx,:),sim%solution)) 
     
     sim%particle(sim%maskv ,:)=xv0(sim%dimx+1:2*sim%dimx,:)/2.0 +  sim%particle(sim%maskv,:)/2.0+  &
              sim%dt*sim%qm*(sim%vxB(sim%particle(sim%maskx,:),sim%particle(sim%maskv,:),t) + E)/2 

     sim%particle(sim%maskx,:)=xv0(1:sim%dimx,:)/2.0 +  sim%particle(sim%maskx,:)/2.0 + &
                                          sim%dt*sim%particle(sim%maskv,:)/2.0

end subroutine heun_generalvp_pif




subroutine rk4_generalvp_pif(sim)
    class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
    sll_real64 :: t
     sll_real64, dimension(sim%dimx,sim%npart_loc) :: k1_xx, k1_vx,k2_xx, &
                    k2_vx,k3_xx, k3_vx,E,k4_xx,k4_vx, xx_up,vx_up
    t=(sim%tstep-1)*sim%dt
! 

!      
!     
!             k1_xx=sim%dt*sim%particle(sim%maskv,:);
!             k1_vx=sim%dt*fx*qm;
!             !Stage 2
!             xx_up=sim%particle(sim%maskx,:)+0.5*k1_xx
!             vx_up=sim%particle(sim%maskv,:)+0.5*k1_vx;
!             
!             call sim%update_weight()
!      sim%rhs=sim%SOLVER%get_rhs_particle(sim%particle(sim%maskxw,:))/sim%npart
!      call sll_collective_globalsum(sll_world_collective, sim%rhs)
!      sim%solution=sim%SOLVER%solve_quasineutral(sim%rhs)
!             
!             
!             
!             k2_xx=sim%dt*vx_up;
!             k2_vx=sim%dt*fx*qm;
!             %Stage 3
!             xx_up=sim%particle(sim%maskx,:)+0.5*k2_xx
!             vx_up=sim%particle(sim%maskv,:)+0.5*k2_vx;
!             rhs(:,tstep+1)=fieldsolver.rhs_particle(xx_up,weight,pcell);
!             solution(:,tstep+1)=fieldsolver.solve_poisson(rhs(:,tstep+1)+rhs_CV-rhs_ION);
!             fx=fieldsolver.eval_gradient(solution(:,tstep+1),xx_up,pcell)+E_ext(xx,t+0.5*dt);
!             CV_new=CV(xx_up,vx_up);
!             weight=forward_weight_CV(weight, weight_const, CV_old, CV_new);CV_old=CV_new;
!             
!             k3_xx=sim%dt*vx_up;
!             k3_vx=sim%dt*fx*qm;
!             %Stage 4
!             xx_up=sim%particle(sim%maskx,:)+k3_xx
!             vx_up=sim%particle(sim%maskv,:)+k3_vx;
!             [rhs(:,tstep+1)]=fieldsolver.rhs_particle(xx_up,weight,pcell);
!             solution(:,tstep+1)=fieldsolver.solve_poisson(rhs(:,tstep+1)+rhs_CV-rhs_ION);
!             fx=fieldsolver.eval_gradient(solution(:,tstep+1),xx_up,pcell)+E_ext(xx,t+sim%dt);
!             CV_new=CV(xx_up,vx_up);
!             weight=forward_weight_CV(weight, weight_const, CV_old, CV_new);CV_old=CV_new;
!             
!             k4_xx=sim%dt*vx_up;
!             k4_vx=sim%dt*fx*qm;
!             sim%particle(sim%maskx,:)=sim%particle(sim%maskx,:)+(k1_xx + 2*k2_xx+2*k3_xx+ k4_xx)/6.0
!             sim%particle(sim%maskv,:)=sim%particle(sim%maskv,:)+(k1_vx + 2*k2_vx+2*k3_vx+ k4_vx)/6;
!             
!             
!             
!             
!             
!             
!             
!             
!             
!             
!             
!             
!             
!             CV_new=CV(xx_up,vx_up);
!             weight=forward_weight_CV(weight, weight_const, CV_old, CV_new);CV_old=CV_new;
!             
!             if (COVAR==0)
!                 [rhs(:,tstep+1)]=fieldsolver.rhs_particle(xx_up,weight,pcell);
!             else
!                 [rhs(:,tstep+1),rhs_cov(:,:,tstep+1)]=fieldsolver.rhs_particle(xx_up,weight,pcell);
!             end
!             solution(:,tstep+1)=fieldsolver.solve_poisson(rhs(:,tstep+1)+rhs_CV-rhs_ION);
!             fx=fieldsolver.eval_gradient(solution(:,tstep+1),xx_up,pcell)+E_ext(xx,t+dt);
!             xx=xx_up;
!             vx=vx_up;
! 
! 


end subroutine rk4_generalvp_pif











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
!    
  c(1,:) = w(1,:) + w(1,:)*(v(2,:)**2 + v(3,:)**2)*(cos(omega) - 1) + v(2,:)*w(3,:)*sin(omega) - v(3,:)*w(2,:)*sin(omega) - v(1,:)*v(2,:)*w(2,:)*(cos(omega) - 1) - v(1,:)*v(3,:)*w(3,:)*(cos(omega) - 1) 
  c(2,:) = w(2,:) + w(2,:)*(v(1,:)**2 + v(3,:)**2)*(cos(omega) - 1) - v(1,:)*w(3,:)*sin(omega) + v(3,:)*w(1,:)*sin(omega) - v(1,:)*v(2,:)*w(1,:)*(cos(omega) - 1) - v(2,:)*v(3,:)*w(3,:)*(cos(omega) - 1) 
  c(3,:) = w(3,:) + w(3,:)*(v(1,:)**2 + v(2,:)**2)*(cos(omega) - 1) + v(1,:)*w(2,:)*sin(omega) - v(2,:)*w(1,:)*sin(omega) - v(1,:)*v(3,:)*w(1,:)*(cos(omega) - 1) - v(2,:)*v(3,:)*w(2,:)*(cos(omega) - 1)  
   
  end function


  function exp_skew_product( v, w) result(c)
  sll_real64, dimension(:,:), intent(in) :: v, w
  sll_real64, dimension(3,size(v,2)) :: c
  sll_real64, dimension(size(v,2)) :: vv
  sll_comp64, dimension(size(v,2)) :: sqrtvv
  
  
  vv=-(v(1,:)**2+v(2,:)**2+v(3,:)**2)
  sqrtvv=sll_i1*sqrt(-vv)
  
     c(1,:) = real((w(1,:)*exp(-sqrtvv)*(v(1,:)**2*exp(sqrtvv)*2.0_f64+v(3,:)**2*exp(sqrtvv*2.0_f64)+v(3,:)**2*exp(sqrtvv*2.0_f64)+v(3,:)**2+v(3,:)**2)*&
     (-1.0_f64/2.0_f64))/vv+1.0_f64/sqrtvv**(3.0_f64)*w(2,:)*exp(-sqrtvv)*(exp(sqrtvv)-1.0_f64)*&
     (v(3,:)**3*exp(sqrtvv)+v(1,:)**2*v(3,:)+v(3,:)**2*v(3,:)+v(3,:)**3+v(1,:)**2*v(3,:)*&
     exp(sqrtvv)+v(3,:)**2*v(3,:)*exp(sqrtvv)-sqrtvv*v(1,:)*v(3,:)+sqrtvv*v(1,:)*v(3,:)*exp(sqrtvv))*(1.0_f64/2.0_f64)+(w(3,:)*&
     exp(-sqrtvv)*(exp(sqrtvv)-1.0_f64)*(sqrtvv*v(3,:)-v(1,:)*v(3,:)+sqrtvv*v(3,:)*exp(sqrtvv)+v(1,:)*v(3,:)*exp(sqrtvv))*(1.0_f64/2.0_f64))/vv)
  
     c(2,:) =real( (w(2,:)*exp(-sqrtvv)*(v(3,:)**2*exp(sqrtvv)*2.0_f64+v(1,:)**2*exp(sqrtvv*2.0_f64)+v(3,:)**2*exp(sqrtvv*2.0_f64)+v(1,:)**2+v(3,:)**2)*&
     (-1.0_f64/2.0_f64))/vv-1.0_f64/sqrtvv**(3.0_f64)*w(1,:)*exp(-sqrtvv)*(exp(sqrtvv)-1.0_f64)*(v(3,:)**3*&
     exp(sqrtvv)+v(1,:)**2*v(3,:)+v(3,:)**2*v(3,:)+v(3,:)**3+v(1,:)**2*v(3,:)*exp(sqrtvv)+v(3,:)**2*v(3,:)*&
     exp(sqrtvv)+sqrtvv*v(1,:)*v(3,:)-sqrtvv*v(1,:)*v(3,:)*exp(sqrtvv))*(1.0_f64/2.0_f64)-&
     (w(3,:)*exp(-sqrtvv)*(exp(sqrtvv)-1.0_f64)*(sqrtvv*v(1,:)+v(3,:)*v(3,:)+sqrtvv*v(1,:)*&
     exp(sqrtvv)-v(3,:)*v(3,:)*exp(sqrtvv))*(1.0_f64/2.0_f64))/vv)
     
     
     c(3,:) = real((w(3,:)*exp(-sqrtvv)*(v(3,:)**2*exp(sqrtvv)*2.0_f64+v(1,:)**2*exp(sqrtvv*2.0_f64)+v(3,:)**2*exp(sqrtvv*2.0_f64)+v(1,:)**2+v(3,:)**2)*(-1.0_f64/2.0_f64))/&
     vv-1.0_f64/sqrtvv**(3.0_f64)*w(2,:)*exp(-sqrtvv)*(exp(sqrtvv)-1.0_f64)*(v(1,:)**3*exp(sqrtvv)+v(1,:)*v(3,:)**2+v(1,:)*v(3,:)**2+v(1,:)**3+v(1,:)&
     *v(3,:)**2*exp(sqrtvv)+v(1,:)*v(3,:)**2*exp(sqrtvv)+sqrtvv*v(3,:)*v(3,:)-sqrtvv*v(3,:)*v(3,:)*exp(sqrtvv))*(1.0_f64/2.0_f64)+&
     1.0_f64/sqrtvv**(3.0_f64)*w(1,:)*exp(-sqrtvv)*(exp(sqrtvv)-1.0_f64)*(v(3,:)**3*exp(sqrtvv)+v(1,:)**2*v(3,:)+v(3,:)*v(3,:)**2+&
     v(3,:)**3+v(1,:)**2*v(3,:)*exp(sqrtvv)+v(3,:)*v(3,:)**2*exp(sqrtvv)-sqrtvv*v(1,:)*v(3,:)+sqrtvv*v(1,:)*v(3,:)*exp(sqrtvv))*(1.0_f64/2.0_f64))
  end function
 
 
subroutine set_symplectic_rungekutta_coeffs_generalvp_pif(sim,rk_order)
  class(sll_simulation_general_vlasov_poisson_pif), intent(inout) :: sim
sll_int32, intent(in) :: rk_order
sll_int32 :: ierr
sll_real64,parameter :: rk4sx=0.17560359597982881702384390448573_f64 
!((2**(1.0/3.0) +2**(-1/3)-1)/6)

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
   CASE DEFAULT
END SELECT

if (sim%coll_rank==0) print *, "Symplectic Runge Kutta of order:", rk_order
! print *, "D_COEFFS:", sim%rk_d
! print *, "C_COEFFS:", sim%rk_c

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
      cv=sqrt(2.0_f64*sll_pi)**(-size(sim%maskv))*exp(-0.5_f64*sum(sim%particle(sim%maskv,:),1)**2)
    CASE (SLL_CONTROLVARIATE_MAXWELLIAN)
      cv=sqrt(2.0_f64*sll_pi)**(-size(sim%maskv))*exp(-0.5_f64*sum(sim%particle(sim%maskv,:),1)**2)
END SELECT

end function control_variate_generalvp_pif



subroutine write_result_generalvp_pif(sim)
  class(sll_simulation_general_vlasov_poisson_pif), intent(in) :: sim       
        integer :: idx,file_id,ierr


        if (sim%coll_rank==0) then

             call sll_new_file_id(file_id, ierr)

            !Write Data File
            !            open(file_id, file = plot_name//"_"//fin//'.dat' )
             open(file_id, file = trim(sim%prefix)//'_result.csv')
            !             open(file_id, file = './'//filename//'.csv')
!             write (file_id, *)  "#Full 1d1v Electrostatic PIC"
!             write (file_id,*)  "#Time steps:", timesteps
!             write (file_id,*)  "#Time stepwidth:", timestepwidth
!             write (file_id,*)  "#Marko particles:", coll_size*nparticles
!             write (file_id,*)  "#Particle Pusher:", particle_pusher
!             write (file_id,*)  "#Finite Elements: 2^(", log(real(mesh_cells,i64))/log(2.0_f64),")"
!             write (file_id,*)  "#Size of MPI Collective: ", coll_size
            write (file_id,*)  ' "time", "kineticenergy", "fieldenergy", "energy_error", "moment_error", "l2potential"'
            do idx=1,sim%tsteps

              write (file_id, '(ES32.16,ES32.16,ES32.16,ES32.16,ES32.16,ES32.16)' ) (idx-1)*sim%dt,sim%kineticenergy(idx), &
                               sim%fieldenergy(idx),sim%energy_error(idx),sim%moment_error(idx),sim%l2potential(idx)
            !                 write (file_id,*) (idx-1)*sim%dt,",",sim%kineticenergy(idx),",",sim%fieldenergy(idx),",",sim%energy_error(idx),",",sim%moment_error(idx),",", sim%l2potential(idx)
            enddo
            close(file_id)
        
          open(file_id, file = trim(sim%prefix)//'_result.gnu')
       
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
        
          write(file_id,*) "set term x11 6"
          write(file_id,*) "set logscale y"
          write(file_id,*) "plot '-' using 1:2 title 'Electrostatic Potential' with lines"
          do idx = 1, sim%tsteps
            write(file_id,*)  (idx-1)*sim%dt, sim%l2potential(idx)
          end do
          write(file_id,*) "e"
          write(file_id,*) "   "
          
          close(file_id)

        
        endif
end subroutine write_result_generalvp_pif








end module sll_m_sim_pif_vp_ndnv_cart
