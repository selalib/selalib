program vp_csl_bloc_unif
!vlasov poisson 1D x 1D solver with block uniform mesh in velocity
#include "sll_working_precision.h"
#include "sll_field_2d.h"
#include "sll_memory.h"
  use sll_constants
  !use distribution_function
  !use sll_diagnostics
  use bgk_mesh_construction
  implicit none

  sll_int :: nc_eta1,nc_eta2,mesh_case,N_x1,N_x2

  sll_real64 ::geom_x(2,2),geom_eta(2,2),alpha_mesh

  sll_real64,dimension(:,:),pointer::x1n_array,x2n_array,x1c_array,x2c_array
  sll_real64,dimension(:,:),pointer::jac_array
  sll_real64,dimension(:,:,:),pointer::integration_points

  sll_real64 :: delta_x1,delta_x2
  sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max,delta_eta1,delta_eta2
  sll_real64 :: x1_min,x1_max,x2_min,x2_max


  mesh_case = 1
  !visu_step = 10
  !test_case = 4
  !rho_case = 2
  !phi_case = 3
 
  alpha_mesh = 1.e-2_f64 !0.1_f64
  
  N_x1 = 64
  N_x2 = 64
  !dt = 0.1_f64
  !nb_step = 600
  
  nc_eta1=N_x1
  nc_eta2=N_x2
  
  !domain for classical vlasov-poisson test case
  
  x1_min = 0._f64
  x1_max = 4._f64*sll_pi
  x2_min = -10._f64
  x2_max = -x2_min

  
  
  eta1_min =  0.0_f64
  eta1_max = 1.0_f64 ! 0.15_f64*x1_max! 1.0_f64
  eta2_min =  0.0_f64
  eta2_max =  1.0_f64
  
  geom_x(1,1)=x1_min
  geom_x(2,1)=x1_max
  geom_x(1,2)=x2_min
  geom_x(2,2)=x2_max
    

  geom_eta(1,1)=eta1_min
  geom_eta(2,1)=eta1_max
  geom_eta(1,2)=eta2_min
  geom_eta(2,2)=eta2_max


  delta_x1 = (x1_max-x1_min)/real(N_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(N_x2,f64)

  delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
  delta_eta2 = (eta2_max-eta2_min)/real(nc_eta2,f64)
  
  
  
  call construct_bgk_mesh(nc_eta1,nc_eta2,mesh_case,&
   &x1n_array,x2n_array,x1c_array,x2c_array,jac_array,integration_points,&
   &geom_x,geom_eta,alpha_mesh,N_x1,N_x2)
  
  

end program vp_csl_bloc_unif
