program VP_DG
#include "sll_working_precision.h"

  use sll_nu_cart_mesh
  use gausslobatto
  use poisson4dg
  use mod_sparse

  implicit none

  !rho,E,phi
  sll_real64,dimension(:),allocatable :: rho,field_e,phi
  !distribution function
  sll_real64,dimension(:,:),allocatable :: dist
  !number of time steb, number of cells un direction x and v
  sll_int32 :: nx,nv,nb_step
  !boudary in direction x and v
  sll_real64 :: x_min,x_max,v_min,v_max
  !mesh
  type(non_unif_cart_mesh) :: mesh
  !coefficients for fluxes
  sll_real64 :: c11,c12,c22
  !number of Gauss-Lobatto points
  sll_int32 :: ng
  !CSC matrix array
  type(t_col) :: matvp,matvm
  !time step and finalt time
  sll_real64 :: dt,tf
  !Gauss-Lobatto
  type(gausslobatto1D) :: gausslob

  !indexes for loops
  sll_int32 :: j,k1,k2

  !definition of geometry and data
  x_min=-1.0d0
  x_max=1.0d0
  v_min=-1.0d0
  v_max=1.0d0
  nx=100
  nv=100
  ng=5
  call init_gausslobatto_1d(ng,gausslob)
  call init_nu_cart_mesh(nx,nv,mesh)
  mesh%d_etat1=(x_max-x_min)/real(nx,8)
  mesh%d_etat2=(v_max-v_min)/real(nv,8)
  call fill_node_nuc_mesh(x_min,v_min,mesh)

  allocate(rho(nx*ng),phi(nx*ng),field_e(nx*ng),dist(nx*ng,nv*ng))

  !definition or time step, delta_t and final time
  dt=0.01d0
  tf=1.0d0
  nb_step=ceiling(tf/dt)

  !flux coefficients
  c22=0.0d0
  c12=0.0d0
  c11=0.0d0
  !build the matrixes for the Poisson-problem
  !((D+F_E).M-¹.(D-F_\Phi)^T - C).\Phi = M.(\rho-1)
  !this is for v>0
  call poisson1d_matrix(gausslob,nx,mesh%jac(:,nv),c11,c12,c22,matvp)
  !this is for v<0
  c12=0.0d0
  call poisson1d_matrix(gausslob,nx,mesh%jac(:,nv),c11,c12,c22,matvm)

  !initialization of distribution f
  !must find a better one
  dist=1.0d0/((x_max-x_min)*(v_max-v_min))

  deallocate(rho,phi,field_e,dist)
  call delete_gausslobatto_1D(gausslob)
  call clear(matvm)
  call clear(matvp)
  call delete_nu_cart_mesh(mesh)

end program VP_DG
