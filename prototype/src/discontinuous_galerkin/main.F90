program VP_DG
#include "sll_working_precision.h"

  use sll_nu_cart_mesh
  use gausslobatto
  use poisson4dg
  use mod_sparse

  implicit none

  !rho,E,phi
  sll_real64,dimension(:),allocatable :: rho,field,phi
  !distribution function
  sll_real64,dimension(:,:),allocatable :: dist
  !number of time steb, number of cells un direction x and v
  sll_int32 :: nx,nv,nb_step
  !mesh
  type(non_unif_cart_mesh) :: mesh
  !coefficients for fluxes
  sll_real64 :: c11,c12,c22
  !number of Gauss-Lobatto points
  sll_int32 :: ng
  !CSC matrix array
  type(t_col) :: mat
  !total number of columns/rows
  sll_int32 :: glob_size
  !time step and finalt time
  sll_real64 :: dt,tf
  !Gauss-Lobatto
  type(gausslobatto1D) :: gausslob

  !indexes for loops
  sll_int32 :: j,k1,k2

  !definition of geometry and data
  nx=100
  nv=100
  ng=5
  call init_gausslobatto_1d(ng,gausslob)
  call init_nu_cart_mesh(nx,nv,mesh)
  mesh%d_etat1=2.0d0/real(nx,8)
  mesh%d_etat2=2.0d0/real(nv,8)
  call fill_node_nuc_mesh(-1.0d0,-1.0d0,mesh)

  !definition or time step, delta_t and final time
  dt=0.01d0
  tf=1.0d0
  nb_step=ceiling(tf/dt)

  glob_size=ng*nx

  !flux coefficients
  c22=0.0d0
  c12=0.0d0
  c11=0.0d0

  !build the matrixes for the Poisson-problem
  !((D+F_E).M-¹.(D-F_\Phi)^T - C).\Phi = M.(\rho-1)

!
!  MUST SOLVE FOR V<0 AND V>0!
!

  !initialization of the matrix
  !the matrix is like the following :
  !
  !  [C1_|__       |]
  !  [ |_C2_|_      ]   where CI is the element I on mesh (cell)
  !  [   |_C3_|     ]
  !  [       ... _|_]
  !  [|         |_CN]
  !
  !And each CI can be writen as
  !
  !  [0..0 x D_{1;1}+x.........D_{1;n}  x....x]
  !  [.  . .  D_{2;1}            .      . 0..0]
  !  [.  . .    .                .      . .  .]
  !  [0..0 .    .                .      . .  .]
  !  [x....x  D_{n;1}.........D_{n;n}+x x 0..0]
  !
  !Note that here it is written in line but we want to use CSC (column)
  !WARNING : the code is written for c12=-c21



end program VP_DG
