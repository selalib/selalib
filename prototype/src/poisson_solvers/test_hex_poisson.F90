program test_hex_poisson

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
  
  use sll_constants
  use sll_hex_meshes
  use hex_poisson

  implicit none

  type(sll_hex_mesh_2d), pointer              :: mesh
  sll_real64, dimension(:),allocatable    :: second_term, rho, sol, phi
  sll_real64, dimension(:,:) ,allocatable :: matrix_poisson
  sll_int32                               :: num_cells, n_points, i 
  sll_int32                               :: ierr
  sll_real64                              :: x, y, erreur = 0._f64

  num_cells = 10

  n_points  = 1 + 3 * num_cells * (num_cells + 1) 

  SLL_ALLOCATE(second_term( n_points),ierr)      ! le b de Ax = b
  SLL_ALLOCATE(rho( n_points),ierr)              ! second terme de l'equation
  SLL_ALLOCATE(phi( n_points),ierr)              ! le x de Ax = b
  SLL_ALLOCATE(sol( n_points),ierr)              ! exact solution
  SLL_ALLOCATE(matrix_poisson( n_points,1+4*num_cells ) , ierr) ! le A de Ax = b
  
  ! initialization
  mesh => new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
  
  do i=1, n_points  
     x = mesh%cartesian_coord(1,i)
     y = mesh%cartesian_coord(2,i)
     rho(i) =  2._f64*(x**2-0.75_f64)*(6._f64*y**2-(x/sqrt(3._f64)-1._f64)**2&
          - ( (x+sqrt(3._f64)*0.5_f64)/sqrt(3._f64) + 0.5_f64)**2) & ! dyyP(x,y)
          + 2._f64*( ( y**2 - (x/sqrt(3._f64)-1._f64)**2 )* &
          (y**2 - ( (x+0.5_f64*sqrt(3._f64))/sqrt(3._f64)+0.5_f64)**2) &
          - 4._f64*x*( (x/sqrt(3._f64)-1._f64)* (y**2 - ( (x+0.5_f64*sqrt(3._f64))/sqrt(3._f64)+0.5_f64)**2) &
          + ( y**2 - (x/sqrt(3._f64)-1._f64)**2 )*( (x+sqrt(3._f64)*0.5_f64)/sqrt(3._f64) + 0.5_f64) )&
          - (x**2-0.75_f64)*( (y**2 - ( (x+0.5_f64*sqrt(3._f64))/sqrt(3._f64)+0.5_f64)**2)/sqrt(3._f64)&
          - 4._f64*(x/sqrt(3._f64)-1._f64)*( (x+sqrt(3._f64)*0.5_f64)/sqrt(3._f64) + 0.5_f64)&
          + ( y**2 - (x/sqrt(3._f64)-1._f64)**2 )/sqrt(3._f64))) ! dxxP(x,y)

     sol(i) = (x**2-0.75_f64) * ( y**2 - (x/sqrt(3._f64)-1._f64)**2 )* &
          ( y**2 - ( (x+0.5_f64*sqrt(3._f64))/sqrt(3._f64)+0.5_f64)**2)
  enddo

  call hex_matrix_poisson( matrix_poisson, mesh )

  call hex_second_terme_poisson( second_term, mesh, rho )
  
  ! solv( matrix_poisson, second_term, phi )  ! solving the linear system


  ! Error between the computed value and the expected value 
  do i=1, n_points  
     erreur = erreur + abs( sol(i) - phi(i))
  enddo

  SLL_DEALLOCATE_ARRAY(second_term,ierr)
  SLL_DEALLOCATE_ARRAY(rho,ierr)
  SLL_DEALLOCATE_ARRAY(sol,ierr)
  SLL_DEALLOCATE_ARRAY(matrix_poisson,ierr)


end program test_hex_poisson
