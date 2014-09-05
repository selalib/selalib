program test_hex_poisson

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
  
  use sll_constants
  use sll_hex_meshes
  use hex_poisson
  use pivotbande

  implicit none

  type(sll_hex_mesh_2d), pointer              :: mesh
  sll_real64, dimension(:),allocatable    :: second_term, rho, sol, phi, phi_end
  sll_real64, dimension(:,:) ,allocatable :: matrix_poisson, l, u
  sll_int32                               :: num_cells, n_points, i, k1, k2 
  sll_int32                               :: ierr, l1,l2, index_tab, global
  sll_real64                              :: x, y, erreur = 0._f64

  num_cells = 60

  n_points  = 1 + 3 * num_cells * (num_cells + 1) 

  SLL_ALLOCATE(second_term( n_points),ierr)      ! le b de Ax = b
  SLL_ALLOCATE(rho( n_points),ierr)              ! second terme de l'equation
  SLL_ALLOCATE(phi( n_points),ierr)              ! le x de Ax = b
  SLL_ALLOCATE(phi_end( n_points),ierr)       
  SLL_ALLOCATE(sol( n_points),ierr)              ! exact solution
  SLL_ALLOCATE(matrix_poisson( n_points,1 + 4*num_cells + 2 ) , ierr) ! le A de Ax = b
  !SLL_ALLOCATE(matrix_poisson( n_points,n_points), ierr)
  SLL_ALLOCATE(l( n_points,n_points), ierr)
  SLL_ALLOCATE(u( n_points,n_points), ierr)



  ! initialization
  mesh => new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
  
  do i=1, n_points  

     x = mesh%cartesian_coord(1,i)
     y = mesh%cartesian_coord(2,i)

     rho(i) = 12._f64*x**2*y**2 - 9._f64*y**2 + 4._f64/3._f64*x**4 + 3._f64 - 3._f64*x**2& ! DyyP
          + 2._f64*y**4 - 8._f64*y**2*x**2 + 10._f64/3._f64*x**4 - 10._f64*x**2 + 3._f64 - 3._f64*y**2 ! DxxP

     sol(i) = (x**2-0.75_f64) * ( y**2 - (x/sqrt(3._f64)-1._f64)**2 )* &
          ( y**2 - (x/sqrt(3._f64)+1._f64)**2)
  enddo

  rho = -rho

  call hex_matrix_poisson( matrix_poisson, mesh )

  call hex_second_terme_poisson( second_term, mesh, rho )
  
  ! call searchband(matrix_poisson,n_points,l1,l2)
  ! print*,l1,l2,n_points
  ! call factolub (matrix_poisson,l,u,n_points,l1,l2)
  ! call solvlub(l,u,phi,second_term,n_points,l1,l2)

  ! need to re-index phi : 

    do global = 1, mesh%num_pts_tot 
       k1 = mesh%hex_coord(1, global)
       k2 = mesh%hex_coord(2, global)
       call index_hex_to_global(mesh, k1, k2, index_tab)
       phi_end(global) = phi(index_tab)
    enddo

  ! Error between the computed value and the expected value 
  do i=1, n_points  
     erreur = erreur + abs( sol(i) - phi_end(i))**2
  enddo

  print*, "erreur", sqrt(erreur/real(num_cells,f64)**2)

  open(unit = 12, file="matrix_hex_poisson.txt", action="write", status="replace")

  do i=1, n_points
        x = mesh%cartesian_coord(1,i)
        y = mesh%cartesian_coord(2,i)
        write(12,*)  x,y,sol(i), phi_end(i)
  enddo

  close(12)

  SLL_DEALLOCATE_ARRAY(second_term,ierr)
  SLL_DEALLOCATE_ARRAY(rho,ierr)
  SLL_DEALLOCATE_ARRAY(phi,ierr)
  SLL_DEALLOCATE_ARRAY(phi_end,ierr)
  SLL_DEALLOCATE_ARRAY(sol,ierr)
  SLL_DEALLOCATE_ARRAY(matrix_poisson,ierr)
  SLL_DEALLOCATE_ARRAY(l,ierr)
  SLL_DEALLOCATE_ARRAY(u,ierr)


end program test_hex_poisson
