program test_hex_poisson

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_hex_poisson, only: &
    compute_hex_fields, &
    hex_matrix_poisson, &
    hex_second_terme_poisson

  use sll_m_hexagonal_meshes, only: &
    new_hex_mesh_2d, &
    sll_hex_mesh_2d

  use sll_m_pivotbande, only: &
    factolub_bande, &
    residue_bande, &
    solvlub_bande

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_hex_mesh_2d), pointer          :: mesh
  sll_real64, dimension(:),allocatable    :: second_term, rho, sol, phi, phi_end
  sll_real64, dimension(:),allocatable    :: uxn, uyn,dxuxn,dyuxn,dxuyn,dyuyn
  sll_real64, dimension(:,:) ,allocatable :: matrix_poisson, l, u
  sll_int32                               :: num_cells, n_points, i, k1, k2
  sll_int32                               :: ierr, l1,l2, index_tab, global
  sll_real64                              :: x, y, erreur = 0._f64
  sll_real64                              :: erreur1, erreur2, erreur3
  sll_real64                              :: t_init,t_inter, t_end, residu

  num_cells = 80


  n_points  = 1 + 3 * num_cells * (num_cells + 1)

  SLL_ALLOCATE(second_term( n_points),ierr)      ! le b de Ax = b
  SLL_ALLOCATE(rho( n_points),ierr)              ! second terme de l'equation
  SLL_ALLOCATE(phi( n_points),ierr)              ! le x de Ax = b
  SLL_ALLOCATE(uxn( n_points),ierr)
  SLL_ALLOCATE(uyn( n_points),ierr)
  SLL_ALLOCATE(dxuxn( n_points),ierr)
  SLL_ALLOCATE(dxuyn( n_points),ierr)
  SLL_ALLOCATE(dyuxn( n_points),ierr)
  SLL_ALLOCATE(dyuyn( n_points),ierr)
  SLL_ALLOCATE(phi_end( n_points),ierr)
  SLL_ALLOCATE(sol( n_points),ierr)              ! exact solution
  SLL_ALLOCATE(matrix_poisson( n_points,1 + 4*num_cells + 2 ) , ierr) ! le A de Ax = b
  SLL_ALLOCATE(l( n_points,1 + 4*num_cells + 2 ) , ierr)
  SLL_ALLOCATE(u( n_points,1 + 4*num_cells + 2 ) , ierr)
  ! SLL_ALLOCATE(matrix_poisson( n_points,n_points), ierr)
  ! SLL_ALLOCATE(l( n_points,n_points), ierr)
  ! SLL_ALLOCATE(u( n_points,n_points), ierr)

  ! initialization
  mesh => new_hex_mesh_2d( num_cells, 0._f64, 0._f64) 
  
  do i=1, n_points  

     x = mesh%cartesian_coord(1,i)
     y = mesh%cartesian_coord(2,i)

     ! rho(i) = -(2._f64*(x**2+y**2)**2 + 6._f64*(1._f64-2._f64*x**2-2._f64*y**2))

     ! sol(i) = x**2*y**4 - 2._f64/3._f64*x**4*y**2 - 0.75_f64 - 1.5_f64*x**2*&
     !      y**2 + x**6/9._f64 - 0.75_f64*(x**4 + y**4) + 1.5_f64*(x**2 + y**2)
     !   = (x**2-0.75_f64) * ( y**2 - (x/sqrt(3._f64)-1._f64)**2 )* &
     !     ( y**2 - (x/sqrt(3._f64)+1._f64)**2)

     rho(i) = -(4._f64*(x**2+y**2)*36._f64**2-4._f64*36._f64)*exp( - (x**2 + y**2)*36._f64 )
     sol(i) = +exp( - (x**2 + y**2)*36._f64 )

  enddo

  ! here rho = - laplacian(sol) , sol ~ phi

  call hex_matrix_poisson( matrix_poisson, mesh,1 )

  call hex_second_terme_poisson( second_term, mesh, rho )
  
  !call searchband(matrix_poisson,n_points,l1,l2)
  l1 = 2*num_cells + 1
  l2 = l1

  call cpu_time(t_init)
  call factolub_bande(matrix_poisson,l,u,n_points,l1,l2)
  call cpu_time(t_inter)
  call solvlub_bande(l,u,phi,second_term,n_points,l1,l2)
  call cpu_time(t_end)
  call residue_bande(matrix_poisson,phi,second_term,l1,l2,n_points,residu)

  print*, "for ",  n_points ," points, the computation time was : ",   t_end - t_init
  print*, " with a computation time for the factorisation, of :",t_inter-t_init
  print*, " residu : " , residu

  !call gauss_seidel_bande(matrix_poisson,phi,second_term,l1,l2,n_points)
  !call gauss_seidel(matrix_poisson,phi,second_term,n_points)

  ! need to re-index phi : 

    do global = 1, mesh%num_pts_tot 
       k1 = mesh%hex_coord(1, global)
       k2 = mesh%hex_coord(2, global)
       call mesh%index_hex_to_global(k1, k2, index_tab)
       phi_end(global) = phi(index_tab)
    enddo

  ! Error between the computed value and the expected value 
  do i=1, n_points  
    erreur = erreur + abs( sol(i) - phi_end(i))**2
  enddo

  print*, "erreur", sqrt(erreur/real(num_cells,f64)**2)


  ! testing computing the field for the guiding center model in a hex. mesh
  ! here : uxn = -dy(phi)  and uyn = +dx(phi)

  call compute_hex_fields(mesh,uxn,uyn,dxuxn,dyuxn,dxuyn,dyuyn,phi_end,1)

  erreur1 = 0._f64
  erreur2 = 0._f64
  do i=1, n_points  
     x = mesh%cartesian_coord(1,i)
     y = mesh%cartesian_coord(2,i)
     erreur1 = erreur1 + abs(-(2._f64*x*36._f64)*exp( - (x**2 + y**2)*36._f64 ) - uyn(i) )**2
     erreur2 = erreur2 + abs(-(2._f64*y*36._f64)*exp( - (x**2 + y**2)*36._f64 ) + uxn(i) )**2
  enddo
  print*, "erreur sur dx(phi)", sqrt(erreur1/real(num_cells,f64)**2)
  print*, "erreur sur dy(phi)", sqrt(erreur2/real(num_cells,f64)**2)


  ! testing computing the derivatives of the field for the guiding center model in a hex. mesh
  ! here : uxn = -dy(phi)  and uyn = +dx(phi)
  ! therefore : DxUx = -DxyUx ; DyUx = -DyyUx ; DxUy = DxxUy ; DyUy = DxyUy 

  erreur1 = 0._f64
  erreur2 = 0._f64
  erreur3 = 0._f64
  do i=1, n_points  
     x = mesh%cartesian_coord(1,i)
     y = mesh%cartesian_coord(2,i)
     erreur1 = erreur1 + abs(-(72._f64-(72._f64*x)**2)*&
          exp( - (x**2 + y**2)*36._f64 ) - dxuyn(i) )**2 
     erreur2 = erreur2 + abs(-(72._f64-(72._f64*y)**2)*&
          exp( - (x**2 + y**2)*36._f64 ) + dyuxn(i) )**2
     erreur3 = erreur3 + abs(+(2._f64*36._f64*2._f64*y*x*36._f64)*&
          exp( - (x**2 + y**2)*36._f64 ) + dxuxn(i) )**2
  enddo
  print*, "erreur sur dxx(phi)", sqrt(erreur1/real(num_cells,f64)**2)
  print*, "erreur sur dyy(phi)", sqrt(erreur2/real(num_cells,f64)**2)
  print*, "erreur sur dxy(phi)", sqrt(erreur3/real(num_cells,f64)**2)



  ! open(unit = 11, file="matrix_hex_poisson.txt", action="write", status="replace")

  ! do i=1, n_points
  !    do j=1,1+l1+l2
  !        write(11,*) i,j,matrix_poisson(i,j)
  !    enddo
  !    write(11,*)
  ! enddo
  
  ! close(11)

  open(unit = 12, file="test_hex_poisson.txt", action="write", status="replace")

  do i=1, n_points
        x = mesh%cartesian_coord(1,i)
        y = mesh%cartesian_coord(2,i)
        write(12,*)  x,y,sol(i), phi_end(i)
  enddo

  close(12)

  SLL_DEALLOCATE_ARRAY(second_term,ierr)
  SLL_DEALLOCATE_ARRAY(rho,ierr)
  SLL_DEALLOCATE_ARRAY(uxn,ierr)
  SLL_DEALLOCATE_ARRAY(uyn,ierr)
  SLL_DEALLOCATE_ARRAY(dxuxn,ierr)
  SLL_DEALLOCATE_ARRAY(dxuyn,ierr)
  SLL_DEALLOCATE_ARRAY(dyuxn,ierr)
  SLL_DEALLOCATE_ARRAY(dyuyn,ierr)
  SLL_DEALLOCATE_ARRAY(phi,ierr)
  SLL_DEALLOCATE_ARRAY(phi_end,ierr)
  SLL_DEALLOCATE_ARRAY(sol,ierr)
  SLL_DEALLOCATE_ARRAY(matrix_poisson,ierr)
  SLL_DEALLOCATE_ARRAY(l,ierr)
  SLL_DEALLOCATE_ARRAY(u,ierr)


end program test_hex_poisson
