program test_hex_hermite

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_constants
  use sll_interpolation_hex_hermite
  use euler_2d_hex
  use sll_hex_meshes
  use hex_poisson
  use pivotbande
  implicit none

  type(sll_hex_mesh_2d), pointer          :: mesh, mesh2
  sll_real64, dimension(:,:), allocatable :: deriv
  sll_real64, dimension(:),allocatable    :: dxuxn,dyuxn,dxuyn,dyuyn
  sll_real64, dimension(:),allocatable    :: dxuxn_center,dyuxn_center,dxuyn_center,dyuyn_center
  sll_real64, dimension(:),allocatable    :: dxuxn_edge,dyuxn_edge,dxuyn_edge,dyuyn_edge
  sll_real64, dimension(:),allocatable    :: second_term, phi,uxn,uyn,phi_interm
  sll_real64, dimension(:),allocatable    :: second_term2, phi2,uxn2,uyn2,phi2_interm,rho2
  sll_real64, dimension(:),allocatable    :: dxuxn2,dyuxn2,dxuyn2,dyuyn2
  sll_real64, dimension(:),allocatable    :: dxuxn_center2,dyuxn_center2,dxuyn_center2,dyuyn_center2
  sll_real64, dimension(:),allocatable    :: dxuxn_edge2,dyuxn_edge2,dxuyn_edge2,dyuyn_edge2
  sll_real64, dimension(:),allocatable    :: second_term_center, phi_center
  sll_real64, dimension(:),allocatable    :: second_term_edge, phi_edge
  sll_real64, dimension(:),allocatable    :: uxn_center,uyn_center
  sll_real64, dimension(:),allocatable    :: uxn_edge,uyn_edge
  sll_real64, dimension(:),allocatable    :: phi_edge_interm, phi_center_interm
  sll_real64, dimension(:,:) ,allocatable :: matrix_poisson, l, u
  sll_real64, dimension(:,:) ,allocatable :: matrix_poisson_center, l_center, u_center
  sll_real64, dimension(:,:) ,allocatable :: matrix_poisson_edge, l_edge, u_edge
  sll_int32                               :: i,j, k1, k2, index_tab
  sll_int32                               :: l1,l2, type  
  sll_int32    :: i1,i2,i3
  sll_int32    :: num_cells, n_points, n_triangle, n_edge, n_points2
  sll_real64   :: center_mesh_x1, center_mesh_x2, radius
  sll_int32    :: nloops,count, ierr, EXTRA_TABLES = 1 ! put 1 for num_method = 15
  sll_real64   :: dt
  sll_real64   :: tmax
  sll_real64   :: t
  sll_real64   :: t_init, t_end
  sll_real64   :: t1,t2,t3,t4,t5,t6
  sll_real64   :: step , aire, h1, h2, f_min, x ,y,xx, yy
  sll_real64   :: r11,r12,r21,r22,det
  sll_int32    :: p = 6!-> degree of the approximation for the derivative 
  ! distribution at time n
  sll_int32    :: num_method = 9
  logical      :: inside
  character(len = 50) :: filename
  character(len = 4)  :: filenum
  sll_real64,dimension(:),allocatable :: rho_tn, rho_center_tn, rho_edge_tn
  ! distribution at time n + 1
  sll_real64,dimension(:),allocatable :: rho_tn1,rho_center_tn1,rho_edge_tn1

  center_mesh_x1 = 0._f64
  center_mesh_x2 = 0._f64

  radius = 14._f64

  call print_method(num_method)
  
  do num_cells = 80,80,40 
     
     t = 0._f64
     tmax  = 100._f64
     dt    = 0.1_f64!*20._f64 !/ real(num_cells,f64)  
     !cfl   = radius * dt / ( radius / real(num_cells,f64)  )
     nloops = 0
     count  = 0

     !*********************************************************
     !             allocation
     !*********************************************************

     n_points   = 1 + 3 * num_cells * (num_cells + 1)  
     n_triangle = 6 * num_cells * num_cells
     n_edge    = 3 * num_cells * ( 3 * num_cells + 1 )

     print*," minimum number of points computed on the mesh: ", n_points

     l1 = 2*num_cells+1
     l2 = l1

     step = radius / real(num_cells,f64)
     aire = step**2*sqrt(3._f64)*0.25_f64

     allocate( deriv(1:6,n_points) )

     SLL_ALLOCATE(rho_tn( n_points),ierr)
     SLL_ALLOCATE(rho_tn1( n_points ),ierr)
     SLL_ALLOCATE(rho_center_tn ( n_triangle),ierr)
     SLL_ALLOCATE(rho_center_tn1( n_triangle),ierr)
     SLL_ALLOCATE(rho_edge_tn ( n_edge),ierr)
     SLL_ALLOCATE(rho_edge_tn1( n_edge),ierr)

     SLL_ALLOCATE(uxn( n_points),ierr)
     SLL_ALLOCATE(uyn( n_points ),ierr)
     SLL_ALLOCATE(uxn_center( n_triangle),ierr)
     SLL_ALLOCATE(uyn_center( n_triangle),ierr)
     SLL_ALLOCATE(uxn_edge( n_edge),ierr)
     SLL_ALLOCATE(uyn_edge( n_edge),ierr)

     SLL_ALLOCATE(dxuxn( n_points),ierr)
     SLL_ALLOCATE(dxuyn( n_points ),ierr)
     SLL_ALLOCATE(dxuxn_center( n_triangle),ierr)
     SLL_ALLOCATE(dxuyn_center( n_triangle),ierr)
     SLL_ALLOCATE(dxuxn_edge( n_edge),ierr)
     SLL_ALLOCATE(dxuyn_edge( n_edge),ierr)
     SLL_ALLOCATE(dyuxn( n_points),ierr)
     SLL_ALLOCATE(dyuyn( n_points ),ierr)
     SLL_ALLOCATE(dyuxn_center( n_triangle),ierr)
     SLL_ALLOCATE(dyuyn_center( n_triangle),ierr)
     SLL_ALLOCATE(dyuxn_edge( n_edge),ierr)
     SLL_ALLOCATE(dyuyn_edge( n_edge),ierr)

     SLL_ALLOCATE(second_term( n_points),ierr)    
     SLL_ALLOCATE(phi( n_points),ierr)           
     SLL_ALLOCATE(phi_interm( n_points),ierr)    
     SLL_ALLOCATE(phi_edge_interm( n_points),ierr)    
     SLL_ALLOCATE(phi_center_interm( n_points),ierr)          
     SLL_ALLOCATE(phi_center( n_triangle),ierr)  
     SLL_ALLOCATE(second_term_center( n_points),ierr)    
     SLL_ALLOCATE(phi_edge(n_edge ),ierr)  
     SLL_ALLOCATE(second_term_edge( n_points),ierr)        
     SLL_ALLOCATE(matrix_poisson( n_points,1 + 4*num_cells + 2 ) , ierr)
     SLL_ALLOCATE(matrix_poisson_center( n_triangle,1 + 4*num_cells + 2), ierr)
     SLL_ALLOCATE(matrix_poisson_edge( n_edge,1 + 4*num_cells + 2 ) , ierr)
     SLL_ALLOCATE(l( n_points,1 + 4*num_cells + 2 ) , ierr)
     SLL_ALLOCATE(u( n_points,1 + 4*num_cells + 2), ierr)
     SLL_ALLOCATE(u_center( n_triangle,1 + 4*num_cells + 2), ierr)
     SLL_ALLOCATE(l_center( n_triangle,1 + 4*num_cells + 2 ) , ierr)
     SLL_ALLOCATE(l_edge( n_edge,1 + 4*num_cells + 2 ) , ierr)
     SLL_ALLOCATE(u_edge( n_edge,1 + 4*num_cells + 2), ierr)
     if  ( num_method == 15 )  then
        n_points2   = 1 + 6 * num_cells * (2*num_cells + 1)  
        SLL_ALLOCATE(second_term2( n_points2 ),ierr) 
        SLL_ALLOCATE(rho2( n_points2),ierr)     
        SLL_ALLOCATE(phi2( n_points2),ierr)           
        SLL_ALLOCATE(phi2_interm( n_points2),ierr)  
        SLL_ALLOCATE(uxn2( n_points2),ierr)
        SLL_ALLOCATE(uyn2( n_points2 ),ierr)
        SLL_ALLOCATE(dxuxn2( n_points),ierr)
        SLL_ALLOCATE(dxuyn2( n_points ),ierr)
        SLL_ALLOCATE(dxuxn_center2( n_triangle),ierr)
        SLL_ALLOCATE(dxuyn_center2( n_triangle),ierr)
        SLL_ALLOCATE(dxuxn_edge2( n_edge),ierr)
        SLL_ALLOCATE(dxuyn_edge2( n_edge),ierr)
        SLL_ALLOCATE(dyuxn2( n_points),ierr)
        SLL_ALLOCATE(dyuyn2( n_points ),ierr)
        SLL_ALLOCATE(dyuxn_center2( n_triangle),ierr)
        SLL_ALLOCATE(dyuyn_center2( n_triangle),ierr)
        SLL_ALLOCATE(dyuxn_edge2( n_edge),ierr)
        SLL_ALLOCATE(dyuyn_edge2( n_edge),ierr)
     endif


     call cpu_time(t_init)
     !*********************************************************
     !  Mesh , distribution function & density initialization   
     !*********************************************************

     mesh => new_hex_mesh_2d( num_cells, center_mesh_x1, center_mesh_x2,&
          radius=radius, EXTRA_TABLES = EXTRA_TABLES ) 

     if  ( num_method == 15 ) then
        mesh2 => new_hex_mesh_2d( num_cells*2, center_mesh_x1, center_mesh_x2,&
             radius=radius, EXTRA_TABLES = EXTRA_TABLES ) 
     endif
        
     
     det = (mesh%r1_x1*mesh%r2_x2 - mesh%r1_x2*mesh%r2_x1)/mesh%delta

     r11 = + mesh%r2_x2/det
     r12 = - mesh%r2_x1/det
     r21 = - mesh%r1_x2/det
     r22 = + mesh%r1_x1/det
     
     call init_distr(rho_tn,rho_center_tn,rho_edge_tn,num_method,mesh)

     if ( num_method /= 15 ) then
        call hex_matrix_poisson( matrix_poisson, mesh,1)
        call factolub_bande(matrix_poisson,l,u,n_points,l1,l2)
        call hex_second_terme_poisson( second_term, mesh, rho_tn )
        call solvlub_bande(l,u,phi_interm,second_term,n_points,l1,l2)

        do i = 1, mesh%num_pts_tot    ! need to re-index phi : 
           k1 = mesh%hex_coord(1, i)
           k2 = mesh%hex_coord(2, i)
           call index_hex_to_global(mesh, k1, k2, index_tab)
           phi(i) = phi_interm(index_tab)
        enddo
        call compute_hex_fields(mesh,uxn,uyn,dxuxn,dyuxn,dxuyn,dyuyn,phi,type=1)
     endif
     
     if ( num_method == 10 ) then
        do i = 1,mesh%num_triangles
           x = mesh%center_cartesian_coord(1,i)
           y = mesh%center_cartesian_coord(2,i)
           call get_cell_vertices_index( x, y, mesh, i1, i2, i3 )
           uxn_center(i) = (uxn(i1)+uxn(i2)+uxn(i3))/3._f64
           uyn_center(i) = (uyn(i1)+uyn(i2)+uyn(i3))/3._f64
        enddo
     endif

     if ( num_method == 15 ) then
        call hex_matrix_poisson( matrix_poisson_edge, mesh2,1)
        call factolub_bande(matrix_poisson_edge,l_edge,u_edge,n_points2,2*l1,2*l2)

        call assemble(rho_tn,rho_edge_tn,rho2,mesh,mesh2)
        call hex_second_terme_poisson( second_term2, mesh2, rho2)
        call solvlub_bande(l_edge,u_edge,phi2_interm,second_term2,n_points2,l1,l2)

        do i = 1, n_points2   ! need to re-index phi : 
           k1 = mesh%hex_coord(1, i)
           k2 = mesh%hex_coord(2, i)
           call index_hex_to_global(mesh, k1, k2, index_tab)
           phi2(i) = phi2_interm(index_tab)
        enddo
        call compute_hex_fields(mesh2,uxn2,uyn2,dxuxn2,dyuxn2,dxuyn2,dyuyn2,phi2,type=1)
        call deassemble(uxn,uxn_edge,uxn2,mesh,mesh2)
        call deassemble(uyn,uyn_edge,uyn2,mesh,mesh2)
     endif

       
     call hex_diagnostics(rho_tn,rho_center_tn,rho_edge_tn,t,mesh,uxn,uyn,nloops, num_method)


     !*********************************************************
     !                          Time loop
     !*********************************************************

     call cpu_time(t3)

     print*,"fin init",t3 - t_init

     do while (t .lt. tmax)

        t = t + dt
        nloops = nloops + 1
        count = count + 1


        !*********************************************************
        !                     interpolation
        !*********************************************************

        call  der_finite_difference( rho_tn, p, step, mesh, deriv )

        do i=1, n_points 

           x = mesh%cartesian_coord(1,i)
           y = mesh%cartesian_coord(2,i)

           !*************************************************
           !       computation of the characteristics
           !*************************************************

           call compute_characteristic_euler_2d_hex( &
                x,y,uxn,uyn,i,xx,yy,dt )

           inside = .true.
           h1 =  xx*r11 + yy*r12
           h2 =  xx*r21 + yy*r22 

           if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
           if ( abs(xx) > (radius-mesh%delta)*sqrt(3._f64)*0.5_f64) inside = .false.

           if ( inside ) then
              call hermite_interpolation(i, xx, yy, rho_tn, rho_center_tn,&
                   rho_edge_tn, rho_tn1, mesh, deriv, aire,& 
                   num_method)
           else 
              rho_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif

        end do ! end of the computation of the mesh points

        if ( num_method == 10 ) then 
           do i = 1,mesh%num_triangles
              x = mesh%center_cartesian_coord(1,i)
              y = mesh%center_cartesian_coord(2,i)
              call compute_characteristic_euler_2d_hex( &
                   x,y,uxn_center,uyn_center,i,xx,yy,dt )


              inside = .true.
              h1 =  xx*r11 + yy*r12
              h2 =  xx*r21 + yy*r22 

              if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
              if ( abs(xx) > (radius-mesh%delta)*sqrt(3._f64)*0.5_f64) inside = .false.
              
              if ( inside ) then
                 call hermite_interpolation(i, x, y, rho_tn, rho_center_tn,&
                      rho_edge_tn, rho_center_tn1, mesh, deriv, aire,& 
                      num_method)
              else 
                 rho_center_tn1(i) = 0._f64 ! dirichlet boundary condition
              endif

           enddo! end of the computation of the center of the triangles
        endif


        if ( num_method == 15 ) then 

           ! besoin de la séparation ici pour avoir rho_tn et rho_edge_tn

           do i=1, n_edge

              x = mesh%edge_center_cartesian_coord(1,i)
              y = mesh%edge_center_cartesian_coord(2,i) 

              call compute_characteristic_euler_2d_hex( &
                   x,y,uxn_edge,uyn_edge,i,xx,yy,dt )

              inside = .true.
              h1 =  xx*r11 + yy*r12
              h2 =  xx*r21 + yy*r22 


              if ( abs(h1) >  radius-mesh%delta .or. abs(h2) >  radius-mesh%delta ) inside = .false.
              if ( abs(xx) > (radius-mesh%delta)*sqrt(3._f64)*0.5_f64) inside = .false.

              if ( inside ) then
                 call hermite_interpolation(i,xx,yy,rho_tn,rho_center_tn&
                      ,rho_edge_tn, rho_edge_tn1,mesh,deriv,aire,&
                      num_method)
              else 
                 rho_edge_tn1(i) = 0._f64 ! dirichlet boundary condition
              endif

           enddo ! end of the computation at the middles of the edges

        endif


        !*********************************************************
        !      computing the solution of the poisson equation
        !*********************************************************

        
        call hex_second_terme_poisson( second_term, mesh, rho_tn )

        call solvlub_bande(l,u,phi_interm,second_term,n_points,l1,l2)

        do i = 1, mesh%num_pts_tot    ! need to re-index phi : 
           k1 = mesh%hex_coord(1, i)
           k2 = mesh%hex_coord(2, i)
           call index_hex_to_global(mesh, k1, k2, index_tab)
           phi(i) = phi_interm(index_tab)
        enddo

        call compute_hex_fields(mesh,uxn,uyn,dxuxn,dyuxn,dxuyn,dyuyn,phi,type=1)

        if ( num_method == 10 ) then
           do i = 1,mesh%num_triangles
              x = mesh%center_cartesian_coord(1,i)
              y = mesh%center_cartesian_coord(2,i)
              call get_cell_vertices_index( x, y, mesh, i1, i2, i3 )
              uxn_center(i) = (uxn(i1)+uxn(i2)+uxn(i3))/3._f64
              uyn_center(i) = (uyn(i1)+uyn(i2)+uyn(i3))/3._f64
           enddo
        endif
        
        if ( num_method == 15 ) then
           call assemble(rho_tn,rho_edge_tn,rho2,mesh,mesh2)
           call hex_second_terme_poisson( second_term2, mesh2, rho2)
           call solvlub_bande(l_edge,u_edge,phi2_interm,second_term2,n_points2,l1,l2)

           do i = 1, n_points2  ! need to re-index phi : 
              k1 = mesh%hex_coord(1, i)
              k2 = mesh%hex_coord(2, i)
              call index_hex_to_global(mesh, k1, k2, index_tab)
              phi2(i) = phi2_interm(index_tab)
           enddo
           call compute_hex_fields(mesh2,uxn2,uyn2,dxuxn2,dyuxn2,dxuyn2,dyuyn2,phi2,type=1)

           call deassemble(uxn,uxn_edge,uxn2,mesh,mesh2)
           call deassemble(uyn,uyn_edge,uyn2,mesh,mesh2)
        endif


        !*********************************************************
        !                  writing diagostics
        !*********************************************************

        call hex_diagnostics(rho_tn,rho_center_tn,rho_edge_tn,t,mesh,uxn,uyn,nloops, num_method)
        if (count == 10.and.nloops<10000) then
           call int2string(nloops,filenum)
           filename  = "center_guide_rho"//trim(filenum)
           call write_field_hex_mesh_xmf(mesh, rho_tn1, trim(filename))
           filename  = "center_guide_phi"//trim(filenum)
           call write_field_hex_mesh_xmf(mesh, phi, trim(filename))
           count = 0
            if ( num_method == 10 ) then

               filename  = "center_guide_rho_center"//trim(filenum)//".dat"
               call write_center(mesh,rho_center_tn1,filename)
            endif
           if ( num_method == 15 ) then 
              filename  = "center_guide_rho_edge"//trim(filenum)
          endif
        endif

        rho_tn = rho_tn1
        rho_center_tn = rho_center_tn1
        rho_edge_tn  = rho_edge_tn1
        
     enddo

     SLL_DEALLOCATE_ARRAY(second_term,ierr)        
     SLL_DEALLOCATE_ARRAY(phi,ierr)        
     SLL_DEALLOCATE_ARRAY(phi_interm,ierr)  
     SLL_DEALLOCATE_ARRAY(phi_edge_interm,ierr) 
     SLL_DEALLOCATE_ARRAY(phi_center_interm,ierr) 
     SLL_DEALLOCATE_ARRAY(second_term_center,ierr)    
     SLL_DEALLOCATE_ARRAY(phi_center,ierr) 
     SLL_DEALLOCATE_ARRAY(second_term_edge,ierr)    
     SLL_DEALLOCATE_ARRAY(phi_edge,ierr)      
     SLL_DEALLOCATE_ARRAY(rho_tn,ierr)
     SLL_DEALLOCATE_ARRAY(rho_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(rho_center_tn,ierr)
     SLL_DEALLOCATE_ARRAY(rho_center_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(rho_edge_tn,ierr)
     SLL_DEALLOCATE_ARRAY(rho_edge_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(uxn,ierr)
     SLL_DEALLOCATE_ARRAY(uyn,ierr)
     SLL_DEALLOCATE_ARRAY(uxn_center,ierr)
     SLL_DEALLOCATE_ARRAY(uyn_center,ierr)
     SLL_DEALLOCATE_ARRAY(uxn_edge,ierr)
     SLL_DEALLOCATE_ARRAY(uyn_edge,ierr)
     deallocate(deriv,matrix_poisson,matrix_poisson_center,matrix_poisson_edge )
     deallocate(l,u,l_center,l_edge,u_center,u_edge)

     call delete_hex_mesh_2d( mesh )

     call cpu_time(t_end)
     print*, "time used =", t_end - t_init

  end do

  
contains

  !*********initialization**************

  subroutine init_distr(f_tn,center_values_tn,edge_values_tn,num_method,mesh)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64,dimension(:)        :: f_tn, center_values_tn, edge_values_tn
    sll_int32  :: num_method
    sll_real64 :: x, y, epsilon = 0.1_f64
    sll_real64 :: rho
    sll_real64 :: r
    sll_int32  :: i

    do i = 1,mesh%num_pts_tot
       x = mesh%cartesian_coord(1,i)
       y = mesh%cartesian_coord(2,i)
       r = sqrt( x**2 + y**2 )

       if ( r <= 8._f64  .and. r >= 5._f64 ) then
          f_tn(i) = (1._f64 + epsilon * cos( 9._f64 * atan2(y,x)) )*&
               exp( -4._f64*(r-6.5_f64)**2)
       else
          f_tn(i) = 0._f64
       endif
    enddo

    if ( num_method == 10 ) then 

       do i = 1,mesh%num_triangles
          x = mesh%center_cartesian_coord(1,i)
          y = mesh%center_cartesian_coord(2,i)
          r = sqrt( x**2 + y**2 )

          if ( r <= 8._f64  .and. r >= 5._f64 ) then
             center_values_tn(i) = (1._f64 + epsilon * cos( 9._f64 * atan2(y,x)) )*&
                  exp( -4._f64*(r-6.5_f64)**2)
          else
             center_values_tn(i) = 0._f64
          endif

       enddo

    endif
    if ( num_method == 15 ) then 

       do i = 1,mesh%num_edges
          x = mesh%edge_center_cartesian_coord(1,i)
          y = mesh%edge_center_cartesian_coord(2,i)
          r = sqrt( x**2 + y**2 )

          if ( r <= 8._f64  .and. r >= 5._f64 ) then
             edge_values_tn(i) = (1._f64 + epsilon *cos( 9._f64 * atan2(y,x)))*&
                  exp( -4._f64*(r-6.5_f64)**2)
          else
             edge_values_tn(i) = 0._f64
          endif

       enddo

    endif


  end subroutine init_distr
  
  !********** diagnostics **************

  subroutine hex_diagnostics(rho,rho_center,rho_edge,t,mesh,uxn,uyn,nloop, num_method)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64,dimension(:) :: rho,rho_center,rho_edge,uxn,uyn
    sll_real64,intent(in)   :: t
    sll_int32 ,intent(in)   :: nloop,num_method
    sll_real64              :: mass,rho_min,norm_l1,norm_l2,norm_linf,energy
    sll_int32               :: i
    character(len = 50)     :: filename
    character(len = 4)      :: filenum

    energy    = 0._f64
    mass      = 0._f64
    rho_min   = rho(1)
    norm_l1   = 0._f64
    norm_l2   = 0._f64
    norm_linf = 0._f64

    call int2string(mesh%num_cells,filenum)
    filename  = "hex_diag_"//trim(filenum)//".dat"
    if (nloop == 0) then
       open(unit = 11, file=filename, action="write", status="replace")
    else
       open(unit = 11, file=filename, action="write", status="old",position = "append")
    endif

    do i = 1,mesh%num_pts_tot
       mass = mass + rho(i)
       norm_l1 = norm_l1 + abs(rho(i))
       norm_l2 = norm_l2 + rho(i)**2
       energy = energy + uxn(i)**2 + uyn(i)**2
       if ( abs(rho(i)) > norm_linf ) norm_linf = rho(i)
       if ( rho(i) < rho_min  ) rho_min  = rho(i)
    enddo

    print*,"diagnostic for t = ",t

     if ( num_method == 10 ) then 
        do i = 1,mesh%num_triangles
           !       mass = mass + rho(i)
           !       norm_l1 = norm_l1 + abs(rho(i))
           !       norm_l2 = norm_l2 + rho(i)**2
           !       energy = energy + uxn(i)**2 + uyn(i)**2
           if ( abs(rho_center(i)) > norm_linf ) norm_linf = rho_center(i)
           if ( rho_center(i) < rho_min  ) rho_min  = rho_center(i)
        enddo
     endif
     
     if ( num_method == 15 ) then 
        do i = 1,mesh%num_edges
           !       mass = mass + rho(i)
           !       norm_l1 = norm_l1 + abs(rho(i))
           !       norm_l2 = norm_l2 + rho(i)**2
           !       energy = energy + uxn(i)**2 + uyn(i)**2
           if ( abs(rho_edge(i)) > norm_linf ) norm_linf = rho_edge(i)
           if ( rho_edge(i) < rho_min  ) rho_min  = rho_edge(i)
        enddo
     endif

    energy  = sqrt(energy * mesh%delta**2)
    mass    = mass * mesh%delta**2
    norm_l1 = norm_l1 * mesh%delta**2
    norm_l2 = sqrt(norm_l2 * mesh%delta**2)

    write(11,*) t,mass,rho_min,norm_l1,norm_l2,norm_linf,energy

    close(11)

  end subroutine hex_diagnostics


  subroutine assemble(rho,rho_edge,rho2,mesh,mesh2)
    type(sll_hex_mesh_2d), pointer :: mesh,mesh2
    sll_real64,dimension(:) :: rho,rho_edge,rho2
    sll_int32               :: i,j,k,i1,i2

    rho2(1) = rho(1)

    rho2(2) = rho_edge(3)
    rho2(3) = rho_edge(2)
    rho2(4) = rho_edge(1)
    rho2(5) = rho_edge(15)
    rho2(6) = rho_edge(17)
    rho2(7) = rho_edge(19)

    ! à optimiser quand ce sera validé

    do i = 1,mesh%num_cells!mesh%num_pts_tot 
       do k = 1,6*i
          j  = 6*i*(2*i-1)+2*k
          i1 = 1+(i-1)*6+k
          rho2(j)   = rho(i1)

       enddo

       do k = 1,i 
          
          ! first edge
          j  = 6*i*(2*i-1)+2*k
          i1 = 1+3*i*(i-1)+k
          i2 = mesh%edge_center_index(1,i1)
          rho2(j+1) = rho_edge(i2)
          ! second edge
          j  = 6*i*(2*i-1)+2*(k+i)
          i1 = 1+3*i*(i-1)+k+  i+1
          i2 = mesh%edge_center_index(3,i1)
          rho2(j+1) = rho_edge(i2)
          ! third edge
          j  = 6*i*(2*i-1)+2*(k+2*i)
          i1 = 1+3*i*(i-1)+k+2*i+1
          i2 = mesh%edge_center_index(2,i1)
          rho2(j+1) = rho_edge(i2)
          ! fourth edge
          j  = 6*i*(2*i-1)+2*(k+3*i)
          i1 = 1+3*i*(i-1)+k+3*i+1
          i2 = mesh%edge_center_index(1,i1)
          rho2(j+1) = rho_edge(i2)
          ! fifth edge
          j  = 6*i*(2*i-1)+2*(k+4*i)
          i1 = 1+3*i*(i-1)+k+4*i
          i2 = mesh%edge_center_index(3,i1)
          rho2(j+1) = rho_edge(i2)
          ! sixth edge
          j  = 6*i*(2*i-1)+2*(k+5*i)
          i1 = 1+3*i*(i-1)+k+5*i
          i2 = mesh%edge_center_index(2,i1)
          rho2(j+1) = rho_edge(i2)

       enddo
       
       if (i<mesh%num_cells) then

          !corners

          j  = 6*i*(2*i-1)+2
          i2 = mesh%edge_center_index(3,i1)
          rho2(j) = rho_edge(i2)

          j  = 6*i*(2*i-1)+2+i
          i2 = mesh%edge_center_index(3,i1)
          rho2(j) = rho_edge(i2)


          do k = 2,i-1 

             ! first edge
             j  = 6*i*(2*i-1)+2*k
             i1 = 1+3*i*(i-1)+k
             i2 = mesh%edge_center_index(3,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(2,i1)
             rho2(j+1) = rho_edge(i2)

             ! second edge
             j  = 6*i*(2*i-1)+(k+i)
             i1 = 1+3*i*(i-1)+k+  i+1
             i2 = mesh%edge_center_index(2,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(1,i1)
             rho2(j+1) = rho_edge(i2)
             ! third edge
             j  = 6*i*(2*i-1)+(k+2*i)
             i1 = 1+3*i*(i-1)+k+2*i+1
             i2 = mesh%edge_center_index(1,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(3,i1)
             rho2(j+1) = rho_edge(i2)
             ! fourth edge
             j  = 6*i*(2*i-1)+(k+3*i)
             i1 = 1+3*i*(i-1)+k+3*i+1
             i2 = mesh%edge_center_index(2,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(3,i1)
             rho2(j+1) = rho_edge(i2)
             ! fifth edge
             j  = 6*i*(2*i-1)+(k+4*i)
             i1 = 1+3*i*(i-1)+k+4*i
             i2 = mesh%edge_center_index(1,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(2,i1)
             rho2(j+1) = rho_edge(i2)
             ! sixth edge
             j  = 6*i*(2*i-1)+2*(k+5*i)
             i1 = 1+3*i*(i-1)+k+5*i
             i2 = mesh%edge_center_index(1,i1)
             rho2(j) = rho_edge(i2)
             i2 = mesh%edge_center_index(3,i1)
             rho2(j+1) = rho_edge(i2)

          enddo
       endif


    enddo


  end subroutine assemble

  subroutine deassemble(rho,rho_edge,rho2,mesh,mesh2)
    type(sll_hex_mesh_2d), pointer :: mesh,mesh2
    sll_real64,dimension(:) :: rho,rho_edge,rho2
    sll_int32               :: i,j,k,i1,i2,j1,j2

    rho(1) = rho2(1)

    do i = 1,mesh%num_cells!mesh2%num_pts_tot 
       i1 = 2*i
       i2 = 2*i + 1 
       ! cas i impair (hexagone de pt milieu)
       do k = 1,i2 
          rho_edge(j) = rho(j1)
       enddo

       ! cas i pair (hexagone mix)
       do k = 1,i 
          rho(j1) = rho2(j)
          rho_edge(j2) = rho(j+1)

       enddo
    enddo


  end subroutine deassemble

  subroutine write_center(mesh,rho,name)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64,dimension(:)        :: rho
    sll_int32                      :: i
    character(len = 50)            :: name
    sll_real64                     :: x,y

    open (unit = 11,file=name, action="write", status="replace")

    do i = 1,mesh%num_triangles
       x = mesh%center_cartesian_coord(1,i)
       y = mesh%center_cartesian_coord(2,i)
       write(11,*) x,y,rho(i)
    enddo

    close (11)

  end subroutine write_center


end program
