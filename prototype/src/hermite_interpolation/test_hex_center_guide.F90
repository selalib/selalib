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

  type(sll_hex_mesh_2d), pointer          :: mesh
  sll_real64, dimension(:,:), allocatable :: deriv
  sll_real64, dimension(:),allocatable    :: second_term, phi,uxn,uyn,phi_interm 
  sll_real64, dimension(:),allocatable    :: second_term_center, phi_center
  sll_real64, dimension(:),allocatable    :: second_term_edge, phi_edge
  sll_real64, dimension(:),allocatable    :: uxn_center,uyn_center
  sll_real64, dimension(:),allocatable    :: uxn_edge,uyn_edge
  sll_real64, dimension(:,:) ,allocatable :: matrix_poisson, l, u
  sll_real64, dimension(:,:) ,allocatable :: matrix_poisson_center, l_center, u_center
  sll_real64, dimension(:,:) ,allocatable :: matrix_poisson_edge, l_edge, u_edge
  sll_int32                               :: i,j, k1, k2, index_tab
  sll_int32                               :: l1,l2, type  
  sll_int32    :: num_cells, n_points, n_triangle, n_edge
  sll_real64   :: center_mesh_x1, center_mesh_x2, radius
  sll_int32    :: nloops,count, ierr, EXTRA_TABLES = 1 ! put 1 for num_method = 15
  sll_real64   :: dt
  sll_real64   :: tmax
  sll_real64   :: t
  sll_real64   :: t_init, t_end
  sll_real64   :: step , aire, h1, h2, f_min, x ,y,xx, yy
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
  
  do num_cells = 40,40,40 

     !*********************************************************
     !             allocation
     !*********************************************************

     n_points   = 1 + 3 * num_cells * (num_cells + 1)  
     n_triangle = 6 * num_cells * num_cells
     n_edge    = 3 * num_cells * ( 3 * num_cells + 1 )

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
     SLL_ALLOCATE(second_term( n_points),ierr)    
     SLL_ALLOCATE(phi( n_points),ierr)           
     SLL_ALLOCATE(phi_interm( n_points),ierr)          
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


     !*********************************************************
     !  Mesh , distribution function & density initialization   
     !*********************************************************

     mesh => new_hex_mesh_2d( num_cells, center_mesh_x1, center_mesh_x2,&
          radius=radius, EXTRA_TABLES = EXTRA_TABLES ) 

     call init_distr(rho_tn,rho_center_tn,rho_edge_tn,num_method,mesh)

     call hex_matrix_poisson( matrix_poisson, mesh)
     ! call hex_matrix_poisson( matrix_poisson_center, mesh,type  )
     ! call hex_matrix_poisson( matrix_poisson_edge, mesh,type  )

     call factolub_bande(matrix_poisson,l,u,n_points,l1,l2)
     !call factolub_bande(matrix_poisson_center,l_center,u_center,n_points,l1,l2)
     !call factolub_bande(matrix_poisson_edge,l_edge,u_edge,n_points,l1,l2)


     tmax  = 100._f64
     dt    = 0.01_f64!*20._f64 !/ real(num_cells,f64)  
     t     = 0._f64
     !cfl   = radius * dt / ( radius / real(num_cells,f64)  )

     !*********************************************************
     !                          Time loop
     !*********************************************************

     nloops = 0
     count = 1

     call cpu_time(t_init)

     do while (t .lt. tmax)

        t = t + dt
        nloops = nloops + 1
        count = count + 1

        !*********************************************************
        !      computing the solution of the poisson equation
        !*********************************************************

        call hex_second_terme_poisson( second_term, mesh, rho_tn )

        call solvlub_bande(l,u,phi_interm,second_term,n_points,l1,l2)

        ! need to re-index phi : 

        do i = 1, mesh%num_pts_tot 
           k1 = mesh%hex_coord(1, i)
           k2 = mesh%hex_coord(2, i)
           call index_hex_to_global(mesh, k1, k2, index_tab)
           phi(i) = phi_interm(index_tab)
        enddo

        !call gauss_seidel_bande(matrix_poisson,phi,second_term,l1,l2,n_points,type)

        ! if ( num_method == 10 ) then
        !    call hex_second_terme_poisson( second_term_center, mesh, rho_center,type )
        !    call gauss_seidel_bande(matrix_poisson_center,phi_center,&
        !         second_term_center,l1,l2,n_triangle)
        ! endif

        ! if ( num_method == 15 ) then
        !    call hex_second_terme_poisson( second_term_edge, mesh, rho_edge,type )
        !    call gauss_seidel_bande(matrix_poisson_edge,phi_edge,second_term_edge,&
        !         l1,l2,n_edge)
        ! endif

        !*********************************************************
        !                  computing the fields
        !*********************************************************

        call compute_hex_fields(mesh,uxn,uyn,phi,type=1)
        ! if ( num_method == 10 ) call compute_hex_fields_center(mesh,uxn_center,uyn_center,phi_center)
        ! if ( num_method == 15 ) call compute_hex_fields_edge(mesh,uxn_edge,uyn_edge,phi_edge)


        !*********************************************************
        !                     interpolation
        !*********************************************************

        call  der_finite_difference( rho_tn, p, step, mesh, deriv )

        do i=1, n_points 

           x = mesh%cartesian_coord(1,i)
           y = mesh%cartesian_coord(2,i)

           ! xx = x*cos(dt) - y*sin(dt);
           ! yy = x*sin(dt) + y*cos(dt);

           xx = x - dt*y
           yy = y + dt*x

           !print*,xx,x*cos(dt) - y*sin(dt),yy,x*sin(dt) + y*cos(dt)

           !*************************************************
           !       computation of the characteristics
           !*************************************************

           ! call compute_characteristic_euler_2d_hex( &
           !      x,y,uxn,uyn,i,xx,yy,dt )
           inside = .true.
           h1 =  xx/sqrt(3.0_f64) + yy
           h2 = -xx/sqrt(3.0_f64) + yy 

           if ( abs(h1) >  radius-1e-15 .or. abs(h2) >  radius-1e-15 ) inside = .false.
           if ( abs(xx) > (radius-1e-15)*sqrt(3._f64)*0.5_f64) inside = .false.

           if ( inside ) then
              call hermite_interpolation(i, xx, yy, rho_tn, rho_center_tn,&
                   rho_edge_tn, rho_tn1, mesh, deriv, aire,& 
                   num_method)
           else 
              rho_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif

        end do ! end of the computation of the mesh points

        !    if ( num_method == 10 ) then 

        !       do i=1, n_triangle

        !          x = mesh%center_cartesian_coord(1,i)
        !          y = mesh%center_cartesian_coord(2,i)

        !          call compute_characteristic_euler_2d_hex( &
        !               x,y,uxn_center,uyn_center,i,xx,yy,dt )

        !          inside = .true.

        !          h1 =  xx/sqrt(3.0_f64) + yy
        !          h2 = -xx/sqrt(3.0_f64) + yy 

        !          if ( h1 >  radius .or. h2 >  radius ) inside = .false.
        !          if ( h1 < -radius .or. h2 < -radius ) inside = .false.
        !          if ( xx  < -radius*sqrt(3._f64)*0.5_f64 .or. xx &
        !               > radius*sqrt(3._f64)*0.5_f64  ) inside = .false.

        !          if ( inside ) then
        !             call hermite_interpolation(i, xx, yy, f_tn, center_values_tn,&
        !                  edge_values_tn, center_values_tn1, mesh, deriv, aire,& 
        !                  num_method)
        !          else 
        !             center_values_tn1(i) = 0._f64 ! dirichlet boundary condition
        !          endif
        !       enddo! end of the computation of the center of the triangles

        !    endif


        !    if ( num_method == 15 ) then 

        !       do i=1, n_edge

        !          x = mesh%edge_center_cartesian_coord(1,i)
        !          y = mesh%edge_center_cartesian_coord(2,i) 

        !          call compute_characteristic_euler_2d_hex( &
        !               x,y,uxn_edge,uyn_edge,i,xx,yy,dt )
        !          inside = .true.
        !          h1 =  xx/sqrt(3.0_f64) + yy
        !          h2 = -xx/sqrt(3.0_f64) + yy 

        !          if ( h1 >  radius .or. h2 >  radius ) inside = .false.
        !          if ( h1 < -radius .or. h2 < -radius ) inside = .false.
        !          if ( xx  < -radius*sqrt(3._f64)*0.5_f64 .or. xx &
        !               > radius*sqrt(3._f64)*0.5_f64  ) inside = .false.

        !          if ( inside ) then
        !             call hermite_interpolation(i, xx, yy, f_tn, center_values_tn,&
        !                  edge_values_tn, edge_values_tn1, mesh, deriv, aire,& 
        !                  num_method)
        !          else 
        !             edge_values_tn1(i) = 0._f64 ! dirichlet boundary condition
        !          endif

        !       enddo ! end of the computation at the middles of the edges

        !    endif

        !    center_values_tn = center_values_tn1
        !    edge_values_tn  = edge_values_tn1
        if (count == 100) then
           call int2string(nloops,filenum)
           filename  = "center_guide_rho"//trim(filenum)
           call write_field_hex_mesh_xmf(mesh, rho_tn1, trim(filename))
           count = 0
        endif
        rho_tn = rho_tn1
        !call diagnostics(mesh,rho_tn)
     enddo

     SLL_DEALLOCATE_ARRAY(second_term,ierr)        
     SLL_DEALLOCATE_ARRAY(phi,ierr)        
     SLL_DEALLOCATE_ARRAY(phi_interm,ierr)  
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

  subroutine init_distr(f_tn,center_values_tn,edge_values_tn,num_method,mesh)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64,dimension(:)        :: f_tn, center_values_tn, edge_values_tn
    sll_int32  :: num_method
    sll_real64 :: x, y, epsilon = 0.0001_f64
    sll_real64 :: rho
    sll_real64 :: r
    sll_int32  :: i

    do i = 1,mesh%num_pts_tot
       x = mesh%cartesian_coord(1,i)
       y = mesh%cartesian_coord(2,i)
       r = sqrt( x**2 + y**2 )

       if ( r <= 8._f64  .and. r >= 5._f64 ) then
          f_tn(i) = (1._f64 + epsilon * cos( 7._f64 * atan2(y,x)) )*&
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
             center_values_tn(i) = (1._f64 + epsilon * cos( 7._f64 * atan2(y,x)) )*&
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
             edge_values_tn(i) = (1._f64 + epsilon *cos( 7._f64 * atan2(y,x)))*&
                  exp( -4._f64*(r-6.5_f64)**2)
          else
             edge_values_tn(i) = 0._f64
          endif

       enddo

    endif


  end subroutine init_distr
  
  subroutine compute_hex_fields(mesh,uxn,uyn,phi,type)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64,dimension(:)        :: uxn, uyn, phi
    sll_int32,          intent(in) :: type
    sll_int32  :: i,h1,h2
    sll_real64 :: phii_2, phii_1, phii1, phii2, phij_2, phij_1, phij1, phij2
    sll_real64 :: uh1, uh2
    sll_real64 :: v1, v2, v3, v4   


    ! v1 = mesh%r1_x1/mesh%delta
    ! v2 = mesh%r1_x2/mesh%delta
    ! v3 = mesh%r2_x1/mesh%delta
    ! v4 = mesh%r2_x2/mesh%delta

    
    if (type==1) then

       do i = 1,mesh%num_pts_tot

          ! h1 = mesh%hex_coord(1,i)
          ! h2 = mesh%hex_coord(2,i)

          ! phii_2 = value_if_inside_phi(h1-2,h2,mesh,phi)
          ! phii_1 = value_if_inside_phi(h1-1,h2,mesh,phi)
          ! phii1  = value_if_inside_phi(h1+1,h2,mesh,phi)
          ! phii2  = value_if_inside_phi(h1+2,h2,mesh,phi)

          ! phij_2 = value_if_inside_phi(h1,h2-2,mesh,phi)
          ! phij_1 = value_if_inside_phi(h1,h2-1,mesh,phi)
          ! phij1  = value_if_inside_phi(h1,h2+1,mesh,phi)
          ! phij2  = value_if_inside_phi(h1,h2+2,mesh,phi)

          ! order 2

          !uh1 = ( phii1 - phii_1 ) / (2._f64)!*mesh%delta)
          !uh2 = ( phij1 - phij_1 ) / (2._f64)!*mesh%delta)

          ! uxn = -(v2*uh1 + v4*uh2)   ! -d(phi)/dy 
          ! uyn = +(v1*uh1 + v3*uh2)   ! +d(phi)/dx

          ! uxn(i) = -( mesh%r1_x2*uh1 + mesh%r2_x2*uh2)   ! -d(phi)/dy 
          ! uyn(i) = +( mesh%r1_x1*uh1 + mesh%r2_x1*uh2)   ! +d(phi)/dx

          uxn(i) = + mesh%cartesian_coord(2,i)   ! +y
          uyn(i) = - mesh%cartesian_coord(1,i)   ! -x

          ! order 4
          ! uxn = - ( phi(nr1i_2) + 8._f64 * ( - phi(nr1i_1) + phi(nr1i1) ) &
          !      - phi(nr1i2) ) / (12._f64*mesh%delta)
          ! uyn = + ( phi(nr2i_2) + 8._f64 * ( - phi(nr2i_1) + phi(nr2i1) ) &
          !      - phi(nr2i2) ) / (12._f64*mesh%delta)

	  ! _Uy[ix][iy]   = +(_phi[ix+1][iy]-_phi[ix-1][iy])/(2.0*dx);

	  ! _dxUx[ix][iy] = -(_phi[ix+1][iy+1]-_phi[ix+1][iy-1]-_phi[ix-1][iy+1]+_phi[ix-1][iy-1])/(4*dx*dy);
	  ! _dyUx[ix][iy] = -(_phi[ix][iy+1]  -2.*_phi[ix][iy]                  +_phi[ix][iy-1]  )/(dy*dy);	

	  ! _dyUy[ix][iy] = +(_phi[ix+1][iy+1]-_phi[ix+1][iy-1]-_phi[ix-1][iy+1]+_phi[ix-1][iy-1])/(4*dx*dy);
	  ! _dxUy[ix][iy] = +(_phi[ix+1][iy]  -2.0*_phi[ix][iy]                 +_phi[ix-1][iy]  )/(dx*dx);
       end do
    endif


  end subroutine compute_hex_fields

  ! subroutine compute_hex_fields_center(mesh,uxn_center,uyn_center,phi_center)
  !   type(sll_hex_mesh_2d), pointer :: mesh
  !   sll_real64,dimension(:)        :: uxn_center, uyn_center, phi_center
  !   sll_int32  :: i,h1,h2
  !   sll_int32  :: nr1i_2,nr1i_1,nr1i1,nr1i2,nr2i_2,nr2i_1,nr2i1,nr2i2
    
  !   do i = 1,mesh%num_triangles

  !      h1 = mesh%hex_coord(1,i)
  !      h2 = mesh%hex_coord(2,i)

  !      nr1i_2 = hex_to_global(mesh,h1-2,h2)
  !      nr1i_1 = hex_to_global(mesh,h1-1,h2)
  !      nr1i1  = hex_to_global(mesh,h1+1,h2)
  !      nr1i2  = hex_to_global(mesh,h1+2,h2) 

  !      nr2i_2 = hex_to_global(mesh,h1,h2-2)
  !      nr2i_1 = hex_to_global(mesh,h1,h2-1)
  !      nr2i1  = hex_to_global(mesh,h1,h2+1)
  !      nr2i2  = hex_to_global(mesh,h1,h2+2) 

  !      ! order 2
  !      uxn_center(i) = - ( phi_center(nr1i1) - phi_center(nr1i_1))&
  !           / (2._f64*mesh%delta)
  !      uyn_center(i) = + ( phi_center(nr2i1) - phi_center(nr2i_1))&
  !           / (2._f64*mesh%delta)
  !      ! order 4
  !      ! uxn_center(i) = - ( phi_center(nr1i_2) + 8._f64 * &
  !      !( - phi_center(nr1i_1) + phi_center(nr1i1) ) - phi_center(nr1i2) ) &
  !      !     / (12._f64*mesh%delta)
  !      ! uyn_center(i) = + ( phi_center(nr2i_2) + 8._f64 * &
  !      !( - phi_center(nr2i_1) + phi_center(nr2i1) ) - phi_center(nr2i2) ) &
  !      !     / (12._f64*mesh%delta)


  !   enddo

  ! end subroutine compute_hex_fields_center

  ! subroutine compute_hex_fields_edge(mesh,uxn_edge,,uyn_edge,phi_edge)
  !   type(sll_hex_mesh_2d), pointer :: mesh
  !   sll_real64,dimension(:)        :: uxn, uyn, phi
  !   sll_int32  :: i,h1,h2
  !   sll_int32  :: nr1i_2,nr1i_1,nr1i1,nr1i2,nr2i_2,nr2i_1,nr2i1,nr2i2
    
  !   do i = 1,mesh%num_edges


  !      h1 = mesh%hex_coord(1,i)
  !      h2 = mesh%hex_coord(2,i)

  !      nr1i_2 = hex_to_global(mesh,h1-2,h2)
  !      nr1i_1 = hex_to_global(mesh,h1-1,h2)
  !      nr1i1  = hex_to_global(mesh,h1+1,h2)
  !      nr1i2  = hex_to_global(mesh,h1+2,h2) 

  !      nr2i_2 = hex_to_global(mesh,h1,h2-2)
  !      nr2i_1 = hex_to_global(mesh,h1,h2-1)
  !      nr2i1  = hex_to_global(mesh,h1,h2+1)
  !      nr2i2  = hex_to_global(mesh,h1,h2+2) 

  !      ! order 2
  !      uxn = - ( phi(nr1i1) - phi(nr1i_1))/ (2._f64*mesh%delta)
  !      uyn = + ( phi(nr2i1) - phi(nr2i_1))/ (2._f64*mesh%delta)
  !      ! order 4
  !      ! uxn = - ( phi(nr1i_2) + 8._f64 * ( - phi(nr1i_1) + phi(nr1i1) ) &
  !      !      - phi(nr1i2) ) / (12._f64*mesh%delta)
  !      ! uyn = + ( phi(nr2i_2) + 8._f64 * ( - phi(nr2i_1) + phi(nr2i1) ) &
  !      !      - phi(nr2i2) ) / (12._f64*mesh%delta)

  !   enddo

  ! end subroutine compute_hex_fields_edge

  function value_if_inside_phi(k1,k2,mesh,rho) result(f)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:)       :: rho
    sll_int32  :: k1, k2, n
    sll_real64 :: f 

    if ( abs(k1) > mesh%num_cells .or. abs(k2) > mesh%num_cells .or. &
         (k1)*(k2)< 0 .and. ( abs(k1) + abs(k2) > mesh%num_cells) ) then
       f = 0._f64 ! null dirichlet boundary condition
    else
       n = hex_to_global(mesh,k1,k2)
       f = rho(n)
    endif

  endfunction value_if_inside_phi

end program
