program test_hex_hermite

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_constants
  use euler_2d_hex
  use sll_hex_meshes
  use hex_poisson
  use pivotbande
  use sll_box_splines
  implicit none

  type(sll_hex_mesh_2d),   pointer        :: mesh
  type(sll_box_spline_2d), pointer        :: spline
  sll_real64, dimension(:),   allocatable :: dxuxn,dyuxn,dxuyn,dyuyn
  sll_real64, dimension(:),   allocatable :: second_term, phi,uxn,uyn,phi_interm
  sll_real64, dimension(:),   allocatable :: second_term2, phi2,uxn2,uyn2,phi2_interm,rho2
  sll_real64, dimension(:),   allocatable :: dxuxn2,dyuxn2,dxuyn2,dyuyn2
  sll_real64, dimension(:,:), allocatable :: matrix_poisson, l, u
  sll_real64, dimension(:),   allocatable :: rho_tn   ! distribution at time n
  sll_real64, dimension(:),   allocatable :: rho_tn1  ! distribution at time n + 1
  sll_real64, dimension(:),   allocatable :: x1_char
  sll_real64, dimension(:),   allocatable :: x2_char

  sll_int32    :: deg = 2
  sll_int32    :: i,j, k1, k2, index_tab
  sll_int32    :: l1,l2, type  
  sll_int32    :: i1,i2,i3
  sll_int32    :: num_cells, n_points, n_triangle, n_points2
  sll_real64   :: center_mesh_x1, center_mesh_x2, radius
  sll_int32    :: nloops,count, ierr, EXTRA_TABLES = 0
  sll_real64   :: dt
  sll_real64   :: tmax
  sll_real64   :: t
  sll_real64   :: t_init, t_end
  sll_real64   :: t1,t2,t3,t4,t5,t6
  sll_real64   :: step , aire, h1, h2, f_min, x ,y,xx, yy
  sll_real64   :: r11,r12,r21,r22,det
  logical      :: inside
  sll_int32    :: p = 6!-> degree of the approximation for the derivative 
  character(len = 50) :: filename
  character(len = 4)  :: filenum

  center_mesh_x1 = 0._f64
  center_mesh_x2 = 0._f64

  radius = 14._f64

  
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

     ! Mesh creation -------------------------
     mesh => new_hex_mesh_2d( num_cells, center_mesh_x1, center_mesh_x2,&
          radius=radius, EXTRA_TABLES = EXTRA_TABLES )         

     n_points   = mesh%num_pts_tot
     n_triangle = mesh%num_triangles

     det = (mesh%r1_x1*mesh%r2_x2 - mesh%r1_x2*mesh%r2_x1)/mesh%delta

     r11 = + mesh%r2_x2/det
     r12 = - mesh%r2_x1/det
     r21 = - mesh%r1_x2/det
     r22 = + mesh%r1_x1/det
     ! ---------------------------------------

     print*," minimum number of points computed on the mesh: ", n_points

     l1 = 2*num_cells+1
     l2 = l1

     step = radius / real(num_cells,f64)
     aire = step**2*sqrt(3._f64)*0.25_f64

     SLL_ALLOCATE(rho_tn( n_points),ierr)
     SLL_ALLOCATE(rho_tn1( n_points ),ierr)
     SLL_ALLOCATE(x1_char( n_points ),ierr)
     SLL_ALLOCATE(x2_char( n_points ),ierr)

     SLL_ALLOCATE(uxn( n_points),ierr)
     SLL_ALLOCATE(uyn( n_points ),ierr)

     SLL_ALLOCATE(dxuxn( n_points),ierr)
     SLL_ALLOCATE(dxuyn( n_points ),ierr)
     SLL_ALLOCATE(dyuxn( n_points),ierr)
     SLL_ALLOCATE(dyuyn( n_points ),ierr)

     SLL_ALLOCATE(second_term( n_points),ierr)    
     SLL_ALLOCATE(phi( n_points),ierr)           
     SLL_ALLOCATE(phi_interm( n_points),ierr)    
     SLL_ALLOCATE(matrix_poisson( n_points,1 + 4*num_cells + 2 ) , ierr)
     SLL_ALLOCATE(l( n_points,1 + 4*num_cells + 2 ) , ierr)
     SLL_ALLOCATE(u( n_points,1 + 4*num_cells + 2), ierr)

     call cpu_time(t_init)

     !*********************************************************
     !  Distribution function & density initialization   
     !*********************************************************

     ! Spline initialization -----------------
     spline => new_box_spline_2d(mesh, SLL_DIRICHLET)
     ! ---------------------------------------
     
     ! Initial distribution ------------------
     call init_distr(rho_tn,mesh)
     ! ---------------------------------------

     ! Poisson solver ------------------------
     call hex_matrix_poisson( matrix_poisson, mesh,1)
     call factolub_bande(matrix_poisson,l,u,n_points,l1,l2)
     call hex_second_terme_poisson( second_term, mesh, rho_tn )
     call solvlub_bande(l,u,phi_interm,second_term,n_points,l1,l2)

     do i = 1, mesh%num_pts_tot 
        k1 = mesh%hex_coord(1, i)
        k2 = mesh%hex_coord(2, i)
        call index_hex_to_global(mesh, k1, k2, index_tab)
        phi(i) = phi_interm(index_tab)
     enddo
     call compute_hex_fields(mesh,uxn,uyn,dxuxn,dyuxn,dxuyn,dyuyn,phi,type=1)
     ! ---------------------------------------     
      
     call hex_diagnostics(rho_tn,t,mesh,uxn,uyn,nloops)


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

        call compute_box_spline_2d( rho_tn, deg, spline )

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
              rho_tn1(i) = hex_interpolate_value(mesh, xx, yy, spline, deg)
           else 
              rho_tn1(i) = 0._f64 ! dirichlet boundary condition
           endif

        end do ! end of the computation of the mesh points


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


        !*********************************************************
        !                  writing diagostics
        !*********************************************************

        call hex_diagnostics(rho_tn,t,mesh,uxn,uyn,nloops)
        if (count == 10.and.nloops<10000) then
           call int2string(nloops,filenum)
           filename  = "center_guide_rho"//trim(filenum)
           call write_field_hex_mesh_xmf(mesh, rho_tn1, trim(filename))
           filename  = "center_guide_phi"//trim(filenum)
           call write_field_hex_mesh_xmf(mesh, phi, trim(filename))
           count = 0
        endif

        rho_tn = rho_tn1
        
     enddo

     SLL_DEALLOCATE_ARRAY(second_term,ierr)        
     SLL_DEALLOCATE_ARRAY(phi,ierr)        
     SLL_DEALLOCATE_ARRAY(phi_interm,ierr)  
     SLL_DEALLOCATE_ARRAY(rho_tn,ierr)
     SLL_DEALLOCATE_ARRAY(rho_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(uxn,ierr)
     SLL_DEALLOCATE_ARRAY(uyn,ierr)
     deallocate(matrix_poisson)
     deallocate(l)
     deallocate(u)

     call delete_hex_mesh_2d( mesh )

     call cpu_time(t_end)
     print*, "time used =", t_end - t_init

  end do

  
contains

  !*********initialization**************

  subroutine init_distr(f_tn,mesh)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:)       :: f_tn
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


  end subroutine init_distr
  
  !********** diagnostics **************

  subroutine hex_diagnostics(rho,t,mesh,uxn,uyn,nloop)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64,dimension(:) :: rho,uxn,uyn
    sll_real64,intent(in)   :: t
    sll_int32 ,intent(in)   :: nloop
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
