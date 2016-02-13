module sll_m_hex_poisson
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_hexagonal_meshes, only: &
    sll_f_new_hex_mesh_2d, &
    sll_t_hex_mesh_2d

  implicit none

  public :: &
    sll_s_compute_hex_fields, &
    sll_s_hex_matrix_poisson, &
    sll_s_hex_second_terme_poisson

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !***********************************************************
  ! self module is used to solve the Poison equation on a hexagonal mesh
  ! by creating the poisson matrix for the hex mesh : matrix_poisson
  ! and the second term of the linear equation      : second_term
  ! therefore we solve with a linear solver : matrix_poisson * X = second_term
  ! we also compute the fields and it derivatives from the potential we get from
  ! the Poisson equation
  ! Precision is of order 4 ( for both poisson solver and fields computation )
  ! for any question : prouveur@math.univ-lyon1.fr
  !***********************************************************

contains

  subroutine sll_s_hex_matrix_poisson( matrix_poisson, mesh,type)
    type(sll_t_hex_mesh_2d), pointer         :: mesh
    sll_real64, dimension(:,:), intent(out):: matrix_poisson

    sll_int32,                   intent(in):: type ! unused parameter atm

    sll_int32                              :: num_cells ! number of hex cells
    sll_int32                              :: global    ! global index

    ! index on the matrix
    sll_int32                              :: index_tab, index_tabij
    sll_int32                              :: index_tabi_1j, index_tabij_1
    sll_int32                              :: index_tabij1, index_tabi1j
    sll_int32                              :: index_tabi1j1, index_tabi_1j_1
    sll_int32                              :: k1, k2, n

    num_cells = mesh%num_cells

    matrix_poisson = 0._f64

    ! indexation of the matrix : we fill row by row with 
    ! j + (i-1)*n(j) the index given by index_tab
    ! therefore we get (i-1,j-1), in hex coordinate, then (i-1,j)
    ! (i,j-1) ; (i,j) ; (i,j+1)  
    ! (i+1,j) ; (i,j+1)

    if (type == 1) then
       n = mesh%num_pts_tot 
    elseif (type == 2) then
       n = mesh%num_triangles 
    elseif (type == 3) then
       n = mesh%num_edges
    endif

    do global = 1, n

       k1 = mesh%hex_coord(1, global)
       k2 = mesh%hex_coord(2, global)

       call mesh%index_hex_to_global(k1, k2, index_tab)

       call mesh%index_hex_to_global(k1-1, k2-1, index_tabi_1j_1)
       call mesh%index_hex_to_global(k1-1, k2  , index_tabi_1j)
       call mesh%index_hex_to_global(k1  , k2-1, index_tabij_1)
       call mesh%index_hex_to_global(k1  , k2+1, index_tabij1 )
       call mesh%index_hex_to_global(k1+1, k2  , index_tabi1j )
       call mesh%index_hex_to_global(k1+1, k2+1, index_tabi1j1 )

       index_tabij = index_tab

       index_tabi_1j_1 = index_tabi_1j_1 - index_tab + 2*num_cells + 2
       index_tabi_1j   = index_tabi_1j   - index_tab + 2*num_cells + 2
       index_tabij_1   = index_tabij_1   - index_tab + 2*num_cells + 2
       index_tabij     = index_tabij     - index_tab + 2*num_cells + 2
       index_tabij1    = index_tabij1    - index_tab + 2*num_cells + 2
       index_tabi1j    = index_tabi1j    - index_tab + 2*num_cells + 2
       index_tabi1j1   = index_tabi1j1   - index_tab + 2*num_cells + 2

       if ( abs(k1) == num_cells .and. abs(k2) == num_cells .or.&
            abs(k1) == num_cells .and. k2 == 0 .or.&
            abs(k2) == num_cells .and. k1 == 0 &
            ) then ! corners

          if (k1 == num_cells.and. k2 == 0) then! corner top right
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          elseif (k1 == num_cells.and. k2 == num_cells) then! corner top
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          else if (k2 == num_cells.and. k1 == 0) then! corner top left
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          else if (k1 == -num_cells.and. k2 == 0) then! corner bottom left
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          else if (-k1== num_cells.and.k2== -num_cells) then!corner bottom
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          else if (k2 == -num_cells.and. k1 == 0) then! corner bottom right
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          endif

       else if ( k1*k2<0 .and. abs(k1)+abs(k2)==num_cells .or.&
            abs(k1) == num_cells .or.  abs(k2) == num_cells&
            ) then ! edges

          if (k1 == num_cells) then! edge top right
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          elseif (k2 == num_cells) then! edge top left             
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          else if (k1*k2<0 .and. -k1+k2==num_cells) then! edge left 
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          else if (k1 == -num_cells) then! edge bottom left        
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          else if (k2 == -num_cells) then! edge bottom right      
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          else if (k1*k2<0 .and. k1-k2==num_cells) then! edge right 
             matrix_poisson(index_tab,index_tabij    ) =  1._f64
          endif

       elseif ( abs(k1) == num_cells-1 .and. abs(k2) == num_cells-1 .or.&
            abs(k1) == num_cells-1 .and. k2 == 0 .or.&
            abs(k2) == num_cells-1 .and. k1 == 0 &
            ) then ! corners

          if (k1 == num_cells-1.and. k2 == 0) then! corner top right
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
          elseif (k1 == num_cells-1.and. k2 == num_cells-1) then! corner top 
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
          else if (k2 == num_cells-1.and. k1 == 0) then! corner top left
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
          else if (k1 == -num_cells+1.and. k2 == 0) then! corner bottom left
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          else if (-k1== num_cells-1.and.k2== -num_cells+1) then!corner bottom
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          else if (k2 == -num_cells+1.and. k1 == 0) then! corner bottom right
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          endif

       else if ( k1*k2<0 .and. abs(k1)+abs(k2)==num_cells-1 .or.&
            abs(k1) == num_cells-1 .or.  abs(k2) == num_cells-1&
            ) then ! edges

          if (k1 == num_cells-1) then! edge top right
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
          elseif (k2 == num_cells-1) then! edge top left              
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
          else if (k1*k2<0 .and. -k1+k2==num_cells-1) then! edge left         
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          else if (k1 == -num_cells+1) then! edge bottom left        
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          else if (k2 == -num_cells+1) then! edge bottom right      
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          else if (k1*k2<0 .and. k1-k2==num_cells-1) then! edge right          
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          endif

       else ! general case

          matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
          matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
          matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
          matrix_poisson(index_tab,index_tabij    ) =  6._f64
          matrix_poisson(index_tab,index_tabij1   ) = -1._f64
          matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
          matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64

       endif
    enddo

  end subroutine sll_s_hex_matrix_poisson

  
  
  subroutine sll_s_hex_second_terme_poisson( second_terme, mesh, rho )
    type(sll_t_hex_mesh_2d), pointer        :: mesh
    sll_real64, dimension(:), intent(out) :: second_terme
    sll_real64, dimension(:), intent(in)  :: rho
    sll_int32                             :: num_cells, i, index_tab, k1, k2
    sll_real64                            :: step
    sll_real64                            :: f1,f2,f3,f4,f5,f6
    sll_int32                             :: global

    num_cells = mesh%num_cells
    step = mesh%delta**2 * 1.5_f64 

    !************************
    ! general case
    !************************
       
       do global = 1, mesh%num_pts_tot 

          k1 = mesh%hex_coord(1, global)
          k2 = mesh%hex_coord(2, global)

          call mesh%index_hex_to_global(k1, k2, index_tab)

          f1 = value_if_inside_rho(k1+1,k2  ,mesh,rho)
          f2 = value_if_inside_rho(k1+1,k2+1,mesh,rho)
          f3 = value_if_inside_rho(k1  ,k2+1,mesh,rho)
          f4 = value_if_inside_rho(k1-1,k2  ,mesh,rho)
          f5 = value_if_inside_rho(k1-1,k2-1,mesh,rho)
          f6 = value_if_inside_rho(k1  ,k2-1,mesh,rho)

          second_terme(index_tab) = ( 0.75_f64*rho(global) + & 
               (f1+f2+f3+f4+f5+f6)/24._f64 ) * step ! ordre 4

          !second_terme(index_tab) = rho(global) * step ! order 2

       enddo


       !************************
       ! Boundaries
       !************************

       ! corners of the hexagon

       call mesh%index_hex_to_global(num_cells , 0, index_tab)
       second_terme(index_tab) = 0._f64
       call mesh%index_hex_to_global(num_cells , num_cells, index_tab)
       second_terme(index_tab) = 0._f64
       call mesh%index_hex_to_global(0 , num_cells, index_tab)
       second_terme(index_tab) = 0._f64
       call mesh%index_hex_to_global(-num_cells , 0, index_tab)
       second_terme(index_tab) = 0._f64
       call mesh%index_hex_to_global(-num_cells , -num_cells, index_tab)
       second_terme(index_tab) = 0._f64
       call mesh%index_hex_to_global(0 , -num_cells, index_tab)
       second_terme(index_tab) = 0._f64


       ! edges of the hexagon  

       do i = 1,num_cells-1   !( 0 and num_cells  are the corners )

          ! top right edge
          call mesh%index_hex_to_global(num_cells, i, index_tab)
          second_terme(index_tab) = 0._f64
          ! top left edge
          call mesh%index_hex_to_global(num_cells - i, num_cells , index_tab)
          second_terme(index_tab) = 0._f64
          ! left edge
          call mesh%index_hex_to_global(- i, num_cells- i , index_tab)
          second_terme(index_tab) = 0._f64
          ! bottom left edge
          call mesh%index_hex_to_global(- num_cells, - i , index_tab)
          second_terme(index_tab) = 0._f64
          ! bottom right edge
          call mesh%index_hex_to_global(- i, - num_cells , index_tab)
          second_terme(index_tab) = 0._f64
          ! right edge
          call mesh%index_hex_to_global(num_cells - i, - i , index_tab)
          second_terme(index_tab) = 0._f64

       enddo

       !************************
       ! Boundary conditions      -> 0 everywhere at the moment 
       !************************

       ! corners of the hexagon

       call mesh%index_hex_to_global(num_cells-1 , 1, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       call mesh%index_hex_to_global(num_cells-1 , num_cells-1, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       call mesh%index_hex_to_global(1 , num_cells-1, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       call mesh%index_hex_to_global(-num_cells+1 , 1, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       call mesh%index_hex_to_global(-num_cells+1 , -num_cells+1, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       call mesh%index_hex_to_global(1 , -num_cells+1, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64


       ! edges of the hexagon

       do i = 2,num_cells-2   !( 1 and num_cells-1 are the corners )

          ! top right edge
          call mesh%index_hex_to_global(num_cells-1, i, index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64
          ! top left edge
          call mesh%index_hex_to_global(num_cells-1 - i, num_cells-1 , index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64
          ! left edge
          call mesh%index_hex_to_global(- i, num_cells-1 - i , index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64
          ! bottom left edge
          call mesh%index_hex_to_global(- num_cells+1, - i , index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64
          ! bottom right edge
          call mesh%index_hex_to_global(- i, - num_cells+1 , index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64
          ! right edge
          call mesh%index_hex_to_global(num_cells-1 - i, - i , index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64

       enddo


  end subroutine sll_s_hex_second_terme_poisson

  !> @details test to check if a point is inside the mesh or not
  function value_if_inside_rho(k1,k2,mesh,rho) result(f)
    type(sll_t_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:)       :: rho
    sll_int32  :: k1, k2, n
    sll_real64 :: f

    if ( abs(k1) > mesh%num_cells .or. abs(k2) > mesh%num_cells .or. &
         (k1)*(k2)< 0 .and. ( abs(k1) + abs(k2) > mesh%num_cells) ) then
       f = 0._f64 ! null dirichlet boundary condition
    else
       n = mesh%hex_to_global(k1,k2)
       f = rho(n)
    endif

  endfunction value_if_inside_rho

! subroutine to compute the fields and its derivatives from the results 
! of solving the poisson equation . Precision is of order 4

subroutine sll_s_compute_hex_fields(mesh,uxn,uyn,dxuxn,dyuxn,dxuyn,dyuyn,phi,type)
    type(sll_t_hex_mesh_2d), pointer :: mesh
    sll_real64,dimension(:)        :: uxn, uyn, phi,dxuxn,dyuxn,dxuyn,dyuyn
    sll_int32,          intent(in) :: type
    sll_int32  :: i,h1,h2
    sll_real64 :: r11,r12,r21,r22,det
    sll_real64 :: phi2j1,phi1j2,phi_2j1,phi_2j_1,phi_1j_2,phi1j_2,phi2j_1,phi_1j2
    sll_real64 :: phi_2j_2,phi2j_2,phi2j2,phi_2j2
    sll_real64 :: phi_3j,phij_3,phij,phi1j1,phi_1j_1,phi1j_1,phi_1j1
    sll_real64 :: phi_2j, phi_1j, phi1j, phi2j, phij_2, phij_1, phij1, phij2
    sll_real64 :: uh1, uh2, uh1h1, uh2h2, uh1h2

    det = (mesh%r1_x1*mesh%r2_x2 - mesh%r1_x2*mesh%r2_x1)/mesh%delta

    r11 = + mesh%r2_x2/det
    r12 = - mesh%r2_x1/det
    r21 = - mesh%r1_x2/det
    r22 = + mesh%r1_x1/det

    if (type==1) then

       do i = 1,mesh%num_pts_tot

          h1 = mesh%hex_coord(1,i)
          h2 = mesh%hex_coord(2,i)

          phi_1j_1 = value_if_inside_phi(h1-1,h2-1,mesh,phi)
          phi1j1   = value_if_inside_phi(h1+1,h2+1,mesh,phi)
          phi1j_1  = value_if_inside_phi(h1+1,h2-1,mesh,phi)
          phi_1j1  = value_if_inside_phi(h1-1,h2+1,mesh,phi)

          phi_2j_2 = value_if_inside_phi(h1-2,h2-2,mesh,phi)
          phi2j2   = value_if_inside_phi(h1+2,h2+2,mesh,phi)
          phi2j_2  = value_if_inside_phi(h1+2,h2-2,mesh,phi)
          phi_2j2  = value_if_inside_phi(h1-2,h2+2,mesh,phi)

          phi_1j_2 = value_if_inside_phi(h1-1,h2-2,mesh,phi)
          phi1j2   = value_if_inside_phi(h1+1,h2+2,mesh,phi)
          phi1j_2  = value_if_inside_phi(h1+1,h2-2,mesh,phi)
          phi_1j2  = value_if_inside_phi(h1-1,h2+2,mesh,phi)

          phi_2j_1 = value_if_inside_phi(h1-2,h2-1,mesh,phi)
          phi2j1   = value_if_inside_phi(h1+2,h2+1,mesh,phi)
          phi2j_1  = value_if_inside_phi(h1+2,h2-1,mesh,phi)
          phi_2j1  = value_if_inside_phi(h1-2,h2+1,mesh,phi)

          phi_3j = value_if_inside_phi(h1-3,h2,mesh,phi)
          phi_2j = value_if_inside_phi(h1-2,h2,mesh,phi)
          phi_1j = value_if_inside_phi(h1-1,h2,mesh,phi)
          phi1j  = value_if_inside_phi(h1+1,h2,mesh,phi)
          phi2j  = value_if_inside_phi(h1+2,h2,mesh,phi)

          phij_3 = value_if_inside_phi(h1,h2-3,mesh,phi)
          phij_2 = value_if_inside_phi(h1,h2-2,mesh,phi)
          phij_1 = value_if_inside_phi(h1,h2-1,mesh,phi)
          phij   = value_if_inside_phi(h1,h2  ,mesh,phi)
          phij1  = value_if_inside_phi(h1,h2+1,mesh,phi)
          phij2  = value_if_inside_phi(h1,h2+2,mesh,phi)


          ! uh1 =  (phii - phii_1)/ (mesh%delta)  ! order 1 - very bad
          ! uh2 =  (phij - phij_1)/ (mesh%delta) 

          ! uh1 = ( phii1 - phii_1 ) / (2._f64*mesh%delta) ! order 2
          ! uh2 = ( phij1 - phij_1 ) / (2._f64*mesh%delta) 

          ! uh1 = ( phii1/3._f64  + phii/2._f64 - phii_1 + phii_2/6._f64 ) / (mesh%delta) ! order 3
          ! uh2 = ( phij1/3._f64  + phij/2._f64 - phij_1 + phij_2/6._f64 ) / (mesh%delta) 

          uh1 = ( phi_2j + 8._f64 * ( phi1j - phi_1j ) - phi2j ) / (12._f64*mesh%delta) ! order 4
          uh2 = ( phij_2 + 8._f64 * ( phij1 - phij_1 ) - phij2 ) / (12._f64*mesh%delta)

          ! uh1 = ( -phii_3/30._f64 + 0.25_f64*phii_2 - phii_1 + phii/3._f64&
          !      + phii1*0.5_f64  - phii2/20._f64 ) / (mesh%delta) ! order 5
          ! uh2 = ( -phij_3/30._f64 + 0.25_f64*phij_2 - phij_1 + phij/3._f64&
          !      + phij1*0.5_f64  - phij2/20._f64 ) / (mesh%delta)

          ! order 4 approximation made with the values in the x and y directions
          !-> not as good
          ! uxn(i) = - ( phiyj_2 + 8._f64 * ( - phiyj_1 + phiyj1 ) &
          !      - phiyj2 ) / (12._f64*mesh%delta)
          ! uyn(i) = + ( phixi_2 + 8._f64 * ( - phixi_1 + phixi1 ) &
          !      - phixi2 ) / (12._f64*sqrt(3._f64)*mesh%delta)

          !approximation made with the values in the r1 and r2 directions

          uxn(i) = - (uh1*r12+uh2*r22)
          uyn(i) = + (uh1*r11+uh2*r21)




          ! order 2 approximation of the second derivatives
          
          ! uh1h1 = ( phi1j  - 2._f64*phij + phi_1j )/mesh%delta**2
          ! uh2h2 = ( phij1  - 2._f64*phij + phij_1 )/mesh%delta**2
          ! uh1h2 = ( phi1j1 - phi1j_1 - phi_1j1 + phi_1j_1 )/(4._f64*mesh%delta**2)

          ! directe approximation of the second order 
          ! dyuxn(i) = -(phi1j1 - 2._f64*phij + phi_1j_1)/(mesh%delta**2) ! -dyy phi
          ! dxuyn(i) = +(phi1j_1- 2._f64*phij + phi_1j1 )/(3._f64*mesh%delta**2) !  dxx phi

          !  order 4 approximation of the second derivatives
          
          uh1h1 = ( - phi_2j + 16._f64*(phi_1j+phi1j) - 30._f64*phij - phi2j)/(12._f64*mesh%delta**2)
          uh2h2 = ( - phij_2 + 16._f64*(phij_1+phij1) - 30._f64*phij - phij2)/(12._f64*mesh%delta**2)
          uh1h2 = ( 64._f64*(phi1j1 - phi1j_1 - phi_1j1 + phi_1j_1 ) +&
               8._f64*(-phi2j1-phi1j2-phi_2j_1-phi_1j_2 + phi_1j2 + phi_2j1 + phi1j_2 + phi2j_1)+&
               phi2j2 - phi2j_2 - phi_2j2 + phi_2j_2)/(144._f64*mesh%delta**2)

          dxuxn(i) = - (uh1h1*r11*r12+uh1h2*(r11*r22+r12*r21)+uh2h2*r21*r22)  ! -dxy
          !dyuxn(i) = - (uh1h1*r12*r12+2._f64*uh1h2*r12*r22+uh2h2*r22*r22)  ! -dyy phi
          dyuyn(i) = + (uh1h1*r11*r12+uh1h2*(r11*r22+r12*r21)+uh2h2*r21*r22) ! dyx
          dxuyn(i) = + (uh1h1*r11*r11+2._f64*uh1h2*r11*r21+uh2h2*r21*r21)  !dxx phi

          ! directe approximation of the fourth order ->more or less the same results
          dyuxn(i) = -( -phi_2j_2 + 16._f64*(phi_1j_1+phi1j1) - 30._f64*phij - phi2j2 )/&
               (12._f64*mesh%delta**2) ! -dyy phi
          !dxuyn(i) = +( -phi_2j2   + 16._f64*(phi1j_1+phi_1j1) - 30._f64*phij - phi2j_2)/&
          !     (3._f64*12._f64*mesh%delta**2) ! +dxx phi

          ! -> by combining the two approach we get the best results


          ! advection circulaire

          ! uxn(i) = + mesh%cartesian_coord(2,i)
          ! uyn(i) = - mesh%cartesian_coord(1,i)

          ! dxuxn(i) = 0
          ! dxuyn(i) = -1

          ! dyuxn(i) = 1
          ! dyuyn(i) = 0

       end do
    endif


  end subroutine sll_s_compute_hex_fields



  function value_if_inside_phi(k1,k2,mesh,phi) result(f)
    type(sll_t_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:)       :: phi
    sll_int32  :: k1, k2, n
    sll_real64 :: f

    if ( abs(k1) > mesh%num_cells .or. abs(k2) > mesh%num_cells .or. &
         (k1)*(k2)< 0 .and. ( abs(k1) + abs(k2) > mesh%num_cells) ) then
       f = 0._f64 ! null dirichlet boundary condition
    else
       n = mesh%hex_to_global(k1,k2)
       f = phi(n)
    endif

  endfunction value_if_inside_phi

end module sll_m_hex_poisson
