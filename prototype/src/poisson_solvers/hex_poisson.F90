module hex_poisson
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_hex_meshes
  implicit none

contains


  subroutine hex_matrix_poisson( matrix_poisson, mesh)
    type(sll_hex_mesh_2d), pointer         :: mesh
    sll_real64, dimension(:,:), intent(out):: matrix_poisson
    sll_int32                              :: num_cells, global
    sll_int32                              :: index_tab, index_tabij
    sll_int32                              :: index_tabi_1j, index_tabij_1 
    sll_int32                              :: index_tabij1, index_tabi1j
    sll_int32                              :: index_tabi1j1, index_tabi_1j_1
    sll_int32                              :: k1, k2

    num_cells = mesh%num_cells

    matrix_poisson = 0._f64

    ! indexation de la matrice on rempli ligne par ligne avec 
    ! j + (i-1)*n(j) l'indice donnée par index_tab
    ! du coup on a (i-1,j-1), en coordonnées hex, puis (i-1,j)
    ! (i,j-1) ; (i,j) ; (i,j+1)  
    ! (i+1,j) ; (i,j+1)

    do global = 1, mesh%num_pts_tot 

       k1 = mesh%hex_coord(1, global)
       k2 = mesh%hex_coord(2, global)

       call index_hex_to_global(mesh, k1, k2, index_tab)

       call index_hex_to_global(mesh, k1-1, k2-1, index_tabi_1j_1)
       call index_hex_to_global(mesh, k1-1, k2  , index_tabi_1j)
       call index_hex_to_global(mesh, k1  , k2-1, index_tabij_1)
       call index_hex_to_global(mesh, k1  , k2+1, index_tabij1 )
       call index_hex_to_global(mesh, k1+1, k2  , index_tabi1j )
       call index_hex_to_global(mesh, k1+1, k2+1, index_tabi1j1 )

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

  end subroutine hex_matrix_poisson

  
  
  subroutine hex_second_terme_poisson( second_terme, mesh, rho )
    type(sll_hex_mesh_2d), pointer         :: mesh
    sll_real64, dimension(:), intent(out)  :: second_terme
    sll_real64, dimension(:), intent(in)   :: rho
    sll_int32                              :: num_cells, i, index_tab, k1, k2
    sll_real64                             :: step
    sll_real64                             :: f1,f2,f3,f4,f5,f6
    sll_int32                              :: global,n1,n2,n3,n4,n5,n6

    num_cells = mesh%num_cells
    step = mesh%delta**2 * 1.5_f64 

    !************************
    ! general case
    !************************


    do global = 1, mesh%num_pts_tot 

       k1 = mesh%hex_coord(1, global)
       k2 = mesh%hex_coord(2, global)

       call index_hex_to_global(mesh, k1, k2, index_tab)
       
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

    call index_hex_to_global(mesh,num_cells , 0, index_tab)
    second_terme(index_tab) = 0._f64
    call index_hex_to_global(mesh,num_cells , num_cells, index_tab)
    second_terme(index_tab) = 0._f64
    call index_hex_to_global(mesh,0 , num_cells, index_tab)
    second_terme(index_tab) = 0._f64
    call index_hex_to_global(mesh,-num_cells , 0, index_tab)
    second_terme(index_tab) = 0._f64
    call index_hex_to_global(mesh,-num_cells , -num_cells, index_tab)
    second_terme(index_tab) = 0._f64
    call index_hex_to_global(mesh,0 , -num_cells, index_tab)
    second_terme(index_tab) = 0._f64


    ! edges of the hexagon  

    do i = 1,num_cells-1   !( 0 and num_cells  are the corners )

       ! top right edge
       call index_hex_to_global(mesh, num_cells, i, index_tab)
       second_terme(index_tab) = 0._f64
       ! top left edge
       call index_hex_to_global(mesh, num_cells - i, num_cells , index_tab)
       second_terme(index_tab) = 0._f64
       ! left edge
       call index_hex_to_global(mesh, - i, num_cells- i , index_tab)
       second_terme(index_tab) = 0._f64
       ! bottom left edge
       call index_hex_to_global(mesh, - num_cells, - i , index_tab)
       second_terme(index_tab) = 0._f64
       ! bottom right edge
       call index_hex_to_global(mesh, - i, - num_cells , index_tab)
       second_terme(index_tab) = 0._f64
       ! right edge
       call index_hex_to_global(mesh, num_cells - i, - i , index_tab)
       second_terme(index_tab) = 0._f64

    enddo

    !************************
    ! Boundary conditions      -> 0 everywhere at the moment 
    !************************

    ! corners of the hexagon

    call index_hex_to_global(mesh,num_cells-1 , 1, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64
    call index_hex_to_global(mesh,num_cells-1 , num_cells-1, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64
    call index_hex_to_global(mesh,1 , num_cells-1, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64
    call index_hex_to_global(mesh,-num_cells+1 , 1, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64
    call index_hex_to_global(mesh,-num_cells+1 , -num_cells+1, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64
    call index_hex_to_global(mesh,1 , -num_cells+1, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64


    ! edges of the hexagon  

    do i = 2,num_cells-2   !( 1 and num_cells-1 are the corners )

       ! top right edge
       call index_hex_to_global(mesh, num_cells-1, i, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       ! top left edge
       call index_hex_to_global(mesh, num_cells-1 - i, num_cells-1 , index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       ! left edge
       call index_hex_to_global(mesh, - i, num_cells-1 - i , index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       ! bottom left edge
       call index_hex_to_global(mesh, - num_cells+1, - i , index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       ! bottom right edge
       call index_hex_to_global(mesh, - i, - num_cells+1 , index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       ! right edge
       call index_hex_to_global(mesh, num_cells-1 - i, - i , index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64

    enddo


  end subroutine hex_second_terme_poisson

  
  function value_if_inside_rho(k1,k2,mesh,rho) result(f)
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

  endfunction value_if_inside_rho

end module hex_poisson
