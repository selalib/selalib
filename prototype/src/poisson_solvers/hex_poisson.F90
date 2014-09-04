module hex_poisson
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_hex_meshes
  implicit none

contains


  subroutine hex_matrix_poisson( matrix_poisson, mesh)
    type(sll_hex_mesh_2d), pointer             :: mesh
    sll_real64, dimension(:,:), intent(out):: matrix_poisson
    sll_int32                              :: num_cells, global
    sll_int32                              :: index_tab, index_tabij
    sll_int32                              :: index_tabi_1j, index_tabij_1 
    sll_int32                              :: index_tabij1, index_tabi1j
    sll_int32                              :: index_tabi1j1, index_tabi_1j_1
    sll_int32                              :: k1, k2

    num_cells = mesh%num_cells

    matrix_poisson = 0._f64


    do global = 1, mesh%num_pts_tot 

       k1 = mesh%hex_coord(1, global)
       k2 = mesh%hex_coord(2, global)

       call index_hex_to_global(mesh, k1, k2, index_tab)

       call index_hex_to_global(mesh, k1-1, k2-1, index_tabi_1j_1)
       call index_hex_to_global(mesh, k1  , k2-1, index_tabi_1j)
       call index_hex_to_global(mesh, k1-1, k2  , index_tabij_1)
       call index_hex_to_global(mesh, k1+1, k2  , index_tabij1 )
       call index_hex_to_global(mesh, k1  , k2+1, index_tabi1j )
       call index_hex_to_global(mesh, k1+1, k2+1, index_tabi1j1 )

       index_tabi_1j_1 = index_tabi_1j_1 - index_tab + 2*num_cells
       index_tabi_1j   = index_tabi_1j   - index_tab + 2*num_cells
       index_tabij_1   = index_tabij_1   - index_tab + 2*num_cells
       index_tabij     = index_tabij     - index_tab + 2*num_cells
       index_tabij1    = index_tabij1    - index_tab + 2*num_cells
       index_tabi1j    = index_tabi1j    - index_tab + 2*num_cells
       index_tabi1j1   = index_tabi1j1   - index_tab + 2*num_cells



       if ( abs(k1) == num_cells .and. abs(k2) == num_cells .or.&
            abs(k1) == num_cells .and. k2 == 0 .or.&
            abs(k2) == num_cells .and. k1 == 0 &
            ) then ! corners

          if (k1 == num_cells.and. k2 == 0) then! corner top right
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) =  0._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j   ) =  0._f64
             matrix_poisson(index_tab,index_tabi1j1  ) =  0._f64
          elseif (k1 == num_cells.and. k2 == num_cells) then! corner top 
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) =  0._f64
             matrix_poisson(index_tab,index_tabi1j   ) =  0._f64
             matrix_poisson(index_tab,index_tabi1j1  ) =  0._f64
          else if (k2 == num_cells.and. k1 == 0) then! corner top left
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) =  0._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) =  0._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) =  0._f64
          else if (k1 == -num_cells.and. k2 == 0) then! corner bottom left
             matrix_poisson(index_tab,index_tabi_1j_1) =  0._f64
             matrix_poisson(index_tab,index_tabi_1j  ) =  0._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) =  0._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          else if (-k1== num_cells.and.k2== -num_cells) then!corner bottom
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) =  0._f64
             matrix_poisson(index_tab,index_tabi1j   ) =  0._f64
             matrix_poisson(index_tab,index_tabi1j1  ) =  0._f64
          else if (k2 == -num_cells.and. k1 == 0) then! corner bottom right
             matrix_poisson(index_tab,index_tabi_1j_1) =  0._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) =  0._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j   ) =  0._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          endif

       else if ( k1*k2<0 .and. abs(k1)+abs(k2)==num_cells .or.&
            abs(k1) == num_cells .or.  abs(k2) == num_cells&
            ) then ! edges

          if (k1 == num_cells) then! edge top right            
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j   ) =  0._f64
             matrix_poisson(index_tab,index_tabi1j1  ) =  0._f64
          elseif (k2 == num_cells) then! edge top left              
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) =  0._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) =  0._f64
          else if (k1*k2<0 .and. -k1+k2==num_cells) then! edge left            
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) =  0._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) =  0._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          else if (k1 == -num_cells) then! edge bottom left            
             matrix_poisson(index_tab,index_tabi_1j_1) =  0._f64
             matrix_poisson(index_tab,index_tabi_1j  ) =  0._f64
             matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          else if (k2 == -num_cells) then! edge bottom right            
             matrix_poisson(index_tab,index_tabi_1j_1) =  0._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) =  0._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
          else if (k1*k2<0 .and. k1-k2==num_cells) then! edge right            
             matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
             matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
             matrix_poisson(index_tab,index_tabij_1  ) =  0._f64
             matrix_poisson(index_tab,index_tabij    ) =  6._f64
             matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             matrix_poisson(index_tab,index_tabi1j   ) =  0._f64
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
    type(sll_hex_mesh_2d), pointer             :: mesh
    sll_real64, dimension(:), intent(out)  :: second_terme
    sll_real64, dimension(:), intent(in)   :: rho
    sll_int32                              :: num_cells, i, index_tab
    sll_real64                             :: step
    sll_int32                              :: k1, k2

    num_cells = mesh%num_cells
    step = mesh%delta**2 * 1.5_f64 

    !************************
    ! general case
    !************************
    second_terme = rho * step

    !************************
    ! Bondary conditions      -> 0 everywhere at the moment 
    !************************

    ! corners of the hexagon

    call index_hex_to_global(mesh,num_cells , 0, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64
    call index_hex_to_global(mesh,num_cells , num_cells, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64
    call index_hex_to_global(mesh,0 , num_cells, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64
    call index_hex_to_global(mesh,-num_cells , 0, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64
    call index_hex_to_global(mesh,-num_cells , -num_cells, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64
    call index_hex_to_global(mesh,0 , -num_cells, index_tab)
    second_terme(index_tab) = second_terme(index_tab) + 0._f64


    ! edges of the hexagon  

    do i = 2,num_cells

       ! top right edge
       k1 = mesh%hex_coord(1 , num_cells)
       k2 = mesh%hex_coord(2 , i)
       call index_hex_to_global(mesh, k1, k2, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       ! top left edge
       k1 = mesh%hex_coord(1, num_cells - i)
       k2 = mesh%hex_coord(2, num_cells)
       call index_hex_to_global(mesh, k1, k2, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       ! left edge
       k1 = mesh%hex_coord(1, - i)
       k2 = mesh%hex_coord(2, num_cells - i)
       call index_hex_to_global(mesh, k1, k2, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       ! bottom left edge
       k1 = mesh%hex_coord(1, - num_cells)
       k2 = mesh%hex_coord(2, - i )
       call index_hex_to_global(mesh, k1, k2, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       ! bottom right edge
       k1 = mesh%hex_coord(1, - i)
       k2 = mesh%hex_coord(2, - num_cells)
       call index_hex_to_global(mesh, k1, k2, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64
       ! right edge
       k1 = mesh%hex_coord(1, num_cells - i)
       k2 = mesh%hex_coord(2,- i )
       call index_hex_to_global(mesh, k1, k2, index_tab)
       second_terme(index_tab) = second_terme(index_tab) + 0._f64

    enddo


  end subroutine hex_second_terme_poisson




end module hex_poisson
