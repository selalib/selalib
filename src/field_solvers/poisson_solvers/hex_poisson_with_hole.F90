module hex_poisson
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_hex_meshes
  implicit none

contains


  subroutine hex_matrix_poisson( matrix_poisson, mesh,n_min, k_min)
    type(sll_hex_mesh_2d), pointer         :: mesh
    sll_real64, dimension(:,:), intent(out):: matrix_poisson
    sll_int32,                   intent(in):: n_min, k_min
    sll_int32                              :: num_cells, global
    sll_int32                              :: index_tab, index_tabij
    sll_int32                              :: index_tabi_1j, index_tabij_1 
    sll_int32                              :: index_tabij1, index_tabi1j
    sll_int32                              :: index_tabi1j1, index_tabi_1j_1
    sll_int32                              :: k1, k2, i

    num_cells = mesh%num_cells

    matrix_poisson = 0._f64

    ! indexation de la matrice on rempli ligne par ligne avec 
    ! j + (i-1)*n(j) l'indice donnée par index_tab
    ! du coup on a (i-1,j-1), en coordonnées hex, puis (i-1,j)
    ! (i,j-1) ; (i,j) ; (i,j+1)  
    ! (i+1,j) ; (i,j+1)

    do global = 1,  mesh%num_pts_tot 

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
          else if (k1== -num_cells.and.k2== -num_cells) then!corner bottom
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
    
    !*******************************************************************
    ! if there is an hexagonal hole 
    !*******************************************************************
    
    if ( n_min > 0 ) then 

       ! the inside of the hole is initialized at zero

       do i = 1,n_min
          k1 = mesh%hex_coord(1, i)
          k2 = mesh%hex_coord(2, i)
          call index_hex_to_global(mesh, k1, k2, index_tab)
          index_tabij = index_tab
          index_tabij = index_tabij - index_tab + 2*num_cells + 2
          matrix_poisson(index_tab,:)            =  0._f64
          matrix_poisson(index_tab,index_tabij ) =  1._f64
       enddo

       do i = n_min+1,n_min+6*k_min

          k1 = mesh%hex_coord(1, i)
          k2 = mesh%hex_coord(2, i)
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

          call index_hex_to_global(mesh, k1, k2, index_tab)

          if ( abs(k1) == k_min .and. abs(k2) == k_min .or.&
               abs(k1) == k_min .and. k2 == 0 .or.&
               abs(k2) == k_min .and. k1 == 0 &
               ) then ! corners

             if (k1 == k_min.and. k2 == 0) then! corner top right
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             elseif (k1 == k_min.and. k2 == k_min) then! corner top
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             else if (k2 == k_min.and. k1 == 0) then! corner top left
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             else if (k1 == -k_min.and. k2 == 0) then! corner bottom left
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             else if (k1== -k_min.and.k2== -k_min) then!corner bottom
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             else if (k2 == -k_min.and. k1 == 0) then! corner bottom right
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             endif

          else if ( k1*k2<0 .and. abs(k1)+abs(k2)==k_min .or.&
               abs(k1) == k_min .or.  abs(k2) == k_min&
               ) then ! edges

             if (k1 == k_min) then! edge top right
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             elseif (k2 == k_min) then! edge top left             
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             else if (k1*k2<0 .and. -k1+k2==k_min) then! edge left 
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             else if (k1 == -k_min) then! edge bottom left        
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             else if (k2 == -k_min) then! edge bottom right      
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             else if (k1*k2<0 .and. k1-k2==k_min) then! edge right 
                matrix_poisson(index_tab,index_tabij    ) =  1._f64
             endif

          elseif ( abs(k1) == k_min+1 .and. abs(k2) == k_min+1 .or.&
               abs(k1) == k_min+1 .and. k2 == 0 .or.&
               abs(k2) == k_min+1 .and. k1 == 0 &
               ) then ! corners

             ! everything is reversed compared to the boundary condition above for the external hexagone

             if (k1 == k_min+1.and. k2 == 0) then! corner top right
                matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
                matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
                matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
             elseif (k1 == k_min+1.and. k2 == k_min+1) then! corner top 
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
                matrix_poisson(index_tab,index_tabij1   ) = -1._f64
                matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
                matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
             else if (k2 == k_min+1.and. k1 == 0) then! corner top left
                matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
                matrix_poisson(index_tab,index_tabij1   ) = -1._f64
                matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
             else if (k1 == -k_min-1.and. k2 == 0) then! corner bottom left
                matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
                matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
                matrix_poisson(index_tab,index_tabij1   ) = -1._f64
             else if (-k1== k_min+1.and.k2== -k_min-1) then!corner bottom
                matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
                matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
                matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
             else if (k2 == -k_min-1.and. k1 == 0) then! corner bottom right
                matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
                matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
                matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
             endif

          else if ( k1*k2<0 .and. abs(k1)+abs(k2)==k_min+1 .or.&
               abs(k1) == k_min+1 .or.  abs(k2) == k_min+1&
               ) then ! edges

             if (k1 == k_min+1) then! edge top right 
                matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
                matrix_poisson(index_tab,index_tabij1   ) = -1._f64
                matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
                matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64
             elseif (k2 == k_min+1) then! edge top left  
                matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
                matrix_poisson(index_tab,index_tabij1   ) = -1._f64
                matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
                matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64       
             else if (k1*k2<0 .and. -k1+k2==k_min+1) then! edge left     
                matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
                matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
                matrix_poisson(index_tab,index_tabij1   ) = -1._f64
                matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64   
             else if (k1 == -k_min-1) then! edge bottom left    
                matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
                matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
                matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
                matrix_poisson(index_tab,index_tabij1   ) = -1._f64   
             else if (k2 == -k_min-1) then! edge bottom right       
                matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
                matrix_poisson(index_tab,index_tabi_1j  ) = -1._f64
                matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
                matrix_poisson(index_tab,index_tabi1j   ) = -1._f64    
             else if (k1*k2<0 .and. k1-k2==k_min+1) then! edge right     
                matrix_poisson(index_tab,index_tabi_1j_1) = -1._f64
                matrix_poisson(index_tab,index_tabij_1  ) = -1._f64
                matrix_poisson(index_tab,index_tabij    ) =  6._f64
                matrix_poisson(index_tab,index_tabi1j   ) = -1._f64
                matrix_poisson(index_tab,index_tabi1j1  ) = -1._f64      
             endif
          endif
       enddo
    endif

  end subroutine hex_matrix_poisson

  
  
  subroutine hex_second_terme_poisson( second_terme, mesh, rho, n_min, k_min)
    type(sll_hex_mesh_2d), pointer        :: mesh
    sll_real64, dimension(:), intent(out) :: second_terme
    sll_real64, dimension(:), intent(in)  :: rho
    sll_int32               , intent(in)  :: n_min, k_min
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

          call index_hex_to_global(mesh, k1, k2, index_tab)

          f1 = value_if_inside_rho(k1+1,k2  ,mesh,rho,k_min)
          f2 = value_if_inside_rho(k1+1,k2+1,mesh,rho,k_min)
          f3 = value_if_inside_rho(k1  ,k2+1,mesh,rho,k_min)
          f4 = value_if_inside_rho(k1-1,k2  ,mesh,rho,k_min)
          f5 = value_if_inside_rho(k1-1,k2-1,mesh,rho,k_min)
          f6 = value_if_inside_rho(k1  ,k2-1,mesh,rho,k_min)

          second_terme(index_tab) = ( 0.75_f64*rho(global) + & 
               (f1+f2+f3+f4+f5+f6)/24._f64 ) * step ! ordre 4

          !second_terme(index_tab) = rho(global) * step ! order 2

       enddo


       !************************
       ! Boundaries
       !************************

       ! corners of the outer hexagon

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


       ! edges of the outer hexagon  

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

       
       ! if there is a hexagonal hole 
       
       if ( n_min > 0 ) then 

          second_terme(1:n_min) = 0._f64

          !************************
          ! inner Boundaries
          !************************

          ! corners of the inner hexagon

          call index_hex_to_global(mesh,k_min , 0, index_tab)
          second_terme(index_tab) = 0._f64
          call index_hex_to_global(mesh,k_min , k_min, index_tab)
          second_terme(index_tab) = 0._f64
          call index_hex_to_global(mesh,0 , k_min, index_tab)
          second_terme(index_tab) = 0._f64
          call index_hex_to_global(mesh,-k_min , 0, index_tab)
          second_terme(index_tab) = 0._f64
          call index_hex_to_global(mesh,-k_min , -k_min, index_tab)
          second_terme(index_tab) = 0._f64
          call index_hex_to_global(mesh,0 , -k_min, index_tab)
          second_terme(index_tab) = 0._f64


          ! edges of the inner hexagon  

          do i = 1,k_min-1   !( 0 and k_min  are the corners )

             ! top right edge
             call index_hex_to_global(mesh, k_min, i, index_tab)
             second_terme(index_tab) = 0._f64
             ! top left edge
             call index_hex_to_global(mesh, k_min - i, k_min , index_tab)
             second_terme(index_tab) = 0._f64
             ! left edge
             call index_hex_to_global(mesh, - i, k_min- i , index_tab)
             second_terme(index_tab) = 0._f64
             ! bottom left edge
             call index_hex_to_global(mesh, - k_min, - i , index_tab)
             second_terme(index_tab) = 0._f64
             ! bottom right edge
             call index_hex_to_global(mesh, - i, - k_min , index_tab)
             second_terme(index_tab) = 0._f64
             ! right edge
             call index_hex_to_global(mesh, k_min - i, - i , index_tab)
             second_terme(index_tab) = 0._f64

          enddo

          !************************
          ! Boundary conditions      -> 0 everywhere at the moment 
          !************************

          ! corners of the hexagon

          call index_hex_to_global(mesh,k_min+1 , 1, index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64
          call index_hex_to_global(mesh,k_min+1 , k_min+1, index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64
          call index_hex_to_global(mesh,1 , k_min+1, index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64
          call index_hex_to_global(mesh,-k_min-1 , 1, index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64
          call index_hex_to_global(mesh,-k_min-1 , -k_min-1, index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64
          call index_hex_to_global(mesh,1 , -k_min-1, index_tab)
          second_terme(index_tab) = second_terme(index_tab) + 0._f64

          ! edges of the hexagon  

          do i = 2,k_min   !( 1 and k_min+1 are the corners )

             ! top right edge
             call index_hex_to_global(mesh, k_min+1, i, index_tab)
             second_terme(index_tab) = second_terme(index_tab) + 0._f64
             ! top left edge
             call index_hex_to_global(mesh, k_min+1 - i, k_min+1 , index_tab)
             second_terme(index_tab) = second_terme(index_tab) + 0._f64
             ! left edge
             call index_hex_to_global(mesh, - i, k_min+1 - i , index_tab)
             second_terme(index_tab) = second_terme(index_tab) + 0._f64
             ! bottom left edge
             call index_hex_to_global(mesh, - k_min-1, - i , index_tab)
             second_terme(index_tab) = second_terme(index_tab) + 0._f64
             ! bottom right edge
             call index_hex_to_global(mesh, - i, - k_min-1 , index_tab)
             second_terme(index_tab) = second_terme(index_tab) + 0._f64
             ! right edge
             call index_hex_to_global(mesh, k_min+1 - i, - i , index_tab)
             second_terme(index_tab) = second_terme(index_tab) + 0._f64

          enddo
       endif

  end subroutine hex_second_terme_poisson


  function value_if_inside_rho(k1,k2,mesh,rho,k_min) result(f)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:)       :: rho
    sll_int32, intent(in)          :: k_min
    sll_int32  :: k1, k2, n
    sll_real64 :: f 
    logical    :: inside, outside

    inside  = .true.
    outside = .false.
    
    ! outer boundary condition
    if ( abs(k1) > mesh%num_cells .or. abs(k2) > mesh%num_cells .or. &
         (k1)*(k2)< 0 .and. ( abs(k1) + abs(k2) > mesh%num_cells) ) outside = .true.
    ! inner boundary condition
    if (k_min == 0 ) then
       inside = .false.
    elseif ( abs(k1) > k_min .or. abs(k2) > k_min .or. &
         (k1)*(k2)< 0 .and. ( abs(k1) + abs(k2) > k_min) ) then
       inside = .false.
    endif

    if ( inside .or. outside ) then
       f = 0._f64 ! null dirichlet boundary condition
    else
       n = hex_to_global(mesh,k1,k2)
       f = rho(n)
    endif

  end function value_if_inside_rho

  subroutine compute_hex_fields(mesh,uxn,uyn,dxuxn,dyuxn,dxuyn,dyuyn,phi,n_min,k_min)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64,dimension(:)        :: uxn, uyn, phi,dxuxn,dyuxn,dxuyn,dyuyn
    sll_int32,          intent(in) :: n_min, k_min
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

    if (n_min>0) then
       uxn  (1:n_min) = 0
       uyn  (1:n_min) = 0
       dxuxn(1:n_min) = 0
       dyuxn(1:n_min) = 0
       dxuyn(1:n_min) = 0
       dyuyn(1:n_min) = 0
    endif

    do i = n_min+1,mesh%num_pts_tot

       h1 = mesh%hex_coord(1,i)
       h2 = mesh%hex_coord(2,i)

       phi_1j_1 = value_if_inside_phi(h1-1,h2-1,mesh,phi,k_min)
       phi1j1   = value_if_inside_phi(h1+1,h2+1,mesh,phi,k_min)
       phi1j_1  = value_if_inside_phi(h1+1,h2-1,mesh,phi,k_min)
       phi_1j1  = value_if_inside_phi(h1-1,h2+1,mesh,phi,k_min)

       phi_2j_2 = value_if_inside_phi(h1-2,h2-2,mesh,phi,k_min)
       phi2j2   = value_if_inside_phi(h1+2,h2+2,mesh,phi,k_min)
       phi2j_2  = value_if_inside_phi(h1+2,h2-2,mesh,phi,k_min)
       phi_2j2  = value_if_inside_phi(h1-2,h2+2,mesh,phi,k_min)

       phi_1j_2 = value_if_inside_phi(h1-1,h2-2,mesh,phi,k_min)
       phi1j2   = value_if_inside_phi(h1+1,h2+2,mesh,phi,k_min)
       phi1j_2  = value_if_inside_phi(h1+1,h2-2,mesh,phi,k_min)
       phi_1j2  = value_if_inside_phi(h1-1,h2+2,mesh,phi,k_min)

       phi_2j_1 = value_if_inside_phi(h1-2,h2-1,mesh,phi,k_min)
       phi2j1   = value_if_inside_phi(h1+2,h2+1,mesh,phi,k_min)
       phi2j_1  = value_if_inside_phi(h1+2,h2-1,mesh,phi,k_min)
       phi_2j1  = value_if_inside_phi(h1-2,h2+1,mesh,phi,k_min)

       phi_3j = value_if_inside_phi(h1-3,h2,mesh,phi,k_min)
       phi_2j = value_if_inside_phi(h1-2,h2,mesh,phi,k_min)
       phi_1j = value_if_inside_phi(h1-1,h2,mesh,phi,k_min)
       phi1j  = value_if_inside_phi(h1+1,h2,mesh,phi,k_min)
       phi2j  = value_if_inside_phi(h1+2,h2,mesh,phi,k_min)

       phij_3 = value_if_inside_phi(h1,h2-3,mesh,phi,k_min)
       phij_2 = value_if_inside_phi(h1,h2-2,mesh,phi,k_min)
       phij_1 = value_if_inside_phi(h1,h2-1,mesh,phi,k_min)
       phij   = value_if_inside_phi(h1,h2  ,mesh,phi,k_min)
       phij1  = value_if_inside_phi(h1,h2+1,mesh,phi,k_min)
       phij2  = value_if_inside_phi(h1,h2+2,mesh,phi,k_min)


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


    end do


  end subroutine compute_hex_fields



  function value_if_inside_phi(k1,k2,mesh,phi,k_min) result(f)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64, dimension(:)       :: phi
    sll_int32, intent(in)          :: k_min
    sll_int32  :: k1, k2, n
    sll_real64 :: f 
    logical    :: inside, outside

    inside  = .true.
    outside = .false.

    ! outer boundary condition
    if ( abs(k1) > mesh%num_cells .or. abs(k2) > mesh%num_cells .or. &
         (k1)*(k2)< 0 .and. ( abs(k1) + abs(k2) > mesh%num_cells) ) outside = .true.
    ! inner boundary condition
    if (k_min == 0 ) then
       inside = .false.
    elseif ( abs(k1) > k_min .or. abs(k2) > k_min .or. &
         (k1)*(k2)< 0 .and. ( abs(k1) + abs(k2) > k_min) ) then
       inside = .false.
    endif

    if ( inside .or. outside ) then
       f = 0._f64 ! null dirichlet boundary condition
    else
       n = hex_to_global(mesh,k1,k2)
       f = phi(n)
    endif

  endfunction value_if_inside_phi

end module hex_poisson
