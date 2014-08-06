module mesh_hex_alt
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none

contains

  subroutine create_hex_mesh(mesh_num, mesh_numh, mesh_coor, n_cell, radius, x0, y0)

    implicit none

    sll_real64, dimension(:,:),intent(out)     :: mesh_num
    sll_int32, dimension(:,:),intent(out)  :: mesh_numh, mesh_coor
    sll_real64, intent(in)                     :: radius, x0, y0
    sll_int32, intent(in)                  :: n_cell 
    sll_int32                              :: i, j, k, nc1, nc2
    sll_real64, dimension(2)                   :: h1, h2, h3, position
    sll_real64, dimension(2)                   :: h1n, h2n, h3n
    sll_real64                                 :: pas

    pas = radius / real(n_cell,f64)

    h1 = (/1,0/)
    h2 = (/0,1/)
    h3 = (/1,1/)
    
    h1n = ( sqrt(3.0_f64)/2.0_f64*h1 + h2*0.5_f64)*pas
    h2n = (-sqrt(3.0_f64)/2.0_f64*h1 + h2*0.5_f64)*pas
    h3n = h2*pas

    mesh_coor = 0

    nc1 = n_cell + 1
    nc2 = n_cell + 2
    
    position                     = (/x0,y0/)
    mesh_num(1:2,1)              = (/x0,y0/)
    mesh_numh(1:2,1)             = (/0,0/)
    mesh_coor(n_cell+1,n_cell+1) = 1 

    k = 1

    do i = 1, n_cell             ! creation of the mesh, hex. by hex.
   
       position = position + h1n ! ( going to the next hexaedre )

       do j = 1, i
          k = k + 1
          mesh_num(1:2,k)           = position
          mesh_numh(1:2,k)          = (/i,j-1/) 
          mesh_coor(nc1+i,n_cell+j) = k 
          position                  = position + h2n  
       enddo

       do j = 1, i
          k = k + 1
          mesh_num(1:2,k)          = position
          mesh_numh(1:2,k)         = (/i-j+1,i/) 
          mesh_coor(nc2+i-j,nc1+i) = k 
          position                 = position - h1n  
       enddo

       do j = 1, i
          k = k + 1
          mesh_num(1:2,k)          = position
          mesh_numh(1:2,k)         = (/-j+1,i-j+1/) 
          mesh_coor(nc2-j,nc2+i-j) = k 
          position                 = position - h3n  
       enddo

       do j = 1, i
          k = k + 1
          mesh_num(1:2,k)             = position
          mesh_numh(1:2,k)            = (/-i,-j+1/) 
          mesh_coor(nc1-i,n_cell-j+2) = k 
          position                    = position - h2n  
       enddo

       do j = 1, i
          k = k + 1
          mesh_num(1:2,k)             = position
          mesh_numh(1:2,k)            = (/-i+j-1,-i/) 
          mesh_coor(n_cell-i+j,nc1-i) = k 
          position                    = position + h1n  
       enddo

       do j = 1, i
          k = k + 1
          mesh_num(1:2,k)                 = position
          mesh_numh(1:2,k)                = (/j-1,-i+j-1/) 
          mesh_coor(n_cell+j,-i+n_cell+j) = k 
          position                        = position + h3n  
       enddo 

    enddo

  end subroutine create_hex_mesh


  subroutine write_mesh_hex(mesh_num, namefile)
    implicit none  
    sll_real64, dimension(:,:),intent(in) :: mesh_num
    character(len = *), intent(in)        :: namefile
    sll_int32                             :: i, n_points 

    n_points = size(mesh_num(1,:))

    open( unit = 11, file = namefile, action = "write")

    do i = 1 , n_points
       write(11,*) mesh_num(1:2,i)
    enddo

    close(11)

  end subroutine write_mesh_hex
  

  subroutine search_tri( x, y, mesh_coor, step, n_cell, radius, s1, s2, s3 )
    implicit none  
    sll_int32, intent(in) :: n_cell 
    sll_real64, intent(in)    :: x, y, radius, step
    sll_int32, dimension(:,:),intent(out) :: mesh_coor
    sll_real64    :: h1, h2, xi
    sll_int32 :: i, j, n1, n2
    sll_int32, intent(out) :: s1, s2, s3
    sll_int32, dimension(2) :: a, b, c, d

    n1 = n_cell + 1
    n2 = n_cell + 2

    !h1 =  x*sqrt(3.0_f64)/2.0_f64 + y*0.5_f64 ! x/sqrt(3.0_f64) + y*0.5_f64 
    !h2 = -x*sqrt(3.0_f64)/2.0_f64 + y*0.5_f64 ! x/sqrt(3.0_f64) + y*0.5_f64 
    h1 =  x/sqrt(3.0_f64) + y 
    h2 = -x/sqrt(3.0_f64) + y  

    !print*, x, y, sqrt(3.0_f64)*0.5_f64*(h1-h2), h1 + h2, h1, h2
    
!print*, x, y, (h1-h2)/sqrt(3.0_f64), h1 + h2, h1, h2
    
    i = floor( (h1+radius) / step ) - n_cell
    j = floor( (h2+radius) / step ) - n_cell 

    a = (/i+1,j+1/)
    b = (/i,j/)
    c = (/i,j+1/)
    d = (/i+1,j/)

    xi = (real(i,f64) - real(j,f64))*step*sqrt(3.0_f64)*0.5_f64

    if ( x >= xi ) then
       s1 = mesh_coor(n1+i,n1+j)
       s2 = mesh_coor(n2+i,n1+j) 
       s3 = mesh_coor(n2+i,n2+j)
    else 
       s1 = mesh_coor(n1+i,n1+j)
       s2 = mesh_coor(n1+i,n2+j)
       s3 = mesh_coor(n2+i,n2+j)
    endif
    

  end subroutine search_tri

end module mesh_hex_alt
