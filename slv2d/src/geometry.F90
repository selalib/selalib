module geometry_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use used_precision

implicit none

type :: geometry
   sll_real64                        :: x0,y0,x1,y1
   sll_real64                        :: dx, dy      
   sll_int32                         :: nx, ny       
   sll_real64, dimension(:), pointer :: xgrid, ygrid  
   character(5)                      :: bc            
end type geometry

interface initialize
   module procedure new_geometry1, new_geometry2
end interface

interface geometry
   module procedure new_geometry
end interface geometry

public :: initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function new_geometry(x0,y0,x1,y1,nx,ny,bc) result(geom)
   
   type(geometry) :: geom
   sll_real64, intent(in) :: x0,y0,x1,y1
   sll_int32, intent(in) :: nx, ny
   character(len=*), intent(in) :: bc 
   sll_int32 :: i    
   sll_int32 :: error 

   geom%x0=x0;geom%x1=x1
   geom%y0=y0;geom%y1=y1
   if ((bc.eq."perxy").or.(bc.eq."perx")) then
      geom%dx=(x1-x0)/nx
   else
      geom%dx=(x1-x0)/(nx-1)
   endif
   if ((bc.eq."perxy").or.(bc.eq."pery")) then
      geom%dy=(y1-y0)/ny
   else
      geom%dy=(y1-y0)/(ny-1)
   endif
   geom%nx=nx
   geom%ny=ny
   geom%bc=bc

   SLL_ALLOCATE(geom%xgrid(nx), error)
   SLL_ALLOCATE(geom%ygrid(ny), error)
   do i=1,nx
      geom%xgrid(i)=geom%x0+(i-1)*geom%dx
   enddo
   do i=1,ny
      geom%ygrid(i)=geom%y0+(i-1)*geom%dy
   enddo

end function new_geometry      

subroutine new_geometry1(geom,x0,y0,nx,ny,dx,dy,error,bc)
   
   type(geometry) :: geom
   sll_real64     :: x0,y0,dx,dy
   sll_int32      :: nx, ny
   character(5)   :: bc 
   sll_int32      :: i    
   sll_int32      :: error

   geom%x0=x0
   geom%y0=y0
   geom%dx=dx
   geom%dy=dy
   geom%nx=nx
   geom%ny=ny
   geom%bc=bc

   SLL_ALLOCATE(geom%xgrid(nx), error)
   SLL_ALLOCATE(geom%ygrid(ny), error)
   do i=1,nx
      geom%xgrid(i)=geom%x0+(i-1)*geom%dx
   enddo
   do i=1,ny
      geom%ygrid(i)=geom%y0+(i-1)*geom%dy
   enddo

 end subroutine new_geometry1

 subroutine new_geometry2(geom,x0,y0,x1,y1,nx,ny,error,bc)
   
   type(geometry), intent(inout) :: geom
   sll_real64, intent(in) :: x0,y0,x1,y1
   sll_int32, intent(in) :: nx, ny
   character(len=*), intent(in) :: bc 
   sll_int32 :: i    
   sll_int32 :: error 

   geom%x0=x0;geom%x1=x1
   geom%y0=y0;geom%y1=y1
   if ((bc.eq."perxy").or.(bc.eq."perx")) then
      geom%dx=(x1-x0)/nx
   else
      geom%dx=(x1-x0)/(nx-1)
   endif
   if ((bc.eq."perxy").or.(bc.eq."pery")) then
      geom%dy=(y1-y0)/ny
   else
      geom%dy=(y1-y0)/(ny-1)
   endif
   geom%nx=nx
   geom%ny=ny
   geom%bc=bc

   SLL_ALLOCATE(geom%xgrid(nx), error)
   SLL_ALLOCATE(geom%ygrid(ny), error)
   do i=1,nx
      geom%xgrid(i)=geom%x0+(i-1)*geom%dx
   enddo
   do i=1,ny
      geom%ygrid(i)=geom%y0+(i-1)*geom%dy
   enddo

end subroutine new_geometry2      

!> compute linear weights of point x,y in mesh geom
subroutine getlinw(this,x,y,inode,jnode,w,error)
   
   type(geometry) :: this
   sll_real64 :: x,y
   sll_int32 :: inode, jnode  
   sll_int32 :: error
   sll_real64, dimension(4) :: w  
   sll_real64 :: aire

   error = 0
   aire = this%dx*this%dy
   ! localize point (x,y) in mesh
   inode = int((x-this%x0)/this%dx)+1
   SLL_ASSERT((inode>=1).or.(inode<=this%nx))

   jnode = int((y-this%y0)/this%dy)+1
   SLL_ASSERT((jnode>=1).or.(jnode<=this%ny))

   ! compute weights
   w(1) = (this%xgrid(inode+1)-x)*(this%ygrid(jnode+1)-y)/aire
   w(2) = (x-this%xgrid(inode))*(this%ygrid(jnode+1)-y)/aire
   w(3) = (x-this%xgrid(inode))*(y-this%ygrid(jnode))/aire
   w(4) = (this%xgrid(inode+1)-x)*(y-this%ygrid(jnode))/aire

end subroutine getlinw

end module geometry_module
