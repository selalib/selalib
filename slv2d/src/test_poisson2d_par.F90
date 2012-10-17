program test_poisson2d_par

#include "selalib.h"

use mpi
use module_mpi 
use geometry_module
use poisson2dpp_seq

implicit none

type(geometry)    :: geom     
type(poisson2dpp) :: poisson 

sll_int32       :: nx, ny      
sll_real64      :: x_min, x_max
sll_real64      :: y_min, y_max
sll_int32       :: error     

sll_real64, dimension(:,:),pointer :: x
sll_real64, dimension(:,:),pointer :: y
sll_real64, dimension(:,:),pointer :: rho
sll_real64, dimension(:,:),pointer :: ex
sll_real64, dimension(:,:),pointer :: ey

sll_int32  :: my_num
sll_int32  :: num_threads
sll_real64 :: tcpu1
sll_real64 :: tcpu2

! initialisation global
tcpu1 = MPI_WTIME()
call initialise_moduleMPI

call MPI_COMM_RANK(MPI_COMM_WORLD, my_num, error)
call MPI_COMM_SIZE(MPI_COMM_WORLD, num_threads, error)
if (my_num == MPI_MASTER) then
   print*,'We are running on ',num_threads, ' processors'
end if

x_min = 0.0
x_max = 2.*sll_pi
y_min = 0.0
y_max = 2.*sll_pi
nx  = 64
ny  = 64


SLL_ALLOCATE(rho(nx,ny),error)
SLL_ALLOCATE(ex(nx,ny),error)
SLL_ALLOCATE(ey(nx,ny),error)

call new_geometry2(geom,x_min,y_min,x_max,y_max,nx,ny,error,"perxy")

call meshgrid(geom%xgrid, geom%ygrid, x, y)
rho = -2_f64 * sin(x) * sin(y)
call plot_solution( x, y, rho, "rho" )

SLL_ALLOCATE(rho(geom%nx,geom%ny),error)
SLL_ALLOCATE(ex(geom%nx,geom%ny),error)
SLL_ALLOCATE(ey(geom%nx,geom%ny),error)

call new(poisson, rho, geom, error)

if (my_num == MPI_MASTER) then
   write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
   write(*,"(2(i3,1x),6(g13.3,1x))") geom%nx, geom%ny, geom%x0, &
                                     geom%x0+(geom%nx)*geom%dx, &
                                     geom%y0, geom%y0+(geom%ny)*geom%dy, &
                                     geom%dx, geom%dy   
endif

call solve(poisson,ex,ey,rho)

call plot_solution( x, y, ex, "ex" )

call plot_solution( x, y, ey, "ey" )

print*, " error ex : ", sum(abs(ex + cos(x)*cos(y)))/(nx*ny)
print*, " error ey : ", sum(abs(ey - sin(x)*sin(y)))/(nx*ny)

tcpu2 = MPI_WTIME()
if (my_num == MPI_MASTER) &
   write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads

call termine_moduleMPI

print*,'PASSED'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_solution( x, y, z, field_name )

   sll_real64, dimension(:,:), intent(in) :: x
   sll_real64, dimension(:,:), intent(in) :: y
   sll_real64, dimension(:,:), intent(in) :: z
   character(len=*), intent(in) :: field_name
   character(len=4) :: prefix = "mesh"
   sll_int32 :: file_id

   call sll_xdmf_open(field_name//".xmf",prefix,size(z,1),size(z,2),file_id,error)
   call sll_xdmf_write_array(prefix,x,'x1',error)
   call sll_xdmf_write_array(prefix,y,'x2',error)
   call sll_xdmf_write_array(prefix,z,field_name,error,file_id,"Node")
   call sll_xdmf_close(file_id,error)

end subroutine plot_solution

subroutine meshgrid(eta1, eta2, x, y)

   sll_real64, dimension(:), intent(in) :: eta1
   sll_real64, dimension(:), intent(in) :: eta2
   sll_real64, dimension(:,:), intent(out), pointer :: x
   sll_real64, dimension(:,:), intent(out), pointer :: y
   sll_int32 :: i, j, error
   
   SLL_ALLOCATE(x(size(eta1),size(eta2)),error)
   SLL_ALLOCATE(y(size(eta1),size(eta2)),error)
   
   do j = 1, size(eta2)
      do i = 1, size(eta2)
         x(i,j) = eta1(i)
         y(i,j) = eta2(j)
      end do
   end do

end subroutine meshgrid

end program test_poisson2d_par
