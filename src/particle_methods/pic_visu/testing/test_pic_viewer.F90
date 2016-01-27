program test_pic_viewer

#include "sll_memory.h"
#include "sll_working_precision.h"

use biot_savart
use sll_m_pic_visu
use sll_m_pic_viewer
use sll_m_constants
use sll_m_cartesian_meshes
use sll_m_common_coordinate_transformations
use sll_m_coordinate_transformation_2d_base
use sll_m_coordinate_transformations_2d
use sll_m_poisson_2d_base
use sll_m_poisson_2d_fft
use sll_m_gnuplot

implicit none

sll_int32                         :: nbpart
sll_int32                         :: iplot
sll_int32                         :: istep
sll_int32                         :: imov
sll_int32                         :: nstep
sll_real64, dimension(:), pointer :: xp
sll_real64, dimension(:), pointer :: yp
sll_real64, dimension(:), pointer :: op
sll_real64, dimension(:), pointer :: up
sll_real64, dimension(:), pointer :: vp
sll_real64                        :: time
sll_real64                        :: t0
sll_real64                        :: t1
sll_real64                        :: tcpu
sll_int32,  parameter             :: nx = 64
sll_int32,  parameter             :: ny = 48
sll_real64, parameter             :: xmin = -3.0_f64
sll_real64, parameter             :: xmax = +3.0_f64
sll_real64, parameter             :: ymin = -2.0_f64
sll_real64, parameter             :: ymax = +2.0_f64
sll_real64                        :: dt
sll_real64                        :: dx
sll_real64                        :: dy
sll_real64                        :: delta
sll_int32                         :: error

class(sll_c_poisson_2d_base),   pointer     :: poisson                              
sll_real64, dimension(:,:),     allocatable :: vgx                                 
sll_real64, dimension(:,:),     allocatable :: vgy                                 
sll_real64, dimension(:,:),     allocatable :: omg                                
                                                                                 
type(sll_t_cartesian_mesh_2d),                  pointer :: mesh
class(sll_c_coordinate_transformation_2d_base), pointer :: tau
type(sll_t_pic_viewer_2d),                      pointer :: viewer

mesh => sll_f_new_cartesian_mesh_2d(nx, ny, xmin, xmax, ymin, ymax)

write(*,"(2f8.3,i4)") mesh%eta1_min, mesh%eta1_max, mesh%num_cells1
write(*,"(2f8.3,i4)") mesh%eta2_min, mesh%eta2_max, mesh%num_cells2

dx = mesh%delta_eta1
dy = mesh%delta_eta2

! "Identity transformation";
tau => sll_f_new_coordinate_transformation_2d_analytic( &
       "identity_transformation",                       &
       mesh,                                            &
       sll_f_identity_x1,                               &
       sll_f_identity_x2,                               &
       sll_f_identity_jac11,                            &
       sll_f_identity_jac12,                            &
       sll_f_identity_jac21,                            &
       sll_f_identity_jac22,                            &
       SLL_NULL_REAL64 )

call tau%write_to_file(SLL_P_IO_MTV)

SLL_ALLOCATE(vgx(nx,ny),error)                                        
SLL_ALLOCATE(vgy(nx,ny),error)                                        
SLL_ALLOCATE(omg(nx,ny),error)                                       
                                                                              
poisson => sll_f_new_poisson_2d_fft_solver(xmin,xmax,nx,ymin,ymax,ny)
                                                                              
nbpart = 1000

call cpu_time(tcpu)
t0 = getRealTimer()
call initialize( nstep, imov, xp, yp, op, delta, dt, nbpart )
SLL_ALLOCATE(up(nbpart), error)
SLL_ALLOCATE(vp(nbpart), error)
iplot = 1
time  = 0.0_f64

viewer => sll_f_new_pic_viewer_2d( mesh, 'pic_viewer' )

call viewer%set_format(SLL_P_IO_XDMF)

do istep = 1, nstep !loop over time
   
  call vitesse(nbpart, xp, yp, op, up, vp, delta, time)
  call deplace(nbpart, xp, yp, up, vp, dt)

  call compute_grid_vorticity( )

  call write_particles_and_field( viewer, xp, yp, op, omg, istep, time )

  time = time + dt

end do      !next time step

call cpu_time(tcpu)
t1 = getRealTimer()
write(*,"(5x,' CPU time = ', G15.3)") tcpu
print*, 'PASSED'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_grid_vorticity( )

sll_real64 :: a1, a2, a3, a4
sll_int32  :: k 
sll_int32  :: i, j 
sll_real64 :: dpx
sll_real64 :: dpy

omg = 0.0_f64    
                
do k=1,nbpart

  i = floor((xp(k)-xmin)/(xmax-xmin)*nx)+1
  j = floor((yp(k)-ymin)/(ymax-ymin)*ny)+1
  dpx = xp(k) - ((i-1)*dx+xmin)
  dpy = yp(k) - ((j-1)*dy+ymin)

  a1  = (dx-dpx) * (dy-dpy) * op(k)
  a2  = (dpx)    * (dy-dpy) * op(k)
  a3  = (dpx)    * (dpy)    * op(k)
  a4  = (dx-dpx) * (dpy)    * op(k)

  omg(i,j)     = omg(i,j)     + a1 
  omg(i+1,j)   = omg(i+1,j)   + a2 
  omg(i+1,j+1) = omg(i+1,j+1) + a3 
  omg(i,j+1)   = omg(i,j+1)   + a4

end do

omg = omg / (dx*dy) / (dx*dy)

end subroutine compute_grid_vorticity


end program test_pic_viewer
