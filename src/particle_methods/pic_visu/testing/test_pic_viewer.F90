program test_pic_viewer

#include "sll_memory.h"
#include "sll_working_precision.h"

use biot_savart
use sll_m_pic_visu
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

class(sll_c_poisson_2d_base),   pointer :: poisson                              
sll_real64, dimension(:,:), allocatable :: vgx                                 
sll_real64, dimension(:,:), allocatable :: vgy                                 
sll_real64, dimension(:,:), allocatable :: omg                                
                                                                                 
type(sll_t_cartesian_mesh_2d), pointer                  :: mesh
class(sll_c_coordinate_transformation_2d_base), pointer :: tau

mesh => sll_f_new_cartesian_mesh_2d(nx, ny, xmin, xmax, ymin, ymax)

write(*,"(3f8.3,i4)") mesh%eta1_min,mesh%eta1_max,mesh%delta_eta1,mesh%num_cells1
write(*,"(3f8.3,i4)") mesh%eta2_min,mesh%eta2_max,mesh%delta_eta2,mesh%num_cells2

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
                                                                              
call poisson%compute_e_from_rho( vgx, vgy, omg )                            
call sll_o_gnuplot_2d(xmin, xmax, nx, ymin, ymax, ny, omg, 'omega', 1, error)


nbpart = 100

call cpu_time(tcpu)
t0 = getRealTimer()
call initialize( nstep, imov, xp, yp, op, delta, dt, nbpart )
SLL_ALLOCATE(up(nbpart), error)
SLL_ALLOCATE(vp(nbpart), error)
iplot   = 0
time = 0.0_f64


do istep = 1, nstep       !loop over time
   
   call vitesse(nbpart, xp, yp, op, up, vp, delta, time)
   call deplace(nbpart, xp, yp, up, vp, dt)

   time = time + dt
   if ( mod(istep, 10) == 0) then
      iplot = iplot + 1
      call sll_o_particles_center_gnuplot( "pic_xy", xp, yp, &
           xmin, xmax, ymin, ymax, iplot, time )
      call sll_o_plot_format_points3d( "pic_xy", xp, yp, op, iplot)
      call sll_s_plot_format_xmdv( "pic_xy", xp, yp, iplot, xmin, xmax, ymin, ymax)
   end if

end do      !next time step

call cpu_time(tcpu)
t1 = getRealTimer()
write(*,"(5x,' CPU time = ', G15.3)") tcpu

end program test_pic_viewer
