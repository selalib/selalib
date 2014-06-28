program particle_sorter
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "particle_representation.h"
  
  use sll_constants, only: sll_pi
  use sll_logical_meshes
  use sll_particle_sort_module
  use sll_particle_group_2d_module
  use gaussian
  use hammersley

#define THERM_SPEED 1._f64
#define num_particles 100000_i32
#define KX    0.5_f64
#define XMIN  0._f64
#define XMAX (2._f64*sll_pi/KX)
#define NC_X  256_i32
#define YMIN  0._f64
#define YMAX  1._f64
#define NC_Y  64_i32
#define ALPHA  0.1_f64
#define QoverM 1._f64


  type(sll_particle_group_2d), pointer  :: group
  type(sll_logical_mesh_2d), pointer    :: m
  type(sll_particle_sorter_2d), pointer :: sorter
!  sll_int32  :: i
  sll_real64 :: x, y, vx, vy, nu, xmin, ymin, rdx, rdy
  sll_int32  :: j
  sll_real32 :: weight, off_x,off_y
  sll_int32  :: ierr, ncx, ic_x,ic_y
  sll_real64 :: tmp1, tmp2

  m => new_logical_mesh_2d( NC_X, NC_Y, XMIN, XMAX, YMIN, YMAX )

  group => new_particle_2d_group( int(num_particles,i32),  &
       int(num_particles,i32), int(num_particles/4,i32), QoverM, m )


  sorter => sll_new_particle_sorter_2d( m )
! the arguments to new_particle_group should be 32bit ints... change

  weight = (m%eta1_max - m%eta1_min) * &
           (m%eta2_max - m%eta2_min)/real(num_particles,f32) 
  rdx = 1._f64/m%delta_eta1
  rdy = 1._f64/m%delta_eta2
  xmin = m%eta1_min
  ymin = m%eta2_min
  ncx  = m%num_cells1

  open(90,file='initialparticles.dat')
  j=1
  do while ( j <= num_particles )
     call random_number(x)
     x = (m%eta1_max - m%eta1_min)*x + m%eta1_min
     call random_number(y)
     y = 2._f64 * y! (2._f64*alpha)*y + 1._f64 - alpha
     if (eval_landau(ALPHA, KX, x) >= y ) then
        y = (m%eta2_max - m%eta2_min)*suite_hamm(j,3) + m%eta2_min
!
        nu = THERM_SPEED*sqrt( -2.0_f64*log(1.0_f64 - &
             (real(j,f64)-0.5_f64)/real(num_particles,f64)) )
        vx = nu * cos(suite_hamm(j,2)*2.0_f64*sll_pi)
        vy = nu * sin(suite_hamm(j,2)*2.0_f64*sll_pi)
        write(90,*) x, y, vx, vy 
!        call set_group_particle_values( group, j, x, y, vx, vy, weight)
        SET_PARTICLE_VALUES(group%p_list(j),x,y,vx,vy,weight,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
        j = j + 1          
     endif
  end do; close(90)

  call sll_sort_particles_2d( sorter, group )

  do j=1,num_particles-1
     if ( group%p_list(j)%ic > group%p_list(j+1)%ic ) then
        print*, 'BAD order of ic:', 'j=', j, 'j+1=', j+1
        stop
     endif
  end do 

!  call sll_delete (sorter) 
  print*, 'PASSED'

contains
    
  function eval_landau(alp, kx, x)
    sll_real64 :: alp, kx, x
    sll_real64 :: eval_landau
    eval_landau = 1._f64 + alp * cos(kx * x)
  end function eval_landau
  
end program particle_sorter
