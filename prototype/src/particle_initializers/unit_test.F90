program initialize_tester
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_constants, only: sll_pi
  use sll_particle_group_2d_module
  use sll_particle_initializers
  use sll_logical_meshes

#define THERM_SPEED 1._f64
#define NUM_PARTICLES 100000_i64
#define GUARD_SIZE    10000_i64
#define PARTICLE_ARRAY_SIZE 150000_i64
#define ALPHA  0.5_f64
#define NC_X 256_i32
#define XMIN 0._f64
#define KX   0.5_f64
#define XMAX 2._f64*sll_pi/KX
#define NC_Y 64_i32
#define YMIN 0._f64
#define YMAX 1._f64

  
  implicit none
  type(sll_particle_group_2d), pointer :: init_group
  type(sll_logical_mesh_2d),   pointer :: m2d
  sll_int64 :: j
  character(5) :: ncx_name, ncy_name

  init_group => new_particle_2d_group( &
       NUM_PARTICLES, &
       PARTICLE_ARRAY_SIZE, &
       GUARD_SIZE )
  
  m2d =>  new_logical_mesh_2d( NC_X, NC_Y, &
       XMIN, XMAX, YMIN, YMAX )

  call sll_initialize_some4Dfunction( THERM_SPEED, &
       ALPHA, KX, m2d, &
       NUM_PARTICLES, init_group )

  write(ncx_name,'(i3)') NC_X
  write(ncy_name,'(i3)') NC_Y
  open(83,file='vit_pos_'//trim(adjustl(ncx_name))//'x'//trim(adjustl(ncy_name))//'.dat')
  do j=1,NUM_PARTICLES
     write(83,*) init_group%p_list(j)%vx, init_group%p_list(j)%vy, &
     init_group%p_list(j)%ic, init_group%p_list(j)%dx, init_group%p_list(j)%dy, &
     init_group%p_list(j)%dx + m2d%delta_eta1 * real( &!  reconstruction of part%X
       mod(init_group%p_list(j)%ic-1,m2d%num_cells1), f64), &
     init_group%p_list(j)%dy + m2d%delta_eta2  &!  reconstruction of part%Y
     * real( int( (init_group%p_list(j)%ic-1)/m2d%num_cells1) )
  enddo
  close(83)
  
  call sll_delete( init_group )
  call delete( m2d )
  print*, "PASSED"

!contains
  
end program initialize_tester
