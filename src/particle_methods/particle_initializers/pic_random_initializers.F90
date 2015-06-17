!**************************************************************
!  Copyright INRIA
!  Authors : 
!     ???
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

module sll_pic_random_initializers

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "particle_representation.h"

  use sll_constants, only: sll_pi
  use sll_simple_pic_4d_group_module
!  use sll_pic_utilities
!  use sll_lt_pic_4d_utilities
  use sll_cartesian_meshes
!  use sll_representation_conversion_module

  implicit none
  

contains
   
  ! initialize a simple_pic group with the landau f0 distribution
  !
  ! this routine takes as arguments a newly created particle group of simple_pic type
  ! and it creates the particle list inside this group
  subroutine sll_pic_4d_random_unweighted_initializer_landau_f0 (   &
      thermal_speed, alpha, k_landau,                               &
      particle_group,                                               &
      space_mesh_2d,                                             &
      number_particles,                                             &
      rand_seed, rank, worldsize                                    &
    )

    sll_real64,                         intent(in)              :: thermal_speed, alpha, k_landau
    type(sll_pic_base_class), pointer,  intent(inout)           :: particle_group
    type(sll_cartesian_mesh_2d),        intent(in)              :: space_mesh_2d
    sll_int32,                          intent(in)              :: number_particles

    sll_int32, dimension(:),            intent(in), optional    :: rand_seed
    sll_int32,                          intent(in), optional    :: rank, worldsize

    sll_int32  :: i_part, ii
    sll_int32  :: ncx, ic_x, ic_y
    sll_int32  :: effective_worldsize
    sll_real64 :: x_min, y_min, x_max, y_max, rdx, rdy
    sll_real32 :: weight
    sll_real64 :: aux_random
    sll_real64 :: x(3)
    sll_real64 :: v(3)
    sll_real64 :: val(1:2)

    SLL_ASSERT( particle_group%dimension_x == 2 )
    SLL_ASSERT( particle_group%dimension_v == 2 )

    rdx = 1._f64/space_mesh_2d%delta_eta1
    rdy = 1._f64/space_mesh_2d%delta_eta2
    x_min = space_mesh_2d%eta1_min
    x_max = space_mesh_2d%eta1_max
    y_min = space_mesh_2d%eta2_min
    y_max = space_mesh_2d%eta2_max

    ncx  = space_mesh_2d%num_cells1

    if ( present(rand_seed) ) then
       call random_seed (put=rand_seed)
    end if

    if( present(worldsize) ) then
      effective_worldsize = worldsize
    else
      effective_worldsize = 1
    end if
    weight = (x_max - x_min) * (y_max - y_min) / real(effective_worldsize * num_particles,f64)
    particle_group%set_common_weight( weight )

    x = 0
    v = 0

    i_part = 1
    ii = 1
    ! Rejection sampling for the function x --> 1 + alpha * cos(k * x)
    ! Each MPI node initializes 'num_particles' particles in phys space and velocity
    do while ( i_part <= num_particles )
      call random_number(aux_random)
      x(1) = (x_max - x_min) * aux_random + x_min
      call random_number(aux_random)
      x(2) = (1._f64 + alpha) * aux_random
      if ( eval_landau(alpha, k, x(1)) >= x(2) ) then
        call random_number(aux_random)
        x(2) = (y_max - y_min) * aux_random + y_min
        call gaussian_deviate_2D(val)
        v(1) = val(1) * thermal_speed
        v(2) = val(2) * thermal_speed

        particle_group%set_x( i_part, x )
        particle_group%set_v( i_part, v )

        !        SET_PARTICLE_VALUES(p_group%p_list(j),x,y,vx,vy,weight,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
        i_part = i_part + 1
       end if
    end do
  end subroutine sll_pic_4d_random_unweighted_initializer_landau_f0


  function eval_landau(alpha, kx, x)
    sll_real64 :: alpha, kx, x
    sll_real64 :: eval_landau
    eval_landau = 1._f64 + alpha * cos(kx * x)
  end function eval_landau

end module sll_pic_random_initializers
