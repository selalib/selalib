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
      particle_group                                                &
  )

    sll_real64, intent(in)                                  :: thermal_speed, alpha, k_landau
    type(sll_pic_base_class), pointer, intent(inout)   :: particle_group



    ! todo: finish rewriting this old function by Sever and Edwin

    sll_real64, intent(in) :: thermal_speed, alpha, k
    type(sll_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32, intent(in)  :: num_particles
    type(sll_particle_group_4d), pointer, intent(inout) :: p_group



    sll_int32  :: j, ii
    sll_int32  :: ncx, ic_x,ic_y
    sll_real64 :: x, y, vx, vy, nu
    sll_real64 :: xmin, ymin, rdx, rdy
    sll_real32 :: weight!  sll_real64 :: weight!
    sll_real32 :: off_x,off_y!  sll_real64 :: off_x,off_y
    sll_real64 :: tmp1, tmp2
    sll_int32, dimension(:), intent(in), optional  :: rand_seed
    sll_int32, optional  :: rank, worldsize
    character(len=8)  :: rank_name
    character(len=40) :: nomfile
    sll_real64 :: yo, val(1:2)


    SLL_ASSERT( particle_group%dimension_x == 2 )
    SLL_ASSERT( particle_group%dimension_v == 2 )




    if ( present(rand_seed) ) then
       call random_seed (put=rand_seed)
    endif

    if( present(worldsize) ) then
       weight = (m2d%eta1_max - m2d%eta1_min) * &
            (m2d%eta2_max - m2d%eta2_min)/real(worldsize*num_particles,f64)
    else
       weight = (m2d%eta1_max - m2d%eta1_min) * &
            (m2d%eta2_max - m2d%eta2_min)/real(num_particles,f64)
    endif

    rdx = 1._f64/m2d%delta_eta1
    rdy = 1._f64/m2d%delta_eta2
    xmin = m2d%eta1_min
    ymin = m2d%eta2_min
    ncx  = m2d%num_cells1

!!$    if(present(rank)) then
!!$       write(rank_name,'(i8)') rank
!!$    else
!!$       rank_name = '00000000'
!!$    end if
!!$    nomfile='initialparts_'//trim(adjustl(rank_name))//'.dat'
!!$    open(90, file=nomfile)
!!$
!!$    write(90,*) '#  POSITIONS in 2d    |||    VELOCITIES in 2d'

    j=1
    ii=1
    !Rejection sampling for the function x --> 1+alpha*cos(k*x)
!Each MPI node initialize 'num_particles' particles in phys space and velocity
    do while ( j <= num_particles )
       call random_number(x)
!!$       x = vandercorput(ii,3,2)! suite_hamm(ii,3)!
!!$       ii = ii+1
       x = (m2d%eta1_max - xmin)*x + xmin
       call random_number(y)
       y = (1._f64+alpha)*y! 2._f64 * y
       if (eval_landau(alpha, k, x) >= y ) then
          call random_number(y)
          y = (m2d%eta2_max - ymin)*y + ymin
!!$          y = (m2d%eta2_max - ymin)*vandercorput(j,5,3) + ymin! suite_hamm(j,2)
          !
!-!          nu = thermal_speed*sqrt( -2.0_f64*log(1.0_f64 - &
!-!               (real(j,f64)-0.5_f64)/real(num_particles,f64)) )
!-!          call random_number(yo)
!-!          vx = nu * cos(yo*2.0_f64*sll_pi)! cos(vandercorput(j,5,2)*2.0_f64*sll_pi)!! yo=suite_hamm(j,5)
!-!          vy = nu * sin(yo*2.0_f64*sll_pi)! sin(vandercorput(j,5,2)*2.0_f64*sll_pi)!
          call gaussian_deviate_2D(val)
          vx = val(1)*thermal_speed
          vy = val(2)*thermal_speed
          !if (j<=50000) write(90,*) x, y, vx, vy
          SET_PARTICLE_VALUES(p_group%p_list(j),x,y,vx,vy,weight,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
          j = j + 1
       endif
    end do
    !close(90)
  end subroutine sll_pic_4d_random_unweighted_initializer_landau_f0


end module sll_pic_random_initializers
