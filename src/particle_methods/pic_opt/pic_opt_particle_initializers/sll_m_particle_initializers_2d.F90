!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
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

module sll_m_particle_initializers_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_particle_representation.h"

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_hammersley, only: &
    sll_f_suite_hamm

  use sll_m_particle_group_2d, only: &
    sll_t_particle_group_2d

  implicit none

  public :: &
    sll_s_initial_random_particles_2d, &
    sll_s_initial_random_particles_2d_kh

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
contains

  subroutine sll_s_initial_random_particles_2d_kh( &
              alpha, kx, &
              m2d,                     &
              num_particles,           &
              p_group,                 &
              rand_seed, rank, worldsize )
    sll_real64, intent(in) :: alpha, kx
    type(sll_t_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32, intent(in)  :: num_particles
    type(sll_t_particle_group_2d), pointer, intent(inout) :: p_group
    sll_int32  :: j
    sll_int32  :: ncx, ic_x, ic_y
    sll_real64 :: x, y, z
    sll_real64 :: xmin, ymin, rdx, rdy
    sll_real32 :: weight
    sll_real32 :: off_x, off_y
    sll_real64 :: tmp1, tmp2
    sll_int32, dimension(:), intent(in), optional  :: rand_seed
    sll_int32, optional  :: rank, worldsize
!!$    character(len=8)  :: rank_name
!!$    character(len=40) :: nomfile

    if ( present(rand_seed) ) then
       call random_seed (put=rand_seed)
    endif
    if( present(worldsize) ) then
       weight = real((1.0 + alpha)*(m2d%eta1_max - m2d%eta1_min) * &
            (m2d%eta2_max - m2d%eta2_min),f32)/real(worldsize*num_particles,f32)
    else
       weight = real((1.0 + alpha)*(m2d%eta1_max - m2d%eta1_min) * &
            (m2d%eta2_max - m2d%eta2_min),f32)/real(num_particles,f32)
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
!!$    write(90,*) '#  POSITIONS in 2d'
    
    j=1
    !Rejection sampling for the function 1 + alp + alp*cos(k*x) + sin(y)
    do while ( j <= num_particles )
       call random_number(x)
       x = (m2d%eta1_max - xmin)*x + xmin
       call random_number(y)
       y = (m2d%eta2_max - ymin)*y + ymin
       call random_number(z)
       z = (2.0_f64 + 2.0_f64*alpha) * z
       if (sll_f_eval_KH(alpha, kx, x, y) >= z ) then
!          write(90,*) x, y
          SET_2DPARTICLE_VALUES(p_group%p_list(j),x,y,weight,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
          j = j + 1          
       endif
    end do
!    close(90)
     return
     SLL_ASSERT(present(rank))

   end subroutine sll_s_initial_random_particles_2d_kh



   subroutine sll_s_initial_random_particles_2d( &
               alpha, k, &
               m2d,                     &
               num_particles,           &
               p_group,                 &
               rand_seed, rank, worldsize )
     sll_real64, intent(in) :: alpha, k
     type(sll_t_cartesian_mesh_2d), intent(in) :: m2d
     sll_int32, intent(in)  :: num_particles
     type(sll_t_particle_group_2d), pointer, intent(inout) :: p_group
     sll_int32  :: j
     sll_int32  :: ncx, ic_x, ic_y
     sll_real64 :: x, y, xmin, ymin, rdx, rdy
     sll_real32 :: weight
     sll_real32 :: off_x, off_y
     sll_real64 :: tmp1, tmp2
     sll_int32, dimension(:), intent(in), optional  :: rand_seed
     sll_int32, optional  :: rank, worldsize

     if ( present(rand_seed) ) then
        call random_seed (put=rand_seed)
     endif

     if( present(worldsize) ) then
        weight = real(m2d%eta1_max-m2d%eta1_min,f32) * &
                 real(m2d%eta2_max-m2d%eta2_min,f32) / &
                 real(worldsize*num_particles,f32)
     else
        weight = real(m2d%eta1_max-m2d%eta1_min,f32) * &
                 real(m2d%eta2_max-m2d%eta2_min,f32) / &
                 real(num_particles,f32)
     endif
     rdx = 1._f64/m2d%delta_eta1
     rdy = 1._f64/m2d%delta_eta2
     xmin = m2d%eta1_min
     ymin = m2d%eta2_min
     ncx  = m2d%num_cells1
 
     j=1
     do while ( j <= num_particles )
        call random_number(x)
        x = (m2d%eta1_max - xmin)*x + xmin
        call random_number(y)
        y = (1._f64+alpha)*y
        if (sll_f_eval_landau(alpha, k, x) >= y ) then
           call random_number(y)
           y = (m2d%eta2_max - ymin)*y + ymin
           SET_2DPARTICLE_VALUES(p_group%p_list(j),x,y,weight,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
           j = j + 1          
        endif
     end do
 
     return
     SLL_ASSERT(present(rank))
   end subroutine sll_s_initial_random_particles_2d


   function sll_f_eval_landau(alp, kx, x)
     sll_real64 :: alp, kx, x
     sll_real64 :: sll_f_eval_landau

     sll_f_eval_landau = 1._f64 + alp * cos(kx * x)
   end function sll_f_eval_landau

   function sll_f_eval_KH(alp, kx, x, y)
     sll_real64 :: x,y, alp, kx
     sll_real64 :: sll_f_eval_KH
     
     sll_f_eval_KH = 1.0_f64 + alp + sin(y) + alp*cos(kx*x)
  end function sll_f_eval_KH


end module sll_m_particle_initializers_2d
