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
!> Initialization of particles in 2d+2v: the Landau damping case
!>
!>\author: S. Hirstoaga

module sll_m_particle_initializers_4d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_particle_representation.h"

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_gaussian, only: &
    sll_s_gaussian_deviate_2d

 use sll_m_hammersley, only: &
       sll_f_suite_hamm

  use sll_m_particle_group_4d, only: &
    sll_t_particle_group_4d

  implicit none

  public :: &
       sll_s_initial_random_particles_4d, &
       sll_s_initial_hammersley_particles_4d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
!   private sll_init_spatial_particle2D, suite_hamm
! !  In the future, we have to build a directory called
! !  random_number_generators 


contains

!> Initialize particles in a 2d 2v phase space using random generators.
!> The Landau damping case with a perturbation in a single direction
!> of the physical space.
  subroutine sll_s_initial_random_particles_4d( &
              thermal_speed, alpha, k, &
              m2d,                     &
              num_particles,           &
              p_group,                 &
              rand_seed, rank, worldsize )
    sll_real64, intent(in) :: thermal_speed, alpha, k
    type(sll_t_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32, intent(in)  :: num_particles
    type(sll_t_particle_group_4d), pointer, intent(inout) :: p_group
    sll_int32  :: j, ii
    sll_int32  :: ncx, ic_x,ic_y
    sll_real64 :: x, y, vx, vy
    sll_real64 :: xmin, ymin, rdx, rdy
    sll_real32 :: weight!  sll_real64 :: weight!
    sll_real32 :: off_x,off_y!  sll_real64 :: off_x,off_y
    sll_real64 :: tmp1, tmp2
    sll_int32, dimension(:), intent(in), optional  :: rand_seed
    sll_int32, optional  :: rank, worldsize
    !character(len=8)  :: rank_name
    !character(len=40) :: nomfile
    sll_real64 :: val(1:2)

    if ( present(rand_seed) ) then
       call random_seed (put=rand_seed)
    endif

    if( present(worldsize) ) then
       weight = real((m2d%eta1_max - m2d%eta1_min) * &
            (m2d%eta2_max - m2d%eta2_min),f32)/real(worldsize*num_particles,f32)
    else
       weight = real((m2d%eta1_max - m2d%eta1_min) * &
            (m2d%eta2_max - m2d%eta2_min),f32)/real(num_particles,f32)
    endif
! Each MPI node initialize 'num_particles' particles in the phase space

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
    do while ( j <= num_particles )
       call random_number(x)
!!$       x = vandercorput(ii,3,2)! suite_hamm(ii,3)!
!!$       ii = ii+1
       x = (m2d%eta1_max - xmin)*x + xmin
       call random_number(y)
       y = (1._f64+alpha)*y
       if (sll_f_eval_landau1d(alpha, k, x) >= y ) then
          call random_number(y)
          y = (m2d%eta2_max - ymin)*y + ymin
!!$          y = (m2d%eta2_max - ymin)*vandercorput(j,5,3) + ymin! suite_hamm(j,2)
          !
!-!          nu = thermal_speed*sqrt( -2.0_f64*log(1.0_f64 - &
!-!               (real(j,f64)-0.5_f64)/real(num_particles,f64)) )
!-!          call random_number(yo)
!-!          vx = nu * cos(yo*2.0_f64*sll_p_pi)! cos(vandercorput(j,5,2)*2.0_f64*sll_p_pi)!! yo=suite_hamm(j,5)
!-!          vy = nu * sin(yo*2.0_f64*sll_p_pi)! sin(vandercorput(j,5,2)*2.0_f64*sll_p_pi)!
          call sll_s_gaussian_deviate_2d(val)
          vx = val(1)*thermal_speed
          vy = val(2)*thermal_speed
          !if (j<=50000) write(90,*) x, y, vx, vy 
          SET_PARTICLE_VALUES(p_group%p_list(j),x,y,vx,vy,weight,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
          j = j + 1          
       endif
    end do
    !close(90)

    return
    SLL_ASSERT(present(rank))
  end subroutine sll_s_initial_random_particles_4d


!> Initialize particles in a 2d 2v phase space using pseudo-random generators.
!> The Landau damping case with a perturbation in a single direction
!> of the physical space.
  subroutine sll_s_initial_hammersley_particles_4d( &
              thermal_speed, alpha, k, &
              m2d,                     &
              num_particles,           &
              p_group,                 &
              rand_seed, rank, worldsize )
    sll_real64, intent(in) :: thermal_speed, alpha, k
    type(sll_t_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32, intent(in)  :: num_particles
    type(sll_t_particle_group_4d), pointer, intent(inout) :: p_group
    sll_int32  :: j, ll
    sll_int32  :: ncx, ic_x,ic_y
    sll_real64 :: x, y, vx, vy
    sll_real64 :: xmin, ymin, rdx, rdy
    sll_real32 :: weight
    sll_real32 :: off_x, off_y
    sll_real64 :: tmp1, tmp2
    sll_int32, dimension(:), intent(in), optional  :: rand_seed
    sll_int32, optional  :: rank, worldsize
    sll_real64, dimension(:), allocatable :: pos_x, pos_y
    sll_int32  :: ierr
    sll_real64 :: nu
    if ( present(rand_seed) ) then
       call random_seed (put=rand_seed)
    endif

    if( present(worldsize) ) then
       weight = real((m2d%eta1_max - m2d%eta1_min) * &
            (m2d%eta2_max - m2d%eta2_min),f32)/real(worldsize*num_particles,f32)
    else
       weight = real((m2d%eta1_max - m2d%eta1_min) * &
            (m2d%eta2_max - m2d%eta2_min),f32)/real(num_particles,f32)
    endif
! Each MPI node initialize 'num_particles' particles in the phase space
    SLL_ALLOCATE(pos_x(1:num_particles), ierr)
    SLL_ALLOCATE(pos_y(1:num_particles), ierr)

    rdx = 1._f64/m2d%delta_eta1
    rdy = 1._f64/m2d%delta_eta2
    xmin = m2d%eta1_min
    ymin = m2d%eta2_min
    ncx  = m2d%num_cells1

    !Rejection sampling for the function x --> 1+alpha*cos(k*x)
    ll = 0
!    do j=1,num_particles+ll
!       do
!          x = sll_f_suite_hamm(j,3)
!          x = (m2d%eta1_max - xmin)*x + xmin
!          call random_number(y)
!          y = y * (1._f64+alpha)
!          if (sll_f_eval_landau1d(alpha, k, x) >= y) then
!             pos_x(j) = x
!             pos_y(j) = sll_f_suite_hamm(j,2) 
!             pos_y(j) = (m2d%eta2_max - ymin) * pos_y(j) + ymin
!             exit
!          else
!             ll = ll + 1
!             exit
!          endif
!       enddo
!    enddo
! -------------------------------------------------
! Need to fix that !!!!!
    j = 1
    do while ( j <= num_particles+ll )
       x =  sll_f_suite_hamm(j,3)
       x = (m2d%eta1_max - xmin) * x + xmin
       call random_number(y)
       y = (1._f64 + alpha) * y
       if (sll_f_eval_landau1d(alpha, k, x) >= y ) then
          pos_x(j-ll) = x
          pos_y(j-ll) = sll_f_suite_hamm(j,2) *(m2d%eta2_max - ymin) +ymin
          j = j + 1
       else
          ll = ll + 1 
          j  = j  + 1
       endif
    enddo
    print*, 'le nb d echecs=', ll
    open(91,file='Particuleinitiale.dat')
    do j=1,num_particles 
       nu = sqrt(-2._f64*log(1._f64-(real(j,f64)-0.5_f64)/real(num_particles,f64)))
       vx = nu * cos(sll_f_suite_hamm(j,2) * 2._f64*sll_p_pi)
       vy = nu * sin(sll_f_suite_hamm(j,2) * 2._f64*sll_p_pi)
       vx = vx*thermal_speed
       vy = vy*thermal_speed
       SET_PARTICLE_VALUES(p_group%p_list(j),pos_x(j),pos_y(j),vx,vy,weight,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
       write(91,*) j, pos_x(j),pos_y(j), vx, vy
    enddo  
    close(91)
    SLL_DEALLOCATE_ARRAY(pos_x, ierr)
    SLL_DEALLOCATE_ARRAY(pos_y, ierr)

    return
    SLL_ASSERT(present(rank))
  end subroutine sll_s_initial_hammersley_particles_4d



!> The Landau damping case with a perturbation in BOTH directions
!> of the physical space.
  subroutine sll_initial_particles_4d_L2d( &
              thermal_speed, alpha, &
              kx, ky, m2d,                     &
              num_particles,           &
              p_group,                 &
              rand_seed, rank, worldsize )
    sll_real64, intent(in) :: thermal_speed, alpha
    sll_real64, intent(in) :: kx, ky
    type(sll_t_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32, intent(in)  :: num_particles
    type(sll_t_particle_group_4d), pointer, intent(inout) :: p_group
    sll_int32  :: j, ii, ll
    sll_int32  :: ncx, ic_x,ic_y
    sll_real64 :: x, y, vx, vy,  z
    sll_real64 :: xmin, ymin, rdx, rdy
    sll_real32 :: weight!  sll_real64 :: weight!
    sll_real32 :: off_x,off_y!  sll_real64 :: off_x,off_y
    sll_real64 :: tmp1, tmp2
    sll_int32, dimension(:), intent(in), optional  :: rand_seed
    sll_int32, optional  :: rank, worldsize
    !character(len=8)  :: rank_name
    !character(len=40) :: nomfile
    sll_real64 :: val(1:2)

    if ( present(rand_seed) ) then
       call random_seed (put=rand_seed)
    endif

    if( present(worldsize) ) then
       weight = real((m2d%eta1_max - m2d%eta1_min) * &
            (m2d%eta2_max - m2d%eta2_min),f32)/real(worldsize*num_particles,f32)
    else
       weight = real((m2d%eta1_max - m2d%eta1_min) * &
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
!!$    write(90,*) '#  POSITIONS in 2d    |||    VELOCITIES in 2d'

    j=1
    ii=1
    !Rejection sampling 
    ll = 0
    do while ( j <= num_particles )
       call random_number(x)
       x = (m2d%eta1_max - xmin)*x + xmin
       call random_number(y)
       y = (m2d%eta2_max - ymin)*y + ymin
       call random_number(z)
       z = (1._f64+alpha)*z
       if (sll_f_eval_landau2d(alpha, kx, x, ky, y) >= z ) then
          call sll_s_gaussian_deviate_2d(val)
          vx = val(1)*thermal_speed
          vy = val(2)*thermal_speed
          !if (j<=50000) write(90,*) x, y, vx, vy 
          SET_PARTICLE_VALUES(p_group%p_list(j),x,y,vx,vy,weight,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
          j = j + 1          
       else
          ll=ll+1
       endif
    end do
    print*, 'nb d echecs', ll
    !close(90)

   return
   SLL_ASSERT(present(rank))

  end subroutine sll_initial_particles_4d_L2d


  function sll_f_eval_landau1d(alp, kx, x)
    sll_real64 :: alp, kx, x
    sll_real64 :: sll_f_eval_landau1d
    sll_f_eval_landau1d = 1._f64 + alp * cos(kx * x)
  end function sll_f_eval_landau1d

  function sll_f_eval_landau2d(alp, kx, x, ky, y)
    sll_real64 :: alp, kx, x, ky, y
    sll_real64 :: sll_f_eval_landau2d
    sll_f_eval_landau2d = 1._f64 + alp * cos(kx * x)*cos(ky * y)
  end function sll_f_eval_landau2d

end module sll_m_particle_initializers_4d
