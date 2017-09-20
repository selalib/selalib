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
!> @brief Initialization of particles in 2d+2v: the Landau damping case
!> @details Rejection sampling for the perturbation function x --> 1 + alpha*cos(k*x)

!> @author P. Navaro and S. Hirstoaga

module sll_m_particle_initializers_6d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_particle_representation.h"

use sll_m_cartesian_meshes, only: sll_t_cartesian_mesh_3d
use sll_m_constants, only: sll_p_pi
use sll_m_gaussian, only: sll_s_gaussian_deviate_2d
use sll_m_hammersley, only: sll_f_suite_hamm
use sll_m_particle_group_6d, only: sll_t_particle_group_6d

implicit none

public :: sll_s_initial_random_particles_6d

private

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
contains

!> @brief 
!> Initialize particles in a 3d+3v phase space using random generators.
!> @details
!> The Landau damping case with a perturbation in a single direction
!> of the physical space.
subroutine sll_s_initial_random_particles_6d( &
            thermal_speed,                    &
            alpha, k,                         &
            mesh,                             &
            num_particles,                    &
            p_group,                          &
            rand_seed,                        &
            rank, worldsize )         

  !> sigma of Gaussian distribution
  sll_real64, intent(in)                        :: thermal_speed
  !> the perturbation's parameters
  sll_real64, intent(in)                        :: alpha
  !> the perturbation's parameters
  sll_real64, intent(in)                        :: k
  !> the mesh of the physical space
  type(sll_t_cartesian_mesh_3d), intent(in)     :: mesh
  !> particles distributed over a single rank
  sll_int32, intent(in)                         :: num_particles
  !> the random seed used by the rank
  sll_int32, dimension(:), intent(in), optional :: rand_seed
  !> when using more than 1 process
  sll_int32, optional                           :: rank
  !> when using more than 1 process
  sll_int32, optional                           :: worldsize
  !> the particle group
  type(sll_t_particle_group_6d), intent(inout)  :: p_group

  sll_int32  :: j
  sll_int32  :: ic_x, ic_y, ic_z
  sll_int32  :: nc_x, nc_y
  sll_real64 :: x, y, z
  sll_real64 :: xmin, ymin, zmin, rdx, rdy, rdz
  sll_real32 :: weight
  sll_real64 :: tmp_x, tmp_y, tmp_z
  sll_real64 :: nu, fdx, phi, theta, rnd
  sll_real32 :: volume

  volume = real((mesh%eta1_max - mesh%eta1_min) * &
                (mesh%eta2_max - mesh%eta2_min) * &
                (mesh%eta3_max - mesh%eta3_min) , f32)

  if ( present(rand_seed) ) then
     call random_seed (put=rand_seed)
  endif

  if( present(worldsize) ) then
     weight = volume / (real(worldsize,f32)*real(num_particles,f32))
  else
     weight = volume / real(num_particles,f32)
  endif

! Each MPI node initialize 'num_particles' particles in the phase space

  rdx  = 1._f64/mesh%delta_eta1
  rdy  = 1._f64/mesh%delta_eta2
  rdz  = 1._f64/mesh%delta_eta3
  xmin = mesh%eta1_min
  ymin = mesh%eta2_min
  zmin = mesh%eta3_min
  nc_x = mesh%num_cells1
  nc_y = mesh%num_cells2

  j=1
  do while ( j <= num_particles )
     call random_number(x)
     x = (mesh%eta1_max - xmin)*x + xmin
     call random_number(fdx)
     fdx = (1._f64+alpha)*fdx
     if (sll_f_eval_landau1d(alpha, k, x) >= fdx ) then
        call random_number(y)
        y = (mesh%eta2_max - ymin)*y + ymin
        call random_number(z)
        z = (mesh%eta3_max - zmin)*z + zmin
        call random_number(rnd)
        nu = thermal_speed*sqrt(-2._f64*log(1.0_f64 - rnd))
        call random_number(theta)
        call random_number(phi)
        p_group%p_list(j)%vx = nu * cos(theta * 2.0_f64*sll_p_pi) * cos(phi * 2.0_f64*sll_p_pi)
        p_group%p_list(j)%vy = nu * sin(theta * 2.0_f64*sll_p_pi) * cos(phi * 2.0_f64*sll_p_pi)
        p_group%p_list(j)%vz = nu * sin(phi   * 2.0_f64*sll_p_pi)
        tmp_x  = (x - xmin) * rdx 
        tmp_y  = (y - ymin) * rdy 
        tmp_z  = (z - zmin) * rdz 
        ic_x   = int(tmp_x)
        ic_y   = int(tmp_y)
        ic_z   = int(tmp_z)
        p_group%p_list(j)%dx = real(tmp_x - real(ic_x,f64),f32)
        p_group%p_list(j)%dy = real(tmp_y - real(ic_y,f64),f32) 
        p_group%p_list(j)%dz = real(tmp_z - real(ic_z,f64),f32) 
        p_group%p_list(j)%ic = 1 + ic_x + ic_y*nc_x + ic_z*(nc_x*nc_y)
        j      = j + 1          
     endif
  end do

  return
  SLL_ASSERT(present(rank))

end subroutine sll_s_initial_random_particles_6d

function sll_f_eval_landau1d(alp, kx, x)

  sll_real64 :: alp, kx, x
  sll_real64 :: sll_f_eval_landau1d

  sll_f_eval_landau1d = 1._f64 + alp * cos(kx * x)

end function sll_f_eval_landau1d

end module sll_m_particle_initializers_6d
