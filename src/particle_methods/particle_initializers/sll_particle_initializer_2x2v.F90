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

module sll_m_particle_initializer

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_module_pic_base
  use gaussian, only: gaussian_deviate_2D
  
  implicit none


  type :: sll_particle_initializer_2x2v
     

   contains
     procedure :: draw_product_density => draw_product_density_2x2v
     
  end type sll_particle_initializer_2x2v
 

contains
  
  subroutine draw_product_density_2x2v (this, particle_group, random_seed)
    class(sll_particle_initializer_2x2v), intent(inout)         :: this
    class(sll_particle_group_base),  intent(inout)              :: particle_group
    !sll_real64, intent(in)                                      :: thermal_velocity
    sll_int32,  intent(in)                                      :: random_seed


    !print*, this%densityx1(0.0)
    

  end subroutine draw_product_density_2x2v


  subroutine sll_particle_initialize_landaux1_2x2v(&
       particle_group, &
       landau_param, &
       xmin, &
       Lx, &
       thermal_velocity, &
       rnd_seed)
    class(sll_particle_group_base),  intent(inout)              :: particle_group
    sll_real64, intent(in)                                      :: landau_param(2)
    sll_real64, intent(in)                                      :: xmin(2)
    sll_real64, intent(in)                                      :: Lx(2)
    sll_real64, intent(in)                                      :: thermal_velocity(2)
    sll_int32,  intent(in)                                      :: rnd_seed(:)

    sll_real64                                                  :: x(3),v(3)
    sll_real64                                                  :: rnd_no
    sll_int32                                                   :: i_part
    

    x = 0.0_f64
    v = 0.0_f64

    call random_seed (put=rnd_seed)

    i_part = 0
    do while( i_part < particle_group%n_particles)
       ! Rejection sampling for x(1)
       call random_number(rnd_no)
       x(1) = Lx(1) * rnd_no + xmin(1)
       call random_number(rnd_no)
       rnd_no = (1.0_f64 + landau_param(1)) * rnd_no
       if ( (1.0_f64 + landau_param(1) * cos(landau_param(2)*x(1))) >= rnd_no ) then  
          ! Draw x(2) from uniform distribution
          i_part = i_part + 1
          call random_number(rnd_no)
          x(2) = Lx(2) * rnd_no + xmin(2)

          ! Draw velocities from 2D Gauss distribution
          call gaussian_deviate_2D(v(1:2))
          v(1:2) = v(1:2) * thermal_velocity

          call particle_group%set_x(i_part, x)
          call particle_group%set_v(i_part, v)
          call particle_group%set_weight(i_part, &
               1.0_f64/real(particle_group%n_total_particles, f64))
       end if
    end do

    

  end subroutine sll_particle_initialize_landaux1_2x2v


end module sll_m_particle_initializer
