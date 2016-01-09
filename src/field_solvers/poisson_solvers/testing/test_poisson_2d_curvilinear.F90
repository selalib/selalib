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

!> @internal [example]
program test_poisson_2d_curvilinear
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann, &
    sll_p_periodic

  use sll_m_poisson_2d_curvilinear, only: &
    sll_f_new_poisson_2d_curvilinear, &
    sll_t_poisson_2d_curvilinear, &
    sll_p_poisson_gauss_legendre, &
    sll_p_poisson_open_knots, &
    sll_p_poisson_periodic_knots

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_poisson_2d_curvilinear), pointer :: poisson
  !type(general_coordinate_elliptic_solver), pointer :: poisson_gen
  sll_int32 :: num_cells1
  sll_int32 :: num_cells2
  character(len=256) :: bc_min1_str
  character(len=256) :: bc_max1_str
  character(len=256) :: bc_min2_str
  character(len=256) :: bc_max2_str
  sll_int32 :: bc_min1
  sll_int32 :: bc_max1
  sll_int32 :: bc_min2
  sll_int32 :: bc_max2
  sll_real64 :: eta1_min
  sll_real64 :: eta1_max
  sll_real64 :: eta2_min
  sll_real64 :: eta2_max
  sll_int32 :: spline_degree1
  sll_int32 :: spline_degree2

  character(len=256) :: filename
  sll_int32 :: IO_stat
  sll_int32 :: params_id
  sll_int32 :: bc_knots1
  sll_int32 :: bc_knots2


  namelist /params/ &
    num_cells1, &
    num_cells2, &
    bc_min1_str, &
    bc_max1_str, &
    bc_min2_str, &
    bc_max2_str, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    spline_degree1, &
    spline_degree2

  
  num_cells1 = 128
  num_cells2 = 128
  bc_min1_str = "SLL_DIRICHLET"
  bc_max1_str = "SLL_DIRICHLET"
  bc_min2_str = "SLL_DIRICHLET"
  bc_max2_str = "SLL_DIRICHLET"
  eta1_min = 0._f64
  eta1_max = 1._f64
  eta2_min = 0._f64
  eta2_max = 1._f64
  spline_degree1 = 3
  spline_degree2 = 3

  
  call get_command_argument(1, filename)

  if (len_trim(filename) .ne. 0)then
    open(unit = params_id, file=trim(filename)//'.nml',IOStat=IO_stat)
    if( IO_stat /= 0 ) then
      print *, '#test_poisson_2d_curvilinear failed to open file ', &
        trim(filename)//'.nml'
      STOP
    end if
    print *,'#initialization with filename:'
    print *,'#',trim(filename)//'.nml'
    read(params_id, params) 
    close(params_id)
  else
    print *,'#initialization with default parameters'    
  endif

  
  bc_min1 = boundary_condition(bc_min1_str)
  bc_max1 = boundary_condition(bc_max1_str)
  bc_min2 = boundary_condition(bc_min2_str)
  bc_max2 = boundary_condition(bc_max2_str)

  if(bc_min1==sll_p_periodic)then
    bc_knots1 = sll_p_poisson_periodic_knots
  else
    bc_knots1 = sll_p_poisson_open_knots      
  endif
  if(bc_min2==sll_p_periodic)then
    bc_knots2 = sll_p_poisson_periodic_knots
  else
    bc_knots2 = sll_p_poisson_open_knots      
  endif
  poisson => sll_f_new_poisson_2d_curvilinear( &
    spline_degree1, &
    spline_degree2, &
    num_cells1, &
    num_cells2, &  
    bc_min1, &
    bc_max1, &
    bc_min2, &
    bc_max2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    quadrature_type1=sll_p_poisson_gauss_legendre, &
    quadrature_type2=sll_p_poisson_gauss_legendre, &
    num_quadrature_points1=spline_degree1+2, &
    num_quadrature_points2=spline_degree1+2, &
    bc_knots_min1 = bc_knots1, &
    bc_knots_max1 = bc_knots1, &
    bc_knots_min2 = bc_knots2, &
    bc_knots_max2 = bc_knots2 )

  
!  poisson => sll_f_new_poisson_2d_curvilinear( &
!    spline_degree1, &
!    spline_degree2, &
!    num_cells1, &
!    num_cells2, &
!    sll_p_poisson_gauss_legendre, &
!    sll_p_poisson_gauss_legendre, &
!    bc_left, &
!    bc_right, &
!    bc_bottom, &
!    bc_top, &
!    eta1_min, &
!    eta1_max, &
!    eta2_min, &
!    eta2_max )
!  
!  print *,poisson%knot_size1
!  print *,poisson%knot_size2

!  poisson_gen => new_general_elliptic_solver( &
!    spline_degree1, &
!    spline_degree2, &
!    num_cells1, &
!    num_cells2, &
!    sll_p_poisson_gauss_legendre, &
!    sll_p_poisson_gauss_legendre, &
!    bc_left, &
!    bc_right, &
!    bc_bottom, &
!    bc_top, &
!    eta1_min, &
!    eta1_max, &
!    eta2_min, &
!    eta2_max )
!
!
!  print *,'#diff of global_spline_indices=',maxval(abs(poisson%global_spline_indices-poisson_gen%global_spline_indices))

  
  print *,'#PASSED'

contains

  function boundary_condition(bc) result(res)
    character(len=256), intent(in) :: bc
    sll_int32 :: res
    select case(bc)
      case ("sll_p_periodic")
        res = sll_p_periodic
      case ("SLL_DIRICHLET")
        res = sll_p_dirichlet
      case ("SLL_NEUMANN")
        res = sll_p_neumann
      case default
        print *,'#boundary condition ',bc
        print *,'#not implemented in file ',__FILE__,'at line ',__LINE__
    end select
  end function boundary_condition
  

end program


