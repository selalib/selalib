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
!> @author
!> Michel Mehrenberger (mehrenbe@math.unistra.fr)
!> @brief 
!> test of initialization of 1d cartesian mesh from namelist
!> and prevent from too long initializations of a simulation

program test_nml_mesh_1d_cart
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_cartesian_mesh_1d

  use sll_m_nml_mesh_1d_cart, only: &
    sll_s_nml_mesh_1d_cart

  use sll_m_nml_mesh_1d_landau_cart, only: &
    sll_s_nml_mesh_1d_landau_cart

  use sll_m_nml_mesh_1d_two_grid_cart, only: &
    sll_s_nml_mesh_1d_two_grid_cart

  use sll_m_nml_mesh_1d_unif_cart, only: &
    sll_s_nml_mesh_1d_unif_cart

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_cartesian_mesh_1d), pointer :: mesh
  type(sll_cartesian_mesh_1d), pointer :: mesh_1
  sll_real64, pointer :: array(:)
  sll_real64, pointer :: array_1(:)
  sll_real64, pointer :: array_2(:)
  !type(sll_cartesian_mesh_1d), pointer :: mesh_landau
  sll_real64 :: err
  sll_real64 :: err_loc
  character(len=256) :: filename
  character(len=256) :: err_msg

  call get_command_argument( 1, filename )
  if (len_trim(filename) == 0) then
     err_msg = 'filename has to be provided'
    SLL_ERROR( 'test_nml_mesh_1d_cart', trim( err_msg ))
  endif

  call sll_s_nml_mesh_1d_landau_cart(trim(filename))
  call sll_s_nml_mesh_1d_landau_cart(trim(filename),clone="_1")
  call sll_s_nml_mesh_1d_landau_cart(trim(filename),clone="_2")

  call sll_s_nml_mesh_1d_unif_cart(trim(filename))
  call sll_s_nml_mesh_1d_unif_cart(trim(filename),clone="_1")
  call sll_s_nml_mesh_1d_unif_cart(trim(filename),clone="_2")

  call sll_s_nml_mesh_1d_two_grid_cart(trim(filename))
  call sll_s_nml_mesh_1d_two_grid_cart(trim(filename),clone="_1")
  call sll_s_nml_mesh_1d_two_grid_cart(trim(filename),clone="_2")

  call sll_s_nml_mesh_1d_cart(trim(filename),array)
  call sll_s_nml_mesh_1d_cart(trim(filename),array_1,clone="_1")
  call sll_s_nml_mesh_1d_cart(trim(filename),array_2,clone="_2")
  
  print *,'#array=',size(array)-1,array(1),array(size(array))
  print *,'#array1=',size(array_1)-1,array_1(1),array_1(size(array_1))
  print *,'#array2=',size(array_2)-1,array_2(1),array_2(size(array_2))

  call sll_s_nml_mesh_1d_landau_cart(trim(filename),mesh)
  call sll_s_nml_mesh_1d_unif_cart(trim(filename),mesh_1,clone="_1")
  
  print *,'#array=',size(array)-1,array(1),array(size(array))
  print *,'#array1=',size(array_1)-1,array_1(1),array_1(size(array_1))
  print *,'#array2=',size(array_2)-1,array_2(1),array_2(size(array_2))
  print *,'#mesh=',mesh%num_cells,mesh%eta_min,mesh%eta_max
  print *,'#mesh_1=',mesh_1%num_cells,mesh_1%eta_min,mesh_1%eta_max

  
  err = 0._f64  
  err_loc = abs(real(size(array)-33,f64))
  err = max(err,err_loc)
  err_loc = abs(real(size(array_1)-33,f64))
  err = max(err,err_loc)
  err_loc = abs(real(size(array_2)-33,f64))
  err = max(err,err_loc)
  err_loc = abs(real(mesh%num_cells-32,f64))
  err = max(err,err_loc)
  err_loc = abs(real(mesh_1%num_cells-32,f64))
  err = max(err,err_loc)
  
  print*,'#err=',err
  if(err<1.e-14)then  
    print *,'#PASSED'
  endif  
      
end program test_nml_mesh_1d_cart
