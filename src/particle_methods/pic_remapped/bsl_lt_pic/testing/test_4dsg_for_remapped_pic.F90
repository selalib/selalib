! \file 
! \brief Testing the sparse interpolation algorithm from [[selalib:src/add_ons/sparse_grid]]
! \author MCP and ALH
!   This code SeLaLib (for Semi-Lagrangian-Library) 
!   is a parallel library for simulating the plasma turbulence 
!   in a tokamak.
! 
!   This software is governed by the CeCILL-B license 
!   under French law and abiding by the rules of distribution 
!   of free software.  You can  use, modify and redistribute 
!   the software under the terms of the CeCILL-B license as 
!   circulated by CEA, CNRS and INRIA at the following URL
!   "http://www.cecill.info". 
! 
! [[elisp:(alh-set-keywords)]] ([[file:~/.emacs::alh-set-keywords]])
! emacs-keywords authors="MCP and ALH" brief="Testing the sparse interpolation algorithm from [[selalib:src/add_ons/sparse_grid]]" dox f90 selalib start=02/11/15

! <<<<test_4dsg_for_remapped_pic>>>> The goal of this test is to reproduce the same function interpolations as in
! [[selalib:src/add_ons/sparse_grid/testing/test_sparse_grid_4d.F90::test_interpolation_4d]] (from directory
! [[selalib:src/add_ons/sparse_grid/testing]]) <<ALH>>

! Compilation commands in [[file:CMakeLists.txt::test_4dsg_for_remapped_pic]]

program test_4dsg_for_remapped_pic

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_m_constants, only : sll_pi
  use sll_m_sparse_grid_4d
  implicit none

  ! Number of hierarchical levels
  sll_int32 :: levels

  ! maximum level in the sparse grid along each direction
  sll_int32,dimension(4) :: levelsini

  ! Order of the sparse grid functions (degree of polynomial)
  sll_int32 :: order
  
  ! 4d domain size
  sll_real64,dimension(4) :: eta_min,eta_max

  ! 4d test point
  sll_real64,dimension(4) :: dx

  ! <<interp>> Interpolation object of type
  ! [[selalib:src/add_ons/sparse_grid/sll_m_sparse_grid_4d.F90::sparse_grid_interpolator_4d]] which extends
  ! [[selalib:src/add_ons/sparse_grid/sll_m_sparse_grid_interpolator.F90::type,%20public%20::%20sparse_grid_interpolator]]
  
  type(sparse_grid_interpolator_4d), target   :: interp

  ! <<f>> Function to interpolate
  sll_real64, dimension(:), allocatable :: f

  ! Function values
  sll_real64 :: finterp,fref,ferr,ferrmax
  
  sll_int32 :: i,j,ierr
  
  ! ------------ end of declarations
  
  ! [[eta_min]], [[eta_max]]
  eta_min(1) = 0.0_f64; eta_max(1) = 2.0_f64*sll_pi;
  eta_min(2) = 0.0_f64; eta_max(2) = 2.0_f64*sll_pi;
  eta_min(3) = 0.0_f64; eta_max(3) = 2.0_f64*sll_pi;
  eta_min(4) = 0.0_f64; eta_max(4) = 2.0_f64*sll_pi;

  ! [[levels]] This is a tunable parameter
  levels = 12;
  levelsini(1)=levels;
  levelsini(2)=levels;
  levelsini(3)=levels;
  levelsini(4)=levels;

  ! [[order]] This is a tunable parameter.
  order = 1;
  
  ! Initialize sparse grid [[selalib:src/add_ons/sparse_grid/sll_m_sparse_grid_4d.F90::subroutine%20initialize_sg4d]] in
  ! [[selalib:src/add_ons/sparse_grid/]]
  
  call interp%initialize(levelsini,order,order+1,0,eta_min,eta_max);

  ! Initialize [[f]].  [[selalib:src/add_ons/sparse_grid/sll_m_sparse_grid_interpolator.F90::size_basis]] is the number of grid
  ! points.

  SLL_ALLOCATE(f(interp%size_basis),ierr);
  do j=1,interp%size_basis

     ! [[testfunction]]
     
     f(j) = testfunction(                    &
          interp%hierarchy(j)%coordinate(1), &
          interp%hierarchy(j)%coordinate(2), &
          interp%hierarchy(j)%coordinate(3), &
          interp%hierarchy(j)%coordinate(4));
  end do

  ! Evaluate a few random points
  ! [[file:~/selalib/src/add_ons/sparse_grid/sll_m_sparse_grid_interpolator.F90::subroutine%20compute_hierarchical_surplus]]

  ferrmax = 0._f64
  call interp%compute_hierarchical_surplus(f)
  do i = 1,10
     dx(1) = randompoint(eta_min(1),eta_max(1));
     dx(2) = randompoint(eta_min(2),eta_max(2));
     dx(3) = randompoint(eta_min(3),eta_max(3));
     dx(4) = randompoint(eta_min(4),eta_max(4));

     ! <<interpolate_value>> [[file:~/selalib/src/add_ons/sparse_grid/sll_m_sparse_grid_4d.F90::function%20interpolate_value]]

     finterp = interp%interpolate_value(f,dx);
     fref = testfunction(dx(1),dx(2),dx(3),dx(4));
     ferr = abs(finterp-fref)
     ferrmax = max(ferr,ferrmax)
     
     ! [[http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html]]
     write(*,"(A,I3,A,G12.5,G12.5,G12.5,G12.5,G12.5)") "point",i," ",dx(1),dx(2),dx(3),dx(4),ferr
  end do

  ! Accuracy depends strongly on the [[tunable]] parameters above.
  
  write(*,"(A,G12.5)") "maximum error = ",ferrmax
  if (ferrmax <= 1e-2) then
     print *,'PASSED'
  else
     print *,'FAILED'
  end if

  DEALLOCATE(f);

  ! ------------ end of program

contains
  
  ! <<testfunction>> An arbitrary 4d function

  function testfunction(x,y,vx,vy) result(f)
    implicit none
    sll_real64 :: f
    sll_real64,intent(in) :: x,y,vx,vy
    f = sin(x) * cos(y) * sin(vx) * cos(vy);
  end function testfunction

  ! <<randompoint>> Returns a random point in a [min,max] interval

  function randompoint(min,max) result(x)
    implicit none
    sll_real64 :: x
    sll_real64,intent(in) :: min,max
    call random_number(x)
    x = min + (max-min)*x
  end function randompoint

end program test_4dsg_for_remapped_pic

! Local Variables:
! mode:F90
! ispell-local-dictionary:"british"
! coding:utf-8
! fill-column:132
! eval:(flyspell-prog-mode)
! eval:(outline-minor-mode)
! End:
! LocalWords: emacs
