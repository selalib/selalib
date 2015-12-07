!**************************************************************
!  Copyright INRIA
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
!
!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!

!> @ingroup maxwell_solvers
!> This module contains common subroutines for  Maxwell solvers
module sll_m_maxwell_solvers_base

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_utilities, only: &
    int2string

  implicit none

  public :: &
    sll_plot_two_fields

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Parent object of all Maxwell solvers
  type :: sll_maxwell_solver

   sll_int32  :: nc_eta1      !< x cells number
   sll_int32  :: nc_eta2      !< y cells number
   sll_int32  :: polarization !< TE or TM
   sll_real64 :: e_0          !< electric conductivity
   sll_real64 :: mu_0         !< magnetic permeability
   sll_real64 :: c            !< speed of light
   sll_real64 :: eta1_min     !< left side 
   sll_real64 :: eta1_max     !< right side
   sll_real64 :: delta_eta1   !< step size
   sll_real64 :: eta2_min     !< bottom side
   sll_real64 :: eta2_max     !< top side
   sll_real64 :: delta_eta2   !< step size

  end type sll_maxwell_solver

contains

!> write files to visualize 2d fields with gnuplot
subroutine sll_plot_two_fields(fname, n1, n2, f1, f2, iplot, time )
character(len=*),             intent(in) :: fname !< output file name
sll_int32,                    intent(in) :: n1    !< size of f1 and f2 first index
sll_int32,                    intent(in) :: n2    !< size of f1 and f2 second index
sll_real64, dimension(n1,n2), intent(in) :: f1    !< first field 2d
sll_real64, dimension(n1,n2), intent(in) :: f2    !< second field 2d
sll_int32,                    intent(in) :: iplot !< plot counter
sll_real64,                   intent(in) :: time  !< step time

integer          :: i, j
character(len=4) :: cplot

call int2string(iplot, cplot)

!write domains
open( 80, file = fname//cplot//".dat" )
   do i=1,n1
      do j=1,n2
         write(80,*) i, j, sngl(f1(i,j)), sngl(f2(i,j))
      end do
      write(80,*) 
   end do
close(80)
   
open( 90, file = fname//'plots.gnu', position="append" )
  if ( iplot == 1 ) then
     rewind(90)
     !write(90,*)"set xr[-0.1:1.1]"
     !write(90,*)"set yr[-0.1:1.1]"
     write(90,*)"set zr[-1.1:1.1]"
     !write(90,*)"set cbrange[-1:1]"
     !write(90,*)"set pm3d"
     write(90,*)"set surf"
     write(90,*)"set term x11"
  end if
write(90,*)"set title 'Time = ",time,"'"
write(90,"(a)",advance='no')"splot '"//fname//cplot//".dat' w lines"
write(90,"(a)",advance='no')",'"//fname//cplot//".dat' u 1:2:4 w lines"
close(90)

end subroutine sll_plot_two_fields

end module sll_m_maxwell_solvers_base
