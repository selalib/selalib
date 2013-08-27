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

module sll_multigrid_base
#include "sll_working_precision.h"
use sll_remapper


type, public  :: sll_multigrid_solver
   type(layout_2D), pointer :: layout
   sll_int32                :: ibdry
   sll_int32                :: jbdry
   sll_real64               :: vbc(4)
   sll_real64               :: phibc(4,20)
   sll_int32                :: comm2d
   sll_real64               :: tolmax
   sll_int32                :: maxcy
   sll_int32                :: kcycle
   sll_int32                :: iprer
   sll_int32                :: ipost
   sll_int32                :: iresw
end type sll_multigrid_solver

end module sll_multigrid_base
