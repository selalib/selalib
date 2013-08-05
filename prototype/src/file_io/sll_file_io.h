#ifndef _SLL_IO_H
#define _SLL_IO_H
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


#define SLL_IO_XDMF    0 
#define SLL_IO_VTK     1
#define SLL_IO_GNUPLOT 2 
#define SLL_IO_MTV     3 

use sll_xml_io
use sll_ascii_io
use sll_binary_io
use sll_gnuplot
use sll_hdf5_io
use sll_xdmf
use sll_plotmtv


#endif
