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

!> \brief
!> This module provides some routines for plotting during PIC simulations with MPI.
!> It extends sll_m_pic_visu with simple MPI functionality.

module sll_m_pic_visu_parallel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_collective_reduce, &
    sll_collective_t, &
    sll_get_collective_rank

  use sll_m_pic_visu, only: &
    compute_df_cic

  use sll_m_utilities, only: &
    int2string

  use sll_m_xdmf, only: &
    sll_xdmf_corect2d_nodes

  use sll_mpi, only: &
    mpi_sum

  implicit none

  public :: &
    distribution_xdmf_coll

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

!>VisIt readable output for particles density
!>Data file format could be XML, HDF5 or Binary (not fully implemented yet)
subroutine distribution_xdmf_coll(plot_name, x, v, w, &
                             xmin, xmax, nx,     &
                             vmin, vmax, nv, iplot, collective,root_rank)
character(len=*), intent(in) :: plot_name
sll_real64, dimension(:), intent(in) :: x
sll_real64, dimension(:), intent(in) :: v
sll_real64, dimension(:), intent(in) :: w
sll_int32, intent(in) :: nx
sll_int32, intent(in) :: nv
sll_int32 :: iplot
sll_real64, dimension(nx,nv) :: df_local,df
sll_real64 :: delta_x, delta_v, xmin, xmax, vmin, vmax
character(len=4) :: fin
type(sll_collective_t), pointer      :: collective
sll_int32, intent(in), optional                :: root_rank

call int2string(iplot, fin)


delta_x = (xmax-xmin)/(nx-1)
delta_v = (vmax-vmin)/(nv-1)

call compute_df_cic(x, v, w, xmin, xmax, nx, vmin, vmax, nv, df_local)

!Mpi reduce
   if (present(root_rank)) then
         !Write the result only to node with root_rank
         call sll_collective_reduce( collective, df_local, nx*nv, MPI_SUM,root_rank, df )

  if ( sll_get_collective_rank( collective ) == root_rank) then
    call sll_xdmf_corect2d_nodes( plot_name//'_'//fin, df, "density", xmin, delta_x, vmin, delta_v)
  endif       
     else
       !Write the result to all nodes, and use parallel io
       
     endif
end subroutine distribution_xdmf_coll



end module sll_m_pic_visu_parallel
