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
      sll_o_collective_reduce, &
      sll_t_collective_t, &
      sll_f_get_collective_rank

   use sll_m_pic_visu, only: &
      sll_s_compute_df_cic

   use sll_m_utilities, only: &
      sll_s_int2string

   use sll_m_xdmf, only: &
      sll_s_xdmf_corect2d_nodes

   use mpi, only: &
      mpi_sum

   implicit none

   public :: &
      sll_s_distribution_xdmf_coll

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

!>VisIt readable output for particles density
!>Data file format could be XML, HDF5 or Binary (not fully implemented yet)
   subroutine sll_s_distribution_xdmf_coll(plot_name, x, v, w, &
                                           xmin, xmax, nx, vmin, vmax, nv, iplot, collective, root_rank)

      character(len=*), intent(in) :: plot_name
      sll_real64, dimension(:), intent(in) :: x
      sll_real64, dimension(:), intent(in) :: v
      sll_real64, dimension(:), intent(in) :: w
      sll_real64, intent(in) :: xmin
      sll_real64, intent(in) :: xmax
      sll_int32, intent(in) :: nx
      sll_real64, intent(in) :: vmin
      sll_real64, intent(in) :: vmax
      sll_int32, intent(in) :: nv
      sll_int32, intent(in) :: iplot
      type(sll_t_collective_t), pointer    :: collective
      sll_int32, optional, intent(in) :: root_rank

      sll_real64, dimension(nx, nv) :: df_local
      sll_real64, dimension(nx, nv) :: df
      sll_real64                   :: delta_x
      sll_real64                   :: delta_v
      character(len=4)             :: fin

      call sll_s_int2string(iplot, fin)

      delta_x = (xmax - xmin)/real(nx - 1, f64)
      delta_v = (vmax - vmin)/real(nv - 1, f64)

      call sll_s_compute_df_cic(x, v, w, xmin, xmax, nx, vmin, vmax, nv, df_local)

      if (present(root_rank)) then

         !Write the result only to node with root_rank
         call sll_o_collective_reduce(collective, df_local, nx*nv, MPI_SUM, root_rank, df)

         if (sll_f_get_collective_rank(collective) == root_rank) then
            call sll_s_xdmf_corect2d_nodes(plot_name//'_'//fin, df, "density", xmin, delta_x, vmin, delta_v)
         end if

      else

         !Write the result to all nodes, and use parallel io

      end if

   end subroutine sll_s_distribution_xdmf_coll

end module sll_m_pic_visu_parallel
