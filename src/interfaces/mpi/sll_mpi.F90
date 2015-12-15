!------------------------------------------------------------------------------
! Selalib
!------------------------------------------------------------------------------
! INTERFACE MODULE: sll_mpi
!
! DESCRIPTION:
!  @ingroup mpi
!! @author Marco Restelli - <marco.restelli@gmail.com>
!! @author Yaman Güçlü    - <yaman.guclu@gmail.com>
!! @brief Interface to MPI.
!!
!! \n
!!
!! @details
!! This is an interface to MPI functions and MPI related utilities.
!! Notice however that MPI specific variables (such as buffers,
!! communicators, tags and so on) must be defined where they are used.
!!
!! One of the main advantages of this module is that it allows an easy
!! switch between the two syntaxes
!! \code
!!   use mpi
!! \endcode
!! and
!! \code
!!   include "mpif.h"
!! \endcode
!! The first syntax should be preferred, but it requires MPI to be
!! compiled with the same compiler used for the application. The
!! second form does not provides the additional checks of fortran 90,
!! but it is the only choice when the 'mpi.mod' file is not available
!! for the chosen compiler.
!!
!! \bug \c mpi_type_create_f90_real has a bug in mpich2, so that one
!! has to use \c mpi_double_precision. The bug is fixed in version 1.3
!! (see https://trac.mcs.anl.gov/projects/mpich2/ticket/1028).
!<----------------------------------------------------------------------
module sll_mpi

!-----------------------------------------------------------------------
! Select here the desired bindings: f77 or f90

  use mpi
  implicit none
!  include "mpif.h"
!  ! These functions are defined in "use mpi" but not in "mpif.h"
!  external :: mpi_init, mpi_init_thread, mpi_initialized, mpi_finalize, &
!    mpi_comm_size, mpi_comm_rank, mpi_comm_split, mpi_comm_free, &
!    mpi_type_create_f90_real, &
!    mpi_barrier, mpi_wait, mpi_waitall 

!-----------------------------------------------------------------------
! Module interface

  public ::              &
    mpi_allgather,       &
    mpi_allgatherv,      &
    mpi_allreduce,       &
    mpi_alltoall,        &
    mpi_alltoallv,       &
    mpi_any_source,      &
    mpi_any_tag,         &
    mpi_barrier,         &
    mpi_bcast,           &
    mpi_cart_coords,     &
    mpi_cart_create,     &
    mpi_cart_get,        &
    mpi_cart_shift,      &
    mpi_character,       &
    mpi_comm_free,       &
    mpi_comm_rank,       &
    mpi_comm_size,       &
    mpi_comm_split,      &
    mpi_comm_world,      &
    mpi_complex,         &
    mpi_dims_create,     &
    mpi_double,          &
    mpi_double_complex,  &
    mpi_double_precision,&
    mpi_finalize,        &
    mpi_gather,          &
    mpi_gatherv,         &
    mpi_get_count,       &
    !mpi_iallreduce,      &
    mpi_in_place,        &
    mpi_info_null,       &
    mpi_init,            &
    mpi_init_thread,     &
    mpi_integer,         &
    mpi_integer8,        &
    mpi_irecv,           &
    mpi_isend,           &
    mpi_land,            &
    mpi_logical,         &
    mpi_lor,             &
    mpi_max,             &
    mpi_proc_null,       &
    mpi_prod,            &
    mpi_real,            &
    mpi_real8,           &
    mpi_recv,            &
    mpi_reduce,          &
    mpi_request_null,    &
    mpi_request_free,    &
    mpi_send,            &
    mpi_sendrecv,        &
    mpi_scatter,         &
    mpi_scatterv,        &
    mpi_source,          &
    mpi_status_ignore,   &
    mpi_status_size,     &
    mpi_success,         &
    mpi_sum,             &
    mpi_thread_funneled, &
    mpi_thread_multiple, &
    mpi_thread_single,   &
    mpi_undefined,       &
    mpi_wait,            &
    mpi_wtime

  private

!-----------------------------------------------------------------------
! These are the subroutines that are not defined in module "mpi"
#ifdef INTEL_MPI
#include "external_intel.F90"
#elif defined(OMPI)
#include "external_openmpi.F90"
#elif defined(MPICH)
#include "external_mpich.F90"
#elif defined(BULLX_MPI)
#include "external_bullx.F90"
#endif

!-----------------------------------------------------------------------
end module sll_mpi
