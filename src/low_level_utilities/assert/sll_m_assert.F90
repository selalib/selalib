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
#ifndef DOXYGEN_SHOULD_SKIP_THIS

module sll_m_assert
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef MPI_VERSION
   use mpi
#endif
   implicit none

   public :: &
      sll_s_assertion

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Instead of using the non-standard subroutine abort() provided by the compiler,
   ! use abort() from the C standard library "stdlib.h"
   interface
      subroutine c_abort() bind(C, name="abort")
      end subroutine
   end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! This routine just takes the information about the assertion error and
   ! writes it on the screen.
   ! This function is only meant to be used by the assert macro. No Doxygen
   ! documentation needed.
   subroutine sll_s_assertion(msg, file, line)
      character(len=*), intent(in) :: msg
      character(len=*), intent(in) :: file
      integer, intent(in) :: line

      write (*, *)
      write (*, '(a)') "ASSERTION FAILURE: condition ( "//trim(msg)//" ) is not satisfied."
      write (*, '(a,i0)') 'Triggered at '//file//':', line

#ifdef MPI_VERSION
      call mpi_abort()
#else
      call c_abort()
#endif

   end subroutine sll_s_assertion

end module sll_m_assert

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
