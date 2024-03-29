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

!> @ingroup file_io
!> @brief
!> Module that contains routines to write data in ASCII format file
!> @details
!> This is an example of how use the sll_m_ascii_io module.
!> More details about this example
!> @snippet file_io/unit_test_ascii.F90 example
module sll_m_ascii_io
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   implicit none

   public :: &
      sll_s_ascii_file_close, &
      sll_s_ascii_file_create, &
      sll_o_ascii_write_array, &
      sll_s_ascii_write_array_1d, &
      sll_s_ascii_write_array_2d, &
      sll_s_ascii_write_array_1d_as_row

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Write nD array in ascii file
!> \param[in] file_id file unit number
!> \param[in] array array
!> \param[out] error error code
   interface sll_o_ascii_write_array
      module procedure sll_s_ascii_write_array_1d
      module procedure sll_s_ascii_write_array_2d
      module procedure sll_ascii_write_array_3d
   end interface

contains

!> Create ASCII file
   subroutine sll_s_ascii_file_create(file_name, file_id, error)

      character(len=*), intent(in)         :: file_name !< file name
      sll_int32, intent(out)               :: error     !< error code
      sll_int32, intent(out)               :: file_id   !< file unit number
      logical                              :: lopen

      error = 0

!       do 100 file_id = 20, 99
!
!          inquire (unit=file_id, opened=lopen)
!          if (lopen) then
!             cycle
!          else
!             open (file_id, status='SCRATCH', err=100)
!             close (file_id, status='DELETE', err=100)
!             goto 200
!          end if
!
! 100      continue
!          error = 1
! 200      continue
!          error = 0

         open(newunit=file_id,FILE=file_name,FORM='FORMATTED',IOSTAT=error)
         rewind (file_id)

         end subroutine sll_s_ascii_file_create

!> Close ASCII file
         subroutine sll_s_ascii_file_close(file_id, error)
            sll_int32, intent(in)  :: file_id !<file unit number
            sll_int32, intent(out) :: error   !<error code

            close (file_id, IOSTAT=error)

         end subroutine sll_s_ascii_file_close

!> Write a 1d array ASCII format
         subroutine sll_ascii_write_array_0d(file_id, array, error)
            sll_int32, intent(in)  :: file_id
            sll_int32, intent(out) :: error
            sll_real64, intent(in)  :: array
            write (file_id, *, IOSTAT=error) array
         end subroutine

!> Write a 1d array ASCII format
         subroutine sll_s_ascii_write_array_1d(file_id, array, error, num_points, array2, array3)
            sll_int32, intent(in)  :: file_id
            sll_int32, intent(out) :: error
            sll_real64, intent(in)  :: array(:)
            sll_int32, intent(in), optional :: num_points
            sll_real64, intent(in), optional  :: array2(:)
            sll_real64, intent(in), optional  :: array3(:)
            sll_int32 :: i
            character(len=30) :: rowfmt

            if (.not. present(num_points)) then
               write (file_id, *, IOSTAT=error) array
            else
               if (size(array) < num_points) then
                  print *, '#bad size for array'
                  print *, '#at line/file', __LINE__, __FILE__
                  stop
               end if
               if (present(array2)) then
                  if (size(array2) < num_points) then
                     print *, '#bad size for array2'
                     print *, '#at line/file', __LINE__, __FILE__
                     stop
                  end if
                  if (present(array3)) then
                     if (size(array3) < num_points) then
                        print *, '#bad size for array3'
                        print *, '#at line/file', __LINE__, __FILE__
                        stop
                     end if
                     do i = 1, num_points
                        write (file_id, *, IOSTAT=error) array(i), array2(i), array3(i)
                     end do
                  else
                     do i = 1, num_points
                        write (file_id, *, IOSTAT=error) array(i), array2(i)
                     end do
                  end if
               else
                  write (rowfmt, '(a,i7,a)') '(', num_points, '(1x,g15.5))'
                  do i = 1, num_points
                     write (file_id, trim(rowfmt), IOSTAT=error) array(i)
                  end do

               end if
            end if
         end subroutine

!> Write a 2d array ASCII format
         subroutine sll_s_ascii_write_array_2d(file_id, array, error)
            sll_int32, intent(in)  :: file_id
            sll_int32, intent(out) :: error
            sll_real64, intent(in)  :: array(:, :)
            write (file_id, "(100(1x,g15.5))", IOSTAT=error) array
         end subroutine

!> Write a 3d array ASCII format
         subroutine sll_ascii_write_array_3d(file_id, array, error)
            sll_int32, intent(in)  :: file_id
            sll_int32, intent(out) :: error
            sll_real64, intent(in)  :: array(:, :, :)
            write (file_id, *, IOSTAT=error) array
          end subroutine sll_ascii_write_array_3d

          
          !> Write a 1d array ASCII format
          subroutine sll_s_ascii_write_array_1d_as_row(file_id,array,num_points)
            sll_int32 , intent(in)  :: file_id
            sll_real64, intent(in)  :: array(:)
            sll_int32, intent(in)   :: num_points
            sll_int32 :: i
            
            do i=1,num_points-1
               write(file_id, '(2g24.16)', advance='no') array(i)
            end do
            write(file_id, '(2g24.16)') array(num_points)
            
          end subroutine sll_s_ascii_write_array_1d_as_row

         end module sll_m_ascii_io

