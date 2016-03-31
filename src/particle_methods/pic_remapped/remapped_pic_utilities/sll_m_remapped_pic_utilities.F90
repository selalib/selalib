!**************************************************************
!  Copyright INRIA
!  Authors : MCP,ALH
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

module sll_m_remapped_pic_utilities

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d

  implicit none

  public :: &
    sll_s_apply_periodic_bc_on_cartesian_mesh_2d, &
    sll_f_is_in_domain_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains



  ! todo: put this in the right module (with the meshes?)
  ! tells whether the given point is in the given domain, with boolean arguments for the domain periodicity
  ! (taken from previous function in_bounds_periodic)
  function sll_f_x__is_in_domain_2d( x, y, mesh, x_periodic, y_periodic ) result(res)

!    use sll_m_cartesian_meshes
    sll_real64,                     intent( in )            :: x, y
    type(sll_t_cartesian_mesh_2d),    intent( in ), pointer   :: mesh
    logical,                        intent( in )            :: x_periodic
    logical,                        intent( in )            :: y_periodic
    logical     :: res

    res = ( x >= mesh%eta1_min )                                                                                    &
          .and.                                                                                                     &
          ( ( x < mesh%eta1_max .and. x_periodic ) .or. ( x <= mesh%eta1_max .and. .not. x_periodic ) )             &
          .and.                                                                                                     &
          ( y >= mesh%eta2_min )                                                                                    &
          .and.                                                                                                     &
          ( ( y < mesh%eta2_max .and. y_periodic ) .or. ( y <= mesh%eta2_max .and. .not. y_periodic) )

  end function sll_f_x_is_in_domain_2d

  ! <<sll_s_apply_periodic_bc_on_cartesian_mesh_2d>>

  ! todo: put this in the right module (with the meshes?)
  subroutine sll_s_apply_periodic_bc_on_cartesian_mesh_2d( mesh, x, y )

!    use sll_m_cartesian_meshes
    ! [[file:../working_precision/sll_m_working_precision.h]]
!    use sll_m_working_precision

    type(sll_t_cartesian_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: x
    sll_real64, intent(inout) :: y

    x = mesh%eta1_min + modulo(x - mesh%eta1_min, mesh%eta1_max - mesh%eta1_min)
    y = mesh%eta2_min + modulo(y - mesh%eta2_min, mesh%eta2_max - mesh%eta2_min)
  end subroutine sll_s_apply_periodic_bc_on_cartesian_mesh_2d

  ! todo: put this in the right module (with the meshes?)
  subroutine sll_s_get_4d_cell_containing_point(eta, grid, i_cell_1, i_cell_2, i_cell_3, i_cell_4, point_is_outside_grid)
    sll_real64, dimension(4),             intent( in )   :: eta     !> point coordinates
    type(sll_t_cartesian_mesh_4d), pointer, intent( in )   :: grid
    sll_int32,                            intent( out)   :: i_cell_1, i_cell_2, i_cell_3, i_cell_4
    logical,                              intent( out )  :: point_is_outside_grid
    sll_real64 ::  tmp

    point_is_outside_grid = .false.

    tmp = ( eta(1) - grid%eta1_min ) / grid%delta_eta1
    i_cell_1 = int( tmp ) + 1
    if( i_cell_1 < 1 .or. i_cell_1 > grid%num_cells1 )then
      point_is_outside_grid = .true.
    end if

    tmp = ( eta(2) - grid%eta2_min ) / grid%delta_eta2
    i_cell_2 = int( tmp ) + 1
    if( i_cell_2 < 1 .or. i_cell_2 > grid%num_cells2 )then
      point_is_outside_grid = .true.
    end if

    tmp = ( eta(3) - grid%eta3_min ) / grid%delta_eta3
    i_cell_3 = int( tmp ) + 1
    if( i_cell_3 < 1 .or. i_cell_3 > grid%num_cells3 )then
      point_is_outside_grid = .true.
    end if

    tmp = ( eta(4) - grid%eta4_min ) / grid%delta_eta4
    i_cell_4 = int( tmp ) + 1
    if( i_cell_4 < 1 .or. i_cell_4 > grid%num_cells4 )then
      point_is_outside_grid = .true.
    end if

  end subroutine sll_s_get_4d_cell_containing_point


  ! ------------------------------------------------------------------------------------------------------------------------
  ! computes the inverse of a matrix with given matrix_size
  !   - here we must have   matrix_size <= 4   and the matrices are 4x4 arrays filled with 0s if matrix_size < 4
  !   - borrows code from David G. Simpson, NASA Goddard Space Flight Center, Greenbelt, Maryland  20771
  subroutine sll_s_get_inverse_matrix_with_given_size(matrix_size, a, inv_a, ok_flag)
    sll_int32,                  intent(in)  :: matrix_size
    sll_real64, dimension(4,4), intent(in)  :: a
    sll_real64, dimension(4,4), intent(out) :: inv_a
    logical,                    intent(out) :: ok_flag

    sll_int32                   :: i,j !,k,l
    sll_real64, parameter       :: epsilon = 1.0d-10
    sll_real64                  :: determinant
    sll_real64                  :: inv_determinant
    sll_real64, dimension(4,4)  :: cofactor

    if( matrix_size == 1 )then

      determinant = a(1,1)

      cofactor(1,1) = 1.0d0

    else if( matrix_size == 2 )then

      determinant =   a(1,1) * a(2,2) - a(1,2) * a(2,1)

      cofactor(1,1) = +a(2,2)
      cofactor(1,2) = -a(2,1)
      cofactor(2,1) = -a(1,2)
      cofactor(2,2) = +a(1,1)

    else if( matrix_size == 3 )then

      determinant =   a(1,1) * a(2,2) * a(3,3)  &
                    - a(1,1) * a(2,3) * a(3,2)  &
                    - a(1,2) * a(2,1) * a(3,3)  &
                    + a(1,2) * a(2,3) * a(3,1)  &
                    + a(1,3) * a(2,1) * a(3,2)  &
                    - a(1,3) * a(2,2) * a(3,1)

      cofactor(1,1) = +(a(2,2) * a(3,3)-a(2,3) * a(3,2))
      cofactor(1,2) = -(a(2,1) * a(3,3)-a(2,3) * a(3,1))
      cofactor(1,3) = +(a(2,1) * a(3,2)-a(2,2) * a(3,1))
      cofactor(2,1) = -(a(1,2) * a(3,3)-a(1,3) * a(3,2))
      cofactor(2,2) = +(a(1,1) * a(3,3)-a(1,3) * a(3,1))
      cofactor(2,3) = -(a(1,1) * a(3,2)-a(1,2) * a(3,1))
      cofactor(3,1) = +(a(1,2) * a(2,3)-a(1,3) * a(2,2))
      cofactor(3,2) = -(a(1,1) * a(2,3)-a(1,3) * a(2,1))
      cofactor(3,3) = +(a(1,1) * a(2,2)-a(1,2) * a(2,1))

    else if( matrix_size == 4 )then

      determinant =  a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)*(a(3,4)*a(4,2)-a(3,2)*a(4,4))+a(2,4)*(a(3,2)*a(4,3)- &
                     a(3,3)*a(4,2)))-a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)*(a(3,4)*a(4,1)-a(3,1)*a(4,4))+ &
                     a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))+a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,2)*(a(3,4)*a(4,1)- &
                     a(3,1)*a(4,4))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))-a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))+ &
                     a(2,2)*(a(3,3)*a(4,1)-a(3,1)*a(4,3))+a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))

      cofactor(1,1) = a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)*(a(3,4)*a(4,2)-a(3,2)*a(4,4))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))
      cofactor(1,2) = a(2,1)*(a(3,4)*a(4,3)-a(3,3)*a(4,4))+a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,3)*a(4,1)-a(3,1)*a(4,3))
      cofactor(1,3) = a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,2)*(a(3,4)*a(4,1)-a(3,1)*a(4,4))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))
      cofactor(1,4) = a(2,1)*(a(3,3)*a(4,2)-a(3,2)*a(4,3))+a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(2,3)*(a(3,2)*a(4,1)-a(3,1)*a(4,2))
      cofactor(2,1) = a(1,2)*(a(3,4)*a(4,3)-a(3,3)*a(4,4))+a(1,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(1,4)*(a(3,3)*a(4,2)-a(3,2)*a(4,3))
      cofactor(2,2) = a(1,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(1,3)*(a(3,4)*a(4,1)-a(3,1)*a(4,4))+a(1,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))
      cofactor(2,3) = a(1,1)*(a(3,4)*a(4,2)-a(3,2)*a(4,4))+a(1,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(1,4)*(a(3,2)*a(4,1)-a(3,1)*a(4,2))
      cofactor(2,4) = a(1,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))+a(1,2)*(a(3,3)*a(4,1)-a(3,1)*a(4,3))+a(1,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))
      cofactor(3,1) = a(1,2)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))+a(1,3)*(a(2,4)*a(4,2)-a(2,2)*a(4,4))+a(1,4)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))
      cofactor(3,2) = a(1,1)*(a(2,4)*a(4,3)-a(2,3)*a(4,4))+a(1,3)*(a(2,1)*a(4,4)-a(2,4)*a(4,1))+a(1,4)*(a(2,3)*a(4,1)-a(2,1)*a(4,3))
      cofactor(3,3) = a(1,1)*(a(2,2)*a(4,4)-a(2,4)*a(4,2))+a(1,2)*(a(2,4)*a(4,1)-a(2,1)*a(4,4))+a(1,4)*(a(2,1)*a(4,2)-a(2,2)*a(4,1))
      cofactor(3,4) = a(1,1)*(a(2,3)*a(4,2)-a(2,2)*a(4,3))+a(1,2)*(a(2,1)*a(4,3)-a(2,3)*a(4,1))+a(1,3)*(a(2,2)*a(4,1)-a(2,1)*a(4,2))
      cofactor(4,1) = a(1,2)*(a(2,4)*a(3,3)-a(2,3)*a(3,4))+a(1,3)*(a(2,2)*a(3,4)-a(2,4)*a(3,2))+a(1,4)*(a(2,3)*a(3,2)-a(2,2)*a(3,3))
      cofactor(4,2) = a(1,1)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))+a(1,3)*(a(2,4)*a(3,1)-a(2,1)*a(3,4))+a(1,4)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      cofactor(4,3) = a(1,1)*(a(2,4)*a(3,2)-a(2,2)*a(3,4))+a(1,2)*(a(2,1)*a(3,4)-a(2,4)*a(3,1))+a(1,4)*(a(2,2)*a(3,1)-a(2,1)*a(3,2))
      cofactor(4,4) = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))+a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

    else
      SLL_ERROR("sll_s_get_inverse_matrix_with_given_size", "incorrect value for matrix_size")
    end if

    if( abs(determinant) .le. epsilon )then
       inv_a = 0.0d0
       print *, "[WARNING --- pbm in routine sll_s_get_inverse_matrix_with_given_size] -- matrix_size, determinant = ", &
                      matrix_size, determinant
       ok_flag = .false.
    else
      inv_determinant = 1./determinant
      inv_a = 0.0d0
      do i = 1, matrix_size
        do j = 1, matrix_size
          inv_a(i,j) = inv_determinant * cofactor(j,i)
        end do
      end do
      ok_flag = .true.
    end if

  end subroutine sll_s_get_inverse_matrix_with_given_size

end module sll_m_remapped_pic_utilities
