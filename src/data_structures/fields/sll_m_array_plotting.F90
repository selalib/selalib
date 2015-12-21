module sll_m_array_plotting
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_4d

  use sll_m_gnuplot, only: &
    sll_o_gnuplot_2d

  implicit none

  public :: &
    sll_p_x1x2, &
    sll_p_x1x3, &
    sll_p_x1x4, &
    sll_p_x2x3, &
    sll_p_x2x4, &
    sll_p_x3x4, &
    sll_s_write_projection_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_int32, parameter :: sll_p_x1x2 = 0
sll_int32, parameter :: sll_p_x1x3 = 1
sll_int32, parameter :: sll_p_x1x4 = 2
sll_int32, parameter :: sll_p_x2x3 = 3
sll_int32, parameter :: sll_p_x2x4 = 4
sll_int32, parameter :: sll_p_x3x4 = 5

contains

!> @brief
!> Write a gnuplot file to display 2d projection of 4d array
!> @details
!> Set the projection, possible values are
!> - sll_p_x1x2
!> - sll_p_x1x3
!> - sll_p_x1x4
!> - sll_p_x2x3
!> - sll_p_x2x4
!> - sll_p_x3x4
!> And set the slice position (integer array with dimension=2)
subroutine sll_s_write_projection_2d( mesh4d,     &
                                array,      &
                                label,      &
                                projection, &
                                slice,      &
                                iplot )

type(sll_t_cartesian_mesh_4d), intent(in) :: mesh4d          !< cartesian mesh
character(len=*),            intent(in) :: label           !< file name
sll_real64,                  intent(in) :: array(:,:,:,:)  !< data array
sll_int32,                   intent(in) :: projection      !< projection plan
sll_int32,                   intent(in) :: slice(2)        !< position of slice
sll_int32,                   intent(in) :: iplot           !< plot number

sll_real64, pointer :: f_split(:,:)

sll_real64 :: eta1_min, eta2_min, eta3_min, eta4_min
sll_real64 :: eta1_max, eta2_max, eta3_max, eta4_max
sll_int32  :: n1, n2, n3, n4
sll_int32  :: error

n1 = mesh4d%num_cells1+1
n2 = mesh4d%num_cells2+1
n3 = mesh4d%num_cells3+1
n4 = mesh4d%num_cells4+1

eta1_min = mesh4d%eta1_min ; eta1_max = mesh4d%eta1_max
eta2_min = mesh4d%eta2_min ; eta2_max = mesh4d%eta2_max
eta3_min = mesh4d%eta3_min ; eta3_max = mesh4d%eta3_max
eta4_min = mesh4d%eta4_min ; eta4_max = mesh4d%eta4_max

select case(projection)

case(sll_p_x1x2)

  SLL_ALLOCATE(f_split(n1,n2),error)
  f_split = array(:,:,slice(1),slice(2))
  call sll_o_gnuplot_2d(eta1_min, eta1_max, n1,       &
                      eta2_min, eta2_max, n2,       &
                      f_split,                      &
                      label, iplot, error)
case(sll_p_x1x3)

  SLL_ALLOCATE(f_split(n1,n3),error)
  f_split = array(:,slice(1),:,slice(2))
  call sll_o_gnuplot_2d(eta1_min, eta1_max, n1,       &
                      eta3_min, eta3_max, n3,       &
                      f_split,                      &
                      label, iplot, error)
case(sll_p_x1x4)

  SLL_ALLOCATE(f_split(n1,n4),error)
  f_split = array(:,slice(1),slice(2),:)
  call sll_o_gnuplot_2d(eta1_min, eta1_max, n1,       &
                      eta4_min, eta4_max, n4,       &
                      f_split,                      &
                      label, iplot, error)
case(sll_p_x2x3)

  SLL_ALLOCATE(f_split(n2,n3),error)
  f_split = array(slice(1),:,:,slice(2))
  call sll_o_gnuplot_2d(eta2_min, eta2_max, n2,       &
                      eta3_min, eta3_max, n3,       &
                      f_split,                      &
                      label, iplot, error)
case(sll_p_x2x4)

  SLL_ALLOCATE(f_split(n2,n4),error)
  f_split = array(slice(1),:,slice(2),:)
  call sll_o_gnuplot_2d(eta2_min, eta2_max, n2,       &
                      eta4_min, eta4_max, n4,       &
                      f_split,                      &
                      label, iplot, error)
case(sll_p_x3x4)

  SLL_ALLOCATE(f_split(n3,n4),error)
  f_split = array(slice(1),slice(2),:,:)
  call sll_o_gnuplot_2d(eta3_min, eta3_max, n3,       &
                      eta4_min, eta4_max, n4,       &
                      f_split,                      &
                      label, iplot, error)

end select


end subroutine sll_s_write_projection_2d

end module sll_m_array_plotting
