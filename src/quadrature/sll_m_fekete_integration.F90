!> @ingroup integration
!> @author Laura S. Mendoza
!> @brief Fekete quadrature rules for a triangle
!> @details
!> This module contains the Fekete quadrature rule adapted to triangles.
!> The main functions are taken from the following site but have been
!> modified to respect Selalib's structure:
!> http://people.sc.fsu.edu/~jburkardt/f_src/triangle_fekete_rule/
module sll_m_fekete_integration

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_m_hexagonal_meshes
  use sll_m_box_splines, only : &
       write_connectivity, &
       boxspline_val_der
  implicit none

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  abstract interface
     !> 2d real function
     function function_2D(x, y)
       use sll_m_working_precision ! can't pass a header file because the
       ! preprocessor prevents double inclusion.
       ! This is very rare.
       sll_real64             :: function_2D
       sll_real64, intent(in) :: x
       sll_real64, intent(in) :: y
     end function function_2D
  end interface
#endif


contains


  !--------------------------------------------------------------------
  !> @brief returns the degree of a Fekete rule for the triangle.
  !> @details This code was first written by John Burkardt and is available
  !> online under the GNU LGPL license.
  !> Reference: Mark Taylor, Beth Wingate, Rachel Vincent,
  !> An Algorithm for Computing Fekete Points in the Triangle,
  !> SIAM Journal on Numerical Analysis,
  !> Volume 38, Number 5, 2000, pages 1707-1720.
  !> @param[IN] rule integer rule of quadrature for the fekete quadrature
  !> @param[OUT] degree integer the polynomial degree of exactness of
  subroutine fekete_degree(rule, degree)
    sll_int32, intent(in)  ::  rule
    sll_int32, intent(out) ::  degree

    if(rule == 1) then
       degree = 3
    else if(rule == 2) then
       degree = 6
    else if(rule == 3) then
       degree = 9
    else if(rule == 4) then
       degree = 12
    else if(rule == 5) then
       degree = 12
    else if(rule == 6) then
       degree = 15
    else if(rule == 7) then
       degree = 18
    else
       degree = -1
       print *, ""
       print *, "In fekete_degree():  Fatal error"
       print *, "                     Illegal RULE = ", rule
       STOP
    end if
  end subroutine fekete_degree


  !----------------------------------------------------------------------
  !> @brief returns the order of a Fekete rule for the triangle.
  !> @details This code was first written by John Burkardt and is available
  !> online under the GNU LGPL license.
  !>  Reference:
  !>   Mark Taylor, Beth Wingate, Rachel Vincent,
  !>   An Algorithm for Computing Fekete Points in the Triangle,
  !>   SIAM Journal on Numerical Analysis,
  !>   Volume 38, Number 5, 2000, pages 1707-1720.
  !> @param[IN] rule integer the index of the rule
  !> @param[OUT] order_num integer the order (number of points) of the rule
  subroutine fekete_order_num(rule, order_num)
    sll_int32, intent(in)  ::  rule
    sll_int32, intent(out) ::  order_num
    sll_int32, allocatable, dimension(:) :: suborder
    sll_int32 ::  suborder_num
    sll_int32 ::  ierr

    call fekete_suborder_num(rule, suborder_num)

    SLL_ALLOCATE(suborder(1:suborder_num), ierr)

    call fekete_suborder(rule, suborder_num, suborder)

    order_num = sum (suborder(1:suborder_num))

    SLL_DEALLOCATE_ARRAY(suborder, ierr)
  end subroutine fekete_order_num

  !---------------------------------------------------------------------------
  !> @brief returns the points and weights of a Fekete rule.
  !> @details This code was first written by John Burkardt and is available
  !> online under the GNU LGPL license.
  !> @param[IN] rule integer the index of the rule
  !> @param[IN] order_num integer the order (number of points of the rule)
  !> @param[OUT] xy real table the points of the rule
  !> @param[OUT] w real table the weigths of the rule
  subroutine fekete_rule(rule, order_num, xy, w)
    sll_int32,  intent(in)  :: rule
    sll_int32,  intent(in)  :: order_num
    sll_real64, intent(out) :: w(order_num)
    sll_real64, intent(out) :: xy(2,order_num)
    sll_int32,  allocatable, dimension (:) :: suborder
    sll_real64, allocatable, dimension (:) :: suborder_w
    sll_real64, allocatable, dimension (:,:) :: suborder_xyz
    sll_int32 ::  suborder_num
    sll_int32 ::  k
    sll_int32 ::  o
    sll_int32 ::  s
    sll_int32 ::  ierr

    !  Get the suborder information.
    call fekete_suborder_num(rule, suborder_num)

    SLL_ALLOCATE(suborder(suborder_num), ierr)
    SLL_ALLOCATE(suborder_xyz(3,suborder_num), ierr)
    SLL_ALLOCATE(suborder_w(suborder_num), ierr)

    call fekete_suborder(rule, suborder_num, suborder)

    call fekete_subrule(rule, suborder_num, suborder_xyz, suborder_w)

    !  Expand the suborder information to a full order rule.
    o = 0
    do s = 1, suborder_num
       if(suborder(s) == 1) then
          o = o + 1
          xy(1:2,o) = suborder_xyz(1:2,s)
          w(o) = suborder_w(s)
       else if(suborder(s) == 3) then
          do k = 1, 3
             o = o + 1
             xy(1,o) = suborder_xyz(wrapping(k,  1,3), s)
             xy(2,o) = suborder_xyz(wrapping(k+1,1,3), s)
             w(o) = suborder_w(s)
          end do
       else if(suborder(s) == 6) then
          do k = 1, 3
             o = o + 1
             xy(1,o) = suborder_xyz(wrapping(k,  1,3), s)
             xy(2,o) = suborder_xyz(wrapping(k+1,1,3), s)
             w(o) = suborder_w(s)
          end do
          do k = 1, 3
             o = o + 1
             xy(1,o) = suborder_xyz(wrapping(k+1,1,3), s)
             xy(2,o) = suborder_xyz(wrapping(k,  1,3), s)
             w(o) = suborder_w(s)
          end do
       else
          print *, ""
          print *, 'FEKETE_RULE - Fatal error!'
          print *, '  Illegal SUBORDER(', s, ') = ', suborder(s)
          STOP
       end if
    end do

    SLL_DEALLOCATE_ARRAY(suborder, ierr)
    SLL_DEALLOCATE_ARRAY(suborder_xyz, ierr)
    SLL_DEALLOCATE_ARRAY(suborder_w, ierr)
  end subroutine fekete_rule

  !---------------------------------------------------------------------------
  !> @brief Re arranges the order of the quadrature point to go from lower-left
  !> to top right.
  !> @param[IN] n integer number of quadrature points
  !> @param[INOUT] xy real table the points of the rule
  !> @param[INOUT] w real table the weigths of the rule
  subroutine rearrange_fekete_rule(n, xy, w)
    sll_int32,  intent(in)  :: n
    sll_real64, intent(out) :: xy(2,n)
    sll_real64, intent(out) :: w(n)
    sll_real64, dimension(3, n) :: xyw_copy
    sll_int32  :: i
    sll_int32  :: low_left
    !sll_real64 :: min_x

    ! we start by making a copy of the points
    xyw_copy(1:2, 1:n) = xy(1:2, 1:n)
    xyw_copy(  3, 1:n) = w(1:n)
    ! now we arrange by the x-coordinate
    do i = 1, n
       low_left = MINLOC(xyw_copy(1,:),1)
       xy(1:2, i) = xyw_copy(1:2, low_left)
       w(i) = xyw_copy(3, low_left)
       ! We put a value that we know will never be the lower left point
       ! as all points are between 0 and 1
       xyw_copy(1:3, low_left) = 10000._f64
    end do

    ! we copy again the points
    xyw_copy(1:2, 1:n) = xy(1:2, 1:n)
    xyw_copy(  3, 1:n) = w(1:n)
    ! now we arrange by the y-coordinate
    do i = 1, n
       low_left = MINLOC(xyw_copy(2,:), 1)
       xy(1:2, i) = xyw_copy(1:2, low_left)
       w(i) = xyw_copy(3, low_left)
       ! We put a value that we know will never be the lower left point
       ! as all points are between 0 and 1
       xyw_copy(1:3, low_left) = 10000._f64
    end do
  end subroutine rearrange_fekete_rule

  !---------------------------------------------------------------------------
  !> @brief returns the number of Fekete rules available.
  !> @details This code was first written by John Burkardt and is available
  !> online under the GNU LGPL license.
  !> @param[OUT] rule_num integer the number of rules available
  subroutine fekete_rule_num(rule_num)
    sll_int32, intent(out) ::  rule_num

    rule_num = 7
  end subroutine fekete_rule_num


  !---------------------------------------------------------------------------
  !> @brief returns the suborders for a Fekete rule.
  !> @details This code was first written by John Burkardt and is available
  !> online under the GNU LGPL license.
  !> @param[IN] rule integer the index of the rule
  !> @param[IN] suborder_num integer the number of suborder of the rule
  !> @param[OUT] suborder integer array the suborders of the rule
  subroutine fekete_suborder(rule, suborder_num, suborder)
    sll_int32, intent(in)  :: suborder_num
    sll_int32, intent(in)  :: rule
    sll_int32, intent(out) :: suborder(suborder_num)

    if(rule == 1) then
       suborder(1:suborder_num) = (/ &
            1, 3, 6 /)
    else if(rule == 2) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 6, 6, 6 /)
    else if(rule == 3) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 3, 6, 6, 6, 6, 6, &
            6, 6 /)
    else if(rule == 4) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 3, 3, 3, 6, 6, 6, &
            6, 6, 6, 6, 6, 6, 6, 6, 6 /)
    else if(rule == 5) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
            3, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
            6  /)
    else if(rule == 6) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
            6, 6, 6, 6, 6, 6, 6, 6  /)
    else if(rule == 7) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
            3, 3, 6, 6, 6, 6, 6, 6, 6, 6, &
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
            6, 6, 6, 6, 6, 6, 6, 6  /)
    else
       print *, ""
       print *,  'FEKETE_SUBORDER - Fatal error!'
       write(*, '(a,i8)') '  Illegal RULE = ', rule
       STOP
    end if
  end subroutine fekete_suborder

  !---------------------------------------------------------------------------
  !> @brief returns the number of suborders for a Fekete rule.
  !> @details This code was first written by John Burkardt and is available
  !> online under the GNU LGPL license.
  !> @param[IN] rule integer the index of the rule
  !> @param[OUTPUT] suborder_num the number of suborder of the rule
  subroutine fekete_suborder_num(rule, suborder_num)
    sll_int32, intent(in)  ::  rule
    sll_int32, intent(out) ::  suborder_num

    if(rule == 1) then
       suborder_num = 3
    else if(rule == 2) then
       suborder_num = 7
    else if(rule == 3) then
       suborder_num = 12
    else if(rule == 4) then
       suborder_num = 19
    else if(rule == 5) then
       suborder_num = 21
    else if(rule == 6) then
       suborder_num = 28
    else if(rule == 7) then
       suborder_num = 38
    else
       suborder_num = -1
       print *, ""
       print *,  'FEKETE_SUBORDER_NUM - Fatal error!'
       write(*, '(a,i8)') '  Illegal RULE = ', rule
       STOP
    end if
  end subroutine fekete_suborder_num

  !---------------------------------------------------------------------------
  !> @brief returns a compressed Fekete rule
  !> @details This code was first written by John Burkardt and is available
  !> online under the GNU LGPL license.
  !> Discussion:
  !>    The listed weights are twice what we want...since we want them
  !>    to sum to 1/2, reflecting the area of a unit triangle.  So we
  !>    simple halve the values before exiting this routine.
  !> @param[IN]  rule integer the index of the rule
  !> @param[IN]  suborder_num the number of suborders
  !> @param[OUT] suborder_xyz(3, suborder_num) real the number of suborders
  !> of the rule
  !> @param[OUT] suborder_w(suborder_num) real the suborder weights
  subroutine fekete_subrule(rule, suborder_num, suborder_xyz, suborder_w)
    sll_int32, intent(in) ::  suborder_num
    sll_int32, intent(in) ::  rule
    sll_real64, intent(out) :: suborder_w(suborder_num)
    sll_real64, intent(out) :: suborder_xyz(3,suborder_num)
    sll_int32 ::  s

    if(rule == 1) then
       suborder_xyz(1:3,1:suborder_num) = reshape((/ &
            1._f64/3._f64,  1._f64/3._f64, 1._f64/3._f64, &
            1.0000000000_f64,  0.0000000000_f64, 0.0000000000_f64, &
            0.0000000000_f64,  0.2763932023_f64, 0.7236067977_f64  &
            /), (/ 3, suborder_num /))

       suborder_w(1:suborder_num) = (/ &
            0.9000000000_f64, &
            0.0333333333_f64, &
            0.1666666667_f64 /)

    else if(rule == 2) then

       suborder_xyz(1:3,1:suborder_num) = reshape((/ &
            1._f64/3._f64,  1._f64/3._f64, 0.3333333334_f64, &
            0.1063354684_f64,  0.1063354684_f64, 0.7873290632_f64, &
            0.5000000000_f64,  0.5000000000_f64, 0.0000000000_f64, &
            1.0000000000_f64,  0.0000000000_f64, 0.0000000000_f64, &
            0.1171809171_f64,  0.3162697959_f64, 0.5665492870_f64, &
            0.0000000000_f64,  0.2655651402_f64, 0.7344348598_f64, &
            0.0000000000_f64,  0.0848854223_f64, 0.9151145777_f64  &
            /), (/ 3, suborder_num /))

       suborder_w(1:suborder_num) = (/ &
            0.2178563571_f64, &
            0.1104193374_f64, &
            0.0358939762_f64, &
            0.0004021278_f64, &
            0.1771348660_f64, &
            0.0272344079_f64, &
            0.0192969460_f64 /)

    else if(rule == 3) then

       suborder_xyz(1:3,1:suborder_num) = reshape((/ &
            1._f64/3._f64,  1._f64/3._f64, 0.3333333334_f64, &
            0.1704318201_f64,  0.1704318201_f64, 0.6591363598_f64, &
            0.0600824712_f64,  0.4699587644_f64, 0.4699587644_f64, &
            0.0489345696_f64,  0.0489345696_f64, 0.9021308608_f64, &
            0.0000000000_f64,  0.0000000000_f64, 1.0000000000_f64, &
            0.1784337588_f64,  0.3252434900_f64, 0.4963227512_f64, &
            0.0588564879_f64,  0.3010242110_f64, 0.6401193011_f64, &
            0.0551758079_f64,  0.1543901944_f64, 0.7904339977_f64, &
            0.0000000000_f64,  0.4173602935_f64, 0.5826397065_f64, &
            0.0000000000_f64,  0.2610371960_f64, 0.7389628040_f64, &
            0.0000000000_f64,  0.1306129092_f64, 0.8693870908_f64, &
            0.0000000000_f64,  0.0402330070_f64, 0.9597669930_f64  &
            /), (/ 3, suborder_num /))

       suborder_w(1:suborder_num) = (/ &
            0.1096011288_f64, &
            0.0767491008_f64, &
            0.0646677819_f64, &
            0.0276211659_f64, &
            0.0013925011_f64, &
            0.0933486453_f64, &
            0.0619010169_f64, &
            0.0437466450_f64, &
            0.0114553907_f64, &
            0.0093115568_f64, &
            0.0078421987_f64, &
            0.0022457501_f64 /)

    else if(rule == 4) then

       suborder_xyz(1:3,1:suborder_num) = reshape((/ &
            1._f64/3._f64,  1._f64/3._f64, 0.3333333334_f64, &
            0.1988883477_f64,  0.4005558262_f64, 0.4005558261_f64, &
            0.2618405201_f64,  0.2618405201_f64, 0.4763189598_f64, &
            0.0807386775_f64,  0.0807386775_f64, 0.8385226450_f64, &
            0.0336975736_f64,  0.0336975736_f64, 0.9326048528_f64, &
            0.0000000000_f64,  0.5000000000_f64, 0.5000000000_f64, &
            0.0000000000_f64,  0.0000000000_f64, 1.0000000000_f64, &
            0.1089969290_f64,  0.3837518758_f64, 0.5072511952_f64, &
            0.1590834479_f64,  0.2454317980_f64, 0.5954847541_f64, &
            0.0887037176_f64,  0.1697134458_f64, 0.7415828366_f64, &
            0.0302317829_f64,  0.4071849276_f64, 0.5625832895_f64, &
            0.0748751152_f64,  0.2874821712_f64, 0.6376427136_f64, &
            0.0250122615_f64,  0.2489279690_f64, 0.7260597695_f64, &
            0.0262645218_f64,  0.1206826354_f64, 0.8530528428_f64, &
            0.0000000000_f64,  0.3753565349_f64, 0.6246434651_f64, &
            0.0000000000_f64,  0.2585450895_f64, 0.7414549105_f64, &
            0.0000000000_f64,  0.1569057655_f64, 0.8430942345_f64, &
            0.0000000000_f64,  0.0768262177_f64, 0.9231737823_f64, &
            0.0000000000_f64,  0.0233450767_f64, 0.9766549233_f64  &
            /), (/ 3, suborder_num /))

       suborder_w(1:suborder_num) = (/ &
            0.0626245179_f64, &
            0.0571359417_f64, &
            0.0545982307_f64, &
            0.0172630326_f64, &
            0.0142519606_f64, &
            0.0030868485_f64, &
            0.0004270742_f64, &
            0.0455876390_f64, &
            0.0496701966_f64, &
            0.0387998322_f64, &
            0.0335323983_f64, &
            0.0268431561_f64, &
            0.0237377452_f64, &
            0.0177255972_f64, &
            0.0043097313_f64, &
            0.0028258057_f64, &
            0.0030994935_f64, &
            0.0023829062_f64, &
            0.0009998683_f64 /)

    else if(rule == 5) then

       suborder_xyz(1:3,1:suborder_num) = reshape((/ &
            1._f64/3._f64,  1._f64/3._f64, 0.3333333334_f64, &
            0.2201371125_f64,  0.3169406831_f64, 0.4629222044_f64, &
            0.2201371125_f64,  0.4629222044_f64, 0.3169406831_f64, &
            0.1877171129_f64,  0.1877171129_f64, 0.6245657742_f64, &
            0.1403402144_f64,  0.4298298928_f64, 0.4298298928_f64, &
            0.0833252778_f64,  0.0833252778_f64, 0.8333494444_f64, &
            0.0664674598_f64,  0.0252297247_f64, 0.9083028155_f64, &
            0.0218884020_f64,  0.4890557990_f64, 0.4890557990_f64, &
            0.0252297247_f64,  0.0664674598_f64, 0.9083028155_f64, &
            0.0000000000_f64,  0.5000000000_f64, 0.5000000000_f64, &
            0.0000000000_f64,  0.0000000000_f64, 1.0000000000_f64, &
            0.1157463404_f64,  0.2842319093_f64, 0.6000217503_f64, &
            0.0672850606_f64,  0.3971764400_f64, 0.5355384994_f64, &
            0.0909839531_f64,  0.1779000668_f64, 0.7311159801_f64, &
            0.0318311633_f64,  0.3025963402_f64, 0.6655724965_f64, &
            0.0273518579_f64,  0.1733665506_f64, 0.7992815915_f64, &
            0.0000000000_f64,  0.3753565349_f64, 0.6246434651_f64, &
            0.0000000000_f64,  0.2585450895_f64, 0.7414549105_f64, &
            0.0000000000_f64,  0.1569057655_f64, 0.8430942345_f64, &
            0.0000000000_f64,  0.0768262177_f64, 0.9231737823_f64, &
            0.0000000000_f64,  0.0233450767_f64, 0.9766549233_f64  &
            /), (/ 3, suborder_num /))

       suborder_w(1:suborder_num) = (/ &
            0.0485965670_f64, &
            0.0602711576_f64, &
            0.0602711576_f64, &
            0.0476929767_f64, &
            0.0453940802_f64, &
            0.0258019417_f64, &
            0.0122004614_f64, &
            0.0230003812_f64, &
            0.0122004614_f64, &
            0.0018106475_f64, &
            -0.0006601747_f64, &
            0.0455413513_f64, &
            0.0334182802_f64, &
            0.0324896773_f64, &
            0.0299402736_f64, &
            0.0233477738_f64, &
            0.0065962854_f64, &
            0.0021485117_f64, &
            0.0034785755_f64, &
            0.0013990566_f64, &
            0.0028825748_f64 /)

    else if(rule == 6) then

       suborder_xyz(1:3,1:suborder_num) = reshape((/ &
            1._f64/3._f64,  1._f64/3._f64, 0.3333333334_f64, &
            0.2379370518_f64,  0.3270403780_f64, 0.4350225702_f64, &
            0.3270403780_f64,  0.2379370518_f64, 0.4350225702_f64, &
            0.1586078048_f64,  0.4206960976_f64, 0.4206960976_f64, &
            0.2260541354_f64,  0.2260541354_f64, 0.5478917292_f64, &
            0.1186657611_f64,  0.1186657611_f64, 0.7626684778_f64, &
            0.0477095725_f64,  0.4761452137_f64, 0.4761452138_f64, &
            0.0531173538_f64,  0.0531173538_f64, 0.8937652924_f64, &
            0.0219495841_f64,  0.0219495841_f64, 0.9561008318_f64, &
            0.0000000000_f64,  0.0000000000_f64, 1.0000000000_f64, &
            0.1585345951_f64,  0.3013819154_f64, 0.5400834895_f64, &
            0.0972525649_f64,  0.3853507643_f64, 0.5173966708_f64, &
            0.0875150140_f64,  0.2749910734_f64, 0.6374939126_f64, &
            0.1339547708_f64,  0.1975591066_f64, 0.6684861226_f64, &
            0.0475622627_f64,  0.3524012205_f64, 0.6000365168_f64, &
            0.0596194677_f64,  0.1978887556_f64, 0.7424917767_f64, &
            0.0534939782_f64,  0.1162464503_f64, 0.8302595715_f64, &
            0.0157189888_f64,  0.4176001732_f64, 0.5666808380_f64, &
            0.0196887324_f64,  0.2844332752_f64, 0.6958779924_f64, &
            0.0180698489_f64,  0.1759511193_f64, 0.8059790318_f64, &
            0.0171941515_f64,  0.0816639421_f64, 0.9011419064_f64, &
            0.0000000000_f64,  0.4493368632_f64, 0.5506631368_f64, &
            0.0000000000_f64,  0.3500847655_f64, 0.6499152345_f64, &
            0.0000000000_f64,  0.2569702891_f64, 0.7430297109_f64, &
            0.0000000000_f64,  0.1738056486_f64, 0.8261943514_f64, &
            0.0000000000_f64,  0.1039958541_f64, 0.8960041459_f64, &
            0.0000000000_f64,  0.0503997335_f64, 0.9496002665_f64, &
            0.0000000000_f64,  0.0152159769_f64, 0.9847840231_f64  &
            /), (/ 3, suborder_num /))

       suborder_w(1:suborder_num) = (/ &
            0.0459710878_f64, &
            0.0346650571_f64, &
            0.0346650571_f64, &
            0.0384470625_f64, &
            0.0386013566_f64, &
            0.0224308157_f64, &
            0.0243531004_f64, &
            0.0094392654_f64, &
            0.0061105652_f64, &
            0.0001283162_f64, &
            0.0305412307_f64, &
            0.0262101254_f64, &
            0.0265367617_f64, &
            0.0269859772_f64, &
            0.0172635676_f64, &
            0.0188795851_f64, &
            0.0158224870_f64, &
            0.0127170850_f64, &
            0.0164489660_f64, &
            0.0120018620_f64, &
            0.0072268907_f64, &
            0.0023599161_f64, &
            0.0017624674_f64, &
            0.0018648017_f64, &
            0.0012975716_f64, &
            0.0018506035_f64, &
            0.0009919379_f64, &
            0.0004893506_f64 /)

    else if(rule == 7) then

       suborder_xyz(1:3,1:suborder_num) = reshape((/ &
            1._f64/3._f64,  1._f64/3._f64, 0.3333333334_f64, &
            0.2515553103_f64,  0.3292984162_f64, 0.4191462735_f64, &
            0.3292984162_f64,  0.2515553103_f64, 0.4191462735_f64, &
            0.1801930996_f64,  0.4099034502_f64, 0.4099034502_f64, &
            0.2438647767_f64,  0.2438647767_f64, 0.5122704466_f64, &
            0.1512564554_f64,  0.1512564554_f64, 0.6974870892_f64, &
            0.0810689493_f64,  0.4594655253_f64, 0.4594655254_f64, &
            0.0832757649_f64,  0.0832757649_f64, 0.8334484702_f64, &
            0.0369065587_f64,  0.0369065587_f64, 0.9261868826_f64, &
            0.0149574850_f64,  0.0149574850_f64, 0.9700850300_f64, &
            0.0000000000_f64,  0.5000000000_f64, 0.5000000000_f64, &
            0.0000000000_f64,  0.0000000000_f64, 1.0000000000_f64, &
            0.1821465920_f64,  0.3095465041_f64, 0.5083069039_f64, &
            0.1246901255_f64,  0.3789288931_f64, 0.4963809814_f64, &
            0.1179441386_f64,  0.2868915642_f64, 0.5951642972_f64, &
            0.1639418454_f64,  0.2204868669_f64, 0.6155712877_f64, &
            0.0742549663_f64,  0.3532533654_f64, 0.5724916683_f64, &
            0.0937816771_f64,  0.2191980979_f64, 0.6870202250_f64, &
            0.0890951387_f64,  0.1446273457_f64, 0.7662775156_f64, &
            0.0409065243_f64,  0.4360543636_f64, 0.5230391121_f64, &
            0.0488675890_f64,  0.2795984854_f64, 0.6715339256_f64, &
            0.0460342127_f64,  0.2034211147_f64, 0.7505446726_f64, &
            0.0420687187_f64,  0.1359040280_f64, 0.8220272533_f64, &
            0.0116377940_f64,  0.4336892286_f64, 0.5546729774_f64, &
            0.0299062187_f64,  0.3585587824_f64, 0.6115349989_f64, &
            0.0132313129_f64,  0.2968103667_f64, 0.6899583204_f64, &
            0.0136098469_f64,  0.2050279257_f64, 0.7813622274_f64, &
            0.0124869684_f64,  0.1232146223_f64, 0.8642984093_f64, &
            0.0365197797_f64,  0.0805854893_f64, 0.8828947310_f64, &
            0.0118637765_f64,  0.0554881302_f64, 0.9326480933_f64, &
            0.0000000000_f64,  0.4154069883_f64, 0.5845930117_f64, &
            0.0000000000_f64,  0.3332475761_f64, 0.6667524239_f64, &
            0.0000000000_f64,  0.2558853572_f64, 0.7441146428_f64, &
            0.0000000000_f64,  0.1855459314_f64, 0.8144540686_f64, &
            0.0000000000_f64,  0.1242528987_f64, 0.8757471013_f64, &
            0.0000000000_f64,  0.0737697111_f64, 0.9262302889_f64, &
            0.0000000000_f64,  0.0355492359_f64, 0.9644507641_f64, &
            0.0000000000_f64,  0.0106941169_f64, 0.9893058831_f64  &
            /), (/ 3, suborder_num /))

       suborder_w(1:suborder_num) = (/ &
            0.0326079297_f64, &
            0.0255331366_f64, &
            0.0255331366_f64, &
            0.0288093886_f64, &
            0.0279490452_f64, &
            0.0174438045_f64, &
            0.0203594338_f64, &
            0.0113349170_f64, &
            0.0046614185_f64, &
            0.0030346239_f64, &
            0.0012508731_f64, &
            0.0000782945_f64, &
            0.0235716330_f64, &
            0.0206304700_f64, &
            0.0204028340_f64, &
            0.0215105697_f64, &
            0.0183482070_f64, &
            0.0174161032_f64, &
            0.0155972434_f64, &
            0.0119269616_f64, &
            0.0147074804_f64, &
            0.0116182830_f64, &
            0.0087639138_f64, &
            0.0098563528_f64, &
            0.0096342355_f64, &
            0.0086477936_f64, &
            0.0083868302_f64, &
            0.0062576643_f64, &
            0.0077839825_f64, &
            0.0031415239_f64, &
            0.0006513246_f64, &
            0.0021137942_f64, &
            0.0004393452_f64, &
            0.0013662119_f64, &
            0.0003331251_f64, &
            0.0011613225_f64, &
            0.0004342867_f64, &
            0.0002031499_f64 /)
    else
       print *, ""
       print *,  'FEKETE_SUBRULE - Fatal error!'
       write(*, '(a,i8)') '  Illegal RULE = ', rule
       STOP
    end if

    !  The listed weights are twice what we want!.
    do s = 1, suborder_num
       suborder_w(s) = 0.5_f64 * suborder_w(s)
    end do
  end subroutine fekete_subrule


  !---------------------------------------------------------
  !> @brief forces an I4 to lie between given limits by wrapping.
  !> @details This code was first written by John Burkardt and is available
  !> online under the GNU LGPL license.
  !>    ILO = 4, IHI = 8
  !>
  !>   I  Value
  !>
  !>   -2     8
  !>   -1     4
  !>    0     5
  !>    1     6
  !>    2     7
  !>    3     8
  !>    4     4
  !>    5     5
  !>    6     6
  !>    7     7
  !>    8     8
  !>    9     4
  !>   10     5
  !>   11     6
  !>   12     7
  !>   13     8
  !>   14     4
  !> @param[IN] ival an integer value
  !> @param[IN] ilo  the desired lower bound
  !> @param[IN] ihi  the desired upper bound
  !> @param[OUT] res a wrapped version of ival
  function wrapping(ival, ilo, ihi) result(res)
    sll_int32, intent(in) ::  ihi
    sll_int32, intent(in) ::  ilo
    sll_int32, intent(in) ::  ival
    sll_int32 ::  jhi
    sll_int32 ::  jlo
    sll_int32 ::  res
    sll_int32 ::  wide

    jlo = min(ilo, ihi)
    jhi = max(ilo, ihi)

    wide = jhi - jlo + 1

    if(wide == 1) then
       res = jlo
    else
       res = jlo + MODULO( ival - jlo, wide)
    end if
  end function wrapping

  !---------------------------------------------------------------------------
  !> @brief maps T3 reference points to physical points.
  !> @details This code was first written by John Burkardt and is available
  !> online under the GNU LGPL license.
  !>    Given the vertices of an order 3 physical triangle and a point
  !>    (XSI,ETA) in the reference triangle, the routine computes the value
  !>    of the corresponding image point (X,Y) in physical space.
  !>
  !>    This routine is also appropriate for an order 4 triangle,
  !>    as long as the fourth node is the centroid of the triangle.
  !>
  !>    This routine may also be appropriate for an order 6
  !>    triangle, if the mapping between reference and physical space
  !>    is linear.  This implies, in particular, that the sides of the
  !>    image triangle are straight and that the "midside" nodes in the
  !>    physical triangle are literally halfway along the sides of
  !>    the physical triangle.
  !>
  !>  Reference Element T3:
  !>
  !>    |
  !>    1  3
  !>    |  |\
  !>    |  | \
  !>    S  |  \
  !>    |  |   \
  !>    |  |    \
  !>    0  1-----2
  !>    |
  !>    +--0--R--1-->
  !>
  !> @param[IN] node_xy(2,3) integer the coordinates of the vertices.
  !> The vertices are assumed to be the images of (0,0), (1,0) and (0,1).
  !> @param[IN] n integer the number of objects to transform
  !> @param[IN] ref(2,n) points in the reference triangle
  !> @param[OUT] phy(2,n) corresponding points in the physical triangle
  subroutine reference_to_physical_t3(node_xy, n, ref, phy)
    sll_int32,  intent(in) ::  n
    sll_real64, intent(in) :: node_xy(2,3)
    sll_real64, intent(in) :: ref(2,n)
    sll_real64, dimension(2,n), intent(out) :: phy
    sll_int32 ::  i

    do i = 1, 2
       phy(i,1:n) = node_xy(i,1) *(1.0_f64 - ref(1,1:n) - ref(2,1:n)) &
            + node_xy(i,2) * ref(1,1:n) &
            + node_xy(i,3) * ref(2,1:n)
    end do
  end subroutine reference_to_physical_t3


  !---------------------------------------------------------
  !> @brief Gives the fekete points coordinates and associated weights
  !> for a certain rule in a given triangle
  !> @details This code was first written by John Burkardt and is available
  !> online under the GNU LGPL license.
  !> @param[IN]  node_xy2 array of dimesion (2,3) containg the coordinates
  !>             of the edges of the triangle
  !> @param[IN] rule integer quadrature rule
  !> @return     xyw array of dimesion (3,n) containg the fekete points
  !>             and weights using the rule number given in parameter.
  !>             xyw(1,:) contains the x coordinates
  !>             of the fekete points, xyw(2, :) the y coordinates and
  !>             xyw(3, :) contains the associated weights.
  function fekete_points_and_weights(node_xy2, rule) result(xyw)
    sll_real64, dimension(2,3), intent(in) :: node_xy2
    sll_int32,  intent(in) :: rule
    !sll_int32  :: order
    sll_int32  :: order_num
    sll_int32  :: ierr
    sll_real64, allocatable :: xy(:,:)
    sll_real64, allocatable :: xy2(:,:)
    sll_real64, allocatable :: w(:)
    sll_real64, allocatable :: xyw(:,:)

    call fekete_order_num(rule, order_num)

    SLL_ALLOCATE(xy(1:2,1:order_num), ierr)
    SLL_ALLOCATE(xy2(1:2,1:order_num), ierr)
    SLL_ALLOCATE(w(1:order_num), ierr)
    SLL_ALLOCATE(xyw(1:3,1:order_num), ierr)

    call fekete_rule(rule, order_num, xy, w)
    call rearrange_fekete_rule(order_num, xy, w)
    xyw(3, :) = w

    call reference_to_physical_t3 (node_xy2, order_num, xy, xy2)
    xyw(1:2, :) = xy2

    SLL_DEALLOCATE_ARRAY(xy,  ierr)
    SLL_DEALLOCATE_ARRAY(xy2, ierr)
    SLL_DEALLOCATE_ARRAY(w,   ierr)

  end function fekete_points_and_weights


  !---------------------------------------------------------------------------
  !> @brief Fekete quadrature rule over a triangle
  !> @details To integrate the function \f$ f(x) \f$
  !> (real-valued and over a triangle) we use the Fekete formula
  !> \f[ \int_{\Omega} f(x)dx \approx \sum_{k=1}^{N} w_k f(x_k) \f]
  !> The only quadrature rule possible for now is 1 (10 points)
  !> @param[in]  f the function to be integrated
  !> @param[in]  pxy array of dimesion (2,3) containg the coordinates
  !>             of the edges of the triangle
  !> @return The value of the integral
  function fekete_integral( f, pxy)
    sll_real64, dimension(2, 3), intent(in) :: pxy
    sll_real64                  :: fekete_integral
    procedure(function_2D)      :: f
    sll_real64, dimension(3,10) :: xyw
    sll_real64, dimension(2)    :: v1
    sll_real64, dimension(2)    :: v2
    sll_real64, dimension(2)    :: v3
    sll_real64 :: a
    sll_real64 :: b
    sll_real64 :: c
    sll_real64 :: p
    sll_real64 :: area
    sll_int32 :: k
    sll_int32 :: N

    ! Initialiting the fekete points and weigths ......
    xyw(:,:) = 0._f64
    xyw = fekete_points_and_weights(pxy, 1)

    N = 10
    fekete_integral = 0._f64
    do k=1,N
       fekete_integral = fekete_integral + f(xyw(1,k), xyw(2,k))*xyw(3,k)
    end do

    ! Computing the area of the triangle
    ! v1 = Vector(p1, p2)
    v1(1) = pxy(1, 2) - pxy(1, 1)
    v1(2) = pxy(2, 2) - pxy(2, 1)
    a = sqrt(v1(1)*v1(1) + v1(2)*v1(2))
    ! v2 = Vector(p1, p3)
    v2(1) = pxy(1, 3) - pxy(1, 1)
    v2(2) = pxy(2, 3) - pxy(2, 1)
    b = sqrt(v2(1)*v2(1) + v2(2)*v2(2))
    ! v3 = Vector(p2, p3)
    v3(1) = pxy(1, 3) - pxy(1, 2)
    v3(2) = pxy(2, 3) - pxy(2, 2)
    c = sqrt(v3(1)*v3(1) + v3(2)*v3(2))
    ! Computing demi-perimeter
    p = 0.5_f64*(a+b+c)
    ! area
    area = sqrt(p*(p-a)*(p-b)*(p-c))

    fekete_integral = fekete_integral * area
  end function fekete_integral

  subroutine triangle_area ( node_xy, area )
    implicit none

    sll_real64 :: area
    sll_real64 :: node_xy(2,3)

    area = 0.5_f64 * ( &
         node_xy(1,1) * ( node_xy(2,2) - node_xy(2,3) ) &
         + node_xy(1,2) * ( node_xy(2,3) - node_xy(2,1) ) &
         + node_xy(1,3) * ( node_xy(2,1) - node_xy(2,2) ) )

    return
  end subroutine triangle_area

  !---------------------------------------------------------------------------
  !> @brief Writes fekete points coordinates of a hex-mesh reference triangle
  !> @details Takes the reference triangle of a hexmesh and computes the
  !> fekete points on it. Then it writes the results in a file following
  !> CAID/Django nomenclature.
  !> Output file : boxsplines_quadrature.txt
  !> @param[in]  rule integer for the fekete quadrature rule
  subroutine write_quadrature(rule)
    sll_int32, intent(in)       :: rule
    sll_int32                   :: out_unit
    character(len=25), parameter :: name = "boxsplines_quadrature.txt"
    sll_real64, dimension(2, 3) :: ref_pts
    sll_real64, dimension(:,:), allocatable :: quad_pw
    sll_int32  :: num_fek
    sll_int32  :: i
    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: w
    sll_real64 :: volume
    sll_int32  :: ierr
    ! Definition of reference triangle, such that:
    !    |
    !    1  3
    !    |  |  \
    !    |  |   \
    !    |  |    \
    !    |  |     \
    !    |  | _____\
    !    0  1      2
    !    |
    !    +--0-----1-->
    ref_pts(:,1) = (/ 0._f64, 0.0_f64 /)
    ! ref_pts(:,2) = (/ 1._f64, 0.0_f64 /)
    ref_pts(:,2) = (/ sqrt(3._f64)*0.5_f64, 0.5_f64 /)
    ref_pts(:,3) = (/ 0._f64, 1.0_f64 /)

    call triangle_area(ref_pts, volume)
    print *, "area triangle = ", volume

    ! Computing fekete points on that triangle
    call fekete_order_num(rule, num_fek)
    SLL_ALLOCATE(quad_pw(1:3, 1:num_fek), ierr)
    quad_pw = fekete_points_and_weights(ref_pts, rule)
    ! For Gaussian quadrature rule:
    ! num_fek = rule + 1
    ! SLL_ALLOCATE(quad_pw(1:3, 1:num_fek), ierr)
    ! quad_pw = gauss_triangle_points_and_weights(ref_pts, rule)

    open( file=name, status="replace", form="formatted", newunit=out_unit )

    write(out_unit, "(i6)") num_fek

    do i=1,num_fek
       x = quad_pw(1,i)
       y = quad_pw(2,i)
       w = quad_pw(3,i) * volume
       write(out_unit, "(2(g25.17,a,1x),(g25.17))") x, ",", y, ",", w
    end do
    close(out_unit)
  end subroutine write_quadrature

  !---------------------------------------------------------------------------
  !> @brief Writes on a file values of boxsplines on fekete points
  !> @details Following CAID structure, we write a file with the values
  !> of the basis function (box splines) on a reference element (triangle)
  !> fekete points. Output for DJANGO.
  !> Output file : boxsplines_basis_values.txt
  !> @param[in] deg integer with degree of splines
  subroutine write_basis_values(deg, rule)
    sll_int32,  intent(in)      :: deg
    sll_int32,  intent(in)      :: rule
    sll_real64, dimension(2, 3) :: ref_pts
    sll_real64, dimension(:, :), allocatable :: quad_pw
    sll_real64, dimension(:, :), allocatable :: disp_vec
    sll_int32,  parameter       :: out_unit=20
    character(len=*), parameter :: name = "boxsplines_basis_values.txt"
    sll_real64  :: x
    sll_real64  :: y
    sll_real64  :: val
    sll_int32   :: ierr
    sll_int32   :: nonZero
    sll_int32   :: ind_nZ
    sll_int32   :: nderiv
    sll_int32   :: idx, idy
    sll_int32   :: num_fek
    sll_int32   :: ind_fek
    ! Definition of reference triangle, such that:
    !    |
    !    1  3
    !    |  |  \
    !    |  |   \
    !    |  |    \
    !    |  |     \
    !    |  | _____\
    !    0  1      2
    !    |
    !    +--0-----1-->
    ref_pts(:,1) = (/ 0._f64, 0.0_f64 /)
    !    ref_pts(:,2) = (/ 1._f64, 0.0_f64 /)
    ref_pts(:,2) = (/ sqrt(3._f64)/2._f64, 0.5_f64 /)
    ref_pts(:,3) = (/ 0._f64, 1.0_f64 /)

    ! Computing fekete points on the reference triangle
    call fekete_order_num ( rule, num_fek )
    SLL_ALLOCATE(quad_pw(1:3, 1:num_fek), ierr)
    quad_pw = fekete_points_and_weights(ref_pts, rule)
    ! ! For Gaussian qudrature:
    ! num_fek = rule + 1
    ! SLL_ALLOCATE(quad_pw(1:3, 1:num_fek), ierr)
    ! quad_pw = gauss_triangle_points_and_weights(ref_pts, rule)

    nonZero = 3*deg*deg !> Number of non null box splines on a cell
    nderiv  = 1 !> Number of derivatives to be computed

    !> The displament vector correspond to the translation
    !> done to obtain the other non null basis functions
    SLL_ALLOCATE(disp_vec(2, nonZero), ierr)
    disp_vec(:,1) = 0._f64
    disp_vec(:,2) = ref_pts(:,1) - ref_pts(:,2)
    disp_vec(:,3) = ref_pts(:,1) - ref_pts(:,3)

    open (unit=out_unit,file=name,action="write",status="replace")

    write(out_unit, "(i6)") deg
    write(out_unit, "(i6)") nderiv

    do ind_nZ = 1, nonZero
       do ind_fek = 1, num_fek
          x = quad_pw(1, ind_fek) + disp_vec(1, ind_nZ)
          y = quad_pw(2, ind_fek) + disp_vec(2, ind_nZ)
          do idx = 0, nderiv
             do idy = 0, nderiv-idx
                val = boxspline_val_der(x, y, deg, idx, idy)
                write(out_unit, "(1(g25.18))", advance='no') val
                if ((idx<nderiv).or.(idy<nderiv-idx))  then
                   write(out_unit, "(1(a,1x))", advance='no') ","
                end if
             end do
          end do
          write(out_unit, *) ""
       end do
    end do

    close(out_unit)
    SLL_DEALLOCATE_ARRAY(disp_vec, ierr)
    SLL_DEALLOCATE_ARRAY(quad_pw, ierr)

  end subroutine write_basis_values

  !> @brief This function is supposed to write all django input files
  !> needed for a Django/Jorek simulation.
  !> @param[in] num_cells integer number of cells in a radius of the hexagonal
  !> mesh
  !> @param[in] deg integer degree of the splines that will be used for the
  !> interpolation
  subroutine write_all_django_files(num_cells, deg, rule, transf)
    sll_int32, intent(in)          :: num_cells
    sll_int32, intent(in)          :: deg
    sll_int32, intent(in)          :: rule
    type(sll_hex_mesh_2d), pointer :: mesh
    character(len=*),  intent(in)  :: transf

    mesh => new_hex_mesh_2d(num_cells, 0._f64, 0._f64, radius = 1._f64)

    call write_caid_files(mesh, transf, deg)
    call write_connectivity(mesh, deg)
    call write_basis_values(deg, rule)
    call write_quadrature(rule)
#ifdef DEBUG
    print*, 'write_all_django_files rule=', rule
#endif

  end subroutine write_all_django_files

end module sll_m_fekete_integration
