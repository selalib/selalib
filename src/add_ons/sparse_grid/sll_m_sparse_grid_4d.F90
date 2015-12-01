!> @ingroup sparse_grid
!> @author Katharina Kormann, IPP 
!> @brief Implementation of a 4D sparse grid with interpolation routines.
!> @details <DETAILED_DESCRIPTION>

module sll_m_sparse_grid_4d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use sll_m_periodic_interpolator_1d
use sll_m_arbitrary_degree_splines
use sll_m_lagrange_interpolator_1d
use sll_m_sparse_grid_interpolator
use sll_m_constants, only: sll_pi

implicit none
private

!> Sparse grid object for 4d with interpolation routines. Note in 4d we have only an implementation of a standard sparse grid with periodic boundary conditions, i.e. modified=0, boundary=0 compared to 2d and 3d).
type, public, extends(sparse_grid_interpolator) :: sparse_grid_interpolator_4d
sll_int32, dimension(:,:,:,:), pointer  :: index !< 4d mapping: for each 4d index l on the sparse grid, \a index gives the index of the first node belonging to this level  

contains
  procedure :: initialize => initialize_sg4d! Initialization routine
  procedure :: interpolate_from_interpolant_value ! Compute the value of the sparse grid interpolant at position eta
  procedure :: interpolate_disp_nconst_in_1d ! Interpolate along one (x)-direction with displacement non-constant in one (v)-direction
  procedure :: interpolate_disp_linnconst_in_1d => interpolate4d_disp_linnconst_in_1d
  procedure :: interpolate_disp_nconst_in_2d ! Interpolate along one (v)-direction with displacement non-constant in all x-directions
  procedure :: interpolate_const_disp

end type sparse_grid_interpolator_4d

contains




!------------------------------------------------------------------------------!
!!!! Interpolation routines !!!!


!> Compute the value of the sparse grid interpolant at position eta
  function interpolate_from_interpolant_value( interpolator,data, eta ) result(val)
    class(sparse_grid_interpolator_4d), intent(inout) :: interpolator !< Sparse grid object
    sll_real64 :: val !< Interpolated value at eta
    sll_real64, dimension(:), intent(in) :: data !< Values of the hierarchical surplus 
    sll_real64, dimension(:), intent(in) :: eta !< Position where to interpolate
    val =  interpolate_from_hierarchical_surplus(&
         interpolator,data,eta)

  end function interpolate_from_interpolant_value

!> Interpolation function for interpolation at (constantly) displaced grid points; displacement only in dimension dim. It is another implementation of the base-class function "interpolate_disp". The advantage is that we can not revisit nodes as we do in the recursive dimension-independently-programmed version.
  subroutine interpolate_const_disp(interpolator,dorder,displacement,data_in, data_out,hiera)
    class(sparse_grid_interpolator_4d), intent(inout) :: interpolator !< Sparse grid object
    sll_real64, dimension(:), intent(inout) :: data_in !< Values of the hierarchical surplus on input.
    sll_real64, dimension(:), intent(out) :: data_out !< Value of the function or the hierarchical surplus (depending on value of \a hiera) for the displaced data points. 
    sll_int32, dimension(:), intent(in) :: dorder !< dorder(1) gives the dimension along which we have the displacement; dorder(2:4) give the remaining dimensions
    sll_real64, intent(in) ::displacement !< Constant displacement along dimension dorder(1)
    logical, intent(in) :: hiera !< If the result should be the hierarchical surplus, define \a hiera = .TRUE.; if the result should be the function values at the data points give \a hiera = .FALSE.

    sll_int32 :: i1,i2,i3,i4,k2,k3,k4,counter,j
    sll_int32, dimension(4) :: l, no,ind_order
    sll_int32, dimension(:,:), allocatable :: ind

  SLL_ALLOCATE(ind(interpolator%max_level+1,4), i1);

  ind_order(dorder(1)) = 0
  no(dorder(1)) = 1
  do i3 = 0,interpolator%levels(dorder(3))
     ind_order(dorder(3)) = i3
     l(dorder(3)) = i3
     no(dorder(3)) = max(2**(i3-1),1);
     do i4 = 0,min(interpolator%max_level - i3,interpolator%levels(dorder(4)))
        ind_order(dorder(4)) = i4
        l(dorder(4)) = i4
        no(dorder(4)) = max(2**(i4-1),1);
        do i2 = 0,min(interpolator%max_level - i3 -i4,interpolator%levels(dorder(2)))
           ind_order(dorder(2)) = i2
           no(dorder(2)) = max(2**(i2-1),1);
           ind(1,dorder(1)) = 0;
           do k2 = 0,no(dorder(2))-1
              ind(ind_order(dorder(2))+1,dorder(2)) = k2;
              do k3 = 0,no(dorder(3))-1
                 ind(ind_order(dorder(3))+1,dorder(3)) = k3;
                 do k4 = 0,no(dorder(4))-1
                    ind(ind_order(dorder(4))+1,dorder(4)) = k4;
                    counter = interpolator%index(&
                         ind_order(1),ind_order(2),ind_order(3),ind_order(4))+&
                         ind(ind_order(1)+1,1)*no(2)*no(3)*no(4)&
                         +ind(ind_order(2)+1,2)*no(3)*no(4)+&
                         ind(ind_order(3)+1,3)*no(4)+ind(ind_order(4)+1,4)
                    ! Evaluate along dorder(1)-stripe
                    call interpolator%interpolate_disp_1d_periodic&
                         (displacement,dorder(1),&
                         min(interpolator%levels(dorder(1)),&
                         interpolator%max_level-ind_order(dorder(2))-&
                         ind_order(dorder(3))-&
                         ind_order(dorder(4))),counter,data_in,data_out,hiera)
                 end do
              end do
           end do
        end do
     end do
  end do

  if (hiera .EQV. .FALSE.) then
     ! Dehierarchization along dimension dorder(1) only
     do j=interpolator%order,2,-1
        call interpolator%dehierarchical_part_order&
             (data_out,&
             interpolator%dim,2,dorder,j)
     end do

     call interpolator%dehierarchical_part(data_out,&
          interpolator%dim,2,dorder)
  end if

end subroutine Interpolate_const_disp


!> Functionality: Interpolates the function values for a displacement on in dimension (periodic b.c. i.e. dimension 1 or 2) where the displacement is allowed to be non-constant in one other dimension (Dirichlet b.c. i.e. dimension 3 or 3).
subroutine interpolate_disp_nconst_in_1d(interpolator,displacement,dorder,data_in, data_out)
  class(sparse_grid_interpolator_4d), intent(inout) :: interpolator !< sparse grid object
  sll_real64, dimension(:), intent(inout) :: data_in !< hierarchical surplus of the present function
  sll_real64, dimension(:), intent(out) :: data_out !< value of the displaced function
  sll_int32, dimension(:), intent(in) :: dorder !< dorder: Ordering of the dimensions. dorder(1) (=1 or 2) gives the dimension where we want to displace, dorder(2) (=3 or 4) gives the dimension of which the displacement is dependent. dorder(3) = 1 or 2 not dorder(1) and dorder(4) = 3 or 4 not dorder(2).
  sll_real64,dimension(:), intent(in) ::displacement !< Vector containing the values of the displacement (in hierarchical order, one dimensional)
  sll_int32 :: i1,i2,i3,i4,k2,k3,k4,counter,j, index_parent,index_parent_old
  sll_int32, dimension(4) :: l, no,ind_order
  sll_int32, dimension(:,:), allocatable :: ind
  sll_real64 :: coordinate_self, coordinate_ancestor, factor,disp

  SLL_ALLOCATE(ind(interpolator%max_level+1,4), i1);

  ! Dehierarchization along dimension dorder(1) only
  do j=interpolator%order,2,-1
     call interpolator%sparse_grid_interpolator%dehierarchical_part_order&
          (data_in,1,1,dorder,j)
  end do

  call interpolator%sparse_grid_interpolator%dehierarchical_part(data_in,1,1,dorder)
  ! Interpolation in dorder(1)/dorder(2)-plane
  ind_order(dorder(1)) = 0
  no(dorder(1)) = 1
  do i3 = 0,interpolator%levels(dorder(3))
     ind_order(dorder(3)) = i3
     l(dorder(3)) = i3
     no(dorder(3)) = max(2**(i3-1),1);
     do i4 = 0,min(interpolator%max_level - i3,interpolator%levels(dorder(4)))
        ind_order(dorder(4)) = i4
        l(dorder(4)) = i4
        no(dorder(4)) = max(2**(i4-1),1);
        do i2 = 0, min(interpolator%max_level - i3 -i4,interpolator%levels(dorder(2)))
           ind_order(dorder(2)) = i2
           no(dorder(2)) = max(2**(i2-1),1);
           ind(1,dorder(1)) = 0;
           do k2 = 0,no(dorder(2))-1
              ind(ind_order(dorder(2))+1,dorder(2)) = k2;
              disp = displacement(2**(i2-1)+k2+1)
              do k3 = 0,no(dorder(3))-1
                 ind(ind_order(dorder(3))+1,dorder(3)) = k3;
                 do k4 = 0,no(dorder(4))-1
                    ind(ind_order(dorder(4))+1,dorder(4)) = k4;
                    counter = interpolator%index(&
                         ind_order(1),ind_order(2),ind_order(3),ind_order(4))+&
                         ind(ind_order(1)+1,1)*no(2)*no(3)*no(4)&
                         +ind(ind_order(2)+1,2)*no(3)*no(4)+&
                         ind(ind_order(3)+1,3)*no(4)+ind(ind_order(4)+1,4)
                    ! Evaluate along dorder(1)-stripe
                    call interpolator%sparse_grid_interpolator%interpolate_disp_1d_periodic_self&
                         (disp,dorder(1),&
                         min(interpolator%levels(dorder(1)),&
                         interpolator%max_level-ind_order(dorder(2))-&
                         ind_order(dorder(3))-&
                         ind_order(dorder(4))),counter,data_in,data_out)
                    ! Evaluate hierarchically along dorder(2) dimension (dorder(1)-stripe-wise)
                    index_parent = max(interpolator%hierarchy(counter)%parent(2*dorder(2)-1),&
                         interpolator%hierarchy(counter)%parent(2*dorder(2)))
                    coordinate_self = interpolator%hierarchy(counter)%coordinate(dorder(2))
                    index_parent_old = counter;
                    do while(index_parent<index_parent_old)
                       coordinate_ancestor = interpolator%hierarchy(index_parent)%&
                            coordinate(dorder(2))
                       call interpolator%basis_function((coordinate_self-coordinate_ancestor)/&
                            interpolator%length(dorder(2))*&
                            2**(interpolator%hierarchy(index_parent)%level(dorder(2))), factor,&
                            interpolator%hierarchy(index_parent)%function_type(dorder(2)))
                       call interpolator%sparse_grid_interpolator%interpolate_disp_1d_periodic_for_neighbor&
                            (disp,factor,&
                            dorder(1),min(interpolator%levels(dorder(1)),&
                            interpolator%max_level-&
                            interpolator%hierarchy(index_parent)%level(dorder(2))-&
                            ind_order(dorder(3))-ind_order(dorder(4))),&
                            min(interpolator%levels(dorder(1)),&
                            interpolator%max_level-ind_order(dorder(2))-&
                            ind_order(dorder(3))-&
                            ind_order(dorder(4))),index_parent,counter,&
                            data_in,data_out)
                       index_parent_old = index_parent;
                       index_parent =  &
                            max(interpolator%hierarchy(index_parent)%parent(2*dorder(2)-1),&
                            interpolator%hierarchy(index_parent)%parent(2*dorder(2)))
                    end do

                 end do
              end do
           end do
        end do
     end do
  end do

  ! Dehierarchization along dimension dorder(3) dorder(4)
  do j=interpolator%order,2,-1
     call interpolator%sparse_grid_interpolator%dehierarchical_part_order&
          (data_out,4,3,dorder,j)
  end do

  call interpolator%sparse_grid_interpolator%dehierarchical_part(data_out,4,3,dorder)


end subroutine Interpolate_disp_nconst_in_1d

!> As \a interpolate_disp_nconst_in_1d but displacement dependent on displacement*coordinate(dorder(2))
subroutine interpolate4d_disp_linnconst_in_1d(interpolator,displacement,dorder,data_in, data_out)
  class(sparse_grid_interpolator_4d), intent(inout) :: interpolator !< Sparse grid object
  sll_real64, dimension(:), intent(inout) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_int32, dimension(:), intent(in) :: dorder
  sll_real64, intent(in) ::displacement
  sll_int32 :: i1,i2,i3,i4,k2,k3,k4,counter,j, index_parent,index_parent_old
  sll_int32, dimension(4) :: l, no,ind_order
  sll_int32, dimension(:,:), allocatable :: ind
  sll_real64 :: coordinate_self, coordinate_ancestor, factor,disp

  SLL_ALLOCATE(ind(interpolator%max_level+1,4), i1);

  ! Dehierarchization along dimension dorder(1) only
  do j=interpolator%order,2,-1
     call interpolator%sparse_grid_interpolator%dehierarchical_part_order&
          (data_in,1,1,dorder,j)
  end do

  call interpolator%sparse_grid_interpolator%dehierarchical_part(data_in,1,1,dorder)

 ! Interpolation in dorder(1)/dorder(2)-plane
  ind_order(dorder(1)) = 0
  no(dorder(1)) = 1
  do i3 = 0,interpolator%levels(dorder(3))
     ind_order(dorder(3)) = i3
     l(dorder(3)) = i3
     no(dorder(3)) = max(2**(i3-1),1);
     do i4 = 0,min(interpolator%max_level - i3,interpolator%levels(dorder(4)))
        ind_order(dorder(4)) = i4
        l(dorder(4)) = i4
        no(dorder(4)) = max(2**(i4-1),1);
        do i2 = 0, min(interpolator%max_level - i3 -i4,interpolator%levels(dorder(2)))
           ind_order(dorder(2)) = i2
           no(dorder(2)) = max(2**(i2-1),1);
           ind(1,dorder(1)) = 0;
           do k2 = 0,no(dorder(2))-1
              ind(ind_order(dorder(2))+1,dorder(2)) = k2;
              do k3 = 0,no(dorder(3))-1
                 ind(ind_order(dorder(3))+1,dorder(3)) = k3;
                 do k4 = 0,no(dorder(4))-1
                    ind(ind_order(dorder(4))+1,dorder(4)) = k4;
                    counter = interpolator%index(&
                         ind_order(1),ind_order(2),ind_order(3),ind_order(4))+&
                         ind(ind_order(1)+1,1)*no(2)*no(3)*no(4)&
                         +ind(ind_order(2)+1,2)*no(3)*no(4)+&
                         ind(ind_order(3)+1,3)*no(4)+ind(ind_order(4)+1,4)
                    disp = displacement*interpolator%hierarchy(counter)%coordinate(dorder(2))
                    ! Evaluate along dorder(1)-stripe
                    call interpolator%sparse_grid_interpolator%interpolate_disp_1d_periodic_self&
                    (disp,dorder(1),&
                    min(interpolator%levels(dorder(1)),interpolator%max_level&
                    -ind_order(dorder(2))-ind_order(dorder(3))-&
                    ind_order(dorder(4))),counter,data_in,data_out)
                    ! Evaluate hierarchically along dorder(2) dimension (dorder(1)-stripe-wise)
                    index_parent = max(interpolator%hierarchy(counter)%parent(2*dorder(2)-1),&
                         interpolator%hierarchy(counter)%parent(2*dorder(2)))
                    coordinate_self = interpolator%hierarchy(counter)%coordinate(dorder(2))
                    index_parent_old = counter;
                    do while(index_parent<index_parent_old)
                       coordinate_ancestor = interpolator%hierarchy(index_parent)%&
                            coordinate(dorder(2))
                       call interpolator%basis_function((coordinate_self-coordinate_ancestor)/&
                            interpolator%length(dorder(2))*&
                            2**(interpolator%hierarchy(index_parent)%level(dorder(2))), factor,&
                            interpolator%hierarchy(index_parent)%function_type(dorder(2)))
                       call interpolator%sparse_grid_interpolator%interpolate_disp_1d_periodic_for_neighbor&
                            (disp,factor,&
                            dorder(1),min(interpolator%levels(dorder(1)),&
                            interpolator%max_level-&
                            interpolator%hierarchy(index_parent)%level(dorder(2))-&
                            ind_order(dorder(3))-ind_order(dorder(4))),&
                            min(interpolator%levels(dorder(1)),&
                            interpolator%max_level-&
                            ind_order(dorder(2))-ind_order(dorder(3))-&
                            ind_order(dorder(4))),index_parent,counter,&
                            data_in,data_out)
                       index_parent_old = index_parent;
                       index_parent =  &
                            max(interpolator%hierarchy(index_parent)%parent(2*dorder(2)-1),&
                            interpolator%hierarchy(index_parent)%parent(2*dorder(2)))
                    end do

                 end do
              end do
           end do
        end do
     end do
  end do

  ! Dehierarchization along dimension dorder(3) dorder(4)
  do j=interpolator%order,2,-1
     call interpolator%sparse_grid_interpolator%dehierarchical_part_order&
          (data_out,4,3,dorder,j)
  end do

  call interpolator%sparse_grid_interpolator%dehierarchical_part(data_out,4,3,dorder)

end subroutine Interpolate4d_disp_linnconst_in_1d



!> As previous function but with displacement displacement*coordinate(dorder(2))
subroutine interpolate_disp_nconst_in_2d(interpolator,displacement,dorder,data_in, data_out)
  class(sparse_grid_interpolator_4d), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_int32, dimension(:), intent(in) :: dorder
  sll_real64, dimension(:), intent(in) ::displacement
  sll_int32 :: i1,i2,i3,i4,k2,k3,k4,counter,j, index_parent, index_old, index_upper,index_upper_old
  sll_int32 :: counter_disp,l23
  sll_int32, dimension(4) :: l, no,ind_order
  sll_int32, dimension(:,:), allocatable :: ind
  sll_real64 :: coordinate_self, coordinate_ancestor, coordinate_self_upper,factor_upper, factor_lower, factor,disp

  SLL_ALLOCATE(ind(interpolator%max_level+1,4), i1);



  ! Dehierarchization along dimension dorder(1) only
  do j=interpolator%order,2,-1
     call interpolator%sparse_grid_interpolator%dehierarchical_part_order&
          (data_in,1,1,dorder,j)
  end do

  call interpolator%sparse_grid_interpolator%dehierarchical_part(&
       data_in,1,1,dorder)

  !print*, data_in

  ! Interpolation in dorder(1)/dorder(2)-plane
  ind_order(dorder(1)) = 0
  no(dorder(1)) = 1
  do i4 = 0,interpolator%levels(dorder(4))
     ind_order(dorder(4)) = i4
     l(dorder(4)) = i4
     no(dorder(4)) = max(2**(i4-1),1);
     counter_disp = 0;
     do l23=0,min(interpolator%max_level-i4, interpolator%levels(dorder(2)))
        do i2 = 0, l23
           ind_order(dorder(2)) = i2
           no(dorder(2)) = max(2**(i2-1),1);
           i3 = l23-i2
           ind_order(dorder(3)) = i3
           l(dorder(3)) = i3
           no(dorder(3)) = max(2**(i3-1),1);
           ind(1,dorder(1)) = 0;
           do k2 = 0,max(2**(ind_order(dorder(2))-1),1)-1
              ind(ind_order(dorder(2))+1,dorder(2)) = k2;
              do k3 = 0,max(2**(ind_order(dorder(3))-1),1)-1
                 counter_disp = counter_disp+1;

                 disp = displacement(counter_disp);
                 ind(ind_order(dorder(3))+1,dorder(3)) = k3;
                 do k4 = 0,max(2**(ind_order(dorder(4))-1),1)-1
                    ind(ind_order(dorder(4))+1,dorder(4)) = k4;
              !      print*, ind_order
              !      print*, k2,k3,k4
              !      print*, ind(ind_order(1)+1,1),ind(ind_order(2)+1,2),ind(ind_order(3)+1,3),ind(ind_order(4)+1,4)
                   ! print*, '------------------------------------------------- '
                    counter = interpolator%index(&
                         ind_order(1),ind_order(2),ind_order(3),ind_order(4))+&
                         ind(ind_order(1)+1,1)*no(2)*no(3)*no(4)&
                         +ind(ind_order(2)+1,2)*no(3)*no(4)+&
                         ind(ind_order(3)+1,3)*no(4)+ind(ind_order(4)+1,4)

                    ! Evaluate along dorder(1)-stripe
                    call interpolator%sparse_grid_interpolator%interpolate_disp_1d_periodic_self&
                         (disp,dorder(1),&
                         min(interpolator%levels(dorder(1)),&
                         interpolator%max_level-ind_order(dorder(2))-&
                         ind_order(dorder(3))-&
                         ind_order(dorder(4))),counter,data_in,data_out)
                    ! Evaluate hierarchically along dorder(2)&dorder(3) dimension (dorder(1)-stripe-wise)
                    index_upper = counter
                    index_upper_old = index_upper+1

                    !print*, counter, counter_disp , disp
                    !print*, counter , interpolator%index(&
                    !     ind_order(1),ind_order(2),ind_order(3),ind_order(4)), ind_order

                    coordinate_self = interpolator%hierarchy(counter)%coordinate(dorder(2))
                    coordinate_self_upper = interpolator%hierarchy(counter)%coordinate(dorder(3))
                    !print*, counter
                    do while (index_upper<index_upper_old)
                       coordinate_ancestor = interpolator%hierarchy(index_upper)%&
                               coordinate(dorder(3))
                       call interpolator%basis_function((coordinate_self_upper-coordinate_ancestor)/&
                            interpolator%length(dorder(3))*&
                            2**(interpolator%hierarchy(index_upper)%level(dorder(3))), &
                            factor_upper,&
                            interpolator%hierarchy(index_upper)%function_type(dorder(3)))

                       if(index_upper == counter) then
                          index_parent = &
                               max(interpolator%hierarchy(index_upper)%parent(2*dorder(2)-1),&
                               interpolator%hierarchy(index_upper)%parent(2*dorder(2)))
                          index_old = index_upper
                       else
                          index_parent = index_upper;
                          index_old = index_parent+1
                       end if

                       do while(index_parent<index_old)
                          coordinate_ancestor = interpolator%hierarchy(index_parent)%&
                               coordinate(dorder(2))
                          call interpolator%basis_function((coordinate_self-coordinate_ancestor)/&
                               interpolator%length(dorder(2))*&
                               2**(interpolator%hierarchy(index_parent)%level(dorder(2))), &
                               factor_lower,&
                               interpolator%hierarchy(index_parent)%function_type(dorder(2)))
                          factor = factor_upper*factor_lower
                          !if(counter==30) then
                          !   print*, index_parent,data_out(30)
                          !end if
                         ! print*, counter, index_parent, factor, factor_upper, factor_lower
                          call interpolator%sparse_grid_interpolator%interpolate_disp_1d_periodic_for_neighbor&
                               (disp,factor,&
                               dorder(1),min(interpolator%levels(dorder(1)),&
                               interpolator%max_level-&
                               interpolator%hierarchy(index_parent)%level(dorder(2))-&
                               interpolator%hierarchy(index_parent)%level(dorder(3))-&
                               ind_order(dorder(4))),&
                               min(interpolator%levels(dorder(1)),&
                               interpolator%max_level-ind_order(dorder(2))-&
                               ind_order(dorder(3))-&
                               ind_order(dorder(4))),index_parent,counter,&
                               data_in,data_out)
                          index_old = index_parent
                          index_parent =  &
                               max(interpolator%hierarchy(index_parent)%parent(2*dorder(2)-1),&
                               interpolator%hierarchy(index_parent)%parent(2*dorder(2)))
                       end do
                       index_upper_old = index_upper;
                       index_upper = &
                            max(interpolator%hierarchy(index_upper)%parent(2*dorder(3)-1),&
                            interpolator%hierarchy(index_upper)%parent(2*dorder(3)));
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

! Hierarchization along dimension dorder(1)-dorder(3)
call interpolator%sparse_grid_interpolator%hierarchical_part(data_out,3,1,dorder)

do j=2,interpolator%order
   call interpolator%sparse_grid_interpolator%hierarchical_part_order&
        (data_out,3,1,dorder,j)
end do


! Dehierarchization of the hierarchical surplus
do j=interpolator%order,2,-1
   call interpolator%sparse_grid_interpolator%dehierarchical_order&
        (data_out,j)
end do

call interpolator%sparse_grid_interpolator%dehierarchical&
     (data_out)


end subroutine Interpolate_disp_nconst_in_2d



!PN DEFINED BUT NOT USED
!! Note dorder should contain: dorder(1) dimension with displacement (1 or 2), dorder(2) dimension where displacement is non-constant (3 or 4), dorder(3) other of 1 or 2, dorder(4) other of 3 or 4
!subroutine displace_disp_nconst1(interpolator, surplus,  data, dorder,displacement)
!  class(sparse_grid_interpolator_4d), intent(inout) :: interpolator
!  sll_real64, dimension(:), intent(in) :: surplus
!  sll_real64, dimension(:), intent(out) :: data
!  sll_int32, dimension(:), intent(in) :: dorder
!  sll_real64, dimension(:), intent(in) :: displacement
!  sll_int32 :: i3, i4, i2, i1, k1,k2,k3,k4, counter
!  sll_int32, dimension(4) :: ind_order,l
!  sll_real64, dimension(2) :: eta
!  sll_int32, dimension(4) :: no
!  sll_int32, dimension(:,:), allocatable :: ind
!
!  SLL_ALLOCATE(ind(interpolator%max_level,4), i1);
!
!
!  ind_order(dorder(1)) = 0
!  do i3 = 0,interpolator%levels(dorder(3))
!     ind_order(dorder(3)) = i3
!     l(dorder(3)) = i3
!     no(dorder(3)) = max(2**(i3-1),1);
!     do i4 = 0,min(interpolator%max_level - i3,interpolator%levels(dorder(4)))
!        ind_order(dorder(4)) = i4
!        l(dorder(4)) = i4
!        no(dorder(4)) = max(2**(i4-1),1);
!        do i2 = 0,min(interpolator%max_level - i3 -i4,interpolator%levels(dorder(2)))
!           ind_order(dorder(2)) = i2
!           no(dorder(2)) = max(2**(i2-1),1);
!           do i1 = 0,min(interpolator%max_level - i3 - i4 -i2,interpolator%levels(dorder(1)))
!              ind_order(dorder(1)) = i1
!              no(dorder(1)) = max(2**(i1-1),1);
!              do k1 = 0,max(2**(ind_order(dorder(1))-1),1)-1
!                 ind(ind_order(dorder(1))+1,dorder(1)) = k1;
!                 do k2 = 0,max(2**(ind_order(dorder(2))-1),1)-1
!                     ind(ind_order(dorder(2))+1,dorder(2)) = k2;
!                    do k3 = 0,max(2**(ind_order(dorder(3))-1),1)-1
!                       ind(ind_order(dorder(3))+1,dorder(3)) = k3;
!                       do k4 = 0,max(2**(ind_order(dorder(4))-1),1)-1
!                          ind(ind_order(dorder(4))+1,dorder(4)) = k4;
!                          counter = interpolator%index(&
!                               ind_order(1),ind_order(2),ind_order(3),ind_order(4))+&
!                               ind(ind_order(1)+1,1)*no(2)*no(3)*no(4)&
!                               +ind(ind_order(2)+1,2)*no(3)*no(4)+&
!                               ind(ind_order(3)+1,3)*no(4)+ind(ind_order(4)+1,4)
!                          eta(1) = interpolator%hierarchy(counter)%coordinate(dorder(1))+&
!                          displacement(2**(i2-1)+k2)
!                          eta(2) = interpolator%hierarchy(counter)%coordinate(dorder(2))
!
!                          l(dorder(2)) = ind_order(dorder(2))
!                          data(counter) =  interpolate_from_2D_hierarchical_surplus( &
!                               interpolator,surplus, eta, dorder, no, ind, l );
!                       end do
!                    end do
!                 end do
!              end do
!           end do
!        end do
!     end do
!  end do
!
!  call interpolator%sparse_grid_interpolator%hierarchical_part(data,2,1,dorder)
!
!
!  call interpolator%sparse_grid_interpolator%dehierarchical(data)
!
!
!end subroutine displace_disp_nconst1


! helper functions

 function interpolate_from_hierarchical_surplus( interpolator,data, eta ) result(val)
   class(sparse_grid_interpolator_4d), intent(inout) :: interpolator
   sll_int32 :: j,l1,l2,l3,level
    sll_real64 :: val
    sll_real64,dimension(:), intent(in) :: data,eta
    sll_real64,dimension(4) :: eta_norm
    sll_real64,dimension(4) :: phi
    sll_int32, dimension(4) :: no,l
    sll_int32,dimension(:,:), allocatable :: ind
    sll_real64 :: scale
    sll_int32 :: index

    SLL_ALLOCATE(ind(0:interpolator%max_level,1:4),j)

    val = 0.0_f64
    ind(0:1,1:4) = 0

    do j=1,4
       eta_norm(j) = (eta(j)-interpolator%eta_min(j))/interpolator%length(j)
       eta_norm(j) = modulo(eta_norm(j),1.0_f64)

       scale = 0.5_f64
       do level = 2, interpolator%max_level
          ind(level,j) = ind(level-1,j)*2
          if (eta_norm(j)> scale*(ind(level,j)+1)) then
             ind(level,j) = ind(level,j)+1
          end if
          scale = scale*0.5_f64
       end do
    end do

    do level = 0, interpolator%max_level
       do l1 = 0, min(level,interpolator%levels(1))
          l(1) = l1
          no(1) = max(2**(l1-1),1)
          do l2=0,min(level-l1,interpolator%levels(2))
             l(2) = l2
             no(2) = max(2**(l2-1),1)
             do l3 =max(0,level-l1-l2-interpolator%levels(4)),min(level-l1-l2,interpolator%levels(3))
                no(3) = max(2**(l3-1),1)
                l(3) = l3
                l(4) = level-l1-l2-l3
                no(4) = max(2**(l(4)-1),1)

                index = interpolator%index(l1,l2,l3,&
                     l(4))+ind(l1,1)*no(2)*no(3)*no(4)&
                     +ind(l2,2)*no(3)*no(4)+ind(l3,3)*no(4)+ind(l(4),4)
                do j=1,4
                   call interpolator%basis_function(real(2**(max(l(j),1)),f64)*eta_norm(j)&
                        -real(2*ind(l(j),j),f64)-1.0_f64, phi(j),&
                        interpolator%hierarchy(index)%function_type(j))
                end do
                val = val + data(index)&
                     *phi(1)*phi(2)*phi(3)*phi(4)
             end do
          end do
       end do
    end do

  end function interpolate_from_hierarchical_surplus




!PN DEFINED BUT NOT USED
! function interpolate_from_2D_hierarchical_surplus( interpolator, surplus,eta, dorder,no_in,ind,l ) result(val)
!   class(sparse_grid_interpolator_4d), intent(inout) :: interpolator
!   sll_real64, dimension(:), intent(in) :: surplus
!    sll_int32 :: j,l1,l2,level
!    sll_real64 :: val
!    sll_real64,dimension(:), intent(in) :: eta
!    sll_int32, dimension(:), intent(in) ::dorder
!    sll_real64,dimension(4) :: eta_norm
!    sll_real64,dimension(4) :: phi
!    sll_int32, dimension(:), intent(inout) :: l
!    sll_int32, dimension(:), intent(in) :: no_in
!    sll_int32, dimension(4) :: no
!    sll_int32,dimension(:,:), intent(inout) :: ind
!    sll_real64 :: scale
!    sll_int32 :: index,maxl2
!
!    val = 0.0_f64
!    ind(1:2,1:4) = 0
!    no(dorder(3)) = no_in(dorder(3));
!    no(dorder(4)) = no_in(dorder(4));
!
!! Note this could be organized more efficiently by reusing data
!    do j=1,2
!       eta_norm(j) = (eta(j)-interpolator%eta_min(dorder(j)))/interpolator%length(dorder(j))
!       eta_norm(j) = modulo(eta_norm(j),1.0_f64)
!       scale = 0.5_f64
!       do level = 0, interpolator%max_level
!          ind(level+1,dorder(j)) = ind(level,dorder(j))*2
!
!          if (eta_norm(j)> scale*(ind(level+1,dorder(j))+1)) then
!             ind(level+1,dorder(j)) = ind(level+1,dorder(j))+1
!          end if
!          scale = scale*0.5_f64
!       end do
!    end do
!
!
!
!    maxl2 = l(dorder(2));
!    do l1 = 0, min(interpolator%max_level-l(dorder(3))-l(dorder(4)),interpolator%levels(dorder(1)))
!       l(dorder(1)) = l1
!       no(dorder(1)) = max(2**(l1-1),1)
!       do l2 = 0,min(maxl2,interpolator%max_level-l(dorder(3))-l(dorder(4))-l1)
!          l(dorder(2)) = l2
!          no(dorder(2)) = max(2**(l2-1),1)
!          index = interpolator%index(l(1),l(2),l(3),l(4))+ind(l(1)+1,1)*no(2)*no(3)*no(4)&
!               +ind(l(2)+1,2)*no(3)*no(4)+ind(l(3)+1,3)*no(4)+ind(l(4)+1,4)
!          do j=1,2
!             call interpolator%basis_function(real(2**(max(l(dorder(j)),1)),f64)*eta_norm(j)&
!                  -real(2*ind(l(dorder(j))+1,dorder(j)),f64)-1.0_f64, phi(j),&
!                  interpolator%hierarchy(index)%function_type(dorder(j)))
!          end do
!          val = val + surplus(index)&
!               *phi(1)*phi(2)
!          !print*, index, val, phi(1)
!
!       end do
!    end do
!
!
!  end function interpolate_from_2D_hierarchical_surplus



!!!! End interpolation routines !!!!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!!!! Initialization routines !!!!


  !> Initialization function. Set up the hierarchy of the sparse grid
  subroutine initialize_sg4d( &
    interpolator, &
    levels, &
    order, &
    interpolation, &
    interpolation_type, &
    eta_min, &
    eta_max)
    class(sparse_grid_interpolator_4d), intent(inout) :: interpolator !< sparse grid object 
    sll_real64,dimension(:), intent(in)           :: eta_min !< \a eta_min defines the lower bound of the domain

    sll_real64,dimension(:), intent(in)           :: eta_max !< \a eta_max defines the upper bound of the domain

    sll_int32, dimension(:), intent(in)           :: levels !> maximum level in the sparse grid along each direction
    sll_int32, intent(in)                         :: order !< \a order of the sparse grid functions (degree of polynomial)
    sll_int32, intent(in)                         :: interpolation !< order of the interpolator (usually order+1 to have the same order)
    sll_int32, intent(in)                         :: interpolation_type !< Choose spline (\a interpolation_type = 0) or Lagrange (\a interpolation_type = 1) interpolation for the 1D interpolators if not traditional sparse grid interpolation is used.
    sll_int32                                     :: i,j, k1,k2,k3,k4,l1,l2,l3,l4,l,counter
    sll_int32                                     :: ierr
    sll_int32, dimension(:) ,allocatable          :: novec,lvec,kvec

    interpolator%dim = 4;
    SLL_ALLOCATE(lvec(interpolator%dim),ierr);
    SLL_ALLOCATE(kvec(interpolator%dim),ierr);
    SLL_ALLOCATE(novec(interpolator%dim),ierr);

    interpolator%max_level = levels(1);
    do l=2,interpolator%dim
       interpolator%max_level = max(levels(l),interpolator%max_level);
    end do  
!    if (interpolator%modified == 1) then
!       interpolator%max_level = interpolator%max_level + 2;
!    end if
    
    interpolator%size_basis = 0;
    do l = 0, interpolator%max_level
       do l1 = 0 , min(l,levels(1))
          do l2 = 0, min(l-l1,levels(2))
             do l3 = max(0,l-l1-l2-levels(4)), min(l-l1-l2,levels(3))
                l4 = l-l1-l2-l3
                   interpolator%size_basis = &
                        interpolator%size_basis + &
                        max(2**(l1-1),1)*max(2**(l2-1),1)*max(2**(l3-1),1)*max(2**(l4-1),1)
             end do
          end do
       end do
    end do

    call  interpolator%initialize_sg( levels, order, interpolation,&
    interpolation_type, eta_min, eta_max);

    SLL_ALLOCATE(interpolator%index(0:interpolator%levels(1),0:interpolator%levels(2),0:interpolator%levels(3),0:interpolator%levels(4)),ierr)

    ! Set the hierarchy of the grid
    counter = 1
    do l = 0, interpolator%max_level
       interpolator%level_mapping(l) = counter;
       do l1 = 0 , min(l,interpolator%levels(1))
          novec(1) = max(2**(l1-1),1)
          lvec(1) = l1;
          do l2 = 0, min(l-l1,interpolator%levels(2))
             novec(2) = max(2**(l2-1),1)
             lvec(2) = l2;
             do l3 = max(0,l-l1-l2-interpolator%levels(4)), min(l-l1-l2,interpolator%levels(3))
                novec(3) = max(2**(l3-1),1)
                lvec(3) = l3
                lvec(4) = l-l1-l2-l3
                l4 = lvec(4)
                novec(4) = max(2**(l4-1),1)
                interpolator%index(l1,l2,l3,l4) = counter
                do k1=0,novec(1)-1
                   kvec(1) = k1
                   do k2=0,novec(2)-1
                      kvec(2) = k2
                      do k3=0,novec(3)-1
                         kvec(3) = k3
                         do k4=0,novec(4)-1
                            kvec(4) = k4;
                            do j=1,interpolator%dim
                               call set_hierarchy_info(interpolator,counter,j,&
                                    lvec,kvec,novec);
                            end do
                            counter = counter +1;
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    interpolator%level_mapping(interpolator%max_level+1) = counter;

    ! Now rescale all the coordinates to the actual mesh size
    do i=1,interpolator%size_basis
       do j=1,interpolator%dim
          interpolator%hierarchy(i)%coordinate(j) = &
               interpolator%hierarchy(i)%coordinate(j)* &
               interpolator%length(j) + &
               interpolator%eta_min(j)
       end do
    end do



  end subroutine initialize_sg4d

!> For a given sparse grid point fill the hierarchy information (4D specific)
subroutine set_hierarchy_info(interpolator,counter,cdim,lvecin,kvecin,novecin)
  class(sparse_grid_interpolator_4d), intent(inout) :: interpolator
  sll_int32 :: ld !< current level along dimension \a dim
  sll_int32 :: kd !< current index within level along dimension \a dim
  sll_int32,intent(in) :: cdim !< dimention to be filled
  sll_int32,intent(in) :: counter !< counter for node
  sll_int32, dimension(:), intent(in) :: lvecin !< level vector
  sll_int32, dimension(:), intent(in) :: kvecin !< vector with level within index
  sll_int32, dimension(:), intent(in) :: novecin !< no. of points along each dimension on current level

  sll_int32, dimension(4) :: lvec,kvec,novec
  sll_int32 :: jj,stride

  do jj=1,interpolator%dim
     lvec(jj) = lvecin(jj);
     kvec(jj) = kvecin(jj);
     novec(jj) = novecin(jj);
  end do
  ld = lvec(cdim);
  kd = kvec(cdim);

  interpolator%hierarchy(counter)%level(cdim) = ld;
  interpolator%hierarchy(counter)%index_on_level(cdim) = kd;

  stride = cdim*2-1
  if (ld==0) then
     interpolator%hierarchy(counter)%coordinate(cdim) = 0.0_f64
     interpolator%hierarchy(counter)%parent(stride) = &
          counter
     interpolator%hierarchy(counter)%parent(stride+1) = &
          counter
     interpolator%hierarchy(counter)%function_type(cdim) = 0
  else
     interpolator%hierarchy(counter)%coordinate(cdim) = &
          1.0_f64/(2.0_f64**ld)+kd*1.0_f64/(2.0_f64**(ld-1))


     lvec(cdim) = lvec(cdim)-1;
     novec(cdim) = max(novec(cdim)/2,1);

     ! This one is actually only a neighboring point not the direct parent.
     kvec(cdim) = modulo((kd+1)/2,max(2**(ld-2),1));
     interpolator%hierarchy(counter)%parent(&
          modulo(kd,2)+stride) = &
          interpolator%hierarchy(&
          interpolator%index(lvec(1),lvec(2),lvec(3),lvec(4))+&
          kvec(1)*novec(2)*novec(3)*novec(4)+&
          kvec(2)*novec(3)*novec(4)+&
          kvec(3)*novec(4)+kvec(4))%parent(stride)
     ! This is the actual parent.
     kvec(cdim) = kd/2;
     interpolator%hierarchy(counter)%parent(&
          modulo(kd+1,2)+stride) = &
          interpolator%index(lvec(1),lvec(2),lvec(3),lvec(4))+&
          kvec(1)*novec(2)*novec(3)*novec(4)+kvec(2)*novec(3)*novec(4)+&
          kvec(3)*novec(4)+kvec(4)
     ! Now tell my parent that I am his child.
     interpolator%hierarchy(&
          interpolator%hierarchy(counter)%parent(&
          modulo(kd+1,2)+stride))%children(modulo(kd,2)+stride) = counter
     if(ld==1) then
        interpolator%hierarchy(&
             interpolator%hierarchy(counter)%parent(&
             modulo(kd+1,2)+stride))%children(modulo(kd+1,2)+stride) = counter
     end if
     if (interpolator%order == 1) then
        interpolator%hierarchy(counter)%&
             function_type(cdim) = -1
     elseif (ld==1 .OR. interpolator%order==2) then
        interpolator%hierarchy(counter)%&
             function_type(cdim) = 1+modulo(kd,2)
     elseif (ld==2 .OR. interpolator%order==3)then !(order==3) then!
        interpolator%hierarchy(counter)%&
             function_type(cdim) = 3 + modulo(kd,4)
     else
        interpolator%hierarchy(counter)%&
             function_type(cdim) = 7 + modulo(kd,8)
     end if
  end if

end subroutine set_hierarchy_info

!!!! End initialization routines !!!!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!PN DEFINED BUT NOT USED
!!> Functions to evaluate fg on sg and sg on fg
!subroutine fg_to_sg(interpolator,fg_values,sg_values)
!sll_real64, dimension(:,:,:,:), intent(in) :: fg_values
!sll_real64, dimension(:), intent(out) :: sg_values
!class(sparse_grid_interpolator_4d), intent(in) :: interpolator
!sll_int32 :: j
!sll_int32, dimension(4) :: fg_ind
!
!do j=1,interpolator%size_basis
!   fg_ind = fg_index(interpolator,j);
!   sg_values(j) = fg_values(fg_ind(1),fg_ind(2),fg_ind(3), fg_ind(4));
!end do
!
!end subroutine fg_to_sg


!> Compute the index of a sparse grid node on level "level" with index "index_on_level" on full grid with of max_level
function fg_index(interpolator,sg_index)  
sll_int32, intent(in) :: sg_index
sll_int32, dimension(4) :: fg_index
class(sparse_grid_interpolator_4d), intent(in) :: interpolator
sll_int32 :: j

do j=1,interpolator%dim
   fg_index(j) = 2**(interpolator%levels(j)-&
        interpolator%hierarchy(sg_index)%level(j))*&
        (1+2*interpolator%hierarchy(sg_index)%index_on_level(j)) + 1;
end do

end function fg_index



! End functions fg_to_sg and sg_to_fg
!------------------------------------------------------------------------------!



end module sll_m_sparse_grid_4d
