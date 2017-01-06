!> @ingroup sparse_grid
!> @author Katharina Kormann, IPP 
!> @brief Implementation of a 3D sparse grid with interpolation routines.
!> @todo Implement the optimized interpolation routines for option boundary=1
!> @details <DETAILED_DESCRIPTION>

module sll_m_sparse_grid_3d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_sparse_grid_interpolator, only: &
    sll_t_sparse_grid_interpolator

  implicit none

  public :: &
    sll_t_sparse_grid_interpolator_3d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Sparse grid object for 3d with interpolation routines.
type, extends(sll_t_sparse_grid_interpolator) :: sll_t_sparse_grid_interpolator_3d
sll_int32, dimension(:,:,:), pointer  :: index !< 3d mapping: for each 3d index l on the sparse grid, \a index gives the index of the first node belonging to this level  

contains
  procedure :: initialize => initialize_sg3d! Initialization routine
  procedure :: interpolate_from_interpolant_value ! Compute the value of the sparse grid interpolant at position eta
  procedure :: interpolate_const_disp
  procedure :: fg_to_sg
  procedure :: SPFFT
  procedure :: ISPFFT
  procedure :: interpolate_array_disp_sgfft

end type sll_t_sparse_grid_interpolator_3d

contains




!------------------------------------------------------------------------------!
!!!! Interpolation routines !!!!


!> Compute the value of the sparse grid interpolant at position \a eta (using standard sparse grid interpolation)
  function interpolate_from_interpolant_value( interpolator,data, eta ) result(val)
    class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator !< sparse grid object
    sll_real64 :: val !< interpolated value
    sll_real64,dimension(:), intent(in) :: data !< Value of hierarchical surplus
    sll_real64,dimension(:), intent(in) :: eta !< Coordinates of the point where to interpolate

    if (interpolator%boundary == 0) then
       val =  interpolate_from_hierarchical_surplus(interpolator,data,eta)
    else
       val = interpolate_from_hierarchical_surplus_boundary(interpolator,data,eta)
    end if
  end function interpolate_from_interpolant_value


!> Interpolation function for interpolation at (constantly) displaced grid
!> points; displacement only in dimension dim. It is another implementation of the
!> base-class function "interpolate_disp". The advantage is that we can not
!> revisit nodes as we do in the recursive dimension-independently-programmed
!> version.
  subroutine interpolate_const_disp(interpolator,dorder,displacement,data_in, data_out,hiera)
    class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator !< sparse grid object
  sll_real64, dimension(:), intent(inout) :: data_in !< hierarchical surplus
  sll_real64, dimension(:), intent(out) :: data_out !< Value of the function or the hierarchical surplus (depending on value of \a hiera) for the displaced data points.
  sll_int32, dimension(:), intent(in) :: dorder !< dorder(1) gives the dimension along which we have the displacement; dorder(2:3) give the remaining dimensions
  sll_real64, intent(in) ::displacement !< Constant displacement along dimension dorder(1)
  logical, intent(in) :: hiera !< If the result should be the hierarchical surplus, define \a hiera = .TRUE.; if the result should be the function values at the data points give \a hiera = .FALSE.

  sll_int32 :: i1,i2,i3,k2,k3,counter,j
  sll_int32, dimension(3) :: l, no,ind_order
  sll_int32, dimension(:,:), allocatable :: ind

  SLL_ALLOCATE(ind(interpolator%max_level+1,3), i1);

  ind_order(dorder(1)) = 0
  no(dorder(1)) = 1
  do i3 = 0,interpolator%levels(dorder(3))
     ind_order(dorder(3)) = i3
     l(dorder(3)) = i3
     no(dorder(3)) = max(2**(i3-1),1);
     do i2 = 0,min(interpolator%max_level -i3,interpolator%levels(dorder(2)))
        ind_order(dorder(2)) = i2
        no(dorder(2)) = max(2**(i2-1),1);
        ind(1,dorder(1)) = 0;
        do k2 = 0,no(dorder(2))-1
           ind(ind_order(dorder(2))+1,dorder(2)) = k2;
           do k3 = 0,no(dorder(3))-1
              ind(ind_order(dorder(3))+1,dorder(3)) = k3;
              counter = interpolator%index(&
                   ind_order(1),ind_order(2),ind_order(3))+&
                   ind(ind_order(1)+1,1)*no(2)*no(3)&
                   +ind(ind_order(2)+1,2)*no(3)+&
                   ind(ind_order(3)+1,3)
              ! Evaluate along dorder(1)-stripe
              call interpolator%interpolate_disp_1d_periodic&
                   (displacement,dorder(1),&
                   min(interpolator%levels(dorder(1)),&
                   interpolator%max_level-ind_order(dorder(2))-&
                   ind_order(dorder(3))),counter,data_in,data_out,hiera)
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




! helper functions
!> Implements \a interpolate_from_interpolant_value for periodic sparse grid
 function interpolate_from_hierarchical_surplus( interpolator,data, eta ) result(val)
    class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator
    sll_int32 :: j,l1,l2,l3,level
    sll_real64 :: val
    sll_real64,dimension(:), intent(in) :: data,eta
    sll_real64,dimension(3) :: eta_norm
    sll_real64,dimension(3) :: phi
    sll_int32, dimension(3) :: no,l
    sll_int32,dimension(:,:), allocatable :: ind
    sll_real64 :: scale
    sll_int32 :: index

    SLL_ALLOCATE(ind(0:interpolator%max_level,1:3),j)

    val = 0.0_f64
    ind(0:1,1:3) = 0

    do j=1,interpolator%dim
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
          do l2=max(0,level-l1-interpolator%levels(3)),min(level-l1,interpolator%levels(2))
             l(2) = l2
             no(2) = max(2**(l2-1),1)
             l(3) = level-l1-l2
             l3 = l(3);
             no(3) = max(2**(l(3)-1),1)

             index = interpolator%index(l1,l2,l3)+ind(l1,1)*no(2)*no(3)&
                  +ind(l2,2)*no(3)+ind(l3,3)
             do j=1,3
                call interpolator%basis_function(real(2**(max(l(j),1)),f64)*eta_norm(j)&
                     -real(2*ind(l(j),j),f64)-1.0_f64, phi(j),&
                     interpolator%hierarchy(index)%function_type(j))
             end do
             val = val + data(index)*phi(1)*phi(2)*phi(3)
          end do
       end do
    end do

  end function interpolate_from_hierarchical_surplus

!> implements interpolation from hierarchical surplus (\a interpolate_from_interpolant_value) non-periodic
 function interpolate_from_hierarchical_surplus_boundary( interpolator,data, eta ) result(val)
    class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator
    sll_int32 :: j,l1,l2, l3,level
    sll_real64 :: val
    sll_real64, dimension(:), intent(in) :: data
    sll_real64,dimension(:), intent(in) :: eta
    sll_real64,dimension(3) :: eta_norm
    sll_real64,dimension(3) :: phi
    sll_int32, dimension(3) :: no,l
    sll_int32,dimension(:,:), allocatable :: ind
    sll_real64 :: scale, phi1a, phi2a, phi3a
    sll_int32 :: index

    SLL_ALLOCATE(ind(0:interpolator%max_level,1:interpolator%dim),j)

    val = 0.0_f64
    ind(0:1,1:interpolator%dim) = 0

    do j=1,interpolator%dim
       eta_norm(j) = (eta(j)-interpolator%eta_min(j))/interpolator%length(j)
 
       scale = 0.5_f64
       do level = 2, interpolator%levels(j)
          ind(level,j) = ind(level-1,j)*2
          if (eta_norm(j)> scale*(ind(level,j)+1)) then
             ind(level,j) = ind(level,j)+1
          end if
          scale = scale*0.5_f64
       end do
    end do



    do level = 0, interpolator%max_level
       do l1 = 0 , min(level, interpolator%levels(1))
          l(1) = l1
          if (l1 == 0) then
             no(1) = 2;
          else
             no(1) = 2**(l1-1);
          end if
          do l2 = max(0,level-l1-interpolator%levels(3)), min(level-l1, interpolator%levels(2))
             l(2) = l2
             if (l2==0) then
                no(2) = 2;
             else
                no(2) = max(2**(l2-1),1);
             end if
             l3 = level-l1-l2;
             if (l3==0) then
                no(3) = 2;
             else
                no(3) = max(2**(l3-1),1);
             end if
             if (l1>0) then
                if (l2>0) then
                   if (l3>0) then
                      index = interpolator%index(l1,l2,l3)+&
                           ind(l1,1)*no(2)*no(3)&
                           +ind(l2,2)*no(3) + ind(l3,3)

                      do j=1,interpolator%dim
                         call interpolator%basis_function(real(2**l(j),f64)*eta_norm(j)&
                              -real(2*ind(l(j),j),f64)-1.0_f64, phi(j),&
                              interpolator%hierarchy(index)%function_type(j))
                      end do
                      val = val + data(index)&
                           *phi(1)*phi(2)*phi(3)
                   else ! l3=0
                      index = interpolator%index(l1,l2,l3)+ &
                           ind(l1,1)*no(2)*no(3)&
                           +ind(l2,2)*no(3) + ind(l3,3)
                      call interpolator%basis_function(real(2**l(1),f64)*eta_norm(1)&
                           -real(2*ind(l(1),1),f64)-1.0_f64, phi(1),&
                           interpolator%hierarchy(index)%function_type(1))
                      call interpolator%basis_function(real(2**l(2),f64)*eta_norm(2)&
                           -real(2*ind(l(2),2),f64)-1.0_f64, phi(2),&
                           interpolator%hierarchy(index)%function_type(2))
                      call interpolator%basis_function(eta_norm(3), phi(3), -1)
                      val = val + data(index)*phi(1)*phi(2)*phi(3);
                      call interpolator%basis_function(eta_norm(3)-1.0_f64, phi(3), -1)
                      val = val + data(index+1)*phi(1)*phi(2)*phi(3);
                   end if
                else !l2 = 0
                   if (l3>0) then
                      index = interpolator%index(l1,l2,l3)+ &
                           ind(l1,1)*no(2)*no(3)&
                           +ind(l2,2)*no(3) + ind(l3,3)
                      call interpolator%basis_function(real(2**l(1),f64)*eta_norm(1)&
                           -real(2*ind(l(1),1),f64)-1.0_f64, phi(1),&
                           interpolator%hierarchy(index)%function_type(1))
                      call interpolator%basis_function(real(2**l(3),f64)*eta_norm(3)&
                           -real(2*ind(l(3),3),f64)-1.0_f64, phi(3),&
                           interpolator%hierarchy(index)%function_type(3))
                      call interpolator%basis_function(eta_norm(2), phi(2), -1)
                      val = val + data(index)*phi(1)*phi(2)*phi(3);
                      call interpolator%basis_function(eta_norm(2)-1.0_f64, phi(2), -1)
                      val = val + data(index+no(3))*phi(1)*phi(2)*phi(3);
                   else !l3=l2=0 
                      index =  interpolator%index(l1,l2,l3)+ &
                           ind(l1,1)*no(2)*no(3)+ind(l2,2)*no(3) + ind(l3,3);
                      call interpolator%basis_function(real(2**l(1),f64)*eta_norm(1)&
                           -real(2*ind(l(1),1),f64)-1.0_f64, phi(1),&
                           interpolator%hierarchy(index)%function_type(1))
                      call interpolator%basis_function(eta_norm(3), phi(3), -1)
                      call interpolator%basis_function(eta_norm(3)-1.0_f64, phi3a, -1)
                      call interpolator%basis_function(eta_norm(2), phi(2), -1)
                      call interpolator%basis_function(eta_norm(2)-1.0_f64, phi2a, -1)
                      val = val + data(index)*phi(1)*phi(2)*phi(3) + &
                           data(index+1)*phi(1)*phi(2)*phi3a + &
                           data(index+2)*phi(1)*phi2a*phi(3) + &
                           data(index+3)*phi(1)*phi2a*phi3a
                   end if
                end if
             else ! l1=0 
                 if (l2>0) then
                   if (l3>0) then
                      index = interpolator%index(l1,l2,l3)+&
                           ind(l1,1)*no(2)*no(3)&
                           +ind(l2,2)*no(3) + ind(l3,3)

                      do j=2,interpolator%dim
                         call interpolator%basis_function(real(2**l(j),f64)*eta_norm(j)&
                              -real(2*ind(l(j),j),f64)-1.0_f64, phi(j),&
                              interpolator%hierarchy(index)%function_type(j))
                      end do
                      call interpolator%basis_function(eta_norm(1), phi(1), -1)
                      val = val + data(index)*phi(1)*phi(2)*phi(3);
                      call interpolator%basis_function(eta_norm(1)-1.0_f64, phi(1), -1)
                      val = val + data(index+no(2)*no(3))*phi(1)*phi(2)*phi(3)
                   else ! l1=l3=0
                      index = interpolator%index(l1,l2,l3)+ &
                           ind(l1,1)*no(2)*no(3)&
                           +ind(l2,2)*no(3) + ind(l3,3)
                      call interpolator%basis_function(eta_norm(1), phi(1), -1)            
                      call interpolator%basis_function(eta_norm(1)-1.0_f64, phi1a, -1)
                      call interpolator%basis_function(real(2**l(2),f64)*eta_norm(2)&
                           -real(2*ind(l(2),2),f64)-1.0_f64, phi(2),&
                           interpolator%hierarchy(index)%function_type(2))
                      call interpolator%basis_function(eta_norm(3), phi(3), -1)
                      call interpolator%basis_function(eta_norm(3)-1.0_f64, phi3a, -1)
                      val = val + data(index)*phi(1)*phi(2)*phi(3)+&
                           data(index+1)*phi(1)*phi(2)*phi3a + &
                           data(index+no(2)*no(3))*phi1a*phi(2)*phi(3) + &
                           data(index+no(2)*no(3)+1)*phi1a*phi(2)*phi3a
                   end if
                else !l1=l2 = 0!AB hier noch nicht ganz ueberarbeitet
                   if (l3>0) then
                      index = interpolator%index(l1,l2,l3)+ &
                           ind(l1,1)*no(2)*no(3)&
                           +ind(l2,2)*no(3) + ind(l3,3)
                      call interpolator%basis_function(eta_norm(1), phi(1), -1)            
                      call interpolator%basis_function(eta_norm(1)-1.0_f64, phi1a, -1)
                      call interpolator%basis_function(real(2**l(3),f64)*eta_norm(3)&
                           -real(2*ind(l(3),3),f64)-1.0_f64, phi(3),&
                           interpolator%hierarchy(index)%function_type(3))
                      call interpolator%basis_function(eta_norm(2), phi(2), -1)
                      call interpolator%basis_function(eta_norm(2)-1.0_f64, phi2a, -1)
                      val = val + data(index)*phi(1)*phi(2)*phi(3)+&
                           data(index+no(3))*phi(1)*phi2a*phi(3) + &
                           data(index+no(2)*no(3))*phi1a*phi(2)*phi(3) + &
                           data(index+no(2)*no(3)+no(3))*phi1a*phi2a*phi(3);
                   else !l3=l2=0 
                      index =  interpolator%index(l1,l2,l3)+ &
                           ind(l1,1)*no(2)*no(3)+ind(l2,2)*no(3) + ind(l3,3);
                      call interpolator%basis_function(eta_norm(1), phi(1), -1)            
                      call interpolator%basis_function(eta_norm(1)-1.0_f64, phi1a, -1)
                      call interpolator%basis_function(eta_norm(3), phi(3), -1)
                      call interpolator%basis_function(eta_norm(3)-1.0_f64, phi3a, -1)
                      call interpolator%basis_function(eta_norm(2), phi(2), -1)
                      call interpolator%basis_function(eta_norm(2)-1.0_f64, phi2a, -1)
                      val = val + data(index)*phi(1)*phi(2)*phi(3) + &
                           data(index+1)*phi(1)*phi(2)*phi3a + &
                           data(index+2)*phi(1)*phi2a*phi(3) + &
                           data(index+3)*phi(1)*phi2a*phi3a + &
                           data(index+4)*phi1a*phi(2)*phi(3) + &
                           data(index+5)*phi1a*phi(2)*phi3a + &
                           data(index+6)*phi1a*phi2a*phi(3) + &
                           data(index+7)*phi1a*phi2a*phi3a
                   end if
                end if
             end if
          end do
       end do
    end do
  end function interpolate_from_hierarchical_surplus_boundary





!!!! End interpolation routines !!!!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!!!! SGFFT routines !!!!!

!>  Compute value at displaced grid points using trigonometric interpolation (based on SG FFT)
subroutine interpolate_array_disp_sgfft(interpolator,dim, displacment_in,data_in,data_out)
  class(sll_t_sparse_grid_interpolator_3d),  intent(inout)       :: interpolator !< Sparse grid object
  sll_int32, intent(in) :: dim !< dimension along which the points should be displaced with \a displacement_in
  sll_real64, intent(in) :: displacment_in !< displacement of the data points along dimension \a dim
  sll_comp64, dimension(:), intent(inout)   :: data_in !< Fourier transformed values on the sparse grid
  sll_real64, dimension(:), intent(out) :: data_out !< Function values on the sparse grid after displacement
  sll_real64:: displacement
  
  displacement = displacment_in*2.0_f64*sll_p_pi/interpolator%length(dim)
  call Displace(interpolator,dim,displacement,data_in);
  call ISPFFT(interpolator,data_in,data_out);
  
end subroutine interpolate_array_disp_sgfft


!!!! End SGFFT routines !!!!!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!!!! Initialization routines !!!!


  !> Initialization function. Set up the hierarchy of the sparse grid
  subroutine initialize_sg3d( &
    interpolator, &
    levels, &
    order, &
    interpolation, &
    interpolation_type, &
    eta_min, &
    eta_max, &
    boundary, &
    modified)

    class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator !< sparse grid object
    sll_real64, dimension(:), intent(in)              :: eta_min !< \a eta_min defines the lower bound of the domain
    sll_real64, dimension(:),  intent(in)             :: eta_max !< \a eta_max defines the upper bound of the domain
    sll_int32, dimension(:),intent(in)                :: levels 
    !< \a levels defines the maximum level in the sparse grid for each direction
    sll_int32, intent(in)                             :: order !< \a order of the sparse grid functions
    sll_int32, intent(in)                             :: interpolation !< order of the interpolator
    sll_int32, intent(in)                             :: interpolation_type !< Choose spline (\a interpolation_type = 0) or Lagrange (\a interpolation_type = 1) interpolation for the 1D interpolators if not traditional sparse grid interpolation is used.
    sll_int32, intent(in)                             :: modified 
    !< \a modified defines if we have a traditional sparse grid for \a modified
    !< = 0 (then the l_1 norm of the levels is bounded by max(\a levels) ) or if the
    !< boundary is sparsified for \a modified = 1 (then the l_1 norm of the levels is
    !< bounded by max(\a levels)+1 )
    sll_int32, intent(in)                             :: boundary 
    !< \a boundary defines the boundary conditions: define 0 for periodic boundary conditions and 1 for zero inflow boundaries

    sll_int32                                     :: i,j, k1,k2,k3,l1,l2,l3,l,counter
    sll_int32                                     :: ierr, no1, no2, no3
    sll_int32, dimension(:) ,allocatable          :: novec,lvec,kvec

    interpolator%dim = 3;
    interpolator%modified = modified;
    interpolator%boundary = boundary;
    SLL_ALLOCATE(lvec(interpolator%dim),ierr);
    SLL_ALLOCATE(kvec(interpolator%dim),ierr);
    SLL_ALLOCATE(novec(interpolator%dim),ierr);

    interpolator%max_level = levels(1);
    do l=2,interpolator%dim
       interpolator%max_level = max(levels(l),interpolator%max_level);
    end do  
    if (interpolator%modified == 1) then
       interpolator%max_level = interpolator%max_level + 2;
    end if
    
    interpolator%size_basis = 0;
    if (interpolator%boundary == 0) then
       do l = 0, interpolator%max_level
          do l1 = 0 , min(l,levels(1))
             do l2 = max(0,l-l1-levels(3)), min(l-l1,levels(2))
                l3 = l-l1-l2
                interpolator%size_basis = &
                     interpolator%size_basis + &
                     max(2**(l1-1),1)*max(2**(l2-1),1)*max(2**(l3-1),1)
             end do
          end do
       end do
    else
        do l = 0, interpolator%max_level
          do l1 = 0, min(l, levels(1))
             do l2 = max(0, l-l1-levels(3)) , min(l-l1, levels(2))
                l3 = l-l1-l2
                if( l1 == 0) then
                   no1 = 2;
                else
                   no1 = 2**(l1-1);
                end if
                if( l2 == 0) then
                   no2 = 2;
                else
                   no2 = 2**(l2-1);
                end if
                if (l3 == 0) then
                   no3 = 2;
                else
                   no3 = 2**(l3-1);
                end if
                interpolator%size_basis = &
                     interpolator%size_basis + no1*no2*no3;
             end do
          end do
       end do
    end if

    call  interpolator%initialize_sg(levels, order, interpolation,&
    interpolation_type, eta_min, eta_max);

#ifdef __PGI
    allocate(interpolator%index(0:interpolator%levels(1),0:interpolator%levels(2),0:interpolator%levels(3)))
#else
    SLL_ALLOCATE(interpolator%index(0:interpolator%levels(1),0:interpolator%levels(2),0:interpolator%levels(3)),ierr)
#endif

    ! Set the hierarchy of the grid
    counter = 1
    do l = 0, interpolator%max_level
       interpolator%level_mapping(l) = counter;
       do l1 = 0 , min(l,interpolator%levels(1))
          novec(1) = max(2**(l1-1),1)
          if (interpolator%boundary == 1 .AND. l1 == 0) then
             novec(1) = 2;
          end if
          lvec(1) = l1;
          do l2 = max(0,l-l1-interpolator%levels(3)), min(l-l1,interpolator%levels(2))
             novec(2) = max(2**(l2-1),1)
             if (interpolator%boundary == 1 .AND. l2 == 0) then
                novec(2) = 2;
             end if
             lvec(2) = l2;
             lvec(3) = l-l1-l2;
             l3 = lvec(3);
             novec(3) = max(2**(l3-1),1)
             if (interpolator%boundary == 1 .AND. l3 == 0) then
                novec(3) = 2;
             end if
             lvec(3) = l3
             interpolator%index(l1,l2,l3) = counter
             do k1=0,novec(1)-1
                kvec(1) = k1
                do k2=0,novec(2)-1
                   kvec(2) = k2
                   do k3=0,novec(3)-1
                      kvec(3) = k3
                      do j=1,interpolator%dim
                         if (interpolator%boundary == 0) then
                            call set_hierarchy_info(interpolator,counter,j,&
                                 lvec,kvec,novec);
                         else
                             call set_hierarchy_info_boundary&
                                 (interpolator,counter,j,lvec,kvec,novec);
                         end if
                      end do
                      counter = counter +1;
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



  end subroutine initialize_sg3d

!> Helfer function for initialization. Setting all the information needed for node \a counter of the sparse grid along dimension \a cdim
!> For a given sparse grid point fill the hierarchy information (3D specific)
subroutine set_hierarchy_info(interpolator,counter,cdim,lvecin,kvecin,novecin)
    class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator
  sll_int32 :: ld !< current level
  sll_int32 :: kd !< current index within level
  sll_int32,intent(in) :: cdim !< current dimension
  sll_int32,intent(in) :: counter !< counter for node  
  sll_int32, dimension(:), intent(in) :: lvecin !< vector of current levels
  sll_int32, dimension(:), intent(in) :: kvecin !< vector of current index within level
  sll_int32, dimension(:), intent(in) :: novecin !< vector with number of points on the current level along each dimension
  
  sll_int32, dimension(3) :: lvec,kvec,novec
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
          interpolator%index(lvec(1),lvec(2),lvec(3))+&
          kvec(1)*novec(2)*novec(3)+&
          kvec(2)*novec(3)+&
          kvec(3))%parent(stride)
     ! This is the actual parent.
     kvec(cdim) = kd/2;
     interpolator%hierarchy(counter)%parent(&
          modulo(kd+1,2)+stride) = &
          interpolator%index(lvec(1),lvec(2),lvec(3))+&
          kvec(1)*novec(2)*novec(3)+kvec(2)*novec(3)+&
          kvec(3)
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

!> Helfer function for initialization. Setting all the information needed for node \a counter of the sparse grid along dimension \a cdim for points at the boundary along dimension \a dim.
subroutine set_hierarchy_info_boundary(interpolator,counter,cdim,lvecin,kvecin,novecin)
  class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator
  sll_int32 :: ld !< current level
  sll_int32 :: kd !< current index within level
  sll_int32,intent(in) :: cdim !< current dimension
  sll_int32,intent(in) :: counter !< counter for node
  sll_int32, dimension(:), intent(in) :: lvecin !< vector of current levels
  sll_int32, dimension(:), intent(in) :: kvecin !< vector of current index within level
  sll_int32, dimension(:), intent(in) :: novecin !< vector with number of points on the current level along each dimension

  sll_int32, dimension(:), allocatable :: lvec,kvec,novec
  sll_int32 :: jj,stride

  SLL_ALLOCATE(lvec(interpolator%dim),jj);
  SLL_ALLOCATE(kvec(interpolator%dim),jj);
  SLL_ALLOCATE(novec(interpolator%dim),jj);

  do jj=1,interpolator%dim
     lvec(jj) = lvecin(jj);
     kvec(jj) = kvecin(jj);
     novec(jj) = novecin(jj);
  end do
  ld = lvec(cdim);
  kd = kvec(cdim);

  interpolator%hierarchy(counter)%level(cdim) = ld;

  stride = cdim*2-1
  if (ld==0) then
     interpolator%hierarchy(counter)%coordinate(cdim) = real(kd,f64);
     interpolator%hierarchy(counter)%parent(stride) = &
          counter
     interpolator%hierarchy(counter)%parent(stride+1) = &
          counter
     interpolator%hierarchy(counter)%function_type(cdim) = -1
  else
     interpolator%hierarchy(counter)%coordinate(cdim) = &
          1.0_f64/(2.0_f64**ld)+kd*1.0_f64/(2.0_f64**(ld-1))


     lvec(cdim) = lvec(cdim)-1;
     novec(cdim) = novec(cdim)/2;

     if(ld==1) then
        kvec(cdim) = kd/2;
        novec(cdim) = 2;
        interpolator%hierarchy(counter)%parent(stride) = &
             interpolator%index(lvec(1),lvec(2),lvec(3))+&
             kvec(1)*novec(2)*novec(3)+kvec(2)*novec(3)+kvec(3)
        if (cdim==1) then
           interpolator%hierarchy(counter)%parent(&
                stride+1) = interpolator%hierarchy(counter)%parent(&
                stride) + novec(2)*novec(3)
        elseif(cdim==2) then
           interpolator%hierarchy(counter)%parent(&
                stride+1) = interpolator%hierarchy(counter)%parent(&
                stride) + novec(3)
        else
           interpolator%hierarchy(counter)%parent(&
             stride+1) = interpolator%hierarchy(counter)%parent(&
             stride) + 1
        end if
        interpolator%hierarchy(interpolator%hierarchy(counter)%parent(stride))%children&
             (stride) = counter
        !interpolator%hierarchy(interpolator%hierarchy(counter)%parent(stride+1))%&
        !     children(stride) = counter
     else
        ! This one is actually only a neighboring point not the direct parent.
        kvec(cdim) =modulo((kd+1)/2,max(2**(ld-2),1));
        kvec(cdim) = (kd+1)/2;
        if (kvec(cdim)<novec(cdim)) then
           interpolator%hierarchy(counter)%parent(&
                modulo(kd,2)+stride) = &
                interpolator%hierarchy(&
                interpolator%index(lvec(1),lvec(2),lvec(3))+&
                kvec(1)*novec(2)*novec(3)+kvec(2)*novec(3)+kvec(3))%&
                parent(stride)
        else
           kvec(cdim) = kvec(cdim)-1;
           interpolator%hierarchy(counter)%parent(&
                modulo(kd,2)+stride) = &
                interpolator%hierarchy(&
                interpolator%index(lvec(1),lvec(2),lvec(3))+&
                kvec(1)*novec(2)*novec(3)+kvec(2)*novec(3)+kvec(3))%&
                parent(stride+1)
        end if
        ! This is the actual parent.
        kvec(cdim) = kd/2;
        interpolator%hierarchy(counter)%parent(&
             modulo(kd+1,2)+stride) = &
             interpolator%index(lvec(1),lvec(2),lvec(3))+&
             kvec(1)*novec(2)*novec(3)+kvec(2)*novec(3)+kvec(3)
        ! Now tell my parent that I am his child.
        interpolator%hierarchy(&
             interpolator%hierarchy(counter)%parent(&
             modulo(kd+1,2)+stride))%children(modulo(kd,2)+stride) = counter
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

end subroutine set_hierarchy_info_boundary


!!!! End initialization routines !!!!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!> Functions to evaluate fg on sg and sg on fg
subroutine fg_to_sg(interpolator,fg_values,sg_values)
sll_real64, dimension(:,:,:), intent(in) :: fg_values
sll_real64, dimension(:), intent(out) :: sg_values
class(sll_t_sparse_grid_interpolator_3d), intent(in) :: interpolator
sll_int32 :: j
sll_int32, dimension(3) :: fg_ind

do j=1,interpolator%size_basis
   fg_ind = fg_index(interpolator,j);
   sg_values(j) = fg_values(fg_ind(1),fg_ind(2),fg_ind(3));
end do

end subroutine fg_to_sg


!> Compute the index of a sparse grid node on level "level" with index "index_on_level" on full grid with of max_level
function fg_index(interpolator,sg_index)  
sll_int32, intent(in) :: sg_index
sll_int32, dimension(3) :: fg_index
class(sll_t_sparse_grid_interpolator_3d), intent(in) :: interpolator
sll_int32 :: j

do j=1,interpolator%dim
   fg_index(j) = 2**(interpolator%levels(j)-&
        interpolator%hierarchy(sg_index)%level(j))*&
        (1+2*interpolator%hierarchy(sg_index)%index_on_level(j)) + 1;
end do

end function fg_index



! End functions fg_to_sg and sg_to_fg
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!!!! SGFFT helper functions (Some clean-up needed) !!!!


subroutine ToHierarchical(interpolator,data_in, data_out)
  class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(in) :: data_in
  sll_comp64, dimension(:), intent(out) :: data_out
  sll_int32 :: i1,i2,i3,j
  
! Maybe possible with index in index to make dimension independent
  do i2 = 0,interpolator%levels(2)
     do i3 = 0, min(interpolator%max_level - i2,interpolator%levels(3))
        do j = interpolator%index(0,i2,i3),&
             interpolator%index(0,i2,i3)+max(2**(i2-1),1)*max(2**(i3-1),1)-1
           call interpolator%ToHierarchical1D(1,&
                min(interpolator%levels(1),interpolator%max_level-i2-i3),&
             j,data_in,data_out)
        end do
     end do
  end do

  do i1 = 0,interpolator%levels(1)
     do i3 = 0, min(interpolator%max_level - i1,interpolator%levels(3))
        do j = interpolator%index(i1,0,i3),&
             interpolator%index(i1,0,i3)+max(2**(i1-1),1)*max(2**(i3-1),1)-1
           call interpolator%ToHierarchical1D_comp &
                (2,&
                min(interpolator%levels(2),interpolator%max_level-i1-i3),j,&
                data_out)
        end do
     end do
  end do

  do i1 = 0,interpolator%levels(1)
     do i2 = 0, min(interpolator%max_level - i1,interpolator%levels(2))
        do j = interpolator%index(i1,i2,0),&
             interpolator%index(i1,i2,0)+max(2**(i1-1),1)*max(2**(i2-1),1)-1
           call interpolator%ToHierarchical1D_comp &
                (3,&
                min(interpolator%levels(3),interpolator%max_level-i1-i2),j,&
                data_out)
        end do
     end do
  end do

end subroutine ToHierarchical


 subroutine ToDehi(interpolator,data_array)
  class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data_array
  sll_int32 :: i1,i2,i3,j

  do i1 = 0,interpolator%levels(1)
     do i2 = 0, min(interpolator%max_level - i1,interpolator%levels(2))
        do j = interpolator%index(i1,i2,0),&
             interpolator%index(i1,i2,0)+max(2**(i1-1),1)*max(2**(i2-1),1)-1
           call interpolator%ToDehi1D &
                (3,&
                min(interpolator%levels(3),interpolator%max_level-i1-i2),j,&
                data_array)
        end do
     end do
  end do

  do i1 = 0,interpolator%levels(1)
     do i3 = 0, min(interpolator%max_level - i1,interpolator%levels(3))
        do j = interpolator%index(i1,0,i3),&
             interpolator%index(i1,0,i3)+max(2**(i1-1),1)*max(2**(i3-1),1)-1
           call interpolator%ToDehi1D &
                (2,&
                min(interpolator%levels(2),interpolator%max_level-i1-i3),j,&
                data_array)
        end do
     end do
  end do

  do i2 = 0,interpolator%levels(2)
     do i3 = 0, min(interpolator%max_level - i2,interpolator%levels(3))
        do j = interpolator%index(0,i2,i3),&
             interpolator%index(0,i2,i3)+max(2**(i2-1),1)*max(2**(i3-1),1)-1
           call interpolator%ToDehi1D(1,&
                min(interpolator%levels(1),interpolator%max_level-i2-i3),&
             j,data_array)
        end do
     end do
  end do



end subroutine ToDehi


subroutine ToHira(interpolator,data_array)
  class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data_array
  sll_int32 :: i1,i2,i3,j


  do i2 = 0,interpolator%levels(2)
     do i3 = 0, min(interpolator%max_level - i2,interpolator%levels(3))
        do j = interpolator%index(0,i2,i3),&
             interpolator%index(0,i2,i3)+max(2**(i2-1),1)*max(2**(i3-1),1)-1
           call interpolator%ToHira1D(1,&
                min(interpolator%levels(1),interpolator%max_level-i2-i3),&
                j,data_array)
        end do
     end do
  end do

  do i1 = 0,interpolator%levels(1)
     do i3 = 0, min(interpolator%max_level - i1,interpolator%levels(3))
        do j = interpolator%index(i1,0,i3),&
             interpolator%index(i1,0,i3)+max(2**(i1-1),1)*max(2**(i3-1),1)-1
           call interpolator%ToHira1D &
                (2,&
                min(interpolator%levels(2),interpolator%max_level-i1-i3),j,&
                data_array)
        end do
     end do
  end do

  do i1 = 0,interpolator%levels(1)
     do i2 = 0, min(interpolator%max_level - i1,interpolator%levels(2))
        do j = interpolator%index(i1,i2,0),&
             interpolator%index(i1,i2,0)+max(2**(i1-1),1)*max(2**(i2-1),1)-1
           call interpolator%ToHira1D &
                (3,&
                min(interpolator%levels(3),interpolator%max_level-i1-i2),j,&
                data_array)
        end do
     end do
  end do

  
end subroutine ToHira

subroutine ToNodal(interpolator,data_in,data_out)
  class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator  
  sll_comp64, dimension(:), intent(inout) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_int32 :: i1,i2,i3,j

  do i1 = 0,interpolator%levels(1)
     do i2 = 0, min(interpolator%max_level - i1,interpolator%levels(2))
        do j = interpolator%index(i1,i2,0),&
             interpolator%index(i1,i2,0)+max(2**(i1-1),1)*max(2**(i2-1),1)-1
           call interpolator%ToNodal1D_comp &
                (3,&
                min(interpolator%levels(3),interpolator%max_level-i1-i2),j,&
                data_in)
        end do
     end do
  end do

  do i1 = 0,interpolator%levels(1)
     do i3 = 0, min(interpolator%max_level - i1,interpolator%levels(3))
        do j = interpolator%index(i1,0,i3),&
             interpolator%index(i1,0,i3)+max(2**(i1-1),1)*max(2**(i3-1),1)-1
           call interpolator%ToNodal1D_comp &
                (2,&
                min(interpolator%levels(2),interpolator%max_level-i1-i3),j,&
                data_in)
        end do
     end do
  end do

  do i2 = 0,interpolator%levels(2)
     do i3 =0, min(interpolator%max_level - i2,interpolator%levels(3))
        do j = interpolator%index(0,i2,i3),&
             interpolator%index(0,i2,i3)+max(2**(i2-1),1)*max(2**(i3-1),1)-1
           call interpolator%ToNodal1D(1,&
                min(interpolator%levels(1),interpolator%max_level-i2-i3),&
                j,data_in,data_out)
        end do
     end do
  end do
 
end subroutine ToNodal


subroutine Displace(interpolator,dim,displacement,data)
  class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data
  sll_real64,intent(in) :: displacement
  sll_int32, intent(in) :: dim
  sll_int32 :: i1,i2,i3,j

  if (dim == 1) then
     do i2 = 0,interpolator%levels(2)
        do i3 = 0 ,min(interpolator%max_level - i2,interpolator%levels(3))
           do j = interpolator%index(0,i2,i3),&
                interpolator%index(0,i2,i3)+max(2**(i2-1),1)*max(2**(i3-1),1)-1
              call interpolator%Displace1D(1,interpolator%levels(1)-i2-i3,&
                   j,displacement,data)
           end do
        end do
     end do
  else if (dim == 2) then
     do i1 = 0,interpolator%levels(1)
        do i3 = 0 ,min(interpolator%max_level - i1,interpolator%levels(3))
           do j = interpolator%index(i1,0,i3),&
                interpolator%index(i1,0,i3)+max(2**(i1-1),1)*max(2**(i3-1),1)-1
              call interpolator%Displace1D(2,interpolator%levels(2)-i1-i3,&
                   j,displacement, data)
           end do
        end do
     end do
  else if (dim == 3) then
     do i1 = 0,interpolator%levels(1)
        do i2 = 0,min(interpolator%max_level - i1,interpolator%levels(2))
           do j = interpolator%index(i1,i2,0),&
                interpolator%index(i1,i2,0)+max(2**(i1-1),1)*max(2**(i2-1),1)-1
              call interpolator%Displace1D(3,interpolator%levels(3)-i1-i2,&
                   j,displacement,data)
           end do
        end do
     end do
  end if


end subroutine DISPLACE

!> Sparse grid FFT
subroutine SPFFT(interpolator,data_in,data_out)
  class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(in) :: data_in
  sll_comp64, dimension(:), intent(out) :: data_out

  call ToHierarchical(interpolator,data_in,data_out)
  call ToDehi(interpolator,data_out)

end subroutine SPFFT

!> Sparse grid inverse FFT
subroutine ISPFFT(interpolator,data_in, data_out)
  class(sll_t_sparse_grid_interpolator_3d), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out

  call ToHira(interpolator,data_in)
  call ToNodal(interpolator,data_in,data_out)

end subroutine ISPFFT



!!!! End SGFFT helper functions !!!!
!------------------------------------------------------------------------------!



end module sll_m_sparse_grid_3d
