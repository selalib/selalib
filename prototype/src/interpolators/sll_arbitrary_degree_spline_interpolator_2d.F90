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

module sll_arbitrary_degree_spline_interpolator_2d_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h" 
use sll_module_interpolators_2d_base
  implicit none

  ! in what follows, the direction '1' is in the contiguous memory direction.
  type, extends(sll_interpolator_2d_base) :: arb_deg_2d_interpolator           
     sll_int32  :: num_pts1
     sll_int32  :: num_pts2
     sll_real64 :: eta1_min
     sll_real64 :: eta1_max
     sll_real64 :: eta2_min
     sll_real64 :: eta2_max
     sll_int32  :: bc_left
     sll_int32  :: bc_right
     sll_int32  :: bc_bottom
     sll_int32  :: bc_top
     sll_int32  :: spline_degree1
     sll_int32  :: spline_degree2
     sll_real64, dimension(:), pointer :: knots1
     sll_real64, dimension(:), pointer :: knots2
     ! some knot-like arrays needed by the spli2d_per routine
     sll_real64, dimension(:), pointer :: t1
     sll_real64, dimension(:), pointer :: t2
     sll_int32  :: size_t1
     sll_int32  :: size_t2 
     sll_int64  :: bc_selector ! this is set in initialization
     sll_real64, dimension(:,:), pointer :: coeff_splines
     sll_int32                  :: size_coeffs1
     sll_int32                  :: size_coeffs2

   contains
    procedure, pass(interpolator) :: initialize=>initialize_ad2d_interpolator
    procedure, pass(interpolator) :: set_coefficients => set_coefficients_ad2d
! better: pre-compute-interpolation-information or something...
    procedure :: compute_interpolants => compute_interpolants_ad2d
   ! procedure,  pass(interpolator) :: compute_spline_coefficients => &
    !     compute_spline_coefficients_ad2d
    !procedure, pass:: compute_spline_coefficients =>compute_spline_coefficients_ad2d
    procedure :: interpolate_value => interpolate_value_ad2d
    procedure :: interpolate_derivative_eta1 => interpolate_derivative1_ad2d
    procedure :: interpolate_derivative_eta2 => interpolate_derivative2_ad2d
    procedure, pass:: interpolate_array => interpolate_array_ad2d
    procedure, pass:: interpolate_array_disp => interpolate_2d_array_disp_ad2d
  end type arb_deg_2d_interpolator

  interface delete
     module procedure delete_arbitrary_degree_2d_interpolator
  end interface delete

contains

  subroutine delete_arbitrary_degree_2d_interpolator( interpolator )
    class(arb_deg_2d_interpolator), intent(inout) :: interpolator
    sll_int32 :: ierr
    SLL_DEALLOCATE(interpolator%knots1,ierr)
    SLL_DEALLOCATE(interpolator%knots2,ierr)
    SLL_DEALLOCATE(interpolator%t1,ierr)
    SLL_DEALLOCATE(interpolator%t2,ierr)
    SLL_DEALLOCATE(interpolator%coeff_splines,ierr)
  end subroutine delete_arbitrary_degree_2d_interpolator

  subroutine initialize_ad2d_interpolator( &
    interpolator, &
    num_pts1, &
    num_pts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    bc_left, &
    bc_right, &
    bc_bottom, &
    bc_top, &
    spline_degree1, &
    spline_degree2 )

    class(arb_deg_2d_interpolator), intent(inout) :: interpolator
    sll_int32, intent(in) :: num_pts1
    sll_int32, intent(in) :: num_pts2
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: bc_left
    sll_int32, intent(in) :: bc_right
    sll_int32, intent(in) :: bc_bottom
    sll_int32, intent(in) :: bc_top
    sll_int32, intent(in) :: spline_degree1
    sll_int32, intent(in) :: spline_degree2
    sll_int32 :: ierr
    sll_int32 :: tmp1
    sll_int32 :: tmp2
    sll_int64 :: bc_selector


    ! do some argument checking...
    if(((bc_left  == SLL_PERIODIC).and.(bc_right.ne. SLL_PERIODIC)).or.&
       ((bc_right == SLL_PERIODIC).and.(bc_left .ne. SLL_PERIODIC)))then
       print *, 'initialize_arbitrary_degree_2d_interpolator, ERROR: ', &
            'if one boundary condition is specified as periodic, then ', &
            'both must be. Error in first direction.'
    end if

    if(((bc_bottom == SLL_PERIODIC).and.(bc_top.ne. SLL_PERIODIC)).or.&
       ((bc_top == SLL_PERIODIC).and.(bc_bottom .ne. SLL_PERIODIC)))then
       print *, 'initialize_arbitrary_degree_2d_interpolator, ERROR: ', &
            'if one boundary condition is specified as periodic, then ', &
            'both must be. Error in second direction.'
    end if

    bc_selector = 0

    if( bc_left == SLL_DIRICHLET ) then
       bc_selector = bc_selector + 1
    end if

    if( bc_left == SLL_NEUMANN ) then
       bc_selector = bc_selector + 2
    end if

    if( bc_left == SLL_HERMITE ) then
       bc_selector = bc_selector + 4
    end if

    if( bc_right == SLL_DIRICHLET ) then
       bc_selector = bc_selector + 8
    end if

    if( bc_right == SLL_NEUMANN ) then
       bc_selector = bc_selector + 16
    end if

   if( bc_right == SLL_HERMITE ) then
       bc_selector = bc_selector + 32
    end if

    if( bc_bottom == SLL_DIRICHLET ) then
       bc_selector = bc_selector + 64
    end if

    if( bc_bottom == SLL_NEUMANN ) then
       bc_selector = bc_selector + 128
    end if

    if( bc_bottom == SLL_HERMITE ) then
       bc_selector = bc_selector + 256
    end if

    if( bc_top == SLL_DIRICHLET ) then
       bc_selector = bc_selector + 512
    end if

    if( bc_top == SLL_NEUMANN ) then
       bc_selector = bc_selector + 1024
    end if

   if( bc_top == SLL_HERMITE ) then
       bc_selector = bc_selector + 2048
    end if

    interpolator%spline_degree1 = spline_degree1
    interpolator%spline_degree2 = spline_degree2
    interpolator%eta1_min = eta1_min
    interpolator%eta1_max = eta1_max
    interpolator%eta2_min = eta2_min
    interpolator%eta2_max = eta2_max
    interpolator%bc_left  = bc_left
    interpolator%bc_right = bc_right
    interpolator%bc_left  = bc_bottom
    interpolator%bc_right = bc_top
    interpolator%bc_selector = bc_selector
    interpolator%num_pts1 = num_pts1
    interpolator%num_pts2 = num_pts2

    select case (bc_selector)
    case (0) ! 1. periodic-periodic
       SLL_ALLOCATE( interpolator%knots1(2*spline_degree1+2),ierr )
       SLL_ALLOCATE( interpolator%knots2(2*spline_degree2+2),ierr )
       tmp1 = num_pts1 + 2*spline_degree1
       tmp2 = num_pts2 + 2*spline_degree2
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(2*spline_degree2+2),ierr )
       tmp1 = num_pts1 + spline_degree1 !- 1
       tmp2 = num_pts2 + 2*spline_degree2
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case (576) ! 3. periodic, dirichlet-bottom, dirichlet-top
       SLL_ALLOCATE( interpolator%knots1(2*spline_degree1+2),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       tmp1 = num_pts1 + 2*spline_degree1
       tmp2 = num_pts2 + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case (585) ! 4. dirichlet in all sides
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       tmp1 = num_pts1 + spline_degree1 !- 1
       tmp2 = num_pts2 + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case default
       print *, 'initialize_ad2d_interpolator: BC combination not implemented.'
    end select

    interpolator%coeff_splines(:,:) = 0.0_f64
    SLL_ALLOCATE( interpolator%t1(num_pts1 + 2*(spline_degree1 + 1)), ierr)
    SLL_ALLOCATE( interpolator%t2(num_pts2 + 2*(spline_degree2 + 1)), ierr)

    interpolator%t1(:) = 0.0_f64
    interpolator%t2(:) = 0.0_f64
    !print*,'SIZE',  num_pts1 + 2*(spline_degree1 + 1)
  end subroutine initialize_ad2d_interpolator

!!$  subroutine compute_interpolants_ad2d( &
!!$    interpolator, &
!!$    data_array )
!!$
!!$    class(arb_deg_2d_interpolator), intent(inout) :: interpolator
!!$    sll_real64, dimension(:,:), intent(in)         :: data_array
!!$
!!$    print *, 'compute interpolants not implemented'
!!$
!!$  end subroutine compute_interpolants_ad2d

  subroutine set_coefficients_ad2d( &
   interpolator, &
   coeffs_1d, &
   coeffs_2d )

   class(arb_deg_2d_interpolator), intent(inout)  :: interpolator
   sll_real64, dimension(:), intent(in), optional :: coeffs_1d
   sll_real64, dimension(:,:), intent(in), optional :: coeffs_2d
   sll_int32 :: sp_deg1
   sll_int32 :: sp_deg2
   sll_int32 :: num_cells1
   sll_int32 :: num_cells2
   sll_int32 :: tmp1, tmp2
   sll_int32 :: i, j
   sll_real64 :: eta1_min, eta1_max
   sll_real64 :: eta2_min, eta2_max
   sll_real64 :: delta1
   sll_real64 :: delta2
   sll_int32  ::  nb_spline_eta1
   sll_int32  ::  nb_spline_eta2
   sll_real64 :: eta1
   sll_real64 :: eta2

   if(present(coeffs_2d)) then
      print *, 'set_coefficients_ad2d(), ERROR: this function has not ', &
           'implmented the option to set the coeffcients provided as a 2d ', &
           'array.'
      stop
   end if

   sp_deg1    = interpolator%spline_degree1
   sp_deg2    = interpolator%spline_degree2
   num_cells1 = interpolator%num_pts1 - 1
   num_cells2 = interpolator%num_pts2 - 1
   eta1_min = interpolator%eta1_min
   eta2_min = interpolator%eta2_min
   eta1_max = interpolator%eta1_max
   eta2_max = interpolator%eta2_max
   delta1 = (eta1_max - eta1_min)/num_cells1
   delta2 = (eta2_max - eta2_min)/num_cells2

   tmp1 = (sp_deg1 + 1)/2
   tmp2 = (sp_deg2 + 1)/2
   print*, tmp1,tmp2
   ! The interpretation and further filling of the spline coefficients array
   ! depends on the boundary conditions.
   select case (interpolator%bc_selector)
   case(0) ! periodic-periodic
      
      interpolator%size_coeffs1=  num_cells1 + sp_deg1
      interpolator%size_coeffs2=  num_cells2 + sp_deg2
      interpolator%size_t1 = 2*sp_deg1 + num_cells1 +1 
      interpolator%size_t2 = 2*sp_deg2 + num_cells2 +1
      ! allocation and definition of knots
      do i = -sp_deg1, num_cells1 + sp_deg1
         interpolator%t1( i + sp_deg1 + 1 ) = eta1_min + i*delta1
      end do
      
      do i = -sp_deg2, num_cells2 + sp_deg2
         interpolator%t2( i + sp_deg2 + 1 ) = eta2_min + i*delta2
      end do


      do i = 1,num_cells1
         do j = 1,num_cells2
            interpolator%coeff_splines(i+tmp1,j+tmp2) = &
                 coeffs_1d( i + num_cells1 *(j-1) )
         end do
      end do
      
      do j=1, tmp2
         interpolator%coeff_splines(:,j) = &
              interpolator%coeff_splines(:,num_cells2 + j)
      end do
      
      if(num_cells2 + tmp2 < num_cells2 + sp_deg2 ) then
         do j = tmp2 + 1, sp_deg2
            interpolator%coeff_splines(:, num_cells2 + j) = &
                 interpolator%coeff_splines(:, j)
         end do
      end if

      do i=1, tmp1
         interpolator%coeff_splines(i,:) = &
              interpolator%coeff_splines(num_cells1 + i,:)
      end do
      
      if (num_cells1 + tmp1 < num_cells1 + sp_deg1 ) then
         do i = tmp1 + 1, sp_deg1
            interpolator%coeff_splines(num_cells1 + i,:) = &
                 interpolator%coeff_splines(i,:)
         end do
      end if
   case (9) ! 2. dirichlet-left, dirichlet-right, periodic
      interpolator%size_coeffs1=  num_cells1 + sp_deg1
      interpolator%size_coeffs2=  num_cells2 + sp_deg2
      interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
      interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
      nb_spline_eta1 = num_cells1 + sp_deg1 - 2
      nb_spline_eta2 = num_cells2
      ! allocation and definition of knots
      do i = - sp_deg2, num_cells2 + sp_deg2
         interpolator%t2( i+ sp_deg2 + 1 ) = eta2_min + i* delta2
      end do
      
      do i = 1, sp_deg1 + 1
         interpolator%t1(i) = eta1_min
      enddo
      eta1 = eta1_min
      do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
         eta1 = eta1 + delta1
         interpolator%t1(i) = eta1
      enddo
       do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
          interpolator%t1(i) = eta1
       enddo
       
       
       do i = 1 ,nb_spline_eta1
          do j = 1,nb_spline_eta2
             interpolator%coeff_splines(i+1,j) = &
                  coeffs_1d(i+nb_spline_eta1*(j-1))
         end do
      end do
     
      
      do j = 1, sp_deg2
         do i = 1,nb_spline_eta1

            interpolator%coeff_splines(i + 1 ,nb_spline_eta2 + j ) = &
                 coeffs_1d(i+nb_spline_eta1*(j-1))
         end do
      end do

      interpolator%coeff_splines(1,:) = 0.0_8
      interpolator%coeff_splines(nb_spline_eta1+2,:) = 0.0_8
      
      
   case(576)!3. periodic, dirichlet-bottom, dirichlet-top
      interpolator%size_coeffs1=  num_cells1 + sp_deg1 
      interpolator%size_coeffs2=  num_cells2 + sp_deg2
      interpolator%size_t1 = 2.0_f64*sp_deg1 + num_cells1 + 1
      interpolator%size_t2 = 2.0_f64*sp_deg2 + num_cells2 + 1
      nb_spline_eta1 = num_cells1
      nb_spline_eta2 = num_cells2 + sp_deg2 - 2
      
      ! allocation and definition of knots
      do i = - sp_deg1, nb_spline_eta1 + sp_deg1

         interpolator%t1( i+ sp_deg1 + 1 ) = eta1_min + i* delta1
      end do
      
      
      do i = 1, sp_deg2 + 1
         interpolator%t2(i) = eta2_min
      enddo
      eta2 = eta2_min
      do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
         eta2 = eta2 + delta2
         interpolator%t2(i) = eta2
      enddo
      do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
         interpolator%t2(i) = eta2
      enddo
       

      
      do i = 1 , nb_spline_eta1
         do j = 1,nb_spline_eta2
            
            interpolator%coeff_splines(i ,j+1) = &
                 coeffs_1d(i+nb_spline_eta1 *(j-1) )
         end do
      end do
      
      do i = 1, sp_deg1
         do j = 1,nb_spline_eta2
            
            interpolator%coeff_splines(nb_spline_eta1 + i ,j+1) = &
                 coeffs_1d(i+nb_spline_eta1 *(j-1) )!nb_spline_eta1 - (tmp1 -i)  + nb_spline_eta1 *(j-1) )
         end do
      end do
      
      interpolator%coeff_splines(:,1) = 0.0_8
      interpolator%coeff_splines(:,nb_spline_eta2+2) = 0.0_8
      
   case(585) ! 4. dirichlet in all sides
      interpolator%size_coeffs1=  num_cells1 + sp_deg1
      interpolator%size_coeffs2=  num_cells2 + sp_deg2
      interpolator%size_t1 = 2.0_f64*sp_deg1 + num_cells1 + 1
      interpolator%size_t2 = 2.0_f64*sp_deg2 + num_cells2 + 1
      nb_spline_eta1 = num_cells1 + sp_deg1 - 2
      nb_spline_eta2 = num_cells2 + sp_deg2 - 2
      
      ! allocation and definition of knots

      do i = 1, sp_deg1 + 1
         interpolator%t1(i) = eta1_min
      enddo
      eta1 = eta1_min
      do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
         eta1 = eta1 + delta1
         interpolator%t1(i) = eta1
      enddo
      do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
         interpolator%t1(i) = eta1
      enddo
       
      do i = 1, sp_deg2 + 1
         interpolator%t2(i) = eta2_min
      enddo
      eta2 = eta2_min
      do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
         eta2 = eta2 + delta2
         interpolator%t2(i) = eta2
      enddo
      do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
         interpolator%t2(i) = eta2
      enddo
      

      interpolator%coeff_splines(:,:) = 0.0_8
      ! allocation coefficient spline
      do i = 1,nb_spline_eta1
         do j = 1,nb_spline_eta2
            
            interpolator%coeff_splines(i+1,j+1) = &
                 coeffs_1d( i + nb_spline_eta1 *(j-1))
         end do
      end do
      
   case default
      print *, 'arbitrary_degree_spline_2d() error: set_spline_coefficients ',&
           'not recognized.'
      stop
   end select
 end subroutine set_coefficients_ad2d

  subroutine compute_interpolants_ad2d( &
    interpolator, &
    data_array, &
    eta1_coords, &
    size_eta1_coords, &
    eta2_coords, &
    size_eta2_coords )

    class(arb_deg_2d_interpolator), intent(inout)  :: interpolator
    sll_real64, dimension(:,:), intent(in)         :: data_array
    sll_real64, dimension(:), intent(in),optional           :: eta1_coords
    sll_real64, dimension(:), intent(in),optional           :: eta2_coords
    sll_int32, intent(in),optional                          :: size_eta1_coords
    sll_int32, intent(in),optional                          :: size_eta2_coords
    sll_int32  :: sz1
    sll_int32  :: sz2
    sll_real64 :: period1
    sll_real64 :: period2
    sll_int32  :: order1
    sll_int32  :: order2


    sz1 = size_eta1_coords
    sz2 = size_eta2_coords
    !PRINT *, 'SZ1 = ', SZ1, 'SZ2 = ', SZ2, 'DATA: ', SIZE(DATA_ARRAY,1), SIZE(DATA_ARRAY,2)
    SLL_ASSERT(sz1 .le. interpolator%num_pts1)
    SLL_ASSERT(sz2 .le. interpolator%num_pts2)
    SLL_ASSERT(size(data_array,1) .ge. sz1)
    SLL_ASSERT(size(data_array,2) .ge. sz2)
    SLL_ASSERT(size(eta1_coords) .ge. sz1)
    SLL_ASSERT(size(eta2_coords) .ge. sz2)

    order1  = interpolator%spline_degree1 + 1
    order2  = interpolator%spline_degree2 + 1
    period1 = interpolator%eta1_max - interpolator%eta1_min
    period2 = interpolator%eta2_max - interpolator%eta2_min


    select case (interpolator%bc_selector)
    case (0) ! periodic-periodic
       interpolator%size_coeffs1 = sz1+1
       interpolator%size_coeffs2 = sz2+1
       interpolator%size_t1 = order1 + sz1 + 1
       interpolator%size_t2 = order2 + sz2 + 1 

       call spli2d_perper( &
            period1, sz1+1, order1, eta1_coords, &
            period2, sz2+1, order2, eta2_coords, &
            data_array, interpolator%coeff_splines(1:sz1+1,1:sz2+1),&
            interpolator%t1(1:order1 + sz1 + 1), &
            interpolator%t2(1:order2 + sz2 + 1) )
       
    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       interpolator%size_coeffs1 = sz1
       interpolator%size_coeffs2 = sz2+1
       interpolator%size_t1 = order1 + sz1
       interpolator%size_t2 = order2 + sz2 + 1
       call spli2d_dirper( sz1, order1, eta1_coords, &
            period2, sz2+1, order2, eta2_coords, &
            data_array, interpolator%coeff_splines(1:sz1,1:sz2+1),&
            interpolator%t1(1:sz1+order1), &
            interpolator%t2(1:sz2+order2+1) )
  
  
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top
       interpolator%size_coeffs1 = sz1+1
       interpolator%size_coeffs2 = sz2
       interpolator%size_t1 = order1 + sz1 + 1
       interpolator%size_t2 = order2 + sz2 
       call spli2d_perdir( period1, sz1+1, order1, eta1_coords, &
            sz2, order2, eta2_coords, &
            data_array, interpolator%coeff_splines(1:sz1+1,1:sz2),&
            interpolator%t1(1:sz1+order1+1), &
            interpolator%t2(1:sz2+order2) )
       
    case (585) ! 4. dirichlet in all sides
       !print*, 'her'
       interpolator%size_coeffs1 = sz1
       interpolator%size_coeffs2 = sz2
       interpolator%size_t1 = order1 + sz1 
       interpolator%size_t2 = order2 + sz2 
 
       call spli2d_custom( sz1, order1, eta1_coords, &
            sz2, order2, eta2_coords, &
            data_array, interpolator%coeff_splines(1:sz1,1:sz2),&
            interpolator%t1(1:sz1+order1), &
            interpolator%t2(1:sz2+order2) )

    end select
  end subroutine compute_interpolants_ad2d


  function interpolate_value_ad2d( &
    interpolator, &
    eta1, &
    eta2 ) result(val)

    class(arb_deg_2d_interpolator), intent(in)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64, intent(in)         :: eta2
    sll_real64                     :: val
    sll_int32 :: size_coeffs1
    sll_int32 :: size_coeffs2
    sll_real64 :: bvalue2d
    sll_real64 :: res1,res2

    size_coeffs1 = interpolator%size_coeffs1
    size_coeffs2 = interpolator%size_coeffs2

    res1 = eta1
    res2 = eta2
    select case (interpolator%bc_selector)
    case (0) ! periodic-periodic
       if ( res1 .ge. interpolator%eta1_max ) then 
          res1 = res1 -interpolator%eta1_max
       end if
       if ( res2 .ge. interpolator%eta2_max ) then 
          res2 = res2 -interpolator%eta2_max
       end if
    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       if ( res2 .ge. interpolator%eta2_max ) then 
          res2 = res2 -interpolator%eta2_max
       end if
       SLL_ASSERT( res1 .lt. interpolator%eta1_min )
       SLL_ASSERT( res1 .gt. interpolator%eta1_max )
  
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top
       if ( res1 .ge. interpolator%eta1_max ) then 
          res1 = res1 -interpolator%eta1_max
       end if
       SLL_ASSERT( res2 .lt. interpolator%eta2_min )
       SLL_ASSERT( res2 .gt. interpolator%eta2_max )
       
    case (585) ! dirichlet-dirichlet 
       SLL_ASSERT( res1 .lt. interpolator%eta1_min )
       SLL_ASSERT( res1 .gt. interpolator%eta1_max )
       SLL_ASSERT( res2 .lt. interpolator%eta2_min )
       SLL_ASSERT( res2 .gt. interpolator%eta2_max )
    end select
          
    val = bvalue2d( &
         res1, &
         res2, &
         size_coeffs1, &
         interpolator%spline_degree1+1, &
         size_coeffs2, &
         interpolator%spline_degree2+1, &
         interpolator%coeff_splines(1:size_coeffs1,1:size_coeffs2), &
         interpolator%t1(1:interpolator%size_t1), &
         interpolator%t2(1:interpolator%size_t2))
  end function interpolate_value_ad2d


  function interpolate_derivative1_ad2d( &
    interpolator, &
    eta1, &
    eta2 ) result(val)

    class(arb_deg_2d_interpolator), intent(in)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64, intent(in)         :: eta2
    sll_real64                     :: val
    sll_int32 :: size_coeffs1
    sll_int32 :: size_coeffs2

    SLL_ASSERT( eta1 .ge. interpolator%eta1_min )
    SLL_ASSERT( eta1 .le. interpolator%eta1_max )
    SLL_ASSERT( eta2 .ge. interpolator%eta2_min )
    SLL_ASSERT( eta2 .le. interpolator%eta2_max )

    size_coeffs1 = interpolator%size_coeffs1
    size_coeffs2 = interpolator%size_coeffs2

    print *, 'interpolate_derivative1d_ad2d: not implemented'

!!$    val = bvalue( &
!!$         eta1, &
!!$         eta2, &
!!$         size_coeffs1, &
!!$         interpolator%spline_degree+1, &
!!$         size_coeffs2, &
!!$         interpolator%spline_degree+1, &
!!$         interpolator%coeffs_splines(1:size_coeffs1,1:size_coeffs2), &
!!$         interpolator%t1, &
!!$         interpolator%t2 )

  end function interpolate_derivative1_ad2d


  function interpolate_derivative2_ad2d( &
    interpolator, &
    eta1, &
    eta2 ) result(val)

    class(arb_deg_2d_interpolator), intent(in)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64, intent(in)         :: eta2
    sll_real64                     :: val
    sll_int32 :: size_coeffs1
    sll_int32 :: size_coeffs2

    SLL_ASSERT( eta1 .ge. interpolator%eta1_min )
    SLL_ASSERT( eta1 .le. interpolator%eta1_max )
    SLL_ASSERT( eta2 .ge. interpolator%eta2_min )
    SLL_ASSERT( eta2 .le. interpolator%eta2_max )

    size_coeffs1 = interpolator%size_coeffs1
    size_coeffs2 = interpolator%size_coeffs2

    print *, 'interpolate_derivative2_ad2d: not implemented'

!!$    val = bvalue( &
!!$         eta1, &
!!$         eta2, &
!!$         size_coeffs1, &
!!$         interpolator%spline_degree+1, &
!!$         size_coeffs2, &
!!$         interpolator%spline_degree+1, &
!!$         interpolator%coeffs_splines(1:size_coeffs1,1:size_coeffs2), &
!!$         interpolator%t1, &
!!$         interpolator%t2 )

  end function interpolate_derivative2_ad2d

  function interpolate_array_ad2d( &
    this, &
    num_points1, &
    num_points2, &
    data_in, &
    eta1, &
    eta2 ) result(res)

    class(arb_deg_2d_interpolator), intent(in)  :: this
    sll_real64,  dimension(:,:), intent(in)         :: eta1
    sll_real64,  dimension(:,:), intent(in)         :: eta2
    sll_real64, dimension(:,:), intent(in)         :: data_in
    sll_int32, intent(in)         :: num_points1
    sll_int32, intent(in)         :: num_points2

    sll_real64, dimension(num_points1,num_points2) :: res
 
    print *, 'interpolate_array_ad2d: not implemented'
  end function interpolate_array_ad2d

    function interpolate_2d_array_disp_ad2d( &
         this,        &
         num_points1, &
         num_points2, &
         data_in,     &
         alpha1,      &
         alpha2) result(res)
      
      class(arb_deg_2d_interpolator), intent(in)    :: this
      sll_int32, intent(in)                          :: num_points1  
      sll_int32, intent(in)                          :: num_points2 
      sll_real64, dimension(:,:), intent(in)         :: data_in
      sll_real64, dimension(:,:), intent(in)         :: alpha1
      sll_real64, dimension(:,:), intent(in)         :: alpha2  
      sll_real64, dimension(num_points1,num_points2) :: res

      print *, 'interpolate_2d_array_disp_ad2d: not implemented.'
    end function interpolate_2d_array_disp_ad2d

  
end module sll_arbitrary_degree_spline_interpolator_2d_module
