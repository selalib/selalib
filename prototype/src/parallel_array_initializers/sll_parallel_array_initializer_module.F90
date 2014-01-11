module sll_parallel_array_initializer_module
  use sll_remapper
  use sll_logical_meshes
  use sll_coordinate_transformation_2d_base_module
#include "sll_working_precision.h"
#include "sll_assert.h"
  implicit none

  ! The parallel array initializer module helps initialize multidimensional
  ! arrays that are distributed among multiple processors. Their distribution
  ! should be described by a layout object (see the documentation for the
  ! sll_remap module. The initializer assumes that the array is connected
  ! to one or more logical meshes and possibly, their corresponding 
  ! transformations. Given the various ways in which the coordinates can
  ! be grouped in terms of their associated logical meshes and transformations,
  ! there are multiple versions of the initialization subroutines. The common
  ! parameter is a user-provided function with the signature described in
  ! the following abstract type.

  abstract interface
     function sll_scalar_initializer_2d( x1, x2, params )
       use sll_working_precision
       sll_real64                                  :: sll_scalar_initializer_2d
       sll_real64, intent(in)                         :: x1
       sll_real64, intent(in)                         :: x2
       sll_real64, dimension(:), intent(in), optional :: params
     end function sll_scalar_initializer_2d
  end interface

  abstract interface
     function sll_scalar_initializer_4d( x1, x2, x3, x4, params )
       use sll_working_precision
       sll_real64                                  :: sll_scalar_initializer_4d
       sll_real64, intent(in)                         :: x1
       sll_real64, intent(in)                         :: x2
       sll_real64, intent(in)                         :: x3
       sll_real64, intent(in)                         :: x4
       sll_real64, dimension(:), intent(in), optional :: params
     end function sll_scalar_initializer_4d
  end interface

  interface sll_4d_parallel_array_initializer_cartesian
    module procedure sll_4d_parallel_array_initializer_cartesian_aux    
    module procedure sll_4d_parallel_array_initializer_cartesian_logical_1d_1d_1d_1d
    !module procedure sll_4d_parallel_array_initializer_cartesian_logical_2d_2d
    module procedure sll_4d_parallel_array_initializer_cartesian_logical_4d  
  end interface 
  
  interface sll_4d_parallel_array_initializer
    module procedure sll_2d_times_2d_parallel_array_initializer
  end interface
  
  interface  sll_2d_parallel_array_initializer_cartesian
    module procedure sll_2d_parallel_array_initializer_cartesian_logical_2d
    module procedure sll_2d_parallel_array_initializer_cartesian_array_1d_1d
  end interface 
contains

  
  subroutine sll_2d_parallel_array_initializer_cartesian_logical_2d( &
       layout, &
       mesh2d, &
       array, &
       func, &
       func_params)

    type(layout_2D), pointer                    :: layout
    type(sll_logical_mesh_2d), pointer          :: mesh2d
    sll_real64, dimension(:,:), intent(out) :: array
    procedure(sll_scalar_initializer_2d)        :: func
    sll_real64, dimension(:), optional          :: func_params


    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: loc_size_x1
    sll_int32  :: loc_size_x2
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta1
    sll_real64 :: eta2
    !sll_real64 :: x1
    !sll_real64 :: x2
    sll_int32, dimension(1:2)  :: gi ! global indices in the distributed array

    if( .not. associated(layout) ) then
       print *, '#sll_2d_parallel_array_initializer_cartesian error: ', &
            '#passed layout is uninitialized.'
    end if

    if( .not. associated(mesh2d) ) then
       print *, '#sll_2d_parallel_array_initializer_cartesian error: ', &
            '#passed mesh2d_eta1_eta2 argument is uninitialized.'
    end if


    call compute_local_sizes( layout, loc_size_x1, loc_size_x2) 

    if( size(array,1) .lt. loc_size_x1 ) then
       print *, '#sll_2d_parallel_array_initializer_cartesian error: ', &
            '#first dimension of passed array is inconsistent with ', &
            '#the size contained in the passed layout.'
    end if

    if( size(array,2) .lt. loc_size_x2 ) then
       print *, '#sll_2d_parallel_array_initializer_cartesian error: ', &
            '#second dimension of passed array is inconsistent with ', &
            '#the size contained in the passed layout.'
    end if




    eta1_min = mesh2d%eta1_min
    eta2_min = mesh2d%eta2_min
    delta1   = mesh2d%delta_eta1
    delta2   = mesh2d%delta_eta2

    ! This initializes a node-centered array. The loop should be repeated
    ! below if cell-centered or if arbitrary positions are specified.

    do j=1,loc_size_x2
      do i=1,loc_size_x1
        gi(:) = local_to_global_2D( layout, (/i,j/) )
        eta1 = eta1_min + real(gi(1)-1,f64)*delta1
        eta2 = eta2_min + real(gi(2)-1,f64)*delta2
        array(i,j) = func(eta1,eta2,func_params)
      end do
    end do

  end subroutine sll_2d_parallel_array_initializer_cartesian_logical_2d

  subroutine sll_2d_parallel_array_initializer_cartesian_array_1d_1d( &
       layout, &
       x1_array, &
       x2_array, &
       array, &
       func, &
       func_params)

    type(layout_2D), pointer                    :: layout
    sll_real64, dimension(:), intent(in) :: x1_array
    sll_real64, dimension(:), intent(in) :: x2_array
    sll_real64, dimension(:,:), intent(out) :: array
    procedure(sll_scalar_initializer_2d)        :: func
    sll_real64, dimension(:), optional          :: func_params


    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: loc_size_x1
    sll_int32  :: loc_size_x2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32, dimension(1:2)  :: gi ! global indices in the distributed array

    if( .not. associated(layout) ) then
       print *, '#sll_2d_parallel_array_initializer_cartesian error: ', &
            '#passed layout is uninitialized.'
    end if

!    if( .not. associated(mesh2d) ) then
!       print *, '#sll_2d_parallel_array_initializer_cartesian error: ', &
!            '#passed mesh2d_eta1_eta2 argument is uninitialized.'
!    end if


    call compute_local_sizes( layout, loc_size_x1, loc_size_x2) 

    if( size(array,1) .lt. loc_size_x1 ) then
       print *, '#sll_2d_parallel_array_initializer_cartesian error: ', &
            '#first dimension of passed array is inconsistent with ', &
            '#the size contained in the passed layout.'
    end if

    if( size(array,2) .lt. loc_size_x2 ) then
       print *, '#sll_2d_parallel_array_initializer_cartesian error: ', &
            '#second dimension of passed array is inconsistent with ', &
            '#the size contained in the passed layout.'
    end if



    do j=1,loc_size_x2
      do i=1,loc_size_x1
        gi(:) = local_to_global_2D( layout, (/i,j/) )
        eta1 = x1_array(gi(1))
        eta2 = x2_array(gi(2))
        array(i,j) = func(eta1,eta2,func_params)
      end do
    end do

  end subroutine sll_2d_parallel_array_initializer_cartesian_array_1d_1d



  subroutine sll_4d_parallel_array_initializer_cartesian_aux( &
       layout, &
       eta1_min, &
       eta2_min, &
       eta3_min, &
       eta4_min, &
       delta1, &
       delta2, &
       delta3, &
       delta4, &
       array, &
       func, &
       func_params) 
    type(layout_4D), pointer                    :: layout
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta3_min
    sll_real64, intent(in) :: eta4_min
    sll_real64, intent(in) :: delta1
    sll_real64, intent(in) :: delta2
    sll_real64, intent(in) :: delta3
    sll_real64, intent(in) :: delta4
    sll_real64, dimension(:,:,:,:), intent(out) :: array
    procedure(sll_scalar_initializer_4d)        :: func
    sll_real64, dimension(:), optional          :: func_params
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: l
    sll_int32  :: loc_size_x1
    sll_int32  :: loc_size_x2
    sll_int32  :: loc_size_x3
    sll_int32  :: loc_size_x4
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: eta3
    sll_real64 :: eta4
    !sll_real64 :: x1
    !sll_real64 :: x2
    !sll_real64 :: x3
    !sll_real64 :: x4
    sll_int32, dimension(1:4)  :: gi ! global indices in the distributed array

    if( .not. associated(layout) ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'passed layout is uninitialized.'
    end if

    call compute_local_sizes( layout, loc_size_x1, loc_size_x2, loc_size_x3, &
         loc_size_x4 ) 

    if( size(array,1) .lt. loc_size_x1 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'first dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    if( size(array,2) .lt. loc_size_x2 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'second dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    if( size(array,3) .lt. loc_size_x3 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'third dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    if( size(array,4) .lt. loc_size_x4 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'fourth dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    do l=1,loc_size_x4
      do k=1,loc_size_x3
        do j=1,loc_size_x2
          do i=1,loc_size_x1
            gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
            eta1 = eta1_min + real(gi(1)-1,f64)*delta1
            eta2 = eta2_min + real(gi(2)-1,f64)*delta2
            eta3 = eta3_min + real(gi(3)-1,f64)*delta3
            eta4 = eta4_min + real(gi(4)-1,f64)*delta4
            array(i,j,k,l) = func(eta1,eta2,eta3,eta4,func_params)
          end do
        end do
      end do
    end do

  end subroutine sll_4d_parallel_array_initializer_cartesian_aux

  subroutine sll_4d_parallel_array_initializer_cartesian_logical_1d_1d_1d_1d(&
    layout, &
    mesh1d_eta1, &
    mesh1d_eta2, &
    mesh1d_eta3, &
    mesh1d_eta4, &
    array, &
    func, &
    func_params) 
    type(layout_4D), pointer                    :: layout
    type(sll_logical_mesh_1d), pointer          :: mesh1d_eta1
    type(sll_logical_mesh_1d), pointer          :: mesh1d_eta2
    type(sll_logical_mesh_1d), pointer          :: mesh1d_eta3
    type(sll_logical_mesh_1d), pointer          :: mesh1d_eta4
    sll_real64, dimension(:,:,:,:), intent(out) :: array
    procedure(sll_scalar_initializer_4d)        :: func
    sll_real64, dimension(:), optional          :: func_params
    sll_real64  :: eta1_min
    sll_real64  :: eta2_min
    sll_real64  :: eta3_min
    sll_real64  :: eta4_min
    sll_real64  :: delta1
    sll_real64  :: delta2
    sll_real64  :: delta3
    sll_real64  :: delta4
    
    eta1_min = mesh1d_eta1%eta_min
    delta1 = mesh1d_eta1%delta_eta
    eta2_min = mesh1d_eta2%eta_min
    delta2 = mesh1d_eta2%delta_eta
    eta3_min = mesh1d_eta3%eta_min
    delta3 = mesh1d_eta3%delta_eta
    eta4_min = mesh1d_eta4%eta_min
    delta4 = mesh1d_eta4%delta_eta
    
    call sll_4d_parallel_array_initializer_cartesian( &
       layout, &
       eta1_min, &
       eta2_min, &
       eta3_min, &
       eta4_min, &
       delta1, &
       delta2, &
       delta3, &
       delta4, &
       array, &
       func, &
       func_params) 

    
  
  end subroutine sll_4d_parallel_array_initializer_cartesian_logical_1d_1d_1d_1d


  subroutine sll_4d_parallel_array_initializer_cartesian_logical_4d( &
       layout, &
       mesh4d, &
       array, &
       func, &
       func_params)

    type(layout_4D), pointer                    :: layout
    type(sll_logical_mesh_4d), pointer          :: mesh4d
    sll_real64, dimension(:,:,:,:), intent(out) :: array
    procedure(sll_scalar_initializer_4d)        :: func
    sll_real64, dimension(:), optional          :: func_params

    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta3_min
    sll_real64 :: eta4_min

    eta1_min = mesh4d%eta1_min
    eta2_min = mesh4d%eta2_min
    eta3_min = mesh4d%eta3_min
    eta4_min = mesh4d%eta4_min
    delta1   = mesh4d%delta_eta1
    delta2   = mesh4d%delta_eta2
    delta3   = mesh4d%delta_eta3
    delta4   = mesh4d%delta_eta4

    call sll_4d_parallel_array_initializer_cartesian( &
       layout, &
       eta1_min, &
       eta2_min, &
       eta3_min, &
       eta4_min, &
       delta1, &
       delta2, &
       delta3, &
       delta4, &
       array, &
       func, &
       func_params) 



  end subroutine sll_4d_parallel_array_initializer_cartesian_logical_4d





  ! it is convenient for example, to separate a 4D mesh into two parts, one can remain 
  ! cartesian while the other may be transformed.

  subroutine sll_2d_times_2d_parallel_array_initializer( &
       layout, &
       mesh2d_eta1_eta2, &
       mesh2d_eta3_eta4, &
       array, &
       func, &
       func_params, &
       transf_x1_x2, &
       transf_x3_x4 )

    type(layout_4D), pointer                    :: layout
    type(sll_logical_mesh_2d), pointer          :: mesh2d_eta1_eta2
    type(sll_logical_mesh_2d), pointer          :: mesh2d_eta3_eta4
    sll_real64, dimension(:,:,:,:), intent(out) :: array
    procedure(sll_scalar_initializer_4d)        :: func
    sll_real64, dimension(:), optional          :: func_params

    class(sll_coordinate_transformation_2d_base), pointer, optional :: &
         transf_x1_x2

    class(sll_coordinate_transformation_2d_base), pointer, optional :: &
         transf_x3_x4

    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: l
    sll_int32  :: loc_size_x1
    sll_int32  :: loc_size_x2
    sll_int32  :: loc_size_x3
    sll_int32  :: loc_size_x4
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta3_min
    sll_real64 :: eta4_min
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: eta3
    sll_real64 :: eta4
    sll_real64 :: x1
    sll_real64 :: x2
    sll_real64 :: x3
    sll_real64 :: x4
    sll_int32  :: case_selector
    sll_int32, dimension(1:4)  :: gi ! global indices in the distributed array

    if( .not. associated(layout) ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'passed layout is uninitialized.'
    end if

    if( .not. associated(mesh2d_eta1_eta2) ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'passed mesh2d_eta1_eta2 argument is uninitialized.'
    end if

    if( .not. associated(mesh2d_eta3_eta4) ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'passed mesh2d_eta3_eta4 argument is uninitialized.'
    end if

    call compute_local_sizes( layout, loc_size_x1, loc_size_x2, loc_size_x3, &
         loc_size_x4 ) 

    if( size(array,1) .lt. loc_size_x1 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'first dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    if( size(array,2) .lt. loc_size_x2 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'second dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    if( size(array,3) .lt. loc_size_x3 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'third dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    if( size(array,4) .lt. loc_size_x4 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'fourth dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    if(present( transf_x1_x2 ) ) then
       if( .not. associated(transf_x1_x2%mesh, mesh2d_eta1_eta2) ) then 
          print *, 'sll_4d_parallel_array_initializer warning: ', &
               'the mesh associated to the transf_x1_x2 transformation ', &
               'is not the same as the mesh2d_eta1_eta2 logical mesh. ', &
               'Unless the parameters of these meshes are the same, ', &
               'bad things will happen.'
       end if
    end if

    if(present( transf_x3_x4 ) ) then
       if( .not. associated(transf_x3_x4%mesh, mesh2d_eta3_eta4) ) then 
          print *, 'sll_4d_parallel_array_initializer warning: ', &
               'the mesh associated to the transf_x3_x4 transformation ', &
               'is not the same as the mesh2d_eta3_eta4 logical mesh. ', &
               'Unless the parameters of these meshes are the same, ', &
               'bad things will happen.'
       end if
    end if

    case_selector = 0

    if( present(transf_x1_x2) ) then
       case_selector = case_selector + 1
    end if

    if( present(transf_x3_x4) ) then
       case_selector = case_selector + 2
    end if
    eta1_min = mesh2d_eta1_eta2%eta1_min
    eta2_min = mesh2d_eta1_eta2%eta2_min
    eta3_min = mesh2d_eta3_eta4%eta1_min
    eta4_min = mesh2d_eta3_eta4%eta2_min
    delta1   = mesh2d_eta1_eta2%delta_eta1
    delta2   = mesh2d_eta1_eta2%delta_eta2
    delta3   = mesh2d_eta3_eta4%delta_eta1
    delta4   = mesh2d_eta3_eta4%delta_eta2

    ! This initializes a node-centered array. The loop should be repeated
    ! below if cell-centered or if arbitrary positions are specified.

    select case (case_selector)

    case (0) ! none of the transformations was provided
       do l=1,loc_size_x4
          do k=1,loc_size_x3
             do j=1,loc_size_x2
                do i=1,loc_size_x1
                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
                   array(i,j,k,l) = func(eta1,eta2,eta3,eta4,func_params)
                end do
             end do
          end do
       end do
    case(1) ! Only the x1,x2 transfomation was provided.
       do l=1,loc_size_x4
          do k=1,loc_size_x3
             do j=1,loc_size_x2
                do i=1,loc_size_x1
                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
                   x1 = transf_x1_x2%x1(eta1,eta2)
                   x2 = transf_x1_x2%x2(eta1,eta2)
                   array(i,j,k,l) = func(x1,x2,eta3,eta4,func_params)
                end do
             end do
          end do
       end do
    case(2) ! Only the x3,x4 transformation was provided
       do l=1,loc_size_x4
          do k=1,loc_size_x3
             do j=1,loc_size_x2
                do i=1,loc_size_x1
                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
                   x3 = transf_x3_x4%x1(eta3,eta4)
                   x4 = transf_x3_x4%x2(eta3,eta4)
                   array(i,j,k,l) = func(eta1,eta2,x3,x4,func_params)
                end do
             end do
          end do
       end do
    case(3)  ! Both transformations were provided.
       do l=1,loc_size_x4
          do k=1,loc_size_x3
             do j=1,loc_size_x2
                do i=1,loc_size_x1
                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
                   x1 = transf_x1_x2%x1(eta1,eta2)
                   x2 = transf_x1_x2%x2(eta1,eta2)
                   x3 = transf_x3_x4%x1(eta3,eta4)
                   x4 = transf_x3_x4%x2(eta3,eta4)
                   array(i,j,k,l) = func(x1,x2,x3,x4,func_params)
                end do
             end do
          end do
       end do
    end select

  end subroutine sll_2d_times_2d_parallel_array_initializer

!  subroutine sll_4d_parallel_array_initializer( &
!       layout, &
!       mesh2d_eta1_eta2, &
!       mesh2d_eta3_eta4, &
!       array, &
!       func, &
!       func_params, &
!       transf_x1_x2, &
!       transf_x3_x4 )
!
!    type(layout_4D), pointer                    :: layout
!    type(sll_logical_mesh_2d), pointer          :: mesh2d_eta1_eta2
!    type(sll_logical_mesh_2d), pointer          :: mesh2d_eta3_eta4
!    sll_real64, dimension(:,:,:,:), intent(out) :: array
!    procedure(sll_scalar_initializer_4d)        :: func
!    sll_real64, dimension(:), optional          :: func_params
!
!    class(sll_coordinate_transformation_2d_base), pointer, optional :: &
!         transf_x1_x2
!
!    class(sll_coordinate_transformation_2d_base), pointer, optional :: &
!         transf_x3_x4
!
!    sll_int32  :: i
!    sll_int32  :: j
!    sll_int32  :: k
!    sll_int32  :: l
!    sll_int32  :: loc_size_x1
!    sll_int32  :: loc_size_x2
!    sll_int32  :: loc_size_x3
!    sll_int32  :: loc_size_x4
!    sll_real64 :: delta1
!    sll_real64 :: delta2
!    sll_real64 :: delta3
!    sll_real64 :: delta4
!    sll_real64 :: eta1_min
!    sll_real64 :: eta2_min
!    sll_real64 :: eta3_min
!    sll_real64 :: eta4_min
!    sll_real64 :: eta1
!    sll_real64 :: eta2
!    sll_real64 :: eta3
!    sll_real64 :: eta4
!    sll_real64 :: x1
!    sll_real64 :: x2
!    sll_real64 :: x3
!    sll_real64 :: x4
!    sll_int32  :: case_selector
!    sll_int32, dimension(1:4)  :: gi ! global indices in the distributed array
!
!    if( .not. associated(layout) ) then
!       print *, 'sll_4d_parallel_array_initializer error: ', &
!            'passed layout is uninitialized.'
!    end if
!
!    if( .not. associated(mesh2d_eta1_eta2) ) then
!       print *, 'sll_4d_parallelsll_parallel_array_initializer_module.F90_array_initializer error: ', &
!            'passed mesh2d_eta1_eta2 argument is uninitialized.'
!    end if
!
!    if( .not. associated(mesh2d_eta3_eta4) ) then
!       print *, 'sll_4d_parallel_array_initializer error: ', &
!            'passed mesh2d_eta3_eta4 argument is uninitialized.'
!    end if
!
!    call compute_local_sizes( layout, loc_size_x1, loc_size_x2, loc_size_x3, &
!         loc_size_x4 ) 
!
!    if( size(array,1) .lt. loc_size_x1 ) then
!       print *, 'sll_4d_parallel_array_initializer error: ', &
!            'first dimension of passed array is inconsistent with ', &
!            'the size contained in the passed layout.'
!    end if
!
!    if( size(array,2) .lt. loc_size_x2 ) then
!       print *, 'sll_4d_parallel_array_initializer error: ', &
!            'second dimension of passed array is inconsistent with ', &
!            'the size contained in the passed layout.'
!    end if
!
!    if( size(array,3) .lt. loc_size_x3 ) then
!       print *, 'sll_4d_parallel_array_initializer error: ', &
!            'third dimension of passed array is inconsistent with ', &
!            'the size contained in the passed layout.'
!    end if
!
!    if( size(array,4) .lt. loc_size_x4 ) then
!       print *, 'sll_4d_parallel_array_initializer error: ', &
!            'fourth dimension of passed array is inconsistent with ', &
!            'the size contained in the passed layout.'
!    end if
!
!    if(present( transf_x1_x2 ) ) then
!       if( .not. associated(transf_x1_x2%mesh, mesh2d_eta1_eta2) ) then 
!          print *, 'sll_4d_parallel_array_initializer warning: ', &
!               'the mesh associated to the transf_x1_x2 transformation ', &
!               'is not the same as the mesh2d_eta1_eta2 logical mesh. ', &
!               'Unless the parameters of these meshes are the same, ', &
!               'bad things will happen.'
!       end if
!    end if
!
!    if(present( transf_x3_x4 ) ) then
!       if( .not. associated(transf_x3_x4%mesh, mesh2d_eta3_eta4) ) then 
!          print *, 'sll_4d_parallel_array_initializer warning: ', &
!               'the mesh associated to the transf_x3_x4 transformation ', &
!               'is not the same as the mesh2d_eta3_eta4 logical mesh. ', &
!               'Unless the parameters of these meshes are the same, ', &
!               'bad things will happen.'
!       end if
!    end if
!
!    case_selector = 0
!
!    if( present(transf_x1_x2) ) then
!       case_selector = case_selector + 1
!    end if
!
!    if( present(transf_x3_x4) ) then
!       case_selector = case_selector + 2
!    end if
!    eta1_min = mesh2d_eta1_eta2%eta1_min
!    eta2_min = mesh2d_eta1_eta2%eta2_min
!    eta3_min = mesh2d_eta3_eta4%eta1_min
!    eta4_min = mesh2d_eta3_eta4%eta2_min
!    delta1   = mesh2d_eta1_eta2%delta_eta1
!    delta2   = mesh2d_eta1_eta2%delta_eta2
!    delta3   = mesh2d_eta3_eta4%delta_eta1
!    delta4   = mesh2d_eta3_eta4%delta_eta2
!
!    ! This initializes a node-centered array. The loop should be repeated
!    ! below if cell-centered or if arbitrary positions are specified.
!
!    select case (case_selector)
!
!    case (0) ! none of the transformations was provided
!       do l=1,loc_size_x4
!          do k=1,loc_size_x3
!             do j=1,loc_size_x2
!                do i=1,loc_size_x1
!                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
!                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
!                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
!                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
!                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
!                   array(i,j,k,l) = func(eta1,eta2,eta3,eta4,func_params)
!                end do
!             end do
!          end do
!       end do
!    case(1) ! Only the x1,x2 transfomation was provided.
!       do l=1,loc_size_x4
!          do k=1,loc_size_x3
!             do j=1,loc_size_x2
!                do i=1,loc_size_x1
!                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
!                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
!                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
!                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
!                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
!                   x1 = transf_x1_x2%x1(eta1,eta2)
!                   x2 = transf_x1_x2%x2(eta1,eta2)
!                   array(i,j,k,l) = func(x1,x2,eta3,eta4,func_params)
!                end do
!             end do
!          end do
!       end do
!    case(2) ! Only the x3,x4 transformation was provided
!       do l=1,loc_size_x4
!          do k=1,loc_size_x3
!             do j=1,loc_size_x2
!                do i=1,loc_size_x1
!                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
!                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
!                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
!                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
!                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
!                   x3 = transf_x3_x4%x1(eta3,eta4)
!                   x4 = transf_x3_x4%x2(eta3,eta4)
!                   array(i,j,k,l) = func(eta1,eta2,x3,x4,func_params)
!                end do
!             end do
!          end do
!       end do
!    case(3)  ! Both transformations were provided.
!       do l=1,loc_size_x4
!          do k=1,loc_size_x3
!             do j=1,loc_size_x2
!                do i=1,loc_size_x1
!                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
!                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
!                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
!                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
!                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
!                   x1 = transf_x1_x2%x1(eta1,eta2)
!                   x2 = transf_x1_x2%x2(eta1,eta2)
!                   x3 = transf_x3_x4%x1(eta3,eta4)
!                   x4 = transf_x3_x4%x2(eta3,eta4)
!                   array(i,j,k,l) = func(x1,x2,x3,x4,func_params)
!                end do
!             end do
!          end do
!       end do
!    end select
!
!  end subroutine sll_4d_parallel_array_initializer

!in case of interpolation in the vx, vy, x, y direction
 subroutine sll_4d_parallel_array_initializer_finite_volume( &
       layout, &
       mesh2d_eta1_eta2, &
       mesh2d_eta3_eta4, &
       array, &
       func, &
       func_params, &
       transf_x1_x2, &
       transf_x3_x4, &
! in case of a mesh with cell subdivisions
       subcells1, &
       subcells2, &
       subcells3, &
       subcells4)

    type(layout_4D), pointer                    :: layout
    type(sll_logical_mesh_2d), pointer          :: mesh2d_eta1_eta2
    type(sll_logical_mesh_2d), pointer          :: mesh2d_eta3_eta4
    sll_real64, dimension(:,:,:,:), intent(out) :: array
    procedure(sll_scalar_initializer_4d)        :: func
    sll_real64, dimension(:), optional          :: func_params
    sll_int32,optional   :: subcells1
    sll_int32,optional   :: subcells2
    sll_int32,optional   :: subcells3
    sll_int32,optional   :: subcells4


    class(sll_coordinate_transformation_2d_base), pointer, optional :: &
         transf_x1_x2

    class(sll_coordinate_transformation_2d_base), pointer, optional :: &
         transf_x3_x4

    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: l
    sll_int32  :: loc_size_x1
    sll_int32  :: loc_size_x2
    sll_int32  :: loc_size_x3
    sll_int32  :: loc_size_x4
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta3_min
    sll_real64 :: eta4_min
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: eta3
    sll_real64 :: eta4
    sll_real64 :: x1
    sll_real64 :: x2
    sll_real64 :: x3
    sll_real64 :: x4
    sll_int32  :: case_selector
    sll_int32    :: sub1,sub2,sub3,sub4
    sll_int32, dimension(1:4)  :: gi ! global indices in the distributed array

    if( .not. associated(layout) ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'passed layout is uninitialized.'
    end if

    if( .not. associated(mesh2d_eta1_eta2) ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'passed mesh2d_eta1_eta2 argument is uninitialized.'
    end if

    if( .not. associated(mesh2d_eta3_eta4) ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'passed mesh2d_eta3_eta4 argument is uninitialized.'
    end if

    call compute_local_sizes( layout, loc_size_x1, loc_size_x2, loc_size_x3, &
         loc_size_x4 ) 

    if( size(array,1) .lt. loc_size_x1 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'first dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    if( size(array,2) .lt. loc_size_x2 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'second dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    if( size(array,3) .lt. loc_size_x3 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'third dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    if( size(array,4) .lt. loc_size_x4 ) then
       print *, 'sll_4d_parallel_array_initializer error: ', &
            'fourth dimension of passed array is inconsistent with ', &
            'the size contained in the passed layout.'
    end if

    if(present( transf_x1_x2 ) ) then
       if( .not. associated(transf_x1_x2%mesh, mesh2d_eta1_eta2) ) then 
          print *, 'sll_4d_parallel_array_initializer warning: ', &
               'the mesh associated to the transf_x1_x2 transformation ', &
               'is not the same as the mesh2d_eta1_eta2 logical mesh. ', &
               'Unless the parameters of these meshes are the same, ', &
               'bad things will happen.'
       end if
    end if

    if(present( transf_x3_x4 ) ) then
       if( .not. associated(transf_x3_x4%mesh, mesh2d_eta3_eta4) ) then 
          print *, 'sll_4d_parallel_array_initializer warning: ', &
               'the mesh associated to the transf_x3_x4 transformation ', &
               'is not the same as the mesh2d_eta3_eta4 logical mesh. ', &
               'Unless the parameters of these meshes are the same, ', &
               'bad things will happen.'
       end if
    end if


   if(.not.present(subcells1  ) ) then
     sub1=1
  else
     sub1=subcells1
  end if

   if(.not.present(subcells2  ) ) then
     sub2=1
  else
     sub2=subcells2
  end if

   if(.not.present(subcells3  ) ) then
     sub3=1
  else
     sub3=subcells3
  end if
    
   if(.not.present(subcells3  ) ) then
     sub4=1
  else
     sub4=subcells4
  end if


    case_selector = 0

    if( present(transf_x1_x2) ) then
       case_selector = case_selector + 1
    end if

    if( present(transf_x3_x4) ) then
       case_selector = case_selector + 2
    end if
    eta1_min = mesh2d_eta1_eta2%eta1_min
    eta2_min = mesh2d_eta1_eta2%eta2_min
    eta3_min = mesh2d_eta3_eta4%eta1_min
    eta4_min = mesh2d_eta3_eta4%eta2_min
    delta1   = mesh2d_eta1_eta2%delta_eta1/sub1
    delta2   = mesh2d_eta1_eta2%delta_eta2/sub2
    delta3   = mesh2d_eta3_eta4%delta_eta1/sub3
    delta4   = mesh2d_eta3_eta4%delta_eta2/sub4

    ! This initializes a node-centered array. The loop should be repeated
    ! below if cell-centered or if arbitrary positions are specified.

    select case (case_selector)

    case (0) ! none of the transformations was provided
       do l=1,loc_size_x4
          do k=1,loc_size_x3
             do j=1,loc_size_x2
                do i=1,loc_size_x1
                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
                   array(i,j,k,l) = func(eta1,eta2,eta3,eta4,func_params)
                end do
             end do
          end do
       end do
    case(1) ! Only the x1,x2 transfomation was provided.
       do l=1,loc_size_x4
          do k=1,loc_size_x3
             do j=1,loc_size_x2
                do i=1,loc_size_x1
                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
                   x1 = transf_x1_x2%x1(eta1,eta2)
                   x2 = transf_x1_x2%x2(eta1,eta2)
                   array(i,j,k,l) = func(x1,x2,eta3,eta4,func_params)
                end do
             end do
          end do
       end do
    case(2) ! Only the x3,x4 transformation was provided
       do l=1,loc_size_x4
          do k=1,loc_size_x3
             do j=1,loc_size_x2
                do i=1,loc_size_x1
                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
                   x3 = transf_x3_x4%x1(eta3,eta4)
                   x4 = transf_x3_x4%x2(eta3,eta4)
                   array(i,j,k,l) = func(eta1,eta2,x3,x4,func_params)
                end do
             end do
          end do
       end do
    case(3)  ! Both transformations were provided.
       do l=1,loc_size_x4
          do k=1,loc_size_x3
             do j=1,loc_size_x2
                do i=1,loc_size_x1
                   gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
                   eta1 = eta1_min + real(gi(1)-1,f64)*delta1
                   eta2 = eta2_min + real(gi(2)-1,f64)*delta2
                   eta3 = eta3_min + real(gi(3)-1,f64)*delta3
                   eta4 = eta4_min + real(gi(4)-1,f64)*delta4
                   x1 = transf_x1_x2%x1(eta1,eta2)
                   x2 = transf_x1_x2%x2(eta1,eta2)
                   x3 = transf_x3_x4%x1(eta3,eta4)
                   x4 = transf_x3_x4%x2(eta3,eta4)
                   array(i,j,k,l) = func(x1,x2,x3,x4,func_params)
                end do
             end do
          end do
       end do
    end select

  end subroutine sll_4d_parallel_array_initializer_finite_volume




!  subroutine sll_4d_parallel_array_initializer_cartesian( &
!       layout, &
!       mesh4d, &
!       array, &
!       func, &
!       func_params)
!
!    type(layout_4D), pointer                    :: layout
!    type(sll_logical_mesh_4d), pointer          :: mesh4d
!    sll_real64, dimension(:,:,:,:), intent(out) :: array
!    procedure(sll_scalar_initializer_4d)        :: func
!    sll_real64, dimension(:), optional          :: func_params
!
!
!    sll_int32  :: i
!    sll_int32  :: j
!    sll_int32  :: k
!    sll_int32  :: l
!    sll_int32  :: loc_size_x1
!    sll_int32  :: loc_size_x2
!    sll_int32  :: loc_size_x3
!    sll_int32  :: loc_size_x4
!    sll_real64 :: delta1
!    sll_real64 :: delta2
!    sll_real64 :: delta3
!    sll_real64 :: delta4
!    sll_real64 :: eta1_min
!    sll_real64 :: eta2_min
!    sll_real64 :: eta3_min
!    sll_real64 :: eta4_min
!    sll_real64 :: eta1
!    sll_real64 :: eta2
!    sll_real64 :: eta3
!    sll_real64 :: eta4
!    sll_int32, dimension(1:4)  :: gi ! global indices in the distributed array
!
!    if( .not. associated(layout) ) then
!       print *, 'sll_4d_parallel_array_initializer_cartesian error: ', &
!            'passed layout is uninitialized.'
!    end if
!
!    if( .not. associated(mesh4d) ) then
!       print *, 'sll_4d_parallel_array_initializer_cartesian error: ', &
!            'passed mesh2d_eta1_eta2 argument is uninitialized.'
!    end if
!
!
!    call compute_local_sizes( layout, loc_size_x1, loc_size_x2, loc_size_x3, &
!         loc_size_x4 ) 
!
!    if( size(array,1) .lt. loc_size_x1 ) then
!       print *, 'sll_4d_parallel_array_initializer_cartesian error: ', &
!            'first dimension of passed array is inconsistent with ', &
!            'the size contained in the passed layout.'
!    end if
!
!    if( size(array,2) .lt. loc_size_x2 ) then
!       print *, 'sll_4d_parallel_array_initializer_cartesian error: ', &
!            'second dimension of passed array is inconsistent with ', &
!            'the size contained in the passed layout.'
!    end if
!    if( size(array,3) .lt. loc_size_x3 ) then
!       print *, 'sll_4d_parallel_array_initializer_cartesian error: ', &
!            'third dimension of passed array is inconsistent with ', &
!            'the size contained in the passed layout.'
!    end if
!
!    if( size(array,4) .lt. loc_size_x4 ) then
!       print *, 'sll_4d_parallel_array_initializer_cartesian error: ', &
!            'fourth dimension of passed array is inconsistent with ', &
!            'the size contained in the passed layout.'
!    end if
!
!
!
!    eta1_min = mesh4d%eta1_min
!    eta2_min = mesh4d%eta2_min
!    eta3_min = mesh4d%eta3_min
!    eta4_min = mesh4d%eta4_min
!    delta1   = mesh4d%delta_eta1
!    delta2   = mesh4d%delta_eta2
!    delta3   = mesh4d%delta_eta3
!    delta4   = mesh4d%delta_eta4
!
!    ! This initializes a node-centered array. The loop should be repeated
!    ! below if cell-centered or if arbitrary positions are specified.
!
!    do l=1,loc_size_x4
!       do k=1,loc_size_x3
!          do j=1,loc_size_x2
!             do i=1,loc_size_x1
!                gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
!                eta1 = eta1_min + real(gi(1)-1,f64)*delta1
!                eta2 = eta2_min + real(gi(2)-1,f64)*delta2
!                eta3 = eta3_min + real(gi(3)-1,f64)*delta3
!                eta4 = eta4_min + real(gi(4)-1,f64)*delta4
!                array(i,j,k,l) = func(eta1,eta2,eta3,eta4,func_params)
!             end do
!          end do
!       end do
!    end do
!
!  end subroutine sll_4d_parallel_array_initializer_cartesian



end module sll_parallel_array_initializer_module
