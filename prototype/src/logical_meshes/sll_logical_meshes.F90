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

!> @file sll_logical_meshes.F90
!> @namespace sll_logical_meshes
!> @brief logical mesh basic types

module sll_logical_meshes
#include "sll_working_precision.h"
#include "sll_memory.h"
  implicit none

  !> @brief 1D logical mesh
  type sll_logical_mesh_1d
     sll_int32  :: num_cells
     sll_real64 :: eta_min
     sll_real64 :: eta_max
     sll_real64 :: delta_eta
  end type sll_logical_mesh_1d


  !> @brief 2D logical mesh
  type sll_logical_mesh_2d
     sll_int32  :: num_cells1 !< number of cells in direction 1
     sll_int32  :: num_cells2 !< number of cells in direction 2
     sll_real64 :: eta1_min   !< minimum value of eta, direction 1
     sll_real64 :: eta1_max   !< maximum value of eta, direction 1
     sll_real64 :: eta2_min   !< minimum value of eta, direction 2
     sll_real64 :: eta2_max   !< maximum value of eta, direction 2
     sll_real64 :: delta_eta1 !< cell spacing, direction 1
     sll_real64 :: delta_eta2 !< cell spacing, direction 2
  end type sll_logical_mesh_2d

  type sll_logical_mesh_2d_ptr
     type(sll_logical_mesh_2d), pointer :: lm
  end type sll_logical_mesh_2d_ptr

  !> @brief 3D logical mesh
  type sll_logical_mesh_3d
     sll_int32  :: num_cells1 !< number of cells in direction 1
     sll_int32  :: num_cells2 !< number of cells in direction 2
     sll_int32  :: num_cells3 !< number of cells in direction 3 
     sll_real64 :: eta1_min   !< minimum value of eta, direction 1
     sll_real64 :: eta1_max   !< maximum value of eta, direction 1
     sll_real64 :: eta2_min   !< minimum value of eta, direction 2
     sll_real64 :: eta2_max   !< maximum value of eta, direction 2
     sll_real64 :: eta3_min   !< minimum value of eta, direction 3
     sll_real64 :: eta3_max   !< maximum value of eta, direction 3
     sll_real64 :: delta_eta1 !< cell spacing, direction 1
     sll_real64 :: delta_eta2 !< cell spacing, direction 2
     sll_real64 :: delta_eta3 !< cell spacing, direction 3
  end type sll_logical_mesh_3d

  !> @brief 4D logical mesh
  type sll_logical_mesh_4d
     sll_int32  :: num_cells1 !< number of cells in direction 1
     sll_int32  :: num_cells2 !< number of cells in direction 2
     sll_int32  :: num_cells3 !< number of cells in direction 3
     sll_int32  :: num_cells4 !< number of cells in direction 4
     sll_real64 :: eta1_min   !< minimum value of eta, direction 1
     sll_real64 :: eta1_max   !< maximum value of eta, direction 1
     sll_real64 :: eta2_min   !< minimum value of eta, direction 2
     sll_real64 :: eta2_max   !< maximum value of eta, direction 2
     sll_real64 :: eta3_min   !< minimum value of eta, direction 3
     sll_real64 :: eta3_max   !< maximum value of eta, direction 3
     sll_real64 :: eta4_min   !< minimum value of eta, direction 4
     sll_real64 :: eta4_max   !< maximum value of eta, direction 4
     sll_real64 :: delta_eta1 !< cell spacing, direction 1
     sll_real64 :: delta_eta2 !< cell spacing, direction 2
     sll_real64 :: delta_eta3 !< cell spacing, direction 3
     sll_real64 :: delta_eta4 !< cell spacing, direction 4
  end type sll_logical_mesh_4d

  ! this should be sll_delete library-wide...
  interface delete
     module procedure delete_logical_mesh_1d
     module procedure delete_logical_mesh_2d
     module procedure delete_logical_mesh_3d
     module procedure delete_logical_mesh_4d
  end interface delete

  interface operator(*)
     module procedure tensor_product_1d_1d
     module procedure tensor_product_2d_2d
  end interface operator(*)

  interface sll_display
     module procedure display_logical_mesh_1d
     module procedure display_logical_mesh_2d
     module procedure display_logical_mesh_3d
     module procedure display_logical_mesh_4d
  end interface sll_display
contains

#define TEST_PRESENCE_AND_ASSIGN_VAL( obj, arg, slot, default_val ) \
  if( present(arg) ) then ; \
    obj%slot = arg; \
  else; \
    obj%slot = default_val; \
end if

  !> @brief allocates the memory space for a new 1D logical mesh on the heap,
  !> initializes it with the given arguments and returns a pointer to the
  !> object.
  !> @param num_cells1 integer denoting the number of cells.
  !> @param eta1_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh.
  !> @param eta1_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh.
  !> return a pointer to the newly allocated object.
  function new_logical_mesh_1d( &
    num_cells, &
    eta_min, &
    eta_max ) result(m)

    type(sll_logical_mesh_1d), pointer :: m
    sll_int32, intent(in)  :: num_cells
    sll_real64, optional, intent(in) :: eta_min
    sll_real64, optional, intent(in) :: eta_max
    !sll_real64 :: delta
    sll_int32 :: ierr
    SLL_ALLOCATE(m, ierr)
    call initialize_logical_mesh_1d( m, num_cells, eta_min, eta_max )
  end function new_logical_mesh_1d


  !> @brief initializes a previously allocated 1D logical mesh object.
  !> @param num_cells1 integer denoting the number of cells.
  !> @param eta1_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh.
  !> @param eta1_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh.
  !> return a pointer to the newly allocated object.
  subroutine initialize_logical_mesh_1d( m, num_cells, eta_min, eta_max )
    type(sll_logical_mesh_1d), pointer :: m
    sll_int32, intent(in)  :: num_cells
    sll_real64, optional, intent(in) :: eta_min
    sll_real64, optional, intent(in) :: eta_max

    m%num_cells = num_cells

    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta_min, eta_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta_max, eta_max, 1.0_f64 )

    m%delta_eta   = (m%eta_max - m%eta_min)/real(num_cells,f64)

    if ( m%eta_max <= m%eta_min) then
       print*,'ERROR, initialize_logical_mesh_1d(): ', &
            'Problem to construct the mesh 1d '
       print*,'because eta_max <= eta_min'
    end if
  end subroutine initialize_logical_mesh_1d


  function tensor_product_1d_1d( m_a, m_b) result(m_c)
    type(sll_logical_mesh_1d), intent(in),  pointer :: m_a
    type(sll_logical_mesh_1d), intent(in),  pointer :: m_b
    type(sll_logical_mesh_2d),              pointer :: m_c

    m_c => new_logical_mesh_2d( &
    m_a%num_cells, &
    m_b%num_cells, &
    m_a%eta_min, &
    m_a%eta_max, &
    m_b%eta_min, &
    m_b%eta_max ) 

  end function tensor_product_1d_1d

  function tensor_product_2d_2d( m_a, m_b) result(m_c)
    type(sll_logical_mesh_2d), intent(in),  pointer :: m_a
    type(sll_logical_mesh_2d), intent(in),  pointer :: m_b
    type(sll_logical_mesh_4d),              pointer :: m_c

    m_c => new_logical_mesh_4d( &
    m_a%num_cells1, &
    m_a%num_cells2, &
    m_b%num_cells1, &
    m_b%num_cells2, &
    m_a%eta1_min,   &
    m_a%eta1_max,   &
    m_a%eta2_min,   &
    m_a%eta2_max,   &
    m_b%eta1_min,   &
    m_b%eta1_max,   &  
    m_b%eta2_min,   &
    m_b%eta2_max ) 

  end function tensor_product_2d_2d

  subroutine initialize_eta1_node_1d( m, eta1_node )
    type(sll_logical_mesh_1d), pointer :: m
    sll_real64, dimension(:), pointer :: eta1_node
    sll_int32  :: num_cells
    sll_real64 :: eta_min
    sll_real64 :: delta_eta
    sll_int32 :: i
    sll_int32 :: ierr
    
    num_cells = m%num_cells
    eta_min = m%eta_min
    delta_eta = m%delta_eta
    SLL_ALLOCATE(eta1_node(num_cells+1), ierr)
    do i=1,num_cells+1
      eta1_node(i) = eta_min+real(i-1,f64)*delta_eta
    enddo    
    
    
  end subroutine initialize_eta1_node_1d

  !> @brief allocates the memory space for a new 2D logical mesh on the heap,
  !> initializes it with the given arguments and returns a pointer to the
  !> object.
  !> @param num_cells1 integer denoting the number of cells, direction 1.
  !> @param num_cells2 integer denoting the number of cells, direction 2.
  !> @param eta1_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh, direction 1.
  !> @param eta1_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh, direction 1.
  !> @param eta2_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh, direction 2.
  !> @param eta2_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh, direction 2.
  !> return a pointer to the newly allocated object.
  function new_logical_mesh_2d( &
    num_cells1, &
    num_cells2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max ) result(m)

    type(sll_logical_mesh_2d), pointer :: m
    sll_int32, intent(in)  :: num_cells1
    sll_int32, intent(in)  :: num_cells2
    sll_real64, optional, intent(in) :: eta1_min
    sll_real64, optional, intent(in) :: eta1_max
    sll_real64, optional, intent(in) :: eta2_min
    sll_real64, optional, intent(in) :: eta2_max
    sll_int32 :: ierr

    SLL_ALLOCATE(m, ierr)
    call initialize_logical_mesh_2d( &
         m, &
         num_cells1, &
         num_cells2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max )

  end function new_logical_mesh_2d

  !> @brief initializes a logical mesh 2D object that has been already 
  !> allocated.
  !> @param num_cells1 integer denoting the number of cells, direction 1.
  !> @param num_cells2 integer denoting the number of cells, direction 2.
  !> @param eta1_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh, direction 1.
  !> @param eta1_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh, direction 1.
  !> @param eta2_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh, direction 2.
  !> @param eta2_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh, direction 2.
  !> return a pointer to the newly allocated object.
  subroutine initialize_logical_mesh_2d( &
    m, & 
    num_cells1, &
    num_cells2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max )

    type(sll_logical_mesh_2d), pointer :: m
    sll_int32, intent(in)  :: num_cells1
    sll_int32, intent(in)  :: num_cells2
    sll_real64, optional, intent(in) :: eta1_min
    sll_real64, optional, intent(in) :: eta1_max
    sll_real64, optional, intent(in) :: eta2_min
    sll_real64, optional, intent(in) :: eta2_max

    m%num_cells1 = num_cells1
    m%num_cells2 = num_cells2

    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta1_min, eta1_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta1_max, eta1_max, 1.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta2_min, eta2_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta2_max, eta2_max, 1.0_f64 )

    m%delta_eta1   = (m%eta1_max - m%eta1_min)/real(num_cells1,f64)
    m%delta_eta2   = (m%eta2_max - m%eta2_min)/real(num_cells2,f64)

    if ( m%eta1_max <= m%eta1_min) then
       print*,'ERROR, initialize_logical_mesh_2d(): ', &
            'Problem to construct the mesh 2d '
       print*,'because eta1_max <= eta1_min'
    end if
    if ( m%eta2_max <= m%eta2_min) then
       print*,'ERROR, initialize_logical_mesh_2d(): ', &
            'Problem to construct the mesh 2d '
       print*,'because eta2_max <= eta2_min'
    end if
  end subroutine initialize_logical_mesh_2d

  !> @brief allocates the memory space for a new 3D logical mesh on the heap,
  !> initializes it with the given arguments and returns a pointer to the
  !> object.
  !> @param num_cells1 integer denoting the number of cells, direction 1.
  !> @param num_cells2 integer denoting the number of cells, direction 2.
  !> @param num_cells3 integer denoting the number of cells, direction 3.
  !> @param eta1_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh, direction 1.
  !> @param eta1_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh, direction 1.
  !> @param eta2_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh, direction 2.
  !> @param eta2_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh, direction 2.
  !> @param eta3_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh, direction 3.
  !> @param eta3_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh, direction 3.
  !> return a pointer to the newly allocated object.
  function new_logical_mesh_3d( &
       num_cells1, &
       num_cells2, &
       num_cells3, &
       eta1_min, &
       eta1_max, &
       eta2_min, &
       eta2_max, &
       eta3_min, &
       eta3_max ) result(m)

    type(sll_logical_mesh_3d), pointer :: m
    sll_int32, intent(in)  :: num_cells1
    sll_int32, intent(in)  :: num_cells2
    sll_int32, intent(in)  :: num_cells3
    sll_real64, optional, intent(in) :: eta1_min
    sll_real64, optional, intent(in) :: eta1_max
    sll_real64, optional, intent(in) :: eta2_min
    sll_real64, optional, intent(in) :: eta2_max
    sll_real64, optional, intent(in) :: eta3_min
    sll_real64, optional, intent(in) :: eta3_max
    !sll_real64 :: delta1
    !sll_real64 :: delta2
    !sll_real64 :: delta3
    sll_int32 :: ierr

    SLL_ALLOCATE(m, ierr)
    m%num_cells1 = num_cells1
    m%num_cells2 = num_cells2
    m%num_cells3 = num_cells3
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta1_min, eta1_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta1_max, eta1_max, 1.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta2_min, eta2_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta2_max, eta2_max, 1.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta3_min, eta3_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta3_max, eta3_max, 1.0_f64 )
    m%delta_eta1   = (m%eta1_max - m%eta1_min)/real(num_cells1,f64)
    m%delta_eta2   = (m%eta2_max - m%eta2_min)/real(num_cells2,f64)
    m%delta_eta3   = (m%eta3_max - m%eta3_min)/real(num_cells3,f64)

    if ( m%eta1_max <= m%eta1_min) then
       print*,'Problem to construct the mesh 3d '
       print*,'because eta1_max <= eta1_min'
    end if
    if ( m%eta2_max <= m%eta2_min) then
       print*,'Problem to construct the mesh 3d '
       print*,'because eta2_max <= eta2_min'
    end if
    if ( m%eta3_max <= m%eta3_min) then
       print*,'Problem to construct the mesh 3d '
       print*,'because eta3_max <= eta3_min'
    end if

  end function new_logical_mesh_3d


  !> @brief allocates the memory space for a new 4D logical mesh on the heap,
  !> initializes it with the given arguments and returns a pointer to the
  !> object.
  !> @param num_cells1 integer denoting the number of cells, direction 1.
  !> @param num_cells2 integer denoting the number of cells, direction 2.
  !> @param num_cells3 integer denoting the number of cells, direction 3.
  !> @param num_cells3 integer denoting the number of cells, direction 4.
  !> @param eta1_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh, direction 1.
  !> @param eta1_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh, direction 1.
  !> @param eta2_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh, direction 2.
  !> @param eta2_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh, direction 2.
  !> @param eta3_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh, direction 3.
  !> @param eta3_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh, direction 3.
  !> @param eta3_min optional double precision value which represents the 
  !> minimum value of the eta1 parameter in the logical mesh, direction 4.
  !> @param eta3_max optional double precision value which represents the 
  !> maximum value of the eta1 parameter in the logical mesh, direction 4.
  !> return a pointer to the newly allocated object.  
  function new_logical_mesh_4d( &
    num_cells1, &
    num_cells2, &
    num_cells3, &
    num_cells4, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    eta3_min, &
    eta3_max, &
    eta4_min, &
    eta4_max ) result(m)
    
    type(sll_logical_mesh_4d), pointer :: m
    sll_int32, intent(in)  :: num_cells1
    sll_int32, intent(in)  :: num_cells2
    sll_int32, intent(in)  :: num_cells3
    sll_int32, intent(in)  :: num_cells4
    sll_real64, optional, intent(in) :: eta1_min
    sll_real64, optional, intent(in) :: eta1_max
    sll_real64, optional, intent(in) :: eta2_min
    sll_real64, optional, intent(in) :: eta2_max
    sll_real64, optional, intent(in) :: eta3_min
    sll_real64, optional, intent(in) :: eta3_max
    sll_real64, optional, intent(in) :: eta4_min
    sll_real64, optional, intent(in) :: eta4_max
    sll_int32 :: ierr

    SLL_ALLOCATE(m, ierr)
    m%num_cells1 = num_cells1
    m%num_cells2 = num_cells2
    m%num_cells3 = num_cells3
    m%num_cells4 = num_cells4
    TEST_PRESENCE_AND_ASSIGN_VAL(m, eta1_min, eta1_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL(m, eta1_max, eta1_max, 1.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL(m, eta2_min, eta2_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL(m, eta2_max, eta2_max, 1.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL(m, eta3_min, eta3_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL(m, eta3_max, eta3_max, 1.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL(m, eta4_min, eta4_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL(m, eta4_max, eta4_max, 1.0_f64 )
    m%delta_eta1   = (m%eta1_max - m%eta1_min)/real(num_cells1,f64)
    m%delta_eta2   = (m%eta2_max - m%eta2_min)/real(num_cells2,f64)
    m%delta_eta3   = (m%eta3_max - m%eta3_min)/real(num_cells3,f64)
    m%delta_eta4   = (m%eta4_max - m%eta4_min)/real(num_cells4,f64)


    if ( m%eta1_max <= m%eta1_min) then
       print*,'Problem to construct the mesh 4d '
       print*,'because eta1_max <= eta1_min'
    end if
    if ( m%eta2_max <= m%eta2_min) then
       print*,'Problem to construct the mesh 4d '
       print*,'because eta2_max <= eta2_min'
    end if
    if ( m%eta3_max <= m%eta3_min) then
       print*,'Problem to construct the mesh 4d '
       print*,'because eta3_max <= eta3_min'
    end if
    if ( m%eta4_max <= m%eta4_min) then
       print*,'Problem to construct the mesh 4d '
       print*,'because eta4_max <= eta4_min'
    end if
  end function new_logical_mesh_4d
  

  !> @brief display contents of a 1D logical mesh. Recommended access through
  !> the generic interface sll_display( mesh ).
  !> @param mesh pointer to a sll_logical_mesh_1d object.
  subroutine display_logical_mesh_1d(mesh)
    type(sll_logical_mesh_1d), pointer :: mesh

    write(*,"(/,(a))") '1D mesh : num_cell eta_min      eta_max       delta_eta'
    write(*,"(10x,(i4,1x),3(g13.3,1x))") mesh%num_cells, &
                                     mesh%eta_min,  &
                                     mesh%eta_max,  &
                                     mesh%delta_eta

  end subroutine display_logical_mesh_1d

  !> @brief display contents of a 2d logical mesh. Recommended access through
  !> the generic interface sll_display( mesh ).
  !> @param mesh pointer to a sll_logical_mesh_2d object.
  subroutine display_logical_mesh_2d(mesh)
    type(sll_logical_mesh_2d), pointer :: mesh

    write(*,"(/,(a))") '2D mesh : num_cell eta_min      eta_max       delta_eta'
    write(*,"(10x,(i4,1x),3(g13.3,1x))") mesh%num_cells1, &
                                         mesh%eta1_min,  &
                                         mesh%eta1_max,  &
                                         mesh%delta_eta1
    write(*,"(10x,(i4,1x),3(g13.3,1x))") mesh%num_cells2, &
                                         mesh%eta2_min,  &
                                         mesh%eta2_max,  &
                                         mesh%delta_eta2
  end subroutine display_logical_mesh_2d

  !> @brief display contents of a 3d logical mesh. Recommended access through
  !> the generic interface sll_display( mesh ).
  !> @param mesh pointer to a sll_logical_mesh_3d object.
  subroutine display_logical_mesh_3d(mesh)
    type(sll_logical_mesh_3d), pointer :: mesh

    write(*,"(/,(a))") '3D mesh : num_cell eta_min      eta_max       delta_eta'
    write(*,"(10x,(i4,1x),3(g13.3,1x))") mesh%num_cells1, &
                                         mesh%eta1_min,  &
                                         mesh%eta1_max,  &
                                         mesh%delta_eta1
    write(*,"(10x,(i4,1x),3(g13.3,1x))") mesh%num_cells2, &
                                         mesh%eta2_min,  &
                                         mesh%eta2_max,  &
                                         mesh%delta_eta2
    write(*,"(10x,(i4,1x),3(g13.3,1x))") mesh%num_cells3, &
                                         mesh%eta3_min,  &
                                         mesh%eta3_max,  &
                                         mesh%delta_eta3
  end subroutine display_logical_mesh_3d
  

  !> @brief display contents of a 4d logical mesh. Recommended access through
  !> the generic interface sll_display( mesh ).
  !> @param mesh pointer to a sll_logical_mesh_4d object.
  subroutine display_logical_mesh_4d(mesh)
    type(sll_logical_mesh_4d), pointer :: mesh

    write(*,"(/,(a))") '4D mesh : num_cell eta_min      eta_max       delta_eta'
    write(*,"(10x,(i4,1x),3(g13.3,1x))") mesh%num_cells1, &
                                         mesh%eta1_min,  &
                                         mesh%eta1_max,  &
                                         mesh%delta_eta1
    write(*,"(10x,(i4,1x),3(g13.3,1x))") mesh%num_cells2, &
                                         mesh%eta2_min,  &
                                         mesh%eta2_max,  &
                                         mesh%delta_eta2
    write(*,"(10x,(i4,1x),3(g13.3,1x))") mesh%num_cells3, &
                                         mesh%eta3_min,  &
                                         mesh%eta3_max,  &
                                         mesh%delta_eta3
    write(*,"(10x,(i4,1x),3(g13.3,1x))") mesh%num_cells4, &
                                         mesh%eta4_min,  &
                                         mesh%eta4_max,  &
                                         mesh%delta_eta4
  end subroutine display_logical_mesh_4d
  
  !> @brief deallocates memory for the 1D logical mesh. Recommended access 
  !> through the generic interface delete( mesh ).
  !> @param mesh pointer to a sll_logical_mesh_1d object.
  subroutine delete_logical_mesh_1d( mesh )
    type(sll_logical_mesh_1d), pointer :: mesh
    sll_int32 :: ierr
    if(.not. associated(mesh))then
       print *, 'delete_logical_mesh_1d, ERROR: passed argument is not ', &
            'associated. Crash imminent...'
    end if
    SLL_DEALLOCATE(mesh, ierr)
  end subroutine delete_logical_mesh_1d

  !> @brief deallocates memory for the 4D logical mesh. Recommended access 
  !> through the generic interface delete( mesh ).
  !> @param mesh pointer to a sll_logical_mesh_4d object.
  subroutine delete_logical_mesh_4d( mesh )
    type(sll_logical_mesh_4d), pointer :: mesh
    sll_int32 :: ierr
    if(.not. associated(mesh))then
       print *, 'delete_logical_mesh_4d, ERROR: passed argument is not ', &
            'associated. Crash imminent...'
    end if
    SLL_DEALLOCATE(mesh, ierr)
  end subroutine delete_logical_mesh_4d

  !> @brief deallocates memory for the 2D logical mesh. Recommended access 
  !> through the generic interface delete( mesh ).
  !> @param mesh pointer to a sll_logical_mesh_2d object.
  subroutine delete_logical_mesh_2d( mesh )
    type(sll_logical_mesh_2d), pointer :: mesh
    sll_int32 :: ierr
    if(.not. associated(mesh))then
       print *, 'delete_logical_mesh_2d, ERROR: passed argument is not ', &
            'associated. Crash imminent...'
    end if
    SLL_DEALLOCATE(mesh, ierr)
  end subroutine delete_logical_mesh_2d

  !> @brief deallocates memory for the 3D logical mesh. Recommended access 
  !> through the generic interface delete( mesh ).
  !> @param mesh pointer to a sll_logical_mesh_3d object.
  subroutine delete_logical_mesh_3d( mesh )
    type(sll_logical_mesh_3d), pointer :: mesh
    sll_int32 :: ierr
    if(.not. associated(mesh))then
       print *, 'delete_logical_mesh_3d, ERROR: passed argument is not ', &
            'associated. Crash imminent...'
    end if
    SLL_DEALLOCATE(mesh, ierr)
  end subroutine delete_logical_mesh_3d
  
#undef TEST_PRESENCE_AND_ASSIGN_VAL

end module sll_logical_meshes
