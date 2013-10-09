module sll_logical_meshes
#include "sll_working_precision.h"
#include "sll_memory.h"
  implicit none

   type sll_logical_mesh_1d
     sll_int32  :: num_cells
     sll_real64 :: eta_min
     sll_real64 :: eta_max
     sll_real64 :: delta_eta
  end type sll_logical_mesh_1d


  type sll_logical_mesh_2d
     sll_int32  :: num_cells1
     sll_int32  :: num_cells2
     sll_real64 :: eta1_min
     sll_real64 :: eta1_max
     sll_real64 :: eta2_min
     sll_real64 :: eta2_max
     sll_real64 :: delta_eta1
     sll_real64 :: delta_eta2
  end type sll_logical_mesh_2d

  type sll_logical_mesh_3d
     sll_int32  :: num_cells1
     sll_int32  :: num_cells2
     sll_int32  :: num_cells3     
     sll_real64 :: eta1_min
     sll_real64 :: eta1_max
     sll_real64 :: eta2_min
     sll_real64 :: eta2_max
     sll_real64 :: eta3_min
     sll_real64 :: eta3_max
     sll_real64 :: delta_eta1
     sll_real64 :: delta_eta2
     sll_real64 :: delta_eta3
  end type sll_logical_mesh_3d



  type sll_logical_mesh_4d
     sll_int32  :: num_cells1
     sll_int32  :: num_cells2
     sll_int32  :: num_cells3
     sll_int32  :: num_cells4
     sll_real64 :: eta1_min
     sll_real64 :: eta1_max
     sll_real64 :: eta2_min
     sll_real64 :: eta2_max
     sll_real64 :: eta3_min
     sll_real64 :: eta3_max
     sll_real64 :: eta4_min
     sll_real64 :: eta4_max
     sll_real64 :: delta_eta1
     sll_real64 :: delta_eta2
     sll_real64 :: delta_eta3
     sll_real64 :: delta_eta4
  end type sll_logical_mesh_4d

  ! this should be sll_delete library-wide...
  interface delete
     module procedure delete_logical_mesh_4d, delete_logical_mesh_2d,delete_logical_mesh_3d
  end interface delete

contains

#define TEST_PRESENCE_AND_ASSIGN_VAL( obj, arg, slot, default_val ) \
  if( present(arg) ) then ; \
    obj%slot = arg; \
  else; \
    obj%slot = default_val; \
end if


  function new_logical_mesh_1d( &
    num_cells, &
    eta_min, &
    eta_max ) result(m)

    type(sll_logical_mesh_1d), pointer :: m
    sll_int32, intent(in)  :: num_cells
    sll_real64, optional, intent(in) :: eta_min
    sll_real64, optional, intent(in) :: eta_max
    sll_real64 :: delta
    sll_int32 :: ierr

    SLL_ALLOCATE(m, ierr)
    m%num_cells = num_cells
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta_min, eta_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta_max, eta_max, 1.0_f64 )
    m%delta_eta   = (m%eta_max - m%eta_min)/real(num_cells,f64)
  end function new_logical_mesh_1d


  subroutine initialize_x1_node_1d( m, x1_node )
    type(sll_logical_mesh_1d), pointer :: m
    sll_real64, dimension(:), pointer :: x1_node
    sll_int32  :: num_cells
    sll_real64 :: eta_min
    sll_real64 :: eta_max
    sll_real64 :: delta_eta
    sll_int32 :: i
    sll_int32 :: ierr
    
    num_cells = m%num_cells
    eta_min = m%eta_min
    delta_eta = m%delta_eta
    SLL_ALLOCATE(x1_node(num_cells+1), ierr)
    do i=1,num_cells+1
      x1_node(i) = eta_min+real(i-1,f64)*delta_eta
    enddo    
    
    
  end subroutine initialize_x1_node_1d



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
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_int32 :: ierr

    SLL_ALLOCATE(m, ierr)
    m%num_cells1 = num_cells1
    m%num_cells2 = num_cells2
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta1_min, eta1_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta1_max, eta1_max, 1.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta2_min, eta2_min, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, eta2_max, eta2_max, 1.0_f64 )
    m%delta_eta1   = (m%eta1_max - m%eta1_min)/real(num_cells1,f64)
    m%delta_eta2   = (m%eta2_max - m%eta2_min)/real(num_cells2,f64)
  end function new_logical_mesh_2d


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
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
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
  end function new_logical_mesh_3d
  
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
  end function new_logical_mesh_4d


  subroutine delete_logical_mesh_4d( mesh )
    type(sll_logical_mesh_4d), pointer :: mesh
    sll_int32 :: ierr
    SLL_DEALLOCATE(mesh, ierr)
  end subroutine delete_logical_mesh_4d

  subroutine delete_logical_mesh_2d( mesh )
    type(sll_logical_mesh_2d), pointer :: mesh
    sll_int32 :: ierr
    SLL_DEALLOCATE(mesh, ierr)
  end subroutine delete_logical_mesh_2d

  subroutine delete_logical_mesh_3d( mesh )
    type(sll_logical_mesh_3d), pointer :: mesh
    sll_int32 :: ierr
    SLL_DEALLOCATE(mesh, ierr)
  end subroutine delete_logical_mesh_3d
  
#undef TEST_PRESENCE_AND_ASSIGN_VAL

end module sll_logical_meshes
