module sll_mesh_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_constants
  use geometry_functions 
  implicit none
  
  integer, parameter :: PERIODIC_MESH = 0, COMPACT_MESH = 1
  !enum, bind(C)
  !   enumerator :: PERIODIC_MESH = 0, COMPACT_MESH = 1
  !end enum

  type, abstract :: mesh_2d
     sll_int32  :: nc_eta1
     sll_int32  :: nc_eta2
     sll_real64 :: delta_eta1
     sll_real64 :: delta_eta2
     sll_int32  :: boundary1_type
     sll_int32  :: boundary2_type
     sll_real64, dimension(:,:), allocatable :: x1_array
     sll_real64, dimension(:,:), allocatable :: x2_array
   contains
     procedure(real_function_2D_int), pass, deferred :: x1_at_node
     procedure(real_function_2D_int), pass, deferred :: x2_at_node
     procedure(real_function_2D_int), pass, deferred :: jac_at_node
     procedure(real_function_2D_real), pass, deferred :: x1
     procedure(real_function_2D_real), pass, deferred :: x2
     procedure(real_function_2D_real), pass, deferred :: jac
  end type mesh_2d
     
  abstract interface
     function real_function_2D_int( this, i1, i2 )
       use sll_working_precision
       import mesh_2d
       class(mesh_2d) :: this
       sll_real64 :: real_function_2D_int
       sll_int32, intent(in)  :: i1
       sll_int32, intent(in)  :: i2
     end function real_function_2D_int
     function real_function_2D_real( this, eta1, eta2 )
       use sll_working_precision
       import mesh_2d
       class(mesh_2d) :: this
       sll_real64 :: real_function_2D_real
       sll_real64, intent(in)  :: eta1
       sll_real64, intent(in)  :: eta2
     end function real_function_2D_real
     function scalar_function_mesh_2d( eta1, eta2 )
       use sll_working_precision
       sll_real64 :: scalar_function_mesh_2d
       sll_real64, intent(in)  :: eta1
       sll_real64, intent(in)  :: eta2
     end function scalar_function_mesh_2d
  end interface

  type, extends(mesh_2d) :: mesh_2d_analytic
     procedure(scalar_function_mesh_2d), pointer, nopass :: x1_func
     procedure(scalar_function_mesh_2d), pointer, nopass :: x2_func
     procedure(scalar_function_mesh_2d), pointer, nopass :: Jacobian_func
   contains
     procedure, pass :: x1_at_node => x1_int
     procedure, pass :: x2_at_node => x2_int
     procedure, pass :: jac_at_node => jac_int
     procedure, pass :: x1 => x1_real
     procedure, pass :: x2 => x2_real
     procedure, pass :: jac => jac_real
  end type mesh_2d_analytic


contains   ! *****************************************************************
  

  subroutine new_mesh_2d_analytic ( this, nc1, nc2, x1, x2, jac, bt1, bt2 )
    class(mesh_2d_analytic)  :: this 
    sll_int32 :: nc1
    sll_int32 :: nc2
    procedure(scalar_function_mesh_2d), pointer :: x1, x2, jac
    sll_int32, optional :: bt1
    sll_int32, optional :: bt2
    ! local variables
    sll_int32 :: ierr
    sll_int32 :: i1, i2
    sll_real64 :: eta1, eta2

    this%nc_eta1 = nc1
    this%nc_eta2 = nc2
    this%delta_eta1 = 1.0_f64 / nc1
    this%delta_eta2 = 1.0_f64 / nc2
    
    this%x1_func => x1
    this%x2_func => x2
    this%Jacobian_func => jac
    if (present(bt1)) then
       this%boundary1_type = bt1
    end if
    if (present(bt2)) then
       this%boundary2_type = bt2
    end if

    SLL_ALLOCATE(this%x1_array(nc1+1,nc2+1),ierr)
    SLL_ALLOCATE(this%x2_array(nc1+1,nc2+1),ierr)
    do i2=1, nc1+1
       do i1 = 1, nc2+1
          eta1 = real(i1-1,f64)/this%nc_eta1
          eta2 = real(i2-1,f64)/this%nc_eta2
          this%x1_array(i1,i2) = x1(eta1,eta2)
          this%x2_array(i1,i2) = x2(eta1,eta2)
       end do
    end do
    ! set sll_kx (this is correct only for tensor product meshes)
    sll_kx = 2*sll_pi/(this%x1_array(nc1+1,1)-this%x1_array(1,1))
  end subroutine new_mesh_2d_analytic

  function x1_int(this,i1,i2)
    class(mesh_2d_analytic)  :: this
    sll_int32 :: i1, i2
    sll_real64 :: x1_int
    ! local variables
    sll_real64 :: eta1, eta2

    eta1 = real(i1-1,f64)/this%nc_eta1
    eta2 = real(i2-1,f64)/this%nc_eta2
    x1_int = this%x1_func(eta1, eta2)
  end function x1_int
  function x2_int(this,i1,i2)
    class(mesh_2d_analytic)  :: this
    sll_int32 :: i1, i2
    sll_real64 :: x2_int
    ! local variables
    sll_real64 :: eta1, eta2

    eta1 = real(i1-1,f64)/this%nc_eta1
    eta2 = real(i2-1,f64)/this%nc_eta2
    x2_int = this%x2_func(eta1, eta2)
  end function x2_int
  function jac_int(this,i1,i2)
    class(mesh_2d_analytic)  :: this
    sll_int32 :: i1, i2
    sll_real64 :: jac_int
    ! local variables
    sll_real64 :: eta1, eta2

    eta1 = real(i1-1,f64)/this%nc_eta1
    eta2 = real(i2-1,f64)/this%nc_eta2
    jac_int = this%Jacobian_func(eta1, eta2)
  end function jac_int
  function x1_real(this,eta1,eta2)
    class(mesh_2d_analytic)  :: this
    sll_real64 :: eta1, eta2
    sll_real64 :: x1_real

    x1_real = this%x1_func(eta1, eta2)
  end function x1_real
  function x2_real(this,eta1,eta2)
    class(mesh_2d_analytic)  :: this
    sll_real64 :: eta1, eta2
    sll_real64 :: x2_real
  
    x2_real = this%x2_func(eta1, eta2)
  end function x2_real
  function jac_real(this,eta1,eta2)
    class(mesh_2d_analytic)  :: this
    sll_real64 :: eta1, eta2
    sll_real64 :: jac_real

    jac_real = this%Jacobian_func(eta1, eta2)
  end function jac_real
  
end module sll_mesh_2d
