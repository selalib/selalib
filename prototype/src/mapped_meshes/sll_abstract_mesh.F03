module m

#define sll_pi 3.14159

  type, abstract :: mesh_2d
    integer :: nc1,nc2
    double precision :: delta1, delta2
  end type mesh_2d

  type, extends(mesh_2d), abstract :: mesh_2d_discrete
!    type(jacobian_matrix_element), dimension(:,:), pointer :: j_matrix
!    type(interpolator_2d) :: interp_2d    

    contains
    
      procedure(geometry_function), deferred, pass(mesh)        :: x1
      procedure(geometry_function), deferred, pass(mesh)        :: x2
      procedure(geometry_function), deferred, pass(mesh)        :: jacobian
      procedure(matrix_geometry_function), deferred, pass(mesh) :: jacobian_matrix
!      procedure :: write_mesh
  end type mesh_2d_discrete

  type, extends(mesh_2d) :: mesh_2d_analytic
!    type(jacobian_matrix_element), dimension(:,:), pointer :: j_matrix
!    type(interpolator_2d) :: interp_2d    
   
      procedure(geometry_function_nopass), pointer, nopass        :: x1
      procedure(geometry_function_nopass), pointer, nopass        :: x2
      procedure(geometry_function_nopass), pointer, nopass        :: jacobian
      procedure(matrix_geometry_function_nopass), pointer, nopass :: jacobian_matrix
!      procedure :: write_mesh
   
      contains
        procedure, pass(mesh) :: new => new_analytic 
  end type mesh_2d_analytic

     
     

!  abstract interface
!     subroutine new_mesh_2d( mesh, n1, n2 )
!       import mesh_2d_discrete
!       class(mesh_2d_discrete)      :: mesh
!       integer, intent(in) :: n1 
!       integer, intent(in) :: n2
!     end subroutine
!  end interface

  abstract interface
     function geometry_function( mesh, eta1, eta2 ) result(res)
       import mesh_2d_discrete
       class(mesh_2d_discrete)               :: mesh
       double precision, intent(in) :: eta1
       double precision, intent(in) :: eta2
       double precision             :: res
     end function 
  end interface

  abstract interface
     function matrix_geometry_function( mesh, eta1, eta2 ) result(res)
       import mesh_2d_discrete
       class(mesh_2d_discrete)               :: mesh
       double precision, intent(in) :: eta1
       double precision, intent(in) :: eta2
       double precision             :: res(2,2)
     end function 
  end interface

  abstract interface
     function geometry_function_nopass( eta1, eta2 ) result(res)
       double precision, intent(in) :: eta1
       double precision, intent(in) :: eta2
       double precision             :: res
     end function 
  end interface

  abstract interface
     function matrix_geometry_function_nopass( eta1, eta2 ) result(res)
       double precision, intent(in) :: eta1
       double precision, intent(in) :: eta2
       double precision             :: res(2,2)
     end function 
  end interface

contains

   subroutine new_analytic(mesh,n1,n2,x1,x2)
     class(mesh_2d_analytic)             :: mesh
     integer                             :: n1,n2
     procedure(geometry_function_nopass) :: x1, x2
     mesh%NC1 = n1
     mesh%NC2 = n2
     mesh%x1 => x1
     mesh%x2 => x2
   end subroutine

end module m







module sll_polar_mesh_2d
  !use m
  
  contains

#define R1 0.1d0
#define R2 1.0d0
  function x1_polar_f( eta1, eta2 )
    double precision             :: x1_polar_f
    double precision, intent(in) :: eta1, eta2
    x1_polar_f = (R1 + (R2-R1)*eta1)*cos(2.0d0*sll_pi*eta2)
  end function x1_polar_f

  function x2_polar_f( eta1, eta2 )
    double precision :: x2_polar_f
    double precision, intent(in) :: eta1, eta2
    x2_polar_f = (R1 + (R2-R1)*eta1)*sin(2.0d0*sll_pi*eta2)
  end function x2_polar_f

  function jacobian_polar_f( eta1, eta2 )
    double precision :: jacobian_polar_f
    double precision, intent(in) :: eta1, eta2
    jacobian_polar_f = 2*sll_pi*(R2-R1)*(R1 + (R2-R1)*eta1)
  end function jacobian_polar_f


  function jacobian_matrix_polar_f( eta1, eta2 ) result(res)
    double precision             :: res(2,2)
    double precision, intent(in) :: eta1, eta2
    double precision :: k

    k = 2.0d0*sll_pi
    res(1,1) = (R2-R1)*cos(k*eta2)
    res(2,1) = (R2-R1)*sin(k*eta2)
    res(1,2) = -(R1+(R2-R1)*eta1)*sin(k*eta2)*k
    res(2,2) = (R1+(R2-R1)*eta1)*cos(k*eta2)*k
  end function jacobian_matrix_polar_f


end module sll_polar_mesh_2d















program test
  use m
  use sll_polar_mesh_2d

  type(mesh_2d_analytic) :: mm
  call mm%new(5,6,x1_polar_f,x2_polar_f)
  print*,mm%x1(0.5d0,0.5d0)
  print*,mm%x2(0.5d0,0.5d0)

end program test
