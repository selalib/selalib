module distribution_function
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"

  use numeric_constants
  use sll_misc_utils   ! for int2string
  implicit none

  type sll_distribution_function_2D_t
     type(field_2D_vec1), pointer :: field
     sll_real64      :: pcharge, pmass
     sll_real64      :: average
     sll_int32       :: plot_counter
     character(32)   :: name
  end type sll_distribution_function_2D_t

  enum, bind(C)
  enumerator :: LANDAU = 0, TWO_STREAM = 1, GAUSSIAN = 2
end enum

contains
  ! intializes some 2D distribution_functions
  function sll_new_distribution_function_2D( mesh_descriptor, name ) 
    type(sll_distribution_function_2D_t), pointer :: &
         sll_new_distribution_function_2D
    type(mesh_descriptor_2D), pointer :: mesh_descriptor
    character(32)   :: name
    sll_int32                         :: ierr
    SLL_ASSERT(associated(mesh_descriptor))
    SLL_ALLOCATE(sll_new_distribution_function_2D, ierr)
    sll_new_distribution_function_2D%field => &
         new_field_2D_vec1( mesh_descriptor )
    sll_new_distribution_function_2D%pcharge = 1.0_f64
    sll_new_distribution_function_2D%pmass = 1.0_f64
    sll_new_distribution_function_2D%plot_counter = 0
    sll_new_distribution_function_2D%name = name
  end function sll_new_distribution_function_2D
  subroutine sll_delete_distribution_function( f )
    type(sll_distribution_function_2D_t), pointer      :: f
    sll_int32 :: ierr
    call delete_field_2D_vec1(f%field)
    SLL_DEALLOCATE(f, ierr)
  end subroutine sll_delete_distribution_function
#define NEW_ACCESS_FUNCTION( func_name, typ, slot ) \
  function func_name( f ) ; \
    typ :: func_name ;\
    type(sll_distribution_function_2D_t), pointer :: f ;\
    func_name = f%field%descriptor%slot ;\
  end function func_name

  NEW_ACCESS_FUNCTION(get_df_nc_eta1,     sll_int32,  nc_eta1)
  NEW_ACCESS_FUNCTION(get_df_eta1_min,   sll_real64, eta1_min)
  NEW_ACCESS_FUNCTION(get_df_eta1_max,   sll_real64, eta1_max)
  NEW_ACCESS_FUNCTION(get_df_delta_eta1, sll_real64, delta_eta1)   
  NEW_ACCESS_FUNCTION(get_df_nc_eta2,     sll_int32,  nc_eta2)
  NEW_ACCESS_FUNCTION(get_df_eta2_min,   sll_real64, eta2_min)
  NEW_ACCESS_FUNCTION(get_df_eta2_max,   sll_real64, eta2_max)
  NEW_ACCESS_FUNCTION(get_df_delta_eta2, sll_real64, delta_eta2)   
  NEW_ACCESS_FUNCTION(get_df_boundary1_type, sll_int32, boundary1_type)
  NEW_ACCESS_FUNCTION(get_df_boundary2_type, sll_int32, boundary2_type)
#undef NEW_ACCESS_FUNCTION

  function get_df_x1 ( f )
    procedure(scalar_function_2D), pointer        :: get_df_x1
    type(sll_distribution_function_2D_t), pointer :: f
    get_df_x1 => f%field%descriptor%geom%x1
  end function get_df_x1

  function get_df_x2 ( f )
    procedure(scalar_function_2D), pointer        :: get_df_x2
    type(sll_distribution_function_2D_t), pointer :: f
    get_df_x2 => f%field%descriptor%geom%x2
  end function get_df_x2

  function get_df_jac11 ( f )
    procedure(scalar_function_2D), pointer        :: get_df_jac11
    type(sll_distribution_function_2D_t), pointer :: f
    get_df_jac11 => f%field%descriptor%geom%Jacobian11
  end function get_df_jac11

  function get_df_jac12 ( f )
    procedure(scalar_function_2D), pointer        :: get_df_jac12
    type(sll_distribution_function_2D_t), pointer :: f
    get_df_jac12 => f%field%descriptor%geom%Jacobian12
  end function get_df_jac12

  function get_df_jac21 ( f )
    procedure(scalar_function_2D), pointer        :: get_df_jac21
    type(sll_distribution_function_2D_t), pointer :: f
    get_df_jac21 => f%field%descriptor%geom%Jacobian21
  end function get_df_jac21

  function get_df_jac22 ( f )
    procedure(scalar_function_2D), pointer        :: get_df_jac22
    type(sll_distribution_function_2D_t), pointer :: f
    get_df_jac22 => f%field%descriptor%geom%Jacobian22
  end function get_df_jac22  

  function get_df_jac ( f )
    procedure(scalar_function_2D), pointer        :: get_df_jac
    type(sll_distribution_function_2D_t), pointer :: f
    get_df_jac => f%field%descriptor%geom%Jacobian
  end function get_df_jac


  function sll_get_df_val( f, i, j )
    sll_real64 :: sll_get_df_val
    type(sll_distribution_function_2D_t), pointer      :: f
    sll_int32 :: i, j
    sll_get_df_val = f%field%data(i,j)
  end function sll_get_df_val

  subroutine sll_set_df_val( f, i, j, val )
    type(sll_distribution_function_2D_t), pointer      :: f
    sll_int32 :: i, j
    sll_real64 :: val
    f%field%data(i,j) = val
  end subroutine sll_set_df_val

  subroutine sll_init_distribution_function_2D( dist_func_2D, test_case, center )
    type(sll_distribution_function_2D_t), pointer      :: dist_func_2D
    sll_int32  :: test_case
    character(4) :: center   ! centering of dist_func_2D, one of ('node' or 'cell')
    ! local variables
    procedure(scalar_function_2D), pointer :: x1, x2, jac
    sll_int32 :: nc_eta1, nc_eta2, i1, i2
    sll_real64 :: delta_eta1, delta_eta2,  eta1_min, eta2_min
    sll_real64 :: x, vx, xx, vv, eps, kx, xi, v0, fval, xoffset, voffset
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: avg, avg_jac
    
    nc_eta1 = get_df_nc_eta1( dist_func_2D ) 
    delta_eta1 = get_df_delta_eta1( dist_func_2D )
    eta1_min = get_df_eta1_min( dist_func_2D )
    nc_eta2 = get_df_nc_eta2( dist_func_2D ) 
    delta_eta2 = get_df_delta_eta2( dist_func_2D )
    eta2_min = get_df_eta2_min( dist_func_2D )
    x1 => get_df_x1 ( dist_func_2D )
    x2 => get_df_x2 ( dist_func_2D )
    jac => get_df_jac ( dist_func_2D )

    if (center=='cell') then ! half cell offset
       eta1_min = eta1_min + 0.5_f64 * delta_eta1
       eta2_min = eta2_min + 0.5_f64 * delta_eta2
    end if

    select case (test_case)
    case (LANDAU)
       eps=0.0_f64
       kx=2*sll_pi/(nc_eta1*delta_eta1)
       do i2 = 1, nc_eta2 + 1
          vx = eta2_min + (i2-0.5_f64) * delta_eta2
          vv = vx*vx
          do i1 = 1, nc_eta1+1
             x = eta1_min + (i1-0.5_f64) * delta_eta1
             xx = x * x
             fval = ( 1 + eps * cos(kx*x) ) / sqrt(2*sll_pi) * exp(-0.5_f64*vv)
             call sll_set_df_val(dist_func_2D, i1, i2, fval)
          end do
       end do
    case (TWO_STREAM)
       xi = 0.90_f64
       eps = 0.001_f64
       !eps = 0.05_f64
       v0 = 2.4_f64
       kx = 2 * sll_pi / (nc_eta1 * delta_eta1)
       do i2 = 1, nc_eta2+1
          vx = eta2_min + (i2-0.5_f64) * delta_eta2
          vv = vx*vx
          do i1=1, nc_eta1 + 1
             x = eta1_min + (i1-0.5_f64) * delta_eta1
             xx = x * x   
             fval=(1+eps*((cos(2*kx*x)+cos(3*kx*x))/1.2_f64+cos(kx*x)))* &
                  (1/sqrt(2*sll_pi))*((2-2*xi)/(3-2*xi))*(1+.5_f64*vv/(1-xi))*exp(-.5_f64*vv)
             ! fval=(1+eps*eps*cos(kx*x)+eps*cos(2*kx*x))*0.5_f64/sqrt(2*pi)*(exp(-.5_f64*(vx-v0)**2)+ exp(-.5_f64*(vx+v0)**2))
             ! fval=(1+eps*cos(kx*x))*0.5_f64/sqrt(2*pi)*(exp(-.5_f64*(vx-v0)**2)+ exp(-.5_f64*(vx+v0)**2))
             ! fval= 1/sqrt (2 * sll_pi) * ( 0.9_f64 * exp (-.5_f64 * vx * vx) + 0.2_f64 * &
             !  exp(-0.5_f64 * (vx - 4.5_f64)*(vx - 4.5_f64)/(0.5_f64 * 0.5_f64))) * (1. + 0.03_f64 * cos (0.3_f64 * x))
             ! fval=(1+eps*cos(kx*x))*1/sqrt(2*pi)*exp(-.5_f64*vv)
             ! fval=exp(-.5_f64*(xx+vv))
             call sll_set_df_val(dist_func_2D, i1, i2, fval)
          end do
       end do
    case (GAUSSIAN)
       xoffset = 1.0_f64
       voffset = 1.0_f64
       avg = 0.0_f64
       avg_jac = 0.0_f64
       eta2 = eta2_min
       do i2 = 1, nc_eta2 
          eta1 = eta1_min
          do i1 = 1, nc_eta1           
             xx = (x1(eta1,eta2) - xoffset )**2
             vv = (x2(eta1,eta2) - voffset )**2
             ! store f * jac for conservative cell centered schemes
             fval =  jac(eta1,eta2) !exp(-0.5_f64*(xx+vv)) * jac(eta1,eta2)
             avg = avg + fval
             avg_jac = avg_jac + jac(eta1,eta2)
             call sll_set_df_val(dist_func_2D, i1, i2, fval)
             eta1 = eta1 +  delta_eta1
          end do
          eta2 = eta2 +  delta_eta2
       end do
       print*, 'dist_func, averages ', avg*delta_eta1*delta_eta2, avg_jac*delta_eta1*delta_eta2
       dist_func_2D%average = 0.0 !avg / avg_jac
       ! subtract average from distribution function
       avg = dist_func_2D%average
       eta2 = eta2_min
       do i2 = 1, nc_eta2 + 1
          eta1 = eta1_min
          do i1 = 1, nc_eta1+1
             fval  =  sll_get_df_val( dist_func_2D, i1, i2 ) - avg * jac(eta1,eta2)
             call sll_set_df_val(dist_func_2D, i1, i2, fval)
             eta1 = eta1 +  delta_eta1
          end do
          eta2 = eta2 +  delta_eta2
       end do
    end select
  end subroutine sll_init_distribution_function_2D

  subroutine write_distribution_function ( f )
    type(sll_distribution_function_2D_t), pointer      :: f
    character(len=4) :: counter
    character(64) :: name
    logical, parameter   :: jacobian = .true.
    call int2string(f%plot_counter,counter)
    name = trim(f%name)//counter
    call write_field_2d_vec1 ( f%field, name, jacobian, f%average )
    f%plot_counter = f%plot_counter + 1
  end subroutine write_distribution_function



end module distribution_function
