module distribution_function
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"

  use numeric_constants
  use advection_field
  use sll_splines
  implicit none

  type sll_distribution_function_2D_t
     type(field_2D_vec1), pointer :: field
     sll_real64 :: pcharge, pmass
  end type sll_distribution_function_2D_t

  enum, bind(C)
  enumerator :: LANDAU = 0, TWO_STREAM = 1, GAUSSIAN = 2
end enum

contains
  ! intializes some 2D distribution_functions
  function sll_new_distribution_function_2D( mesh_descriptor ) 
    type(sll_distribution_function_2D_t), pointer :: &
         sll_new_distribution_function_2D
    type(mesh_descriptor_2D), pointer :: mesh_descriptor
    sll_int32                         :: ierr
    SLL_ASSERT(associated(mesh_descriptor))
    SLL_ALLOCATE(sll_new_distribution_function_2D, ierr)
    sll_new_distribution_function_2D%field => &
         new_field_2D_vec1( mesh_descriptor )
    sll_new_distribution_function_2D%pcharge = 1.0_f64
    sll_new_distribution_function_2D%pmass = 1.0_f64
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

#undef NEW_ACCESS_FUNCTION

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

  subroutine sll_init_distribution_function_2D( dist_func_2D, test_case )
    type(sll_distribution_function_2D_t), pointer      :: dist_func_2D
    sll_int32  :: test_case
    ! local variables
    sll_int32 :: nc_eta1, nc_eta2, i1, i2
    sll_real64 :: delta_eta1, delta_eta2,  eta1_min, eta2_min
    sll_real64 :: vx, x, x2, v2, eps, kx, xi, v0, fval, xoffset, vxoffset

    nc_eta1 = get_df_nc_eta1( dist_func_2D ) 
    delta_eta1 = get_df_delta_eta1( dist_func_2D )
    eta1_min = get_df_eta1_min( dist_func_2D )
    nc_eta2 = get_df_nc_eta2( dist_func_2D ) 
    delta_eta2 = get_df_delta_eta2( dist_func_2D )
    eta2_min = get_df_eta2_min( dist_func_2D )

    select case (test_case)
    case (LANDAU)
       eps=0.0_f64
       kx=2*sll_pi/(nc_eta1*delta_eta1)
       do i2 = 1, nc_eta2 + 1
          vx = eta2_min + (i2-1) * delta_eta2
          v2 = vx*vx
          do i1 = 1, nc_eta1+1
             x = eta1_min + (i1-1) * delta_eta1
             x2 = x * x
             fval = ( 1 + eps * cos(kx*x) ) / sqrt(2*sll_pi) * exp(-0.5_f64*v2)
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
          vx = eta2_min + (i2-1) * delta_eta2
          v2 = vx*vx
          do i1=1, nc_eta1 + 1
             x = eta1_min + (i1-1) * delta_eta1
             x2 = x * x   
             fval=(1+eps*((cos(2*kx*x)+cos(3*kx*x))/1.2_f64+cos(kx*x)))* &
                  (1/sqrt(2*sll_pi))*((2-2*xi)/(3-2*xi))*(1+.5_f64*v2/(1-xi))*exp(-.5_f64*v2)
             ! fval=(1+eps*eps*cos(kx*x)+eps*cos(2*kx*x))*0.5_f64/sqrt(2*pi)*(exp(-.5_f64*(vx-v0)**2)+ exp(-.5_f64*(vx+v0)**2))
             ! fval=(1+eps*cos(kx*x))*0.5_f64/sqrt(2*pi)*(exp(-.5_f64*(vx-v0)**2)+ exp(-.5_f64*(vx+v0)**2))
             ! fval= 1/sqrt (2 * sll_pi) * ( 0.9_f64 * exp (-.5_f64 * vx * vx) + 0.2_f64 * &
             !  exp(-0.5_f64 * (vx - 4.5_f64)*(vx - 4.5_f64)/(0.5_f64 * 0.5_f64))) * (1. + 0.03_f64 * cos (0.3_f64 * x))
             ! fval=(1+eps*cos(kx*x))*1/sqrt(2*pi)*exp(-.5_f64*v2)
             ! fval=exp(-.5_f64*(x2+v2))
             call sll_set_df_val(dist_func_2D, i1, i2, fval)
          end do
       end do
    case (GAUSSIAN)
       xoffset = 1.0
       vxoffset = 1.0
       do i2 = 1, nc_eta2 + 1
          vx = eta2_min + (i2-1) * delta_eta2
          v2 = (vx - vxoffset )**2
          do i1 = 1, nc_eta1+1
             x = eta1_min + (i1-1) * delta_eta1
             x2 = (x - xoffset )**2
             fval = exp(-0.5_f64*(x2+v2))
             call sll_set_df_val(dist_func_2D, i1, i2, fval)
          end do
       end do
    end select
  end subroutine sll_init_distribution_function_2D

  subroutine write_distribution_function ( f )
    type(sll_distribution_function_2D_t), pointer      :: f
    call write_field_2d_vec1 ( f%field )
  end subroutine write_distribution_function



end module distribution_function
