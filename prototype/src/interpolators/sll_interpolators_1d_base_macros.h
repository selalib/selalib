#ifndef _sll_interpolators_1d_base_macros_h_
#define _sll_interpolators_1d_base_macros_h_

#include "sll_utilities.h"

  ! NOTE: THIS IS EXPERIMENTAL...
  ! In case that a developer wants to write an interpolator using this
  ! base class, but at the same time does not want to implement all the 
  ! functions, here we wish to provide some functions with the required
  ! signatures whose only side-effect is to announce that they should not be 
  ! called with a subsequent stopping of the program.
  !
  ! Unfortunately, each interpolator implementation (module) must write its
  ! own functions of this kind, the reason being that the interpolator
  ! type needs to be specifically of the child class that is being defined
  ! in the module. Here we offer a macro-based solution. The idea is that
  ! each child module, when in need of any/all of the null functions, can
  ! call the specific macro and have an implementation for it. This would 
  ! avoid a lot of repetitive work. If we had the full capabilities of a
  ! C preprocessor, a single macro would be needed, since we could build
  ! the specific function names on the fly. Unfortunately, we do not have 
  ! the ## operator and thus we can't build new names of functions. Therefore
  ! we are stuck with the solution of offering multiple macros.

#define DEFINE_NULL_INTERP_ONE_ARG_MSG( child_class, func_name )      \
  function func_name( interpolator, eta1 ) result(val);               \
    class(child_class), intent(inout) :: interpolator;                \
    sll_real64, intent(in) :: eta1;                                   \
    sll_real64 :: val;                                                \
    val = 0.0_F64;                                                    \
    print *, STRNG(func_name), ': this function is not meant to ',    \
         'be called. Need an actual implementation.';                 \
    stop;                                                             \
  end function func_name

#define DEFINE_NULL_INTERP_1D_ARRAY_MSG( child_class, func_name )     \
  subroutine func_name( interpolator, data_array );                   \
    class(child_class), intent(inout) :: interpolator;                \
    sll_real64, dimension(:), intent(in) :: data_array;               \
    print *, STRNG(func_name), ': this function is not meant to ',    \
         'be called. Need an actual implementation.';                 \
    stop;                                                             \
  end subroutine func_name

#define DEFINE_NULL_INTERP_1D_ARRAY_SUB( child_class, func_name)      \
  subroutine func_name(interpolator, num_pts, vals_to_interpolate, output_array); \
    class(child_class), intent(in) :: interpolator;                   \
    sll_int32, intent(in)                       :: num_pts;           \
    sll_real64, dimension(:), intent(in)        :: vals_to_interpolate; \
    sll_real64, dimension(:), intent(out)       :: output_array;      \
    print *, STRNG(func_name), ': this function is not meant to ',    \
         'be called. Need an actual implementation.';                 \
    stop;                                                             \
  end subroutine func_name

#define DEFINE_NULL_INTERP_1D_POINTER_SUB( child_class, func_name )   \
  subroutine func_name(interpolator, num_pts, vals_to_interpolate, output); \
    class(child_class), intent(in) :: interpolator;                   \
    sll_int32, intent(in)                       :: num_pts;           \
    sll_real64, dimension(:), pointer           :: vals_to_interpolate; \
    sll_real64, dimension(:), pointer           :: output;            \
    print *, STRNG(func_name), ': this function is not meant to ',    \
         'be called. Need an actual implementation.';                 \
    stop;                                                             \
  end subroutine func_name

#define DEFINE_NULL_INTERP_1D_ARRAY( child_class, func_name )         \
  function func_name(this, num_points, data, coordinates)             \
       result(res);                                                   \
    class(child_class), intent(in)     :: this;                       \
    sll_int32, intent(in)  :: num_points;                             \
    sll_real64, dimension(:), intent(in) :: data;                     \
    sll_real64, dimension(:), intent(in) :: coordinates;              \
    sll_real64, dimension(num_points)    :: res;                      \
    res = 0.0_F64;                                                    \      
    print *, STRNG(func_name), ': this function is not meant to ',    \
         'be called. Need an actual implementation.';                 \
    stop;                                                             \
  end function func_name

#define DEFINE_NULL_RECONSTRUCT_1D_ARRAY( child_class, func_name )    \
  function func_name(this, num_points, data) result(res);             \
    class(child_class), intent(in)     :: this;                       \
    sll_int32, intent(in)     :: num_points;                          \
    sll_real64, dimension(:), intent(in)  :: data;                    \
    sll_real64, dimension(num_points)     :: res;                     \
    res = 0.0_F64;                                                    \
    print *, STRNG(func_name), ': this function is not meant to ',    \
         'be called. Need an actual implementation.';                 \
    stop;                                                             \
  end function func_name

#define DEFINE_NULL_INTERPOLATE_1D_DISP( child_class, func_name )     \
  function func_name( this, num_points, data, alpha ) result(res);    \
    class(child_class), intent(in) :: this;                           \
    sll_int32, intent(in)     :: num_points;                          \
    sll_real64, dimension(num_points)     :: res;                     \
    sll_real64, dimension(:), intent(in)  :: data;                    \
    sll_real64, intent(in) :: alpha;                                  \
    res = 0.0_F64;                                                    \
    print *, STRNG(func_name), ': this function is not meant to ',    \
         'be called. Need an actual implementation.';                 \
    stop;                                                             \
  end function func_name


#endif
