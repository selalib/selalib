#ifndef field_1d_h_
#define field_1d_h_
  use sll_scalar_field_1d
! Why this with macros? We explore hiding the access to an array behind a
! defined interface but at the same time we don't want to be penalized by
! function calls. Access functions like these are expected to live within
! inner loops and thus would have an appreciable overhead. 

! get macros
#define GET_FIELD_MESH_MESH( f )           f%mesh
#define GET_FIELD_ETA1_MIN( f )            f%mesh%eta1_min
#define GET_FIELD_ETA1_MAX( f )            f%mesh%eta1_max
#define GET_FIELD_NC_ETA1( f )             f%mesh%nc_eta1
#define GET_FIELD_DELTA_ETA1( f )          f%mesh%delta_eta1
#define GET_FIELD_BOUNDARY1( f )           f%mesh%boundary1_type

! There are no SET_FIELD() macros since one is not supposed to set those 
! values outside of the initialization.
#define FIELD_1D_AT_I( f, index1 )         f%data(index1)

! The intent of this macro is to pass all the data to functions
! that require it, eg FFT
! It works for all types of fields
#define FIELD_DATA(f)   f%data


! The following macro should call some interpolation with splines or something
#define GET_FIELD_1D_VALUE_AT_X( f, x )      
#define SET_FIELD_1D_VALUE_AT_I( f, index, val )  f%data(index) = val

#define GET_MESH_NC_ETA1( m )     m%nc_eta1

! NOTE: if the DEBUG flag is not set, the assertions get expanded into nothing.
! For this reason, the semicolon has to be part of the macro expansion. If
! not, we would end up in this case with dangling semicolons, yielding
! compilation errors. This means that there are situations in which we want
! the semicolons to be part of the macro and others in which they are not.
! So much for consistency...
!
! At the same time, SLL_ASSERT is an exceptional case in which it is the only
! macro that is meant to expand into nothing. All other macros could be 
! written without the semicolons in their definition, thus maintaining some
! resemblance of uniformity in the syntax.

#endif