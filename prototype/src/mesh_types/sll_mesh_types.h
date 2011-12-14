#ifndef _mesh_types_h_
#define _mesh_types_h_

#include "sll_assert.h"
  use sll_mesh_types
! Why this with macros? We explore hiding the access to an array behind a
! defined interface but at the same time we don't want to be penalized by
! function calls. Access functions like these are expected to live within
! inner loops and thus would have an appreciable overhead. Here we also 
! hide some protection against incorrect array indexing which should not
! represent any additional cost if the DEBUG flag is turned off.

! get macros
#define GET_FIELD_MESH_DESCRIPTOR( f )     f%descriptor
#define GET_FIELD_ETA1_MIN( f )            f%descriptor%eta1_min
#define GET_FIELD_ETA1_MAX( f )            f%descriptor%eta1_max
#define GET_FIELD_NC_ETA1( f )             f%descriptor%nc_eta1
#define GET_FIELD_DELTA_ETA1( f )          f%descriptor%delta_eta1
#define GET_FIELD_ETA2_MIN( f )            f%descriptor%eta2_min
#define GET_FIELD_ETA2_MAX( f )            f%descriptor%eta2_max
#define GET_FIELD_NC_ETA2( f )             f%descriptor%nc_eta2
#define GET_FIELD_DELTA_ETA2( f )          f%descriptor%delta_eta2

! There are no SET_FIELD() macros since one is not supposed to set those 
! values outside of the initialization.
#define FIELD_1D_AT_I( f, index ) \
     f%data(index)
#define FIELD_2D_AT_I( f, index1, index2 ) \
     f%data(index1,index2)
#define FIELD_2D_AT_I_V1( f, index1, index2 ) \
     f%data(index1,index2)%v1
#define FIELD_2D_AT_I_V2( f, index1, index2 ) \
     f%data(index1,index2)%v2

! The intent of this macro is to pass all the data to functions that require it, eg FFT
! It works for all types of fields
#define FIELD_DATA(f)   f%data


! The following macro should call some interpolation with splines or something
#define GET_FIELD_1D_VALUE_AT_X( f, x )      
#define SET_FIELD_1D_VALUE_AT_I( f, index, val )  f%data(index) = val

#define GET_MESH_RMIN( m )        m%rmin
#define GET_MESH_RMAX( m )        m%rmax
#define GET_MESH_ZMIN( m )        m%zmin
#define GET_MESH_ZMAX( m )        m%zmax
#define GET_MESH_NCR( m )         m%ncr
#define GET_MESH_NCTHETA( m )     m%nctheta
#define GET_MESH_NCZ( m )         m%ncz
#define GET_MESH_DELTA_R( m )     m%delta_r
#define GET_MESH_DELTA_THETA( m ) m%delta_theta
#define GET_MESH_DELTA_Z( m )     m%delta_z


#define GET_MESH_NC_ETA1( m )     m%nc_eta1
#define GET_MESH_NC_ETA2( m )     m%nc_eta2

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

#define GET_MESH_VALUE( m,i,j,k )                                \
   SLL_ASSERT( (i .ge. 1) .and. (i .le. GET_MESH_NCR(m)+1))      \
   SLL_ASSERT( (j .ge. 1) .and. (j .le. GET_MESH_NCTHETA(m)+1))  \
   SLL_ASSERT( (k .ge. 1) .and. (k .le. GET_MESH_NCZ(m)+1))      \
   m%data(i,j,k)

! set macros
#define SET_MESH_RMIN( m, val )    m%rmin = val
#define SET_MESH_RMAX( m, val )    m%rmax = val
#define SET_MESH_ZMIN( m, val )    m%zmin = val
#define SET_MESH_ZMAX( m, val )    m%zmax = val
#define SET_MESH_NCR( m, val )     m%ncr = val
#define SET_MESH_NCTHETA( m,val )  m%nctheta = val
#define SET_MESH_NCZ( m, val )     m%ncz = val
! note: the values of delta can't be set, they are defined and written
! when the instance is created.

#define SET_MESH_VALUE( m, i, j, k, val )                        \
   SLL_ASSERT( (i .ge. 1) .and. (i .le. GET_MESH_NCR(m)+1))      \
   SLL_ASSERT( (j .ge. 1) .and. (j .le. GET_MESH_NCTHETA(m)+1))  \
   SLL_ASSERT( (k .ge. 1) .and. (k .le. GET_MESH_NCZ(m)+1))      \
   m%data(i,j,k) = val

#endif
