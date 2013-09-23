#ifndef _sll_fftw_h
#define _sll_fftw_h

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

#ifdef FFTW_F2003

#define fftw_plan type(C_PTR) 
#define fftw_comp complex(C_DOUBLE_COMPLEX)
#define fftw_int  integer(C_SIZE_T)

#define FFTW_ALLOCATE(array,array_size,sz_array,p_array)     \
sz_array = int(array_size,C_SIZE_T);                         \
p_array = fftw_alloc_complex(sz_array);                      \
call c_f_pointer(p_array, array, [array_size])               \

#define NEW_FFTW_PLAN_1D(plan,n,in,out,direction)            \
plan = fftw_plan_dft_1d(n,in,out,direction,FFTW_PATIENT)

#define NEW_FFTW_PLAN_R2C_1D(plan,n,in,out)                  \
plan = fftw_plan_dft_r2c_1d(n,in,out,FFTW_PATIENT)

#define NEW_FFTW_PLAN_C2R_1D(plan,n,in,out)                  \
plan = fftw_plan_dft_c2r_1d(n,in,out,FFTW_PATIENT)

#define NEW_FFTW_PLAN_R2HC_1D(plan,n,in,out)                 \
plan = fftw_plan_r2r_1d(n,in,out,FFTW_R2HC,FFTW_PATIENT)

#define NEW_FFTW_PLAN_HC2R_1D(plan,n,in,out)                 \
plan = fftw_plan_r2r_1d(n,in,out,FFTW_HC2R,FFTW_PATIENT)

#define NEW_FFTW_PLAN_2D(plan,n1,n2,in,out,direction)        \
plan = fftw_plan_dft_2d(n2,n1,in,out,direction,FFTW_PATIENT)

#define NEW_FFTW_PLAN_R2C_2D(plan,n1,n2,in,out)              \
plan = fftw_plan_dft_r2c_2d(n2,n1,in,out,FFTW_PATIENT)

#define NEW_FFTW_PLAN_C2R_2D(plan,n1,n2,in,out)              \
plan = fftw_plan_dft_c2r_2d(n2,n1,in,out,FFTW_PATIENT)

#else

#define FFTW_ALLOCATE(array,array_size,sz_array,p_array)     \
SLL_ALLOCATE(array(array_size), error)

#define NEW_FFTW_PLAN_1D(plan,n,in,out,direction)            \
call dfftw_plan_dft_1d(plan,n,in,out,direction,FFTW_PATIENT)

#define NEW_FFTW_PLAN_R2C_1D(plan,n,in,out)                  \
call dfftw_plan_dft_r2c_1d(plan,n,in,out,FFTW_PATIENT)

#define NEW_FFTW_PLAN_C2R_1D(plan,n,in,out)                  \
call dfftw_plan_dft_c2r_1d(plan,n,in,out,FFTW_PATIENT)

#define NEW_FFTW_PLAN_R2HC_1D(plan,n,in,out)                 \
call dfftw_plan_r2r_1d(plan,n,in,out,FFTW_R2HC,FFTW_PATIENT)

#define NEW_FFTW_PLAN_HC2R_1D(plan,n,in,out)                 \
call dfftw_plan_r2r_1d(plan,n,in,out,FFTW_HC2R,FFTW_PATIENT)

#define NEW_FFTW_PLAN_2D(plan,n1,n2,in,out,direction)        \
call dfftw_plan_dft_2d(plan,n1,n2,in,out,direction,FFTW_PATIENT)

#define NEW_FFTW_PLAN_R2C_2D(plan,n1,n2,in,out)              \
call dfftw_plan_dft_r2c_2d(plan,n1,n2,in,out,FFTW_PATIENT)

#define NEW_FFTW_PLAN_C2R_2D(plan,n1,n2,in,out)              \
call dfftw_plan_dft_c2r_2d(plan,n1,n2,in,out,FFTW_PATIENT)

#define fftw_plan  sll_int64
#define fftw_comp  sll_comp64
#define fftw_int   sll_int32

#define fftw_plan_dft_1d dfftw_plan_dft_1d
#define fftw_plan_dft_c2r_1d dfftw_plan_dft_c2r_1d
#define fftw_plan_dft_r2c_1d dfftw_plan_dft_r2c_1d
#define fftw_plan_dft_2d dfftw_plan_dft_2d
#define fftw_plan_dft_c2r_2d dfftw_plan_dft_c2r_2d
#define fftw_plan_dft_r2c_2d dfftw_plan_dft_r2c_2d

#define fftw_execute_dft dfftw_execute_dft
#define fftw_execute_dft_r2c dfftw_execute_dft_r2c
#define fftw_execute_dft_c2r dfftw_execute_dft_c2r
#define fftw_execute_r2r dfftw_execute_r2r

#define fftw_destroy_plan    dfftw_destroy_plan

#endif

#endif
