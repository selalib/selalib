#ifndef _SLL_FFT_H_
#define _SLL_FFT_H_
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


use sll_fft

#define GET_MODE0(mode,data) \
      mode = cmplx(data(1),0.0_f64,kind=f64)
#define SET_MODE0(new_value,data) \
      data(1) = real(new_value,kind=f64)
#define GET_MODE_COMPLEX_1D(mode,array,k) \
	  mode=array(k)
#define SET_MODE_COMPLEX_1D(array,new_value,k) \
	  array(k)=nex_value
#define GET_MODE_COMPLEX_2D(mode,array,k,l) \
	  mode=array(k,l)
#define SET_MODE_COMPLEX_2D(array,new_value,k,l) \
	  array(k,l)=nex_value
#define GET_MODE_COMPLEX_3D(mode,array,k,l,m) \
	  mode=array(k,l,m)
#define SET_MODE_COMPLEX_3D(array,new_value,k,l,m) \
	  array(k,l,m)=nex_value

#ifdef SLLFFT

#define GET_MODE_N_2(mode,data,n) \
        mode = cmplx(data(2),0.0_f64,kind=f64)

#define GET_MODE_GT_N_2(mode,data,k,n) \
        mode = cmplx( data(2*(n-k)+1) , -data(2*(n-k)+2),kind=f64)

#define GET_MODE_LT_N_2(mode,data,k,n) \
        mode = cmplx( data(2*k+1) , data(2*k+2) ,kind=f64)

#define SET_MODE_N_2(new_value,data,n) \
        data(2) = real(new_value,kind=f64)

#define SET_MODE_GT_N_2(new_value,data,k,n) \
        data(2*(n-k)+1) = real(new_value,kind=f64); \
        data(2*(n-k)+2) = -dimag(new_value)

#define SET_MODE_LT_N_2(new_value,data,k,n) \
        data(2*k+1) = real(new_value,kind=f64); \
        data(2*k+2) = dimag(new_value)

#endif

#ifdef FFTPACK

#define GET_MODE_N_2(mode,data,n) \
        mode = cmplx(data(n),0.0_f64,kind=f64)
#define GET_MODE_GT_N_2(mode,data,k,n) \
        mode = cmplx( data(2*(n-k)) , -data(2*(n-k)+1) ,kind=f64)
#define GET_MODE_LT_N_2(mode,data,k,n) \
        mode = cmplx( data(2*k) , data(2*k+1) ,kind=f64)     
#define SET_MODE_N_2(new_value,data,n) \
        data(n) = real(new_value,kind=f64)
#define SET_MODE_GT_N_2(new_value,data,k,n) \
        data(2*(n-k)) = real(new_value,kind=f64); \
        data(2*(n-k)+1) = -dimag(new_value)
#define SET_MODE_LT_N_2(new_value,data,k,n) \
        data(2*k) = real(new_value,kind=f64); \
        data(2*k+1) = dimag(new_value)

#endif

#ifdef FFTW
#define GET_MODE_N_2(mode,data,n) \
        mode = cmplx(data(n/2+1),0.0_f64,kind=f64)
#define GET_MODE_GT_N_2(mode,data,k,n) \
        mode = cmplx( data(n-k+1) , -data(k+1) ,kind=f64)
#define GET_MODE_LT_N_2(mode,data,k,n) \
        mode = cmplx( data(k+1) , data(n-k+1) ,kind=f64)  
#define SET_MODE_N_2(new_value,data,n) \
        data(n/2+1) = real(new_value,kind=f64)
#define SET_MODE_GT_N_2(new_value,data,k,n) \
        data(n-k+1) = real(new_value,kind=f64); \
        data(k+1) = -dimag(new_value)
#define SET_MODE_LT_N_2(new_value,data,k,n) \
        data(k+1) = real(new_value,kind=f64); \
        data(n-k+1) = dimag(new_value)
#endif

#endif
