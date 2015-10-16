#ifndef _assert_h_
#define _assert_h_

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

#if (defined ( GFORTRAN ) || defined ( G95 ) || defined(MPIF90))
# define STRNG(x) "x"
#else
# define STRNG(x) #x
#endif

    ! Note the semicolon when the SLL_ASSERT() macro gets expanded. We need
    ! that this be a part of the macro since we also want to be able to
    ! use this macro call within other macros. If the expansion yields 
    ! nothing (second case) then we don't want dangling semicolons...
#ifdef DEBUG
# define SLL_ASSERT(x) if ( .not. (x) ) then;          \
      call sll_assertion( STRNG(x), __FILE__, __LINE__ ); \
   end if;
#else
# define SLL_ASSERT(x) 
#endif


   ! The intent of the SLL_PRIV macro is to enforce the interface 'walls' for
   ! a given type, preventing unintended access to the given type internals.
   ! However, while debugging it may also be desired to permit direct access.
   ! the SLL_PRIV macro serves as a switchable 'private' declaration. The
   ! tricky part is that it should be used always as the last attribute in a
   ! declaration, since if it were to follow a comma, a compiler would not 
   ! react happily.
#ifdef DEBUG
# define SLL_PRIV
#else
# define SLL_PRIV ,private
#endif



use sll_m_assert


#endif
