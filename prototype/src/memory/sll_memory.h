#ifndef _sll_memory_h_
#define _sll_memory_h_

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

  ! **************************************************************************
  ! Tentative Selalib's basic memory allocator.
  !
  ! This service is provided by two files:
  ! - sll_memory.h:   the present header file that exposes the allocation,
  !                   initialization and clearing macros, and
  ! - sll_memory.F90: which implements the error testing function and
  !                   other related functionalities.
  ! 
  ! Usage:
  ! The directive
  ! #include "sll_memory.h" will make available macros:
  ! SLL_ALLOCATE
  ! SLL_DEALLOCATE
  ! SLL_INIT_ARRAY
  ! SLL_CLEAR_ALLOCATE
  ! and will provide any additional functionality included in the 
  ! sll_memory.F90 module. This directive should be put in the 'use' section 
  ! of the subprogram that wants to use the allocator.
  !
  ! Examples: 
  !
  ! To allocate the simple array 'a' of 1000 elements:
  !
  ! SLL_ALLOCATE(a(1000),err)
  !
  ! To allocate the 2D array 'a' with non-default indexing:
  !
  ! SLL_ALLOCATE(a(5:55,1:20),err),
  !
  ! etc. The syntax comes directly from the native Fortran allocate function.
  !
  ! The deallocation macro also requires the integer error variable:
  !
  ! SLL_DEALLOCATE(a, err)
  !
  ! To initialize an array:
  !
  ! SLL_INIT_ARRAY(a,value)
  ! 
  ! To allocate and initialize an array to 0:
  !
  ! SLL_CLEAR_ALLOCATE(a(1:1000,4:500), err)
  !
  ! Any file that uses this include file will need to be compiled such that
  ! there is no upper limit to the length of a source code line. In gfortran,
  ! for instance, this behavior is obtained with the flag  
  ! -ffree-line-length-none
  !
  ! *************************************************************************
 
use sll_memory

#define SLL_ALLOCATE(array_name_and_lims, error_var)   \
  allocate(array_name_and_lims, stat=error_var);      \
  call test_error_code(error_var, 'Memory allocation Failure.', __FILE__, \
  __LINE__);

#define SLL_DEALLOCATE(ptr, error_var)  \
  deallocate(ptr, stat=error_var);      \
  call test_error_code(error_var, 'Error in memory deallocation.', __FILE__,\
  __LINE__); \
  nullify(ptr);

#define SLL_DEALLOCATE_ARRAY(array, error_var) \
  deallocate(array, stat=error_var);           \
  call test_error_code(error_var, 'Failed array deallocation: ', __FILE__, \
  __LINE__ ); 

#define SLL_INIT_ARRAY(arry,val) arry = val;

#define SLL_CLEAR_ALLOCATE(arry_name_and_lims, error_var) \
  SLL_ALLOCATE(arry_name_and_lims, error_var)             \
  SLL_INIT_ARRAY(arry_name_and_lims, 0.0) 

  ! **************************************************************************
  ! IMPLEMENTATION NOTES FOR sll_memory.h:
  !
  ! This is an exploratory implementation of Selalib's basic memory allocator.
  ! There are only two things that we ask of the memory allocator:
  ! 1. To allocate the memory requested, and
  ! 2. In case of failure, to stop the program execution with a message 
  !    informing where was the call made.
  ! Thus, we are leaving the initialization of the memory or any other
  ! behavior, like more detailed logging, to other facilities built on top of
  ! this basic allocator. 
  !
  ! SLL_ALLOCATE has virtually the same syntax as the native Fortran 
  ! 'allocate' function, except for the mandatory integer error variable that 
  ! must be passed as an argument. The caller is responsible for providing 
  ! this variable, as with the native function.
  !
  ! The macro is kind of a wrapper around 'allocate' with error checking.
  ! It is written in a way that can be dealt with by the Fortran
  ! preprocessor. The success of the macro relies on the preprocessor 
  ! requirement that argument parentheses balance, thus, a comma within
  ! those parentheses will not end the argument. This is why a call like
  ! SLL_ALLOCATE(array(1:100,1:200),err) has two commas, but only the 
  ! second acts as an argument separator. This gives a convenient
  ! flexibility in this case.
  !
  ! Possible pitfall: Fortran has a limit on the length of a line of 132
  !                   characters. The length of some of the above macros
  !                   will exceed this limit. Any file that includes this
  !                   header will need to be compiled with a flag to 
  !                   allow for much longer lines. In gfortran the flag is:
  !                   -ffree-line-length-none
  !
  ! In case that this way of allocating memory is found not satisfactory,
  ! the file sll_memory.F90 has some elements of an alternative implementation,
  ! which also uses macros to handle most of the redundancies. However, these 
  ! macros, while still manageable by fpp, will also yield very long lines 
  ! that may confuse a Fortran compiler if one does not use flags like:
  ! -ffree-line-length-none
  ! The alternative implementation consists of allocating
  ! subroutines, each specialized for a particular type and dimension of an
  ! array. In order to pass along the __FILE__ and __LINE__ information to
  ! the functions conveniently, some macro call has to be made as well, unless
  ! one goes for the truly inconvenient and ugly solution of passing those
  ! manually.
  !
  ! **************************************************************************

#endif
 
