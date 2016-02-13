! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/doc/build/html/doxygen/html/defgroup memorys.html 
! The following lines will be read by doxygen to generate documentation:


!> @defgroup memory sll_memory 
!> @author Edwin Chacon-Golcher
!> @brief 
!> Implements the error testing function and allocation related functionalities.
!> @details
!> Provides the header file that exposes the allocation,
!> initialization and clearing macros, and other related functionalities.
!> 
!> The following directive will make available macros:
!> @code
!> #include "sll_m_memory.h" 
!> @endcode
!>
!> - SLL_ALLOCATE
!> - SLL_DEALLOCATE
!> - SLL_INIT_ARRAY
!> - SLL_CLEAR_ALLOCATE
!>
!> and will provide any additional functionality included in the 
!> sll_m_memory.F90 module. This directive should be put in the 'use' section 
!> of the subprogram that wants to use the allocator.
!>
!> Examples: 
!>
!> To allocate the simple array 'a' of 1000 elements:
!>@code
!> SLL_ALLOCATE(a(1000),err)
!>@endcode
!> To allocate the 2D array 'a' with non-default indexing:
!>@code
!> SLL_ALLOCATE(a(5:55,1:20),err),
!>@endcode
!> etc. The syntax comes directly from the native Fortran allocate function.
!>
!> The deallocation macro also requires the integer error variable:
!>@code
!> SLL_DEALLOCATE(a, err)
!>@endcode
!> To initialize an array:
!>@code
!> SLL_INIT_ARRAY(a,value)
!>@endcode 
!> To allocate and initialize an array to 0:
!>@code
!> SLL_CLEAR_ALLOCATE(a(1:1000,4:500), err)
!>@endcode
!> Any file that uses this include file will need to be compiled such that
!> there is no upper limit to the length of a source code line. In gfortran,
!> for instance, this behavior is obtained with the flag  
!> -ffree-line-length-none
!>
