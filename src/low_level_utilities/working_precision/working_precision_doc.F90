! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/documentation/build/html/doxygen/html/defgroup working_precisions.html 
! The following lines will be read by doxygen to generate documentation:


!> @defgroup working_precision sll_working_precision 
!> @brief 
!> Define the kind type parameter for intern type data.
!> @author Edwin Chacon-Golcher
!> @details
!>
!> <b> Headers file available </b>
!>  - sll_m_working_precision.h
!>
!> <b> How to use this module: </b> \n
!>
!> If you want to select the kind parameter \a n, you need to write \n
!> \code
!> real(kind=n) :: var1
!> real*n       :: var2 
!> \endcode
!> The two entries \a var1 and \a var2 are equivalents.
!> You can also define the constant like this \a 23.455_n.
!>
!> First, call the module \a sll_woring_precision like that
!> \code #include "sll_m_working_precision.h" \endcode
!> Now, you can use the types :
!> \code
!> _sll_int32  :: i        !integer simple precision
!> _sll_int64  :: N        !integer double precision
!> _sll_real32 :: theta    !real simple precision
!> _sll_real64 :: my_pi    !real double precision
!>
!> my_pi = 3.1415926535897932384626433_f64
!> theta = 2.0*real(N,f32)
!> \endcode
