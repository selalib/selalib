!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: sll_fft
!
! DESCRIPTION:
!> @file fft_documentation.F90
!> @namespace sll_fft
!> @author EDWIN C. GOLCHER & SAMUEL DE SANTIS
!> @brief Interface around fftpack, fftw and the interne selalib fft.
!> @details 
!>
!> 
!>
!> \section how How to use sll_fft module?
!>
!> The first thing is to add the line \code use sll_fft \endcode
!> The sll_fft module can use internal or external librarie 
!> 
!> 1. Declare a fft plan
!> \code type(sll_fft_plan), pointer :: p \endcode
!> 2. Initialize the plan
!> \code p => fft_new_plan(size,in,out,direction,flags) \endcode
!> The arrays in and out can be real and/or complex, 1d or 2d. The size is only a power of two (radix-2).
!> \warning For complex to real and real to complex transform, there is no direction flag.
!>          \code p => fft_new_plan(size,in,out,flags) \endcode
!>
!> \a direction can take two values : FFT_FORWARD and FFT_INVERSE
!>
!> \a flags optional argument that can be : FFT_NORMALIZE
!>                                FFT_ONLY_FIRST_DIRECTION, FFT_ONLY_SECOND_DIRECTION (2d case only)
!> You can combine flags with '+'.
!>
!> 3. Execute the plan
!> \code call fft_apply_plan(p,in,out) \endcode
!> 4. Delete the plan
!> \code call fft_delete_plan(p) \endcode
!>
!>
!> \section sum Summary:
!>
!> 1D
!> <table border="1">
!> <tr>
!> <th> size problem </th>
!> <th> type of in </th>
!> <th> type of out </th>
!> <th> size of in </th>
!> <th> size of out </th>
!> <th> direction </th>
!> <th> flags </th>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> real </td>
!> <td> real </td>
!> <td> n </td>
!> <td> n </td>
!> <td> FFT_FORWARD <br /> FFT_INVERSE </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> complex </td>
!> <td> complex </td>
!> <td> n </td>
!> <td> n </td>
!> <td> FFT_FORWARD <br /> FFT_INVERSE </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> real </td>
!> <td> complex </td>
!> <td> n </td>
!> <td> n/2 </td>
!> <td> ----- </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> complex </td>
!> <td> real </td>
!> <td> n/2 </td>
!> <td> n </td>
!> <td> ----- </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> </table>
!>
!>
!> 2D
!> <table border="1">
!> <tr>
!> <th> size problem </th>
!> <th> type of in </th>
!> <th> type of out </th>
!> <th> size of in </th>
!> <th> size of out </th>
!> <th> direction present </th>
!> <th> flags </th>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> real </td>
!> <td> real </td>
!> <td> n,m </td>
!> <td> n,m </td>
!> <td> FFT_FORWARD <br /> FFT_INVERSE </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> complex </td>
!> <td> complex </td>
!> <td> n,m </td>
!> <td> n,m </td>
!> <td> FFT_FORWARD <br /> FFT_INVERSE </td>
!> <td> FFT_NORMALIZE <br /> FFT_ONLY_FIRST_DIRECTION <br /> FFT_ONLY_SECOND_DIRECTION </td>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> real </td>
!> <td> complex </td>
!> <td> n,m </td>
!> <td> n/2,m </td>
!> <td> ----- </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> complex </td>
!> <td> real </td>
!> <td> n/2,m </td>
!> <td> n,m </td>
!> <td> ----- </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> </table>
!>
!> \section example Examples:
!>
!> In-place transform
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_comp64, dimension(0,n-1) :: in
!> type(sll_fft_plan), pointer  :: p
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,in,FFT_FORWARD,FFT_NORMALIZE)
!> call fft_apply_plan(p,in,in)
!> call fft_delete_plan(p)
!> \endcode
!>
!> Two-dimensional transform  
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_int32, parameter :: m = 2**3
!> sll_comp64, dimension(n/2,m) :: in
!> sll_real64, dimension(n,m) :: out
!> type(sll_fft_plan), pointer  :: p
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_INVERSE)
!> call fft_apply_plan(p,in,out)
!> call fft_delete_plan(p)
!> \endcode
!>
!> Transform in one direction
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_int32, parameter :: m = 2**3
!> sll_comp64, dimension(n,m) :: in
!> sll_comp64, dimension(n,m) :: out
!> type(sll_fft_plan), pointer  :: p
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_FORWARD,FFT_ONLY_FIRST_DIRECTION)
!> call fft_apply_plan(p,in,out)
!> call fft_delete_plan(p)
!> \endcode
!>
!> \section acc Access the mode
!> 
!> To get the value of a mode call the function fft_get_mode(plan,data,num_mode)
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_comp64, dimension(0,n-1) :: in
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: mode
!> sll_int32 :: k = 3
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,in,FFT_FORWARD,FFT_NORMALIZE)
!> call fft_apply_plan(p,in,in)
!> mode = fft_get_mode(plan,in,k)
!> call fft_delete_plan(p)
!> \endcode
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_int32, parameter :: m = 2**3
!> sll_comp64, dimension(n/2,m) :: in
!> sll_real64, dimension(n,m) :: out
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: mode
!> sll_int32 :: k=4, l=2
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_INVERSE)
!> call fft_apply_plan(p,in,out)
!> mode = fft_get_mode(plan,out,k,l)
!> call fft_delete_plan(p)
!> \endcode
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_real64, dimension(0,n-1) :: in
!> sll_real64, dimension(0,n-1) :: out
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: mode
!> sll_int32 :: k = 3
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_FORWARD,FFT_NORMALIZE)
!> call fft_apply_plan(p,in,out)
!> mode = fft_get_mode(plan,out,k)
!> call fft_delete_plan(p)
!> \endcode
!>
!> To set a mode call the subroutine fft_set_mode(plan,data,new_value,num_mode)
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_comp64, dimension(0,n-1) :: in
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: new_value
!> sll_int32 :: k = 3
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,in,FFT_FORWARD,FFT_NORMALIZE)
!> call fft_apply_plan(p,in,in)
!> new_value = complex(5.0_f64,3.2_f64)
!> call fft_set_mode(plan,in,new_value,k)
!> call fft_delete_plan(p)
!> \endcode
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_int32, parameter :: m = 2**3
!> sll_comp64, dimension(n/2,m) :: in
!> sll_real64, dimension(n,m) :: out
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: new_value
!> sll_int32 :: k=4, l=2
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_INVERSE)
!> call fft_apply_plan(p,in,out)
!> new_value = complex(5.0_f64,3.2_f64)
!> call fft_set_mode(plan,out,new_value,k,l)
!> call fft_delete_plan(p)
!> \endcode
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_real64, dimension(0,n-1) :: in
!> sll_real64, dimension(0,n-1) :: out
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: new_value
!> sll_int32 :: k = 0
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_FORWARD,FFT_NORMALIZE)
!> call fft_apply_plan(p,in,out)
!> new_value = complex(5.0_f64,3.2_f64)
!> call fft_set_mode(plan,out,new_value,k)
!> call fft_delete_plan(p)
!> \endcode
!>
!> \section what What sll_fft really computes
!>
!> The forward (FFT_FORWARD) DFT of a 1d complex array x of size n computes an array X, where:
!>
!> \f[ X_k = \sum_{i=0}^{n-1} x_i e^{-2\pi i j k/n}. \f]
!>
!> The backward (FFT_INVERSE) DFT computes:
!>
!> \f[ x_i = \sum_{k=0}^{n-1} X_k e^{2\pi k j i/n}. \f]
!>
! For the real transform, we have
! \f$ (x_0,x_1,\dots,x_{n-1}) \rightarrow
!     (r_0,r_{n/2},r_1,i_1,\dots,r_{n/2-1},i_{n/2-1})\f$
! which must be interpreted as the complex array
! \f[ \begin{pmatrix} r_0 &,& 0
!                     \\ r_1 &,& i_1
!                     \\ \vdots  & & \vdots 
!                     \\ r_{n/2-1} &,& i_{n/2-1}
!                     \\ r_{n/2} &,& 0
!                     \\ r_{n/2-1} &,& -i_{n/2-1}
!                     \\ \vdots    & & \vdots
!                     \\ r_1 &,& -i_1 
! \end{pmatrix}\f] 
! \warning Note that ffw use \f$(r_0,r_1,\dots,r_{n/2-1},r_{n/2},i_{n/2-1},\dots,i_1)\f$
!          convention whereas fftpack use \f$(r_0,r_1,i_1,\dots,r_{n/2-1},i_{n/2-1},r_{n/2})\f$
! 
!
! By example if, the input data is \f$(x_0,x_1,x_2,x_3)\f$ the output is \f$(X_0,X_2,X_1,X_3)\f$.
! Thus, sll_get_index(1) returns 2 (cause data[1]=X_2) and sll_get_mode(1) returns X_1.
!------------------------------------------------------------------------------
