!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: sll_m_fft
!
! DESCRIPTION:
!> @defgroup fft sll_fft
!> @author Edwin Chacon-Golcher, Samuel De Santis, Pierre Navaro, Katharina Kormann
!> @brief Interface around fftpack, fftw and the selalib fft.
!> @details 
!> 
!>
!> \section how How to use sll_m_fft module?
!>
!> The first thing is to add the line \code use sll_m_fft \endcode
!> The sll_m_fft module can use internal or external librarie s
!> 
!> 1. Declare a fft plan
!> \code type(sll_t_fft) :: p \endcode
!> 2. Initialize the plan
!> \code call sll_s_fft_init_c2c_1d( p, size, in, out, direction, normalized, aligned, optimization ) \endcode
!>
!> The arrays in and out can be real and/or complex, 1d or 2d. Change c2c_1d accordingly.
!> \warning For complex to real and real to complex transform, there is no direction flag.
!>          \code call sll_s_fft_init_r2c_1d( p, size, in, out, normalized, aligned, optimization ) \endcode
!>
!> \a direction can take two values : FFT_FORWARD and FFT_BACKWARD
!>
!> \a normalized is a boolean optional argument: If true, the transformed data is normalized by the (product) of the element length.  [default: false]
!> \a aligned is a boolean optional argument (only used by FFTW): If true, FFTW assumes the the data is aligned (use \a fft_alloc for in/out data allocation in this case).  [default: false]         
!> \a optimization is an optional argument (only used by FFTW): With this argument, you can set how rigorous the initialization of the planner will be. It can take the values FFT_ESTIMATE, FFT_WISDOM_ONLY, FFT_MEASURE, FFT_PATIENT, FFT_EXHAUSTIVE. Note that \a in and \a out are overwritten for the last three options. FFT_WISDOM_ONLY only works if wisdom is available. [default: FFT_ESTIMATE]                      
!>
!> 3. Execute the plan
!> \code call sll_s_fft_exec_c2c_1d( p, in, out ) \endcode
!> 4. Delete the plan
!> \code call sll_s_fft_free( p ) \endcode
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
!> <th> normalized </th>
!> <th> aligned </th>
!> <th> optimization </th>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> real </td>
!> <td> real </td>
!> <td> n </td>
!> <td> n </td>
!> <td> FFT_FORWARD <br /> FFT_BACKWARD </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> FFT_ESTIMATE <br /> FFT_WISDOM_ONLY <br />  FFT_MEASURE <br /> FFT_PATIENT <br /> FFT_EXHAUSTIVE  </td>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> complex </td>
!> <td> complex </td>
!> <td> n </td>
!> <td> n </td>
!> <td> FFT_FORWARD <br /> FFT_BACKWARD </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> FFT_ESTIMATE <br /> FFT_WISDOM_ONLY <br />  FFT_MEASURE <br /> FFT_PATIENT <br /> FFT_EXHAUSTIVE  </td>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> real </td>
!> <td> complex </td>
!> <td> n </td>
!> <td> n/2+1 </td>
!> <td> ----- </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> FFT_ESTIMATE <br /> FFT_WISDOM_ONLY <br />  FFT_MEASURE <br /> FFT_PATIENT <br /> FFT_EXHAUSTIVE  </td>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> complex </td>
!> <td> real </td>
!> <td> n/2+1 </td>
!> <td> n </td>
!> <td> ----- </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> FFT_ESTIMATE <br /> FFT_WISDOM_ONLY <br />  FFT_MEASURE <br /> FFT_PATIENT <br /> FFT_EXHAUSTIVE  </td>
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
!> <th> normalized </th>
!> <th> aligned </th>
!> <th> optimization </th>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> real </td>
!> <td> real </td>
!> <td> n,m </td>
!> <td> n,m </td>
!> <td> FFT_FORWARD <br /> FFT_BACKWARD </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> FFT_ESTIMATE <br /> FFT_WISDOM_ONLY <br />  FFT_MEASURE <br /> FFT_PATIENT <br /> FFT_EXHAUSTIVE  </td>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> complex </td>
!> <td> complex </td>
!> <td> n,m </td>
!> <td> n,m </td>
!> <td> FFT_FORWARD <br /> FFT_BACKWARD </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> FFT_ESTIMATE <br /> FFT_WISDOM_ONLY <br />  FFT_MEASURE <br /> FFT_PATIENT <br /> FFT_EXHAUSTIVE  </td>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> real </td>
!> <td> complex </td>
!> <td> n,m </td>
!> <td> n/2,m </td>
!> <td> ----- </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> FFT_ESTIMATE <br /> FFT_WISDOM_ONLY <br />  FFT_MEASURE <br /> FFT_PATIENT <br /> FFT_EXHAUSTIVE  </td>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> complex </td>
!> <td> real </td>
!> <td> n/2,m </td>
!> <td> n,m </td>
!> <td> ----- </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> .TRUE. <br /> .FALSE. </td>
!> <td> FFT_ESTIMATE <br /> FFT_WISDOM_ONLY <br />  FFT_MEASURE <br /> FFT_PATIENT <br /> FFT_EXHAUSTIVE  </td>
!> </tr>
!> </table>
!>
!> \section example Examples:
!>
!> In-place transform
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_comp64, dimension(0,n-1) :: in
!> type(sll_t_fft)  :: p
!>
!> !** INIT DATA **
!>
!> call sll_s_fft_init_c2c_1d( p, n, in, in, FFT_FORWARD, normalized = .TRUE. )
!> call sll_s_fft_exec_c2c_1d( p, in, in )
!> call sll_s_fft_free( p )
!> \endcode
!>
!> Two-dimensional transform  
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_int32, parameter :: m = 2**3
!> sll_comp64, dimension(n/2,m) :: in
!> sll_real64, dimension(n,m) :: out
!> type(sll_t_fft)  :: p
!>
!> !** INIT DATA **
!>
!> call sll_s_fft_init_c2r_2d( n, in, out )
!> call sll_s_fft_exec_c2r_2d( p, in, out )
!> call sll_s_fft_free( p )
!> \endcode
!>
!>
!> \section what What sll_m_fft really computes
!>
!> The forward (FFT_FORWARD) DFT of a 1d complex array x of size n computes an array X, where:
!>
!> \f[ X_k = \sum_{i=0}^{n-1} x_i e^{-2\pi i j k/n}. \f]
!>
!> The backward (FFT_BACKWARD) DFT computes:
!>
!> \f[ x_i = \sum_{k=0}^{n-1} X_k e^{2\pi k j i/n}. \f]
!>
!> For the real transform, we have
!> \f$ (x_0,x_1,\dots,x_{n-1}) \rightarrow
!>     (r_0,r_{n/2},r_1,i_1,\dots,r_{n/2-1},i_{n/2-1})\f$
!> which must be interpreted as the complex array
!> \f[ \begin{pmatrix} r_0 &,& 0
!>                     \\ r_1 &,& i_1
!>                     \\ \vdots  & & \vdots 
!>                     \\ r_{n/2-1} &,& i_{n/2-1}
!>                     \\ r_{n/2} &,& 0
!>                     \\ r_{n/2-1} &,& -i_{n/2-1}
!>                     \\ \vdots    & & \vdots
!>                     \\ r_1 &,& -i_1 
!> \end{pmatrix}\f] 
!> \warning Note that ffw uses \f$(r_0,r_1,\dots,r_{n/2-1},r_{n/2},i_{n/2-1},\dots,i_1)\f$
!>          convention whereas fftpack uses \f$(r_0,r_1,i_1,\dots,r_{n/2-1},i_{n/2-1},r_{n/2})\f$. Through the interface, FFTW ordering is enforced also for FFTPACK. 
!> For r2c and r2c, the lower half (plus one element) is stored (non-negative frequencies). In higher dimensions, the first dimension is reduced to n/2+1 whereas the others remain n.
!> 
!>
!------------------------------------------------------------------------------
