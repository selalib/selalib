!> @defgroup fftpack5 external 
!> @brief 
!> FFTPACK5 is a FORTRAN77 library which computes the Fast Fourier Transform, by Paul Swarztrauber and Dick Valent.
!> @details
!> Note: An apparent indexing problem in the 2D complex codes 
!> CFFT2B/CFFT2F/CFFT2I and ZFFT2B/ZFFT2F/ZFFT2I was reported to the authors on 10 May 2010. A fix has been promised.
!> 
!> Special features include:
!> 
!> real or complex data can be handled;
!> separate routines for forward analysis (data => Fourier coefficients) and backward analysis (Fourier coefficients => data);
!> sine and cosine transform routines;
!> quarter wave sine and cosine transform routines;
!> the amount of data is NOT required to be a power of 2.
!> Routines in the library come in groups of three:
!> 
!> an initialization routine;
!> the forward computational routine (data to FFT coefficients);
!> the backward computational routine (FFT coefficients to data).
!> 
!>
!> <b> How to use it </b>
!> - Link with   <code>-lfftpack5</code>
