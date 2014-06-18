IF(Fortran_COMPILER MATCHES IBM)
   SET(CMAKE_Fortran_COMPILER xlf)
ENDIF(Fortran_COMPILER MATCHES IBM)


SET(DFFTPACK_SRCS zfftb.f cfftb1.f zfftf.f cfftf1.f zffti.f cffti1.f
                  dcosqb.f cosqb1.f dcosqf.f cosqf1.f dcosqi.f dcost.f
                  dcosti.f ezfft1.f dzfftb.f dzfftf.f dzffti.f passb.f
                  passb2.f passb3.f passb4.f passb5.f passf.f passf2.f
                  passf3.f passf4.f passf5.f radb2.f radb3.f radb4.f
                  radb5.f radbg.f radf2.f radf3.f radf4.f radf5.f radfg.f
                  dfftb.f rfftb1.f dfftf.f rfftf1.f dffti.f rffti1.f
                  dsinqb.f dsinqf.f dsinqi.f dsint.f sint1.f dsinti.f)

ADD_LIBRARY(dfftpack STATIC ${DFFTPACK_SRCS})
SET(FFTPACK_LIBRARIES dfftpack)
INSTALL(TARGETS dfftpack ARCHIVE DESTINATION lib)
