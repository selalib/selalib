IF(Fortran_COMPILER MATCHES IBM)
   SET(CMAKE_Fortran_COMPILER xlf)
ENDIF(Fortran_COMPILER MATCHES IBM)

#SET(CMAKE_Fortran_FLAGS_DEBUG ${CMAKE_Fortran_FLAGS_RELEASE})

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


IF(Fortran_COMPILER_NAME MATCHES "gfortran")
   set(CMAKE_Fortran_FLAGS_DEBUG "-g -w -cpp -fdefault-real-8 ")
   set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -w -cpp -fdefault-real-8 ")
ELSEIF(Fortran_COMPILER_NAME MATCHES "ifort")
   SET(CMAKE_Fortran_FLAGS_RELEASE "-r8 -nowarn -O3 -xHost -ip -openmp")
   SET(CMAKE_Fortran_FLAGS_DEBUG "-g -r8 -nowarn")
ELSEIF(Fortran_COMPILER MATCHES IBM)
   SET(CMAKE_Fortran_FLAGS_RELEASE "-qfixed -qrealsize=8")
   SET(CMAKE_Fortran_FLAGS_DEBUG "-qfixed -qrealsize=8")
ELSEIF(Fortran_COMPILER_NAME MATCHES "pgf")
   SET(CMAKE_Fortran_FLAGS_DEBUG "-g -r8 ")
   SET(CMAKE_Fortran_FLAGS_RELEASE "-fast -r8 ")
ENDIF()

# New version of FFTPACK (this will substitute old 'dfftpack')
ADD_LIBRARY(fftpack5 STATIC fftpack5.f)

# Every external library should be installed (otherwise 'package' does not work)
INSTALL(TARGETS fftpack5 ARCHIVE DESTINATION lib)

# CTest
IF(BUILD_TESTING)
  ADD_EXECUTABLE(test_fftpack5 fftpack5_prb.f)
  TARGET_LINK_LIBRARIES(test_fftpack5 fftpack5)
ENDIF(BUILD_TESTING)
