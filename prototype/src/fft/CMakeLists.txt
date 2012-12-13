set(DEPENDENCIES
	sll_utilities
	sll_timer
	sll_constants
	sll_memory
	sll_assertion
	sll_working_precision
)

IF(FFTW_FOUND AND FFTW_ENABLED)
  add_library(fftw_module STATIC fftw_module.F90)
  LIST(APPEND DEPENDENCIES fftw_module)
  target_link_libraries(fftw_module ${FFTW_LIBRARIES})
  include_directories(${FFTW_INCLUDE_DIRS})
ENDIF()

IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTPACK")
  set(FFT_SRC sll_fft_mod_through_fftpack.F90)
  set(DEPENDENCIES ${DEPENDENCIES} fftpack_module)
  add_definitions(-D_DEFAULTFFTLIB=FFTPACK_MOD)
ELSEIF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW")
  set(FFT_SRC sll_fft_mod_through_legacy_fftw.F90)
  add_definitions(-D_DEFAULTFFTLIB=FFTW_MOD)
ELSE()
  set(FFT_SRC sll_fft_selalib_implementation.F90)
  add_definitions(-D_DEFAULTFFTLIB=SLLFFT_MOD)
ENDIF()

add_library(sll_fft STATIC
	${FFT_SRC}
	fft_utils.F90

)
add_dependencies(sll_fft ${DEPENDENCIES})

IF(FFTW_FOUND AND FFTW_ENABLED)
	target_link_libraries(sll_fft fftw_module)
ENDIF()
IF(FFTPACK_ENABLED)
	target_link_libraries(sll_fft fftpack_module)
ENDIF()

add_executable(test_fft unit_test.F90)
target_link_libraries(test_fft
	sll_fft
	${DEPENDENCIES}
)


IF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTPACK")
add_executable( test_macro_fftpack test_macro_fftpack.F90)
target_link_libraries( test_macro_fftpack
        sll_fft
	sll_utilities 
	sll_working_precision
	sll_constants
	sll_low_level_file_io
  	sll_utilities
	sll_memory
        sll_assertion
	fftpack_module
	${LAPACK_LIBRARIES}
        sll_timer
  )
ELSEIF(${FFT_DEFAULT_LIBRARY} STREQUAL "FFTW")
add_executable( test_macro_fftw test_macro_fftw.F90)
target_link_libraries( test_macro_fftw
	sll_constants
        sll_fft
	sll_constants
	sll_low_level_file_io
  	sll_utilities
        sll_assertion
	sll_memory
	fftpack_module
	${LAPACK_LIBRARIES}
        sll_timer
  )
ELSE()
add_executable( test_macro_fftselalib test_macro_fftselalib.F90)
target_link_libraries( test_macro_fftselalib
	sll_constants
        sll_fft
	sll_constants
	sll_low_level_file_io
  	sll_utilities
	sll_memory
        sll_assertion
	sll_constants
	fftpack_module
	${LAPACK_LIBRARIES}
        sll_timer
  )
ENDIF()

