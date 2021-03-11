set(XDMF_SRC sll_xdmf_io.F90 sll_m_xml.F90 sll_m_xdmf_light_serial.F90
             sll_m_io_utilities.F90)

# Object-oriented wrapper around Pierre's "sll_hdf5_io_serial"
set(XDMF_SRC ${XDMF_SRC} sll_m_hdf5_serial.F90)
set(XDMF_DEP sll_file_io)

add_library(sll_xdmf_io STATIC ${XDMF_SRC})
target_link_libraries(sll_xdmf_io ${XDMF_DEP})

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
