# XML file creation (using xml type)
add_executable(test_xml unit_test_xml.F90)
target_link_libraries(test_xml sll_xdmf_io)

set(XML_REFERENCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ref1.xml)
add_test(NAME xml COMMAND test_xml ${XML_REFERENCE_FILE})
set_tests_properties(xml PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 1)
