module sll_xdmf_io

   use sll_m_io_utilities, only: &
      sll_s_remove_file, &  ! delete an existing file
      sll_s_read_file, &  ! read file into string
      sll_f_check_equal_files, &  ! verify if 2 files are identical
      sll_f_check_empty_file, &  ! verify if a file is totally empty
      sll_s_ints_to_string, &  ! write integer array into a string
      sll_s_split_path             ! split path string into dirname + filename

   use sll_m_xml, only: &
      sll_t_xml_document, &        ! XML document tree, only for output
      sll_t_xml_element            ! most basic object in an XML tree

   use sll_m_xdmf_light_serial, only: &
      sll_t_xdmf_file              ! helper: XML doc creation in XDMF database

#ifndef NOHDF5
   use sll_m_hdf5_serial, only: &
      sll_t_hdf5_serial            ! OO wrapper to Pierre's sll_m_hdf5_io_serial
#endif /* NOHDF5 */

end module sll_xdmf_io
