#ifndef _SLL_IO_H
#define _SLL_IO_H

#define SLL_IO_XDMF    0 
#define SLL_IO_VTK     1
#define SLL_IO_GNUPLOT 2 

use sll_m_xml_io
use sll_m_ascii_io
use sll_m_binary_io
use sll_m_gnuplot
use sll_m_gnuplot_parallel
#ifndef NOHDF5
use sll_m_hdf5_io_parallel
use sll_m_hdf5_io_serial
#endif
use sll_m_xdmf_parallel

#endif
