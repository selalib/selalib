# ~/sllrzg/src/particle_methods/pic_remapped/particle_lbf/CMakeLists.txt
# ======================================================================
# 
# _entry points_ [[particle_lbf_library]]
# 
# [[elisp:(alh-set-keywords)]] ([[file:~/.emacs::alh-set-keywords]])
# emacs-keywords make start=08/06/16

# <<<<particle_lbf_library>>>> [[elisp:(startup-persistent-process "rox-filer")]]

ADD_LIBRARY( sll_particle_lbf_4d_marker STATIC
  sll_m_pic_lbfr_4d_marker.F90)

TARGET_LINK_LIBRARIES( sll_particle_lbf_4d_marker
  sll_working_precision
  sll_memory
  )

ADD_LIBRARY( sll_particle_lbf STATIC
  sll_m_pic_lbfr_4d_group.F90
  sll_m_pic_lbfr_4d_utilities.F90
  )

# Selalib compilation primer: [[https://gitlab.mpcdf.mpg.de/clapp/selalib]]
# Configuring: [[elisp:(compile "mkdir -p ${SELALIB}/build && cd ${SELALIB}/build && cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-4.8 -DMPI_Fortran_COMPILER=mpifort -DCMAKE_VERBOSE_MAKEFILE=TRUE -DOPENMP_ENABLED=OFF ${SELALIB}" t)]]
# Compiling: [[elisp:(compile "cd ${SELALIB}/build && OMPI_FC=/usr/bin/gfortran-4.8 make")]]
#
# Configuration cache: [[selalib:build/CMakeCache.txt]]
# Restart from scratch: [[shell:rm -rf ${SELALIB}/build]]

TARGET_LINK_LIBRARIES( sll_particle_lbf
  sll_working_precision
  sll_memory
  sll_assert
  sll_meshes
  sll_pic_interface
  sll_particle_lbf_4d_marker
  sll_remapped_pic
  sll_remapped_pic_utilities
  sll_initial_density_parameters
  # sll_pic_accumulators
  sll_sparse_grid
  sll_pic_opt2d_accumulators
  sll_timer # needed?
  sll_file_io # needed?
  sobol
  )

# Ctest
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# eval:(flyspell-prog-mode)
# eval:(outline-minor-mode)
# End:
# LocalWords: emacs
