# PIC-LBFR tests cmakefile
# ==========================
# 
# MCP and ALH
#   This code SeLaLib (for Semi-Lagrangian-Library) 
#   is a parallel library for simulating the plasma turbulence 
#   in a tokamak.
# 
#   This software is governed by the CeCILL-B license 
#   under French law and abiding by the rules of distribution 
#   of free software.  You can  use, modify and redistribute 
#   the software under the terms of the CeCILL-B license as 
#   circulated by CEA, CNRS and INRIA at the following URL
#   "http://www.cecill.info". 
# 
# [[elisp:(alh-set-keywords)]] ([[file:~/.emacs::alh-set-keywords]])
# emacs-keywords authors="MCP and ALH" brief="PIC-LBFR tests cmakefile" make selalib start=02/11/15

# <<<<test_pic_lbfr_4d>>>> [[file:test_pic_lbfr_4d.F90]]

ADD_EXECUTABLE( test_pic_lbfr_4d test_pic_lbfr_4d.F90 )
TARGET_LINK_LIBRARIES( test_pic_lbfr_4d
  sll_timer
  sll_assert
  sll_pic_lbfr
  sll_initial_density_parameters
  sll_particle_sampling_interface
  sll_remapped_pic
  )

ADD_TEST( NAME pic_lbfr_4d COMMAND test_pic_lbfr_4d )
SET_TESTS_PROPERTIES( pic_lbfr_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

# <<test_4dsg_for_pic_lbfr>> [[file:test_4dsg_for_pic_lbfr.F90]]
# Configuration from [[selalib:src/simulations/parallel/pic_vp_2d2v_cart_remapped/README]]
# [[elisp:(compile "cd ${SELALIB}/build && make test_4dsg_for_pic_lbfr && ${SELALIB}/build/bin/test_4dsg_for_pic_lbfr")]]

ADD_EXECUTABLE(test_4dsg_for_pic_lbfr test_4dsg_for_pic_lbfr.F90)
TARGET_LINK_LIBRARIES(test_4dsg_for_pic_lbfr sll_sparse_grid sll_interpolators sll_periodic_interpolation)
ADD_TEST(NAME 4dsg_for_pic_lbfr  COMMAND test_4dsg_for_pic_lbfr)
SET_TESTS_PROPERTIES(4dsg_for_pic_lbfr PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# eval:(flyspell-prog-mode)
# eval:(outline-minor-mode)
# End:
# LocalWords: emacs
