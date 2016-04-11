# BSL-LT-PIC tests cmakefile
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
# emacs-keywords authors="MCP and ALH" brief="BSL-LT-PIC tests cmakefile" make selalib start=02/11/15

# <<<<test_bsl_lt_pic_4d>>>> [[file:test_bsl_lt_pic_4d.F90]]

ADD_EXECUTABLE( test_bsl_lt_pic_4d test_bsl_lt_pic_4d.F90 )
TARGET_LINK_LIBRARIES( test_bsl_lt_pic_4d
  sll_timer
  sll_assert
  sll_remapped_pic
  sll_bsl_lt_pic
  )

ADD_TEST( NAME bsl_lt_pic_4d COMMAND test_bsl_lt_pic_4d )
SET_TESTS_PROPERTIES( bsl_lt_pic_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

# <<test_4dsg_for_remapped_pic>> [[file:test_4dsg_for_remapped_pic.F90]]
# Configuration from [[selalib:src/simulations/parallel/pic_vp_2d2v_cart_remapped/README]]
# [[elisp:(compile "cd ${SELALIB}/build && make test_4dsg_for_remapped_pic && ${SELALIB}/build/bin/test_4dsg_for_remapped_pic")]]

ADD_EXECUTABLE(test_4dsg_for_remapped_pic test_4dsg_for_remapped_pic.F90)
TARGET_LINK_LIBRARIES(test_4dsg_for_remapped_pic sll_sparse_grid sll_interpolators sll_periodic_interpolation)
ADD_TEST(NAME 4dsg_for_remapped_pic  COMMAND test_4dsg_for_remapped_pic)
SET_TESTS_PROPERTIES(4dsg_for_remapped_pic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# eval:(flyspell-prog-mode)
# eval:(outline-minor-mode)
# End:
# LocalWords: emacs
