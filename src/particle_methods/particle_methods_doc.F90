! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/doc/build/html/doxygen/html/defgroup simulations.html 
! The following lines will be read by doxygen to generate documentation:








  ! [[shell:header\alh 'particle_methods_doc.F90']] (cf [[file:~/alh/bin/headeralh]])
  ! headeralh authors=- default=0 f90 start=25/06/15

  ! To rebuild the documentation: [[elisp:(compile "cd ${SELALIB}/build && make doc")]].

  ! This file creates the page [[selalib:doc/build/html/doxygen/html/group__particle__methods.html]].
  
  ! The Doxygen documentation index is [[selalib:doc/build/html/doxygen/html/index.html]].  The Doxygen configuration
  ! file is [[selalib:doc/doxygen/Doxyfile]].

  !> @defgroup particle_methods sll_particle_methods
  !> @brief Centralized location for all particle methods
  !! @authors Martin Campos Pinto   - <campos@ann.jussieu.fr>
  !! @authors Antoine Le Hyaric     - <Antoine.Le_Hyaric@upmc.fr>

  ! This "brief" appears in the Doxygen list of "Libraries" at [[selalib:doc/build/html/doxygen/html/modules.html]].
  
  
  !> @details
  !! In this directory we are implementing an abstract class of particle method and two classes which 
  !! extend it: 
  !! - one is a "simple" PIC method based on the one implemented by Sever and Edwin, but with no optimization
  !! - the other one is a remapped PIC method which follows from the BSL and LTPIC approach
  !! These classes (and the corresponding deposition and remapping schemes) are used in a simulation 
  !! implemented in the module sll_m_sim_4d_vp_generic_pic_cartesian

  ! this is directory [[selalib:src/particle_methods]]
  
  ! Local Variables:
  ! mode:F90
  ! ispell-local-dictionary:"british"
  ! coding:utf-8
  ! eval:(flyspell-prog-mode)
  ! eval:(outline-minor-mode)
  ! End:
  ! LocalWords: headeralh Doxygen selalib html Centralized src doxygen Doxyfile
