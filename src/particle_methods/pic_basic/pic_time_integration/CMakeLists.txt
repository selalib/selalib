SET( PIC_TIME_SRCS 
  sll_m_time_propagator_base.F90
  sll_m_operator_splitting_pic_vp_2d2v.F90 
  sll_m_time_propagator_pic_vm_1d2v_hs.F90
  sll_m_time_propagator_pic_vm_1d2v_momentum.F90
  sll_m_time_propagator_pic_vm_1d2v_zigsub.F90
  sll_m_time_propagator_pic_vm_1d2v_cef.F90
  sll_m_time_propagator_pic_vm_1d2v_boris.F90
  sll_m_time_propagator_pic_vm_1d2v_subcyc.F90
  sll_m_time_propagator_pic_vm_1d2v_subcyci.F90
  sll_m_time_propagator_pic_vm_1d2v_helper.F90
  sll_m_time_propagator_pic_vm_1d2v_disgradE.F90
  sll_m_time_propagator_pic_vm_1d2v_disgradEC.F90
  sll_m_time_propagator_pic_vm_1d2v_disgradEC_sub.F90
  sll_m_time_propagator_pic_vm_1d2v_hs_trafo.F90
  sll_m_time_propagator_pic_vm_1d2v_trafo_helper.F90
  sll_m_time_propagator_pic_vm_1d2v_disgradE_trafo.F90
  sll_m_time_propagator_pic_vm_1d2v_disgradEC_trafo.F90
  sll_m_time_propagator_pic_vm_2d3v_hs.F90
  sll_m_time_propagator_pic_vm_3d3v_hs.F90
  sll_m_time_propagator_pic_vm_3d3v_disgradEC.F90
  sll_m_time_propagator_pic_vm_3d3v_helper.F90
  sll_m_particle_mass_1d_base.F90
  sll_m_linear_operator_particle_mass_1d.F90
  sll_m_linear_operator_particle_mass_cl_1d.F90
  sll_m_linear_operator_particle_mass_smooth_1d.F90
  sll_m_particle_mass_3d_base.F90
  sll_m_linear_operator_particle_mass_3d_diag.F90
  sll_m_linear_operator_particle_mass_3d_od.F90
  sll_m_linear_operator_particle_mass_cl_3d_diag.F90
  sll_m_linear_operator_particle_mass_cl_3d_od.F90
  sll_m_linear_operator_schur_ev_1d.F90
  sll_m_linear_operator_schur_ev_3d.F90
  sll_m_linear_operator_schur_phiv_1d.F90
  sll_m_linear_operator_schur_phiv_3d.F90
  sll_m_linear_operator_ecsim.F90
  sll_m_linear_operator_ecsim_eb.F90
  sll_m_linear_operator_ecsim_ev.F90
  sll_m_time_propagator_pic_vm_1d2v_ecsim.F90
  sll_m_time_propagator_pic_vm_1d2v_ecsim2o.F90  
  sll_m_time_propagator_pic_vm_3d3v_disgradE.F90
  sll_m_time_propagator_pic_vm_3d3v_trafo_helper.F90
  sll_m_time_propagator_pic_vm_3d3v_disgradE_trafo.F90
  sll_m_time_propagator_pic_vm_3d3v_hs_trafo.F90
  sll_m_time_propagator_pic_vm_3d3v_cef_trafo.F90
  sll_m_time_propagator_pic_vm_3d3v_disgradEC_trafo.F90
  sll_m_time_propagator_pic_vm_3d3v_cl_helper.F90
  sll_m_time_propagator_pic_vm_3d3v_disgradE_trunc.F90
  sll_m_time_propagator_pic_vm_3d3v_cef.F90)

SET( PIC_TIME_DEPS
  sll_memory
  sll_assert
  sll_errors
  sll_splitting_methods
  sll_particle_groups
  sll_pic_weights
  sll_particle_mesh_coupling
  sll_pic_field_solvers
  sll_poisson_solvers
  sll_maxwell_solvers
  sll_collective
  sll_mapping
  sll_filter
  sll_linear_solvers)



ADD_LIBRARY(sll_pic_time_integration STATIC ${PIC_TIME_SRCS})
TARGET_LINK_LIBRARIES( sll_pic_time_integration ${PIC_TIME_DEPS})



## Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF()
