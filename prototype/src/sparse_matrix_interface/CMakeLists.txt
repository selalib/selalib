##########################################################################
# spm directory
##########################################################################
set(spm_files getdia.F90 getu.F90 getl.F90 lsol.F90 remove_zeros.F90  basic_cg.F90 gauss_seidel.F90 spm_def.F90 spm_tools.F90 spm_basic.F90)
set(spm_mods basic_cg_mod.mod gauss_seidel_mod.mod spm_def.mod spm.mod)
#add_library(spm STATIC ${spm_files})
add_library(spm SHARED ${spm_files})
install(TARGETS spm DESTINATION lib)
install(FILES ${spm_mods} DESTINATION include)
##########################################################################
