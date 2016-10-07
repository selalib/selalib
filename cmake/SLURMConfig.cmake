if( NOT POE_FOUND )

  find_program( SLURM_EXECUTABLE
    NAMES srun
    DOC "Open-source workload manager." )

  if( SLURM_EXECUTABLE )
    set( SLURM_FOUND "YES" )
  endif()

  mark_as_advanced( SLURM_FOUND SLURM_EXECUTABLE )

  if( SLURM_FOUND )
    set( MPIEXEC ${SLURM_EXECUTABLE} )
    set( MPIEXEC_NUMPROC_FLAG "-n" )
  endif( SLURM_FOUND )

endif( NOT POE_FOUND )
