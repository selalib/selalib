# Control threads, processors and the parallel environment.
# They have C linkage, and do not throw exceptions.
control_routines = { \
'omp_get_active_level',  # Number of active parallel regions
'omp_get_ancestor_thread_num',  # Ancestor thread ID
'omp_get_cancellation',  # Whether cancellation support is enabled
'omp_get_default_device',  # Get the default device for target regions
'omp_get_dynamic',  # Dynamic teams setting
'omp_get_level',  # Number of parallel regions
'omp_get_max_active_levels',  # Maximum number of active regions
'omp_get_max_task_priority',  # Maximum task priority value that can be set
'omp_get_max_threads',  # Maximum number of threads of parallel region
'omp_get_nested',  # Nested parallel regions
'omp_get_num_devices',  # Number of target devices
'omp_get_num_procs',  # Number of processors online
'omp_get_num_teams',  # Number of teams
'omp_get_num_threads',  # Size of the active team
'omp_get_proc_bind',  # Whether theads may be moved between CPUs
'omp_get_schedule',  # Obtain the runtime scheduling method
'omp_get_team_num',  # Get team number
'omp_get_team_size',  # Number of threads in a team
'omp_get_thread_limit',  # Maximum number of threads
'omp_get_thread_num',  # Current thread ID
'omp_in_parallel',  # Whether a parallel region is active
'omp_in_final',  # Whether in final or included task region
'omp_is_initial_device',  # Whether executing on the host device
'omp_set_default_device',  # Set the default device for target regions
'omp_set_dynamic',  # Enable/disable dynamic teams
'omp_set_max_active_levels',  # Limits the number of active parallel regions
'omp_set_nested',  # Enable/disable nested parallel regions
'omp_set_num_threads',  # Set upper team size limit
'omp_set_schedule',  # Set the runtime scheduling method 
}

# Initialize, set, test, unset and destroy simple and nested locks.
lock_routines = { \
'omp_init_lock',  # Initialize simple lock
'omp_set_lock',  # Wait for and set simple lock
'omp_test_lock',  # Test and set simple lock if available
'omp_unset_lock',  # Unset simple lock
'omp_destroy_lock',  # Destroy simple lock
'omp_init_nest_lock',  # Initialize nested lock
'omp_set_nest_lock',  # Wait for and set simple lock
'omp_test_nest_lock',  # Test and set nested lock if available
'omp_unset_nest_lock',  # Unset nested lock
'omp_destroy_nest_lock',  # Destroy nested lock 
}

# Portable, thread-based, wall clock timer.
timer_routines = { \
'omp_get_wtick',  # Get timer precision.
'omp_get_wtime',  # Elapsed wall clock time. 
}

all_routines = set.union( control_routines, lock_routines, timer_routines )
