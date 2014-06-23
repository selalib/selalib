#!/bin/bash
#@ class            = clallmds
#@ job_name         = temp
#@ total_tasks      = 16
#@ node             = 2
#@ node_usage       = not_shared
#@ wall_clock_limit = 00:30:00
#@ job_type         = mpich
#@ environment      = COPY_ALL 
#@ queue
