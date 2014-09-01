# @ job_name = job_simple
# @ job_type = BLUEGENE
# Fichier sortie standard du travail
# @ output = $(job_name).$(jobid)
# Fichier erreur standard du travail
# @ error = $(output)
# Temps elapsed maximum demande
# @ wall_clock_limit = 1:00:00
# Taille bloc d'execution
# @ bg_size = 64
# @ queue

# Copy executable and input file to TMPDIR
# Warning: if you need to transfer important volumes
# of data, please use a multi-step job
cp my_code $TMPDIR
cp data.in $TMPDIR
cd $TMPDIR

#Run job with 32 processes by compute node (2 processes by core,
# maximum allowed 4 proc/core) for a total of 2048 processes
runjob --ranks-per-node 16 --np 16 : ./my_code my_args

# Copy output file to submission directory
# Warning: if you need to transfer important volumes
# of data, please use a multi-step job
# $LOADL_STEP_INITDIR is the submission directory
cp data.out $LOADL_STEP_INITDIR/
