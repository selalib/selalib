# Voici un exemple de soumission d'un job pour exécuter
# un code sur 16 cœurs avec 8 processus MPI et 
# 4 threads OpenMP par processus MPI, soit un total de 
# 32 threads avec 2 threads par cœur physique de la
# machine et 8 processus MPI par nœud de calcul. 
# La soumission se fait, en supposant que le fichier #
# de soumission s'appelle job.ll, via la commande :
# llsubmit job.ll


# @ job_name = selalib
# @ job_type = BLUEGENE
# Fichier sortie standard du travail
# @ output = $(job_name).$(jobid)
# Fichier erreur standard du travail
# @ error = $(output)
# Temps elapsed maximum demande
# @ wall_clock_limit = 1:00:00
# Taille bloc d'execution
# en nombre de noeuds de calcul (16 coeurs/noeud)
# @ bg_size = 64
# @ queue
  
runjob --ranks-per-node 8 --envs "OMP_NUM_THREADS=4" --np 8 : ${workdir}/build/bin/test_2d_vp_cartesian ${home}/selalib/prototype/src/simulation/vpsim2d_cartesian_keen
