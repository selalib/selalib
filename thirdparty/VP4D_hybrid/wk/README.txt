#-----------------------------------------------------------
# HOW TO RUN VP4D ?
#-----------------------------------------------------------
llrun -f ~/interactive.sh
export SLL_SIMU_DIR=/gpfshome/mds/staff/vgrandgirard/SELALIB-simulations_git/VP4D_hybrid
mpirun -np 4 ${SLL_SIMU_DIR}/build/bin/VP4D_hybrid.exe < ${SLL_SIMU_DIR}/wk/sim4d_qns_general_input.txt > VP4D_res.out
