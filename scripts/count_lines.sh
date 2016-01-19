#/bin/bash

#Please put your name when you are the original author of the code present in these files.
#Optimisation, documentation, refactoring don't make you the owner.
#Several names could be set to one file.
#I already put some names but i can make mistake or forget something.
#Pierre
#
#EC : Edwin Chacon-Golcher
#AB : Aurore Back
#PN : Pierre Navaro
#PH : Philippe Helluy
#AR : Ahmed Ratnani
#ES : Eric Sonnendrucker
#CS : Christophe Steiner
#YG : Yaman Guclu
#KK : Katharina Kormann
#KR : Klaus Reuter
#MC : Martin Campos-Pinto
#AL : Antoine Le Hyaric
#LM : Laura Mendoza
#MM : Michel Mehrenberger
#NC : Nicolas Crouseilles
#AH : Adnane Hamiaz
#SH : Sever Hirstoaga
#CP : Charles Prouveur

AB=$(git ls-files ${1} add_ons/multipatch/ | grep F90| xargs cat | wc -l)
AB=$((AB))

KK1=$(git ls-files add_ons/sparse_grid | grep F90| xargs cat | wc -l)
KK2=$(git ls-files interfaces/fft | grep F90| xargs cat | wc -l)
KK3=$(git ls-files particle_methods/pic_basic | grep F90| xargs cat | wc -l)
KK4=$(git ls-files time_integration/splitting_methods | grep F90| xargs cat | wc -l)

KK=$((KK1+KK2+KK3+KK4))

EC00=$(git ls-files data_structures/boundary_condition_descriptors | grep F90| xargs cat | wc -l)
EC01=$(git ls-files data_structures/distribution_function | grep F90| xargs cat | wc -l)
EC02=$(git ls-files data_structures/fields | grep F90| xargs cat | wc -l)
EC03=$(git ls-files initialization_tools/function_initialization | grep F90| xargs cat | wc -l)
EC04=$(git ls-files interpolation/interpolators | grep F90| xargs cat | wc -l)
EC05=$(git ls-files interpolation/lagrange_interpolation | grep F90| xargs cat | wc -l)
EC06=$(git ls-files field_solvers/poisson_solvers_parallel | grep F90| xargs cat | wc -l)
EC07=$(git ls-files linear_solvers | grep F90| xargs cat | wc -l)
EC08=$(git ls-files parallelization/collective | grep F90| xargs cat | wc -l)
EC09=$(git ls-files low_level_utilities/assert | grep F90| xargs cat | wc -l)
EC10=$(git ls-files low_level_utilities/constants | grep F90| xargs cat | wc -l)
EC11=$(git ls-files low_level_utilities/memory | grep F90| xargs cat | wc -l)
EC12=$(git ls-files low_level_utilities/timer | grep F90| xargs cat | wc -l)
EC13=$(git ls-files low_level_utilities/utilities | grep F90| xargs cat | wc -l)
EC14=$(git ls-files low_level_utilities/working_precision | grep F90| xargs cat | wc -l)
EC15=$(git ls-files mesh/coordinate_transformations | grep F90| xargs cat | wc -l)
EC16=$(git ls-files mesh/mesh_calculus | grep F90| xargs cat | wc -l)
EC17=$(git ls-files parallelization/parallel_utilities | grep F90| xargs cat | wc -l)
EC18=$(git ls-files parallelization/point_to_point_communications | grep F90| xargs cat | wc -l)
EC19=$(git ls-files parallelization/remap | grep F90| xargs cat | wc -l)

EC=$((EC00+EC01+EC02+EC03+EC04+EC05+EC06+EC07+EC08+EC09+EC10+EC11+EC12+EC13+EC14+EC15+EC16+EC17+EC18+EC19))

YG0=$(git ls-files data_structures/vector_space | grep F90| xargs cat | wc -l)
YG1=$(git ls-files time_integration/ode_integrators | grep F90| xargs cat | wc -l)
YG2=$(git ls-files io/xdmf_io | grep F90| xargs cat | wc -l)
YG3=$(git ls-files io/xdmf_io_parallel | grep F90| xargs cat | wc -l)
YG4=$(git ls-files low_level_utilities/errors | grep F90| xargs cat | wc -l)
YG5=$(git ls-files parallelization/parallel_array_utilities | grep F90| xargs cat | wc -l)
YG6=$(git ls-files time_integration/ode_solvers | grep F90| xargs cat | wc -l)

YG=$((YG0+YG1+YG2+YG3+YG4+YG5+YG6))

PH=$(git ls-files field_solvers/general_coordinate_elliptic_solvers | grep F90| xargs cat | wc -l)

PH=$((PH))

CS=$(git ls-files field_solvers/gyroaverage | grep F90| xargs cat | wc -l)

CS=$((CS))

PN0=$(git ls-files field_solvers/maxwell_solvers | grep F90| xargs cat | wc -l)
PN1=$(git ls-files field_solvers/maxwell_solvers_parallel | grep F90| xargs cat | wc -l)
PN2=$(git ls-files io/file_io | grep F90| xargs cat | wc -l)
PN3=$(git ls-files io/file_io_parallel | grep F90| xargs cat | wc -l)
PN4=$(git ls-files particle_methods/pic_visu | grep F90| xargs cat | wc -l)

PN=$((PN0+PN1+PN2+PN3+PN4))

MM0=$(git ls-files interpolation/hermite_interpolation | grep F90| xargs cat | wc -l)
MM1=$(git ls-files field_solvers/quasi_neutral_solvers_parallel | grep F90| xargs cat | wc -l)
MM2=$(git ls-files field_solvers/quasi_neutral_solvers | grep F90| xargs cat | wc -l)
MM3=$(git ls-files initialization_tools/namelist_initialization | grep F90| xargs cat | wc -l)
MM4=$(git ls-files interpolation/periodic_interpolation | grep F90| xargs cat | wc -l)
MM5=$(git ls-files parallelization/reduction | grep F90| xargs cat | wc -l)
MM6=$(git ls-files semi_lagrangian/advection | grep F90| xargs cat | wc -l)
MM7=$(git ls-files semi_lagrangian/fcisl | grep F90| xargs cat | wc -l)
MM8=$(git ls-files time_integration/characteristics | grep F90| xargs cat | wc -l)

MM=$((MM0+MM1+MM2+MM3+MM4+MM5+MM6+MM7+MM8))

AR=$(git ls-files interfaces/sparse_matrix_manager | grep F90| xargs cat | wc -l)

LM0=$(git ls-files mesh/meshes | grep F90| xargs cat | wc -l)
LM1=$(git ls-files quadrature | grep F90| xargs cat | wc -l)
LM=$((LM0+LM1))

X=$(git ls-files field_solvers/poisson_solvers | grep F90| xargs cat | wc -l)
PN=$((PN+X/3))
AH=$((AH+X/3))
MM=$((MM+X/3))

KR=$(git ls-files parallelization/decomposition | grep F90| xargs cat | wc -l)

JA0=$(git ls-files particle_methods/pic_1d | grep F90| xargs cat | wc -l)
JA1=$(git ls-files particle_methods/pif | grep F90| xargs cat | wc -l)
JA=$((JA0+JA1))

SH0=$(git ls-files particle_methods/pic_2d_standard | grep F90| xargs cat | wc -l)
SH1=$(git ls-files particle_methods/random_deviate_generators | grep F90| xargs cat | wc -l)
SH=$((SH0+SH1))

X=$(git ls-files particle_methods/pic_remapped | grep F90| xargs cat | wc -l)
MC=$((X/2))
AL=$((X/2))

ES=$(git ls-files splines | grep F90| xargs cat | wc -l)

SL=$((AB+AR+KK+EC+YG+CS+PN+MM+AR+LM+KR+H+JA+MC+ES+AB+AR))
echo 'SL='$((SL))
echo 'SL='$(git ls-files . | grep F90| xargs cat | wc -l)

echo 'EC='$((EC*1000/SL))
echo 'YG='$((YG*1000/SL))
echo 'CS='$((CS*1000/SL))
echo 'PN='$((PN*1000/SL))
echo 'MM='$((MM*1000/SL))
echo 'AR='$((AR*1000/SL))
echo 'LM='$((LM*1000/SL))
echo 'KR='$((KR*1000/SL))
echo 'SH='$((SH*1000/SL))
echo 'JA='$((JA*1000/SL))
echo 'MC='$((MC*1000/SL))
echo 'AL='$((AL*1000/SL))
echo 'ES='$((ES*1000/SL))
echo 'AB='$((AB*1000/SL))
echo 'AR='$((AR*1000/SL))
echo 'KK='$((KK*1000/SL))
