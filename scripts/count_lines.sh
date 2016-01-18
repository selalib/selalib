#/bin/bash

#Please put your name when you are the original author of the code present in these files.
#Optimisation, documentation, refactoring don't make you the owner.
#Several names could be set to one file.
#I already put some names but i can make mistake or forget something.
#Pierre
#
#ECG : Edwin Chacon-Golcher
#AB  : Aurore Back
#PN  : Pierre Navaro
#PH  : Philippe Helluy
#AR  : Ahmed Ratnani
#ES  : Eric Sonnendrucker
#CS  : Christophe Steiner
#YG  : Yaman Guclu
#KK  : Katharina Kormann
#KR  : Klaus Reuter
#MCP : Martin Campos-Pinto
#ALH : Antoine Le Hyaric
#LM  : Laura Mendoza
#MM  : Michel Mehrenberger
#NC  : Nicolas Crouseilles
#AH  : Adnane Hamiaz
#SH  : Sever Hirstoaga
#CP  : Charles Prouveur

AB=$(git ls-files -z ${1} add_ons/multipatch/ | xargs -0 cat | wc -l)
echo 'AB='$((AB))

KK1=$(git ls-files -z ${1} add_ons/sparse_grid | xargs -0 cat | wc -l)
KK2=$(git ls-files -z ${1} interfaces/fft | xargs -0 cat | wc -l)
KK3=$(git ls-files -z ${1} particle_methods/pic_basic | xargs -0 cat | wc -l)
KK4=$(git ls-files -z ${1} time_integration/splitting_methods | xargs -0 cat | wc -l)
echo 'KK='$((KK1+KK2+KK3+KK4))

ECG00=$(git ls-files -z ${1} data_structures/boundary_condition_descriptors | xargs -0 cat | wc -l)
ECG01=$(git ls-files -z ${1} data_structures/distribution_function | xargs -0 cat | wc -l)
ECG02=$(git ls-files -z ${1} data_structures/fields | xargs -0 cat | wc -l)
ECG03=$(git ls-files -z ${1} initialization_tools/function_initialization | xargs -0 cat | wc -l)
ECG04=$(git ls-files -z ${1} interpolation/interpolators | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} interpolation/lagrange_interpolation | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} field_solvers/poisson_solvers_parallel | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} linear_solvers | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} parallelization/collective | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} low_level_utilities/assert | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} low_level_utilities/constants | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} low_level_utilities/memory | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} low_level_utilities/timer | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} low_level_utilities/utilities | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} low_level_utilities/working_precision | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} mesh/coordinate_transformations | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} mesh/mesh_calculus | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} parallelization/parallel_utilities | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} parallelization/point_to_point_communications | xargs -0 cat | wc -l)
ECG00=$(git ls-files -z ${1} parallelization/remap | xargs -0 cat | wc -l)
echo 'ECG'$ECG

YG=0
YG+=$(git ls-files -z ${1} data_structures/vector_space | xargs -0 cat | wc -l)
YG+=$(git ls-files -z ${1} time_integration/ode_integrators | xargs -0 cat | wc -l)
YG+=$(git ls-files -z ${1} io/xdmf_io | xargs -0 cat | wc -l)
YG+=$(git ls-files -z ${1} io/xdmf_io_parallel | xargs -0 cat | wc -l)
YG+=$(git ls-files -z ${1} low_level_utilities/errors | xargs -0 cat | wc -l)
YG+=$(git ls-files -z ${1} parallelization/parallel_array_utilities | xargs -0 cat | wc -l)
YG+=$(git ls-files -z ${1} time_integration/ode_solvers | xargs -0 cat | wc -l)
echo 'YG='$YG

PH=$(git ls-files -z ${1} field_solvers/general_coordinate_elliptic_solvers | xargs -0 cat | wc -l)
echo 'PH='$PH

CS=$(git ls-files -z ${1} field_solvers/gyroaverage | xargs -0 cat | wc -l)
echo 'CS='$CS

PN=$(git ls-files -z ${1} field_solvers/maxwell_solvers | xargs -0 cat | wc -l)
PN=$(git ls-files -z ${1} field_solvers/maxwell_solvers_parallel | xargs -0 cat | wc -l)
PN=$(git ls-files -z ${1} io/file_io | xargs -0 cat | wc -l)
PN=$(git ls-files -z ${1} io/file_io_parallel | xargs -0 cat | wc -l)
PN=$(git ls-files -z ${1} particle_methods/pic_visu | xargs -0 cat | wc -l)
echo 'PN='$PN

MM=$(git ls-files -z ${1} interpolation/hermite_interpolation | xargs -0 cat | wc -l)
MM=$(git ls-files -z ${1} field_solvers/quasi_neutral_solvers_parallel | xargs -0 cat | wc -l)
MM=$(git ls-files -z ${1} field_solvers/quasi_neutral_solvers | xargs -0 cat | wc -l)
MM=$(git ls-files -z ${1} initialization_tools/namelist_initialization | xargs -0 cat | wc -l)
MM=$(git ls-files -z ${1} interpolation/periodic_interpolation | xargs -0 cat | wc -l)
MM=$(git ls-files -z ${1} parallelization/reduction | xargs -0 cat | wc -l)
MM=$(git ls-files -z ${1} semi_lagrangian/advection | xargs -0 cat | wc -l)
MM=$(git ls-files -z ${1} semi_lagrangian/fcisl | xargs -0 cat | wc -l)
MM=$(git ls-files -z ${1} time_integration/characteristics | xargs -0 cat | wc -l)
echo 'MM='$MM

AR=$(git ls-files -z ${1} interfaces/sparse_matrix_manager | xargs -0 cat | wc -l)
echo 'AR='$AR

LM=$(git ls-files -z ${1} mesh/meshes | xargs -0 cat | wc -l)
LM=$(git ls-files -z ${1} quadrature | xargs -0 cat | wc -l)
echo 'LM='$LM

X=$(git ls-files -z ${1} field_solvers/poisson_solvers | xargs -0 cat | wc -l)
PN+=PN/3
AH+=X/3
MM+=X/3
exit

echo 'KR'
KR=$(git ls-files -z ${1} parallelization/decomposition | xargs -0 cat | wc -l)

echo 'JA'
JA=$(git ls-files -z ${1} particle_methods/pic_1d | xargs -0 cat | wc -l)
JA=$(git ls-files -z ${1} particle_methods/pif | xargs -0 cat | wc -l)

echo 'SH'
SH=$(git ls-files -z ${1} particle_methods/pic_2d_standard | xargs -0 cat | wc -l)
SH=$(git ls-files -z ${1} particle_methods/random_deviate_generators | xargs -0 cat | wc -l)

echo 'MCP+ALH' 
MCP=$(git ls-files -z ${1} particle_methods/pic_remapped | xargs -0 cat | wc -l)

echo 'ECG+LM+MM+ES'
ES=$(git ls-files -z ${1} splines | xargs -0 cat | wc -l)
