Modular library for the gyrokinetic simulation model by a semi-Lagrangian method.

As part of french action Fusion, Calvi INRIA projetct developed in collaboration with CEA Cadarache GYSELA simulation code for gyrokinetic simulation of plasma turbulence in Tokamaks. Development and testing of new numerical methods is generally done on different simplified models before its integration into GYSELA. No specification is implemented for this aim, which involves several rewriting code lines before the final integration, which is cumbersome and inefficient. SeLaLib is an API for components of basic numerical implementation also in GYSELA, a framework for parallel single streamline built in order to improve reliability and simplify the work of development.

Sign in to the INRIA GForge and Subscribe to the project selalib by sending a email
to navaro@math.unistra.fr

The selalib-build script runs 

     1. architecture reservation, OS deployement, OS configuration on pipol platform
     2. git pull of the sources from gforge 
     3. compilation
     4. test and visualisation on a dashboard
     5. packages generation
     6. packages upload on the gforge

selalib compilation, testing and packaging
-----------

mkdir build
cd build
cmake <the path of this directory>
make Experimental
(the test result goes to http://irma-webhpc.u-strasbg.fr/index.php?project=SeLaLib)
To display results from outside of IRMA network you need to connect your VPN client.

For package generation
cpack -G <the cpack generator for your system : DEB, RPM, ... >
