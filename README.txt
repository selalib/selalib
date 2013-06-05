Modular library for the gyrokinetic simulation model by a semi-Lagrangian method.

As part of french action Fusion, Calvi INRIA projetct developed in collaboration with CEA Cadarache GYSELA simulation code for gyrokinetic simulation of plasma turbulence in Tokamaks. Development and testing of new numerical methods is generally done on different simplified models before its integration into GYSELA. No specification is implemented for this aim, which involves several rewriting code lines before the final integration, which is cumbersome and inefficient. SeLaLib is an API for components of basic numerical implementation also in GYSELA, a framework for parallel single streamline built in order to improve reliability and simplify the work of development.

Sign in to the INRIA GForge and Subscribe to the project selalib by sending a email
to navaro@math.unistra.fr

To develop in Selalib please read :
   - GitQuickstart.txt
   - CMakeQuickstart.txt


The selalib-build script runs 
------------------------------------------

     1. architecture reservation, OS deployement, OS configuration on pipol platform
     2. git pull of the sources from gforge 
     3. compilation
     4. test and visualisation on a dashboard
     5. packages generation
     6. packages upload on the gforge

selalib compilation, installation
------------------------------------------

mkdir build 
cd build
cmake -DCMAKE_BUILD_TYPE=Release                                    \
      -DSLL_PACKAGE=1                                               \
      -DCMAKE_INSTALL_PREFIX=<path where selalib will be installed> \
      <the path of this directory>/prototype/src 
make install

selalib compilation, testing and packaging
------------------------------------------

mkdir build
cd build
cmake <the path of this directory>
make Experimental
(the test result goes to http://cdash.inria.fr/CDash/index.php?project=Selalib)
To display results from outside of IRMA network you need to connect your VPN client.

For package generation
cpack -G <the cpack generator for your system : DEB, RPM, ... >

AUTHORS

 - Aurore Back
 - Morgane Bergot
 - Raphael Blanchard
 - Edwin Chacon-Golcher
 - Samuel Desantis
 - Aliou Diouf
 - Emmanuel Frenod
 - Adnane Hamiaz
 - Philippe Helluy
 - Sever Hirstoaga
 - Eric Madaule
 - Michel Mehrenberger
 - Pierre Navaro
 - Laurent Navoret
 - Eric Sonnendrucker
